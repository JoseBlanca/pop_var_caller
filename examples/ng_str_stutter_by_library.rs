//! **Where does a tract start to stutter, and how much does the answer move between libraries?** —
//! the survey that sets ng's per-period copy floors.
//!
//! `doc/devel/ng/spec/parameter_prepass_ssr.md` §5.1 defines an STR locus as one that is **likely to
//! stutter**, not merely one that contains a repeat, and puts the copy floor where a period's tracts
//! start behaving that way. Two things follow, and this program measures both:
//!
//! 1. **Where each period crosses**, per (period, reference repeat count) — which §5's three pooled
//!    repeat-count bands cannot give, since they are pooled across periods.
//! 2. **How far that crossing moves between libraries**, which is the axis that actually drives
//!    stutter. PCR amplification stutters more than a PCR-free preparation, so the floor a library
//!    deserves is a property of the library. Every number ng has today comes from a **single**
//!    library per species, so the one axis that matters is the one nothing has varied.
//!
//! ## The two criteria, and why a period needs both
//!
//! - **`off_ref_share`** — reads differing from the reference tract length, over all reads. *Does it
//!   move at all.* An over-estimate of stutter, because a locus's own alleles differ from the
//!   reference too, and the only criterion available for **period 1**.
//! - **`guard_share`** — of the reads that differ, the share differing by something other than a
//!   whole number of copies. *Is the movement the kind this noise model can express.* The sharper
//!   test, because a stratum failing it produces a confident slippage rate that is mostly ordinary
//!   indel. **Identically zero at period 1** — every integer is a multiple of one — so it locates a
//!   period's crossing on its own curve and its level is not comparable between periods.
//!
//! ## Read length is a confound, and this program cannot see it — join to the manifest
//!
//! A 100 bp library cannot span a tract a 150 bp library spans, so its long-repeat strata are empty
//! for purely geometric reasons. `benchmarks/ssr_tomato1/scripts/rick_sample_manifest.sh` says the
//! same thing and refuses to mix lengths without stratifying: *"the complete-vs-partial frontier IS
//! a read-length effect … you cannot ignore it"*.
//!
//! **This program does not have the read length.** A locus observation carries the tract the read
//! showed, not the read — `SequenceObservation::bases` is "allele content, in read coordinates" —
//! so the per-group `mean_tract_bases_per_read` below is the mean *observed tract*, which is a
//! property of the loci walked and not of the library's chemistry. **Join the output to that
//! manifest on the file or the `@RG` id to stratify by read length**; do not read this column as a
//! substitute for it.
//!
//! The **floors** question is the least exposed part of the confound, since the low-repeat strata a
//! floor sits in are spanned by every read length. It is the long-repeat end that a mixed survey
//! would misread.
//!
//! ## Without `--min-copies` the answer is censored at the floors it is trying to measure
//!
//! Region typing applies `MinCopies` **while it types**, so a walk at the defaults never emits a
//! tract below them and the survey's curves start exactly at `[6, 4, 4, 3, 3, 3]`. That can say
//! whether a floor should be **raised**; it cannot say whether one could be lowered, and it cannot
//! see the shape of the curve on the side the floor is meant to cut off. **Pass `--min-copies 2`
//! for a survey meant to place floors** — ng made `MinCopies` a knob for exactly this question,
//! where production hardcodes it in a `const fn`.
//!
//! ```text
//! ng_str_stutter_by_library [--contigs a,b] [--regions r.bed] [--min-copies 2] \
//!     <reference.fa> <sample.cram> [more...]
//! ```
//!
//! **Restrict the walk.** The cost is one region-typing pass plus a read fetch per sample, so it
//! scales with samples; a few hundred thousand loci settle these curves, which is one tomato
//! chromosome. `--contigs SL4.0ch01` is the recommended shape for a wide survey.
//!
//! Output: the read-group table as `#rg` lines, the per-stratum TSV, and — per read group — the
//! copy floor each period's own data implies.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, ReadWitness, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::{NameOrigin, build_read_groups};
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::MinCopies;
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegionConfig, TypedRegionIterator,
};
use pop_var_caller::ng::types::{Bp, ContigId, ReadGroupId};
use pop_var_caller::regions::ContigBounds;

use std::collections::BTreeMap;

/// Offsets are recorded over `±OFFSET_HALF_RANGE`, the ends saturating (arch §2.1).
const OFFSET_HALF_RANGE: i64 = 4;

/// Above this share of the reads that differ from the reference, the stratum is one the STR noise
/// model does not describe (spec §5, §5.1).
const GUARD_SHARE_LIMIT: f64 = 0.10;

/// How often a read must differ from the reference length before a period-1 tract is worth the STR
/// path. **Soft, and the one number here that is chosen rather than derived** — the guard share is
/// vacuous at period 1, so monos need a level criterion, and this is where the survey's own curve
/// should be read to place it.
const MONO_OFF_REFERENCE_LIMIT: f64 = 0.01;

/// The fewest reads a stratum needs before its shares are worth reading. A floor against noise, not
/// a precision target.
const MIN_READS_TO_JUDGE: u64 = 500;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Stratum {
    period: u8,
    repeats: u32,
}

#[derive(Default, Clone)]
struct StratumCounts {
    loci: u64,
    reads: u64,
    reads_off_reference: u64,
    reads_not_whole_repeat: u64,
    reads_in_end_buckets: u64,
}

impl StratumCounts {
    fn off_reference_share(&self) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            self.reads_off_reference as f64 / self.reads as f64
        }
    }
    fn guard_share(&self) -> f64 {
        if self.reads_off_reference == 0 {
            0.0
        } else {
            self.reads_not_whole_repeat as f64 / self.reads_off_reference as f64
        }
    }
}

/// What one read group contributed, beside its strata — the fields an analysis needs to tell
/// geometry from chemistry.
#[derive(Default, Clone)]
struct GroupCounts {
    loci: u64,
    reads: u64,
    read_bases: u64,
}

/// The stratum a locus belongs to, from the **reference** tract alone.
fn stratum_of(locus: &SampleLocusObservations) -> Option<Stratum> {
    let LocusKind::Ssr(detail) = &locus.kind else {
        return None;
    };
    let period = detail.motif.period();
    if period == 0 || !locus.reference_bases.len().is_multiple_of(period) {
        return None;
    }
    Some(Stratum {
        period: period as u8,
        repeats: (locus.reference_bases.len() / period) as u32,
    })
}

/// Fold one locus into the per-(read group × stratum) counts.
fn observe_locus(
    strata: &mut BTreeMap<(ReadGroupId, Stratum), StratumCounts>,
    groups: &mut BTreeMap<ReadGroupId, GroupCounts>,
    locus: &SampleLocusObservations,
) {
    let Some(stratum) = stratum_of(locus) else {
        return;
    };
    let period = i64::from(stratum.period);
    let reference_len = locus.reference_bases.len() as i64;

    let mut seen_groups: BTreeMap<ReadGroupId, ()> = BTreeMap::new();
    for obs in &locus.observations {
        if obs.read_witness != ReadWitness::Complete {
            continue;
        }
        let reads = u64::from(obs.num_obs);
        let entry = strata.entry((obs.read_group, stratum)).or_default();
        entry.reads += reads;
        let group = groups.entry(obs.read_group).or_default();
        group.reads += reads;
        group.read_bases += reads * obs.bases.len() as u64;
        seen_groups.insert(obs.read_group, ());

        let difference = obs.bases.len() as i64 - reference_len;
        if difference == 0 {
            continue;
        }
        entry.reads_off_reference += reads;
        if !difference.rem_euclid(period).eq(&0) {
            entry.reads_not_whole_repeat += reads;
            continue;
        }
        let offset = difference / period;
        if offset.abs() >= OFFSET_HALF_RANGE {
            entry.reads_in_end_buckets += reads;
        }
    }
    for (group, ()) in seen_groups {
        strata.entry((group, stratum)).or_default().loci += 1;
        groups.entry(group).or_default().loci += 1;
    }
}

/// **The copy floor a period's own data implies**: scanning upward, the first repeat count that
/// passes and whose neighbour above also passes.
///
/// **Two consecutive passes, and neither a first-crossing nor a last-failure rule.** A
/// first-crossing rule picks a single dip — tomato's dinucleotides pass at 5 and fail again at 6.
/// A last-failure rule ("one above the highest failing stratum") is worse and was this function's
/// first version: the high-repeat strata are thin, so one noisy 46-locus stratum at 9 copies pushed
/// the trinucleotide floor to 11 when the curve settles at 7. The question is where the **low** end
/// stops misbehaving, so the scan runs from the bottom and asks for the behaviour to persist.
fn implied_floor(
    strata: &BTreeMap<(ReadGroupId, Stratum), StratumCounts>,
    group: ReadGroupId,
    period: u8,
) -> Option<u32> {
    let rows: Vec<(u32, &StratumCounts)> = strata
        .iter()
        .filter(|((g, s), c)| *g == group && s.period == period && c.reads >= MIN_READS_TO_JUDGE)
        .map(|((_, s), c)| (s.repeats, c))
        .collect();
    if rows.is_empty() {
        return None;
    }
    let passes = |counts: &StratumCounts| -> bool {
        if period == 1 {
            // The guard share is vacuous at period 1 — every integer is a multiple of one — so the
            // level criterion is the only one available.
            counts.off_reference_share() >= MONO_OFF_REFERENCE_LIMIT
        } else {
            counts.guard_share() <= GUARD_SHARE_LIMIT
        }
    };
    for (index, (repeats, counts)) in rows.iter().enumerate() {
        if !passes(counts) {
            continue;
        }
        let neighbour_holds = match rows.get(index + 1) {
            Some((next, next_counts)) if *next == repeats + 1 => passes(next_counts),
            // Nothing judgeable above it: the curve has run out of loci, not out of good behaviour.
            _ => true,
        };
        if neighbour_holds {
            return Some(*repeats);
        }
    }
    // Every judgeable stratum misbehaves — the floor is above what this walk can see.
    rows.last().map(|(repeats, _)| repeats + 1)
}

fn origin_label(origin: NameOrigin) -> &'static str {
    match origin {
        NameOrigin::Declared => "declared",
        NameOrigin::Synthesized => "synthesized",
    }
}

fn run(
    fasta: &Path,
    alignments: &[PathBuf],
    contig_filter: &[String],
    regions_bed: Option<&Path>,
    min_copies: Option<u32>,
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs: ContigList = info.contig_list();
    // One reference for the whole cohort, and so one copy of its bases — the per-file repository
    // this replaces cost ~752 MiB of resident tomato genome per open CRAM (ng_ssr_cohort_stutter).
    let reference = OpenReference::new(info);

    let read_groups = build_read_groups(alignments)?;
    let samples: Vec<SampleReads> = read_groups
        .read_groups_per_sample()
        .iter()
        .map(|entry| {
            SampleReads::open(
                entry,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                true,
            )
        })
        .collect::<Result<_, _>>()?;
    eprintln!(
        "  {} sample(s), {} read group(s)",
        samples.len(),
        read_groups.iter().count()
    );

    // **Lowering the floor is what lets the survey see below it.** At the defaults the walk emits
    // nothing under `[6, 4, 4, 3, 3, 3]`, so the curves start where the answer is meant to be
    // decided and the measurement is censored at exactly the wrong place.
    let mut walk_config = TypedRegionConfig::default();
    if let Some(copies) = min_copies {
        walk_config.criteria.min_copies = MinCopies::uniform(copies);
        eprintln!("  typing from {copies} copies at every period, below ng's defaults");
    }
    let bundle_threshold = Bp(walk_config.criteria.bundle_threshold);
    #[expect(
        clippy::arc_with_non_send_sync,
        reason = "RawRefSeq is implemented for Arc only; this walk is single-threaded"
    )]
    let shared_reference = Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));
    // One generator per sample: a generator holds a reader positioned for one sample's files, so a
    // shared one would answer every sample out of the first one's reads.
    let mut generators = samples
        .iter()
        .map(|_| {
            SsrGenerator::with_default_aligner(
                Arc::clone(&shared_reference),
                {
                    let reference = Arc::clone(&shared_reference);
                    move || Arc::clone(&reference)
                },
                SsrGeneratorConfig::default(),
                bundle_threshold,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut strata: BTreeMap<(ReadGroupId, Stratum), StratumCounts> = BTreeMap::new();
    let mut groups: BTreeMap<ReadGroupId, GroupCounts> = BTreeMap::new();

    let bed_spans = regions_bed
        .map(|bed| {
            let bounds: Vec<ContigBounds<'_>> = contigs
                .entries
                .iter()
                .map(|e| ContigBounds {
                    name: &e.name,
                    length: e.length as u32,
                })
                .collect();
            GenomeRegions::from_bed_path(bed, &bounds)
        })
        .transpose()?;

    let wanted_contig = |contig: ContigId| {
        contig_filter.is_empty()
            || contigs
                .entries
                .get(contig.0 as usize)
                .is_some_and(|e| contig_filter.iter().any(|n| n == &e.name))
    };

    let mut walks: Vec<(String, TypedRegionIterator<WindowedRefSeq>)> = Vec::new();
    match bed_spans {
        Some(spans) => {
            let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
            walks.push((
                "BED".to_string(),
                TypedRegionIterator::over_regions(walk_reference, spans, walk_config.clone())?,
            ));
        }
        None => {
            for (index, entry) in contigs.entries.iter().enumerate() {
                if !wanted_contig(ContigId(index as u32)) {
                    continue;
                }
                let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
                walks.push((
                    entry.name.clone(),
                    TypedRegionIterator::over_contig(
                        walk_reference,
                        ContigId(index as u32),
                        walk_config.clone(),
                    )?,
                ));
            }
        }
    }

    for (label, mut walk) in walks {
        eprintln!("  walking {label}");
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            if !wanted_contig(region.region.contig) {
                continue;
            }
            for (sample, generator) in samples.iter().zip(generators.iter_mut()) {
                generator.begin_segment(region.region);
                while let Some(locus) = generator.next_locus(segment, sample)? {
                    observe_locus(&mut strata, &mut groups, &locus);
                }
            }
        }
    }

    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());

    // The read-group table, once. **`mean_tract_bases_per_read` is not a read length** — a locus
    // observation carries the tract the read showed, not the read — so stratifying by read length
    // means joining this to `rick_sample_manifest.sh`'s output on the file or the `@RG` id.
    for (id, group) in read_groups.iter() {
        let counts = groups.get(&id).cloned().unwrap_or_default();
        writeln!(
            out,
            "#rg\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            id.get(),
            group.id,
            group.sample,
            group.library.value,
            origin_label(group.library.origin),
            group.experiment.value,
            origin_label(group.experiment.origin),
            group.platform.as_deref().unwrap_or(""),
            counts.loci,
            counts.reads,
            if counts.reads > 0 {
                counts.read_bases as f64 / counts.reads as f64
            } else {
                0.0
            },
        )?;
    }
    writeln!(
        out,
        "#rg_columns\tread_group\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\t\
         experiment_origin\tplatform\tloci\treads\tmean_tract_bases_per_read"
    )?;

    // The per-read-group floors, which are the answer this survey exists for.
    writeln!(
        out,
        "#floor_columns\tread_group\tperiod\timplied_floor\tcriterion"
    )?;
    for (id, _) in read_groups.iter() {
        for period in 1..=6u8 {
            if let Some(floor) = implied_floor(&strata, id, period) {
                writeln!(
                    out,
                    "#floor\t{}\t{}\t{}\t{}",
                    id.get(),
                    period,
                    floor,
                    if period == 1 {
                        "off_reference_share"
                    } else {
                        "guard_share"
                    }
                )?;
            }
        }
    }

    writeln!(
        out,
        "read_group\tperiod\trepeats\tloci\treads\toff_ref_reads\toff_ref_share\t\
         not_whole_reads\tguard_share\tend_bucket_reads"
    )?;
    for ((group, stratum), counts) in &strata {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{:.5}\t{}",
            group.get(),
            stratum.period,
            stratum.repeats,
            counts.loci,
            counts.reads,
            counts.reads_off_reference,
            counts.off_reference_share(),
            counts.reads_not_whole_repeat,
            counts.guard_share(),
            counts.reads_in_end_buckets,
        )?;
    }
    out.flush()?;

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut contig_filter: Vec<String> = Vec::new();
    let mut regions_bed: Option<PathBuf> = None;
    let mut min_copies: Option<u32> = None;
    let mut rest = std::env::args().skip(1);
    while let Some(arg) = rest.next() {
        match arg.as_str() {
            "--regions" => match rest.next() {
                Some(path) => regions_bed = Some(PathBuf::from(path)),
                None => {
                    eprintln!("error: --regions needs a BED path");
                    return ExitCode::from(2);
                }
            },
            "--contigs" => match rest.next() {
                Some(list) => {
                    contig_filter = list.split(',').map(|s| s.trim().to_string()).collect()
                }
                None => {
                    eprintln!("error: --contigs needs a comma-separated list");
                    return ExitCode::from(2);
                }
            },
            "--min-copies" => match rest.next().and_then(|v| v.parse().ok()) {
                Some(copies) => min_copies = Some(copies),
                None => {
                    eprintln!("error: --min-copies needs a whole number of copies");
                    return ExitCode::from(2);
                }
            },
            _ => positional.push(arg),
        }
    }
    if positional.len() < 2 {
        eprintln!(
            "usage: ng_str_stutter_by_library [--contigs a,b] [--regions r.bed] <reference.fa> \
             <sample.bam|cram> [sample ...]\n\
             measures, per (read group, motif period, reference repeat count), how often a read \
             differs from the reference tract length and how much of that difference is a whole \
             number of copies — the two curves ng's per-period copy floors are placed on.\n\
             Restrict the walk: the cost scales with samples, and one chromosome settles these \
             curves."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&positional[0]);
    let alignments: Vec<PathBuf> = positional[1..].iter().map(PathBuf::from).collect();

    match run(
        &fasta,
        &alignments,
        &contig_filter,
        regions_bed.as_deref(),
        min_copies,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
