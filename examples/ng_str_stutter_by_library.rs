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
//! ## Lowering `--min-copies` does not extend the curve downward — measured, 2026-08-07
//!
//! Region typing applies `MinCopies` **while it types**, so a walk at the defaults never emits a
//! tract below `[6, 4, 4, 3, 3, 3]` and this survey's curves start exactly where the floor is meant
//! to be decided. An earlier version of this file told the reader to pass `--min-copies 2` for that
//! reason, and `scripts/ng_str_library_survey.sh` defaulted to it. **It measures nothing.**
//!
//! The floor is read **twice**, and the second reading is the one nobody accounted for: `prefilter`
//! applies it before bundling (`segment_criteria.rs:601`) and `classify` applies it again
//! (`:985`). So it decides not only *what is admitted as a locus* but *what counts as a
//! neighbouring repeat*. Two copies of a mononucleotide is any `AA`, which occurs every few bases,
//! so at a uniform floor of two every real tract has a neighbour inside the bundle threshold, no
//! tract has clean flanks, and region typing emits `SsrBundle` — which this survey does not count,
//! because a bundle names no locus. Over a 2 Mb slice of tomato SL4.0ch01:
//!
//! | copy floors | STR loci | bundles | bundle bp |
//! |---|---:|---:|---:|
//! | `[6,4,4,3,3,3]` — ng's defaults | 6,237 | 1,943 | 74,289 |
//! | uniform 4 | 7,434 | 11,720 | 675,372 |
//! | uniform 3 | 848 | 7,950 | 1,623,636 |
//! | uniform 2 | **0** | 1 | 1,177,849 |
//! | `[2,4,4,3,3,3]` | **0** | 18 | 1,599,157 |
//! | `[6,2,4,3,3,3]` | 1,618 | 13,913 | 1,017,906 |
//!
//! **And the damage is not confined to the period that moved**, which is what makes a one-period
//! sweep no safer than a uniform one: lowering period 2 alone takes period-1 tracts at six copies
//! from 2,678 loci to 225 — 92% gone at a period whose own floor never changed. Two settings
//! therefore measure two different populations of loci rather than two ends of one curve.
//!
//! ## One step down is usable, and the criterion that places the floors survives it
//!
//! **`--min-copies 5,3,3,3,3,3` is what the archive survey runs**, and it is a step rather than a
//! plunge for the reason above. One step still costs loci — about half of every stratum shared with
//! a default walk — but the two criteria are not equally affected, and the sharper one is the one
//! that holds:
//!
//! - **The guard share is essentially unmoved.** Dinucleotides at 4 repeats read 0.345 in a
//!   default walk and 0.344 in a swept one; at 20 bp of flank, 0.431 against 0.408. So the
//!   criterion that decides periods 2 to 6 can be read off a swept walk directly.
//! - **The off-reference share is not**, coming back ~30% low, because the tracts that survive are
//!   the isolated ones and isolated tracts sit in cleaner context. That is the only criterion
//!   **mononucleotides** have, so the period-1 floor rests on softer evidence than the rest — and
//!   the bias understates stutter, which pushes a floor *up*, the direction that fails silently.
//!
//! **And the answer the survey exists for is robust to all of it.** Dinucleotides at 3 repeats read
//! a guard share of 63% against a one-in-ten threshold, unchanged across every flank width and both
//! copy-floor settings tried, on thousands of loci. Six-fold over is not a marginal call, so no
//! correction reaches it. See `doc/devel/ng/spec/parameter_prepass_ssr.md` §5.1.
//!
//! ```text
//! ng_str_stutter_by_library [--contigs a,b] [--regions r.bed] [--min-copies 5,3,3,3,3,3] \
//!     [--bundle-threshold 20] <reference.fa> <sample.cram> [more...]
//! ```
//!
//! `--bundle-threshold` is the clean sequence a tract needs either side to be nameable as a locus,
//! and it is also the aligner's anchor. **ng's default is 20 bp and this survey does not need to
//! set it**; the flag exists so the knob can be swept, which is how 20 was chosen.
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
use pop_var_caller::ng::region_typing::segment_criteria::{MAX_MOTIF_LEN, MinCopies};
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegionConfig, TypedRegionCounts, TypedRegionIterator,
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

/// Parse `--min-copies`: one number for every period, or a table of six.
///
/// **The table form is what the floors question actually needs.** Lowering every period at once
/// admits two-base homopolymers, which carpet the genome, so every tract acquires a neighbour
/// inside the bundle threshold and no locus survives (see [`EmptySurvey`]). Lowering **one** period
/// and leaving the rest at ng's defaults keeps the walk intact and is how a period's own curve is
/// extended below its floor.
fn parse_min_copies(value: &str) -> Option<MinCopies> {
    let fields: Vec<&str> = value.split(',').map(str::trim).collect();
    match fields.as_slice() {
        [one] => one.parse().ok().map(MinCopies::uniform),
        table if table.len() == MAX_MOTIF_LEN => {
            let mut by_period = [0u32; MAX_MOTIF_LEN];
            for (slot, field) in by_period.iter_mut().zip(table) {
                *slot = field.parse().ok()?;
            }
            // Periods wider than the table keep production's floor: this survey measures 1 to 6.
            Some(MinCopies::new(by_period, 3))
        }
        _ => None,
    }
}

/// Add one walk's typing tally into the run's. Field by field, because the walks are one per contig
/// or one per BED and the survey's question is about the whole run.
fn add_typing_counts(total: &mut TypedRegionCounts, walk: &TypedRegionCounts) {
    total.spans += walk.spans;
    total.ssr_loci += walk.ssr_loci;
    total.ssr_bundles += walk.ssr_bundles;
    total.ssr_bundle_bp += walk.ssr_bundle_bp;
    total.generic += walk.generic;
    total.satellites += walk.satellites;
    total.satellite_bp += walk.satellite_bp;
    total.repeat_bp_with_no_locus += walk.repeat_bp_with_no_locus;
    total.rejected_by_reason.copy_floor += walk.rejected_by_reason.copy_floor;
    total.rejected_by_reason.purity += walk.rejected_by_reason.purity;
    total.rejected_by_reason.compound += walk.rejected_by_reason.compound;
    total.rejected_by_reason.no_clean_trim += walk.rejected_by_reason.no_clean_trim;
    total.rejected_by_reason.flank_clamped += walk.rejected_by_reason.flank_clamped;
}

/// **A walk that typed nothing this survey can measure.** It is an error and not an empty table:
/// the driver that runs this over an archive treats a written result as a finished file, so a
/// silent empty success is skipped on every later pass and the library is never surveyed.
///
/// The message names the mechanism from the typing tally, because the two ways to type nothing look
/// identical from the outside and want opposite responses. Repeats that became **bundles** mean the
/// copy floor was set so low that neighbouring tracts cluster and no locus can be named — lower
/// floors then measure *less*, which is the trap `--min-copies` walks into. No repeats **at all**
/// means the walk covered sequence that has none, which is a region or reference problem.
#[derive(Debug)]
struct EmptySurvey {
    min_copies: Option<MinCopies>,
    typing: TypedRegionCounts,
}

impl std::fmt::Display for EmptySurvey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let counts = &self.typing;
        write!(
            f,
            "the walk produced no stratum rows — nothing was measured. Region typing over \
             {spans} span(s) emitted {loci} STR locus/loci, {bundles} bundle(s) covering \
             {bundle_bp} bp, {satellites} satellite(s) and {generic} generic region(s); \
             {no_locus} bp of repeat yielded no locus (copy floor {floor}, purity {purity}, \
             compound {compound}, no clean trim {trim}, flank clamped {clamped}).",
            spans = counts.spans,
            loci = counts.ssr_loci,
            bundles = counts.ssr_bundles,
            bundle_bp = counts.ssr_bundle_bp,
            satellites = counts.satellites,
            generic = counts.generic,
            no_locus = counts.repeat_bp_with_no_locus,
            floor = counts.rejected_by_reason.copy_floor,
            purity = counts.rejected_by_reason.purity,
            compound = counts.rejected_by_reason.compound,
            trim = counts.rejected_by_reason.no_clean_trim,
            clamped = counts.rejected_by_reason.flank_clamped,
        )?;
        if counts.ssr_loci == 0 && counts.ssr_bundles > 0 {
            write!(
                f,
                " Every repeat became a bundle rather than a locus: at this copy floor \
                 neighbouring tracts sit within the bundle threshold of each other, so no tract \
                 has the clean flanks a locus needs."
            )?;
            if let Some(floors) = self.min_copies {
                let table: Vec<String> = (1..=6u8)
                    .map(|p| floors.for_period(p).to_string())
                    .collect();
                write!(
                    f,
                    " --min-copies was {}: lower one period at a time and leave the rest at ng's \
                     defaults (`--min-copies 6,2,4,3,3,3`), because lowering every period past the \
                     point where tracts start touching measures less, not more.",
                    table.join(",")
                )?;
            }
        } else if counts.ssr_loci > 0 {
            write!(
                f,
                " Loci were typed but no read witnessed one completely — check that the \
                 alignments cover the walked regions."
            )?;
        }
        Ok(())
    }
}

impl std::error::Error for EmptySurvey {}

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
    min_copies: Option<MinCopies>,
    bundle_threshold_override: Option<u64>,
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

    // **Lowering the floor is what lets the survey see below it, and lowering it at every period at
    // once is what stops the survey seeing anything.** At the defaults the walk emits nothing under
    // `[6, 4, 4, 3, 3, 3]`, so the curves start where the answer is meant to be decided. But a floor
    // low enough to admit two-copy tracts admits a two-base homopolymer, which occurs every few
    // bases, so every real tract acquires a neighbour inside the bundle threshold and the whole
    // walk becomes one bundle carrying no locus at all — measured below on tomato.
    let mut walk_config = TypedRegionConfig::default();
    if let Some(floors) = min_copies {
        walk_config.criteria.min_copies = floors;
        eprintln!(
            "  typing from {} copies at periods 1-6",
            (1..=6u8)
                .map(|p| floors.for_period(p).to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
    }
    // **The other lever on the same problem** (§5.1). A tract needs this many clean bases either
    // side before it can be named a locus, so it is the radius in which a neighbouring repeat
    // spoils a tract. Narrowing it keeps more tracts nameable at a lowered copy floor, which is the
    // selection that biases a swept walk. It is not free: `flank_bp` is held at or below this, so
    // it is also the anchor the STR aligner places a read against.
    if let Some(bases) = bundle_threshold_override {
        walk_config.criteria.bundle_threshold = bases;
        eprintln!("  a tract needs {bases} clean bases either side, against ng's default of 30");
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
                {
                    // **The anchor moves with the radius, because the two are one number in
                    // practice.** `flank_bp` is the sequence the STR aligner places a read
                    // against, and it must not reach past the repeat-free radius region typing
                    // guarantees — so narrowing the radius to keep more tracts nameable also
                    // shortens every read's anchor. That is the cost of this lever.
                    SsrGeneratorConfig {
                        flank_bp: bundle_threshold,
                        ..SsrGeneratorConfig::default()
                    }
                },
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

    // **What region typing produced, not only what this survey could use.** The two differ, and
    // when they differ the survey measures nothing while reporting success: only `SsrSegment` names
    // a locus, so a walk that types everything as `SsrBundle` yields zero strata and looks exactly
    // like a walk over a genome with no repeats in it. Keeping the tally is what lets the failure
    // below say *which* it was.
    let mut typing = TypedRegionCounts::default();
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
        add_typing_counts(&mut typing, walk.counts());
    }

    // **An empty survey must not look like a finished one.** A resumable driver treats a written
    // result as done, so a run that typed nothing usable would be skipped forever
    // (`scripts/ng_str_library_survey.sh`). Refuse before anything is written, and name the cause
    // from the tally rather than leaving the reader to guess.
    if strata.is_empty() {
        return Err(Box::new(EmptySurvey { min_copies, typing }));
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
            "#rg\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
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
            group.file.display(),
        )?;
    }
    // **`file` is last and it is what makes two runs joinable.** The numeric `read_group` is minted
    // per run and means nothing across them; `(file, rg_id)` is the stable identity, because the SAM
    // specification makes `@RG ID` unique within its file. A survey run in batches merges on that
    // pair — which `scripts/ng_str_library_survey.sh` does.
    writeln!(
        out,
        "#rg_columns\tread_group\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\t\
         experiment_origin\tplatform\tloci\treads\tmean_tract_bases_per_read\tfile"
    )?;

    // **What region typing gave, beside what the reads said** — after the `#rg` lines, so a merge
    // that keys this block by file has already read the file name. A library whose strata are thin
    // because its tracts bundled is a different finding from one whose library is shallow, and this
    // line is what separates them.
    writeln!(
        out,
        "#typing\t{}\t{}\t{}\t{}\t{}\t{}",
        typing.spans,
        typing.ssr_loci,
        typing.ssr_bundles,
        typing.ssr_bundle_bp,
        typing.satellites,
        typing.repeat_bp_with_no_locus,
    )?;
    writeln!(
        out,
        "#typing_columns\tspans\tssr_loci\tssr_bundles\tssr_bundle_bp\tsatellites\t\
         repeat_bp_with_no_locus"
    )?;

    // **The settings this walk actually ran under, so a merge can refuse to mix two of them.**
    // The driver stamps what it was *asked* for, which cannot cover the copy floors and the bundle
    // radius when they come from the binary's own defaults — so a rebuild part way through a survey
    // changes the walk and leaves nothing to say so. That nearly happened on 2026-08-07, when the
    // radius default moved from 20 to 15 mid-survey. Emitting it per file makes the mix detectable
    // after the fact, whoever caused it.
    writeln!(
        out,
        "#config\tbundle_threshold={}\tmin_copies={}",
        walk_config.criteria.bundle_threshold,
        (1..=6u8)
            .map(|p| walk_config.criteria.min_copies.for_period(p).to_string())
            .collect::<Vec<_>>()
            .join(","),
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
    let mut min_copies: Option<MinCopies> = None;
    let mut bundle_threshold: Option<u64> = None;
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
            "--bundle-threshold" => match rest.next().and_then(|v| v.parse::<u64>().ok()) {
                Some(bases) if bases >= 1 => bundle_threshold = Some(bases),
                _ => {
                    eprintln!(
                        "error: --bundle-threshold needs a whole number of bases, at least 1 \
                         (ng's default is 30)"
                    );
                    return ExitCode::from(2);
                }
            },
            "--min-copies" => match rest.next().as_deref().and_then(parse_min_copies) {
                Some(floors) => min_copies = Some(floors),
                None => {
                    eprintln!(
                        "error: --min-copies needs one whole number for every period, or a \
                         comma-separated table of six — `--min-copies 6,2,4,3,3,3` lowers period 2 \
                         alone and leaves the rest at ng's defaults"
                    );
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
        bundle_threshold,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
