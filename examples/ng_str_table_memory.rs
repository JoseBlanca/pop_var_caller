//! **What the STR path's per-locus table costs on real data** — the measurement three open
//! questions in `doc/devel/ng/arch/parameter_prepass_ssr.md` all wait on.
//!
//! ng step 4's STR accumulator holds **one entry per locus**, keyed by
//! `(read group, motif period, reference repeat count)` and holding how many of that locus's reads
//! fell at each whole-repeat offset from the reference tract length. Entries of identical shape are
//! counted together, so the table's size is the number of *distinct locus shapes* — and that grows
//! with a locus's depth, which is why a per-locus read cap exists. Nothing has measured where
//! between its two bounds real data sits: a few hundred shapes a stratum at three reads a locus,
//! or one entry per locus at 300×.
//!
//! This walks region typing and the STR locus generator over real alignments — the same two steps
//! the real walk uses — and builds that table at several read caps at once. It answers:
//!
//! 1. **`MAX_LOCUS_READS`** — the entry count against the cap, so the knee is visible rather than
//!    argued (arch §2.1).
//! 2. **`ALLELE_OFFSET_LIMIT`** — how far a sample's tract lengths actually sit from the reference
//!    length, per stratum. That is the width the fit's genotype frequencies span, and the one the
//!    measurements showed decides the answer (arch §2.1, spec §8.1).
//! 3. **`OFFSET_HALF_RANGE`** — how much of the read mass the saturating end buckets absorb at ±4,
//!    which is the number that says whether four is comfortable.
//!
//! It also reports the **guard-bucket share** per stratum against `GUARD_SHARE_LIMIT`, which is
//! spec §5's threshold measured at per-stratum grain rather than in three repeat-count bands.
//!
//! ```text
//! ng_str_table_memory [--contigs a,b] [--regions r.bed] [--caps 4,8,12,20,0] \
//!     <reference.fa> <sample.bam|cram>
//! ```
//!
//! `--caps` takes the per-locus read caps to sweep; `0` means no cap, which is the upper bound.
//! `--regions` and `--contigs` restrict the walk, and **`--regions` is what makes a run over
//! region-restricted inputs finish** — typing every tract in a reference no read reaches is pure
//! waste (see `ng_ssr_cohort_stutter`, which this borrows its walk from).

use std::collections::BTreeMap;
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
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::GenomeRegion;
use pop_var_caller::ng::types::{Bp, ContigId, ReadGroupId};
use pop_var_caller::regions::ContigBounds;

#[path = "shared/catalog_regions.rs"]
mod catalog_regions;

// ---------------------------------------------------------------------------
// The types the architecture doc specifies, as far as this measurement needs them
// ---------------------------------------------------------------------------

/// Offsets are recorded over `±OFFSET_HALF_RANGE`, the ends saturating (arch §2.1).
const OFFSET_HALF_RANGE: i32 = 4;
/// Nine offset buckets, plus one for reads whose length is not a whole number of copies.
const OFFSET_BUCKETS: usize = (2 * OFFSET_HALF_RANGE + 1) as usize;
const GUARD_BUCKET: usize = OFFSET_BUCKETS;
const TOTAL_BUCKETS: usize = OFFSET_BUCKETS + 1;

/// Above this share of the reads that differ from the reference length, the stratum is one this
/// noise model does not describe (spec §5).
const GUARD_SHARE_LIMIT: f64 = 0.10;

/// One group of loci that gets its own fitted parameters: a motif period and the **reference**
/// tract's repeat count. Both are properties of the reference, so every sample strata-fies
/// identically.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Stratum {
    period: u8,
    repeats: u32,
}

/// One locus's reads laid out across the buckets — the table's key. `counts.iter().sum() == depth`
/// always, which is the "no bucket is charged a negative number of reads" invariant.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct LocusShape {
    counts: [u8; TOTAL_BUCKETS],
}

impl LocusShape {
    fn depth(&self) -> u32 {
        self.counts.iter().map(|&c| u32::from(c)).sum()
    }
}

/// One stratum's evidence, for one read group.
#[derive(Default)]
struct StratumTable {
    /// The distinct shapes and how many loci had each.
    loci_by_shape: BTreeMap<LocusShape, u32>,
    loci: u64,
    reads: u64,
    reads_in_end_buckets: u64,
    reads_off_reference: u64,
    reads_not_whole_repeat: u64,
}

/// What one entry costs: the key, the count, and a `BTreeMap` node's share of its overhead.
///
/// The key is 10 bytes and the count 4, so the payload is 14. A `BTreeMap` stores entries in
/// nodes of up to 11, with two pointers and a length per node, so the amortised overhead is small
/// — but it is not zero and it is not worth pretending to know exactly. **Sixteen bytes of
/// overhead an entry is the figure used here, and it is an assumption, not a measurement**: the
/// entry count either side of it is what this program actually establishes.
const BYTES_PER_ENTRY: u64 = 14 + 16;

// ---------------------------------------------------------------------------
// One locus → one shape
// ---------------------------------------------------------------------------

/// The stratum a locus belongs to, from the **reference** tract alone, or `None` where the
/// reference tract is not a whole number of motif copies.
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

/// A tiny reproducible generator, seeded from the locus's position so a region-sharded walk and a
/// single-threaded one keep the same reads (spec §4.1).
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in `0..bound`, by rejection so the result is not skewed by the modulo.
    fn below(&mut self, bound: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % bound;
        loop {
            let drawn = self.next_u64();
            if drawn < zone {
                return drawn % bound;
            }
        }
    }
}

/// Draw `cap` of the reads a shape holds, without replacement — the subsample the design specifies,
/// which is exact rather than approximate because thinning a locus's reads uniformly leaves the
/// bucket counts distributed as they would be at the lower depth.
fn subsample(counts: &[u32; TOTAL_BUCKETS], cap: u32, seed: u64) -> [u8; TOTAL_BUCKETS] {
    let depth: u32 = counts.iter().sum();
    if cap == 0 || depth <= cap {
        let mut out = [0u8; TOTAL_BUCKETS];
        for (slot, &count) in out.iter_mut().zip(counts) {
            *slot = count.min(u32::from(u8::MAX)) as u8;
        }
        return out;
    }
    let mut remaining = *counts;
    let mut left = depth;
    let mut kept = [0u32; TOTAL_BUCKETS];
    let mut rng = SplitMix64(seed);
    for _ in 0..cap {
        let mut pick = rng.below(u64::from(left));
        for (bucket, count) in remaining.iter_mut().enumerate() {
            if pick < u64::from(*count) {
                *count -= 1;
                kept[bucket] += 1;
                break;
            }
            pick -= u64::from(*count);
        }
        left -= 1;
    }
    let mut out = [0u8; TOTAL_BUCKETS];
    for (slot, &count) in out.iter_mut().zip(&kept) {
        *slot = count.min(u32::from(u8::MAX)) as u8;
    }
    out
}

/// Everything one run accumulates, at one read cap.
#[derive(Default)]
struct Tables {
    cap: u32,
    by_stratum: BTreeMap<(ReadGroupId, Stratum), StratumTable>,
}

impl Tables {
    fn entries(&self) -> u64 {
        self.by_stratum
            .values()
            .map(|t| t.loci_by_shape.len() as u64)
            .sum()
    }
    fn loci(&self) -> u64 {
        self.by_stratum.values().map(|t| t.loci).sum()
    }
}

/// How far this sample's own tract lengths sit from the reference length: the locus's **modal**
/// observed whole-repeat offset, tallied. This is what `ALLELE_OFFSET_LIMIT` has to cover — not
/// the recorded range, which the harness measured to matter far less.
#[derive(Default)]
struct AlleleOffsets {
    /// Index `i` is offset `i - MAX_TRACKED`, saturating at the ends.
    counts: Vec<u64>,
    beyond_low: u64,
    beyond_high: u64,
    loci: u64,
}

const MAX_TRACKED_ALLELE_OFFSET: i32 = 20;

impl AlleleOffsets {
    fn new() -> Self {
        Self {
            counts: vec![0; (2 * MAX_TRACKED_ALLELE_OFFSET + 1) as usize],
            ..Default::default()
        }
    }
    fn observe(&mut self, offset: i32) {
        self.loci += 1;
        if offset < -MAX_TRACKED_ALLELE_OFFSET {
            self.beyond_low += 1;
        } else if offset > MAX_TRACKED_ALLELE_OFFSET {
            self.beyond_high += 1;
        } else {
            self.counts[(offset + MAX_TRACKED_ALLELE_OFFSET) as usize] += 1;
        }
    }
    /// The narrowest symmetric limit holding at least `share` of the loci.
    fn limit_covering(&self, share: f64) -> i32 {
        let wanted = (self.loci as f64 * share).ceil() as u64;
        let mut held = self.counts[MAX_TRACKED_ALLELE_OFFSET as usize];
        for limit in 1..=MAX_TRACKED_ALLELE_OFFSET {
            if held >= wanted {
                return limit - 1;
            }
            held += self.counts[(MAX_TRACKED_ALLELE_OFFSET - limit) as usize]
                + self.counts[(MAX_TRACKED_ALLELE_OFFSET + limit) as usize];
        }
        MAX_TRACKED_ALLELE_OFFSET
    }
}

/// Fold one locus into every cap's table, and into the allele-offset tally.
fn observe_locus(
    tables: &mut [Tables],
    alleles: &mut AlleleOffsets,
    locus: &SampleLocusObservations,
    skipped_no_whole_reference: &mut u64,
) {
    let Some(stratum) = stratum_of(locus) else {
        *skipped_no_whole_reference += 1;
        return;
    };
    let period = usize::from(stratum.period);
    let reference_len = locus.reference_bases.len() as i64;

    // One shape per read group that covered the locus: chemistry is per read group, and each
    // entry's own distribution is correctly specified, so the split costs precision and not
    // correctness (spec §4.1).
    let mut per_group: BTreeMap<ReadGroupId, [u32; TOTAL_BUCKETS]> = BTreeMap::new();
    // The modal observed offset needs every group's reads together — it is a property of the
    // locus, not of a chemistry.
    let mut whole_offsets: BTreeMap<i32, u32> = BTreeMap::new();

    for obs in &locus.observations {
        if obs.read_witness != ReadWitness::Complete {
            continue;
        }
        let difference = obs.bases.len() as i64 - reference_len;
        let buckets = per_group
            .entry(obs.read_group)
            .or_insert([0; TOTAL_BUCKETS]);
        if difference % period as i64 != 0 {
            buckets[GUARD_BUCKET] += obs.num_obs;
            continue;
        }
        let offset = (difference / period as i64) as i32;
        buckets
            [(offset.clamp(-OFFSET_HALF_RANGE, OFFSET_HALF_RANGE) + OFFSET_HALF_RANGE) as usize] +=
            obs.num_obs;
        *whole_offsets.entry(offset).or_default() += obs.num_obs;
    }
    if per_group.is_empty() {
        return;
    }
    if let Some((&offset, _)) = whole_offsets
        .iter()
        .max_by_key(|&(offset, &count)| (count, -offset))
    {
        alleles.observe(offset);
    }

    let seed = u64::from(locus.region.contig.0) << 40 ^ locus.region.start.get();
    for (group, counts) in per_group {
        let depth: u32 = counts.iter().sum();
        let off_reference: u32 = depth - counts[OFFSET_HALF_RANGE as usize];
        let in_ends = counts[0] + counts[OFFSET_BUCKETS - 1];
        for table in tables.iter_mut() {
            let shape = LocusShape {
                counts: subsample(&counts, table.cap, seed),
            };
            if shape.depth() == 0 {
                continue;
            }
            let entry = table.by_stratum.entry((group, stratum)).or_default();
            *entry.loci_by_shape.entry(shape).or_default() += 1;
            entry.loci += 1;
            entry.reads += u64::from(depth);
            entry.reads_off_reference += u64::from(off_reference);
            entry.reads_in_end_buckets += u64::from(in_ends);
            entry.reads_not_whole_repeat += u64::from(counts[GUARD_BUCKET]);
        }
    }
}

// ---------------------------------------------------------------------------
// The walk — borrowed from ng_ssr_cohort_stutter
// ---------------------------------------------------------------------------

fn run(
    fasta: &Path,
    alignments: &Path,
    contig_filter: &[String],
    regions_bed: Option<&Path>,
    caps: &[u32],
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs: ContigList = info.contig_list();
    // **The typed regions come from the catalog beside the reference**, checked against what
    // the pass just reported. No catalog, no run: the error names the command that writes one.
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    let reference = OpenReference::new(info);

    let inputs = [alignments.to_path_buf()];
    let read_groups = build_read_groups(&inputs)?;
    let sample_entries = read_groups.read_groups_per_sample();
    let sample = SampleReads::open(
        &sample_entries[0],
        &read_groups,
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;
    eprintln!(
        "  sample {} with {} read group(s)",
        sample.sample_name(),
        sample_entries[0].read_groups.len()
    );

    let walk_config = TypedRegionConfig::default();
    let bundle_threshold = Bp(walk_config.criteria.bundle_threshold);
    #[expect(
        clippy::arc_with_non_send_sync,
        reason = "RawRefSeq is implemented for Arc only; this walk is single-threaded"
    )]
    let walk_reference_shared = Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));
    let mut generator = SsrGenerator::with_default_aligner(
        Arc::clone(&walk_reference_shared),
        {
            let reference = Arc::clone(&walk_reference_shared);
            move || Arc::clone(&reference)
        },
        SsrGeneratorConfig::default(),
        bundle_threshold,
    )?;

    let mut tables: Vec<Tables> = caps
        .iter()
        .map(|&cap| Tables {
            cap,
            ..Default::default()
        })
        .collect();
    let mut alleles = AlleleOffsets::new();
    let mut skipped_no_whole_reference = 0u64;
    let mut typed_segments = 0u64;

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

    // The stretches to ask the catalog for, labelled. The spans are collected rather than the
    // readers, because a reader borrows the catalog.
    let mut walks: Vec<(String, Vec<GenomeRegion>)> = Vec::new();
    match bed_spans {
        Some(spans) => {
            eprintln!("  {} BED spans", spans.iter().count());
            walks.push((
                "BED".to_string(),
                spans.iter().collect::<Vec<GenomeRegion>>(),
            ));
        }
        None => {
            for (index, entry) in contigs.entries.iter().enumerate() {
                if !wanted_contig(ContigId(index as u32)) {
                    continue;
                }
                walks.push((
                    entry.name.clone(),
                    vec![catalog_regions::whole_contig(
                        ContigId(index as u32),
                        entry.length,
                    )],
                ));
            }
        }
    }

    let criteria = StrRepeatCriteria::from(&walk_config);
    for (label, spans) in walks {
        eprintln!("  walking {label}");
        let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(&spans))?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            if !wanted_contig(region.region.contig) {
                continue;
            }
            typed_segments += 1;
            generator.begin_segment(region.region);
            while let Some(locus) = generator.next_locus(segment, &sample)? {
                observe_locus(
                    &mut tables,
                    &mut alleles,
                    &locus,
                    &mut skipped_no_whole_reference,
                );
            }
        }
        eprintln!(
            "    {} typed segments so far, {} loci entered, {} entries at cap {}",
            typed_segments,
            tables[0].loci(),
            tables[0].entries(),
            tables[0].cap
        );
    }

    report(
        &tables,
        &alleles,
        typed_segments,
        skipped_no_whole_reference,
    );

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

fn report(tables: &[Tables], alleles: &AlleleOffsets, segments: u64, skipped: u64) {
    println!("\n## What the per-locus table costs\n");
    println!(
        "{segments} typed STR segments; {skipped} loci skipped because the reference tract is not \
         a whole number of motif copies."
    );
    println!(
        "\nOne entry per distinct locus shape, per (read group, period, reference repeat count).\n\
         `cap 0` is no cap at all — the upper bound, one entry per locus wherever depth is high.\n"
    );
    println!(
        "{:>10} {:>14} {:>14} {:>12} {:>12} {:>14}",
        "read cap", "loci entered", "entries", "entries/locus", "strata", "table (MB)"
    );
    for table in tables {
        let entries = table.entries();
        let loci = table.loci();
        println!(
            "{:>10} {:>14} {:>14} {:>12.3} {:>12} {:>14.2}",
            if table.cap == 0 {
                "none".to_string()
            } else {
                table.cap.to_string()
            },
            loci,
            entries,
            if loci > 0 {
                entries as f64 / loci as f64
            } else {
                0.0
            },
            table.by_stratum.len(),
            (entries * BYTES_PER_ENTRY) as f64 / (1024.0 * 1024.0),
        );
    }
    println!(
        "\n`table (MB)` assumes {BYTES_PER_ENTRY} bytes an entry — a 10-byte key, a 4-byte count \
         and 16 bytes of\nmap overhead, the last of which is an assumption. The entry count is what \
         is measured."
    );

    // Where the entries are: the strata carrying most of them.
    if let Some(table) = tables.last() {
        println!("\n### The ten strata holding the most entries, at the widest cap\n");
        println!(
            "{:>8} {:>10} {:>12} {:>12} {:>14} {:>14} {:>12}",
            "period", "repeats", "loci", "entries", "reads", "end buckets", "guard"
        );
        let mut rows: Vec<_> = table.by_stratum.iter().collect();
        rows.sort_by_key(|(_, t)| std::cmp::Reverse(t.loci_by_shape.len()));
        for ((_, stratum), t) in rows.iter().take(10) {
            println!(
                "{:>8} {:>10} {:>12} {:>12} {:>14} {:>13.3}% {:>11.3}%",
                stratum.period,
                stratum.repeats,
                t.loci,
                t.loci_by_shape.len(),
                t.reads,
                if t.reads > 0 {
                    100.0 * t.reads_in_end_buckets as f64 / t.reads as f64
                } else {
                    0.0
                },
                if t.reads_off_reference > 0 {
                    100.0 * t.reads_not_whole_repeat as f64 / t.reads_off_reference as f64
                } else {
                    0.0
                },
            );
        }

        // **Every stratum, as a table a spreadsheet can take.** This is what decides where the
        // STR path's copy floors go: an STR locus is one that is *likely to stutter*, not merely
        // one that contains a repeat, so the floor is the repeat count at which a period's tracts
        // start behaving that way. `off_ref` is reads differing from the reference length as a
        // share of all reads — a slippage proxy, and an over-estimate, since a locus's own alleles
        // differ from the reference too. `guard` is the share of *those differing reads* that
        // differ by something other than a whole number of copies, which is the half of the
        // question the model can be wrong about rather than merely small.
        println!("\n### Every stratum\n");
        println!(
            "period\trepeats\tloci\treads\toff_ref_reads\toff_ref_share\tnot_whole\tguard_share\tend_bucket_share"
        );
        let mut all: Vec<_> = table.by_stratum.iter().collect();
        all.sort_by_key(|((_, stratum), _)| (stratum.period, stratum.repeats));
        for ((_, stratum), t) in all {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{:.5}\t{:.5}",
                stratum.period,
                stratum.repeats,
                t.loci,
                t.reads,
                t.reads_off_reference,
                if t.reads > 0 {
                    t.reads_off_reference as f64 / t.reads as f64
                } else {
                    0.0
                },
                t.reads_not_whole_repeat,
                if t.reads_off_reference > 0 {
                    t.reads_not_whole_repeat as f64 / t.reads_off_reference as f64
                } else {
                    0.0
                },
                if t.reads > 0 {
                    t.reads_in_end_buckets as f64 / t.reads as f64
                } else {
                    0.0
                },
            );
        }

        // The two constants this run is here to settle, over the whole table.
        let reads: u64 = table.by_stratum.values().map(|t| t.reads).sum();
        let in_ends: u64 = table
            .by_stratum
            .values()
            .map(|t| t.reads_in_end_buckets)
            .sum();
        let off_reference: u64 = table
            .by_stratum
            .values()
            .map(|t| t.reads_off_reference)
            .sum();
        let guard: u64 = table
            .by_stratum
            .values()
            .map(|t| t.reads_not_whole_repeat)
            .sum();
        println!("\n### The recorded offset range, at ±{OFFSET_HALF_RANGE}\n");
        println!(
            "  {:.4}% of reads land in a saturating end bucket ({in_ends} of {reads}).",
            if reads > 0 {
                100.0 * in_ends as f64 / reads as f64
            } else {
                0.0
            }
        );
        println!(
            "  {:.3}% of the reads that differ from the reference length differ by something \
             that is not\n  a whole number of copies ({guard} of {off_reference}) — the guard \
             bucket, against a\n  threshold of {:.0}%.",
            if off_reference > 0 {
                100.0 * guard as f64 / off_reference as f64
            } else {
                0.0
            },
            100.0 * GUARD_SHARE_LIMIT
        );
        let above: usize = table
            .by_stratum
            .values()
            .filter(|t| {
                t.reads_off_reference > 0
                    && t.reads_not_whole_repeat as f64 / t.reads_off_reference as f64
                        > GUARD_SHARE_LIMIT
            })
            .count();
        println!(
            "  {above} of {} strata sit above that threshold — the ones this noise model does \
             not describe.",
            table.by_stratum.len()
        );
    }

    println!("\n### How far this sample's tract lengths sit from the reference length\n");
    println!(
        "The locus's modal observed whole-repeat offset, over {} loci. **This is the width the\n\
         fit's allele support has to cover** — not the recorded offset range.\n",
        alleles.loci
    );
    println!("{:>10} {:>16} {:>12}", "offset", "loci", "share");
    for offset in -8..=8 {
        let index = (offset + MAX_TRACKED_ALLELE_OFFSET) as usize;
        let count = alleles.counts[index];
        if count == 0 {
            continue;
        }
        println!(
            "{:>10} {:>16} {:>11.3}%",
            offset,
            count,
            100.0 * count as f64 / alleles.loci.max(1) as f64
        );
    }
    let far: u64 = alleles.beyond_low
        + alleles.beyond_high
        + (0..(MAX_TRACKED_ALLELE_OFFSET - 8))
            .map(|i| {
                alleles.counts[i as usize]
                    + alleles.counts[(2 * MAX_TRACKED_ALLELE_OFFSET - i) as usize]
            })
            .sum::<u64>();
    println!(
        "{:>10} {:>16} {:>11.3}%",
        "beyond ±8",
        far,
        100.0 * far as f64 / alleles.loci.max(1) as f64
    );
    println!(
        "\n  narrowest symmetric limit holding 99% of loci:    ±{}",
        alleles.limit_covering(0.99)
    );
    println!(
        "  narrowest holding 99.9%:                          ±{}",
        alleles.limit_covering(0.999)
    );
    println!(
        "  narrowest holding 99.99%:                         ±{}",
        alleles.limit_covering(0.9999)
    );
}

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut contig_filter: Vec<String> = Vec::new();
    let mut regions_bed: Option<PathBuf> = None;
    let mut caps: Vec<u32> = vec![4, 8, 12, 20, 0];
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
            "--caps" => match rest.next() {
                Some(list) => {
                    caps = list
                        .split(',')
                        .filter_map(|s| s.trim().parse().ok())
                        .collect();
                }
                None => {
                    eprintln!("error: --caps needs a comma-separated list");
                    return ExitCode::from(2);
                }
            },
            _ => positional.push(arg),
        }
    }
    if positional.len() != 2 || caps.is_empty() {
        eprintln!(
            "usage: ng_str_table_memory [--contigs a,b] [--regions r.bed] [--caps 4,8,12,20,0] \
             <reference.fa> <sample.bam|cram>\n\
             measures the STR path's per-locus table on real reads: entries against the read cap, \
             how far tract lengths sit from the reference, and what the end buckets absorb."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&positional[0]);
    let alignments = PathBuf::from(&positional[1]);
    eprintln!("caps: {caps:?}");

    match run(
        &fasta,
        &alignments,
        &contig_filter,
        regions_bed.as_deref(),
        &caps,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
