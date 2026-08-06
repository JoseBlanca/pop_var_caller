//! **SKETCH 1 — what the generator emits, and what the parameter pre-pass reads.**
//!
//! Throwaway. Built for `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md` §4 and
//! deleted after the decision. Forked from `ng_generic_walk_probe`, which walks the same
//! pipeline and drops every locus; this one hands each locus to a stand-in for the
//! parameter pre-pass instead.
//!
//! ```text
//! cargo run --release --example ng_prepass_sketch -- <reference.fa> <sample.bam|cram> [contig]
//! ```
//!
//! # The four arms
//!
//! - `PVC_SKETCH_ARM=A` — today's shipped generator, one owned locus per covered base,
//!   consumed record-shaped. **The baseline**, measured here rather than cited so that all
//!   the arms are one experiment.
//! - `PVC_SKETCH_ARM=B` — the walk fills a block of parallel arrays the caller owns and
//!   reuses; the pre-pass reads columns. No locus is ever materialised.
//! - `PVC_SKETCH_ARM=C` — the same block, but the pre-pass is handed one locus at a time
//!   through a view over a buffer **it** owns and refills per locus. Every line of the
//!   pre-pass stays record-shaped.
//! - `PVC_SKETCH_ARM=Cb` — arm C's consumer reading **borrowed** slices straight out of the
//!   block, with no per-locus refill. This is the shape production built and rejected
//!   (`BlockColumns<'a>`, test-only, because a consumer holding a locus across calls is
//!   self-referential). It is here only to price arm C's copy: `C − Cb` is what refilling
//!   the consumer's own buffer costs.
//!
//! Every arm runs the **same** accumulation arithmetic — [`Accumulators::add_site`] — over
//! the same reduction of a locus to `(depth, alt reads)` per read group and whole. Only the
//! ten lines that walk a locus differ. `summary=` is printed by all four and must agree.
//!
//! # Knobs
//!
//! Everything `ng_generic_walk_probe` reads, plus:
//!
//! - `PVC_SKETCH_ARM` — `A` (default), `B`, `C`, `Cb`.
//! - `PVC_SKETCH_BLOCK_KB` — the block's payload budget in KiB (default 256). Block size
//!   is the known memory lever in this project, so it is swept, not chosen.
//! - `PVC_SKETCH_SAMPLE_ONE_IN` — the random-locus sampler's rate (default 29,000, which
//!   is 100,000 loci out of a human genome at 30×). `0` switches the sampler off.
//! - `PVC_SKETCH_NO_FIT=1` — skip the end-of-run parameter fit.
//! - `PVC_SKETCH_DIGEST=1` — also digest every field of every locus, which pins the arms
//!   to the same payload and not merely to the same summary. Off for timed runs, because
//!   it is per-observation work no real consumer does.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::rc::Rc;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::block::{LocusBlock, LocusSink};
use pop_var_caller::ng::locus_generation::pileup::{
    PileupGenerator, PileupGeneratorConfig, PileupGeneratorCounts,
};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, ReadWitness, SampleLocusObservations,
    SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{
    RegionKind, TypedRegion, TypedRegionConfig, TypedRegionError, TypedRegionIterator,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::{reference_check_from_env, reference_check_label};

// ---------------------------------------------------------------------------------------
// The parameter pre-pass stand-in
// ---------------------------------------------------------------------------------------
//
// Faithful to `doc/devel/ng/spec/parameter_prepass{,_generic,_census_sites}.md` in what it
// accumulates, not in what it estimates:
//
// - two accumulated objects, not three: a `(depth, alt reads)` histogram keyed by
//   `(read group, ploidy)`, and one keyed by `(contig, 100 kb window, ploidy)`. The
//   whole-sample histogram is the second summed over its windows and is not built;
// - depth bins are exact to eight and widen geometrically above, the table is ragged
//   (`alt ≤ depth`), and `depth_sums` is kept **per cell**, not per bin;
// - a site deeper than the ladder's top is scaled down, its alt count scaled by the same
//   ratio and rounded **stochastically, keyed on the position** so that every arm and every
//   shard draws the same number;
// - depth counts complete witnesses only; `reads_without_observation` does not enter it
//   and `reads_discarded_by_cap` does not skip the locus;
// - base qualities and mapping qualities are not read: the model has one error rate per
//   read group;
// - the census keeps positions by a **hash threshold**, a pure function of the position —
//   not a reservoir — so the kept set is identical in every arm.

/// The depth ladder. Exact integers to eight, widening geometrically above; the table is
/// ragged because a bin covering depths 100–124 holds cells up to `alt = 124`.
struct DepthBinEdges {
    /// Inclusive upper depth of each bin.
    upper: Vec<u32>,
    /// Where each bin's row starts in the flat cell table.
    row_start: Vec<u32>,
    cells: usize,
    max_depth: u32,
}

impl DepthBinEdges {
    fn new() -> Self {
        let upper: Vec<u32> = vec![
            0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 19, 24, 31, 39, 49, 59, 79, 99, 124,
        ];
        let mut row_start = Vec::with_capacity(upper.len());
        let mut at = 0u32;
        for bound in &upper {
            row_start.push(at);
            at += bound + 1;
        }
        Self {
            max_depth: *upper.last().expect("the ladder is not empty"),
            cells: at as usize,
            upper,
            row_start,
        }
    }

    /// The bin a depth falls in. Depths are scaled to `max_depth` before this is asked.
    #[inline]
    fn bin_of(&self, depth: u32) -> usize {
        // Twenty-one bins and the mass at the bottom, so a forward scan beats a search.
        for (bin, bound) in self.upper.iter().enumerate() {
            if depth <= *bound {
                return bin;
            }
        }
        self.upper.len() - 1
    }

    /// The flat cell index of `(depth, alt)`.
    #[inline]
    fn cell_of(&self, depth: u32, alt: u32) -> usize {
        let bin = self.bin_of(depth);
        self.row_start[bin] as usize + alt.min(self.upper[bin]) as usize
    }
}

/// One `(depth bin × alt count)` table, with the exact depths summed per cell.
struct DepthAltHistogram {
    counts: Vec<u32>,
    depth_sums: Vec<u64>,
}

impl DepthAltHistogram {
    fn new(cells: usize) -> Self {
        Self {
            counts: vec![0; cells],
            depth_sums: vec![0; cells],
        }
    }
}

/// Everything the walk is folded into, plus the loci the census kept.
struct Accumulators {
    edges: Rc<DepthBinEdges>,
    by_read_group: BTreeMap<(u32, u8), DepthAltHistogram>,
    by_window: BTreeMap<(u32, u32, u8), DepthAltHistogram>,
    total_loci: u64,
    total_covered_positions: u64,
    observations: u64,
    /// Σ `num_obs` over every observation — one per **read at a position**, which is the
    /// denominator the shipped price table calls a *contributor visit*.
    visits: u64,
    complete_observations: u64,
    reads_without_observation: u64,
    loci_subsampled: u64,
    sites_scaled_down: u64,
    /// Reduction scratch, one entry per read group present at the locus.
    groups: Vec<(u32, u32, u32)>,
    /// One in this many positions is kept. Zero switches the census off.
    sample_one_in: u64,
    kept: Vec<SampleLocusObservations>,
    /// Bytes the census copied out of a block. Zero on the record path, which moves.
    kept_bytes: u64,
    /// Optional payload digest, over every field of every locus.
    digest_on: bool,
    digest: u64,
}

/// 100 kb, as the spec's windowed histogram uses.
const WINDOW_BP: u32 = 100_000;
/// Diploid everywhere in the stand-in; the real thing resolves ploidy per region.
const PLOIDY: u8 = 2;

/// SplitMix64 — the census's position hash, and the stochastic rounder's seed.
#[inline]
fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Accumulators {
    fn new(sample_one_in: u64, digest_on: bool) -> Self {
        Self {
            edges: Rc::new(DepthBinEdges::new()),
            by_read_group: BTreeMap::new(),
            by_window: BTreeMap::new(),
            total_loci: 0,
            total_covered_positions: 0,
            observations: 0,
            visits: 0,
            complete_observations: 0,
            reads_without_observation: 0,
            loci_subsampled: 0,
            sites_scaled_down: 0,
            groups: Vec::new(),
            sample_one_in,
            kept: Vec::new(),
            kept_bytes: 0,
            digest_on,
            digest: 0xcbf2_9ce4_8422_2325,
        }
    }

    /// Scale a site deeper than the ladder down to its top, the alt count with it, rounding
    /// stochastically on a draw keyed to the position. Deterministic rounding is not a
    /// proportional subsample and biases the error rate.
    #[inline]
    fn scale_to_ladder(&self, contig: u32, start: u32, depth: u32, alt: u32) -> (u32, u32) {
        let max = self.edges.max_depth;
        if depth <= max {
            return (depth, alt);
        }
        let scaled = f64::from(alt) * f64::from(max) / f64::from(depth);
        let floor = scaled.floor();
        let draw = (mix64(
            (u64::from(contig) << 40) ^ (u64::from(start) << 8) ^ u64::from(alt),
        ) >> 11) as f64
            / (1u64 << 53) as f64;
        let alt_scaled = floor as u32 + u32::from(draw < scaled - floor);
        (max, alt_scaled.min(max))
    }

    /// **The one place a site enters the accumulators**, shared by every arm so that three
    /// implementations cannot disagree about anything but how they read a locus.
    ///
    /// `groups` holds `(read group, depth, alt reads)` once per read group that covered the
    /// site; `whole` is the same pair summed over them.
    fn add_site(&mut self, contig: u32, start: u32, span: u32, whole: (u32, u32)) {
        let cells = self.edges.cells;

        // Read-group histogram: a site enters once per read group that covered it, at that
        // group's own depth.
        for index in 0..self.groups.len() {
            let (group, depth, alt) = self.groups[index];
            if depth == 0 {
                continue;
            }
            let (depth, alt) = self.scale_to_ladder(contig, start, depth, alt);
            let cell = self.edges.cell_of(depth, alt);
            let histogram = self
                .by_read_group
                .entry((group, PLOIDY))
                .or_insert_with(|| DepthAltHistogram::new(cells));
            histogram.counts[cell] += 1;
            histogram.depth_sums[cell] += u64::from(depth);
        }

        // Windowed histogram: the site once, whole.
        if whole.0 > 0 {
            let (depth, alt) = self.scale_to_ladder(contig, start, whole.0, whole.1);
            if depth != whole.0 {
                self.sites_scaled_down += 1;
            }
            let cell = self.edges.cell_of(depth, alt);
            let histogram = self
                .by_window
                .entry((contig, start / WINDOW_BP, PLOIDY))
                .or_insert_with(|| DepthAltHistogram::new(cells));
            histogram.counts[cell] += 1;
            histogram.depth_sums[cell] += u64::from(depth);
        }

        self.total_loci += 1;
        self.total_covered_positions += u64::from(span);
    }

    /// Whether the census keeps the position — a hash threshold, so every arm and every
    /// sample of a cohort selects the identical set.
    #[inline]
    fn census_keeps(&self, contig: u32, start: u32) -> bool {
        self.sample_one_in != 0
            && mix64((u64::from(contig) << 32) | u64::from(start)) % self.sample_one_in == 0
    }

    #[inline]
    fn digest_u64(&mut self, value: u64) {
        self.digest = (self.digest ^ value).wrapping_mul(0x100_0000_01b3);
    }

    fn summary(&self) -> String {
        // Integer digests over the cell tables, in the `BTreeMap`'s fixed key order — the
        // gate that says three arms accumulated the same thing.
        let mut rg = 0xcbf2_9ce4_8422_2325u64;
        for ((group, ploidy), histogram) in &self.by_read_group {
            rg = fold(rg, u64::from(*group));
            rg = fold(rg, u64::from(*ploidy));
            for (cell, count) in histogram.counts.iter().enumerate() {
                if *count != 0 {
                    rg = fold(rg, cell as u64);
                    rg = fold(rg, u64::from(*count));
                    rg = fold(rg, histogram.depth_sums[cell]);
                }
            }
        }
        let mut win = 0xcbf2_9ce4_8422_2325u64;
        for ((contig, window, ploidy), histogram) in &self.by_window {
            win = fold(win, u64::from(*contig));
            win = fold(win, u64::from(*window));
            win = fold(win, u64::from(*ploidy));
            for (cell, count) in histogram.counts.iter().enumerate() {
                if *count != 0 {
                    win = fold(win, cell as u64);
                    win = fold(win, u64::from(*count));
                    win = fold(win, histogram.depth_sums[cell]);
                }
            }
        }
        let mut kept = 0xcbf2_9ce4_8422_2325u64;
        for locus in &self.kept {
            kept = fold(kept, locus.region.contig.0.into());
            kept = fold(kept, locus.region.start.get());
            kept = fold(kept, locus.observations.len() as u64);
            for observation in &locus.observations {
                for byte in observation.bases.iter() {
                    kept = fold(kept, u64::from(*byte));
                }
                kept = fold(kept, u64::from(observation.num_obs));
                kept = fold(kept, observation.q_sum.to_bits());
                kept = fold(kept, observation.chain_ids.len() as u64);
            }
        }
        format!(
            "summary loci={} covered={} visits={} observations={} complete={} rwo={} subsampled={} \
             scaled={} groups={} windows={} rg_digest={rg:016x} win_digest={win:016x} \
             kept={} kept_digest={kept:016x} payload_digest={:016x}",
            self.total_loci,
            self.total_covered_positions,
            self.visits,
            self.observations,
            self.complete_observations,
            self.reads_without_observation,
            self.loci_subsampled,
            self.sites_scaled_down,
            self.by_read_group.len(),
            self.by_window.len(),
            self.kept.len(),
            self.digest,
        )
    }

    /// The end-of-run fit: one error rate per read group, on the spec's ladder, by a
    /// profile scan over a small grid of genotype frequencies.
    ///
    /// Real arithmetic, run once, so it cancels exactly when two locus counts of one arm
    /// are differenced. It is here because a stand-in whose accumulators nothing ever reads
    /// can accumulate the wrong thing and never notice.
    fn fit(&self) -> String {
        let mut out = String::new();
        for ((group, ploidy), histogram) in &self.by_read_group {
            let (phred, het, hom_alt, score) = fit_error_rate(&self.edges, histogram);
            out.push_str(&format!(
                " fit rg={group} ploidy={ploidy} phred={phred:.2} het={het:.5} \
                 hom_alt={hom_alt:.5} lnL={score:.6}\n"
            ));
        }
        out
    }
}

#[inline]
fn fold(state: u64, value: u64) -> u64 {
    (state ^ value).wrapping_mul(0x100_0000_01b3)
}

/// The ladder: Phred 10 to 50 in quarter-decibel steps, 161 rungs about 6 % apart.
const LADDER_MIN_PHRED: f64 = 10.0;
const LADDER_MAX_PHRED: f64 = 50.0;
const LADDER_STEP_PHRED: f64 = 0.25;

/// Maximise the diploid site likelihood over `(ε, het, hom-alt)`.
///
/// ```text
///                  P
///   P(site)  =  Σ    π_j · p_j^k · (1 − p_j)^(n−k),
///                 j=0
///   p_j = (j/P)·(1 − ε/3) + (1 − j/P)·ε
/// ```
///
/// The `3` is load-bearing: writing `1 − ε` charges the alternative copy three times too
/// much reversion. Ties resolve to the **lower** error rate, stated so two implementations
/// cannot differ.
fn fit_error_rate(
    edges: &DepthBinEdges,
    histogram: &DepthAltHistogram,
) -> (f64, f64, f64, f64) {
    let hets: Vec<f64> = (0..11).map(|i| 1e-4 * 10f64.powf(i as f64 * 0.3)).collect();
    let hom_alts: Vec<f64> = (0..6).map(|i| 1e-5 * 10f64.powf(i as f64 * 0.6)).collect();
    // The cells that carry any evidence, with the mean depth **of the cell**.
    let mut cells: Vec<(f64, f64, f64)> = Vec::new();
    for bin in 0..edges.upper.len() {
        let row = edges.row_start[bin] as usize;
        for alt in 0..=edges.upper[bin] {
            let cell = row + alt as usize;
            let count = histogram.counts[cell];
            if count == 0 {
                continue;
            }
            let depth = histogram.depth_sums[cell] as f64 / f64::from(count);
            cells.push((depth, f64::from(alt), f64::from(count)));
        }
    }
    let rungs = ((LADDER_MAX_PHRED - LADDER_MIN_PHRED) / LADDER_STEP_PHRED) as usize + 1;
    let mut best = (LADDER_MIN_PHRED, hets[0], hom_alts[0], f64::NEG_INFINITY);
    for rung in 0..rungs {
        let phred = LADDER_MIN_PHRED + LADDER_STEP_PHRED * rung as f64;
        let epsilon = 10f64.powf(-phred / 10.0);
        for het in &hets {
            for hom_alt in &hom_alts {
                let hom_ref = 1.0 - het - hom_alt;
                if hom_ref <= 0.0 {
                    continue;
                }
                let p = [epsilon, 0.5 + epsilon / 3.0 - epsilon / 2.0, 1.0 - epsilon / 3.0];
                let pi = [hom_ref, *het, *hom_alt];
                let mut score = 0.0;
                for (depth, alt, count) in &cells {
                    let mut probability = 0.0;
                    for genotype in 0..3 {
                        probability += pi[genotype]
                            * p[genotype].powf(*alt)
                            * (1.0 - p[genotype]).powf(depth - alt);
                    }
                    score += count * probability.max(f64::MIN_POSITIVE).ln();
                }
                // Strictly greater, so a tie keeps the lower rung — the lower error rate.
                if score > best.3 {
                    best = (phred, *het, *hom_alt, score);
                }
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------------------
// Arm A and arm C's consumer: record-shaped
// ---------------------------------------------------------------------------------------

/// One observation, however it is stored. Arms A, C and Cb all reduce through this; arm B
/// reads the block's arrays directly, which is what makes it the columnar consumer.
struct ObsRef<'a> {
    bases: &'a [u8],
    complete: bool,
    read_group: u32,
    num_obs: u32,
}

/// The reduction: a locus to `(depth, alt reads)` whole and per read group.
///
/// **Complete witnesses only.** A read spanning part of the locus witnessed neither the
/// reference allele nor an alternative one at the positions it missed.
#[inline]
fn reduce_observation(
    accumulators: &mut Accumulators,
    reference_bases: &[u8],
    observation: ObsRef<'_>,
    whole: &mut (u32, u32),
) {
    accumulators.observations += 1;
    accumulators.visits += u64::from(observation.num_obs);
    if !observation.complete {
        return;
    }
    accumulators.complete_observations += 1;
    let alt = observation.bases != reference_bases;
    let depth = observation.num_obs;
    let alt_reads = if alt { depth } else { 0 };
    whole.0 += depth;
    whole.1 += alt_reads;
    match accumulators
        .groups
        .iter_mut()
        .find(|entry| entry.0 == observation.read_group)
    {
        Some(entry) => {
            entry.1 += depth;
            entry.2 += alt_reads;
        }
        None => accumulators
            .groups
            .push((observation.read_group, depth, alt_reads)),
    }
}

/// Arm A: an owned record. The census **moves** a kept locus rather than copying it, which
/// is the record path's honest cost.
fn consume_record(accumulators: &mut Accumulators, locus: SampleLocusObservations) {
    accumulators.groups.clear();
    let mut whole = (0u32, 0u32);
    accumulators.reads_without_observation += u64::from(locus.reads_without_observation);
    if locus.reads_discarded_by_cap > 0 {
        accumulators.loci_subsampled += 1;
    }
    for observation in &locus.observations {
        reduce_observation(
            accumulators,
            &locus.reference_bases,
            ObsRef {
                bases: &observation.bases,
                complete: observation.read_witness == ReadWitness::Complete,
                read_group: observation.read_group.0,
                num_obs: observation.num_obs,
            },
            &mut whole,
        );
    }
    if accumulators.digest_on {
        digest_record(accumulators, &locus);
    }
    let contig = locus.region.contig.0;
    let start = locus.region.start.get() as u32;
    accumulators.add_site(contig, start, locus.region.len() as u32, whole);
    if accumulators.census_keeps(contig, start) {
        accumulators.kept.push(locus);
    }
}

fn digest_record(accumulators: &mut Accumulators, locus: &SampleLocusObservations) {
    accumulators.digest_u64(u64::from(locus.region.contig.0));
    accumulators.digest_u64(locus.region.start.get());
    accumulators.digest_u64(locus.region.end.get());
    for byte in locus.reference_bases.iter() {
        accumulators.digest_u64(u64::from(*byte));
    }
    accumulators.digest_u64(u64::from(locus.reads_without_observation));
    accumulators.digest_u64(u64::from(locus.reads_discarded_by_cap));
    for observation in &locus.observations {
        for byte in observation.bases.iter() {
            accumulators.digest_u64(u64::from(*byte));
        }
        match &observation.read_witness {
            ReadWitness::Complete => accumulators.digest_u64(0),
            ReadWitness::Partial { positions } => {
                accumulators.digest_u64(1);
                for (from, to) in positions.runs() {
                    accumulators.digest_u64(u64::from(from));
                    accumulators.digest_u64(u64::from(to));
                }
            }
        }
        accumulators.digest_u64(u64::from(observation.read_group.0));
        accumulators.digest_u64(u64::from(observation.num_obs));
        accumulators.digest_u64(u64::from(observation.num_fwd));
        accumulators.digest_u64(observation.q_sum.to_bits());
        accumulators.digest_u64(u64::from(observation.mapq_sum));
        accumulators.digest_u64(observation.mapq_sum_sq);
        accumulators.digest_u64(u64::from(observation.placed_left));
        for id in &observation.chain_ids {
            accumulators.digest_u64(*id);
        }
    }
}

// ---------------------------------------------------------------------------------------
// Arm B: the columnar consumer
// ---------------------------------------------------------------------------------------

fn consume_block_columns(accumulators: &mut Accumulators, block: &LocusBlock) {
    for locus in 0..block.len() {
        accumulators.groups.clear();
        let mut whole = (0u32, 0u32);
        accumulators.reads_without_observation +=
            u64::from(block.reads_without_observation[locus]);
        if block.reads_discarded_by_cap[locus] > 0 {
            accumulators.loci_subsampled += 1;
        }
        let reference_bases = block.reference_bases_of(locus);
        for observation in block.observations_of(locus) {
            accumulators.observations += 1;
            accumulators.visits += u64::from(block.obs_num_obs[observation]);
            if block.obs_witness_kind[observation] != 0 {
                continue;
            }
            accumulators.complete_observations += 1;
            let depth = block.obs_num_obs[observation];
            let alt_reads = if block.bases_of(observation) != reference_bases {
                depth
            } else {
                0
            };
            whole.0 += depth;
            whole.1 += alt_reads;
            let group = block.obs_read_group[observation];
            match accumulators.groups.iter_mut().find(|entry| entry.0 == group) {
                Some(entry) => {
                    entry.1 += depth;
                    entry.2 += alt_reads;
                }
                None => accumulators.groups.push((group, depth, alt_reads)),
            }
        }
        if accumulators.digest_on {
            digest_block_locus(accumulators, block, locus);
        }
        let contig = block.contig[locus];
        let start = block.start[locus];
        let span = block.end[locus] - block.start[locus] + 1;
        accumulators.add_site(contig, start, span, whole);
        if accumulators.census_keeps(contig, start) {
            let (owned, bytes) = block.copy_out(locus);
            accumulators.kept_bytes += bytes as u64;
            accumulators.kept.push(owned);
        }
    }
}

fn digest_block_locus(accumulators: &mut Accumulators, block: &LocusBlock, locus: usize) {
    accumulators.digest_u64(u64::from(block.contig[locus]));
    accumulators.digest_u64(u64::from(block.start[locus]));
    accumulators.digest_u64(u64::from(block.end[locus]));
    for byte in block.reference_bases_of(locus) {
        accumulators.digest_u64(u64::from(*byte));
    }
    accumulators.digest_u64(u64::from(block.reads_without_observation[locus]));
    accumulators.digest_u64(u64::from(block.reads_discarded_by_cap[locus]));
    for observation in block.observations_of(locus) {
        for byte in block.bases_of(observation) {
            accumulators.digest_u64(u64::from(*byte));
        }
        if block.obs_witness_kind[observation] == 0 {
            accumulators.digest_u64(0);
        } else {
            accumulators.digest_u64(1);
            for (from, to) in block.runs_of(observation) {
                accumulators.digest_u64(u64::from(*from));
                accumulators.digest_u64(u64::from(*to));
            }
        }
        accumulators.digest_u64(u64::from(block.obs_read_group[observation]));
        accumulators.digest_u64(u64::from(block.obs_num_obs[observation]));
        accumulators.digest_u64(u64::from(block.obs_num_fwd[observation]));
        accumulators.digest_u64(block.obs_q_sum[observation].to_bits());
        accumulators.digest_u64(u64::from(block.obs_mapq_sum[observation]));
        accumulators.digest_u64(block.obs_mapq_sum_sq[observation]);
        accumulators.digest_u64(u64::from(block.obs_placed_left[observation]));
        for id in block.chain_ids_of(observation) {
            accumulators.digest_u64(*id);
        }
    }
}

// ---------------------------------------------------------------------------------------
// Arm C: a view over a buffer the consumer owns and refills per locus
// ---------------------------------------------------------------------------------------

/// One observation in the consumer's own buffer. Never freed, only cleared and refilled.
#[derive(Default)]
struct ScratchObservation {
    bases: Vec<u8>,
    complete: bool,
    runs: Vec<(u16, u16)>,
    read_group: u32,
    num_obs: u32,
    num_fwd: u32,
    q_sum: f64,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    placed_left: u32,
    chain_ids: Vec<u64>,
}

/// The buffer the pre-pass owns, refilled once per locus.
///
/// **This is the pattern that survived in production** (`LocusObservations<'a>`,
/// `src/paralog/locus_score.rs:62`) rather than the borrow into the decoder's block that
/// did not (`BlockColumns<'a>`, `src/psp/reader.rs`, test-only). A consumer may hold this
/// across calls because it owns it.
#[derive(Default)]
struct LocusScratch {
    contig: u32,
    start: u32,
    end: u32,
    reference_bases: Vec<u8>,
    reads_without_observation: u32,
    reads_discarded_by_cap: u32,
    observations: Vec<ScratchObservation>,
    used: usize,
    bytes_copied: u64,
}

impl LocusScratch {
    /// Refill from locus `i` of `block`. Everything is copied; nothing borrows the block.
    fn refill_from(&mut self, block: &LocusBlock, i: usize) {
        self.contig = block.contig[i];
        self.start = block.start[i];
        self.end = block.end[i];
        self.reads_without_observation = block.reads_without_observation[i];
        self.reads_discarded_by_cap = block.reads_discarded_by_cap[i];
        self.reference_bases.clear();
        self.reference_bases
            .extend_from_slice(block.reference_bases_of(i));
        let mut copied = block.reference_bases_of(i).len();
        self.used = 0;
        for j in block.observations_of(i) {
            if self.used == self.observations.len() {
                self.observations.push(ScratchObservation::default());
            }
            let slot = &mut self.observations[self.used];
            slot.bases.clear();
            slot.bases.extend_from_slice(block.bases_of(j));
            slot.complete = block.obs_witness_kind[j] == 0;
            slot.runs.clear();
            slot.runs.extend_from_slice(block.runs_of(j));
            slot.read_group = block.obs_read_group[j];
            slot.num_obs = block.obs_num_obs[j];
            slot.num_fwd = block.obs_num_fwd[j];
            slot.q_sum = block.obs_q_sum[j];
            slot.mapq_sum = block.obs_mapq_sum[j];
            slot.mapq_sum_sq = block.obs_mapq_sum_sq[j];
            slot.placed_left = block.obs_placed_left[j];
            slot.chain_ids.clear();
            slot.chain_ids.extend_from_slice(block.chain_ids_of(j));
            copied +=
                slot.bases.len() + slot.runs.len() * 4 + slot.chain_ids.len() * 8 + 41;
            self.used += 1;
        }
        self.bytes_copied += copied as u64;
    }
}

/// Arm C's consumer: **record-shaped, line for line the same as arm A's**, over a locus the
/// pre-pass refilled into its own buffer.
fn consume_block_via_scratch(
    accumulators: &mut Accumulators,
    block: &LocusBlock,
    scratch: &mut LocusScratch,
) {
    for i in 0..block.len() {
        scratch.refill_from(block, i);
        accumulators.groups.clear();
        let mut whole = (0u32, 0u32);
        accumulators.reads_without_observation += u64::from(scratch.reads_without_observation);
        if scratch.reads_discarded_by_cap > 0 {
            accumulators.loci_subsampled += 1;
        }
        for observation in &scratch.observations[..scratch.used] {
            reduce_observation(
                accumulators,
                &scratch.reference_bases,
                ObsRef {
                    bases: &observation.bases,
                    complete: observation.complete,
                    read_group: observation.read_group,
                    num_obs: observation.num_obs,
                },
                &mut whole,
            );
        }
        if accumulators.digest_on {
            digest_block_locus(accumulators, block, i);
        }
        let span = scratch.end - scratch.start + 1;
        accumulators.add_site(scratch.contig, scratch.start, span, whole);
        if accumulators.census_keeps(scratch.contig, scratch.start) {
            let (owned, bytes) = block.copy_out(i);
            accumulators.kept_bytes += bytes as u64;
            accumulators.kept.push(owned);
        }
    }
}

/// Arm Cb: arm C's consumer over slices **borrowed straight out of the block**, with no
/// per-locus refill. Not a candidate — this is the self-referential shape production
/// rejected — and measured only so that `C − Cb` prices the refill.
fn consume_block_borrowed(accumulators: &mut Accumulators, block: &LocusBlock) {
    for i in 0..block.len() {
        accumulators.groups.clear();
        let mut whole = (0u32, 0u32);
        accumulators.reads_without_observation +=
            u64::from(block.reads_without_observation[i]);
        if block.reads_discarded_by_cap[i] > 0 {
            accumulators.loci_subsampled += 1;
        }
        let reference_bases = block.reference_bases_of(i);
        for j in block.observations_of(i) {
            reduce_observation(
                accumulators,
                reference_bases,
                ObsRef {
                    bases: block.bases_of(j),
                    complete: block.obs_witness_kind[j] == 0,
                    read_group: block.obs_read_group[j],
                    num_obs: block.obs_num_obs[j],
                },
                &mut whole,
            );
        }
        if accumulators.digest_on {
            digest_block_locus(accumulators, block, i);
        }
        let span = block.end[i] - block.start[i] + 1;
        accumulators.add_site(block.contig[i], block.start[i], span, whole);
        if accumulators.census_keeps(block.contig[i], block.start[i]) {
            let (owned, bytes) = block.copy_out(i);
            accumulators.kept_bytes += bytes as u64;
            accumulators.kept.push(owned);
        }
    }
}

// ---------------------------------------------------------------------------------------
// The driver — the probe's walk, unchanged
// ---------------------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Arm {
    Records,
    Columns,
    Scratch,
    Borrowed,
}

impl Arm {
    fn from_env() -> Result<Self, String> {
        match std::env::var("PVC_SKETCH_ARM").as_deref() {
            Err(_) | Ok("A") => Ok(Arm::Records),
            Ok("B") => Ok(Arm::Columns),
            Ok("C") => Ok(Arm::Scratch),
            Ok("Cb") => Ok(Arm::Borrowed),
            Ok(other) => Err(format!("PVC_SKETCH_ARM must be A, B, C or Cb, not {other:?}")),
        }
    }

    fn label(self) -> &'static str {
        match self {
            Arm::Records => "A_records",
            Arm::Columns => "B_columns",
            Arm::Scratch => "C_scratch_view",
            Arm::Borrowed => "Cb_borrowed_view",
        }
    }

    fn wants_block(self) -> bool {
        self != Arm::Records
    }
}

fn contig_is_selected(name: &str, filter: Option<&str>) -> bool {
    filter.is_none_or(|wanted| name == wanted)
}

fn split_generic(
    item: Result<TypedRegion, TypedRegionError>,
    chunk_bp: Option<u64>,
) -> Vec<Result<TypedRegion, TypedRegionError>> {
    let Ok(region) = &item else {
        return vec![item];
    };
    let Some(chunk) = chunk_bp else {
        return vec![item];
    };
    if region.kind != RegionKind::Generic {
        return vec![item];
    }
    let mut out = Vec::new();
    let mut at = region.region.start.get();
    let end = region.region.end.get();
    while at <= end {
        let stop = (at + chunk - 1).min(end);
        out.push(Ok(TypedRegion {
            region: GenomeRegion {
                contig: region.region.contig,
                start: Position(at),
                end: Position(stop),
            },
            kind: RegionKind::Generic,
        }));
        at = stop + 1;
    }
    out
}

fn parse_env_u64(name: &str) -> Result<Option<u64>, String> {
    match std::env::var(name) {
        Err(_) => Ok(None),
        Ok(value) => value
            .parse::<u64>()
            .map(Some)
            .map_err(|_| format!("{name} must be a whole number, not {value:?}")),
    }
}

fn parse_env_flag(name: &str) -> bool {
    std::env::var(name).is_ok_and(|value| value == "1")
}

struct Run {
    arm: Arm,
    block_bytes: usize,
    max_loci: Option<u64>,
    whole_contig: bool,
    chunk_bp: Option<u64>,
    sample_one_in: u64,
    digest_on: bool,
    fit_on: bool,
}

#[allow(clippy::too_many_lines)]
fn walk(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: Option<&str>,
    run: Run,
    cache: &Arc<ReferenceInfoCache>,
) -> Result<(), Box<dyn std::error::Error>> {
    let reference_check = reference_check_from_env()?;
    let (info, verify) =
        read_reference_verifying_or_creating_fai(cache, fasta.to_path_buf(), reference_check)?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;
    // The same preparer the probe builds, so the two walks are the same pipeline.
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
    let reference = OpenReference::new(info);
    let sample =
        SampleReads::open_only_sample(bams, &reference, ReadFilterConfig::default(), true)?;

    #[allow(clippy::arc_with_non_send_sync)]
    let walk_reference = Arc::new(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
    let make_reference = {
        let fasta = fasta.to_path_buf();
        let contigs = contigs.clone();
        let index = index.clone();
        move || WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
    };
    let mut config = PileupGeneratorConfig::default();
    if let Some(cap) = parse_env_u64("PVC_PROBE_MAX_ACTIVE_READS")? {
        config.max_active_reads = cap as u32;
    }
    if let Some(span) = parse_env_u64("PVC_PROBE_MAX_RECORD_SPAN")? {
        config.max_record_span = span as u32;
    }
    let generator = PileupGenerator::new(walk_reference, make_reference, preparer, config)?;
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    let typed: Box<dyn Iterator<Item = Result<TypedRegion, TypedRegionError>>> =
        if run.whole_contig {
            let mut whole = Vec::new();
            for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
                if !contig_is_selected(&entry.name, contig_filter) {
                    continue;
                }
                whole.push(Ok(TypedRegion {
                    region: GenomeRegion {
                        contig: ContigId(index_of_contig as u32),
                        start: Position(1),
                        end: Position(entry.length),
                    },
                    kind: RegionKind::Generic,
                }));
            }
            Box::new(whole.into_iter())
        } else {
            let walk_config = TypedRegionConfig::default();
            let mut walks = Vec::new();
            for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
                if !contig_is_selected(&entry.name, contig_filter) {
                    continue;
                }
                walks.push(TypedRegionIterator::over_contig(
                    WindowedRefSeq::with_shared_index(
                        fasta.to_path_buf(),
                        contigs.clone(),
                        index.clone(),
                    ),
                    ContigId(index_of_contig as u32),
                    walk_config.clone(),
                )?);
            }
            Box::new(walks.into_iter().flatten())
        };
    let chunk_bp = run.chunk_bp;
    let regions = typed.flat_map(move |item| split_generic(item, chunk_bp));

    let mut accumulators = Accumulators::new(run.sample_one_in, run.digest_on);
    let mut stream = SampleLocusObservationsIterator::new(regions, sample, generators);
    let mut blocks = 0u64;
    let mut block_payload_high_water = 0usize;
    let mut scratch = LocusScratch::default();
    let started = Instant::now();

    if run.arm.wants_block() {
        stream.install_block_sink(LocusSink::columns(run.block_bytes));
        loop {
            let more = stream.fill_block()?;
            let Some(sink) = stream.block_sink_mut() else {
                break;
            };
            let block = sink.block();
            if !block.is_empty() {
                blocks += 1;
                block_payload_high_water = block_payload_high_water.max(block.payload_bytes());
                match run.arm {
                    Arm::Columns => consume_block_columns(&mut accumulators, block),
                    Arm::Scratch => {
                        consume_block_via_scratch(&mut accumulators, block, &mut scratch)
                    }
                    Arm::Borrowed => consume_block_borrowed(&mut accumulators, block),
                    Arm::Records => unreachable!("the record arm does not fill blocks"),
                }
            }
            // PANIC-FREE: the sink was just read through the same `Option`.
            stream
                .block_sink_mut()
                .expect("the sink is installed")
                .block_mut()
                .clear();
            if !more {
                break;
            }
            if run.max_loci.is_some_and(|max| accumulators.total_loci >= max) {
                break;
            }
        }
    } else {
        for locus in &mut stream {
            let locus = locus?;
            consume_record(&mut accumulators, locus);
            if run.max_loci.is_some_and(|max| accumulators.total_loci >= max) {
                break;
            }
        }
    }
    let seconds = started.elapsed().as_secs_f64();

    let counts = stream.counts();
    let generic = match stream.generators().generic_counts() {
        Some(pop_var_caller::ng::locus_generation::GeneratorCounts::Pileup(generic)) => *generic,
        _ => return Err("the generic slot reported no counts".into()),
    };
    let PileupGeneratorCounts {
        reads_admitted,
        records_outside_region,
        mate_overlap_positions,
        column_depth_high_water,
        ..
    } = generic;

    println!("arm={}", run.arm.label());
    println!("reference_check={}", reference_check_label(reference_check));
    println!("seconds={seconds:.3}");
    println!("block_budget_bytes={}", run.block_bytes);
    println!("blocks={blocks}");
    println!("block_payload_high_water={block_payload_high_water}");
    println!("regions_in={}", counts.regions_in);
    println!("regions_handled={}", counts.regions_handled);
    println!("loci_emitted={}", counts.loci_emitted);
    println!("reads_admitted={reads_admitted}");
    println!("records_outside_region={records_outside_region}");
    println!("mate_overlap_positions={mate_overlap_positions}");
    println!("column_depth_high_water={column_depth_high_water}");
    println!("sample_one_in={}", run.sample_one_in);
    println!("census_kept={}", accumulators.kept.len());
    println!("census_bytes_copied={}", accumulators.kept_bytes);
    println!("scratch_bytes_copied={}", scratch.bytes_copied);
    println!("{}", accumulators.summary());
    if run.fit_on {
        print!("{}", accumulators.fit());
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: {} <reference.fa> <sample.bam|cram> [contig]",
            args.first().map_or("ng_prepass_sketch", |s| s.as_str())
        );
        return ExitCode::FAILURE;
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig = args.get(3).map(String::as_str);

    let run = match (|| -> Result<Run, String> {
        Ok(Run {
            arm: Arm::from_env()?,
            block_bytes: parse_env_u64("PVC_SKETCH_BLOCK_KB")?.unwrap_or(256) as usize * 1024,
            max_loci: parse_env_u64("PVC_PROBE_MAX_LOCI")?,
            whole_contig: parse_env_flag("PVC_PROBE_WHOLE_CONTIG"),
            chunk_bp: parse_env_u64("PVC_GENERIC_REGION_CHUNK_BP")?,
            sample_one_in: parse_env_u64("PVC_SKETCH_SAMPLE_ONE_IN")?.unwrap_or(29_000),
            digest_on: parse_env_flag("PVC_SKETCH_DIGEST"),
            fit_on: !parse_env_flag("PVC_SKETCH_NO_FIT"),
        })
    })() {
        Ok(run) => run,
        Err(message) => {
            eprintln!("{message}");
            return ExitCode::FAILURE;
        }
    };

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let cache = Arc::new(ReferenceInfoCache::new());
    match walk(&fasta, &[bam], contig, run, &cache) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
