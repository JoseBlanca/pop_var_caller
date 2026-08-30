//! The STR locus generator — one microsatellite tract → one locus.
//!
//! The first [`LocusGenerator`](super::LocusGenerator): it consumes an `SsrSegment`,
//! fetches the reads over the tract, aligns each to read off its repeat, and tallies the
//! answers into a [`SampleLocusObservations`](super::SampleLocusObservations). A port of
//! production `src/ssr/pileup/`, adapted at two seams (the split of coordinates from bases,
//! and a widened admission gate). See `doc/devel/ng/spec/locus_generation_ssr.md` (design)
//! and `doc/devel/ng/arch/locus_generation_ssr.md` (types & interfaces).
//!
//! This module lands across the STR generator plan: the config, counts, cap constant, the
//! working input ([`SsrLocus`] + its margin fetch), the reservoir cap, the per-locus read
//! fetch, the read-region CIGAR mapping, the per-read classify pipeline (delimited by a chosen
//! [`RepeatDelimiter`] — algorithm 4, the unit-slip aligner, by default; algorithm 3, the
//! flat-gap production-parity port, for the bake-off and the parity oracle), the tally, and — the
//! public surface — [`SsrGenerator`], the [`LocusGenerator`](super::LocusGenerator) that turns one
//! `SsrSegment` into one locus.

use super::{
    GeneratorCounts, LocusGenerationError, LocusGenerator, LocusKind, ReadWitness,
    SampleLocusObservations, SsrDetail,
};
#[cfg(test)]
use crate::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use crate::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner;
use crate::ng::alignment::{
    BestPathAligner, PerQualityEmission, RepeatContext, RepeatSpan, StutterModel,
};
use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::input::cursor::{CursorCounts, CursorError};
use crate::ng::read::input::sample_cursor::SampleCursor;
use crate::ng::read::input::{SampleIdentity, SampleReads};
use crate::ng::ref_seq::{ContigTable, EvictableRefSeq, RawRefSeq, RefSeq, RefSeqError};
use crate::ng::region_typing::segment_criteria::SsrSegment;
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};

/// An STR locus ready to align against: the segment plus the reference bases the aligner
/// aligns the reads to.
///
/// The ng counterpart of production's `Locus`, **split** so the coordinates come from region
/// typing (`SsrSegment`) and the bases are fetched here (spec §2). Fetching them makes the
/// port an *adaptation, not a lift* — the most likely place for the port to go subtly wrong,
/// which is why the flank lengths are always **measured** from the clamped span, never
/// assumed to be `flank_bp`.
#[derive(Debug, Clone, PartialEq)]
pub struct SsrLocus {
    /// The tract's coordinates, motif and purity — from region typing.
    pub segment: SsrSegment,
    /// The tract plus its query margin, canonical `{A,C,G,T,N}` bases, **clamped at contig
    /// ends** — so this may be shorter than `2 * flank_bp + tract`, and each flank must be
    /// measured, never assumed (spec §2).
    pub tract_with_margin_bases: Box<[u8]>,
    /// 1-based position of `tract_with_margin_bases[0]`.
    pub margin_start: Position,
}

impl SsrLocus {
    /// Fetch the tract ± `flank_bp` for `segment` into an [`SsrLocus`], clamped at the contig
    /// ends, using the reused `buffer` as scratch (spec §2).
    ///
    /// `contig` is the tract's contig (from the region handed to `begin_segment`); its length
    /// is read from `reference`'s contig table to clamp the right margin, and `fetch_into`
    /// validates the final window. The bases are canonical — the aligner compares against
    /// them, and production upper-cases the same way.
    pub fn fetch<R: RefSeq + ContigTable>(
        reference: &R,
        contig: ContigId,
        segment: SsrSegment,
        flank_bp: Bp,
        buffer: &mut Vec<u8>,
    ) -> Result<Self, RefSeqError> {
        let flank = flank_bp.get();
        let tract_start = segment.start(); // 1-based inclusive
        let tract_end = segment.end();
        // The left margin cannot run below base 1; the right cannot run past the contig.
        let margin_start = tract_start.saturating_sub(flank).max(1);
        let margin_end = match reference.contigs().entries.get(contig.get() as usize) {
            // Clamp the right margin to the contig — but only when the tract actually fits
            // inside it. A tract reaching past the contig end is a broken segment: leaving
            // the window unclamped makes `fetch_into` reject it as out-of-bounds, rather than
            // fetching a window that is missing the tract (which would underflow the derived
            // right flank). An unknown contig (`None`) is likewise left for `fetch_into`.
            Some(entry) if tract_end <= entry.length => (tract_end + flank).min(entry.length),
            _ => tract_end + flank,
        };
        let length = (margin_end + 1).saturating_sub(margin_start);

        reference.fetch_into(contig, margin_start, length, buffer)?;
        Ok(Self {
            segment,
            tract_with_margin_bases: buffer.as_slice().into(),
            margin_start: Position(margin_start),
        })
    }

    /// The left flank's length in bases — **measured** from the clamp (`tract_start −
    /// margin_start`), so a tract near the contig start reports a short flank, not `flank_bp`.
    pub fn left_flank_len(&self) -> usize {
        (self.segment.start() - self.margin_start.get()) as usize
    }

    /// The right flank's length in bases — measured as what remains after the left flank and
    /// the tract, so a tract near the contig end reports a short flank.
    pub fn right_flank_len(&self) -> usize {
        self.tract_with_margin_bases.len()
            - self.left_flank_len()
            - self.segment.tract_len() as usize
    }
}

/// The STR generator's knobs — owned and taken at construction
/// ([shared config discipline](super); spec §4). Each generator owns its own knobs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrGeneratorConfig {
    /// The flank fetched either side of the tract — the aligner's anchor and the read
    /// query margin. Must be `<= bundle_threshold`, the radius region typing guarantees
    /// repeat-free (checked by [`Self::check_flank_within`]); equal by default (spec §4).
    pub flank_bp: Bp,
    /// Reads kept per locus, reservoir-sampled. `None` = no cap (spec §4).
    pub max_reads_per_locus: Option<u32>,
}

/// ng's own per-locus read cap — **not** production's `MAX_READS_PER_LOCUS`, so the two can
/// diverge. Starts at 1000 (matching production) but is never-measured and soft: to be set
/// by experiment (spec §4).
pub const DEFAULT_SSR_MAX_READS_PER_LOCUS: u32 = 1000;

impl Default for SsrGeneratorConfig {
    /// The flank equals the bundle threshold's default (spec §4, "equal by default"), and
    /// the cap starts at [`DEFAULT_SSR_MAX_READS_PER_LOCUS`].
    fn default() -> Self {
        Self {
            flank_bp: Bp(crate::ng::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            max_reads_per_locus: Some(DEFAULT_SSR_MAX_READS_PER_LOCUS),
        }
    }
}

impl SsrGeneratorConfig {
    /// Check the cross-config invariant `flank_bp <= bundle_threshold`.
    ///
    /// A wider flank than the radius region typing guarantees repeat-free would let the read
    /// query hit a neighbouring repeat, leaving the aligner's anchor no longer clean (spec
    /// §4). It is a relation between two configs, so no newtype can hold it — the generator's
    /// constructor calls this. `bundle_threshold` is [`SsrSegmentCriteria::bundle_threshold`]
    /// ([`crate::ng::region_typing::segment_criteria`]).
    pub fn check_flank_within(&self, bundle_threshold: Bp) -> Result<(), SsrGeneratorConfigError> {
        if self.flank_bp.get() > bundle_threshold.get() {
            return Err(SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: self.flank_bp.get(),
                bundle_threshold: bundle_threshold.get(),
            });
        }
        Ok(())
    }
}

/// A malformed STR generator configuration. `#[non_exhaustive]`; raised at construction so a
/// bad knob combination never reaches a locus.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SsrGeneratorConfigError {
    /// The flank is wider than region typing's clean-flank guarantee — the read query could
    /// reach a neighbouring repeat and spoil the aligner's anchor (spec §4).
    #[error(
        "flank_bp ({flank_bp}) exceeds bundle_threshold ({bundle_threshold}): the read query \
         would reach past the repeat-free radius region typing guarantees"
    )]
    FlankExceedsBundleThreshold {
        flank_bp: u64,
        bundle_threshold: u64,
    },
}

/// Run-level STR counts, reported alongside the shared
/// [`LocusCounts`](super::LocusCounts). The locus records *that* reads yielded nothing
/// ([`reads_without_observation`](super::SampleLocusObservations::reads_without_observation));
/// **why** is this generator's to report, because the reasons are specific to how it reads a
/// tract and mean nothing to a pileup (spec §4).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SsrGeneratorCounts {
    /// Reads fetched over the tract-plus-margin query span.
    pub reads_fetched: u64,
    /// Reads the reservoir cap dropped.
    pub reads_discarded_by_cap: u64,
    /// Observations pinning the tract length (both borders seen).
    pub observations_complete: u64,
    /// Observations proving a lower bound (one border, ran off the other end).
    pub observations_partial: u64,
    /// Reads that reached the aligner and anchored no border (`read_preparation_ssr.md` §4).
    pub no_border_anchored: u64,
    /// Reads dropped for low base quality.
    pub low_quality: u64,
    /// Reads whose allele stayed flank-truncated even after widening.
    pub window_truncated: u64,
    /// Reads that anchored a flank but crossed **no** tract position — they overlap the
    /// locus *window* and lie outside the repeat the locus is. Logged rather than merely
    /// dropped, because it is a large population and the number is the evidence that the
    /// fetch window is wider than the locus: **6,704 reads against 7,085 real partials** on
    /// tomato chr01 of `SRR7279503` when this counter was added.
    pub outside_tract: u64,
}

impl SsrGeneratorCounts {
    /// Reads that reached the aligner and yielded nothing, **by every reason there is**.
    ///
    /// # It lives here because the sum is what drifts
    ///
    /// Two tools print this total and each summed the reasons itself. When C0 added a
    /// fourth, one of them was updated and the other was not — and it under-reported by
    /// 6,704 reads of ~9,265 on one tomato chromosome without a warning, because adding a
    /// field to a struct does not break a `+` chain. A reviewer confirmed the shape rather
    /// than the instance: adding a *fifth* reason and changing nothing else left
    /// `clippy --lib --examples --all-features -- -D warnings` clean, with both tools now
    /// silently short (Milestone C review, F5).
    ///
    /// **The destructure is the guard, and it has to name every field with no `..`.** A
    /// reason added to the struct then fails to compile *here*, once, with
    /// `error[E0027]: pattern does not mention field`, which is the question being asked:
    /// does this new reason mean the read yielded no observation? The observation counters
    /// are named and discarded for the same reason — so a field added to *them* also stops
    /// here rather than being silently swept into a no-observation total.
    pub fn reads_without_observation(&self) -> u64 {
        let Self {
            reads_fetched: _,
            reads_discarded_by_cap: _,
            observations_complete: _,
            observations_partial: _,
            no_border_anchored,
            low_quality,
            window_truncated,
            outside_tract,
        } = self;
        no_border_anchored + low_quality + window_truncated + outside_tract
    }
}

// ---------------------------------------------------------------------
// The per-locus read cap — a faithful port of production's reservoir sampler
// (src/ssr/pileup/fetch_reads.rs), keyed to ng's own seed and cap constant.
// ---------------------------------------------------------------------

/// A tiny deterministic PRNG (SplitMix64) — seeded per locus so the depth-cap subsample is
/// reproducible and thread-count-invariant, with no external RNG whose stream could shift.
/// Ported verbatim from production; the constants are load-bearing for byte parity.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7615);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}

/// Deterministic per-locus reservoir seed from a contig name and a **0-based** tract start:
/// FNV-1a over the name bytes, folded with the start. Ported verbatim from production so the
/// kept set matches byte-for-byte (the parity oracle depends on it).
///
/// **The trap (spec §4):** the seed is over the contig **name** and the **0-based** start.
/// ng speaks `ContigId` and 1-based positions, so seeding from the id or the 1-based start
/// silently produces a *different* kept set — deterministic, so the parity test fails looking
/// like an aligner bug. Callers seed through [`seed_for_segment`], which does the conversion.
fn locus_seed(chrom: &str, start_0based: u32) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    for &b in chrom.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= start_0based as u64;
    h.wrapping_mul(FNV_PRIME)
}

/// The reservoir seed for an STR segment — the one place the seed trap is discharged: the
/// contig **name** and the tract's **0-based** start (`start() - 1`, since `SsrSegment` is
/// 1-based).
///
/// The `as u32` matches production's seed domain (its `Locus.start` is `u32`); per-contig
/// positions are far below `2^32` (the largest chromosome is ~250 Mb), so the cast never
/// truncates in practice and parity holds. `start() - 1` cannot underflow: `SsrSegment::new`
/// enforces `1 <= start`.
pub fn seed_for_segment(segment: &SsrSegment) -> u64 {
    locus_seed(segment.chrom(), (segment.start() - 1) as u32)
}

/// Reservoir sampler (Algorithm R) — an effectively-uniform sample of up to `capacity` items
/// from a stream of unknown length, in one pass with `O(capacity)` memory. Ported verbatim
/// from production: the eviction index `next_u64() % seen` carries a modulo bias bounded by
/// `seen / 2^64` (negligible at any real depth), accepted deliberately because the draw is
/// deterministic and thread-count-invariant — an unbiased reduction would change the kept set
/// and break that. The caller must `offer` admitted reads in a fixed total order —
/// `SampleReads`' merge order (spec §4).
pub struct Reservoir<T> {
    capacity: usize,
    held: Vec<T>,
    /// Admitted items offered so far (the `i` of Algorithm R).
    seen: u64,
    rng: SplitMix64,
}

impl<T> Reservoir<T> {
    /// A reservoir of `capacity` items seeded by `seed` (from [`seed_for_segment`]).
    pub fn new(capacity: usize, seed: u64) -> Self {
        Self {
            capacity,
            held: Vec::with_capacity(capacity),
            seen: 0,
            rng: SplitMix64::new(seed),
        }
    }

    /// Offer one admitted item. Keeps the first `capacity`; for the `i`-th item
    /// (`i > capacity`) keeps it with probability `capacity / i`, evicting one held item
    /// uniformly at random if kept.
    pub fn offer(&mut self, item: T) {
        self.seen += 1;
        if self.held.len() < self.capacity {
            self.held.push(item);
        } else {
            let j = (self.rng.next_u64() % self.seen) as usize;
            if j < self.capacity {
                self.held[j] = item;
            }
        }
    }

    /// The admitted depth — total items offered (the reservoir sees only admitted reads).
    pub fn seen(&self) -> u64 {
        self.seen
    }

    /// Consume the reservoir, yielding the sampled items (≤ `capacity`).
    pub fn into_held(self) -> Vec<T> {
        self.held
    }
}

// ---------------------------------------------------------------------
// The per-locus read fetch (D1): fetch over the tract+margin span, depth-cap.
// ---------------------------------------------------------------------

/// The reads kept for one locus after depth-capping, and how many were fetched — the caller
/// records `reads_discarded_by_cap = fetched - kept.len()`.
pub struct CappedReads {
    pub kept: Vec<AlignedRead>,
    pub fetched: u64,
}

/// Fetch the reads over a locus's query span (tract + margin) and depth-cap them with the
/// reservoir, seeded per locus (spec §2, §4).
///
/// Admission is **relevance** — overlap with `query_span`, which `SampleReads` already applies
/// — **not** spanning: production's `reaches_locus` gate is deliberately not ported, because a
/// read that runs off mid-tract is exactly what a partial observation is made of (spec §2). The
/// cap sits between fetch and align; `None` keeps every read.
///
/// **Pointed at a cursor, not handed a query — D3.** The cursor covers `query_span`'s
/// chromosome and stays open across every locus on it, so consecutive loci that share reads
/// decode them once (`spec/alignment_cursor.md` §1). Making the cursor and moving it to the
/// right chromosome is the caller's job; this only points it at the span.
pub fn fetch_capped_reads<R: RawRefSeq>(
    reads: &mut SampleCursor<R>,
    query_span: GenomeRegion,
    seed: u64,
    max_reads_per_locus: Option<u32>,
) -> Result<CappedReads, LocusGenerationError> {
    // **The region these failures are attributed to is the span queried**, this
    // function knowing no other: it is handed a tract-plus-margin span, not the
    // segment it came from. Every other attachment in this module names the
    // segment's own region (`LocusGenerationError`'s doc).
    reads
        .move_to_region(query_span)
        .map_err(|source| LocusGenerationError::OpenReadQuery {
            region: query_span,
            source: source.into(),
        })?;
    let read_failed = |source: CursorError| LocusGenerationError::Reads {
        region: query_span,
        source: source.into(),
    };
    let mut fetched = 0u64;
    let kept = match max_reads_per_locus {
        Some(cap) => {
            let mut reservoir = Reservoir::new(cap as usize, seed);
            while let Some(read) = reads.next_read() {
                reservoir.offer(read.map_err(read_failed)?);
                fetched += 1;
            }
            reservoir.into_held()
        }
        None => {
            let mut all = Vec::new();
            while let Some(read) = reads.next_read() {
                all.push(read.map_err(read_failed)?);
                fetched += 1;
            }
            all
        }
    };
    Ok(CappedReads { kept, fetched })
}

// ---------------------------------------------------------------------
// D2a — the read-region mapping: CIGAR → the read-coordinate slice covering the locus
// window, plus the region quality gate. Ported from production
// (`src/ssr/pileup/footprint.rs`, `alignment.rs`), adapted to take window coordinates
// directly and a plain scratch buffer. Consumed by the classify pipeline (D2b).
// ---------------------------------------------------------------------
mod read_region {
    use crate::pileup::walker::CigarOp;
    use std::ops::Range;

    /// The lower-quartile base-quality floor a delimited tract must clear (production's
    /// `MIN_REGION_Q1`).
    pub(super) const MIN_REGION_Q1: u8 = 15;

    /// A read's reference footprint from its CIGAR — the aligned span (0-based, half-open)
    /// and the soft-clip on each end (the optimistic long-allele reach).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) struct Footprint {
        pub(super) ref_start: u32,
        pub(super) ref_end: u32,
        pub(super) leading_clip: u32,
        pub(super) trailing_clip: u32,
    }

    /// First soft-clip length at an end (skipping a hard-clip), or `0` if the end is aligned.
    fn end_soft_clip(op: &CigarOp) -> Option<u32> {
        match op {
            CigarOp::HardClip(_) => None,
            CigarOp::SoftClip(n) => Some(*n),
            _ => Some(0),
        }
    }

    /// The read's reference footprint from its CIGAR and 1-based mapping position.
    pub(super) fn read_footprint(cigar: &[CigarOp], pos: u64) -> Footprint {
        let ref_start = pos.saturating_sub(1) as u32;
        let ref_span: u32 = cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n)
                | CigarOp::Deletion(n)
                | CigarOp::Skip(n)
                | CigarOp::SeqMatch(n)
                | CigarOp::SeqMismatch(n) => *n,
                _ => 0,
            })
            .sum();
        Footprint {
            ref_start,
            ref_end: ref_start + ref_span,
            leading_clip: cigar.iter().find_map(end_soft_clip).unwrap_or(0),
            trailing_clip: cigar.iter().rev().find_map(end_soft_clip).unwrap_or(0),
        }
    }

    /// Read coordinate of a 0-based reference position `target` inside the aligned span. A
    /// dual-cursor CIGAR walk; a `target` inside a deletion maps to the read position the
    /// deletion sits at.
    fn ref_to_read(cigar: &[CigarOp], ref_start: u32, leading_clip: u32, target: u32) -> usize {
        let mut ref_cur = ref_start;
        let mut read_cur = leading_clip as usize;
        for op in cigar {
            match op {
                CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                    if target < ref_cur + n {
                        return read_cur + (target - ref_cur) as usize;
                    }
                    ref_cur += n;
                    read_cur += *n as usize;
                }
                CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                    if target < ref_cur + n {
                        return read_cur;
                    }
                    ref_cur += n;
                }
                CigarOp::Insertion(n) => read_cur += *n as usize,
                CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
            }
        }
        read_cur
    }

    /// The read-coordinate span covering the locus window `[window_start, window_start +
    /// window_len)` (0-based reference), where the whole soft-clip is included on a side the
    /// window opens past the aligned span — that is where a long allele's extra tract lives,
    /// and the aligner realigns within the grabbed bases (spec §2).
    pub(super) fn extract_region(
        cigar: &[CigarOp],
        fp: Footprint,
        read_len: usize,
        window_start: u32,
        window_len: u32,
    ) -> Range<usize> {
        let w_end = window_start + window_len;
        let r_start = if window_start < fp.ref_start {
            0 // window opens left of the alignment → take the full leading clip
        } else {
            ref_to_read(cigar, fp.ref_start, fp.leading_clip, window_start)
        };
        let r_end = if w_end > fp.ref_end {
            read_len // window closes right of the alignment → take the full trailing clip
        } else {
            ref_to_read(cigar, fp.ref_start, fp.leading_clip, w_end)
        };
        // Clamp into `[0, read_len]` — defense-in-depth against a length-inconsistent CIGAR.
        let r_start = r_start.min(read_len);
        r_start..r_end.min(read_len).max(r_start)
    }

    /// Whether a delimited tract looks truncated by the extraction window rather than by the
    /// read genuinely ending — the long-allele recovery trigger. A side is window-bounded when
    /// [`extract_region`] did not reach the read edge there; the tract is suspicious when a
    /// window-bounded side carries fewer flank bytes than the locus declares. `tract` is
    /// region-relative.
    pub(super) fn flank_truncated(
        region: &Range<usize>,
        tract: &Range<usize>,
        read_len: usize,
        left_flank_len: usize,
        right_flank_len: usize,
    ) -> bool {
        let region_len = region.end - region.start;
        let left_flank_bytes = tract.start;
        let right_flank_bytes = region_len - tract.end;
        let left_window_bounded = region.start > 0;
        let right_window_bounded = region.end < read_len;
        (left_window_bounded && left_flank_bytes < left_flank_len)
            || (right_window_bounded && right_flank_bytes < right_flank_len)
    }

    /// Widen a read-coordinate region by one full reference flank each side, clamped to the
    /// read — in read, not reference, coordinates, sidestepping the unreliable CIGAR of a
    /// mis-aligned long-allele read.
    pub(super) fn widen_region(
        region: Range<usize>,
        read_len: usize,
        left_flank_len: usize,
        right_flank_len: usize,
    ) -> Range<usize> {
        let start = region.start.saturating_sub(left_flank_len);
        let end = (region.end + right_flank_len).min(read_len);
        start..end
    }

    /// Whether the tract's base qualities clear the gate: the nearest-rank lower quartile (the
    /// element at sorted index `⌊(n-1)/4⌋`) is at least `threshold`. `buffer` is reused
    /// scratch. An empty tract passes vacuously.
    pub(super) fn passes_quality_gate(quals: &[u8], threshold: u8, buffer: &mut Vec<u8>) -> bool {
        if quals.is_empty() {
            return true;
        }
        buffer.clear();
        buffer.extend_from_slice(quals);
        let k = (buffer.len() - 1) / 4;
        let (_, q1, _) = buffer.select_nth_unstable(k);
        *q1 >= threshold
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// A window fully inside an all-Match read maps to the same read offsets as reference
        /// offsets (shifted by the mapping position).
        #[test]
        fn extract_region_maps_a_window_inside_an_all_match_read() {
            // read at pos 11 (0-based 10), 30M; window [15, 25) (0-based).
            let cigar = vec![CigarOp::Match(30)];
            let fp = read_footprint(&cigar, 11);
            let region = extract_region(&cigar, fp, 30, 15, 10);
            // ref 15 → read 5 (15 - 10), ref 25 → read 15.
            assert_eq!(region, 5..15);
        }

        /// An internal deletion shifts the window's right edge left in read coordinates.
        #[test]
        fn extract_region_maps_across_an_internal_deletion() {
            // pos 11 (0-based 10): 10M 4D 16M — ref span 30, read len 26. Window [15, 28).
            let cigar = vec![CigarOp::Match(10), CigarOp::Deletion(4), CigarOp::Match(16)];
            let fp = read_footprint(&cigar, 11);
            assert_eq!(fp.ref_end, 40);
            let region = extract_region(&cigar, fp, 26, 15, 13);
            // ref 15 → read 5; ref 28 is past the deletion (ref 20..24 deleted), read = 10 + (28-24) = 14.
            assert_eq!(region, 5..14);
        }

        /// A window opening left of the alignment takes the full leading soft-clip.
        #[test]
        fn extract_region_takes_the_leading_softclip_when_the_window_opens_left() {
            // pos 21 (0-based 20): 5S 20M. Window [18, 30) opens left of ref_start 20.
            let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(20)];
            let fp = read_footprint(&cigar, 21);
            assert_eq!(fp.leading_clip, 5);
            let region = extract_region(&cigar, fp, 25, 18, 12);
            assert_eq!(region.start, 0, "the full leading clip is grabbed");
        }

        #[test]
        fn flank_truncated_flags_a_window_bounded_short_flank_only() {
            // region [3, 40) within a 50-base read (both sides window-bounded), tract [5, 32)
            // → left flank 5, right flank region_len(37) - 32 = 5.
            assert!(
                flank_truncated(&(3..40), &(5..32), 50, 6, 6),
                "a flank of 5 below the declared 6 is truncated"
            );
            assert!(
                !flank_truncated(&(3..40), &(5..32), 50, 5, 5),
                "flanks exactly the declared length are not truncated"
            );
            // A region reaching both read edges is never flagged (genuine allele ≥ read length).
            assert!(!flank_truncated(&(0..50), &(0..45), 50, 6, 6));
        }

        #[test]
        fn widen_region_extends_a_flank_each_side_clamped_to_the_read() {
            assert_eq!(widen_region(20..40, 100, 10, 10), 10..50);
            // Clamp at both ends.
            assert_eq!(widen_region(5..95, 100, 10, 10), 0..100);
        }

        #[test]
        fn the_quality_gate_keys_on_the_nearest_rank_lower_quartile() {
            let mut buffer = Vec::new();
            assert!(
                passes_quality_gate(&[], MIN_REGION_Q1, &mut buffer),
                "empty passes"
            );
            // len 1 → k = 0 (the lowest); 15 passes, 14 fails.
            assert!(passes_quality_gate(&[15], MIN_REGION_Q1, &mut buffer));
            assert!(!passes_quality_gate(&[14], MIN_REGION_Q1, &mut buffer));
            // len 8 → k = ⌊7/4⌋ = 1 (the 2nd-lowest): one base below the floor still passes.
            assert!(passes_quality_gate(
                &[10, 20, 20, 20, 20, 20, 20, 20],
                MIN_REGION_Q1,
                &mut buffer
            ));
            assert!(!passes_quality_gate(
                &[10, 10, 20, 20, 20, 20, 20, 20],
                MIN_REGION_Q1,
                &mut buffer
            ));
        }
    }
}

// ---------------------------------------------------------------------
// D2b — classify one read against a locus: delimit (the chosen `RepeatDelimiter` — algorithm 3 or
// 4), recover a long allele by widening, gate the tract quality, and map the result
// to an observation or a no-observation reason. The per-read pipeline ported from
// production's `classify_read`, over the new delimiter. Consumed by the tally + generator.
// ---------------------------------------------------------------------
mod classify {
    use super::read_region::{
        MIN_REGION_Q1, extract_region, flank_truncated, passes_quality_gate, read_footprint,
        widen_region,
    };
    use super::{RepeatDelimiter, SsrLocus};
    use crate::ng::alignment::{
        ReadBases, RepeatContext, RepeatGeometry, RepeatSpan, StutterModel,
    };
    use crate::ng::locus_generation::{LocusLen, ReadWitness};
    use crate::ng::read::aligned_read::AlignedRead;
    use crate::ng::types::Bp;
    use crate::pileup::walker::CigarOp;
    use std::ops::Range;

    /// Why a read yielded no usable observation — the tally increments the matching
    /// `SsrGeneratorCounts` reason (and `reads_without_observation`).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) enum NoObservationReason {
        /// Neither flank anchored (the read lies wholly inside the repeat), or a malformed read.
        NoBorderAnchored,
        /// The delimited tract failed the base-quality gate.
        LowQuality,
        /// The allele stayed flank-truncated even after widening.
        WindowTruncated,
        /// **The read is not in this locus.** It anchored a flank and crossed no tract
        /// position at all, so the "lower bound" it would supply is *at least zero* — no
        /// evidence about the repeat's length, on a read that never entered the repeat.
        ///
        /// It arises because the fetch queries the tract **plus its margin**, so a read
        /// overlapping only the flank is delimited too; the aligner then reports
        /// [`RepeatSpan::FromLeft`](crate::ng::alignment::RepeatSpan::FromLeft) or
        /// `FromRight` with an empty span, the rejected side having fallen back to the
        /// read's own edge. Measured on tomato chr01 of `SRR7279503`: **6,704 such reads
        /// against 7,085 genuine partials**, at a median window overlap of 16 bases against
        /// a 30-base flank.
        ///
        /// Those bases belong to the SNP/indel path, which analyses them; the STR path
        /// discards them and counts them (owner, 2026-07-30). Before that decision they
        /// became observations with **empty bases** and a witness covering zero locus
        /// positions — 3,180 rows of the STR dump — contributing to `num_obs` while adding
        /// nothing to depth anywhere along the locus.
        OutsideTract,
    }

    /// What one read contributes to a locus: an observed tract (complete or partial), or a
    /// reason it observed nothing.
    #[derive(Debug, Clone, PartialEq)]
    pub(super) enum Classified {
        Observed {
            bases: Box<[u8]>,
            read_witness: ReadWitness,
            /// The read's base-quality error mass over the tract, `Σ ln(P_err)`, in log-error
            /// space (freebayes' `q_sum` convention). The **BQ** support moment the spec names
            /// (spec §3, "strand/BQ/MAPQ moments") — computed here because the tract base
            /// qualities are already sliced, a free by-product; the tally folds it into
            /// [`SequenceObservation::q_sum`](crate::ng::locus_generation::SequenceObservation). MAPQ
            /// is carried separately, off the [`AlignedRead`], in the tally. **Soft** — filled,
            /// unconsumed today.
            q_sum: f64,
        },
        NoObservation(NoObservationReason),
    }

    /// Classify one read against `locus` using `aligner` (the chosen [`RepeatDelimiter`]). `stutter`
    /// feeds the context — algorithm 3 (flat-gap) ignores it, algorithm 4 (unit-slip) prices its
    /// slips from it; `align_scratch` (the aligner's own scratch) and `qual_buffer` are reused
    /// across reads.
    ///
    /// The pipeline mirrors production: extract the read's slice over the locus window, delimit
    /// it, and — if a *complete* tract looks truncated by the window (a mapper-collapsed long
    /// allele) — widen the slice and re-delimit, giving up as `WindowTruncated` if it stays
    /// truncated. A complete tract is quality-gated; **partials are kept as lower bounds**
    /// without the gate (spec §3, the new behaviour production discards).
    pub(super) fn classify_read<A: RepeatDelimiter>(
        read: &AlignedRead,
        locus: &SsrLocus,
        aligner: &A,
        stutter: &StutterModel,
        align_scratch: &mut A::Scratch,
        qual_buffer: &mut Vec<u8>,
    ) -> Classified {
        // A malformed record (qual length ≠ seq length) cannot be sliced safely.
        if read.qual.len() != read.seq.len() {
            return Classified::NoObservation(NoObservationReason::NoBorderAnchored);
        }
        let read_len = read.seq.len();
        let reference = &*locus.tract_with_margin_bases;
        let left = locus.left_flank_len();
        let right = locus.right_flank_len();
        let geometry = RepeatGeometry {
            left_flank_len: Bp(left as u64),
            right_flank_len: Bp(right as u64),
            motif: locus.segment.motif(),
        };
        let context = RepeatContext {
            geometry: &geometry,
            stutter,
        };

        let fp = read_footprint(&read.cigar, read.pos);
        let window_start = (locus.margin_start.get() - 1) as u32; // 0-based
        let window_len = reference.len() as u32;
        let region = extract_region(&read.cigar, fp, read_len, window_start, window_len);

        // **A read that copies the reference across the whole window has nothing to align.**
        // Its repeat count is the reference's, and the tract sits where the reference's tract
        // sits — so the delimiter would spend a full dynamic-programming pass rediscovering a
        // span this function already knows. On tomato at three reads a position that is 4 reads
        // in every 9 reaching the aligner (57,507 of 132,069 over 24 spans of chromosome 1).
        //
        // The base-quality gate below still runs: skipping the aligner skips *where the tract
        // is*, not *whether the read is good enough to be counted*, which is a separate question
        // and is what turns the span into an observation.
        if let Some(tract) = tract_needing_no_alignment(read, &region, locus, reference) {
            return complete_or_low_quality(read, &region, &tract, qual_buffer);
        }

        let span = delimit(aligner, read, &region, reference, context, align_scratch);
        match span {
            RepeatSpan::Between(tract)
                if flank_truncated(&region, &to_usize(&tract), read_len, left, right) =>
            {
                // A complete tract whose flanks look eaten by the window: widen and retry.
                let wide = widen_region(region, read_len, left, right);
                let wspan = delimit(aligner, read, &wide, reference, context, align_scratch);
                match wspan {
                    RepeatSpan::Between(wtract)
                        if !flank_truncated(&wide, &to_usize(&wtract), read_len, left, right) =>
                    {
                        complete_or_low_quality(read, &wide, &wtract, qual_buffer)
                    }
                    _ => Classified::NoObservation(NoObservationReason::WindowTruncated),
                }
            }
            RepeatSpan::Between(tract) => {
                complete_or_low_quality(read, &region, &tract, qual_buffer)
            }
            // **An empty one-sided span is a read outside the locus, not a lower bound of
            // zero** — and `partial` is where that is decided, in one place.
            //
            // This arm used to hold a `if tract.is_empty()` guard of its own, added at C0
            // when `ReadWitness::from_left`/`from_right` were infallible and nothing else
            // could reject the case. C2 gave them an `Option`, so `partial` has to handle
            // "no position covered" anyway and answers exactly this — which made the guard
            // **unfalsifiable**: the Milestone C review deleted it outright and the whole
            // suite stayed green, because no input can tell the two paths apart. Two
            // spellings of one decision, one of them untestable, is worse than one.
            RepeatSpan::FromLeft(tract) => {
                partial(read, &region, &tract, AnchoredBorder::Left, locus)
            }
            RepeatSpan::FromRight(tract) => {
                partial(read, &region, &tract, AnchoredBorder::Right, locus)
            }
            RepeatSpan::Unanchored => {
                Classified::NoObservation(NoObservationReason::NoBorderAnchored)
            }
        }
    }

    /// The tract's span inside `region`, for a read the aligner has nothing to decide about.
    ///
    /// **The condition is that the read reproduces the reference across the whole window.** Then
    /// the read carries the reference's repeat count — a different count would have to appear in
    /// the CIGAR as an insertion or a deletion — and the tract sits at the reference's own
    /// offsets, so the answer is arithmetic rather than a search.
    ///
    /// Four things are required, and each rules out a way the CIGAR could be lying:
    ///
    /// - **No insertion or deletion anywhere in the read.** One outside the window still shifts
    ///   every reference-to-read offset after it, and `region` is derived from those offsets.
    /// - **No clip at either end.** A soft clip is the one place a long allele's extra units can
    ///   hide from the aligned span — production's own gate admits a clipped read unconditionally
    ///   for exactly this reason (`src/ssr/pileup/footprint.rs:210-222`) — so a clipped read goes
    ///   to the aligner however clean the rest of its CIGAR looks.
    /// - **The aligned span brackets the whole window**, or the read is a partial observation and
    ///   belongs on the aligner's path, which is what decides how far into the tract it reached.
    /// - **The bases match the reference over the whole window**, not merely their count. The
    ///   delimiter scores emissions by base quality, so *equal length* alone would not let this
    ///   function predict its answer; *equal bases* does, because then the best path is the
    ///   diagonal at any quality — every alternative either mismatches a base or pays a slip's
    ///   open cost, and both are strictly worse than a run of exact matches.
    ///
    /// The last is why a `SeqMismatch` in the CIGAR is not enough to disqualify a read on its own
    /// and is not tested for: what matters is the bases, and they are compared directly.
    ///
    /// **A fifth condition is about the locus rather than the read, and it was found by
    /// measurement rather than by argument.** Where the reference's own repeat run continues past
    /// the boundary region typing drew, the delimiter answers a *lower bound* and not a length —
    /// on a tomato run it returned `FromLeft(15..42)` for a 12-base `A` tract whose right flank is
    /// more `A`, because the run reaches the window edge and might go further. That is the right
    /// answer and this function must not overrule it, so a locus whose run is extendable in the
    /// reference is left to the aligner entirely. It is checked against the reference and so is
    /// the same verdict for every read at that locus.
    fn tract_needing_no_alignment(
        read: &AlignedRead,
        region: &Range<usize>,
        locus: &SsrLocus,
        reference: &[u8],
    ) -> Option<Range<u64>> {
        if reference_run_escapes_the_tract(locus, reference) {
            return None;
        }
        let straight = read.cigar.iter().all(|op| {
            matches!(
                op,
                CigarOp::Match(_) | CigarOp::SeqMatch(_) | CigarOp::SeqMismatch(_)
            )
        });
        if !straight {
            return None;
        }
        // With no indel and no clip, `extract_region` returns exactly the window's width when the
        // aligned span brackets it, and something shorter when it does not — so this one length
        // test covers the bracketing requirement.
        if region.len() != reference.len() {
            return None;
        }
        if &read.seq[region.clone()] != reference {
            return None;
        }
        let left = locus.left_flank_len() as u64;
        Some(left..left + locus.segment.tract_len())
    }

    /// Does the reference's own repeat run carry on past the tract region typing drew?
    ///
    /// A run of period `p` continues one base to the left when that base repeats the one `p`
    /// further along, and one base to the right when it repeats the one `p` back. Either way the
    /// delimiter sees a run that does not end inside the window, so it answers a lower bound
    /// rather than a length, and [`tract_needing_no_alignment`] must stand aside.
    ///
    /// A flank shorter than the period cannot be tested and is treated as escaping, which costs
    /// an alignment and never a wrong answer.
    fn reference_run_escapes_the_tract(locus: &SsrLocus, reference: &[u8]) -> bool {
        let period = locus.segment.motif().as_bytes().len();
        let left = locus.left_flank_len();
        let tract_end = left + locus.segment.tract_len() as usize;
        if period == 0 || left < period || reference.len() < tract_end + period {
            return true;
        }
        reference[left - 1] == reference[left - 1 + period]
            || reference[tract_end] == reference[tract_end - period]
    }

    /// Align the read's `region` slice against the reference frame.
    fn delimit<A: RepeatDelimiter>(
        aligner: &A,
        read: &AlignedRead,
        region: &Range<usize>,
        reference: &[u8],
        context: RepeatContext<'_>,
        scratch: &mut A::Scratch,
    ) -> RepeatSpan {
        let bases = ReadBases::try_new(&read.seq[region.clone()], &read.qual[region.clone()])
            .expect("seq and qual are the same slice, hence equal length");
        aligner.align(bases, reference, context, scratch)
    }

    fn to_usize(range: &Range<u64>) -> Range<usize> {
        range.start as usize..range.end as usize
    }

    /// A complete tract, if it clears the base-quality gate; else `LowQuality`. `tract` is
    /// relative to the `region` slice.
    fn complete_or_low_quality(
        read: &AlignedRead,
        region: &Range<usize>,
        tract: &Range<u64>,
        qual_buffer: &mut Vec<u8>,
    ) -> Classified {
        let region_seq = &read.seq[region.clone()];
        let region_qual = &read.qual[region.clone()];
        let tract = to_usize(tract);
        let tract_qual = &region_qual[tract.clone()];
        if passes_quality_gate(tract_qual, MIN_REGION_Q1, qual_buffer) {
            Classified::Observed {
                bases: region_seq[tract].into(),
                read_witness: ReadWitness::Complete,
                q_sum: ln_p_err_sum(tract_qual),
            }
        } else {
            Classified::NoObservation(NoObservationReason::LowQuality)
        }
    }

    /// `Σ ln(P_err)` over base qualities — the per-read base-quality error mass on the tract, in
    /// log-error space (freebayes' `q_sum` convention). Negative and monotone: higher quality →
    /// less error mass. The BQ support moment (spec §3); there is no BAQ on this path and MAPQ is
    /// carried separately, so this is the base-quality term alone.
    fn ln_p_err_sum(quals: &[u8]) -> f64 {
        const LN10_OVER_10: f64 = std::f64::consts::LN_10 / 10.0;
        quals.iter().map(|&q| -(q as f64) * LN10_OVER_10).sum()
    }

    /// Which border of the tract a partial read held — the half of `RepeatSpan`'s answer that
    /// decides where the witnessed run sits inside the locus.
    ///
    /// A plain marker, **not** a `ReadWitness` constructor passed by value: since the reshape
    /// (spec §6) the run is a struct variant, which cannot be used as a function, and building it
    /// needs the locus length as well as the reach.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum AnchoredBorder {
        Left,
        Right,
    }

    /// A partial (lower-bound) observation: the tract bases the read showed, as the run of locus
    /// positions it witnessed, anchored at whichever border held.
    fn partial(
        read: &AlignedRead,
        region: &Range<usize>,
        tract: &Range<u64>,
        border: AnchoredBorder,
        locus: &SsrLocus,
    ) -> Classified {
        // **The span has to be ordered, and only a `debug_assert` upstream says so.**
        // `classify`'s own ordering check (`ssr_best_path_flat_gap.rs`) is debug-only on a
        // struct with public fields — kept public for the delimiter parity harness — and
        // this repo has recorded twice that debug-only guards compile out of the build it
        // actually runs. An inverted span underflows `tract.end - tract.start` below: a
        // panic in debug, and in release a wrap that clamps to `u16::MAX` and then dies on
        // the slice index, naming neither the cause nor the read. The `tract.is_empty()`
        // guard removed from `classify_read` used to swallow this shape too, since
        // `Range::is_empty` is `start >= end` (Milestone C review).
        debug_assert!(
            tract.start <= tract.end,
            "inverted tract span {tract:?} reached `partial`",
        );
        let region_seq = &read.seq[region.clone()];
        let region_qual = &read.qual[region.clone()];
        let tract = to_usize(tract);
        // `reach` is the observed tract length in **read** coordinates, which diverge from locus
        // positions under stutter — an expanded allele reaches further in read bases than the
        // reference tract has positions.
        //
        // # The convention that turns a read length into reference positions (owner, 2026-07-31)
        //
        // A witness is a set of **reference** positions, so `reach` has to be placed on the
        // reference before it can be stored, and inside a repeat that placement needs a rule:
        // if a read shows 12 repeat bases where the reference tract has 10 positions, *which*
        // two are the extra copies is not determined by the sequence.
        //
        // **The rule is: lay the read's repeat down from the border it anchored.** A read that
        // held the left flank starts its repeat at the tract's left border; one that held the
        // right flank ends its repeat at the right border. Bases beyond the tract's far border
        // are then extra copies — an insertion — rather than positions of the tract. It is the
        // same rule indel left-alignment uses, with the anchored side choosing the direction.
        //
        // Two consequences, and the clamp in `ReadWitness::from_left` / `from_right` is exactly
        // this rule implemented:
        //
        // 1. a reach **shorter** than the tract covers that many positions from the anchored
        //    border, which is what the constructors build directly;
        // 2. a reach **at or past** the tract length covers the tract end to end — every
        //    reference position — with the surplus falling outside as inserted copies. That is
        //    what the clamp produces, and it is correct rather than a saturation artefact.
        //
        // **Covering every position is still not a measurement.** The read anchored one border
        // and then ran out; the allele can continue past what it showed, so the evidence stays a
        // lower bound. That is why this path builds a *partial* witness in both cases and only a
        // read holding **both** flanks is `Complete` — and why the dumps spell case 2
        // `partial:both` rather than folding it into `partial:left` (spec §8; 2,530 of 6,216
        // partial observations on chr01 of tomato SRR7279503).
        let reach = (tract.end - tract.start).min(u16::MAX as usize) as u16;
        let locus_len = LocusLen::from_positions(locus.segment.tract_len());
        // **The read is outside the locus, and this is where that is decided.** The
        // constructors answer `None` when the clamped run covers no position, which is
        // exactly the read that anchored a flank and crossed no tract position: the
        // delimiter's rejected side falls back to the read's own edge, so an alignment that
        // never entered the repeat comes back as `tract_start == tract_end`. The bound it
        // would supply is *at least zero* — no evidence — and it arrives in bulk, because
        // the fetch queries the tract plus its margin so every read clipping the flank is
        // delimited. Those bases are the SNP/indel path's to analyse (C0; see
        // `NoObservationReason::OutsideTract` for the measurement).
        //
        // A second `tract.is_empty()` guard stood in `classify_read` from C0 until the
        // Milestone C review showed no input could distinguish the two — C2's `Option` had
        // made this branch answer first. One decision, one place.
        let Some(read_witness) = (match border {
            AnchoredBorder::Left => ReadWitness::from_left(reach, locus_len),
            AnchoredBorder::Right => ReadWitness::from_right(reach, locus_len),
        }) else {
            return Classified::NoObservation(NoObservationReason::OutsideTract);
        };
        Classified::Observed {
            bases: region_seq[tract.clone()].into(),
            read_witness,
            // Partials are kept without the quality gate (spec §3), but still carry the BQ moment.
            q_sum: ln_p_err_sum(&region_qual[tract]),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::ng::alignment::PerQualityEmission;
        use crate::ng::alignment::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
        use crate::ng::region_typing::segment_criteria::{Motif, SsrSegment};
        use crate::ng::types::{Position, ReadGroupId};
        use crate::pileup::walker::CigarOp;

        // Reference frame: 6-base flanks around a CACACA tract → "GGGGGGCACACATTTTTT".
        const FRAME: &[u8] = b"GGGGGGCACACATTTTTT";

        fn locus() -> SsrLocus {
            SsrLocus {
                segment: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"CA").unwrap(), 1.0)
                    .unwrap(),
                tract_with_margin_bases: FRAME.into(),
                margin_start: Position(1),
            }
        }

        fn read(seq: &[u8], qual_value: u8) -> AlignedRead {
            AlignedRead {
                qname: b"r".to_vec(),
                flag: 0,
                ref_id: 0,
                pos: 1,
                mapq: 60,
                cigar: vec![CigarOp::Match(seq.len() as u32)],
                seq: seq.to_vec(),
                qual: vec![qual_value; seq.len()],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                read_group: ReadGroupId(0),
            }
        }

        fn classify(read: &AlignedRead) -> Classified {
            let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
            let stutter = StutterModel::hipstr_shipped();
            let mut scratch = ViterbiScratch::new();
            let mut qual_buffer = Vec::new();
            classify_read(
                read,
                &locus(),
                &aligner,
                &stutter,
                &mut scratch,
                &mut qual_buffer,
            )
        }

        /// A read spanning the whole frame measures the tract exactly — a complete observation
        /// of the reference allele.
        #[test]
        fn a_spanning_read_yields_a_complete_observation() {
            match classify(&read(FRAME, 40)) {
                Classified::Observed {
                    bases,
                    read_witness: ReadWitness::Complete,
                    q_sum,
                } => {
                    assert_eq!(&*bases, b"CACACA");
                    // Six Q40 bases: Σ ln(P_err) = 6 · 40 · (−ln10/10), all negative.
                    assert!(q_sum < 0.0, "the BQ error mass is negative, got {q_sum}");
                }
                other => panic!("expected a complete observation, got {other:?}"),
            }
        }

        /// A read whose tract quality is below the gate is `LowQuality`, not an observation.
        #[test]
        fn a_low_quality_tract_is_rejected() {
            // Whole read at quality 5 (< MIN_REGION_Q1 = 15).
            let classified = classify(&read(FRAME, 5));
            assert_eq!(
                classified,
                Classified::NoObservation(NoObservationReason::LowQuality)
            );
        }

        /// A read that anchors the left flank but ends inside the tract is a partial (left)
        /// observation — a lower bound production would have discarded.
        #[test]
        fn a_read_running_off_the_right_is_a_left_partial() {
            // left flank + 4 tract bases, then the read ends (no right flank).
            let classified = classify(&read(b"GGGGGGCACA", 40));
            match classified {
                Classified::Observed {
                    bases,
                    read_witness: ReadWitness::Partial { positions },
                    ..
                } => {
                    assert_eq!(&*bases, b"CACA");
                    // Flush with the left border: the run starts at locus position 0 and
                    // covers four. Half-open, so it ends at 4.
                    assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(0, 4)]);
                }
                other => panic!("expected a left partial, got {other:?}"),
            }
        }

        /// A read that anchors the right flank but begins inside the tract is a partial
        /// (right) — the mirror of the left case, which a left/right swap would fail.
        ///
        /// **The offset is the assertion that matters here.** Since the reshape the side is not
        /// a variant but a derivation from where the run sits, so `offset_in_locus == 2` on a
        /// 6-base tract is what says "flush with the right border" — and it is the one number a
        /// mint that forgot the locus length could not produce.
        #[test]
        fn a_read_running_off_the_left_is_a_right_partial() {
            // 4 tract bases + full right flank, mapped at the tract's 3rd base (0-based 8).
            let mut read = read(b"CACATTTTTT", 40);
            read.pos = 9;
            match classify(&read) {
                Classified::Observed {
                    bases,
                    read_witness: ReadWitness::Partial { positions },
                    ..
                } => {
                    assert_eq!(&*bases, b"CACA");
                    // The tract is "CACACA" — 6 positions — so a 4-position run flush with the
                    // right border starts at 2 and ends at the locus's own end.
                    assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(6 - 4, 6)]);
                }
                other => panic!("expected a right partial, got {other:?}"),
            }
        }

        /// **A read that stops at the tract's edge is outside the locus, not a partial of
        /// length zero.**
        ///
        /// It covers the left flank and none of `CACACA`. The delimiter anchors the left
        /// flank and, having no right anchor, falls the far end back to the read's own edge
        /// — which is the same offset, so the span is empty. Before this was classified it
        /// became an observation with **empty bases** and a witness covering zero positions,
        /// and it was not a corner: 6,704 reads against 7,085 genuine partials on tomato
        /// chr01 of `SRR7279503`, at a median window overlap of 16 bases against a 30-base
        /// flank. Those bases belong to the SNP/indel path (owner, 2026-07-30).
        ///
        /// Both borders are asserted because the two arrive by mirror-image routes — the
        /// left case through `tract_start == region_len`, the right through
        /// `tract_end == 0`. **They now share one decision point** (`partial`'s handling of
        /// the constructors' `None`), so neither can be fixed without the other; asserting
        /// both is what keeps that true if the routes ever separate again. The mutation that
        /// discriminates is on the constructor: make `from_left` saturate a zero-length run
        /// to one position and both cases here fail, along with
        /// `a_constructor_asked_for_no_positions_answers_none`.
        #[test]
        fn a_read_covering_only_a_flank_is_outside_the_tract() {
            // The left flank alone: the read ends where the tract begins.
            assert_eq!(
                classify(&read(b"GGGGGG", 40)),
                Classified::NoObservation(NoObservationReason::OutsideTract),
                "a read that stops at the tract's first base witnessed no tract position",
            );

            // The right flank alone, mapped at the base just past the tract.
            let mut right_only = read(b"TTTTTT", 40);
            right_only.pos = 13;
            assert_eq!(
                classify(&right_only),
                Classified::NoObservation(NoObservationReason::OutsideTract),
                "and the mirror case, which arrives through the other end of the span",
            );
        }

        /// A read wholly inside the tract anchors neither flank — no per-read fact, so no
        /// observation.
        #[test]
        fn a_read_wholly_inside_the_tract_is_unanchored() {
            let mut read = read(b"CACACA", 40); // exactly the tract, no flanks
            read.pos = 7;
            assert_eq!(
                classify(&read),
                Classified::NoObservation(NoObservationReason::NoBorderAnchored)
            );
        }

        /// A long-allele read the mapper laid down all-Match pushes its far flank out of the
        /// ref-sized window; the first delimit sees a truncated flank, and widening the read
        /// slice recovers the full long allele (spec §2, the long-allele window recovery).
        #[test]
        fn a_window_truncated_long_allele_is_recovered_by_widening() {
            // 5 CA units (10 bp tract) instead of the reference's 3, mapped 22M.
            match classify(&read(b"GGGGGGCACACACACATTTTTT", 40)) {
                Classified::Observed {
                    bases,
                    read_witness: ReadWitness::Complete,
                    ..
                } => assert_eq!(&*bases, b"CACACACACA"),
                other => panic!("expected the recovered long allele, got {other:?}"),
            }
        }
    }
}

// ---------------------------------------------------------------------
// D2c — the tally: fold each kept read's `Classified` into the locus's observed sequences,
// deduping by `(bases, read_witness)`, accumulating the support moments, and counting the
// no-observation reasons. A port of production's `tally` (`src/ssr/pileup/locus_tally.rs`),
// extended with partial observations and the strand/BQ/MAPQ moments the shared type carries
// (spec §3, §6). Consumed by the generator (D3).
// ---------------------------------------------------------------------
mod tally {
    use super::SsrGeneratorCounts;
    use super::classify::{Classified, NoObservationReason};
    use crate::bam::alignment_input::FLAG_REVERSE_STRAND;
    use crate::ng::locus_generation::{ReadWitness, SequenceObservation};
    use crate::ng::read::aligned_read::AlignedRead;
    use crate::ng::types::ReadGroupId;
    use std::collections::HashMap;

    /// The per-locus tally the generator folds onto the `SampleLocusObservations`: the deduped
    /// observed sequences and how many reads reached the aligner but yielded nothing.
    /// (`reads_fetched` / `reads_discarded_by_cap` are the caller's to set — they come from the
    /// cap, not the outcomes.)
    pub(super) struct SsrTally {
        pub(super) observations: Vec<SequenceObservation>,
        pub(super) reads_without_observation: u32,
    }

    /// Accumulated support for one distinct `(bases, read_witness, read_group)` bucket — the moments summed
    /// as reads fold in, materialised into a [`SequenceObservation`] at the end.
    #[derive(Default)]
    struct Support {
        num_obs: u32,
        num_fwd: u32,
        q_sum: f64,
        mapq_sum: u32,
        mapq_sum_sq: u64,
        placed_left: u32,
    }

    /// Fold each kept read's classification into the locus tally.
    ///
    /// Observations dedup by **`(bases, read_witness, read_group)`** — a `Complete` and a
    /// partial of the same bases are different evidence and stay separate observations, two identical
    /// partials from one read group merge, and the **same allele seen from two read groups is two
    /// observations** (spec §3, §6). Each bucket accumulates the strand (`num_fwd` off the reverse-strand
    /// flag), BQ (`q_sum`, off the read's tract error mass), MAPQ and `placed_left` moments;
    /// `chain_ids` stays empty because the STR path does not phase. `observations` is sorted
    /// by `(bases, witness, read_group)`, so — like production's `tally` — the bases, the counts
    /// and the integer moments are independent of the order reads were folded. `q_sum` is the one
    /// exception: it sums in `f64`, which is commutative but not associative, so a bucket's
    /// `q_sum` can differ in its low bits under a different fold order. That is immaterial while
    /// `q_sum` is soft and unconsumed, and the parity oracle checks only bytes and counts (spec
    /// §6).
    ///
    /// `locus_start` is the tract's 1-based anchor, needed only for `placed_left` — the count of
    /// supporting reads that began strictly left of it, which is production's own definition
    /// (`open_record.rs`'s `alignment_start < rec_pos`).
    ///
    /// `counts` carries the **run-level** totals (complete/partial observations and the four
    /// no-observation reasons); the returned `reads_without_observation` is this locus's own
    /// total.
    pub(super) fn tally<'a>(
        reads_and_outcomes: impl IntoIterator<Item = (&'a AlignedRead, Classified)>,
        locus_start: u64,
        counts: &mut SsrGeneratorCounts,
    ) -> SsrTally {
        let mut buckets: HashMap<(Box<[u8]>, ReadWitness, ReadGroupId), Support> = HashMap::new();
        let mut reads_without_observation = 0u32;
        for (read, outcome) in reads_and_outcomes {
            match outcome {
                Classified::Observed {
                    bases,
                    read_witness,
                    q_sum,
                } => {
                    match read_witness {
                        ReadWitness::Complete => counts.observations_complete += 1,
                        ReadWitness::Partial { .. } => counts.observations_partial += 1,
                    }
                    let support = buckets
                        .entry((bases, read_witness, read.read_group))
                        .or_default();
                    support.num_obs += 1;
                    if read.flag & FLAG_REVERSE_STRAND == 0 {
                        support.num_fwd += 1;
                    }
                    support.q_sum += q_sum;
                    let mapq = u32::from(read.mapq);
                    support.mapq_sum += mapq;
                    support.mapq_sum_sq += u64::from(mapq) * u64::from(mapq);
                    // Production's rule verbatim: strictly left of the anchor, not at it.
                    // `placed_start` — the "exactly at it" half — is deliberately not carried.
                    support.placed_left += u32::from(read.pos < locus_start);
                }
                Classified::NoObservation(reason) => {
                    reads_without_observation += 1;
                    match reason {
                        NoObservationReason::NoBorderAnchored => counts.no_border_anchored += 1,
                        NoObservationReason::LowQuality => counts.low_quality += 1,
                        NoObservationReason::WindowTruncated => counts.window_truncated += 1,
                        NoObservationReason::OutsideTract => counts.outside_tract += 1,
                    }
                }
            }
        }

        let mut observations: Vec<SequenceObservation> = buckets
            .into_iter()
            .map(
                |((bases, read_witness, read_group), support)| SequenceObservation {
                    bases,
                    read_witness,
                    read_group,
                    num_obs: support.num_obs,
                    num_fwd: support.num_fwd,
                    q_sum: crate::ng::types::SummedLogError::from_nats(support.q_sum),
                    mapq_sum: support.mapq_sum,
                    mapq_sum_sq: support.mapq_sum_sq,
                    placed_left: support.placed_left,
                    // The STR path does not phase — there is no chain id to fold (spec §3).
                    chain_ids: Vec::new(),
                },
            )
            .collect();
        // The read group joins the sort key because it joined the bucket key: without it two
        // observations differing only by group would tie, and `HashMap` iteration order is seeded per
        // process — the output would be non-deterministic run to run on any multi-group sample.
        observations.sort_unstable_by(|a, b| {
            a.bases
                .cmp(&b.bases)
                .then_with(|| a.read_witness.sort_key().cmp(&b.read_witness.sort_key()))
                .then_with(|| a.read_group.cmp(&b.read_group))
        });

        SsrTally {
            observations,
            reads_without_observation,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::ng::locus_generation::LocusLen;
        use crate::ng::types::ReadGroupId;
        use crate::pileup::walker::CigarOp;

        /// An `AlignedRead` with a given strand flag and MAPQ — the only fields the tally reads
        /// off the read; the sequence and qualities live in the `Classified` handed alongside.
        fn read(flag: u16, mapq: u8) -> AlignedRead {
            AlignedRead {
                qname: b"r".to_vec(),
                flag,
                ref_id: 0,
                pos: 1,
                mapq,
                cigar: vec![CigarOp::Match(6)],
                seq: b"CACACA".to_vec(),
                qual: vec![40; 6],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                read_group: ReadGroupId(0),
            }
        }

        /// The same, in a named read group and starting at `pos` — the two fields the
        /// `(…, read_group)` bucket key and `placed_left` read.
        fn read_in_group(group: u32, pos: u64) -> AlignedRead {
            AlignedRead {
                pos,
                read_group: ReadGroupId(group),
                ..read(0, 60)
            }
        }

        fn observed(bases: &[u8], read_witness: ReadWitness, q_sum: f64) -> Classified {
            Classified::Observed {
                bases: bases.into(),
                read_witness,
                q_sum,
            }
        }

        /// **The read group is part of the identity**: one allele seen from two read groups is
        /// two observations, and their `num_obs` sum to what a single-group tally would have reported
        /// (spec §6). This is the check that the split is *computed* and not defaulted — with
        /// `read_group` left at a constant both reads would merge into one observation of two.
        #[test]
        fn one_allele_from_two_read_groups_is_two_observations_that_sum_back() {
            let rg0 = read_in_group(0, 10);
            let rg1 = read_in_group(1, 10);
            let outcomes = vec![
                (&rg0, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (&rg0, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (&rg1, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let observations = tally(outcomes, 10, &mut counts).observations;

            assert_eq!(
                observations.len(),
                2,
                "one observation per (allele, read group)"
            );
            assert!(
                observations
                    .iter()
                    .all(|observation| &*observation.bases == b"CACACA"),
                "both observations are the same allele"
            );
            assert_eq!(
                observations
                    .iter()
                    .map(|observation| observation.read_group)
                    .collect::<Vec<_>>(),
                vec![ReadGroupId(0), ReadGroupId(1)],
                "and they are ordered by group, which is what keeps the output deterministic"
            );
            assert_eq!(
                observations
                    .iter()
                    .map(|observation| observation.num_obs)
                    .sum::<u32>(),
                3,
                "collapsing the group axis recovers the single-group total exactly"
            );
        }

        /// **The group tie-break is asserted on the sort key, not on a lucky fold order.**
        ///
        /// The obvious form of this test — fold two groups and check the observation order — rides on
        /// `HashMap` iteration, which is seeded per process: with the `then_with` deleted it
        /// passes roughly six runs in ten, so the regression it exists to catch would reach a
        /// green CI most of the time. Enough observations that a wrong order cannot be a coin flip is
        /// the fix: six groups over two alleles, whose sorted order is fully determined.
        #[test]
        fn observations_sort_by_group_within_an_allele_deterministically() {
            let reads: Vec<AlignedRead> = (0..6).map(|g| read_in_group(g, 10)).collect();
            let mut outcomes = Vec::new();
            for (i, r) in reads.iter().enumerate() {
                // Alternate the alleles so group order and fold order disagree.
                let bases: &[u8] = if i % 2 == 0 { b"AA" } else { b"CC" };
                outcomes.push((r, observed(bases, ReadWitness::Complete, -1.0)));
            }
            let mut counts = SsrGeneratorCounts::default();
            let observations = tally(outcomes, 10, &mut counts).observations;

            let order: Vec<(&[u8], u32)> = observations
                .iter()
                .map(|observation| (observation.bases.as_ref(), observation.read_group.get()))
                .collect();
            assert_eq!(
                order,
                vec![
                    (b"AA".as_ref(), 0),
                    (b"AA".as_ref(), 2),
                    (b"AA".as_ref(), 4),
                    (b"CC".as_ref(), 1),
                    (b"CC".as_ref(), 3),
                    (b"CC".as_ref(), 5),
                ],
                "bases first, then read group — ascending, and never fold order"
            );
        }

        /// **An expanded allele merges the two sides into one observation, and that is a real
        /// behaviour change the reshape brought.**
        ///
        /// `reach` is measured in *read* bases; on an allele longer than the reference tract it
        /// exceeds the locus length, so `from_left` and `from_right` both clamp to the whole
        /// locus and produce the **same run**. Two reads anchored at opposite borders with the
        /// same bases then share a bucket key and merge — where `PartialLeft(n)` and
        /// `PartialRight(n)` kept them as two observations.
        ///
        /// It is arguably the right answer (identical constraints are one observation) but it is not
        /// the pre-reshape answer, it is invisible on any fixture whose reads are exact
        /// reference slices, and the plan's stated equivalence
        /// `PartialRight(n) ⇔ Partial { len - n, n }` silently stops holding at `n = len`.
        /// Pinned here so it is a decision on the record rather than a surprise in a dump.
        #[test]
        fn an_expanded_allele_merges_the_two_sides_into_one_observation() {
            let locus_len = LocusLen::from_positions(6);
            let left =
                ReadWitness::from_left(9, locus_len).expect("a run covering at least one position");
            let right = ReadWitness::from_right(9, locus_len)
                .expect("a run covering at least one position");
            assert_eq!(
                left, right,
                "the two sides denote the same run once saturated"
            );

            let r = read(0, 60);
            let outcomes = vec![
                (&r, observed(b"CACACACACA", left, -1.0)),
                (&r, observed(b"CACACACACA", right, -1.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let observations = tally(outcomes, 1, &mut counts).observations;

            assert_eq!(
                observations.len(),
                1,
                "one observation, not two — the sides are indistinguishable once the run saturates"
            );
            assert_eq!(observations[0].num_obs, 2, "and both reads support it");
            assert_eq!(
                counts.observations_partial, 2,
                "the run-level per-read tally is unaffected: two reads, two partials"
            );
        }

        /// `placed_left` counts supporting reads that began **strictly left** of the locus
        /// anchor — production's own rule (`alignment_start < rec_pos`), not `<=`.
        ///
        /// The read starting exactly on the anchor is the discriminating case: it is what
        /// separates `placed_left` from `placed_start`, the sibling field ng deliberately does
        /// not carry, and an implementation using `<=` passes every other assertion here.
        #[test]
        fn placed_left_counts_only_reads_starting_strictly_left_of_the_anchor() {
            let before = read_in_group(0, 9);
            let on = read_in_group(0, 10);
            let after = read_in_group(0, 11);
            let outcomes = vec![
                (&before, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (&on, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (&after, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let observations = tally(outcomes, 10, &mut counts).observations;

            assert_eq!(
                observations.len(),
                1,
                "one allele, one group — one observation"
            );
            assert_eq!(observations[0].num_obs, 3);
            assert_eq!(
                observations[0].placed_left, 1,
                "only the read starting at 9 is left of the anchor at 10"
            );
        }

        /// A `Complete` and a partial run of the **same** bases are different evidence, so they
        /// stay as two separate observations (spec §3) — the property the `(bases, read_witness, read_group)` dedup
        /// key rests on.
        #[test]
        fn a_complete_and_a_partial_of_the_same_bases_stay_separate() {
            let fwd = read(0, 60);
            let outcomes = vec![
                (&fwd, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (
                    &fwd,
                    observed(
                        b"CACACA",
                        ReadWitness::from_left(6, LocusLen::from_positions(6))
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            assert_eq!(result.observations.len(), 2);
            assert_eq!(counts.observations_complete, 1);
            assert_eq!(counts.observations_partial, 1);
        }

        /// Two identical partials (same bases, same witness) are the identical constraint, so
        /// they merge into one observation with `num_obs == 2` (spec §3).
        #[test]
        fn two_identical_partials_merge_into_one_count() {
            let fwd = read(0, 60);
            let outcomes = vec![
                (
                    &fwd,
                    observed(
                        b"CACA",
                        ReadWitness::from_left(4, LocusLen::from_positions(6))
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
                (
                    &fwd,
                    observed(
                        b"CACA",
                        ReadWitness::from_left(4, LocusLen::from_positions(6))
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            assert_eq!(result.observations.len(), 1);
            assert_eq!(result.observations[0].num_obs, 2);
        }

        /// The strand, BQ and MAPQ moments accumulate across the reads folded into a bucket:
        /// `num_fwd` counts the forward-strand reads, `q_sum` sums the per-read BQ mass,
        /// `mapq_sum` / `mapq_sum_sq` sum MAPQ and MAPQ². `chain_ids` stays empty.
        #[test]
        fn the_support_moments_accumulate() {
            let fwd = read(0, 60);
            let rev = read(FLAG_REVERSE_STRAND, 30);
            let outcomes = vec![
                (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
                (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            let obs = &result.observations[0];
            assert_eq!(obs.num_obs, 2);
            assert_eq!(obs.num_fwd, 1, "one forward, one reverse");
            assert_eq!(obs.q_sum.nats(), -5.0, "the BQ masses sum");
            assert_eq!(obs.mapq_sum, 90, "60 + 30");
            assert_eq!(obs.mapq_sum_sq, 60 * 60 + 30 * 30);
            assert!(obs.chain_ids.is_empty(), "the STR path does not phase");
        }

        /// Each no-observation reason lands in its own run-level counter, and every such read is
        /// counted in the locus's `reads_without_observation` total.
        #[test]
        fn no_observation_reasons_are_counted_by_reason_and_in_total() {
            let r = read(0, 60);
            let outcomes = vec![
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::NoBorderAnchored),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::WindowTruncated),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
                // C0's reason, and the largest of the four on real data — its tally arm went
                // untested until the Milestone C review said so.
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::OutsideTract),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            assert!(result.observations.is_empty());
            assert_eq!(result.reads_without_observation, 5);
            assert_eq!(counts.no_border_anchored, 1);
            assert_eq!(counts.low_quality, 2);
            assert_eq!(counts.window_truncated, 1);
            assert_eq!(counts.outside_tract, 1);
        }

        /// `observations` is sorted by bytes, then by witness — so the record is identical
        /// regardless of the order reads folded in (production's order-independence, extended to
        /// the witness tie-break).
        #[test]
        fn observations_are_sorted_by_bases_then_witness() {
            let r = read(0, 60);
            let outcomes = vec![
                (&r, observed(b"GG", ReadWitness::Complete, -1.0)),
                (
                    &r,
                    observed(
                        b"AA",
                        ReadWitness::from_right(2, LocusLen::from_positions(6))
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
                (&r, observed(b"AA", ReadWitness::Complete, -1.0)),
                (
                    &r,
                    observed(
                        b"AA",
                        ReadWitness::from_left(2, LocusLen::from_positions(6))
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            let order: Vec<(&[u8], ReadWitness)> = result
                .observations
                .iter()
                .map(|o| (o.bases.as_ref(), o.read_witness.clone()))
                .collect();
            assert_eq!(
                order,
                vec![
                    (b"AA".as_ref(), ReadWitness::Complete),
                    (
                        b"AA".as_ref(),
                        ReadWitness::from_left(2, LocusLen::from_positions(6))
                            .expect("a run covering at least one position")
                    ),
                    (
                        b"AA".as_ref(),
                        ReadWitness::from_right(2, LocusLen::from_positions(6))
                            .expect("a run covering at least one position")
                    ),
                    (b"GG".as_ref(), ReadWitness::Complete),
                ]
            );
        }

        /// **Which of the tie-break's two components dominates** — the claim
        /// [`ReadWitness::sort_key`]'s doc makes, and which the test above cannot check.
        ///
        /// `observations_are_sorted_by_bases_then_witness` builds both of its partials with the
        /// **same** `positions_covered` (`from_left(2, len 6)` = `{0,2}`,
        /// `from_right(2, len 6)` = `{4,2}`), so exchanging the two components of the sort key
        /// maps them to `(1,2,0)` and `(1,2,4)` and the asserted order survives by accident.
        /// The review made exactly that mutation and the whole suite stayed green.
        ///
        /// Here the long run is the left-flush one and the short run the right-flush one, so
        /// "offset outranks length" and "shortest first" disagree — which is the only shape
        /// that can fail. It matters because on the STR path the two are different genetic
        /// constraints, a prefix and a suffix, and their emission order is what a cohort merge
        /// reads.
        #[test]
        fn tally_orders_two_partials_of_one_sequence_by_offset_before_length() {
            let r = read(0, 60);
            let len = LocusLen::from_positions(6);
            let outcomes = vec![
                (
                    &r,
                    observed(
                        b"AA",
                        ReadWitness::from_right(2, len)
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
                (
                    &r,
                    observed(
                        b"AA",
                        ReadWitness::from_left(4, len)
                            .expect("a run covering at least one position"),
                        -1.0,
                    ),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, &mut counts);
            let order: Vec<ReadWitness> = result
                .observations
                .iter()
                .map(|o| o.read_witness.clone())
                .collect();
            assert_eq!(
                order,
                vec![
                    ReadWitness::from_left(4, len).expect("a run covering at least one position"),
                    ReadWitness::from_right(2, len).expect("a run covering at least one position"),
                ],
                "a left-flush run must precede a right-flush one whatever their lengths",
            );
        }

        /// The integer moments of a **multi-read bucket** are order-independent: the same reads
        /// (differing strand and MAPQ) folded into one `(bases, witness)` bucket in two orders
        /// give the same `num_obs` / `num_fwd` / `mapq_sum` / `mapq_sum_sq`. This is the case the
        /// singleton buckets of `tally_is_order_independent` never exercise, and it isolates
        /// `q_sum` as the sole order-sensitive field (its f64 sum is not asserted here).
        #[test]
        fn a_multi_read_bucket_accumulates_integer_moments_order_independently() {
            let fwd = read(0, 60);
            let rev = read(FLAG_REVERSE_STRAND, 30);
            let moments = |outcomes: Vec<(&AlignedRead, Classified)>| {
                let mut counts = SsrGeneratorCounts::default();
                let obs = tally(outcomes, 1, &mut counts).observations;
                assert_eq!(obs.len(), 1, "all reads share one bucket");
                let o = &obs[0];
                (o.num_obs, o.num_fwd, o.mapq_sum, o.mapq_sum_sq)
            };
            let forward_first = moments(vec![
                (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
                (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
            ]);
            let reverse_first = moments(vec![
                (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
                (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
            ]);
            assert_eq!(forward_first, reverse_first);
            assert_eq!(forward_first, (2, 1, 90, 60 * 60 + 30 * 30));
        }

        /// The tally is order-independent: the same multiset of outcomes folded in two different
        /// orders yields the same observed sequences and the same counts. (Its buckets are
        /// singletons, so no `q_sum` f64 reordering occurs — the multi-read integer case is
        /// `a_multi_read_bucket_accumulates_integer_moments_order_independently`.)
        #[test]
        fn tally_is_order_independent() {
            let r = read(0, 60);
            let run = |outcomes: Vec<(&AlignedRead, Classified)>| {
                let mut counts = SsrGeneratorCounts::default();
                let result = tally(outcomes, 1, &mut counts);
                (result.observations, counts)
            };
            let a = run(vec![
                (&r, observed(b"CACACA", ReadWitness::Complete, -1.0)),
                (&r, observed(b"CACA", ReadWitness::Complete, -1.0)),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
            ]);
            let b = run(vec![
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
                (&r, observed(b"CACA", ReadWitness::Complete, -1.0)),
                (&r, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            ]);
            assert_eq!(a, b);
        }
    }
}

// ---------------------------------------------------------------------
// D3 — the generator: one `SsrSegment` → one `SampleLocusObservations`. The four steps
// (build the locus, fetch+cap the reads, classify each, tally) run inside the first
// `next_locus`; the second returns `None`. One locus per segment, including zero coverage
// (spec §2, arch §1, §2).
// ---------------------------------------------------------------------

/// A tract delimiter the generator can drive — any alignment-module aligner that measures a read's
/// repeat: [algorithm 3](crate::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner) (the
/// flat-gap port of production's `delimit_read`, kept as the byte-parity oracle and baseline),
/// [algorithm 4](crate::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner) (the unit-slip
/// model), or [algorithm 4u](crate::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner) (unit-slip
/// hardened by the delimiter bake-off — **the recommended default**). Any [`BestPathAligner`] whose
/// `Output` is a [`RepeatSpan`] over a [`RepeatContext`] qualifies.
///
/// It is a **trait alias** (a supertrait bundle with a blanket impl), the one bound the generator
/// and [`classify::classify_read`] repeat — so the aligner is a **type parameter**, chosen per
/// generator and monomorphised. The choice is thus **static dispatch**: `align` is a direct call in
/// the per-read loop, never a `dyn` indirection (a bake-off runs millions of reads). The only
/// dynamic dispatch is the one-per-run `Box<dyn LocusGenerator>` the dispatcher already uses to hold
/// a generator — outside the hot loop.
pub trait RepeatDelimiter:
    BestPathAligner<Output = RepeatSpan> + for<'a> BestPathAligner<Context<'a> = RepeatContext<'a>>
{
}

impl<A> RepeatDelimiter for A where
    A: BestPathAligner<Output = RepeatSpan>
        + for<'a> BestPathAligner<Context<'a> = RepeatContext<'a>>
{
}

/// The STR locus generator: turns one microsatellite tract into one locus.
///
/// Holds its own accessors and reusable scratch, the "a generator holds its own accessors"
/// convention ([`LocusGenerator`](super::LocusGenerator); spec §2). The delimiter `A` is a **type
/// parameter** ([`RepeatDelimiter`]) so the algorithm is chosen per generator and dispatched
/// statically — [`with_default_aligner`](Self::with_default_aligner) builds the recommended
/// algorithm 4u (the unit-robust aligner — algorithm 4 hardened by the delimiter bake-off), while
/// [`new`](Self::new) takes any aligner, which is how a bake-off swaps in algorithm 3 (the flat-gap
/// port) and how the parity oracle pins it to unit-slip / flat-gap (spec §6).
///
/// **The reference seam (the Arc gap).** Two reference handles, because they serve two needs the
/// current [`RawRefSeq`] design keeps apart:
/// - `reference` — the margin fetch ([`SsrLocus::fetch`], [`RefSeq`] + [`ContigTable`]), a
///   persistent accessor the generator holds.
/// - `make_reference` — a **factory** the cursor needs ([`SampleReads::cursor`]): each file's
///   cursor owns its own reader, and a `RawRefSeq` is a stateful reader that is neither `Clone`
///   nor shareable by `&`, so the cursor takes `FnMut() -> R` rather than a borrow.
///
/// These are the same reference logically; a future `Arc`-shared reference would collapse them to
/// one field with a cheap per-file view. Holding both is the working stopgap the spec flags (spec
/// §8) — **and D3 took most of its cost away**: the factory is called once per file per
/// *chromosome* now, where it used to run at every locus, so a file-backed `R` whose factory
/// reloads is no longer the per-locus tax it was.
pub struct SsrGenerator<R: RawRefSeq + EvictableRefSeq, A: RepeatDelimiter> {
    /// The reference the margin fetch reads the flanks from.
    reference: R,
    /// Builds a fresh read reference for the cursor's mismatch-fraction filter — one per file,
    /// per chromosome (see the type doc).
    ///
    /// **Boxed, which is what keeps it off this type's parameter list**, the same trade the
    /// generic generator took at D2 (arch §3.6). It is called at chromosome boundaries and
    /// nowhere else, so the indirection is not on any path that matters.
    ///
    /// **Not `+ Send`, and a real caller is why.** Review proposed it as free insurance for the
    /// fan-out. It is not free: `Arc<T>` implements [`RawRefSeq`], and `ng_ssr_cohort_stutter`
    /// uses that to hand every file a clone of one accessor — but `WindowedRefSeq` holds a
    /// `RefCell`, so it is `!Sync`, so `Arc<WindowedRefSeq>` is `!Send`, so the closure
    /// returning it is `!Send`. Requiring `Send` here rules out a caller that exists. The
    /// fan-out does not need it either: a worker takes its own generator and its own accessor,
    /// and both generators are already `!Send` through other fields.
    make_reference: Box<dyn FnMut() -> R>,
    /// The sample's cursor and the chromosome it covers — **D3.** One per chromosome, minted by
    /// the first locus on each and rebuilt at every boundary, because nothing in a cursor
    /// survives a chromosome change (`spec/alignment_cursor.md` §4).
    ///
    /// `None` before the first locus. It is *not* cleared per segment: the whole point is that
    /// it outlives them.
    cursor: Option<(ContigId, SampleCursor<R>)>,
    /// What the cursors of *already-retired* chromosomes did — see
    /// [`cursor_counts`](Self::cursor_counts).
    retired_cursor_counts: CursorCounts,
    /// The sample this generator opened its first cursor for. **One sample per generator** —
    /// a second one is refused rather than answered out of the first one's files. See
    /// [`cursor_for`](Self::cursor_for) for what that mistake looks like when it is not caught.
    ///
    /// `None` until the first cursor opens, so a generator is unclaimed until it reads
    /// something.
    sample: Option<SampleIdentity>,
    /// The tract delimiter — the chosen [`RepeatDelimiter`] (algorithm 3 or 4).
    aligner: A,
    /// Reused alignment matrices (the aligner's own scratch type), so it does not reallocate per
    /// read.
    align_scratch: A::Scratch,
    /// The stutter model feeding the aligner's context. Algorithm 3 (flat-gap) ignores it;
    /// algorithm 4 (unit-slip) prices its slips from it. The context carries it either way (spec
    /// §2, arch §5).
    stutter: StutterModel,
    /// Reused scratch for the margin fetch.
    margin_buffer: Vec<u8>,
    /// Reused scratch for the quality gate.
    qual_buffer: Vec<u8>,
    config: SsrGeneratorConfig,
    counts: SsrGeneratorCounts,
    /// The region of the segment begun, from [`begin_segment`](LocusGenerator::begin_segment) —
    /// its `ContigId` is what the margin fetch and the read query key on (`SsrSegment` carries a
    /// contig *name*, not an id).
    current_region: Option<GenomeRegion>,
    /// Latched once the segment's single locus has been produced, so the second `next_locus`
    /// returns `None` (spec §2).
    produced: bool,
}

impl<R> SsrGenerator<R, SsrUnitRobustAligner<PerQualityEmission>>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq,
{
    /// A generator with the **recommended** delimiter — algorithm 4u, the *unit-robust* aligner over
    /// production's [`PerQualityEmission`] table. It is algorithm 4 (unit-slip) hardened by the
    /// delimiter bake-off with a narrow junction guard (a sequencing error near a flank/tract
    /// boundary can no longer slide the boundary) and an evidence-based anchor test (a complete is
    /// only reported when the flank was actually matched, else a lower-bound partial), the demotion
    /// capped so a partial is never lost. On HG002 it removes ~4,500 fabricated exact-lengths per
    /// chr20–22 that unit-slip reported and preserves every real partial (`ng_ssr_gain_loss`).
    ///
    /// The default for real use; a bake-off, or the production byte-parity oracle, passes a specific
    /// aligner to [`new`](Self::new) — parity is checked against unit-slip / the flat-gap port, not
    /// against this, because unit-robust deliberately departs from production's measurement.
    pub fn with_default_aligner(
        reference: R,
        make_reference: impl FnMut() -> R + 'static,
        config: SsrGeneratorConfig,
        bundle_threshold: Bp,
    ) -> Result<Self, SsrGeneratorConfigError> {
        Self::new(
            reference,
            make_reference,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            config,
            bundle_threshold,
        )
    }
}

impl<R, A> SsrGenerator<R, A>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq,
    A: RepeatDelimiter,
{
    /// Build a generator over `reference` (the margin fetch), `make_reference` (the read-query
    /// factory) and `aligner` (the chosen delimiter — algorithm 3 or 4), with `config` checked
    /// against `bundle_threshold`, the region-typing radius the flank must stay within (spec §4).
    /// Fails if `flank_bp > bundle_threshold`. For the recommended default, use
    /// [`with_default_aligner`](Self::with_default_aligner).
    pub fn new(
        reference: R,
        make_reference: impl FnMut() -> R + 'static,
        aligner: A,
        config: SsrGeneratorConfig,
        bundle_threshold: Bp,
    ) -> Result<Self, SsrGeneratorConfigError> {
        config.check_flank_within(bundle_threshold)?;
        Ok(Self {
            reference,
            make_reference: Box::new(make_reference),
            cursor: None,
            retired_cursor_counts: CursorCounts::default(),
            sample: None,
            aligner,
            align_scratch: A::Scratch::default(),
            stutter: StutterModel::hipstr_shipped(),
            margin_buffer: Vec::new(),
            qual_buffer: Vec::new(),
            config,
            counts: SsrGeneratorCounts::default(),
            current_region: None,
            produced: false,
        })
    }

    /// The running STR counts — accumulated across every segment, readable at any point.
    pub fn counts(&self) -> &SsrGeneratorCounts {
        &self.counts
    }

    /// **What the cursors did**, summed over every chromosome walked so far — see
    /// [`PileupGenerator::cursor_counts`](super::pileup::PileupGenerator::cursor_counts) for the
    /// argument, which is the same one.
    ///
    /// In short: a generator that mints a cursor per locus and one that keeps a cursor per
    /// chromosome emit **identical loci**, so nothing in the output can tell them apart, and the
    /// D2 review proved the feature could be switched off with the whole suite green. What
    /// moves is only what the reader avoided (`spec/alignment_cursor.md` §11.5).
    pub fn cursor_counts(&self) -> CursorCounts {
        let mut total = self.retired_cursor_counts;
        if let Some((_, cursor)) = &self.cursor {
            total += cursor.counts();
        }
        total
    }

    /// This sample's cursor for `contig`, minting one if the last locus was on another
    /// chromosome — **D3, and the whole of it.**
    ///
    /// The cursor is what replaces a read query per locus. It stays open and positioned, so two
    /// STR loci close enough to share reads decode them once
    /// (`spec/alignment_cursor.md` §1). A chromosome change is the one thing it cannot absorb:
    /// nothing in a cursor survives it — the kept reads are useless and, on CRAM, so are the
    /// reference bases — so the boundary mints a new one rather than repositioning the old (§4).
    ///
    /// **The STR generator's own regions are far apart**, so this reuses far less than the
    /// generic generator's tiling walk does — and the ratio the Checkpoint C audit measured
    /// collapses with stride (23× dense, 3.8× at a 13× stride, 0.5× walked backwards). So the
    /// gain here is **measured, not assumed**: `ng_ssr_loci_dump` on chromosome 21 of HG002 30×
    /// goes from **17.3 s to 12.5 s of user CPU, −28 %**, with the dump byte-identical. Wall time
    /// does not move, because that tool's wall is the background whole-genome reference MD5.
    ///
    /// The two generators hold **separate** cursors on purpose: their regions interleave, and
    /// sharing one would tie their lifetimes together (spec §12).
    ///
    /// # One sample per generator, and it is checked here
    ///
    /// A cursor is opened for **one sample's files** and kept. Handing this generator a second
    /// sample would answer that sample out of the first one's files — every individual in a
    /// cohort reporting the first one's reads, with no error and output of exactly the right
    /// shape. So the sample is remembered alongside the chromosome and a foreign one is refused
    /// ([`ForeignSample`](LocusGenerationError::ForeignSample)).
    ///
    /// **Refused rather than re-opened**, and the arithmetic is why. A cohort tool walks
    /// region-major — every sample is asked about one repeat before anything moves to the next
    /// repeat — so re-opening on a change of sample would open one reader per sample *per
    /// locus*, which is worse than the per-locus query this replaced. A generator per sample
    /// opens one reader per sample per chromosome, which is what the design assumes anyway
    /// (`arch/alignment_cursor.md` §2.4 counts open files as `files × generators × workers`).
    fn cursor_for(
        &mut self,
        contig: ContigId,
        reads: &SampleReads,
        region: GenomeRegion,
    ) -> Result<&mut SampleCursor<R>, LocusGenerationError> {
        let sample = reads.identity();
        if self.sample.as_ref().is_some_and(|opened| *opened != sample) {
            return Err(LocusGenerationError::ForeignSample { region });
        }
        // Let both readers drop what this walk has gone past. A reference reader is a
        // sliding window over the FASTA that only ever grows unless it is told to shrink, and
        // nobody was telling it: a forward walk of a chromosome ended up holding one byte for
        // every base it had passed. The margin is one flank, which is the furthest back either
        // reader looks from a repeat's own start; being generous costs a re-read and never an
        // answer, because releasing is a hint.
        let keep_from = region.start.get().saturating_sub(self.config.flank_bp.0);
        self.reference.evict_before(keep_from);
        if let Some((_, cursor)) = &self.cursor {
            cursor.evict_reference_before(keep_from);
        }

        if self.cursor.as_ref().is_none_or(|(on, _)| *on != contig) {
            // Taken before the new one is built, so a run never holds two chromosomes'
            // cursors — and their kept reads — at once. The tallies are taken with it: they
            // die with the cursor, and they are what says whether it did anything.
            if let Some((_, retiring)) = self.cursor.take() {
                self.retired_cursor_counts += retiring.counts();
            }
            let cursor = reads
                .cursor(contig, &mut self.make_reference)
                .map_err(|source| LocusGenerationError::OpenReadQuery { region, source })?;
            self.cursor = Some((contig, cursor));
            // Adopted only once a cursor exists, so a failed open leaves the generator
            // unclaimed and a retry with the same sample still works.
            self.sample = Some(sample);
        }
        Ok(&mut self.cursor.as_mut().expect("the cursor was just ensured").1)
    }

    /// Delimit every kept read for `segment` with this generator's aligner and return the per-read
    /// outcomes — the view *behind* [`next_locus`](LocusGenerator::next_locus)'s tally, exposed for a
    /// **delimiter bake-off** (run two generators, one per aligner, and compare read by read).
    ///
    /// It reuses the exact fetch + `classify_read` pipeline `next_locus` uses, so the measurements
    /// match the real calling path — but it does **not** tally observations or touch the run-level
    /// [`counts`](Self::counts), and it does not consume the one-locus latch, so it is a read-only
    /// diagnostic sitting beside `next_locus`. `begin_segment` must have set the region first, the
    /// same contract as `next_locus`.
    pub fn delimit_segment_reads(
        &mut self,
        segment: &SsrSegment,
        reads: &SampleReads,
    ) -> Result<SegmentDelimitations, LocusGenerationError> {
        let region = self
            .current_region
            .expect("begin_segment is called before delimit_segment_reads");
        let contig = region.contig;

        // (1)+(2) mirror next_locus: build the locus (tract ± flank) and fetch+cap its reads. No
        // counts side effects — this path is a diagnostic, not part of the tally.
        let locus = SsrLocus::fetch(
            &self.reference,
            contig,
            segment.clone(),
            self.config.flank_bp,
            &mut self.margin_buffer,
        )
        .map_err(|source| LocusGenerationError::Reference { region, source })?;
        let margin_start = locus.margin_start.get();
        let margin_len = locus.tract_with_margin_bases.len() as u64;
        let query_span = GenomeRegion {
            contig,
            start: Position(margin_start),
            end: Position(margin_start + margin_len - 1),
        };
        // Read out before the cursor is borrowed: `cursor_for` takes `&mut self`, and a config
        // field read in the same call would hold a shared borrow across it.
        let max_reads_per_locus = self.config.max_reads_per_locus;
        let seed = seed_for_segment(segment);
        let capped = fetch_capped_reads(
            self.cursor_for(contig, reads, region)?,
            query_span,
            seed,
            max_reads_per_locus,
        )?;

        // (3) Classify each kept read, keeping the per-read outcome instead of tallying it.
        let mut delimited = Vec::with_capacity(capped.kept.len());
        for read in &capped.kept {
            let observation = match classify::classify_read(
                read,
                &locus,
                &self.aligner,
                &self.stutter,
                &mut self.align_scratch,
                &mut self.qual_buffer,
            ) {
                classify::Classified::Observed {
                    bases,
                    read_witness,
                    ..
                } => Some((read_witness, bases.into_vec())),
                classify::Classified::NoObservation(_) => None,
            };
            delimited.push(ReadDelimitation {
                qname: read.qname.clone(),
                read_seq: read.seq.clone(),
                observation,
            });
        }

        let left = locus.left_flank_len();
        let tract_len = segment.tract_len() as usize;
        let bases = &locus.tract_with_margin_bases;
        Ok(SegmentDelimitations {
            reference_tract: bases[left..left + tract_len].to_vec(),
            left_flank: bases[..left].to_vec(),
            right_flank: bases[left + tract_len..].to_vec(),
            reads: delimited,
        })
    }
}

/// One kept read's delimitation under a generator's aligner — the per-read view *behind* the tally,
/// exposed for a delimiter bake-off ([`SsrGenerator::delimit_segment_reads`]). Not on the normal
/// calling path.
#[derive(Debug, Clone)]
pub struct ReadDelimitation {
    /// The read name — matches the same read across two generators (two aligners).
    pub qname: Vec<u8>,
    /// The read's sequence — the raw evidence a bake-off eyeballs to judge which aligner is right.
    pub read_seq: Vec<u8>,
    /// The measured tract bases and how the read covered the tract, or `None` for a no-observation
    /// (anchored no border, or quality-gated out).
    pub observation: Option<(ReadWitness, Vec<u8>)>,
}

/// Every kept read's delimitation for one segment, plus the locus reference context a bake-off reads
/// the outcomes against. Produced by [`SsrGenerator::delimit_segment_reads`].
#[derive(Debug, Clone)]
pub struct SegmentDelimitations {
    /// The reference tract (REF), tract bases only, no flanks.
    pub reference_tract: Vec<u8>,
    pub left_flank: Vec<u8>,
    pub right_flank: Vec<u8>,
    /// One entry per kept read, in the reservoir's kept order (identical across aligners for the
    /// same config, so two generators' vectors align by index as well as by `qname`).
    pub reads: Vec<ReadDelimitation>,
}

impl<R, A> LocusGenerator<SsrSegment> for SsrGenerator<R, A>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq,
    A: RepeatDelimiter,
{
    fn begin_segment(&mut self, region: GenomeRegion) {
        self.current_region = Some(region);
        self.produced = false;
    }

    /// The STR counts, reachable through a boxed generator — the same need the generic
    /// generator's nine had. This tool-driven path reads them off the concrete type
    /// ([`SsrGenerator::counts`]); a dispatcher-driven one cannot, and now does not have
    /// to.
    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        Some(GeneratorCounts::Ssr(SsrGenerator::counts(self)))
    }

    fn next_locus(
        &mut self,
        segment: &SsrSegment,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        if self.produced {
            return Ok(None);
        }
        self.produced = true;
        let region = self
            .current_region
            .expect("begin_segment is called before next_locus (the LocusGenerator contract)");
        let contig = region.contig;

        // The margin fetch and the read query key on `region.contig` (an id), the reservoir seed
        // on `segment.chrom()` (a name); the `LocusGenerator` contract leaves the region/segment
        // pairing unenforced, so a mis-pairing would silently fetch the wrong contig. Guard it in
        // debug/test builds against the reference's own name for that id.
        debug_assert_eq!(
            self.reference
                .contigs()
                .entries
                .get(contig.get() as usize)
                .map(|entry| entry.name.as_str()),
            Some(segment.chrom()),
            "the region's contig and the segment must name the same contig"
        );

        // (1) Build the locus — fetch the tract ± flank, clamped at contig ends.
        let locus = SsrLocus::fetch(
            &self.reference,
            contig,
            segment.clone(),
            self.config.flank_bp,
            &mut self.margin_buffer,
        )
        .map_err(|source| LocusGenerationError::Reference { region, source })?;

        // (2) Fetch the reads over the tract-plus-margin query span, admitting on relevance
        // (overlap, which `SampleReads` applies) — not spanning — and depth-cap them.
        let margin_start = locus.margin_start.get();
        // `margin_len >= 1` on a successful `fetch` (the tract is non-empty), so the `- 1` for
        // the inclusive end cannot underflow.
        let margin_len = locus.tract_with_margin_bases.len() as u64;
        let query_span = GenomeRegion {
            contig,
            start: Position(margin_start),
            end: Position(margin_start + margin_len - 1),
        };
        // Read out before the cursor is borrowed: `cursor_for` takes `&mut self`, and a config
        // field read in the same call would hold a shared borrow across it.
        let max_reads_per_locus = self.config.max_reads_per_locus;
        let seed = seed_for_segment(segment);
        let capped = fetch_capped_reads(
            self.cursor_for(contig, reads, region)?,
            query_span,
            seed,
            max_reads_per_locus,
        )?;
        // Per-locus discards cannot approach 2^32 (the cap and any real depth are far below), so
        // the `as u32` on the locus field below is exact; the run-level counter keeps the u64.
        let reads_discarded_by_cap = capped.fetched - capped.kept.len() as u64;
        self.counts.reads_fetched += capped.fetched;
        self.counts.reads_discarded_by_cap += reads_discarded_by_cap;

        // (3) Classify each kept read against the locus (delimit, recover, gate).
        let mut outcomes = Vec::with_capacity(capped.kept.len());
        for read in &capped.kept {
            outcomes.push(classify::classify_read(
                read,
                &locus,
                &self.aligner,
                &self.stutter,
                &mut self.align_scratch,
                &mut self.qual_buffer,
            ));
        }

        // (4) Tally the outcomes into the locus's observed sequences + the run-level counts.
        let tallied = tally::tally(
            capped.kept.iter().zip(outcomes),
            segment.start(),
            &mut self.counts,
        );

        // Assemble the output: the tract coordinates, the tract bases only (flanks go in the
        // kind), and the flanks split out of the fetched margin (spec §3).
        let left = locus.left_flank_len();
        let tract_len = segment.tract_len() as usize;
        let bases = &locus.tract_with_margin_bases;
        Ok(Some(SampleLocusObservations {
            region: GenomeRegion {
                contig,
                start: Position(segment.start()),
                end: Position(segment.end()),
            },
            reference_bases: bases[left..left + tract_len].into(),
            observations: tallied.observations,
            reads_without_observation: tallied.reads_without_observation,
            reads_discarded_by_cap: reads_discarded_by_cap as u32,
            kind: LocusKind::Ssr(SsrDetail {
                motif: segment.motif(),
                left_flank: bases[..left].into(),
                right_flank: bases[left + tract_len..].into(),
            }),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::LocusLen;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::region_typing::segment_criteria::Motif;
    use crate::ng::types::ReadGroupId;
    use std::collections::HashSet;

    /// A 100-base contig `chr1` with a known repeating pattern, so a fetched span can be
    /// checked byte-for-byte.
    fn reference_100() -> (InMemoryRefSeq, Vec<u8>) {
        let bases: Vec<u8> = b"ACGT".iter().cycle().take(100).copied().collect();
        (
            InMemoryRefSeq::from_named_contigs(vec![("chr1".to_string(), bases.clone())]),
            bases,
        )
    }

    fn tract(start: u64, end: u64) -> SsrSegment {
        SsrSegment::new("chr1".into(), start, end, Motif::new(b"AC").unwrap(), 1.0).unwrap()
    }

    /// A tract with room either side: both flanks are exactly `flank_bp`, the margin span is
    /// `tract + 2·flank_bp`, and the fetched bytes are exactly that slice of the reference.
    #[test]
    fn a_mid_contig_tract_fetches_equal_flanks_and_the_exact_bytes() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [40, 49] (10 bases), flank 10 → margin [30, 59].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(40, 49), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(30));
        assert_eq!(locus.left_flank_len(), 10);
        assert_eq!(locus.right_flank_len(), 10);
        assert_eq!(locus.tract_with_margin_bases.len(), 30); // 10 + 10 + 10
        // margin [30, 59] 1-based == bases[29..59] 0-based.
        assert_eq!(&*locus.tract_with_margin_bases, &bases[29..59]);
    }

    /// A tract near the contig **start**: the left margin clamps to base 1, so the left flank
    /// is short while the right is full — the flanks are unequal, measured not assumed.
    #[test]
    fn a_tract_near_the_contig_start_has_a_short_left_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [5, 14], flank 10 → margin_start clamps to 1, margin_end 24.
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(5, 14), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(1));
        assert_eq!(locus.left_flank_len(), 4, "clamped: 5 - 1");
        assert_eq!(locus.right_flank_len(), 10, "full");
        assert_ne!(locus.left_flank_len(), locus.right_flank_len());
        assert_eq!(&*locus.tract_with_margin_bases, &bases[0..24]);
    }

    /// A tract near the contig **end**: the right margin clamps to the contig length, so the
    /// right flank is short — again unequal.
    #[test]
    fn a_tract_near_the_contig_end_has_a_short_right_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [92, 96], flank 10 → margin [82, 100].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(92, 96), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(82));
        assert_eq!(locus.left_flank_len(), 10, "full");
        assert_eq!(locus.right_flank_len(), 4, "clamped: 100 - 96");
        assert_eq!(locus.tract_with_margin_bases.len(), 19); // 10 + 5 + 4
        assert_eq!(&*locus.tract_with_margin_bases, &bases[81..100]);
    }

    /// A tract reaching past the contig end is a broken segment: `fetch` leaves the window
    /// unclamped so `fetch_into` rejects it, rather than returning a locus whose derived
    /// right flank would underflow.
    #[test]
    fn a_tract_past_the_contig_end_is_rejected_not_silently_truncated() {
        let (reference, _bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [98, 110] runs past the end.
        assert!(
            SsrLocus::fetch(&reference, ContigId(0), tract(98, 110), Bp(10), &mut buffer).is_err()
        );
    }

    #[test]
    fn default_flank_equals_the_bundle_threshold_default_and_caps_at_1000() {
        let config = SsrGeneratorConfig::default();
        assert_eq!(
            config.flank_bp,
            Bp(crate::ng::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            "flank equals the bundle threshold by default (spec §4)"
        );
        assert_eq!(
            config.max_reads_per_locus,
            Some(DEFAULT_SSR_MAX_READS_PER_LOCUS)
        );
    }

    #[test]
    fn the_flank_check_accepts_within_and_at_the_threshold() {
        // Equal is allowed (the default case), and strictly narrower.
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(20),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
    }

    #[test]
    fn the_flank_check_rejects_a_flank_wider_than_the_bundle_threshold() {
        let error = SsrGeneratorConfig {
            flank_bp: Bp(50),
            max_reads_per_locus: None,
        }
        .check_flank_within(Bp(30))
        .expect_err("50 > 30 must be refused");
        assert!(matches!(
            error,
            SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: 50,
                bundle_threshold: 30,
            }
        ));
    }

    // --- the reservoir port, oracle'd against production's own tests ---------

    #[test]
    fn keeps_everything_when_offered_at_most_capacity() {
        let mut r = Reservoir::new(5, locus_seed("chr1", 100));
        for x in [10u32, 20, 30] {
            r.offer(x);
        }
        assert_eq!(r.seen(), 3);
        assert_eq!(r.into_held(), vec![10, 20, 30]); // first-K kept, in order
    }

    #[test]
    fn caps_at_capacity_and_counts_all_offers() {
        let mut r = Reservoir::new(10, locus_seed("chr1", 100));
        for x in 1..=100u32 {
            r.offer(x);
        }
        assert_eq!(r.seen(), 100);
        assert_eq!(r.into_held().len(), 10);
    }

    #[test]
    fn reservoir_is_deterministic_for_a_fixed_seed_and_order() {
        let run = || {
            let mut r = Reservoir::new(10, locus_seed("chr7", 4242));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_eq!(run(), run());
    }

    #[test]
    fn reservoir_keeps_a_deterministic_subset_far_past_capacity() {
        let run = || {
            let mut r = Reservoir::new(8, locus_seed("chrX", 7));
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            (r.seen(), r.into_held())
        };
        let (seen, held) = run();
        assert_eq!(seen, 10_000);
        assert_eq!(held.len(), 8);
        assert!(
            held.iter().all(|x| (1..=10_000).contains(x)),
            "kept set is a subset of the stream"
        );
        assert_eq!(run().1, held, "kept set is identical across runs");
    }

    #[test]
    fn different_loci_sample_differently() {
        let sample = |chrom, start| {
            let mut r = Reservoir::new(10, locus_seed(chrom, start));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_ne!(sample("chr1", 100), sample("chr1", 101));
        assert_ne!(sample("chr1", 100), sample("chr2", 100));
    }

    #[test]
    fn every_item_can_be_selected_no_structural_exclusion() {
        let mut covered = HashSet::new();
        for seed in 0..500u64 {
            let mut r = Reservoir::new(10, seed);
            for x in 1..=100u32 {
                r.offer(x);
            }
            covered.extend(r.into_held());
        }
        assert_eq!(covered.len(), 100);
    }

    #[test]
    fn locus_seed_is_deterministic_and_distinguishes_loci() {
        assert_eq!(locus_seed("chr1", 100), locus_seed("chr1", 100));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr1", 101));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr2", 100));
    }

    /// The seed trap: [`seed_for_segment`] seeds from the contig **name** and the **0-based**
    /// start (`start - 1`), not the 1-based start — feeding the 1-based start would produce a
    /// different, deterministic kept set that fails parity looking like an aligner bug (spec §4).
    #[test]
    fn seed_for_segment_uses_the_contig_name_and_the_zero_based_start() {
        let segment = tract(101, 110); // 1-based start 101 → 0-based 100
        assert_eq!(seed_for_segment(&segment), locus_seed("chr1", 100));
        assert_ne!(
            seed_for_segment(&segment),
            locus_seed("chr1", 101),
            "the 1-based start is the trap the conversion avoids"
        );
    }

    /// **The parity oracle for the port**: ng's seed and reservoir must produce output
    /// identical to frozen production (`src/ssr/pileup/fetch_reads.rs`), byte for byte. The
    /// self-consistency tests above would survive a drifted constant; this one would not —
    /// it is what makes "byte-faithful port" a checked claim rather than an asserted one.
    /// (Calling production as a test-only oracle mirrors region typing's `build_loci`
    /// differential; ng does not depend on production at run time.)
    #[test]
    fn ng_seed_and_reservoir_match_frozen_production_byte_for_byte() {
        use crate::ssr::pileup::fetch_reads as production;

        for (chrom, start) in [
            ("chr1", 0u32),
            ("chr1", 100),
            ("chrX", 4242),
            ("scaffold_7", 999_999),
        ] {
            assert_eq!(
                locus_seed(chrom, start),
                production::locus_seed(chrom, start),
                "seed for ({chrom}, {start})"
            );
        }

        // The eviction branch, far past capacity — the kept set is where a drifted PRNG
        // constant would show.
        let kept = |seed| {
            let mut r = Reservoir::new(8, seed);
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            r.into_held()
        };
        let prod_kept = |seed| {
            let mut r = production::Reservoir::new(8, seed);
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            r.into_held()
        };
        let seed = locus_seed("chrX", 7);
        assert_eq!(
            kept(seed),
            prod_kept(seed),
            "the kept set must be identical"
        );
    }

    // --- D1: the per-locus read fetch + cap ---------------------------------

    fn span(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A `RawRefSeq` over the fixture contigs (all `A`s), for the cursor's
    /// mismatch-fraction filter — matches the fixture reference the reads are opened against.
    /// Named (`chr1` / `chr2`) so a generator's contig-name invariant holds against the
    /// `SsrSegment`'s `chrom()`; the fetch itself keys on `ContigId`, so the names are only for
    /// that check.
    fn fixture_ref_bases() -> InMemoryRefSeq {
        use crate::ng::read::input::test_fixtures::matching_contigs;
        InMemoryRefSeq::from_named_contigs(
            matching_contigs()
                .iter()
                .map(|(name, len, _)| ((*name).to_string(), vec![b'A'; *len]))
                .collect(),
        )
    }

    /// A cursor over the fixture sample's first contig — what
    /// [`fetch_capped_reads`] is pointed at since D3, in place of the sample plus a
    /// per-locus query.
    fn a_cursor(reads: &SampleReads) -> SampleCursor<InMemoryRefSeq> {
        reads
            .cursor(ContigId(0), fixture_ref_bases)
            .expect("the fixture sample opens a cursor")
    }

    /// Open a `SampleReads` over one indexed BAM of `records`, against the fixture reference.
    fn sample_reads_with(
        records: &[noodles_sam::alignment::RecordBuf],
    ) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::test_fixtures::{
            fixture_reference, header, indexed_bam, matching_contigs,
        };
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// The cap keeps `max_reads_per_locus` of the fetched reads and counts them all; `None`
    /// keeps everything (spec §4).
    #[test]
    fn fetch_caps_reads_and_counts_the_fetched() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..5)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 10, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        let mut cursor = a_cursor(&reads);
        let capped =
            fetch_capped_reads(&mut cursor, span(1, 100), 42, Some(3)).expect("fetch succeeds");
        assert_eq!(capped.fetched, 5);
        assert_eq!(capped.kept.len(), 3, "capped to 3 of the 5 fetched");

        let mut cursor = a_cursor(&reads);
        let uncapped =
            fetch_capped_reads(&mut cursor, span(1, 100), 42, None).expect("fetch succeeds");
        assert_eq!(uncapped.fetched, 5);
        assert_eq!(uncapped.kept.len(), 5, "no cap keeps every read");
    }

    /// The kept subset is a deterministic function of the seed (the thread-invariance the seed
    /// buys, spec §4).
    #[test]
    fn the_capped_kept_set_is_deterministic_for_a_seed() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..20)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 4, 20))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let kept_positions = |seed| {
            let mut positions: Vec<u64> =
                fetch_capped_reads(&mut a_cursor(&reads), span(1, 100), seed, Some(5))
                    .expect("fetch succeeds")
                    .kept
                    .iter()
                    .map(|read| read.pos)
                    .collect();
            positions.sort_unstable();
            positions
        };
        assert_eq!(
            kept_positions(7),
            kept_positions(7),
            "same seed → same kept set"
        );
    }

    // --- D3: the generator ---------------------------------------------------

    /// A generator over the all-`A` fixture reference (`chr1` = 100 bp, `chr2` = 200 bp), with
    /// `flank_bp = 10` (within a bundle threshold of 10) and no cap. The reference and the
    /// read-query factory are the same all-`A` bases, so fixture reads (also all-`A`) clear the
    /// mismatch filter — the tract content is not real STR sequence, which is fine for the
    /// structural checks D3 owns (real tract parity is E's).
    fn ssr_generator() -> SsrGenerator<InMemoryRefSeq, SsrUnitSlipAligner<PerQualityEmission>> {
        // Pinned to unit-slip (algorithm 4), not `with_default_aligner` (now unit-robust): these
        // structural checks assert byte-exact observation baselines established against the
        // production-derived delimiter, so they must not move when the recommended default changes.
        SsrGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrUnitSlipAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold")
    }

    /// **The cursor is kept across loci rather than minted for each — D3, and the test that can
    /// tell.**
    ///
    /// A generator that opens a query per locus and one that keeps a cursor per chromosome emit
    /// **identical loci**; that is the correctness requirement, and it is why the D2 review could
    /// switch the whole feature off with 1,557 tests green. The only observable is what the
    /// reader avoided (`spec/alignment_cursor.md` §11.5), which is what `cursor_counts` exists
    /// for.
    ///
    /// Asserted on `regions_reusing`, not on `reads_decoded`: a fresh cursor has no last region
    /// to compare against, so its first move always jumps — three loci through one cursor jump
    /// once and reuse twice, and three cursors jump three times.
    #[test]
    fn the_cursor_is_kept_across_loci_rather_than_minted_for_each() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("r0", 0, 20, 30),
            read_named_with_length("r1", 0, 40, 30),
            read_named_with_length("r2", 0, 55, 30),
        ]);
        let mut generator = ssr_generator();

        for (start, end) in [(30u64, 39u64), (45, 54), (60, 69)] {
            generator.begin_segment(span(start, end));
            generator
                .next_locus(&tract(start, end), &reads)
                .expect("no fetch error")
                .expect("one locus per segment");
        }

        let counts = generator.cursor_counts();
        assert_eq!(
            (counts.regions_jumping, counts.regions_reusing),
            (1, 2),
            "three ascending loci through one cursor jump once — at the first, which has \
             nothing to reuse — and reuse twice: {counts:?}",
        );
        assert_eq!(
            (counts.reads_decoded, counts.reads_replayed),
            (3, 4),
            "the three reads are decoded once each and handed back four more times — locus 2 \
             replays r0 and r1, locus 3 replays r1 and r2. That ratio *is* the feature: \
             {counts:?}",
        );
        assert_eq!(
            counts.reads_evicted, 1,
            "and r0 ends before the third locus begins, so it is dropped — the only witness \
             that 'the kept set is bounded' is an observation rather than an argument: \
             {counts:?}",
        );
    }

    /// **Every retired chromosome's tallies survive, not just the last one's.**
    ///
    /// With a single retirement, accumulating and *overwriting* are numerically identical — so
    /// `retired += counts` mutated to `retired = counts` passed, and so did dropping four of the
    /// five fields. Both matter: `regions_reusing` and `regions_jumping` are the only observable
    /// that says whether the cursor is being kept at all, and a fold that loses them would leave
    /// the feature detectable on the first chromosome and undetectable on every other.
    ///
    /// chr1 → chr2 → chr1 is two retirements, which is what tells the two apart.
    #[test]
    fn returning_to_a_chromosome_keeps_every_retired_cursors_tallies() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("chr1-a", 0, 20, 30),
            read_named_with_length("chr2-a", 1, 20, 30),
            read_named_with_length("chr2-b", 1, 25, 30),
        ]);
        let mut generator = ssr_generator();

        let mut decoded = Vec::new();
        for (contig, name) in [(0u32, "chr1"), (1, "chr2"), (0, "chr1")] {
            let segment =
                SsrSegment::new(name.into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0).unwrap();
            generator.begin_segment(GenomeRegion {
                contig: ContigId(contig),
                start: Position(30),
                end: Position(39),
            });
            generator
                .next_locus(&segment, &reads)
                .expect("no fetch error")
                .expect("one locus per segment");
            decoded.push(generator.cursor_counts().reads_decoded);
        }

        assert_eq!(
            decoded,
            vec![1, 3, 4],
            "chr1 decodes one read, chr2 two more, and returning to chr1 decodes its one \
             again — a fold that overwrote instead of accumulating would report [1, 2, 1]",
        );
    }

    /// **A locus on a new chromosome mints a cursor, and the old one's tallies survive it.**
    ///
    /// A cursor covers one chromosome and refuses a region on any other (§4), so the generator
    /// has to notice the boundary. Nothing before D3 could see this: a query per locus carried
    /// its contig in the query.
    #[test]
    fn a_locus_on_a_new_chromosome_mints_a_cursor_and_keeps_the_old_tallies() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("on-chr1", 0, 20, 30),
            read_named_with_length("on-chr2", 1, 20, 30),
        ]);
        let mut generator = ssr_generator();

        generator.begin_segment(span(30, 39));
        let first = generator
            .next_locus(&tract(30, 39), &reads)
            .expect("no fetch error")
            .expect("one locus per segment");
        let after_first = generator.cursor_counts();

        // The same coordinates on chr2 — `SsrSegment` carries a contig *name*, the region a
        // `ContigId`, and the generator's own debug assertion pins that they agree.
        let on_chr2 = SsrSegment::new("chr2".into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0)
            .expect("a valid segment");
        generator.begin_segment(GenomeRegion {
            contig: ContigId(1),
            start: Position(30),
            end: Position(39),
        });
        let second = generator
            .next_locus(&on_chr2, &reads)
            .expect("chr2 must be reachable — a cursor stuck on chr1 would refuse it")
            .expect("one locus per segment");
        let after_second = generator.cursor_counts();

        // Not asserted on `region.contig`: a locus's region is built from `begin_segment`'s,
        // never from the cursor, so it re-states the test's own input and cannot fail under any
        // mutation here (review). What is load-bearing is that chr2's *read* arrived.
        assert!(
            !first.observations.is_empty(),
            "chr1's read must reach its locus, or the comparison below is between two nothings",
        );
        assert!(
            !second.observations.is_empty(),
            "chr2's read must reach the locus, not merely fail to error",
        );
        assert!(
            after_first.reads_decoded > 0,
            "the first chromosome must decode something, or this test cannot fail: \
             {after_first:?}",
        );
        assert!(
            after_second.reads_decoded > after_first.reads_decoded,
            "the retiring cursor's tallies must be taken before it is dropped, not replaced: \
             {after_first:?} then {after_second:?}",
        );
    }

    /// **A second sample through one generator is refused, not answered.**
    ///
    /// This is the shape of a real defect, caught by review before it shipped. A generator
    /// opens a reader for one sample's files and keeps it for a whole chromosome, but the
    /// `LocusGenerator` trait hands it a `&SampleReads` afresh on every call. Keyed on the
    /// chromosome alone, it answered the second sample out of the **first** sample's files —
    /// so a cohort tool asking 51 plants about one repeat would have got one plant's reads 51
    /// times, under 51 names, with no error and rows of exactly the right shape.
    ///
    /// The fixture is built so the wrong answer is unmistakable: sample A has a read on the
    /// tract, and sample B's nearest read is 45 bases away — outside the queried span
    /// `[20,49]` — so *any* observation from B is A's.
    #[test]
    fn a_second_sample_through_one_generator_is_refused() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 25, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 0, 80, 30)]);
        let mut generator = ssr_generator();
        let seg = tract(30, 39);

        generator.begin_segment(span(30, 39));
        let a = generator
            .next_locus(&seg, &sample_a)
            .expect("no fetch error")
            .expect("one locus per segment");
        assert_eq!(
            a.observations.len(),
            1,
            "sample A's read must reach the tract, or this test cannot fail",
        );

        generator.begin_segment(span(30, 39));
        let refused = generator.next_locus(&seg, &sample_b);

        match refused {
            Err(LocusGenerationError::ForeignSample { region }) => {
                assert_eq!(
                    region,
                    span(30, 39),
                    "the refusal names the region asked for"
                );
            }
            Ok(Some(locus)) => panic!(
                "sample B was answered instead of refused, and with sample A's read — B has \
                 nothing within 45 bases of this tract: {:?}",
                locus.observations,
            ),
            other => panic!("expected a refusal naming the foreign sample, got {other:?}"),
        }
    }

    /// **A generator per sample is the shape that works** — the other half of the test above,
    /// because "refuses everything" would pass that one on its own.
    #[test]
    fn a_generator_per_sample_gives_each_sample_its_own_reads() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 25, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 0, 80, 30)]);
        let seg = tract(30, 39);

        let observations_for = |reads: &SampleReads| {
            let mut generator = ssr_generator();
            generator.begin_segment(span(30, 39));
            generator
                .next_locus(&seg, reads)
                .expect("no fetch error")
                .expect("one locus per segment")
                .observations
                .len()
        };

        assert_eq!(
            observations_for(&sample_a),
            1,
            "sample A's read covers the tract",
        );
        assert_eq!(
            observations_for(&sample_b),
            0,
            "sample B's read is 45 bases away and must not appear — which is the whole point",
        );
    }

    /// A reference that records what it was told to release, delegating everything else.
    ///
    /// The fixture reference is in memory and holds no window, so "how much is resident" reads
    /// zero however badly the release is wired. Asking *what was released, and when* is
    /// falsifiable where asking *how much is held* is not.
    struct ReleaseSpy {
        inner: InMemoryRefSeq,
        released: std::sync::Arc<std::sync::Mutex<Vec<u64>>>,
    }

    impl EvictableRefSeq for ReleaseSpy {
        fn evict_before(&self, pos: u64) {
            self.released.lock().expect("no panic holds this").push(pos);
        }
    }

    impl RawRefSeq for ReleaseSpy {
        fn fetch_raw_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            self.inner.fetch_raw_into(contig, start_1based, length, dst)
        }
    }

    impl RefSeq for ReleaseSpy {
        fn fetch_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            self.inner.fetch_into(contig, start_1based, length, dst)
        }
    }

    impl ContigTable for ReleaseSpy {
        fn contigs(&self) -> &crate::fasta::ContigList {
            self.inner.contigs()
        }
    }

    /// **Every repeat tells the reference readers what they may release, and the mark moves
    /// forward.**
    ///
    /// A reference reader is a sliding window over the FASTA. It grows as it is asked for more
    /// and shrinks only when told to — so a walk that never tells it holds one byte for every
    /// base it has passed, about 250 MB on human chromosome 1. `ref_seq`'s own tests measure
    /// that directly; this pins the two things *this* generator is responsible for: that it
    /// releases at all, and that the mark advances with the walk.
    ///
    /// Reviewers found this exact call site deletable with the whole suite green, which is why
    /// it has a test of its own rather than sharing the generic generator's.
    #[test]
    fn every_repeat_releases_the_reference_behind_it() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        // Kept clear of the 100-base fixture contig's end: the last read runs 50..79.
        let records: Vec<_> = (0..4)
            .map(|index| read_named_with_length(&format!("r{index}"), 0, 5 + index * 15, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        // Two recorders: the margin fetch's reader, and the one the read filter holds inside
        // the cursor.
        let released = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let released_by_cursor = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let spy = ReleaseSpy {
            inner: fixture_ref_bases(),
            released: std::sync::Arc::clone(&released),
        };
        let make_reference = {
            let released_by_cursor = std::sync::Arc::clone(&released_by_cursor);
            move || ReleaseSpy {
                inner: fixture_ref_bases(),
                released: std::sync::Arc::clone(&released_by_cursor),
            }
        };
        let mut generator = SsrGenerator::new(
            spy,
            make_reference,
            SsrUnitSlipAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold");

        for index in 0..4u64 {
            let start = 20 + index * 10;
            generator.begin_segment(span(start, start + 9));
            generator
                .next_locus(&tract(start, start + 9), &reads)
                .expect("no fetch error");
        }

        let marks = released.lock().expect("no panic holds this").clone();
        assert_eq!(
            marks.len(),
            4,
            "one release per repeat, so a walk of any length holds one repeat's worth: \
             {marks:?}",
        );
        assert!(
            marks.windows(2).all(|pair| pair[1] > pair[0]),
            "the mark must move forward with the walk, or it releases the same nothing every \
             time: {marks:?}",
        );
        // The first repeat starts at 20 and the margin is one flank, 10.
        assert_eq!(
            marks[0], 10,
            "the mark is the repeat's start less one flank, which is the furthest back either \
             reader looks from a repeat's own start: {marks:?}",
        );

        // The cursor's reader has none to release at the first repeat, which opens it.
        let by_cursor = released_by_cursor
            .lock()
            .expect("no panic holds this")
            .clone();
        assert_eq!(
            by_cursor,
            marks[1..].to_vec(),
            "the cursor's reader is released too, from the repeat after the one that opens it: \
             {by_cursor:?}",
        );
    }

    /// `new` refuses a flank wider than the bundle threshold — the cross-config check at
    /// construction (spec §4).
    #[test]
    fn new_rejects_a_flank_wider_than_the_bundle_threshold() {
        // `SsrGenerator` is not `Debug` (it holds a closure and the aligner), so match on the
        // `Result` directly rather than `expect_err`.
        let result = SsrGenerator::with_default_aligner(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrGeneratorConfig {
                flank_bp: Bp(50),
                max_reads_per_locus: None,
            },
            Bp(10),
        );
        assert!(matches!(
            result,
            Err(SsrGeneratorConfigError::FlankExceedsBundleThreshold { .. })
        ));
    }

    /// A tract no read covers still yields **one** locus — present, with an empty observation
    /// table and zeroed drop counts — and the second `next_locus` returns `None`. "We looked and
    /// saw nothing" ≠ "we never looked" (spec §2). The locus carries the tract coordinates, the
    /// tract bases only, and the flanks split out of the fetched margin.
    #[test]
    fn a_zero_coverage_tract_yields_one_empty_but_present_locus() {
        // One read on chr1 at 65..94 — clear of the tract's query span (30..59), and long
        // enough (30 bp) to clear the min-read-length filter, so it is a genuine miss rather
        // than a filtered read.
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r0", 0, 65, 30)]);
        let mut generator = ssr_generator();
        let segment = tract(40, 49); // chr1, motif AC, 10 bp

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus per segment, even at zero coverage");

        assert_eq!(locus.region, span(40, 49));
        assert_eq!(locus.region.len(), 10);
        assert_eq!(
            &*locus.reference_bases,
            &b"AAAAAAAAAA"[..],
            "the tract only"
        );
        assert!(locus.observations.is_empty(), "no read covered it");
        assert_eq!(locus.reads_without_observation, 0);
        assert_eq!(locus.reads_discarded_by_cap, 0);
        match locus.kind {
            LocusKind::Ssr(SsrDetail {
                motif,
                left_flank,
                right_flank,
            }) => {
                assert_eq!(motif, Motif::new(b"AC").unwrap());
                assert_eq!(&*left_flank, &b"AAAAAAAAAA"[..], "10 bp left flank");
                assert_eq!(&*right_flank, &b"AAAAAAAAAA"[..], "10 bp right flank");
            }
            other => panic!("expected an Ssr locus, got {other:?}"),
        }

        assert!(
            generator
                .next_locus(&segment, &reads)
                .expect("no error")
                .is_none(),
            "the second poll ends the segment"
        );
        assert_eq!(
            generator.counts().reads_fetched,
            0,
            "no read reached the tract"
        );
    }

    /// A covered tract wires fetch → classify → tally: every fetched read is accounted for,
    /// either as an observation's support or in `reads_without_observation` (no cap, so nothing
    /// is discarded). The tract content here is not real STR sequence, so *which* bucket a read
    /// lands in is E's concern; this pins the conservation the four steps must uphold.
    #[test]
    fn a_covered_tract_accounts_every_fetched_read() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        // Four 30 bp reads (clearing the min-read-length filter) overlapping the query span
        // (30..59) on chr1.
        let records: Vec<_> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 30 + i * 3, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        let fetched = generator.counts().reads_fetched;
        assert_eq!(fetched, 4, "all four overlapping reads reached the tract");
        let supported: u32 = locus.observations.iter().map(|obs| obs.num_obs).sum();
        assert_eq!(
            supported as u64 + locus.reads_without_observation as u64,
            fetched,
            "every fetched read is either an observation or a no-observation"
        );
        assert_eq!(locus.reads_discarded_by_cap, 0, "no cap");
    }

    /// **`placed_left` is counted against the tract anchor, and this is the only test that says
    /// which coordinate that is.**
    ///
    /// The fold itself (`tally`) is handed `locus_start` as a bare `u64`, and its unit test
    /// supplies that argument itself — so the unit test pins the comparison (`<`, not `<=`) and
    /// nothing about *which* coordinate the generator passes. A 0-based/1-based slip, or the
    /// margin start, or a widened region's start would all be silent: `placed_left` is a QUAL
    /// bias term, so a wrong anchor is a wrong QUAL on every call with no panic.
    ///
    /// The fixture separates the two plausible anchors. The tract is `[40, 49]` and `flank_bp`
    /// is 10, so the margin starts at 30 — and the reads start at 30, 33, 36 and 39, every one
    /// of them **left of the tract anchor** and **not left of the margin anchor**. So the
    /// correct anchor counts them all and the margin anchor counts none.
    #[test]
    fn placed_left_is_counted_against_the_tract_anchor_not_the_margin() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 30 + i * 3, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        let supporting: u32 = locus.observations.iter().map(|obs| obs.num_obs).sum();
        let placed_left: u32 = locus.observations.iter().map(|obs| obs.placed_left).sum();

        assert!(
            supporting > 0,
            "the fixture must produce observations, or this test cannot discriminate"
        );
        assert_eq!(
            placed_left, supporting,
            "every supporting read started left of the tract anchor at 40; counting against \
             the margin start (30) would give 0"
        );
    }

    /// The mirror, so the anchor is pinned from **both** sides: reads starting at or after the
    /// tract anchor contribute nothing to `placed_left`. Together with the test above this
    /// brackets the anchor — one case fails if it is too far left, the other if it is too far
    /// right.
    #[test]
    fn a_read_starting_on_the_tract_anchor_is_not_placed_left() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        // Starting exactly on the anchor (40) and after it (43) — 30 bp each, so they still
        // cover the tract and clear the min-length filter.
        let records: Vec<_> = [40u64, 43]
            .iter()
            .enumerate()
            .map(|(i, start)| read_named_with_length(&format!("r{i}"), 0, *start as usize, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        let placed_left: u32 = locus.observations.iter().map(|obs| obs.placed_left).sum();
        assert_eq!(
            placed_left, 0,
            "strictly left, so the read starting *on* the anchor does not count — this is what \
             separates `placed_left` from the `placed_start` ng deliberately does not carry"
        );
    }

    /// **The length the mint clamps against is the length the emitted locus reports.**
    ///
    /// They come from different expressions — the mint uses `segment.tract_len()`, because it runs
    /// per read before the locus exists, while every consumer asks the locus
    /// (`SampleLocusObservations::locus_len`, over `region`). They agree by construction, since the
    /// region *is* `segment.start()..=segment.end()`, but nothing said so and a divergence would be
    /// silent: a run clamped against one length and tested for flushness against another reports
    /// the wrong side.
    #[test]
    fn the_mints_locus_length_is_the_emitted_regions() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records = vec![read_named_with_length("r0", 0, 30, 30)];
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        assert_eq!(
            locus.locus_len(),
            LocusLen::from_positions(segment.tract_len()),
            "the consumer's length and the mint's are the same quantity"
        );
    }

    /// The cap wires through `next_locus`: with `max_reads_per_locus = 2` over four overlapping
    /// reads, the locus reports two discarded, and the run-level counter agrees (the `as u32`
    /// path and `counts.reads_discarded_by_cap`).
    #[test]
    fn a_capped_locus_reports_the_discarded_reads() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 30 + i * 3, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = SsrGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrUnitSlipAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: Some(2),
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold");
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no error")
            .expect("one locus");

        assert_eq!(generator.counts().reads_fetched, 4);
        assert_eq!(locus.reads_discarded_by_cap, 2, "4 fetched, 2 kept");
        assert_eq!(generator.counts().reads_discarded_by_cap, 2);
    }

    /// A tract near the contig **end** splits a short right flank end-to-end: the fetched margin
    /// is clamped, so `SsrDetail.right_flank` is shorter than the left, and `reference_bases` is
    /// still the tract only. This is the contig-end-clamp path through the generator's
    /// flank-split (proven not to misalign; pinned here as a regression anchor).
    #[test]
    fn a_near_end_tract_splits_a_short_right_flank() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        // A read clear of the tract's query span (82..100) — zero coverage, so the locus is
        // present-but-empty and only the flank split is under test.
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r0", 0, 1, 30)]);
        let mut generator = ssr_generator(); // flank_bp = 10
        // chr1 is 100 bp; tract [92, 96] → margin [82, 100], right flank clamped to 100 - 96 = 4.
        let segment = tract(92, 96);

        generator.begin_segment(span(92, 96));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no error")
            .expect("one locus");

        assert_eq!(locus.region, span(92, 96));
        assert_eq!(
            &*locus.reference_bases,
            &b"AAAAA"[..],
            "the 5 bp tract only"
        );
        match locus.kind {
            LocusKind::Ssr(SsrDetail {
                left_flank,
                right_flank,
                ..
            }) => {
                assert_eq!(left_flank.len(), 10, "full left flank");
                assert_eq!(right_flank.len(), 4, "clamped right flank (100 - 96)");
            }
            other => panic!("expected an Ssr locus, got {other:?}"),
        }
    }

    // --- E2: the parity oracle — the port anchor -----------------------------

    /// **The north-star parity oracle** (spec §6, §9.3): every COMPLETE observation ng produces
    /// must match production's `SsrLocusObs.observed` byte for byte, in bases and count, **with the
    /// cap disabled** on a fixture shallower than any cap.
    ///
    /// Both sides run over the **same** reads and the **same** reference frame. ng runs its real
    /// per-read pipeline (`classify::classify_read`, the `SsrFlatGapAligner` over production's
    /// `PerQualityEmission` table) and its real `tally`; production runs its real delimiter
    /// (`delimit_read` over `HmmModel`) and its real `tally`. The two aligners are *different
    /// algorithms* — this is exactly the claim under test: that on the complete class of clean reads
    /// they delimit the identical tract. The generator's fetch/cap are not exercised here because,
    /// with the cap disabled, they only select reads and never change an observation's bases; that
    /// wiring is covered by D3 and the dump-tool fixture (E1). The production side reproduces
    /// production's *clean-read* classify path (extract → delimit → the Region tract), which for
    /// these Q40, fully-spanned reads is exactly what production's `classify_read` does — asserted,
    /// not assumed, via `flank_truncated` and the quality floor.
    ///
    /// Partial observations are ng's new behaviour with **no oracle**: this checks only that they
    /// *exist* (ng's classify keeps as a partial the read production's delimiter drops as
    /// border-off-end), not their bytes (spec §6).
    ///
    /// **Scope.** The fixture is all-Match reads of the reference allele and a one-unit *contraction*
    /// (`CACA`), so it pins parity on the complete class over the delimiters' contraction and
    /// exact-match behaviour. An *expansion* allele or a soft-clipped read cannot be added without
    /// leaving this regime: both drive the tract past the reference-sized window, tripping the
    /// long-allele **widening** recovery — a branch this production-side reproduction deliberately
    /// does not model (it asserts `!flank_truncated`). Widening-path and soft-clip parity are not
    /// covered here.
    #[test]
    fn ng_complete_observations_match_frozen_production_byte_for_byte() {
        use crate::ng::alignment::ssr_best_path_flat_gap::{
            SsrFlatGapAligner, ViterbiScratch as NgViterbiScratch,
        };
        use crate::ng::alignment::{PerQualityEmission, StutterModel};
        use crate::ng::locus_generation::ReadWitness;
        use crate::ng::read::aligned_read::AlignedRead;
        use crate::pileup::walker::CigarOp;
        // Frozen production oracle (called test-only, as the reservoir parity test does; ng does not
        // depend on production at run time).
        use crate::ssr::pileup::alignment::{
            Delimited, HmmModel, ViterbiScratch as ProdViterbiScratch, delimit_read,
        };
        use crate::ssr::pileup::footprint::{extract_region, flank_truncated, read_footprint};
        use crate::ssr::pileup::locus_tally::{QcCounts, ReadObs, tally as production_tally};
        use crate::ssr::types::{Locus, Motif as ProductionMotif};

        // A 6 bp G-flank + CACACA + 6 bp T-flank — production's own delimiter fixture frame.
        const FRAME: &[u8] = b"GGGGGGCACACATTTTTT";

        /// A reference-frame read: `seq` mapped all-Match at 1-based `pos`, Q40.
        fn mapped(seq: &[u8], pos: u64) -> AlignedRead {
            AlignedRead {
                qname: b"r".to_vec(),
                flag: 0,
                ref_id: 0,
                pos,
                mapq: 60,
                cigar: vec![CigarOp::Match(seq.len() as u32)],
                seq: seq.to_vec(),
                qual: vec![40u8; seq.len()],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                read_group: ReadGroupId(0),
            }
        }

        // Three complete reads of the reference allele (CACACA) and two of a clean shorter allele
        // (CACA, one CA unit deleted) — two distinct observed sequences, so the parity is not a
        // single-allele triviality. Plus one partially-covering read (left flank + two CA units, no
        // right flank) — the read ng keeps as a partial and production drops as border-off-end.
        let reads = [
            mapped(FRAME, 1),
            mapped(FRAME, 1),
            mapped(FRAME, 1),
            mapped(b"GGGGGGCACATTTTTT", 1),
            mapped(b"GGGGGGCACATTTTTT", 1),
            mapped(b"GGGGGGCACA", 1),
        ];

        // --- ng: real classify + real tally ---------------------------------
        let ng_locus = SsrLocus {
            segment: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"CA").unwrap(), 1.0)
                .unwrap(),
            tract_with_margin_bases: FRAME.into(),
            margin_start: Position(1),
        };
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let mut ng_scratch = NgViterbiScratch::new();
        let mut qual_buffer = Vec::new();
        let mut counts = SsrGeneratorCounts::default();
        let ng_outcomes: Vec<_> = reads
            .iter()
            .map(|read| {
                classify::classify_read(
                    read,
                    &ng_locus,
                    &aligner,
                    &stutter,
                    &mut ng_scratch,
                    &mut qual_buffer,
                )
            })
            .collect();
        let ng = tally::tally(
            reads.iter().zip(ng_outcomes),
            ng_locus.segment.start(),
            &mut counts,
        );
        let ng_complete: Vec<(Vec<u8>, u32)> = ng
            .observations
            .iter()
            .filter(|obs| obs.read_witness == ReadWitness::Complete)
            .map(|obs| (obs.bases.to_vec(), obs.num_obs))
            .collect();

        // --- production: real delimit_read + real tally ---------------------
        let production_locus = Locus::new(
            "chr1".into(),
            6, // 0-based tract start
            12,
            ProductionMotif::new(b"CA").unwrap(),
            1.0,
            FRAME.into(),
            0, // ref_bytes_start (0-based)
        )
        .unwrap();
        // The clean-read glue below reproduces production's `classify_read`
        // (src/ssr/pileup/driver.rs:195) — private, so it cannot be called directly. Only its
        // control flow is copied; the delimiter (`delimit_read`) and the `tally` are the real frozen
        // production code. Both simplifications the reproduction makes are **self-checking**, so it
        // cannot silently drift from production: the widening branch is asserted away
        // (`!flank_truncated`) and the quality gate is asserted to pass. If production's
        // `classify_read` grows a step before its Region→Sequence path, update here too.
        let model = HmmModel::new();
        let mut production_scratch = ProdViterbiScratch::new();
        let outcomes: Vec<ReadObs> = reads
            .iter()
            .map(|read| {
                let footprint = read_footprint(&read.cigar, read.pos);
                let region =
                    extract_region(&read.cigar, footprint, read.seq.len(), &production_locus);
                match delimit_read(
                    &read.seq[region.clone()],
                    &read.qual[region.clone()],
                    &production_locus,
                    &model,
                    &mut production_scratch,
                ) {
                    Delimited::Region(tract) => {
                        // The fixture reads span cleanly: a complete tract with full flanks, so
                        // production's classify_read takes its Region→Sequence path, not widening.
                        assert!(
                            !flank_truncated(
                                &region,
                                &tract,
                                read.seq.len(),
                                production_locus.left_flank().len(),
                                production_locus.right_flank().len(),
                            ),
                            "the fixture read is not window-truncated"
                        );
                        // Production gates on the tract's lower-quartile base quality
                        // (MIN_REGION_Q1 = 15); every tract base here is Q40, so its
                        // `sequence_or_low_quality` returns `Sequence`. Asserted so that skipping
                        // the gate cannot mask a divergence.
                        assert!(
                            read.qual[region.clone()][tract.clone()]
                                .iter()
                                .all(|&q| q >= 15),
                            "the fixture tract clears production's quality floor"
                        );
                        ReadObs::Sequence(read.seq[region][tract].into())
                    }
                    Delimited::BorderOffEnd => ReadObs::BorderOffEnd,
                }
            })
            .collect();
        let qc = QcCounts {
            depth: reads.len() as u32,
            n_filtered: 0,
            mapped_reads: reads.len() as u32,
        };
        let production = production_tally(&production_locus, &outcomes, qc);
        let production_observed: Vec<(Vec<u8>, u32)> = production
            .observed
            .iter()
            .map(|(bases, count)| (bases.to_vec(), *count))
            .collect();

        // --- the parity assertion -------------------------------------------
        assert_eq!(
            ng_complete,
            vec![(b"CACA".to_vec(), 2), (b"CACACA".to_vec(), 3)],
            "ng's two complete alleles, sorted by bytes with their counts"
        );
        assert_eq!(
            ng_complete, production_observed,
            "ng's complete observations match production's SsrLocusObs.observed byte for byte"
        );

        // ng's classify keeps a partial where production's delimiter drops the same read.
        assert_eq!(
            counts.observations_partial, 1,
            "ng's classify keeps the partially-covering read as a partial"
        );
        assert_eq!(
            production.n_border_off_end, 1,
            "production's delimiter drops that same read as border-off-end — no oracle for its bytes"
        );
    }
}
