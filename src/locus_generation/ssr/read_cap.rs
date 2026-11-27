//! The per-locus read fetch and its depth cap.
//!
//! [`fetch_capped_reads`] fetches the reads over a locus's tract plus margin and, when a cap is
//! set, keeps at most that many with a reservoir sampler seeded from the locus itself — so the
//! kept set is the same at any thread count and in any run. The seed, the random-number
//! generator and the sampler are copied from the caller this one replaced, and
//! `ng_seed_and_reservoir_match_frozen_production_byte_for_byte` pins the kept sets that caller
//! produced.

use crate::locus_generation::LocusGenerationError;
use crate::read::aligned_read::AlignedRead;
use crate::read::input::cursor::CursorError;
use crate::read::input::sample_cursor::SampleCursor;
use crate::ref_seq::RawRefSeq;
use crate::region_typing::segment_criteria::SsrSegment;
use crate::types::GenomeRegion;

/// A tiny deterministic random-number generator — seeded per locus so the depth-cap subsample
/// is reproducible and thread-count-invariant, with no external generator whose stream could
/// shift.
///
/// **SplitMix64 in shape, but not the standard SplitMix64.** Its increment,
/// `0x9e3779b97f4a7615`, differs in one hex digit from the standard constant
/// `0x9e37_79b9_7f4a_7c15` that [`crate::splitmix`] uses. The value is still odd, so the
/// sequence is sound; it is kept as copied from the caller this one replaced because changing
/// it changes which reads the cap keeps, and
/// `ng_seed_and_reservoir_match_frozen_production_byte_for_byte` pins those reads. Do not
/// replace it with the shared function.
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
pub(super) fn seed_for_segment(segment: &SsrSegment) -> u64 {
    seed_for_repeat(segment.chrom(), segment.start())
}

/// The reservoir seed for any repeat starting at 1-based `start` on the contig named `chrom` —
/// [`seed_for_segment`]'s arithmetic, for an STR bundle's hull, which carries no contig name of
/// its own. The same trap and the same discharge: the **name**, and the **0-based** start.
pub(super) fn seed_for_repeat(chrom: &str, start: u64) -> u64 {
    locus_seed(chrom, (start - 1) as u32)
}

/// Reservoir sampler (Algorithm R) — an effectively-uniform sample of up to `capacity` items
/// from a stream of unknown length, in one pass with `O(capacity)` memory. Ported verbatim
/// from production: the eviction index `next_u64() % seen` carries a modulo bias bounded by
/// `seen / 2^64` (negligible at any real depth), accepted deliberately because the draw is
/// deterministic and thread-count-invariant — an unbiased reduction would change the kept set
/// and break that. The caller must `offer` admitted reads in a fixed total order —
/// `SampleReads`' merge order (spec §4).
pub(super) struct Reservoir<T> {
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

    #[cfg(test)]
    /// The admitted depth — total items offered (the reservoir sees only admitted reads).
    pub fn seen(&self) -> u64 {
        self.seen
    }

    /// Consume the reservoir, yielding the sampled items (≤ `capacity`).
    pub fn into_held(self) -> Vec<T> {
        self.held
    }
}

/// The reads kept for one locus after depth-capping, and how many were fetched — the caller
/// records `reads_discarded_by_cap = fetched - kept.len()`.
pub(super) struct CappedReads {
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
/// **Pointed at a cursor, not handed a query.** The cursor covers `query_span`'s
/// chromosome and stays open across every locus on it, so consecutive loci that share reads
/// decode them once (`spec/alignment_cursor.md` §1). Making the cursor and moving it to the
/// right chromosome is the caller's job; this only points it at the span.
pub(super) fn fetch_capped_reads<R: RawRefSeq>(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Motif;
    use std::collections::HashSet;

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
        // 1-based start 101 → 0-based 100.
        let segment =
            SsrSegment::new("chr1".into(), 101, 110, Motif::new(b"AC").unwrap(), 1.0).unwrap();
        assert_eq!(seed_for_segment(&segment), locus_seed("chr1", 100));
        assert_ne!(
            seed_for_segment(&segment),
            locus_seed("chr1", 101),
            "the 1-based start is the trap the conversion avoids"
        );
    }

    /// **The parity oracle for the port**: ng's seed and reservoir must produce output
    /// identical, byte for byte, to the deleted caller's `src/ssr/pileup/fetch_reads.rs` as of
    /// commit `d9e7b076`. The self-consistency tests above would survive a drifted constant;
    /// this one would not — it is what makes "byte-faithful port" a checked claim rather than an
    /// asserted one.
    ///
    /// **That caller's answers are frozen here**, since it was deleted on 2026-09-13: its
    /// `locus_seed` for the four loci and the kept set its `Reservoir` held, recorded at
    /// commit `d9e7b076`.
    #[test]
    fn ng_seed_and_reservoir_match_frozen_production_byte_for_byte() {
        for (chrom, start, production_seed) in [
            ("chr1", 0u32, 0x17ca_5576_8bd0_a6ad_u64),
            ("chr1", 100, 0x17ca_b176_8bd1_4301),
            ("chrX", 4242, 0x16ab_7176_8aea_2282),
            ("scaffold_7", 999_999, 0xb4f8_cd00_1980_491e),
        ] {
            assert_eq!(
                locus_seed(chrom, start),
                production_seed,
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
        // Production's `Reservoir` with the same capacity, seed and offers kept these.
        const PRODUCTION_KEPT: [u32; 8] = [9794, 3095, 2292, 9898, 13, 764, 4816, 3616];
        let seed = locus_seed("chrX", 7);
        assert_eq!(
            seed, 0x169b_de76_8acf_abb9,
            "production's seed for (chrX, 7)"
        );
        assert_eq!(
            kept(seed),
            PRODUCTION_KEPT,
            "the kept set must be identical"
        );
    }
}
