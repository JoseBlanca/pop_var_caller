//! Where a read's tract length is recorded: one bucket per whole-repeat offset from the
//! reference tract's length, with the two ends absorbing everything beyond them.
//!
//! **Its own file so that the bucket cannot be built wrongly by the three files that use
//! it.** Rust privacy reaches a module's descendants, so a private field declared in
//! `ssr/mod.rs` is still reachable from `ssr/locus_offsets.rs`, `ssr/stratum_table.rs` and
//! `ssr/slippage.rs` — which are exactly the modules that will index an entry's counts with
//! it. As a **sibling** of those three, the invariant holds against all of them:
//! `OffsetBucket(200)` is `error[E0603]` there, where inside `ssr/mod.rs`'s subtree it
//! compiled and panicked at run time with an out-of-bounds index.

use std::fmt;

use super::WholeRepeatOffset;

/// The offsets an entry records either side of the reference length — `±OFFSET_HALF_RANGE`,
/// with the two end buckets absorbing everything beyond them.
///
/// **Measured to matter far less than it looks, and the reason is how the ends are scored.**
/// An end bucket's probability is the sum over every offset it absorbs, never the
/// probability of sitting exactly on the edge; with that rule a range of **±1**, against
/// loci whose alleles reach ±3, still returns the slippage level to within 0.05% and both
/// shares to within 0.002 (`spec/parameter_prepass_ssr.md` §4.1). What a narrow range costs
/// is the heterozygosity that falls out of the fitted genotype frequencies — 1.5% at ±1 —
/// which this path does not emit.
///
/// **Four is comfortable on real data**: the saturating end buckets took 0.89% of reads
/// across GRCh38's typed tracts and 0.14% across tomato's 138 million (spec §4.1).
///
/// **The width that decides the answer is [`ALLELE_OFFSET_LIMIT`](super::ALLELE_OFFSET_LIMIT),
/// not this one**, and it is the wider of the two — a relation
/// [`super`] holds by compile-time assertion rather than by prose.
pub const OFFSET_HALF_RANGE: i8 = 4;

/// How many buckets an entry holds: one per offset in `-OFFSET_HALF_RANGE..=OFFSET_HALF_RANGE`.
///
/// Derived rather than restated, so the two cannot drift apart.
pub const OFFSET_BUCKETS: usize = (2 * OFFSET_HALF_RANGE + 1) as usize;

/// The half-range must leave room for a bucket either side of the origin. A zero or
/// negative value would make [`OFFSET_BUCKETS`] one or astronomically large, and neither is
/// a table anything can index.
const _: () = assert!(
    OFFSET_HALF_RANGE > 0,
    "the recorded offset range needs at least one bucket either side of the reference length"
);

/// Which bucket of an entry an offset falls in: `0 ..= 2·OFFSET_HALF_RANGE`, with index
/// [`OFFSET_HALF_RANGE`] the reference length.
///
/// **The two end buckets saturate**, so bucket 0 means "at least `OFFSET_HALF_RANGE`
/// repeats short of the reference" rather than "exactly that many short". A scoring rule
/// that forgets this and plugs in the edge is not a probability distribution at all — its
/// buckets sum to 0.9488 at ±1 — and rescaled to sum to one it still costs +33% of the
/// slippage level where 30 in 100 slipped reads take a second step
/// (`spec/parameter_prepass_ssr.md` §4.1).
///
/// **Its value is private and [`bucket_of`] is the only way to make one**, because an
/// out-of-range bucket would index past an entry's counts. Any `i8` is a legal offset; only
/// `0..OFFSET_BUCKETS` is a legal bucket.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct OffsetBucket(u8);

impl OffsetBucket {
    /// The bucket's index into an entry's counts, always inside `0..OFFSET_BUCKETS`.
    #[inline]
    #[must_use]
    pub fn index(self) -> usize {
        usize::from(self.0)
    }

    /// Whether this is one of the two saturating ends — the buckets that hold more than one
    /// offset, and so the ones a scoring rule must reach by summing rather than by plugging
    /// in a single value.
    #[inline]
    #[must_use]
    pub fn is_saturating_end(self) -> bool {
        self.0 == 0 || self.index() == OFFSET_BUCKETS - 1
    }
}

impl fmt::Display for OffsetBucket {
    /// The offset the bucket sits at, signed, with the two ends marked as the open-ended
    /// ranges they are: `≤-4`, `-3` … `+3`, `≥+4`. A report that printed a bare `-4` for the
    /// low end would claim an exact offset the bucket does not hold.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let offset = self.0 as i8 - OFFSET_HALF_RANGE;
        match self.0 {
            0 => write!(f, "≤{offset:+}"),
            _ if self.is_saturating_end() => write!(f, "≥{offset:+}"),
            _ => write!(f, "{offset:+}"),
        }
    }
}

/// Which bucket an offset is recorded in, clamping to the saturating ends.
///
/// Total by construction: every `i8` maps to a bucket, and offsets beyond the range land in
/// the end that absorbs them.
#[inline]
#[must_use]
pub fn bucket_of(offset: WholeRepeatOffset) -> OffsetBucket {
    let clamped = offset.get().clamp(-OFFSET_HALF_RANGE, OFFSET_HALF_RANGE);
    // Distance from the low end of the clamped range, which is `0 ..= 2·OFFSET_HALF_RANGE`
    // by construction and needs no cast to prove it: `abs_diff` on two `i8`s returns `u8`.
    OffsetBucket(clamped.abs_diff(-OFFSET_HALF_RANGE))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two constants describe one range and a mismatch between them would size every
    /// entry wrongly, so this reads the derived one against the definition rather than
    /// against a literal.
    #[test]
    fn the_bucket_count_is_the_offset_range_it_is_derived_from() {
        assert_eq!(OFFSET_BUCKETS, 9);
        assert_eq!(
            OFFSET_BUCKETS,
            (-OFFSET_HALF_RANGE..=OFFSET_HALF_RANGE).count()
        );
    }

    /// Inside the range each offset gets its own bucket and the mapping rises with the
    /// offset — the reference length sits in the middle, at index `OFFSET_HALF_RANGE`.
    #[test]
    fn every_offset_inside_the_range_has_its_own_bucket_in_order() {
        let buckets: Vec<usize> = (-OFFSET_HALF_RANGE..=OFFSET_HALF_RANGE)
            .map(|offset| bucket_of(WholeRepeatOffset(offset)).index())
            .collect();

        assert_eq!(buckets, (0..OFFSET_BUCKETS).collect::<Vec<_>>());
        assert_eq!(
            bucket_of(WholeRepeatOffset(0)).index(),
            usize::try_from(OFFSET_HALF_RANGE).unwrap(),
            "the reference length is the middle bucket"
        );
    }

    /// **The end buckets absorb rather than reject**, which is the property the scoring rule
    /// has to know about: bucket 0 means "at least four repeats short", not "exactly four
    /// short". A read 40 repeats short is a real observation — a locus carrying a long
    /// deletion — and it has to land somewhere.
    #[test]
    fn offsets_past_the_range_saturate_into_the_end_buckets() {
        for far_short in [-5, -12, -40, i8::MIN] {
            assert_eq!(bucket_of(WholeRepeatOffset(far_short)).index(), 0);
        }
        for far_long in [5, 12, 40, i8::MAX] {
            assert_eq!(
                bucket_of(WholeRepeatOffset(far_long)).index(),
                OFFSET_BUCKETS - 1
            );
        }
    }

    /// **The two properties this file claims over a whole domain, checked over the whole
    /// domain rather than at sampled points**: every `i8` lands in a bucket that can index an
    /// entry's counts, and the mapping never runs backwards as the offset rises. 256 inputs,
    /// so the claim is settled rather than sampled — and from Milestone B on this index
    /// subscripts a fixed-size array, where leaving the range is an out-of-bounds panic.
    #[test]
    fn every_i8_offset_maps_into_the_bucket_range_without_going_backwards() {
        let mut previous = 0;
        for raw in i8::MIN..=i8::MAX {
            let index = bucket_of(WholeRepeatOffset(raw)).index();
            assert!(
                index < OFFSET_BUCKETS,
                "offset {raw} indexed past the counts, at {index}"
            );
            assert!(
                index >= previous,
                "offset {raw} went backwards, {previous} to {index}"
            );
            previous = index;
        }
    }

    /// Only the two ends hold more than one offset, and a scoring rule that treats an
    /// interior bucket as saturating — or an end bucket as exact — is wrong in the direction
    /// measured at +33% of the slippage level.
    #[test]
    fn exactly_the_two_ends_are_saturating() {
        let saturating: Vec<usize> = (-OFFSET_HALF_RANGE..=OFFSET_HALF_RANGE)
            .map(|offset| bucket_of(WholeRepeatOffset(offset)))
            .filter(|bucket| bucket.is_saturating_end())
            .map(OffsetBucket::index)
            .collect();

        assert_eq!(saturating, vec![0, OFFSET_BUCKETS - 1]);
    }

    /// A bucket renders as what it holds, and the two ends say they are open-ended: a report
    /// column printing a bare `-4` beside a `-3` would claim an exact offset for a bucket
    /// that holds every offset below it.
    #[test]
    fn the_end_buckets_render_as_ranges_and_the_interior_as_exact_offsets() {
        assert_eq!(bucket_of(WholeRepeatOffset(-40)).to_string(), "≤-4");
        assert_eq!(bucket_of(WholeRepeatOffset(-3)).to_string(), "-3");
        assert_eq!(bucket_of(WholeRepeatOffset(0)).to_string(), "+0");
        assert_eq!(bucket_of(WholeRepeatOffset(40)).to_string(), "≥+4");
    }
}
