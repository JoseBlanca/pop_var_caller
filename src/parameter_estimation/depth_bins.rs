//! The depth ladder a census writes each position's depth on: which depths share a code.
//!
//! **One code per depth up to 124 reads, then ten widening rungs up to about 1,500.** A census
//! keeps every position's allele counts exactly up to the per-position cap of 124 reads and thins
//! them to it above, so the ladder is exact wherever the record is and coarse only where the
//! record is already approximate. [`CENSUS_DEPTH_BIN_COUNT`] has what the exact region is worth
//! and what it costs; [`CENSUS_REACH_BINS`] has why the ladder reaches past the cap at all.
//!
//! **The rungs above the cap come from an older, twenty-bin ladder, and that is why three of its
//! constants are still here.** Before the census route existed, a run kept a histogram of
//! (depth, alternative-read count) cells whose memory grew with the square of the depth, and its
//! ladder was exact only to 8 reads, then widened in eleven geometric steps to the cap of 124
//! ([research note](../../../../doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md)
//! §4.3 has why twenty bins and that cap). That route is gone. The census ladder's ten rungs above
//! the cap continue at the same growth ratio, so [`EXACT_DEPTH_LIMIT`], [`MAX_BINNED_DEPTH`] and
//! [`DEPTH_BIN_COUNT`] are the named inputs to that ratio. **Changing any of them moves the
//! rungs, which changes the ladder's digest (`DepthLadderDigest`) and makes every census already
//! on disk unreadable.**
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2,
//! `doc/devel/ng/spec/parameter_prepass_joint_records.md` §4.1.

use std::ops::RangeInclusive;

use crate::float;

/// The deepest depth the old twenty-bin ladder gave a bin to itself (`0..=8`), and the base the
/// growth ratio of the census ladder's upper rungs is measured from.
///
/// **Fixed, not a knob**: one of the three inputs to [`widening_ratio`], so changing it moves
/// every rung above the cap (module doc). It also bounds `bin_for`'s fast path, which is correct
/// for any value at or below the cap because the census ladder has a bin per depth there.
const EXACT_DEPTH_LIMIT: u32 = 8;

/// The per-position depth cap: the deepest depth with a code of its own. A deeper position has
/// its allele counts thinned to this many reads; its depth is still written, on a rung above.
///
/// **Fixed, not a knob.** It is one of the three inputs to [`widening_ratio`] and the top of the
/// exact region.
///
/// **What losing the reads above costs is not measured.** The argument is that it costs
/// nothing — at 124 reads a heterozygote shows about 62 alternative reads and a
/// homozygous-reference site about 0.12, so the genotype is already certain and more
/// depth cannot make it more so — but no world in either research harness reaches the
/// cap, and research §4.6 lists this as the one depth mechanism shipping without a
/// measurement behind it. It fires on samples above ~124×, which is HG002 and never
/// tomato.
const MAX_BINNED_DEPTH: u32 = 124;

/// How many bins the old twenty-bin ladder had: nine exact ones (`0..=8`) and eleven widening to
/// the cap. **Kept only as an input to [`widening_ratio`]**: eleven steps from 8 to 124 fix the
/// ratio the census ladder's upper rungs grow by.
const DEPTH_BIN_COUNT: usize = 20;

/// How many rungs the **census** ladder has above the cap: ten, at the growth ratio the
/// twenty-bin ladder's widening bins were built with, carrying the top from 124 to about
/// 1,500 reads a position.
///
/// **The ceiling is what needs them, not the resolution.** A position a sample carries two
/// copies of holds about twice that sample's median depth, so the signal dies wherever twice
/// the median runs past what the encoding can say. At a top of 124 that is 76 reads a
/// position for a doubled position to stop reading deeper than an ordinary one, and 98 for
/// the two to be written identically — inside the depth range this caller commits to
/// (`spec/parameter_prepass_joint_records.md` §4.1).
const CENSUS_REACH_BINS: usize = 10;

/// How many bins the **census** ladder has: one per depth up to the cap of 124, then the
/// [ten reach rungs](CENSUS_REACH_BINS) above it.
///
/// **Every depth the counts are exact at has a bin to itself, and coarseness starts where the
/// counts stop being exact.** A position deeper than the cap has its allele counts thinned to
/// it, so above 124 the record is approximate whatever the ladder does; below it the record is
/// exact and the ladder was the only thing making it otherwise.
///
/// **What the coarse rungs cost was measured on contamination, which is the parameter that
/// reads a depth and a disagreeing-read count against each other.** On a drawn panel of 40
/// samples at 30 reads a position with one sample contaminated at 3%, the widening ladder — at
/// 30 reads its code stood for "between 29 and 36" — returned **0.0120**, and one bin per depth
/// returns **0.0263**, against a floor that stays at 0.0000 for all 39 clean samples either way
/// (`reports/census_depth_resolution_2026-08-16.md`). The count of disagreeing reads was always
/// exact, so a 3% contamination moving a read share by 1.5% sat inside a ±12% uncertainty that
/// came from nothing but how the depth was written down.
///
/// **It is not free, and the price is a third of a byte a position.** The code no longer fits in
/// five bits: 135 rungs plus the never-walked sentinel needs eight, which takes the depth array
/// from 1.25 MB to 2.0 MB per sample per read group at two million kept positions — 79 MB to
/// 126 MB across a cohort of 63. Beside the same sample's list of disagreeing reads (2.7 MB at
/// 30 reads a position) it is a few percent of the census.
pub(crate) const CENSUS_DEPTH_BIN_COUNT: usize = MAX_BINNED_DEPTH as usize + 1 + CENSUS_REACH_BINS;

/// The index of the first bin above [`EXACT_DEPTH_LIMIT`], where `bin_for`'s search starts.
const FIRST_WIDENING_BIN: usize = EXACT_DEPTH_LIMIT as usize + 1;

// The ratio needs at least one widening step and a cap above its base, or `widening_ratio`
// divides by zero steps or raises a ratio of one; and every bin index has to fit the `u16` a
// `DepthBin` holds. Checked at build time because these are consts a later edit can change.
const _: () = assert!(
    EXACT_DEPTH_LIMIT >= 1
        && DEPTH_BIN_COUNT > FIRST_WIDENING_BIN
        && MAX_BINNED_DEPTH > EXACT_DEPTH_LIMIT
        && CENSUS_DEPTH_BIN_COUNT - 1 <= u16::MAX as usize,
    "the depth ladder needs a widening step, a cap above its exact region, and u16 bin indices"
);

/// Which rung of the census ladder a position's depth falls in — **an index, not a depth**.
/// `u16` because the 135 rungs fit with room to spare.
///
/// Only `0..DepthBinEdges::bin_count()` is a legal value. The field is public because a stored
/// depth code is built from one and matched on; the accessors that take a bin check the range
/// rather than trusting it — see `DepthBinEdges::depth_range`.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct DepthBin(pub u16);

impl DepthBin {
    #[inline]
    #[must_use]
    pub fn get(self) -> u16 {
        self.0
    }
}

/// The census depth ladder: where one rung ends and the next begins.
///
/// Built by [`DepthBinEdges::for_census`], the only constructor, so every ladder a run holds is
/// the same one. Named for depths and not made generic over "any binned quantity" on purpose: a
/// repeat count and a depth are both `u32`, and edges that accepted either would let the two be
/// transposed silently.
#[derive(Debug)]
pub struct DepthBinEdges {
    /// The deepest depth in each bin, ascending and strictly increasing. `bin_tops[b]`
    /// is the top of bin `b`, so the ladder is stated once and everything else is
    /// derived from it.
    bin_tops: Vec<u32>,
}

/// How much wider each widening bin of the old twenty-bin ladder was than the one below it: the
/// ratio that climbs from `exact_limit` to `cap` in exactly `widening_bins` steps. The census
/// ladder's rungs above the cap grow by this same ratio ([`DepthBinEdges::for_census`]).
fn widening_ratio(exact_limit: u32, cap: u32, widening_bins: usize) -> f64 {
    float::powf(
        f64::from(cap) / f64::from(exact_limit),
        1.0 / widening_bins as f64,
    )
}

impl DepthBinEdges {
    /// **The ladder the census stores its depth codes on: one bin per depth to the cap of
    /// 124, then ten widening rungs to about 1,500 reads a position.**
    ///
    /// Exact wherever the record itself is exact, and coarse only above the cap, where a
    /// position's allele counts have been thinned and the record is approximate anyway. What
    /// the coarse version cost, and what the extra three bits buy, is on
    /// [`CENSUS_DEPTH_BIN_COUNT`]. The ten rungs above the cap grow by the twenty-bin ladder's
    /// ratio (module doc).
    #[must_use]
    pub fn for_census() -> Self {
        let mut bin_tops: Vec<u32> = (0..=MAX_BINNED_DEPTH).collect();
        let ratio = widening_ratio(
            EXACT_DEPTH_LIMIT,
            MAX_BINNED_DEPTH,
            DEPTH_BIN_COUNT - EXACT_DEPTH_LIMIT as usize - 1,
        );
        for step in 1..=CENSUS_REACH_BINS {
            let top =
                (f64::from(MAX_BINNED_DEPTH) * float::powi(ratio, step as i32)).round() as u32;
            // PANIC-FREE: the vector above always holds the cap.
            let previous = *bin_tops.last().expect("the ladder has at least one bin");
            bin_tops.push(top.max(previous + 1));
        }
        Self { bin_tops }
    }

    /// Which bin a depth falls in. **Total over every `u32`**: a depth above the ladder's top
    /// answers the last bin.
    #[must_use]
    pub fn bin_for(&self, depth: u32) -> DepthBin {
        // The shallowest depths answer without a search. This runs once per covered position
        // over hundreds of millions of them, and at three reads a position 97 in 100 land here.
        // A single `partition_point` over the whole ladder would give the same answer for
        // every depth — `bin_tops[d] == d` below the limit — at more comparisons.
        if depth <= EXACT_DEPTH_LIMIT {
            return DepthBin(depth as u16);
        }
        let widening = self.bin_tops[FIRST_WIDENING_BIN..].partition_point(|&top| top < depth);
        DepthBin((FIRST_WIDENING_BIN + widening).min(self.bin_count() - 1) as u16)
    }

    /// The depths this bin holds.
    ///
    /// A range rather than one endpoint, and `RangeInclusive` rather than a
    /// `(u32, u32)` pair, because a consumer sizes a denominator from it where an off-by-one
    /// silently mis-sizes the answer — inclusivity belongs in the type.
    #[must_use]
    pub fn depth_range(&self, bin: DepthBin) -> RangeInclusive<u32> {
        let index = bin.get() as usize;
        assert!(
            index < self.bin_count(),
            "depth bin {index} is outside a ladder of {} bins",
            self.bin_count()
        );
        let first = if index == 0 {
            0
        } else {
            self.bin_tops[index - 1] + 1
        };
        first..=self.bin_tops[index]
    }

    /// How many bins the ladder has.
    #[must_use]
    pub fn bin_count(&self) -> usize {
        self.bin_tops.len()
    }

    // **`recorded_depths` was here and is gone — 2026-08-14.** It answered `depth_range`
    // except at the ladder's top bin, where it answered the cap and nothing else, and that
    // was right only while the ladder's top and the per-position depth cap were the same
    // number: a position deeper than 124 was thinned to exactly 124 before it was written
    // down, so the top bin held one depth rather than the twenty-seven it spans. It was also
    // *wrong* for a position that genuinely held 100 reads, which was never thinned and yet
    // read back as 124.
    //
    // The census now stores the position's true depth
    // (`spec/parameter_prepass_joint_records.md` §2.2), so a code means its bin's whole range
    // and `depth_range` is the answer. What a consumer needs for the **allele counts'**
    // denominator is a different question, because those are still thinned: it is
    // `DepthCap::denominator_for`, in `joint/census.rs`, where the cap lives.

    /// Every bin, in order.
    ///
    /// It exists so that a caller walking the ladder does not write `DepthBin(bin as
    /// u16)` by hand: the cast is safe only because the const assertion at the top of
    /// this file pins `CENSUS_DEPTH_BIN_COUNT - 1` inside a `u16`, and that guarantee lives
    /// here rather than at the call sites.
    pub fn bins(&self) -> impl Iterator<Item = DepthBin> + '_ {
        // PANIC-FREE: the const assertion above pins the ladder's last index inside a
        // `u16`, and `for_census` is the only constructor.
        (0..self.bin_count()).map(|bin| DepthBin(bin as u16))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The deepest depth the census ladder has a bin for.
    fn ceiling(edges: &DepthBinEdges) -> u32 {
        *edges
            .depth_range(DepthBin(edges.bin_count() as u16 - 1))
            .end()
    }

    /// **The census ladder, pinned rung by rung.** One bin per depth to the cap, then ten rungs
    /// grown at the twenty-bin ladder's ratio. These literals are the oracle against a moved
    /// constant: a change to any of the three inputs to [`widening_ratio`] moves the upper
    /// rungs, changes the ladder's digest, and makes every census on disk unreadable.
    #[test]
    fn the_census_ladder_is_one_bin_per_depth_to_the_cap_then_ten_measured_rungs() {
        let census = DepthBinEdges::for_census();
        let tops: Vec<u32> = census
            .bins()
            .map(|bin| *census.depth_range(bin).end())
            .collect();

        assert_eq!(census.bin_count(), CENSUS_DEPTH_BIN_COUNT);
        assert_eq!(census.bin_count(), 135);
        assert_eq!(
            tops[..=MAX_BINNED_DEPTH as usize],
            (0..=MAX_BINNED_DEPTH).collect::<Vec<u32>>()[..],
            "the census ladder holds one bin per depth up to the cap"
        );
        assert_eq!(
            tops[MAX_BINNED_DEPTH as usize + 1..],
            [159, 204, 262, 336, 431, 553, 709, 910, 1168, 1498],
            "the ten rungs above the cap, at the ratio the twenty-bin ladder's widening bins used"
        );
    }

    /// **The half of this the extension exists for**, and the half a ladder that
    /// saturated early would still pass: a position a sample carries two copies of holds
    /// about twice that sample's median depth, so twice a depth has to land at least one
    /// rung above it or the two are written identically and the signal is gone. At the
    /// twenty-bin ladder that failed from 98 reads a position upwards.
    #[test]
    fn a_doubled_depth_lands_at_least_one_rung_higher_all_the_way_to_the_ceiling() {
        let census = DepthBinEdges::for_census();
        let ceiling = ceiling(&census);
        assert!(
            (1_400..=1_600).contains(&ceiling),
            "the census ladder tops out at {ceiling}, not near 1,500"
        );
        for depth in 4..=(ceiling / 2) {
            assert!(
                census.bin_for(2 * depth) > census.bin_for(depth),
                "{depth} reads a position and {} land in the same bin, so a doubled \
                 position cannot be told from an ordinary one there",
                2 * depth
            );
        }
    }

    /// `bin_for` is **monotone and total** over every depth the ladder covers: a
    /// deeper position never lands in an earlier bin, every bin is reachable, and the bin
    /// a depth answers is the bin whose range contains it. The fast path below
    /// [`EXACT_DEPTH_LIMIT`] and the search above it have to agree, and this is what checks
    /// they do.
    #[test]
    fn bin_for_is_monotone_and_total_over_every_depth_the_ladder_covers() {
        let edges = DepthBinEdges::for_census();
        let mut previous = edges.bin_for(0);
        let mut seen = vec![false; edges.bin_count()];
        seen[previous.get() as usize] = true;

        for depth in 0..=ceiling(&edges) {
            let bin = edges.bin_for(depth);
            assert!(
                bin >= previous,
                "depth {depth} went backwards to bin {bin:?}"
            );
            assert!(
                edges.depth_range(bin).contains(&depth),
                "depth {depth} is outside the range of the bin it answers, {:?}",
                edges.depth_range(bin)
            );
            seen[bin.get() as usize] = true;
            previous = bin;
        }

        assert!(
            seen.iter().all(|&reached| reached),
            "every bin is reachable"
        );
        assert_eq!(previous, DepthBin(edges.bin_count() as u16 - 1));
    }

    /// The bins **partition** the depths: consecutive ranges with no gap and no
    /// overlap, from 0 to the ceiling. Stated separately from `bin_for` because a rule
    /// that agrees with itself is not the same as a rule that covers everything —
    /// this is what says no depth is lost between two bins.
    #[test]
    fn the_bin_ranges_partition_the_depths_from_zero_to_the_ceiling() {
        let edges = DepthBinEdges::for_census();
        let mut expected_first = 0;

        for bin in 0..edges.bin_count() {
            let range = edges.depth_range(DepthBin(bin as u16));
            assert_eq!(
                *range.start(),
                expected_first,
                "gap or overlap at bin {bin}"
            );
            assert!(range.start() <= range.end(), "bin {bin} is inverted");
            expected_first = range.end() + 1;
        }

        assert_eq!(expected_first - 1, ceiling(&edges));
    }

    /// `bins()` yields every bin once, in order, and nothing else — it exists so that a
    /// caller walking the ladder does not spell the `u16` cast by hand, and a walk that
    /// skipped or repeated a bin would leave a rung out of the ladder's digest with nothing
    /// to say so.
    #[test]
    fn bins_yields_every_bin_once_in_order() {
        let edges = DepthBinEdges::for_census();
        let walked: Vec<DepthBin> = edges.bins().collect();

        assert_eq!(walked.len(), edges.bin_count());
        assert!(walked.windows(2).all(|pair| pair[0] < pair[1]));
        assert_eq!(walked.first().copied(), Some(DepthBin(0)));
        assert_eq!(
            walked.last().copied(),
            Some(DepthBin(edges.bin_count() as u16 - 1))
        );
    }

    /// A depth above the ladder's top answers the last bin rather than panicking.
    #[test]
    fn a_depth_above_the_ceiling_answers_the_last_bin() {
        let edges = DepthBinEdges::for_census();
        let last = DepthBin(edges.bin_count() as u16 - 1);

        assert_eq!(edges.bin_for(ceiling(&edges)), last);
        assert_eq!(edges.bin_for(ceiling(&edges) + 1), last);
        assert_eq!(edges.bin_for(u32::MAX), last);
    }

    /// **A bin index one past the end must fail, not answer**: a stored code past the ladder
    /// would otherwise read as some depth range.
    #[test]
    #[should_panic(expected = "outside a ladder of 135 bins")]
    fn depth_range_rejects_a_bin_past_the_end_of_the_ladder() {
        let edges = DepthBinEdges::for_census();
        let _ = edges.depth_range(DepthBin(edges.bin_count() as u16));
    }

    proptest::proptest! {
        /// `bin_for`'s doc claims totality over every `u32`, and the index it returns
        /// is stored as a depth code — so "always in range" is the invariant that matters.
        #[test]
        fn bin_for_answers_an_existing_bin_for_every_depth(depth in proptest::num::u32::ANY) {
            let edges = DepthBinEdges::for_census();
            let bin = edges.bin_for(depth);

            proptest::prop_assert!((bin.get() as usize) < edges.bin_count());
            let inside = edges.depth_range(bin).contains(&depth);
            proptest::prop_assert!(inside || depth > ceiling(&edges));
        }
    }
}
