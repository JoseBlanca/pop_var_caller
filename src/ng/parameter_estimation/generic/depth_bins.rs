//! The depth ladder: which depths share a bin.
//!
//! A site's depth is not kept exactly. Keeping every depth to 100 needs 5,151 cells
//! — `101 × 102 / 2`, since a site's alternative count runs from 0 to its depth —
//! which at eight bytes a cell is 330 MB per tomato sample, and the design rejects
//! that (`spec/parameter_prepass_generic.md` §9). Depth 100 and depth 105 say almost
//! the same thing about a genotype, while depth 2 and depth 3 say very different
//! things at three reads a plant, so the bins are one-per-depth at the bottom and
//! widen going up.
//!
//! **Twenty bins: exact integers to 8, then eleven geometrically widening bins to a
//! cap of 124.** That is not a memory choice, and the earlier draft that called it one
//! was wrong. Across twenty worlds — error-rate ratios of 1 and 4, mean depths 3 to
//! 60, even and 90/10 read splits — this ladder biases the fitted error rate by
//! **0.054 rungs** and each genotype frequency by **0.3%**, where sixteen bins at the
//! same cap costs **0.55 rungs and 1.8%**, and a cap of 300 at sixteen bins costs
//! **1.04 rungs and 8.0%**. Nothing downstream would show any of it
//! ([research note](../../../../doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md)
//! §4.3).
//!
//! **The bin count and the cap are one decision, not two.** A cap buys reach out of
//! the same twenty bins that buy resolution: raising it from 124 to 300 at a fixed bin
//! count doubles the error-rate bias and quadruples the one in the
//! homozygous-non-reference rate — measured on data where no site is deeper than 125,
//! so the extra reach is spent on depths nothing occupies and paid for out of the
//! depths everything occupies.
//!
//! **Where a ladder can hurt is 10 to 30 reads a site**, which is where an ordinary
//! whole-genome run sits. At three reads, 97 sites in 100 are at depth 6 or below and
//! never binned at all; at 60 the genotype is certain whatever the exact depth. Any
//! replacement must be checked in that band — one checked at tomato's own depth would
//! pass whatever it did.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2,
//! `doc/devel/ng/spec/parameter_prepass_generic.md` §4.

use std::ops::RangeInclusive;

/// How many depths get a bin to themselves, starting at zero: depths `0..=8` are
/// never merged with anything.
///
/// Eight rather than four or sixteen, and measured: an exact region of 4 costs 0.15
/// to 0.21 rungs even at tomato's three reads a site, because a three-read sample's
/// Poisson tail reaches depth 15; an exact region of 16 spends so many of the twenty
/// bins down here that the widening ones above become coarse enough to cost 0.98
/// rungs (research note §4.3).
pub const EXACT_DEPTH_LIMIT: u32 = 8;

/// The deepest a site is entered at. A deeper site is subsampled down to it
/// (Milestone C2), never rescaled.
///
/// Losing the reads above costs nothing: at 124 reads a heterozygote shows about 62
/// alternative reads and a homozygous-reference site about 0.12, so the genotype is
/// already certain and more depth cannot make it more so.
pub const MAX_BINNED_DEPTH: u32 = 124;

/// How many bins the ladder has: nine exact ones (`0..=8`) and eleven widening.
pub const DEPTH_BIN_COUNT: usize = 20;

/// Which depth bin a site's depth falls in — an index into a histogram's rows, not a
/// depth. `u16` because twenty bins fit with room to spare and the value is stored
/// once per cell key.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct DepthBin(pub u16);

impl DepthBin {
    #[inline]
    pub fn get(self) -> u16 {
        self.0
    }
}

/// The binning rule: where one bin ends and the next begins, plus where each bin's row
/// starts in a histogram's flat vector of cells.
///
/// **Built once per run and shared by every histogram in step 4**, so two accumulators
/// cannot drift apart and their cells stay comparable. Histograms hold it by `Arc`
/// rather than copying it, which is what lets `merge` prove two tables are binned the
/// same way by pointer identity instead of comparing lengths and hoping.
///
/// Named for depths and not made generic over "any binned quantity" on purpose: a
/// repeat count and a depth are both `u32`, and edges that accepted either would let
/// the two be transposed silently.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DepthBinEdges {
    /// The deepest depth in each bin, ascending. `bin_tops[b]` is the top of bin `b`,
    /// so the ladder is stated once and everything else is derived from it.
    bin_tops: Vec<u32>,
    /// Where bin `b`'s row starts in a histogram's flat cell vector. One longer than
    /// `bin_tops`, so the last entry is the total cell count.
    row_starts: Vec<usize>,
}

impl Default for DepthBinEdges {
    fn default() -> Self {
        Self::new()
    }
}

impl DepthBinEdges {
    /// The adopted ladder: exact integers to [`EXACT_DEPTH_LIMIT`], then geometrically
    /// widening bins to [`MAX_BINNED_DEPTH`], [`DEPTH_BIN_COUNT`] bins in all.
    ///
    /// The widening bins' tops are `EXACT_DEPTH_LIMIT · r^k` rounded, where `r` is the
    /// ratio that reaches the cap in exactly the number of bins left over. That gives
    /// tops of 10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124.
    pub fn new() -> Self {
        let widening_bins = DEPTH_BIN_COUNT - EXACT_DEPTH_LIMIT as usize - 1;
        let ratio = (f64::from(MAX_BINNED_DEPTH) / f64::from(EXACT_DEPTH_LIMIT))
            .powf(1.0 / widening_bins as f64);

        let mut bin_tops: Vec<u32> = (0..=EXACT_DEPTH_LIMIT).collect();
        for widening in 1..=widening_bins {
            let top = (f64::from(EXACT_DEPTH_LIMIT) * ratio.powi(widening as i32)).round() as u32;
            // Forced strictly increasing and capped, so a rounding collision at the
            // bottom of the widening region cannot produce an empty bin and the last
            // bin cannot overshoot the cap.
            let previous = *bin_tops.last().expect("the exact bins are never empty");
            bin_tops.push(top.max(previous + 1).min(MAX_BINNED_DEPTH));
        }

        // A bin's row must be as wide as its deepest site's alternative count, which
        // runs 0..=top — so `top + 1` cells.
        let mut row_starts = Vec::with_capacity(bin_tops.len() + 1);
        let mut cells = 0usize;
        for &top in &bin_tops {
            row_starts.push(cells);
            cells += top as usize + 1;
        }
        row_starts.push(cells);

        Self {
            bin_tops,
            row_starts,
        }
    }

    /// Which bin a depth falls in. **Total over every `u32`**: a depth above the cap
    /// answers the last bin, because a site deeper than the cap has been subsampled
    /// down to it before it reaches a histogram (Milestone C2) and a panic here would
    /// be a second, later guard on the same invariant.
    pub fn bin_for(&self, depth: u32) -> DepthBin {
        if depth <= EXACT_DEPTH_LIMIT {
            return DepthBin(depth as u16);
        }
        let widening =
            self.bin_tops[EXACT_DEPTH_LIMIT as usize + 1..].partition_point(|&top| top < depth);
        DepthBin((EXACT_DEPTH_LIMIT as usize + 1 + widening).min(self.bin_count() - 1) as u16)
    }

    /// Where this bin's row starts in a histogram's flat cell vector.
    pub fn row_start(&self, bin: DepthBin) -> usize {
        self.row_starts[bin.get() as usize]
    }

    /// The depths this bin holds.
    ///
    /// A range rather than one endpoint, and `RangeInclusive` rather than a
    /// `(u32, u32)` pair, because the sole consumer is a row width where an off-by-one
    /// silently mis-sizes the table — inclusivity belongs in the type.
    pub fn depth_range(&self, bin: DepthBin) -> RangeInclusive<u32> {
        let index = bin.get() as usize;
        let first = if index == 0 {
            0
        } else {
            self.bin_tops[index - 1] + 1
        };
        first..=self.bin_tops[index]
    }

    /// How many cells a histogram binned this way holds — the flat vector's length.
    pub fn cell_count(&self) -> usize {
        *self
            .row_starts
            .last()
            .expect("row_starts always carries its trailing total")
    }

    /// How many bins the ladder has.
    pub fn bin_count(&self) -> usize {
        self.bin_tops.len()
    }

    /// The deepest depth the ladder has a bin for. A site above it is subsampled down
    /// to it; nothing else in step 4 states this number, so the two cannot disagree.
    pub fn max_depth(&self) -> u32 {
        *self
            .bin_tops
            .last()
            .expect("the ladder always has at least one bin")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The measured ladder, pinned depth by depth.** These eleven numbers are the
    /// ladder the research note adopted, and the reason this step is its own commit:
    /// the same cap on sixteen bins biases the fitted error rate by 0.55 rungs and the
    /// homozygous-non-reference rate by 1.8%, against 0.05 and 0.3% here — and nothing
    /// downstream would show either. A `git bisect` over a moved parameter has to be
    /// able to land on one commit.
    #[test]
    fn the_widening_bins_top_out_at_the_measured_depths() {
        let edges = DepthBinEdges::new();
        let widening_tops: Vec<u32> = (EXACT_DEPTH_LIMIT as usize + 1..edges.bin_count())
            .map(|bin| *edges.depth_range(DepthBin(bin as u16)).end())
            .collect();

        assert_eq!(
            widening_tops,
            vec![10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124],
            "the adopted ladder (research note §4.3)"
        );
    }

    /// Depths 0 to 8 each get a bin to themselves — the part of the ladder that
    /// carries tomato's three-reads-a-site cohort, where 97 sites in 100 sit at depth
    /// 6 or below and are never merged with anything.
    #[test]
    fn the_bottom_of_the_ladder_is_one_bin_per_depth() {
        let edges = DepthBinEdges::new();
        for depth in 0..=EXACT_DEPTH_LIMIT {
            assert_eq!(edges.bin_for(depth), DepthBin(depth as u16));
            assert_eq!(edges.depth_range(DepthBin(depth as u16)), depth..=depth);
        }
        assert_eq!(edges.bin_count(), DEPTH_BIN_COUNT);
    }

    /// **583 cells**, because a bin's row must be as wide as its deepest site's
    /// alternative count: 45 in the nine exact bins and 538 in the eleven widening
    /// ones. It is the number `spec/parameter_prepass_generic.md` §9 prices the
    /// accumulator's memory against, so it is stated here rather than left to be
    /// counted.
    #[test]
    fn the_ladder_holds_583_cells() {
        let edges = DepthBinEdges::new();

        let exact_cells: usize = (0..=EXACT_DEPTH_LIMIT).map(|d| d as usize + 1).sum();
        assert_eq!(exact_cells, 45);
        assert_eq!(edges.cell_count(), 583);
        assert_eq!(edges.cell_count() - exact_cells, 538);
    }

    /// `bin_for` is **monotone and total** over every depth the ladder covers: a
    /// deeper site never lands in an earlier bin, every bin is reachable, and the bin
    /// a depth answers is the bin whose range contains it. A ladder that skipped a bin
    /// would leave a row of the histogram permanently empty, which no fit would
    /// report.
    #[test]
    fn bin_for_is_monotone_and_total_over_every_depth_the_ladder_covers() {
        let edges = DepthBinEdges::new();
        let mut previous = edges.bin_for(0);
        let mut seen = vec![false; edges.bin_count()];
        seen[previous.get() as usize] = true;

        for depth in 0..=MAX_BINNED_DEPTH {
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
    /// overlap, from 0 to the cap. Stated separately from `bin_for` because a rule
    /// that agrees with itself is not the same as a rule that covers everything —
    /// this is what says no depth is lost between two bins.
    #[test]
    fn the_bin_ranges_partition_the_depths_from_zero_to_the_cap() {
        let edges = DepthBinEdges::new();
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

        assert_eq!(expected_first - 1, MAX_BINNED_DEPTH);
        assert_eq!(edges.max_depth(), MAX_BINNED_DEPTH);
    }

    /// Each bin's row is exactly as wide as its deepest site's alternative count
    /// allows, and the rows sit end to end with no gap. An off-by-one here would have
    /// one bin's cells reading into the next bin's row — a wrong count with no symptom.
    #[test]
    fn each_row_is_as_wide_as_its_bins_deepest_alternative_count() {
        let edges = DepthBinEdges::new();

        for bin in 0..edges.bin_count() {
            let bin = DepthBin(bin as u16);
            let width = *edges.depth_range(bin).end() as usize + 1;
            let next_start = edges.row_start(bin) + width;
            let expected = if bin.get() as usize + 1 == edges.bin_count() {
                edges.cell_count()
            } else {
                edges.row_start(DepthBin(bin.get() + 1))
            };
            assert_eq!(next_start, expected, "row width wrong at bin {bin:?}");
        }
        assert_eq!(edges.row_start(DepthBin(0)), 0);
    }

    /// A depth above the cap answers the last bin rather than panicking. The cap is
    /// enforced by the subsampling that runs before a site reaches a histogram
    /// (Milestone C2); a panic here would be a second guard on the same invariant, and
    /// the one that fires would be the less informative of the two.
    #[test]
    fn a_depth_above_the_cap_answers_the_last_bin() {
        let edges = DepthBinEdges::new();
        let last = DepthBin(edges.bin_count() as u16 - 1);

        assert_eq!(edges.bin_for(MAX_BINNED_DEPTH), last);
        assert_eq!(edges.bin_for(MAX_BINNED_DEPTH + 1), last);
        assert_eq!(edges.bin_for(u32::MAX), last);
    }
}
