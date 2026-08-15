//! The depth ladder: which depths share a bin.
//!
//! A site's depth is not kept exactly. Keeping every depth to 100 needs 5,151 cells
//! per window — `101 × 102 / 2`, since a site's alternative count runs from 0 to its
//! depth — and at eight bytes a cell across a tomato sample's ~8,000 windows that is
//! the 330 MB the design rejects (`spec/parameter_prepass_generic.md` §9). Depth 100
//! and depth 105 say almost the same thing about a genotype, while depth 2 and depth 3
//! say very different things at three reads a plant, so the bins are one-per-depth at
//! the bottom and widen going up.
//!
//! **Twenty bins: exact integers to 8, then eleven geometrically widening bins to a
//! cap of 124.** That is not a memory choice, and the earlier draft of the
//! architecture that called it one was wrong. Across twenty worlds — error-rate ratios
//! of 1 and 4, mean depths 3 to 60, even and 90/10 read splits — this ladder biases
//! the fitted error rate by **0.054 rungs** of the quarter-Phred error-rate ladder and
//! each genotype frequency by **0.3%**, where sixteen bins at the same cap costs
//! **0.55 rungs and 1.8%**. Nothing downstream would show either
//! ([research note](../../../../doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md)
//! §4.3).
//!
//! **The bin count and the cap are one decision, not two.** A cap buys reach out of
//! the same twenty bins that buy resolution. Holding the bin count at twenty and
//! raising the cap from 124 to 300 takes the error-rate bias from 0.054 to 0.190 rungs
//! and the homozygous-non-reference rate from 0.30% to 0.88% — 3.5-fold and 2.9-fold,
//! measured on data where no site is deeper than 125, so the extra reach is spent on
//! depths nothing occupies and paid for out of the depths everything occupies. (At
//! sixteen bins the same change costs more still, 0.55 → 1.04 rungs and 1.8% → 8.0%,
//! which is the comparison the architecture's §2.2 quotes.)
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

/// The deepest depth that keeps a bin to itself: depths `0..=8` are never merged with
/// anything, which is nine bins.
///
/// **Fixed, not a knob**, and measured. An exact region reaching only 4 costs 0.15 to
/// 0.21 rungs of the error-rate ladder even at tomato's three reads a site, because a
/// three-read sample's Poisson tail reaches depth 15. One reaching 16 spends so many
/// of the twenty bins down here that the widening ones above go coarse enough to cost
/// 0.98 rungs (research note §4.3).
const EXACT_DEPTH_LIMIT: u32 = 8;

/// The deepest depth the ladder has a bin for. A deeper site is subsampled down to it
/// (Milestone C2), never rescaled.
///
/// **Fixed, not a knob**, and one decision with [`DEPTH_BIN_COUNT`] rather than an
/// independent knob — the module doc has what raising it costs.
///
/// **Private, and read through [`DepthBinEdges::max_depth`].** Two statements of this
/// number is the shape of a bug the architecture records as having already shipped once
/// in the design: an earlier draft set the subsampling cap to 300 while the ladder
/// ended at 124, which left depths 125–300 with no bin and made every cell count in two
/// documents wrong.
///
/// **What losing the reads above costs is not measured.** The argument is that it costs
/// nothing — at 124 reads a heterozygote shows about 62 alternative reads and a
/// homozygous-reference site about 0.12, so the genotype is already certain and more
/// depth cannot make it more so — but no world in either research harness reaches the
/// cap, and research §4.6 lists this as the one depth mechanism shipping without a
/// measurement behind it. It fires on samples above ~124×, which is HG002 and never
/// tomato.
const MAX_BINNED_DEPTH: u32 = 124;

/// How many bins the ladder has: nine exact ones (`0..=8`) and eleven widening.
///
/// **Fixed, not a knob.** Twenty against sixteen at the same cap is the difference
/// between 0.054 rungs and 0.55, for 30% more cells — a tenfold change in the answer
/// from a choice an earlier draft of the architecture left open on the premise that it
/// bought only memory.
const DEPTH_BIN_COUNT: usize = 20;

/// How many bins the **census** ladder has: the twenty above, plus ten more rungs at the
/// same growth ratio, carrying the top from 124 to about 1,500 reads a position.
///
/// **A longer ladder, not a different one.** The first twenty rungs are the twenty-bin
/// ladder's own tops, appended to rather than recomputed, so every edge of the shorter
/// ladder is an edge of this one and a census code maps to a histogram bin by collapsing
/// everything above 124. Two independently-derived ladders would let the two routes bin
/// differently, which is the thing the shared ladder exists to prevent.
///
/// **The ceiling is what needed the extension, not the resolution.** A position a sample
/// carries two copies of holds about twice that sample's median depth, so the signal dies
/// wherever twice the median runs past what the encoding can say. At a top of 124 that is
/// 76 reads a position for a doubled position to stop reading deeper than an ordinary one,
/// and 98 for the two to be written identically — inside the depth range this caller
/// commits to (`spec/parameter_prepass_joint_records.md` §4.1).
///
/// **It is free.** The stored code is five bits, which holds 32 values; twenty bins plus the
/// never-walked sentinel used 21 of them and left eleven spare. Thirty bins plus the
/// sentinel is 31.
pub(crate) const CENSUS_DEPTH_BIN_COUNT: usize = 30;

/// The index of the first widening bin — the one just above the exact region. Named
/// once, because the offset is otherwise re-derived at each site that needs it, and two
/// derivations that drift apart return a bin from the wrong region.
const FIRST_WIDENING_BIN: usize = EXACT_DEPTH_LIMIT as usize + 1;

// The ladder needs an exact region, at least one widening bin above it, and a cap above
// that region — otherwise `ladder_tops` divides by zero bins, or raises a ratio of
// infinity, or cannot reach its own cap. Checked at build time for the reason
// `generic/mod.rs` gives for the error-rate ladder's constants: these are consts a later
// edit can change, and every accessor here rests on them.
const _: () = assert!(
    EXACT_DEPTH_LIMIT >= 1
        && DEPTH_BIN_COUNT > FIRST_WIDENING_BIN
        && MAX_BINNED_DEPTH > EXACT_DEPTH_LIMIT
        && DEPTH_BIN_COUNT - 1 <= u16::MAX as usize,
    "the depth ladder needs an exact region, at least one widening bin, and a cap above it"
);

/// Which depth bin a site's depth falls in — **an index into a histogram's rows, not a
/// depth**. `u16` because twenty bins fit with room to spare and the value is stored
/// once per cell key.
///
/// Only `0..DepthBinEdges::bin_count()` is a legal value. The field is public because
/// the cell key is built from one and matched on; the accessors that take a bin check
/// the range rather than trusting it — see [`DepthBinEdges::row_start`].
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct DepthBin(pub u16);

impl DepthBin {
    #[inline]
    #[must_use]
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
/// **Deliberately neither `PartialEq` nor `Clone`, and both for the same reason.** Two
/// ladders are "the same rule" only when they are the same `Arc`, which is what
/// [`DepthAltHistogram::merge`](super::histogram::DepthAltHistogram::merge) checks. A
/// derived `==` would compare two vectors both derived from the same compile-time
/// constants, so it would answer `true` for any two ladders that exist at all — a check
/// that cannot fail, sitting one keystroke away from the check that must not be skipped.
/// And a `Clone` would let `Arc::new((*edges).clone())` hand a second worker a ladder
/// identical in every observable way and yet not the shared one, so that every later
/// merge against it panics at run time. Without the derive that line does not compile,
/// which is the earlier failure. **Clone the handle instead — `Arc::clone(&edges)`.**
/// There is nothing a separate ladder could be: [`DepthBinEdges::new`] takes no
/// arguments, so every ladder a run can build is byte-identical to every other.
///
/// Named for depths and not made generic over "any binned quantity" on purpose: a
/// repeat count and a depth are both `u32`, and edges that accepted either would let
/// the two be transposed silently.
#[derive(Debug)]
pub struct DepthBinEdges {
    /// The deepest depth in each bin, ascending and strictly increasing. `bin_tops[b]`
    /// is the top of bin `b`, so the ladder is stated once and everything else is
    /// derived from it.
    bin_tops: Vec<u32>,
    /// Where bin `b`'s row starts in a histogram's flat cell vector. One longer than
    /// `bin_tops`: the trailing entry is the total cell count, which is what
    /// [`DepthBinEdges::cell_count`] reads and what
    /// [`DepthBinEdges::row_start`] must refuse to return.
    row_starts: Vec<usize>,
}

impl Default for DepthBinEdges {
    fn default() -> Self {
        Self::new()
    }
}

/// The bin tops for a given ladder shape: strictly increasing, and ending exactly at
/// the cap.
///
/// **Parameterised, and private, so the clamp below is reachable from a test.** Under
/// the adopted constants it never fires — the raw geometric tops already land strictly
/// increasing and hit 124 exactly — so without a way to drive other shapes, the guard
/// against a degenerate ladder would be the one branch in this file that nothing
/// exercises, on the parameter this whole step was isolated to protect.
///
/// **Strict increase is guaranteed by construction, and the cap is what gets checked.**
/// Clamping to the cap *before* forcing each top above the last is what makes the
/// difference: two bins sharing a top would give the second an inverted range that
/// `bin_for` can never answer — a bin that exists in `bin_count()` and owns cells in
/// the flat vector, but that no site can land in — and applying the clamps the other
/// way round produces exactly that. Since `previous + 1` always exceeds `previous`, the
/// sequence cannot repeat; what it can do instead is walk past the cap, which is the
/// one thing left to assert. A ladder that ends above its cap leaves the depths between
/// its real top and the stated cap without a bin.
/// How much wider each widening bin is than the one below it: the ratio that climbs from
/// the exact region's top to the cap in exactly the bins left over.
///
/// **Named once because two users must not derive it separately.** The census ladder is the
/// twenty-bin one with ten more rungs *at this same ratio* ([`DepthBinEdges::for_census`]),
/// and a second derivation that drifted in the last decimal would round a shared rung
/// differently and break the property that makes the two ladders one.
fn widening_ratio(exact_limit: u32, cap: u32, widening_bins: usize) -> f64 {
    (f64::from(cap) / f64::from(exact_limit)).powf(1.0 / widening_bins as f64)
}

fn ladder_tops(exact_limit: u32, cap: u32, bin_count: usize) -> Vec<u32> {
    let widening_bins = bin_count - exact_limit as usize - 1;
    let ratio = widening_ratio(exact_limit, cap, widening_bins);

    let widening = (1..=widening_bins).scan(exact_limit, |previous, step| {
        let top = (f64::from(exact_limit) * ratio.powi(step as i32)).round() as u32;
        *previous = top.min(cap).max(*previous + 1);
        Some(*previous)
    });
    let tops: Vec<u32> = (0..=exact_limit).chain(widening).collect();

    // PANIC-FREE: `(0..=exact_limit)` yields at least one element for every `u32`, so
    // the chained vector is never empty.
    let highest = *tops
        .last()
        .expect("the exact region always contributes at least one bin");
    assert!(
        highest == cap,
        "a ladder of {bin_count} bins over depths 0..={cap}, exact to {exact_limit}, \
         tops out at {highest} instead of its cap, so the depths between would have no \
         bin — there are more bins here than the cap has room for"
    );
    tops
}

impl DepthBinEdges {
    /// The adopted ladder: exact integers to 8, then eleven geometrically widening bins
    /// to a cap of 124 — twenty bins, 583 cells.
    ///
    /// The widening bins' tops are `8 · r^k` rounded, where `r` is the ratio that
    /// reaches the cap in exactly the bins left over, giving 10, 13, 17, 22, 28, 36,
    /// 46, 59, 75, 97, 124.
    #[must_use]
    pub fn new() -> Self {
        Self::from_tops(ladder_tops(
            EXACT_DEPTH_LIMIT,
            MAX_BINNED_DEPTH,
            DEPTH_BIN_COUNT,
        ))
    }

    /// The same ladder with ten more rungs on top — thirty bins, topping out near 1,500
    /// reads a position. **This is the one the census stores its depth codes on.**
    ///
    /// **The twenty-bin ladder's rungs are taken, not recomputed.** The extension appends to
    /// them at the ratio they were built with, so every edge of the shorter ladder is an edge
    /// of this one by construction and not by a rounding coincidence — which is what lets a
    /// census code be turned into a histogram bin by collapsing everything above 124, and
    /// what stops the two routes binning differently. Rebuilding a thirty-rung geometric
    /// ladder from scratch over the wider range would move the shared rungs by a base or two
    /// and quietly break that.
    ///
    /// Why the reach is needed rather than the resolution: [`CENSUS_DEPTH_BIN_COUNT`].
    #[must_use]
    pub fn for_census() -> Self {
        let mut bin_tops = ladder_tops(EXACT_DEPTH_LIMIT, MAX_BINNED_DEPTH, DEPTH_BIN_COUNT);
        let ratio = widening_ratio(
            EXACT_DEPTH_LIMIT,
            MAX_BINNED_DEPTH,
            DEPTH_BIN_COUNT - EXACT_DEPTH_LIMIT as usize - 1,
        );
        for step in 1..=(CENSUS_DEPTH_BIN_COUNT - DEPTH_BIN_COUNT) {
            let top = (f64::from(MAX_BINNED_DEPTH) * ratio.powi(step as i32)).round() as u32;
            // PANIC-FREE: `ladder_tops` never returns an empty vector.
            let previous = *bin_tops.last().expect("the ladder has at least one bin");
            bin_tops.push(top.max(previous + 1));
        }
        Self::from_tops(bin_tops)
    }

    /// Both ladders' shared tail: the row offsets a histogram indexes by.
    fn from_tops(bin_tops: Vec<u32>) -> Self {
        // A bin's row must be as wide as its deepest site's alternative count, which
        // runs 0..=top — so `top + 1` cells. The leading 0 is bin 0's start; the
        // trailing entry is the total, which `cell_count` reads.
        let row_starts: Vec<usize> = std::iter::once(0)
            .chain(bin_tops.iter().scan(0usize, |cells, &top| {
                *cells += top as usize + 1;
                Some(*cells)
            }))
            .collect();

        Self {
            bin_tops,
            row_starts,
        }
    }

    /// Which bin a depth falls in. **Total over every `u32`**: a depth above the cap
    /// answers the last bin, because a site deeper than the cap has been subsampled
    /// down to it before it reaches a histogram (Milestone C2), and a panic here would
    /// be a second, later guard on the same invariant.
    #[must_use]
    pub fn bin_for(&self, depth: u32) -> DepthBin {
        // The exact region answers without a search, and that is a deliberate fast
        // path rather than a shortcut that happens to work: it is where 97 sites in
        // 100 of a three-read cohort land, and this runs once per covered position
        // over hundreds of millions of them. A single `partition_point` over the whole
        // ladder would give the same answer for every depth — `bin_tops[d] == d` below
        // the limit — at four or five comparisons instead of one.
        if depth <= EXACT_DEPTH_LIMIT {
            return DepthBin(depth as u16);
        }
        let widening = self.bin_tops[FIRST_WIDENING_BIN..].partition_point(|&top| top < depth);
        DepthBin((FIRST_WIDENING_BIN + widening).min(self.bin_count() - 1) as u16)
    }

    /// Where this bin's row starts in a histogram's flat cell vector.
    ///
    /// **The range check is not ceremony.** `row_starts` carries a trailing total, so
    /// without it a bin index one past the end would read that total and hand back
    /// `cell_count()` — a number shaped exactly like a row offset, one past the end of
    /// the vector the caller is about to index. In a module whose whole difficulty is
    /// that a wrong number here has no symptom, that is the worst available failure.
    #[must_use]
    pub fn row_start(&self, bin: DepthBin) -> usize {
        assert!(
            (bin.get() as usize) < self.bin_count(),
            "depth bin {} is outside a ladder of {} bins",
            bin.get(),
            self.bin_count()
        );
        self.row_starts[bin.get() as usize]
    }

    /// The depths this bin holds.
    ///
    /// A range rather than one endpoint, and `RangeInclusive` rather than a
    /// `(u32, u32)` pair, because the sole consumer is a row width where an off-by-one
    /// silently mis-sizes the table — inclusivity belongs in the type.
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

    /// How many cells a histogram binned this way holds — the flat vector's length.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        // PANIC-FREE: `new` pushes the leading 0 before anything else, so `row_starts`
        // is never empty.
        *self
            .row_starts
            .last()
            .expect("row_starts always carries its leading zero")
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
    /// this file pins `DEPTH_BIN_COUNT - 1` inside a `u16`, and that guarantee lives
    /// here rather than at the call sites.
    pub fn bins(&self) -> impl Iterator<Item = DepthBin> + '_ {
        // PANIC-FREE: the const assertion above pins the ladder's last index inside a
        // `u16`, and `new` is the only constructor.
        (0..self.bin_count()).map(|bin| DepthBin(bin as u16))
    }

    /// The deepest depth the ladder has a bin for. A site above it is subsampled down
    /// to it before it is entered.
    ///
    /// **The one number a consumer should read for that cap**, which is why the
    /// constant behind it is private: the subsampling target and the ladder's own top
    /// cannot then disagree.
    #[must_use]
    pub fn max_depth(&self) -> u32 {
        // PANIC-FREE: `ladder_tops` asserts the ladder is non-empty and ends at its cap
        // before returning, and `new` is its only caller.
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
    ///
    /// **This test and [`the_ladder_holds_583_cells`] are the whole oracle against a
    /// moved constant.** Every other test in this file is written in terms of
    /// `EXACT_DEPTH_LIMIT`, `MAX_BINNED_DEPTH` and `DEPTH_BIN_COUNT`, so its
    /// expectations move with an edit to them: those check that whatever ladder is
    /// configured is *internally consistent*, which is a different question from
    /// whether it is the ladder that was measured. Only the literals here and in the
    /// cell count carry the measurement.
    #[test]
    fn the_widening_bins_top_out_at_the_measured_depths() {
        let edges = DepthBinEdges::new();
        let widening_tops: Vec<u32> = (FIRST_WIDENING_BIN..edges.bin_count())
            .map(|bin| *edges.depth_range(DepthBin(bin as u16)).end())
            .collect();

        assert_eq!(
            widening_tops,
            vec![10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124],
            "the adopted ladder (research note §4.3)"
        );
        assert_eq!(edges.bin_count(), 20);
        assert_eq!(edges.max_depth(), 124);
    }

    // ---- the census ladder: the same rungs, and ten more ---------------------

    /// **The property that makes the two ladders one**: every depth at which the
    /// twenty-bin ladder changes bin is a depth at which the thirty-bin one changes bin
    /// too, so a census code turns into a histogram bin by collapsing everything above
    /// 124 and the two routes cannot bin differently.
    ///
    /// Asserted on the edges themselves rather than on `bin_for`, because it is the
    /// edges that a stored code indexes.
    #[test]
    fn every_edge_of_the_twenty_bin_ladder_is_an_edge_of_the_census_one() {
        let histogram = DepthBinEdges::new();
        let census = DepthBinEdges::for_census();

        let tops_of = |edges: &DepthBinEdges| -> Vec<u32> {
            edges
                .bins()
                .map(|bin| *edges.depth_range(bin).end())
                .collect()
        };
        let short = tops_of(&histogram);
        let long = tops_of(&census);

        assert_eq!(
            long[..short.len()],
            short[..],
            "the census ladder's first twenty rungs are the histogram ladder's own"
        );
        assert_eq!(census.bin_count(), CENSUS_DEPTH_BIN_COUNT);
        assert_eq!(
            long[short.len()..],
            [159, 204, 262, 336, 431, 553, 709, 910, 1168, 1498],
            "the ten rungs the extension adds, at the same ratio as the eleven below them"
        );
        assert_eq!(
            histogram.cell_count(),
            583,
            "extending the census ladder must not move the histogram route's cell table"
        );
    }

    /// Below the twenty-bin ladder's cap the two agree on every depth, so a run that
    /// records on one and reads on the other is only wrong above 124.
    #[test]
    fn the_two_ladders_answer_alike_at_every_depth_up_to_the_shorter_ones_cap() {
        let histogram = DepthBinEdges::new();
        let census = DepthBinEdges::for_census();
        for depth in 0..=MAX_BINNED_DEPTH {
            assert_eq!(
                histogram.bin_for(depth),
                census.bin_for(depth),
                "the two ladders disagree at depth {depth}"
            );
        }
    }

    /// **The half of this the extension exists for**, and the half a ladder that
    /// saturated early would still pass: a position a sample carries two copies of holds
    /// about twice that sample's median depth, so twice a depth has to land at least one
    /// rung above it or the two are written identically and the signal is gone. At the
    /// twenty-bin ladder that failed from 98 reads a position upwards.
    #[test]
    fn a_doubled_depth_lands_at_least_one_rung_higher_all_the_way_to_the_ceiling() {
        let census = DepthBinEdges::for_census();
        let ceiling = census.max_depth();
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

    /// **583 cells**, because a bin's row must be as wide as its deepest site's
    /// alternative count: 45 in the nine exact bins and 538 in the eleven widening
    /// ones. It is the number `spec/parameter_prepass_generic.md` §9 prices the
    /// accumulator's memory against, so it is stated here rather than left to be
    /// counted. The second half of the oracle — see the test above.
    #[test]
    fn the_ladder_holds_583_cells() {
        let edges = DepthBinEdges::new();

        let exact_cells: usize = (0..=EXACT_DEPTH_LIMIT).map(|d| d as usize + 1).sum();
        assert_eq!(exact_cells, 45);
        assert_eq!(edges.cell_count(), 583);
        assert_eq!(edges.cell_count() - exact_cells, 538);
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
        assert_eq!(edges.max_depth(), 124, "the measured cap, as a literal");
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

    /// `bins()` yields every bin once, in order, and nothing else — it exists so that a
    /// caller walking the ladder does not spell the `u16` cast by hand, and a walk that
    /// skipped or repeated a bin would leave a histogram row unreported with nothing to
    /// say so.
    #[test]
    fn bins_yields_every_bin_once_in_order() {
        let edges = DepthBinEdges::new();
        let walked: Vec<DepthBin> = edges.bins().collect();

        assert_eq!(walked.len(), edges.bin_count());
        assert!(walked.windows(2).all(|pair| pair[0] < pair[1]));
        assert_eq!(walked.first().copied(), Some(DepthBin(0)));
        assert_eq!(
            walked.last().copied(),
            Some(DepthBin(edges.bin_count() as u16 - 1))
        );
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

    /// **A bin index one past the end must fail, not answer.** `row_starts` carries a
    /// trailing total, so before the range check this returned `cell_count()` — a
    /// number shaped exactly like a row offset — where `depth_range` panicked on the
    /// identical input. Milestone B indexes a flat cell vector at
    /// `row_start(bin) + alt_reads`, so an off-by-one bin would have addressed another
    /// bin's cells rather than failing.
    #[test]
    #[should_panic(expected = "outside a ladder of 20 bins")]
    fn row_start_rejects_a_bin_past_the_end_of_the_ladder() {
        let edges = DepthBinEdges::new();
        let _ = edges.row_start(DepthBin(edges.bin_count() as u16));
    }

    /// The sibling accessor rejects the same input, so the two agree about what an
    /// out-of-range bin means.
    #[test]
    #[should_panic(expected = "outside a ladder of 20 bins")]
    fn depth_range_rejects_a_bin_past_the_end_of_the_ladder() {
        let edges = DepthBinEdges::new();
        let _ = edges.depth_range(DepthBin(edges.bin_count() as u16));
    }

    /// **A cap too low for the bin count is rejected, not quietly flattened.** This is
    /// the branch the adopted constants never exercise, which is the whole reason
    /// `ladder_tops` takes its shape as parameters. At `exact = 8, cap = 12` there is
    /// room for four more bins and the ladder is asked for eleven, so it walks to 19 —
    /// and depths 13 to 19 would sit above a cap that says 12.
    #[test]
    #[should_panic(expected = "tops out at 19 instead of its cap")]
    fn a_cap_too_low_for_the_bin_count_is_rejected() {
        let _ = ladder_tops(8, 12, 20);
    }

    /// The same guard from the other direction: **too many bins for the cap**, which is
    /// the plausible mis-edit — raising the bin count for resolution without noticing
    /// that the widening region has only 116 depths to spread over.
    #[test]
    #[should_panic(expected = "more bins here than the cap has room for")]
    fn more_bins_than_the_cap_has_room_for_is_rejected() {
        let _ = ladder_tops(8, 124, 200);
    }

    /// **Strict increase is by construction, not by assertion** — clamping to the cap
    /// before forcing each top above the last is what buys it, and applying the two
    /// the other way round is what produced repeated tops and inverted ranges. Asserted
    /// over the adopted shape and over one where the geometric step rounds below a
    /// whole depth a bin, so the `+1` clamp is doing all the work.
    #[test]
    fn every_ladder_shape_is_strictly_increasing() {
        for (exact, cap, bins) in [
            (EXACT_DEPTH_LIMIT, MAX_BINNED_DEPTH, DEPTH_BIN_COUNT),
            (8, 124, 60),
            (1, 124, 20),
        ] {
            let tops = ladder_tops(exact, cap, bins);
            assert!(
                tops.windows(2).all(|pair| pair[0] < pair[1]),
                "exact={exact} cap={cap} bins={bins} gave {tops:?}"
            );
            assert_eq!(tops.last().copied(), Some(cap));
            assert_eq!(tops.len(), bins);
        }
    }

    /// `Default` is the constructor a caller reaches through `#[derive(Default)]` on a
    /// containing struct, and nothing otherwise pins it to `new()`. Compared through
    /// the observable ladder rather than by `==`, because `DepthBinEdges` deliberately
    /// does not derive `PartialEq`.
    #[test]
    fn default_edges_are_the_adopted_ladder() {
        let (default, new) = (DepthBinEdges::default(), DepthBinEdges::new());

        assert_eq!(default.bin_count(), new.bin_count());
        assert_eq!(default.cell_count(), new.cell_count());
        assert_eq!(default.max_depth(), new.max_depth());
        for bin in 0..new.bin_count() {
            let bin = DepthBin(bin as u16);
            assert_eq!(default.depth_range(bin), new.depth_range(bin));
            assert_eq!(default.row_start(bin), new.row_start(bin));
        }
    }

    proptest::proptest! {
        /// `bin_for`'s doc claims totality over every `u32`, and the index it returns
        /// addresses a table — so "always in range" is the invariant that matters.
        /// Three sampled points could satisfy a change to the search that still left a
        /// gap; this cannot.
        #[test]
        fn bin_for_answers_an_existing_bin_for_every_depth(depth in proptest::num::u32::ANY) {
            let edges = DepthBinEdges::new();
            let bin = edges.bin_for(depth);

            proptest::prop_assert!((bin.get() as usize) < edges.bin_count());
            proptest::prop_assert!(edges.row_start(bin) < edges.cell_count());
        }
    }
}
