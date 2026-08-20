//! How the slippage level rises with repeat count, and where an emitted level came from.
//!
//! **The problem this exists for.** The four slippage numbers are fitted separately in every
//! (motif period, repeat count) cell, so neighbouring cells are free to disagree, and a cell too
//! thin to fit takes a neighbour's value whole — which costs 15 to 25% of the level per repeat
//! count borrowed across. On the tomato bench run 65 of 71 cells got no answer from their own
//! tracts at all.
//!
//! **What this module holds.** The *level* — how often a read reports a tract length other than
//! its allele's — becomes a curve in repeat count, fitted once per motif period and evaluated at
//! every cell. The other three numbers are untouched and keep their per-cell fit
//! (`doc/devel/ng/spec/str_slippage_level_curve.md` §1.2).
//!
//! **The family, and why it has a shape number instead of a fixed shape.** Two cohorts fitted
//! with nothing linking their cells prefer *opposite* shapes over the same repeat counts: over
//! 8 to 12 repeats, tomato's homopolymers are predicted best by an exponential (12.4% held-out
//! against a straight line's 33.6%) and HG002's by a straight line (8.0% against 31.2%). And both
//! fixed shapes produce impossible numbers a few repeat counts outside the range they saw — an
//! exponential fitted on HG002's 8-to-12 cells says the level at 30 repeats is 205, where the cell
//! at 30 repeats fits 0.120. So the shape is fitted rather than chosen
//! (`doc/devel/ng/reports/str_slippage_shape_2026-08-20.md` §4).
//!
//! Design: `doc/devel/ng/spec/str_slippage_level_curve.md`. Build order:
//! `doc/devel/ng/impl_plan/str_slippage_level_curve.md`.

use std::fmt;

// ---------------------------------------------------------------------
// The shape number
// ---------------------------------------------------------------------

/// How the slippage level compounds as repeat count rises.
///
/// **At 0 each extra repeat multiplies the level by a fixed factor; at 1 each extra repeat adds
/// a fixed amount.** In between it is the power the level is raised to before a straight line is
/// fitted through it:
///
/// ```text
/// level ^ rise_shape  =  intercept + slope · repeat_count      for rise_shape > 0
/// log(level)          =  intercept + slope · repeat_count      for rise_shape = 0
/// ```
///
/// **1 is production's own shape** — `baseline + slope · units`, clamped, in
/// `src/ssr/cohort/em.rs`'s `candidate_level` — so the family contains what production does
/// rather than replacing it with something unrelated.
///
/// Fitted per motif period and shared by every slippage group at that period, because a
/// curvature needs the whole span of repeat counts visible at once where a level and a slope do
/// not (spec §3).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct RiseShape(f64);

impl RiseShape {
    /// The multiplying end: each extra repeat multiplies the level by a fixed factor.
    pub const MULTIPLYING: Self = Self(0.0);

    /// The adding end: each extra repeat adds a fixed amount, which is production's shape.
    pub const ADDING: Self = Self(1.0);

    /// `None` outside `[0, 1]` or on a value that is not a number.
    ///
    /// **The range is closed at both ends and not widened**: past 1 the curve is convex in
    /// repeat count, which would say slippage accelerates without limit, and no cell of either
    /// cohort asks for it — HG002's homopolymers land exactly at 1 and its dinucleotides at 0.8.
    pub fn new(value: f64) -> Option<Self> {
        (value.is_finite() && (0.0..=1.0).contains(&value)).then_some(Self(value))
    }

    pub fn get(self) -> f64 {
        self.0
    }

    /// Whether this is the multiplying end, where the straight line is fitted to the logarithm.
    ///
    /// **The two branches are genuinely different arithmetic**, not a limit taken numerically:
    /// the power transform `(level^s − 1) / s` is undefined at `s = 0`, and its limit there is
    /// `log(level)`.
    pub fn is_multiplying(self) -> bool {
        self.0 == 0.0
    }
}

impl fmt::Display for RiseShape {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.2}", self.0)
    }
}

// ---------------------------------------------------------------------
// What the curve is fitted from, and what it is
// ---------------------------------------------------------------------

/// One cell's own answer, as the curve fit sees it.
///
/// **Only cells that were fitted from their own tracts are here.** A cell that borrowed, or that
/// was refused, has nothing to say about the shape and would pull hardest if it were let in
/// (spec §4.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedCell {
    /// The reference tract's repeat count — the axis the curve runs along.
    pub repeats: u64,
    /// The cell's own fitted slippage level, in `(0, 1)`.
    pub level: f64,
    /// **How many of the cell's reads the fitted level says slipped** — `level × reads crossing`.
    /// This is what sets how precisely the cell determines its own level, and it is the weight
    /// the curve fit uses (spec §4.2).
    ///
    /// **Not the count of reads sitting off the reference tract length**, which is much larger
    /// and at a polymorphic tract is mostly genuine allele length: at HG002's 30-repeat
    /// homopolymer cell, 60 reads in 100 that cross the tract report a length other than the
    /// reference's and the fit attributes 12 of those 60 to slippage.
    pub slipped_reads: f64,
}

impl FittedCell {
    /// How precisely this cell determines its own level, as a share of the level itself.
    ///
    /// `1 / sqrt(slipped reads)` — the best case, since the level is a proportion over the reads
    /// that moved. The real estimator sums over the genotype and is noisier than that, so this
    /// **understates** the cell's error and therefore leans the blend of spec §7 toward the cell.
    ///
    /// **Fewer than one slipped read is read as one**, so the answer is a share of the level and
    /// never infinite. Such a cell determines nothing and does not feed a curve (spec §4.1); the
    /// clamp exists so that weighing it against a curve is arithmetic rather than a special case.
    pub fn relative_standard_error(&self) -> f64 {
        1.0 / self.slipped_reads.max(1.0).sqrt()
    }
}

/// How the slippage level rises with repeat count, for one slippage group at one motif period.
///
/// **Held flat outside the repeat counts it was fitted over, never continued** (spec §6); the
/// evaluation itself arrives in step A2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageCurve {
    /// How the level compounds; shared by every slippage group at this motif period.
    pub rise_shape: RiseShape,
    /// The straight line fitted to `level ^ rise_shape` against repeat count.
    pub intercept: f64,
    /// The same line's slope. **Positive by construction**: a fit returning a level that falls
    /// with repeat count is refused rather than emitted (spec §9).
    pub slope: f64,
    /// The lowest repeat count of a cell that fed this curve.
    pub fitted_from: u64,
    /// The highest. Beyond `fitted_from ..= fitted_to` the level is held at the nearer end.
    pub fitted_to: u64,
    /// **How far the curve landed from a cell it had not seen** — the median relative error over
    /// leaving each contributing cell out in turn. Spec §7 reads it as the curve's own precision
    /// when weighing the curve against a cell's own answer.
    pub held_out_error: f64,
    /// How many cells stood behind it. A curve through four cells and one through twenty-three
    /// are both curves, and a consumer must be able to tell them apart.
    pub cells: usize,
}

// ---------------------------------------------------------------------
// Where an emitted level came from
// ---------------------------------------------------------------------

/// Where a cell's emitted slippage level came from.
///
/// **After this change a level fitted from 8,000 slipped reads and one interpolated across a gap
/// look the same in the number alone**, and the mechanism that used to mark the second — the
/// borrowing rule — no longer sets the level. This is its replacement (spec §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LevelSource {
    /// The cell's own fit, whole. Its period had no curve.
    Cell,
    /// The curve, whole. The cell had no fit of its own.
    Curve,
    /// Both, weighted by how precisely each determines the level.
    ///
    /// `curve_weight` is the share the curve carried, in `[0, 1]`. It runs about 0.93 at a cell
    /// with 40 slipped reads behind it and about 0.06 at one with 8,000, at a curve whose
    /// held-out error is 4.4% (spec §7).
    Blend { curve_weight: f64 },
}

impl LevelSource {
    /// The share the curve carried, whichever variant this is — 1 for [`LevelSource::Curve`] and
    /// 0 for [`LevelSource::Cell`].
    ///
    /// **The three variants are one formula with the weight at its two ends and in between**, so
    /// a consumer that only wants the number should not have to match on which end it is.
    pub fn curve_weight(self) -> f64 {
        match self {
            Self::Cell => 0.0,
            Self::Curve => 1.0,
            Self::Blend { curve_weight } => curve_weight,
        }
    }
}

/// Whether a cell sat inside the repeat counts its curve was fitted over.
///
/// **Held flat is wrong in a known direction and that is why it is recorded**: the level
/// genuinely keeps rising, so a cell above the fitted range is under-stated (spec §6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveReach {
    /// Inside `fitted_from ..= fitted_to`.
    Inside,
    /// Below the lowest cell that fed the curve; the level is that cell's.
    BelowFitted,
    /// Above the highest; the level is that cell's.
    AboveFitted,
}

// ---------------------------------------------------------------------
// The knobs, and where each number comes from
// ---------------------------------------------------------------------

/// How many rungs the shape number is searched over, spanning `[0, 1]` inclusive.
///
/// **21, which is steps of 0.05.** The two cohorts' periods land at 0.00, 0.80 and 1.00, so the
/// grid has to resolve 0.05 to distinguish period 2's answer from the adding end; nothing
/// measured asks for finer, and each extra rung costs one weighted least-squares fit per cell per
/// slippage group.
pub const RISE_SHAPE_RUNGS: usize = 21;

/// How few contributing cells leave a period without a curve.
///
/// **4, and it is arithmetic rather than a measurement — the smallest count at which leaving one
/// cell out still leaves a line and a spare.** It is soft, and spec §11 records the measurement
/// that would settle it along with the reason to expect it is too low: HG002's period 3 has
/// exactly four cells and its best rung predicts a held-out cell only to 31%, against 3.8% at
/// period 2's twenty.
pub const MIN_CELLS_FOR_A_CURVE: usize = 4;

/// How far a cell may sit from the curve, in the two errors combined, before the curve is taken
/// to be wrong about that cell rather than the cell unlucky.
///
/// **2.5 combined errors.** The case it exists for is measured: at HG002's 9-repeat homopolymer
/// the curve is 55% high against a cell whose own sampling error is 1.7%, a gap of **9.3**
/// combined errors, which no sampling noise produces. *The value is a conventional outlier knee
/// and is soft* — spec §7.2 records that without it the blend was already within 5.8% there
/// rather than 55%, because a cell with 3,520 slipped reads outweighs the curve on its own.
pub const DISAGREEMENT_KNEE: f64 = 2.5;

/// What a run may change about how the level's curve is fitted.
///
/// **The weight each cell carries is deliberately not here.** It is the cell's slipped reads, and
/// the measurement says the choice barely matters: the winning family's held-out error moves from
/// 5.13% unweighted to 4.39% weighted by slipped reads, and the ranking of families is identical
/// under all four weights tried. A knob nobody sets is how a measured decision gets changed by
/// accident.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageCurveConfig {
    /// Rungs of the shape-number grid; see [`RISE_SHAPE_RUNGS`].
    pub rise_shape_rungs: usize,
    /// Below this many contributing cells a period gets no curve; see [`MIN_CELLS_FOR_A_CURVE`].
    pub min_cells_for_a_curve: usize,
    /// Where the curve stops being trusted about a cell; see [`DISAGREEMENT_KNEE`].
    pub disagreement_knee: f64,
    /// **`false` keeps every cell's own fitted level and draws no curve at all**, which is the
    /// arm the parity oracle runs: nothing here touches how a cell is fitted, so a cell's own
    /// level moving with the curve off is a defect in the plumbing.
    pub draw_curves: bool,
}

impl SlippageCurveConfig {
    /// The shape numbers the fit will try, evenly spaced over `[0, 1]` inclusive.
    ///
    /// **The grid lives here rather than inside the fit** so that a test asserting which shapes
    /// are reachable is asserting the same list the fit walks. One rung degenerates to the
    /// multiplying end alone, which is a configuration no caller should reach but which must not
    /// divide by zero.
    pub fn rise_shape_grid(&self) -> Vec<RiseShape> {
        let rungs = self.rise_shape_rungs.max(1);
        if rungs == 1 {
            return vec![RiseShape::MULTIPLYING];
        }
        (0..rungs)
            .map(|rung| {
                RiseShape::new(rung as f64 / (rungs - 1) as f64)
                    .expect("an evenly spaced rung of [0, 1] is a rise shape")
            })
            .collect()
    }
}

impl Default for SlippageCurveConfig {
    fn default() -> Self {
        Self {
            rise_shape_rungs: RISE_SHAPE_RUNGS,
            min_cells_for_a_curve: MIN_CELLS_FOR_A_CURVE,
            disagreement_knee: DISAGREEMENT_KNEE,
            draw_curves: true,
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_rise_shape_outside_zero_to_one_is_refused() {
        assert!(RiseShape::new(0.0).is_some());
        assert!(RiseShape::new(1.0).is_some());
        assert!(RiseShape::new(0.8).is_some());
        assert!(RiseShape::new(-0.01).is_none());
        assert!(RiseShape::new(1.01).is_none());
        assert!(RiseShape::new(f64::NAN).is_none());
        assert!(RiseShape::new(f64::INFINITY).is_none());
    }

    /// The two ends are the two arithmetics, and only exactly zero takes the logarithm branch.
    #[test]
    fn only_the_multiplying_end_takes_the_logarithm_branch() {
        assert!(RiseShape::MULTIPLYING.is_multiplying());
        assert!(!RiseShape::ADDING.is_multiplying());
        assert!(
            !RiseShape::new(0.05)
                .expect("inside the range")
                .is_multiplying()
        );
        assert_eq!(RiseShape::MULTIPLYING.get(), 0.0);
        assert_eq!(RiseShape::ADDING.get(), 1.0);
    }

    /// A cell's own precision is set by the reads that moved, not by the reads that crossed.
    #[test]
    fn a_cells_own_error_falls_as_the_square_root_of_its_slipped_reads() {
        let thin = FittedCell {
            repeats: 9,
            level: 0.0067,
            slipped_reads: 100.0,
        };
        let fat = FittedCell {
            repeats: 9,
            level: 0.0067,
            slipped_reads: 10_000.0,
        };
        assert!((thin.relative_standard_error() - 0.1).abs() < 1e-12);
        assert!((fat.relative_standard_error() - 0.01).abs() < 1e-12);
    }

    /// A cell the fit says nothing slipped in must not divide by zero.
    #[test]
    fn a_cell_with_no_slipped_reads_gets_an_error_of_one_rather_than_infinity() {
        let empty = FittedCell {
            repeats: 8,
            level: 0.0,
            slipped_reads: 0.0,
        };
        assert!(empty.relative_standard_error().is_finite());
        assert!((empty.relative_standard_error() - 1.0).abs() < 1e-12);
    }

    /// The three variants are one weight, so a consumer never matches on which end it is.
    #[test]
    fn every_level_source_reports_the_share_the_curve_carried() {
        assert_eq!(LevelSource::Cell.curve_weight(), 0.0);
        assert_eq!(LevelSource::Curve.curve_weight(), 1.0);
        assert_eq!(
            LevelSource::Blend { curve_weight: 0.93 }.curve_weight(),
            0.93
        );
    }

    #[test]
    fn the_default_config_is_the_measured_constants_with_curves_on() {
        let config = SlippageCurveConfig::default();
        assert_eq!(config.rise_shape_rungs, RISE_SHAPE_RUNGS);
        assert_eq!(config.min_cells_for_a_curve, MIN_CELLS_FOR_A_CURVE);
        assert_eq!(config.disagreement_knee, DISAGREEMENT_KNEE);
        assert!(config.draw_curves);
    }

    /// The grid must be able to name the answers both cohorts gave, or the fit cannot return
    /// them: 21 rungs over `[0, 1]` are steps of 0.05, and 0.80 is one of them.
    #[test]
    fn the_shape_grid_resolves_the_answers_both_cohorts_gave() {
        let grid = SlippageCurveConfig::default().rise_shape_grid();
        assert_eq!(grid.len(), RISE_SHAPE_RUNGS);
        for wanted in [0.0, 0.8, 1.0] {
            assert!(
                grid.iter().any(|shape| (shape.get() - wanted).abs() < 1e-9),
                "the grid cannot name {wanted}"
            );
        }
        assert_eq!(grid.first().copied(), Some(RiseShape::MULTIPLYING));
        assert_eq!(grid.last().copied(), Some(RiseShape::ADDING));
    }

    /// A grid of one rung must not divide by zero when spacing the rungs.
    #[test]
    fn a_single_rung_grid_is_the_multiplying_end_alone() {
        let config = SlippageCurveConfig {
            rise_shape_rungs: 1,
            ..SlippageCurveConfig::default()
        };
        assert_eq!(config.rise_shape_grid(), vec![RiseShape::MULTIPLYING]);
        let none = SlippageCurveConfig {
            rise_shape_rungs: 0,
            ..SlippageCurveConfig::default()
        };
        assert_eq!(none.rise_shape_grid(), vec![RiseShape::MULTIPLYING]);
    }
}
