//! A slippage-level curve as an older parameters file carries it, and where an emitted level came
//! from.
//!
//! **Nothing in a run draws these curves now.** The repeat-tract fit before the penalised one —
//! the climb — fitted each (motif period, repeat count) cell on its own, then drew the slippage
//! *level*, how often a read reports a tract length other than its allele's, as a curve across
//! repeat count and blended each cell's own level with it
//! (`doc/devel/ng/spec/str_slippage_level_curve.md`, which every "spec §" below cites). The climb
//! was deleted on 2026-10-16, and the code that fitted and blended the curves later the same day.
//! A parameters file the climb wrote still records, for each level, the curve it came from or was
//! blended with; this module holds what reading and calling such a file needs:
//!
//! - [`SlippageCurve`] and its shape number [`RiseShape`]: the curve as written, read back at any
//!   repeat count, held flat outside the counts it was fitted over;
//! - [`LevelSource`], [`PenalisedPoint`] and [`CurveReach`]: where an emitted level came from —
//!   the cell, the curve, a blend of the two, or a point on the penalised fit's curve — and whether
//!   the cell sat inside the curve's range;
//! - [`LEVEL_FLOOR`] and [`LEVEL_CEILING`], the bounds every level read off a curve is held inside.

use std::fmt;

use crate::float;

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
// The curve
// ---------------------------------------------------------------------

/// How the slippage level rises with repeat count, for one slippage group at one motif period.
///
/// **Held flat outside the repeat counts it was fitted over, never continued** — see
/// [`SlippageCurve::level_at`] and spec §6.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageCurve {
    /// How the level compounds; shared by every slippage group at this motif period.
    pub rise_shape: RiseShape,
    /// The straight line fitted to `level ^ rise_shape` against repeat count.
    pub intercept: f64,
    /// The same line's slope. **Positive by construction**: the climb's curve fit refused a level
    /// that falls with repeat count rather than writing it (spec §9).
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

impl SlippageCurve {
    /// The slippage level this curve gives at `repeats`.
    ///
    /// **Inside the repeat counts the curve was fitted over it is the fitted line, read back
    /// through the shape number. Outside, it is held at the nearer fitted end.**
    ///
    /// Holding flat rather than continuing the line is the whole of spec §6, and the reason is
    /// that continuing it produces numbers that are not probabilities: an exponential fitted on
    /// HG002's 8-to-12 homopolymer cells says 205 at 30 repeats, where the cell at 30 repeats
    /// fits 0.120. **Held flat is wrong in a known direction** — slippage genuinely keeps rising,
    /// so a tract longer than anything fitted is under-stated — and [`SlippageCurve::reach`] is
    /// what lets a consumer see it.
    pub fn level_at(&self, repeats: u64) -> f64 {
        assert!(
            self.fitted_from <= self.fitted_to,
            "a curve fitted from {} to {} repeats has no range to hold at its ends",
            self.fitted_from,
            self.fitted_to
        );
        let inside = repeats.clamp(self.fitted_from, self.fitted_to);
        self.level_on_the_line(inside as f64)
    }

    /// Whether `repeats` sat inside the repeat counts this curve was fitted over.
    pub fn reach(&self, repeats: u64) -> CurveReach {
        if repeats < self.fitted_from {
            CurveReach::BelowFitted
        } else if repeats > self.fitted_to {
            CurveReach::AboveFitted
        } else {
            CurveReach::Inside
        }
    }

    /// The fitted line read back into a level, with no regard for the fitted range.
    ///
    /// **Clamped into `(0, 1)` at both ends**, because the line is fitted in a transformed space
    /// where nothing stops it crossing either boundary: at the multiplying end the exponential
    /// passes 1 a few repeat counts above its range, and at the adding end the line passes 0 a
    /// few below — HG002's homopolymer line goes negative below 7.4 repeats. The clamp is not
    /// what makes extrapolation safe; holding at the fitted ends is. It is here so that no
    /// caller can be handed a number that is not a probability.
    fn level_on_the_line(&self, repeats: f64) -> f64 {
        let line = self.intercept + self.slope * repeats;
        let level = if self.rise_shape.is_multiplying() {
            float::exp(line)
        } else {
            // `level ^ s = line`, so `level = line ^ (1/s)`; a line at or below zero has no
            // real root and means the level has fallen off the bottom of the family.
            if line <= 0.0 {
                0.0
            } else {
                float::powf(line, 1.0 / self.rise_shape.get())
            }
        };
        if level.is_finite() {
            level.clamp(LEVEL_FLOOR, LEVEL_CEILING)
        } else {
            LEVEL_CEILING
        }
    }
}

/// The smallest level a curve may report.
///
/// **Not zero.** A level of exactly zero says a read can never misread the tract, which the
/// emission model then treats as certainty; and the blend of spec §7 works on the logarithm,
/// which has no value there. One read in ten million is far below anything either cohort
/// measures — the thinnest cell fitted is 3.7 in 1,000 — so the floor cannot bind on real data
/// and exists to keep the arithmetic total.
pub const LEVEL_FLOOR: f64 = 1e-7;

/// The largest level a curve may report.
///
/// **Not one.** At a level of one every read has misread, and the emission model divides by the
/// share of reads that did not. HG002's longest homopolymer cell fits 0.120, so a curve reaching
/// this has been asked for a repeat count far outside anything fitted.
pub const LEVEL_CEILING: f64 = 0.999;

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
    /// `curve_weight` is the share the curve carried, in `[0, 1]`. It runs about 0.81 at a cell
    /// with 40 slipped reads behind it and about 0.02 at one with 8,000, at a curve whose
    /// held-out error is 7.7% — HG002's homopolymers at the corrected ±8 recording window
    /// (spec §7, whose own figures are the ±4 ones).
    Blend { curve_weight: f64 },
    /// **A point on its motif length's penalised curve, fitted with this cell's own reads.**
    ///
    /// The penalised fit fits every repeat count of a motif length together, pulling each
    /// count's level toward its neighbours'. The value is neither the cell's fit alone nor a curve
    /// borrowed over a gap — the fit only reports cells with reads — so it is a variant of its own
    /// rather than [`LevelSource::Cell`] (`str_fit_to_calling.md` §2.3).
    PenalisedCurve(PenalisedPoint),
}

impl LevelSource {
    /// The share the parametric curve carried — 1 for [`LevelSource::Curve`], 0 for
    /// [`LevelSource::Cell`], the weight for a blend — and `None` for
    /// [`LevelSource::PenalisedCurve`], which no parametric curve was weighed against.
    ///
    /// **The first three variants are one formula with the weight at its two ends and in
    /// between**, so a consumer that only wants that number should not have to match on which end
    /// it is. A penalised point is not on that scale: it is neither borrowed nor the cell's fit
    /// alone, so it has no weight to report, and a consumer deciding whether a number is the
    /// cell's own must say what it does with one.
    pub fn curve_weight(self) -> Option<f64> {
        match self {
            Self::Cell => Some(0.0),
            Self::Curve => Some(1.0),
            Self::Blend { curve_weight } => Some(curve_weight),
            Self::PenalisedCurve(_) => None,
        }
    }
}

/// **What the parameters file records about one point on a penalised curve.**
///
/// The penalised fit draws one curve per motif length through every repeat count with reads; how
/// stiff the curve is was chosen from the data, and at some repeat counts the fit let the curve
/// bend more sharply because the reads insisted on a kink. Both facts travel with each number so
/// that a reader of the file can tell a value held close to its neighbours from one the reads set
/// on their own.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PenalisedPoint {
    /// The curve's stiffness strength: larger holds each repeat count's value closer to its
    /// neighbours'.
    pub strength: f64,
    /// Whether a bend was eased at this repeat count: the reads insisted on a kink here.
    pub eased_here: bool,
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

    /// The three variants are one weight, so a consumer never matches on which end it is.
    #[test]
    fn every_level_source_reports_the_share_the_curve_carried() {
        assert_eq!(LevelSource::Cell.curve_weight(), Some(0.0));
        assert_eq!(LevelSource::Curve.curve_weight(), Some(1.0));
        assert_eq!(
            LevelSource::Blend { curve_weight: 0.93 }.curve_weight(),
            Some(0.93)
        );
    }

    /// A point on a penalised curve was never weighed against a parametric curve, so it reports
    /// no weight rather than one of the two ends.
    #[test]
    fn a_penalised_point_reports_no_curve_weight() {
        let point = PenalisedPoint {
            strength: 12.5,
            eased_here: true,
        };
        assert_eq!(LevelSource::PenalisedCurve(point).curve_weight(), None);
    }

    /// A curve built at the adding end is a straight line in the level itself, and reading it
    /// back must return the line.
    #[test]
    fn the_adding_end_reads_back_as_a_straight_line_in_the_level() {
        // level = 0.005 * repeats - 0.035, which is 0.005 at 8 repeats and 0.115 at 30.
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        assert!((curve.level_at(8) - 0.005).abs() < 1e-12);
        assert!((curve.level_at(20) - 0.065).abs() < 1e-12);
        assert!((curve.level_at(30) - 0.115).abs() < 1e-12);
    }

    /// A curve built at the multiplying end multiplies by a fixed factor each repeat.
    #[test]
    fn the_multiplying_end_reads_back_as_a_fixed_factor_a_repeat() {
        // log(level) = -6.5 + 0.4 * repeats, so each repeat multiplies by exp(0.4) = 1.4918.
        let curve = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        let at_eight = curve.level_at(8);
        let at_nine = curve.level_at(9);
        assert!((float::exp(-6.5 + 0.4 * 8.0) - at_eight).abs() < 1e-15);
        assert!((at_nine / at_eight - float::exp(0.4)).abs() < 1e-12);
    }

    /// Beyond the repeat counts the curve saw, the level is the nearer fitted end's — never the
    /// line continued. This is the property that stops an exponential reporting 205.
    #[test]
    fn outside_the_fitted_range_the_level_is_held_at_the_nearer_end() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        assert_eq!(curve.level_at(7), curve.level_at(8));
        assert_eq!(curve.level_at(1), curve.level_at(8));
        assert_eq!(curve.level_at(13), curve.level_at(12));
        assert_eq!(curve.level_at(150), curve.level_at(12));

        // Continued rather than held, this curve reaches exp(-6.5 + 0.4*30) = 244.7 at 30
        // repeats — the shape of the failure that made holding flat the rule — where held it
        // stays at its 12-repeat value of 0.1827.
        assert_eq!(curve.level_on_the_line(30.0), LEVEL_CEILING);
        assert!((curve.level_at(30) - float::exp(-6.5 + 0.4 * 12.0)).abs() < 1e-15);
        assert!((curve.level_at(30) - 0.1827).abs() < 1e-4);
    }

    /// The property every consumer of this curve relies on: a longer tract never comes back
    /// less slippery than a shorter one. It has to hold at both ends of the shape range and in
    /// between, because a file may carry any shape number the climb's grid could pick.
    #[test]
    fn a_positive_slope_never_lets_a_longer_tract_come_back_less_slippery() {
        // Every shape number a file may carry, in steps of 0.05 — the grid the climb searched.
        let shapes = (0..=20).map(|step| RiseShape::new(step as f64 / 20.0).expect("in [0, 1]"));
        for shape in shapes {
            let (low, high) = (0.004_f64, 0.12_f64);
            let (low_t, high_t) = if shape.is_multiplying() {
                (float::ln(low), float::ln(high))
            } else {
                (
                    float::powf(low, shape.get()),
                    float::powf(high, shape.get()),
                )
            };
            let slope = (high_t - low_t) / (30.0 - 8.0);
            let curve = SlippageCurve {
                rise_shape: shape,
                intercept: low_t - slope * 8.0,
                slope,
                fitted_from: 8,
                fitted_to: 30,
                held_out_error: 0.05,
                cells: 23,
            };
            let mut previous = 0.0;
            for repeats in 1..=60 {
                let level = curve.level_at(repeats);
                assert!(
                    level >= previous - 1e-15,
                    "shape {shape} fell from {previous} to {level} at {repeats} repeats"
                );
                assert!((LEVEL_FLOOR..=LEVEL_CEILING).contains(&level));
                previous = level;
            }
            assert!((curve.level_at(8) - low).abs() < 1e-12, "shape {shape}");
            assert!((curve.level_at(30) - high).abs() < 1e-12, "shape {shape}");
        }
    }

    #[test]
    #[should_panic(expected = "has no range to hold at its ends")]
    fn a_curve_whose_range_runs_backwards_says_so_rather_than_clamping() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: 0.0,
            slope: 0.005,
            fitted_from: 30,
            fitted_to: 8,
            held_out_error: 0.044,
            cells: 23,
        };
        let _ = curve.level_at(12);
    }

    #[test]
    fn a_cell_knows_whether_it_sat_inside_the_curves_range() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        assert_eq!(curve.reach(7), CurveReach::BelowFitted);
        assert_eq!(curve.reach(8), CurveReach::Inside);
        assert_eq!(curve.reach(30), CurveReach::Inside);
        assert_eq!(curve.reach(31), CurveReach::AboveFitted);
    }

    /// A line fitted at the adding end goes negative below its range; the level may not.
    #[test]
    fn a_level_is_never_reported_outside_the_open_unit_interval() {
        let adding = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        // -0.035 + 0.005*5 = -0.010, a level a line is happy to produce and a caller cannot use.
        assert_eq!(adding.level_on_the_line(5.0), LEVEL_FLOOR);
        assert!(adding.level_on_the_line(1000.0) <= LEVEL_CEILING);

        let multiplying = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        assert_eq!(multiplying.level_on_the_line(100.0), LEVEL_CEILING);
        assert!(multiplying.level_on_the_line(-1000.0) >= LEVEL_FLOOR);
    }

    /// A shape number strictly between the two ends inverts through a root, and the round trip
    /// has to land back on the value the line was built from.
    #[test]
    fn a_shape_number_between_the_ends_round_trips_through_its_root() {
        let shape = RiseShape::new(0.8).expect("inside the range");
        // Build the line so that it passes exactly through 0.03 at 12 repeats and 0.06 at 25.
        let (low, high) = (float::powf(0.03, 0.8), float::powf(0.06, 0.8));
        let slope = (high - low) / (25.0 - 12.0);
        let curve = SlippageCurve {
            rise_shape: shape,
            intercept: low - slope * 12.0,
            slope,
            fitted_from: 12,
            fitted_to: 25,
            held_out_error: 0.038,
            cells: 20,
        };
        assert!((curve.level_at(12) - 0.03).abs() < 1e-12);
        assert!((curve.level_at(25) - 0.06).abs() < 1e-12);
        // and it sits between the two ends: below a straight line, above an exponential.
        let midpoint = curve.level_at(18);
        assert!(midpoint > 0.03 && midpoint < 0.06);
    }
}
