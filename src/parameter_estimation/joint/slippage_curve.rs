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
//! - [`SlippageCurve`] and its shape number [`RiseShape`]: the curve as written, carried beside
//!   the level so a reader can see which curve stood behind it. **Nothing in a run evaluates the
//!   curve**: reading a file copies the level, the curve and whether the cell sat inside the
//!   curve's range exactly as written, and calling uses the level;
//! - [`LevelSource`], [`PenalisedPoint`] and [`CurveReach`]: where an emitted level came from —
//!   the cell, the curve, a blend of the two, or a point on the penalised fit's curve — and whether
//!   the cell sat inside the curve's range;
//! - [`LEVEL_CEILING`], the largest level the climb's curves reported.

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
/// **The climb held it flat outside the repeat counts it was fitted over, never continued** (spec
/// §6): a level for a cell outside that range was the nearer end's, and the file records which
/// side of the range the cell sat on ([`CurveReach`]).
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
}

/// **The largest level the climb's curves reported**, and so the largest level an older file can
/// carry from one.
///
/// **Not one.** At a level of one every read has misread, and the emission model divides by the
/// share of reads that did not. HG002's longest homopolymer cell fits 0.120, so a curve reaching
/// this had been asked for a repeat count far outside anything fitted.
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
}
