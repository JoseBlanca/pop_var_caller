//! A stratum's two slippage *shares* as an older parameters file carries them: the curve each was
//! read off, and where an emitted share came from.
//!
//! **What the two shares are.** When a read crosses a repeat tract, the copying steps before
//! sequencing sometimes add or drop whole repeat units, so the read reports a tract of the wrong
//! length. Of the reads that slipped, the **direction split** is the share that came back
//! *shorter* rather than longer, and the **fall-off** is how much rarer a two-unit slip is than a
//! one-unit slip. The third number, the slippage *level*, is how often a read slips at all; its
//! curve is in [`slippage_curve`](super::slippage_curve).
//!
//! **Nothing in a run draws these curves now.** The climb, the repeat-tract fit before the
//! penalised one, drew each share across repeat count per motif period — one of three shapes on
//! the logit scale, `log(p / (1 − p))`, chosen by leaving each stratum out in turn — and blended
//! each stratum's own share with it (`doc/devel/ng/spec/str_slippage_level_curve.md` §5.1). The
//! climb was deleted on 2026-10-16, and the code that fitted the curves later the same day. A
//! parameters file the climb wrote still records the curve each share came from; this module holds
//! what reading such a file needs: [`ShareCurve`] with its [`ShareShape`] and
//! [`ShareCurveSource`], and [`ShareSource`]. **Nothing in a run evaluates the curve**: reading a
//! file copies each share, its curve and whether the stratum sat inside the curve's range exactly
//! as written, and calling uses the share.

use super::slippage_curve::{CurveReach, PenalisedPoint};

// ---------------------------------------------------------------------
// The shape a share follows
// ---------------------------------------------------------------------

/// Which shape one share follows across repeat count, on the logit scale.
///
/// **The climb fitted all three and chose one by measurement**, per motif period and per share,
/// because no single shape is right everywhere. Over the ten cases both cohorts can speak to —
/// four motif periods on HG002, one on tomato, two shares each — [`ShareShape::Flat`] wins six,
/// [`ShareShape::Sloping`] two and [`ShareShape::Turning`] two
/// (`doc/devel/ng/reports/str_slippage_share_families_2026-08-20.md`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ShareShape {
    /// One value for the whole motif period — the share does not move with repeat count.
    ///
    /// **This is the honest answer wherever there is no trend to fit**, and it is a real winner
    /// rather than a fallback: HG002's tetranucleotide direction split is predicted better by the
    /// weighted mean than by either bending shape.
    Flat,
    /// The share moves one way across repeat count and never turns back.
    Sloping,
    /// The share turns once — falls and then rises, or rises and then falls.
    ///
    /// **Its worst failure is a period with barely more strata than it has coefficients**, where
    /// every leave-one-out fit passes exactly through the strata that remain and the score has no
    /// residual left to notice a bad shape with. It is allowed at four strata anyway, for two
    /// reasons. Outside the fitted repeat range the curve is held flat, so it cannot run away at
    /// the ends. And at the one real period where it happens — HG002's trinucleotide direction
    /// split, whose four strata swing 0.89, 0.33, 0.36, 0.89 across 6 to 9 repeats — it predicts
    /// a stratum it did not see to 0.21 logit units, better than eight of the other nine cases
    /// manage with more strata behind them.
    Turning,
}

// ---------------------------------------------------------------------
// The curve
// ---------------------------------------------------------------------

/// How one share moves across repeat count, for one slippage group at one motif period.
///
/// **Held flat outside the repeat counts it was fitted over, never continued** — the same rule
/// the level's curve follows, and for the same reason: a shape fitted over 6 to 9 repeats says
/// nothing about 30.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShareCurve {
    /// Which shape won, or [`ShareShape::Flat`] where no shape could be chosen.
    pub shape: ShareShape,
    /// The fitted curve on the logit scale, read at `repeats`:
    /// `intercept + slope · (repeats − centre) + bend · (repeats − centre)²`.
    ///
    /// `slope` is zero for a flat curve and `bend` for anything but a turning one.
    pub intercept: f64,
    /// See [`ShareCurve::intercept`].
    pub slope: f64,
    /// See [`ShareCurve::intercept`].
    pub bend: f64,
    /// The repeat count the other two are measured from — the weighted mean of the strata that
    /// fed the curve. **Centring is arithmetic, not design**: a quadratic in raw repeat count
    /// puts a term of 900 beside a term of 1 and loses precision solving for them.
    pub centre: f64,
    /// The lowest repeat count of a stratum that fed this curve.
    ///
    /// **A curve that does not depend on repeat count spans the whole axis.** The two bottom
    /// rungs of the ladder — another period's mean, and the built-in default — are one value
    /// everywhere, so there is no range to hold them at the ends of and
    /// [`ShareCurve::reach`] reports every repeat count as inside.
    pub fitted_from: u64,
    /// The highest. Beyond `fitted_from ..= fitted_to` the share is held at the nearer end.
    pub fitted_to: u64,
    /// **How far the curve landed from a stratum it had not seen**, in logit units — the median
    /// over leaving each contributing stratum out in turn, refitting, and predicting it.
    ///
    /// **Logit units, because that is what it was weighed against**: a stratum's own error — a
    /// share `p` measured on `S` slipped reads has a logit standard error of
    /// `1 / sqrt(S · p · (1 − p))` — is in the same units, so the two combined without a
    /// conversion. For scale, the winning shape scores between 0.167 and 0.842 across the ten
    /// cases both cohorts can speak to, where a well-measured stratum's own error is 0.033.
    pub held_out_error: f64,
    /// How many strata stood behind it. A curve through four strata and one through twenty-three
    /// are both curves, and a consumer must be able to tell them apart.
    pub strata: usize,
    /// Which rung of the fallback ladder produced this curve.
    pub source: ShareCurveSource,
}

/// Which rung of the fallback ladder a share curve came from.
///
/// **The ladder exists because these numbers are a prior rather than a measurement.** A stratum
/// with no share at all cannot be scored by the read likelihood, so a coarse answer that says it
/// is coarse beats no answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShareCurveSource {
    /// Fitted through this motif period's own strata, with the shape chosen by leaving each of
    /// them out in turn. This is the intended case.
    ThisPeriod,
    /// This period had strata but too few to choose a shape between, so it is their weighted mean
    /// and nothing more. [`ShareCurve::strata`] says how few.
    ThisPeriodUnscored,
    /// This motif period had no fitted stratum at all, so the curve is the weighted mean over
    /// every other period of the same run.
    ///
    /// **Crossing motif periods is known to be worse than staying inside one** — the direction
    /// split runs 1.4-fold across tomato's homopolymers and 4.9-fold across its dinucleotides —
    /// which is why this rung is recorded rather than blended into the one above.
    OtherPeriods,
    /// The run fitted no stratum anywhere, so the curve is a built-in constant: the climb used a
    /// direction split of 0.6 and a fall-off of 0.5. **The crudest number in the parameter set**,
    /// and the only one that is not this run's own data in any form.
    BuiltInDefault,
}

impl ShareCurve {
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

// ---------------------------------------------------------------------
// Where an emitted share came from
// ---------------------------------------------------------------------

/// Where one emitted share came from.
///
/// **The three are one formula with the curve's weight at its two ends and in between**, so a
/// consumer that only wants the number need not match on which end it is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShareSource {
    /// The stratum's own fit, whole. Its period had no curve at all.
    Stratum,
    /// The curve, whole. The stratum had no fit of its own.
    Curve,
    /// Both, weighted by how precisely each determines the share.
    Blend { curve_weight: f64 },
    /// **A point on its motif length's penalised curve, fitted with this stratum's own reads** —
    /// neither the stratum's fit alone nor a borrowed curve, as for
    /// [`LevelSource::PenalisedCurve`](super::slippage_curve::LevelSource::PenalisedCurve).
    PenalisedCurve(PenalisedPoint),
}

impl ShareSource {
    /// The share the parametric curve carried — 1 for [`ShareSource::Curve`], 0 for
    /// [`ShareSource::Stratum`], the weight for a blend — and `None` for
    /// [`ShareSource::PenalisedCurve`], which no parametric curve was weighed against.
    pub fn curve_weight(self) -> Option<f64> {
        match self {
            Self::Stratum => Some(0.0),
            Self::Curve => Some(1.0),
            Self::Blend { curve_weight } => Some(curve_weight),
            Self::PenalisedCurve(_) => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three parametric variants report the weight the curve carried; a point on a penalised
    /// curve, never weighed against one, reports none.
    #[test]
    fn every_share_source_reports_the_curve_weight_or_none_for_a_penalised_point() {
        assert_eq!(ShareSource::Stratum.curve_weight(), Some(0.0));
        assert_eq!(ShareSource::Curve.curve_weight(), Some(1.0));
        assert_eq!(
            ShareSource::Blend { curve_weight: 0.6 }.curve_weight(),
            Some(0.6)
        );
        let point = PenalisedPoint {
            strength: 3.0,
            eased_here: false,
        };
        assert_eq!(ShareSource::PenalisedCurve(point).curve_weight(), None);
    }

    /// A turning curve fitted over 8 to 12 repeats.
    fn a_turning_curve() -> ShareCurve {
        ShareCurve {
            shape: ShareShape::Turning,
            intercept: 0.5,
            slope: 0.0,
            bend: 0.6,
            centre: 10.0,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.21,
            strata: 5,
            source: ShareCurveSource::ThisPeriod,
        }
    }

    /// A stratum outside the repeat counts its curve was fitted over is told apart from one inside.
    #[test]
    fn a_stratum_knows_whether_it_sat_inside_the_curves_range() {
        let curve = a_turning_curve();
        assert_eq!(curve.reach(3), CurveReach::BelowFitted);
        assert_eq!(curve.reach(8), CurveReach::Inside);
        assert_eq!(curve.reach(12), CurveReach::Inside);
        assert_eq!(curve.reach(60), CurveReach::AboveFitted);
    }
}
