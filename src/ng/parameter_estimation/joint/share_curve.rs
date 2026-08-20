//! How a stratum's two slippage *shares* move across repeat count.
//!
//! **What the two shares are.** When a read crosses a repeat tract, the copying steps before
//! sequencing sometimes add or drop whole repeat units, so the read reports a tract of the wrong
//! length. Of the reads that slipped, the **direction split** is the share that came back
//! *shorter* rather than longer, and the **fall-off** is how much rarer a two-unit slip is than a
//! one-unit slip. Both are proportions; the third number, the slippage *level*, is how often a
//! read slips at all and lives in [`slippage_curve`](super::slippage_curve).
//!
//! **The problem this exists for.** A share was either measured on the stratum's own 4,000
//! slipped reads or copied whole from the nearest stratum that had them. Measured on both
//! cohorts, the best-measured stratum of a motif period reaches 4,000 at **one period out of the
//! twelve** — HG002's homopolymers, with 8,840 — so 69 of HG002's strata and every one of
//! tomato's got nothing at all. A curve fitted from *every* stratum, each weighted by how
//! precisely it measured its own share, has no gate to fail
//! (`doc/devel/ng/spec/str_slippage_level_curve.md` §5.1).
//!
//! **Why the level's curve cannot simply be reused.** That family is monotone and refuses a
//! falling fit, because every measurement says slippage rises with repeat count. The shares do no
//! such thing: HG002's dinucleotide direction split falls from 0.80 at 6 repeats to 0.22 at 10
//! and climbs to 0.98 by 23. So the shape is chosen from three candidates by measurement, per
//! motif period and per share — and all three win somewhere
//! (`doc/devel/ng/reports/str_slippage_share_families_2026-08-20.md`).
//!
//! **Everything here happens on the logit scale**, `log(p / (1 − p))`. A share is a proportion,
//! and on that scale a fitted curve may take any value at all while the share it maps back to
//! still lies strictly between 0 and 1 — so unlike the level, no fitted share ever has to be
//! clamped into being a probability.
//!
//! **A curve always comes back.** Where a period has too few strata to choose a shape it gets a
//! flat one; where it has none it gets the run's own average over the other periods; where the
//! run fitted nothing anywhere it gets a built-in default. These numbers are a prior the read
//! likelihood consults, not a measurement anyone reports, so refusing to answer is worse than
//! answering coarsely — and [`ShareCurve::source`] says which of the four happened.
//!
//! Design: `doc/devel/ng/spec/str_slippage_level_curve.md` §5.1. Build order:
//! `doc/devel/ng/impl_plan/str_slippage_level_curve.md`, Milestone E.

use super::slippage_curve::CurveReach;

// ---------------------------------------------------------------------
// The shape a share follows
// ---------------------------------------------------------------------

/// Which shape one share follows across repeat count, on the logit scale.
///
/// **All three are fitted and one is chosen by measurement**, per motif period and per share,
/// because no single shape is right everywhere. Over the ten cases both cohorts can speak to —
/// four motif periods on HG002, one on tomato, two shares each — [`ShareShape::Flat`] wins six,
/// [`ShareShape::Sloping`] two and [`ShareShape::Turning`] two
/// (`tests/share_curve_on_real_cells.rs` asserts every one of them).
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

impl ShareShape {
    /// The shapes the fit tries, simplest first. **The order is the tie rule**: an equal
    /// held-out score keeps the shape already held, so the simplest shape wins a tie.
    pub const ALL: [Self; 3] = [Self::Flat, Self::Sloping, Self::Turning];

    /// How many coefficients this shape fits. A shape cannot be fitted from fewer strata than
    /// this, and cannot be *scored* from fewer than one more.
    pub fn coefficients(self) -> usize {
        match self {
            Self::Flat => 1,
            Self::Sloping => 2,
            Self::Turning => 3,
        }
    }
}

// ---------------------------------------------------------------------
// One stratum's own answer, and how precisely it holds it
// ---------------------------------------------------------------------

/// One stratum's own fitted share, as the curve fit sees it.
///
/// **Only a stratum that fitted this share from its own tracts belongs here.** A share copied
/// from a neighbour, or supplied by a curve, would make the next fit a fit to the previous fit's
/// output — the circularity the design forbids in so many words: *a stratum feeds its curve only
/// through its own fit, never a blended one*.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedShare {
    /// The reference tract's repeat count — the axis the curve runs along.
    pub repeats: u64,
    /// The share itself, a proportion. Values at or beyond the ends are pulled just inside them
    /// before the logarithm; see [`SHARE_FLOOR`].
    pub share: f64,
    /// **How many of the stratum's reads its fitted level says slipped** — `level × reads
    /// crossing`. The two shares are proportions *over the reads that slipped*, so this and not
    /// the read count is what sets how precisely the stratum holds them.
    pub slipped_reads: f64,
}

impl FittedShare {
    /// Whether this stratum says anything about where the curve runs.
    ///
    /// **A stratum with no slipped reads behind it does not**, and it is dropped from the fit and
    /// from the score alike — scoring a curve on how well it predicts a share nothing measured
    /// would charge it for the fit's own noise, and that score is what the blend reads as the
    /// curve's precision.
    pub fn is_usable(&self) -> bool {
        self.share.is_finite() && self.slipped_reads > 0.0
    }

    /// The share pulled inside `(0, 1)`, so its logarithm exists.
    pub fn bounded(&self) -> f64 {
        self.share.clamp(SHARE_FLOOR, SHARE_CEILING)
    }

    /// This stratum's share on the logit scale.
    pub fn logit(&self) -> f64 {
        let share = self.bounded();
        (share / (1.0 - share)).ln()
    }

    /// **How precisely this stratum holds its own share, in the units the curve is fitted in.**
    ///
    /// A proportion `p` measured on `S` slipped reads has variance `p(1 − p) / S`; carried onto
    /// the logit scale that is `1 / (S · p · (1 − p))`, so the standard error is the square root
    /// of it. **This is the number that weights the fit** — inverse variance of the quantity
    /// being fitted, which is the logit rather than the share.
    ///
    /// **Fewer than one slipped read is read as one**, so a stratum that measured nothing has a
    /// wide error rather than an infinite one and weighing it against a curve stays arithmetic.
    pub fn logit_standard_error(&self) -> f64 {
        let share = self.bounded();
        1.0 / (self.slipped_reads.max(1.0) * share * (1.0 - share)).sqrt()
    }

    /// The weight this stratum carries in the fit: one over its squared logit standard error,
    /// which is `slipped reads × p × (1 − p)`.
    ///
    /// **A thin stratum contributes without swamping anything.** At HG002's homopolymers a
    /// stratum with 10 slipped reads at a split of 0.6 carries a weight of 2.4 against a
    /// well-measured stratum's 2,100 — about one part in 900 — which is what makes "every
    /// stratum contributes" a statement about weight rather than about admission.
    pub fn weight(&self) -> f64 {
        let error = self.logit_standard_error();
        1.0 / (error * error)
    }

    /// How precisely the stratum holds the share **as a share of the share itself** —
    /// `sqrt((1 − p) / (p · S))`.
    ///
    /// **This is not what weights the fit, and it is here because it is where the floor of 4,000
    /// slipped reads came from.** At the values `parameter_prepass_ssr.md` §3 measures, holding a
    /// direction split of 0.17 to 6% of itself takes 1,356 slipped reads where the architecture
    /// says "about 1,400", and holding a fall-off of 0.065 to the same takes 3,996 where it says
    /// "about 4,000".
    ///
    /// **It disagrees with [`FittedShare::logit_standard_error`] by a factor of `1 / (1 − p)`**,
    /// which is nothing for a share near zero and everything for one near one: at HG002's
    /// dinucleotides a stratum splitting 0.983 on 1,940 slipped reads would carry 173 times the
    /// weight of one splitting 0.223 on 2,262 under this measure, and one twelfth of it under the
    /// other. The fit runs on the logit, so the logit's error is the one that weights it.
    pub fn relative_standard_error(&self) -> f64 {
        let share = self.bounded();
        ((1.0 - share) / (share * self.slipped_reads.max(1.0))).sqrt()
    }
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
    /// **Logit units, because that is what it will be weighed against**: a stratum's own error
    /// ([`FittedShare::logit_standard_error`]) is in the same units, so the two combine without a
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
    /// The run fitted no stratum anywhere, so the curve is a built-in constant; see
    /// [`DEFAULT_SHORTER_SHARE`] and [`DEFAULT_FALL_OFF`]. **The crudest number in the parameter
    /// set**, and the only one that is not this run's own data in any form.
    BuiltInDefault,
}

impl ShareCurve {
    /// The share this curve gives at `repeats`.
    ///
    /// **Inside the repeat counts it was fitted over it is the fitted curve; outside, it is held
    /// at the nearer fitted end.** A turning shape continued past its range runs to 0 or 1 within
    /// a few repeat counts, and neither is a share any stratum measured.
    pub fn share_at(&self, repeats: u64) -> f64 {
        let inside = repeats.clamp(self.fitted_from, self.fitted_to);
        self.share_on_the_curve(inside as f64)
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

    /// The fitted curve read at any repeat count, with no regard for the fitted range.
    ///
    /// **Used by the score and not by a consumer.** Leaving out the lowest or highest stratum
    /// puts it outside what remains, so scoring it the way a deployed curve behaves would predict
    /// it with its neighbour's value — which says nothing about the shape, and the end strata are
    /// exactly where a shape shows. This is the same rule the level's curve is scored under.
    fn share_on_the_curve(&self, repeats: f64) -> f64 {
        let from_centre = repeats - self.centre;
        let logit =
            self.intercept + self.slope * from_centre + self.bend * from_centre * from_centre;
        let share = 1.0 / (1.0 + (-logit).exp());
        if share.is_finite() {
            share.clamp(SHARE_FLOOR, SHARE_CEILING)
        } else if logit > 0.0 {
            SHARE_CEILING
        } else {
            SHARE_FLOOR
        }
    }
}

/// Why a set of strata produced no curve of a given shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoShareCurve {
    /// Fewer usable strata than this shape has coefficients.
    TooFewStrata { strata: usize, needs: usize },
    /// The strata do not spread over enough distinct repeat counts to determine the shape — two
    /// points cannot pin a curve that turns.
    TooFewRepeatCounts { repeat_counts: usize, needs: usize },
    /// The normal equations had no solution, or one that is not a number.
    NoSolution,
}

// ---------------------------------------------------------------------
// Fitting one shape
// ---------------------------------------------------------------------

/// Fit one shape through a period's strata, weighted by how precisely each holds its own share.
///
/// **Weighted least squares of the logit of the share on repeat count**, with each stratum
/// weighted by `slipped reads × p × (1 − p)` — the inverse variance of the quantity being
/// fitted. The `held_out_error` of the returned curve is a placeholder that
/// [`choose_share_shape`] fills; a caller fitting one shape directly gets zero and must not read
/// it as a claim.
pub fn fit_share_curve(
    strata: &[FittedShare],
    shape: ShareShape,
) -> Result<ShareCurve, NoShareCurve> {
    let usable: Vec<&FittedShare> = strata
        .iter()
        .filter(|stratum| stratum.is_usable())
        .collect();
    let coefficients = shape.coefficients();
    if usable.len() < coefficients {
        return Err(NoShareCurve::TooFewStrata {
            strata: usable.len(),
            needs: coefficients,
        });
    }

    let mut repeat_counts: Vec<u64> = usable.iter().map(|stratum| stratum.repeats).collect();
    repeat_counts.sort_unstable();
    repeat_counts.dedup();
    if repeat_counts.len() < coefficients {
        return Err(NoShareCurve::TooFewRepeatCounts {
            repeat_counts: repeat_counts.len(),
            needs: coefficients,
        });
    }

    let total_weight: f64 = usable.iter().map(|stratum| stratum.weight()).sum();
    if !(total_weight.is_finite() && total_weight > 0.0) {
        return Err(NoShareCurve::NoSolution);
    }
    let centre: f64 = usable
        .iter()
        .map(|stratum| stratum.weight() * stratum.repeats as f64)
        .sum::<f64>()
        / total_weight;

    // The normal equations, `(Xᵀ W X) c = Xᵀ W y`, built as one augmented matrix. The design
    // matrix's columns are 1, (repeats − centre) and its square, as far as the shape goes.
    let mut augmented = [[0.0_f64; 4]; 3];
    for stratum in &usable {
        let weight = stratum.weight();
        let from_centre = stratum.repeats as f64 - centre;
        let terms = [1.0, from_centre, from_centre * from_centre];
        let target = stratum.logit();
        for row in 0..coefficients {
            for column in 0..coefficients {
                augmented[row][column] += weight * terms[row] * terms[column];
            }
            augmented[row][coefficients] += weight * terms[row] * target;
        }
    }
    let Some(solution) = solve(&mut augmented, coefficients) else {
        return Err(NoShareCurve::NoSolution);
    };
    if !solution.iter().all(|value| value.is_finite()) {
        return Err(NoShareCurve::NoSolution);
    }

    Ok(ShareCurve {
        shape,
        intercept: solution[0],
        slope: solution[1],
        bend: solution[2],
        centre,
        fitted_from: *repeat_counts.first().expect("one repeat count at least"),
        fitted_to: *repeat_counts.last().expect("one repeat count at least"),
        held_out_error: 0.0,
        strata: usable.len(),
        source: ShareCurveSource::ThisPeriod,
    })
}

/// Gaussian elimination with partial pivoting on an augmented `size × (size + 1)` matrix.
///
/// Returns the coefficients, padded with zeros to three. `None` when the matrix is singular —
/// which here means the strata do not determine the shape, not that anything went wrong.
fn solve(augmented: &mut [[f64; 4]; 3], size: usize) -> Option<[f64; 3]> {
    for step in 0..size {
        let pivot = (step..size)
            .max_by(|left, right| {
                augmented[*left][step]
                    .abs()
                    .total_cmp(&augmented[*right][step].abs())
            })
            .expect("a non-empty range");
        if augmented[pivot][step].abs() < 1e-12 {
            return None;
        }
        augmented.swap(step, pivot);
        let pivot_row = augmented[step];
        for row in augmented.iter_mut().take(size).skip(step + 1) {
            let factor = row[step] / pivot_row[step];
            for (column, value) in row.iter_mut().enumerate().take(size + 1).skip(step) {
                *value -= factor * pivot_row[column];
            }
        }
    }
    let mut solution = [0.0_f64; 3];
    for row in (0..size).rev() {
        let mut value = augmented[row][size];
        for column in (row + 1)..size {
            value -= augmented[row][column] * solution[column];
        }
        solution[row] = value / augmented[row][row];
    }
    Some(solution)
}

// ---------------------------------------------------------------------
// Choosing the shape
// ---------------------------------------------------------------------

/// Choose a shape for one period's strata and return the curve fitted at it.
///
/// **The shape is chosen on strata it did not see.** Each shape is fitted, then every stratum is
/// left out in turn, the shape refitted without it, and the stratum predicted; the shape with the
/// lowest median error in logit units wins. Fit quality on the strata a shape *did* see decides
/// nothing — a turning shape can always follow what it saw more closely than a flat one.
///
/// **Ties go to the simplest shape**, since the three fail in the same direction rather than in
/// opposite ones: the shapes are tried simplest first and an equal score does not displace the
/// shape already held.
///
/// `None` when no shape could be fitted and scored at all — fewer strata than
/// `min_strata_for_a_curve`, or every leave-one-out fit refusing.
pub fn choose_share_shape(strata: &[FittedShare], config: &ShareCurveConfig) -> Option<ShareCurve> {
    // **Counted over the strata that say something**, which is the set the fit and the score both
    // use, so a period cannot reach the floor on strata that are then dropped.
    let usable = strata.iter().filter(|stratum| stratum.is_usable()).count();
    if usable < config.min_strata_for_a_curve {
        return None;
    }
    let mut best: Option<(f64, ShareCurve)> = None;
    for shape in ShareShape::ALL {
        let Ok(mut curve) = fit_share_curve(strata, shape) else {
            continue;
        };
        let Some(held_out_error) = held_out_error_of(strata, shape) else {
            continue;
        };
        curve.held_out_error = held_out_error.max(SHARE_CURVE_ERROR_FLOOR);
        if best
            .as_ref()
            .is_none_or(|(best_error, _)| held_out_error < *best_error)
        {
            best = Some((held_out_error, curve));
        }
    }
    best.map(|(_, curve)| curve)
}

/// The median error, in logit units, of leaving each stratum out in turn and predicting it.
///
/// `None` when not one stratum could be predicted — every leave-one-out fit refused, which
/// happens when a shape has exactly as many strata as coefficients.
fn held_out_error_of(strata: &[FittedShare], shape: ShareShape) -> Option<f64> {
    let mut errors: Vec<f64> = Vec::with_capacity(strata.len());
    for held in 0..strata.len() {
        if !strata[held].is_usable() {
            continue;
        }
        let rest: Vec<FittedShare> = strata
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != held)
            .map(|(_, stratum)| *stratum)
            .collect();
        let Ok(curve) = fit_share_curve(&rest, shape) else {
            continue;
        };
        let predicted = curve.share_on_the_curve(strata[held].repeats as f64);
        let predicted_logit = (predicted / (1.0 - predicted)).ln();
        let error = (predicted_logit - strata[held].logit()).abs();
        if error.is_finite() {
            errors.push(error);
        }
    }
    if errors.is_empty() {
        return None;
    }
    errors.sort_by(f64::total_cmp);
    let middle = errors.len() / 2;
    Some(if errors.len().is_multiple_of(2) {
        (errors[middle - 1] + errors[middle]) / 2.0
    } else {
        errors[middle]
    })
}

// ---------------------------------------------------------------------
// The ladder: a curve always comes back
// ---------------------------------------------------------------------

/// The curve one motif period's strata get for one share — **always one, however little there is
/// to go on**.
///
/// Four rungs, in order, each recorded in [`ShareCurve::source`]:
///
/// 1. **This period's strata, shape chosen by measurement** — enough of them to leave one out.
/// 2. **This period's strata, flat and unscored** — one to three of them. Their weighted mean,
///    credited with [`UNSCORED_SHARE_CURVE_ERROR`] so it carries less weight than any curve that
///    could be scored.
/// 3. **The run's other periods, flat** — this period has nothing. Worse than staying inside a
///    period and recorded as such.
/// 4. **`fallback`** — the run fitted nothing anywhere.
///
/// `elsewhere` is every fitted stratum of the same run at *other* motif periods, and `fallback`
/// the built-in constant for this share ([`DEFAULT_SHORTER_SHARE`] or [`DEFAULT_FALL_OFF`]).
pub fn share_curve_for_a_period(
    at_this_period: &[FittedShare],
    elsewhere: &[FittedShare],
    fallback: f64,
    config: &ShareCurveConfig,
) -> ShareCurve {
    if let Some(curve) = choose_share_shape(at_this_period, config) {
        return curve;
    }
    if let Ok(mut flat) = fit_share_curve(at_this_period, ShareShape::Flat) {
        flat.held_out_error = UNSCORED_SHARE_CURVE_ERROR;
        flat.source = ShareCurveSource::ThisPeriodUnscored;
        return flat;
    }
    if let Ok(mut flat) = fit_share_curve(elsewhere, ShareShape::Flat) {
        flat.held_out_error = UNSCORED_SHARE_CURVE_ERROR;
        flat.source = ShareCurveSource::OtherPeriods;
        // A curve borrowed from other periods says nothing about *this* period's repeat counts,
        // so it is flat everywhere rather than held at the ends of a range it never saw.
        flat.fitted_from = u64::MIN;
        flat.fitted_to = u64::MAX;
        return flat;
    }
    let bounded = fallback.clamp(SHARE_FLOOR, SHARE_CEILING);
    ShareCurve {
        shape: ShareShape::Flat,
        intercept: (bounded / (1.0 - bounded)).ln(),
        slope: 0.0,
        bend: 0.0,
        centre: 0.0,
        fitted_from: u64::MIN,
        fitted_to: u64::MAX,
        held_out_error: UNSCORED_SHARE_CURVE_ERROR,
        strata: 0,
        source: ShareCurveSource::BuiltInDefault,
    }
}

// ---------------------------------------------------------------------
// The knobs, and where each number comes from
// ---------------------------------------------------------------------

/// How few contributing strata leave a period without a chosen shape.
///
/// **4, the same floor the level's curve uses, and for the same arithmetic reason** — the
/// smallest count at which leaving one stratum out still leaves a fit and a spare. It is also the
/// smallest at which the turning shape can be scored at all.
pub const MIN_STRATA_FOR_A_SHARE_CURVE: usize = 4;

/// The smallest share a curve or a stratum may report, so its logit exists.
///
/// **A share of exactly 0 or 1 is a real fitted value, not a defect**: HG002's one fitted
/// pentanucleotide stratum splits 1.000 — every one of its 39 slipped reads came back short. On
/// the logit scale that is infinity, so it is pulled one part in ten thousand inside the boundary,
/// which is far tighter than any share is measured and leaves the stratum with almost no weight —
/// which is the right answer for a share pinned against its boundary.
pub const SHARE_FLOOR: f64 = 1e-4;

/// The largest share a curve or a stratum may report; see [`SHARE_FLOOR`].
pub const SHARE_CEILING: f64 = 1.0 - SHARE_FLOOR;

/// The smallest held-out error a share curve is credited with, in logit units.
///
/// **A curve fitted through strata that happen to lie exactly on it scores zero, and zero error
/// is infinite weight** — it would outweigh any stratum however well measured. The smallest
/// error either cohort produces is 0.167 logit units, so one part in a thousand is 167 times
/// below anything measured: the floor cannot bind on a real fit and exists to keep a drawn or
/// degenerate one from swamping the arithmetic.
pub const SHARE_CURVE_ERROR_FLOOR: f64 = 1e-3;

/// What a curve nobody could score is credited with, in logit units.
///
/// **One logit unit — a factor of `e` in the odds.** It is wider than every scored curve either
/// cohort produced (0.167 to 0.842 across the ten cases), so a curve that could not be scored
/// always carries less weight than one that could, without being so wide that it stops
/// contributing at a stratum that measured nothing itself.
pub const UNSCORED_SHARE_CURVE_ERROR: f64 = 1.0;

/// The direction split a run gets when it fitted no stratum anywhere.
///
/// **0.6 — of the reads that slipped, three in five came back short.** Both cohorts' own values,
/// pooled over every fitted stratum and weighted, sit either side of it: 0.641 on HG002 and 0.592
/// on tomato. It is reached only when nothing in the run was fitted, which on the range this
/// caller has to work over means a single sample too shallow for any stratum to clear the fit
/// floor.
pub const DEFAULT_SHORTER_SHARE: f64 = 0.6;

/// The fall-off a run gets when it fitted no stratum anywhere.
///
/// **0.5 — a two-unit slip is half as likely as a one-unit slip.** The two cohorts' pooled
/// weighted values are 0.405 on HG002 and 0.689 on tomato, and this sits between them. *It is the
/// crudest number in the parameter set and the two cohorts do not agree well about it*, which is
/// the price of answering at all where nothing was measured.
pub const DEFAULT_FALL_OFF: f64 = 0.5;

/// What a run may change about how a share's curve is fitted.
///
/// **The weight each stratum carries is deliberately not here.** It is the inverse variance of
/// the logit, which is what a weighted fit's weight has to be; a knob nobody sets is how a
/// settled decision gets changed by accident.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShareCurveConfig {
    /// Below this many contributing strata a period gets no chosen shape, only a flat mean; see
    /// [`MIN_STRATA_FOR_A_SHARE_CURVE`].
    pub min_strata_for_a_curve: usize,
}

impl Default for ShareCurveConfig {
    fn default() -> Self {
        Self {
            min_strata_for_a_curve: MIN_STRATA_FOR_A_SHARE_CURVE,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A stratum at `repeats` whose share is `share`, measured on `slipped` slipped reads.
    fn stratum(repeats: u64, share: f64, slipped: f64) -> FittedShare {
        FittedShare {
            repeats,
            share,
            slipped_reads: slipped,
        }
    }

    fn logit(share: f64) -> f64 {
        (share / (1.0 - share)).ln()
    }

    fn expit(value: f64) -> f64 {
        1.0 / (1.0 + (-value).exp())
    }

    /// Strata lying exactly on a logit line, drawn over eight repeat counts.
    fn drawn_on_a_line(intercept: f64, slope: f64) -> Vec<FittedShare> {
        (6..14)
            .map(|repeats| stratum(repeats, expit(intercept + slope * repeats as f64), 2_000.0))
            .collect()
    }

    // -----------------------------------------------------------------
    // Recovering a shape the strata were drawn with
    // -----------------------------------------------------------------

    /// A share that does not move with repeat count comes back flat, at the value it was drawn
    /// at — and the two bending shapes do not win by fitting the strata they saw more closely.
    #[test]
    fn a_share_that_does_not_move_comes_back_flat() {
        let strata: Vec<FittedShare> = (6..14).map(|n| stratum(n, 0.62, 2_000.0)).collect();
        let curve = choose_share_shape(&strata, &ShareCurveConfig::default())
            .expect("eight strata are enough to choose between shapes");
        assert_eq!(curve.shape, ShareShape::Flat);
        assert_eq!(curve.source, ShareCurveSource::ThisPeriod);
        assert!(
            (curve.share_at(9) - 0.62).abs() < 1e-9,
            "the flat curve came back at {}",
            curve.share_at(9)
        );
    }

    /// A share drawn on a logit line is recovered as a sloping one, and the line's own
    /// coefficients come back.
    #[test]
    fn a_share_drawn_on_a_logit_line_is_recovered() {
        let strata = drawn_on_a_line(-2.0, 0.25);
        let curve = choose_share_shape(&strata, &ShareCurveConfig::default())
            .expect("eight strata over eight repeat counts");
        assert_eq!(curve.shape, ShareShape::Sloping);
        assert!(
            (curve.slope - 0.25).abs() < 1e-9,
            "the slope came back as {}",
            curve.slope
        );
        for repeats in 6..14 {
            let drawn = expit(-2.0 + 0.25 * repeats as f64);
            assert!(
                (curve.share_at(repeats) - drawn).abs() < 1e-9,
                "at {repeats} repeats the curve says {} and the strata were drawn at {drawn}",
                curve.share_at(repeats)
            );
        }
    }

    /// A share that falls and then rises — the shape HG002's dinucleotide direction split has —
    /// needs the turning shape, and neither simpler one can follow it.
    #[test]
    fn a_share_that_turns_once_needs_the_turning_shape() {
        let strata: Vec<FittedShare> = (6..16)
            .map(|repeats| {
                let from_centre = repeats as f64 - 11.0;
                stratum(
                    repeats,
                    expit(-1.0 + 0.08 * from_centre * from_centre),
                    2_000.0,
                )
            })
            .collect();
        let curve = choose_share_shape(&strata, &ShareCurveConfig::default())
            .expect("ten strata over ten repeat counts");
        assert_eq!(curve.shape, ShareShape::Turning);
        assert!(
            (curve.bend - 0.08).abs() < 1e-9,
            "the bend came back as {}",
            curve.bend
        );
        // The turn is where it was drawn: the lowest share sits at the middle of the range.
        let lowest = (6..16)
            .min_by(|left, right| curve.share_at(*left).total_cmp(&curve.share_at(*right)))
            .expect("a non-empty range");
        assert_eq!(lowest, 11);
    }

    /// **The tie rule.** Strata that lie exactly on a flat line are fitted equally well by all
    /// three shapes, and the simplest one is kept.
    #[test]
    fn an_exact_tie_between_shapes_goes_to_the_simplest() {
        let strata: Vec<FittedShare> = (6..12).map(|n| stratum(n, 0.4, 5_000.0)).collect();
        for shape in ShareShape::ALL {
            let held_out = held_out_error_of(&strata, shape).expect("every shape is scorable");
            assert!(held_out < 1e-9, "{shape:?} scored {held_out}");
        }
        let curve = choose_share_shape(&strata, &ShareCurveConfig::default()).expect("six strata");
        assert_eq!(curve.shape, ShareShape::Flat);
    }

    // -----------------------------------------------------------------
    // What a curve gives outside what it saw, and at the boundaries
    // -----------------------------------------------------------------

    /// Outside the repeat counts a curve was fitted over, the share is held at the nearer end
    /// rather than continued — a turning shape continued a few counts past its range runs to 0
    /// or 1, and neither is a share any stratum measured.
    #[test]
    fn outside_the_fitted_range_the_share_is_held_at_the_nearer_end() {
        let strata: Vec<FittedShare> = (8..13)
            .map(|repeats| {
                let from_centre = repeats as f64 - 10.0;
                stratum(
                    repeats,
                    expit(0.5 + 0.6 * from_centre * from_centre),
                    3_000.0,
                )
            })
            .collect();
        let curve = fit_share_curve(&strata, ShareShape::Turning).expect("five strata");
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 12));

        assert_eq!(curve.share_at(3), curve.share_at(8));
        assert_eq!(curve.share_at(60), curve.share_at(12));
        assert_eq!(curve.reach(3), CurveReach::BelowFitted);
        assert_eq!(curve.reach(10), CurveReach::Inside);
        assert_eq!(curve.reach(60), CurveReach::AboveFitted);

        // Continued instead of held, the same curve leaves the range of anything measured.
        assert!(curve.share_on_the_curve(60.0) > SHARE_CEILING - 1e-12);
    }

    /// Every share a curve reports is a proportion, at every repeat count, for the shape most
    /// able to leave the range.
    #[test]
    fn an_emitted_share_is_always_a_proportion() {
        let strata: Vec<FittedShare> = (6..12)
            .map(|repeats| {
                let from_centre = repeats as f64 - 9.0;
                stratum(repeats, expit(-3.0 * from_centre), 1_000.0)
            })
            .collect();
        for shape in ShareShape::ALL {
            let Ok(curve) = fit_share_curve(&strata, shape) else {
                continue;
            };
            for repeats in 1..200 {
                let share = curve.share_at(repeats);
                assert!(
                    (0.0..=1.0).contains(&share) && share.is_finite(),
                    "{shape:?} gave {share} at {repeats} repeats"
                );
            }
        }
    }

    /// **A share of exactly 1 is a real fitted value, not a defect** — HG002's one fitted
    /// pentanucleotide stratum splits 1.000, every one of its 39 slipped reads coming back
    /// short. It must stay finite, and it must carry almost no weight, because a share pinned
    /// against its boundary says nothing about where the curve runs.
    #[test]
    fn a_share_pinned_at_one_stays_finite_and_carries_almost_no_weight() {
        let pinned = stratum(5, 1.0, 39.0);
        assert!(pinned.logit().is_finite());
        assert!(pinned.logit_standard_error().is_finite());

        let ordinary = stratum(5, 0.6, 39.0);
        assert!(
            pinned.weight() * 100.0 < ordinary.weight(),
            "pinned {} against ordinary {}",
            pinned.weight(),
            ordinary.weight()
        );

        // And a stratum at zero is the same case at the other end.
        assert!(stratum(5, 0.0, 39.0).logit().is_finite());
    }

    // -----------------------------------------------------------------
    // The weight
    // -----------------------------------------------------------------

    /// **What "every stratum contributes" means.** A stratum with ten slipped reads moves the
    /// curve by about one part in nine hundred against a well-measured one, so admitting it
    /// costs nothing and refusing it gains nothing.
    #[test]
    fn a_thin_stratum_moves_the_curve_by_about_one_part_in_nine_hundred() {
        let thin = stratum(9, 0.6, 10.0);
        let full = stratum(9, 0.6, 8_000.0);
        let ratio = full.weight() / thin.weight();
        assert!(
            (780.0..=1_000.0).contains(&ratio),
            "a well-measured stratum outweighs a thin one {ratio:.0} to one"
        );
    }

    /// A stratum pulls the curve in proportion to its weight: moving the share of a
    /// well-measured stratum moves the fitted curve, and moving a thin one's barely does.
    #[test]
    fn a_stratum_pulls_the_curve_in_proportion_to_its_slipped_reads() {
        let base: Vec<FittedShare> = (6..12).map(|n| stratum(n, 0.5, 2_000.0)).collect();
        let flat = fit_share_curve(&base, ShareShape::Flat).expect("six strata");

        let mut with_a_thin_outlier = base.clone();
        with_a_thin_outlier.push(stratum(12, 0.95, 5.0));
        let thin_pull = fit_share_curve(&with_a_thin_outlier, ShareShape::Flat).expect("seven");

        let mut with_a_full_outlier = base.clone();
        with_a_full_outlier.push(stratum(12, 0.95, 5_000.0));
        let full_pull = fit_share_curve(&with_a_full_outlier, ShareShape::Flat).expect("seven");

        let moved_by_thin = (thin_pull.intercept - flat.intercept).abs();
        let moved_by_full = (full_pull.intercept - flat.intercept).abs();
        assert!(
            moved_by_full > moved_by_thin * 50.0,
            "a thin outlier moved the curve by {moved_by_thin:.4} and a full one by \
             {moved_by_full:.4}"
        );
    }

    /// **Where the floor of 4,000 slipped reads came from**, reproduced from the precision model
    /// the design gives: at the values `parameter_prepass_ssr.md` §3 measures, holding a
    /// direction split of 0.17 to 6% of itself takes about 1,400 slipped reads and holding a
    /// fall-off of 0.065 to the same takes about 4,000.
    #[test]
    fn the_relative_precision_model_is_where_four_thousand_slipped_reads_came_from() {
        let reads_to_hold = |share: f64, to_within: f64| {
            let mut slipped = 1.0;
            while stratum(9, share, slipped).relative_standard_error() > to_within {
                slipped += 1.0;
            }
            slipped
        };
        assert_eq!(reads_to_hold(0.17, 0.06), 1_357.0);
        assert_eq!(reads_to_hold(0.065, 0.06), 3_996.0);
    }

    /// **The two precisions rank the same two strata in opposite orders**, which is why the fit
    /// uses the logit's and not the share's: dividing by `p` makes a share near certainty look
    /// best measured, where on the scale the curve is fitted on it is worst.
    #[test]
    fn the_two_precisions_rank_a_near_certain_share_in_opposite_orders() {
        let near_certain = stratum(23, 0.983, 1_940.0);
        let middling = stratum(10, 0.223, 2_262.0);

        let by_relative = |share: &FittedShare| {
            let error = share.relative_standard_error();
            1.0 / (error * error)
        };
        let relative_ratio = by_relative(&near_certain) / by_relative(&middling);
        let logit_ratio = near_certain.weight() / middling.weight();

        assert!(
            (170.0..=176.0).contains(&relative_ratio),
            "by the share's own relative error the near-certain stratum carries \
             {relative_ratio:.0} times the weight"
        );
        assert!(
            (0.075..=0.09).contains(&logit_ratio),
            "on the logit scale it carries {logit_ratio:.3} of it"
        );
    }

    // -----------------------------------------------------------------
    // The score
    // -----------------------------------------------------------------

    /// **Fit quality on the strata a shape saw decides nothing, and this is why.** Over strata
    /// scattered around a flat truth, the turning shape follows what it saw more closely than
    /// either simpler one — it has three coefficients to their one — and yet it predicts a
    /// stratum it did not see worse. The score is the held-out one for exactly this reason.
    #[test]
    fn the_shape_that_fits_what_it_saw_best_is_not_the_shape_that_wins() {
        // A flat truth of 0.5 with a fixed scatter, so the test draws the same strata every run.
        let scatter = [0.30, -0.22, 0.26, -0.35, 0.14, -0.27, 0.31, -0.17];
        let strata: Vec<FittedShare> = scatter
            .iter()
            .enumerate()
            .map(|(index, offset)| stratum(6 + index as u64, expit(*offset), 2_000.0))
            .collect();

        let residual_on_what_it_saw = |shape: ShareShape| {
            let curve = fit_share_curve(&strata, shape).expect("eight strata");
            let mut residuals: Vec<f64> = strata
                .iter()
                .map(|held| (logit(curve.share_at(held.repeats)) - held.logit()).abs())
                .collect();
            residuals.sort_by(f64::total_cmp);
            residuals[residuals.len() / 2]
        };

        let flat_own = residual_on_what_it_saw(ShareShape::Flat);
        let turning_own = residual_on_what_it_saw(ShareShape::Turning);
        assert!(
            turning_own < flat_own,
            "the turning shape should follow what it saw more closely: {turning_own:.3} \
             against {flat_own:.3}"
        );

        let flat_held_out = held_out_error_of(&strata, ShareShape::Flat).expect("scorable");
        let turning_held_out = held_out_error_of(&strata, ShareShape::Turning).expect("scorable");
        assert!(
            turning_held_out > flat_held_out,
            "and predict a stratum it did not see worse: {turning_held_out:.3} against \
             {flat_held_out:.3}"
        );

        let curve = choose_share_shape(&strata, &ShareCurveConfig::default()).expect("eight");
        assert_eq!(curve.shape, ShareShape::Flat);
    }

    /// A shape with exactly as many strata as it has coefficients cannot be scored at all — each
    /// leave-one-out fit is refused — so it does not win by default.
    #[test]
    fn a_shape_with_no_room_to_leave_one_out_is_not_scored() {
        let strata: Vec<FittedShare> = (6..9)
            .map(|repeats| stratum(repeats, expit(-1.0 + 0.3 * repeats as f64), 500.0))
            .collect();
        assert_eq!(strata.len(), 3);
        assert!(held_out_error_of(&strata, ShareShape::Turning).is_none());
        assert!(held_out_error_of(&strata, ShareShape::Sloping).is_some());
    }

    /// A curve fitted through strata that lie exactly on it is credited with the error floor
    /// rather than with zero, since zero error is infinite weight.
    #[test]
    fn a_curve_that_fits_perfectly_is_credited_with_the_error_floor() {
        let strata = drawn_on_a_line(-2.0, 0.25);
        let curve = choose_share_shape(&strata, &ShareCurveConfig::default()).expect("eight");
        assert_eq!(curve.held_out_error, SHARE_CURVE_ERROR_FLOOR);
    }

    // -----------------------------------------------------------------
    // Refusals inside the fit
    // -----------------------------------------------------------------

    /// Two repeat counts cannot pin a shape that turns, however many strata sit at them.
    #[test]
    fn two_repeat_counts_cannot_pin_a_shape_that_turns() {
        let strata = vec![
            stratum(8, 0.4, 1_000.0),
            stratum(8, 0.45, 1_000.0),
            stratum(12, 0.6, 1_000.0),
            stratum(12, 0.65, 1_000.0),
        ];
        assert_eq!(
            fit_share_curve(&strata, ShareShape::Turning),
            Err(NoShareCurve::TooFewRepeatCounts {
                repeat_counts: 2,
                needs: 3
            })
        );
        assert!(fit_share_curve(&strata, ShareShape::Sloping).is_ok());
    }

    /// A stratum with no slipped reads behind it has nothing to say and is dropped before the
    /// fit, rather than counted as evidence at whatever share it happens to hold.
    #[test]
    fn a_stratum_with_no_slipped_reads_is_dropped_before_the_fit() {
        let strata = vec![
            stratum(8, 0.4, 1_000.0),
            stratum(9, 0.5, 0.0),
            stratum(10, 0.6, 1_000.0),
        ];
        let curve = fit_share_curve(&strata, ShareShape::Sloping).expect("two usable strata");
        assert_eq!(curve.strata, 2);
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 10));
    }

    // -----------------------------------------------------------------
    // The ladder — a curve always comes back
    // -----------------------------------------------------------------

    /// The intended rung: enough strata at this period to choose a shape between.
    #[test]
    fn a_period_with_enough_strata_gets_its_own_shape() {
        let here = drawn_on_a_line(-2.0, 0.25);
        let curve = share_curve_for_a_period(
            &here,
            &[],
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        assert_eq!(curve.source, ShareCurveSource::ThisPeriod);
        assert_eq!(curve.shape, ShareShape::Sloping);
    }

    /// **A period with one to three strata still gets a curve** — their weighted mean, credited
    /// with a wider error than any curve that could be scored, so a stratum's own answer
    /// outweighs it more easily.
    #[test]
    fn a_period_with_too_few_strata_gets_their_flat_mean() {
        let here = vec![stratum(8, 0.4, 1_000.0), stratum(10, 0.6, 1_000.0)];
        let curve = share_curve_for_a_period(
            &here,
            &[stratum(20, 0.9, 9_000.0)],
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        assert_eq!(curve.source, ShareCurveSource::ThisPeriodUnscored);
        assert_eq!(curve.shape, ShareShape::Flat);
        assert_eq!(curve.strata, 2);
        assert_eq!(curve.held_out_error, UNSCORED_SHARE_CURVE_ERROR);
        // The mean of the two logits, which for these two is the middle share.
        assert!((curve.share_at(9) - 0.5).abs() < 1e-9);
    }

    /// **A period with nothing takes the run's other periods**, flat, and says so — crossing
    /// motif periods is worse than staying inside one and is recorded rather than hidden.
    #[test]
    fn a_period_with_no_strata_takes_the_runs_other_periods() {
        let elsewhere = vec![stratum(8, 0.7, 4_000.0), stratum(12, 0.7, 4_000.0)];
        let curve = share_curve_for_a_period(
            &[],
            &elsewhere,
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        assert_eq!(curve.source, ShareCurveSource::OtherPeriods);
        assert!((curve.share_at(30) - 0.7).abs() < 1e-9);
        // It does not depend on repeat count, so no repeat count is outside it.
        assert_eq!(curve.reach(1), CurveReach::Inside);
        assert_eq!(curve.reach(10_000), CurveReach::Inside);
    }

    /// **A run that fitted nothing anywhere still answers.** This is the bottom of the range
    /// this caller has to work over — one sample too shallow for any stratum to clear the fit
    /// floor — and the read likelihood cannot score a stratum with no shares at all.
    #[test]
    fn a_run_that_fitted_nothing_gets_the_built_in_default() {
        for (fallback, expected) in [(DEFAULT_SHORTER_SHARE, 0.6), (DEFAULT_FALL_OFF, 0.5)] {
            let curve = share_curve_for_a_period(&[], &[], fallback, &ShareCurveConfig::default());
            assert_eq!(curve.source, ShareCurveSource::BuiltInDefault);
            assert_eq!(curve.shape, ShareShape::Flat);
            assert_eq!(curve.strata, 0);
            assert!((curve.share_at(9) - expected).abs() < 1e-9);
            assert!((curve.share_at(50) - expected).abs() < 1e-9);
        }
    }

    /// The four rungs are ordered by how much they know, and the error each is credited with
    /// keeps that order: a scored curve is trusted more than an unscored one.
    #[test]
    fn a_scored_curve_is_credited_with_less_error_than_an_unscored_one() {
        let scored = share_curve_for_a_period(
            &drawn_on_a_line(-2.0, 0.25),
            &[],
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        let unscored = share_curve_for_a_period(
            &[stratum(8, 0.4, 1_000.0)],
            &[],
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        assert!(scored.held_out_error < unscored.held_out_error);
    }

    /// **A stratum with no slipped reads is dropped from the score as well as from the fit.**
    /// Counting it would charge the curve for missing a share nothing measured, and that score is
    /// what the blend reads as the curve's precision.
    #[test]
    fn a_stratum_with_nothing_behind_it_is_left_out_of_the_score_too() {
        let measured: Vec<FittedShare> = (6..12)
            .map(|repeats| stratum(repeats, expit(-1.0 + 0.2 * repeats as f64), 2_000.0))
            .collect();
        let mut with_an_empty_one = measured.clone();
        with_an_empty_one.push(stratum(12, 0.01, 0.0));

        let scored = held_out_error_of(&measured, ShareShape::Sloping).expect("six strata");
        let scored_with_the_empty_one =
            held_out_error_of(&with_an_empty_one, ShareShape::Sloping).expect("six usable");
        assert_eq!(scored, scored_with_the_empty_one);
    }

    /// A period reaches the floor on strata that say something, not on rows in a table: four
    /// strata of which one is empty get the flat mean of the other three, and say so.
    #[test]
    fn a_period_reaches_the_floor_on_strata_that_say_something() {
        let strata = vec![
            stratum(8, 0.40, 1_000.0),
            stratum(9, 0.45, 1_000.0),
            stratum(10, 0.50, 1_000.0),
            stratum(11, 0.55, 0.0),
        ];
        assert!(choose_share_shape(&strata, &ShareCurveConfig::default()).is_none());
        let curve = share_curve_for_a_period(
            &strata,
            &[],
            DEFAULT_SHORTER_SHARE,
            &ShareCurveConfig::default(),
        );
        assert_eq!(curve.source, ShareCurveSource::ThisPeriodUnscored);
        assert_eq!(curve.strata, 3);
    }

    /// The default configuration is the measured constants.
    #[test]
    fn the_default_config_is_the_measured_constants() {
        let config = ShareCurveConfig::default();
        assert_eq!(config.min_strata_for_a_curve, MIN_STRATA_FOR_A_SHARE_CURVE);
        assert_eq!(MIN_STRATA_FOR_A_SHARE_CURVE, 4);
    }
}
