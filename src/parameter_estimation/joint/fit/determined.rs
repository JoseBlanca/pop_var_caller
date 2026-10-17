//! **What the data determine about the frequency density and the duplicated class**, and the verdict
//! on each pair of numbers the data cannot tell apart.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_fit.md` §2.1.3. Build order:
//! `doc/devel/ng/impl_plan/snp_fit_convergence.md`, step C3.
//!
//! # The trade
//!
//! A position segregating at one in a million looks invariant in any cohort we can sequence. So the
//! fit can trade a larger invariant share for more of the Beta's mass near zero — a smaller first
//! shape `a` — at almost no cost in likelihood: the reads fix **how many positions carry an allele at
//! a frequency they can see**, not how many carry one at all. The duplicated class has the same kind
//! of pair: a duplication carried by one plant in a thousand has no carrier among a few samples, so
//! the class's share trades against its carrier Beta's first shape.
//!
//! In each pair, `σ` is the share the trade moves and `S` the part of it the reads see — the share
//! above [`STATED_FREQUENCY`]:
//!
//! | pair | the two numbers | `σ` | `S` |
//! |---|---|---|---|
//! | frequency density | invariant share, shape `a` | `1 − p_invariant − p_fixed_alt` | `σ · P(f > 0.01)`, `f ~ Beta(a, b)` |
//! | duplicated class | duplicated share, carrier shape `a` | the duplicated share `δ` | `δ · P(c > 0.01)`, `c ~ Beta(carrier_a, carrier_b)` |
//!
//! **The test asks what halving or doubling `σ` costs with `S` held** and every other number left to
//! refit. Under two log-likelihood units ([`SEPARATION_LOG_LIKELIHOOD`]) the data cannot tell those
//! apart, and the pair's two numbers are **not separately determined**: they are reported without
//! standard errors and with that mark, and the determined share `S` is the number to quote. With `C`
//! the covariance of the cohort's eight numbers and `∇` a gradient over them:
//!
//! ```text
//! Var(x) = ∇xᵀ C ∇x,   Cov(x, y) = ∇xᵀ C ∇y
//! V_trade    = Var(σ) − Cov(σ, S)² / Var(S)       how far σ can move with S held
//! trade_cost = ½ (ln 2)² σ² / V_trade              (V_trade ≤ 10⁻¹² Var(σ): no trade, cost +∞)
//! ```
//!
//! The pair is also marked, with no trade cost, when a number it is computed from has no standard
//! error — dropped from the inversion as not identified, the density's at one sample and the
//! duplicated class's at one and two, or kept in it with an error wider than its range, as every
//! number of the density's group at two samples and 3 reads — or when `Var(S)` is not a positive
//! finite number.
//!
//! # How each number is computed
//!
//! `C` is [`cohort_covariance`](super::standard_errors::cohort_covariance), from the information the
//! winning start's final pass summed — the matrix the errors come from. `σ` is linear in its group,
//! so its gradient is exact. `S`, the expected heterozygosity and the mean alternative frequency
//! are linear in the shares — the invariant and fixed non-reference shares, the duplicated share —
//! with slopes that depend on the shapes alone, so their slopes along the shares are exact too;
//! along the shapes they are central differences, a step of [`DIFFERENCE_STEP`] of each shape's
//! value, through the same functions that compute their values. A difference would fail along a
//! share: one on its lower bound of 10⁻¹² gives a step below the last place of `σ`. The errors are
//! `√Var` from the same `C` (the delta method), and none where a number of the group has no
//! covariance or where `√Var` is wider than 1, the whole range of each of them.

use std::f64::consts::LN_2;

use crate::genetics::regularised_incomplete_beta;

use super::information::{COHORT_PARAMETERS, cohort};
use super::standard_errors::CohortCovariance;
use super::{DuplicatedPositions, FrequencyDensity};

/// **The frequency a share must lie above for the reads to be said to see it**: 1%, for both pairs.
/// The owner's choice for the density (2026-10-16), and soft — at 59 samples one copy among 118
/// chromosomes is 0.85%. A stated constant rather than one copy in the panel, so that the number
/// means the same at every cohort size and two runs can be compared. The duplicated class uses the
/// same constant: a carrier share is a frequency in the panel as `f` is in the population.
pub const STATED_FREQUENCY: f64 = 0.01;

/// **The log-likelihood that halving or doubling `σ` must cost, `S` held, for the pair's two
/// numbers to count as separately determined**: two units, about the half-width of a 95%
/// likelihood interval for one number (1.92). Soft.
pub const SEPARATION_LOG_LIKELIHOOD: f64 = 2.0;

/// The central difference's step, as a share of the number it is taken on.
const DIFFERENCE_STEP: f64 = 1e-5;

/// **The share of `σ`'s own variance at or below which what is left of it, `S` held, is rounding**:
/// 10⁻¹², ten thousand units in the last place of the variance it is the difference from. A trade
/// that small would cost `½ (ln 2)² / 10⁻¹² ≈ 2.4 × 10¹¹` units times `σ² / Var(σ)`, past the line
/// of two unless `σ`'s standard error were over 300,000 times `σ` itself.
const NO_TRADE_SHARE: f64 = 1e-12;

/// **One number as the fit reports it**: its value, and its standard error where it has one — none
/// when the number is marked not separately determined, or when the fit gives it none.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ReportedNumber {
    pub value: f64,
    pub standard_error: Option<f64>,
}

/// **What the data determine about one of the two pairs**: the share above [`STATED_FREQUENCY`], what
/// trading `σ` against it costs, and the verdict.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DeterminedShare {
    /// `S`: the share above [`STATED_FREQUENCY`] — of positions segregating there, for the density;
    /// of positions duplicated in more than that share of the panel, for the duplicated class.
    pub share_above: ReportedNumber,
    /// The log-likelihood lost halving or doubling `σ` with `S` held: `+∞` when there is no trade;
    /// `None` when a number of the group has no standard error or `Var(S)` is not a positive finite
    /// number, so the cost cannot be computed.
    pub trade_cost: Option<f64>,
    /// `false`: the pair's two numbers are written without standard errors and with the mark.
    pub separately_determined: bool,
}

/// **The frequency density as the fit reports it**: its four numbers, the two integrals the
/// genotype prior's seed is built from, and what the data determine about the trade between the
/// invariant share and the first shape.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DensityReport {
    /// One of the pair: no error when the pair is not separately determined.
    pub invariant_share: ReportedNumber,
    pub fixed_non_reference_share: ReportedNumber,
    /// One of the pair: no error when the pair is not separately determined.
    pub segregating_shape_a: ReportedNumber,
    pub segregating_shape_b: ReportedNumber,
    /// `Hexp`, [`FrequencyDensity::expected_heterozygosity`].
    pub expected_heterozygosity: ReportedNumber,
    /// `E[f]`, [`FrequencyDensity::expected_alternative_frequency`].
    pub expected_alternative_frequency: ReportedNumber,
    pub determined: DeterminedShare,
}

/// **The duplicated class as the fit reports it**: its three numbers, and what the data determine
/// about the trade between its share and the carrier Beta's first shape.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct DuplicatedReport {
    /// One of the pair: no error when the pair is not separately determined.
    pub duplicated_share: ReportedNumber,
    /// One of the pair: no error when the pair is not separately determined.
    pub carrier_shape_a: ReportedNumber,
    pub carrier_shape_b: ReportedNumber,
    pub determined: DeterminedShare,
}

impl DensityReport {
    /// The report of `density`, from the cohort's standard errors `errors` (by slot of
    /// [`cohort`], `None` where a number has none) and `covariance`,
    /// both from the information the fit's errors come from.
    pub(super) fn of(
        density: &FrequencyDensity,
        errors: &[Option<f64>; COHORT_PARAMETERS],
        covariance: &CohortCovariance,
    ) -> Self {
        let slots = [
            cohort::P_INVARIANT,
            cohort::P_FIXED_ALT,
            cohort::DENSITY_A,
            cohort::DENSITY_B,
        ];
        let group = Group::of(
            slots,
            [
                density.p_invariant,
                density.p_fixed_alt,
                density.a,
                density.b,
            ],
            errors,
            covariance,
        );
        let as_density = |values: &[f64; 4]| FrequencyDensity {
            p_invariant: values[0],
            p_fixed_alt: values[1],
            a: values[2],
            b: values[3],
        };
        let segregating = |values: &[f64; 4]| 1.0 - values[0] - values[1];
        // **The three derived numbers are linear in the two shares**, with slopes set by the
        // shapes alone: a segregating position's chance of lying above 1%, its heterozygosity and
        // its mean frequency. Those slopes are exact; only the shapes are differenced.
        let (a, b) = (density.a, density.b);
        let seen = share_of_beta_above(a, b);
        let per_segregating_position = as_density(&[0.0, 0.0, a, b]);
        let heterozygosity_slope = -per_segregating_position.expected_heterozygosity();
        let mean_frequency = per_segregating_position.expected_alternative_frequency();
        let determined = group.determined_share(
            [-1.0, -1.0, 0.0, 0.0],
            segregating,
            |values| segregating(values) * share_of_beta_above(values[2], values[3]),
            [Some(-seen), Some(-seen), None, None],
        );
        let heterozygosity = |values: &[f64; 4]| as_density(values).expected_heterozygosity();
        let alternative = |values: &[f64; 4]| as_density(values).expected_alternative_frequency();
        let pair_error = |slot: usize| errors[slot].filter(|_| determined.separately_determined);
        Self {
            invariant_share: ReportedNumber {
                value: density.p_invariant,
                standard_error: pair_error(cohort::P_INVARIANT),
            },
            fixed_non_reference_share: ReportedNumber {
                value: density.p_fixed_alt,
                standard_error: errors[cohort::P_FIXED_ALT],
            },
            segregating_shape_a: ReportedNumber {
                value: density.a,
                standard_error: pair_error(cohort::DENSITY_A),
            },
            segregating_shape_b: ReportedNumber {
                value: density.b,
                standard_error: errors[cohort::DENSITY_B],
            },
            expected_heterozygosity: group.reported(
                heterozygosity,
                [
                    Some(heterozygosity_slope),
                    Some(heterozygosity_slope),
                    None,
                    None,
                ],
            ),
            expected_alternative_frequency: group.reported(
                alternative,
                [
                    Some(-mean_frequency),
                    Some(1.0 - mean_frequency),
                    None,
                    None,
                ],
            ),
            determined,
        }
    }
}

impl DuplicatedReport {
    /// The report of `duplicated`, from the cohort's standard errors and covariance as
    /// [`DensityReport::of`] takes them.
    pub(super) fn of(
        duplicated: &DuplicatedPositions,
        errors: &[Option<f64>; COHORT_PARAMETERS],
        covariance: &CohortCovariance,
    ) -> Self {
        let slots = [
            cohort::DUPLICATED_SHARE,
            cohort::CARRIER_A,
            cohort::CARRIER_B,
        ];
        let group = Group::of(
            slots,
            [duplicated.share, duplicated.carrier_a, duplicated.carrier_b],
            errors,
            covariance,
        );
        let determined = group.determined_share(
            [1.0, 0.0, 0.0],
            |values| values[0],
            |values| values[0] * share_of_beta_above(values[1], values[2]),
            [
                Some(share_of_beta_above(
                    duplicated.carrier_a,
                    duplicated.carrier_b,
                )),
                None,
                None,
            ],
        );
        let pair_error = |slot: usize| errors[slot].filter(|_| determined.separately_determined);
        Self {
            duplicated_share: ReportedNumber {
                value: duplicated.share,
                standard_error: pair_error(cohort::DUPLICATED_SHARE),
            },
            carrier_shape_a: ReportedNumber {
                value: duplicated.carrier_a,
                standard_error: pair_error(cohort::CARRIER_A),
            },
            carrier_shape_b: ReportedNumber {
                value: duplicated.carrier_b,
                standard_error: errors[cohort::CARRIER_B],
            },
            determined,
        }
    }
}

/// The share of `Beta(a, b)` above [`STATED_FREQUENCY`].
fn share_of_beta_above(a: f64, b: f64) -> f64 {
    1.0 - regularised_incomplete_beta(STATED_FREQUENCY, a, b)
}

/// **The numbers a pair's quantities are computed from**: their values, their covariance with
/// each other, `None` when any entry is missing — a number of the group dropped from the inversion
/// — and whether every one of them has a standard error.
struct Group<const N: usize> {
    values: [f64; N],
    covariance: Option<[[f64; N]; N]>,
    /// `false` when a number of the group has no standard error: dropped from the inversion, or
    /// kept in it with an error wider than its range. The pair is then marked, with no trade cost.
    every_number_has_an_error: bool,
}

impl<const N: usize> Group<N> {
    fn of(
        slots: [usize; N],
        values: [f64; N],
        errors: &[Option<f64>; COHORT_PARAMETERS],
        covariance: &CohortCovariance,
    ) -> Self {
        let every_number_has_an_error = slots.iter().all(|&slot| errors[slot].is_some());
        let mut within = [[0.0; N]; N];
        for (row, &slot_of_row) in within.iter_mut().zip(&slots) {
            for (entry, &slot_of_column) in row.iter_mut().zip(&slots) {
                match covariance[slot_of_row][slot_of_column] {
                    Some(covariance) => *entry = covariance,
                    None => {
                        return Self {
                            values,
                            covariance: None,
                            every_number_has_an_error,
                        };
                    }
                }
            }
        }
        Self {
            values,
            covariance: Some(within),
            every_number_has_an_error,
        }
    }

    /// **The gradient of `quantity` at the group's values**: `exact[k]`, the slope along the `k`th
    /// number, where the caller knows it — where the quantity is linear in that number — and a
    /// central difference where it is `None`. A share on its lower bound of 10⁻¹² must have its
    /// slope given: the difference's step there is below the last place of the numbers it moves.
    fn gradient(&self, quantity: impl Fn(&[f64; N]) -> f64, exact: [Option<f64>; N]) -> [f64; N] {
        std::array::from_fn(|k| {
            if let Some(slope) = exact[k] {
                return slope;
            }
            let step = DIFFERENCE_STEP * self.values[k].abs();
            let mut above = self.values;
            above[k] += step;
            let mut below = self.values;
            below[k] -= step;
            (quantity(&above) - quantity(&below)) / (2.0 * step)
        })
    }

    /// `leftᵀ C right`: the covariance of two quantities with these gradients; `None` without `C`.
    fn covariance_of(&self, left: &[f64; N], right: &[f64; N]) -> Option<f64> {
        let covariance = self.covariance.as_ref()?;
        Some(
            left.iter()
                .zip(covariance)
                .map(|(l, row)| l * row.iter().zip(right).map(|(c, r)| c * r).sum::<f64>())
                .sum(),
        )
    }

    /// A quantity of the group as the fit reports it: its value, and the delta method's error, from
    /// a gradient that takes the slopes in `exact` as [`Self::gradient`] does.
    fn reported(
        &self,
        quantity: impl Fn(&[f64; N]) -> f64,
        exact: [Option<f64>; N],
    ) -> ReportedNumber {
        let gradient = self.gradient(&quantity, exact);
        ReportedNumber {
            value: quantity(&self.values),
            standard_error: self
                .covariance_of(&gradient, &gradient)
                .and_then(error_from_variance),
        }
    }

    /// **The verdict on the pair** whose traded share `σ` has the exact gradient `share_gradient`
    /// and value `share`, and whose determined share is `share_above`, with the slopes of
    /// `share_above` known exactly given in `above_exact`.
    fn determined_share(
        &self,
        share_gradient: [f64; N],
        share: impl Fn(&[f64; N]) -> f64,
        share_above: impl Fn(&[f64; N]) -> f64,
        above_exact: [Option<f64>; N],
    ) -> DeterminedShare {
        let above_gradient = self.gradient(&share_above, above_exact);
        let variance_above = self.covariance_of(&above_gradient, &above_gradient);
        // **A number of the group with no error leaves no cost to compute.** Its covariance, where
        // it kept one, describes a curvature too shallow to place it anywhere in its range; a cost
        // from it would be a verdict on numbers the data do not determine — at two samples, `+∞`
        // with every number of the density's group so placed (spec §2.1.3, amended 2026-10-17).
        let trade_cost = match (
            self.every_number_has_an_error,
            self.covariance_of(&share_gradient, &share_gradient),
            variance_above,
            self.covariance_of(&share_gradient, &above_gradient),
        ) {
            (true, Some(variance), Some(variance_above), Some(together))
                if variance_above > 0.0 && variance_above.is_finite() =>
            {
                trade_cost_of(
                    share(&self.values),
                    variance,
                    variance - together * together / variance_above,
                )
            }
            _ => None,
        };
        DeterminedShare {
            share_above: ReportedNumber {
                value: share_above(&self.values),
                standard_error: variance_above.and_then(error_from_variance),
            },
            trade_cost,
            separately_determined: is_separately_determined(trade_cost),
        }
    }
}

/// **What halving or doubling the share `σ` costs** when it can move by a variance of `trade_variance`
/// with the determined share held, out of the `variance` it has alone: `½ (ln 2)² σ² / V_trade`.
///
/// **A `V_trade` at or below [`NO_TRADE_SHARE`] of `σ`'s own variance is no trade** and costs
/// `+∞`, never a negative number or a finite one. There `σ` and the determined share move together
/// — on the 4 oracle samples their correlation is 1.0000 — and `V_trade` is the cancellation of two
/// equal numbers: exactly zero, or a unit in the last place either side, which as a cost would be
/// about 10¹⁶. One state, so one answer. `None` when `V_trade` is not a number.
fn trade_cost_of(share: f64, variance: f64, trade_variance: f64) -> Option<f64> {
    if trade_variance.is_nan() {
        None
    } else if trade_variance <= NO_TRADE_SHARE * variance {
        Some(f64::INFINITY)
    } else {
        Some(0.5 * LN_2 * LN_2 * share * share / trade_variance)
    }
}

/// **The verdict**: separately determined when trading costs at least
/// [`SEPARATION_LOG_LIKELIHOOD`]; not when it costs less, or when the cost cannot be computed.
fn is_separately_determined(trade_cost: Option<f64>) -> bool {
    trade_cost.is_some_and(|cost| cost >= SEPARATION_LOG_LIKELIHOOD)
}

/// **The standard error a variance gives a derived number**: its square root when it is a positive
/// finite number no wider than 1. Every number derived here — the share above
/// [`STATED_FREQUENCY`], `Hexp`, `E[f]` — lies in [0, 1], and an error wider than that whole range
/// says nothing about where in it the number lies; the cohort's own numbers are withheld by the same
/// rule ([`StandardError::WiderThanItsRange`](super::standard_errors::StandardError::WiderThanItsRange)). Measured
/// at two samples and 3 reads: the share above 1% came out at 0.073 ± 11.0.
fn error_from_variance(variance: f64) -> Option<f64> {
    (variance > 0.0 && variance.is_finite())
        .then(|| variance.sqrt())
        .filter(|&error| error <= DERIVED_NUMBER_RANGE)
}

/// The width of the range every derived number lies in: each is a share or a frequency in [0, 1].
const DERIVED_NUMBER_RANGE: f64 = 1.0;

#[cfg(test)]
mod tests {
    use super::*;

    /// A covariance from a dense 8 × 8 matrix, with the slots `dropped` says `None` across their
    /// row and column.
    fn covariance_of(dense: [[f64; 8]; 8], dropped: &[usize]) -> CohortCovariance {
        std::array::from_fn(|j| {
            std::array::from_fn(|k| {
                (!dropped.contains(&j) && !dropped.contains(&k)).then_some(dense[j][k])
            })
        })
    }

    /// A diagonal covariance of the given variances.
    fn diagonal(variances: [f64; 8]) -> [[f64; 8]; 8] {
        std::array::from_fn(|j| std::array::from_fn(|k| if j == k { variances[j] } else { 0.0 }))
    }

    fn duplicated(share: f64, carrier_a: f64, carrier_b: f64) -> DuplicatedPositions {
        DuplicatedPositions {
            share,
            carrier_a,
            carrier_b,
        }
    }

    /// **`V_trade` at or below 10⁻¹² of `σ`'s variance is no trade, never a negative or finite
    /// cost**: zero, a negative value, a unit in the last place of the variance above zero (the 4
    /// oracle samples' case, which cost 4.66 × 10¹⁶ before) and 10⁻¹² of it all cost `+∞`; just
    /// above that, and an ordinary value, cost `½ (ln 2)² σ² / V_trade`; not a number has no cost.
    #[test]
    fn a_trade_variance_at_rounding_level_is_no_trade() {
        let variance = 1e-3;
        for no_trade in [
            0.0,
            -1e-20,
            variance * f64::EPSILON,
            NO_TRADE_SHARE * variance,
        ] {
            assert_eq!(
                trade_cost_of(0.1, variance, no_trade),
                Some(f64::INFINITY),
                "{no_trade:e}"
            );
        }
        for trade in [(NO_TRADE_SHARE * variance).next_up(), 1e-4] {
            let cost = trade_cost_of(0.1, variance, trade).expect("a cost");
            assert!(cost.is_finite(), "{trade:e}");
            assert!(
                (cost - 0.5 * LN_2 * LN_2 * 0.01 / trade).abs() <= 1e-15 * cost,
                "{cost}"
            );
        }
        assert_eq!(trade_cost_of(0.1, variance, f64::NAN), None);
    }

    /// **Two units is the line**: a cost just under two marks the pair, two does not, and so do
    /// `+∞` and no cost at all, the other way.
    #[test]
    fn a_trade_cost_just_under_two_marks_the_pair_and_two_does_not() {
        assert!(!is_separately_determined(Some(
            SEPARATION_LOG_LIKELIHOOD.next_down()
        )));
        assert!(is_separately_determined(Some(SEPARATION_LOG_LIKELIHOOD)));
        assert!(is_separately_determined(Some(f64::INFINITY)));
        assert!(!is_separately_determined(None));
        assert!(!is_separately_determined(Some(f64::NAN)));
    }

    /// **A share that moves only with what the reads see has no trade.** Every carrier above 1%
    /// (`Beta(30, 60)`) and only the duplicated share uncertain — the carrier shapes known to
    /// 10⁻¹⁵, so that they still have errors: `S` is `σ` times a constant, so `V_trade` is zero up
    /// to rounding — the cost is `+∞` — and the pair is separately determined, its two numbers
    /// keeping their errors.
    #[test]
    fn a_share_that_moves_only_with_what_is_seen_has_no_trade() {
        let mut variances = [1e-6; 8];
        variances[cohort::DUPLICATED_SHARE] = 4e-6;
        variances[cohort::CARRIER_A] = 1e-30;
        variances[cohort::CARRIER_B] = 1e-30;
        let errors: [Option<f64>; 8] = variances.map(error_from_variance);
        let report = DuplicatedReport::of(
            &duplicated(0.01, 30.0, 60.0),
            &errors,
            &covariance_of(diagonal(variances), &[]),
        );
        assert_eq!(report.determined.trade_cost, Some(f64::INFINITY));
        assert!(report.determined.separately_determined);
        assert_error_near(report.duplicated_share, 2e-3, 1e-15, "the duplicated share");
        assert!((report.determined.share_above.value - 0.01).abs() < 1e-12);
    }

    /// **A number of the group with no covariance marks the pair, and the derived numbers carry no
    /// error**: the density's second shape dropped (as at one sample). No trade cost, the invariant
    /// share and the first shape written without their errors, `S`, `Hexp` and `E[f]` without
    /// theirs; the fixed non-reference share keeps its own, as it has one.
    #[test]
    fn a_number_without_a_covariance_marks_the_pair() {
        let errors = [Some(1e-3); 8];
        let report = DensityReport::of(
            &FrequencyDensity {
                p_invariant: 0.9,
                p_fixed_alt: 0.001,
                a: 0.3,
                b: 4.0,
            },
            &errors,
            &covariance_of(diagonal([1e-6; 8]), &[cohort::DENSITY_B]),
        );
        assert_eq!(report.determined.trade_cost, None);
        assert!(!report.determined.separately_determined);
        assert_eq!(report.invariant_share.standard_error, None);
        assert_eq!(report.segregating_shape_a.standard_error, None);
        assert_eq!(report.fixed_non_reference_share.standard_error, Some(1e-3));
        assert_eq!(report.determined.share_above.standard_error, None);
        assert_eq!(report.expected_heterozygosity.standard_error, None);
        assert_eq!(report.expected_alternative_frequency.standard_error, None);
        assert!((report.invariant_share.value - 0.9).abs() < 1e-15);
    }

    /// **A number of the group kept in the inversion but with no standard error marks the pair,
    /// with no trade cost** (spec §2.1.3, amended 2026-10-17): the duplicated class with every
    /// carrier above 1%, which with every error present has no trade (the test above), its carrier
    /// shape `b` given no error — as when it came out wider than its range. Its covariance is still
    /// there, so `S` keeps its error.
    #[test]
    fn a_number_with_a_covariance_but_no_error_marks_the_pair() {
        let mut variances = [1e-6; 8];
        variances[cohort::DUPLICATED_SHARE] = 4e-6;
        variances[cohort::CARRIER_A] = 1e-30;
        variances[cohort::CARRIER_B] = 1e-30;
        let mut errors: [Option<f64>; 8] = variances.map(error_from_variance);
        errors[cohort::CARRIER_B] = None;
        let report = DuplicatedReport::of(
            &duplicated(0.01, 30.0, 60.0),
            &errors,
            &covariance_of(diagonal(variances), &[]),
        );
        assert_eq!(report.determined.trade_cost, None);
        assert!(!report.determined.separately_determined);
        assert_eq!(report.duplicated_share.standard_error, None);
        assert_eq!(report.carrier_shape_a.standard_error, None);
        assert!(report.determined.share_above.standard_error.is_some());
    }

    /// **A derived number's error wider than 1, the whole range of each, is withheld**: an error
    /// of exactly 1 is kept, the next number up is not. Through the report: the invariant share's
    /// variance at 100 (an error of 10, as at two samples) leaves `S`, `Hexp` and `E[f]` with none,
    /// while at 10⁻⁶ they keep theirs.
    #[test]
    fn a_derived_error_wider_than_its_range_is_withheld() {
        assert_eq!(error_from_variance(1.0), Some(1.0));
        let just_over_one = 1.0_f64.next_up();
        assert_eq!(error_from_variance(just_over_one * just_over_one), None);
        assert_eq!(error_from_variance(0.25), Some(0.5));
        let density = FrequencyDensity {
            p_invariant: 0.9,
            p_fixed_alt: 0.001,
            a: 0.6,
            b: 2.2,
        };
        for (variance, kept) in [(100.0, false), (1e-6, true)] {
            let mut variances = [1e-6; 8];
            variances[cohort::P_INVARIANT] = variance;
            let report = DensityReport::of(
                &density,
                &[Some(1e-3); 8],
                &covariance_of(diagonal(variances), &[]),
            );
            for derived in [
                report.determined.share_above,
                report.expected_heterozygosity,
                report.expected_alternative_frequency,
            ] {
                assert_eq!(
                    derived.standard_error.is_some(),
                    kept,
                    "{variance}: {derived:?}"
                );
            }
        }
    }

    /// **A fixed non-reference share on its lower bound gets exact slopes**, not a central
    /// difference whose step, 10⁻⁵ of 10⁻¹², is below the last place of the segregating share. Only
    /// that share is uncertain (variance 10⁻⁶), so each derived number's error is its slope along
    /// it times 10⁻³. The slopes are `∂S/∂p_fixed_alt = −P(f > 0.01)`, `∂Hexp/∂p_fixed_alt =
    /// −2ab / ((a+b)(a+b+1))` and `∂E[f]/∂p_fixed_alt = 1 − a/(a+b)`; the code review measured the
    /// difference giving −0.0867 against −0.0533 for `S` in the first case, −0.694 against −0.893 in
    /// the second, and 0.694 against 0.786 for `E[f]` there. The exact `P(f > 0.01)` values are
    /// SciPy's `betainc` (`tmp/review_CD/py/pfa_bound.py`).
    #[test]
    fn a_fixed_non_reference_share_on_its_bound_gets_exact_slopes() {
        let cases = [
            // (p_invariant, a, b, −∂S/∂p_fixed_alt)
            (0.84, 0.02, 4.28, 0.05334440746162261),
            (0.9, 0.6, 2.2, 0.893035940028032),
            (0.99, 0.3, 14.5, 0.3988729496860478),
        ];
        for (p_invariant, a, b, share_seen) in cases {
            let mut variances = [0.0; 8];
            variances[cohort::P_FIXED_ALT] = 1e-6;
            let report = DensityReport::of(
                &FrequencyDensity {
                    p_invariant,
                    p_fixed_alt: 1e-12,
                    a,
                    b,
                },
                &[Some(1e-3); 8],
                &covariance_of(diagonal(variances), &[]),
            );
            let what = format!("at p_invariant {p_invariant}, a {a}, b {b}");
            assert_error_near(
                report.determined.share_above,
                share_seen * 1e-3,
                1e-9,
                &format!("S {what}"),
            );
            assert_error_near(
                report.expected_heterozygosity,
                2.0 * a * b / ((a + b) * (a + b + 1.0)) * 1e-3,
                1e-12,
                &format!("Hexp {what}"),
            );
            assert_error_near(
                report.expected_alternative_frequency,
                (1.0 - a / (a + b)) * 1e-3,
                1e-12,
                &format!("E[f] {what}"),
            );
        }
    }

    /// The returned fit of one of the stopping-rule report's reference runs, run to a decrement of
    /// 10⁻⁴ (`tmp/snp_stop/runs/<run>.judge.tsv`, the start-0 rows the prototype wrote after the
    /// final pass): the cohort's eight values and errors, and their covariance from the prototype's
    /// `cohort_covariance`, which is this crate's but for one thing — it left a dropped column's
    /// entries at the zero its solve gave them, where this crate has `None`, as here.
    struct ReferenceFit {
        values: [f64; 8],
        errors: [Option<f64>; 8],
        covariance: CohortCovariance,
    }

    /// One oracle accession, `s1_d1e-4`: the density's two shapes and the duplicated class dropped.
    fn one_oracle_sample() -> ReferenceFit {
        ReferenceFit {
            values: [
                0.00871559894373236,
                0.9935775412379461,
                0.0007574635009272126,
                0.31140552500589974,
                14.567269519264244,
                0.00037139638322625614,
                1.1347797046209376,
                20.3458279004514,
            ],
            errors: [
                Some(0.0009576881),
                Some(0.0003897345),
                Some(2.080833e-05),
                None,
                None,
                None,
                None,
                None,
            ],
            covariance: {
                let mut covariance = [[None; 8]; 8];
                let kept = [
                    [9.1716656745e-07, 3.8109600763e-08, 1.0362793683e-10],
                    [3.8109600763e-08, 1.5189295242e-07, 3.2137879714e-10],
                    [1.0362793683e-10, 3.2137879714e-10, 4.3298650849e-10],
                ];
                for (j, row) in kept.iter().enumerate() {
                    for (k, &entry) in row.iter().enumerate() {
                        covariance[j][k] = Some(entry);
                    }
                }
                covariance
            },
        }
    }

    /// The 4 oracle accessions, `o4_d1e-4`: nothing dropped; the density's second shape and the
    /// duplicated share have errors wider than their ranges, so no error but a covariance.
    fn four_oracle_samples() -> ReferenceFit {
        ReferenceFit {
            values: [
                0.005154554164246633,
                0.9983384064044813,
                5.8896057925424846e-05,
                15.439289432714222,
                50.0,
                0.013922414984305343,
                0.02,
                2.0735990491108494,
            ],
            errors: [
                Some(0.0003122813),
                Some(0.0003361555),
                Some(1.326464e-05),
                Some(36.77845),
                None,
                None,
                Some(0.2772658),
                Some(0.5503759),
            ],
            covariance: covariance_of(
                [
                    [
                        9.7519580308e-08,
                        -5.239418544e-10,
                        1.4308185177e-11,
                        2.5638992775e-06,
                        9.8134672785e-06,
                        2.6437587412e-07,
                        -4.4691148187e-07,
                        -2.8282917069e-06,
                    ],
                    [
                        -5.239418544e-10,
                        1.1300051417e-07,
                        -3.089206791e-09,
                        0.0041043751841,
                        0.012802170841,
                        -5.7476780831e-05,
                        8.498849815e-05,
                        0.00015871775705,
                    ],
                    [
                        1.4308185177e-11,
                        -3.089206791e-09,
                        1.7595068956e-10,
                        0.00014361052418,
                        0.00045070474433,
                        2.1106893454e-06,
                        -3.1209829894e-06,
                        -5.8484110569e-06,
                    ],
                    [
                        2.5638992775e-06,
                        0.0041043751841,
                        0.00014361052418,
                        1352.6544786,
                        4235.8377362,
                        0.33302831538,
                        -0.4958177385,
                        -1.202247844,
                    ],
                    [
                        9.8134672785e-06,
                        0.012802170841,
                        0.00045070474433,
                        4235.8377362,
                        13269.853122,
                        1.0196849288,
                        -1.518943799,
                        -3.7370165893,
                    ],
                    [
                        2.6437587412e-07,
                        -5.7476780831e-05,
                        2.1106893454e-06,
                        0.33302831538,
                        1.0196849288,
                        0.035167304994,
                        -0.051993944706,
                        -0.097047839975,
                    ],
                    [
                        -4.4691148187e-07,
                        8.498849815e-05,
                        -3.1209829894e-06,
                        -0.4958177385,
                        -1.518943799,
                        -0.051993944706,
                        0.076876313341,
                        0.14380761028,
                    ],
                    [
                        -2.8282917069e-06,
                        0.00015871775705,
                        -5.8484110569e-06,
                        -1.202247844,
                        -3.7370165893,
                        -0.097047839975,
                        0.14380761028,
                        0.30291363416,
                    ],
                ],
                &[],
            ),
        }
    }

    /// The 59 tomato2 accessions, `t59_d1e-4`: nothing dropped, both first shapes on their bound.
    fn fifty_nine_samples() -> ReferenceFit {
        ReferenceFit {
            values: [
                0.010591692025235237,
                0.84161686120661,
                2.9594621109687785e-05,
                0.02,
                4.280399457742142,
                0.015427236567361826,
                0.02,
                1.4076589942132522,
            ],
            errors: [
                Some(0.0001401734),
                Some(0.06779098),
                Some(4.05028e-06),
                Some(0.008922428),
                Some(0.06050911),
                Some(0.02293727),
                Some(0.03094622),
                Some(0.04663929),
            ],
            covariance: covariance_of(
                [
                    [
                        1.9648578244e-08,
                        -5.6675895024e-08,
                        -9.0329319521e-12,
                        -1.3688453354e-08,
                        -1.8206569827e-07,
                        -2.7789628197e-07,
                        3.8043001191e-07,
                        3.1789088553e-07,
                    ],
                    [
                        -5.6675895024e-08,
                        0.0045956174407,
                        4.6838585419e-09,
                        0.00060392963452,
                        0.0025189466003,
                        0.0001068586123,
                        -0.00014212305617,
                        -7.9009864329e-05,
                    ],
                    [
                        -9.0329319521e-12,
                        4.6838585419e-09,
                        1.64047648e-11,
                        9.3692368395e-10,
                        4.7107988058e-09,
                        1.5940702116e-08,
                        -2.1491348499e-08,
                        -1.9915024976e-08,
                    ],
                    [
                        -1.3688453354e-08,
                        0.00060392963452,
                        9.3692368395e-10,
                        7.9609718746e-05,
                        0.00033503096978,
                        2.4484372855e-05,
                        -3.2773974333e-05,
                        -2.3925085456e-05,
                    ],
                    [
                        -1.8206569827e-07,
                        0.0025189466003,
                        4.7107988058e-09,
                        0.00033503096978,
                        0.0036613520292,
                        3.7159880624e-05,
                        -5.0398520573e-05,
                        -4.8033441296e-05,
                    ],
                    [
                        -2.7789628197e-07,
                        0.0001068586123,
                        1.5940702116e-08,
                        2.4484372855e-05,
                        3.7159880624e-05,
                        0.0005261181728,
                        -0.00070968050205,
                        -0.00066762356601,
                    ],
                    [
                        3.8043001191e-07,
                        -0.00014212305617,
                        -2.1491348499e-08,
                        -3.2773974333e-05,
                        -5.0398520573e-05,
                        -0.00070968050205,
                        0.00095766857101,
                        0.00091047058239,
                    ],
                    [
                        3.1789088553e-07,
                        -7.9009864329e-05,
                        -1.9915024976e-08,
                        -2.3925085456e-05,
                        -4.8033441296e-05,
                        -0.00066762356601,
                        0.00091047058239,
                        0.0021752232511,
                    ],
                ],
                &[],
            ),
        }
    }

    impl ReferenceFit {
        fn reports(&self) -> (DensityReport, DuplicatedReport) {
            let v = self.values;
            (
                DensityReport::of(
                    &FrequencyDensity {
                        p_invariant: v[1],
                        p_fixed_alt: v[2],
                        a: v[3],
                        b: v[4],
                    },
                    &self.errors,
                    &self.covariance,
                ),
                DuplicatedReport::of(
                    &duplicated(v[5], v[6], v[7]),
                    &self.errors,
                    &self.covariance,
                ),
            )
        }
    }

    /// `got` within `relative` of `expected`.
    fn assert_near(got: f64, expected: f64, relative: f64, what: &str) {
        assert!(
            (got - expected).abs() <= relative * expected.abs(),
            "{what}: {got:e} against {expected:e}"
        );
    }

    /// `number`'s error within `relative` of `expected`.
    fn assert_error_near(number: ReportedNumber, expected: f64, relative: f64, what: &str) {
        assert_near(number.standard_error.expect(what), expected, relative, what);
    }

    /// **The measured verdicts, reproduced from the reference runs' recorded values and
    /// covariances** (spec §2.1.3's table). The expected numbers are computed the way this module
    /// does — central differences at a step of 10⁻⁵ of each value, SciPy's `betainc` — by
    /// `tmp/impl_CD/verdicts_reference.py`, and agree with `tmp/snp_plan/dup_trade.py`'s
    /// forward differences to the figures the spec quotes. At one accession both pairs are marked,
    /// with no trade cost and no error on any derived number.
    #[test]
    fn at_one_oracle_sample_both_pairs_are_marked() {
        let (density, duplicated) = one_oracle_sample().reports();
        for determined in [density.determined, duplicated.determined] {
            assert!(!determined.separately_determined);
            assert_eq!(determined.trade_cost, None);
            assert_eq!(determined.share_above.standard_error, None);
        }
        assert_near(
            density.determined.share_above.value,
            0.002326684888570621,
            1e-12,
            "the share segregating above 1%",
        );
        assert_near(
            duplicated.determined.share_above.value,
            0.00031947552675612575,
            1e-12,
            "the share duplicated in more than 1% of the panel",
        );
        assert_eq!(density.invariant_share.standard_error, None);
        assert_eq!(density.segregating_shape_a.standard_error, None);
        assert_eq!(
            density.fixed_non_reference_share.standard_error,
            Some(2.080833e-05)
        );
        assert_eq!(density.expected_heterozygosity.standard_error, None);
        assert_eq!(density.expected_alternative_frequency.standard_error, None);
        assert_eq!(duplicated.duplicated_share.standard_error, None);
        assert_eq!(duplicated.carrier_shape_a.standard_error, None);
    }

    /// **At 4 oracle accessions both pairs are marked, each because a number of its group has no
    /// standard error** (spec §2.1.3, amended 2026-10-17): the density's second shape `b` (on its
    /// bound of 50, its error 115 against a range of 49.98) and the duplicated share. Neither pair
    /// has a trade cost; the derived numbers keep their errors, each well inside [0, 1].
    ///
    /// **The trades themselves are as measured before the amendment**: given an error for the
    /// number that has none, the density has no trade — every segregating position is above 1%
    /// (`a` = 15.4), so `σ` is the determined share and `V_trade` comes out at rounding level, the
    /// cost `+∞` — and the duplicated pair costs 0.0055 units (`carrier_a` on its bound).
    #[test]
    fn at_four_oracle_samples_both_pairs_are_marked_as_a_number_of_each_has_no_error() {
        let reference = four_oracle_samples();
        let (density, duplicated) = reference.reports();
        assert!(!density.determined.separately_determined);
        assert_eq!(density.determined.trade_cost, None);
        assert_eq!(density.invariant_share.standard_error, None);
        assert_eq!(density.segregating_shape_a.standard_error, None);
        assert_eq!(density.segregating_shape_b.standard_error, None);
        assert_eq!(
            density.fixed_non_reference_share.standard_error,
            Some(1.326464e-05)
        );
        assert_near(
            density.determined.share_above.value,
            0.0016026975375932264,
            1e-12,
            "S",
        );
        assert_error_near(
            density.determined.share_above,
            0.00032710556595231063,
            1e-6,
            "S",
        );
        assert_error_near(
            density.expected_heterozygosity,
            0.00010859828950885404,
            1e-6,
            "Hexp",
        );
        assert_error_near(
            density.expected_alternative_frequency,
            8.362336563428386e-05,
            1e-6,
            "E[f]",
        );

        assert!(!duplicated.determined.separately_determined);
        assert_eq!(duplicated.determined.trade_cost, None);
        assert_near(
            duplicated.determined.share_above.value,
            0.0009618499283254206,
            1e-12,
            "the duplicated S",
        );
        assert_error_near(
            duplicated.determined.share_above,
            0.0001594938578793182,
            1e-6,
            "the duplicated S",
        );
        assert_eq!(duplicated.duplicated_share.standard_error, None);
        assert_eq!(duplicated.carrier_shape_a.standard_error, None);
        assert_eq!(duplicated.carrier_shape_b.standard_error, Some(0.5503759));

        let mut every_error = reference.errors;
        every_error[cohort::DENSITY_B] = Some(1.0);
        every_error[cohort::DUPLICATED_SHARE] = Some(1e-3);
        let (density, duplicated) = ReferenceFit {
            errors: every_error,
            ..reference
        }
        .reports();
        assert_eq!(density.determined.trade_cost, Some(f64::INFINITY));
        assert!(density.determined.separately_determined);
        assert_near(
            duplicated.determined.trade_cost.expect("a cost"),
            0.005494469013440156,
            1e-4,
            "the duplicated pair's trade cost",
        );
        assert!(!duplicated.determined.separately_determined);
    }

    /// **At 59 tomato2 accessions both pairs are marked**: the density's trade costs 1.35 units, the
    /// duplicated class's 0.116. The numbers are those of the spec's example tables: `S` 8.447 × 10⁻³
    /// ± 2.04 × 10⁻⁴, `Hexp` 1.1895 × 10⁻³ ± 3.03 × 10⁻⁵, `E[f]` 7.661 × 10⁻⁴ ± 2.02 × 10⁻⁵, the
    /// duplicated `S` 1.212 × 10⁻³ ± 3.4 × 10⁻⁵; the marked numbers keep their values and lose their
    /// errors, the others keep theirs.
    #[test]
    fn at_fifty_nine_samples_both_pairs_are_marked() {
        let (density, duplicated) = fifty_nine_samples().reports();
        assert!(!density.determined.separately_determined);
        assert_near(
            density.determined.trade_cost.expect("a cost"),
            1.346728287322403,
            1e-4,
            "the density's trade cost",
        );
        assert_near(
            density.determined.share_above.value,
            0.008446973378575733,
            1e-12,
            "S",
        );
        assert_error_near(
            density.determined.share_above,
            0.00020398488887385483,
            1e-6,
            "S",
        );
        assert_near(
            density.expected_heterozygosity.value,
            0.0011894732518292063,
            1e-12,
            "Hexp",
        );
        assert_error_near(
            density.expected_heterozygosity,
            3.0347029828571595e-05,
            1e-6,
            "Hexp",
        );
        assert_near(
            density.expected_alternative_frequency.value,
            0.000766054318532408,
            1e-12,
            "E[f]",
        );
        assert_error_near(
            density.expected_alternative_frequency,
            2.021977341875278e-05,
            1e-6,
            "E[f]",
        );
        assert_eq!(density.invariant_share.value, 0.84161686120661);
        assert_eq!(density.invariant_share.standard_error, None);
        assert_eq!(density.segregating_shape_a.value, 0.02);
        assert_eq!(density.segregating_shape_a.standard_error, None);
        assert_eq!(density.segregating_shape_b.standard_error, Some(0.06050911));
        assert_eq!(
            density.fixed_non_reference_share.standard_error,
            Some(4.05028e-06)
        );

        assert!(!duplicated.determined.separately_determined);
        assert_near(
            duplicated.determined.trade_cost.expect("a cost"),
            0.11570630999060695,
            1e-4,
            "the duplicated pair's trade cost",
        );
        assert_near(
            duplicated.determined.share_above.value,
            0.0012121951526186507,
            1e-12,
            "the duplicated S",
        );
        assert_error_near(
            duplicated.determined.share_above,
            3.378019850602306e-05,
            1e-6,
            "the duplicated S",
        );
        assert_eq!(duplicated.duplicated_share.standard_error, None);
        assert_eq!(duplicated.carrier_shape_a.standard_error, None);
        assert_eq!(duplicated.carrier_shape_b.standard_error, Some(0.04663929));
    }
}
