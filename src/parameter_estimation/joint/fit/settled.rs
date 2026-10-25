//! **When the SNP/indel fit stops**: once the step still left to the likelihood's maximum, as a
//! Newton step estimates it, is worth less than a thousandth of a log-likelihood unit.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §2 (as amended at plan step B1, and again on
//! 2026-10-16) and §3.3. Build order: `doc/devel/ng/impl_plan/snp_fit_convergence.md`, steps A1–A2.
//!
//! # The rule
//!
//! Once a cycle of the fit has gained less than
//! [`JointFitConfig::log_likelihood_stillness`](super::JointFitConfig::log_likelihood_stillness) of
//! log-likelihood a position, every later cycle's first pass — a *judging pass* — also sums the
//! information and each parameter's score, the slope of the whole log-likelihood. Solved against
//! each other they give the Newton step to the maximum within the parameters' bounds
//! ([`newton_step`](super::standard_errors::newton_step)): how far each parameter still is from where
//! the likelihood peaks, **however slowly the fit happens to be moving towards it**, and the
//! log-likelihood that step would gain under the quadratic model — the *Newton decrement*.
//!
//! **A start has converged when the decrement is at least zero and below [`DECREMENT_THRESHOLD`]**
//! ([`has_converged`]), and stops at once at the parameters that pass judged; the final pass follows
//! there. One number bounds every parameter that has a standard error and every quantity computed
//! from those alone: the step moves any of them by at most `√(2λ)` of its own standard error; at
//! 0.01 that is 0.14 errors. The parameters with none — not identified, or with an error wider
//! than their range — are profiled out of the decrement (plan step B4,
//! [`newton_step`](super::standard_errors::newton_step)), and the number bounds nothing about them.
//! On every cohort measured (1, 4 and 59 tomato samples) every determined quantity stopped within a
//! tenth of its error of a fit run to 10⁻⁴; the decrement understates what is left by two and a half
//! to six times near the maximum, so the threshold was chosen from that measurement and not from the
//! bound. A decrement that is not a number is not convergence: it means the bounded solve reached
//! its round limit short of the quadratic model's maximum. Nor is a negative one, which the solve
//! cannot give (`fit_precision.md` §2, amendment of 2026-10-17) and the test refuses all the same.
//!
//! **The per-parameter test survives for reporting only.** A parameter is *settled* when its distance
//! is below [`SETTLED_FRACTION`] of its standard error; a parameter without an error or a distance is
//! settled by definition. The progress line, and a start that did not converge, say how many are not
//! settled and which is furthest ([`not_settled`], [`furthest`]). It stopped the fit until
//! 2026-10-16, and stops nothing now: with 185 parameters each a tenth of an error from the maximum,
//! the remaining step can still be worth almost a unit.
//!
//! # Agreement between starts
//!
//! Spec §3.4, plan step B2. The fit runs its starting points one after another. Once one has
//! converged, **a later start stops as soon as it is heading where that one arrived**: when, on one of
//! its judging passes, its *projected endpoint* — each parameter's value plus its Newton distance — lies
//! within [`AGREEMENT_FRACTION`] of the earlier answer's standard error of that answer, on every
//! parameter the earlier answer gives an error ([`agrees`]). It is recorded as having agreed
//! ([`StartEnding::Agreed`]) and does not compete to be the fit's answer. A start whose judging pass
//! finds it converged is not stopped by agreement: it converges and competes.
//!
//! The earlier answer is the best-scoring start that has converged so far ([`EarlierAnswer`]); a start
//! stopped at the pass limit is never one, since its errors need not be the maximum's.

use super::standard_errors::StandardError;

/// The share of its own standard error a parameter's distance to the maximum must fall below for the
/// parameter to count as settled — the spec's starting value, kept at checkpoint A (owner,
/// 2026-09-28).
pub(in crate::parameter_estimation::joint) const SETTLED_FRACTION: f64 = 0.1;

/// **The log-likelihood the remaining Newton step must be worth less than for a start to have
/// converged** — the Newton decrement's threshold: 0.001 from 2026-10-16, 0.01 since 2026-10-25
/// (owner), from the measurements in `fit_precision.md` §2. Since the starts explore on a thinned
/// census and the winner finishes on all of it, 0.01 stops every number the file gives an error
/// within 0.055 of that error of a fit run to 10⁻⁴ (1, 8 and 62 tomato samples), for 24 to 30% less
/// of the fit's arithmetic than 0.001.
pub(in crate::parameter_estimation::joint) const DECREMENT_THRESHOLD: f64 = 0.01;

/// **Whether a start whose remaining Newton step is worth `decrement` has converged**, at
/// `threshold`: when the decrement is at least zero and below the threshold. A decrement that is not
/// a number comes from a box solve that reached its round limit, and is not convergence. Nor is a
/// negative one: the solve climbs the quadratic model from staying put and cannot give one, so a
/// negative value would mean a defect in it, never a start that has arrived (`fit_precision.md` §2,
/// traps).
pub(super) fn has_converged(decrement: f64, threshold: f64) -> bool {
    decrement >= 0.0 && decrement < threshold
}

/// **How many parameters are not yet settled**: those whose distance to the maximum (`distances`, a
/// Newton step's) is not below `fraction` of their standard error (`errors`), both in the order the
/// fit's vector lists them. A parameter missing either is settled.
pub(super) fn not_settled(
    errors: &[StandardError],
    distances: &[Option<f64>],
    fraction: f64,
) -> usize {
    // Both come from one pass's sums laid out by the same parameters (`by_coordinate` and
    // `newton_step`), so their lengths agree unless one was laid out for another fit.
    assert_eq!(
        errors.len(),
        distances.len(),
        "one error and one distance a parameter"
    );
    errors
        .iter()
        .zip(distances)
        .filter(|(error, distance)| match (error.value(), distance) {
            // Not below the fraction, a distance that is not a number included.
            (Some(error), Some(distance)) => !matches!(
                distance.abs().partial_cmp(&(fraction * error)),
                Some(std::cmp::Ordering::Less)
            ),
            _ => false,
        })
        .count()
}

/// **The parameter furthest from the maximum in units of its own error**: its position in the order
/// the fit's vector lists them and its distance over its error, among those with both; `None` when no
/// parameter has both. What a fit that has not settled prints to say which parameter holds it.
pub(super) fn furthest(
    errors: &[StandardError],
    distances: &[Option<f64>],
) -> Option<(usize, f64)> {
    errors
        .iter()
        .zip(distances)
        .enumerate()
        .filter_map(|(j, (error, distance))| Some((j, distance.as_ref()?.abs() / error.value()?)))
        .max_by(|(_, left), (_, right)| left.total_cmp(right))
}

/// **How close a later start's projected endpoint must come to an earlier start's answer**, in units
/// of that answer's standard errors, for the later start to stop as having agreed (spec §3.4) — the
/// spec's starting value.
pub(super) const AGREEMENT_FRACTION: f64 = 0.5;

/// **How one start of the SNP/indel fit ended** (spec §3.5).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum StartEnding {
    /// The remaining Newton step was worth less than the decrement threshold ([`has_converged`]).
    Converged,
    /// It ran out of passes first.
    AtTheLimit,
    /// It was heading where an earlier start had converged ([`agrees`]), and stopped there: it adds
    /// nothing new, and does not compete to be the fit's answer. `with_start` counts from one.
    Agreed { with_start: usize },
}

/// **Why a start that did not converge had not** (spec §3.5): at the values it returned, how many
/// parameters were not yet within the settled fraction of their standard error of the maximum, of
/// how many, and the one furthest from it, named as the run's log names it, with its distance in
/// units of its own standard error.
#[derive(Clone, Debug, PartialEq)]
pub struct FurthestFromSettled {
    /// How many parameters were not yet within the settled fraction of their error.
    pub not_settled: usize,
    /// Of how many parameters the fit carries.
    pub parameters: usize,
    /// The furthest, as the run's log names it: `invariant share`, `homozygote excess of LA1589`,
    /// `error rate at ordinary positions of LA1589's library 2`.
    pub parameter: String,
    /// Its distance to the maximum over its standard error.
    pub errors: f64,
}

/// **An earlier start's answer, as a later start is judged against it**: which start it was, counted
/// from one, its log-likelihood, and its values and standard errors at the parameters it returned, in
/// the order the fit's vector lists them.
pub(super) struct EarlierAnswer {
    pub(super) start: usize,
    pub(super) log_likelihood: f64,
    pub(super) values: Vec<f64>,
    pub(super) errors: Vec<StandardError>,
}

impl EarlierAnswer {
    /// The answer start `start` (counted from one) returned: its parameters, and the errors from the
    /// information its final pass summed.
    pub(super) fn of(
        start: usize,
        outcome: &super::StartOutcome,
        group_index: &[Vec<usize>],
    ) -> Self {
        Self {
            start,
            log_likelihood: outcome.statistics.log_likelihood,
            values: outcome.parameters.values(),
            errors: super::StandardErrors::of(&outcome.information)
                .by_coordinate(&outcome.parameters, group_index),
        }
    }
}

/// **Whether a start is heading where an earlier one arrived**: its projected endpoint — each
/// parameter's value plus its distance to the maximum (`distances`, a Newton step's; a parameter
/// without one stays where it is) — within `fraction` of the earlier answer's standard error of that
/// answer, on every parameter the earlier answer gives an error. The endpoint, not the value: early
/// in a start the value is far from anything (spec §3.4).
///
/// **An earlier answer that gives no parameter an error is agreed with by nothing**: there is no
/// yardstick to be within. At a `fraction` of zero, below it or not a number, no start agrees.
pub(super) fn agrees(
    values: &[f64],
    distances: &[Option<f64>],
    earlier: &EarlierAnswer,
    fraction: f64,
) -> bool {
    // Every list is the fit's vector in its own order, from one fit: the same parameters.
    assert!(
        distances.len() == values.len()
            && earlier.values.len() == values.len()
            && earlier.errors.len() == values.len(),
        "two starts of one fit list the same parameters"
    );
    let mut judged = 0;
    let within = values
        .iter()
        .zip(distances)
        .zip(earlier.values.iter().zip(&earlier.errors))
        .all(|((value, distance), (answer, error))| {
            let Some(error) = error.value() else {
                return true;
            };
            judged += 1;
            let endpoint = value + distance.unwrap_or(0.0);
            (endpoint - answer).abs() < fraction * error
        });
    within && judged > 0
}

/// **What a start's log line adds about how far from settled it stopped** (spec §3.5): nothing for a
/// start that converged; for one that did not, at the values it returned, how many parameters were
/// still further than `fraction` of a standard error from where the likelihood peaks and the furthest
/// — or that every one was within it, for a start that ran out of passes before a cycle judged it —
/// or, with no parameter to judge, that too. Begins with `"; "` when not empty.
pub(super) fn describe_short_of_settled(
    ended: StartEnding,
    furthest: Option<&FurthestFromSettled>,
    fraction: f64,
) -> String {
    if ended == StartEnding::Converged {
        return String::new();
    }
    match furthest {
        None => "; no parameter has both a standard error and a distance to the maximum, so how \
                 far it stopped from settled cannot be said"
            .to_owned(),
        Some(furthest) if furthest.not_settled == 0 => format!(
            "; at the values it returned every one of its {} parameter(s) was within {fraction} \
             standard errors of where the likelihood peaks (furthest, the {}, at {:.2})",
            furthest.parameters, furthest.parameter, furthest.errors
        ),
        Some(furthest) => format!(
            "; at the values it returned {} of {} parameter(s) still more than {fraction} standard \
             errors from where the likelihood peaks; furthest, the {}, at {:.2} standard errors",
            furthest.not_settled, furthest.parameters, furthest.parameter, furthest.errors
        ),
    }
}

/// **Whether a start that ended so, scoring `score`, becomes the answer later starts are judged
/// against**, in place of one scoring `current`: only a start that converged — an answer found at the
/// pass limit has errors that need not be the maximum's (spec §3.4) — and only above the current one,
/// so later starts are judged against the best converged answer.
pub(super) fn is_the_new_yardstick(ended: StartEnding, score: f64, current: Option<f64>) -> bool {
    ended == StartEnding::Converged && current.is_none_or(|current| score > current)
}

/// **Whether a start that ended so, scoring `score`, becomes the fit's answer** in place of one
/// scoring `current`: the best log-likelihood wins, among the starts that did not agree. A start that
/// agreed stopped part-way on purpose, heading where another had already arrived.
pub(super) fn is_the_new_best(ended: StartEnding, score: f64, current: Option<f64>) -> bool {
    !matches!(ended, StartEnding::Agreed { .. }) && current.is_none_or(|current| score > current)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The verdict at the decrements `fit_precision.md` §3.6 item 4 names: −2,959 (measured on the
    /// subsampled oracle, a bounded solve ended at its round limit) and NaN are not convergence; 0 and
    /// just below the threshold are; the threshold itself is not.
    #[test]
    fn a_start_has_converged_below_the_threshold_and_never_at_a_negative_or_nan_decrement() {
        assert!(!has_converged(-2_959.0, DECREMENT_THRESHOLD));
        assert!(!has_converged(f64::NAN, DECREMENT_THRESHOLD));
        assert!(has_converged(0.0, DECREMENT_THRESHOLD));
        assert!(has_converged(
            DECREMENT_THRESHOLD.next_down(),
            DECREMENT_THRESHOLD
        ));
        assert!(!has_converged(DECREMENT_THRESHOLD, DECREMENT_THRESHOLD));
        assert_eq!(DECREMENT_THRESHOLD, 0.01);
    }

    #[test]
    fn a_parameter_is_settled_below_its_share_of_its_error() {
        let errors = [
            StandardError::Estimated(1.0),
            StandardError::Estimated(0.01),
        ];
        // 0.09 of the first's error and 0.5 of the second's: one not settled.
        assert_eq!(
            not_settled(&errors, &[Some(0.09), Some(-0.005)], SETTLED_FRACTION),
            1
        );
        // A step downwards counts by its size.
        assert_eq!(
            not_settled(&errors, &[Some(-0.2), Some(0.0005)], SETTLED_FRACTION),
            1
        );
        assert_eq!(
            not_settled(&errors, &[Some(0.0), Some(-0.0009)], SETTLED_FRACTION),
            0
        );
        // Exactly at the fraction is not below it.
        assert_eq!(not_settled(&errors, &[Some(0.25), Some(0.0)], 0.25), 1);
    }

    #[test]
    fn a_parameter_without_an_error_or_a_distance_is_settled() {
        let errors = [
            StandardError::NoInformation,
            StandardError::HeldFixed,
            StandardError::NotIdentified,
            StandardError::WiderThanItsRange(80.0),
            StandardError::Estimated(1.0),
        ];
        assert_eq!(
            not_settled(
                &errors,
                &[Some(9.0), Some(9.0), Some(9.0), Some(9.0), None],
                SETTLED_FRACTION
            ),
            0
        );
    }

    #[test]
    fn at_a_fraction_of_zero_nothing_with_an_error_and_a_distance_settles() {
        let errors = [StandardError::Estimated(1.0), StandardError::Estimated(1.0)];
        assert_eq!(not_settled(&errors, &[Some(0.0), Some(1e-300)], 0.0), 2);
    }

    #[test]
    fn a_distance_that_is_not_a_number_is_not_settled() {
        let errors = [StandardError::Estimated(1.0)];
        assert_eq!(not_settled(&errors, &[Some(f64::NAN)], SETTLED_FRACTION), 1);
    }

    #[test]
    fn the_furthest_parameter_is_the_largest_distance_in_errors() {
        let errors = [
            StandardError::Estimated(1.0),
            StandardError::Estimated(0.25),
            StandardError::NoInformation,
            StandardError::Estimated(2.0),
        ];
        // 0.5, 0.25 and 0.375 errors; the third has no error and is not a candidate however far.
        assert_eq!(
            furthest(&errors, &[Some(-0.5), Some(0.0625), Some(90.0), Some(0.75)]),
            Some((0, 0.5))
        );
        assert_eq!(
            furthest(&errors, &[None, Some(-0.1875), None, Some(1.0)]),
            Some((1, 0.75))
        );
        assert_eq!(furthest(&errors[2..3], &[Some(1.0)]), None);
    }

    fn an_earlier_answer() -> EarlierAnswer {
        EarlierAnswer {
            start: 1,
            log_likelihood: -100.0,
            values: vec![1.0, 10.0, 0.5],
            errors: vec![
                StandardError::Estimated(0.1),
                StandardError::Estimated(2.0),
                StandardError::NotIdentified,
            ],
        }
    }

    #[test]
    fn a_start_agrees_when_its_endpoint_is_within_half_an_error_on_every_parameter() {
        let earlier = an_earlier_answer();
        // Endpoints 1.04 and 10.5: 0.4 and 0.25 errors away. The third has no error in the earlier
        // answer and is not judged, however far.
        assert!(agrees(
            &[1.0, 11.0, 0.9],
            &[Some(0.04), Some(-0.5), Some(5.0)],
            &earlier,
            AGREEMENT_FRACTION
        ));
        // 1.06 is 0.6 errors away.
        assert!(!agrees(
            &[1.0, 11.0, 0.9],
            &[Some(0.06), Some(-0.5), None],
            &earlier,
            AGREEMENT_FRACTION
        ));
    }

    #[test]
    fn the_endpoint_and_not_the_value_is_compared() {
        let earlier = an_earlier_answer();
        // The value is 5 errors off, the endpoint on the answer.
        assert!(agrees(
            &[1.5, 10.0, 0.5],
            &[Some(-0.5), None, None],
            &earlier,
            AGREEMENT_FRACTION
        ));
        // The value is on the answer, the endpoint 5 errors off.
        assert!(!agrees(
            &[1.0, 10.0, 0.5],
            &[Some(0.5), None, None],
            &earlier,
            AGREEMENT_FRACTION
        ));
        // Without a distance a parameter's endpoint is its value.
        assert!(agrees(
            &[1.0, 10.0, 0.5],
            &[None, None, None],
            &earlier,
            AGREEMENT_FRACTION
        ));
    }

    /// With no parameter given an error by the earlier answer there is nothing to agree on, at any
    /// fraction; and an endpoint that is not a number agrees with nothing.
    #[test]
    fn nothing_agrees_with_an_answer_that_gives_no_error_or_at_a_nan_endpoint() {
        let mut earlier = an_earlier_answer();
        earlier.errors = vec![StandardError::NotIdentified; 3];
        assert!(!agrees(
            &[1.0, 10.0, 0.5],
            &[None, None, None],
            &earlier,
            f64::INFINITY
        ));
        let earlier = an_earlier_answer();
        assert!(!agrees(
            &[1.0, 10.0, 0.5],
            &[Some(f64::NAN), None, None],
            &earlier,
            AGREEMENT_FRACTION
        ));
    }

    #[test]
    #[should_panic(expected = "list the same parameters")]
    fn a_distance_list_of_another_length_is_refused() {
        agrees(
            &[1.0, 10.0, 0.5],
            &[None, None],
            &an_earlier_answer(),
            AGREEMENT_FRACTION,
        );
    }

    /// Only a converged start becomes the yardstick, and only above the current one — the best
    /// converged answer, never one at the pass limit, however well it scores.
    #[test]
    fn the_yardstick_is_the_best_converged_start() {
        assert!(is_the_new_yardstick(StartEnding::Converged, -10.0, None));
        assert!(is_the_new_yardstick(
            StartEnding::Converged,
            -10.0,
            Some(-11.0)
        ));
        assert!(!is_the_new_yardstick(
            StartEnding::Converged,
            -12.0,
            Some(-11.0)
        ));
        assert!(!is_the_new_yardstick(StartEnding::AtTheLimit, -1.0, None));
        assert!(!is_the_new_yardstick(
            StartEnding::AtTheLimit,
            -1.0,
            Some(-11.0)
        ));
        assert!(!is_the_new_yardstick(
            StartEnding::Agreed { with_start: 1 },
            -1.0,
            None
        ));
    }

    /// The best log-likelihood wins among the starts that did not agree; one that agreed does not
    /// compete however well it scores.
    #[test]
    fn the_answer_is_the_best_start_that_did_not_agree() {
        assert!(is_the_new_best(StartEnding::AtTheLimit, -1.0, Some(-10.0)));
        assert!(is_the_new_best(StartEnding::Converged, -1.0, None));
        assert!(!is_the_new_best(StartEnding::Converged, -12.0, Some(-10.0)));
        assert!(!is_the_new_best(
            StartEnding::Agreed { with_start: 1 },
            -1.0,
            Some(-10.0)
        ));
        assert!(!is_the_new_best(
            StartEnding::Agreed { with_start: 1 },
            -1.0,
            None
        ));
    }

    /// The start's log line says nothing more for a start that converged, and for one that did not,
    /// the count, the name and the distance — or that every parameter was within the fraction, or
    /// that none could be judged.
    #[test]
    fn a_start_that_did_not_converge_says_how_far_it_stopped_from_settled() {
        let furthest = FurthestFromSettled {
            not_settled: 4,
            parameters: 17,
            parameter: "allele-frequency shape a".to_owned(),
            errors: 0.2344,
        };
        assert_eq!(
            describe_short_of_settled(StartEnding::Converged, Some(&furthest), 0.1),
            ""
        );
        let at_the_limit = describe_short_of_settled(StartEnding::AtTheLimit, Some(&furthest), 0.1);
        for part in [
            "4 of 17 parameter(s) still more than 0.1 standard errors",
            "the allele-frequency shape a, at 0.23 standard errors",
        ] {
            assert!(at_the_limit.contains(part), "{at_the_limit}");
        }
        let agreed =
            describe_short_of_settled(StartEnding::Agreed { with_start: 1 }, Some(&furthest), 0.1);
        assert_eq!(agreed, at_the_limit);
        let settled = FurthestFromSettled {
            not_settled: 0,
            errors: 0.04,
            ..furthest
        };
        let never_judged = describe_short_of_settled(StartEnding::AtTheLimit, Some(&settled), 0.1);
        assert!(
            never_judged.contains("every one of its 17 parameter(s) was within 0.1"),
            "{never_judged}"
        );
        let nothing = describe_short_of_settled(StartEnding::AtTheLimit, None, 0.1);
        assert!(nothing.contains("cannot be said"), "{nothing}");
    }

    #[test]
    fn at_a_fraction_of_zero_no_start_agrees() {
        let earlier = an_earlier_answer();
        assert!(!agrees(
            &[1.0, 10.0, 0.5],
            &[None, None, None],
            &earlier,
            0.0
        ));
    }
}
