//! **When the SNP/indel fit stops**: once every parameter is within a tenth of its own standard error
//! of the likelihood's maximum, as a Newton step estimates the distance.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §2 (as amended at plan step B1) and §3.3. Build
//! order: `doc/devel/implementation_plans/fit_precision.md`, step B1.
//!
//! # The rule
//!
//! Once a cycle of the fit has gained less than
//! [`JointFitConfig::log_likelihood_stillness`](super::JointFitConfig::log_likelihood_stillness) of
//! log-likelihood a position, every later cycle's first pass also sums the information and each
//! parameter's score — the slope of the whole log-likelihood. Solved against each other they give
//! the Newton step to the maximum ([`newton_step`](super::standard_errors::newton_step)): how far each
//! parameter still is from where the likelihood peaks, **however slowly the fit happens to be moving
//! towards it**. That is what the first version of this rule lacked: it projected the distance from
//! the fit's last two moves, and on the four-accession oracle cohort it stopped every start 74
//! log-likelihood units short of a fit run to 1,000 passes, on a plateau the fit crosses slowly.
//!
//! **A parameter is settled when its distance is below [`SETTLED_FRACTION`] of its standard error**,
//! both on the parameter's own scale and both from the same pass: a distance of 10⁻⁹ in a quantity
//! known to ±0.05 is nothing, whatever the quantity's size. A parameter without a standard error or
//! without a distance — no position's likelihood depends on it, the fit holds it fixed, the data
//! cannot tell it from others, or its error is wider than its whole range — is settled by definition:
//! nothing the fit does to it matters to the likelihood. **The fit has converged when every parameter
//! is settled**, and stops at the end of that cycle.

use super::standard_errors::StandardError;

/// The share of its own standard error a parameter's distance to the maximum must fall below for the
/// parameter to count as settled — the spec's starting value, kept at checkpoint A (owner,
/// 2026-09-28).
pub(super) const SETTLED_FRACTION: f64 = 0.1;

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
