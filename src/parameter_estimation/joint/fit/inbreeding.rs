//! **The fit's homozygote excess, offered to calling as each sample's inbreeding coefficient.**
//!
//! It lives beside [`HomozygoteExcess`], in the fit that measures it, because it reads nothing
//! else. The `estimate-parameters` driver calls it through `crate::run`, which re-exports it, and
//! the fit's own tests call it without importing the driver.

use std::collections::BTreeMap;

use super::HomozygoteExcess;
use crate::parameter_estimation::Estimate;
use crate::types::InbreedingF;

/// **The inbreeding coefficients this fit measured, by sample name** — the middle rung of
/// [`DeclaredInbreeding::of_each_sample_over`](crate::calling::parameters_file::DeclaredInbreeding::of_each_sample_over)'s
/// ladder.
///
/// The quantity is each sample's **homozygote excess**: how much less heterozygous it is than the
/// allele frequencies this same fit produced predict. **That is circular**, and it is the only
/// fitted source there is since the runs-of-homozygosity estimator was removed
/// (`impl_plan/remove_histogram_route.md`).
///
/// **Keyed by name, never by position**, because the run's sample order is the run's and a list
/// joined by order sends one plant's coefficient to another with nothing to notice.
///
/// It takes the fit's `hom_excess` map rather than the whole fit, because that map is all it
/// reads — and a function that took the fit would need one built to be tested at all.
///
/// # A sample is left out rather than coerced, in one case
///
/// The excess is a fraction in `[0, 1]` and a coefficient is one in `[0, 1)`, so an excess of
/// exactly one has no coefficient to become. Such a sample is **absent from this map** and falls to
/// the rung below, which is the honest handling: a maximisation that returns its own bound has
/// reported where it stopped looking rather than what it measured. *The shipped fit does not return
/// exactly one, so this guards the type's contract rather than an expected case.* It moves an excess
/// two ways, and neither lands on 1: a golden-section search, which reports the middle of a bracket
/// strictly inside `[0, 1]`, and a Newton step, whose target stops half that search's stopping
/// width short of 1 (`Coordinate::newton_range` in `joint::fit`). Before that target was kept inside, it was
/// clamped onto 1 itself, and on 59 tomato samples two nearly fully inbred samples came out of here
/// with no coefficient and were written `defaulted` at 0.
#[must_use]
pub fn fitted_inbreeding_of(
    hom_excess: &BTreeMap<String, Estimate<HomozygoteExcess>>,
) -> BTreeMap<Box<str>, Estimate<InbreedingF>> {
    hom_excess
        .iter()
        .filter_map(|(sample, excess)| {
            InbreedingF::try_new(excess.value.get())
                .ok()
                .map(|coefficient| {
                    (
                        sample.as_str().into(),
                        Estimate {
                            value: coefficient,
                            // **The fit's own warrant travels with the number.** It is
                            // `FittedHere` where the cohort could identify the excess and
                            // `Defaulted` at a single sample, where it comes out zero whatever
                            // the truth — and a reader of the file has to be able to tell those
                            // apart.
                            provenance: excess.provenance,
                            observations: excess.observations,
                            // The coefficient is the excess itself, so its error is the excess's.
                            standard_error: excess.standard_error,
                        },
                    )
                })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    //! **What the fit's homozygote excess becomes when it is offered as a coefficient.**

    use super::*;
    use crate::parameter_estimation::Provenance;

    fn excess(value: f64, provenance: Provenance) -> Estimate<HomozygoteExcess> {
        Estimate {
            value: HomozygoteExcess::try_new(value).expect("a fraction in [0, 1]"),
            provenance,
            observations: 1_806,
            standard_error: None,
        }
    }

    /// **The fit's own warrant travels with its number, and the two warrants it produces mean
    /// different things to a reader.** `FittedHere` is *this cohort identified it*; `Defaulted` is
    /// what a single-sample fit returns, where the excess comes out zero whatever the truth
    /// because one genome's totals cannot identify it. A conversion that stamped everything
    /// `FittedHere` would tell a one-sample run it had measured its plant.
    #[test]
    fn the_fits_warrant_and_evidence_count_travel_with_the_coefficient() {
        let hom_excess: BTreeMap<String, Estimate<HomozygoteExcess>> = [
            (
                "identified".to_string(),
                excess(0.78, Provenance::FittedHere),
            ),
            ("one_sample".to_string(), excess(0.0, Provenance::Defaulted)),
        ]
        .into_iter()
        .collect();
        let mut hom_excess = hom_excess;
        hom_excess
            .get_mut("identified")
            .expect("inserted above")
            .standard_error = Some(0.031);

        let coefficients = fitted_inbreeding_of(&hom_excess);

        assert_eq!(coefficients.len(), 2);
        let identified = &coefficients["identified"];
        assert!((identified.value.get() - 0.78).abs() < 1e-15);
        assert_eq!(identified.provenance, Provenance::FittedHere);
        assert_eq!(identified.observations, 1_806);
        // The coefficient is the excess, so it keeps the excess's error, and an excess with none
        // stays without one.
        assert_eq!(identified.standard_error, Some(0.031));
        assert_eq!(coefficients["one_sample"].provenance, Provenance::Defaulted);
        assert_eq!(coefficients["one_sample"].standard_error, None);
    }

    /// **An excess of exactly one has no coefficient to become, so that sample is left out** and
    /// falls to the rung below rather than being coerced to something just under one — which would
    /// be the same number to every consumer while claiming to be a measurement.
    ///
    /// The excess is a fraction in `[0, 1]` and a coefficient one in `[0, 1)`: at `F = 1` every
    /// genotype is homozygous by construction, which is why the coefficient type excludes it.
    /// **Its neighbour keeps its own number**, which is the property that matters — one
    /// unconvertible sample must not cost the cohort its other coefficients.
    #[test]
    fn a_sample_whose_excess_is_exactly_one_is_left_out_and_its_neighbour_is_not() {
        let hom_excess: BTreeMap<String, Estimate<HomozygoteExcess>> = [
            (
                "at_the_ceiling".to_string(),
                excess(1.0, Provenance::FittedHere),
            ),
            ("ordinary".to_string(), excess(0.42, Provenance::FittedHere)),
        ]
        .into_iter()
        .collect();

        let coefficients = fitted_inbreeding_of(&hom_excess);

        assert!(
            !coefficients.contains_key("at_the_ceiling"),
            "an excess of one is not a coefficient: {coefficients:?}"
        );
        assert!((coefficients["ordinary"].value.get() - 0.42).abs() < 1e-15);
    }

    /// **Just below the ceiling still converts**, so the exclusion above is the type's boundary
    /// and not a range this function narrowed on its own.
    #[test]
    fn an_excess_just_below_one_still_becomes_a_coefficient() {
        let hom_excess: BTreeMap<String, Estimate<HomozygoteExcess>> = [(
            "nearly".to_string(),
            excess(1.0 - 1e-9, Provenance::FittedHere),
        )]
        .into_iter()
        .collect();

        assert_eq!(fitted_inbreeding_of(&hom_excess).len(), 1);
    }
}
