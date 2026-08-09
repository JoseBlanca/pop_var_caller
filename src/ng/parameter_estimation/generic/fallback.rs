//! What each parameter falls back to when its own data will not carry it — and which of
//! them are allowed to fall back at all.
//!
//! **The three fitting failures are not interchangeable, and only one of the four
//! parameters may be guessed** (`arch/parameter_prepass_generic.md` §5.4):
//!
//! - **The error rate has a ladder**: fitted from this read group's own sites, borrowed
//!   from the sample's other groups, supplied by the run, or defaulted to a stated
//!   constant. Chemistry varies far less between sequencing runs than biology does between
//!   samples, so a constant is defensible here.
//! - **The genotype frequencies have no rung at all.** A sample has one heterozygosity, so
//!   there is no sibling to borrow it from, and no constant is worth inventing because it
//!   is the biology. Too few sites is an error, raised by the coupled fit.
//! - **The inbreeding coefficient has one rung and it is not a default.** It may be
//!   *supplied* — which is what [`InbreedingMode::Supplied`] carries, and it changes what
//!   is accumulated rather than only what is reported — and otherwise it is fitted or it
//!   fails. It is the parameter that differs most between an outcrosser and a selfing
//!   landrace, and a cohort's diversity divides by `1 − F`, so a wrong constant would be
//!   amplified rather than absorbed.
//!
//! Every value carries its [`Provenance`], because a consumer that treats a measurement and
//! a guess alike is the failure that record exists to prevent.

use std::collections::BTreeMap;

use super::DEFAULT_ERROR_RATE;
use super::accumulators::InbreedingMode;
use super::read_group_error_rate::ReadGroupErrorRateFit;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::types::{ErrorRate, InbreedingF, ReadGroupId};

/// Settle each read group's error rate, walking the ladder for any group whose own sites
/// are too few.
///
/// `fitted` is what the coupled fit arrived at, one entry per read group with a table;
/// `reads_of_group` is how many read observations each group contributed, which is the unit
/// an error rate's warrant is counted in (§2.4). `supplied` is what the run was told, and
/// `floor` is ordinarily [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT).
///
/// **Every group in `reads_of_group` gets an answer**, including one with no fit at all —
/// that is the case the ladder's lower rungs exist for, and a group silently absent from
/// the result is a sample whose reads are scored against a rate nobody chose.
///
/// **The borrowed value is the plain mean of the qualifying groups' rates**, not a
/// read-weighted one. What is being borrowed is a property of chemistry that this group's
/// own data could not resolve, and weighting by the *lenders'* yields would say the
/// deepest library's chemistry is the most likely to be this one's.
///
/// **⚠ The order is the design's, and what it turns on is what *supplied* means.** `arch`
/// §5.4 and the plan both list *fitted here → borrowed → supplied → defaulted*, so a thin
/// group in a sample with one good group takes the borrowed rate even when the run supplied
/// one.
///
/// Read as a **fallback**, that order is principled and matches the rest of the ladder —
/// measurement over assertion, the same reason a fit outranks a borrow. A borrowed rate is
/// still a number measured from this sample's own reads; a supplied one is an assertion
/// from outside the data. And the docs do phrase it that way: rung 3 is "Supplied, **if**
/// the run was given one", and the default applies "when none can be fitted and **none was
/// supplied**".
///
/// Read as an **override**, it is wrong, and arch §5.4 supplies the counterargument itself:
/// borrowing is "a compromise and is marked as one", because chemistry differs between
/// libraries. An operator who supplies a rate *because* this library is unlike its siblings
/// — a different platform in the same sample — gets the neighbour's number instead, and
/// `Provenance::Borrowed` lets a consumer notice without letting anyone fix it.
///
/// **Settled (owner, 2026-08-09), and F1 is what settled it**: the configuration field is
/// named `fallback_error_rates`, so by this comment's own rule the order stands as
/// implemented. `GenericEstimationConfig::fallback_error_rates` is now a production source
/// for the third rung, exercised by `a_supplied_rate_is_used_where_no_group_can_fit_or_lend`.
///
/// # Panics
///
/// If a computed mean is not a probability, which it cannot be: it is the mean of values
/// already inside `[0, 1]`.
#[must_use]
pub fn resolve_error_rates(
    fitted: &BTreeMap<ReadGroupId, ReadGroupErrorRateFit>,
    reads_of_group: &BTreeMap<ReadGroupId, u64>,
    supplied: &BTreeMap<ReadGroupId, ErrorRate>,
    floor: u64,
) -> BTreeMap<ReadGroupId, Estimate<ErrorRate>> {
    // Which groups' own fits stand on enough sites to be worth keeping — computed once,
    // because every borrowing group asks the same question about the same set.
    let groups_above_the_floor: Vec<(ReadGroupId, ErrorRate)> = fitted
        .iter()
        .filter(|(_, fit)| fit.sites >= floor)
        .map(|(&group, fit)| (group, fit.error_rate))
        .collect();

    reads_of_group
        .keys()
        .map(|&group| {
            let own = fitted.get(&group).filter(|fit| fit.sites >= floor);
            if let Some(fit) = own {
                return (
                    group,
                    Estimate {
                        value: fit.error_rate,
                        provenance: Provenance::FittedHere,
                        observations: reads_of_group[&group],
                    },
                );
            }

            // **The lenders are the *other* groups**, which matters only when this group
            // qualifies too — and it cannot, since it is here. Filtering anyway, so that
            // the rule reads as what it is rather than as a consequence of the branch
            // above.
            let lenders_rates: Vec<f64> = groups_above_the_floor
                .iter()
                .filter(|&&(lender, _)| lender != group)
                .map(|&(_, rate)| rate.get())
                .collect();
            if !lenders_rates.is_empty() {
                let mean = lenders_rates.iter().sum::<f64>() / lenders_rates.len() as f64;
                return (
                    group,
                    Estimate {
                        value: ErrorRate::try_new(mean)
                            .expect("the mean of probabilities is a probability"),
                        provenance: Provenance::Borrowed,
                        // **The lenders' reads, not this group's.** The warrant belongs to
                        // the evidence the number came from, and this group's own reads are
                        // precisely the ones that were not enough.
                        observations: groups_above_the_floor
                            .iter()
                            .filter(|&&(lender, _)| lender != group)
                            .filter_map(|(lender, _)| reads_of_group.get(lender))
                            .sum(),
                    },
                );
            }

            if let Some(&rate) = supplied.get(&group) {
                return (
                    group,
                    Estimate {
                        value: rate,
                        provenance: Provenance::Supplied,
                        observations: 0,
                    },
                );
            }

            (
                group,
                Estimate {
                    value: ErrorRate::try_new(DEFAULT_ERROR_RATE)
                        .expect("the default error rate is a probability"),
                    provenance: Provenance::Defaulted,
                    observations: 0,
                },
            )
        })
        .collect()
}

/// Take the inbreeding coefficient a run was handed, if it was handed one.
///
/// **Its own function rather than a branch inside the fit**, because a supplied `F` changes
/// what the accumulator *stores*: [`InbreedingMode::Supplied`] drops the window key
/// altogether, collapsing the windowed table from about 37 MB to a few kB, so by the time
/// this is asked there is no chain to walk and the runs model cannot be run even if a
/// caller wanted it.
///
/// `covered_positions` says how much genome the number is being applied to, which is the only
/// quantity a consumer can compare across samples.
///
/// **It is not quite what the fitted value's warrant would have been, and on a mixed-ploidy
/// genome it is larger.** A fitted `F` is warranted by the **diploid** windows' covered
/// positions, because that is the only part of the genome its chain walked; this counts every
/// ploidy, since a supplied coefficient was not fitted from any of them. On a fixture whose
/// two arms are equal that is 20,002 positions against the 10,001 a fit would have reported.
/// A consumer comparing the two across samples is comparing different denominators.
#[must_use]
pub fn take_supplied_inbreeding(
    mode: InbreedingMode,
    covered_positions: u64,
) -> Option<Estimate<InbreedingF>> {
    match mode {
        InbreedingMode::Fitted => None,
        InbreedingMode::Supplied(value) => Some(Estimate {
            value,
            provenance: Provenance::Supplied,
            observations: covered_positions,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::generic::MIN_SITES_TO_FIT;
    use crate::ng::types::LogProb;

    fn rate(value: f64) -> ErrorRate {
        ErrorRate::try_new(value).expect("a probability")
    }

    fn fit(error_rate: f64, sites: u64) -> ReadGroupErrorRateFit {
        ReadGroupErrorRateFit {
            error_rate: rate(error_rate),
            rung: 80,
            log_likelihood: LogProb(-1.0),
            argmax_at_ladder_end: false,
            sites,
        }
    }

    fn reads(entries: &[(u32, u64)]) -> BTreeMap<ReadGroupId, u64> {
        entries
            .iter()
            .map(|&(group, reads)| (ReadGroupId(group), reads))
            .collect()
    }

    /// **A group with enough sites of its own keeps its own rate**, and the warrant it
    /// carries is its own reads.
    #[test]
    fn a_group_that_stands_on_its_own_sites_is_fitted_here() {
        let fitted = BTreeMap::from([(ReadGroupId(1), fit(0.001, MIN_SITES_TO_FIT))]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 4_000_000)]),
            &BTreeMap::new(),
            MIN_SITES_TO_FIT,
        );

        assert_eq!(resolved[&ReadGroupId(1)].provenance, Provenance::FittedHere);
        assert_eq!(resolved[&ReadGroupId(1)].value, rate(0.001));
        assert_eq!(resolved[&ReadGroupId(1)].observations, 4_000_000);
    }

    /// **A thin group borrows the mean of the sample's other groups**, and it is the mean
    /// of the *others* — a group that borrowed from itself would keep the number the floor
    /// just rejected.
    ///
    /// Two lenders at different rates, so a rule that took one of them rather than their
    /// mean lands on neither. **The two lenders sit exactly *at* the floor**, and that is
    /// what makes its comparison reachable: `>=` accepts them and `>` does not. A group one
    /// site below the floor — which the thin group is — is rejected by both, so it says
    /// nothing about the operator.
    #[test]
    fn a_thin_group_borrows_the_mean_of_the_others() {
        let fitted = BTreeMap::from([
            (ReadGroupId(1), fit(0.001, MIN_SITES_TO_FIT)),
            (ReadGroupId(2), fit(0.003, MIN_SITES_TO_FIT * 2)),
            (ReadGroupId(3), fit(0.05, MIN_SITES_TO_FIT - 1)),
        ]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 4_000_000), (2, 8_000_000), (3, 900)]),
            &BTreeMap::new(),
            MIN_SITES_TO_FIT,
        );

        let borrowed = &resolved[&ReadGroupId(3)];
        assert_eq!(borrowed.provenance, Provenance::Borrowed);
        assert!(
            (borrowed.value.get() - 0.002).abs() < 1e-15,
            "the mean of 0.001 and 0.003, not either of them: {}",
            borrowed.value.get()
        );
        assert_eq!(
            borrowed.observations, 12_000_000,
            "the warrant is the lenders' reads, not the 900 that were too few"
        );
        assert_eq!(resolved[&ReadGroupId(1)].provenance, Provenance::FittedHere);
        assert_eq!(resolved[&ReadGroupId(2)].provenance, Provenance::FittedHere);
    }

    /// **With every group thin there is nothing to borrow**, so the ladder falls through to
    /// what the run supplied, and to the default where it supplied nothing. Both groups are
    /// in the same sample, so this also says the two lower rungs are decided per group.
    #[test]
    fn with_every_group_thin_the_ladder_falls_through_to_supplied_then_defaulted() {
        let fitted = BTreeMap::from([
            (ReadGroupId(1), fit(0.05, 10)),
            (ReadGroupId(2), fit(0.05, 10)),
        ]);
        let supplied = BTreeMap::from([(ReadGroupId(1), rate(0.002))]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 100), (2, 100)]),
            &supplied,
            MIN_SITES_TO_FIT,
        );

        assert_eq!(resolved[&ReadGroupId(1)].provenance, Provenance::Supplied);
        assert_eq!(resolved[&ReadGroupId(1)].value, rate(0.002));
        assert_eq!(resolved[&ReadGroupId(2)].provenance, Provenance::Defaulted);
        assert!((resolved[&ReadGroupId(2)].value.get() - DEFAULT_ERROR_RATE).abs() < 1e-15);
        for group in [1, 2] {
            assert_eq!(
                resolved[&ReadGroupId(group)].observations,
                0,
                "a value nothing in this sample stood behind carries no observations"
            );
        }
    }

    /// **A group with no fit at all still gets an answer.** It is the case the lower rungs
    /// exist for, and a group missing from the result would be a set of reads scored
    /// against a rate nobody chose.
    #[test]
    fn a_group_with_no_fit_at_all_is_still_answered() {
        let fitted = BTreeMap::from([(ReadGroupId(1), fit(0.001, MIN_SITES_TO_FIT))]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 4_000_000), (2, 5)]),
            &BTreeMap::new(),
            MIN_SITES_TO_FIT,
        );

        assert_eq!(resolved.len(), 2);
        assert_eq!(resolved[&ReadGroupId(2)].provenance, Provenance::Borrowed);
        assert_eq!(resolved[&ReadGroupId(2)].value, rate(0.001));
    }

    /// **A supplied rate does not displace a group's own fit.** `fitted here` is the top
    /// rung, so a run that supplied a rate for a group holding enough sites of its own is
    /// told what was measured rather than what it asked for — the same ordering question as
    /// `borrowing_outranks_a_supplied_rate`, asked one rung higher.
    #[test]
    fn a_supplied_rate_does_not_displace_a_group_fitted_on_its_own_sites() {
        let fitted = BTreeMap::from([(ReadGroupId(1), fit(0.001, MIN_SITES_TO_FIT))]);
        let supplied = BTreeMap::from([(ReadGroupId(1), rate(0.02))]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 4_000_000)]),
            &supplied,
            MIN_SITES_TO_FIT,
        );

        assert_eq!(resolved[&ReadGroupId(1)].provenance, Provenance::FittedHere);
        assert_eq!(resolved[&ReadGroupId(1)].value, rate(0.001));
        assert_eq!(
            resolved[&ReadGroupId(1)].observations,
            4_000_000,
            "the warrant is the group's own reads, not the supplied value's nothing"
        );
    }

    /// The ladder's order is the design's: **borrowed outranks supplied**. Recorded on a
    /// test rather than only in prose, so that reordering it is a visible decision.
    #[test]
    fn borrowing_outranks_a_supplied_rate() {
        let fitted = BTreeMap::from([
            (ReadGroupId(1), fit(0.001, MIN_SITES_TO_FIT)),
            (ReadGroupId(2), fit(0.05, 10)),
        ]);
        let supplied = BTreeMap::from([(ReadGroupId(2), rate(0.02))]);

        let resolved = resolve_error_rates(
            &fitted,
            &reads(&[(1, 4_000_000), (2, 100)]),
            &supplied,
            MIN_SITES_TO_FIT,
        );

        assert_eq!(resolved[&ReadGroupId(2)].provenance, Provenance::Borrowed);
        assert_eq!(resolved[&ReadGroupId(2)].value, rate(0.001));
    }

    /// A run handed an inbreeding coefficient reports it as **supplied**, with the genome
    /// it applies to behind it; a run fitting one has nothing to report here.
    #[test]
    fn a_supplied_inbreeding_coefficient_is_marked_as_one() {
        let value = InbreedingF::try_new(0.4).expect("a fraction");

        let supplied = take_supplied_inbreeding(InbreedingMode::Supplied(value), 800_000_000)
            .expect("a supplied mode yields a value");

        assert_eq!(supplied.value, value);
        assert_eq!(supplied.provenance, Provenance::Supplied);
        assert_eq!(supplied.observations, 800_000_000);
        assert!(take_supplied_inbreeding(InbreedingMode::Fitted, 800_000_000).is_none());
    }
}
