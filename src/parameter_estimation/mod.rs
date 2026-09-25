//! ng step 4 — the parameters the caller runs on, measured from the cohort's own loci before
//! anything is called.
//!
//! **One route, and it is [`joint`]**: every parameter fitted once, over every sample at the same
//! bounded set of positions — the *census* — so that a position's own allele frequency in the
//! population is a quantity the fit can weigh a genotype against. What comes out is a per-library
//! error rate and mismapped-position rate, the population's allele-frequency curve, each sample's
//! departure from Hardy–Weinberg proportions, each library's contamination, and the repeat-tract
//! slippage numbers per stratum.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass.md` (the shared framing),
//! `parameter_prepass_joint_fit.md` (what the route is), `parameter_prepass_joint_loci.md` (which
//! positions) and `parameter_prepass_joint_records.md` (what is recorded at each).
//!
//! # There used to be a second route, and it was removed on 2026-09-11
//!
//! The **per-sample whole-genome histogram** route folded each sample's every covered position
//! into histograms and fitted from those, one sample at a time. It was never run by any shipped
//! command, and the join that would have handed its results to the caller had no caller outside
//! its own tests. What only it could produce was an inbreeding coefficient read off *which windows
//! of the genome lie in a run of homozygosity* — a local quantity a census cannot support, since a
//! census window holds about a quarter of one heterozygote at the shipped budget.
//!
//! **The coefficient a caller reads is now the per-sample departure from Hardy–Weinberg**, which
//! [`joint::fit`] measures from the census because it is an average over positions rather than a
//! local quantity. That estimator is **circular** — it is measured against a population
//! expectation the same fit produced — and `joint::census_moments`'s output says so rather than
//! hiding it. Removing the runs estimator gave that up knowingly
//! (`impl_plan/remove_histogram_route.md`, and the measurement behind the decision is
//! `doc/devel/reports/ng_census_inbreeding_budget_2026-09-11.md`).
//!
//! # What sits here rather than in [`joint`]
//!
//! Three things the census route needs that are not the fit, and that calling reads too — so they
//! are shared vocabulary rather than a route's internals:
//!
//! - [`depth_bins`] — the depth ladder a census codes its depths on, and the ladder a run's cells
//!   are keyed by. Two samples binned under different edges hold codes that mean different
//!   depths, which is why the ladder's digest is one of the settings a fit refuses on.
//! - [`calibration`] — what a library's own base qualities *claimed*, summed as the loci go past.
//!   A base-quality calibration is fitted from that and the measured error rate together.
//! - [`repeat_strata`] — how a repeat tract's stratum is named: a motif period, the reference's
//!   repeat count, and the library and ploidy that make one fitted set of slippage numbers.
//!
//! And what **every** parameter this step emits carries: where the number came from
//! ([`Provenance`]) and how much data stood behind it ([`Estimate`]).

pub mod calibration;
pub mod depth_bins;
pub mod joint;
pub(crate) mod progress;
pub mod repeat_strata;

/// Where a parameter came from.
///
/// **Not an error condition.** A rate fitted on 80,000 reads and one borrowed from a
/// neighbouring read group are both usable, and the consumer has to be able to tell
/// them apart — a consumer that treats all four alike is the failure this exists to
/// prevent, because a defaulted error rate is a guess and a fitted one is a
/// measurement.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Provenance {
    /// Fitted from this sample's own sites, at this grain.
    FittedHere,
    /// Too little data here, so the mean of the sample's other read groups was taken.
    /// Chemistry differs between libraries, which is the whole reason for the
    /// read-group grain, so this is a compromise and is marked as one.
    Borrowed,
    /// Nothing could be fitted and nothing was supplied, so a stated constant was used.
    Defaulted,
    /// The run was given this value rather than fitting it.
    Supplied,
}

impl Provenance {
    /// The weaker of two warrants — what a value derived from both is entitled to claim.
    ///
    /// **Consumers combine provenances; they do not branch on them** (`spec/read_likelihoods.md`
    /// §4.4). A score resting on one fitted parameter and one borrowed parameter is a borrowed
    /// score, and saying otherwise would launder the weaker of the two — the failure this enum's
    /// own documentation exists to prevent.
    ///
    /// **The order is the ladder this step has always stated**: *fitted here, borrowed from the
    /// sample's other read groups, supplied, defaulted*. So [`Self::FittedHere`] is the strongest
    /// and
    /// [`Self::Defaulted`] the weakest, and **[`Self::Supplied`] sits below [`Self::Borrowed`]**
    /// — a number the run was handed says nothing about this data, where a borrowed one is at
    /// least a measurement of a neighbouring grain.
    ///
    /// *(That last placement is the ladder's, not this method's invention. If a run's supplied
    /// values should outrank a borrowed fit, the ladder is where to change it and this follows.)*
    #[must_use]
    pub fn weaker_of(self, other: Self) -> Self {
        if other.strength() < self.strength() {
            other
        } else {
            self
        }
    }

    /// Where this warrant sits on the ladder — higher is better founded. Private, because the
    /// number is an implementation of [`Self::weaker_of`] and not a quantity to expose.
    fn strength(self) -> u8 {
        match self {
            Self::FittedHere => 3,
            Self::Borrowed => 2,
            Self::Supplied => 1,
            Self::Defaulted => 0,
        }
    }
}

/// A fitted number with its warrant: what it is, where it came from, and how much data
/// stood behind it.
///
/// Generic over the quantity, so an `Estimate<ErrorRate>` and an `Estimate<InbreedingF>`
/// stay unmixable — the warrant travels without erasing which quantity it is a warrant
/// for.
///
/// **No uncertainty interval.** These are priors; a caller mixes them into a genotype
/// prior rather than reporting them, and an interval on a prior is not a quantity any
/// consumer in the design reads.
#[derive(Clone, PartialEq, Debug)]
pub struct Estimate<T> {
    pub value: T,
    pub provenance: Provenance,
    /// Reads for a per-read rate, sites for a per-site one. The unit follows the
    /// quantity, which is why it is not named in the type.
    pub observations: u64,
}

/// The per-base error rate used when none could be fitted and none was supplied.
///
/// **Soft, and the only defaulted parameter of this step.** Chemistry varies far less between runs
/// than biology does between samples, so a stated constant is defensible here in a way it is not
/// for the population's diversity or for inbreeding — which is why a run without those says so
/// rather than inventing them.
///
/// A library that takes this rate also takes a defaulted base-quality calibration, because the two
/// come out of one pass over one set of reads: `RunParameters::assemble` refuses a fitted rate with
/// no accumulator total behind it, and says why.
pub const DEFAULT_ERROR_RATE: f64 = 0.001;

#[cfg(test)]
mod tests {
    use super::*;

    /// The four provenances are distinct values, not a scale — `Borrowed` is not
    /// "better" than `Defaulted` in any ordering the type imposes, and deliberately so:
    /// what a consumer should do with each is the consumer's decision.
    #[test]
    fn the_four_provenances_are_distinct() {
        let all = [
            Provenance::FittedHere,
            Provenance::Borrowed,
            Provenance::Defaulted,
            Provenance::Supplied,
        ];
        for (index, one) in all.iter().enumerate() {
            for other in &all[index + 1..] {
                assert_ne!(one, other);
            }
        }
    }

    /// An `Estimate` carries the quantity in its type, so the warrant cannot be moved
    /// from one parameter to another by accident.
    #[test]
    fn an_estimate_carries_its_quantity_provenance_and_observation_count() {
        let error_rate = Estimate {
            value: crate::types::ErrorRate::try_new(0.001).unwrap(),
            provenance: Provenance::FittedHere,
            observations: 80_000,
        };

        assert_eq!(error_rate.value.get(), 0.001);
        assert_eq!(error_rate.provenance, Provenance::FittedHere);
        assert_eq!(error_rate.observations, 80_000);
    }
}
