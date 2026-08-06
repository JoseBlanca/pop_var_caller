//! ng step 4 — the parameters the caller runs on, measured from the sample's own
//! loci before anything is called.
//!
//! Four numbers per sample come out of the SNP/indel path: a per-read-group error
//! rate, the sample's heterozygosity, its homozygous-non-reference rate, and its
//! inbreeding coefficient. They are measured from **every** covered position,
//! including the overwhelming majority that show no alternative allele at all —
//! which is what separates this step from production's estimator. Production writes
//! the pure-reference columns; it is production's *heterozygosity accumulator* that
//! never looks at them (`spec/parameter_prepass.md` §2.1), so the loss is in the
//! estimator rather than in the data — and what is lost is the strongest evidence
//! there is about the error rate.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_generic.md` (the design and its
//! rationale), `doc/devel/ng/spec/parameter_prepass.md` (the shared framing), and
//! `doc/devel/ng/arch/parameter_prepass_generic.md` (types and interfaces).
//!
//! Two sub-units, split so that the shaping of data and the mathematics on it never
//! live in one file:
//!
//! - [`fitting`] — the mathematics. Knows nothing about markers, loci or windows: it
//!   is given a table of numbers and returns the values that best explain them. A
//!   folder rather than a file because it is the one genuine swappable seam — one
//!   trait, an implementation on this path and a second on the STR path.
//! - [`generic`] — the SNP/indel path: the two accumulators, the cell table, the
//!   vocabulary they are keyed on, and what each of the four numbers is fitted from.
//!
//! A third sub-unit for the STR path joins them later, which is why the path-specific
//! vocabulary is not here: an error-rate ladder in per-base probabilities and a window
//! size for runs of homozygosity are the SNP/indel path's, and they live in [`generic`]
//! where the STR path will not inherit them. What this file does hold is what **every**
//! parameter step 4 emits carries — where the number came from, and how much data stood
//! behind it — and the step's error type.

pub mod fitting;
pub mod generic;

use crate::ng::parameter_estimation::generic::{MIN_SITES_TO_FIT, MIN_WINDOWS_TO_FIT_INBREEDING};
use crate::ng::types::{DomainError, Ploidy};

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

/// What went wrong while estimating a sample's parameters.
///
/// `#[non_exhaustive]` because the STR path and the two censuses add their own
/// conditions as they land.
///
/// **The three fitting failures are not interchangeable, and only one of the four
/// parameters may be guessed.** An error rate has a ladder of fallbacks — fitted here,
/// borrowed from the sample's other read groups, supplied, defaulted — because
/// chemistry varies far less between runs than biology does between samples. Neither
/// the genotype frequencies nor the inbreeding coefficient has any such rung.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ParameterEstimationError {
    /// Too few sites to fit this sample's genotype frequencies at this ploidy.
    ///
    /// **Not recoverable here.** There is no sibling to borrow from — a sample has one
    /// heterozygosity — and no constant worth inventing, because this is the biology.
    #[error(
        "sample {sample}: {sites} sites at ploidy {ploidy} is too few to fit genotype \
         frequencies (need {MIN_SITES_TO_FIT}); supply them or drop the sample"
    )]
    GenotypeFrequenciesNotFittable {
        sample: String,
        ploidy: Ploidy,
        sites: u64,
    },

    /// `F` was to be fitted and the runs model had too few windows to run on.
    ///
    /// **Deliberately has no default.** Inbreeding is the parameter that differs most
    /// between an outcrosser and a selfing landrace, so any constant would be wrong for
    /// half the runs — and the cohort's diversity divides by `1 − F`, so a wrong one is
    /// amplified rather than absorbed.
    #[error(
        "sample {sample}: {windows} usable windows is too few to fit the inbreeding \
         coefficient (need {MIN_WINDOWS_TO_FIT_INBREEDING}); supply one instead"
    )]
    InbreedingNotFittable { sample: String, windows: usize },

    /// Every starting point emptied one of the runs model's two states, so no
    /// separation between them was found.
    ///
    /// **This is not `F` = 0**, and the distinction is the reason the variant exists.
    /// An outcrossing genome and a search that failed leave identical fitted values —
    /// an empty inside state and its frequencies at their starting guesses — and only
    /// the scores across starting points tell them apart. Returning zero here is the
    /// one way this estimator produces a confidently wrong number rather than a visible
    /// failure.
    #[error(
        "sample {sample}: the runs model found no second state from any of {starts} \
         starting points — this is a search that failed, not an inbreeding coefficient \
         of zero; widen the state separations or supply F"
    )]
    InbreedingStatesNotSeparated { sample: String, starts: usize },

    /// A domain invariant was violated on the way — a rate outside `[0, 1]`, a ploidy
    /// of zero.
    #[error(transparent)]
    Domain(#[from] DomainError),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    /// **Every message names the sample and the number that was too small**, next to
    /// the floor it fell short of. A parameter-estimation failure is read out of a log
    /// on a cohort run of hundreds of samples, so a message that omits which sample, or
    /// omits how far short it fell, sends the reader back to the data to find out.
    #[test]
    fn each_fitting_failure_names_the_sample_and_the_number_that_was_too_small() {
        let frequencies = ParameterEstimationError::GenotypeFrequenciesNotFittable {
            sample: "SL_landrace_07".to_string(),
            ploidy: ploidy(2),
            sites: 812,
        };
        let message = frequencies.to_string();
        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(message.contains("812"), "{message}");
        assert!(
            message.contains("10000"),
            "the floor it fell short of: {message}"
        );

        let inbreeding = ParameterEstimationError::InbreedingNotFittable {
            sample: "HG002_chr20".to_string(),
            windows: 1_200,
        };
        let message = inbreeding.to_string();
        assert!(message.contains("HG002_chr20"), "{message}");
        assert!(message.contains("1200"), "{message}");
        assert!(
            message.contains("3000"),
            "the floor it fell short of: {message}"
        );
    }

    /// The one message that has to say what it is **not**. An outcrosser and a failed
    /// search leave the same fitted values, so a reader who takes this for `F = 0` has
    /// been told a confident wrong number — which is exactly what the variant exists to
    /// prevent.
    #[test]
    fn the_unseparated_states_message_says_it_is_not_an_inbreeding_coefficient_of_zero() {
        let unseparated = ParameterEstimationError::InbreedingStatesNotSeparated {
            sample: "SL_landrace_07".to_string(),
            starts: 9,
        };
        let message = unseparated.to_string();

        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(
            message.contains('9'),
            "how many starts were tried: {message}"
        );
        assert!(
            message.contains("not an inbreeding coefficient"),
            "the reader must not take this for F = 0: {message}"
        );
    }

    /// A domain violation reaching a fit is reported transparently rather than being
    /// re-worded, so the newtype's own message — which names the quantity and the value
    /// — is what the reader sees.
    #[test]
    fn a_domain_violation_passes_through_with_its_own_message() {
        let rejected = ParameterEstimationError::from(DomainError::InbreedingF(1.5));

        assert_eq!(
            rejected.to_string(),
            DomainError::InbreedingF(1.5).to_string()
        );
        assert!(rejected.to_string().contains("1.5"));
    }

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
            value: crate::ng::types::ErrorRate::try_new(0.001).unwrap(),
            provenance: Provenance::FittedHere,
            observations: 80_000,
        };

        assert_eq!(error_rate.value.get(), 0.001);
        assert_eq!(error_rate.provenance, Provenance::FittedHere);
        assert_eq!(error_rate.observations, 80_000);
    }
}
