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
    ///
    /// `floor` travels on the variant rather than being read from a constant, so that
    /// the message quotes the floor the *caller* applied. This error is
    /// `#[non_exhaustive]` and shared with the STR path, whose floors are its own; a
    /// message hard-wired to the SNP/indel constant would quote the wrong number at the
    /// wrong grain the first time the STR path raises it.
    #[error(
        "sample {sample}: {sites} sites at ploidy {ploidy} is too few to fit genotype \
         frequencies (need {floor}); supply them or drop the sample"
    )]
    GenotypeFrequenciesNotFittable {
        sample: String,
        ploidy: Ploidy,
        sites: u64,
        floor: u64,
    },

    /// `F` was to be fitted and the runs model had too few windows to run on.
    ///
    /// **Deliberately has no default.** Inbreeding is the parameter that differs most
    /// between an outcrosser and a selfing landrace, so any constant would be wrong for
    /// half the runs — and the cohort's diversity divides by `1 − F`, so a wrong one is
    /// amplified rather than absorbed.
    #[error(
        "sample {sample}: {windows} usable windows is too few to fit the inbreeding \
         coefficient (need {floor}); supply one instead"
    )]
    InbreedingNotFittable {
        sample: String,
        windows: usize,
        floor: usize,
    },

    /// A sample's genotype frequencies were built with the wrong number of entries for
    /// their ploidy, or with entries that do not sum to one.
    ///
    /// **This is our own arithmetic being broken, not a data condition** — which is why
    /// the fits construct through the checked door and `.expect()`. It is here rather
    /// than in [`DomainError`] because the invariant is about a *set* of frequencies
    /// against a ploidy, which no single constrained scalar can express: entry `k` of
    /// the set is the rate for `k` non-reference copies, so a set one entry short hands
    /// back the wrong dosage's rate under the right name.
    #[error(
        "a ploidy-{ploidy} sample needs one genotype frequency per dosage 0..={ploidy}, \
         got {entries} summing to {total}"
    )]
    GenotypeFrequenciesOffSimplex {
        ploidy: Ploidy,
        entries: usize,
        total: f64,
    },

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

    /// The runs model's starting points scored alike and answered differently, so the
    /// data did not determine an inbreeding coefficient.
    ///
    /// **Also not `F` = 0, and a different failure from
    /// [`Self::InbreedingStatesNotSeparated`].** There the search never found a second
    /// state; here it found one from every start and found a *different* one each time,
    /// with nothing to choose between them — which is what a chain reading sampling noise
    /// looks like from the outside. Measured on genomes drawn with no runs at all, where
    /// the two states came out at 0.31 and 0.62 of each other, well inside the ratio the
    /// other check refuses on, and nine starts returned `F` from 0.0003 to 0.8497 while
    /// scoring **within 0.91 nats of one another**
    /// (`doc/devel/ng/research/inbreeding_resolution_2026-08-09.md` §1, §4, §6).
    ///
    /// **Only the tied starts count**, which is what keeps this from refusing a fit the
    /// data did determine: on a genome with real runs and a floor of spurious
    /// heterozygotes five times the real rate, six starts agree on `F` = 0.3157 and three
    /// collapse to zero — and those three score 1,473 nats worse, so they have been
    /// rejected rather than disagreed with.
    #[error(
        "sample {sample}: {tied_starts} of the runs model's {starts} starting points \
         scored alike and disagreed about the inbreeding coefficient by {spread:.4}, more \
         than the {threshold} that leaves it identified — this is a genome the model could \
         not read, not an inbreeding coefficient of zero; supply F instead"
    )]
    InbreedingStartsDisagree {
        sample: String,
        /// How many of them could not be told apart by score — the ones the spread is
        /// over.
        tied_starts: usize,
        /// How many were tried in all.
        starts: usize,
        /// The largest fitted `F` minus the smallest, across the tied starting points.
        spread: f64,
        /// What it was compared against, so the message does not have to be read against
        /// a constant the reader has to go and find.
        threshold: f64,
    },

    /// A constrained scalar rejected its value while a named fit was running — a rate
    /// outside `[0, 1]`, a ploidy of zero.
    ///
    /// **Not transparent, and not `#[from]`, and both for the same reason.** The inner
    /// [`DomainError`] names the quantity and the offending value, which is half of what a
    /// reader needs; what it cannot know is *whose* data and *which* fit produced it, and
    /// on a cohort run of hundreds of samples that is the half that locates the fault. A
    /// transparent variant forwards the inner message unchanged and drops both. A `#[from]`
    /// conversion is worse than the missing context alone: it makes `?` silently mint this
    /// variant at every one of the five constructors that can raise a `DomainError`, so the
    /// sample and the fit would have to be *remembered* to be attached, and the compiler
    /// would never ask. Constructing it by hand is what forces each site to say where it
    /// was.
    #[error("sample {sample}: {fit} rejected a value — {source}")]
    Domain {
        sample: String,
        /// Which fit was running, in the words the emitted summary uses — "the error-rate
        /// scan", "the runs model". A reader who has the sample and the fit can find the
        /// site; one who has only the quantity cannot.
        fit: &'static str,
        source: DomainError,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::generic::MIN_SITES_TO_FIT;
    use crate::ng::parameter_estimation::generic::runs::MIN_WINDOWS_TO_FIT_INBREEDING;

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
            floor: MIN_SITES_TO_FIT,
        };
        let message = frequencies.to_string();
        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(message.contains("812"), "{message}");
        assert!(
            message.contains("ploidy 2"),
            "the ploidy renders: {message}"
        );
        // Against the constant's own rendering rather than a literal: `contains("10000")`
        // is also true of a floor of 100,000, so a floor raised tenfold would have left
        // this assertion green.
        assert!(
            message.contains(&MIN_SITES_TO_FIT.to_string()),
            "the floor it fell short of: {message}"
        );

        let inbreeding = ParameterEstimationError::InbreedingNotFittable {
            sample: "HG002_chr20".to_string(),
            windows: 1_200,
            floor: MIN_WINDOWS_TO_FIT_INBREEDING,
        };
        let message = inbreeding.to_string();
        assert!(message.contains("HG002_chr20"), "{message}");
        assert!(message.contains("1200"), "{message}");
        assert!(
            message.contains(&MIN_WINDOWS_TO_FIT_INBREEDING.to_string()),
            "the floor it fell short of: {message}"
        );
    }

    /// A set of genotype frequencies can be wrong for its ploidy in two ways, and the
    /// message has to say which. Entry `k` is the rate for `k` non-reference copies, so
    /// a set one entry short is not merely incomplete — it hands back the wrong dosage's
    /// rate under the right name.
    #[test]
    fn the_off_simplex_message_names_the_ploidy_the_entry_count_and_the_total() {
        let malformed = ParameterEstimationError::GenotypeFrequenciesOffSimplex {
            ploidy: ploidy(2),
            entries: 1,
            total: 1.0,
        };
        let message = malformed.to_string();

        assert!(message.contains("ploidy-2"), "{message}");
        assert!(message.contains("0..=2"), "the dosages expected: {message}");
        assert!(message.contains("got 1"), "{message}");
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

    /// A domain violation carries **three** things a reader needs, and the inner error has
    /// only one of them. The quantity and the offending value come from the newtype; the
    /// sample and the fit have to be attached here, because on a cohort run of hundreds of
    /// samples those are what say where to look.
    #[test]
    fn a_domain_violation_names_the_sample_and_the_fit_as_well_as_the_quantity() {
        let rejected = ParameterEstimationError::Domain {
            sample: "SL_landrace_07".to_string(),
            fit: "the runs model",
            source: DomainError::InbreedingF(1.5),
        };
        let message = rejected.to_string();

        assert!(message.contains("SL_landrace_07"), "the sample: {message}");
        assert!(message.contains("the runs model"), "the fit: {message}");
        assert!(message.contains("1.5"), "the offending value: {message}");
        assert!(
            message.contains("inbreeding coefficient"),
            "the quantity, in the newtype's own words: {message}"
        );
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
