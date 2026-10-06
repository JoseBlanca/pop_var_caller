//! **Fitting a cohort of censuses** — both halves, in the order the estimator needs them.
//!
//! # Two halves, one after the other
//!
//! The **generic half**
//! ([`fit_jointly`](crate::parameter_estimation::joint::fit::fit_jointly)) reads the
//! ordinary positions: the per-read-group noise
//! rates, the population's allele-frequency density, each sample's departure from Hardy–Weinberg
//! proportions, and each library's contamination.
//!
//! The **repeat-tract half** reads the kept tracts, one stratum at a time, and fits the slippage
//! numbers. It consumes exactly one thing from the generic half — each sample's homozygote
//! excess, which weights a genotype drawn from a locus's length frequencies — and hands nothing
//! back (`parameter_prepass_joint_records.md` §6.2). That is why the generic records can be
//! dropped before a single tract record is read.
//!
//! # Why this needs the reference and the catalog
//!
//! **A census stores a tract by its index within its stratum and nothing else** — no coordinate,
//! no stratum — so the order has to be rebuilt from the same kept-loci object the writer was
//! given ([`strata_of_kept_loci`]).
//! Rebuilding it means choosing the selection again, which is a function of the seed, the
//! reference, the analysed ground and the catalog.
//!
//! **The rebuild is checked rather than trusted.** Every census carries a digest of the loci it
//! was written against, and a selection rebuilt from another reference or another catalog would
//! index one stratum's tracts by another's — a wrong answer with no symptom. So the rebuilt loci
//! are digested and compared before a tract is read, and a mismatch is refused.

use std::collections::BTreeMap;

use crate::calling::parameters_file::{
    CensusIdentity, ParametersFile, ReadsBehindEachCalibration, SnpIndelFitStart, StartOutcome,
};
use crate::calling::run_parameters::RunParameters;
use crate::parameter_estimation::calibration::MintedReadErrors;
use crate::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, RecordingTerms,
};
use crate::parameter_estimation::joint::contamination::ContaminationEstimate;
use crate::parameter_estimation::joint::fit::{
    HomozygoteExcess, JointFit, JointFitConfig, JointFitError, StartEnding, StartRecord,
};
use crate::parameter_estimation::joint::loci::ReferenceDigest;
use crate::parameter_estimation::joint::loci::{CensusLoci, CensusLociDigester};
use crate::parameter_estimation::joint::sample_order::{SAMPLE_ORDER_SEED, sample_order};
use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::parameter_estimation::joint::ssr_fit::{
    self, SsrFitConfig, StratumEvidence, StratumOutcome, StratumSubstitutionCounts, gather_strata,
    strata_of_kept_loci,
};
use crate::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::parameter_estimation::repeat_strata::{RepeatCount, Stratum as SsrStratum, StratumKey};
use crate::parameter_estimation::{Estimate, Provenance};
use crate::read::input::read_groups::ReadGroups;
use crate::repeat_catalog::StrRepeatCriteria;
use crate::run::census_freshness::THE_COMMAND_THAT_REBUILDS_A_CENSUS;
use crate::types::{ContigId, ErrorRate, InbreedingF, Ploidy, ReadGroupId, SsrPeriod};

/// What a fit over a cohort of censuses produced.
#[derive(Debug)]
#[must_use]
pub struct CohortFit {
    /// The generic half: noise rates, the allele-frequency density, contamination, and each
    /// sample's homozygote excess.
    pub generic: JointFit,
    /// The repeat-tract half, one outcome a stratum — fitted, furnished from its period's curve,
    /// or refused for want of tracts.
    pub strata: Vec<StratumOutcome>,
    /// How many kept tracts the strata were rebuilt over, which is what says whether the
    /// repeat-tract half had anything to read.
    pub tracts: usize,
    /// Each stratum's bases compared and bases mismatching, one entry a stratum the tract half
    /// gathered, in the order it gathered them.
    ///
    /// **Measured from the evidence, not fitted**: a stratum whose fit was refused still has
    /// these counts, and its substitution rate is written from them. They are the only part of
    /// the evidence the parameters file needs, which is why the evidence itself can be dropped
    /// as soon as the fit returns.
    pub substitution_counts: Vec<StratumSubstitutionCounts>,
}

/// Why a cohort could not be fitted.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CohortFitError {
    /// The selection rebuilt here is not the one the censuses were written against.
    ///
    /// **The third of spec §4.2's three causes, and the only one a cheap read cannot reach**: it
    /// needs the run's reference read and its selection rebuilt. **`estimate-parameters` makes
    /// that comparison itself and reports it as rows of
    /// [`CensusesToRegenerate`](super::CensusesToRegenerate)** (plan step C5), so what reaches
    /// here from a command is nothing, and what reaches here at all is a caller that fitted
    /// without comparing.
    ///
    /// **Two different faults end here, and they have opposite fixes**, so the message offers
    /// both — in the order of what they cost:
    ///
    /// - **the run was pointed at another reference or catalog** from the ones the psps were
    ///   walked against. Fitting with the right ones regenerates nothing. Regenerating instead
    ///   would rebuild every census against the wrong reference — a quarter of an hour a sample at
    ///   whole-genome scale (spec §8) — and the next fit, with the right reference, would be
    ///   refused again;
    /// - **this build chooses census positions differently** from the one that wrote them — its
    ///   selection seed, budget or cap changed. Only regenerating fixes that.
    ///
    /// **Which of the two it is, this variant cannot say** — what is compared is a digest of the
    /// kept positions, and it records no reason — **and `estimate-parameters` no longer reaches
    /// it.** Before fitting, that command compares its reference and its catalog with the psp
    /// headers, refusing either in its own words, and then each census's twelve recorded settings
    /// with the ones its own plan records under, naming every stale sample and the first setting
    /// that differs (plan step C5 of `psp_census_pair.md`). This is the backstop behind those, and
    /// it is reached today only by this module's own tests — `estimate-parameters` is the one
    /// caller there is. For a caller that fits without comparing, both faults remain possible,
    /// which is why the message keeps both ways out.
    #[error(
        "the census positions chosen from this reference and catalog are not the ones these \
         censuses were written against, and fitting them anyway would read one stratum's tracts \
         as another's; if the psps were walked against another reference or catalog, fit with \
         those, which regenerates nothing; if not, this build chooses census positions \
         differently from the one that wrote them, so regenerate them with \
         {THE_COMMAND_THAT_REBUILDS_A_CENSUS} and fit again"
    )]
    AnotherSelection,

    /// The cohort holds no recording terms to check the selection against, which means it holds
    /// no samples.
    #[error("a fit needs at least one sample")]
    NoSamples,

    /// The generic half failed.
    #[error("fitting the ordinary positions")]
    Generic {
        /// What the estimator said.
        #[source]
        source: Box<JointFitError>,
    },

    /// A census section could not be read while the tracts were gathered.
    #[error("reading a sample's repeat-tract evidence")]
    Tracts {
        /// What the reader said.
        #[source]
        source: Box<CensusError>,
    },
}

/// **Refuse a cohort whose censuses were written against another selection than `loci`**: the
/// digest every census records of the ordinary positions it kept, against the digest of `loci`'s.
///
/// # Errors
///
/// [`CohortFitError::AnotherSelection`] when they differ; [`CohortFitError::NoSamples`] for a cohort
/// of none.
pub fn refuse_another_selection(
    cohort: &CohortCensusEvidence,
    loci: &CensusLoci,
) -> Result<(), CohortFitError> {
    // The digest is over the kept ordinary positions, in the order the writer digested them —
    // which is the order `CensusWriter` holds them in, straight from the selection.
    let recorded = cohort
        .terms()
        .ok_or(CohortFitError::NoSamples)?
        .kept_loci
        .clone();
    let mut digester = CensusLociDigester::new();
    for (index, position) in loci.generic().iter().enumerate() {
        digester.observe(index, *position);
    }
    if digester.finish() != recorded {
        return Err(CohortFitError::AnotherSelection);
    }
    Ok(())
}

/// **A cohort's repeat-tract evidence**: every stratum with each sample's reads at its tracts, and
/// how many tracts the selection kept.
pub struct CohortTractStrata {
    /// One entry a stratum, its samples indexed by the cohort's own sample order.
    pub strata: Vec<StratumEvidence>,
    /// The repeat tracts the selection kept, over every stratum.
    pub tracts: usize,
}

/// **Every repeat-tract stratum of a cohort, with each sample's reads at its tracts** — what the
/// tract half of [`fit_a_cohort`] fits.
///
/// # Errors
///
/// [`CohortFitError::Tracts`] when a census's tract sections cannot be read.
pub fn the_tract_strata_of_a_cohort(
    cohort: &mut CohortCensusEvidence,
    loci: &CensusLoci,
    contig_of: &dyn Fn(&str) -> Option<ContigId>,
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
) -> Result<CohortTractStrata, CohortFitError> {
    let kept = strata_of_kept_loci(loci, contig_of);
    let strata = gather_strata(cohort, &kept, slippage_group_of).map_err(|source| {
        CohortFitError::Tracts {
            source: Box::new(source),
        }
    })?;
    Ok(CohortTractStrata {
        strata,
        tracts: kept.len(),
    })
}

/// **Fit a cohort of censuses, both halves.**
///
/// `loci` is the selection rebuilt from this run's reference and catalog; it is checked against
/// the digest every census carries before a tract is read. `contig_of` turns a contig's name
/// into the identifier the records use, and `slippage_group_of` says which slippage group each
/// read group is in — a declaration of the run's, never estimated.
///
/// # Errors
///
/// [`CohortFitError::AnotherSelection`] when the rebuilt loci are not the ones the censuses were
/// written against, and the two halves' own failures otherwise.
pub fn fit_a_cohort(
    cohort: &mut CohortCensusEvidence,
    loci: &CensusLoci,
    contig_of: &dyn Fn(&str) -> Option<ContigId>,
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
    generic: &JointFitConfig,
    tracts: &SsrFitConfig,
) -> Result<CohortFit, CohortFitError> {
    // **Checked before anything is fitted**, because the selection decides how every tract record
    // is indexed and a wrong one produces numbers rather than a failure.
    //
    // The digest is over the kept ordinary positions, in the order the writer digested them —
    // which is the order `CensusWriter` holds them in, straight from the selection.
    refuse_another_selection(cohort, loci)?;

    let fit = crate::parameter_estimation::joint::fit::fit_jointly(cohort, generic).map_err(
        |source| CohortFitError::Generic {
            source: Box::new(source),
        },
    )?;

    // **The one number the tract half takes from the generic one**, per sample and in the
    // cohort's own sample order.
    let homozygote_excess: Vec<f64> = cohort
        .sample_names()
        .map(|name| {
            fit.hom_excess
                .get(name)
                .map_or(0.0, |estimate| estimate.value.get())
        })
        .collect();

    let CohortTractStrata {
        strata: evidence,
        tracts: kept_tracts,
    } = the_tract_strata_of_a_cohort(cohort, loci, contig_of, slippage_group_of)?;
    // **A large cohort's strata are each read from a subset of its samples**, the first in a
    // fixed order drawn from the samples' names (`fit_precision.md` §4.4); the order is indexed as
    // the evidence's samples and `homozygote_excess` are, the cohort's own order.
    let names: Vec<&str> = cohort.sample_names().collect();
    let order = sample_order(&names, SAMPLE_ORDER_SEED);
    let outcomes =
        ssr_fit::fit_strata_on_sample_subsets(&evidence, &homozygote_excess, &order, tracts);
    // **The evidence is the largest thing the tract half holds and nothing after the fit needs
    // more of it than these two counts a stratum**, so it goes here rather than when the
    // parameters file has been written.
    let substitution_counts = evidence
        .iter()
        .map(|evidence| evidence.substitution_counts())
        .collect();
    drop(evidence);

    Ok(CohortFit {
        generic: fit,
        strata: outcomes,
        tracts: kept_tracts,
        substitution_counts,
    })
}

/// **Every read group in one slippage group**, which is what a cohort too thin to fit them apart
/// can afford.
///
/// One slippage group per read group is the specified grain; pooling is a run's own declaration
/// and is recorded as such, never estimated
/// (`doc/devel/ng/spec/parameters_file.md` — a read group's slippage group is the run's
/// declaration and is not estimated at all).
#[must_use]
pub fn every_read_group_pooled(cohort: &CohortCensusEvidence) -> BTreeMap<ReadGroupId, u32> {
    cohort
        .read_groups()
        .iter()
        .map(|group| (*group, 0))
        .collect()
}

// ---------------------------------------------------------------------
// Turning the fit into the numbers a calling run scores with
// ---------------------------------------------------------------------

/// **The numbers a fit produced, in the shape a calling run reads.**
///
/// `RunParameters::assemble` takes nine groups of numbers, and this is where a cohort of
/// censuses supplies them. **It gathers rather than fits**: every value it hands on comes from
/// [`fit_a_cohort`] unchanged, except the genotype prior's seed, which
/// [`RunParameters::seed_from_moments`] solves in closed form from the fit's two moments.
///
/// # The base-quality calibration, and where its two halves come from
///
/// A library's calibration is fitted from two numbers together: the error rate the run measured,
/// and what the library's own base qualities claimed — Σ over reads of `ln P(this read is
/// wrong)`, with the count it ran over. **`assemble` refuses one without the other**, and says
/// why: they come from one pass over one set of reads.
///
/// The rate comes from this fit. The totals come from the census, which accumulates them as the
/// loci go past — the same accumulator the per-sample calibration pre-pass uses, with its unit
/// unchanged. **A read group whose census accumulated nothing is left out of both maps**, so the
/// pair stays whole and that library takes the defaulted calibration rather than a fitted rate
/// with no evidence behind it.
///
/// `inbreeding` is one coefficient a sample **in the run's own sample order**, already resolved
/// against the ladder in
/// [`DeclaredInbreeding::of_each_sample_over`](crate::calling::parameters_file::DeclaredInbreeding::of_each_sample_over)
/// — so it may be what the operator stated, what this fit measured, or the default, and it carries
/// which. **It arrives as the warranted estimates rather than as bare values so that the numbers
/// this assembles and the numbers the parameters file writes cannot be two different lists**: an
/// earlier version took the values here and the estimates in `parameters_file_of`, which let a
/// caller hand the run one coefficient and the file another with nothing to notice.
#[must_use]
pub fn parameters_from_the_fit(
    fit: &CohortFit,
    cohort: &CohortCensusEvidence,
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
    inbreeding: &[Estimate<InbreedingF>],
    ploidy: Ploidy,
) -> RunParameters {
    // **The clean class's rate is the sequencing error rate.** The fit models two: how often a
    // read misreads a base at an ordinary position, and how often a read disagrees with the
    // reference at a mismapped one. Only the first is chemistry.
    let error_rate_by_read_group = fit.generic.sequencing_error_rates();

    // Per sample, per read group — and a read group belongs to one sample, so flattening cannot
    // collide.
    let contamination_by_read_group: BTreeMap<ReadGroupId, ContaminationEstimate> = fit
        .generic
        .contamination
        .values()
        .flat_map(|of_sample| of_sample.iter().cloned())
        .collect();

    // **One rate per (read group, stratum, ploidy)**, from the counts the evidence measured
    // rather than out of the fit: a stratum whose fit was refused still measured a substitution
    // rate, and a run that dropped it would score its tracts against nothing.
    let mut ssr_substitution_rate: BTreeMap<StratumKey, Estimate<ErrorRate>> = BTreeMap::new();
    for counts in &fit.substitution_counts {
        let Some(rate) = counts.substitution_rate() else {
            continue;
        };
        let Ok(rate) = ErrorRate::try_new(rate) else {
            continue;
        };
        // **The two `Stratum` types are different and the conversion is here.** The census
        // names a stratum by a period in bases and a reference repeat count as plain numbers;
        // the calling key names it by the checked types. A stratum the census holds but the
        // checked types reject is skipped rather than coerced.
        let (Ok(period), Ok(repeats)) = (
            SsrPeriod::try_new(counts.stratum.period as usize),
            u32::try_from(counts.stratum.reference_repeats),
        ) else {
            continue;
        };
        let repeats = RepeatCount(repeats);
        // **A counted rate's error is the binomial one**, √(p(1 − p)/n) over the bases compared:
        // each base compared either matched or did not. **None where no base mismatched, or
        // every one did**: there the formula gives zero, which would claim the rate known exactly
        // from a count that found only one outcome.
        let standard_error = (counts.mismatching_bases > 0
            && counts.mismatching_bases < counts.bases_compared)
            .then(|| (rate.get() * (1.0 - rate.get()) / counts.bases_compared as f64).sqrt());
        for group in cohort.read_groups() {
            ssr_substitution_rate.insert(
                StratumKey {
                    read_group: *group,
                    stratum: SsrStratum::new(period, repeats),
                    ploidy,
                },
                Estimate {
                    value: rate,
                    provenance: Provenance::FittedHere,
                    observations: counts.bases_compared,
                    standard_error,
                },
            );
        }
    }

    // **What each library's own base qualities claimed**, which every census carries beside who
    // its read groups are. A calibration is fitted from this and the measured rate together;
    // `assemble` refuses one without the other, and says why.
    //
    // **A read group nothing scored gets empty totals**, which is what they are: the census sums
    // them at ordinary positions only (`calibration::minted_error_by_read_group`), so a sample whose
    // walk wrote repeat-tract sections and no ordinary one has none, while the fit still lists its
    // read group with a defaulted rate from no observation. Empty totals have no mean, and
    // `ReadGroupCalibration::from_fitted_rate` gives such a library the defaulted multiplier of
    // one. Until 2026-10-06 a cohort holding such a sample stopped here with a panic. **Only such a
    // read group**: a fitted rate with no totals behind it is a lost read-group axis, and
    // `RunParameters::assemble` still refuses it.
    let census_totals: BTreeMap<ReadGroupId, MintedReadErrors> = cohort
        .samples()
        .iter()
        .flat_map(|sample| {
            sample
                .minted_read_errors()
                .iter()
                .map(|(group, totals)| (*group, *totals))
        })
        .collect();
    let minted_by_read_group: BTreeMap<ReadGroupId, MintedReadErrors> = error_rate_by_read_group
        .iter()
        .filter_map(|(group, rate)| match census_totals.get(group) {
            Some(totals) => Some((*group, *totals)),
            None if rate.provenance == Provenance::Defaulted && rate.observations == 0 => {
                Some((*group, MintedReadErrors::default()))
            }
            None => None,
        })
        .collect();

    RunParameters::assemble(
        &error_rate_by_read_group,
        &minted_by_read_group,
        &contamination_by_read_group,
        SequencingBatches::all_together_over(cohort.read_groups().len(), cohort.len()),
        inbreeding.iter().map(|estimate| estimate.value).collect(),
        RunParameters::seed_from_moments(
            fit.generic.fitted_alternative_frequency(),
            fit.generic.fitted_diversity(),
        ),
        StratumFits::over(&fit.strata, slippage_group_of.clone()),
        ssr_substitution_rate,
        ploidy,
    )
}

/// **The inbreeding coefficients this fit measured, by sample name** — the middle rung of
/// [`DeclaredInbreeding::of_each_sample_over`](crate::calling::parameters_file::DeclaredInbreeding::of_each_sample_over)'s
/// ladder.
///
/// The quantity is each sample's **homozygote excess**: how much less heterozygous it is than the
/// allele frequencies this same fit produced predict. **That is circular and the output says so**
/// — `joint::census_moments` carries the warning — and it is the only fitted source there is since
/// the runs-of-homozygosity estimator was removed
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
/// reported where it stopped looking rather than what it measured. *The shipped maximisation cannot
/// return exactly one — it is a golden-section search that reports the midpoint of a bracket
/// strictly inside `[0, 1]` — so this guards the type's contract rather than an observed case.*
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

/// **The parameters file a fit over censuses writes.**
///
/// Everything the file needs that is not a fitted number comes from the cohort itself: the
/// read-group table is built from what the censuses declare
/// ([`read_groups_of`](super::census_cohort::read_groups_of)), and the census identity is the
/// recording terms every one of them agreed on — so a calling run handed this file can tell
/// whether its own evidence was recorded the same way.
///
/// **The calibration rows carry no `observations`.** Their warrant and standard error are the
/// calibration's own, but the count this fit keeps beside a library's rate is the census positions
/// it read, and the file counts a multiplier's evidence in reads: writing it would put a position
/// count under the name of a read count. `ReadsBehindEachCalibration::no_count_in_reads` is what
/// states no count.
///
/// `snp_indel_fit_starts` is the fit's record of how each of its starts ended, written into the
/// file's `[fitted_from]` (`doc/devel/ng/spec/fit_precision.md` §5.2).
#[must_use]
pub fn parameters_file_of(
    parameters: &RunParameters,
    read_groups: &ReadGroups,
    inbreeding: &[Estimate<InbreedingF>],
    reference: &ReferenceDigest,
    terms: &RecordingTerms,
    repeat_routing: &StrRepeatCriteria,
    snp_indel_fit_starts: &[StartRecord],
) -> ParametersFile {
    ParametersFile::of_run(
        parameters,
        read_groups,
        &ReadsBehindEachCalibration::no_count_in_reads(read_groups.len()),
        inbreeding,
        reference,
        CensusIdentity::of(terms),
        repeat_routing,
    )
    .with_snp_indel_fit_starts(Some(the_starts_as_written(snp_indel_fit_starts)))
}

/// **How each start of the SNP/indel fit ended, in the file's words**: its number, how it ended,
/// its passes, and for one that agreed, which earlier start.
fn the_starts_as_written(starts: &[StartRecord]) -> Vec<SnpIndelFitStart> {
    starts
        .iter()
        .map(|record| {
            let (ended, agreed_with_start) = match record.ended {
                StartEnding::Converged => (StartOutcome::Converged, None),
                StartEnding::AtTheLimit => (StartOutcome::AtThePassLimit, None),
                StartEnding::Agreed { with_start } => (
                    StartOutcome::AgreedWithAnEarlierStart,
                    Some(u32::try_from(with_start).expect("a fit's starts are few")),
                ),
            };
            SnpIndelFitStart {
                start: u32::try_from(record.number).expect("a fit's starts are few"),
                ended,
                passes: record.passes,
                agreed_with_start,
            }
        })
        .collect()
}

#[cfg(test)]
mod fitted_inbreeding_tests {
    //! **What the fit's homozygote excess becomes when it is offered as a coefficient.**

    use super::*;

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

#[cfg(test)]
mod tests {
    //! **`parameter_prepass_runs.md` plan step C4**: a cohort of censuses is fitted, both halves,
    //! and the selection it is fitted against is the one it was written against.

    use super::*;
    use crate::cli::generate_psps::{GeneratePspsArgs, run_generate_psps};
    use crate::cli::test_fixtures::a_varying_cohort_on_disk;
    use crate::run::census_cohort::every_census_in_the_cohorts_psps;
    use crate::run::psp_caller::OpenPspCohort;
    use crate::run::test_fixtures::a_census_plan_over_selecting;
    use std::path::PathBuf;

    /// **The varying fixture cohort**, walked, each psp carrying its own census.
    ///
    /// **Not the plain on-disk cohort**: that one's reference is all `A`, so every base is a
    /// homopolymer, the whole genome routes to the repeat path, and the selection keeps no tract
    /// at all — measured, 0 strata over 0 tracts, which would leave the repeat-tract half of the
    /// fit untested. This one carries a deliberate ten-copy `GT` tract.
    pub(super) fn a_fitted_cohorts_inputs()
    -> (crate::cli::test_fixtures::AVaryingCohort, Vec<PathBuf>) {
        use crate::region_typing::DEFAULT_MAX_STR_LEN;
        use crate::region_typing::segment_criteria::{
            DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
        };

        let cohort = a_varying_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        let walk = GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: false,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
            max_reads_per_position: crate::locus_generation::pileup::DEFAULT_MAX_SNP_COLUMN_DEPTH,
        };
        run_generate_psps(&walk).expect("the cohort walks into psps");
        let mut walked: Vec<PathBuf> = std::fs::read_dir(&psps)
            .expect("the walk made the directory")
            .map(|entry| entry.expect("an entry").path())
            .filter(|path| path.extension().is_some_and(|it| it == "psp"))
            .collect();
        walked.sort();
        assert_eq!(
            walked.len(),
            2,
            "one psp a sample, each carrying its census"
        );
        (cohort, walked)
    }

    /// The cohort's censuses, read out of its psps whole, where the command reads them a half at a time so that it can judge each census in between.
    ///
    /// **The cohort is dropped here and the fit below still reads**, because a census read this
    /// way holds a path and an offset rather than the psp's descriptor
    /// (`SampleCensusEvidence::backed`). The command keeps its cohort alive because it wants the
    /// paths and the settings, not because the evidence needs it.
    pub(super) fn the_censuses_in(psps: &[PathBuf]) -> CohortCensusEvidence {
        let cohort = OpenPspCohort::open(psps).expect("the walk wrote one cohort");
        every_census_in_the_cohorts_psps(&cohort).expect("each psp carries its census")
    }

    fn a_generic_config() -> JointFitConfig {
        JointFitConfig::default()
    }

    /// **Both halves run over a cohort read from census files, and both have something to
    /// read.**
    ///
    /// The numbers themselves are not asserted: two samples over 600 bases is not a population,
    /// and a figure from it would be a fact about the fixture. What is asserted is that the
    /// evidence reached the estimator — the tract half over the fixture's own tract, the generic
    /// half over its positions — and that the same cohort fitted twice gives the same answers,
    /// which is the property a run depends on and the one an unstable fit would break.
    #[test]
    fn a_cohort_of_censuses_is_fitted_both_halves() {
        let (cohort, psps) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = the_censuses_in(&psps);

        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open);

        let fitted = fit_a_cohort(
            &mut open,
            &plan.loci,
            &contig_of,
            &pooled,
            &a_generic_config(),
            &SsrFitConfig::default(),
        );

        let fit = fitted.expect("a cohort of two censuses fits");
        assert!(
            fit.tracts > 0,
            "the repeat-tract half read no tract, so this fixture tests only half the fit",
        );
        assert!(
            !fit.strata.is_empty(),
            "a stratum's outcome a stratum, even where it was refused for want of tracts",
        );
        assert!(
            fit.generic.noisy_share.is_finite(),
            "the generic half returned a number rather than a NaN: {}",
            fit.generic.noisy_share,
        );

        // The same cohort, fitted again from the same files.
        let mut again = the_censuses_in(&psps);
        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&again);
        let twice = fit_a_cohort(
            &mut again,
            &plan.loci,
            &contig_of,
            &pooled,
            &a_generic_config(),
            &SsrFitConfig::default(),
        )
        .expect("it fits the second time too");

        assert_eq!(
            twice.generic.noisy_share, fit.generic.noisy_share,
            "one cohort fitted twice gives one answer",
        );
        assert_eq!(twice.tracts, fit.tracts);

        // **What outlives the evidence**: one entry a stratum, in the fit's order, and at least
        // one stratum whose reads were compared at all — otherwise the parameters file's
        // substitution-rate loop has nothing to read on this fixture.
        assert_eq!(
            fit.substitution_counts
                .iter()
                .map(|counts| counts.stratum)
                .collect::<Vec<_>>(),
            fit.strata
                .iter()
                .map(StratumOutcome::stratum)
                .collect::<Vec<_>>(),
            "one substitution count a stratum the tract half gathered, in its order",
        );
        assert!(
            fit.substitution_counts
                .iter()
                .any(|counts| counts.bases_compared > 0),
            "the fixture's tracts compared no base, so no substitution rate is tested end to end",
        );
        assert_eq!(twice.substitution_counts, fit.substitution_counts);
    }

    /// **A cohort fitted against another selection is refused before a tract is read.**
    ///
    /// A census stores a tract by its index within its stratum, so a selection rebuilt from
    /// another seed indexes one stratum's tracts as another's — numbers rather than a failure.
    #[test]
    fn a_cohort_fitted_against_another_selection_is_refused() {
        let (cohort, psps) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = the_censuses_in(&psps);

        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open);

        // A selection holding one position the census's does not is enough: the digest is over
        // the kept positions in order.
        let shifted = crate::parameter_estimation::joint::loci::CensusLoci::from_parts(
            plan.loci
                .generic()
                .iter()
                .skip(1)
                .copied()
                .collect::<Vec<_>>(),
            plan.loci.ssr().clone(),
            plan.loci.ssr_stratum_counts().clone(),
        );

        let error = fit_a_cohort(
            &mut open,
            &shifted,
            &contig_of,
            &pooled,
            &a_generic_config(),
            &SsrFitConfig::default(),
        )
        .expect_err("these censuses were written against another selection");

        assert!(
            matches!(error, CohortFitError::AnotherSelection),
            "{error:?}"
        );
        // **And it says what to do about it** (spec §4.2, third row). No command reaches this
        // message: `estimate-parameters` compares the recorded settings itself and reports them as
        // rows of its regeneration report (plan step C5). What it serves is a caller that fits
        // without comparing, which is why it still offers both ways out — such a caller has ruled
        // out neither.
        let said = error.to_string();
        let free = said
            .find("if the psps were walked against another reference or catalog, fit with those")
            .unwrap_or_else(|| panic!("the fix that regenerates nothing is offered: {said}"));
        let costly = said
            .find(&format!(
                "so regenerate them with {THE_COMMAND_THAT_REBUILDS_A_CENSUS} and fit again"
            ))
            .unwrap_or_else(|| panic!("and the command, as an instruction: {said}"));
        // **The free check first.** A wrong reference ends here as surely as a changed build, and
        // a person who regenerated in that case would spend a quarter of an hour a sample and be
        // refused again by the next fit.
        assert!(
            free < costly,
            "the fix that regenerates nothing comes before the one that costs hours: {said}",
        );
    }
}

#[cfg(test)]
mod writing_the_parameters_file {
    //! **Plan step C5**: the fit's numbers, assembled and written as the file a calling run
    //! takes.
    //!
    //! # What made this run at all
    //!
    //! `RunParameters::assemble` refuses a read group that has a fitted error rate and no minted
    //! read-error total, and says why: the two come from one pass over one set of reads. The plan
    //! had assumed the pair would fall back to a defaulted calibration; it does not, it panics,
    //! and these tests found that on their first run. **The owner's answer, 2026-09-05, was to
    //! bring Milestone E forward**: the census now accumulates the totals as its loci go past, so
    //! both halves are present and the calibration is fitted rather than defaulted.

    use super::tests::the_censuses_in;
    use super::*;
    use crate::calling::parameters_file::Warrant;
    use crate::run::census_cohort::read_groups_of;
    use crate::run::census_fit::tests::a_fitted_cohorts_inputs;

    /// **The first parameters file this tree produces from data**, and what it may and may not
    /// claim.
    ///
    /// The population's numbers are fitted. The base-quality calibration is not, and the file
    /// says `defaulted` against it — because a read group's calibration needs its minted
    /// read-error total as well as its rate, and no part of a census carries one.
    #[test]
    fn the_fit_writes_a_parameters_file_that_says_what_it_fitted() {
        let (parameters, file) = a_fitted_cohorts_parameters();

        let toml = file.to_toml();
        assert!(
            toml.contains("read_groups"),
            "the file names the run's read groups: {toml:.400}",
        );
        assert!(
            toml.contains("fitted from reads"),
            "and says how many of its groups of numbers were fitted at all: {toml:.400}",
        );
        for row in &file.fitted_from.read_groups {
            assert!(
                !row.declared_id.is_empty(),
                "every read group carries the @RG ID it was declared under, which is what a \
                 calling run checks the file by",
            );
        }
        assert_eq!(
            file.fitted_from.read_groups.len(),
            parameters.read_group_count(),
            "one row a library",
        );
    }

    /// **The base-quality calibration is fitted, and the file says which library's it is.**
    ///
    /// Both halves are present: the rate from this fit, the minted totals from the census. A
    /// library whose census accumulated nothing is left out of both maps and takes the defaulted
    /// calibration — never a fitted rate with no evidence behind it.
    #[test]
    fn the_base_quality_calibration_says_where_it_came_from() {
        let (_parameters, file) = a_fitted_cohorts_parameters();

        assert!(!file.base_quality_calibration.by_read_group.is_empty());
        for row in &file.base_quality_calibration.by_read_group {
            assert!(
                matches!(
                    row.error_probability_multiplier.warrant,
                    Warrant::FittedHere | Warrant::Defaulted,
                ),
                "a calibration is either fitted from this library's own evidence or honestly \
                 defaulted, never anything else: {:?}",
                row.error_probability_multiplier.warrant,
            );
        }
        assert!(
            file.base_quality_calibration
                .by_read_group
                .iter()
                .any(|row| row.error_probability_multiplier.warrant == Warrant::FittedHere),
            "at least one library of the fixture has both halves, so at least one calibration \
             is fitted — otherwise this route still cannot produce one and the change did \
             nothing",
        );
    }

    /// **The file names the census these numbers were fitted from**, term by term — which is
    /// what lets a calling run tell whether its own evidence was recorded the same way.
    #[test]
    fn the_file_names_the_census_it_was_fitted_from() {
        let (_parameters, file) = a_fitted_cohorts_parameters();

        assert!(
            !file.fitted_from.census.terms.is_empty(),
            "a run that fitted from a census names its terms, unlike a defaults or direct run",
        );
    }

    /// **A stratum that compared no base gets no substitution rate, and one that did gets its
    /// own**, for every read group, counted over the bases it compared.
    ///
    /// The fixture's one stratum compares 460 bases and finds none disagreeing, so on its own it
    /// reaches neither case: its rate is zero and no stratum lacks one. Three strata are added to
    /// the fit's counts before the file is assembled — one with nothing compared, which must be
    /// skipped rather than written as a fitted zero, one with 3 disagreeing in 4,000, and one whose
    /// 5 bases all disagree. The fixture's own zero rate and the third stratum's rate of one carry
    /// no standard error; the stratum with 3 in 4,000 carries its binomial error.
    #[test]
    fn a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own() {
        use crate::parameter_estimation::joint::census::Stratum as CensusStratum;

        let nothing = CensusStratum {
            period: 3,
            reference_repeats: 7,
        };
        let some = CensusStratum {
            period: 4,
            reference_repeats: 6,
        };
        let all = CensusStratum {
            period: 5,
            reference_repeats: 5,
        };
        let (parameters, _) = a_fitted_cohorts_parameters_with(&[
            StratumSubstitutionCounts {
                stratum: nothing,
                bases_compared: 0,
                mismatching_bases: 0,
            },
            StratumSubstitutionCounts {
                stratum: some,
                bases_compared: 4_000,
                mismatching_bases: 3,
            },
            StratumSubstitutionCounts {
                stratum: all,
                bases_compared: 5,
                mismatching_bases: 5,
            },
        ]);
        let rates: Vec<_> = parameters.ssr_substitution_rate().collect();
        let of_period = |period: u8| {
            rates
                .iter()
                .filter(|(key, _)| key.stratum.period.get() == period)
                .collect::<Vec<_>>()
        };

        assert!(
            of_period(3).is_empty(),
            "a stratum that compared no base has no rate, not a fitted zero"
        );
        let fixtures_own = of_period(2);
        assert!(
            !fixtures_own.is_empty(),
            "the fixture's own stratum is written, which is what the count below is measured by"
        );
        // **A count that found no mismatch gives a rate of zero and no error**, not an error of
        // zero, which would claim the rate known exactly.
        for (_, estimate) in &fixtures_own {
            assert_eq!(estimate.value.get(), 0.0);
            assert!(estimate.observations > 0);
            assert_eq!(estimate.standard_error, None);
        }
        // Nor does one where every base compared mismatched.
        let every_base = of_period(5);
        assert_eq!(every_base.len(), fixtures_own.len());
        for (_, estimate) in every_base {
            assert_eq!(estimate.value.get(), 1.0);
            assert_eq!(estimate.standard_error, None);
        }
        let added = of_period(4);
        assert_eq!(
            added.len(),
            fixtures_own.len(),
            "one rate a read group, as the fixture's own stratum has"
        );
        for (key, estimate) in added {
            assert_eq!(key.stratum.repeats, RepeatCount(6));
            assert_eq!(estimate.value.get(), 3.0 / 4_000.0);
            assert_eq!(estimate.observations, 4_000);
            assert_eq!(estimate.provenance, Provenance::FittedHere);
            // **The binomial error of 3 mismatches in 4,000 bases**, √(p(1 − p)/n): 0.000433,
            // about the rate over √3 — one over the square root of the mismatches, as a count's
            // error is.
            let rate: f64 = 3.0 / 4_000.0;
            let error = estimate
                .standard_error
                .expect("a counted rate has an error");
            assert_eq!(error, (rate * (1.0 - rate) / 4_000.0).sqrt());
            assert!((error - 0.000_433).abs() < 1e-6, "{error}");
        }
    }

    /// Fit the fixture cohort and assemble its parameters file.
    fn a_fitted_cohorts_parameters() -> (RunParameters, ParametersFile) {
        a_fitted_cohorts_parameters_with(&[])
    }

    /// Fit the fixture cohort, add `extra` to the substitution counts the fit measured, and
    /// assemble its parameters file.
    fn a_fitted_cohorts_parameters_with(
        extra: &[StratumSubstitutionCounts],
    ) -> (RunParameters, ParametersFile) {
        use crate::types::InbreedingF;

        let (cohort, psps) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = crate::run::test_fixtures::a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = the_censuses_in(&psps);
        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open);
        let mut fit = fit_a_cohort(
            &mut open,
            &plan.loci,
            &contig_of,
            &pooled,
            &JointFitConfig::default(),
            &SsrFitConfig::default(),
        )
        .expect("the cohort fits");
        fit.substitution_counts.extend_from_slice(extra);

        let read_groups = read_groups_of(&open, &psps);
        // **Each sample's read groups name that sample's psp** — the column exists only so a
        // message can name a file, and it is indexed by position, so a table built against
        // another order would name the wrong file in every message about it.
        assert_eq!(
            read_groups
                .get(read_groups.read_groups_per_sample()[0].read_groups[0])
                .file
                .as_ref(),
            psps[0].as_path(),
            "the first sample's read groups name the first sample's psp",
        );
        // **A stated coefficient, which is what this fixture is about** — the resolution between
        // stated and fitted is `DeclaredInbreeding`'s and is tested there.
        let declared = InbreedingF::try_new(0.0).expect("zero is a coefficient");
        let stated: Vec<Estimate<InbreedingF>> = (0..open.len())
            .map(|_| Estimate {
                value: declared,
                provenance: Provenance::Supplied,
                observations: 0,
                standard_error: None,
            })
            .collect();

        let parameters = parameters_from_the_fit(
            &fit,
            &open,
            &pooled,
            &stated,
            crate::types::Ploidy::try_new(2).expect("diploid"),
        );
        let terms = open
            .terms()
            .expect("a cohort of one or more samples records terms")
            .clone();
        let file = parameters_file_of(
            &parameters,
            &read_groups,
            &stated,
            &plan.terms.reference,
            &terms,
            &plan.terms.ssr_criteria,
            &fit.generic.starts,
        );
        (parameters, file)
    }
}

#[cfg(test)]
mod a_sample_with_repeat_tracts_only {
    //! **A sample whose census holds repeat tracts and no ordinary position, through to the file**
    //! (plan `fit_precision.md` step E2, the owner at checkpoint A′).

    use super::*;
    use crate::calling::likelihood::ReadGroupCalibration;
    use crate::calling::parameters_file::{DeclaredInbreeding, Warrant};
    use crate::parameter_estimation::joint::census::{
        NamedReadGroup, SampleCensusEvidence, Section, SectionKey, SsrEvidence,
        Stratum as CensusStratum,
    };
    use crate::parameter_estimation::joint::fit::bench_fixtures::{as_cohort, draw_cohort};
    use crate::parameter_estimation::joint::fit::{FrequencyDensity, StartingPoint, fit_jointly};
    use crate::repeat_catalog::StrRepeatCriteria;

    /// **Its read group and its sample are written `defaulted`, from nothing, with no error — and
    /// the run gets that far.** Four drawn samples and a fifth with one repeat-tract section. The
    /// census sums what a library's qualities claimed at ordinary positions only, so the fifth
    /// read group has no totals, while the drawn four are given some. Before 2026-10-06 the run
    /// stopped assembling at that read group with a panic; now it takes the defaulted multiplier of
    /// one, beside the four fitted ones. Its sample's inbreeding coefficient is the default,
    /// `defaulted`, where it was a fitted number no read informed. The drawn samples keep fitted
    /// coefficients with errors.
    #[test]
    fn its_read_group_and_its_coefficient_are_written_defaulted() {
        let mut drawn = draw_cohort(
            4,
            2_000,
            6.0,
            (0.003, 0.06, 0.02),
            FrequencyDensity {
                p_invariant: 0.90,
                p_fixed_alt: 0.01,
                a: 0.7,
                b: 2.5,
            },
            0.2,
            0x5EED_0E20_0000_0002,
        );
        // **Real claimed-error totals for the four drawn read groups** — a thousand reads at Q30
        // each — so that their rates make fitted multipliers beside the fifth's defaulted one. Drawn
        // samples are built without any.
        drawn.samples = std::mem::take(&mut drawn.samples)
            .into_iter()
            .enumerate()
            .map(|(s, sample)| {
                let group = ReadGroupId(u32::try_from(s).expect("four samples"));
                sample.with_minted_read_errors(BTreeMap::from([(
                    group,
                    MintedReadErrors::of_observation(-6.907_755 * 1_000.0, 1_000),
                )]))
            })
            .collect();
        let terms = as_cohort(&drawn.samples)
            .terms()
            .expect("a drawn cohort records terms")
            .clone();
        drawn.samples.push(SampleCensusEvidence::resident(
            "tracts_only".to_string(),
            terms.clone(),
            NamedReadGroup::drawn_for("tracts_only", [ReadGroupId(4)]),
            BTreeMap::new(),
            BTreeMap::from([(
                SectionKey::Ssr(
                    ReadGroupId(4),
                    CensusStratum {
                        period: 2,
                        reference_repeats: 6,
                    },
                ),
                Section::Ssr(SsrEvidence::never_walked(0)),
            )]),
        ));
        let mut cohort = as_cohort(&drawn.samples);
        let generic = fit_jointly(
            &mut cohort,
            &JointFitConfig {
                quadrature_nodes: 8,
                max_passes: 40,
                duplicated_positions: false,
                estimate_contamination: false,
                starting_points: vec![StartingPoint::spanning_the_class_separation()[1]],
                ..JointFitConfig::default()
            },
        )
        .expect("four drawn samples identify the cohort's numbers");
        let fit = CohortFit {
            generic,
            strata: Vec::new(),
            tracts: 0,
            substitution_counts: Vec::new(),
        };
        let read_groups = ReadGroups::of_lanes(&[
            ("rg0", "s0", "lib0"),
            ("rg1", "s1", "lib1"),
            ("rg2", "s2", "lib2"),
            ("rg3", "s3", "lib3"),
            ("rg4", "tracts_only", "lib4"),
        ]);
        let inbreeding = DeclaredInbreeding::nothing_said()
            .of_each_sample_over(&read_groups, &fitted_inbreeding_of(&fit.generic.hom_excess));
        let pooled: BTreeMap<ReadGroupId, u32> =
            (0..5).map(|group| (ReadGroupId(group), 0)).collect();
        let parameters = parameters_from_the_fit(
            &fit,
            &cohort,
            &pooled,
            &inbreeding,
            crate::types::Ploidy::try_new(2).expect("diploid"),
        );
        assert_eq!(
            parameters.calibration_by_read_group()[4],
            ReadGroupCalibration::defaulted()
        );
        for group in 0..4 {
            let fitted = parameters.calibration_by_read_group()[group];
            assert_eq!(
                fitted.provenance,
                Provenance::FittedHere,
                "read group {group}"
            );
            assert!(fitted.scale_standard_error.is_some(), "read group {group}");
        }

        let file = parameters_file_of(
            &parameters,
            &read_groups,
            &inbreeding,
            &ReferenceDigest([7; 16]),
            &terms,
            &StrRepeatCriteria::default(),
            &fit.generic.starts,
        );
        let multiplier =
            &file.base_quality_calibration.by_read_group[4].error_probability_multiplier;
        assert_eq!(
            (
                multiplier.value,
                multiplier.warrant,
                multiplier.observations,
                multiplier.standard_error
            ),
            (1.0, Warrant::Defaulted, None, None)
        );
        let coefficient = |sample: &str| {
            file.inbreeding
                .by_sample
                .iter()
                .find(|row| row.sample == sample)
                .expect("every sample has a row")
                .inbreeding_coefficient
        };
        let silent = coefficient("tracts_only");
        assert_eq!(
            (
                silent.value,
                silent.warrant,
                silent.observations,
                silent.standard_error
            ),
            (0.0, Warrant::Defaulted, None, None)
        );
        for sample in ["s0", "s1", "s2", "s3"] {
            let fitted = coefficient(sample);
            assert_eq!(fitted.warrant, Warrant::FittedHere, "{sample}");
            assert!(fitted.standard_error.is_some(), "{sample} has an error");
        }
        let read_back = ParametersFile::from_toml(&file.to_toml()).expect("the file reads back");
        assert_eq!(read_back, file);
        read_back.validate().expect("and passes its own checks");
    }
}

#[cfg(test)]
mod the_starts_in_the_file {
    use super::*;

    /// **Each way a start can end is written as the file spells it**, with the start it agreed
    /// with, numbered as the fit numbers its starts.
    #[test]
    fn each_ending_is_written_with_its_number_and_passes() {
        let record = |number: usize, ended: StartEnding, passes: u32| StartRecord {
            number,
            ended,
            passes,
            log_likelihood: -1.0,
            furthest_from_settled: None,
        };
        let written = the_starts_as_written(&[
            record(1, StartEnding::Converged, 41),
            record(2, StartEnding::AtTheLimit, 200),
            record(3, StartEnding::Agreed { with_start: 1 }, 12),
        ]);
        assert_eq!(
            written,
            vec![
                SnpIndelFitStart {
                    start: 1,
                    ended: StartOutcome::Converged,
                    passes: 41,
                    agreed_with_start: None,
                },
                SnpIndelFitStart {
                    start: 2,
                    ended: StartOutcome::AtThePassLimit,
                    passes: 200,
                    agreed_with_start: None,
                },
                SnpIndelFitStart {
                    start: 3,
                    ended: StartOutcome::AgreedWithAnEarlierStart,
                    passes: 12,
                    agreed_with_start: Some(1),
                },
            ]
        );
    }
}
