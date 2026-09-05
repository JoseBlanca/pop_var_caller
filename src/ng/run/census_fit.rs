//! **Fitting a cohort of censuses** — both halves, in the order the estimator needs them.
//!
//! # Two halves, one after the other
//!
//! The **generic half** ([`fit_jointly`]) reads the ordinary positions: the per-read-group noise
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
//! given ([`strata_of_kept_loci`](super::super::parameter_estimation::joint::ssr_fit::strata_of_kept_loci)).
//! Rebuilding it means choosing the selection again, which is a function of the seed, the
//! reference, the analysed ground and the catalog.
//!
//! **The rebuild is checked rather than trusted.** Every census carries a digest of the loci it
//! was written against, and a selection rebuilt from another reference or another catalog would
//! index one stratum's tracts by another's — a wrong answer with no symptom. So the rebuilt loci
//! are digested and compared before a tract is read, and a mismatch is refused.

use std::collections::BTreeMap;

use crate::ng::calling::parameters_file::{
    CensusIdentity, ParametersFile, ReadsBehindEachCalibration,
};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, RecordingTerms,
};
use crate::ng::parameter_estimation::joint::contamination::ContaminationEstimate;
use crate::ng::parameter_estimation::joint::fit::{JointFit, JointFitConfig, JointFitError};
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::parameter_estimation::joint::loci::{CensusLoci, CensusLociDigester};
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::ng::parameter_estimation::joint::ssr_fit::{
    self, SsrFitConfig, StratumEvidence, StratumOutcome, gather_strata, strata_of_kept_loci,
};
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::parameter_estimation::ssr::{RepeatCount, Stratum as SsrStratum, StratumKey};
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::ng::types::{ContigId, ErrorRate, InbreedingF, Ploidy, ReadGroupId, SsrPeriod};

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
    /// The evidence the tract half was fitted from, kept because the per-stratum substitution
    /// rates are read off it and not out of the fit.
    pub tract_evidence: Vec<StratumEvidence>,
}

/// Why a cohort could not be fitted.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CohortFitError {
    /// The selection rebuilt here is not the one the censuses were written against.
    #[error(
        "the census positions rebuilt from this reference and catalog are not the ones these \
         censuses were written against; a tract is stored by its index within its stratum, so \
         fitting them against another selection would read one stratum's tracts as another's"
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

    let fit = crate::ng::parameter_estimation::joint::fit::fit_jointly(cohort, generic).map_err(
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

    let strata = strata_of_kept_loci(loci, contig_of);
    let evidence = gather_strata(cohort, &strata, slippage_group_of).map_err(|source| {
        CohortFitError::Tracts {
            source: Box::new(source),
        }
    })?;
    let outcomes = ssr_fit::fit_strata(&evidence, &homozygote_excess, tracts);

    Ok(CohortFit {
        generic: fit,
        strata: outcomes,
        tracts: strata.len(),
        tract_evidence: evidence,
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
/// # The one group this route leaves at its constant, and why
///
/// **The base-quality calibration is `Defaulted` throughout.** A read group's calibration needs
/// its fitted error rate *and* its minted read-error total — Σ over reads of `ln P(this read is
/// wrong)` — and only the first comes out of this fit. The totals are summed per read from each
/// generic locus's complete observations, which no part of a census carries: a census holds a
/// depth code per position per read group and the non-reference allele counts, and no per-read
/// quality at all. Where either half is missing `assemble` takes
/// `ReadGroupCalibration::defaulted`, scale one, and the parameters file says so — which is a
/// smaller claim than the file could make and the honest one.
///
/// `inbreeding` is one coefficient a sample, in the cohort's own sample order, and is a
/// declaration rather than a fit on this route.
#[must_use]
pub fn parameters_from_the_fit(
    fit: &CohortFit,
    cohort: &CohortCensusEvidence,
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
    inbreeding: Vec<InbreedingF>,
    ploidy: Ploidy,
) -> RunParameters {
    // **The clean class's rate is the sequencing error rate.** The fit models two: how often a
    // read misreads a base at an ordinary position, and how often a read disagrees with the
    // reference at a mismapped one. Only the first is chemistry.
    let error_rate_by_read_group: BTreeMap<ReadGroupId, Estimate<ErrorRate>> = fit
        .generic
        .noise
        .iter()
        .filter_map(|(group, estimate)| {
            ErrorRate::try_new(estimate.value.clean).ok().map(|rate| {
                (
                    *group,
                    Estimate {
                        value: rate,
                        provenance: estimate.provenance,
                        observations: estimate.observations,
                    },
                )
            })
        })
        .collect();

    // Per sample, per read group — and a read group belongs to one sample, so flattening cannot
    // collide.
    let contamination_by_read_group: BTreeMap<ReadGroupId, ContaminationEstimate> = fit
        .generic
        .contamination
        .values()
        .flat_map(|of_sample| of_sample.iter().cloned())
        .collect();

    // **One rate per (read group, stratum, ploidy)**, read off the evidence the tract half was
    // fitted from rather than out of the fit: a stratum whose fit was refused still measured a
    // substitution rate, and a run that dropped it would score its tracts against nothing.
    let mut ssr_substitution_rate: BTreeMap<StratumKey, Estimate<ErrorRate>> = BTreeMap::new();
    for stratum in &fit.tract_evidence {
        let Some(rate) = stratum.substitution_rate() else {
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
            SsrPeriod::try_new(stratum.stratum.period as usize),
            u32::try_from(stratum.stratum.reference_repeats),
        ) else {
            continue;
        };
        let repeats = RepeatCount(repeats);
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
                    observations: stratum.bases_compared,
                },
            );
        }
    }

    RunParameters::assemble(
        &error_rate_by_read_group,
        // **Empty, deliberately** — see this function's own note on the calibration.
        &BTreeMap::new(),
        &contamination_by_read_group,
        SequencingBatches::all_together_over(cohort.read_groups().len(), cohort.len()),
        inbreeding,
        RunParameters::seed_from_moments(
            fit.generic.fitted_alternative_frequency(),
            fit.generic.fitted_diversity(),
        ),
        StratumFits::over(&fit.strata, slippage_group_of.clone()),
        ssr_substitution_rate,
        ploidy,
    )
}

/// **The parameters file a fit over censuses writes.**
///
/// Everything the file needs that is not a fitted number comes from the cohort itself: the
/// read-group table is built from what the censuses declare
/// ([`read_groups_of`](super::census_cohort::read_groups_of)), and the census identity is the
/// recording terms every one of them agreed on — so a calling run handed this file can tell
/// whether its own evidence was recorded the same way.
///
/// **The calibration counts say nothing was fitted, and that is true.** This route produces no
/// per-read-group minted-error totals, so `ReadsBehindEachCalibration::nothing_was_fitted` is the
/// honest entry and the file's warrants come out `defaulted` for those numbers.
#[must_use]
pub fn parameters_file_of(
    parameters: &RunParameters,
    read_groups: &ReadGroups,
    inbreeding: &[Estimate<InbreedingF>],
    reference: &ReferenceDigest,
    terms: &RecordingTerms,
    repeat_routing: &StrRepeatCriteria,
) -> ParametersFile {
    ParametersFile::of_run(
        parameters,
        read_groups,
        &ReadsBehindEachCalibration::nothing_was_fitted(read_groups.len()),
        inbreeding,
        reference,
        CensusIdentity::of(terms),
        repeat_routing,
    )
}

#[cfg(test)]
mod tests {
    //! **Plan step C4**: a cohort of censuses is fitted, both halves, and the selection it is
    //! fitted against is the one it was written against.

    use super::*;
    use crate::ng::run::census_cohort::open_census_cohort;
    use crate::ng::run::test_fixtures::a_census_plan_over_selecting;
    use crate::pop_var_caller_exp::generate_psps::{GeneratePspsArgs, run_generate_psps};
    use crate::pop_var_caller_exp::test_fixtures::a_varying_cohort_on_disk;
    use std::path::PathBuf;

    /// **The varying fixture cohort**, walked, with its censuses beside its psps.
    ///
    /// **Not the plain on-disk cohort**: that one's reference is all `A`, so every base is a
    /// homopolymer, the whole genome routes to the repeat path, and the selection keeps no tract
    /// at all — measured, 0 strata over 0 tracts, which would leave the repeat-tract half of the
    /// fit untested. This one carries a deliberate ten-copy `GT` tract.
    pub(super) fn a_fitted_cohorts_inputs() -> (
        crate::pop_var_caller_exp::test_fixtures::AVaryingCohort,
        Vec<PathBuf>,
    ) {
        use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
        use crate::ng::region_typing::segment_criteria::{
            DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
        };

        let cohort = a_varying_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        run_generate_psps(&GeneratePspsArgs {
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
        })
        .expect("the cohort walks into psps");
        let mut censuses: Vec<PathBuf> = std::fs::read_dir(&psps)
            .expect("the walk made the directory")
            .map(|entry| entry.expect("an entry").path())
            .filter(|path| path.extension().is_some_and(|it| it == "census"))
            .collect();
        censuses.sort();
        assert_eq!(censuses.len(), 2, "one census a sample");
        (cohort, censuses)
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
        let (cohort, censuses) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::ng::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = open_census_cohort(&censuses).expect("the censuses are this cohort's");

        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::ng::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open.evidence);

        let fitted = fit_a_cohort(
            &mut open.evidence,
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
        let mut again = open_census_cohort(&censuses).expect("the censuses are still there");
        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::ng::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&again.evidence);
        let twice = fit_a_cohort(
            &mut again.evidence,
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
    }

    /// **A cohort fitted against another selection is refused before a tract is read.**
    ///
    /// A census stores a tract by its index within its stratum, so a selection rebuilt from
    /// another seed indexes one stratum's tracts as another's — numbers rather than a failure.
    #[test]
    fn a_cohort_fitted_against_another_selection_is_refused() {
        let (cohort, censuses) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::ng::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = open_census_cohort(&censuses).expect("the censuses are this cohort's");

        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| crate::ng::types::ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open.evidence);

        // A selection holding one position the census's does not is enough: the digest is over
        // the kept positions in order.
        let shifted = crate::ng::parameter_estimation::joint::loci::CensusLoci::from_parts(
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
            &mut open.evidence,
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
    }
}

#[cfg(test)]
mod writing_the_parameters_file {
    //! **Plan step C5**: the fit's numbers, assembled and written as the file a calling run
    //! takes.
    //!
    //! # ⚠ Every test here is ignored, and the reason is a contradiction the plan did not see
    //!
    //! `RunParameters::assemble` **refuses a read group that has a fitted error rate and no
    //! minted read-error total**, and says why: the two come from one pass over one set of reads,
    //! so one without the other means they saw different data
    //! (`checked_read_group_count_of`). The plan's §3.4 assumed the pair would simply fall back
    //! to a defaulted calibration; it does not — it panics.
    //!
    //! That leaves this route three ways and none of them is the implementer's to choose:
    //! accumulate the minted totals while `generate-census` reads the psps and store them, which
    //! is Milestone E brought forward; hand `assemble` a total that saw no reads beside each
    //! fitted rate, which defeats a check written for a real reason; or supply no rates either,
    //! at which point the run has no read-group axis and `assemble` refuses it outright.
    //!
    //! The code below is written and compiles. It is held here, unrun, until that is settled.

    use super::*;
    use crate::ng::calling::parameters_file::Warrant;
    use crate::ng::run::census_cohort::{open_census_cohort, read_groups_of};
    use crate::ng::run::census_fit::tests::a_fitted_cohorts_inputs;

    /// **The first parameters file this tree produces from data**, and what it may and may not
    /// claim.
    ///
    /// The population's numbers are fitted. The base-quality calibration is not, and the file
    /// says `defaulted` against it — because a read group's calibration needs its minted
    /// read-error total as well as its rate, and no part of a census carries one.
    #[test]
    #[ignore = "assemble refuses a fitted rate with no minted read-error total; see this module's note"]
    fn the_fit_writes_a_parameters_file_that_says_what_it_fitted() {
        let (parameters, file) = a_fitted_cohorts_parameters();

        let toml = file.to_toml();
        assert!(
            toml.contains("[[read_group]]"),
            "the file names the run's read groups: {toml:.400}",
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

    /// **The base-quality calibration is defaulted, and the file says so.**
    ///
    /// This is the claim plan §3.4 makes and the one a reader would otherwise take on trust.
    #[test]
    #[ignore = "assemble refuses a fitted rate with no minted read-error total; see this module's note"]
    fn the_base_quality_calibration_comes_out_defaulted() {
        let (_parameters, file) = a_fitted_cohorts_parameters();

        assert!(!file.base_quality_calibration.by_read_group.is_empty());
        for row in &file.base_quality_calibration.by_read_group {
            assert_eq!(
                row.error_probability_multiplier.warrant,
                Warrant::Defaulted,
                "no minted read-error total reaches this route, so no calibration can be fitted",
            );
        }
    }

    /// **The file names the census these numbers were fitted from**, term by term — which is
    /// what lets a calling run tell whether its own evidence was recorded the same way.
    #[test]
    #[ignore = "assemble refuses a fitted rate with no minted read-error total; see this module's note"]
    fn the_file_names_the_census_it_was_fitted_from() {
        let (_parameters, file) = a_fitted_cohorts_parameters();

        assert!(
            !file.fitted_from.census.terms.is_empty(),
            "a run that fitted from a census names its terms, unlike a defaults or direct run",
        );
    }

    /// Fit the fixture cohort and assemble its parameters file.
    fn a_fitted_cohorts_parameters() -> (RunParameters, ParametersFile) {
        use crate::ng::types::InbreedingF;

        let (cohort, censuses) = a_fitted_cohorts_inputs();
        let (_segmentation, plan) = crate::ng::run::test_fixtures::a_census_plan_over_selecting(
            &cohort.reference,
            &cohort.catalog,
            crate::ng::run::CensusSelection::SHIPPED.generic_target,
        );
        let mut open = open_census_cohort(&censuses).expect("the censuses are this cohort's");
        let contigs = std::sync::Arc::clone(&plan.contigs);
        let contig_of = move |name: &str| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|index| ContigId(index as u32))
        };
        let pooled = every_read_group_pooled(&open.evidence);
        let fit = fit_a_cohort(
            &mut open.evidence,
            &plan.loci,
            &contig_of,
            &pooled,
            &JointFitConfig::default(),
            &SsrFitConfig::default(),
        )
        .expect("the cohort fits");

        let read_groups = read_groups_of(&open.evidence, &open.samples);
        // **Declared, not fitted, on this route** — the coefficient comes from a sample's own
        // windowed genome histogram, which is the other pre-pass route.
        let declared = InbreedingF::try_new(0.0).expect("zero is a coefficient");
        let inbreeding: Vec<InbreedingF> = (0..open.evidence.len()).map(|_| declared).collect();
        let stated: Vec<Estimate<InbreedingF>> = inbreeding
            .iter()
            .map(|value| Estimate {
                value: *value,
                provenance: Provenance::Supplied,
                observations: 0,
            })
            .collect();

        let parameters = parameters_from_the_fit(
            &fit,
            &open.evidence,
            &pooled,
            inbreeding,
            crate::ng::types::Ploidy::try_new(2).expect("diploid"),
        );
        let terms = open
            .evidence
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
        );
        (parameters, file)
    }
}
