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

use crate::ng::parameter_estimation::joint::census::{CensusError, CohortCensusEvidence};
use crate::ng::parameter_estimation::joint::fit::{JointFit, JointFitConfig, JointFitError};
use crate::ng::parameter_estimation::joint::loci::{CensusLoci, CensusLociDigester};
use crate::ng::parameter_estimation::joint::ssr_fit::{
    self, SsrFitConfig, StratumOutcome, gather_strata, strata_of_kept_loci,
};
use crate::ng::types::{ContigId, ReadGroupId};

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
    fn a_fitted_cohorts_inputs() -> (
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
