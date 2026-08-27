//! **Everything the parameter pre-pass froze, gathered once for a run** — and the four rules
//! that turn its outputs into what calling reads.
//!
//! [`FrozenParameters`] borrows; something has to own. This module is that owner, and it is
//! the one place the pre-pass's shapes meet calling's: the pre-pass reports per sample and per
//! read group in maps keyed by id, calling reads dense slices indexed by id, and the
//! translation is exactly the kind of join that is silent when it goes wrong.
//!
//! # The four rules, and each has a failure attached
//!
//! **A read group with no usable rate gets a scale of one and says so.** The calibration is a
//! multiplier on each read's own reported error probability, and *no fit* is not the same claim
//! as *a fitted zero*: a zero scale charges every read of the library the floor, which is
//! maximal confidence about every base drawn from a number that says the fit found no errors at
//! all. [`ReadGroupCalibration::from_fitted_rate`] refuses that, and this takes
//! [`ReadGroupCalibration::defaulted`] — scale one, provenance `Defaulted` — whenever it does.
//!
//! **Contamination is absent or measured, never a fitted zero.** Where *no* read group
//! identified a fraction, the run is [`FrozenParameters::uncontaminated`] and the read
//! likelihood computes its plain formula — which is the simple case for that model, not the
//! weak one (`doc/devel/ng/spec/read_likelihoods.md` §3.6). Where **some** did, every read group
//! needs an entry, and a group that identified nothing gets a view whose fraction is zero and
//! whose evidence counts are zero: [`ContaminationView::was_measured`] is what tells that apart
//! from a group measured and found clean, and only the counts can.
//!
//! **The inbreeding coefficients arrive in the run's sample order and are stored as they
//! arrive.** The pre-pass keys them by sample *name*; the run owns the order, and mapping one to
//! the other belongs where the run is assembled, not here. What this checks is that there is at
//! least one.
//!
//! **The prior's seed is built once**, by [`RunParameters::seed_from_moments`], and not per
//! locus: what varies per locus is how the seed is spread across that locus's alleles
//! (`doc/devel/ng/spec/calling_priors.md` §2.3).
//!
//! # The read-group axis has to be dense, and this is where that is checked
//!
//! [`crate::ng::calling::likelihood::ReadGroupParameters::calibration_of`] indexes by
//! `read_group.get() as usize`, so the run's read-group ids are `0..n` with nothing missing. The
//! pre-pass's maps are keyed by id and could carry any set.
//!
//! **What a gap actually costs, measured rather than assumed.** The dense vectors are built over
//! `0..count` by *keyed lookup*, so nothing slides: the missing id's slot holds
//! [`ReadGroupCalibration::defaulted`], and the vectors come out shorter than the highest id, so
//! the **highest read group is dropped entirely**. Its symptom is a panic in `calibration_of` at
//! whichever locus first carries one of that library's reads — which names the read group and not
//! the run, arrives after the pre-pass is long finished, and is a whole run's work thrown away.
//! So the ids are checked here, once, where the dense vectors are built: not because a gap is
//! silent, but because failing at assembly is the difference between a message about the run and
//! a message about a locus.

use std::collections::BTreeMap;

use super::genotype_prior::SpectrumSeed;
use super::genotype_prior::seed_generic::seed_from_population_moments;
use super::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT;
use super::run_report::{
    ContaminationUsed, ReadGroupContamination, RunParameterReport, SequencingBatchingUsed,
};
use super::{ContaminationView, FrozenParameters, ReadGroupCalibration};
use crate::ng::parameter_estimation::Estimate;
use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
use crate::ng::parameter_estimation::joint::contamination::{
    ContaminationEstimate, ContaminationSource,
};
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::parameter_estimation::ssr::StratumKey;
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::types::{
    ErrorRate, ExpectedAlternativeFrequency, ExpectedHeterozygosity, InbreedingF, Ploidy,
    ReadGroupId,
};

// **What stood here until 2026-08-27: `FittedFrequencySpectrum`.**
//
// It owned the vector of `2N + 1` allele-count class weights the genotype prior's search borrowed,
// projected from the joint fit's density at the panel the run had. **The search is deleted**
// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5) and both of the seed's numbers are now
// integrals of the density itself, so nothing evaluates it into classes on the way and there is
// no vector to own.
//
// **⚑ One thing it computed is still owed, and this is where the debt is recorded.** Alongside the
// weights it carried **how many census positions came out variable across the panel** — one minus
// the two end classes, times the positions the density was fitted on, which is *not* the share
// that segregates in the population and differed from it 6.6-fold at one individual on a
// tomato-like density. Spec §6.2 needs that quantity beside the two moments and re-sources it from
// the fit's own per-position posteriors, as a **soft count** — the sum over positions of the
// probability that the position segregates in this panel. **Nothing computes it today**; it is
// Milestone C's step C2.

/// **The pre-pass's outputs, owned for the run** — what [`FrozenParameters`] borrows.
///
/// Assembled once, before any locus is called, and never written afterwards: every error rate,
/// contamination fraction and inbreeding coefficient arrives fitted and leaves unchanged
/// (`doc/devel/ng/spec/calling_em_loop.md` §5).
#[derive(Debug)]
pub struct RunParameters {
    calibration_by_read_group: Vec<ReadGroupCalibration>,
    contamination_by_read_group: Vec<ContaminationView>,
    /// Who was sequenced beside whom, as the run was told — the grouping the contaminating
    /// population is drawn from. **Owned even where no contamination was fitted**, because it
    /// is a fact about the run rather than about the fit; what the view does with it then is
    /// drop it, since a run with no mixture has nothing to read it
    /// ([`FrozenParameters::uncontaminated`]).
    sequencing_batches: SequencingBatches,
    inbreeding_coefficient_by_sample: Vec<InbreedingF>,
    prior_seed: SpectrumSeed,
    ssr_slippage_fits: StratumFits,
    ssr_substitution_rate: BTreeMap<StratumKey, Estimate<ErrorRate>>,
    ploidy: Ploidy,
}

impl RunParameters {
    /// **The prior's seed for this run**, built from the two moments of the population the
    /// pre-pass fitted.
    ///
    /// A named step of its own rather than an argument, because it is the one number here that is
    /// *derived* rather than gathered. See [`seed_from_population_moments`] for what it does with
    /// each moment and what it falls back to when one is missing.
    ///
    /// **Both moments come off the same fitted density**, in closed form and with no panel in
    /// either: `FrequencyDensity::expected_alternative_frequency` and
    /// [`JointFit::fitted_diversity`](crate::ng::parameter_estimation::joint::fit::JointFit::fitted_diversity).
    /// **The panel's inbreeding coefficient is no longer an input**, because nothing in the seed
    /// reads a panel any more (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §2).
    #[must_use]
    pub fn seed_from_moments(
        expected_frequency: Option<ExpectedAlternativeFrequency>,
        diversity: Option<ExpectedHeterozygosity>,
    ) -> SpectrumSeed {
        seed_from_population_moments(expected_frequency, diversity)
    }

    /// **Gather one run's frozen parameters.**
    ///
    /// Three of the four maps are the pre-pass's own per-read-group shapes, keyed by
    /// [`ReadGroupId`]; the fourth, the repeat-tract substitution rates, is keyed by its own
    /// `StratumKey`. The inbreeding coefficients arrive already in the run's sample order,
    /// because the run owns that order and the pre-pass keys by sample name.
    ///
    /// # Panics
    ///
    /// Held in release, because each is a run-assembly bug whose symptom is a wrong genotype
    /// rather than a crash (`doc/devel/ng/spec/calling_em_loop.md` §8):
    ///
    /// - **the run has at least one read group and at least one sample.** Every read of a run
    ///   belongs to a read group, so an empty axis is a run whose read groups went missing
    ///   rather than a run with none.
    /// - **the read-group ids are `0..n` with nothing missing.** The dense vectors are built
    ///   over `0..count` by keyed lookup, so a gap does not misattribute anything — it drops the
    ///   highest read group, and `ReadGroupParameters::calibration_of` panics at whichever locus
    ///   first carries one of its reads. Refusing here turns that into a message about the run.
    /// - **every read group with a contamination estimate is one of the run's.** Nothing else
    ///   checks that map's keys, and an estimate for a read group past the axis is dropped in
    ///   silence — a contaminated library left uncorrected.
    /// - **every read group that has a minted-error total has a fitted rate, or neither.**
    ///   The two come from the same pass over the same reads; one without the other means the
    ///   accumulator and the fit saw different data.
    /// - **the batching covers this run's read groups and this run's samples.** Both axes, since
    ///   a run of one library per sample cannot tell them apart.
    ///
    /// **`sequencing_batches` is the one argument here that is *declared* rather than fitted** —
    /// who ran beside whom, as the run was told, which is the population a contaminant's
    /// genotype is drawn against. Its default is one batch holding everything, which is what a
    /// run that says nothing gets and what every benchmark cohort here has.
    #[allow(
        clippy::too_many_arguments,
        reason = "these are the run's inputs to calling, and the list is the point: a bundle \
                  would be a second type naming the same nine things, and the constructor exists \
                  so that a run cannot forget one of them"
    )]
    #[must_use]
    pub fn assemble(
        error_rate_by_read_group: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        minted_by_read_group: &BTreeMap<ReadGroupId, MintedReadErrors>,
        contamination_by_read_group: &BTreeMap<ReadGroupId, ContaminationEstimate>,
        sequencing_batches: SequencingBatches,
        inbreeding_coefficient_by_sample: Vec<InbreedingF>,
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: StratumFits,
        ssr_substitution_rate: BTreeMap<StratumKey, Estimate<ErrorRate>>,
        ploidy: Ploidy,
    ) -> Self {
        assert!(
            !inbreeding_coefficient_by_sample.is_empty(),
            "every sample of the run carries an inbreeding coefficient and a run has at least \
             one sample, so an empty list is a run whose sample order went missing"
        );
        let read_group_count =
            checked_read_group_count_of(error_rate_by_read_group, minted_by_read_group);
        assert!(
            read_group_count > 0,
            "every read of a run belongs to a read group and a run has at least one, so a run \
             whose calibration inputs name none is one whose read-group axis went missing"
        );

        let calibration_by_read_group: Vec<ReadGroupCalibration> = (0..read_group_count)
            .map(|group| {
                let group = ReadGroupId(group as u32);
                // **Both halves are there by construction** — the count above refuses a rate
                // without its accumulator total, or the reverse, so this `zip` never takes its
                // empty arm. What is left for `from_fitted_rate` to refuse is a **zero** rate,
                // which would give a zero scale and charge every read of the library the error
                // floor; that takes the honest `Defaulted`, scale one, marked as such.
                error_rate_by_read_group
                    .get(&group)
                    .zip(minted_by_read_group.get(&group))
                    .and_then(|(rate, minted)| {
                        ReadGroupCalibration::from_fitted_rate(rate, *minted)
                    })
                    .unwrap_or_else(ReadGroupCalibration::defaulted)
            })
            .collect();

        for group in contamination_by_read_group.keys() {
            assert!(
                (group.get() as usize) < read_group_count,
                "read group {} has a contamination estimate and is not one of the run's {} read \
                 groups, so the contamination fit and the calibration fit were run over \
                 different read-group sets — and an estimate past the axis is dropped in \
                 silence, which leaves a contaminated library uncorrected",
                group.get(),
                read_group_count
            );
        }
        let views: Vec<Option<ContaminationView>> = (0..read_group_count)
            .map(|group| {
                contamination_by_read_group
                    .get(&ReadGroupId(group as u32))
                    .and_then(ContaminationView::of_estimate)
            })
            .collect();
        // **Absent everywhere, or a view for every group.** `FrozenParameters` carries either
        // an empty list — the uncontaminated run, whose read likelihood computes the plain
        // formula — or one entry per read group. A group that identified nothing inside an
        // otherwise-contaminated run gets a view of zero fraction and zero evidence, because
        // there is no correction to make for what could not be measured, and
        // `ContaminationView::was_measured` is what still tells that apart from a group
        // measured and found clean.
        let contamination_by_read_group = if views.iter().all(Option::is_none) {
            Vec::new()
        } else {
            views
                .into_iter()
                .map(|view| view.unwrap_or(UNMEASURED_READ_GROUP))
                .collect()
        };

        // **The batching is checked against both axes here**, where the axes are minted, rather
        // than at the first locus that carries a read. It is declared by the user and the
        // parameters are fitted, so the two can disagree about how many libraries a run has
        // without either being wrong on its own.
        assert_eq!(
            sequencing_batches.read_group_count(),
            calibration_by_read_group.len(),
            "the declared sequencing batching covers {} read groups and the run has {}; a \
             batching minted over a different read-group set would score some library against \
             the neighbours of another",
            sequencing_batches.read_group_count(),
            calibration_by_read_group.len()
        );
        assert_eq!(
            sequencing_batches.sample_count(),
            inbreeding_coefficient_by_sample.len(),
            "the declared sequencing batching covers {} samples and the run has {}; the \
             sample-keyed batching is read by the run's own sample index",
            sequencing_batches.sample_count(),
            inbreeding_coefficient_by_sample.len()
        );
        Self {
            calibration_by_read_group,
            contamination_by_read_group,
            sequencing_batches,
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ssr_substitution_rate,
            ploidy,
        }
    }

    /// What calling borrows for the whole run.
    #[must_use]
    pub fn view(&self) -> FrozenParameters<'_> {
        if self.contamination_by_read_group.is_empty() {
            FrozenParameters::uncontaminated(
                &self.calibration_by_read_group,
                &self.inbreeding_coefficient_by_sample,
                self.prior_seed,
                &self.ssr_slippage_fits,
                &self.ssr_substitution_rate,
                self.ploidy,
            )
        } else {
            FrozenParameters::new(
                &self.calibration_by_read_group,
                &self.contamination_by_read_group,
                &self.sequencing_batches,
                &self.inbreeding_coefficient_by_sample,
                self.prior_seed,
                &self.ssr_slippage_fits,
                &self.ssr_substitution_rate,
                self.ploidy,
            )
        }
    }

    /// How many read groups the run has.
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.calibration_by_read_group.len()
    }

    /// **What this run scored its reads under, in a form an output can print** — the
    /// contamination fraction each read group was corrected for, the batching those fractions
    /// were drawn against, and the repeat-tract constant nothing measured.
    ///
    /// **A genotype computed at a fraction of 3 in 100 and one computed at zero are
    /// indistinguishable in the call**, which is why spec §3.6 requires the run to state what it
    /// used. This is the route from the parameters to that statement; it reads them and computes
    /// nothing.
    ///
    /// **The grain is one row per read group, naming the sample it belongs to.** The spec asks
    /// for the fraction per sample and the fit produces one per read group; a sample's read
    /// groups can carry genuinely different fractions, so a per-sample line would have to pick
    /// one or average them and would erase the claim that they differ. A sample sequenced once,
    /// in one lane — every sample of both benchmark cohorts here — gets exactly one row.
    ///
    /// `read_groups` is the run's read-group table, which is where the sample names, the read
    /// groups' own `@RG ID`s and the library names live; the parameters carry none of them.
    ///
    /// # Panics
    ///
    /// If `read_groups` does not describe this run — a different read-group count, or a
    /// different sample count. Both are the same failure as the batching check at assembly: two
    /// tables minted from different inputs, joined positionally. Here the symptom would be a
    /// fraction reported under another read group's name, which is worse than a crash because it
    /// looks like an answer.
    #[must_use]
    pub fn report(&self, read_groups: &ReadGroups) -> RunParameterReport {
        assert_eq!(
            read_groups.len(),
            self.calibration_by_read_group.len(),
            "the read-group table names {} libraries and the run's parameters cover {}; the two \
             were minted from different inputs, and joining them positionally would report one \
             library's contamination fraction under another's name",
            read_groups.len(),
            self.calibration_by_read_group.len()
        );
        assert_eq!(
            read_groups.read_groups_per_sample().len(),
            self.inbreeding_coefficient_by_sample.len(),
            "the read-group table names {} samples and the run's parameters cover {}; the sample \
             index each row carries is an index into the run's sample order, which is this \
             table's own first-seen order",
            read_groups.read_groups_per_sample().len(),
            self.inbreeding_coefficient_by_sample.len()
        );

        // **Absent, not a fitted zero.** The gathered list is empty exactly where no read group
        // identified a fraction, which is the same condition `view` reads to hand calling the
        // uncontaminated parameters — so the two cannot drift apart.
        let contamination = if self.contamination_by_read_group.is_empty() {
            ContaminationUsed::NoneFitted
        } else {
            // **Walked by sample and then by that sample's read groups**, so the rows come out in
            // the run's sample order however the read-group ids are spread across the samples. A
            // walk over the read-group axis instead would come out in id order, which is the same
            // list only where the two orders agree — and they agree in every fixture that gives
            // each sample one library.
            let mut rows = Vec::with_capacity(self.contamination_by_read_group.len());
            for (sample, of_sample) in read_groups.read_groups_per_sample().iter().enumerate() {
                for &read_group in &of_sample.read_groups {
                    let declared = read_groups.get(read_group);
                    rows.push(ReadGroupContamination {
                        sample,
                        sample_name: of_sample.sample.clone(),
                        read_group,
                        read_group_name: declared.id.clone(),
                        library: declared.library.clone(),
                        estimate: self.contamination_by_read_group[read_group.get() as usize],
                    });
                }
            }
            ContaminationUsed::PerReadGroup(rows)
        };

        let sequencing_batching = if self.sequencing_batches.is_default() {
            SequencingBatchingUsed::DefaultedToOneBatch
        } else {
            SequencingBatchingUsed::Declared {
                batches: self.sequencing_batches.batch_count(),
            }
        };

        RunParameterReport::new(contamination, sequencing_batching, DEFAULT_OUTLIER_WEIGHT)
    }

    /// Who was sequenced beside whom, as the run was told.
    ///
    /// **The one thing here that is declared rather than fitted**, and the one an output has to
    /// be able to report alongside a contamination fraction: two runs under different batchings
    /// produce frequencies that are not comparable, and
    /// [`SequencingBatches::is_default`] is the only thing that can tell a declared batching
    /// from an assumed one (`doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
    #[must_use]
    pub fn sequencing_batches(&self) -> &SequencingBatches {
        &self.sequencing_batches
    }
}

/// A read group that identified no contamination, inside a run where some group did.
///
/// **Zero evidence rather than a zero fraction with evidence behind it**, which is the whole of
/// the distinction: `was_measured` reads the counts, and a group that touched no marker was not
/// measured whatever its fraction says.
const UNMEASURED_READ_GROUP: ContaminationView = ContaminationView {
    fraction: 0.0,
    markers_with_reads: 0,
    reads_on_markers: 0,
    // **`source` is meaningless here, and no variant says so.** Nothing was fitted for this read
    // group, so neither *this library's reads* nor *the whole sample's reads* is true — and that
    // field is the one a consumer reads to answer whether two libraries of one plant may be said
    // to differ. **The counts are what must be read first**, through
    // `ContaminationView::was_measured`; a `NotMeasured` variant on `ContaminationSource` is what
    // would make this unrepresentable, and it belongs to that type's owner.
    source: ContaminationSource::TheWholeSamplesReads,
};

/// How many read groups the run has — **and the two refusals that make the dense build safe**,
/// which is why the name says `checked`.
///
/// # Panics
///
/// On a gap in the ids, and on a rate without a minted total or the reverse. Both are
/// run-assembly bugs: the first drops the highest read group and defers the failure to a locus,
/// the second says the accumulator and the fit saw different reads.
fn checked_read_group_count_of(
    error_rate_by_read_group: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    minted_by_read_group: &BTreeMap<ReadGroupId, MintedReadErrors>,
) -> usize {
    let named: std::collections::BTreeSet<ReadGroupId> = error_rate_by_read_group
        .keys()
        .chain(minted_by_read_group.keys())
        .copied()
        .collect();
    for (position, group) in named.iter().enumerate() {
        assert_eq!(
            group.get() as usize,
            position,
            "the run's read-group ids are 0..n with nothing missing, because the dense vectors \
             are built over 0..count and calling indexes them by the id itself — so a gap makes \
             the vectors shorter than the highest id, the highest read group gets no calibration \
             at all, and `calibration_of` panics at whichever locus first carries one of its \
             reads. The ids given were {:?}",
            named.iter().map(|group| group.get()).collect::<Vec<_>>()
        );
        assert_eq!(
            error_rate_by_read_group.contains_key(group),
            minted_by_read_group.contains_key(group),
            "read group {} has {} — the fitted rate and the accumulator's total come from one \
             pass over one set of reads, so one without the other means they saw different data",
            group.get(),
            if error_rate_by_read_group.contains_key(group) {
                "a fitted rate and no minted-error total"
            } else {
                "a minted-error total and no fitted rate"
            }
        );
    }
    named.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::genotype_prior::SeedRegime;
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::contamination::NotIdentifiedReason;
    use crate::ng::parameter_estimation::joint::fit::FrequencyDensity;
    use crate::ng::parameter_estimation::ssr::{RepeatCount, Stratum};
    use crate::ng::read::input::read_groups::NameOrigin;
    use crate::ng::types::BatchId;
    use crate::ng::types::SsrPeriod;
    use std::collections::BTreeSet;

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
    }
    /// The diversity that density itself implies — **the number the seed's total is now solved
    /// from**, so a test that projects the density must hand the projection this and not a
    /// convenient constant (`doc/devel/ng/spec/ordinary_site_seed.md` §3).
    ///
    /// On this fixture it is 6.06 differences per 10,000 bases: 4.0 in 1,000 positions segregate
    /// and a Beta(0.20, 1.00) population is heterozygous at 0.152 of those.
    fn its_own_diversity(density: &Estimate<FrequencyDensity>) -> Option<ExpectedHeterozygosity> {
        Some(
            ExpectedHeterozygosity::try_new(density.value.expected_heterozygosity())
                .expect("a fitted density's heterozygosity is a probability"),
        )
    }

    /// The mean alternative-allele frequency that density itself implies — **the seed's other
    /// input**, and the one that replaced the panel-fitted search on 2026-08-27
    /// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §2).
    ///
    /// On this fixture it is 1.67 in 1,000: 4.0 in 1,000 positions segregate at a Beta(0.20,
    /// 1.00) mean of 0.167, plus the 1.0 in 1,000 fixed for a non-reference base.
    fn its_own_frequency(
        density: &Estimate<FrequencyDensity>,
    ) -> Option<ExpectedAlternativeFrequency> {
        Some(
            ExpectedAlternativeFrequency::try_new(density.value.expected_alternative_frequency())
                .expect("a fitted density's mean frequency is a probability"),
        )
    }

    /// A density whose four parameters are all different, at a diversity near this project's
    /// tomato cohort's — 6 differences per 10,000 bases.
    fn a_fitted_density() -> Estimate<FrequencyDensity> {
        Estimate {
            value: FrequencyDensity {
                p_invariant: 0.9950,
                p_fixed_alt: 0.0010,
                a: 0.20,
                b: 1.00,
            },
            provenance: Provenance::FittedHere,
            observations: 250_000,
        }
    }

    /// **The seam that was missing**: a fitted density produces both of the seed's numbers, and
    /// the run seeds from them rather than from a species-range constant.
    ///
    /// **The panel it was fitted on no longer enters.** Both numbers are integrals of the curve,
    /// so a run of one sample and a run of a thousand get the same seed off the same density
    /// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §2).
    ///
    /// **⚠ What this can and cannot see, because a review measured it.** Both assertions below
    /// take their truth from the same two accessors that supplied the seed's inputs, so with
    /// respect to *those functions* they are identities: whatever
    /// `expected_alternative_frequency` returns, the seed reports it back. What they do catch is
    /// everything between — the two moments reaching the builder in the right order, the identity
    /// that solves the total, and the pair coming back unswapped, which the first assertion alone
    /// could not see. **The test that can see a defect in either moment function is
    /// [`the_seed_is_what_the_identity_gives_from_the_four_fitted_numbers`] below**, which starts
    /// from the four literal fitted numbers.
    #[test]
    fn a_fitted_density_seeds_the_run_from_its_own_moments() {
        let density = a_fitted_density();

        let seed = RunParameters::seed_from_moments(
            its_own_frequency(&density),
            its_own_diversity(&density),
        );
        assert!(
            matches!(seed.regime(), SeedRegime::FittedCurve),
            "both moments arrived, so the run is on its own measurement rather than a fallback: \
             {:?}",
            seed.regime()
        );
        assert!(
            seed.alpha_alt_total() > 0.0 && seed.alpha_ref() > 0.0,
            "the pair is ({}, {})",
            seed.alpha_ref(),
            seed.alpha_alt_total()
        );
        // **The seed reproduces the density's own heterozygosity**, which is the whole claim: a
        // Dirichlet(`α_ref`, `α_alt`) makes a diploid heterozygous
        // `2 α_ref α_alt / (A (A + 1))` of the time, with `A` the pair's total.
        let total = seed.alpha_ref() + seed.alpha_alt_total();
        let implied = 2.0 * seed.alpha_ref() * seed.alpha_alt_total() / (total * (total + 1.0));
        let theta = density.value.expected_heterozygosity();
        assert!(
            (implied / theta - 1.0).abs() < 1e-12,
            "the seed implies {implied:e} where the density's own heterozygosity is {theta:e}"
        );
        // **And its expected frequency is the density's own**, which is the number that used to
        // come out of a search over the panel's allele-count classes.
        let frequency = seed.alpha_alt_total() / total;
        let densitys_own = density.value.expected_alternative_frequency();
        assert!(
            (frequency / densitys_own - 1.0).abs() < 1e-12,
            "the seed's expected frequency is {frequency:e} where the density's is {densitys_own:e}"
        );
    }

    /// **The seed is the pair a hand calculation gives, to the last few bits** — on the fixture
    /// density, whose four numbers make both moments short enough to write down.
    ///
    /// **The density's own two moments**: 4.0 in 1,000 positions segregate, at a `Beta(0.20,
    /// 1.00)` whose mean is `0.2/1.2` and whose `E[2f(1−f)]` is `2·0.2·1.0/(1.2·2.2)`, plus 1.0
    /// in 1,000 positions fixed for a non-reference base. So the mean frequency is
    /// `0.001 + 0.004 · 0.2/1.2` and the heterozygosity is `0.004 · 2 · 0.2/(1.2 · 2.2)`.
    ///
    /// **Then `ordinary_site_seed.md` §3's identity**: `t = θ / (2 f (1 − f))`, `A = t/(1 − t)`,
    /// and the pair is `(A(1 − f), A f)`.
    ///
    /// **This is what makes the seam's arithmetic checkable rather than merely self-consistent.**
    /// The two assertions in the test above compare the seed against the density through the
    /// library's own accessors; this one writes the whole chain out from the four fitted numbers
    /// and compares against that.
    ///
    /// **It runs on two densities, and the second is there because the first has `b = 1`** — where
    /// `a/(a+b)` and `a/(a+1)` are the same number, so a formula that lost the `b` would pass. The
    /// second is the human-like shape, `Beta(0.35, 1.20)`, whose `a · b` is also not one, so
    /// deleting the Beta's shape from the heterozygosity fails there too.
    #[test]
    fn the_seed_is_what_the_identity_gives_from_the_four_fitted_numbers() {
        for (p_invariant, p_fixed_alt, a, b) in
            [(0.9950, 0.0010, 0.20, 1.00), (0.9949, 0.0004, 0.35, 1.20)]
        {
            let density = Estimate {
                value: FrequencyDensity {
                    p_invariant,
                    p_fixed_alt,
                    a,
                    b,
                },
                provenance: Provenance::FittedHere,
                observations: 250_000,
            };
            let segregating = 1.0 - p_invariant - p_fixed_alt;
            let frequency = p_fixed_alt + segregating * a / (a + b);
            let theta = segregating * 2.0 * a * b / ((a + b) * (a + b + 1.0));
            let share_of_ceiling = theta / (2.0 * frequency * (1.0 - frequency));
            let total = share_of_ceiling / (1.0 - share_of_ceiling);

            let seed = RunParameters::seed_from_moments(
                its_own_frequency(&density),
                its_own_diversity(&density),
            );
            assert!(
                (seed.alpha_ref() / (total * (1.0 - frequency)) - 1.0).abs() < 1e-12,
                "on Beta({a}, {b}) the reference concentration is {} where the identity gives {}",
                seed.alpha_ref(),
                total * (1.0 - frequency)
            );
            assert!(
                (seed.alpha_alt_total() / (total * frequency) - 1.0).abs() < 1e-12,
                "on Beta({a}, {b}) the alternative concentration is {} where the identity gives {}",
                seed.alpha_alt_total(),
                total * frequency
            );
        }
    }

    fn outbred() -> InbreedingF {
        InbreedingF::try_new(0.0).expect("an outbred sample")
    }

    fn human_like_seed() -> SpectrumSeed {
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape)
    }

    fn no_strata() -> StratumFits {
        StratumFits::over(&[], BTreeMap::new())
    }

    /// The default batching for the run these two maps describe, over `samples` samples — one
    /// batch holding all of it, which is what a run that declared nothing gets.
    ///
    /// **The read-group count is taken from the maps rather than restated**, so that a fixture
    /// which grows a read group cannot leave a batching behind describing the old run — which
    /// `assemble` would then refuse, for a reason that has nothing to do with what the test is
    /// about. Both counts are floored at one, because two of the tests here hand `assemble` an
    /// empty axis on purpose and expect *its* refusal rather than this fixture's.
    fn one_batch_for(
        rates: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        totals: &BTreeMap<ReadGroupId, MintedReadErrors>,
        samples: usize,
    ) -> SequencingBatches {
        let libraries = rates
            .keys()
            .chain(totals.keys())
            .map(|group| group.get() as usize + 1)
            .max()
            .unwrap_or(1);
        one_batch(libraries, samples.max(1))
    }

    /// The default batching for a run of `libraries` read groups over `samples` samples.
    ///
    /// **The two counts are separate arguments** because they are separate axes: the batching's
    /// two views are keyed by different things and are different lengths whenever a sample has
    /// more than one library.
    fn one_batch(libraries: usize, samples: usize) -> SequencingBatches {
        let names: Vec<(String, String)> = (0..libraries)
            .map(|library| {
                (
                    format!("rg{library}"),
                    format!("s{}", library.min(samples.saturating_sub(1))),
                )
            })
            .collect();
        let borrowed: Vec<(&str, &str)> = names
            .iter()
            .map(|(id, sample)| (id.as_str(), sample.as_str()))
            .collect();
        let groups = crate::ng::read::input::read_groups::ReadGroups::of_libraries(&borrowed);
        assert_eq!(
            groups.read_groups_per_sample().len(),
            samples,
            "the fixture asked for {samples} samples over {libraries} libraries"
        );
        SequencingBatches::all_together(&groups)
    }

    /// A fitted per-base error rate with the warrant the fit gave it.
    fn fitted_rate(rate: f64, provenance: Provenance) -> Estimate<ErrorRate> {
        Estimate {
            value: ErrorRate::try_new(rate).expect("a legal error rate"),
            provenance,
            observations: 1_000,
        }
    }

    /// The accumulator's total for a read group whose reads averaged `error` per read.
    fn minted(error: f64, reads: u32) -> MintedReadErrors {
        MintedReadErrors::of_observation(error.ln() * f64::from(reads), reads)
    }

    /// A contamination fraction this read group's own reads produced.
    fn estimated(alpha: f64, markers_with_reads: u64) -> ContaminationEstimate {
        ContaminationEstimate::Estimated {
            alpha,
            source: ContaminationSource::ThisReadGroupsReads,
            panel_markers: 10_000,
            markers_with_reads,
            reads_on_markers: markers_with_reads * 3,
            leverage: 1.0,
        }
    }

    fn one_read_group(
        rates: &[(u32, Estimate<ErrorRate>)],
        totals: &[(u32, MintedReadErrors)],
    ) -> (
        BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        BTreeMap<ReadGroupId, MintedReadErrors>,
    ) {
        (
            rates
                .iter()
                .map(|(group, rate)| (ReadGroupId(*group), rate.clone()))
                .collect(),
            totals
                .iter()
                .map(|(group, total)| (ReadGroupId(*group), *total))
                .collect(),
        )
    }

    /// **A read group the fit measured gets a scale, and it carries the *rate's* warrant.**
    ///
    /// The scale is the fitted rate over the reads' own mean reported error, so both halves have
    /// to be there. A rate borrowed from a sibling read group makes a *borrowed* calibration,
    /// and stamping it `FittedHere` would launder it.
    #[test]
    fn a_measured_read_group_gets_a_scale_and_keeps_the_rates_warrant() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::Borrowed))],
            &[(0, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        let calibration = view.calibration_by_read_group();
        assert_eq!(calibration.len(), 1);
        // The reads' mean reported error is 0.004 and the fit says 0.002, so every read's own
        // probability is halved. **Not to the bit**: the accumulator sums each read's log error
        // in fixed point, precisely so that merging shards in different orders gives the same
        // denominator, and the price of that determinism is a quantised mean — measured here at
        // 0.500000000282.
        assert!(
            (calibration[0].scale - 0.5).abs() < 1e-9,
            "the fitted rate over the reads' own mean: {}",
            calibration[0].scale
        );
        assert_eq!(
            calibration[0].provenance,
            Provenance::Borrowed,
            "a borrowed rate makes a borrowed calibration"
        );
    }

    /// **A read group with a total and no usable rate gets scale one, and says `Defaulted`.**
    ///
    /// *No fit* and *a fitted zero* are different claims and only one is safe to multiply by: a
    /// zero rate would give a zero scale, which charges every read of the library the floor —
    /// maximal confidence about every base, from a number saying the fit found no errors at all.
    #[test]
    fn a_read_group_the_fit_could_not_measure_gets_scale_one_and_says_so() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.0, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        let calibration = view.calibration_by_read_group();
        assert_eq!(calibration[0].scale, 1.0);
        assert_eq!(calibration[0].provenance, Provenance::Defaulted);
    }

    /// **Where no read group identified a fraction, the run is uncontaminated** — absent, not a
    /// fitted zero. The read likelihood then computes its plain formula, which is the simple
    /// case for that model rather than the weak one.
    #[test]
    fn a_run_where_nothing_was_identified_is_uncontaminated() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let contamination = BTreeMap::from([(
            ReadGroupId(0),
            ContaminationEstimate::NotIdentified {
                reason: NotIdentifiedReason::NoPanel,
            },
        )]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        assert!(run.view().contamination_is_absent());
    }

    /// **Where some read group identified a fraction, every group gets a view — and a group
    /// that identified nothing is distinguishable from one measured and found clean.**
    ///
    /// Both come back near zero: the search keeps zero when it has nothing to go on, and a clean
    /// library is genuinely near zero. Only the evidence counts tell them apart, which is what
    /// `was_measured` reads.
    #[test]
    fn an_unmeasured_read_group_is_told_apart_from_one_measured_and_clean() {
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
                (2, fitted_rate(0.002, Provenance::FittedHere)),
                (3, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[
                (0, minted(0.004, 500)),
                (1, minted(0.004, 500)),
                (2, minted(0.004, 500)),
                (3, minted(0.004, 500)),
            ],
        );
        let contamination = BTreeMap::from([
            (ReadGroupId(0), estimated(0.03, 4_000)),
            // Measured, and found clean.
            (ReadGroupId(1), estimated(0.000_01, 4_000)),
            // Measured and refused.
            (
                ReadGroupId(2),
                ContaminationEstimate::NotIdentified {
                    reason: NotIdentifiedReason::TooFewMarkers,
                },
            ),
            // **And read group 3 has no entry at all** — which is what stops "one view per read
            // group" and "one view per estimate" being the same number on this fixture.
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        assert!(!view.contamination_is_absent());
        let views = view.contamination_by_read_group();
        assert_eq!(
            views.len(),
            4,
            "one entry per read group, not per estimate — and read group 3 has no estimate"
        );
        assert!((views[0].fraction - 0.03).abs() < 1e-12);
        assert!(
            views[1].was_measured() && views[1].fraction < 1e-4,
            "read group 1 was measured and came back clean"
        );
        for unmeasured in [2, 3] {
            assert!(
                !views[unmeasured].was_measured() && views[unmeasured].fraction == 0.0,
                "read group {unmeasured} was not measured, and its zero is not a measurement"
            );
            // **The one field of an unmeasured view that is not a fact about it.** Nothing was
            // fitted, so no `source` is true; the value is pinned so that a change to it is a
            // decision rather than a drift, and `was_measured` stays the gate a consumer reads
            // first.
            assert_eq!(
                views[unmeasured].source,
                ContaminationSource::TheWholeSamplesReads
            );
        }
    }

    /// **A declared batching minted over a different set of libraries is refused where the run
    /// is assembled.** The batching comes from the user and the read-group axis from the fit, so
    /// the two can disagree without either being wrong on its own — and what a run would get
    /// instead is some library scored against the neighbours of another.
    #[test]
    #[should_panic(expected = "batching covers 2 read groups and the run has 1")]
    fn a_batching_over_another_runs_libraries_is_refused_at_assembly() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch(2, 2),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// **And the other axis.** A run of one library per sample cannot tell the two apart, which
    /// is why both are checked: the sample-keyed batching is read by the run's own sample index.
    #[test]
    #[should_panic(expected = "batching covers 1 samples and the run has 2")]
    fn a_batching_over_another_runs_samples_is_refused_at_assembly() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch(1, 1),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// **The batching travels to calling and back out**, because a run that reports a
    /// contamination fraction has to be able to say what population it was drawn against —
    /// and a defaulted batching and a declared one holding every library are the same dense
    /// value.
    #[test]
    fn the_run_keeps_the_batching_it_was_assembled_with() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch(1, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        assert!(run.sequencing_batches().is_default());
        assert_eq!(run.sequencing_batches().batch_count(), 1);
    }

    /// The inbreeding coefficients are the run's own order, stored as they arrive: the pre-pass
    /// keys them by sample name, and mapping that onto the run's order belongs where the run is
    /// assembled.
    #[test]
    fn the_inbreeding_coefficients_keep_the_order_they_arrive_in() {
        // **One library per sample, because the batching this run carries covers both axes** —
        // and a run of three samples reading one read group is one where two of them were never
        // sequenced.
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
                (2, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[
                (0, minted(0.004, 500)),
                (1, minted(0.004, 500)),
                (2, minted(0.004, 500)),
            ],
        );
        // **Three distinct values, because a palindrome cannot see a reversal** — and reversal
        // is the wrong implementation a map-keyed source most easily produces. `[0.0, 0.9, 0.0]`
        // is bit-identical to its own reverse, and the suite stayed green under one.
        let arriving = [0.1, 0.9, 0.5];
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, arriving.len()),
            arriving
                .iter()
                .map(|f| InbreedingF::try_new(*f).expect("a legal coefficient"))
                .collect(),
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        let stored: Vec<f64> = view
            .inbreeding_coefficient_by_sample()
            .iter()
            .map(|f| f.get())
            .collect();
        assert_eq!(stored.len(), 3);
        for (sample, expected) in arriving.iter().enumerate() {
            assert!(
                (stored[sample] - expected).abs() < 1e-12,
                "sample {sample} arrived at {expected} and is stored at {}",
                stored[sample]
            );
        }
    }

    /// **The substitution rate is looked up at the candidate's repeat count, not the reference
    /// tract's** — the same rule the slippage lookup states, and for the same reason: a read's
    /// chance of mismatching is a property of the tract it was copied from.
    #[test]
    fn the_substitution_rate_is_keyed_by_the_candidates_repeat_count() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let period = SsrPeriod::try_new(2).expect("a dinucleotide");
        let substitution = BTreeMap::from([
            (
                StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: Stratum::new(period, RepeatCount(6)),
                    ploidy: diploid(),
                },
                fitted_rate(0.001, Provenance::FittedHere),
            ),
            (
                StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: Stratum::new(period, RepeatCount(12)),
                    ploidy: diploid(),
                },
                fitted_rate(0.004, Provenance::FittedHere),
            ),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            substitution,
            diploid(),
        );

        let view = run.view();
        let at_six = view
            .ssr_substitution_rate_at(ReadGroupId(0), period, RepeatCount(6))
            .expect("the fit holds this stratum");
        let at_twelve = view
            .ssr_substitution_rate_at(ReadGroupId(0), period, RepeatCount(12))
            .expect("the fit holds this stratum");
        assert!((at_six.value.get() - 0.001).abs() < 1e-12);
        assert!(
            (at_twelve.value.get() - 0.004).abs() < 1e-12,
            "a candidate twice as long is a different stratum, and this one mismatches four \
             times as often"
        );
        assert!(
            view.ssr_substitution_rate_at(ReadGroupId(0), period, RepeatCount(30))
                .is_none(),
            "a candidate several repeats from any fitted stratum is an ordinary absence"
        );
    }

    /// The reverse of `a_rate_without_its_accumulator_total_is_refused`, and **the likelier
    /// direction**: the accumulator runs over every read of every library, so a total with no
    /// rate beside it is what a fit that declined to model one library produces. Measured:
    /// dropping the minted map from the id union leaves the suite green and the read group
    /// silently absent from the run.
    #[test]
    #[should_panic(expected = "saw different data")]
    fn a_minted_total_without_its_fitted_rate_is_refused() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500)), (1, minted(0.004, 500))],
        );
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// **A contamination estimate for a read group the run does not have is refused**, because
    /// nothing else looks at that map's keys and an estimate past the axis is dropped in
    /// silence — which leaves a contaminated library uncorrected, and its fraction charged
    /// nowhere.
    #[test]
    #[should_panic(expected = "not one of the run's")]
    fn a_contamination_estimate_off_the_read_group_axis_is_refused() {
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[(0, minted(0.004, 500)), (1, minted(0.004, 500))],
        );
        let contamination = BTreeMap::from([
            (ReadGroupId(1), estimated(0.03, 4_000)),
            (ReadGroupId(2), estimated(0.09, 4_000)),
        ]);
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// **A read group whose accumulator saw no read gets scale one too**, and for the third of
    /// the three reasons `from_fitted_rate` names: there is no average to divide by. The
    /// module's other fixture covers the zero-rate reason; this one covers the zero-reads one,
    /// which is what a library present in the header and absent from the data looks like.
    #[test]
    fn a_read_group_whose_accumulator_saw_no_read_gets_scale_one() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, MintedReadErrors::default())],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        assert_eq!(view.calibration_by_read_group()[0].scale, 1.0);
        assert_eq!(
            view.calibration_by_read_group()[0].provenance,
            Provenance::Defaulted
        );
    }

    /// **A fit noisier than the reads reported gives a scale above one**, which every other
    /// fixture here cannot see: they all put the fit *below* the reported mean, where a
    /// transposed division gives the same 0.5 the test asserts.
    #[test]
    fn a_fit_noisier_than_the_reads_reported_scales_above_one() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.008, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let view = run.view();
        assert!(
            (view.calibration_by_read_group()[0].scale - 2.0).abs() < 1e-8,
            "0.008 over a reads' mean of 0.004 doubles every read's own probability: {}",
            view.calibration_by_read_group()[0].scale
        );
    }

    /// **The lookup key carries the run's ploidy**, because the pre-pass fits each ploidy's loci
    /// apart — and a run that is not diploid must find its own strata rather than looking for
    /// somebody else's. Every other fixture in this module is diploid, where a hard-coded two is
    /// indistinguishable from reading the run's.
    #[test]
    fn the_substitution_rate_is_keyed_by_the_runs_ploidy() {
        let haploid = Ploidy::try_new(1).expect("a haploid");
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let period = SsrPeriod::try_new(2).expect("a dinucleotide");
        let substitution = BTreeMap::from([(
            StratumKey {
                read_group: ReadGroupId(0),
                stratum: Stratum::new(period, RepeatCount(6)),
                ploidy: haploid,
            },
            fitted_rate(0.001, Provenance::FittedHere),
        )]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            substitution,
            haploid,
        );
        let view = run.view();
        let at_six = view
            .ssr_substitution_rate_at(ReadGroupId(0), period, RepeatCount(6))
            .expect("a haploid run finds its own stratum");
        assert!((at_six.value.get() - 0.001).abs() < 1e-12);
    }

    /// **And it carries the read group**, which one library alone cannot say: two libraries at
    /// the same stratum mismatch at different rates, because chemistry belongs to the library.
    #[test]
    fn the_substitution_rate_is_keyed_by_the_read_group() {
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[(0, minted(0.004, 500)), (1, minted(0.004, 500))],
        );
        let period = SsrPeriod::try_new(2).expect("a dinucleotide");
        let substitution = BTreeMap::from([
            (
                StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: Stratum::new(period, RepeatCount(6)),
                    ploidy: diploid(),
                },
                fitted_rate(0.001, Provenance::FittedHere),
            ),
            (
                StratumKey {
                    read_group: ReadGroupId(1),
                    stratum: Stratum::new(period, RepeatCount(6)),
                    ploidy: diploid(),
                },
                fitted_rate(0.006, Provenance::FittedHere),
            ),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            substitution,
            diploid(),
        );
        let view = run.view();
        let first = view
            .ssr_substitution_rate_at(ReadGroupId(0), period, RepeatCount(6))
            .expect("the fit holds this library's stratum");
        let second = view
            .ssr_substitution_rate_at(ReadGroupId(1), period, RepeatCount(6))
            .expect("and the other library's");
        assert!((first.value.get() - 0.001).abs() < 1e-12);
        assert!(
            (second.value.get() - 0.006).abs() < 1e-12,
            "the second library mismatches six times as often at the same stratum"
        );
    }

    /// The run's read-group count is the calibration axis's length — not the contamination
    /// vector's, which is empty on an uncontaminated run.
    #[test]
    fn the_read_group_count_is_the_calibration_axis() {
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[(0, minted(0.004, 500)), (1, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        assert!(run.view().contamination_is_absent());
        assert_eq!(run.read_group_count(), 2);
    }

    /// A gap in the read-group ids is refused: calling indexes its per-group slices by the id
    /// itself, so read group 2's calibration would be read where read group 1's belongs.
    #[test]
    #[should_panic(expected = "0..n with nothing missing")]
    fn a_gap_in_the_read_group_ids_is_refused() {
        let (rates, totals) = one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (2, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[(0, minted(0.004, 500)), (2, minted(0.004, 500))],
        );
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// A fitted rate with no accumulator total, or the reverse, means the fit and the
    /// accumulator saw different reads — refused rather than defaulted, because *defaulted* is
    /// an answer about the data and this is a statement about the run.
    #[test]
    #[should_panic(expected = "saw different data")]
    fn a_rate_without_its_accumulator_total_is_refused() {
        let (rates, totals) =
            one_read_group(&[(0, fitted_rate(0.002, Provenance::FittedHere))], &[]);
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// A run has at least one sample, so no coefficients is a run whose sample order went
    /// missing rather than a run with nobody in it.
    #[test]
    #[should_panic(expected = "sample order went missing")]
    fn a_run_with_no_samples_is_refused() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let _ = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            Vec::new(),
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    /// Every read of a run belongs to a read group, so a run whose calibration inputs name none
    /// is one whose read-group axis went missing.
    #[test]
    #[should_panic(expected = "read-group axis went missing")]
    fn a_run_with_no_read_groups_is_refused() {
        let _ = RunParameters::assemble(
            &BTreeMap::new(),
            &BTreeMap::new(),
            &BTreeMap::new(),
            one_batch(1, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
    }

    // ── E2b: what the run says it used ────────────────────────────────────────────────────

    /// **A run of three read groups over two samples, where the sample order and the read-group
    /// order are deliberately different.**
    ///
    /// Read groups are minted in the order given, samples in first-seen order: `rgA` and `rgC`
    /// both name `s2`, so sample 0 is `s2` holding read groups 0 and 2, and sample 1 is `s1`
    /// holding read group 1. **A report walked over the read-group axis would come out
    /// `rgA, rgB, rgC`; walked over the samples it comes out `rgA, rgC, rgB`** — which is the
    /// accident every contamination fixture in *this* file shared, each giving one sample one
    /// read group so that the two walks were the same list. (The loop's own fixtures include one
    /// sample with two read groups — `summarise_condition`'s
    /// `a_sample_with_two_libraries_reads_its_own_samples_batch` — but there the two are adjacent
    /// identifiers, so the two walks still agree.)
    fn two_samples_one_of_them_two_read_groups() -> ReadGroups {
        ReadGroups::of_libraries(&[("rgA", "s2"), ("rgB", "s1"), ("rgC", "s2")])
    }

    /// The same run's calibration inputs — three read groups, all measured, which this step's
    /// tests do not read but `assemble` requires.
    fn three_calibrated_read_groups() -> (
        BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        BTreeMap<ReadGroupId, MintedReadErrors>,
    ) {
        one_read_group(
            &[
                (0, fitted_rate(0.002, Provenance::FittedHere)),
                (1, fitted_rate(0.002, Provenance::FittedHere)),
                (2, fitted_rate(0.002, Provenance::FittedHere)),
            ],
            &[
                (0, minted(0.004, 500)),
                (1, minted(0.004, 500)),
                (2, minted(0.004, 500)),
            ],
        )
    }

    /// **The run says what contamination it used, one row per read group, listed under its
    /// sample** — spec §3.6's requirement, at the grain the fit produces.
    ///
    /// The three fractions are all different, so a row that read another read group's view is a
    /// different number rather than the same one; and the sample a row names is checked against
    /// the run's sample order, which is not the read-group order here.
    #[test]
    fn every_read_group_reports_its_own_fraction_under_its_own_sample() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let contamination = BTreeMap::from([
            (ReadGroupId(0), estimated(0.03, 4_000)),
            (ReadGroupId(1), estimated(0.07, 2_000)),
            (ReadGroupId(2), estimated(0.05, 1_000)),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        let rows = match report.contamination() {
            ContaminationUsed::PerReadGroup(rows) => rows,
            ContaminationUsed::NoneFitted => panic!("three read groups were measured"),
        };
        assert_eq!(
            rows.len(),
            3,
            "one row per read group, none summarised away"
        );

        // **Sample order, not read-group order**: s2's two read groups first, then s1's.
        let named: Vec<(usize, &str, u32, f64)> = rows
            .iter()
            .map(|row| {
                (
                    row.sample,
                    row.sample_name.as_ref(),
                    row.read_group.get(),
                    row.estimate.fraction,
                )
            })
            .collect();
        assert_eq!(
            named,
            vec![(0, "s2", 0, 0.03), (0, "s2", 2, 0.05), (1, "s1", 1, 0.07),]
        );

        // The two evidence counts and the source travel beside each fraction, unchanged.
        assert_eq!(rows[0].estimate.markers_with_reads, 4_000);
        assert_eq!(rows[0].estimate.reads_on_markers, 12_000);
        assert_eq!(rows[1].estimate.markers_with_reads, 1_000);
        assert_eq!(rows[1].estimate.reads_on_markers, 3_000);
        assert_eq!(rows[2].estimate.markers_with_reads, 2_000);
        assert_eq!(rows[2].estimate.reads_on_markers, 6_000);
        for row in rows {
            assert_eq!(
                row.estimate.source,
                ContaminationSource::ThisReadGroupsReads
            );
            assert!(row.was_measured());
        }
        // Every read group of this run is its own library, so the two names agree here — which
        // is exactly why `a_librarys_two_lanes_are_two_rows_that_name_one_library` exists.
        assert_eq!(rows[0].read_group_name.as_ref(), "rgA");
        assert_eq!(rows[1].read_group_name.as_ref(), "rgC");
        assert_eq!(rows[2].read_group_name.as_ref(), "rgB");
    }

    /// **One sample's two read groups carry two different fractions, and the report keeps
    /// both.**
    ///
    /// This is the whole reason the row is a read group rather than a sample: a neighbouring
    /// library hopping its index contaminates the run it is on and not the plant,
    /// so a per-sample line would have to pick one of 0.06 and 0.0008 or average them, and each
    /// of the three answers says something the fit did not.
    #[test]
    fn a_samples_two_read_groups_keep_their_two_fractions() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let contamination = BTreeMap::from([
            (ReadGroupId(0), estimated(0.06, 4_000)),
            (ReadGroupId(1), estimated(0.01, 4_000)),
            (ReadGroupId(2), estimated(0.0008, 4_000)),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        let of_sample_zero: Vec<f64> = report
            .contamination()
            .rows()
            .iter()
            .filter(|row| row.sample == 0)
            .map(|row| row.estimate.fraction)
            .collect();
        assert_eq!(of_sample_zero, vec![0.06, 0.0008]);
    }

    /// **A library sequenced over two lanes is two rows naming one library, and only the read
    /// group's own name tells them apart.**
    ///
    /// `@RG LB` is a grouping key rather than an identity — the lanes of one preparation share
    /// it — so a row that carried the library name alone would show the same name twice against
    /// two different fractions, and a reader could not say which lane either belonged to.
    ///
    /// **Two lanes of one preparation really can differ in this number**, which is why the two
    /// rows are not a redundancy: index hopping happens on a flowcell, so a library run beside a
    /// contaminated neighbour on one lane and alone on another picks up stray reads on the first.
    ///
    /// Every other fixture here is built by a helper that names each read group's library after
    /// the read group, so on those the two names agree and neither can be told from the other.
    #[test]
    fn a_librarys_two_lanes_are_two_rows_that_name_one_library() {
        // One plant, one library preparation, two lanes; a second plant beside it.
        let groups = ReadGroups::of_lanes(&[
            ("lane1", "plant", "prep1"),
            ("lane2", "plant", "prep1"),
            ("lane3", "other", "prep2"),
        ]);
        let (rates, totals) = three_calibrated_read_groups();
        let contamination = BTreeMap::from([
            (ReadGroupId(0), estimated(0.04, 4_000)),
            (ReadGroupId(1), estimated(0.0002, 4_000)),
            (ReadGroupId(2), estimated(0.01, 4_000)),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        let named: Vec<(&str, &str, &str, f64)> = report
            .contamination()
            .rows()
            .iter()
            .map(|row| {
                (
                    row.sample_name.as_ref(),
                    row.read_group_name.as_ref(),
                    row.library.value.as_ref(),
                    row.estimate.fraction,
                )
            })
            .collect();
        assert_eq!(
            named,
            vec![
                ("plant", "lane1", "prep1", 0.04),
                ("plant", "lane2", "prep1", 0.0002),
                ("other", "lane3", "prep2", 0.01),
            ]
        );

        // **And whether that library name is the file's or ours travels with it.** `of_lanes`
        // builds runs whose headers declared `@RG LB`, so every row here says `Declared`; a
        // report that dropped the origin would say the same thing about a name this pipeline
        // invented, which is a claim about the run that nobody made.
        for row in report.contamination().rows() {
            assert_eq!(row.library.origin, NameOrigin::Declared);
        }
    }

    /// **A library that could not be measured and one measured and found clean report the same
    /// fraction and different rows**, which is the distinction spec §3.6 says must survive to
    /// the output: only the evidence counts tell them apart.
    #[test]
    fn a_library_nothing_could_be_measured_on_is_not_a_library_measured_and_clean() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let contamination = BTreeMap::from([
            (ReadGroupId(0), estimated(0.03, 4_000)),
            // Measured, and found clean.
            (ReadGroupId(1), estimated(0.0, 4_000)),
            // Measured and refused — and read group 2 gets the same fraction from a different
            // absence.
            (
                ReadGroupId(2),
                ContaminationEstimate::NotIdentified {
                    reason: NotIdentifiedReason::TooFewMarkers,
                },
            ),
        ]);
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &contamination,
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        let rows = report.contamination().rows();
        let clean = rows
            .iter()
            .find(|row| row.read_group == ReadGroupId(1))
            .expect("the library that was measured and found clean");
        let unmeasured = rows
            .iter()
            .find(|row| row.read_group == ReadGroupId(2))
            .expect("the library nothing could be measured on");
        assert_eq!(clean.estimate.fraction, unmeasured.estimate.fraction);
        assert!(clean.was_measured());
        assert!(!unmeasured.was_measured());
        assert_eq!(unmeasured.estimate.markers_with_reads, 0);
        assert_eq!(unmeasured.estimate.reads_on_markers, 0);
        // **And the one field of an unmeasured row that is true of nothing.** No variant of
        // `ContaminationSource` says *nothing was fitted here*, so the row carries a value an
        // output would print as a claim; it is pinned so that changing it is a decision rather
        // than a drift, and `was_measured` stays the gate a consumer reads first.
        assert_eq!(
            unmeasured.estimate.source,
            ContaminationSource::TheWholeSamplesReads,
            "meaningless where nothing was measured — read the counts before the source"
        );
    }

    /// **A run whose fit identified nothing anywhere says so, rather than reporting three zero
    /// fractions.** *Absent, not a fitted zero*: what this fixture exercises is a fit that
    /// emitted no estimate for any read group — which is what a one-sample run always gets, since
    /// contamination is a comparison between samples — and the report says *nothing was
    /// estimable* rather than three zeroes, because those are different claims.
    #[test]
    fn a_run_that_fitted_no_contamination_reports_none_rather_than_zeroes() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        assert_eq!(report.contamination(), &ContaminationUsed::NoneFitted);
        assert!(report.contamination().rows().is_empty());
    }

    /// **The batching those fractions were drawn against travels with them**, and a declared
    /// batching of one batch is not the same claim as a defaulted one — the dense per-read-group
    /// view the caller holds cannot tell them apart, because they are the same values.
    #[test]
    fn the_report_says_whether_the_batching_was_declared_or_assumed() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let assemble = |batches: SequencingBatches| {
            RunParameters::assemble(
                &rates,
                &totals,
                &BTreeMap::from([(ReadGroupId(0), estimated(0.03, 4_000))]),
                batches,
                vec![outbred(), outbred()],
                human_like_seed(),
                no_strata(),
                BTreeMap::new(),
                diploid(),
            )
        };

        let defaulted = assemble(SequencingBatches::all_together(&groups));
        assert_eq!(
            defaulted.report(&groups).sequencing_batching(),
            SequencingBatchingUsed::DefaultedToOneBatch
        );

        // s2's two read groups ran together and s1's ran apart, which is a partition no sample
        // straddles — the one shape a declaration of two batches can take on this run.
        let two = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(0), ReadGroupId(2)]),
                BTreeSet::from([ReadGroupId(1)]),
            ],
        )
        .expect("a partition no sample straddles");
        assert_eq!(
            assemble(two).report(&groups).sequencing_batching(),
            SequencingBatchingUsed::Declared { batches: 2 }
        );

        // **And the case this test's own headline is about, which the two above do not reach.**
        // Two batches are told from a defaulted one by their count alone; a run that *declares*
        // one batch holding every read group carries the same batch identifiers as a run that
        // declared nothing, so `is_default` is the only thing that separates them — and a
        // frequency somebody chose to draw over the whole cohort is a different claim from one
        // drawn over it because nobody said otherwise.
        let one = SequencingBatches::declared(
            &groups,
            &[BTreeSet::from([
                ReadGroupId(0),
                ReadGroupId(1),
                ReadGroupId(2),
            ])],
        )
        .expect("one batch holding every read group is a partition");
        assert_eq!(one.of_each_read_group().0, [BatchId(0); 3]);
        assert_eq!(
            SequencingBatches::all_together(&groups)
                .of_each_read_group()
                .0,
            [BatchId(0); 3],
            "the two batchings carry identical identifiers, which is why the report cannot read \
             them off the dense view"
        );
        assert_eq!(
            assemble(one).report(&groups).sequencing_batching(),
            SequencingBatchingUsed::Declared { batches: 1 },
            "declared, not defaulted, though every read group is in one batch either way"
        );
    }

    /// **The run states the repeat-tract outlier weight and states that nothing fitted it.**
    ///
    /// It is one run-wide constant, so it cannot enter a tract's per-`(read group, candidate)`
    /// warrant without marking every tract of every run as defaulted; reported once per run it
    /// says the same true thing and erases nothing. The value is the read likelihood's own
    /// constant rather than a second spelling of 0.01.
    #[test]
    fn the_run_states_the_outlier_weight_it_inherited() {
        let groups = two_samples_one_of_them_two_read_groups();
        let (rates, totals) = three_calibrated_read_groups();
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            SequencingBatches::all_together(&groups),
            vec![outbred(), outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );

        let report = run.report(&groups);
        assert_eq!(
            report.inherited_repeat_tract_outlier_weight(),
            DEFAULT_OUTLIER_WEIGHT
        );
        assert_eq!(report.inherited_repeat_tract_outlier_weight(), 0.01);
    }

    /// **A read-group table that does not describe this run is refused**, because the join is
    /// positional: the sample names and the library names come from the table and the fractions
    /// from the fit, so two tables minted from different inputs would report one library's
    /// fraction under another's name — an answer rather than a crash.
    #[test]
    #[should_panic(expected = "names 3 libraries and the run's parameters cover 1")]
    fn a_report_over_another_runs_read_groups_is_refused() {
        let (rates, totals) = one_read_group(
            &[(0, fitted_rate(0.002, Provenance::FittedHere))],
            &[(0, minted(0.004, 500))],
        );
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch_for(&rates, &totals, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let _ = run.report(&two_samples_one_of_them_two_read_groups());
    }

    /// **And the other axis**, which a run of one library per sample cannot tell from the
    /// first: the sample index each row carries indexes the run's own sample order.
    #[test]
    #[should_panic(expected = "names 2 samples and the run's parameters cover 1")]
    fn a_report_over_another_runs_samples_is_refused() {
        let (rates, totals) = three_calibrated_read_groups();
        let groups = two_samples_one_of_them_two_read_groups();
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            one_batch(3, 1),
            vec![outbred()],
            human_like_seed(),
            no_strata(),
            BTreeMap::new(),
            diploid(),
        );
        let _ = run.report(&groups);
    }
}
