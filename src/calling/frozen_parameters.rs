//! [`FrozenParameters`]: everything the parameter pre-pass fitted, borrowed once for the whole
//! run.

use crate::calling::genotype_prior::SpectrumSeed;
use crate::calling::likelihood::ssr::RepeatTractOutlierWeight;
use crate::calling::{ContaminationView, ReadGroupCalibration};
use crate::misread_context::MisreadTable;
use crate::parameter_estimation::Estimate;
use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::parameter_estimation::joint::stratum_fits::{AlleleLengthVariation, StratumFits};
use crate::parameter_estimation::repeat_strata::{RepeatCount, Stratum as SsrStratum, StratumKey};
use crate::types::{
    BatchId, BatchOfEachReadGroup, BatchOfEachSample, ErrorRate, InbreedingF, Ploidy, ReadGroupId,
    SsrPeriod,
};
use std::collections::BTreeMap;

/// Everything the parameter pre-pass froze, borrowed for the whole run.
///
/// **Assembled once per run and never written during calling.** Every error rate,
/// contamination fraction and inbreeding coefficient arrives fitted and leaves unchanged;
/// what a locus may move is its own allele frequencies
/// (`doc/devel/ng/spec/calling_em_loop.md` §5). Each field is another document's to define;
/// this type only gathers them so that the loop borrows them once, at the boundary where
/// calling begins.
///
/// **Two different axes are in here and they are not interchangeable**, which is why every
/// field says which one it is on. The calibration and the contamination views are per
/// **read group** — a sample sequenced from two libraries has two of each. The inbreeding
/// coefficients are per **sample**, in the run's sample order, the same order
/// [`LocusEvidence`] and [`LocusInference::per_sample`] use.
///
/// [`LocusEvidence`]: super::LocusEvidence
/// [`LocusInference::per_sample`]: super::LocusInference::per_sample
#[derive(Copy, Clone, Debug)]
pub struct FrozenParameters<'a> {
    calibration_by_read_group: &'a [ReadGroupCalibration],
    contamination_by_read_group: &'a [ContaminationView],
    /// **Which sequencing batch each read group ran in** — which row of the contaminant
    /// frequency table this library's reads are scored against.
    ///
    /// **Empty exactly where [`Self::contamination_by_read_group`] is**, and that pairing is
    /// this type's own rule rather than an accident: the batching's only consumer is the
    /// mixture's second half, and a run the fit found no contamination in has no second half
    /// to build. `ContaminationMixture::new` refuses one without the other for the same reason.
    batch_of_each_read_group: &'a [BatchId],
    /// **Which sequencing batch each sample ran in**, in the run's sample order — which batch
    /// row a sample's expected allele copies are added into when the frequencies are built.
    ///
    /// Two views of one partition, and they are different lengths whenever a sample has more
    /// than one library ([`BatchOfEachReadGroup`] says why they are two types).
    batch_of_each_sample: &'a [BatchId],
    /// How many batches the run declares — the number of rows the contaminant frequency table
    /// has. **One under the default batching**, which is one batch holding the whole run.
    batch_count: usize,
    inbreeding_coefficient_by_sample: &'a [InbreedingF],
    prior_seed: SpectrumSeed,
    ssr_slippage_fits: &'a StratumFits,
    /// **The fourth fitted number a repeat tract's scoring context needs, and the one
    /// [`StratumFits`] does not carry.** The slippage lookup holds the level, the direction
    /// split and the fall-off; the per-base substitution rate is fitted alongside them, per
    /// `(read group, stratum)`, and the pre-pass emits it as a map of its own
    /// (`parameter_estimation::ssr::substitution_rates`).
    ///
    /// **Never the SNP/indel path's ε and never a read's own summed quality**
    /// (`doc/devel/ng/spec/read_likelihoods.md` §4.3): a read's error probability is a
    /// per-*read* number, and this term needs a per-*base* rate applied once for each of the
    /// tract's twenty or forty bases. Using the first as the second overcharges by the tract's
    /// length.
    ssr_substitution_rate: &'a BTreeMap<StratumKey, Estimate<ErrorRate>>,
    ploidy: Ploidy,
    /// **How often a read at a repeat tract came from somewhere the model cannot explain** —
    /// the one number here that no fit produces, so it carries its warrant rather than being a
    /// bare share (`doc/devel/ng/spec/parameters_file.md` §3.8).
    ///
    /// **Defaulted unless a constructor is told otherwise** ([`Self::with_repeat_tract_outlier_weight`]),
    /// and the value a run gets that way is the inherited constant it would have scored under
    /// anyway. The one caller with another value is
    /// [`RunParameters::view`], which passes on whatever
    /// its own parameters hold.
    ///
    /// [`RunParameters::view`]: super::run_parameters::RunParameters::view
    repeat_tract_outlier_weight: RepeatTractOutlierWeight,
    /// **Each read group's misread table after runs** (`str_context_error_fit.md` §6); a read
    /// group absent from it scores under [`MisreadTable::defaults`].
    misread_tables: &'a BTreeMap<ReadGroupId, MisreadTable>,
}

/// No read group's own misread table: what a view built without one points at.
static NO_MISREAD_TABLES: BTreeMap<ReadGroupId, MisreadTable> = BTreeMap::new();

impl<'a> FrozenParameters<'a> {
    /// Gather the run's frozen parameters for a run the fit **did** measure contamination
    /// in, with what can be checked between them checked.
    ///
    /// - `calibration_by_read_group` — one entry per read group, indexed by
    ///   [`ReadGroupId`].
    /// - `contamination_by_read_group` — one entry per read group as well.
    /// - `batching` — who was sequenced beside whom, as the run was told. **It arrives with
    ///   the contamination views rather than beside them**, because it has exactly one
    ///   consumer: the frequency the contaminant's genotype is drawn against is the frequency
    ///   among the samples that ran together, so a run with no fitted fraction has nothing to
    ///   read it. The default — one batch holding the whole run — is what a run that declared
    ///   nothing gets, and it costs that run nothing: the table is one row and every read group
    ///   reads the cohort frequency.
    /// - `inbreeding_coefficient_by_sample` — one entry per sample, **in the run's sample
    ///   order**.
    /// - `prior_seed` — the genotype prior's starting concentration for the run: how many
    ///   chromosomes' worth of belief each allele carries before any read is looked at
    ///   ([`genotype_prior`]).
    /// - `ssr_slippage_fits` — the fitted slippage numbers, how often a read gains or loses
    ///   a whole repeat, looked up by read group and by **stratum**, the group of tracts
    ///   that share a motif length and a repeat count. A run with no repeat tracts supplies
    ///   `StratumFits::over(&[], BTreeMap::new())`, a gather over no outcomes rather than
    ///   nothing at all, so the lookup answers *no such stratum* rather than being absent.
    /// - `ploidy` — how many copies of the genome a sample carries.
    ///
    /// # Panics
    ///
    /// If `contamination_by_read_group` is empty. **No contamination fitted anywhere is
    /// spelled [`Self::uncontaminated`]** — one named way to say it, so that a caller
    /// reaches the decision rather than the shortest thing that compiles. That is the same
    /// rule [`ContaminationMixture::uncontaminated`] already carries for the per-locus half
    /// of the same mixture, and the reason is that *absent* and *fitted zero* are different
    /// claims (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
    ///
    /// If `contamination_by_read_group` is not one entry per read group. The two are read
    /// off the same axis and a run whose calibration covers ten groups and whose
    /// contamination covers four is a caller bug the row could only find lazily — at
    /// whichever locus first held a read from the fifth, or never.
    ///
    /// If `batching` does not cover this run's read groups, or does not cover its samples.
    /// **Both axes, because a run of one library per sample cannot tell them apart** — which is
    /// every sample of every benchmark cohort here. The batching is declared by the user and
    /// the parameters are fitted, so the two can disagree about how many libraries a run has
    /// without either being wrong on its own; what a run would get instead is one library
    /// scored against the neighbours of another, or a sample dropped out of every batch's
    /// copies.
    ///
    /// If either per-axis list is empty. A run has at least one sample and at least one
    /// read group, so an empty list is that axis going missing — and the symptom would
    /// otherwise be an out-of-range index at whichever locus first carried a read, naming
    /// neither the axis nor the locus.
    ///
    /// [`genotype_prior`]: super::genotype_prior
    /// [`ContaminationMixture::uncontaminated`]: super::ContaminationMixture::uncontaminated
    #[allow(
        clippy::too_many_arguments,
        reason = "the run's frozen parameters, and the list is the point: every one of them is \
                  another document's to define, and the constructor exists so that a run cannot \
                  forget one"
    )]
    #[must_use]
    pub fn new(
        calibration_by_read_group: &'a [ReadGroupCalibration],
        contamination_by_read_group: &'a [ContaminationView],
        batching: &'a SequencingBatches,
        inbreeding_coefficient_by_sample: &'a [InbreedingF],
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: &'a StratumFits,
        ssr_substitution_rate: &'a BTreeMap<StratumKey, Estimate<ErrorRate>>,
        ploidy: Ploidy,
    ) -> Self {
        assert!(
            !contamination_by_read_group.is_empty(),
            "no contamination fitted anywhere is spelled `FrozenParameters::uncontaminated`, \
             not an empty list — one named way to say it, so that a caller reaches the \
             decision rather than the shortest thing that compiles"
        );
        assert_eq!(
            contamination_by_read_group.len(),
            calibration_by_read_group.len(),
            "contamination is fitted per read group: the run supplied {} calibrations and \
             {} contamination views",
            calibration_by_read_group.len(),
            contamination_by_read_group.len()
        );
        assert_eq!(
            batching.read_group_count(),
            calibration_by_read_group.len(),
            "the batching covers {} read groups and the run supplied {} calibrations, so the \
             batching and the parameters describe different runs — and the mixture would score \
             some library against a batch that is not the one it ran in",
            batching.read_group_count(),
            calibration_by_read_group.len()
        );
        assert_eq!(
            batching.sample_count(),
            inbreeding_coefficient_by_sample.len(),
            "the batching covers {} samples and the run has {}; the sample-keyed batching is \
             read by the run's own sample index, so a shorter one drops the last samples out \
             of every batch's copies and a longer one belongs to another run",
            batching.sample_count(),
            inbreeding_coefficient_by_sample.len()
        );
        let BatchOfEachReadGroup(batch_of_each_read_group) = batching.of_each_read_group();
        let BatchOfEachSample(batch_of_each_sample) = batching.of_each_sample();
        Self::gather(
            calibration_by_read_group,
            contamination_by_read_group,
            batch_of_each_read_group,
            batch_of_each_sample,
            batching.batch_count(),
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ssr_substitution_rate,
            ploidy,
        )
    }

    /// The run's frozen parameters where **no contamination was fitted anywhere** — a single
    /// sample, which has no panel to measure a contaminant against, or a fit that identified
    /// none (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
    ///
    /// **Absent, not a fitted zero**, and this is the one named way to say it. A consumer
    /// that reached for `contamination.get(group).map(|v| v.fraction).unwrap_or(0.0)` would
    /// turn *not estimable* into *estimated and found clean*, which are different claims
    /// about the sample.
    ///
    /// # Panics
    ///
    /// As [`Self::new`], on an empty per-axis list.
    #[must_use]
    pub fn uncontaminated(
        calibration_by_read_group: &'a [ReadGroupCalibration],
        inbreeding_coefficient_by_sample: &'a [InbreedingF],
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: &'a StratumFits,
        ssr_substitution_rate: &'a BTreeMap<StratumKey, Estimate<ErrorRate>>,
        ploidy: Ploidy,
    ) -> Self {
        Self::gather(
            calibration_by_read_group,
            &[],
            // **No batching either, and for the same reason the views are absent**: the
            // grouping's only consumer is the mixture's second half, and there is no mixture.
            // A batching carried here would be a value nothing reads, which is how a
            // *declared* batching and a defaulted one come to be indistinguishable.
            &[],
            &[],
            0,
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ssr_substitution_rate,
            ploidy,
        )
    }

    /// The checks both constructors share, and the only place the fields are written.
    #[allow(
        clippy::too_many_arguments,
        reason = "the private gather behind two public constructors; its list is theirs"
    )]
    fn gather(
        calibration_by_read_group: &'a [ReadGroupCalibration],
        contamination_by_read_group: &'a [ContaminationView],
        batch_of_each_read_group: &'a [BatchId],
        batch_of_each_sample: &'a [BatchId],
        batch_count: usize,
        inbreeding_coefficient_by_sample: &'a [InbreedingF],
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: &'a StratumFits,
        ssr_substitution_rate: &'a BTreeMap<StratumKey, Estimate<ErrorRate>>,
        ploidy: Ploidy,
    ) -> Self {
        assert!(
            !calibration_by_read_group.is_empty(),
            "every read of the run belongs to a read group and a run has at least one, so \
             an empty calibration list is a run whose read-group axis went missing"
        );
        assert!(
            !inbreeding_coefficient_by_sample.is_empty(),
            "every sample of the run carries an inbreeding coefficient and a run has at \
             least one sample, so an empty list is a run whose sample order went missing"
        );
        // **The batching is present exactly where the contamination views are**, and this is
        // what says so. **Held in debug only, and that is not an accident**: this function is
        // private and its two callers are the only doors — [`Self::new`], which always passes
        // a batching, and [`Self::uncontaminated`], which always passes none — so no test can
        // reach either check, and a release check no test can reach is one the suite cannot
        // keep honest. What it guards is a *third* caller being added without the pairing.
        debug_assert_eq!(
            contamination_by_read_group.is_empty(),
            batch_of_each_read_group.is_empty(),
            "the run carries {} contamination views and a batching over {} read groups; the \
             batching's only consumer is the mixture's second half, so a run has both or \
             neither",
            contamination_by_read_group.len(),
            batch_of_each_read_group.len()
        );
        debug_assert_eq!(
            batch_of_each_read_group.is_empty(),
            batch_count == 0,
            "a run with a batching has at least one batch and a run without one has none, and \
             this has {} read groups batched into {batch_count}",
            batch_of_each_read_group.len()
        );
        Self {
            calibration_by_read_group,
            contamination_by_read_group,
            batch_of_each_read_group,
            batch_of_each_sample,
            batch_count,
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ssr_substitution_rate,
            ploidy,
            repeat_tract_outlier_weight: RepeatTractOutlierWeight::defaulted(),
            misread_tables: &NO_MISREAD_TABLES,
        }
    }

    /// One calibration per read group, indexed by
    /// [`ReadGroupId`].
    #[inline]
    #[must_use]
    pub fn calibration_by_read_group(&self) -> &'a [ReadGroupCalibration] {
        self.calibration_by_read_group
    }

    /// **Score repeat tracts under a supplied outlier weight rather than the inherited one.**
    ///
    /// **A builder rather than an eleventh argument to the two constructors**, because no fit
    /// produces this number and almost every call site wants the compiled-in value. The
    /// exception is [`RunParameters::view`], whose two
    /// arms pass on what a parameters file supplied; the other **25** call sites in this
    /// repository — 27 in all, counted with
    /// `grep -rn "FrozenParameters::new(\|FrozenParameters::uncontaminated("` — would each
    /// have gained an argument naming the same constant.
    ///
    /// **A caller that omits it is not silently wrong**: it gets the inherited constant, which
    /// is what the run would have scored under.
    ///
    /// [`RunParameters::view`]: super::run_parameters::RunParameters::view
    #[inline]
    #[must_use]
    pub fn with_repeat_tract_outlier_weight(mut self, weight: RepeatTractOutlierWeight) -> Self {
        self.repeat_tract_outlier_weight = weight;
        self
    }

    /// **The same view with each read group's own misread table** — what
    /// [`RunParameters::view`] passes on.
    ///
    /// [`RunParameters::view`]: super::run_parameters::RunParameters::view
    #[inline]
    #[must_use]
    pub fn with_misread_tables(mut self, tables: &'a BTreeMap<ReadGroupId, MisreadTable>) -> Self {
        self.misread_tables = tables;
        self
    }

    /// **The misread table a read group's reads are scored under at repeat tracts**: its own,
    /// or the default where the run holds none for it.
    #[inline]
    #[must_use]
    pub fn misread_table(&self, read_group: ReadGroupId) -> &'a MisreadTable {
        self.misread_tables
            .get(&read_group)
            .unwrap_or(&*crate::misread_context::DEFAULT_MISREAD_TABLE)
    }

    /// **The share of a repeat tract's reads charged to none of its candidate alleles**, with
    /// its warrant — what
    /// [`TractScoringFits`] reads once a
    /// locus and hands the scoring row.
    ///
    /// [`TractScoringFits`]: super::inference::repeat_tract_parameters::TractScoringFits
    #[inline]
    #[must_use]
    pub fn repeat_tract_outlier_weight(&self) -> RepeatTractOutlierWeight {
        self.repeat_tract_outlier_weight
    }

    /// One contamination view per read group, or empty where the fit identified none —
    /// which [`Self::contamination_is_absent`] is the name for.
    #[inline]
    #[must_use]
    pub fn contamination_by_read_group(&self) -> &'a [ContaminationView] {
        self.contamination_by_read_group
    }

    /// Whether the fit identified no contamination anywhere, in which case the read
    /// likelihood computes the plain formula.
    ///
    /// **The branch that decides which formula runs deserves a name** — see
    /// [`ContaminationMixture::is_absent`], which exists for this reason. Asking it as an
    /// emptiness test announces at the call site only that a list is empty.
    ///
    /// [`ContaminationMixture::is_absent`]: super::ContaminationMixture::is_absent
    #[inline]
    #[must_use]
    pub fn contamination_is_absent(&self) -> bool {
        self.contamination_by_read_group.is_empty()
    }

    /// Which sequencing batch each read group ran in — **empty on an uncontaminated run**,
    /// where there is no mixture to read it.
    ///
    /// The typed wrapper rather than the slice, because the sample-keyed batching is the same
    /// slice type and means something else: the two agree in length at one library per sample,
    /// which is every sample of every benchmark cohort here, so a transposition passes every
    /// shape check and comes back as a wrong contaminant frequency
    /// ([`BatchOfEachReadGroup`]).
    #[inline]
    #[must_use]
    pub fn batch_of_each_read_group(&self) -> BatchOfEachReadGroup<'a> {
        BatchOfEachReadGroup(self.batch_of_each_read_group)
    }

    /// Which sequencing batch each sample ran in, in the run's sample order — **empty on an
    /// uncontaminated run**.
    #[inline]
    #[must_use]
    pub fn batch_of_each_sample(&self) -> BatchOfEachSample<'a> {
        BatchOfEachSample(self.batch_of_each_sample)
    }

    /// The sequencing batch **one sample** ran in, by its index in the run's sample order.
    ///
    /// **It takes a `usize` where `Self::batch_of_read_group` takes a [`ReadGroupId`], and
    /// that is the point.** The two batchings are the same slice type over different axes, and
    /// at one library per sample — every sample of every benchmark cohort here — they are the
    /// same length too, so an index handed to the wrong one passes every shape check and comes
    /// back as a wrong contaminant frequency. Asking through these two makes the transposition
    /// a type error rather than a number.
    ///
    /// # Panics
    ///
    /// On a sample past the run's own count, and on an uncontaminated run, which carries no
    /// batching at all — asking one which neighbours a contaminant it does not have was drawn
    /// from.
    #[inline]
    #[must_use]
    pub fn batch_of_sample(&self, run_sample: usize) -> BatchId {
        *self
            .batch_of_each_sample
            .get(run_sample)
            .unwrap_or_else(|| {
                panic!(
                    "sample {run_sample} has no sequencing batch; this run's batching covers {} \
                 samples",
                    self.batch_of_each_sample.len()
                )
            })
    }

    #[cfg(test)]
    /// The sequencing batch **one read group** ran in. See [`Self::batch_of_sample`] for why the
    /// two are separate calls.
    ///
    /// # Panics
    ///
    /// On a read group the batching does not cover, and on an uncontaminated run.
    #[inline]
    #[must_use]
    pub fn batch_of_read_group(&self, read_group: ReadGroupId) -> BatchId {
        let at = read_group.get() as usize;
        *self.batch_of_each_read_group.get(at).unwrap_or_else(|| {
            panic!(
                "read group {at} has no sequencing batch; this run's batching covers {} read \
                 groups",
                self.batch_of_each_read_group.len()
            )
        })
    }

    /// How many sequencing batches the run declares — **zero on an uncontaminated run**, one
    /// under the default batching, and the number of rows the per-locus contaminant frequency
    /// table has.
    #[inline]
    #[must_use]
    pub fn batch_count(&self) -> usize {
        self.batch_count
    }

    /// One inbreeding coefficient per sample, in the run's sample order.
    #[inline]
    #[must_use]
    pub fn inbreeding_coefficient_by_sample(&self) -> &'a [InbreedingF] {
        self.inbreeding_coefficient_by_sample
    }

    /// How many samples the run has.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.inbreeding_coefficient_by_sample.len()
    }

    /// How many read groups the run has.
    #[inline]
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.calibration_by_read_group.len()
    }

    /// The genotype prior's starting concentration for the run.
    #[inline]
    #[must_use]
    pub fn prior_seed(&self) -> SpectrumSeed {
        self.prior_seed
    }

    /// **The fitted per-base substitution rate for one read group at one candidate's
    /// stratum**, or nothing where the fit has no such stratum — which is ordinary, since a
    /// candidate several repeats from its reference tract's length lands there on perfectly
    /// good data and a caller owes an answer for it.
    ///
    /// **Keyed by the *candidate's* repeat count, never the reference tract's**, which is the
    /// rule [`StratumFits::at`] states and for the same reason: a read's chance of mismatching
    /// is a property of the tract it was copied from, and that is the candidate allele
    /// (`doc/devel/ng/spec/read_likelihoods.md` §4.4). The argument is named for the candidate
    /// so that handing over the locus's own tract length is a mistake somebody has to type on
    /// purpose — the same device that lookup uses.
    ///
    /// The key carries the run's ploidy because the pre-pass fits each ploidy's loci apart.
    #[inline]
    #[must_use]
    pub fn ssr_substitution_rate_at(
        &self,
        read_group: ReadGroupId,
        period: SsrPeriod,
        candidate_repeats: RepeatCount,
    ) -> Option<&'a Estimate<ErrorRate>> {
        self.ssr_substitution_rate.get(&StratumKey {
            read_group,
            stratum: SsrStratum::new(period, candidate_repeats),
            ploidy: self.ploidy,
        })
    }

    /// The fitted slippage numbers, looked up by read group and stratum.
    #[inline]
    #[must_use]
    pub fn ssr_slippage_fits(&self) -> &'a StratumFits {
        self.ssr_slippage_fits
    }

    /// **What this tract's genotype prior is seeded from: a shape and a strength**, with the
    /// rung of the tract ladder they came from beside them
    /// (`doc/devel/ng/spec/population_diversity.md` §4).
    ///
    /// **Keyed by the tract's own reference repeat count, where its two neighbours here are
    /// keyed by the candidate's — and the three sit within a screen of each other in the
    /// tract's parameter assembly.** [`Self::ssr_substitution_rate_at`] and
    /// [`StratumFits::at`] answer *how does a read of **this candidate** go wrong*, which is a
    /// property of the tract a read was copied from. This answers *which lengths can **this
    /// tract** be*, which is one question per locus: the spectrum is indexed by whole-repeat
    /// offset from the reference tract length, so passing a candidate's count would re-centre
    /// the shape on that candidate and flatten the prior. The argument is named for the tract
    /// for the same reason the others are named for the candidate.
    ///
    /// **It always answers**, at one of four rungs — the stratum's own fit, the stratum's own
    /// shares held at the stated concentration, its motif period's pooled tracts, or a flat
    /// shape at the stated concentration — so a repeat tract is never left without a prior.
    /// Which rung it was is on the returned value and belongs in the run's output.
    #[inline]
    #[must_use]
    pub fn ssr_allele_length_variation_at(
        &self,
        period: SsrPeriod,
        reference_repeats: RepeatCount,
    ) -> AlleleLengthVariation<'a> {
        self.ssr_slippage_fits
            .allele_length_variation_at(period.get(), u64::from(reference_repeats.get()))
    }

    /// How many copies of the genome a sample carries.
    #[inline]
    #[must_use]
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calling::tests::{
        NO_SUBSTITUTION_RATES, diploid_ploidy, frozen_parameters, neutral_seed, no_strata,
        one_batch, one_read_group, outbred_samples,
    };
    use crate::parameter_estimation::joint::contamination::ContaminationSource;

    /// A run whose fit produced an allele length variation at one dinucleotide stratum and at one
    /// trinucleotide stratum, with **different shapes and different concentrations** — so that a
    /// lookup answering from the wrong stratum, or from the wrong period, is a different number
    /// rather than the same one.
    fn strata_with_allele_length_variations() -> StratumFits {
        use crate::parameter_estimation::joint::census::Stratum;
        use crate::parameter_estimation::joint::share_curve::ShareSource;
        use crate::parameter_estimation::joint::slippage_curve::LevelSource;
        use crate::parameter_estimation::joint::ssr_fit::{
            LevelProvenance, ShareProvenance, SharesProvenance, Slippage, StratumFit,
            StratumOutcome,
        };

        let level = LevelProvenance {
            source: LevelSource::Cell,
            curve: None,
            reach: None,
            slipped_reads: Some(400.0),
        };
        let share = ShareProvenance {
            source: ShareSource::Stratum,
            curve: None,
            reach: None,
        };
        let one =
            |period: u8, repeats: u64, allele_length_variation: Vec<f64>, concentration: f64| {
                StratumOutcome::Fitted(Box::new(StratumFit {
                    standard_errors: None,
                    stratum: Stratum {
                        period,
                        reference_repeats: repeats,
                    },
                    slippage: vec![Some(Slippage {
                        level: 0.05,
                        shorter_share: 0.83,
                        fall_off: 0.25,
                    })],
                    allele_length_variation,
                    concentration: Some(concentration),
                    samples_fitted_on: None,
                    tracts_of_its_own: 40,
                    reads_crossing: 400,
                    level_provenance: vec![Some(level)],
                    shares_provenance: vec![Some(SharesProvenance {
                        slipped_reads: Some(400.0),
                        shorter_share: share,
                        fall_off: share,
                    })],
                }))
            };
        StratumFits::over(
            &[
                one(2, 10, vec![0.6, 0.3, 0.1], 4.0),
                one(3, 10, vec![0.1, 0.3, 0.6], 25.0),
            ],
            BTreeMap::from([(ReadGroupId(0), 0)]),
        )
    }

    /// **The run's frozen parameters answer a tract's prior shape from the tract's own stratum**,
    /// and the accessor is keyed by the tract's reference repeat count where its two neighbours
    /// on this type are keyed by the candidate's.
    ///
    /// **Three fixtures' worth of coincidence is removed on purpose.** The two strata carry
    /// different motif periods, different shapes and different concentrations, and the repeat
    /// count asked for is not one either stratum would answer by accident — so a lookup off by
    /// one repeat, or reading the wrong period, is a visibly different answer.
    #[test]
    fn the_run_answers_a_tracts_prior_shape_from_its_own_stratum() {
        use crate::parameter_estimation::joint::stratum_fits::AlleleLengthVariationRung;
        use crate::parameter_estimation::repeat_strata::RepeatCount;

        let strata = strata_with_allele_length_variations();
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let batching = SequencingBatches::all_together_over(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let dinucleotide = parameters.ssr_allele_length_variation_at(
            SsrPeriod::try_new(2).expect("a dinucleotide"),
            RepeatCount(10),
        );
        assert_eq!(
            dinucleotide.rung(),
            AlleleLengthVariationRung::StratumsOwnFit
        );
        assert_eq!(dinucleotide.fitted_weights(), Some(&[0.6, 0.3, 0.1][..]));
        assert_eq!(dinucleotide.concentration(), 4.0);

        // Same repeat count, other motif period: a different stratum and a different answer.
        let trinucleotide = parameters.ssr_allele_length_variation_at(
            SsrPeriod::try_new(3).expect("a trinucleotide"),
            RepeatCount(10),
        );
        assert_eq!(trinucleotide.fitted_weights(), Some(&[0.1, 0.3, 0.6][..]));
        assert_eq!(trinucleotide.concentration(), 25.0);

        // A repeat count neither stratum holds falls down the ladder rather than answering.
        let elsewhere = parameters.ssr_allele_length_variation_at(
            SsrPeriod::try_new(2).expect("a dinucleotide"),
            RepeatCount(11),
        );
        assert_eq!(elsewhere.rung(), AlleleLengthVariationRung::StatedFlat);
        // The median is taken on the between-locus share 1 / (1 + c): 4 and 25 have shares 1/5
        // and 1/26, whose mean is 31/260, a concentration of (1 − 31/260) / (31/260) = 229/31 —
        // not the mean of the two concentrations, 14.5.
        let expected = 229.0 / 31.0;
        assert!(
            (elsewhere.concentration() - expected).abs() < 1e-12,
            "the run fitted 4.0 and 25.0, whose median share converts back to 229/31 ≈ 7.39; \
             got {}",
            elsewhere.concentration()
        );
    }

    fn measured_contamination() -> ContaminationView {
        ContaminationView {
            fraction: 0.03,
            markers_with_reads: 400,
            reads_on_markers: 1_200,
            source: ContaminationSource::ThisReadGroupsReads,
        }
    }

    /// Contamination is fitted per read group, so a run whose calibration covers one group
    /// and whose contamination covers two has one of the two indexed by something else.
    ///
    /// **The check has to be at construction, not at the first read.** The row looks a
    /// fraction up by read group id, so a mismatched pair is found at whichever locus first
    /// carries a read from the group past the end — or never, if no read ever comes from
    /// one, in which case every genotype of the run is scored under somebody else's
    /// contamination.
    #[test]
    #[should_panic(expected = "contamination is fitted per read group")]
    fn frozen_parameters_refuse_a_contamination_list_of_another_read_group_count() {
        let calibration = one_read_group();
        let contamination = vec![measured_contamination(), measured_contamination()];
        let inbreeding = outbred_samples(2);
        let strata = no_strata();
        // A batching that matches neither list, so that what refuses this is the check the
        // test is about rather than a fixture that happened to disagree first.
        let _ = frozen_parameters(
            &calibration,
            &contamination,
            &one_batch(2, 2),
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );
    }

    /// No contamination anywhere has a **named** door, and it is a different claim from a
    /// fitted zero: at one sample there is no panel to be surprised by, so the fraction is
    /// not estimable at all rather than estimated and found clean.
    #[test]
    fn a_run_with_no_contamination_fitted_says_so_by_name() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let parameters = FrozenParameters::uncontaminated(
            &calibration,
            &inbreeding,
            neutral_seed(),
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid_ploidy(),
        );

        assert!(parameters.contamination_is_absent());
        assert!(parameters.contamination_by_read_group().is_empty());
        assert_eq!(parameters.read_group_count(), 1);
        assert_eq!(parameters.sample_count(), 1);
        assert_eq!(parameters.ploidy().get(), 2);
        assert_eq!(parameters.prior_seed().alpha_ref(), 1.0);
        assert_eq!(parameters.ssr_slippage_fits().strata(), 0);
    }

    /// The shortest thing that compiles is refused, so a caller reaches the decision.
    ///
    /// **The two spellings are not equivalent and that is the whole point.** A consumer
    /// handed an empty list writes `contamination.get(group).map(|v| v.fraction)
    /// .unwrap_or(0.0)` and has silently turned *not estimable* into *estimated and found
    /// clean* — which are different claims about the sample. The sibling half of the same
    /// mixture, `ContaminationMixture::new`, refuses the empty spelling for this reason,
    /// and this is the frozen half spelled to match.
    #[test]
    #[should_panic(expected = "is spelled `FrozenParameters::uncontaminated`")]
    fn an_empty_contamination_list_is_refused_in_favour_of_the_named_constructor() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let _ = FrozenParameters::new(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            neutral_seed(),
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid_ploidy(),
        );
    }

    /// The two axes the type calls "not interchangeable" are told apart by a fixture where
    /// they differ: one library, three samples.
    ///
    /// **Every fixture with one of each passes a swapped implementation.** `sample_count()`
    /// is what the evidence check compares against, so `read_group_count()` and
    /// `sample_count()` reading each other's field turns the run-order guard into a
    /// read-group-count guard — invisible on a panel of 63 accessions with one library
    /// each, and wrong on any run where a sample is sequenced twice.
    #[test]
    fn read_group_count_and_sample_count_are_two_different_axes() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(3);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        assert_eq!(parameters.read_group_count(), 1);
        assert_eq!(parameters.sample_count(), 3);
        assert_eq!(parameters.calibration_by_read_group().len(), 1);
        assert_eq!(parameters.inbreeding_coefficient_by_sample().len(), 3);
    }

    /// A run whose read-group axis went missing is refused at construction, rather than
    /// producing an out-of-range index at whichever locus first carries a read.
    #[test]
    #[should_panic(expected = "read-group axis went missing")]
    fn frozen_parameters_refuse_a_run_with_no_read_groups() {
        let inbreeding = outbred_samples(2);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let _ = frozen_parameters(
            &[],
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );
    }

    /// **A batching minted over a different set of libraries is refused where the run is
    /// assembled**, not at whichever locus first carries a read from the library it misses.
    ///
    /// The batching is declared by the user and the parameters are fitted, so the two can
    /// disagree about how many libraries a run has without either being wrong on its own — and
    /// the symptom without this is a mixture scoring some library against the neighbours of
    /// another.
    #[test]
    #[should_panic(expected = "the batching covers 2 read groups")]
    fn a_batching_over_another_runs_read_groups_is_refused() {
        let calibration = one_read_group();
        let contamination = vec![measured_contamination()];
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let two_libraries = one_batch(2, 2);
        let _ = FrozenParameters::new(
            &calibration,
            &contamination,
            &two_libraries,
            &inbreeding,
            neutral_seed(),
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid_ploidy(),
        );
    }

    /// **And the other axis**, which is the one a run of one library per sample cannot tell
    /// apart from the first: the sample-keyed batching is read by the run's own sample index, so
    /// a shorter one drops the last samples out of every batch's copies.
    #[test]
    #[should_panic(expected = "the batching covers 1 samples")]
    fn a_batching_over_another_runs_samples_is_refused() {
        let calibration = one_read_group();
        let contamination = vec![measured_contamination()];
        let inbreeding = outbred_samples(2);
        let strata = no_strata();
        let one_library_one_sample = one_batch(1, 1);
        let _ = FrozenParameters::new(
            &calibration,
            &contamination,
            &one_library_one_sample,
            &inbreeding,
            neutral_seed(),
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid_ploidy(),
        );
    }

    /// A run whose sample order went missing is refused, rather than producing a locus with
    /// no calls in it.
    #[test]
    #[should_panic(expected = "an empty list is a run whose sample order went missing")]
    fn frozen_parameters_refuse_a_run_with_no_samples() {
        let calibration = one_read_group();
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let _ = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &[],
            &strata,
            &NO_SUBSTITUTION_RATES,
        );
    }
}
