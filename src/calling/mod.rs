//! Turning a locus's reads into a call: which alleles are in play, how probable each
//! sample's reads are under each candidate genotype, how likely each genotype is before
//! the reads are looked at, and what comes out. Steps 6 to 9 of the ng caller.
//!
//! **One folder rather than four, and a dependency decided it.** As four sibling
//! folders at the top of `src/`, the three calling architecture documents each had
//! to state by hand that the prior and the likelihood must never import from the
//! inference step — a rule the tree made necessary and could not enforce. Under one
//! folder the shared vocabulary sits here, in `calling/mod.rs`, one level above all
//! four sub-modules, so each imports downward and the rule disappears
//! (`doc/devel/ng/arch/module_layout.md`, principle 1b).
//!
//! Steps 10 and 11 are deliberately elsewhere: phasing and locus filtering read what
//! calling produces rather than taking part in producing it.
//!
//! ## What is here
//!
//! The vocabulary that needs nothing from the sub-modules: [`CandidateAlleles`], the
//! alleles a locus is called over; [`ExpectedAlleleCopies`], the fractional allele
//! counts the loop feeds back to itself; and [`LocusInference`] with
//! [`SampleGenotypeCall`], what a locus produces. Beside them, the three the calling
//! loop takes and gives back: [`LocusEvidence`], one locus's reads per sample and the
//! one place the SNP/indel and repeat-tract paths meet; [`FrozenParameters`], everything
//! the parameter pre-pass fitted, gathered so one borrow crosses the seam; and
//! [`CallingScratch`], every buffer a locus's calling fills, allocated once per worker.
//! Of the four sub-modules — the candidate
//! step, the likelihood, the genotype prior and the inference loop — three are here:
//! [`allele_candidates`], step 6, which narrows the merge's allele table to the sequences
//! worth calling over (`doc/devel/ng/impl_plan/candidate_alleles.md`);
//! [`likelihood`], step 7, how probable this sample's reads are given each genotype
//! (`doc/devel/ng/impl_plan/calling_read_likelihoods.md`); and
//! [`genotype_prior`], step 8, how likely each genotype is before any read is looked at
//! (`doc/devel/ng/impl_plan/calling_prior.md`). **The fourth is [`inference`], step 9, the
//! loop that consumes all three** — it holds the seam every way of handling a cohort crosses
//! and the configuration of the three nested loops, two of which ship switched off
//! (`doc/devel/ng/impl_plan/calling_loop.md`).
//!
//! Beside this file rather than inside any one sub-module: [`genotype_table`], which
//! says which genotypes a locus's alleles make and holds the three flat tables the
//! genotype prior and the read likelihood both read them from. It sits here for the
//! same reason the types above do — two of the four sub-modules consume it, so it
//! belongs one level above them — and its three types are re-exported here, so that
//! all of calling's shared vocabulary is named at one depth. Beside it again, and
//! test-only: `genotype_table_parity`, which checks that table value for value against
//! production's `GenotypeShape`, frozen in `testdata/` at promotion step C5. It is a file of its
//! own rather than a block inside [`genotype_table`]'s own tests so that the table's tests stay
//! free of the comparison.
//!
//! The scalars other steps also name — `AlleleId`, `Phred`, `Genotype` — are not here
//! but in [`crate::types`], with the rest of ng's shared vocabulary.

#[cfg(test)]
mod genotype_table_parity;
#[cfg(test)]
mod loop_parity;
#[cfg(test)]
mod quality_parity;

pub mod allele_candidates;
pub mod evidence_shaping;
pub mod genotype_prior;
pub mod genotype_table;
pub mod inference;
pub mod likelihood;
pub mod parameters_file;
pub mod quality;
pub mod run_parameters;

mod candidate_alleles;
mod frozen_parameters;
mod locus_evidence;
mod scratch;

pub use candidate_alleles::{CandidateAlleles, ExpectedAlleleCopies};
pub use frozen_parameters::FrozenParameters;
pub use genotype_table::{GenotypeIdx, GenotypeTable, GenotypeTableView};
pub use likelihood::generic::{
    ERROR_SPREAD_BASES, NO_ERROR_SPREAD, ReferenceBesideLocus, WrongReadCharge,
    WrongReadChargeTable, allele_is_compatible_with_partial, fill_wrong_read_charges,
    partial_row_fits_allele,
};
pub use likelihood::{
    ContaminationMixture, ContaminationView, GenericObservation, GenericRowScratch,
    GenericSampleEvidence, MIN_BASE_ERROR, MIN_CONTAMINANT_FREQUENCY, ReadGroupCalibration,
    ReadGroupParameters, SsrRowScratch, SsrSampleEvidence, fill_batch_allele_copies,
    fill_own_batch_contaminant_row, fill_shared_contaminant_allele_frequencies,
};
use locus_evidence::locus_kind_word;
pub use locus_evidence::{GenericLocusSample, LocusEvidence, SsrLocusShape};
// `EmissionCost` is named only by tests and by doc links in `inference`, which the unused-import
// lint does not count.
#[allow(unused_imports)]
pub(crate) use scratch::EmissionCost;
pub use scratch::{CallingScratch, UNWRITTEN_SCRATCH_VALUE};
pub(crate) use scratch::{CohortSummingBuffers, SampleScoringBuffers};

use crate::calling::quality::ArtifactTestCounts;
use crate::locus_generation::LocusKind;
use crate::parameter_estimation::Provenance;
use crate::parameter_estimation::joint::stratum_fits::AlleleLengthVariationRung;
#[cfg(test)]
use crate::types::Ploidy;
use crate::types::{GenomeRegion, Genotype, GenotypeRows, Phred};

/// One sample's call at one locus.
///
/// **Two outcomes, and the second is not a low-confidence version of the first.** A sample
/// the candidate step declared uncallable has no genotype at all: the locus is called over a
/// set of alleles that cannot represent what the sample carries, so it is set aside before
/// the first pass, scored against nothing, and emitted as missing
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0, §9). Emission must not write that as a
/// genotype with a poor quality beside it, and an enum is what stops it: there is no
/// quality to read.
///
/// **A sample with no reads at the locus is [`Self::Called`], not [`Self::Missing`]**, and
/// the distinction is the one a reader is most likely to get backwards. Such a sample scores
/// every genotype alike, so the prior decides it alone and a genotype comes out — which is
/// the right answer rather than a special case
/// (`doc/devel/ng/spec/calling_em_loop.md` §7). [`Self::Missing`] is not *missing data*; it
/// is *the caller declined to invent a genotype over a set that cannot hold this sample's
/// allele*.
///
/// **Why the type changed shape here.** `doc/devel/ng/arch/calling_em_loop.md` §2 sketches
/// this as a struct of a genotype and a quality, which cannot express the missing case, and
/// [`Genotype::new`] refuses an empty multiset by design — *an empty multiset is not a
/// haploid call, it is a sample with no genome*. The spec's §5.0 records that the ruling had
/// a producer and no carrier and that the loop's own plan adds one; this is it.
#[derive(Clone, PartialEq, Debug)]
pub enum SampleGenotypeCall {
    /// The caller reached a genotype for this sample.
    Called {
        /// Which alleles this sample carries, one per copy of its genome.
        genotype: Genotype,
        /// How sure the caller is of that genotype — the posterior-derived `GQ` of the
        /// loop's last pass, how much of the sample's genotype probability the winning
        /// genotype took. Step 13's quality model **refines** this number; it does not
        /// replace it (`doc/devel/ng/arch/calling_em_loop.md` §2).
        genotype_quality: Phred,
        /// **Whether this sample's own reads said nothing about which genotype it has** —
        /// its genotype likelihoods were flat, so every genotype was equally probable
        /// under them and the prior decided alone.
        ///
        /// **It is what turns a called sample into a `./.`** at emission
        /// (`doc/devel/ng/spec/vcf_output.md` §7.1), and it is recorded **here**, by
        /// whoever scored the sample, because nothing downstream can recover it: the
        /// likelihoods live in per-sample scratch the loop overwrites at the next locus.
        /// `vcf::assemble`'s own module note named this as the one input it could not be
        /// given; this is the answer to it.
        ///
        /// **It must be the likelihood and not the posterior.** A sample with no reads is
        /// scored by the loop and comes back with a genotype, because the prior decides it
        /// alone — and where the fitted frequency is low that posterior is sharply peaked,
        /// so no threshold on the confidence would catch this sample. The likelihood is
        /// what the reads said; the posterior is what the reads said plus what the cohort
        /// assumed.
        ///
        /// **True is common and is not a defect.** On six tomato accessions over 400 kb,
        /// every locus is called on every sample of the run and most samples cover most
        /// loci — but a sample that covered nothing at a locus is scored there all the
        /// same, and this is what says so.
        reads_were_uninformative: bool,
        /// The sample's likelihood and posterior rows over every genotype of the table the locus
        /// was called with — the `PL` and `GP` the VCF writes (`vcf_output.md` §7.2). Taken in
        /// the final pass, where both exist for one sample at a time.
        genotype_rows: GenotypeRows,
    },
    /// The candidate step ruled this sample uncallable at this locus, so it took no part in
    /// the loop and emission writes its `GT` as missing.
    ///
    /// **This happens on the SNP/indel path only.** A repeat tract sets no sample aside
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), and [`LocusInference::new`] refuses
    /// a locus that carries one there.
    Missing,
}

#[cfg(test)]
impl SampleGenotypeCall {
    /// **Test fixtures only**: a called sample's rows replaced by flat ones — every genotype equally
    /// likely and equally probable — over `allele_count` alleles at the genotype's ploidy, so a
    /// fixture built without a calling pass still carries rows the VCF assembly can index.
    #[must_use]
    pub fn with_flat_rows(self, allele_count: usize) -> Self {
        match self {
            Self::Called {
                genotype,
                genotype_quality,
                reads_were_uninformative,
                genotype_rows: _,
            } => {
                let ploidy = Ploidy::try_new(
                    u8::try_from(genotype.alleles().len()).expect("a fixture's ploidy fits a u8"),
                )
                .expect("a fixture's genotype names at least one copy");
                let genotypes = GenotypeTable::build(ploidy, allele_count).genotype_count();
                Self::Called {
                    genotype,
                    genotype_quality,
                    reads_were_uninformative,
                    genotype_rows: GenotypeRows {
                        ploidy,
                        log_likelihoods: vec![0.0; genotypes].into(),
                        posteriors: vec![1.0 / genotypes as f64; genotypes].into(),
                    },
                }
            }
            Self::Missing => Self::Missing,
        }
    }
}

/// **Test fixtures only**: rows with no entries, for a call that is never assembled into a record.
#[cfg(test)]
#[must_use]
pub fn placeholder_genotype_rows() -> GenotypeRows {
    GenotypeRows {
        ploidy: Ploidy::try_new(2).expect("two copies is a ploidy"),
        log_likelihoods: Box::new([]),
        posteriors: Box::new([]),
    }
}

impl SampleGenotypeCall {
    /// The alleles called, or `None` where the sample was set aside.
    #[inline]
    #[must_use]
    pub fn genotype(&self) -> Option<&Genotype> {
        match self {
            Self::Called { genotype, .. } => Some(genotype),
            Self::Missing => None,
        }
    }

    #[cfg(test)]
    /// How sure the caller is, or `None` where the sample was set aside — which carries no
    /// quality, because no genotype was scored to have one.
    #[inline]
    #[must_use]
    pub fn score_best_genotype(&self) -> Option<Phred> {
        match self {
            Self::Called {
                genotype_quality, ..
            } => Some(*genotype_quality),
            Self::Missing => None,
        }
    }

    /// Whether this sample's `GT` is written as missing.
    #[inline]
    #[must_use]
    pub fn is_missing(&self) -> bool {
        matches!(self, Self::Missing)
    }
}

/// **What a repeat tract's call rested on** — the four things about a tract's parameters that
/// the genotypes cannot be asked.
///
/// A tract is scored under two fitted numbers per `(read group, candidate)` pair — how often
/// the copying steps before sequencing add or drop a whole repeat, and how often a base is
/// misread — and under a prior shape saying which lengths are plausible. Each of the three can
/// fall back to a stated constant, and **a call resting on measurements and one resting on
/// constants are different claims about the same reads**
/// (`doc/devel/ng/spec/population_diversity.md` §1). Nothing in the called genotype says which
/// it was, so it travels here.
///
/// # The counts are per `(read group, candidate)` cell, and the table covers the whole run
///
/// A tract's parameter table is `read groups × candidates`, over **every** read group of the
/// run rather than the ones whose reads reached this tract — the read likelihood's own table is
/// indexed by read-group identifier directly. So on a run of many libraries a tract can report
/// defaulted cells on account of a library that sent it nothing. That is the conservative
/// direction and it is what [`Self::scoring_cells`] counts against.
///
/// # Why the fields are read through accessors
///
/// Three of the six are counts of the same type, and two of them are *shares of* another: a
/// record claiming seven defaulted cells out of six, or more cells defaulted by an unknown
/// library than defaulted at all, is arithmetic nonsense that would print as an answer.
/// [`Self::new`] refuses those, and public fields would let a caller build one around it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct RepeatTractProvenance {
    /// **Which rung of the tract ladder this tract's prior shape came from** — the stratum's own
    /// fitted allele length variation, its own shares at the stated concentration, its motif
    /// period's pooled tracts, or a flat shape at the stated concentration
    /// (`doc/devel/ng/spec/population_diversity.md` §4.4).
    ///
    /// Every rung answers, so this is never *no shape*; what it says is how well founded the
    /// shape was. The SNP/indel path has a ladder of its own with different rungs, resolved once
    /// per run rather than per locus, and none of its rungs can appear here — which is why this
    /// whole value is absent there rather than this field being optional.
    ///
    /// *(It replaced a `seed_diversity_unreachable` flag, whose whole subject was a failure that
    /// no longer exists: the prior used to scale a constructed geometric shape to reproduce a
    /// measured gene diversity, which was impossible above a ceiling the shape itself set.
    /// Seeding from the fit asserts no such scaling — `population_diversity.md` §4.2.)*
    allele_length_variation_rung: AlleleLengthVariationRung,
    /// How many `(read group, candidate)` cells this tract was scored over — the denominator the
    /// three counts below are read against.
    scoring_cells: usize,
    /// How many of them took the shipped stutter row because the slippage fit had no numbers for
    /// them.
    ///
    /// **Ordinary rather than an error**: a candidate several repeats from its tract's reference
    /// length lands in no fitted stratum on perfectly good data.
    cells_with_no_fitted_slippage: usize,
    /// Of those, how many were defaulted because **the slippage fit does not describe this run's
    /// read groups** — a library the pre-pass never saw, or a slippage map naming more libraries
    /// than the fit was run over.
    ///
    /// **This absence means the parameters and the reads came from different runs**, where the
    /// other one is a candidate sitting off the fitted range on perfectly good data. That is why
    /// it is counted apart rather than folded in: anything above zero is a run to look at, not a
    /// tract.
    ///
    /// **It sees the slippage fit only, and the substitution rates can fail the same way
    /// unseen.** Those are looked up in a plain map keyed by `(read group, stratum, ploidy)`
    /// whose absence carries no reason
    /// ([`FrozenParameters::ssr_substitution_rate_at`]), so a rate map fitted over a different
    /// set of libraries lands in the count below, indistinguishable from a stratum that was
    /// simply never fitted. Splitting that count the same way means giving that lookup a typed
    /// absence, which is a change to the parameter side and is recorded as open rather than made
    /// here.
    cells_whose_read_group_the_fit_does_not_describe: usize,
    /// How many cells found no fitted substitution rate and took the stated constant, which is
    /// defined as the SNP/indel path's own default so that a run cannot default its two error
    /// parameters to two different guesses.
    ///
    /// **One count, two causes, and nothing here separates them** — see the field above.
    cells_with_no_fitted_substitution_rate: usize,
    /// **Whether the mixture's third term was built** — how common each reachable tract length
    /// is in the contaminating population.
    ///
    /// True exactly where the run's parameter fit found a contamination fraction, since that is
    /// what decides whether a tract gets the three-term form or the two
    /// (`doc/devel/ng/spec/read_likelihoods.md` §4.5.1). It is a run-wide condition read at the
    /// tract, and it is here rather than only in the run's report because a reader holding one
    /// locus's record should not have to fetch the run's to know whether its reads were shared
    /// out with a contaminant.
    contaminant_term_was_built: bool,
}

impl RepeatTractProvenance {
    /// One repeat tract's record, with its three counts checked against each other.
    ///
    /// # Panics
    ///
    /// Held in release, because each is a record that would print as an answer
    /// (`doc/devel/ng/spec/calling_em_loop.md` §8):
    ///
    /// - **more cells with no fitted slippage than cells.** The count is over the same table its
    ///   denominator measures, so a larger one means the two came from different tracts.
    /// - **more cells defaulted by a library the fit does not describe than defaulted at all.**
    ///   The first is documented as a share of the second and is counted inside its own arm; a
    ///   larger one means the two counters were swapped, which is the one mistake a reader of
    ///   these numbers could not detect.
    /// - **more cells with no fitted substitution rate than cells**, for the first reason again
    ///   on the other parameter.
    #[must_use]
    pub fn new(
        allele_length_variation_rung: AlleleLengthVariationRung,
        scoring_cells: usize,
        cells_with_no_fitted_slippage: usize,
        cells_whose_read_group_the_fit_does_not_describe: usize,
        cells_with_no_fitted_substitution_rate: usize,
        contaminant_term_was_built: bool,
    ) -> Self {
        assert!(
            cells_with_no_fitted_slippage <= scoring_cells,
            "a tract cannot have {cells_with_no_fitted_slippage} of its {scoring_cells} cells \
             fall back on the slippage fit: the count and its denominator are over one table, \
             so a larger count came from a different tract"
        );
        assert!(
            cells_whose_read_group_the_fit_does_not_describe <= cells_with_no_fitted_slippage,
            "{cells_whose_read_group_the_fit_does_not_describe} cells were defaulted by a \
             library the fit does not describe and only {cells_with_no_fitted_slippage} were \
             defaulted at all: the first is a share of the second, so this is the two counters \
             swapped"
        );
        assert!(
            cells_with_no_fitted_substitution_rate <= scoring_cells,
            "a tract cannot have {cells_with_no_fitted_substitution_rate} of its \
             {scoring_cells} cells fall back on the substitution rate: the count and its \
             denominator are over one table"
        );
        Self {
            allele_length_variation_rung,
            scoring_cells,
            cells_with_no_fitted_slippage,
            cells_whose_read_group_the_fit_does_not_describe,
            cells_with_no_fitted_substitution_rate,
            contaminant_term_was_built,
        }
    }

    /// Which rung of the tract ladder this tract's prior shape came from — see
    /// [`Self::allele_length_variation_rung`].
    #[inline]
    #[must_use]
    pub fn allele_length_variation_rung(&self) -> AlleleLengthVariationRung {
        self.allele_length_variation_rung
    }

    /// How many `(read group, candidate)` cells this tract was scored over.
    #[inline]
    #[must_use]
    pub fn scoring_cells(&self) -> usize {
        self.scoring_cells
    }

    /// How many of them took the shipped stutter row — see
    /// [`Self::cells_with_no_fitted_slippage`].
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_slippage(&self) -> usize {
        self.cells_with_no_fitted_slippage
    }

    /// Of those, how many because the slippage fit does not describe this run's read groups —
    /// see [`Self::cells_whose_read_group_the_fit_does_not_describe`].
    #[inline]
    #[must_use]
    pub fn cells_whose_read_group_the_fit_does_not_describe(&self) -> usize {
        self.cells_whose_read_group_the_fit_does_not_describe
    }

    /// How many took the stated substitution rate — see
    /// [`Self::cells_with_no_fitted_substitution_rate`].
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_substitution_rate(&self) -> usize {
        self.cells_with_no_fitted_substitution_rate
    }

    /// Whether the mixture's third term was built — see [`Self::contaminant_term_was_built`].
    #[inline]
    #[must_use]
    pub fn contaminant_term_was_built(&self) -> bool {
        self.contaminant_term_was_built
    }
}

/// What calling produces at one locus: the alleles it settled on, every sample's call,
/// and the evidence for how it got there.
///
/// **How the answer was reached travels with the answer, because nothing downstream can
/// reconstruct it.** A genotype from a loop that ran out of passes is a different claim
/// from one that settled; a parameter that was defaulted is a guess where a fitted one
/// is a measurement; a prior seed that could not reach its target is a prior that was
/// bent to fit. Each of those is a field here, and each would otherwise be invisible by
/// the time the record is written (`doc/devel/ng/spec/calling_em_loop.md` §6,
/// `doc/devel/ng/arch/calling_priors.md` §5).
///
/// Plain data, and six of its ten fields are public because every consumer reads
/// them. **The four that are not are private for two different reasons.**
/// [`Self::alleles`], [`Self::cohort_expected_copies`] and
/// [`Self::artifact_test_counts`] are all indexed by, or parallel to, the allele table: a
/// public allele table would let a consumer widen it with [`CandidateAlleles::admit`]
/// against unchanged copies, or leave an artifact summary naming an allele the locus no
/// longer holds — breaking, after construction, the pairings [`Self::new`] checks at it.
/// [`Self::site_quality`] is private for a different rule, which its own comment gives:
/// **there is one quality field, and nothing may read it between the worker that writes
/// the baseline and the stage that overwrites it with the corrected value**
/// (`doc/devel/ng/spec/calling_quality.md` §3.5).
#[derive(Clone, PartialEq, Debug)]
pub struct LocusInference {
    /// Where on the reference this locus is.
    pub region: GenomeRegion,
    /// The alleles it was finally called over — after any discovery round admitted
    /// alleles the first pass missed, and after the prune dropped the ones no sample
    /// carried. Read through [`Self::alleles`].
    alleles: CandidateAlleles,
    /// One call per sample, **in the run's sample order** — the same order every
    /// per-sample slice the loop was given is indexed by. That one order for the whole
    /// run is what makes the loop's fixed-order sum over samples reproducible
    /// (`doc/devel/ng/spec/calling_em_loop.md` §8).
    pub per_sample: Vec<SampleGenotypeCall>,
    /// The cohort's expected allele copies, one entry per allele of [`Self::alleles`].
    ///
    /// **Handed on rather than left to be recomputed**: deriving it downstream from the
    /// called genotypes gives a different number, because a call has already thrown away
    /// the uncertainty these counts carry, and site filtering and emission both read it
    /// (`doc/devel/ng/spec/calling_em_loop.md` §9). Read through
    /// [`Self::cohort_expected_copies`].
    cohort_expected_copies: ExpectedAlleleCopies,
    /// Whether the allele frequencies stopped moving, or the loop simply ran out of
    /// passes.
    ///
    /// **`false` is emitted, never dropped and never fatal.** One hard locus must not
    /// kill a cohort run — production retired its non-convergence error for exactly this
    /// reason — but the flag has to reach the output, because a genotype from a loop that
    /// did not settle is a weaker claim than one from a loop that did, and nothing
    /// downstream can tell them apart otherwise
    /// (`doc/devel/ng/spec/calling_em_loop.md` §6).
    pub converged: bool,
    /// How many passes the frequency loop took. At least one.
    ///
    /// **One pass is a legitimate answer, not a defect.** The loop starts with an E-step
    /// that uses no prior at all — the reads alone — and only then begins iterating
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3), so the first pass already has a
    /// previous estimate to be compared against and can settle against it.
    ///
    /// Carried as an instrument, not as a diagnostic: the pass cap and the convergence
    /// threshold are both inherited from production and neither has been measured on
    /// this caller's range — a thousand samples, or three reads a position — so the
    /// distribution of this number across a real panel is what will set them
    /// (`doc/devel/ng/spec/calling_em_loop.md` §12, question 4).
    pub passes: u32,
    /// The weakest warrant among every parameter that entered this locus's arithmetic.
    ///
    /// The scoring contexts propagate it rather than branching on it: each carries the
    /// weakest [`Provenance`] of the parameters that reached it, and the loop copies the
    /// weakest of those onto the locus (`doc/devel/ng/arch/read_likelihoods.md` §1.4). A
    /// consumer that treats a defaulted error rate like a fitted one is the failure this
    /// exists to prevent.
    ///
    /// **Which warrant is weakest is settled**, and this field is the first consumer of the
    /// answer: [`Provenance::weaker_of`] combines two along the ladder
    /// `parameter_estimation`'s own documentation states — fitted here, borrowed from the
    /// sample's other read groups, supplied, defaulted — so a value the run *supplied* ranks
    /// below a borrowed fit, on the ground that a handed-in number says nothing about this
    /// data where a borrowed one is at least a measurement of a neighbouring grain.
    ///
    /// **The two paths fold two different lists, because they read two different
    /// parameters.** A SNP/indel locus folds the calibrations of the read groups whose reads
    /// reached it. A repeat tract folds the stutter and substitution warrants of its own
    /// `(read group, candidate)` parameter table, which `inference::repeat_tract_parameters`
    /// gathers — the calibration scale never enters a tract's likelihood, so charging a tract
    /// for it would report a worse warrant than the call has.
    ///
    /// **What is in neither is the prior's own shape**, which carries no provenance on the
    /// SNP/indel path at all; at a tract, the rung its shape came from travels on
    /// [`Self::repeat_tract`] rather than through this fold, because that is a
    /// statement about the prior where every rung here is a statement about the reads.
    pub weakest_provenance: Provenance,
    /// **What a repeat tract's call rested on, beyond the warrant above** — `None` at a SNP or
    /// indel, `Some` at every repeat tract.
    ///
    /// See [`RepeatTractProvenance`] for the four things it carries and why each of them cannot
    /// be recovered from the call.
    pub repeat_tract: Option<RepeatTractProvenance>,
    /// **How unlikely it is that no sample here carries a copy of any non-reference
    /// allele** — the site quality, on the Phred scale, as the worker computed it and
    /// *before* the artifact correction that the first output stage applies
    /// (`doc/devel/ng/spec/calling_quality.md` §5).
    ///
    /// **One quality field, written twice, and that is a rule with a shipped defect behind
    /// it.** Production keeps its corrected value nowhere and recomputes it at VCF-encode
    /// time; for sixteen days its `--min-qual` gate compared the engine's baseline while
    /// the *corrected* number went into the `QUAL` column, so sites were emitted `PASS`
    /// carrying a written quality of zero — 40 false positives at 30× on GIAB HG002 and 64
    /// at 50×, against 14 and 14 once both read one function (§3.5). ng carries one field:
    /// the worker writes the baseline here, the correction stage overwrites it in place,
    /// and there is never a second quality for anything to read by mistake.
    ///
    /// **Which is why this field is private and has no public reader.** §3.5's other half
    /// is that nothing between the worker and that stage may read it, and visibility is
    /// what enforces it: [`Self::uncorrected_site_quality`] is `pub(crate)`, so no consumer
    /// outside this crate can see the uncorrected number at all. The public accessor
    /// arrives with the stage that makes the value public-worthy. The check that the value
    /// came from the function that owns the ceiling is [`Self::new`]'s.
    site_quality: Phred,
    /// The nine pooled read counts the artifact correction consumes — or `None` at a locus
    /// that gives the two tests nothing to test.
    ///
    /// **Built here rather than downstream because one of the nine needs the calls**: how
    /// many alternative-allele reads the called genotypes lead you to expect. The other
    /// eight are the evidence's, and the evidence is released with the locus
    /// (`doc/devel/ng/spec/calling_quality.md` §3.3).
    ///
    /// **`None` is the ordinary outcome at more than one built locus in four**, not an
    /// error: the merge builds a locus when some sample's non-reference reads *pooled*
    /// reach its rule, so a locus whose candidate table came back as the reference alone is
    /// a first-class result — measured at 27.4% of built loci on the 63-accession tomato
    /// panel and 27.3% on HG002 at 30×
    /// ([`allele_candidates::SelectionVerdict::Selected`]). Both tests compare an
    /// alternative allele against the reference, so a locus with no alternative, or one
    /// whose alternatives drew no read at all, has nothing for them to weigh; production
    /// hands back its baseline unchanged in exactly those two cases
    /// ([`qual_refine.rs:79`](../../../src/vcf/qual_refine.rs)), and `None` is that same
    /// answer as a type rather than as an early return.
    artifact_test_counts: Option<ArtifactTestCounts>,
}

impl LocusInference {
    /// One locus's finished call, with everything checkable about it checked.
    ///
    /// # Panics
    ///
    /// If `cohort_expected_copies` is not one entry per allele of `alleles`, in either
    /// direction. [`ExpectedAlleleCopies`] is already built against *an* allele table, so
    /// what this catches is the residue: copies built against a **different** table that
    /// happened to be a different width. **Wider is the direction the pipeline
    /// produces** — the final prune shrinks the allele table, so copies not re-cut
    /// alongside it are left long, and every consumer indexes by [`AlleleId`], so the
    /// trailing entries would ride along unread. Both are in scope exactly here, which is
    /// why the check is here.
    ///
    /// If `per_sample` names no sample's call. A cohort has at least one sample and every
    /// sample gets a call, so an empty one is a locus that lost its calls rather than a
    /// locus with nothing to say.
    ///
    /// If `passes` is zero. Every locus takes at least one pass, so zero is a counter
    /// that was never incremented, and it would distort the distribution the pass cap is
    /// to be set from.
    ///
    /// If `repeat_tract` is set on a SNP/indel locus, **or missing at a repeat tract or
    /// bundle**. The
    /// tract ladder belongs to the repeat-tract prior; a generic locus seeds from a *frequency*
    /// spectrum, which has a ladder of its own with different rungs
    /// (`doc/devel/ng/spec/population_diversity.md` §3.4, §4.4), so a rung set there means one
    /// path's ladder was wired onto the other. The other direction is the newer half of the
    /// check: a tract whose record went missing reports no rung and no defaulted-cell counts,
    /// which reads downstream as *this call rested on nothing worth stating* rather than as a
    /// dropped field.
    ///
    /// If any call is [`SampleGenotypeCall::Missing`] at a repeat tract or bundle. **A
    /// repeat tract sets no sample aside** — a discovery round there can put back a length
    /// the cap cut, so no sample is locked out of the locus
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1, §9). This is the mirror of the seed
    /// marker's check above: one ruling belongs to each path, and each is refused on the
    /// other. What it catches is a sample silently dropped from a tract's output — spec §9
    /// says such a sample has no expected copies at all, so a `Missing` wrongly emitted
    /// there also removes that sample from the cohort's expected-copies denominator, and
    /// the locus's allele frequencies come out wrong with a well-formed record beside them.
    ///
    /// If a called genotype names an allele this locus was not called over. `Genotype::new`
    /// checks only that the multiset is non-empty, and the final prune renumbers every id
    /// above the one it drops ([`CandidateAlleles::admit`]), so a call not remapped
    /// alongside it arrives stale. An out-of-range id reaches the VCF's `GT` column as an
    /// index past the `ALT` list — an unparseable record, the same failure class
    /// [`CandidateAlleles::new`] refuses an empty reference allele for. **It catches the
    /// out-of-range half only:** an id that stays in range after a renumber names a
    /// different allele silently, which is an argument for the prune returning its
    /// remapping rather than for a wider check here.
    ///
    /// If `region` runs backwards. [`GenomeRegion`] is plain data with no constructor of
    /// its own, and its documentation says a caller that requires `start <= end` says so
    /// itself; a called locus is such a caller, since the region reaches the `POS` column
    /// and the writer's span arithmetic.
    ///
    /// If `site_quality` is above [`quality::MAX_SITE_QUALITY`]. That ceiling belongs to
    /// [`quality::score_uncorrected_site_quality`], which is the only thing that may fill
    /// this field, so a value above it did not come from there — and the field is the one
    /// the correction stage later overwrites in place
    /// (`doc/devel/ng/spec/calling_quality.md` §3.5, §5.3). The lower end needs no check:
    /// [`Phred`] is finite and non-negative by construction.
    ///
    /// If an artifact summary names the reference as its primary alternative, or names an
    /// allele this locus was not called over. Both tests weigh one alternative against the
    /// reference, so a locus with nothing to weigh carries **no summary at all** rather
    /// than one pointing at allele 0; and an id past the table is the same stale-after-a-
    /// renumber failure the genotype check above catches, arriving by the other door.
    ///
    /// **What is deliberately *not* checked: `converged` against `passes`.** A locus that
    /// hit the cap should report the cap, but the cap is run configuration this type does
    /// not see. Nor is `converged` required to imply more than one pass: the loop's first
    /// pass compares against a reads-only estimate made before it
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3), so settling on the first pass is a
    /// real outcome and not a comparison against nothing.
    #[allow(
        clippy::too_many_arguments,
        reason = "the architecture fixes this type's fields as a flat list \
                  (arch/calling_em_loop.md §2, and calling_quality.md §10's two additions); \
                  grouping them here to satisfy the lint would be a design change, not a \
                  refactor"
    )]
    pub fn new(
        region: GenomeRegion,
        alleles: CandidateAlleles,
        per_sample: Vec<SampleGenotypeCall>,
        cohort_expected_copies: ExpectedAlleleCopies,
        converged: bool,
        passes: u32,
        weakest_provenance: Provenance,
        repeat_tract: Option<RepeatTractProvenance>,
        site_quality: Phred,
        artifact_test_counts: Option<ArtifactTestCounts>,
    ) -> Self {
        assert_eq!(
            cohort_expected_copies.len(),
            alleles.len(),
            "a locus's expected allele copies run parallel to the alleles it was called \
             over: one entry per allele, reference first"
        );
        assert!(
            !per_sample.is_empty(),
            "a locus carries one call per sample and a cohort has at least one sample, so \
             a locus naming no call has lost them rather than having none to make"
        );
        assert!(
            passes > 0,
            "every locus takes at least one pass of the frequency loop, so a pass count \
             of zero is a counter that was never incremented"
        );
        // **A repeat *bundle* is on the demanding side of this, deliberately and unreachably.**
        // Written against `Generic` rather than against `Ssr(_)`, so a bundle would have to carry
        // a record — which it could not honestly build, having no single reference repeat count
        // and so no one allele length variation. Nothing can reach it: `tract_candidates` refuses a
        // bundle several frames earlier. When bundle scoring lands, this is one of the places
        // that has to answer for it, which is why it is stated rather than left to be discovered.
        assert_eq!(
            repeat_tract.is_some(),
            !matches!(alleles.kind(), LocusKind::Generic),
            "a {} locus {}: the tract ladder's rung and the defaulted-cell \
             counts belong to the repeat-tract path, and a SNP/indel locus is seeded from the \
             population curve's two moments, whose own ladder has different rungs",
            locus_kind_word(alleles.kind()),
            if repeat_tract.is_some() {
                "carries a repeat-tract record"
            } else {
                "carries no repeat-tract record"
            }
        );
        assert!(
            matches!(alleles.kind(), LocusKind::Generic)
                || per_sample.iter().all(|call| !call.is_missing()),
            "a repeat tract sets no sample aside, so a missing call at a {} locus is the \
             SNP/indel path's ruling wired onto the wrong path",
            locus_kind_word(alleles.kind())
        );
        assert!(
            per_sample
                .iter()
                .all(|call| call.genotype().is_none_or(|genotype| genotype
                    .alleles()
                    .iter()
                    .all(|id| alleles.bases_of(*id).is_some()))),
            "a called genotype names an allele this locus was not called over: the locus \
             holds {} alleles, and the final prune renumbers every id above the one it \
             drops, so a call not remapped alongside it is stale",
            alleles.len()
        );
        assert!(
            region.start <= region.end,
            "a called locus covers a stretch of reference, so its region cannot run \
             backwards: {region}"
        );
        assert!(
            site_quality.get() <= quality::MAX_SITE_QUALITY,
            "a site quality of {} is above the ceiling {} that \
             quality::score_uncorrected_site_quality caps at, so this number did not come \
             from it — and the field it is going into is the one the correction stage \
             overwrites in place",
            site_quality.get(),
            quality::MAX_SITE_QUALITY
        );
        if let Some(counts) = artifact_test_counts {
            assert!(
                !counts.primary_alternative.is_reference(),
                "the artifact tests compare an alternative allele against the reference, so \
                 a summary whose primary alternative *is* the reference has no comparison to \
                 make — a locus with nothing to test carries no summary at all rather than \
                 one naming allele 0"
            );
            assert!(
                alleles.bases_of(counts.primary_alternative).is_some(),
                "the artifact summary's primary alternative is allele {} and this locus \
                 holds {} alleles: the summary was pooled against a different allele table, \
                 or the final prune renumbered the table without re-cutting it",
                counts.primary_alternative.get(),
                alleles.len()
            );
        }
        Self {
            region,
            alleles,
            per_sample,
            cohort_expected_copies,
            converged,
            passes,
            weakest_provenance,
            repeat_tract,
            site_quality,
            artifact_test_counts,
        }
    }

    /// The alleles this locus was called over.
    #[inline]
    pub fn alleles(&self) -> &CandidateAlleles {
        &self.alleles
    }

    /// The cohort's expected allele copies, one entry per allele of [`Self::alleles`].
    #[inline]
    pub fn cohort_expected_copies(&self) -> &ExpectedAlleleCopies {
        &self.cohort_expected_copies
    }

    /// The site quality **as the worker wrote it**, before the artifact correction.
    ///
    /// **`pub(crate)`, and that is the enforcement of
    /// `doc/devel/ng/spec/calling_quality.md` §3.5 rather than a default choice of
    /// visibility**: the rule is that nothing between the worker and the correction stage
    /// may read this number, and the two that legitimately do — that stage, and this
    /// crate's tests — are both inside the crate. When the stage lands it overwrites the
    /// field in place and publishes a reader for the corrected value; there is no moment at
    /// which two qualities exist.
    ///
    /// **That stage now exists and this is read by it** (2026-09-01):
    /// [`evidence_for_output`](crate::run::records::evidence_for_output) takes this as the
    /// baseline, applies `correct_site_quality`, and the record carries the corrected number.
    /// The `dead_code` waiver this carried until then is spent.
    #[inline]
    #[must_use]
    pub(crate) fn uncorrected_site_quality(&self) -> Phred {
        self.site_quality
    }

    /// The nine pooled counts the artifact correction reads, or `None` where this locus
    /// gave its two tests nothing to weigh — see [`Self::artifact_test_counts`].
    #[inline]
    #[must_use]
    pub fn artifact_test_counts(&self) -> Option<ArtifactTestCounts> {
        self.artifact_test_counts
    }
}

#[cfg(test)]
mod tests;
