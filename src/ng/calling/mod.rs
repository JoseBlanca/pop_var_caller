//! Turning a locus's reads into a call: which alleles are in play, how probable each
//! sample's reads are under each candidate genotype, how likely each genotype is before
//! the reads are looked at, and what comes out. Steps 6 to 9 of the ng caller.
//!
//! **One folder rather than four, and a dependency decided it.** As four sibling
//! folders at the top of `src/ng/`, the three calling architecture documents each had
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
//! test-only: `genotype_table_parity`, which checks that table against production's
//! `GenotypeShape` value for value. It is a file of its own rather than a block inside
//! [`genotype_table`]'s own tests so that ng's single `use crate::var_calling::` — an
//! oracle, never a dependency of anything shipped — sits in one greppable place.
//!
//! The scalars other steps also name — `AlleleId`, `Phred`, `Genotype` — are not here
//! but in [`crate::ng::types`], with the rest of ng's shared vocabulary.

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
pub mod run_report;

pub use genotype_table::{GenotypeIdx, GenotypeTable, GenotypeTableView};
/// Re-exported for a different reason from [`genotype_table`]'s three, and the reason is
/// worth stating: **without a `use` of it that has to compile, deleting `pub mod
/// likelihood;` orphans the whole module silently.** The crate still builds, clippy is
/// still clean, and `cargo test --lib ng::calling::likelihood` reports `0 passed; ok` —
/// a green run naming a module that is no longer compiled. With this line the same
/// deletion is `error[E0432]: unresolved import`.
pub use likelihood::generic::{
    ERROR_SPREAD_BASES, ErrorSpreadTable, NO_ERROR_SPREAD, allele_is_compatible_with_partial,
    fill_error_spreads, genotype_log_likelihood_row,
};
pub use likelihood::{
    ContaminationMixture, ContaminationView, GenericEvidenceBuffer, GenericObservation,
    GenericRowScratch, GenericSampleEvidence, MIN_BASE_ERROR, MIN_CONTAMINANT_FREQUENCY,
    ReadGroupCalibration, ReadGroupParameters, SsrRowScratch, SsrSampleEvidence,
    fill_batch_allele_copies, fill_contaminant_allele_frequencies,
};

use crate::ng::calling::genotype_prior::SpectrumSeed;
use crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight;
use crate::ng::calling::quality::{ArtifactTestCounts, SiteQualityBuffers};
use crate::ng::locus_generation::{LocusKind, SsrDetail};
use crate::ng::parameter_estimation::Estimate;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::ng::parameter_estimation::joint::stratum_fits::{
    LengthSpectrum, LengthSpectrumRung, StratumFits,
};
use crate::ng::parameter_estimation::ssr::{RepeatCount, Stratum as SsrStratum, StratumKey};
use crate::ng::types::{
    AlleleId, BatchId, BatchOfEachReadGroup, BatchOfEachSample, ErrorRate, GenomeRegion, Genotype,
    InbreedingF, LogProb, Phred, Ploidy, ReadGroupId, SsrPeriod,
};
use std::collections::BTreeMap;
use std::num::NonZeroU32;

/// The alleles one locus is called over: the reference allele, and every alternative
/// the candidate step admitted — each stored as the bases it spells.
///
/// **The reference is allele 0 and is always present**, which is why there is no way
/// to build one of these without it. Every downstream branch tests against the
/// reference — REF against ALT in the VCF, the homozygous-reference genotype in the
/// prior — so a table that could lose it, or hold it somewhere other than the front,
/// would make every one of those tests conditional. That invariant is also why the
/// alleles are private where the architecture sketched a public `Vec`: a public one
/// admits `clear()`, `insert(0, …)` and `swap_remove(0)`, and nothing in the type
/// would notice.
///
/// **Owned rather than borrowed, because the table changes while a locus is called.**
/// A discovery round admits an allele the reads suggested and the first pass missed;
/// a final prune drops the ones no sample carries. Neither is possible against a
/// borrowed slice of the observations the loop was handed, and both are why the
/// alleles are stored here rather than pointed at
/// (`doc/devel/ng/arch/calling_em_loop.md` §2).
///
/// **The prune is not built here**, deliberately. Which alleles to drop is a policy
/// the calling loop owns, and dropping allele *k* renumbers every id above it — so it
/// cannot be `Vec::retain` handed out free, and every [`AlleleId`] minted before it
/// goes stale. The method that does it has to return the remapping, and it lands with
/// the step that needs one.
///
/// [`Self::kind`] routes two things. Which **row builder** scores the reads — that is,
/// which model turns a read into a per-allele likelihood, and a repeat tract's model
/// is not a SNP's. And which **seed** the genotype prior starts from, the prior's
/// opening guess at how the alleles are distributed before any read is read. A locus's
/// kind and the evidence handed to the loop must agree, and disagreeing is a caller
/// bug rather than a data condition.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CandidateAlleles {
    alleles: Vec<Box<[u8]>>,
    kind: LocusKind,
}

impl CandidateAlleles {
    /// A locus with only its reference allele, which is where every locus starts —
    /// the candidate step then admits what it found.
    ///
    /// # Panics
    ///
    /// On a reference allele spelling no bases. A locus's reference is the bases it
    /// covers — one for a SNP, several for an indel, the whole tract for a repeat —
    /// and an empty one would reach the VCF's `REF` column, which no writer in this
    /// repository rejects, as an unparseable record rather than a crash.
    pub fn new(reference: Box<[u8]>, kind: LocusKind) -> Self {
        assert!(
            !reference.is_empty(),
            "a locus's reference allele is the bases it spells — one for a SNP, several \
             for an indel, the whole tract for a repeat — and never nothing"
        );
        Self {
            alleles: vec![reference],
            kind,
        }
    }

    /// The bases of the reference allele. Infallible: allele 0 exists by construction.
    #[inline]
    pub fn reference(&self) -> &[u8] {
        &self.alleles[usize::from(AlleleId::REFERENCE.get())]
    }

    /// What kind of locus this is, and so which read model and which prior seed it
    /// gets. Read-only: it is fixed when the locus is generated, and a table whose
    /// kind could be changed afterwards could be re-routed to the wrong read model
    /// after its alleles were chosen by the right one.
    #[inline]
    pub fn kind(&self) -> &LocusKind {
        &self.kind
    }

    /// Every allele in id order, reference first — what a scorer walks.
    ///
    /// Yields the bases rather than the boxes that hold them, so the storage stays
    /// this type's business.
    #[inline]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &[u8]> {
        self.alleles.iter().map(Box::as_ref)
    }

    /// How many alleles this locus is called over — the reference plus every
    /// alternative admitted so far.
    #[inline]
    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    /// Never true: the reference is always present. Here because the compiler's lints
    /// ask for it beside [`Self::len`], and because saying so at the call site is the
    /// clearest statement of the invariant this type exists to hold.
    #[inline]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The bases one id names, or `None` if this table has no such allele.
    ///
    /// **Looked up rather than indexed, and that is the promise [`AlleleId`] makes.**
    /// An id carries no locus, so an id minted at one locus is a legal `u16` at the
    /// next: indexing would panic on a narrow table and, on a wider one, hand back a
    /// real but wrong allele without complaint. `AlleleId`'s own documentation says an
    /// out-of-range id is caught when the table is read, and this is where that
    /// happens.
    #[inline]
    pub fn bases_of(&self, id: AlleleId) -> Option<&[u8]> {
        self.alleles.get(usize::from(id.get())).map(Box::as_ref)
    }

    /// Admit an alternative allele and return the id that now names it.
    ///
    /// This is what the candidate step and a discovery round both call — the
    /// architecture's word, because discovered sequences enter the table through the
    /// same door as selected ones (`doc/devel/ng/arch/calling_em_loop.md` §6.2). The
    /// id is the table's previous length, so ids stay dense and in admission order,
    /// and the reference keeps index 0 however many rounds run.
    ///
    /// # Panics
    ///
    /// If the table already holds `u16::MAX + 1` alleles, since an [`AlleleId`] could
    /// not then name the new one.
    ///
    /// **ng has no allele cap of its own yet** — the candidate step (step 6) is where
    /// one will live — so this check is the only thing standing between a pathological
    /// locus and an id that wraps onto `AlleleId::REFERENCE`, which would silently make
    /// a discovered alternative *be* the reference for every downstream test. It is
    /// checked rather than trusted for that reason. Production's own caps are not an
    /// argument here: they belong to a caller ng shares no code with.
    pub fn admit(&mut self, allele: Box<[u8]>) -> AlleleId {
        let id = u16::try_from(self.alleles.len())
            .expect("a locus cannot hold more alleles than an AlleleId can name");
        self.alleles.push(allele);
        AlleleId(id)
    }
}

/// How many copies of each allele **the cohort** carries, summed over every sample's
/// own genotype probabilities — the quantity the loop feeds back to itself.
///
/// **The cohort's sum, not one sample's.** Both quantities exist, they have the same
/// shape, and the prior's leave-one-out term is a subtraction of the second from the
/// first (`doc/devel/ng/arch/calling_priors.md` §3), so confusing them is a wrong
/// prior rather than a compile error. A sample's own expected copies are a bare
/// `&[f64]` handed to the prior's builder; this type is the cohort total that travels
/// in the locus's output.
///
/// **Fractional, and never a call.** A sample that is 70% likely to be heterozygous
/// contributes 0.7 of a copy, not one or zero. That is what lets the loop work at
/// three reads a position, where no single sample's genotype is certain but the
/// cohort's allele frequencies still are; a version that rounded to called genotypes
/// first would throw away exactly the uncertainty that makes the low-coverage case
/// work (`doc/devel/ng/spec/calling_em_loop.md` §1.3).
///
/// **Parallel to the locus's [`CandidateAlleles`]**: entry *i* is allele *i*, so entry
/// 0 is the reference. [`Self::new`] takes the table for that reason — a copies vector
/// of the wrong length is not a value this type can hold, because a short one hands
/// every consumer that indexes by [`AlleleId`] a different allele's count.
///
/// It travels in the locus's output rather than being recomputed downstream, because
/// recomputing it from the called genotypes gives a *different* number — a call has
/// already discarded the uncertainty these counts still carry
/// (`doc/devel/ng/spec/calling_em_loop.md` §9).
#[derive(Clone, PartialEq, Debug)]
pub struct ExpectedAlleleCopies(Vec<f64>);

impl ExpectedAlleleCopies {
    /// # Panics
    ///
    /// If `copies` is not one entry per allele of `alleles`. The two are parallel by
    /// definition and this is the one place both are in scope, so it is where the
    /// pairing stops being a promise. A short vector is worse than a crash: every
    /// consumer indexing by [`AlleleId`] would read a different allele's count.
    ///
    /// If any entry is negative, infinite or a `NaN`. These are counts of genome
    /// copies; none of the three is a low-coverage answer, they are arithmetic that
    /// went wrong upstream. `NaN` is the one that must not pass silently — the loop
    /// stops when successive copy vectors stop moving, and a `NaN` is never equal to
    /// itself, so the locus would run to its pass cap and be emitted as unconverged
    /// with nothing saying why.
    pub fn new(copies: Vec<f64>, alleles: &CandidateAlleles) -> Self {
        assert_eq!(
            copies.len(),
            alleles.len(),
            "expected allele copies run parallel to the locus's allele table: one entry \
             per allele, reference first"
        );
        assert!(
            copies.iter().all(|c| c.is_finite() && *c >= 0.0),
            "expected allele copies are counts of genome copies: every entry must be \
             finite and at or above zero, got {copies:?}"
        );
        Self(copies)
    }

    /// The copies, allele by allele in id order, reference first.
    #[inline]
    pub fn copies(&self) -> &[f64] {
        &self.0
    }

    /// How many alleles these copies run parallel to.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Never true: these are parallel to an allele table, which always holds at least
    /// the reference. Here for the reason [`CandidateAlleles::is_empty`] is.
    #[inline]
    pub fn is_empty(&self) -> bool {
        false
    }
}

/// One run sample's evidence at one SNP/indel locus, together with the candidate step's
/// ruling on whether the sample can be genotyped there at all.
///
/// **Two facts that must travel together, because separating them is a silent join.** The
/// evidence says what the sample's reads showed; the ruling says whether the alleles those
/// reads earned survived candidate selection's cap. A run keeps both per sample and in the
/// run's own sample order, so holding them in one entry is what stops a later step pairing
/// sample *i*'s reads with sample *j*'s ruling.
///
/// **Why the ruling is not on [`GenericSampleEvidence`] itself**, where
/// `doc/devel/ng/arch/read_likelihoods.md` §2.1 sketches it. A sample the cap ruled
/// uncallable never reaches the read likelihood: it leaves the loop before the first pass
/// and is scored against no genotype at all
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0). A field on the evidence view would
/// therefore be one the row builder is handed and never reads, and the shipped
/// `GenericSampleEvidence` does not carry it.
#[derive(Copy, Clone, Debug)]
pub struct GenericLocusSample<'a> {
    /// What this sample's reads showed, as the SNP/indel row consumes it. A sample with no
    /// coverage at the locus gets [`GenericSampleEvidence::empty`], which scores every
    /// genotype alike and leaves the prior to decide — the right answer rather than a
    /// special case (`doc/devel/ng/spec/calling_em_loop.md` §7).
    pub evidence: GenericSampleEvidence<'a>,
    /// **Whether the allele cap cut a sequence this sample's own reads had earned**, copied
    /// from `UnmatchedSupport::genotype_must_be_missing`
    /// ([`allele_candidates::UnmatchedSupport::genotype_must_be_missing`]).
    ///
    /// True means the locus is no longer called over something this sample carries, so
    /// every genotype the caller could form for it is wrong. Such a sample is set aside
    /// before the first pass: it is scored against nothing, contributes nothing to the
    /// cohort's expected copies, and its call is emitted as
    /// [`SampleGenotypeCall::Missing`] (`doc/devel/ng/spec/calling_em_loop.md` §5.0).
    ///
    /// **Nothing downstream could derive it.** The pooled error mass of the cut sequences
    /// is the same number under every genotype and cancels, so the read likelihood cannot
    /// tell a sample that lost a real allele from one that had a handful of error reads
    /// dropped — and the second is nearly every sample at nearly every locus.
    pub genotype_must_be_missing: bool,
}

impl GenericLocusSample<'_> {
    /// Whether the locus can be called on this sample at all — the negation of
    /// [`genotype_must_be_missing`](Self::genotype_must_be_missing), and **the one spelling
    /// of it**.
    ///
    /// The predicate is asked in three places — when the scratch's rows are counted, when they
    /// are claimed, and when the artifact summary decides whose reads to pool — and a fourth
    /// site that wrote it out again would be invisible to anyone changing what *callable*
    /// means.
    #[inline]
    #[must_use]
    pub fn is_callable(&self) -> bool {
        !self.genotype_must_be_missing
    }
}

/// One locus's evidence, per sample, in whichever shape its path uses.
///
/// **This enum is the only place the two calling paths meet.** Below it everything is
/// path-pure: the discriminant chooses the row builder, and it is the same discriminant
/// that chose the candidates (`doc/devel/ng/arch/calling_em_loop.md` §2).
///
/// **Every variant's `per_sample` is one entry per sample of the run, in the run's sample
/// order** — not one entry per *covering* sample. That is the order every per-sample slice
/// in [`FrozenParameters`] is indexed by and the order [`LocusInference::per_sample`] comes
/// back in, and holding one order for the whole run is what makes the loop's fixed-order
/// sum over samples reproducible (`doc/devel/ng/spec/calling_em_loop.md` §8). The merge's
/// own `CohortObservation::per_sample` holds only the covering samples and is a different
/// list; converting between them is the input edge's work, not this type's.
#[derive(Copy, Clone, Debug)]
pub enum LocusEvidence<'a> {
    /// A SNP/indel locus: what each sample showed, and which samples the candidate cap
    /// already ruled uncallable.
    Generic {
        /// Where on the reference this locus is.
        region: GenomeRegion,
        /// One entry per run sample, in run order.
        per_sample: &'a [GenericLocusSample<'a>],
    },
    /// A repeat tract, or a repeat bundle: what each sample showed, and the tract's motif
    /// and flanks.
    ///
    /// **No sample is set aside here**, and the absence is structural rather than an
    /// oversight. A tract's discovery round can put back a length the cap cut, so a sample
    /// that lost one is not locked out of the locus for the rest of its calling
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1).
    Ssr {
        /// Where on the reference this locus is.
        region: GenomeRegion,
        /// One entry per run sample, in run order.
        per_sample: &'a [SsrSampleEvidence<'a>],
        /// The tract's repeat unit and its two flanks — what the repeat read model scores
        /// against.
        detail: &'a SsrDetail,
        /// **How many whole repeats each candidate carries**, parallel to the locus's
        /// [`CandidateAlleles`] and in the same order.
        ///
        /// **It is not derivable from the bases, which is the whole reason it travels.** An
        /// interrupted tract — one whose repeat is broken by a substitution — holds fewer whole
        /// repeats than its length suggests, and the slippage fit is keyed by the count rather
        /// than by the length (`doc/devel/ng/spec/read_likelihoods.md` §4.4). Two candidates of
        /// equal length can therefore be drawn from different strata.
        ///
        /// **Its producer is repeat-tract candidate selection, which is unwritten**
        /// ([`candidate_alleles_ssr.md`](../../../doc/devel/ng/impl_plan/candidate_alleles_ssr.md)),
        /// so until that lands every caller supplies it — and a test that supplies one must say
        /// so, because a supplied candidate set read as a selected one is a claim about a step
        /// that does not exist.
        candidate_repeat_counts: &'a [NonZeroU32],
    },
}

impl<'a> LocusEvidence<'a> {
    /// A SNP/indel locus's evidence.
    ///
    /// # Panics
    ///
    /// If `per_sample` names no sample. A run has at least one sample and every sample gets
    /// an entry — an empty one is a locus whose evidence was lost on the way in, not a
    /// locus nobody covered, which is one entry per sample of
    /// [`GenericSampleEvidence::empty`].
    #[must_use]
    pub fn generic(region: GenomeRegion, per_sample: &'a [GenericLocusSample<'a>]) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the SNP/indel evidence at {region} names no sample: a locus carries one entry \
             per sample of the run and a run has at least one sample, so an empty list is \
             evidence that went missing rather than a locus nobody covered — a sample with \
             no reads gets GenericSampleEvidence::empty()"
        );
        Self::Generic { region, per_sample }
    }

    /// A repeat tract's or repeat bundle's evidence.
    ///
    /// # Panics
    ///
    /// If `per_sample` names no sample, for the reason [`Self::generic`] gives.
    #[must_use]
    pub fn ssr(
        region: GenomeRegion,
        per_sample: &'a [SsrSampleEvidence<'a>],
        detail: &'a SsrDetail,
        candidate_repeat_counts: &'a [NonZeroU32],
    ) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the repeat-tract evidence at {region} names no sample: a locus carries one \
             entry per sample of the run and a run has at least one sample, so an empty \
             list is evidence that went missing rather than a locus nobody covered"
        );
        assert!(
            !candidate_repeat_counts.is_empty(),
            "the repeat tract at {region} supplied no candidate repeat counts: a tract is \
             called over at least its reference length, and the count is what keys the \
             slippage fit. Their number is checked against the candidate table by \
             `assert_matches_locus_and_run`, which is where the two meet"
        );
        Self::Ssr {
            region,
            per_sample,
            detail,
            candidate_repeat_counts,
        }
    }

    /// Where on the reference this locus is.
    #[inline]
    #[must_use]
    pub fn region(&self) -> GenomeRegion {
        match self {
            Self::Generic { region, .. } | Self::Ssr { region, .. } => *region,
        }
    }

    /// How many samples the run has — the length of every per-sample list at this locus.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        match self {
            Self::Generic { per_sample, .. } => per_sample.len(),
            Self::Ssr { per_sample, .. } => per_sample.len(),
        }
    }

    /// **How many of the run's samples this locus can actually be called on** — the rest were
    /// ruled uncallable by the candidate step, for having earned a sequence the allele cap cut
    /// (`doc/devel/ng/spec/candidate_alleles.md` §4.1).
    ///
    /// **It counts the run's samples and not the locus's covering ones**, which is the whole of
    /// what a caller has to get right here. A sample that covered nothing is **callable**: its
    /// evidence is empty, an empty sum is zero, every genotype scores alike and the prior
    /// decides alone. So this reaches zero only where **every sample of the run** covered the
    /// locus and every one of them lost a sequence its own reads had earned — which in a large
    /// cohort at low coverage is close to never, because some sample almost always covers
    /// nothing.
    ///
    /// **A locus can still have none.** A run must count such a locus and carry on rather than
    /// die on it (owner's ruling, 2026-09-01), so whoever drives the calling asks this before
    /// handing the locus to a genotyper — and
    /// [`LocusGenotyper::call_locus`](inference::LocusGenotyper::call_locus) keeps its
    /// precondition that there is somebody to call, since its scratch cannot be prepared for
    /// no rows.
    ///
    /// **A repeat tract sets no sample aside**, so on that path this is always the run's whole
    /// sample count (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1).
    #[inline]
    #[must_use]
    pub fn callable_sample_count(&self) -> usize {
        match self {
            Self::Generic { per_sample, .. } => per_sample
                .iter()
                .filter(|sample| sample.is_callable())
                .count(),
            Self::Ssr { per_sample, .. } => per_sample.len(),
        }
    }

    /// Which calling path this evidence is on, in the words the panic messages use.
    #[inline]
    fn path_word(&self) -> &'static str {
        match self {
            Self::Generic { .. } => "SNP/indel",
            Self::Ssr { .. } => "repeat-tract",
        }
    }

    /// Check this evidence against the locus it claims to describe and the run it claims to
    /// come from — **the ordering contract, made a runtime fact rather than a promise.**
    ///
    /// Three things are checked, and each is a caller bug whose symptom is a wrong genotype
    /// rather than a crash:
    ///
    /// - **the path.** A repeat tract's evidence handed to a locus whose alleles say
    ///   SNP/indel would be scored by the wrong read model, which is a different likelihood
    ///   at every sample with nothing failing
    ///   (`doc/devel/ng/arch/calling_em_loop.md` §2). A repeat **bundle** goes to the
    ///   repeat path with a tract, which is how every other consumer of
    ///   [`LocusKind`] groups the two. *(What happens to a bundle after that is a
    ///   different question and the calling loop has no answer yet: the tract's candidate
    ///   table refuses one by name, because nothing scores a bundle. No producer emits one
    ///   into calling.)*
    /// - **the cohort.** One run-wide sample order indexes this evidence, every per-sample
    ///   slice of `parameters`, and the calls that come back. Two lists of different lengths
    ///   are two different orders, and a positional join between them silently pairs one
    ///   sample's reads with another's inbreeding coefficient.
    /// - **the repeat counts, at a tract only.** They run parallel to the candidate table,
    ///   and a table with more candidates than counts would pair one candidate's bases with
    ///   another candidate's stratum.
    ///
    /// **The locus's genotype table is not checked here**, because this method never sees
    /// it. That check has its own home, at the one point the shape is fixed:
    /// [`CallingScratch::prepare_for_locus`].
    ///
    /// # Panics
    ///
    /// On any of the three disagreements, in release as well as debug: each is a comparison
    /// of two integers or two discriminants, against a defect that would otherwise reach the
    /// output as a genotype.
    pub fn assert_matches_locus_and_run(
        &self,
        alleles: &CandidateAlleles,
        parameters: &FrozenParameters<'_>,
    ) {
        let paths_agree = matches!(
            (self, alleles.kind()),
            (Self::Generic { .. }, LocusKind::Generic)
                | (Self::Ssr { .. }, LocusKind::Ssr(_) | LocusKind::SsrBundle)
        );
        assert!(
            paths_agree,
            "the evidence at {} is on the {} path and its allele table is a {} locus, so \
             the two belong to different loci — the discriminant that chose the candidates \
             is the one that chooses the read model",
            self.region(),
            self.path_word(),
            locus_kind_word(alleles.kind())
        );
        // **The join a repeat tract adds, and it is the one no type enforces.** Its repeat
        // counts run
        // parallel to the candidate table, and a table with more candidates than counts would
        // pair one candidate's bases with another's stratum — a wrong slippage model per
        // candidate, with a well-formed genotype coming back.
        if let Self::Ssr {
            candidate_repeat_counts,
            ..
        } = self
        {
            assert_eq!(
                candidate_repeat_counts.len(),
                alleles.len(),
                "the repeat tract at {} is called over {} candidates and {} repeat counts \
                 arrived, so one of the two belongs to a different locus",
                self.region(),
                alleles.len(),
                candidate_repeat_counts.len()
            );
        }
        assert_eq!(
            self.sample_count(),
            parameters.sample_count(),
            "the evidence at {} covers {} samples and the run's frozen parameters cover {}, \
             so one of them is indexed by a different sample order",
            self.region(),
            self.sample_count(),
            parameters.sample_count()
        );
    }
}

/// What kind of locus this is, in one word, for a panic message.
///
/// **Not `{:?}` on the kind**, which for a repeat tract prints both flanks as decimal byte
/// arrays — a screenful of digits for a fact the discriminant alone carries.
fn locus_kind_word(kind: &LocusKind) -> &'static str {
    match kind {
        LocusKind::Generic => "SNP/indel",
        LocusKind::Ssr(_) => "repeat-tract",
        LocusKind::SsrBundle => "repeat-bundle",
    }
}

/// Everything the parameter pre-pass froze, borrowed for the whole run.
///
/// **Assembled once per run and never written during calling.** Every error rate,
/// contamination fraction and inbreeding coefficient arrives fitted and leaves unchanged;
/// what a locus may move is its own allele frequencies, and — at a repeat tract with the
/// re-fit switched on — its own slippage numbers
/// (`doc/devel/ng/spec/calling_em_loop.md` §5). Each field is another document's to define;
/// this type only gathers them so that the loop borrows them once, at the boundary where
/// calling begins.
///
/// **Two different axes are in here and they are not interchangeable**, which is why every
/// field says which one it is on. The calibration and the contamination views are per
/// **read group** — a sample sequenced from two libraries has two of each. The inbreeding
/// coefficients are per **sample**, in the run's sample order, the same order
/// [`LocusEvidence`] and [`LocusInference::per_sample`] use.
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
    /// [`RunParameters::view`](run_parameters::RunParameters::view), which passes on whatever
    /// its own parameters hold.
    repeat_tract_outlier_weight: RepeatTractOutlierWeight,
}

impl<'a> FrozenParameters<'a> {
    /// Gather the run's frozen parameters for a run the fit **did** measure contamination
    /// in, with what can be checked between them checked.
    ///
    /// - `calibration_by_read_group` — one entry per read group, indexed by
    ///   [`ReadGroupId`](crate::ng::types::ReadGroupId).
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
        }
    }

    /// One calibration per read group, indexed by
    /// [`ReadGroupId`](crate::ng::types::ReadGroupId).
    #[inline]
    #[must_use]
    pub fn calibration_by_read_group(&self) -> &'a [ReadGroupCalibration] {
        self.calibration_by_read_group
    }

    /// **Score repeat tracts under a supplied outlier weight rather than the inherited one.**
    ///
    /// **A builder rather than an eleventh argument to the two constructors**, because no fit
    /// produces this number and almost every call site wants the compiled-in value. The
    /// exception is [`RunParameters::view`](run_parameters::RunParameters::view), whose two
    /// arms pass on what a parameters file supplied; the other **25** call sites in this
    /// repository — 27 in all, counted with
    /// `grep -rn "FrozenParameters::new(\|FrozenParameters::uncontaminated("` — would each
    /// have gained an argument naming the same constant.
    ///
    /// **A caller that omits it is not silently wrong**: it gets the inherited constant, which
    /// is what the run would have scored under.
    #[inline]
    #[must_use]
    pub fn with_repeat_tract_outlier_weight(mut self, weight: RepeatTractOutlierWeight) -> Self {
        self.repeat_tract_outlier_weight = weight;
        self
    }

    /// **The share of a repeat tract's reads charged to none of its candidate alleles**, with
    /// its warrant — what
    /// [`TractScoringFits`](inference::repeat_tract_parameters::TractScoringFits) reads once a
    /// locus and hands the scoring row.
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
    /// ([`BatchOfEachReadGroup`](crate::ng::types::BatchOfEachReadGroup)).
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
    /// **It takes a `usize` where [`Self::batch_of_read_group`] takes a [`ReadGroupId`], and
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
    /// **It always answers**, at one of three rungs — the stratum's own fit, its motif period's
    /// pooled tracts, or a flat shape at a stated concentration — so a repeat tract is never
    /// left without a prior. Which rung it was is on the returned value and belongs in the
    /// run's output.
    #[inline]
    #[must_use]
    pub fn ssr_length_spectrum_at(
        &self,
        period: SsrPeriod,
        reference_repeats: RepeatCount,
    ) -> LengthSpectrum<'a> {
        self.ssr_slippage_fits
            .length_spectrum_at(period.get(), u64::from(reference_repeats.get()))
    }

    /// How many copies of the genome a sample carries.
    #[inline]
    #[must_use]
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }
}

/// What every scratch slot holds between [`CallingScratch::prepare_for_locus`] and the pass
/// that writes it: **not a number**, so a slot some pass forgot to fill fails the next
/// arithmetic check it is handed to rather than reaching a genotype as a plausible zero.
///
/// Zero is exactly the mistake to avoid: an expected copy count of zero, a log-prior of
/// zero and a concentration of zero are all legal values a pass could have written, so a
/// zero-filled buffer and a correctly-filled one are indistinguishable. It is the same
/// argument `SsrRowScratch::prepare_emissions` makes for its own fill value.
pub const UNWRITTEN_SCRATCH_VALUE: f64 = f64::NAN;

/// **What building one locus's genotype-likelihood table cost** — the instrument
/// `doc/devel/ng/spec/calling_em_loop.md` §13's test 5 reads.
///
/// **The property it exists to pin is that the expensive part does not grow with the pass
/// count.** A version that recomputed it per pass would give identical genotypes, only slower,
/// so nothing but a counter can tell the two apart (§2, §8).
///
/// # The build has two halves, and this has to say which it is counting
///
/// **With contamination on at an ordinary site, a caller may no longer cache a whole row across
/// iterations**: `q(o)`, the contaminating population's frequency for the allele an observation
/// shows, is the locus's own number and moves with the loop
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.6, corrected 2026-08-24). **A repeat tract's third
/// term does not move** — it is the fit's length spectrum for the tract's stratum, frozen before
/// calling starts (§4.5.1) — so a contaminated tract caches its whole row like an uncontaminated
/// one.
///
/// What it may still cache — and what the first three fields count — is the **emission**: the
/// answer to how one copy of one allele produced one observed sequence. It reads no frequency,
/// so it stays computed once per `(sample, observation, candidate)` per locus (§6.1). What runs
/// again is the per-genotype **assembly** of the row from those emissions, which the last two
/// fields count: one multiply and one add inside a logarithm the row was taking anyway.
///
/// **⚠ "Expensive" is the spec's word for the general case, and on the SNP/indel path with no
/// partial reads it is the wrong way round.** Counted on this module's own three-sample
/// contaminated fixture — 3 rows × 3 observations × 3 alleles, diploid, so 6 genotypes and 9
/// assemblies — the assemblies cost about 486 multiply-adds and 54 logarithms against the
/// emission fill's 9 charged-error calls. What the emission genuinely is expensive for is the
/// per-`(partial read, candidate)` byte comparison, and, on the repeat-tract path, an alignment
/// per `(observation, candidate)`. **What the split buys where the emission is cheap is not
/// arithmetic but the count**: `emission_evaluations` cannot grow with the pass count, so the
/// invariant §13's test 5 pins stays a statement about the model rather than about a fixture.
///
/// **An uncontaminated run assembles once too**, because with no frequency in the formula the
/// assembled row is the same value at every pass — so it keeps exactly the behaviour it had
/// before the split.
///
/// **`emission_evaluations` is accumulated per row fill from the shapes that fill was
/// handed, not counted inside the emission.** The repeat-tract path has an emission seam a
/// counting model could wrap; the SNP/indel path has none — its per-read term is inside the
/// fill. So this counts what the driver *asked for*, which is what catches a rebuild,
/// and does not catch a fill that evaluates the wrong number of terms. Saying which
/// is which is the point of the field names.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub(crate) struct EmissionCost {
    /// **The expensive half.** How many times the locus's frequency-free emissions were
    /// computed for the whole cohort — **one**, whatever the pass count and whether or not the
    /// run fitted a contamination fraction.
    pub(crate) emission_builds: u64,
    /// How many times one sample's emissions were filled — one per row per emission build.
    pub(crate) emission_row_fills: u64,
    /// Summed over those fills, `observations × candidates`: the emission evaluations the
    /// fills were asked for. **`Σ_s`, not a three-way product** — a fixture whose samples
    /// hold equal observation counts is the one shape that hides the difference (§13 test 5).
    pub(crate) emission_evaluations: u64,
    /// **The cheap half.** How many times the whole `rows × genotypes` table was assembled
    /// from those emissions.
    ///
    /// **One wherever the assembled row reads no allele frequency**, whatever the pass count —
    /// which is every uncontaminated locus, **and every repeat tract**, whose third mixture term
    /// is frozen before calling starts (spec §4.5.1). At a **contaminated ordinary site** it is
    /// one for the initialisation pass, one at the head of each pass, and one more against the
    /// settled frequencies before the final pass.
    pub(crate) table_assemblies: u64,
    /// How many rows were assembled, counted one at a time — `rows × table_assemblies` on a
    /// whole table, and less where an assembly stopped short.
    pub(crate) row_assemblies: u64,
}

/// Every buffer one locus's calling fills — **allocated once per worker and reused at every
/// locus**, so calling a locus costs no allocation of its own.
///
/// **The reason is measured, not stylistic.** Production lifted exactly these buffers out of
/// its own iteration after a profile put the allocator's own self-time at about one cycle in
/// six (`src/var_calling/posterior_engine.rs`; `doc/devel/ng/spec/calling_em_loop.md` §8).
///
/// **Candidate selection's buffers live here too** — the per-allele running totals it
/// accumulates and its survivor list — rather than standing alone as a second per-worker
/// object. The same worker selects a locus's alleles and then calls it, so a second
/// allocation would buy nothing (`doc/devel/ng/arch/candidate_alleles.md` §2.4). Nothing
/// about their shape changed on the way in.
///
/// **Sized by [`prepare_for_locus`](Self::prepare_for_locus), one call per locus**, and
/// every buffer comes back holding [`UNWRITTEN_SCRATCH_VALUE`]. Until that call the scratch
/// has no shape, and every accessor refuses rather than handing back an empty slice: a fold
/// over a zero-length buffer writes nothing and sums to a plausible `0.0`, which is the
/// failure the fill value exists to prevent, reached by the door the fill cannot cover. The
/// two row scratches and the candidate-selection buffers own their own sizing.
///
/// **Amortised, not allocation-free.** A buffer grows when a locus is wider than every locus
/// this worker has met and never shrinks back, so a worker ends up holding its widest
/// locus's shape for its lifetime.
///
/// `SsrEmissionScratch` is the repeat-tract emission model's own working memory — the
/// associated `Scratch` of whichever [`SsrEmissionModel`](likelihood::ssr_emission::SsrEmissionModel)
/// the run scores tracts with. **It has no default**, so every construction names the model
/// it is for, the way the `SsrRowScratch` it wraps already does: that seam exists for a
/// bake-off between two emission models, and a defaulted parameter would let a run switch
/// models with no call site changing.
#[derive(Default, Debug)]
pub struct CallingScratch<SsrEmissionScratch> {
    /// Candidate selection's per-allele running totals and its survivor list.
    candidate_selection: allele_candidates::SelectionScratch,
    /// `samples × genotypes`, sample-major: how probable each sample's reads are under each
    /// candidate genotype — the **genotype** likelihood, written `Lg` in
    /// `doc/devel/ng/spec/read_likelihoods.md` §1, as against `Lr`, which is one read
    /// against one allele.
    ///
    /// **Assembled from the emissions below, and how often depends on whether the assembled row
    /// reads an allele frequency.** With no fraction fitted it does not, so the row is the same
    /// value at every pass and is written once per locus. With one fitted, an **ordinary site's**
    /// row holds `q(o)` — the cohort's own frequency for the allele an observation shows — which
    /// moves with the loop, so it is written again at every pass; a **repeat tract's** third term
    /// is the fit's length spectrum for its stratum, frozen before calling starts, so its row is
    /// written once too (`doc/devel/ng/spec/read_likelihoods.md` §3.6, §4.5.1, §6.1). What is
    /// written once in every case is the emissions.
    genotype_likelihoods: Vec<LogProb>,
    /// One entry per genotype: the sample being scored, its log-prior.
    prior_row: Vec<LogProb>,
    /// One entry per genotype: the sample being scored, its posterior over genotypes.
    posterior_row: Vec<f64>,
    /// One entry per allele: the locus's **seed concentration** — one positive number per
    /// allele, read as chromosomes the prior behaves as though it had already seen
    /// ([`genotype_prior`]'s module doc). Built once per locus.
    seed_concentration: Vec<f64>,
    /// One entry per allele: the sample being scored, its own concentration — the seed plus
    /// what the *other* samples showed here.
    sample_concentration: Vec<f64>,
    /// One entry per allele: the genotype prior's own working space, which
    /// `PriorRow::new` borrows mutably while it reads the concentration. **One per allele,
    /// not one per genotype** — the buffer beside it, `prior_row`, is the other, and at a
    /// diploid biallelic locus they are 2 and 3, so each is a legal length for the other.
    /// `PriorRow::new` calls this same buffer `per_allele_scratch`.
    prior_per_allele_workspace: Vec<f64>,
    /// One entry per allele: the **cohort's** expected allele copies as this pass has them,
    /// summed over every sample.
    cohort_expected_copies: Vec<f64>,
    /// One entry per allele: what the previous pass's cohort copies were — the convergence
    /// comparison.
    previous_cohort_expected_copies: Vec<f64>,
    /// `samples × alleles`, sample-major: each **sample's own** expected allele copies. The
    /// prior's leave-one-out term is the cohort's total minus this.
    per_sample_expected_copies: Vec<f64>,
    /// `samples × (ploidy + 1)`, sample-major: the site quality's first step, each sample's
    /// log-likelihood of carrying exactly `c` non-reference copies
    /// (`doc/devel/ng/spec/calling_quality.md` §5.2).
    copy_count_log_likelihoods: Vec<f64>,
    /// The site quality's count axis, and the buffer its fold alternates with. Both are
    /// `ploidy + samples × ploidy + 1` long — **padded by `ploidy` on the left**, so a
    /// copy-count tap reads `[padding − copies]` without a bounds test in the quadratic inner
    /// loop.
    allele_count_distribution: Vec<f64>,
    /// The other half of that alternation.
    allele_count_distribution_next: Vec<f64>,
    /// `samples × ploidy + 1`: the fold's result back in the log domain, then the
    /// unnormalised log-posterior over the cohort's allele count once the prior is applied.
    log_allele_count_distribution: Vec<f64>,
    /// One entry per allele: how many reads that allele drew, pooled over the samples the
    /// locus was called on — the walk the artifact summary picks its **primary
    /// alternative** from (`doc/devel/ng/spec/calling_quality.md` §6.3).
    ///
    /// **The one buffer here with no [`UNWRITTEN_SCRATCH_VALUE`] sentinel**, and it is
    /// counts rather than the sentinel's type that decide it: reads are whole, so this is a
    /// `u64` and a `NaN` is not expressible in it. What replaces the sentinel is that the
    /// only thing that reads this buffer is the same call that zeroes and fills it, in that
    /// order, so a stale entry from the previous locus cannot be read. A zero here is a
    /// real answer — an allele no read reached — which is why nothing downstream may treat
    /// it as *unwritten*.
    pooled_allele_reads: Vec<u64>,
    /// `genotypes × alleles`, genotype-major: how far each allele's own error mass is spread
    /// across the locus's others under each genotype — the SNP/indel row's
    /// `fill_error_spreads` table.
    ///
    /// **Filled once per locus rather than once per sample**: how many alleles a wrong read
    /// could have shown is a property of the candidate sequences and of the genotype being
    /// scored, and of nothing a sample showed.
    error_spreads: Vec<f64>,
    /// One entry per **scratch row**: which sample of the run that row holds.
    ///
    /// **The scratch is sized for the samples the locus is called on, not for the run**, so
    /// this is the map back. A sample the candidate step ruled uncallable takes no part in
    /// the loop at all (`doc/devel/ng/spec/calling_em_loop.md` §5.0) and has no row here —
    /// which is what makes the cohort's expected copies, the convergence denominator and the
    /// site quality's count axis all run over the same cohort without any of them being told
    /// to skip anything (§9).
    run_sample_of_each_row: Vec<usize>,
    /// One entry per scratch row: that sample's inbreeding coefficient, in run order.
    ///
    /// A compacted copy of [`FrozenParameters::inbreeding_coefficient_by_sample`], because
    /// the loop walks one coefficient per row and the run's slice has an entry for every
    /// sample including the ones with no row.
    inbreeding_coefficient_by_row: Vec<InbreedingF>,
    /// What the last `prepare_for_locus` reset and the table build has counted since — spec
    /// §13 test 5's instrument.
    emission_cost: EmissionCost,
    /// **One SNP/indel emission cache per scratch row**, and one per row rather than one
    /// reused, because with contamination on the loop reads them again after the next sample
    /// has been scored: the emissions are filled once per locus and the rows are assembled from
    /// them once per pass (`doc/devel/ng/spec/read_likelihoods.md` §6.1).
    ///
    /// **Grown to the locus's row count and never shrunk**, like every other buffer here —
    /// and, like the repeat-tract row's, each one sizes itself per sample, which is why they
    /// are outside [`buffer_fingerprints`](Self::buffer_fingerprints).
    generic_rows: Vec<GenericRowScratch>,
    /// `batches × alleles`, batch-major: how many copies of each allele the samples of each
    /// sequencing batch are expected to carry, this locus, at the loop's current estimate —
    /// the first half of spec §3.6's `q(o)`, and the half that does not depend on which
    /// sample is being scored. Refilled **once per pass**.
    ///
    /// Empty wherever nothing reads it — a run the fit found no contamination in, and **every
    /// repeat tract of any run**, whose mixture reads no per-batch frequency.
    batch_allele_copies: Vec<f64>,
    /// `batches × alleles`, batch-major: the contaminant allele frequencies **one sample** is
    /// scored against — the copies above turned into a distribution with that sample's own
    /// copies taken out of its own batch. Refilled **once per sample per pass**, because the
    /// subtraction is that sample's.
    ///
    /// Empty wherever nothing reads it — see [`batch_allele_copies`](Self::batch_allele_copies).
    contaminant_allele_frequencies: Vec<f64>,
    /// `run samples × alleles`, sample-major: the expected copies above, scattered from the
    /// scratch rows back onto the **run's** sample axis, with zero at every sample this locus
    /// has no row for.
    ///
    /// **It exists because the batching is the run's and the rows are the locus's.** A batch
    /// whose samples were all ruled uncallable here has no row at all, and summing over rows
    /// would leave that batch's row of the copy table unwritten — which the copy fill refuses,
    /// by name, as a batching that does not describe the run. Scattering onto the run's axis
    /// gives such a batch a row of zeros instead, which is what the M-step already does with an
    /// uncallable sample: it contributes nothing.
    ///
    /// Empty wherever nothing reads it — see [`batch_allele_copies`](Self::batch_allele_copies).
    expected_copies_by_run_sample: Vec<f64>,
    /// How many sequencing batches the two tables above are sized for. **Zero means the
    /// contaminant tables were not prepared**, which is what an uncontaminated locus leaves
    /// behind and what every contaminant accessor refuses.
    batch_count: usize,
    /// How many of the run's samples [`expected_copies_by_run_sample`](Self::expected_copies_by_run_sample)
    /// is sized for. Zero alongside `batch_count`.
    run_sample_count: usize,
    /// The repeat-tract row's own scratch, including the emission model's.
    ssr_row: SsrRowScratch<SsrEmissionScratch>,
    /// **The tract's fitted scoring parameters, one cell per `(read group, candidate)`** —
    /// gathered once per repeat tract and read by every row of it
    /// ([`inference::repeat_tract_parameters`](super::calling::inference::repeat_tract_parameters)).
    ///
    /// **Held here rather than built per locus** for the reason every buffer in this type is:
    /// its seven vectors are cleared and refilled at each tract, so a worker allocates them once
    /// for a whole run. What cannot be held here is the *contexts* built from it — they borrow
    /// these vectors, so a struct owning both would be self-referential — and those are one
    /// allocation per tract, which is the cost the repeat path pays and the SNP/indel path does
    /// not.
    tract_fits: crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits,
    /// How many **rows** the buffers above are sized for. Zero means never prepared.
    ///
    /// **Rows, not the run's samples, and the difference is this type's whole shape**: a
    /// locus is prepared for the samples it is *called on*, which is the run's samples minus
    /// the ones the candidate step ruled uncallable
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0). `LocusEvidence::sample_count` and
    /// `FrozenParameters::sample_count` both answer the other question, and the two numbers
    /// appear within a few lines of each other in the final pass.
    row_count: usize,
    /// How many genotypes they are sized for.
    genotype_count: usize,
    /// How many alleles they are sized for.
    allele_count: usize,
}

impl<SsrEmissionScratch> CallingScratch<SsrEmissionScratch> {
    /// Size every buffer for one locus and fill it with [`UNWRITTEN_SCRATCH_VALUE`], so
    /// nothing survives from the last one.
    ///
    /// **The shape comes from the locus's own allele table and genotype table, and from
    /// nowhere else.** Taking the two objects rather than two integers is what stops the
    /// allele count and the genotype count being swapped: at a diploid biallelic locus they
    /// are 2 and 3, and every buffer would still be a legal length.
    ///
    /// # Panics
    ///
    /// If the genotype table does not index the locus's alleles. This is the first of the
    /// three caller bugs `doc/devel/ng/spec/calling_em_loop.md` §8 names as assertions, and
    /// this is the one point where the locus's shape is fixed, so it is where the two are
    /// compared. A table built before a discovery round admitted an allele is exactly how
    /// the two come apart; with a table one allele narrow, every per-allele buffer is sized
    /// for the old count and the locus is called over a set that silently excludes the
    /// allele just admitted.
    ///
    /// If `row_count` is zero. Every locus is called on at least one sample, so a locus
    /// prepared for no rows is a run whose sample order went missing rather than a locus
    /// nobody covered.
    ///
    /// If `rows × genotypes` or `rows × alleles` overflows a `usize`.
    /// `row_count` is **one row per sample this locus is called on** — not the run's sample
    /// count, which at a locus with an uncallable sample is larger.
    pub fn prepare_for_locus(
        &mut self,
        row_count: usize,
        alleles: &CandidateAlleles,
        genotypes: &GenotypeTableView<'_>,
    ) {
        assert_eq!(
            genotypes.allele_count(),
            alleles.len(),
            "the genotype table indexes {} alleles and this locus is called over {}, so the \
             table was built for a different allele set — a discovery round that admitted \
             an allele needs its table rebuilt with it",
            genotypes.allele_count(),
            alleles.len()
        );
        assert!(
            row_count > 0,
            "every locus is called on at least one sample, so a locus prepared for no rows is \
             a run whose sample order went missing"
        );
        let genotype_count = genotypes.genotype_count();
        let allele_count = genotypes.allele_count();
        let table_len = row_count.checked_mul(genotype_count).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows over {genotype_count} genotypes needs a \
                 genotype-likelihood table longer than a usize can index"
            )
        });
        let copies_len = row_count.checked_mul(allele_count).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows over {allele_count} alleles needs a \
                 per-sample copies table longer than a usize can index"
            )
        });

        let unwritten = LogProb(UNWRITTEN_SCRATCH_VALUE);
        resize_and_fill(&mut self.genotype_likelihoods, table_len, unwritten);
        resize_and_fill(&mut self.prior_row, genotype_count, unwritten);
        resize_and_fill(
            &mut self.posterior_row,
            genotype_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.seed_concentration,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.sample_concentration,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.prior_per_allele_workspace,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.cohort_expected_copies,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.previous_cohort_expected_copies,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.per_sample_expected_copies,
            copies_len,
            UNWRITTEN_SCRATCH_VALUE,
        );

        // The site quality's four buffers. **Sized from the genotype table's ploidy**, which
        // is the one place this locus's ploidy is stated, so they cannot disagree with the
        // table the fold walks (`doc/devel/ng/spec/calling_quality.md` §5.2, §9).
        let ploidy = usize::from(genotypes.ploidy().get());
        let copy_count_table_len = row_count.checked_mul(ploidy + 1).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows at ploidy {ploidy} needs a copy-count log-likelihood \
                 table longer than a usize can index"
            )
        });
        let largest_count = row_count.checked_mul(ploidy).unwrap_or_else(|| {
            panic!(
                "a cohort of {row_count} rows at ploidy {ploidy} carries more \
                 chromosomes than a usize can count"
            )
        });
        resize_and_fill(
            &mut self.copy_count_log_likelihoods,
            copy_count_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        let padded_axis = ploidy + largest_count + 1;
        resize_and_fill(
            &mut self.allele_count_distribution,
            padded_axis,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.allele_count_distribution_next,
            padded_axis,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.log_allele_count_distribution,
            largest_count + 1,
            UNWRITTEN_SCRATCH_VALUE,
        );

        // The artifact summary's per-allele read totals. **Zero rather than a sentinel** —
        // see the field's own comment: its only reader zeroes it in the same call.
        resize_and_fill(&mut self.pooled_allele_reads, allele_count, 0);
        let spread_len = genotype_count.checked_mul(allele_count).unwrap_or_else(|| {
            panic!(
                "a locus of {genotype_count} genotypes over {allele_count} alleles needs an \
                 error-spread table longer than a usize can index"
            )
        });
        resize_and_fill(&mut self.error_spreads, spread_len, UNWRITTEN_SCRATCH_VALUE);
        // **The two per-row maps are cleared, not resized.** Their length is what says how
        // many rows the caller has claimed so far, and it claims them one at a time after this
        // call — so a resize here would be a length nothing yet means.
        self.run_sample_of_each_row.clear();
        self.inbreeding_coefficient_by_row.clear();
        self.emission_cost = EmissionCost::default();

        // **One emission cache per row**, grown and never shrunk. Each sizes itself when it is
        // filled, so nothing here says how wide a sample is.
        if self.generic_rows.len() < row_count {
            self.generic_rows
                .resize_with(row_count, GenericRowScratch::default);
        }

        // **The contaminant tables are un-sized here and re-sized by the driver**, and only on a
        // run that fitted a fraction. Clearing rather than leaving them is what stops a locus
        // called without contamination from reading the last contaminated locus's frequencies:
        // every accessor refuses a zero batch count, and `resize` alone would have left a table
        // of the right length holding another locus's numbers.
        self.batch_count = 0;
        self.run_sample_count = 0;
        self.batch_allele_copies.clear();
        self.contaminant_allele_frequencies.clear();
        self.expected_copies_by_run_sample.clear();

        self.row_count = row_count;
        self.genotype_count = genotype_count;
        self.allele_count = allele_count;
    }

    /// How many **rows** the buffers are currently sized for — one per sample this locus is
    /// called on, which at a locus with an uncallable sample is fewer than the run's samples.
    #[inline]
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.row_count
    }

    /// How many genotypes the buffers are currently sized for.
    #[inline]
    #[must_use]
    pub fn genotype_count(&self) -> usize {
        self.genotype_count
    }

    /// How many alleles the buffers are currently sized for.
    #[inline]
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// One sample's row of the genotype-likelihood table — one entry per candidate genotype.
    ///
    /// **The one spelling of the index.** The table is flat and sample-major, so a caller
    /// slicing it itself would, on a slip, read a different sample's likelihoods and score
    /// that sample's reads onto this one's genotype with nothing failing.
    ///
    /// Until this locus's pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub fn sample_genotype_likelihoods(&self, sample: usize) -> &[LogProb] {
        let range = self.genotype_row_range(sample);
        &self.genotype_likelihoods[range]
    }

    /// One sample's row of the genotype-likelihood table, to fill.
    ///
    /// # Panics
    ///
    /// As [`sample_genotype_likelihoods`](Self::sample_genotype_likelihoods).
    #[inline]
    pub fn sample_genotype_likelihoods_mut(&mut self, sample: usize) -> &mut [LogProb] {
        let range = self.genotype_row_range(sample);
        &mut self.genotype_likelihoods[range]
    }

    /// One sample's own expected allele copies — one entry per allele. The cohort's sum is
    /// [`cohort_expected_copies`](Self::cohort_expected_copies), and the prior's
    /// leave-one-out term is the difference of the two.
    ///
    /// Until this locus's pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub fn sample_expected_copies(&self, sample: usize) -> &[f64] {
        let range = self.allele_row_range(sample);
        &self.per_sample_expected_copies[range]
    }

    /// One sample's own expected allele copies, to fill.
    ///
    /// # Panics
    ///
    /// As [`sample_expected_copies`](Self::sample_expected_copies).
    #[inline]
    pub fn sample_expected_copies_mut(&mut self, sample: usize) -> &mut [f64] {
        let range = self.allele_row_range(sample);
        &mut self.per_sample_expected_copies[range]
    }

    /// The **cohort's** expected allele copies as this pass has them, summed over every
    /// sample. Until this pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn cohort_expected_copies(&self) -> &[f64] {
        self.assert_prepared();
        &self.cohort_expected_copies
    }

    /// The cohort's expected allele copies, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn cohort_expected_copies_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.cohort_expected_copies
    }

    /// What the previous pass's cohort expected copies were — the convergence comparison.
    /// Before the first pass has advanced, every entry is [`UNWRITTEN_SCRATCH_VALUE`], and
    /// every comparison against it is therefore false.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn previous_cohort_expected_copies(&self) -> &[f64] {
        self.assert_prepared();
        &self.previous_cohort_expected_copies
    }

    /// Make this pass's cohort copies the previous pass's, and hand back a buffer to write
    /// the next pass into.
    ///
    /// **A swap rather than a copy**, so there is one spelling of *which buffer is now
    /// which*. The returned buffer arrives holding [`UNWRITTEN_SCRATCH_VALUE`] in every
    /// entry: a pass that leaves an allele unwritten fails the next check it reaches,
    /// rather than reading the pass-before-last's value and letting the convergence test
    /// compare a number no pass wrote.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn advance_cohort_expected_copies(&mut self) -> &mut [f64] {
        self.assert_prepared();
        std::mem::swap(
            &mut self.cohort_expected_copies,
            &mut self.previous_cohort_expected_copies,
        );
        self.cohort_expected_copies.fill(UNWRITTEN_SCRATCH_VALUE);
        &mut self.cohort_expected_copies
    }

    /// The locus's seed concentration. Until this locus fills it, every entry is
    /// [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn seed_concentration(&self) -> &[f64] {
        self.assert_prepared();
        &self.seed_concentration
    }

    /// The locus's seed concentration, to fill once per locus.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn seed_concentration_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.seed_concentration
    }

    /// The sample being scored, its own concentration.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn sample_concentration(&self) -> &[f64] {
        self.assert_prepared();
        &self.sample_concentration
    }

    /// The sample being scored, its own concentration, to fill per sample per pass.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn sample_concentration_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.sample_concentration
    }

    /// The sample being scored, its log-prior over every candidate genotype.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn prior_row(&self) -> &[LogProb] {
        self.assert_prepared();
        &self.prior_row
    }

    /// The sample's log-prior row, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn prior_row_mut(&mut self) -> &mut [LogProb] {
        self.assert_prepared();
        &mut self.prior_row
    }

    /// The sample being scored, its posterior over every candidate genotype.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn posterior_row(&self) -> &[f64] {
        self.assert_prepared();
        &self.posterior_row
    }

    /// The sample's posterior row, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn posterior_row_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.posterior_row
    }

    /// The genotype prior's per-allele working space.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn prior_per_allele_workspace(&self) -> &[f64] {
        self.assert_prepared();
        &self.prior_per_allele_workspace
    }

    /// The genotype prior's per-allele working space, to lend it.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn prior_per_allele_workspace_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.prior_per_allele_workspace
    }

    /// Candidate selection's buffers. **Not sized by
    /// [`prepare_for_locus`](Self::prepare_for_locus)** — it owns its own, and selection
    /// runs before the locus's shape is known.
    #[inline]
    pub fn candidate_selection_mut(&mut self) -> &mut allele_candidates::SelectionScratch {
        &mut self.candidate_selection
    }

    /// The three buffers the SNP/indel row builder needs, for one scratch row.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, and on a row past the ones it was prepared for.
    #[must_use]
    pub(crate) fn generic_row_buffers_mut(&mut self, row: usize) -> GenericRowBuffers<'_> {
        let range = self.genotype_row_range(row);
        GenericRowBuffers {
            error_spreads: &self.error_spreads,
            contaminant_allele_frequencies: &self.contaminant_allele_frequencies,
            row_scratch: &mut self.generic_rows[row],
            genotype_likelihoods: &mut self.genotype_likelihoods[range],
        }
    }

    /// The SNP/indel row's per-allele error-spread table, to fill once per locus.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub(crate) fn error_spreads_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.error_spreads
    }

    /// Claim a scratch row for one run sample, with its inbreeding coefficient.
    ///
    /// Called once per sample the locus is **called on**, in the run's sample order, and
    /// **after** [`prepare_for_locus`](Self::prepare_for_locus), which clears the previous
    /// locus's map. Claiming first throws the claims away and leaves the locus with none,
    /// which the map's first reader then refuses by name — measured, "0 of them were
    /// claimed", at the first locus.
    pub(crate) fn claim_row_for(&mut self, run_sample: usize, inbreeding: InbreedingF) {
        self.run_sample_of_each_row.push(run_sample);
        self.inbreeding_coefficient_by_row.push(inbreeding);
    }

    /// Which run sample each scratch row holds, in row order.
    ///
    /// # Panics
    ///
    /// If the claimed rows are not one per prepared row. **The check is here rather than at
    /// the claim**, because a caller that prepared for one shape and claimed for another has
    /// a map whose every entry is a different sample's — and the symptom would be one
    /// sample's reads scored against another's likelihoods rather than an index out of range.
    #[inline]
    #[must_use]
    pub(crate) fn run_sample_of_each_row(&self) -> &[usize] {
        self.assert_rows_claimed();
        &self.run_sample_of_each_row
    }

    /// The claimed rows' inbreeding coefficients, one per row — what the loop walks.
    ///
    /// # Panics
    ///
    /// As [`run_sample_of_each_row`](Self::run_sample_of_each_row).
    #[inline]
    #[must_use]
    pub(crate) fn inbreeding_coefficient_by_row(&self) -> &[InbreedingF] {
        self.assert_rows_claimed();
        &self.inbreeding_coefficient_by_row
    }

    /// Refuse a scratch whose per-row map is not one entry per prepared row.
    #[inline]
    fn assert_rows_claimed(&self) {
        self.assert_prepared();
        assert_eq!(
            self.run_sample_of_each_row.len(),
            self.row_count,
            "this scratch was prepared for {} rows and {} of them were claimed, so the map \
             from rows back to the run's samples describes a different locus",
            self.row_count,
            self.run_sample_of_each_row.len()
        );
    }

    /// What this locus's table build has cost so far — the instrument spec §13's test 5 reads.
    ///
    /// **Nothing in a run reads it**, and that is what it is for: the cost it records is
    /// invisible in the output, because emissions recomputed every pass give the same
    /// genotypes.
    /// Step D2 of the calling loop's plan is what asserts on it.
    #[inline]
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "an instrument: what it measures — a table built once rather than once \
                      per pass — is invisible in a run's output, so its readers are step D2's \
                      tests and the bakeoffs plan's benches"
        )
    )]
    pub(crate) fn emission_cost(&self) -> EmissionCost {
        self.emission_cost
    }

    /// Record one sample's emission fill over `observations` observations and `candidates`
    /// candidates.
    pub(crate) fn charge_emission_row_fill(&mut self, observations: usize, candidates: usize) {
        self.emission_cost.emission_row_fills += 1;
        self.emission_cost.emission_evaluations += (observations * candidates) as u64;
    }

    /// Record that the locus's frequency-free emissions were computed once more.
    pub(crate) fn charge_emission_build(&mut self) {
        self.emission_cost.emission_builds += 1;
    }

    /// Record that a fresh assembly of the whole `rows × genotypes` table has begun.
    ///
    /// **The rows are charged one at a time by [`charge_row_assembly`](Self::charge_row_assembly)
    /// rather than added up here**, so that the two counts can disagree. Charged from the row
    /// count, `row_assemblies` is `table_assemblies × rows` by construction and can only repeat
    /// what the other says — measured: an assembly that stops one row short of the table then
    /// leaves both counts untouched, and only the genotypes it did not write say so.
    pub(crate) fn charge_table_assembly(&mut self) {
        self.emission_cost.table_assemblies += 1;
    }

    /// Record one row of the current assembly.
    pub(crate) fn charge_row_assembly(&mut self) {
        self.emission_cost.row_assemblies += 1;
    }

    /// **Where the twenty per-locus buffers' bytes are and how many of them there are** —
    /// the cheap half of the loop's zero-allocation invariant
    /// (`doc/devel/ng/spec/calling_em_loop.md` §13's test 7).
    ///
    /// A `Vec` that **grew** during the loop moves its bytes, so its pointer changes; one
    /// refilled in place does not.
    ///
    /// **⚠ Two of them exchange pointers with every pass, and that is by design**:
    /// [`advance_cohort_expected_copies`](Self::advance_cohort_expected_copies) `mem::swap`s the
    /// cohort's expected copies with the previous pass's rather than copying either. So this
    /// list's *order* depends on the parity of the pass count, and a caller comparing two runs
    /// must compare it as a set — sorted — or compare only runs whose pass counts have the same
    /// parity. Measured: a 2-pass run and a 7-pass run over identical evidence differ in exactly
    /// those two entries, both still `(pointer, length)` pairs the other list holds.
    ///
    /// **What this cannot see is a buffer that was reallocated and came back to the same
    /// address**, which is what a freed block of the same size usually gets — so a `Vec` cloned
    /// once a row leaves this list unchanged after an even number of clones. It cannot see a
    /// temporary allocated and dropped inside a pass either, which leaves no trace anywhere.
    /// **Both halves are counted for real**, by `tests/ng_calling_loop_allocation.rs`, which
    /// installs `dhat`'s counting allocator and reads `total_blocks` across runs at different
    /// pass counts. Measured: reallocating the contaminant frequency table once a row once a
    /// pass leaves every fingerprint here identical and takes that counter from 8 blocks to 24.
    /// It lives in a test binary of its own because a global allocator counts the whole process,
    /// and the lib suite runs its tests in parallel.
    ///
    /// **The row scratches are deliberately not here**, and their absence is the invariant
    /// rather than a gap: each `GenericRowScratch` sizes itself per *sample*, inside the
    /// emission build, so it legitimately grows within a locus when a wider sample arrives. The
    /// emission build happens once, outside the frequency loop, so a pass still allocates
    /// nothing — which is what §13's test 7 claims. Fingerprinting them would make this test
    /// fail on correct code. **The three contaminant buffers above are here**, because they are
    /// refilled inside the loop and sized outside it, which is exactly the property this
    /// measures.
    ///
    /// Test-only, because the pointers are an implementation detail that no run should read.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn buffer_fingerprints(&self) -> Vec<(usize, usize)> {
        fn of<T>(buffer: &[T]) -> (usize, usize) {
            (buffer.as_ptr() as usize, buffer.len())
        }
        vec![
            of(&self.genotype_likelihoods),
            of(&self.prior_row),
            of(&self.posterior_row),
            of(&self.seed_concentration),
            of(&self.sample_concentration),
            of(&self.prior_per_allele_workspace),
            of(&self.cohort_expected_copies),
            of(&self.previous_cohort_expected_copies),
            of(&self.per_sample_expected_copies),
            of(&self.copy_count_log_likelihoods),
            of(&self.allele_count_distribution),
            of(&self.allele_count_distribution_next),
            of(&self.log_allele_count_distribution),
            of(&self.pooled_allele_reads),
            of(&self.error_spreads),
            of(&self.run_sample_of_each_row),
            of(&self.inbreeding_coefficient_by_row),
            of(&self.batch_allele_copies),
            of(&self.contaminant_allele_frequencies),
            of(&self.expected_copies_by_run_sample),
        ]
    }

    /// **Size this locus's two contaminant tables** — one row per sequencing batch, and the
    /// run's whole sample axis for the copies those rows are summed from.
    ///
    /// **Called only where the locus's own mixture reads a per-batch frequency**, which is an
    /// **ordinary site** in a run that fitted a fraction — a repeat tract's third term is the
    /// fit's frozen length spectrum and reads none. And called after
    /// [`prepare_for_locus`](Self::prepare_for_locus), which un-sizes them. That order is what
    /// makes a locus unable to read another's frequencies: the tables come back empty and every
    /// accessor below refuses them.
    ///
    /// **The run's sample axis and not the locus's rows.** The batching is declared over the
    /// run, so a batch whose samples were all ruled uncallable at this locus still has a row in
    /// the table; summing over rows would leave it unwritten, which
    /// [`fill_batch_allele_copies`] refuses by name. Scattering onto the run's axis gives it a
    /// row of zeros — the same thing the M-step already does with an uncallable sample.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch; on a run with no batches or no samples, which is a run whose
    /// batching went missing rather than a run without one; and on a table longer than a
    /// `usize` can index.
    pub(crate) fn prepare_contaminant_tables(&mut self, batch_count: usize, run_samples: usize) {
        self.assert_prepared();
        assert!(
            batch_count > 0,
            "a contaminated run has at least one sequencing batch, the default being one that \
             holds all of it, so a locus prepared for none is a run whose batching went missing"
        );
        // **The precondition the scatter actually has**, which is not `run_samples >=
        // row_count`: the copies are written at `run_sample × alleles`, so what must hold is
        // that every claimed row names a sample the run has. A one-row locus whose row names
        // sample 5 of a two-sample run passes the count comparison and panics later on a slice
        // range that names neither the row nor the sample.
        self.assert_rows_claimed();
        let highest_named = self
            .run_sample_of_each_row
            .iter()
            .copied()
            .max()
            .expect("a prepared locus has claimed at least one row");
        assert!(
            highest_named < run_samples,
            "this locus's rows name sample {highest_named} and the run has {run_samples}; a row \
             is one of the run's samples, so a row naming one the run does not have is not this \
             run's locus"
        );
        let batch_table_len = batch_count
            .checked_mul(self.allele_count)
            .unwrap_or_else(|| {
                panic!(
                    "a run of {batch_count} sequencing batches over {} alleles needs a contaminant \
                 table longer than a usize can index",
                    self.allele_count
                )
            });
        let run_table_len = run_samples
            .checked_mul(self.allele_count)
            .unwrap_or_else(|| {
                panic!(
                    "a run of {run_samples} samples over {} alleles needs a copies table longer \
                 than a usize can index",
                    self.allele_count
                )
            });
        // **All three take the unwritten sentinel here**, because none of them is read before
        // something writes it: `batch_copy_buffers_mut` zeroes the run-keyed copies in the same
        // call that scatters into them, and the two batch tables are written whole. The zero
        // that matters is that one, and its reason lives beside it.
        resize_and_fill(
            &mut self.expected_copies_by_run_sample,
            run_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.batch_allele_copies,
            batch_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.contaminant_allele_frequencies,
            batch_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        self.batch_count = batch_count;
        self.run_sample_count = run_samples;
    }

    /// **Scatter this pass's per-row expected copies back onto the run's sample axis**, then
    /// hand back that axis and the batch table to sum it into.
    ///
    /// The scatter is here rather than at the call site because it is the join the two axes
    /// meet at: the row map is this scratch's and so is the run-keyed buffer, and a caller
    /// writing the loop itself would be writing this scratch's row arithmetic a second time.
    ///
    /// # Panics
    ///
    /// On a scratch whose contaminant tables were not prepared.
    pub(crate) fn batch_copy_buffers_mut(&mut self) -> BatchCopyBuffers<'_> {
        self.assert_contaminant_tables_prepared();
        self.assert_rows_claimed();
        // **Zeroed here rather than at sizing, and the zero is an answer rather than a slot
        // nobody filled.** Every sample this locus has no row for keeps it, and zero copies is
        // exactly what an uncallable sample contributes to the cohort's expected copies too — so
        // a batch made entirely of samples this locus set aside gets a row of zeros and takes
        // the frequency's no-evidence fallback, rather than a row nobody wrote.
        self.expected_copies_by_run_sample.fill(0.0);
        let allele_count = self.allele_count;
        for (row, &run_sample) in self.run_sample_of_each_row.iter().enumerate() {
            let from = row * allele_count..(row + 1) * allele_count;
            let onto = run_sample * allele_count..(run_sample + 1) * allele_count;
            self.expected_copies_by_run_sample[onto]
                .copy_from_slice(&self.per_sample_expected_copies[from]);
        }
        BatchCopyBuffers {
            expected_copies_by_run_sample: &self.expected_copies_by_run_sample,
            allele_count,
            batch_allele_copies: &mut self.batch_allele_copies,
        }
    }

    /// The three buffers one **sample's** contaminant frequencies are built between.
    ///
    /// # Panics
    ///
    /// On a scratch whose contaminant tables were not prepared, and on a row past what this
    /// locus was prepared for.
    pub(crate) fn contaminant_frequency_buffers_mut(
        &mut self,
        row: usize,
    ) -> ContaminantFrequencyBuffers<'_> {
        self.assert_contaminant_tables_prepared();
        let range = self.allele_row_range(row);
        ContaminantFrequencyBuffers {
            batch_allele_copies: &self.batch_allele_copies,
            own_expected_copies: &self.per_sample_expected_copies[range],
            allele_count: self.allele_count,
            contaminant_allele_frequencies: &mut self.contaminant_allele_frequencies,
        }
    }

    /// How many sequencing batches this locus's contaminant tables hold — **zero where they
    /// were never prepared**, which is every locus of an uncontaminated run.
    ///
    /// Test-only: a run reads the tables, never their shape.
    #[cfg(test)]
    #[inline]
    #[must_use]
    pub(crate) fn contaminant_batch_count(&self) -> usize {
        self.batch_count
    }

    /// The contaminant allele frequencies as the last row's fill left them, batch-major.
    ///
    /// Test-only: in a run they are read through
    /// [`generic_row_buffers_mut`](Self::generic_row_buffers_mut), one row at a time and while
    /// that row is being assembled, so nothing outside a test ever sees them afterwards.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn contaminant_allele_frequencies(&self) -> &[f64] {
        &self.contaminant_allele_frequencies
    }

    /// Each sequencing batch's expected allele copies at this locus, batch-major, as the last
    /// fill left them. Test-only, for [`contaminant_allele_frequencies`](Self::contaminant_allele_frequencies)'s
    /// reason.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn batch_allele_copies(&self) -> &[f64] {
        &self.batch_allele_copies
    }

    /// Refuse a scratch whose contaminant tables were never sized for this locus.
    ///
    /// **It says which table went missing, which is the whole of what it adds.** The fills do
    /// refuse an empty table on their own — `fill_batch_allele_copies` compares the copies
    /// against `samples × alleles` before anything else — but they refuse it as a shape
    /// mismatch between two slices, at a locus, with nothing in the message about the one call
    /// that was skipped. This names it.
    #[inline]
    fn assert_contaminant_tables_prepared(&self) {
        self.assert_prepared();
        assert!(
            self.batch_count > 0,
            "this locus's contaminant tables were not prepared: `prepare_contaminant_tables` \
             runs after `prepare_for_locus` on a run that fitted a contamination fraction, and \
             nothing else sizes them"
        );
    }

    /// One scratch row's SNP/indel emission cache, which owns its own sizing.
    ///
    /// # Panics
    ///
    /// On a row past what [`prepare_for_locus`](Self::prepare_for_locus) sized this scratch
    /// for.
    #[inline]
    pub fn generic_row_mut(&mut self, row: usize) -> &mut GenericRowScratch {
        assert!(
            row < self.row_count,
            "row {row} is past the {} this scratch was prepared for",
            self.row_count
        );
        &mut self.generic_rows[row]
    }

    /// The repeat-tract row's own scratch, which owns its own sizing.
    #[inline]
    pub fn ssr_row_mut(&mut self) -> &mut SsrRowScratch<SsrEmissionScratch> {
        &mut self.ssr_row
    }

    /// **The tract's fitted scoring parameters, to gather into once per repeat tract.**
    ///
    /// Cleared and refilled by
    /// [`TractScoringFits::gather_for_locus`](super::calling::inference::repeat_tract_parameters::TractScoringFits::gather_for_locus),
    /// which is the only thing that writes it.
    #[inline]
    pub(crate) fn tract_fits_mut(
        &mut self,
    ) -> &mut crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits {
        &mut self.tract_fits
    }

    /// The same fits to read — **what the locus's warrant is folded from**, once the rows have
    /// been scored.
    ///
    /// **Separate from [`tract_locus_buffers_mut`](Self::tract_locus_buffers_mut) because the
    /// warrant is read after the scoring rather than during it**, and by then nothing else of
    /// this scratch is borrowed.
    #[inline]
    pub(crate) fn tract_fits(
        &self,
    ) -> &crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits {
        &self.tract_fits
    }

    /// **Everything a repeat tract's rows are scored from, borrowed from one scratch at once.**
    ///
    /// **It exists because the contexts the row takes borrow the fits, and the row writes into
    /// this same scratch.** The contexts are built once per tract, from `fits`, and every row
    /// then needs them alive while `ssr_row` and the likelihood table are written — three
    /// borrows of three different fields, which no sequence of one-field accessors can hold at
    /// once. **The borrows are disjoint by construction**: every field below is a different
    /// field of this scratch.
    ///
    /// **The whole table rather than one row**, unlike the SNP/indel path's, because the tract's
    /// row scratch is one buffer reused across rows rather than one per row — so the loop over
    /// rows lives inside this borrow instead of outside it.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, and on a row map that is not one entry per prepared row —
    /// both of which say the locus's shape was never fixed. **The row map is checked here
    /// rather than where it is read**, because the walk over rows reads it inside this borrow,
    /// where no accessor of this scratch can be called.
    #[inline]
    pub(crate) fn tract_locus_buffers_mut(&mut self) -> TractLocusBuffers<'_, SsrEmissionScratch> {
        self.assert_rows_claimed();
        TractLocusBuffers {
            fits: &self.tract_fits,
            run_sample_of_each_row: &self.run_sample_of_each_row,
            row_scratch: &mut self.ssr_row,
            genotype_likelihoods: &mut self.genotype_likelihoods,
            genotype_count: self.genotype_count,
        }
    }

    /// Every buffer scoring **one sample** at this locus reads or writes, handed out
    /// together.
    ///
    /// **It exists because four of them have to be live at once and they are fields of one
    /// scratch.** Scoring a sample builds its concentration from three buffers into a
    /// fourth, its log-prior from that concentration and a working buffer into a fifth, its
    /// posterior from that prior and its likelihood row into a sixth, and finally its own
    /// expected copies from the posterior. Reaching for those one accessor at a time does
    /// not compile — each borrows the whole scratch — and the way round it that needs no
    /// new type is to copy a buffer, which would put an allocation back into a pass whose
    /// whole shape exists to have none.
    ///
    /// **The borrows are disjoint by construction**: every borrowed field below is a
    /// different field of this scratch, so nothing here can alias.
    ///
    /// **`_mut` because it is the mutable half of five accessor pairs at once**, not the
    /// shared half of a sixth — five of the eight borrows it hands out are `&mut`, and the
    /// suffix is the only thing on the name that says a caller needs a mutable scratch.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub(crate) fn sample_scoring_buffers_mut(&mut self, sample: usize) -> SampleScoringBuffers<'_> {
        let genotype_row = self.genotype_row_range(sample);
        let allele_row = self.allele_row_range(sample);
        SampleScoringBuffers {
            sample,
            seed_concentration: &self.seed_concentration,
            cohort_expected_copies: &self.cohort_expected_copies,
            genotype_likelihoods: &self.genotype_likelihoods[genotype_row],
            sample_concentration: &mut self.sample_concentration,
            prior_per_allele_workspace: &mut self.prior_per_allele_workspace,
            prior_row: &mut self.prior_row,
            posterior_row: &mut self.posterior_row,
            sample_expected_copies: &mut self.per_sample_expected_copies[allele_row],
        }
    }

    /// The two buffers the **M-step** reads and writes: every sample's own expected allele
    /// copies, and the cohort row their sum goes into.
    ///
    /// **A second bundle for the same reason as the first** — the per-sample copies and the
    /// cohort's are two fields of one scratch, and reaching for them one accessor at a time
    /// does not compile. Measured while B1 was in review: an M-step written against the
    /// per-buffer accessors gives two `error[E0502]`.
    ///
    /// **The whole per-sample table, not one row.** The M-step's defining property is that it
    /// adds the samples up *in the run's fixed order*, so it takes the table and walks it
    /// itself rather than being handed a row at a time by a caller whose loop could reorder.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub(crate) fn cohort_summing_buffers_mut(&mut self) -> CohortSummingBuffers<'_> {
        self.assert_prepared();
        CohortSummingBuffers {
            sample_count: self.row_count,
            per_sample_expected_copies: &self.per_sample_expected_copies,
            cohort_expected_copies: &mut self.cohort_expected_copies,
        }
    }

    /// The artifact summary's per-allele read totals, to zero and refill.
    ///
    /// **Handed out mutably and never read back through a shared accessor**, because the
    /// one caller that fills it is the one caller that reads it, in the same call: the
    /// final pass pools every called sample's reads onto the alleles, picks the primary
    /// alternative from the totals, and is done with them
    /// (`doc/devel/ng/spec/calling_quality.md` §6.3). There is nothing for a later stage to
    /// come back for.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub(crate) fn pooled_allele_reads_mut(&mut self) -> &mut [u64] {
        self.assert_prepared();
        &mut self.pooled_allele_reads
    }

    /// The five buffers the **site quality** reads and writes: the likelihood table the loop
    /// built, and the four the fold needs.
    ///
    /// **A third bundle, for the third reason of the same kind** — five fields of one
    /// scratch, and the fold swaps two of them, which needs both to be borrowed at once from
    /// the same call.
    ///
    /// **Called once per locus, after the loop has stopped.** The likelihood table it hands
    /// over is the one the loop leaves behind, which is what makes the site quality cost no
    /// emission evaluations of its own (`doc/devel/ng/spec/calling_quality.md` §3.2) — and
    /// where a fraction was fitted, the one assembled against the frequencies the loop settled
    /// on, so the quality and the genotypes beside it come from the same table.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn site_quality_buffers_mut(&mut self) -> SiteQualityBuffers<'_> {
        self.assert_prepared();
        SiteQualityBuffers {
            sample_count: self.row_count,
            genotype_likelihoods: &self.genotype_likelihoods,
            copy_count_log_likelihoods: &mut self.copy_count_log_likelihoods,
            allele_count_distribution: &mut self.allele_count_distribution,
            allele_count_distribution_next: &mut self.allele_count_distribution_next,
            log_allele_count_distribution: &mut self.log_allele_count_distribution,
        }
    }

    /// Refuse a scratch that was never sized for a locus.
    ///
    /// # Panics
    ///
    /// If [`prepare_for_locus`](Self::prepare_for_locus) has not been called. Every buffer
    /// of a freshly-made scratch is empty, so without this a pass folding into
    /// [`cohort_expected_copies_mut`](Self::cohort_expected_copies_mut) would run zero
    /// iterations, write nothing, and leave the cohort's copies summing to a plausible
    /// `0.0` — which the fill value cannot catch, because an unprepared buffer has no slots
    /// to fill.
    #[inline]
    fn assert_prepared(&self) {
        assert!(
            self.row_count > 0,
            "this scratch has not been prepared for a locus: call prepare_for_locus first"
        );
    }

    /// Where one sample's row of the genotype-likelihood table starts and ends.
    #[inline]
    fn genotype_row_range(&self, sample: usize) -> std::ops::Range<usize> {
        self.row_range(sample, self.genotype_count, "genotype-likelihood")
    }

    /// Where one sample's row of the per-sample expected copies starts and ends.
    #[inline]
    fn allele_row_range(&self, sample: usize) -> std::ops::Range<usize> {
        self.row_range(sample, self.allele_count, "per-sample expected copies")
    }

    /// Where one sample's row of a flat sample-major table starts and ends.
    ///
    /// **Private, and reached only through the two typed wrappers**, so no caller picks the
    /// row width by hand — the two widths are the genotype count and the allele count, and
    /// at a diploid biallelic locus they are 3 and 2, so a wrong one is a legal range over
    /// the wrong table.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given — a caller that walked one
    /// sample too far would otherwise read the next sample's row, or, at the last sample,
    /// panic somewhere that says nothing about which sample it was.
    #[inline]
    fn row_range(&self, sample: usize, width: usize, table: &str) -> std::ops::Range<usize> {
        self.assert_prepared();
        assert!(
            sample < self.row_count,
            "sample {sample} is past the {} this scratch was prepared for, indexing the \
             {table} table",
            self.row_count
        );
        let start = sample * width;
        start..start + width
    }
}

/// Every buffer one sample's turn through a pass touches, borrowed from one
/// [`CallingScratch`] at once.
///
/// **Every run-time construction goes through
/// [`CallingScratch::sample_scoring_buffers_mut`]**, which is where the flat sample-major
/// tables are sliced, so the row arithmetic has one spelling. See that method for why the
/// buffers travel together rather than one at a time.
///
/// **The fields are public because the tests build one by hand**, deliberately: the only way
/// to hand the scorer a row of the wrong width — the mis-shape its release-held checks exist
/// for — is to reach past the accessor, which cannot produce one. Nothing in a run should.
/// The shapes are not this type's to enforce: they are checked where they are used, in
/// `score_one_sample`, `fill_sample_concentration` and `PriorRow::new`.
///
/// Of the eight buffers, three are read-only — what the locus's setup and the *previous* pass
/// left — and four are this sample's share of what this pass produces. The eighth is both:
/// this sample's expected copies are read as the leave-one-out term the prior needs, and then
/// overwritten with what this pass's posterior implies. The `sample` index beside them is not
/// a buffer; it is there so a panic can say which sample was being scored.
#[derive(Debug)]
pub(crate) struct SampleScoringBuffers<'a> {
    /// Which sample of the run these buffers belong to — carried so a panic raised while
    /// scoring can name it. At a thousand samples that is the first thing the reader of a
    /// message wants, and no other field records it.
    pub sample: usize,
    /// The locus's seed concentration — one entry per allele, the same for every sample,
    /// and what the prior falls back to where the cohort says nothing.
    pub seed_concentration: &'a [f64],
    /// The cohort's expected allele copies as the previous pass left them, this sample's
    /// own contribution included.
    pub cohort_expected_copies: &'a [f64],
    /// This sample's row of the genotype-likelihood table — one entry per candidate
    /// genotype, assembled from emissions the locus computed once.
    pub genotype_likelihoods: &'a [LogProb],
    /// One entry per allele, to fill: this sample's own concentration — the seed plus what
    /// the **other** samples showed here.
    pub sample_concentration: &'a mut [f64],
    /// One entry per allele of working space, to hand the genotype prior. Its contents on
    /// entry are ignored and on exit unspecified.
    pub prior_per_allele_workspace: &'a mut [f64],
    /// One entry per candidate genotype, to fill: this sample's log-prior.
    pub prior_row: &'a mut [LogProb],
    /// One entry per candidate genotype, to fill: this sample's posterior over genotypes.
    pub posterior_row: &'a mut [f64],
    /// One entry per allele: this sample's own expected allele copies — **read first** as
    /// the leave-one-out term, then **overwritten** with this pass's answer.
    pub sample_expected_copies: &'a mut [f64],
}

/// The three buffers the **SNP/indel row builder** reads and writes, borrowed from one
/// [`CallingScratch`] in one call.
///
/// **A bundle for the reason the other three are**: they are three fields of one scratch, and
/// reaching for them an accessor at a time does not compile. The error spreads are shared by
/// every sample and filled once per locus; the other two are this row's.
#[derive(Debug)]
pub(crate) struct GenericRowBuffers<'a> {
    /// The locus's per-allele error-spread table, filled once before the first row.
    pub(crate) error_spreads: &'a [f64],
    /// The contaminant allele frequencies **this row** is scored against, batch-major — empty
    /// on a run the fit found no contamination in, where the mixture is absent and the row
    /// computes the plain formula.
    pub(crate) contaminant_allele_frequencies: &'a [f64],
    /// This row's own emission cache, which it sizes itself.
    pub(crate) row_scratch: &'a mut GenericRowScratch,
    /// This sample's row of the genotype-likelihood table, to fill.
    pub(crate) genotype_likelihoods: &'a mut [LogProb],
}

/// **Everything one repeat tract's rows are scored from**, borrowed from one
/// [`CallingScratch`] at once — [`CallingScratch::tract_locus_buffers_mut`] is the only way to
/// make one.
///
/// **It hands out the whole likelihood table rather than one row**, which is the difference from
/// [`GenericRowBuffers`]: a tract's row scratch is one reused buffer rather than one per row, so
/// the walk over rows happens inside this borrow.
#[derive(Debug)]
pub(crate) struct TractLocusBuffers<'a, SsrEmissionScratch> {
    /// This tract's fitted scoring parameters — what the contexts are built from, and what the
    /// locus's warrant is folded from.
    pub(crate) fits: &'a crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits,
    /// Which sample of the run each scratch row holds — **the map the walk over rows needs and
    /// cannot ask for separately**, since asking would be a second borrow of this scratch while
    /// the three below are held.
    ///
    /// **At a repeat tract this is provably the run's sample order unchanged**, because a tract
    /// rules no sample out — so no test of the tract path exercises the difference between a
    /// row and the sample it holds. It is carried rather than assumed so that the day a tract
    /// can set a sample aside, the walk that reads it needs no change.
    pub(crate) run_sample_of_each_row: &'a [usize],
    /// The repeat-tract row's working memory, reused across rows.
    pub(crate) row_scratch: &'a mut SsrRowScratch<SsrEmissionScratch>,
    /// `rows × genotypes`, row-major: the whole genotype-likelihood table, to fill.
    pub(crate) genotype_likelihoods: &'a mut [LogProb],
    /// How many candidate genotypes this locus has — the table's stride.
    pub(crate) genotype_count: usize,
}

/// The buffers each **sequencing batch's** expected allele copies are summed between, borrowed
/// from one [`CallingScratch`] at once — [`CallingScratch::batch_copy_buffers_mut`] is the only
/// way to make one, because it is what performs the scatter the first field's name implies.
#[derive(Debug)]
pub(crate) struct BatchCopyBuffers<'a> {
    /// This pass's expected allele copies on the **run's** sample axis, sample-major, zero at
    /// every sample this locus has no row for.
    pub(crate) expected_copies_by_run_sample: &'a [f64],
    /// How many alleles this locus is called over — the stride of both tables.
    pub(crate) allele_count: usize,
    /// One row per sequencing batch, to fill.
    pub(crate) batch_allele_copies: &'a mut [f64],
}

/// The buffers **one sample's** contaminant allele frequencies are built between, borrowed from
/// one [`CallingScratch`] at once.
#[derive(Debug)]
pub(crate) struct ContaminantFrequencyBuffers<'a> {
    /// Every batch's expected allele copies at this locus, batch-major — this sample's own
    /// among them, which is what the fill takes back out.
    pub(crate) batch_allele_copies: &'a [f64],
    /// This sample's own expected allele copies.
    pub(crate) own_expected_copies: &'a [f64],
    /// How many alleles this locus is called over — the stride of both tables.
    pub(crate) allele_count: usize,
    /// One row per sequencing batch, to fill: the frequencies **this** sample is scored
    /// against.
    pub(crate) contaminant_allele_frequencies: &'a mut [f64],
}

/// The two buffers the **M-step** works between, borrowed from one [`CallingScratch`] at once.
///
/// **The M-step** is the half of one pass that turns every sample's genotype probabilities back
/// into the cohort-wide summary the next pass conditions on: it adds up how many copies of each
/// allele the samples are expected to carry. Its counterpart is the E-step, which goes the other
/// way — [`SampleScoringBuffers`] is that half's bundle
/// (`doc/devel/ng/spec/calling_em_loop.md` §2).
///
/// **Made only by [`CallingScratch::cohort_summing_buffers_mut`]** in a run; the tests build one
/// by hand to hand the sum a mis-shaped table, which the accessor cannot produce.
#[derive(Debug)]
pub(crate) struct CohortSummingBuffers<'a> {
    /// How many samples the table below holds — the run's whole sample count, which is what
    /// the scratch was prepared for.
    ///
    /// **Deliberately redundant against the table's own length, because that redundancy is
    /// what makes the shape check a check.** Derive the count from the length instead and any
    /// table divides evenly by the allele count, so a two-row table presented as three samples
    /// is accepted and the cohort comes back short by a sample with nothing raised. Measured
    /// during B2's review: with the count derived, the suite stays green and the
    /// wrong-size test starts accepting that table and returning `[2.0, 2.0]`.
    pub sample_count: usize,
    /// `samples × alleles`, sample-major: each sample's own expected allele copies, as the
    /// E-step of this pass left them. Read whole, so that the order of the sum is this
    /// function's and not its caller's.
    pub per_sample_expected_copies: &'a [f64],
    /// One entry per allele, to fill: the cohort's expected allele copies, the sum of the
    /// rows above.
    pub cohort_expected_copies: &'a mut [f64],
}

/// Resize a scratch buffer and overwrite every entry, including the ones that were already
/// there.
///
/// **`Vec::resize` alone is the bug this exists to avoid**: it leaves the leading entries
/// untouched, so a locus with the same shape as the last one would silently reuse the last
/// one's values — a wrong prior or a wrong likelihood at every sample, with nothing failing.
fn resize_and_fill<T: Copy>(buffer: &mut Vec<T>, len: usize, fill: T) {
    buffer.clear();
    buffer.resize(len, fill);
}

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
    },
    /// The candidate step ruled this sample uncallable at this locus, so it took no part in
    /// the loop and emission writes its `GT` as missing.
    ///
    /// **This happens on the SNP/indel path only.** A repeat tract sets no sample aside
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), and [`LocusInference::new`] refuses
    /// a locus that carries one there.
    Missing,
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
    /// fitted length spectrum, its motif period's pooled tracts, or a flat shape at a stated
    /// concentration (`doc/devel/ng/spec/population_diversity.md` §4.4).
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
    length_spectrum_rung: LengthSpectrumRung,
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
        length_spectrum_rung: LengthSpectrumRung,
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
            length_spectrum_rung,
            scoring_cells,
            cells_with_no_fitted_slippage,
            cells_whose_read_group_the_fit_does_not_describe,
            cells_with_no_fitted_substitution_rate,
            contaminant_term_was_built,
        }
    }

    /// Which rung of the tract ladder this tract's prior shape came from — see
    /// [`Self::length_spectrum_rung`].
    #[inline]
    #[must_use]
    pub fn length_spectrum_rung(&self) -> LengthSpectrumRung {
        self.length_spectrum_rung
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
        // and so no one length spectrum. Nothing can reach it: `tract_candidates` refuses a
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
    /// [`evidence_for_output`](crate::ng::run::records::evidence_for_output) takes this as the
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
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Motif, Position};

    fn generic_reference() -> CandidateAlleles {
        CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic)
    }

    /// A table cannot be built without its reference allele, and the reference stays
    /// at index 0 however many alternatives are admitted afterwards.
    ///
    /// This is the invariant every downstream branch rests on — REF against ALT in the
    /// VCF, the homozygous-reference genotype in the prior — and the reason the
    /// `alleles` field is private rather than the `Vec` the architecture sketches.
    #[test]
    fn the_reference_allele_is_id_zero_and_stays_there() {
        let mut candidates = generic_reference();
        assert_eq!(candidates.reference(), b"A");
        assert_eq!(
            candidates.bases_of(AlleleId::REFERENCE),
            Some(b"A".as_slice())
        );
        assert_eq!(candidates.len(), 1);

        let first_alternate = candidates.admit(Box::from(b"AT".as_slice()));
        let second_alternate = candidates.admit(Box::from(b"C".as_slice()));

        // The reference did not move, and the alternates got the ids that name them.
        assert_eq!(candidates.reference(), b"A");
        assert_eq!(candidates.len(), 3);
        assert_eq!(first_alternate, AlleleId(1));
        assert_eq!(second_alternate, AlleleId(2));
        assert_eq!(candidates.bases_of(first_alternate), Some(b"AT".as_slice()));
        assert_eq!(candidates.bases_of(second_alternate), Some(b"C".as_slice()));
        assert!(AlleleId::REFERENCE.is_reference());
        assert!(!first_alternate.is_reference());
    }

    /// A reference allele spelling no bases is not a locus without variation; it is a
    /// locus whose reference went missing, and it would reach the VCF's `REF` column
    /// as an unparseable record rather than a crash.
    #[test]
    #[should_panic(expected = "never nothing")]
    fn a_reference_allele_that_spells_no_bases_is_refused() {
        let _ = CandidateAlleles::new(Box::from(b"".as_slice()), LocusKind::Generic);
    }

    /// The order a scorer walks is the order `bases_of` resolves and the order the
    /// expected copies run in — checked by resolving every id the walk offers.
    #[test]
    fn the_walk_is_in_id_order_with_the_reference_first() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"AT".as_slice()));
        candidates.admit(Box::from(b"C".as_slice()));

        let spelled: Vec<&[u8]> = candidates.iter().collect();
        assert_eq!(
            spelled,
            vec![b"A".as_slice(), b"AT".as_slice(), b"C".as_slice()]
        );
        assert_eq!(candidates.iter().len(), candidates.len());

        for (index, allele) in candidates.iter().enumerate() {
            let id = AlleleId(u16::try_from(index).expect("a small fixture"));
            assert_eq!(candidates.bases_of(id), Some(allele));
        }
    }

    /// An id from another locus is a legal `u16` here and names nothing — which is
    /// what `AlleleId`'s doc comment promises is caught when the table is read.
    ///
    /// **The multi-allele table is the load-bearing half.** At one allele, an
    /// implementation that indexed instead of looking up can only panic; it takes a
    /// wider table for the worse failure — a foreign id landing on a real but wrong
    /// allele — to be possible at all, and that is the one the checked lookup exists
    /// to stop.
    #[test]
    fn an_id_this_table_does_not_hold_resolves_to_nothing() {
        let candidates = generic_reference();
        assert_eq!(candidates.bases_of(AlleleId(1)), None);
        assert_eq!(candidates.bases_of(AlleleId(u16::MAX)), None);

        let mut wider = generic_reference();
        wider.admit(Box::from(b"AT".as_slice()));
        wider.admit(Box::from(b"C".as_slice()));
        assert_eq!(wider.bases_of(AlleleId(2)), Some(b"C".as_slice())); // the last real id
        assert_eq!(wider.bases_of(AlleleId(3)), None); // one past the end
        assert_eq!(wider.bases_of(AlleleId(u16::MAX)), None);
    }

    /// A table wider than a byte, where an id above 255 must still resolve to its own
    /// allele.
    ///
    /// **This is the width no other fixture reaches, and one arithmetic slip needs
    /// it.** Narrowing the id inside the lookup — `id.get() as u8` — is invisible
    /// below 256 alleles and, above it, hands back a real but *wrong* allele and
    /// reports a present one for an id past the end. That is the failure the checked
    /// lookup exists to stop, and it is the reason `AlleleId` is a `u16` rather than a
    /// `u8`. Each allele spells its own id in bases, so a wrong answer names itself.
    #[test]
    fn an_id_above_a_bytes_worth_resolves_to_its_own_allele() {
        const WIDTH: u16 = 300;

        let mut candidates = generic_reference();
        for id in 1..WIDTH {
            candidates.admit(format!("A{id}").into_bytes().into_boxed_slice());
        }
        assert_eq!(candidates.len(), usize::from(WIDTH));

        for id in [1u16, 255, 256, 257, WIDTH - 1] {
            assert_eq!(
                candidates.bases_of(AlleleId(id)),
                Some(format!("A{id}").as_bytes()),
                "allele {id} did not resolve to itself"
            );
        }
        assert_eq!(candidates.bases_of(AlleleId(WIDTH)), None);
    }

    /// The table is one id short of holding an id for every `u16`, so the next
    /// admission has no id left to mint and says so.
    ///
    /// **Without this the guard cannot be told from its absence.** Spelling the
    /// conversion `len() as u16` instead wraps at 65,536 to `AlleleId(0)`, which is
    /// the reference — so a discovered alternative would silently *become* the
    /// reference for every downstream test, with nothing panicking anywhere.
    #[test]
    #[should_panic(expected = "more alleles than an AlleleId can name")]
    fn admitting_past_the_widest_table_an_id_can_name_is_refused() {
        let mut candidates = generic_reference();
        for _ in 1..=u16::MAX {
            candidates.admit(Box::from(b"A".as_slice()));
        }
        assert_eq!(candidates.len(), 65_536);
        candidates.admit(Box::from(b"T".as_slice()));
    }

    /// The kind travels with the table, because it is what routes the read model and
    /// the prior's seed. All three variants, and the payload-carrying one in full:
    /// for a repeat tract the flanks **are** the read model's anchor, so a dropped or
    /// truncated one is a wrong likelihood at every repeat locus with nothing
    /// panicking.
    #[test]
    fn the_locus_kind_travels_with_the_allele_table() {
        assert_eq!(generic_reference().kind(), &LocusKind::Generic);

        let bundle = CandidateAlleles::new(Box::from(b"CAG".as_slice()), LocusKind::SsrBundle);
        assert_eq!(bundle.kind(), &LocusKind::SsrBundle);
        assert_eq!(bundle.reference(), b"CAG");

        let motif = Motif::new(b"AT").expect("a dinucleotide motif");
        let tract = CandidateAlleles::new(
            Box::from(b"ATAT".as_slice()),
            LocusKind::Ssr(SsrDetail {
                motif,
                left_flank: Box::from(b"CCCGGG".as_slice()),
                right_flank: Box::from(b"TTTAAA".as_slice()),
            }),
        );
        match tract.kind() {
            LocusKind::Ssr(detail) => {
                assert_eq!(detail.motif, motif);
                assert_eq!(&*detail.left_flank, b"CCCGGG");
                assert_eq!(&*detail.right_flank, b"TTTAAA");
            }
            other => panic!("the kind came back as {other:?}"),
        }
    }

    /// Expected copies are fractional. A cohort where every sample is uncertain has no
    /// whole number anywhere in it, and a constructor that rounded or clamped would
    /// destroy exactly the signal the low-coverage case runs on.
    ///
    /// The fixture straddles 1.0 in both directions and includes a value below 0.5, so
    /// a clamp at either end of `[0, 1]` and a round to the nearest whole number all
    /// show. Checked bit for bit, not merely equal.
    #[test]
    fn expected_copies_are_fractional_and_kept_as_given() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"T".as_slice()));

        let copies = ExpectedAlleleCopies::new(vec![1.4, 0.6], &candidates);
        assert_eq!(copies.copies(), [1.4, 0.6]);
        assert_eq!(copies.copies()[0].to_bits(), 1.4f64.to_bits());
        assert_eq!(copies.copies()[1].to_bits(), 0.6f64.to_bits());
        assert_eq!(copies.len(), candidates.len());

        // The smallest legal locus: one allele, and a diploid cohort of one carrying
        // two copies of it. The immediate neighbour of the length the constructor
        // refuses.
        let monomorphic = ExpectedAlleleCopies::new(vec![2.0], &generic_reference());
        assert_eq!(monomorphic.copies(), [2.0]);
    }

    /// Copies of a length the allele table does not have are not a locus without
    /// variation; every consumer indexing by `AlleleId` would read a different
    /// allele's count, or index past the end.
    #[test]
    #[should_panic(expected = "run parallel to the locus's allele table")]
    fn expected_copies_cannot_be_built_against_a_table_of_another_width() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"T".as_slice()));
        let _ = ExpectedAlleleCopies::new(vec![2.0], &candidates);
    }

    /// A `NaN` copy count is the one that must not pass silently: it is never equal to
    /// itself, so the loop's stopping test could never be satisfied and the locus
    /// would run to its pass cap and be emitted as unconverged, naming nothing.
    #[test]
    #[should_panic(expected = "finite and at or above zero")]
    fn expected_copies_reject_a_count_that_is_not_a_number() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"T".as_slice()));
        let _ = ExpectedAlleleCopies::new(vec![1.0, f64::NAN], &candidates);
    }

    /// A negative count is not a low-coverage answer — a genome cannot carry fewer
    /// than none of an allele — it is arithmetic that went wrong upstream.
    #[test]
    #[should_panic(expected = "finite and at or above zero")]
    fn expected_copies_reject_a_negative_count() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"T".as_slice()));
        let _ = ExpectedAlleleCopies::new(vec![1.0, -0.25], &candidates);
    }

    proptest::proptest! {
        /// Over arbitrary sequences of admitted alleles: every id `admit` mints
        /// resolves back to the bases that were admitted, ids are dense from 1, the
        /// reference never moves off index 0 however many discovery rounds run, and
        /// the first id past the end resolves to nothing.
        ///
        /// The point tests above walk a fixed three-allele fixture, which cannot
        /// separate a dense id scheme from any scheme that happens to agree at widths
        /// one to three — and the table's width is the one dimension a discovery round
        /// moves.
        #[test]
        fn every_id_admit_mints_resolves_back_to_the_allele_that_was_admitted(
            alleles in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(vec![b'A', b'C', b'G', b'T']),
                    1..8,
                ),
                0..40,
            )
        ) {
            let mut candidates = generic_reference();
            let mut minted = Vec::new();
            for allele in &alleles {
                minted.push(candidates.admit(allele.clone().into_boxed_slice()));
            }

            proptest::prop_assert_eq!(candidates.reference(), b"A");
            proptest::prop_assert_eq!(candidates.len(), alleles.len() + 1);

            for (offset, (id, allele)) in minted.iter().zip(&alleles).enumerate() {
                let expected_id = AlleleId(u16::try_from(offset + 1).expect("under 40"));
                proptest::prop_assert_eq!(*id, expected_id);
                proptest::prop_assert_eq!(candidates.bases_of(*id), Some(allele.as_slice()));
                proptest::prop_assert!(!id.is_reference());
            }

            let past_the_end = AlleleId(u16::try_from(alleles.len() + 1).expect("under 41"));
            proptest::prop_assert_eq!(candidates.bases_of(past_the_end), None);
        }
    }
    /// A site quality standing in for one the worker computed, where what the test is
    /// about is something else.
    ///
    /// **Deliberately not zero.** Zero is a real answer — a locus at which nobody carries
    /// anything — so a test whose fixture wrote zero would pass against an
    /// implementation that lost the value on the way in.
    fn a_worker_written_site_quality() -> Phred {
        Phred::try_new(37.0).expect("a legal quality, and below the site-quality ceiling")
    }

    /// A run whose repeat-tract substitution rates were never fitted — the empty map, which is
    /// what a fixture calling no tract needs, and what
    /// `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
    ///
    /// A `static` rather than a function, so that a call site can borrow it for as long as the
    /// parameters live: `BTreeMap::new` is a `const fn`, and a temporary would be freed at the
    /// end of the statement that built the view.
    static NO_SUBSTITUTION_RATES: std::collections::BTreeMap<
        crate::ng::parameter_estimation::ssr::StratumKey,
        crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
    > = std::collections::BTreeMap::new();

    /// **What a repeat tract's parameters rested on, for a fixture that is testing something
    /// else** — a tract scored over six cells, one of which found no fitted slippage and none
    /// of which was defaulted by a library the fit had never seen.
    ///
    /// **The four counts are all different from one another**, so that the record is a different
    /// value under any transposition — which matters not here, where `LocusInference::new` moves
    /// the record whole and never reads inside it, but for anything that later does. **The
    /// fixture that really discriminates the counts is the driver's**,
    /// `a_tract_reports_how_many_of_its_cells_fell_back_and_why`, which reads them off a gather.
    fn a_tract_record(rung: LengthSpectrumRung) -> RepeatTractProvenance {
        RepeatTractProvenance::new(rung, 6, 1, 0, 2, false)
    }

    fn diploid_call(first: u16, second: u16, quality: f32) -> SampleGenotypeCall {
        SampleGenotypeCall::Called {
            genotype: Genotype::new(vec![AlleleId(first), AlleleId(second)]),
            genotype_quality: Phred::try_new(quality).expect("a legal quality"),
            // These fixtures are about what a called locus holds, not about what the reads
            // said, so they use the case where the reads said something.
            reads_were_uninformative: false,
        }
    }

    fn two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
        let mut alleles = generic_reference();
        alleles.admit(Box::from(b"T".as_slice()));
        let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
        (alleles, copies)
    }

    /// A repeat tract with two candidate lengths — the path the gene-diversity seed
    /// marker belongs to, and the only one on which a locus may carry it.
    fn str_two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
        let mut alleles = CandidateAlleles::new(
            Box::from(b"ATAT".as_slice()),
            LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").expect("a dinucleotide motif"),
                left_flank: Box::from(b"CCCGGG".as_slice()),
                right_flank: Box::from(b"TTTAAA".as_slice()),
            }),
        );
        alleles.admit(Box::from(b"ATATAT".as_slice()));
        let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
        (alleles, copies)
    }

    fn region() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(3),
            start: Position(940),
            end: Position(940),
        }
    }

    /// A settled locus keeps everything a consumer needs and cannot reconstruct: the
    /// alleles finally called over, one call per sample in run order, the cohort's
    /// expected copies, and the four pieces of evidence for how the answer was reached.
    #[test]
    fn a_locus_carries_its_calls_and_the_evidence_for_how_they_were_reached() {
        let (alleles, copies) = two_allele_locus();
        let inference = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 40.0), diploid_call(0, 1, 25.0)],
            copies,
            true,
            4,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );

        assert_eq!(inference.region, region());
        assert_eq!(inference.alleles().len(), 2);
        assert_eq!(inference.cohort_expected_copies().copies(), [2.6, 1.4]);
        assert_eq!(inference.passes, 4);
        assert!(inference.converged);
        assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
        assert_eq!(inference.repeat_tract, None);

        // The per-sample calls are a sequence, in the run's sample order — the second
        // sample is the heterozygote, and swapping them would be a different record.
        assert_eq!(inference.per_sample.len(), 2);
        assert_eq!(
            inference.per_sample[0]
                .genotype()
                .expect("a called sample")
                .alleles(),
            [AlleleId(0), AlleleId(0)]
        );
        assert_eq!(
            inference.per_sample[1]
                .genotype()
                .expect("a called sample")
                .alleles(),
            [AlleleId(0), AlleleId(1)]
        );
        assert_eq!(
            inference.per_sample[1]
                .score_best_genotype()
                .expect("a called sample")
                .get(),
            25.0
        );
    }

    /// A locus that ran out of passes is emitted with the flag set, never dropped and
    /// never fatal — one hard locus must not kill a cohort run. The flag is the whole
    /// point: a genotype from a loop that did not settle is a weaker claim, and nothing
    /// downstream can tell it from a settled one otherwise.
    ///
    /// Built on the **repeat** path, because the tract ladder's rung this also carries
    /// cannot arise on the SNP/indel one — a fixture that set it there would be pinning a
    /// state its own field documents as impossible.
    #[test]
    fn a_locus_that_ran_out_of_passes_is_emitted_with_the_flag_set() {
        let (alleles, copies) = str_two_allele_locus();
        let capped = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 3.0)],
            copies,
            false,
            50,
            Provenance::Defaulted,
            Some(a_tract_record(LengthSpectrumRung::PeriodsPooledTracts)),
            a_worker_written_site_quality(),
            None,
        );

        assert!(!capped.converged, "the loop hit its cap");
        assert_eq!(capped.passes, 50);
        // The call is still there. Nothing about not converging removes it.
        assert_eq!(capped.per_sample.len(), 1);
        assert_eq!(
            capped.per_sample[0]
                .genotype()
                .expect("a called sample")
                .alleles(),
            [AlleleId(0), AlleleId(1)]
        );
        // And the two other warrants travel independently of it.
        assert_eq!(capped.weakest_provenance, Provenance::Defaulted);
        assert_eq!(
            capped.repeat_tract,
            Some(a_tract_record(LengthSpectrumRung::PeriodsPooledTracts))
        );
    }

    /// One pass is a legitimate answer and must be constructible, because the loop's
    /// first pass already has something to settle against: the reads-only estimate made
    /// before iteration begins. A bound that drifted from `> 0` to `> 1` would reject
    /// real single-pass loci at run time — and a cohort of one at three reads a position,
    /// this caller's hardest committed case, is where a loop settles fastest.
    #[test]
    fn a_locus_that_settled_on_its_first_pass_is_a_locus() {
        let (alleles, copies) = two_allele_locus();
        let quick = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 60.0)],
            copies,
            true,
            1,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
        assert_eq!(quick.passes, 1);
        assert!(quick.converged);
    }

    /// Copies of a width the allele table does not have, in **both** directions.
    ///
    /// Narrower is the obvious slip. **Wider is the one the pipeline produces**: the
    /// final prune shrinks the allele table, so a copies vector that was not re-cut
    /// alongside it is left long, and since every consumer indexes by `AlleleId` the
    /// trailing entries would ride along unread — a pruned-away allele's expected copies
    /// travelling in the output with nothing to notice.
    #[test]
    fn a_locus_cannot_carry_copies_of_a_width_its_alleles_do_not_have() {
        let build = |allele_count: usize, copy_count: usize| {
            let mut alleles = generic_reference();
            for _ in 1..allele_count {
                alleles.admit(Box::from(b"T".as_slice()));
            }
            let mut copies_alleles = generic_reference();
            for _ in 1..copy_count {
                copies_alleles.admit(Box::from(b"T".as_slice()));
            }
            let copies = ExpectedAlleleCopies::new(vec![1.0; copy_count], &copies_alleles);
            LocusInference::new(
                region(),
                alleles,
                vec![diploid_call(0, 0, 30.0)],
                copies,
                true,
                2,
                Provenance::FittedHere,
                None,
                a_worker_written_site_quality(),
                None,
            )
        };

        // Copies narrower than the alleles.
        let narrower = std::panic::catch_unwind(|| build(3, 2));
        assert!(narrower.is_err(), "3 alleles cannot carry 2 copy counts");

        // Copies wider than the alleles — what an un-re-cut vector looks like after the
        // prune, and the direction no fixture reached before.
        let wider = std::panic::catch_unwind(|| build(2, 3));
        assert!(wider.is_err(), "2 alleles cannot carry 3 copy counts");
    }

    /// Zero passes is not a locus that converged instantly — the loop always runs its
    /// body at least once — it is a counter that was never incremented, and it would
    /// distort the distribution the pass cap is to be set from.
    #[test]
    #[should_panic(expected = "at least one pass")]
    fn a_locus_cannot_report_no_passes_at_all() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 30.0)],
            copies,
            true,
            0,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
    }

    /// A locus naming no call at all has lost them: a cohort has at least one sample and
    /// every sample is called, so an empty list is a dropped result rather than a locus
    /// with nothing to say.
    #[test]
    #[should_panic(expected = "has lost them")]
    fn a_locus_cannot_name_no_sample_at_all() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            Vec::new(),
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
    }

    /// The tract ladder belongs to the repeat prior. A SNP/indel locus seeds from a
    /// *frequency* spectrum, whose own ladder has different rungs, so a locus carrying a
    /// tract record has had one path's ladder wired onto the other — which is exactly what an
    /// implementation slip in the seed's routing would look like.
    #[test]
    #[should_panic(expected = "SNP/indel locus carries a repeat-tract record")]
    fn a_snp_locus_cannot_carry_what_a_repeat_tract_rested_on() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 30.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            Some(a_tract_record(LengthSpectrumRung::PeriodsPooledTracts)),
            a_worker_written_site_quality(),
            None,
        );
    }

    /// **And the other direction, which is the newer half of the same check.** A tract whose
    /// record went missing reports no rung and no defaulted-cell counts, which downstream reads
    /// as *this call rested on nothing worth stating* rather than as a dropped field — so the
    /// record is refused as absent at a tract exactly as it is refused as present at a SNP.
    ///
    /// This direction is what a driver that built the record on the wrong branch would produce,
    /// and the fixture is a tract in every other respect.
    #[test]
    #[should_panic(expected = "repeat-tract locus carries no repeat-tract record")]
    fn a_repeat_tract_cannot_be_called_without_saying_what_it_rested_on() {
        let (alleles, copies) = str_two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 30.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
    }

    /// **A tract's record cannot claim more fallbacks than it has cells**, in either parameter.
    ///
    /// The count and its denominator are read off one table, so a larger count is two tracts'
    /// numbers in one record — and it would print as an answer rather than fail.
    #[test]
    #[should_panic(expected = "fall back on the slippage fit")]
    fn a_tract_record_cannot_default_more_slippage_cells_than_it_has() {
        let _ = RepeatTractProvenance::new(LengthSpectrumRung::StratumsOwnFit, 6, 7, 0, 0, false);
    }

    /// The same on the other parameter, checked apart because a run of one candidate per library
    /// cannot tell the two denominators from each other.
    #[test]
    #[should_panic(expected = "fall back on the substitution rate")]
    fn a_tract_record_cannot_default_more_substitution_cells_than_it_has() {
        let _ = RepeatTractProvenance::new(LengthSpectrumRung::StratumsOwnFit, 6, 0, 0, 7, false);
    }

    /// **And the share cannot exceed what it is a share of.** The cells defaulted by a library
    /// the fit never named are counted inside the arm that counts every defaulted cell, so a
    /// larger one is the two counters swapped — the one mistake a reader of these numbers could
    /// not detect, since both are plausible counts of the same table.
    #[test]
    #[should_panic(expected = "the two counters swapped")]
    fn a_tract_record_cannot_blame_more_cells_on_an_unknown_library_than_it_defaulted() {
        let _ = RepeatTractProvenance::new(LengthSpectrumRung::StratumsOwnFit, 6, 2, 3, 0, false);
    }

    /// A region that runs backwards is not a locus. `GenomeRegion` holds no constructor
    /// of its own and says a caller needing `start <= end` must say so; a called locus
    /// needs it, because the region reaches the output's position column.
    #[test]
    #[should_panic(expected = "cannot run backwards")]
    fn a_locus_cannot_cover_a_region_that_runs_backwards() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            GenomeRegion {
                contig: ContigId(3),
                start: Position(1_100),
                end: Position(940),
            },
            alleles,
            vec![diploid_call(0, 0, 30.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
    }

    /// A cohort of one is in scope, and is where this caller is weakest, so it must be
    /// representable: one sample, one call, and the cohort's expected copies are that
    /// sample's own, since it is the whole cohort.
    #[test]
    fn a_cohort_of_one_sample_is_a_locus_like_any_other() {
        let (alleles, copies) = two_allele_locus();
        let single = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 12.0)],
            copies,
            true,
            2,
            Provenance::Borrowed,
            None,
            a_worker_written_site_quality(),
            None,
        );
        assert_eq!(single.per_sample.len(), 1);
        assert_eq!(single.cohort_expected_copies().copies(), [2.6, 1.4]);
    }

    /// The nine pooled counts as the final pass builds them at a biallelic locus whose
    /// alternative drew 6 of the 8 reads — a summary of the right shape, for the two checks
    /// this type makes on it.
    fn artifact_counts_naming(primary_alternative: AlleleId) -> ArtifactTestCounts {
        ArtifactTestCounts {
            primary_alternative,
            reference_reads: 2.0,
            reference_forward_reads: 1.0,
            reference_placed_left_reads: 1.0,
            alternative_reads: 6.0,
            alternative_forward_reads: 5.0,
            alternative_placed_left_reads: 2.0,
            total_reads: 8.0,
            genotype_expected_alternative_reads: 4.0,
        }
    }

    /// **Both of the quality plan's fields travel onto the record**, and the site quality is
    /// read back through the one accessor there is.
    ///
    /// `doc/devel/ng/spec/calling_quality.md` §3.5 makes that single accessor the design: the
    /// worker writes the baseline, the artifact-correction stage overwrites the same field
    /// with the corrected value, and **there is never a second quality field** for a gate and
    /// a column to disagree about. Production shipped the other shape for sixteen days and
    /// emitted 40 false positives at 30× on GIAB HG002 carrying a written `QUAL` of 0.
    #[test]
    fn the_site_quality_and_the_artifact_summary_travel_onto_the_record() {
        let (alleles, copies) = two_allele_locus();
        let inference = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 25.0)],
            copies,
            true,
            3,
            Provenance::FittedHere,
            None,
            Phred::try_new(431.5).expect("a legal quality"),
            Some(artifact_counts_naming(AlleleId(1))),
        );

        assert_eq!(inference.uncorrected_site_quality().get(), 431.5);
        let counts = inference
            .artifact_test_counts()
            .expect("this locus carries a summary");
        assert_eq!(counts.primary_alternative, AlleleId(1));
        assert_eq!(counts.total_reads, 8.0);
    }

    /// A locus with no summary carries `None` rather than a summary of zeroes — the two are
    /// different claims, and only one of them can be told from *the tests found no bias*.
    #[test]
    fn a_locus_with_nothing_for_the_artifact_tests_to_weigh_carries_no_summary() {
        let (alleles, copies) = two_allele_locus();
        let inference = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 25.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
        assert!(inference.artifact_test_counts().is_none());
    }

    /// A site quality above the ceiling did not come from the function that owns the ceiling.
    ///
    /// `quality::score_uncorrected_site_quality` caps at `MAX_SITE_QUALITY`, which is also
    /// its answer to a probability of exactly zero, so nothing it returns can be above it.
    /// The field this is going into is the one the correction stage overwrites in place, and
    /// a number from somewhere else in it is the shape §3.5 exists to prevent.
    #[test]
    #[should_panic(expected = "above the ceiling")]
    fn a_site_quality_above_the_ceiling_is_refused() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 25.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            Phred::try_new(quality::MAX_SITE_QUALITY + 1.0).expect("a legal Phred"),
            None,
        );
    }

    /// **The ceiling itself is accepted**, and the boundary is not academic: it is exactly
    /// what `quality::score_uncorrected_site_quality` returns for a probability of no variant
    /// of zero — a locus whose reads exclude the reference outright. A `<`/`<=` slip in the
    /// check above would panic the worker on the loci the caller is most confident about.
    #[test]
    fn a_site_quality_at_the_ceiling_is_accepted() {
        let (alleles, copies) = two_allele_locus();
        let inference = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(1, 1, 25.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            Phred::try_new(quality::MAX_SITE_QUALITY).expect("a legal Phred"),
            None,
        );
        assert_eq!(
            inference.uncorrected_site_quality().get(),
            quality::MAX_SITE_QUALITY
        );
    }

    /// An artifact summary that names the **reference** as its primary alternative is
    /// refused: both tests weigh one alternative against the reference, so a locus with
    /// nothing to weigh carries no summary at all rather than one pointing at allele 0.
    #[test]
    #[should_panic(expected = "has no comparison to make")]
    fn an_artifact_summary_naming_the_reference_as_its_alternative_is_refused() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 25.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            Some(artifact_counts_naming(AlleleId::REFERENCE)),
        );
    }

    /// An artifact summary naming an allele this locus was not called over is refused, for
    /// the reason a stale genotype is: the final prune renumbers every id above the one it
    /// drops, so a summary not re-cut alongside the table points somewhere else.
    #[test]
    #[should_panic(expected = "pooled against a different allele table")]
    fn an_artifact_summary_naming_an_allele_the_locus_lacks_is_refused() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles, // two alleles: ids 0 and 1
            vec![diploid_call(0, 1, 25.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            Some(artifact_counts_naming(AlleleId(4))),
        );
    }

    // ────────────────────────────────────────────────────────────────────────────────
    // The three types the calling loop takes and gives back, and the missing genotype.
    // ────────────────────────────────────────────────────────────────────────────────

    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
    use std::collections::BTreeMap;

    /// Two candidates' repeat counts, **different from each other**, for a tract fixture whose
    /// candidate table holds two lengths. Distinct so that a lookup reading one candidate's
    /// count for another is a different stratum rather than the same one.
    fn two_repeat_counts() -> Vec<NonZeroU32> {
        vec![
            NonZeroU32::new(6).expect("six repeats"),
            NonZeroU32::new(7).expect("seven repeats"),
        ]
    }

    fn ssr_detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }
    }

    fn no_strata() -> StratumFits {
        StratumFits::over(&[], BTreeMap::new())
    }

    /// A run whose fit produced a length spectrum at one dinucleotide stratum and at one
    /// trinucleotide stratum, with **different shapes and different concentrations** — so that a
    /// lookup answering from the wrong stratum, or from the wrong period, is a different number
    /// rather than the same one.
    fn strata_with_length_spectra() -> StratumFits {
        use crate::ng::parameter_estimation::joint::census::Stratum;
        use crate::ng::parameter_estimation::joint::share_curve::ShareSource;
        use crate::ng::parameter_estimation::joint::slippage_curve::LevelSource;
        use crate::ng::parameter_estimation::joint::ssr_fit::{
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
        let one = |period: u8, repeats: u64, length_spectrum: Vec<f64>, concentration: f64| {
            StratumOutcome::Fitted(Box::new(StratumFit {
                stratum: Stratum {
                    period,
                    reference_repeats: repeats,
                },
                slippage: vec![Some(Slippage {
                    level: 0.05,
                    shorter_share: 0.83,
                    fall_off: 0.25,
                })],
                length_spectrum,
                concentration,
                log_likelihood_a_tract: -1.5,
                tracts_fitted: 40,
                borrowed: Vec::new(),
                converged: true,
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
        use crate::ng::parameter_estimation::joint::stratum_fits::LengthSpectrumRung;
        use crate::ng::parameter_estimation::ssr::RepeatCount;

        let strata = strata_with_length_spectra();
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

        let dinucleotide = parameters.ssr_length_spectrum_at(
            SsrPeriod::try_new(2).expect("a dinucleotide"),
            RepeatCount(10),
        );
        assert_eq!(dinucleotide.rung(), LengthSpectrumRung::StratumsOwnFit);
        assert_eq!(dinucleotide.fitted_weights(), Some(&[0.6, 0.3, 0.1][..]));
        assert_eq!(dinucleotide.concentration(), 4.0);

        // Same repeat count, other motif period: a different stratum and a different answer.
        let trinucleotide = parameters.ssr_length_spectrum_at(
            SsrPeriod::try_new(3).expect("a trinucleotide"),
            RepeatCount(10),
        );
        assert_eq!(trinucleotide.fitted_weights(), Some(&[0.1, 0.3, 0.6][..]));
        assert_eq!(trinucleotide.concentration(), 25.0);

        // A repeat count neither stratum holds falls down the ladder rather than answering.
        let elsewhere = parameters.ssr_length_spectrum_at(
            SsrPeriod::try_new(2).expect("a dinucleotide"),
            RepeatCount(11),
        );
        assert_eq!(elsewhere.rung(), LengthSpectrumRung::StatedFlat);
        assert_eq!(
            elsewhere.concentration(),
            14.5,
            "the run fitted 4.0 and 25.0, whose median is their mean"
        );
    }

    fn diploid_ploidy() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
    }

    fn outbred_samples(samples: usize) -> Vec<InbreedingF> {
        vec![InbreedingF::try_new(0.0).expect("a legal coefficient"); samples]
    }

    fn one_read_group() -> Vec<ReadGroupCalibration> {
        vec![ReadGroupCalibration::defaulted()]
    }

    fn measured_contamination() -> ContaminationView {
        ContaminationView {
            fraction: 0.03,
            markers_with_reads: 400,
            reads_on_markers: 1_200,
            source: ContaminationSource::ThisReadGroupsReads,
        }
    }

    fn neutral_seed() -> SpectrumSeed {
        SpectrumSeed::new(1.0, 1e-3, genotype_prior::SeedRegime::NeutralShape)
    }

    /// The frozen parameters for a run of diploids, routed to whichever constructor the
    /// contamination list calls for — an empty one is the run where the fit identified
    /// none, which has its own named door.
    ///
    /// **The contaminated branch takes the batching as an argument**, because it is borrowed
    /// for as long as the parameters are and a helper that built one inside would hand back a
    /// borrow of something already dropped. A run with no mixture has nothing to read it, so
    /// the uncontaminated branch never looks at it.
    fn frozen_parameters<'a>(
        calibration: &'a [ReadGroupCalibration],
        contamination: &'a [ContaminationView],
        batching: &'a SequencingBatches,
        inbreeding: &'a [InbreedingF],
        strata: &'a StratumFits,
        substitution: &'a std::collections::BTreeMap<
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
        >,
    ) -> FrozenParameters<'a> {
        if contamination.is_empty() {
            FrozenParameters::uncontaminated(
                calibration,
                inbreeding,
                neutral_seed(),
                strata,
                substitution,
                diploid_ploidy(),
            )
        } else {
            FrozenParameters::new(
                calibration,
                contamination,
                batching,
                inbreeding,
                neutral_seed(),
                strata,
                substitution,
                diploid_ploidy(),
            )
        }
    }

    /// The default batching — one batch holding a run of `libraries` read groups over
    /// `samples` samples, which is what a run that declared nothing gets.
    pub(crate) fn one_batch(libraries: usize, samples: usize) -> SequencingBatches {
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

    /// A two-allele diploid locus: its allele table and the genotype table over it. Three
    /// genotypes against two alleles, so a scratch that swapped the two counts is caught.
    fn biallelic_locus() -> (CandidateAlleles, std::sync::Arc<GenotypeTable>) {
        let (alleles, _) = two_allele_locus();
        let table = GenotypeTable::build(diploid_ploidy(), alleles.len());
        (alleles, table)
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

    /// Evidence and the allele table have to be on the same path: the discriminant that
    /// chose the candidates is the one that chooses the read model.
    ///
    /// **Handing a repeat tract's alleles a SNP/indel's evidence would not crash** — the
    /// generic row would score the tract's sequences as though they were substitutions,
    /// giving a different likelihood at every sample and a plausible genotype at the end of
    /// it.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn evidence_on_the_wrong_path_for_its_allele_table_is_refused() {
        let (tract_alleles, _) = str_two_allele_locus();
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
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

        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&tract_alleles, &parameters);
    }

    /// The evidence and the run's per-sample parameters are indexed by one sample order, so
    /// two lists of different lengths are two different orders.
    #[test]
    #[should_panic(expected = "indexed by a different sample order")]
    fn evidence_covering_a_different_sample_count_from_the_run_is_refused() {
        let (alleles, _) = two_allele_locus();
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

        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&alleles, &parameters);
    }

    /// Both paths pass their own check, and each carries the locus's region and the run's
    /// sample count.
    ///
    /// **A repeat tract's evidence carries no missing-genotype flag at all**, and the
    /// absence is the design: a tract's discovery round can put back a length the cap cut,
    /// so no sample is set aside there.
    #[test]
    fn evidence_that_matches_its_locus_and_its_run_is_accepted_on_both_paths() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(2);
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

        let (generic_alleles, _) = two_allele_locus();
        let generic_samples = [
            GenericLocusSample {
                evidence: GenericSampleEvidence::empty(),
                genotype_must_be_missing: false,
            },
            GenericLocusSample {
                evidence: GenericSampleEvidence::empty(),
                genotype_must_be_missing: true,
            },
        ];
        let generic = LocusEvidence::generic(region(), &generic_samples);
        generic.assert_matches_locus_and_run(&generic_alleles, &parameters);
        assert_eq!(generic.region(), region());
        assert_eq!(generic.sample_count(), 2);
        match generic {
            LocusEvidence::Generic { per_sample, .. } => {
                assert!(!per_sample[0].genotype_must_be_missing);
                assert!(per_sample[1].genotype_must_be_missing);
            }
            LocusEvidence::Ssr { .. } => panic!("built on the SNP/indel path"),
        }

        let (tract_alleles, _) = str_two_allele_locus();
        let detail = ssr_detail();
        let tract_samples = [
            SsrSampleEvidence::new(&[], &detail),
            SsrSampleEvidence::new(&[], &detail),
        ];
        let counts = two_repeat_counts();
        let tract = LocusEvidence::ssr(region(), &tract_samples, &detail, &counts);
        tract.assert_matches_locus_and_run(&tract_alleles, &parameters);
        assert_eq!(tract.region(), region());
        assert_eq!(tract.sample_count(), 2);
    }

    /// A locus whose evidence list is empty is evidence that went missing, not a locus
    /// nobody covered — a sample with no reads gets an empty *entry*, and the prior decides
    /// its genotype alone.
    #[test]
    #[should_panic(expected = "evidence that went missing")]
    fn evidence_naming_no_sample_at_all_is_refused() {
        let _ = LocusEvidence::generic(region(), &[]);
    }

    /// Each sample reads and writes its own row of the two flat tables, and the rows do not
    /// overlap.
    ///
    /// **The flat tables are where a slip is silent.** A caller slicing `lg_table` itself
    /// with the allele count instead of the genotype count would read a window straddling
    /// two samples' rows: every entry a real log-likelihood, none of them this sample's.
    /// The fixture writes each sample's index into its own row so a wrong window names the
    /// sample it actually came from.
    #[test]
    fn each_sample_reads_and_writes_its_own_row_of_the_scratch_tables() {
        let (alleles, table) = biallelic_locus();
        let view = table.view();
        assert_eq!(view.genotype_count(), 3);
        assert_eq!(view.allele_count(), 2);

        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(3, &alleles, &view);
        assert_eq!(scratch.row_count(), 3);
        assert_eq!(scratch.genotype_count(), 3);
        assert_eq!(scratch.allele_count(), 2);

        for sample in 0..3 {
            for (genotype, slot) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .enumerate()
            {
                *slot = LogProb(-(sample as f64) - 0.1 * genotype as f64);
            }
            for (allele, slot) in scratch
                .sample_expected_copies_mut(sample)
                .iter_mut()
                .enumerate()
            {
                *slot = sample as f64 + 0.5 * allele as f64;
            }
        }

        for sample in 0..3 {
            let row = scratch.sample_genotype_likelihoods(sample);
            assert_eq!(row.len(), 3);
            assert_eq!(row[0].get(), -(sample as f64));
            assert_eq!(row[2].get(), -(sample as f64) - 0.2);

            let copies = scratch.sample_expected_copies(sample);
            assert_eq!(copies.len(), 2);
            assert_eq!(copies, [sample as f64, sample as f64 + 0.5]);
        }
    }

    /// **[`CallingScratch::sample_scoring_buffers_mut`] hands back that sample's two rows**
    /// of the flat sample-major tables — its likelihood row and its expected-copies row —
    /// and not another sample's.
    ///
    /// **The scoring path's own tests cannot see this**, which is why the check lives here.
    /// Measured: with the method's two `sample` arguments replaced by `0`, so every sample
    /// is handed sample 0's rows, the whole of `ng::calling` still passes. The test that
    /// looks like the one to catch it — the three-sample scoring test, whose `assert_ne!`
    /// is commented *"sample 1 was scored on sample 0's likelihood row"* — passes under
    /// that mutation, because it also redirects the expected-copies row, which the scorer
    /// **overwrites** on every call; so each sample gets a different leave-one-out term and
    /// therefore a different posterior, from the same reads.
    ///
    /// The two rows are 3 and 2 entries wide here, so a wrong row is a legal slice of the
    /// right table — the hazard [`row_range`](CallingScratch::row_range)'s own doc names.
    #[test]
    fn the_scoring_buffers_hand_back_that_samples_own_rows() {
        let (alleles, table) = biallelic_locus();
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(3, &alleles, &view);
        for sample in 0..3 {
            for (genotype, slot) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .enumerate()
            {
                *slot = LogProb(10.0 * sample as f64 + genotype as f64);
            }
            scratch
                .sample_expected_copies_mut(sample)
                .copy_from_slice(&[sample as f64, 2.0 - sample as f64]);
        }

        let buffers = scratch.sample_scoring_buffers_mut(1);
        assert_eq!(buffers.sample, 1);
        assert_eq!(
            buffers.genotype_likelihoods,
            &[LogProb(10.0), LogProb(11.0), LogProb(12.0)],
            "sample 1's likelihood row"
        );
        assert_eq!(
            buffers.sample_expected_copies,
            &[1.0, 1.0],
            "sample 1's own expected copies"
        );
    }

    /// **The buffers a repeat tract's rows are scored from refuse a scratch whose rows were
    /// never claimed**, and refuse it at the door rather than inside the walk.
    ///
    /// The walk over a tract's rows runs inside that one borrow, where no accessor of this
    /// scratch can be called — so the map from a row back to the run's sample has to be checked
    /// on the way in. Without the check, a scratch prepared for two rows and claimed for one
    /// reads row 1 of a map one entry long, which is an index panic naming a slice where the
    /// reader needs the cohort.
    #[test]
    #[should_panic(expected = "of them were claimed")]
    fn the_tracts_locus_buffers_refuse_a_scratch_whose_rows_were_never_claimed() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &table.view());
        scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
        let _ = scratch.tract_locus_buffers_mut();
    }

    /// A sample past the count the scratch was prepared for is refused by name, rather than
    /// reading the next sample's row.
    #[test]
    #[should_panic(expected = "sample 2 is past the 2 this scratch was prepared for")]
    fn a_sample_past_the_prepared_count_is_refused() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &table.view());
        let _ = scratch.sample_genotype_likelihoods(2);
    }

    /// **Preparing a locus overwrites every entry, even when the shape has not changed.**
    ///
    /// This is the failure `Vec::resize` alone would leave: at two loci of the same shape it
    /// keeps the leading entries, so the second locus would be scored against the first
    /// one's likelihoods and priors with nothing failing. Every buffer comes back `NaN`, so
    /// a value that survived is one that reaches an arithmetic check rather than a genotype.
    #[test]
    fn preparing_a_locus_overwrites_the_previous_locus_of_the_same_shape() {
        let (alleles, table) = biallelic_locus();
        let view = table.view();

        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &view);
        scratch.sample_genotype_likelihoods_mut(0)[0] = LogProb(-7.0);
        scratch.sample_genotype_likelihoods_mut(1)[2] = LogProb(-9.0);
        scratch.sample_expected_copies_mut(0)[1] = 1.75;
        scratch.cohort_expected_copies_mut()[0] = 3.5;
        scratch.seed_concentration_mut()[0] = 1.0;

        scratch.prepare_for_locus(2, &alleles, &view);

        assert!(
            scratch
                .sample_genotype_likelihoods(0)
                .iter()
                .all(|entry| entry.get().is_nan())
        );
        assert!(
            scratch
                .sample_genotype_likelihoods(1)
                .iter()
                .all(|entry| entry.get().is_nan())
        );
        assert!(
            scratch
                .sample_expected_copies(0)
                .iter()
                .all(|copy| copy.is_nan())
        );
        assert!(
            scratch
                .cohort_expected_copies()
                .iter()
                .all(|copy| copy.is_nan())
        );
        assert!(scratch.seed_concentration().iter().all(|a| a.is_nan()));
    }

    /// Advancing hands the previous pass's copies to the convergence test and gives back a
    /// buffer for the next pass — a swap, so no allele-length copy is paid per pass.
    #[test]
    fn advancing_makes_the_current_copies_the_previous_ones() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &table.view());

        scratch
            .cohort_expected_copies_mut()
            .copy_from_slice(&[1.25, 0.75]);
        let next = scratch.advance_cohort_expected_copies();
        // Handed back unwritten, not holding the pass before last's real numbers: a pass
        // that skipped an allele would otherwise leave a value no pass wrote for the
        // convergence test to compare.
        assert!(next.iter().all(|copy| copy.is_nan()), "{next:?}");
        next.copy_from_slice(&[1.5, 0.5]);

        assert_eq!(scratch.previous_cohort_expected_copies(), [1.25, 0.75]);
        assert_eq!(scratch.cohort_expected_copies(), [1.5, 0.5]);
    }

    /// A locus prepared for no rows at all is a run whose sample order went missing. **Rows
    /// rather than samples**: the scratch is sized for the samples the locus is *called on*,
    /// so zero here is a locus nobody can be called at rather than a cohort of nobody.
    #[test]
    #[should_panic(expected = "a locus prepared for no rows")]
    fn a_scratch_prepared_for_no_rows_is_refused() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(0, &alleles, &table.view());
    }

    /// A missing genotype has no genotype and no quality — the two are absent together,
    /// which is what an enum buys over a struct with an optional field.
    ///
    /// **A quality beside a missing genotype would be the failure.** Emission has to write
    /// this sample's `GT` as missing because no genotype was scored for it, not because a
    /// scored one came out weak, and the two must not be conflated in the output.
    #[test]
    fn a_missing_call_carries_neither_a_genotype_nor_a_quality() {
        let missing = SampleGenotypeCall::Missing;
        assert!(missing.is_missing());
        assert_eq!(missing.genotype(), None);
        assert_eq!(missing.score_best_genotype(), None);

        let called = diploid_call(0, 1, 25.0);
        assert!(!called.is_missing());
        assert_eq!(
            called.genotype().expect("a called sample").alleles(),
            [AlleleId(0), AlleleId(1)]
        );
        assert_eq!(
            called.score_best_genotype().expect("a called sample").get(),
            25.0
        );
    }

    /// A locus can hand back a missing call beside called ones, in the run's sample order —
    /// which is what a cohort whose cap cut one sample's earned allele produces.
    #[test]
    fn a_locus_carries_a_missing_call_beside_called_ones_in_run_order() {
        let (alleles, copies) = two_allele_locus();
        let inference = LocusInference::new(
            region(),
            alleles,
            vec![
                diploid_call(0, 0, 40.0),
                SampleGenotypeCall::Missing,
                diploid_call(0, 1, 25.0),
            ],
            copies,
            true,
            3,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );

        assert!(!inference.per_sample[0].is_missing());
        assert!(inference.per_sample[1].is_missing());
        assert!(!inference.per_sample[2].is_missing());
    }

    /// A repeat **bundle**'s allele table takes repeat-tract evidence, exactly as a tract's
    /// does — which is how every other consumer of `LocusKind` groups the two
    /// (`src/ng/run/cohort_merge/close.rs`).
    ///
    /// **No fixture reached this cell before, and two mutations lived there.** Deleting
    /// `| LocusKind::SsrBundle` from the accepting arm left the whole suite green, and so
    /// did widening the SNP/indel arm to accept a bundle — the second of which routes a
    /// repeat bundle to the SNP/indel row: a plausible genotype at every sample, nothing
    /// failing.
    #[test]
    fn ssr_evidence_against_a_bundle_allele_table_is_accepted() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
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

        let bundle = CandidateAlleles::new(Box::from(b"CAGCAG".as_slice()), LocusKind::SsrBundle);
        let detail = ssr_detail();
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
        // One repeat count, because a bundle built from a reference alone is called over one
        // allele — the counts run parallel to the candidate table, which the same method checks.
        let counts = [NonZeroU32::new(6).expect("six repeats")];
        let evidence = LocusEvidence::ssr(region(), &per_sample, &detail, &counts);
        evidence.assert_matches_locus_and_run(&bundle, &parameters);
    }

    /// The other half of the same cell: SNP/indel evidence at a repeat bundle is refused.
    ///
    /// This is the mutation that bites — an accepting arm widened to admit a bundle scores
    /// a repeat locus with the SNP/indel read model, and nothing anywhere fails.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn generic_evidence_against_a_bundle_allele_table_is_refused() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
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

        let bundle = CandidateAlleles::new(Box::from(b"CAGCAG".as_slice()), LocusKind::SsrBundle);
        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&bundle, &parameters);
    }

    /// The repeat constructor's own empty-list guard, which its SNP/indel twin's test does
    /// not reach.
    ///
    /// **Measured: this was the one assertion in the module no test noticed.** Downgrading
    /// all sixteen of the implementation's release-held checks to `debug_assert!` and
    /// running under `--release` failed a test for fifteen of them; this was the sixteenth.
    /// The two constructors are a copy-paste pair, which is the shape in which one half
    /// drifts from the other.
    #[test]
    #[should_panic(expected = "evidence that went missing")]
    fn ssr_evidence_naming_no_sample_at_all_is_refused() {
        let detail = ssr_detail();
        let _ = LocusEvidence::ssr(region(), &[], &detail, &two_repeat_counts());
    }

    /// **A tract with no candidate repeat counts at all is refused**, and it is a different
    /// failure from the empty sample list beside it: a tract is called over at least its
    /// reference length, so a run that supplied none lost the counts on the way in.
    ///
    /// **Refused here rather than left to the parameter assembly**, which would otherwise ask
    /// for a table of `read groups × 0` cells and refuse it naming the table rather than the
    /// locus.
    #[test]
    #[should_panic(expected = "supplied no candidate repeat counts")]
    fn ssr_evidence_with_no_candidate_repeat_counts_is_refused() {
        let detail = ssr_detail();
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
        let _ = LocusEvidence::ssr(region(), &per_sample, &detail, &[]);
    }

    /// **Repeat counts that are not one per candidate are refused where the two still name each
    /// other**, rather than at whichever genotype first read the wrong cell.
    ///
    /// The counts run parallel to the candidate table, so a table with more candidates than
    /// counts would pair one candidate's bases with another candidate's stratum — a wrong
    /// slippage model per candidate, with a well-formed genotype coming back and nothing
    /// failing.
    #[test]
    #[should_panic(expected = "called over 2 candidates and 1 repeat counts")]
    fn a_tract_whose_repeat_counts_are_not_one_per_candidate_is_refused() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let batching = SequencingBatches::all_together_over(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let detail = ssr_detail();
        let mut alleles = CandidateAlleles::new(
            Box::from(b"ATATAT".as_slice()),
            LocusKind::Ssr(ssr_detail()),
        );
        alleles.admit(Box::from(b"ATATATAT".as_slice()));
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
        // Two candidates, one count.
        let counts = [NonZeroU32::new(6).expect("six repeats")];
        let evidence = LocusEvidence::ssr(region(), &per_sample, &detail, &counts);
        evidence.assert_matches_locus_and_run(&alleles, &parameters);
    }

    /// The third class the constructor's own documentation refuses, and the only one that
    /// had no test.
    ///
    /// **Infinity is the arithmetic-gone-wrong shape a log-domain sum produces**, and it is
    /// the one that survives a check weakened from `is_finite()` to `!is_nan()`. It would
    /// ride into the locus's output, and the convergence comparison of two infinities is
    /// itself not a number — the failure the not-a-number check exists to stop, arriving
    /// through the untested door.
    #[test]
    #[should_panic(expected = "finite and at or above zero")]
    fn expected_copies_reject_an_infinite_count() {
        let mut candidates = generic_reference();
        candidates.admit(Box::from(b"T".as_slice()));
        let _ = ExpectedAlleleCopies::new(vec![1.0, f64::INFINITY], &candidates);
    }

    /// A scratch that was never sized for a locus is refused by name, rather than handing a
    /// pass an empty buffer to fold nothing into.
    ///
    /// **The unwritten-value fill cannot cover this door**, and that is why the check
    /// exists: an unprepared buffer has no slots to fill, so a fold over it runs zero
    /// iterations, writes nothing, and leaves the cohort's expected copies summing to a
    /// plausible `0.0` — which is exactly the value the fill exists to keep out of a
    /// genotype. Step B2's oracle is a bitwise comparison of those sums, and an all-zero
    /// sum passes against another all-zero sum.
    #[test]
    #[should_panic(expected = "has not been prepared for a locus")]
    fn an_unprepared_scratch_is_refused() {
        let mut scratch = CallingScratch::<()>::default();
        for slot in scratch.cohort_expected_copies_mut() {
            *slot = 99.0;
        }
    }

    /// A genotype table built for a different allele set is refused where the locus's shape
    /// is fixed, rather than at the far end of the locus.
    ///
    /// **This is the first of the three caller bugs the spec names as assertions**, and a
    /// discovery round admitting an allele is exactly how the two come apart. With a table
    /// one allele narrow, every per-allele buffer is sized for the old count, and the first
    /// thing to notice is the copies-versus-alleles check when the locus is finished —
    /// thousands of arithmetic operations later, naming the output rather than the table.
    #[test]
    #[should_panic(expected = "the table was built for a different allele set")]
    fn a_genotype_table_built_for_another_allele_set_is_refused() {
        let (mut alleles, table) = biallelic_locus();
        // A discovery round admits a third allele; the table is still the two-allele one.
        alleles.admit(Box::from(b"G".as_slice()));

        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &table.view());
    }

    /// The previous pass's copies are poisoned by the next locus too, not only the current
    /// pass's.
    ///
    /// **Without this, locus *n*'s first convergence comparison runs against locus *n−1*'s
    /// final allele copies** — a locus declared settled on the previous locus's numbers,
    /// with nothing failing. The fixture writes the previous-pass buffer between the two
    /// loci, which is the one state the same-shape test does not reach.
    #[test]
    fn preparing_a_locus_poisons_the_previous_passs_copies_as_well() {
        let (alleles, table) = biallelic_locus();
        let view = table.view();

        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &view);
        scratch
            .cohort_expected_copies_mut()
            .copy_from_slice(&[1.25, 0.75]);
        scratch
            .advance_cohort_expected_copies()
            .copy_from_slice(&[1.5, 0.5]);
        assert_eq!(scratch.previous_cohort_expected_copies(), [1.25, 0.75]);

        // The next locus, same shape.
        scratch.prepare_for_locus(2, &alleles, &view);
        assert!(
            scratch
                .previous_cohort_expected_copies()
                .iter()
                .all(|copy| copy.is_nan()),
            "the previous pass's copies survived into the next locus: {:?}",
            scratch.previous_cohort_expected_copies()
        );
    }

    /// The four buffers no accessor reached are sized on the axis they belong to.
    ///
    /// **Three genotypes against two alleles is what makes this discriminating.** A buffer
    /// sized on the wrong axis is still a legal length, so at a locus where the two counts
    /// were equal every wrong sizing would pass; here the prior and posterior rows must be
    /// 3 and the two per-allele buffers must be 2.
    #[test]
    fn the_prior_and_posterior_buffers_are_sized_on_their_own_axes() {
        let (alleles, table) = biallelic_locus();
        let view = table.view();

        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &view);

        assert_eq!(scratch.prior_row().len(), view.genotype_count());
        assert_eq!(scratch.posterior_row().len(), view.genotype_count());
        assert_eq!(scratch.sample_concentration().len(), view.allele_count());
        assert_eq!(
            scratch.prior_per_allele_workspace().len(),
            view.allele_count()
        );
        assert_ne!(view.genotype_count(), view.allele_count());

        // And they arrive unwritten, like every other buffer.
        assert!(scratch.prior_row().iter().all(|entry| entry.get().is_nan()));
        assert!(scratch.posterior_row().iter().all(|entry| entry.is_nan()));
        assert!(scratch.sample_concentration().iter().all(|a| a.is_nan()));
        assert!(
            scratch
                .prior_per_allele_workspace()
                .iter()
                .all(|a| a.is_nan())
        );

        // Written, then overwritten by the next locus of the same shape.
        scratch.prior_row_mut()[0] = LogProb(-1.5);
        scratch.posterior_row_mut()[0] = 0.25;
        scratch.sample_concentration_mut()[0] = 1.0;
        scratch.prior_per_allele_workspace_mut()[0] = 2.0;
        scratch.prepare_for_locus(2, &alleles, &view);
        assert!(scratch.prior_row()[0].get().is_nan());
        assert!(scratch.posterior_row()[0].is_nan());
        assert!(scratch.sample_concentration()[0].is_nan());
        assert!(scratch.prior_per_allele_workspace()[0].is_nan());
    }

    /// The configuration a run actually uses — the shipped repeat-tract emission model's
    /// own scratch — is built here and nowhere else in the module's tests.
    ///
    /// **Every other fixture uses `()` as the model scratch**, which compiles whether or
    /// not the shipped model's `Scratch` still satisfies what the type needs. This also
    /// gives the three sub-scratch accessors their only call site, and pins that a scratch
    /// can cross a thread boundary, which "allocated once per worker" requires and nothing
    /// else checked.
    #[test]
    fn the_shipped_emission_scratch_builds_and_can_cross_a_worker_boundary() {
        fn assert_send<T: Send>() {}
        assert_send::<CallingScratch<likelihood::ssr_emission::StutterSubstitutionScratch>>();

        let (alleles, table) = biallelic_locus();
        let mut scratch =
            CallingScratch::<likelihood::ssr_emission::StutterSubstitutionScratch>::default();
        scratch.prepare_for_locus(2, &alleles, &table.view());
        assert_eq!(scratch.row_count(), 2);

        // The three sub-scratches are reachable, and are the worker's rather than a locus's.
        let _ = scratch.candidate_selection_mut();
        let _ = scratch.generic_row_mut(0);
        let _ = scratch.ssr_row_mut();
    }

    /// **The contaminant tables refuse to be read before they are sized**, and that refusal is
    /// what the sizing's own separateness costs.
    ///
    /// `prepare_for_locus` un-sizes them and only a run that fitted a fraction re-sizes them, so
    /// a driver that forgot the second call would otherwise sum every batch's copies into a
    /// zero-length table: the fill writes nothing, returns having done nothing wrong, and the
    /// mixture is built over an empty slice.
    #[test]
    #[should_panic(expected = "contaminant tables were not prepared")]
    fn the_contaminant_tables_refuse_to_be_read_before_they_are_sized() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &table.view());
        scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
        let _ = scratch.batch_copy_buffers_mut();
    }

    /// A locus prepared for no sequencing batches is a run whose batching went missing rather
    /// than a run without one — the default is one batch holding all of it.
    #[test]
    #[should_panic(expected = "at least one sequencing batch")]
    fn contaminant_tables_for_no_batches_are_refused() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &table.view());
        scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
        scratch.prepare_contaminant_tables(0, 1);
    }

    /// **A row naming a sample the run does not have is not this run's locus**, and the copies
    /// are scattered onto the run's axis by that name — so it would put one sample's copies
    /// outside the table they are indexed in.
    ///
    /// **The row count is not what decides it.** A one-row locus whose row names sample 5 of a
    /// two-sample run has fewer rows than the run has samples and is still wrong, which is what
    /// this fixture is: one row, a run of two, and the row names sample 5.
    #[test]
    #[should_panic(expected = "rows name sample 5 and the run has 2")]
    fn a_row_naming_a_sample_the_run_does_not_have_is_refused() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &table.view());
        scratch.claim_row_for(5, InbreedingF::try_new(0.0).expect("an outbred sample"));
        scratch.prepare_contaminant_tables(1, 2);
    }

    /// One emission cache per scratch row, and a row past the locus's own count is refused —
    /// which is what a driver walking one row too far would otherwise fill, leaving the last
    /// row's cache holding another sample's observations.
    #[test]
    #[should_panic(expected = "row 2 is past the 2 this scratch was prepared for")]
    fn an_emission_cache_past_the_prepared_rows_is_a_caller_bug() {
        let (alleles, table) = biallelic_locus();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &table.view());
        let _ = scratch.generic_row_mut(2);
    }

    /// **The contaminant tables are un-sized by every `prepare_for_locus`**, so a worker that
    /// called a contaminated locus and then an uncontaminated one cannot read the first one's
    /// frequencies at the second.
    ///
    /// Sizing them and then preparing another locus is the whole of the fixture: `resize` alone
    /// would have left a table of the right length holding the last locus's numbers, and only
    /// the count coming back to zero says the tables were released rather than reused.
    #[test]
    fn preparing_a_locus_releases_the_last_ones_contaminant_tables() {
        let (alleles, table) = biallelic_locus();
        let outbred = InbreedingF::try_new(0.0).expect("an outbred sample");
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &table.view());
        scratch.claim_row_for(0, outbred);
        scratch.prepare_contaminant_tables(3, 4);
        assert_eq!(scratch.contaminant_batch_count(), 3);

        scratch.prepare_for_locus(1, &alleles, &table.view());
        assert_eq!(
            scratch.contaminant_batch_count(),
            0,
            "the next locus has to size them again, and an uncontaminated one never will"
        );
    }

    /// A repeat tract sets no sample aside, so a missing call at one is the SNP/indel path's
    /// ruling wired onto the wrong path.
    ///
    /// **The mirror of the tract-ladder rung's check**, which refuses the repeat path's
    /// ruling at a SNP/indel locus. What this one catches is a sample silently dropped from
    /// a tract's output: such a sample has no expected copies at all, so it also leaves the
    /// cohort's expected-copies denominator, and the locus's allele frequencies come out
    /// wrong with a well-formed record beside them.
    #[test]
    #[should_panic(expected = "a repeat tract sets no sample aside")]
    fn a_repeat_tract_locus_cannot_carry_a_missing_call() {
        let (alleles, copies) = str_two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 1, 20.0), SampleGenotypeCall::Missing],
            copies,
            true,
            2,
            Provenance::FittedHere,
            Some(a_tract_record(LengthSpectrumRung::StratumsOwnFit)),
            a_worker_written_site_quality(),
            None,
        );
    }

    /// A call naming an allele the locus was not called over is refused.
    ///
    /// **The prune is what produces one.** Dropping an allele renumbers every id above it,
    /// so a call minted before the prune and not remapped alongside it names an id past the
    /// end — which reaches the `GT` column as an index past the `ALT` list, an unparseable
    /// record rather than a crash.
    #[test]
    #[should_panic(expected = "names an allele this locus was not called over")]
    fn a_call_naming_an_allele_the_locus_lost_is_refused() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles, // two alleles: ids 0 and 1
            vec![diploid_call(0, 2, 20.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        );
    }

    /// The two "never true" emptiness answers are the ones their documentation promises.
    ///
    /// Both types hold at least one entry by construction — the reference allele, and the
    /// copy count beside it — so these exist to say so at a call site rather than to be a
    /// question. Nothing called either before.
    #[test]
    fn the_allele_table_and_its_copies_are_never_empty() {
        let candidates = generic_reference();
        assert!(!candidates.is_empty());
        assert_eq!(candidates.len(), 1);

        let copies = ExpectedAlleleCopies::new(vec![2.0], &candidates);
        assert!(!copies.is_empty());
        assert_eq!(copies.len(), 1);
    }
}
