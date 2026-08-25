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

pub mod allele_candidates;
pub mod genotype_prior;
pub mod genotype_table;
pub mod inference;
pub mod likelihood;

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
use crate::ng::locus_generation::{LocusKind, SsrDetail};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::types::{AlleleId, GenomeRegion, Genotype, InbreedingF, LogProb, Phred, Ploidy};

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
    ) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the repeat-tract evidence at {region} names no sample: a locus carries one \
             entry per sample of the run and a run has at least one sample, so an empty \
             list is evidence that went missing rather than a locus nobody covered"
        );
        Self::Ssr {
            region,
            per_sample,
            detail,
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
    /// Two things are checked, and each is a caller bug whose symptom is a wrong genotype
    /// rather than a crash:
    ///
    /// - **the path.** A repeat tract's evidence handed to a locus whose alleles say
    ///   SNP/indel would be scored by the wrong read model, which is a different likelihood
    ///   at every sample with nothing failing
    ///   (`doc/devel/ng/arch/calling_em_loop.md` §2). A repeat **bundle** goes to the
    ///   repeat path with a tract, which is how every other consumer of
    ///   [`LocusKind`] groups the two.
    /// - **the cohort.** One run-wide sample order indexes this evidence, every per-sample
    ///   slice of `parameters`, and the calls that come back. Two lists of different lengths
    ///   are two different orders, and a positional join between them silently pairs one
    ///   sample's reads with another's inbreeding coefficient.
    ///
    /// **The locus's genotype table is not checked here**, because this method never sees
    /// it. That check has its own home, at the one point the shape is fixed:
    /// [`CallingScratch::prepare_for_locus`].
    ///
    /// # Panics
    ///
    /// On either disagreement, in release as well as debug: both are comparisons of two
    /// integers or two discriminants, against a defect that would otherwise reach the
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
    inbreeding_coefficient_by_sample: &'a [InbreedingF],
    prior_seed: SpectrumSeed,
    ssr_slippage_fits: &'a StratumFits,
    ploidy: Ploidy,
}

impl<'a> FrozenParameters<'a> {
    /// Gather the run's frozen parameters for a run the fit **did** measure contamination
    /// in, with what can be checked between them checked.
    ///
    /// - `calibration_by_read_group` — one entry per read group, indexed by
    ///   [`ReadGroupId`](crate::ng::types::ReadGroupId).
    /// - `contamination_by_read_group` — one entry per read group as well.
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
    /// If either per-axis list is empty. A run has at least one sample and at least one
    /// read group, so an empty list is that axis going missing — and the symptom would
    /// otherwise be an out-of-range index at whichever locus first carried a read, naming
    /// neither the axis nor the locus.
    #[must_use]
    pub fn new(
        calibration_by_read_group: &'a [ReadGroupCalibration],
        contamination_by_read_group: &'a [ContaminationView],
        inbreeding_coefficient_by_sample: &'a [InbreedingF],
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: &'a StratumFits,
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
        Self::gather(
            calibration_by_read_group,
            contamination_by_read_group,
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
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
        ploidy: Ploidy,
    ) -> Self {
        Self::gather(
            calibration_by_read_group,
            &[],
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ploidy,
        )
    }

    /// The checks both constructors share, and the only place the fields are written.
    fn gather(
        calibration_by_read_group: &'a [ReadGroupCalibration],
        contamination_by_read_group: &'a [ContaminationView],
        inbreeding_coefficient_by_sample: &'a [InbreedingF],
        prior_seed: SpectrumSeed,
        ssr_slippage_fits: &'a StratumFits,
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
        Self {
            calibration_by_read_group,
            contamination_by_read_group,
            inbreeding_coefficient_by_sample,
            prior_seed,
            ssr_slippage_fits,
            ploidy,
        }
    }

    /// One calibration per read group, indexed by
    /// [`ReadGroupId`](crate::ng::types::ReadGroupId).
    #[inline]
    #[must_use]
    pub fn calibration_by_read_group(&self) -> &'a [ReadGroupCalibration] {
        self.calibration_by_read_group
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

    /// The fitted slippage numbers, looked up by read group and stratum.
    #[inline]
    #[must_use]
    pub fn ssr_slippage_fits(&self) -> &'a StratumFits {
        self.ssr_slippage_fits
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
    /// against one allele. Built once per set of slippage numbers and read by every pass.
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
    /// The SNP/indel row's own scratch.
    generic_row: GenericRowScratch,
    /// The repeat-tract row's own scratch, including the emission model's.
    ssr_row: SsrRowScratch<SsrEmissionScratch>,
    /// How many samples the buffers above are sized for. Zero means never prepared.
    sample_count: usize,
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
    /// If `sample_count` is zero. A cohort has at least one sample, so a locus prepared for
    /// none is a run whose sample order went missing rather than a locus nobody covered.
    ///
    /// If `samples × genotypes` or `samples × alleles` overflows a `usize`.
    pub fn prepare_for_locus(
        &mut self,
        sample_count: usize,
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
            sample_count > 0,
            "a cohort has at least one sample, so a locus prepared for none is a run whose \
             sample order went missing"
        );
        let genotype_count = genotypes.genotype_count();
        let allele_count = genotypes.allele_count();
        let table_len = sample_count.checked_mul(genotype_count).unwrap_or_else(|| {
            panic!(
                "a locus of {sample_count} samples over {genotype_count} genotypes needs a \
                 genotype-likelihood table longer than a usize can index"
            )
        });
        let copies_len = sample_count.checked_mul(allele_count).unwrap_or_else(|| {
            panic!(
                "a locus of {sample_count} samples over {allele_count} alleles needs a \
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

        self.sample_count = sample_count;
        self.genotype_count = genotype_count;
        self.allele_count = allele_count;
    }

    /// How many samples the buffers are currently sized for.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.sample_count
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

    /// The SNP/indel row's own scratch, which owns its own sizing.
    #[inline]
    pub fn generic_row_mut(&mut self) -> &mut GenericRowScratch {
        &mut self.generic_row
    }

    /// The repeat-tract row's own scratch, which owns its own sizing.
    #[inline]
    pub fn ssr_row_mut(&mut self) -> &mut SsrRowScratch<SsrEmissionScratch> {
        &mut self.ssr_row
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
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "its only caller is `score_one_sample`, which step D1 of the \
                      calling-loop plan is the first to call from outside a test"
        )
    )]
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
            self.sample_count > 0,
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
            sample < self.sample_count,
            "sample {sample} is past the {} this scratch was prepared for, indexing the \
             {table} table",
            self.sample_count
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
    /// genotype, built once per set of slippage numbers and read by every pass.
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
    pub fn genotype_quality(&self) -> Option<Phred> {
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
/// Plain data, and six of its eight fields are public because every consumer reads
/// them. The two that are not are [`Self::alleles`] and [`Self::cohort_expected_copies`]:
/// they are one entry per allele of each other, and a public allele table would let a
/// consumer widen it with [`CandidateAlleles::admit`] against unchanged copies — breaking,
/// after construction, the pairing [`Self::new`] checks at it. They are read through
/// accessors instead.
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
    /// **Which warrant is weakest is not yet decided.** [`Provenance`] defines no
    /// ordering — `FittedHere`, `Borrowed`, `Defaulted` and `Supplied` are four names,
    /// not a scale — and in particular where a value the run *supplied* sits against one
    /// the caller *fitted* is an open question. The step that first has to compare two
    /// of them is the one that must settle it.
    pub weakest_provenance: Provenance,
    /// Set when the STR prior's seed could not be scaled to reproduce the measured gene
    /// diversity at this locus — the geometry has a ceiling and the measurement was
    /// above it.
    ///
    /// The loop uses the ceiling and marks the locus rather than silently rescaling, so
    /// that a call resting on a bent prior is distinguishable from one that is not.
    /// **The ceiling is provisional**, not a settled rule: it stands until the STR
    /// prior's own open question about what to do instead is answered
    /// (`doc/devel/ng/arch/calling_priors.md` §5).
    ///
    /// Never set on the SNP/indel path, which seeds from a different quantity — and
    /// [`Self::new`] refuses a locus that sets it there.
    pub seed_diversity_unreachable: bool,
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
    /// If `seed_diversity_unreachable` is set on a SNP/indel locus. That marker belongs
    /// to the STR prior's seed; a generic locus seeds from a different quantity
    /// (`doc/devel/ng/arch/calling_priors.md` §5) and can never raise it, so setting it
    /// there means the marker was wired to the wrong path.
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
    /// **What is deliberately *not* checked: `converged` against `passes`.** A locus that
    /// hit the cap should report the cap, but the cap is run configuration this type does
    /// not see. Nor is `converged` required to imply more than one pass: the loop's first
    /// pass compares against a reads-only estimate made before it
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3), so settling on the first pass is a
    /// real outcome and not a comparison against nothing.
    #[allow(
        clippy::too_many_arguments,
        reason = "the architecture fixes this type's eight fields as a flat list \
                  (arch/calling_em_loop.md §2); grouping them here to satisfy the lint \
                  would be a design change, not a refactor"
    )]
    pub fn new(
        region: GenomeRegion,
        alleles: CandidateAlleles,
        per_sample: Vec<SampleGenotypeCall>,
        cohort_expected_copies: ExpectedAlleleCopies,
        converged: bool,
        passes: u32,
        weakest_provenance: Provenance,
        seed_diversity_unreachable: bool,
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
        assert!(
            !(seed_diversity_unreachable && matches!(alleles.kind(), LocusKind::Generic)),
            "the gene-diversity seed marker belongs to the STR prior: a SNP/indel locus \
             seeds from a different quantity and can never raise it"
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
        Self {
            region,
            alleles,
            per_sample,
            cohort_expected_copies,
            converged,
            passes,
            weakest_provenance,
            seed_diversity_unreachable,
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
    fn diploid_call(first: u16, second: u16, quality: f32) -> SampleGenotypeCall {
        SampleGenotypeCall::Called {
            genotype: Genotype::new(vec![AlleleId(first), AlleleId(second)]),
            genotype_quality: Phred::try_new(quality).expect("a legal quality"),
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
            false,
        );

        assert_eq!(inference.region, region());
        assert_eq!(inference.alleles().len(), 2);
        assert_eq!(inference.cohort_expected_copies().copies(), [2.6, 1.4]);
        assert_eq!(inference.passes, 4);
        assert!(inference.converged);
        assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
        assert!(!inference.seed_diversity_unreachable);

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
                .genotype_quality()
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
    /// Built on the **repeat** path, because the seed marker this also carries cannot
    /// arise on the SNP/indel one — a fixture that set it there would be pinning a state
    /// its own field documents as impossible.
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
            true,
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
        assert!(capped.seed_diversity_unreachable);
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
            false,
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
                false,
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
            false,
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
            false,
        );
    }

    /// The gene-diversity seed marker belongs to the repeat prior. A SNP/indel locus
    /// seeds from a different quantity and can never raise it, so a locus that sets it
    /// there has had the marker wired onto the wrong path — which is exactly what an
    /// implementation slip in the seed's routing would look like.
    #[test]
    #[should_panic(expected = "belongs to the STR prior")]
    fn a_snp_locus_cannot_carry_the_repeat_seed_marker() {
        let (alleles, copies) = two_allele_locus();
        let _ = LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 30.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            true,
        );
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
            false,
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
            false,
        );
        assert_eq!(single.per_sample.len(), 1);
        assert_eq!(single.cohort_expected_copies().copies(), [2.6, 1.4]);
    }

    // ────────────────────────────────────────────────────────────────────────────────
    // The three types the calling loop takes and gives back, and the missing genotype.
    // ────────────────────────────────────────────────────────────────────────────────

    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
    use std::collections::BTreeMap;

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
    fn frozen_parameters<'a>(
        calibration: &'a [ReadGroupCalibration],
        contamination: &'a [ContaminationView],
        inbreeding: &'a [InbreedingF],
        strata: &'a StratumFits,
    ) -> FrozenParameters<'a> {
        if contamination.is_empty() {
            FrozenParameters::uncontaminated(
                calibration,
                inbreeding,
                neutral_seed(),
                strata,
                diploid_ploidy(),
            )
        } else {
            FrozenParameters::new(
                calibration,
                contamination,
                inbreeding,
                neutral_seed(),
                strata,
                diploid_ploidy(),
            )
        }
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
        let _ = frozen_parameters(&calibration, &contamination, &inbreeding, &strata);
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
        let _ = FrozenParameters::new(
            &calibration,
            &[],
            &inbreeding,
            neutral_seed(),
            &strata,
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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

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
        let _ = frozen_parameters(&[], &[], &inbreeding, &strata);
    }

    /// A run whose sample order went missing is refused, rather than producing a locus with
    /// no calls in it.
    #[test]
    #[should_panic(expected = "an empty list is a run whose sample order went missing")]
    fn frozen_parameters_refuse_a_run_with_no_samples() {
        let calibration = one_read_group();
        let strata = no_strata();
        let _ = frozen_parameters(&calibration, &[], &[], &strata);
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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

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
        let tract = LocusEvidence::ssr(region(), &tract_samples, &detail);
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
        assert_eq!(scratch.sample_count(), 3);
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

    /// A locus prepared for no samples at all is a run whose sample order went missing.
    #[test]
    #[should_panic(expected = "a locus prepared for none")]
    fn a_scratch_prepared_for_no_samples_is_refused() {
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
        assert_eq!(missing.genotype_quality(), None);

        let called = diploid_call(0, 1, 25.0);
        assert!(!called.is_missing());
        assert_eq!(
            called.genotype().expect("a called sample").alleles(),
            [AlleleId(0), AlleleId(1)]
        );
        assert_eq!(
            called.genotype_quality().expect("a called sample").get(),
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
            false,
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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

        let bundle = CandidateAlleles::new(Box::from(b"CAGCAG".as_slice()), LocusKind::SsrBundle);
        let detail = ssr_detail();
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
        let evidence = LocusEvidence::ssr(region(), &per_sample, &detail);
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
        let parameters = frozen_parameters(&calibration, &[], &inbreeding, &strata);

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
        let _ = LocusEvidence::ssr(region(), &[], &detail);
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
        assert_eq!(scratch.sample_count(), 2);

        // The three sub-scratches are reachable, and are the worker's rather than a locus's.
        let _ = scratch.candidate_selection_mut();
        let _ = scratch.generic_row_mut();
        let _ = scratch.ssr_row_mut();
    }

    /// A repeat tract sets no sample aside, so a missing call at one is the SNP/indel path's
    /// ruling wired onto the wrong path.
    ///
    /// **The mirror of the gene-diversity marker's check**, which refuses the repeat path's
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
            false,
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
            false,
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
