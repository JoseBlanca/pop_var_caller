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
//! [`SampleGenotypeCall`], what a locus produces. Of the four sub-modules — the candidate
//! step, the likelihood, the genotype prior and the inference loop — two are here:
//! [`genotype_prior`], step 8, how likely each genotype is before any read is looked at
//! (`doc/devel/ng/impl_plan/calling_prior.md`), and [`likelihood`], step 7, how probable
//! this sample's reads are given each genotype
//! (`doc/devel/ng/impl_plan/calling_read_likelihoods.md`). The other two, and the shared
//! types that borrow from them, arrive with their own plans
//! (`doc/devel/ng/impl_plan/calling_foundations.md`).
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

pub mod genotype_prior;
pub mod genotype_table;
pub mod likelihood;

pub use genotype_table::{GenotypeIdx, GenotypeTable, GenotypeTableView};
/// Re-exported for a different reason from [`genotype_table`]'s three, and the reason is
/// worth stating: **without a `use` of it that has to compile, deleting `pub mod
/// likelihood;` orphans the whole module silently.** The crate still builds, clippy is
/// still clean, and `cargo test --lib ng::calling::likelihood` reports `0 passed; ok` —
/// a green run naming a module that is no longer compiled. With this line the same
/// deletion is `error[E0432]: unresolved import`.
pub use likelihood::generic::{
    ERROR_SPREAD_BASES, LOG_ERROR_SPREAD, LogErrorSpreadTable, NO_LOG_ERROR_SPREAD,
    fill_log_error_spreads, genotype_log_likelihood_row,
};
pub use likelihood::{
    ContaminationMixture, ContaminationView, GenericEvidenceBuffer, GenericObservation,
    GenericSampleEvidence, MAX_BASE_ERROR, MIN_BASE_ERROR, MIN_CONTAMINANT_FREQUENCY,
    ReadGroupCalibration, SsrRowScratch, SsrSampleEvidence, fill_contaminant_allele_frequencies,
};

use crate::ng::locus_generation::LocusKind;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::types::{AlleleId, GenomeRegion, Genotype, Phred};

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

/// One sample's call at one locus: the alleles it carries, and how sure the caller is.
///
/// [`Self::genotype_quality`] is the posterior-derived `GQ` of the loop's last pass —
/// how much of the sample's genotype probability the winning genotype took. Step 13's
/// quality model **refines** this number; it does not replace it
/// (`doc/devel/ng/arch/calling_em_loop.md` §2).
#[derive(Clone, PartialEq, Debug)]
pub struct SampleGenotypeCall {
    /// Which alleles this sample carries, one per copy of its genome.
    pub genotype: Genotype,
    /// How sure the caller is of that genotype.
    pub genotype_quality: Phred,
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
    use crate::ng::locus_generation::SsrDetail;
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
        SampleGenotypeCall {
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
            inference.per_sample[0].genotype.alleles(),
            [AlleleId(0), AlleleId(0)]
        );
        assert_eq!(
            inference.per_sample[1].genotype.alleles(),
            [AlleleId(0), AlleleId(1)]
        );
        assert_eq!(inference.per_sample[1].genotype_quality.get(), 25.0);
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
            capped.per_sample[0].genotype.alleles(),
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
}
