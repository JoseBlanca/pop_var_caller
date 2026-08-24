//! Step 6 — **choosing the short list of sequences a locus is called over.**
//!
//! The cohort merge collects every sequence any sample's reads showed over a stretch of
//! genome and unifies them into one table; it narrows nothing
//! (`doc/devel/ng/spec/cohort_merge.md` §4.2). Everything downstream is defined over the
//! narrowed list — the read likelihood scores each observation against each candidate,
//! the genotype prior lays its mass over the genotypes those candidates make, and the
//! VCF's `ALT` column is what survived. This module does the narrowing
//! (`doc/devel/ng/spec/candidate_alleles.md`, and `doc/devel/ng/arch/candidate_alleles.md`
//! for the shapes).
//!
//! **A folder rather than a file, and no trait.** Two paths — the ordinary SNP/indel one
//! in `generic.rs` and the repeat tract's — take different evidence and return different
//! extras, and which runs is decided by the locus's kind rather than by a swappable
//! recipe. Two functions, not two impls of one seam. What both share lives here: the
//! config, the verdict, the leftover, the remapping and the ranking.
//!
//! **The rule, in one line.** An alternative survives if *some single sample's* reads
//! lent it at least `max(2 reads, 5 in 100 of that sample's reads at the locus)`, and a
//! locus is called over at most six alleles counting the reference. **No term of the bar
//! reads the cohort** — one sample reaching it admits the sequence for everyone —
//! because otherwise a sample's candidate list would depend on who else is in the run
//! (spec §3.2).

use crate::ng::calling::CandidateAlleles;
use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare, MinAltReads};
use crate::ng::types::AlleleId;

/// **The support one sample must lend one sequence for it to be called over, and the cap
/// on how many sequences a locus is called over at all** — the two halves of the
/// narrowing, and the run's only knobs on it (spec §3, §4).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct CandidateSelectionConfig {
    /// `max(floor, ceil(share × that sample's reads at the locus))` non-reference reads,
    /// asked of **each sample separately against that sample's own reads**; one sample
    /// reaching it admits the sequence for the whole cohort (spec §3.2).
    ///
    /// **The merge's own type, reused rather than copied**, so the rule that decides
    /// whether a locus is built and the rule that decides which of its alleles are
    /// called over cannot drift apart — a sweep of one is a sweep of both. Only the
    /// number differs; see [`DEFAULT_MIN_ALLELE_SUPPORT`].
    pub min_allele_support: MinAltReads,
    /// How many alleles a locus may be called over, **counting the reference**. Above it
    /// the list is cut to the best-ranked and the locus is still called; it is never
    /// refused (spec §4.1).
    pub max_candidate_alleles: MaxCandidateAlleles,
}

impl CandidateSelectionConfig {
    /// [`DEFAULT_MIN_ALLELE_SUPPORT`] and [`DEFAULT_MAX_CANDIDATE_ALLELES`] — both soft,
    /// and both carrying their source in their own documentation.
    pub const DEFAULT: Self = Self {
        min_allele_support: DEFAULT_MIN_ALLELE_SUPPORT,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES,
    };
}

impl Default for CandidateSelectionConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// **How many alleles a locus may be called over, counting the reference** — at least
/// two, because a cap below that could not admit a single alternative.
///
/// **A cap of 0 or 1 is refusal under another name**, and spec §4.1 rules refusal out: a
/// locus carrying two obvious variants and six noise sequences must keep the two, not lose
/// them with the six. At either value the reference is the only survivor and every
/// alternative becomes a truncation, so the type refuses them rather than leaving a
/// later step to discover it.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MaxCandidateAlleles(u16);

impl MaxCandidateAlleles {
    /// The default cap, [`DEFAULT_MAX_CANDIDATE_ALLELES`].
    pub const DEFAULT: Self = DEFAULT_MAX_CANDIDATE_ALLELES;

    /// The smallest cap that is a cap rather than a refusal: the reference and one
    /// alternative.
    pub const SMALLEST: u16 = 2;

    /// The cap, or `None` below [`SMALLEST`](Self::SMALLEST). Refusing rather than
    /// clamping, for the reason [`MinAltReadShare::new`] gives about its own range: a
    /// mistyped cap is a run whose output looks ordinary.
    pub fn new(alleles: u16) -> Option<Self> {
        (alleles >= Self::SMALLEST).then_some(Self(alleles))
    }

    /// The same cap, for a `const` that has to name one — and it **panics** where
    /// [`new`](Self::new) returns `None`, for the reason
    /// [`MinAltReadShare::new_or_panic`] gives.
    pub const fn new_or_panic(alleles: u16) -> Self {
        assert!(
            alleles >= Self::SMALLEST,
            "a locus is called over at least the reference and one alternative"
        );
        Self(alleles)
    }

    /// The cap, counting the reference.
    #[inline]
    pub const fn get(self) -> u16 {
        self.0
    }

    /// How many *alternatives* fit under it — one fewer, and never negative because the
    /// constructors refuse a cap below two.
    #[inline]
    pub const fn alternatives(self) -> u16 {
        self.0 - 1
    }
}

impl Default for MaxCandidateAlleles {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// **Two reads, or 5 in 100 of that sample's reads at the locus, whichever is more.**
///
/// **The floor is the merge's own** ([`MinAltObs::DEFAULT`], production's number) **and
/// should stay at 2**: measured against the GIAB trio's v4.2.1 truth set over 572 kb on
/// 2026-08-24, **at 30×** raising it from 2 to 3 loses five true alternative alleles,
/// where raising the share to 10 in 100 loses two for the same reduction in table size —
/// 1,539 alternatives kept against 1,601. The floor is the expensive knob (spec §3.3).
///
/// **The share is 5 in 100 where the merge's keep rule uses 2**, because the
/// allele-level question tolerates a stricter share at depth. On the same trio **at
/// 300×** it cuts the merge's 15,474 alternatives to 2,308, where a bar of 2 reads alone
/// keeps 10,793 — and it loses the same two true alleles the 2-in-100 share loses.
///
/// **It is inert below 41 compared reads a sample**, which is the arithmetic rather than
/// a measurement: `ceil(0.05 × 40) = 2` is the floor, and 41 is the first count at which
/// the share asks for more. So a tomato-depth run — about 11 compared reads a sample at a
/// locus — sees the identical rule it would have seen at 2 in 100. What was measured on
/// that panel is the neighbouring comparison: turning the share off entirely against the
/// merge's 2 in 100 moves 4 loci in 53,935 (spec §3.3).
///
/// **Soft.** Measured on one human trio over 572 kb (spec §11, Q3); what would move it is
/// the same scoring on a second high-depth cohort.
pub const DEFAULT_MIN_ALLELE_SUPPORT: MinAltReads = MinAltReads {
    floor: MinAltObs::DEFAULT,
    share: MinAltReadShare::new_or_panic(0.05),
};

/// **Six alleles including the reference** — production's `DEFAULT_MAX_ALLELES_PER_RECORD`
/// (`src/var_calling/per_group_merger.rs`), inherited and declared inherited.
///
/// **It is not GATK's number, though production's comment says it is.** GATK's
/// `--max-alternate-alleles` defaults to 6 *alternates*
/// (`DEFAULT_MAX_ALTERNATE_ALLELES = 6`, documented "Maximum number of alternate alleles
/// to genotype"), so GATK genotypes over seven alleles where ng genotypes over six.
/// Production's constant counts a record's whole allele set — `enforce_max_alleles`
/// compares `unified.alleles.len()` and protects the reference from pruning — and its own
/// doc comment nonetheless claims to match GATK's. **ng is the tighter of the two by one
/// allele**, which is 21 genotypes against 28 at diploid.
///
/// **Measured to be a safety valve rather than a working part at the cohort sizes we
/// have** (spec §4.2, 2026-08-24): it binds at 23 of 53,935 tomato loci — one in 2,300 —
/// and at none of the GIAB trio's 4,177 loci at 30× or 7,478 at 300×. What it guards
/// against grows with the cohort, which is why it is here: holding the tomato allele
/// table fixed and asking the bar of 1, 4, 16 and 63 samples gives 0, 0, 3 and 23 loci
/// above six alleles.
///
/// **Those counts were taken with the merge's 2-in-100 share, not the 5 in 100 this
/// module ships**, which spec §4.2 states in its own header and a reader of this constant
/// would otherwise not know. The direction is safe: at tomato depth the two shares are
/// provably the same rule (see [`DEFAULT_MIN_ALLELE_SUPPORT`]), and everywhere else 5 in 100
/// admits fewer alternatives, so the cap binds no more often than these numbers say.
///
/// **Soft, and never measured at its own value.** Whether it becomes load-bearing past a
/// few hundred samples is an extrapolation from that table, not a measurement
/// (spec §11, Q2).
pub const DEFAULT_MAX_CANDIDATE_ALLELES: MaxCandidateAlleles = MaxCandidateAlleles::new_or_panic(6);

/// **What selection did at one locus, beyond the list itself** — whether the list is
/// everything that cleared the bar, or was cut to fit the cap.
///
/// **There is no depth variant, and its absence is a decision** (spec §6.2). The
/// architecture once sketched `Admission { Ok, LowDepth, NotPeriodic, TooManyAlleles }`
/// and it was never built. `TooManyAlleles` named a refusal, and spec §4.1 chose
/// truncation instead. `LowDepth` would re-ask the merge's keep rule with a different
/// denominator — production's version of that is a sum over the cohort, measured refusing
/// 98.6% of repeat tracts at one sample sequenced to 5× against 0.2% at 300×. **Depth is
/// asked once, upstream, per sample.**
///
/// `#[non_exhaustive]` because the repeat-tract path adds the third variant and a fourth
/// is not ruled out; a match on this in another crate must keep a wildcard arm.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SelectionVerdict {
    /// Everything that cleared the bar is in the list.
    ///
    /// **Including when the list is the reference alone**, which is a first-class outcome
    /// and not an error: the merge builds a locus when some sample's non-reference reads
    /// *pooled* reach its rule, and two reads split one and one across two alternatives
    /// clear that while clearing neither allele bar. Measured at more than one built locus
    /// in four, and the fraction is the same on both benchmarks — 27.4% on the 63-accession
    /// tomato panel, 27.3% on the GIAB trio at 30× and 28.0% at 300× (spec §6.2).
    Selected,
    /// The cap bound, and `dropped` alternatives were cut — the lowest-ranked first, by
    /// the ranking `ranks_above` will define at step B2: the largest share of one sample's
    /// compared reads, then how many samples cleared the bar, then the cohort's read
    /// total, then the bases. The reference is never among them (spec §4.1).
    Truncated {
        /// How many alternatives the cap removed. Not how many were dropped in total: an
        /// alternative that failed the bar was never a candidate for the cap.
        dropped: u16,
    },
    /// **Repeat tracts only, and never minted by the ordinary path.** A tract whose reads
    /// do not agree on a repeat unit is not called as a tract
    /// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7); `generic.rs` has no way to reach
    /// this and no reason to.
    NotPeriodic,
}

/// **One sample's reads whose sequence selection dropped, and the error mass they carry** —
/// the pool the SNP/indel read likelihood scores reads against no candidate with
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.3's `q_sum_other`). **Nothing upstream
/// produces it, because nothing upstream drops anything**, so selection owes it.
///
/// **The count is not decoration: it is what makes truncation defensible.** The mass is
/// the same under every genotype and cancels in genotyping, so without the count a sample
/// whose true allele was cut is scored confidently against a set that does not contain it,
/// with nothing per-sample saying so. The count is what lets a later step no-call that
/// sample. It costs one integer per sample per locus. **Truncation and this count are one
/// decision, not two: drop the count and refusing the locus becomes the correct policy
/// again** (spec §5).
///
/// **What is not in here** (spec §5.1): partial reads, which say the sample carries *at
/// least* this rather than what it carries and are scored on their own axis; reads that
/// produced no observation and reads removed as evidence, neither of which carries a
/// quality sum and neither of which was ever in the allele table; and the reference's own
/// reads, since the reference is never dropped.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct UnmatchedSupport {
    /// How many of this sample's reads showed a sequence selection dropped.
    pub num_reads: u32,
    /// Σ `ln P(error)` over those reads — **summed straight from the merge's own per-row
    /// `q_sum`**, never re-derived from a count and a rate. Zero, not negative, where
    /// nothing was dropped.
    pub q_sum: f64,
    /// **Of those reads, the ones on a sequence *this sample itself* showed convincingly —
    /// it cleared the support rule for this sample — which the cap then cut.** Non-zero
    /// means this sample cannot be genotyped here; see
    /// [`genotype_must_be_missing`](Self::genotype_must_be_missing).
    ///
    /// **Why this is separate from [`num_reads`](Self::num_reads), and it is the whole
    /// point of the field.** Sequences are dropped two ways, and only one of them says
    /// anything about a sample. The support rule drops sequences almost nobody showed,
    /// which are overwhelmingly sequencing error — on the GIAB trio at 300× that is
    /// 13,166 of the merge's 15,474 alternatives (spec §3.3). Every sample has a few
    /// error reads at nearly every locus, so a rule keyed on `num_reads` would emit a
    /// missing genotype almost everywhere. The cap is the other way, and it only ever
    /// cuts sequences that already cleared the rule for *somebody* — but not necessarily
    /// for the sample whose reads are being counted here, which is why "cleared it for
    /// this sample" is the condition and not "the cap cut it".
    pub earned_reads_cut_by_the_cap: u32,
}

impl UnmatchedSupport {
    /// **Whether this sample's genotype must be emitted as missing at this locus.**
    ///
    /// True when the cap removed a sequence this sample's own reads had earned. The
    /// sample carries something the locus is no longer called over, so every genotype the
    /// caller can form for it is wrong — and the read likelihood cannot say so on its
    /// own, because the pooled error mass is identical under every genotype and cancels
    /// (`doc/devel/ng/spec/read_likelihoods.md` §3.3). **Emitting a made-up genotype
    /// instead is the behaviour this rule exists to prevent** (owner's decision,
    /// 2026-08-24; spec §5).
    ///
    /// **This is what makes truncation defensible at all.** Without it the honest policy
    /// would be to refuse the whole locus, which is what HipSTR does above 1,000
    /// haplotypes and what the existing repeat-tract caller does above 24 candidates —
    /// and refusing costs the other samples a locus they were called at perfectly well.
    #[inline]
    pub fn genotype_must_be_missing(&self) -> bool {
        self.earned_reads_cut_by_the_cap > 0
    }
}

/// **For each allele of the merge's table, in that table's own index order: the id it now
/// has among the candidates, or nothing where selection dropped it.**
///
/// **The evidence builder cannot be written without this.** [`CandidateAlleles`] ids are
/// dense and in admission order, while a `SupportedAllele::allele` indexes the *merge's*
/// table; after narrowing the two are different numbers and nothing else records the
/// correspondence (arch §2.3).
/// **The ids it hands out are dense, ascending and gapless**, because they are
/// [`CandidateAlleles`]' own — that table's `admit` returns its previous length, so the
/// *n*th allele admitted is always `AlleleId(n)`. [`admit`](Self::admit) asserts it rather
/// than trusting it, which is what stops two merge alleles being recorded onto one id: a
/// collision there is not an out-of-range id and no bounds check would see it, and the
/// evidence hand-off would then re-key two different sequences' reads onto one candidate
/// and the read likelihood would score them as one.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct AlleleRemap {
    /// One entry per allele of the merge's table, in that table's index order.
    to_candidate: Box<[Option<AlleleId>]>,
    /// How many have been admitted so far — the id the next admission must carry, and
    /// what makes [`num_admitted`](Self::num_admitted) a read rather than a scan.
    ///
    /// **A `u32` where an [`AlleleId`] is a `u16`**, so that the count can exceed what an
    /// id can name without wrapping. It cannot get there in a run — the cap is six — but
    /// at exactly 65,536 admissions a `u16` counter would wrap in release and panic in
    /// debug, and a check that behaves differently in the two profiles is worse than the
    /// case it guards.
    num_admitted: u32,
}

impl AlleleRemap {
    /// A remapping over a merge table of `table_len` alleles, with every allele dropped.
    /// Selection fills the survivors in as it admits them.
    #[inline]
    pub fn with_all_dropped(table_len: usize) -> Self {
        Self {
            to_candidate: vec![None; table_len].into_boxed_slice(),
            num_admitted: 0,
        }
    }

    /// The candidate id for one of the merge table's alleles, or `None` where it was
    /// dropped.
    ///
    /// **`None` means absent and is never a sentinel id**, which is the whole reason this
    /// is an `Option` rather than a reserved number: a sentinel would be a legal
    /// [`AlleleId`] and would index a real allele.
    ///
    /// # Panics
    ///
    /// On a `table_index` outside the merge table this was built for. That is a caller
    /// bug — a support row naming an allele the table does not hold — and spec §8 makes
    /// it an assertion held in release rather than a value.
    #[inline]
    pub fn candidate_for(&self, table_index: usize) -> Option<AlleleId> {
        assert!(
            table_index < self.to_candidate.len(),
            "a support row named allele {table_index} of a merge table holding {}",
            self.to_candidate.len()
        );
        self.to_candidate[table_index]
    }

    /// Record that the merge table's `table_index` survived as `candidate`.
    ///
    /// # Panics
    ///
    /// Three caller bugs, all of them silent if they were let through:
    ///
    /// - an out-of-range `table_index`;
    /// - an allele admitted twice, which would leave the first id naming nothing;
    /// - a `candidate` that is not the next dense id. **This is the one a bounds check
    ///   cannot see.** Two merge alleles recorded onto one [`AlleleId`] are both in range
    ///   and both written once, and the evidence hand-off would then re-key two different
    ///   sequences' reads onto one candidate — two alleles scored as one, with a genotype
    ///   that looks ordinary coming out. Since every id comes from
    ///   [`CandidateAlleles::admit`], which returns the table's previous length, the next
    ///   dense id is what a correct caller always passes.
    #[inline]
    pub fn admit(&mut self, table_index: usize, candidate: AlleleId) {
        assert!(
            table_index < self.to_candidate.len(),
            "cannot admit allele {table_index} of a merge table holding {}",
            self.to_candidate.len()
        );
        assert!(
            self.to_candidate[table_index].is_none(),
            "allele {table_index} of the merge table was admitted twice, first as {:?} and now as {candidate:?}",
            self.to_candidate[table_index]
        );
        assert_eq!(
            u32::from(candidate.get()),
            self.num_admitted,
            "candidate ids are dense and in admission order, so allele {table_index} of the \
             merge table must be admitted as AlleleId({}) and not {candidate:?}",
            self.num_admitted
        );
        self.to_candidate[table_index] = Some(candidate);
        self.num_admitted += 1;
    }

    /// How many alleles the merge table held — the length this remapping is indexed over.
    #[inline]
    pub fn table_len(&self) -> usize {
        self.to_candidate.len()
    }

    /// How many of them survived — which is also how long the candidate table is, since
    /// the two are filled together.
    #[inline]
    pub fn num_admitted(&self) -> usize {
        self.num_admitted as usize
    }
}

/// **Everything selection produces at one locus**: the narrowed table, what it did, the
/// per-sample leftover, and the correspondence between the merge's allele indices and the
/// new ids.
/// **The fields are private and [`new`](Self::new) is the only way to build one**, so the
/// two invariants below cannot be written around. That is a departure from arch §2.4's
/// `pub` fields, taken at Checkpoint A: a bypassable check on a value whose defect is a
/// wrong genotype rather than a crash is not a check.
#[derive(Clone, PartialEq, Debug)]
pub struct LocusSelection {
    /// The surviving sequences, the reference at [`AlleleId::REFERENCE`].
    alleles: CandidateAlleles,
    /// What selection did beyond the list — everything that cleared the bar, or a list cut
    /// to fit the cap. See [`SelectionVerdict`].
    verdict: SelectionVerdict,
    /// **Parallel to `CohortObservation::per_sample`** — same length, same order, so entry
    /// `i` is the leftover of the sample that `per_sample[i]` describes.
    ///
    /// **Not indexed by the run's sample order**, and the distinction is not cosmetic:
    /// `per_sample` holds only the samples that *covered* the locus, each naming its own
    /// sample index, so a run of 63 accessions can produce a locus whose `per_sample` has
    /// 4 entries. Indexing this by the run's order would put four leftovers at four
    /// scattered positions and 59 zeroed rows between them, and a zeroed row is
    /// indistinguishable from a covering sample that dropped nothing.
    unmatched: Vec<UnmatchedSupport>,
    /// The merge table's indices mapped onto the surviving ids — see [`AlleleRemap`].
    remap: AlleleRemap,
}

impl LocusSelection {
    /// **The one door**, and the only place the two parallelism invariants are checked.
    ///
    /// `covering_samples` is `observation.per_sample.len()` — how many samples had reads
    /// over the locus, which is what `unmatched` must be as long as. It is passed rather
    /// than the observation itself so that this module need not import the merge's
    /// assembled locus type.
    ///
    /// **The fields stay public because arch §2.4 declares them so**, which means this
    /// constructor can be bypassed by a struct literal. It is still worth having: it is
    /// what `select_generic` and `select_ssr` call, so the checks run on every value a run
    /// produces, and a test that writes a literal is a test that has said it wants to.
    /// Making the fields private would make this the *only* door and is raised at
    /// Checkpoint A.
    ///
    /// # Panics
    ///
    /// If `unmatched` is not one entry per covering sample, or if the remapping admitted a
    /// different number of alleles from the table's length. Both are caller bugs whose
    /// symptom is a wrong genotype rather than a crash: the first shifts every sample's
    /// leftover onto its neighbour, and the second means an admitted allele has no bases
    /// or a table entry has no evidence.
    pub fn new(
        alleles: CandidateAlleles,
        verdict: SelectionVerdict,
        unmatched: Vec<UnmatchedSupport>,
        remap: AlleleRemap,
        covering_samples: usize,
    ) -> Self {
        assert_eq!(
            unmatched.len(),
            covering_samples,
            "the leftover runs parallel to the locus's covering samples, one entry each"
        );
        assert_eq!(
            remap.num_admitted(),
            alleles.len(),
            "every admitted allele must be in the candidate table and every candidate must \
             have been admitted"
        );
        Self {
            alleles,
            verdict,
            unmatched,
            remap,
        }
    }

    /// The surviving sequences, the reference at [`AlleleId::REFERENCE`].
    #[inline]
    pub fn alleles(&self) -> &CandidateAlleles {
        &self.alleles
    }

    /// What selection did beyond the list.
    #[inline]
    pub fn verdict(&self) -> SelectionVerdict {
        self.verdict
    }

    /// One leftover per covering sample, in `CohortObservation::per_sample`'s order.
    #[inline]
    pub fn unmatched(&self) -> &[UnmatchedSupport] {
        &self.unmatched
    }

    /// The merge table's indices mapped onto the surviving ids.
    #[inline]
    pub fn remap(&self) -> &AlleleRemap {
        &self.remap
    }

    /// The four parts, by value — what the calling loop's input edge needs, since it takes
    /// ownership of the table and builds its evidence views from the remapping.
    #[inline]
    pub fn into_parts(
        self,
    ) -> (
        CandidateAlleles,
        SelectionVerdict,
        Vec<UnmatchedSupport>,
        AlleleRemap,
    ) {
        (self.alleles, self.verdict, self.unmatched, self.remap)
    }

    /// How many alternatives survived — the number the genotype prior divides its
    /// alternative concentration by (`doc/devel/ng/spec/calling_priors.md` §4, which
    /// spells the same quantity `alternative_allele_count`).
    ///
    /// Zero at a locus that selected down to the reference alone, which is legal and
    /// happens at more than one built locus in four on both benchmarks (see
    /// [`SelectionVerdict::Selected`]). **The subtraction cannot underflow**:
    /// [`CandidateAlleles`] has one constructor, which pushes the reference, and one
    /// mutator, which only pushes, so its length is at least one by construction.
    #[inline]
    pub fn alternative_allele_count(&self) -> usize {
        self.alleles.len() - 1
    }
}

/// **The fold's buffers, one set per worker** — cleared and refilled at every locus, so a
/// locus costs no allocation for the fold itself. What a locus *does* allocate is its
/// output: the surviving table, and the [`AlleleRemap`] over the merge's table (arch §1,
/// and §6's open item on whether that remapping should be a bitset plus a prefix sum).
///
/// **It stands alone here and becomes a field of `CallingScratch` when that type exists**
/// (`doc/devel/ng/impl_plan/calling_loop.md` A1, which builds it). The same worker runs
/// selection and then the calling loop on the same locus, so a second per-worker
/// allocation would buy nothing; arch §2.4 records the reasoning. Nothing about the shape
/// changes when it moves.
///
/// **Not `Clone`.** Two workers sharing one set of buffers by copy would each pay the
/// allocation this type exists to avoid, and a cloned scratch carrying the previous
/// locus's fold is the one state that must never be read.
#[derive(Default, Debug)]
pub struct SelectionScratch {
    /// One entry per allele of the merge's table, in that table's index order.
    per_allele: Vec<AlleleSummary>,
    /// Merge table indices, ordered by the cap's ranking. Only the alternatives that
    /// cleared the bar ever enter it.
    ranked_table_indices: Vec<u32>,
}

impl SelectionScratch {
    /// Empty buffers. Cheap: the first locus sizes them and every later one reuses them.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty the buffers and make room for a table of `table_len` alleles, without
    /// releasing what earlier loci already reserved.
    ///
    /// **`clear` and not `truncate`, and a fresh `resize` rather than a reused row**, so
    /// that no locus can read a value another locus left behind — the failure this whole
    /// scratch shape invites, and the one nothing downstream would notice. Three wrong
    /// versions are caught by `resetting_the_scratch_leaves_no_value_from_an_earlier_locus`
    /// (a `resize` with no `clear`; a `clear` that resizes to the old length; a `truncate`
    /// in place of the `clear`) and a fourth — resetting rows in place and growing without
    /// ever shrinking — only by
    /// `resetting_the_scratch_to_a_smaller_table_shrinks_it`, which is why both are here.
    ///
    /// **The buffers are destructured rather than named through `self`**, so that a field
    /// added later — the repeat-tract path commits to adding one
    /// (`doc/devel/ng/arch/candidate_alleles_ssr.md` §5) — fails to compile here instead
    /// of silently carrying the previous locus's values into the next.
    pub fn reset_for(&mut self, table_len: usize) {
        let Self {
            per_allele,
            ranked_table_indices,
        } = self;
        per_allele.clear();
        per_allele.resize(table_len, AlleleSummary::default());
        ranked_table_indices.clear();
        ranked_table_indices.reserve(table_len);
    }

    /// How many alleles the buffers are currently sized for.
    #[inline]
    pub fn table_len(&self) -> usize {
        self.per_allele.len()
    }
}

/// **One allele's fold across every covering sample** — what the bar and the cap read.
///
/// Private to the module: it is the shape of a computation, not part of what selection
/// hands back. The fold that fills it is step B1 and the ranking that reads it is B2.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
struct AlleleSummary {
    /// **The largest share of one sample's compared reads this allele took**, maximised
    /// over samples — the cap's first ranking key (spec §4.1). Not a cohort share: a
    /// cohort total would truncate the private alleles first at scale, which is the one
    /// thing the ranking exists to avoid.
    best_within_sample_share: f64,
    /// How many samples' reads reached the bar. **One is enough to admit the allele**
    /// (spec §3.2); the count is the cap's tie-break, where it can only reorder and never
    /// exclude.
    samples_clearing_the_bar: u32,
    /// This allele's reads over the whole cohort — the ranking's third key, and the only
    /// place a cohort sum appears in this module.
    cohort_reads: u64,
}

impl AlleleSummary {
    /// **Whether this allele survives the bar** — whether *some single sample's* reads
    /// reached it, which is the admission rule of spec §3 in one line.
    ///
    /// Derived rather than stored. A first draft carried a `reached_the_bar` flag beside
    /// the count, and the two can only ever agree: the fold raises the count in the same
    /// branch that would have set the flag, so a stored flag is a second copy of one fact
    /// and the only thing it can do is disagree with the first.
    ///
    /// **The reference is not asked this question.** It is admitted before any sample's
    /// evidence is read and is exempt from both the bar and the cap (spec §6.1), so C1
    /// seeds it structurally rather than reading a flag that says it passed.
    #[inline]
    #[allow(
        dead_code,
        reason = "its shipping caller is the admission pass of step C1; the test below is \
                  its only caller today, which `expect` cannot express — that would be \
                  unfulfilled in the test build and satisfied in the library one"
    )]
    fn cleared_the_bar(self) -> bool {
        self.samples_clearing_the_bar > 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::LocusKind;

    /// The two numbers themselves, and the floor's **coupling** to the merge's constant
    /// rather than to the digit 2 — the doc comment says the floor *is* the merge's own,
    /// and only the second assertion holds that.
    #[test]
    fn the_default_bar_is_two_reads_or_five_in_a_hundred() {
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.floor.get(), 2);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.floor, MinAltObs::DEFAULT);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.share.get(), 0.05);
    }

    /// The two ends of the committed depth range, as spec §3 states them: at 3 compared
    /// reads the rule asks 2, and at 300 it asks 15.
    ///
    /// The third count is there because the first two cannot see the rounding: `0.05 ×
    /// 300` is exactly 15, so rounding the share *down* would answer 15 as well. At 301
    /// the share is 15.05, and up and down are 16 and 15.
    #[test]
    fn the_floor_decides_at_three_reads_and_the_share_at_three_hundred() {
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(3), 2);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(300), 15);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(301), 16);
    }

    /// The share is stricter than the merge's own at depth and **indistinguishable from
    /// it below 41 compared reads** — the claim [`DEFAULT_MIN_ALLELE_SUPPORT`]'s
    /// documentation makes, held against the merge's constant rather than against a
    /// number retyped here.
    ///
    /// **40 and 41 are the fixture, and stopping short of them is what makes this test
    /// vacuous.** With the equality arm ending at 20 compared reads, a share of 10 in 100
    /// also passes — `ceil(0.10 × 20) = 2` is still the floor — so the test would admit a
    /// bar twice as strict as the one shipped, which spec §3.3 measures losing two more
    /// true alleles at 300×. Carrying the pair to 40 and 41 pins the share to **more than
    /// 2/41 — about 0.0488 — and no more than 0.05**, which is the narrowest window the
    /// rule's own integer arithmetic can express.
    #[test]
    fn the_allele_share_binds_only_above_forty_compared_reads() {
        for compared_reads in [1_u32, 3, 11, 20, 40] {
            assert_eq!(
                DEFAULT_MIN_ALLELE_SUPPORT.required_of(compared_reads),
                MinAltReads::DEFAULT.required_of(compared_reads),
                "at {compared_reads} compared reads the floor decides for both rules"
            );
        }
        assert!(
            DEFAULT_MIN_ALLELE_SUPPORT.required_of(41) > MinAltReads::DEFAULT.required_of(41),
            "41 compared reads is where the allele rule first asks for more than the merge's"
        );
        assert!(
            DEFAULT_MIN_ALLELE_SUPPORT.required_of(300) > MinAltReads::DEFAULT.required_of(300),
            "at 300 compared reads the allele rule must ask for more than the merge's"
        );
    }

    #[test]
    fn the_cap_default_is_six_and_the_config_carries_it() {
        assert_eq!(DEFAULT_MAX_CANDIDATE_ALLELES.get(), 6);
        assert_eq!(
            DEFAULT_MAX_CANDIDATE_ALLELES.alternatives(),
            5,
            "six counting the reference is five alternatives — the number GATK's own \
             default of six names, which is why ng's cap is the tighter of the two"
        );
        assert_eq!(
            CandidateSelectionConfig::default().max_candidate_alleles,
            DEFAULT_MAX_CANDIDATE_ALLELES
        );
    }

    /// **A cap below two is refused rather than clamped**, because at 0 or 1 the reference
    /// is the only survivor and every alternative becomes a truncation — refusal under
    /// another name, which spec §4.1 rules out.
    #[test]
    fn a_cap_that_cannot_hold_one_alternative_is_refused() {
        assert!(MaxCandidateAlleles::new(0).is_none(), "no allele at all");
        assert!(
            MaxCandidateAlleles::new(1).is_none(),
            "the reference and nothing else"
        );
        assert_eq!(
            MaxCandidateAlleles::new(2).map(MaxCandidateAlleles::get),
            Some(2),
            "the reference and one alternative is the smallest cap that is a cap"
        );
        assert_eq!(
            MaxCandidateAlleles::new(2).map(MaxCandidateAlleles::alternatives),
            Some(1)
        );
    }

    #[test]
    #[should_panic(expected = "at least the reference and one alternative")]
    fn a_const_cap_below_two_fails_the_build() {
        let _ = MaxCandidateAlleles::new_or_panic(1);
    }

    /// **The config a run gets by default is the two announced constants, and its support
    /// rule is the allele one rather than the merge's.**
    ///
    /// Both names are in scope in this file and both are a `MinAltReads`, so writing
    /// `MinAltReads::DEFAULT` here is a one-token slip that nothing else would catch: the
    /// type is right, the floor is right, and only the share moves — from 5 in 100 to 2 in
    /// 100. It is invisible at tomato depth, where the two are the same rule, and on the
    /// GIAB trio at 300× it is the difference between keeping 2,308 alternatives and
    /// keeping 5,596 (spec §3.3), each of which the genotype prior divides its
    /// concentration by. The last assertion is what makes the test about the *number* and
    /// not merely about the type.
    #[test]
    fn the_default_config_is_the_two_announced_constants_and_not_the_merges_rule() {
        let config = CandidateSelectionConfig::default();
        assert_eq!(config, CandidateSelectionConfig::DEFAULT);
        assert_eq!(config.min_allele_support, DEFAULT_MIN_ALLELE_SUPPORT);
        assert_eq!(config.max_candidate_alleles, DEFAULT_MAX_CANDIDATE_ALLELES);
        assert_ne!(
            config.min_allele_support,
            MinAltReads::DEFAULT,
            "the allele rule and the merge's keep rule share a type, not a share"
        );
    }

    // ---- the output vocabulary (A2) ----

    /// A remapping over five merge alleles with the middle one dropped — the shape the
    /// evidence hand-off of arch §3.2 walks. Reference, then alternatives 1, 3 and 4
    /// admitted in table order, so their dense ids are 1, 2 and 3.
    fn remap_with_a_hole() -> AlleleRemap {
        let mut remap = AlleleRemap::with_all_dropped(5);
        remap.admit(0, AlleleId::REFERENCE);
        remap.admit(1, AlleleId(1));
        remap.admit(3, AlleleId(2));
        remap.admit(4, AlleleId(3));
        remap
    }

    /// **The dropped allele returns `None`, and every survivor returns its dense id** —
    /// the correspondence nothing else in the caller records.
    ///
    /// The hole is in the *middle* on purpose: a remapping that simply subtracted a
    /// constant, or that returned the table index unchanged, would agree with this one on
    /// alleles 0 and 1 and disagree on 3 and 4.
    #[test]
    fn the_remapping_answers_none_for_a_dropped_allele_and_a_dense_id_for_a_survivor() {
        let remap = remap_with_a_hole();
        assert_eq!(remap.candidate_for(0), Some(AlleleId::REFERENCE));
        assert_eq!(remap.candidate_for(1), Some(AlleleId(1)));
        assert_eq!(remap.candidate_for(2), None, "the dropped one");
        assert_eq!(remap.candidate_for(3), Some(AlleleId(2)));
        assert_eq!(remap.candidate_for(4), Some(AlleleId(3)));
        assert_eq!(remap.table_len(), 5);
        assert_eq!(remap.num_admitted(), 4);
    }

    /// **The ids the remapping hands out are dense and start at the reference**, which is
    /// what lets a candidate id index [`CandidateAlleles`] directly. Written as a property
    /// over the survivors rather than as four literals, so it still means something if the
    /// fixture changes.
    #[test]
    fn the_admitted_ids_are_dense_from_the_reference_with_no_gaps() {
        let remap = remap_with_a_hole();
        let mut admitted: Vec<u16> = (0..remap.table_len())
            .filter_map(|table_index| remap.candidate_for(table_index))
            .map(AlleleId::get)
            .collect();
        admitted.sort_unstable();
        let dense: Vec<u16> = (0..admitted.len() as u16).collect();
        assert_eq!(admitted, dense);
    }

    /// A support row naming an allele outside the merge's table is a caller bug, and
    /// spec §8 makes it an assertion rather than a value that flows on.
    #[test]
    #[should_panic(expected = "named allele 5 of a merge table holding 5")]
    fn a_row_naming_an_allele_the_table_does_not_hold_is_a_caller_bug() {
        let _ = remap_with_a_hole().candidate_for(5);
    }

    /// **Admitting one merge allele twice is refused**, because the second write would
    /// leave the first candidate id pointing at nothing and one merge allele silently
    /// sharing another's evidence.
    #[test]
    #[should_panic(expected = "admitted twice")]
    fn one_merge_allele_cannot_be_admitted_twice() {
        let mut remap = AlleleRemap::with_all_dropped(3);
        remap.admit(0, AlleleId::REFERENCE);
        remap.admit(1, AlleleId(1));
        remap.admit(1, AlleleId(2));
    }

    /// **Two merge alleles cannot be recorded onto one candidate id.** A bounds check
    /// cannot see this — both indices are in range and each is written once — and the
    /// consequence is worse than the off-by-one the double-admission check catches: the
    /// evidence hand-off would re-key two different sequences' reads onto one candidate,
    /// and the read likelihood would score two alleles as one.
    #[test]
    #[should_panic(expected = "dense and in admission order")]
    fn two_merge_alleles_cannot_share_one_candidate_id() {
        let mut remap = AlleleRemap::with_all_dropped(3);
        remap.admit(0, AlleleId::REFERENCE);
        remap.admit(1, AlleleId(1));
        remap.admit(2, AlleleId(1));
    }

    /// An id that skips ahead leaves a candidate table entry with no merge allele behind
    /// it, which is the same defect seen from the other side.
    #[test]
    #[should_panic(expected = "dense and in admission order")]
    fn a_candidate_id_cannot_skip_ahead_of_the_table() {
        let mut remap = AlleleRemap::with_all_dropped(3);
        remap.admit(0, AlleleId::REFERENCE);
        remap.admit(1, AlleleId(7));
    }

    /// The write-side range check, whose absence would let a survivor be recorded past the
    /// end of the merge's table.
    #[test]
    #[should_panic(expected = "cannot admit allele 3 of a merge table holding 3")]
    fn a_survivor_cannot_be_recorded_past_the_end_of_the_table() {
        let mut remap = AlleleRemap::with_all_dropped(3);
        remap.admit(3, AlleleId::REFERENCE);
    }

    /// **A locus that dropped nothing has a zero pool, with no branch taken to produce
    /// it** — the property spec §12 asks for, and the reason [`UnmatchedSupport`] derives
    /// `Default` rather than being built by a constructor that could round to zero.
    #[test]
    fn a_sample_with_nothing_dropped_has_a_zero_pool() {
        let empty = UnmatchedSupport::default();
        assert_eq!(empty.num_reads, 0);
        assert_eq!(empty.q_sum, 0.0);
        assert!(
            empty.q_sum.is_sign_positive(),
            "an empty pool is zero, not negative zero, so a later sum cannot inherit a sign"
        );
        assert_eq!(empty.earned_reads_cut_by_the_cap, 0);
        assert!(
            !empty.genotype_must_be_missing(),
            "a sample that lost nothing keeps its genotype"
        );
    }

    /// **Error reads in the pool must not cost a sample its genotype**, and this is the
    /// assertion that separates the two ways a sequence leaves the table.
    ///
    /// The fixture is the ordinary case, not a corner: a sample with 40 reads of pooled
    /// error mass, all of it on sequences the support rule rejected. On the GIAB trio at
    /// 300× the support rule rejects 13,166 of 15,474 alternatives, so nearly every sample
    /// at nearly every locus has a pool like this — a rule keyed on the pool's size would
    /// emit a missing genotype almost everywhere.
    #[test]
    fn a_pool_of_rejected_error_reads_does_not_cost_a_sample_its_genotype() {
        let error_only = UnmatchedSupport {
            num_reads: 40,
            q_sum: -123.5,
            earned_reads_cut_by_the_cap: 0,
        };
        assert!(!error_only.genotype_must_be_missing());
    }

    /// **A sample that lost a sequence it had earned cannot be genotyped here**, however
    /// little of its depth that sequence took — one read is enough, because the sequence
    /// cleared the support rule for this sample before the cap removed it, and every
    /// genotype the caller can now form for that sample is over a set that does not
    /// contain what it carries.
    #[test]
    fn a_sample_whose_earned_sequence_the_cap_cut_is_emitted_as_missing() {
        let truncated = UnmatchedSupport {
            num_reads: 40,
            q_sum: -123.5,
            earned_reads_cut_by_the_cap: 1,
        };
        assert!(truncated.genotype_must_be_missing());
    }

    /// The verdict's two ordinary values are distinguishable, and `Truncated` carries how
    /// many alternatives the *cap* cut — never how many the bar rejected.
    #[test]
    fn the_verdict_separates_a_full_list_from_a_truncated_one() {
        assert_ne!(
            SelectionVerdict::Selected,
            SelectionVerdict::Truncated { dropped: 0 }
        );
        assert_eq!(
            SelectionVerdict::Truncated { dropped: 3 },
            SelectionVerdict::Truncated { dropped: 3 }
        );
        assert_ne!(
            SelectionVerdict::Truncated { dropped: 3 },
            SelectionVerdict::Truncated { dropped: 4 }
        );
    }

    /// **The scratch is emptied, not merely resized**, so that no locus can read a value
    /// an earlier locus left behind.
    ///
    /// The fixture writes a non-default summary into the buffer, resets to a *larger*
    /// table, and requires every entry to be the default — a `resize` without a `clear`
    /// would grow the buffer and leave the written row in place at index 0, which is the
    /// mistake this method exists to prevent and which no downstream test could see.
    #[test]
    fn resetting_the_scratch_leaves_no_value_from_an_earlier_locus() {
        let mut scratch = SelectionScratch::new();
        scratch.reset_for(2);
        scratch.per_allele[0] = AlleleSummary {
            best_within_sample_share: 0.9,
            samples_clearing_the_bar: 7,
            cohort_reads: 42,
        };
        scratch.ranked_table_indices.push(1);

        scratch.reset_for(4);

        assert_eq!(scratch.table_len(), 4);
        assert!(
            scratch
                .per_allele
                .iter()
                .all(|summary| *summary == AlleleSummary::default()),
            "a reset locus must not see the previous locus's fold"
        );
        assert!(scratch.ranked_table_indices.is_empty());
    }

    /// Resetting to a *smaller* table shrinks the visible buffer.
    ///
    /// **This is not redundant with the test above, and the reason is narrower than it
    /// looks.** The grow test catches the three obvious wrong versions — a `resize` with
    /// no `clear`, a `clear` that resizes to the old length, a `truncate` in place of the
    /// `clear`. What this one catches *alone* is a fourth: a `reset_for` that overwrites
    /// each row in place and grows when it must, but never shrinks, leaving `table_len()`
    /// over-reporting and the fold walking alleles the merge's table does not hold.
    #[test]
    fn resetting_the_scratch_to_a_smaller_table_shrinks_it() {
        let mut scratch = SelectionScratch::new();
        scratch.reset_for(9);
        scratch.reset_for(2);
        assert_eq!(scratch.table_len(), 2);
    }

    /// One allele's summary clears the bar exactly when some single sample's reads reached
    /// it — the admission rule of spec §3.2, and the reason the flag that used to sit
    /// beside this count was removed.
    #[test]
    fn an_allele_clears_the_bar_when_one_sample_reached_it() {
        let unreached = AlleleSummary::default();
        assert!(!unreached.cleared_the_bar());
        let one_sample = AlleleSummary {
            samples_clearing_the_bar: 1,
            ..AlleleSummary::default()
        };
        assert!(
            one_sample.cleared_the_bar(),
            "one sample suffices — the cohort never enters the bar (spec §3.2)"
        );
    }

    /// A `LocusSelection` over `alleles` sequences, with one leftover per covering sample
    /// and a remapping that admitted every one of them from a merge table of the same size.
    fn selection_over(alleles: usize, covering_samples: usize) -> LocusSelection {
        let mut table = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        let mut remap = AlleleRemap::with_all_dropped(alleles);
        remap.admit(0, AlleleId::REFERENCE);
        for alternative in 1..alleles {
            let id = table.admit(Box::from(b"C".as_slice()));
            remap.admit(alternative, id);
        }
        LocusSelection::new(
            table,
            SelectionVerdict::Selected,
            vec![UnmatchedSupport::default(); covering_samples],
            remap,
            covering_samples,
        )
    }

    /// **The alternative count is the table's length less the reference**, which is the
    /// number the genotype prior divides its concentration by — so an answer one too high
    /// dilutes every real allele at every locus, and nothing panics.
    ///
    /// Both ends are asserted because they fail differently: at a locus that selected down
    /// to the reference alone the right answer is zero, and that is more than one built
    /// locus in four on both benchmarks.
    #[test]
    fn the_alternative_count_excludes_the_reference() {
        assert_eq!(selection_over(1, 2).alternative_allele_count(), 0);
        assert_eq!(selection_over(2, 2).alternative_allele_count(), 1);
        assert_eq!(selection_over(6, 2).alternative_allele_count(), 5);
    }

    /// The leftover runs parallel to the locus's covering samples, and a length that does
    /// not match is refused rather than silently shifting every sample's pool onto its
    /// neighbour.
    #[test]
    #[should_panic(expected = "runs parallel to the locus's covering samples")]
    fn a_leftover_shorter_than_the_covering_samples_is_refused() {
        let mut remap = AlleleRemap::with_all_dropped(1);
        remap.admit(0, AlleleId::REFERENCE);
        let _ = LocusSelection::new(
            CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic),
            SelectionVerdict::Selected,
            vec![UnmatchedSupport::default(); 2],
            remap,
            3,
        );
    }

    /// A remapping that admitted fewer alleles than the table holds means a candidate has
    /// no evidence behind it, which is the other half of the same invariant.
    #[test]
    #[should_panic(expected = "every admitted allele must be in the candidate table")]
    fn a_table_wider_than_the_remapping_admitted_is_refused() {
        let mut table = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        table.admit(Box::from(b"C".as_slice()));
        let mut remap = AlleleRemap::with_all_dropped(2);
        remap.admit(0, AlleleId::REFERENCE);
        let _ = LocusSelection::new(
            table,
            SelectionVerdict::Selected,
            vec![UnmatchedSupport::default()],
            remap,
            1,
        );
    }
}
