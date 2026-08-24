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
//! lent it at least `max(2 reads, 10 in 100 of that sample's reads at the locus)`, and a
//! locus is called over at most six alleles counting the reference. **No term of the bar
//! reads the cohort** — one sample reaching it admits the sequence for everyone —
//! because otherwise a sample's candidate list would depend on who else is in the run
//! (spec §3.2).
//!
//! **"The bar" and "the admission rule" are the same thing here, and this file uses both.**
//! The identifiers say *bar* — `cleared_the_bar`, `samples_clearing_the_bar` — and the prose
//! says whichever reads better in the sentence; the spec does the same (§3 is "The admission
//! rule" and its paragraphs say "the bar"). **The cap is a different thing**: the bar decides
//! which sequences are *worth* calling over, one sample at a time, and the cap decides how
//! many a locus may be called over at all. Which of the two dropped an allele is what decides
//! whether a sample keeps its genotype (§5), so the two words are never interchangeable.

pub mod generic;

use crate::ng::calling::CandidateAlleles;
use crate::ng::run::cohort_merge::build::SupportedAllele;
use crate::ng::run::cohort_merge::build::{CohortObservation, SampleSupport};
use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare, MinAltReads};
use crate::ng::types::{AlleleId, GenomeRegion};
use std::cmp::Ordering;

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

/// **Two reads, or 10 in 100 of that sample's reads at the locus, whichever is more.**
///
/// **The floor is the merge's own** ([`MinAltObs::DEFAULT`], production's number) **and
/// should stay at 2**: measured against the GIAB trio's v4.2.1 truth set over 572 kb on
/// 2026-08-24, **at 30×** raising it from 2 to 3 loses five true alternative alleles,
/// where raising the share to 10 in 100 loses two for the same reduction in table size —
/// 1,539 alternatives kept against 1,601. The floor is the expensive knob (spec §3.3).
///
/// **The share is 10 in 100 where the merge's keep rule uses 2** (owner's decision,
/// 2026-08-24), **and it is set against the recall measurement rather than by it.** On
/// this trio **at 300×** it cuts the merge's 15,474 alternatives to 1,273 where a bar of
/// 2 reads alone keeps 10,793, and **it costs two true alleles that 5 in 100 keeps** —
/// `chr1:193718424` `T→C` at 6 of one sample's 107 compared reads, and `chr1:120579074`
/// `C→A` at 4 of 42.
///
/// **The reason it is 10 and not 5 is the axis recall does not measure.** An allele one
/// sample shows at a twentieth of its reads is far likelier to be error than variation,
/// and every one admitted is a column in every genotype table, at every sample, for the
/// life of the locus — memory and wall time spent on candidates that are mostly not real.
/// Two true alleles in 913, at loci whose local depth is 107 and 42 reads rather than the
/// run's nominal 300, is the price accepted for that. **Provisional in the strict sense:
/// it is to be re-decided once the calling loop and emission exist and the cost can be
/// measured** rather than reasoned about — recall is only one side of it, and nothing here
/// yet weighs the other.
///
/// **It is inert below 21 compared reads a sample**, which is the arithmetic rather than
/// a measurement: `ceil(0.10 × 20) = 2` is the floor, and 21 is the first count at which
/// the share asks for more. So a tomato-depth run — about 11 compared reads a sample at a
/// locus — sees the identical rule it would have seen at 2 in 100, and **the whole of this
/// decision lands at depth and nowhere else.** What was measured on that panel is the
/// neighbouring comparison: turning the share off entirely against the merge's 2 in 100
/// moves 4 loci in 53,935 (spec §3.3).
///
/// **Soft, and the softest constant in this module.** Measured for recall on one human trio
/// over 572 kb (spec §11, Q3); the candidate-count argument that decided it is not measured
/// at all yet.
pub const DEFAULT_MIN_ALLELE_SUPPORT: MinAltReads = MinAltReads {
    floor: MinAltObs::DEFAULT,
    share: MinAltReadShare::new_or_panic(0.10),
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
/// **Those counts were taken with the merge's 2-in-100 share, not the 10 in 100 this
/// module ships**, which spec §4.2 states in its own header and a reader of this constant
/// would otherwise not know. The direction is safe: at tomato depth the two shares are
/// provably the same rule (see [`DEFAULT_MIN_ALLELE_SUPPORT`]), and everywhere else 10 in 100
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
    /// The cap bound, and `dropped` alternatives were cut — the lowest-ranked first, by the
    /// ranking [`compare_best_first`] defines: the largest share of one sample's compared
    /// reads, then how many samples cleared the bar, then the cohort's read total, then the
    /// bases. The reference is never among them (spec §4.1).
    Truncated {
        /// How many alternatives the cap removed. Not how many were dropped in total: an
        /// alternative that failed the bar was never a candidate for the cap.
        ///
        /// **A `u32`, because a `u16` could not count what the cap can cut** (owner's
        /// decision, 2026-08-24). A review built a locus of 70,001 alternatives that all
        /// cleared the bar, and narrowing 69,995 into a `u16` panics — so the one input on
        /// which this step would have refused a locus is the one spec §4.1's whole argument
        /// says must be truncated instead. Nothing upstream bounds a locus's allele table at
        /// 65,536: [`CandidateAlleles::admit`]'s refusal at that width guards the *candidate*
        /// table, which the cap holds at six. It takes 65,536 distinct sequences at one
        /// position, which neither benchmark approaches; the two extra bytes buy the
        /// guarantee rather than a measured case.
        dropped: u32,
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
    /// cuts sequences that already cleared the bar for *somebody* — but not necessarily
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
    /// **The fields are private, so this is the only door** (the owner's ruling at
    /// Checkpoint A, departing from arch §2.4's `pub` fields): a bypassable check on a value
    /// whose defect is a wrong genotype rather than a crash is not a check. Every value a run
    /// produces comes through here, tests included.
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
    /// Merge table indices. Only the alternatives that cleared the bar ever enter it.
    ///
    /// **It holds two orders in turn, and the second is not decoration.** While the cap
    /// chooses, the indices are ordered by [`compare_best_first`]; once it has, the survivors
    /// are sorted back into the merge table's own index order, because that is the order they
    /// are admitted to [`CandidateAlleles`] in (arch §3.1). So the ranking decides *which*
    /// alleles survive and the merge decides what order they appear in — which is the order
    /// that reaches the VCF's `ALT` column.
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

    /// **The cap's first ranking key for one allele of the last locus folded** — the largest
    /// share of one sample's compared reads that allele took, over the samples that cleared
    /// the bar for it (spec §4.1).
    ///
    /// **Zero means no sample cleared the bar for it**, which is not the same as no sample
    /// showing it: those reads are in [`cohort_reads_of`](Self::cohort_reads_of). Allele 0 is
    /// the reference, which is folded like any other but exempt from the bar (spec §6.1).
    ///
    /// **This exists for `examples/ng_candidate_selection_probe.rs` and nothing in the
    /// pipeline calls it** (step D1). The measurement reports the ranking's own keys, and the
    /// only other way to report them is to recompute them — which is the duplicate rule D1
    /// exists to delete, and which had drifted from this one twice over. Two scalars are
    /// exported rather than the fold's own type, so the shape of the computation stays
    /// private (arch §2.4).
    ///
    /// **It reads a reused buffer**, so it means nothing until a fold has filled it and it
    /// means the previous locus after the next [`reset_for`](Self::reset_for). Ask it between
    /// a [`select_generic`](generic::select_generic) on a locus and the next call on any
    /// other.
    ///
    /// # Panics
    ///
    /// If `table_index` is not an allele of the locus last folded — a measurement reading a
    /// buffer sized for a different locus, which would otherwise answer with a neighbour's
    /// number.
    #[inline]
    pub fn best_within_sample_share_of(&self, table_index: usize) -> f64 {
        self.per_allele[table_index].best_within_sample_share
    }

    /// **One allele's reads across the cohort, over every covering sample** — including the
    /// samples that did not clear the bar for it, which is what separates it from
    /// [`best_within_sample_share_of`](Self::best_within_sample_share_of).
    ///
    /// The ranking's third key, and production's only one. Same buffer, same lifetime and
    /// same panic as the accessor above.
    ///
    /// # Panics
    ///
    /// If `table_index` is not an allele of the locus last folded.
    #[inline]
    pub fn cohort_reads_of(&self, table_index: usize) -> u64 {
        self.per_allele[table_index].cohort_reads
    }
}

/// **One allele's fold across every covering sample** — what the bar and the cap read.
///
/// Private to the module: it is the shape of a computation, not part of what selection
/// hands back. The fold that fills it is step B1 and the ranking that reads it is B2.
///
/// **Two of its three fields are readable from outside, one scalar at a time** —
/// [`SelectionScratch::best_within_sample_share_of`] and
/// [`SelectionScratch::cohort_reads_of`], which exist for
/// `examples/ng_candidate_selection_probe.rs` (step D1) and are the reason the measurement
/// no longer carries its own copy of the ranking's keys. **The type itself stays private**,
/// so nothing outside can hold, build or match on a summary.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
struct AlleleSummary {
    /// **The largest share of one sample's compared reads this allele took, over the samples
    /// that cleared the admission rule for it** — the cap's first ranking key (spec §4.1).
    /// Not a cohort share: a cohort total would truncate the private alleles first at scale,
    /// which is the one thing the ranking exists to avoid.
    ///
    /// **Which samples the maximum runs over is the whole of it, and the wider reading was
    /// wrong** (owner's decision, 2026-08-24; spec §4.1). Maximised over *every* sample, a
    /// sample with one compared read contributes a share of 1.0 to whatever that read landed
    /// on — and because this is the ranking's first key, no later key can overturn it. The
    /// case a review built and this now refuses: an allele 40 samples carry at 150 of their
    /// 300 reads each, 6,000 reads across the cohort, cut by the cap in favour of an allele
    /// one sample carries, because a *second* sample's lone read landed on it — after which
    /// the leftover's second count emits all 40 carriers as a missing genotype.
    ///
    /// **Zero where no sample cleared the rule**, which is not a share of nothing but an
    /// allele the cap never ranks: only alternatives that cleared the rule enter
    /// [`SelectionScratch::ranked_table_indices`].
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
    fn cleared_the_bar(self) -> bool {
        self.samples_clearing_the_bar > 0
    }
}

/// **How many of one sample's reads at this locus were read off whole and compared with
/// the reference** — the denominator every question in this module divides by
/// (spec §1.3, §3).
///
/// It is the sum of that sample's rows, over alleles **and** over read groups, because the
/// merge admits only complete observations onto alleles: a read that reached a row is a
/// read that spanned the locus and was compared. **The reference's own rows are in it** —
/// the bar is a share of the sample's reads, not of its non-reference reads, which is what
/// makes a heterozygote's alternative sit near a half rather than near one.
///
/// **Three of the merge's other counts are deliberately not in it** (spec §5.1), and each
/// would make the bar ask a sample for more than its comparable reads can answer:
/// [`partials`](SampleSupport::partials), which say the sample carries *at least* this
/// rather than what it carries and are scored on their own axis;
/// [`reads_without_observation`](SampleSupport::reads_without_observation), which showed
/// nothing; and [`reads_removed_as_evidence`](SampleSupport::reads_removed_as_evidence),
/// which the merge withheld from the table. None of the three ever reached an allele, so
/// none of them can be in a numerator either, and a denominator they entered alone would
/// be a bar that rises with a sample's *unusable* depth.
///
/// Saturating, like the merge's own sums, which needs four billion reads of one sample at
/// one locus to bind.
fn compared_reads_of(sample: &SampleSupport) -> u32 {
    sample
        .supported
        .iter()
        .fold(0, |total, row| total.saturating_add(row.support.num_reads))
}

/// **What selection dropped for one sample, and what it costs that sample** (spec §5).
///
/// Every allele the admission rule or the cap removed keeps its reads' error mass in the
/// arithmetic: the SNP/indel genotype likelihood carries a term for reads matching no candidate
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.3's `q_sum_other`). **Nothing upstream produces it,
/// because nothing upstream drops anything** — selection creates the pool, so selection owes it.
/// It costs no new producer: the merge already stores `q_sum` per `(allele, read group)`, and this
/// is a sum over the rows whose allele is no longer in the table.
///
/// **The mass is summed straight from the merge's rows and never re-derived from a count and a
/// rate**, which is why this walks the rows rather than multiplying anything.
///
/// **What is not in it** (spec §5.1): partial reads, which say the sample carries *at least* this
/// and are scored on their own axis; reads that produced no observation and reads removed as
/// evidence, neither of which carries a quality sum and neither of which was ever in the allele
/// table; and the reference's own reads, since the reference is never dropped and so always has a
/// candidate id.
///
/// # The second count, and why it is not the first one narrowed
///
/// [`UnmatchedSupport::earned_reads_cut_by_the_cap`] is this sample's reads on an allele **this
/// sample's own reads earned** and the **cap** then cut. Non-zero means the sample carries
/// something the locus is no longer called over, so every genotype the caller can form for it is
/// wrong, and emission writes a missing genotype (spec §4.1, §5; the owner's decision of
/// 2026-08-24).
///
/// **Keying it on the pool instead would no-call almost everybody.** The rule drops alleles almost
/// nobody showed — 13,166 of 15,474 alternatives on the GIAB trio at 300× (spec §3.3) — which are
/// overwhelmingly sequencing error, and every sample carries a few error reads at nearly every
/// locus, so a sample's pool is almost always non-zero and says nothing. The cap is the other way,
/// and it only ever cuts alleles that already cleared the bar for *somebody* — but not necessarily
/// for the sample being counted, which is why "cleared it for this sample" is the condition.
///
/// **And asking that needs no list of what the cap cut.** A sample that cleared the bar for an
/// allele is, by construction, a sample that put that allele among the cap's candidates — so an
/// allele this sample earned and the remapping no longer holds can only have been cut by the cap.
/// The test is therefore `dropped && this sample reached the bar`, asked with the same pooled
/// reads and the same denominator the fold used, which is what [`one_run_per_allele`] and
/// [`compared_reads_of`] exist to guarantee.
///
/// # Panics
///
/// On a dropped allele whose quality mass is not finite, and on a sample's rows out of ascending
/// allele order (see [`one_run_per_allele`]). Both are caller bugs held in release, which spec §8
/// names — **this is the third of the three it lists, and this step is where it becomes
/// reachable**, because nothing before C3 read `q_sum` at all. A non-finite mass is not a crash
/// waiting to happen: it flows into the pool, the pool into every genotype's data likelihood, and
/// a non-finite likelihood prefers no genotype over any other, so the locus comes out called with
/// nothing chosen and nothing failed.
fn leftover_of(
    sample: &SampleSupport,
    locus: GenomeRegion,
    remap: &AlleleRemap,
    min_allele_support: MinAltReads,
) -> UnmatchedSupport {
    let compared_reads = compared_reads_of(sample);
    let mut leftover = UnmatchedSupport::default();
    for rows in one_run_per_allele(sample, locus) {
        let allele = rows[0].allele;
        if remap.candidate_for(allele).is_some() {
            continue;
        }
        let pooled_reads = rows.iter().fold(0_u32, |total, row| {
            total.saturating_add(row.support.num_reads)
        });
        leftover.num_reads = leftover.num_reads.saturating_add(pooled_reads);
        let mass = rows.iter().map(|row| row.support.q_sum).sum::<f64>();
        assert!(
            mass.is_finite(),
            "sample {}'s rows for allele {allele} at {locus} sum to a quality mass of {mass}, \
             which is not a number the arithmetic can carry: it flows into the pool, the pool \
             into every genotype's data likelihood, and a non-finite likelihood prefers no \
             genotype at all — so the locus comes out called with nothing chosen and nothing \
             failed",
            sample.sample
        );
        leftover.q_sum += mass;
        if min_allele_support.reached_by(pooled_reads, compared_reads) {
            leftover.earned_reads_cut_by_the_cap = leftover
                .earned_reads_cut_by_the_cap
                .saturating_add(pooled_reads);
        }
    }
    leftover
}

/// **One sample's support rows, grouped so each allele's read-group rows arrive together** —
/// the one place this module decides what "a sample's reads for this allele" means, and it is
/// shared so that the two walks over these rows cannot come to disagree.
///
/// The merge writes the rows in ascending `(allele, read group)` order, so one allele's rows are
/// one contiguous run and `chunk_by` groups exactly the read groups. **A read is a read whichever
/// lane produced it**, and asking a rule of each row separately would be a stricter rule applied
/// to exactly the samples that carry more than one read group (`doc/devel/ng/spec/read_groups.md`
/// §1 — 157 of 1,707 in a surveyed tomato archive carry more than one).
///
/// **Both walks of a locus go through here**: [`summarise_alleles`], which asks the admission
/// rule, and the leftover, which sums what selection dropped and asks the same bar again per
/// sample. If those two pooled differently, a sample with two libraries could clear the bar in
/// one walk and not the other, and its genotype would be emitted as missing for no reason.
///
/// # Panics
///
/// On rows out of ascending allele order. Out of order, one allele's rows split into two runs and
/// each asks the bar with part of the sample's reads, so a sequence the sample really did earn
/// quietly fails.
///
/// **The check runs as the iterator is advanced, not when it is built**, so a caller that stops
/// early can walk past a disorder it never reaches — `.next()` alone on rows ordered 1, 0, 1
/// returns without panicking. Both callers here exhaust it, which is what makes the guarantee
/// hold today; a caller that does not must not treat this as a validation.
fn one_run_per_allele(
    sample: &SampleSupport,
    locus: GenomeRegion,
) -> impl Iterator<Item = &[SupportedAllele]> {
    let mut previous_allele: Option<usize> = None;
    sample
        .supported
        .chunk_by(|left, right| left.allele == right.allele)
        .inspect(move |rows| {
            let allele = rows[0].allele;
            assert!(
                previous_allele.is_none_or(|previous| previous < allele),
                "sample {}'s rows must be in ascending allele order, and allele {allele} \
                 follows {previous_allele:?} at {locus}: out of order, one allele's rows split \
                 into two runs and each asks the bar with part of the sample's reads, so a \
                 sequence the sample really did earn quietly fails",
                sample.sample
            );
            previous_allele = Some(allele);
        })
}

/// **Fold one locus's rows into one summary per allele, asking the bar of each sample
/// separately** — the pass that decides which sequences are worth calling over
/// (spec §3; arch §3.1).
///
/// One sample at a time: its denominator is [`compared_reads_of`], and then each allele it
/// showed — **its read-group rows pooled** — is asked whether it reached the bar against
/// that denominator. **One sample reaching it admits the sequence for the whole cohort**
/// (spec §3.2), so the count of samples that did is not a term of the bar; it is only the
/// cap's first tie-break, where it can reorder and never exclude.
///
/// **A sample that did not reach the rule lends the allele nothing** — not a share, not a
/// count — and only its reads reach the cohort total (owner's decision, 2026-08-24; spec
/// §4.1). The share and the count are raised in the same branch for that reason;
/// [`AlleleSummary::best_within_sample_share`] carries what the other reading cost.
///
/// **The read-group rows are pooled here, and this is the one place pooling is right**
/// (arch §3.1). The merge keys its rows on `(allele, read group)` because a read likelihood
/// may fold reads into one term only when every one of them would get the same number, and
/// two lanes have different error rates (`doc/devel/ng/spec/read_likelihoods.md` §2.3). The
/// bar counts reads, and a read is a read whichever lane produced it. Asking it of each row
/// separately would be a **stricter rule applied to exactly the samples that carry more
/// than one read group**, which a surveyed tomato archive found in 157 of 1,707 samples
/// — 133 with two libraries, 20 with three, and four with 7, 16, 16 and 42
/// (`doc/devel/ng/spec/read_groups.md` §1). A sample showing an allele 3 reads from one lane
/// and 2 from another would then be asked to reach a bar of 5 twice over, and would fail it
/// with the 5 reads that clear it.
///
/// **Pooled here rather than through [`SampleSupport::pooled_support_for`]**, which arch §3.1
/// nominates. That method answers one allele, so reaching a sample's distinct alleles through
/// it means scanning the whole locus table per sample — rows the sample never showed
/// included — where grouping the rows walks each sample's own rows once. It also rebuilds all
/// six of the merge's quality moments where the bar reads one of them.
///
/// **The reference is folded like every other allele, the bar included.** Its read total, its
/// within-sample share and its count of samples reaching the bar are all filled, because they
/// are honest facts about it and cost nothing to keep. **Nothing downstream reads that last
/// one**: the reference is exempt from the bar, and step C1 seeds it into the candidate table
/// structurally, before any sample's evidence is read (spec §6.1). It is recorded rather than
/// suppressed so that no reader has to wonder which of the three fields is missing for
/// allele 0.
///
/// # Panics
///
/// On three merge bugs, each an assertion held in release because its symptom would be a
/// wrong candidate list rather than a crash, which is the convention spec §8 sets for this
/// step and names two of these three in:
///
/// - a support row naming an allele the locus's table does not hold;
/// - a sample's rows out of ascending allele order, or a locus's samples out of ascending
///   sample order;
/// - **a sample that has rows and no reads on any of them.** A sample with *no* rows is a
///   different thing and is legitimate — it covers the locus and every one of its reads
///   stopped inside it, so all it has is partials — and that one is stepped over.
fn summarise_alleles(
    observation: &CohortObservation,
    min_allele_support: MinAltReads,
    scratch: &mut SelectionScratch,
) {
    let table_len = observation.alleles.len();
    let locus = observation.region;
    scratch.reset_for(table_len);
    let mut previous_sample: Option<usize> = None;
    for sample in &observation.per_sample {
        assert!(
            previous_sample.is_none_or(|previous| previous < sample.sample),
            "the locus's covering samples must be in ascending sample order, and sample {} \
             follows {previous_sample:?} at {locus}: one sample listed twice is folded twice, \
             which lifts its allele's cohort total and lets one sample clear the bar as two",
            sample.sample
        );
        previous_sample = Some(sample.sample);
        if sample.supported.is_empty() {
            // A covering sample whose reads all stopped inside the locus: it has partials
            // and nothing else, and partials count toward no bar (spec §5.1). Stepped over,
            // and the samples after it are still folded.
            continue;
        }
        let compared_reads = compared_reads_of(sample);
        assert!(
            compared_reads > 0,
            "sample {} has {} support rows and no reads on any of them at {locus}: the merge \
             writes no row for a pair a sample showed no reads for, and a zero denominator \
             would divide the bar's share by nothing",
            sample.sample,
            sample.supported.len()
        );
        for rows in one_run_per_allele(sample, locus) {
            let allele = rows[0].allele;
            let pooled_reads = rows.iter().fold(0_u32, |total, row| {
                total.saturating_add(row.support.num_reads)
            });
            let summary = scratch.per_allele.get_mut(allele).unwrap_or_else(|| {
                panic!(
                    "sample {}'s support row named allele {allele} of a locus whose table \
                     holds {table_len}, at {locus}",
                    sample.sample
                )
            });
            summary.cohort_reads = summary.cohort_reads.saturating_add(u64::from(pooled_reads));
            if min_allele_support.reached_by(pooled_reads, compared_reads) {
                summary.samples_clearing_the_bar =
                    summary.samples_clearing_the_bar.saturating_add(1);
                // The share and the count are raised in the same branch, so they describe the
                // same set of samples — see `best_within_sample_share`'s own documentation for
                // what happened when they did not.
                let share = f64::from(pooled_reads) / f64::from(compared_reads);
                summary.best_within_sample_share = summary.best_within_sample_share.max(share);
            }
        }
    }
}

/// **One alternative as the cap's ranking reads it** — the fold [`summarise_alleles`] filled
/// for it, and the bases that break the last tie.
///
/// **The two travel together so that a call site cannot pair one allele's summary with
/// another's bases**, which is the whole reason the type exists. The comparison used to take
/// four positional arguments — two summaries and two base slices — and a review swapped the
/// two slices at a call site: it compiled, `clippy` was silent, and the test still passed,
/// because the bases only decide when all three numeric keys tie. The mis-pairing is
/// therefore invisible at exactly the loci where the ranking does its work, and the caller
/// this is built for is worse than the test — step C2 sorts a buffer of table indices, so
/// every argument is an index expression.
#[derive(Clone, Copy, Debug)]
struct RankedAlternative<'bases> {
    /// What the fold recorded for this allele across the cohort's samples.
    summary: AlleleSummary,
    /// The allele's own sequence, from `CohortObservation::alleles`.
    bases: &'bases [u8],
}

/// **Order two alternatives best-first for the cap** — the largest share of one sample's
/// compared reads, then how many samples cleared the rule, then the cohort's read total,
/// then the bases (spec §4.1; arch §2.5).
///
/// **The name says the direction because the return cannot.** `Ordering::Less` means `left`
/// belongs earlier in a best-first list, so the three ranking keys are compared
/// right-against-left and only the bases, the last, are compared the way they read. Sorting
/// with it gives best-first, and `min_by` gives the best allele — **`max_by` gives the
/// worst**, which is the one idiom that compiles and means the opposite of what it looks
/// like. *(Arch §2.5 names this `ranks_above`; renamed here because a third-person `-s` verb
/// is the shape Rust reserves for `-> bool`, and answering it with `Less` for "yes" is a trap
/// a doc comment can only paper over.)*
///
/// **The first key is a share of one sample's reads, not of the cohort's, and that is the
/// design.** Production ranks by the cohort's raw read total
/// (`enforce_max_alleles`, `src/var_calling/per_group_merger.rs`, a stable sort on
/// `Reverse(cohort_count)`), and at scale that truncates the *private* alleles first: an
/// allele one sample really carries, heterozygous at 30×, scores 15 reads, where a systematic
/// mismapping artefact at 1 read in 100 across 800 samples at 30× scores 240. As
/// within-sample shares the same two are 0.5 and 0.01. bcftools divides each sample's quality
/// sums by that sample's total before adding across the cohort (`bam2bcf.c`), and HipSTR pools
/// `count/sample_reads` rather than counts (`HaplotypeGenerator.cpp`), for the same reason.
///
/// **The tie-break order is what makes the ranking degrade across the depth range without a
/// branch.** At 300 reads a sample the shares separate and the first key decides. At about 3
/// reads a position every admitted allele takes two of them, so the shares all sit near two
/// thirds, the first key ties, and how many samples cleared the rule is the only signal there
/// is (spec §4.1). Production's key is this ranking's third.
///
/// **In a cohort of mixed depth the first key does not compare amounts of evidence**, which
/// that argument states one depth at a time and never says outright: a homozygous alternative
/// scores 1.0 whatever its depth and a heterozygote about 0.5 whatever its depth, so a sample
/// sequenced at 3 reads outranks every heterozygous sample sequenced at 300 — three agreeing
/// reads beat 150. **Kept, and now stated** (owner's decision, 2026-08-24; spec §4.1): three
/// agreeing reads are evidence, and the alternative is a share shrunk toward a half by the
/// sample's depth, which would change what the key means everywhere to fix a case the depth
/// range makes rare.
///
/// **The bases are the tie-break that cannot tie**, which is why they are here rather than a
/// stable sort of the rows. The merge's table keys its alleles by their own bytes
/// (`AlleleTable`, `src/ng/run/cohort_merge/build.rs`), so two entries always differ and the
/// order can never fall through to however the samples were walked — which spec §8 requires,
/// since the output must be byte-identical at any worker count. Production's ranking rests on
/// a stable sort of equal keys, which is deterministic but arbitrary.
///
/// **Shares compare with [`f64::total_cmp`]**, a total order, so there is no `NaN` branch and
/// no partial-order footgun. Nothing here can be handed a `NaN`: [`summarise_alleles`] asserts
/// a non-zero denominator for every sample that has rows, so every share is one `u32` over
/// another.
fn compare_best_first(left: RankedAlternative<'_>, right: RankedAlternative<'_>) -> Ordering {
    right
        .summary
        .best_within_sample_share
        .total_cmp(&left.summary.best_within_sample_share)
        .then(
            right
                .summary
                .samples_clearing_the_bar
                .cmp(&left.summary.samples_clearing_the_bar),
        )
        .then(right.summary.cohort_reads.cmp(&left.summary.cohort_reads))
        .then(left.bases.cmp(right.bases))
}

/// **Hand-built cohort loci, shared by this module and by [`generic`]'s tests.**
///
/// One home rather than a copy per file, for the reason `cohort_merge`'s own `fixtures`
/// module gives: a locus fixture is fiddly enough that two copies drift, and the two
/// callers here test the same rule from either side of it — the fold in this file, the
/// whole narrowing in [`generic`].
#[cfg(test)]
pub(super) mod fixtures {
    use super::*;
    use crate::ng::locus_generation::WitnessedLocusPositions;
    use crate::ng::run::cohort_merge::build::{AlleleSupport, PartialObservation, SupportedAllele};
    use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};
    use std::num::NonZeroU32;

    /// A support rule of `floor` reads or `share` of a sample's compared reads.
    ///
    /// Written out rather than taken from [`DEFAULT_MIN_ALLELE_SUPPORT`] because the
    /// shipped share of 10 in 100 is inert below 21 compared reads, so a fixture of a
    /// handful of reads built on it could not tell a right denominator from a wrong one:
    /// the floor would decide either way. **Raising the share is not by itself enough to
    /// make it decide** — the fixture also has to be deep enough that
    /// `ceil(share × compared reads)` exceeds the floor, which is what
    /// `the_share_refuses_what_the_floor_would_admit` is for.
    pub(super) fn support_rule_of(floor: u32, share: f64) -> MinAltReads {
        MinAltReads {
            floor: MinAltObs(NonZeroU32::new(floor).expect("a floor of at least one read")),
            share: MinAltReadShare::new(share).expect("a share that is a fraction of one"),
        }
    }

    /// One allele's row from the sample's **first** read group — the shape of every sample
    /// that carries one library, which is most samples of most runs.
    ///
    /// `q_sum` is carried although no assertion in this step reads it: it is the field
    /// step C3 sums into the leftover, so the fixtures hold a plausible error mass from
    /// the start rather than gaining one when C3 arrives.
    pub(super) fn row(allele: usize, num_reads: u32, q_sum: f64) -> SupportedAllele {
        row_from_group(allele, ReadGroupId(0), num_reads, q_sum)
    }

    /// One `(allele, read group)` row, naming the group — for the fixtures about two lanes
    /// of one sample. The group is a [`ReadGroupId`] rather than a bare number so that it
    /// cannot be transposed with `num_reads`, which is the neighbouring argument and also
    /// counts something.
    pub(super) fn row_from_group(
        allele: usize,
        read_group: ReadGroupId,
        num_reads: u32,
        q_sum: f64,
    ) -> SupportedAllele {
        SupportedAllele {
            allele,
            read_group,
            support: AlleleSupport {
                num_reads,
                q_sum,
                ..AlleleSupport::default()
            },
        }
    }

    /// One covering sample showing `rows` and nothing on the merge's other axes.
    pub(super) fn sample_showing(sample: usize, rows: Vec<SupportedAllele>) -> SampleSupport {
        SampleSupport {
            sample,
            supported: rows,
            partials: Vec::new(),
            reads_without_observation: 0,
            reads_removed_as_evidence: 0,
            reads_composed_across_records: 0,
        }
    }

    /// A covering sample whose reads **all stopped inside the locus**: no support rows at
    /// all, `num_reads` partial ones. The merge builds this — `per_sample` holds the
    /// samples that covered the span, not the ones that spanned it — and it is the only
    /// input that reaches the fold with a denominator of zero.
    pub(super) fn sample_with_only_partials(sample: usize, num_reads: u32) -> SampleSupport {
        let mut only_partials = sample_showing(sample, Vec::new());
        only_partials.partials = vec![partial_of(num_reads)];
        only_partials
    }

    /// `num_reads` reads that stopped inside the locus — the axis the denominator must not
    /// read.
    pub(super) fn partial_of(num_reads: u32) -> PartialObservation {
        PartialObservation {
            witnessed_in_locus: WitnessedLocusPositions::one_run_from_offset_and_length(0, 2)
                .expect("a two-position witness"),
            read_group: ReadGroupId(0),
            bases: Box::from(b"AC".as_slice()),
            num_reads,
            q_sum: -3.0,
        }
    }

    /// A cohort locus over `alleles`, the reference first, covered by `per_sample`.
    pub(super) fn locus_of(alleles: &[&[u8]], per_sample: Vec<SampleSupport>) -> CohortObservation {
        CohortObservation {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(100),
                end: Position(100),
            },
            alleles: alleles.iter().map(|bases| Box::from(*bases)).collect(),
            per_sample,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fixtures::*;
    use super::*;
    use proptest::prelude::*;

    use crate::ng::locus_generation::LocusKind;
    use crate::ng::types::ReadGroupId;

    /// The two numbers themselves, and the floor's **coupling** to the merge's constant
    /// rather than to the digit 2 — the doc comment says the floor *is* the merge's own,
    /// and only the second assertion holds that.
    #[test]
    fn the_default_bar_is_two_reads_or_ten_in_a_hundred() {
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.floor.get(), 2);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.floor, MinAltObs::DEFAULT);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.share.get(), 0.10);
    }

    /// The two ends of the committed depth range, as spec §3 states them: at 3 compared
    /// reads the rule asks 2, and at 300 it asks 30.
    ///
    /// The third count is there because the first two cannot see the rounding: `0.10 ×
    /// 300` is exactly 30, so rounding the share *down* would answer 30 as well. At 301
    /// the share is 30.1, and up and down are 31 and 30.
    #[test]
    fn the_floor_decides_at_three_reads_and_the_share_at_three_hundred() {
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(3), 2);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(300), 30);
        assert_eq!(DEFAULT_MIN_ALLELE_SUPPORT.required_of(301), 31);
    }

    /// The share is stricter than the merge's own at depth and **indistinguishable from
    /// it below 21 compared reads** — the claim [`DEFAULT_MIN_ALLELE_SUPPORT`]'s
    /// documentation makes, held against the merge's constant rather than against a
    /// number retyped here.
    ///
    /// **20 and 21 are the fixture, and they are what stops this test being vacuous.**
    /// Equal up to 20 and strictly greater at 21 pins the share to **more than 2/21 —
    /// about 0.0952 — and no more than 2/20, which is 0.10**, the narrowest window the
    /// rule's own integer arithmetic can express. Both neighbours the project has argued
    /// over fail it: at 5 in 100, `ceil(0.05 × 21) = 2` is still the floor and the strict
    /// arm fails; at 20 in 100, `ceil(0.20 × 20) = 4` and the equality arm fails.
    #[test]
    fn the_allele_share_binds_only_above_twenty_compared_reads() {
        for compared_reads in [1_u32, 3, 11, 20] {
            assert_eq!(
                DEFAULT_MIN_ALLELE_SUPPORT.required_of(compared_reads),
                MinAltReads::DEFAULT.required_of(compared_reads),
                "at {compared_reads} compared reads the floor decides for both rules"
            );
        }
        assert!(
            DEFAULT_MIN_ALLELE_SUPPORT.required_of(21) > MinAltReads::DEFAULT.required_of(21),
            "21 compared reads is where the allele rule first asks for more than the merge's"
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
    /// type is right, the floor is right, and only the share moves — from 10 in 100 to 2 in
    /// 100. It is invisible at tomato depth, where the two are the same rule, and on the
    /// GIAB trio at 300× it is the difference between keeping 1,273 alternatives and
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

    /// The fold's summaries for `observation` under `min_allele_support`. It stays here
    /// rather than in [`fixtures`] because it reaches into the scratch's private buffer,
    /// which is this file's business and not a locus fixture.
    fn summaries_of(
        observation: &CohortObservation,
        min_allele_support: MinAltReads,
    ) -> Vec<AlleleSummary> {
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, min_allele_support, &mut scratch);
        scratch.per_allele.clone()
    }

    /// **The plan's oracle for this step**, and the four things it separates.
    ///
    /// One sample, whose compared reads are 9 by an independent count — 4 reference reads
    /// and 5 on one alternative — while the merge's other three axes hold 6 partial reads,
    /// 5 that showed nothing and 3 removed as evidence. The rule is 2 reads or half of the
    /// sample's compared reads, so it asks for 5, and the alternative has exactly 5. Every
    /// wrong denominator this step can produce moves the answer:
    ///
    /// - counting the **partials** makes it 15 compared reads and a bar of 8, so the
    ///   alternative fails;
    /// - counting the **silent reads** makes it 14 and a bar of 7, and it fails;
    /// - counting the **reads removed as evidence** makes it 12 and a bar of 6, and it
    ///   fails;
    /// - dropping the **reference's own rows** makes it 5 and a bar of 3, which the
    ///   alternative passes — so the count is asserted directly as well, since that one
    ///   error is invisible in the verdict.
    ///
    /// And the alternative's 5 reads are 3 from one read group and 2 from another, so a
    /// fold that took the larger row instead of the sum would ask 5 of 3 and fail.
    #[test]
    fn the_denominator_is_the_samples_compared_reads_and_nothing_else() {
        let mut sample = sample_showing(
            0,
            vec![
                row(0, 4, -4.0),
                row_from_group(1, ReadGroupId(0), 3, -3.0),
                row_from_group(1, ReadGroupId(1), 2, -2.0),
            ],
        );
        sample.partials = vec![partial_of(6)];
        sample.reads_without_observation = 5;
        sample.reads_removed_as_evidence = 3;

        assert_eq!(
            compared_reads_of(&sample),
            9,
            "four reference reads and five on the alternative, pooled over both read groups"
        );

        let observation = locus_of(&[b"A", b"C"], vec![sample]);
        let summaries = summaries_of(&observation, support_rule_of(2, 0.5));
        assert_eq!(
            summaries[1].samples_clearing_the_bar, 1,
            "half of nine compared reads is five, and the alternative has exactly five"
        );
        assert_eq!(
            summaries[1].best_within_sample_share,
            5.0 / 9.0,
            "the share is of the sample's compared reads, the reference's included"
        );
    }

    /// **The share half of the rule, in the regime where it is the half that decides** —
    /// which no other fixture here enters, because below 21 compared reads the floor
    /// decides for any share up to the shipped 10 in 100.
    ///
    /// One sample with 100 compared reads shows 3 on the alternative, against 2 reads or 5
    /// in 100: the floor would admit it and the share asks for 5, so it is refused; at 5
    /// reads it is admitted. **The 5 in 100 here is the fixture's own and not the shipped
    /// number** ([`DEFAULT_MIN_ALLELE_SUPPORT`] is 10 in 100), for the reason
    /// [`support_rule_of`] gives: a rule written out is a rule a reader can check against
    /// the counts beside it. **Without this pair the whole share term could be deleted and
    /// every other test here would still pass**, and a fold applying the floor alone admits
    /// sequencing error as a candidate allele at 300× — the depth at which the share is the
    /// only half of the rule doing any work (spec §3.3: the shipped 10-in-100 share cuts the
    /// GIAB trio's 15,474 alternatives to 1,273 where a count-only bar keeps 10,793).
    ///
    /// It also pins the **denominator** the share is taken of: asked against the allele's
    /// own 3 reads instead of the sample's 100, `ceil(0.05 × 3) = 1` and the alternative
    /// would clear.
    #[test]
    fn the_share_refuses_what_the_floor_would_admit() {
        let refused = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 97, -97.0), row(1, 3, -3.0)])],
        );
        assert_eq!(
            summaries_of(&refused, support_rule_of(2, 0.05))[1].samples_clearing_the_bar,
            0,
            "three reads clear a floor of two and fail a share of five in a hundred of 100"
        );

        let admitted = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 95, -95.0), row(1, 5, -5.0)])],
        );
        assert_eq!(
            summaries_of(&admitted, support_rule_of(2, 0.05))[1].samples_clearing_the_bar,
            1,
            "and five reads is exactly the share, so this is a boundary and not a refusal \
             of everything at depth"
        );
    }

    /// **A sample's two read groups sum rather than the larger winning**, in the direction
    /// the oracle above cannot show: here each row on its own clears a rule of 2 reads, so
    /// a fold asking the rule once per row would count the sample twice and report an
    /// allele two samples cleared where one did.
    ///
    /// That count is the cap's first tie-break (spec §4.1), so the wrong answer is a
    /// different allele kept at a truncated locus, with nothing failing.
    #[test]
    fn one_samples_two_read_groups_are_one_sample() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 5, -5.0),
                    row_from_group(1, ReadGroupId(0), 3, -3.0),
                    row_from_group(1, ReadGroupId(1), 2, -2.0),
                ],
            )],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(
            summaries[1].samples_clearing_the_bar, 1,
            "one sample showed the allele, from two lanes"
        );
        assert_eq!(
            summaries[1].cohort_reads, 5,
            "and its five reads are five reads whichever lane produced them"
        );
    }

    /// **The count is a count and not a flag** (spec §4.1). Two samples clear the rule on
    /// the alternative and the summary says two.
    ///
    /// **This is the number the cap breaks its first tie on, and at three reads a position
    /// it is the only key that separates two alleles at all** — every admitted allele there
    /// has a share near two thirds, so the share key ties and this one decides (spec §4.1).
    /// A count stuck at one would drop the ranking through to the cohort read total, which
    /// is production's ranking and the one spec §4.1 exists not to be: at a thousand
    /// samples it truncates the private alleles first.
    #[test]
    fn every_sample_that_cleared_the_rule_is_counted() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 1, -1.0), row(1, 3, -3.0)]),
                sample_showing(1, vec![row(0, 1, -1.0), row(1, 4, -4.0)]),
            ],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(summaries[1].samples_clearing_the_bar, 2);
        assert_eq!(summaries[1].cohort_reads, 7);
    }

    /// **The share is the largest of one sample's, never the cohort's** (spec §4.1): the
    /// allele takes 3 of 4 reads in one sample and 1 of 40 in another, and the summary
    /// keeps 0.75 — where a cohort share would be 4 in 44, about 0.09, and would rank a
    /// private allele below a systematic artefact at scale.
    ///
    /// The sample of 40 reads is also what makes the rule's independence visible: at 2
    /// reads or a tenth it asks 4 of that sample and gets 1, so it does not clear, and the
    /// count stays at the one sample that did.
    ///
    /// **The larger share arrives first here**, which is what separates a maximum from a
    /// last-wins assignment; `the_largest_share_wins_when_it_arrives_last` is the other
    /// direction, and both are needed.
    #[test]
    fn the_share_is_one_samples_own_and_the_rule_is_asked_of_each_sample_alone() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 1, -1.0), row(1, 3, -3.0)]),
                sample_showing(1, vec![row(0, 39, -39.0), row(1, 1, -1.0)]),
            ],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.10));
        assert_eq!(summaries[1].best_within_sample_share, 0.75);
        assert_eq!(
            summaries[1].samples_clearing_the_bar, 1,
            "the second sample's one read in forty reaches neither half of the rule"
        );
        assert_eq!(
            summaries[1].cohort_reads, 4,
            "the cohort total is the one place a sum over samples appears"
        );
    }

    /// **The share is maximised in both directions**, and the two fixtures fail different
    /// wrong folds: with the larger share first, a last-wins assignment is caught; with it
    /// last, a first-wins one is.
    ///
    /// A first-wins fold would make the cap's first ranking key depend on the order the
    /// samples happen to be walked in, so which allele survives a truncated locus would
    /// change when the cohort is re-ordered — and spec §8 requires the output to be
    /// byte-identical at any worker count.
    #[test]
    fn the_largest_share_wins_when_it_arrives_last() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 2, -2.0), row(1, 2, -2.0)]),
                sample_showing(1, vec![row(0, 1, -1.0), row(1, 3, -3.0)]),
            ],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(
            summaries[1].samples_clearing_the_bar, 2,
            "both samples must clear the rule, or only one share reaches the maximum and the \
             test cannot tell a maximum from an assignment"
        );
        assert_eq!(summaries[1].best_within_sample_share, 0.75);
    }

    /// **A sample that did not clear the admission rule lends the allele no share** (owner's
    /// decision, 2026-08-24; spec §4.1), and this is the case that decided it.
    ///
    /// One sample carries the allele at 20 of its 100 reads and clears the rule. A second
    /// sample has a single compared read, and it landed there — a share of 1.0, and the floor
    /// of 2 refuses it. Maximised over every sample, that one read would set the allele's
    /// first ranking key to 1.0 and **no later key could overturn it**, because the share is
    /// compared first: at a binding cap it displaces alleles that dozens of samples earned,
    /// and the leftover's second count then emits every one of those samples as missing.
    ///
    /// The cohort read total still counts the refused sample's read, which is the honest
    /// place for it: it is a read on the allele, it is just not evidence about a sample.
    #[test]
    fn a_sample_that_did_not_clear_the_rule_lends_no_share() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 80, -80.0), row(1, 20, -20.0)]),
                sample_showing(1, vec![row(1, 1, -1.0)]),
            ],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.05));
        assert_eq!(
            summaries[1].samples_clearing_the_bar, 1,
            "one compared read reaches neither half of the rule"
        );
        assert_eq!(
            summaries[1].best_within_sample_share, 0.20,
            "the carrier's share, not the single read's 1.0"
        );
        assert_eq!(
            summaries[1].cohort_reads, 21,
            "the refused sample's read is still a read on the allele"
        );
    }

    /// **Adding a sample that shows only reference reads changes nothing about the
    /// alternatives** — spec §3.2's principle, that no term of the rule reads the cohort,
    /// as a test. It is the one that fails first if a cohort term ever creeps into the
    /// denominator or the share.
    #[test]
    fn a_sample_showing_only_reference_reads_changes_no_alternatives_summary() {
        let carrier = || sample_showing(0, vec![row(0, 1, -1.0), row(1, 2, -2.0)]);
        let alone = summaries_of(
            &locus_of(&[b"A", b"C"], vec![carrier()]),
            support_rule_of(2, 0.5),
        );
        let with_a_bystander = summaries_of(
            &locus_of(
                &[b"A", b"C"],
                vec![carrier(), sample_showing(1, vec![row(0, 60, -60.0)])],
            ),
            support_rule_of(2, 0.5),
        );
        assert_eq!(alone[1], with_a_bystander[1]);
    }

    /// **What the fold records for the reference**, pinned rather than left to a reader of
    /// the doc comment. It is folded like every other allele: its reads are totalled, its
    /// within-sample share is kept, and the rule is asked of it too.
    ///
    /// **Nothing downstream reads that last answer** — step C1 seeds the reference into the
    /// candidate table structurally, before any sample's evidence is read (spec §6.1) — so
    /// this test exists to stop a C1 that loops `cleared_the_bar()` over the whole table
    /// from either double-seeding the reference or dropping it, neither of which any other
    /// test here would notice.
    #[test]
    fn the_reference_row_is_folded_like_every_other_allele() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)])],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.5));
        assert_eq!(summaries[0].cohort_reads, 3);
        assert_eq!(summaries[0].best_within_sample_share, 0.75);
        assert_eq!(
            summaries[0].samples_clearing_the_bar, 1,
            "the rule is asked of the reference too; C1 does not read the answer"
        );
    }

    /// A support row naming an allele the locus's table does not hold is a merge bug, and
    /// it is refused rather than folded into a neighbouring allele's summary.
    #[test]
    #[should_panic(expected = "named allele 2 of a locus whose table holds 2")]
    fn a_row_naming_an_allele_the_locus_does_not_hold_is_refused() {
        summaries_of(
            &locus_of(
                &[b"A", b"C"],
                vec![sample_showing(7, vec![row(0, 1, -1.0), row(2, 3, -3.0)])],
            ),
            support_rule_of(2, 0.0),
        );
    }

    /// **Rows out of ascending allele order are refused**, because pooling the read groups
    /// reads them as contiguous runs: split into two runs, one allele asks the rule twice
    /// with part of the sample's reads each time, and a sequence the sample earned quietly
    /// fails. Here 3 reads and 2 reads would each be asked against a bar of 5.
    #[test]
    #[should_panic(expected = "must be in ascending allele order")]
    fn rows_out_of_allele_order_are_refused() {
        summaries_of(
            &locus_of(
                &[b"A", b"C"],
                vec![sample_showing(
                    4,
                    vec![
                        row_from_group(1, ReadGroupId(0), 3, -3.0),
                        row(0, 5, -5.0),
                        row_from_group(1, ReadGroupId(1), 2, -2.0),
                    ],
                )],
            ),
            support_rule_of(2, 0.5),
        );
    }

    /// **A sample listed twice is refused**, which is the same failure shape as rows out of
    /// allele order one level up: folded twice, one sample's evidence lifts its allele's
    /// cohort total and clears the rule as two samples, and the cap's first tie-break moves
    /// with nothing failing.
    #[test]
    #[should_panic(expected = "must be in ascending sample order")]
    fn a_sample_listed_twice_is_refused() {
        let twice = || sample_showing(3, vec![row(0, 1, -1.0), row(1, 3, -3.0)]);
        summaries_of(
            &locus_of(&[b"A", b"C"], vec![twice(), twice()]),
            support_rule_of(2, 0.0),
        );
    }

    /// **A sample with rows and no reads on any of them is refused** — spec §8 names it
    /// among the caller bugs this step asserts on, beside the out-of-range allele index.
    /// The merge cannot build it: it writes no row for a pair a sample showed no reads for.
    ///
    /// The assertion is what makes the two other assertions reachable for such a sample. A
    /// guard that skipped it instead would step over the ordering and range checks as well,
    /// so a zero-read row naming an allele outside the table would pass unnoticed.
    #[test]
    #[should_panic(expected = "no reads on any of them")]
    fn a_sample_with_rows_and_no_reads_is_refused() {
        summaries_of(
            &locus_of(
                &[b"A", b"C"],
                vec![sample_showing(2, vec![row(0, 0, 0.0), row(1, 0, 0.0)])],
            ),
            support_rule_of(2, 0.5),
        );
    }

    /// **A covering sample with no rows at all is legitimate, and is the merge's only
    /// realisable zero denominator**: it covered the locus and every one of its reads
    /// stopped inside it, so all it has is partials — which count toward no rule and enter
    /// no denominator (spec §5.1).
    #[test]
    fn a_covering_sample_whose_reads_all_stopped_inside_the_locus_has_no_compared_reads() {
        assert_eq!(compared_reads_of(&sample_with_only_partials(0, 7)), 0);
    }

    /// That sample is **stepped over and is not a stop**: the samples listed after it are
    /// still folded, which is what separates the `continue` from a `return`.
    #[test]
    fn a_sample_with_only_partial_reads_is_stepped_over_and_the_next_sample_is_folded() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_with_only_partials(0, 7),
                sample_showing(1, vec![row(0, 1, -1.0), row(1, 3, -3.0)]),
            ],
        );
        let summaries = summaries_of(&observation, support_rule_of(2, 0.5));
        assert_eq!(summaries[1].samples_clearing_the_bar, 1);
        assert_eq!(summaries[1].best_within_sample_share, 0.75);
        assert_eq!(summaries[1].cohort_reads, 3);
    }

    /// A locus whose table is the reference alone gets one summary, not none and not two.
    /// **It is a first-class outcome rather than a corner** — selection reaches it at more
    /// than one built locus in four on both benchmarks (spec §6.2) — and the fold has to
    /// size its buffers for it before C1 can return it.
    #[test]
    fn a_reference_only_table_gets_one_summary() {
        let observation = locus_of(&[b"A"], vec![sample_showing(0, vec![row(0, 3, -3.0)])]);
        let summaries = summaries_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].cohort_reads, 3);
    }

    /// A locus no sample covers leaves every summary at its default — and this is the
    /// fixture that catches a reset moved inside the sample loop, which would leave the
    /// *previous* locus's summaries standing for a locus that has no samples to overwrite
    /// them.
    #[test]
    fn a_locus_no_sample_covers_leaves_every_summary_at_its_default() {
        let mut scratch = SelectionScratch::new();
        let covered = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 1, -1.0), row(1, 4, -4.0)])],
        );
        summarise_alleles(&covered, support_rule_of(2, 0.0), &mut scratch);
        let uncovered = locus_of(&[b"A", b"C"], Vec::new());
        summarise_alleles(&uncovered, support_rule_of(2, 0.0), &mut scratch);
        assert_eq!(scratch.table_len(), 2);
        assert!(
            scratch
                .per_allele
                .iter()
                .all(|summary| *summary == AlleleSummary::default())
        );
    }

    /// The fold resets its buffers, so a second locus carries nothing of the first —
    /// neither a summary from an allele the second locus does not have, nor a count added
    /// on top of the first locus's.
    #[test]
    fn folding_a_second_locus_carries_nothing_from_the_first() {
        let mut scratch = SelectionScratch::new();
        let wide = locus_of(
            &[b"A", b"C", b"G"],
            vec![sample_showing(
                0,
                vec![row(0, 1, -1.0), row(1, 4, -4.0), row(2, 4, -4.0)],
            )],
        );
        summarise_alleles(&wide, support_rule_of(2, 0.0), &mut scratch);
        assert_eq!(scratch.table_len(), 3);

        let narrow = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 1, -1.0), row(1, 4, -4.0)])],
        );
        summarise_alleles(&narrow, support_rule_of(2, 0.0), &mut scratch);
        assert_eq!(scratch.table_len(), 2, "the third allele is gone");
        assert_eq!(
            scratch.per_allele[1].cohort_reads, 4,
            "and the second locus's four reads are four, not eight"
        );
        assert_eq!(scratch.per_allele[1].samples_clearing_the_bar, 1);
    }

    // ---- the cap's ranking (step B2) ---------------------------------------------------

    /// One alternative for the ranking: its three folded numbers, spelled in the order the
    /// ranking reads them, and its bases.
    fn alternative(
        best_within_sample_share: f64,
        samples_clearing_the_bar: u32,
        cohort_reads: u64,
        bases: &[u8],
    ) -> RankedAlternative<'_> {
        RankedAlternative {
            summary: AlleleSummary {
                best_within_sample_share,
                samples_clearing_the_bar,
                cohort_reads,
            },
            bases,
        }
    }

    /// `alternatives` sorted best-first by [`compare_best_first`], as their bases, so a
    /// failure names the alleles rather than their positions.
    fn ranked_bases<'bases>(alternatives: &[RankedAlternative<'bases>]) -> Vec<&'bases [u8]> {
        let mut ranked = alternatives.to_vec();
        ranked.sort_unstable_by(|left, right| compare_best_first(*left, *right));
        ranked.into_iter().map(|entry| entry.bases).collect()
    }

    /// **Every pairwise fixture below sets the bases *against* the expected answer**, and
    /// that is deliberate rather than decorative. The bases are the last key, so a fixture
    /// whose winner also sorts first alphabetically is passed by a comparator that consults
    /// nothing else: a review replaced the whole function body with `left.bases.cmp(right.bases)`
    /// and seven of the eight tests here still passed. With the bases opposed, that comparator
    /// fails every one of them.
    ///
    /// The one exception is the test about the bases themselves, where they have to decide.
    #[test]
    fn the_better_ranked_allele_compares_less() {
        let better = alternative(0.6, 1, 10, b"G");
        let worse = alternative(0.2, 1, 10, b"C");
        assert_eq!(compare_best_first(better, worse), Ordering::Less);
        assert_eq!(compare_best_first(worse, better), Ordering::Greater);
    }

    /// **The first key: the largest share one sample gave the allele** — and here it decides
    /// *against* all three of the other keys at once. The allele one sample showed at 3 of
    /// its 4 reads outranks the one 40 samples showed at a fiftieth of theirs, which has
    /// eighty times the cohort reads and forty times the carriers.
    ///
    /// **That inversion is the reason the key is a within-sample share** (spec §4.1). At a
    /// thousand samples, ranking by the cohort total keeps a systematic mismapping artefact
    /// present at 1 read in 100 everywhere and truncates the private allele a single carrier
    /// really has.
    #[test]
    fn the_largest_within_sample_share_decides_first_against_every_other_key() {
        let private = alternative(0.75, 1, 3, b"G");
        let widespread = alternative(0.02, 40, 240, b"C");
        assert_eq!(
            ranked_bases(&[widespread, private]),
            vec![b"G".as_slice(), b"C".as_slice()]
        );
    }

    /// **The second key, in the regime spec §4.1 says it is the only signal there is.** At
    /// about 3 reads a sample every admitted allele takes two of them, so all the shares are
    /// 2/3 and the first key ties; how many samples cleared the rule then decides.
    ///
    /// The remaining two keys are set against the answer — the allele two samples cleared
    /// has 4 cohort reads where the other has 9, and its bases sort later — so a ranking
    /// that skipped this key and fell through to the cohort total, which is production's
    /// key, would return the other order.
    #[test]
    fn how_many_samples_cleared_the_rule_decides_when_the_shares_tie() {
        let two_carriers = alternative(2.0 / 3.0, 2, 4, b"G");
        let one_carrier = alternative(2.0 / 3.0, 1, 9, b"C");
        assert_eq!(
            ranked_bases(&[one_carrier, two_carriers]),
            vec![b"G".as_slice(), b"C".as_slice()]
        );
    }

    /// **The third key**, reached only when the share and the sample count both tie: the
    /// cohort's read total, which is production's *first* key and this ranking's last
    /// numeric one. The bases are set against it, so deleting the key does not leave the
    /// fourth one quietly producing the same answer.
    #[test]
    fn the_cohort_read_total_decides_when_the_share_and_the_sample_count_both_tie() {
        let deeper = alternative(0.5, 2, 30, b"G");
        let shallower = alternative(0.5, 2, 12, b"C");
        assert_eq!(
            ranked_bases(&[shallower, deeper]),
            vec![b"G".as_slice(), b"C".as_slice()]
        );
    }

    /// **The fourth key, which cannot tie.** The merge's table keys its alleles by their own
    /// bytes, so two entries always differ, no two ever compare `Equal`, and the order cannot
    /// fall through to however the samples happened to be walked. Spec §8 requires the output
    /// to be byte-identical at any worker count, and this is what closes it.
    #[test]
    fn the_bases_decide_when_all_three_numbers_tie() {
        let earlier = alternative(0.5, 2, 20, b"AC");
        let later = alternative(0.5, 2, 20, b"AG");
        assert_eq!(
            compare_best_first(earlier, later),
            Ordering::Less,
            "the lexicographically smaller sequence ranks above"
        );
        assert_eq!(
            ranked_bases(&[later, earlier]),
            vec![b"AC".as_slice(), b"AG".as_slice()]
        );
    }

    /// **The whole ranking on one table, and the same table fed in reverse order** — the
    /// plan's oracle for this step. Six alternatives make five adjacent pairs, and each of
    /// the four keys decides at least one of them.
    ///
    /// Reading the expected list against the fixture: `TT` wins on the share alone (0.9).
    /// `AC` and `AG` both sit at 0.5, `AC` cleared the rule in 3 samples against 2 — **and
    /// `AG` carries the larger cohort total, 30 against 8, so a ranking that consulted the
    /// totals before the sample count would swap them.** `AG` then beats `GG` on the share
    /// again (0.5 against 0.4). `GG` beats `CC` on the cohort total with both other numbers
    /// tied. And `CA` and `CC` tie on all three, so the bases decide and `CA` — the
    /// lexicographically smaller — comes first.
    ///
    /// **That one number is what makes this the key-*order* oracle it claims to be.** A
    /// review swapped the second and third keys in the comparator and this test still passed,
    /// because `AC` and `AG` then had equal cohort totals and the third key was a no-op on
    /// the only pair the second key was meant to decide.
    #[test]
    fn every_key_decides_a_pair_and_the_row_order_does_not_matter() {
        let table = [
            alternative(0.9, 1, 5, b"TT"),
            alternative(0.5, 3, 8, b"AC"),
            alternative(0.5, 2, 30, b"AG"),
            alternative(0.4, 2, 30, b"GG"),
            alternative(0.4, 2, 12, b"CC"),
            alternative(0.4, 2, 12, b"CA"),
        ];
        let expected = vec![
            b"TT".as_slice(),
            b"AC".as_slice(),
            b"AG".as_slice(),
            b"GG".as_slice(),
            b"CA".as_slice(),
            b"CC".as_slice(),
        ];
        assert_eq!(ranked_bases(&table), expected);

        let mut reversed = table;
        reversed.reverse();
        assert_eq!(
            ranked_bases(&reversed),
            expected,
            "the ranking must not inherit the order the rows arrived in (spec §8)"
        );
    }

    /// **At 300 reads a sample the first key decides on its own** — the deep end of the
    /// committed range, against the 3-read fixture above. A heterozygote showing 151 of one
    /// sample's 300 reads outranks an error-level allele at 16 in 300, although the second
    /// allele was seen in 3 samples against 1, carries 250 cohort reads against 151, and
    /// sorts first by bases. All three lower keys point the other way, so only the share can
    /// produce this answer.
    #[test]
    fn at_three_hundred_reads_a_sample_the_share_alone_decides() {
        let heterozygous = alternative(151.0 / 300.0, 1, 151, b"G");
        let error_level = alternative(16.0 / 300.0, 3, 250, b"C");
        assert_eq!(
            ranked_bases(&[error_level, heterozygous]),
            vec![b"G".as_slice(), b"C".as_slice()]
        );
    }

    /// **At one sample this ranking and production's are the same ranking**, which is worth
    /// pinning because it bounds what spec §4.1's argument buys: every allele's share has the
    /// same denominator — that one sample's compared reads — so ordering by share and
    /// ordering by the cohort read total agree exactly, and the second and third keys never
    /// come into it.
    ///
    /// **The argument for the within-sample share is a cohort-size argument**, and at the
    /// thin end of the range it is worth nothing. It starts paying as soon as two samples of
    /// different depth are present, which is what the next test shows.
    #[test]
    fn at_one_sample_the_share_ranking_and_the_cohort_total_ranking_agree() {
        let half = alternative(0.5, 1, 50, b"AA");
        let three_tenths = alternative(0.3, 1, 30, b"CC");
        let a_fifth = alternative(0.2, 1, 20, b"GG");
        assert_eq!(
            ranked_bases(&[a_fifth, half, three_tenths]),
            vec![b"AA".as_slice(), b"CC".as_slice(), b"GG".as_slice()]
        );
        let mut by_cohort_reads = [a_fifth, half, three_tenths];
        by_cohort_reads.sort_unstable_by(|left, right| {
            right.summary.cohort_reads.cmp(&left.summary.cohort_reads)
        });
        assert_eq!(
            by_cohort_reads
                .iter()
                .map(|entry| entry.bases)
                .collect::<Vec<_>>(),
            vec![b"AA".as_slice(), b"CC".as_slice(), b"GG".as_slice()],
            "at one sample the two rankings cannot disagree"
        );
    }

    /// **A cohort of mixed depth, which spec §4.1's range argument does not cover** — its two
    /// halves each describe a cohort at one depth, and here the two regimes meet inside one
    /// comparison.
    ///
    /// A sample sequenced at 3 reads, homozygous for its allele, scores a share of 1.0 on
    /// three reads in the whole cohort. A sample sequenced at 300, heterozygous, scores 0.5
    /// on 150. **The shallow sample's allele outranks the deep one**, and at a binding cap it
    /// is the last thing cut. Recorded as a test rather than argued about: three agreeing
    /// reads are genuinely evidence, and whether they should outrank 150 is the owner's call
    /// at Checkpoint B — but the behaviour should not be a surprise to whoever reads a
    /// truncated locus.
    #[test]
    fn a_shallow_homozygote_outranks_a_deep_heterozygote_on_the_first_key() {
        let three_of_three = alternative(1.0, 1, 3, b"G");
        let one_hundred_fifty_of_three_hundred = alternative(0.5, 1, 150, b"C");
        assert_eq!(
            ranked_bases(&[one_hundred_fifty_of_three_hundred, three_of_three]),
            vec![b"G".as_slice(), b"C".as_slice()]
        );
    }

    /// **A table wider than the cap**, which every other fixture here is not: the largest is
    /// exactly six rows and the default cap is six *counting the reference*, so the cap first
    /// bites at seven. What the cap will keep is the first five of this list, and what it
    /// cuts is the last two — the ordering C2 will `truncate` on.
    #[test]
    fn a_table_above_the_cap_ranks_its_survivors_into_the_front() {
        let table = [
            alternative(0.10, 1, 6, b"AAAAAA"),
            alternative(0.90, 3, 60, b"CCCCCC"),
            alternative(0.05, 1, 3, b"GGGGGG"),
            alternative(0.70, 2, 40, b"TTTTTT"),
            alternative(0.60, 5, 90, b"ACACAC"),
            alternative(0.80, 1, 20, b"AGAGAG"),
            alternative(0.65, 1, 15, b"ATATAT"),
        ];
        let ranked = ranked_bases(&table);
        assert_eq!(
            &ranked[..usize::from(DEFAULT_MAX_CANDIDATE_ALLELES.alternatives())],
            &[
                b"CCCCCC".as_slice(),
                b"AGAGAG".as_slice(),
                b"TTTTTT".as_slice(),
                b"ATATAT".as_slice(),
                b"ACACAC".as_slice(),
            ],
            "the five alternatives the default cap leaves room for, best first"
        );
        assert_eq!(
            &ranked[usize::from(DEFAULT_MAX_CANDIDATE_ALLELES.alternatives())..],
            &[b"AAAAAA".as_slice(), b"GGGGGG".as_slice()],
            "and the two the cap would cut, still ordered among themselves"
        );
    }

    /// Two shares that are equal to the bit fall through to the next key rather than holding
    /// their input order.
    ///
    /// **No fixture here can separate [`f64::total_cmp`] from `partial_cmp`, and saying so is
    /// more use than a test that pretends otherwise:** the two differ only on `NaN`, and a
    /// `NaN` share cannot reach this function, because [`summarise_alleles`] asserts a
    /// non-zero denominator for every sample that has rows. Substituting
    /// `partial_cmp(..).unwrap()` leaves every test in this module green. `total_cmp` is
    /// still the right call — it is the spelling that cannot panic if that assertion is ever
    /// relaxed — but its value here is a guarantee, not a behaviour any fixture observes.
    #[test]
    fn two_shares_equal_to_the_bit_fall_through_to_the_next_key() {
        let from_a_quarter = alternative(2.0 / 4.0, 3, 10, b"G");
        let from_a_half = alternative(1.0 / 2.0, 1, 10, b"C");
        assert_eq!(
            from_a_quarter.summary.best_within_sample_share,
            from_a_half.summary.best_within_sample_share
        );
        assert_eq!(
            compare_best_first(from_a_quarter, from_a_half),
            Ordering::Less,
            "the shares tie, so the sample count decides"
        );
    }

    proptest! {
        /// **The comparison is a strict weak ordering**, which `sort_unstable_by` requires and
        /// does not check: given a comparator that is not one, the standard library is
        /// entitled to leave the slice in any order at all, and spec §8 requires the output to
        /// be byte-identical at any worker count.
        ///
        /// Asserted over triples drawn across the shares, sample counts and read totals a
        /// locus can produce, plus two-byte alleles so the last key both ties and decides:
        /// the relation is asymmetric, and it is transitive on both `Less` and on the
        /// equivalence `Equal` induces.
        #[test]
        fn the_ranking_is_a_strict_weak_ordering(
            shares in prop::collection::vec(0.0_f64..=1.0, 3),
            counts in prop::collection::vec(0_u32..8, 3),
            reads in prop::collection::vec(0_u64..40, 3),
            bases in prop::collection::vec("[ACGT]{2}", 3),
        ) {
            let entries: Vec<RankedAlternative<'_>> = (0..3)
                .map(|i| alternative(shares[i], counts[i], reads[i], bases[i].as_bytes()))
                .collect();
            for left in &entries {
                for right in &entries {
                    prop_assert_eq!(
                        compare_best_first(*left, *right),
                        compare_best_first(*right, *left).reverse(),
                        "the comparison must be asymmetric"
                    );
                }
            }
            let (a, b, c) = (entries[0], entries[1], entries[2]);
            if compare_best_first(a, b) == Ordering::Less
                && compare_best_first(b, c) == Ordering::Less
            {
                prop_assert_eq!(compare_best_first(a, c), Ordering::Less);
            }
            if compare_best_first(a, b) == Ordering::Equal
                && compare_best_first(b, c) == Ordering::Equal
            {
                prop_assert_eq!(compare_best_first(a, c), Ordering::Equal);
            }
        }
    }
}
