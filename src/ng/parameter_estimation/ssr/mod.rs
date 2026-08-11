//! The STR path: how often a read at a repeat tract shows a length its allele does not
//! have, and the four numbers that describe it.
//!
//! A read at an ordinary site can only be misread — one base for another. A read at a
//! repeat tract can also **slip**, showing a whole motif copy more or fewer than the
//! allele it was drawn from, and no per-base error rate has anywhere to put that. So this
//! path carries the generic path's noise model and adds slippage to it: how often a read
//! slips at all, which way it slips, how far it slips when it does, and the per-base
//! substitution rate that changes a tract's composition at fixed length.
//!
//! **A stratum is a motif period and a reference repeat count** — one group of loci that
//! gets its own fitted numbers, because how much a tract slips depends on how many copies
//! it holds: 9 reads in 10,000 below four repeats against 2 in 100 at six or more, a
//! twenty-two-fold spread inside one dataset (`spec/parameter_prepass_ssr.md` §4). **The
//! fits are filed per `(read group, stratum)`**, the read group being the other half of the
//! key rather than part of the stratum, because slippage is a property of the chemistry
//! rather than of the individual.
//!
//! **Nothing here calls a genotype first.** The genotype is summed over rather than
//! chosen, which is what separates this from production's stutter pre-pass: that one pools
//! reads from loci that passed a confident-genotype gate, and measured against HG002's
//! assembly truth it reports gains as *more* common than losses where the truth is 3.4
//! times the other way (`spec/parameter_prepass.md` §2.2).
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_ssr.md` (what is fitted and why),
//! `doc/devel/ng/arch/parameter_prepass_ssr.md` (types and interfaces), on the shared
//! framing of `doc/devel/ng/spec/parameter_prepass.md`.
//!
//! Five files, split so that the shaping of data and the mathematics on it never live
//! together:
//!
//! - this one — the vocabulary a stratum is keyed on, the widths the fit works inside, and
//!   (from Milestone A5) what a fit emits;
//! - [`offset_bucket`] — where a read's tract length is recorded, and the only place that
//!   can build one of those buckets;
//! - [`locus_offsets`] — one STR locus reduced to one table entry;
//! - [`stratum_table`] — the sparse table of locus shapes, per stratum;
//! - [`slippage`] — the noise model and the search that fits it.
//!
//! **No trait over the accumulator.** Nothing generic drives it and the walk knows which
//! object it is filling; `fitting/` is the one genuine swappable seam in step 4.

pub mod locus_offsets;
pub mod offset_bucket;
pub mod slippage;
pub mod stratum_table;

use std::collections::BTreeMap;
use std::fmt;

use smallvec::SmallVec;

use crate::ng::parameter_estimation::Estimate;
use crate::ng::parameter_estimation::ssr::slippage::SlippageModel;
use crate::ng::types::{DomainError, ErrorRate, ReadGroupId, SsrPeriod};

pub use offset_bucket::{OFFSET_BUCKETS, OFFSET_HALF_RANGE, OffsetBucket, bucket_of};

/// How many whole motif copies a tract holds.
///
/// **The reference tract's count, never the sample's.** It is a pure function of the
/// reference, which is what makes every sample file a locus under the same stratum and so
/// lets a cohort compare one sample's stutter with another's
/// (`arch/parameter_prepass_ssr.md` §2.1). A sample whose alleles differ from the
/// reference does not move between strata; its reads land at an offset instead
/// ([`WholeRepeatOffset`]).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct RepeatCount(pub u32);

impl RepeatCount {
    #[inline]
    #[must_use]
    pub fn get(self) -> u32 {
        self.0
    }
}

impl fmt::Display for RepeatCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// One group of loci that gets its own fitted stutter parameters: a motif period and a
/// reference repeat count (`spec/parameter_prepass_ssr.md` §4).
///
/// **Ordered by `(period, repeats)`**, so walking a map of strata visits each period's
/// repeat counts in ascending order — which is what the monotonicity rule needs: slippage
/// genuinely rises with repeat count, so a fitted sequence that dips in the middle is
/// reporting noise in one stratum rather than a fact about repeats, and finding that means
/// comparing each stratum with the one before it (§4.3).
///
/// **The reference's own catalog already speaks this concept and spells it differently**:
/// [`repeat_catalog::strata`](crate::ng::repeat_catalog::strata) keys its counts and its
/// per-stratum sample on a raw `(u8, u64)` pair, so a driver that asks the catalog how many
/// loci a stratum holds converts at that seam. The two orderings agree — both are period
/// first, then repeat count — and nothing yet observes that they do, because nothing
/// converts between them until Milestone C.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Stratum {
    pub period: SsrPeriod,
    pub repeats: RepeatCount,
}

impl Stratum {
    #[must_use]
    pub fn new(period: SsrPeriod, repeats: RepeatCount) -> Self {
        Self { period, repeats }
    }
}

impl fmt::Display for Stratum {
    /// "period 2, 6 repeats" — the words the emitted summary and every error message use,
    /// so neither has to spell the pair out again.
    ///
    /// Destructured rather than read field by field, so that a third field added to
    /// `Stratum` is a compile error here rather than a rendering that silently names two
    /// thirds of the key.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self { period, repeats } = self;
        write!(f, "period {period}, {repeats} repeats")
    }
}

/// How far a read's tract sits from the **reference** tract length, in whole motif copies.
/// Negative is shorter than the reference, positive is longer.
///
/// **The origin is the reference tract's length and not each locus's most common observed
/// length**, which is what an earlier draft of the design used. Measured, the modal origin
/// returns a slippage level 50% to 408% high and destroys the direction asymmetry the model
/// exists to carry — a split of 0.48 where the truth is 0.17, a 1.1-fold imbalance where
/// the truth is 4.9-fold (`spec/parameter_prepass_ssr.md` §4.1). The reason is that
/// centring on the mode makes the origin a function of the reads, so a fit treating it as a
/// property of the locus is answering a question about a quantity that moved.
///
/// Unconstrained: any offset is a legal thing to observe. What is bounded is the bucket it
/// is *recorded* in ([`OffsetBucket`]) and the offsets a fitted allele may sit at
/// ([`ALLELE_OFFSET_LIMIT`]).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct WholeRepeatOffset(pub i8);

impl WholeRepeatOffset {
    #[inline]
    #[must_use]
    pub fn get(self) -> i8 {
        self.0
    }
}

impl fmt::Display for WholeRepeatOffset {
    /// Signed, always — `+2` and `-2` are different observations and a bare `2` in a
    /// report would not say which.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:+}", self.0)
    }
}

/// How far from the reference tract length the fit may place an allele — the support of a
/// stratum's genotype frequencies, and **a different width from the recorded range, wider
/// than it**.
///
/// It is what lets a saturating end bucket be explained by a distant *allele* rather than
/// by a far *slip*, which is why this width and not [`OFFSET_HALF_RANGE`] decides the
/// answer. Too narrow and the fit charges a real long allele to slippage; too wide and a
/// thin stratum is fitting frequencies for lengths no locus carries, at `A(A+1)/2` numbers
/// for `A` lengths.
///
/// **Six comes from the measured distribution of allele lengths and is a threshold to
/// clear rather than a number to tune.** A locus's most common observed length against the
/// reference, on HG002 at 300×: 88.9% sit exactly at the reference, ±4 holds 99%, ±12 holds
/// 99.9%, ±19 holds 99.99%; tomato is tighter, 95.7% at zero with ±1 holding 99%
/// (`spec/parameter_prepass_ssr.md` §8, item 1). What a locus outside the support costs is
/// nothing and then everything: leaving 2.5% of loci out costs 0.1% of the slippage level,
/// 7.9% costs 2.5%, and **19.3% costs +499% with the direction asymmetry destroyed**. Six
/// leaves about one human locus in 200 outside — a fifth of the way to the row that is
/// already free.
///
/// **It clips at the low end, so the number of allele lengths is a per-stratum quantity**:
/// an allele cannot be shorter than nothing, so a stratum at 4 repeats reaches only −4.
/// [`allele_support`] is where that happens.
pub const ALLELE_OFFSET_LIMIT: i8 = 6;

/// The relation the design turns on, held by the compiler rather than by prose: the fit may
/// place an allele **further** from the reference than an entry records offsets. That is
/// what lets a saturating end bucket be attributed to a distant allele instead of to a far
/// slip; invert it and every locus whose allele sits past the recorded range has its reads
/// explained the only way left, as slippage.
const _: () = assert!(
    ALLELE_OFFSET_LIMIT > OFFSET_HALF_RANGE,
    "the fit must be able to place an allele beyond the offsets an entry records"
);

/// Reads entered from one locus. A deeper locus is entered from a **random subsample** of
/// its reads down to this, seeded from the locus's position so a region-sharded walk and a
/// single-threaded one keep the same reads and merging stays exact.
///
/// **Not the memory knob it looks like.** Measured on HG002 at 300× over the GIAB
/// tandem-repeat set, the uncapped table is 0.43 entries a locus — deep data deduplicates,
/// because most loci at a clean tract are "every read at the reference length" and what
/// separates two entries is mostly their depth (`spec/parameter_prepass_ssr.md` §4.1).
///
/// **Nor is it a correctness limit**: the scoring rule is exactly unbiased at every depth
/// tried, to 45 reads a locus. What it does decide is the width of an entry's counters — a
/// `u8` bucket count wraps silently above 255 — so the cap and that width are one decision.
/// Twelve is a low starting value whose only cost is the precision of the reads it drops.
///
/// Consumed from Milestone C, where a locus first becomes an entry.
pub const MAX_LOCUS_READS: u32 = 12;

/// Above this share of the reads that differ from the reference tract length, a stratum is
/// one this noise model does not describe: what moved the reads is ordinary indel rather
/// than repeat slippage, and a slippage rate fitted there is mostly mis-modelled indel
/// however many loci stood behind it (`spec/parameter_prepass_ssr.md` §5).
///
/// **One in ten, and the bands it separates are 0.9% against 33.8% and 58.5%** — a factor
/// of three either way, so a stratum crossing it is unambiguous. *Soft*: three bands of one
/// dataset, and the number moves if the per-stratum distribution turns out to be continuous
/// rather than the two clumps that table suggests.
///
/// Consumed from Milestone B, where the table starts counting the two denominators.
pub const GUARD_SHARE_LIMIT: f64 = 0.10;

/// The allele lengths a stratum's fit may place mass on, as offsets from the reference
/// tract length, in ascending order.
///
/// `±ALLELE_OFFSET_LIMIT` around the reference, **clipped at the low end because an allele
/// cannot be shorter than nothing**: a tract of 4 reference copies reaches −4, so it has 11
/// lengths where a tract of 6 or more has the full 13. The genotype count follows —
/// `A(A+1)/2` unordered pairs, so 66 against 91 — which is why the fit is handed a stratum's
/// support rather than assuming one (`arch/parameter_prepass_ssr.md` §2.1).
#[must_use]
pub fn allele_support(reference_repeats: RepeatCount) -> Vec<WholeRepeatOffset> {
    // How far below the reference an allele could reach at this stratum, before the fit's
    // own limit applies: a tract of `n` copies can lose at most `n` of them.
    let reach_down = i8::try_from(reference_repeats.get()).unwrap_or(i8::MAX);
    // Parenthesised deliberately: `-ALLELE_OFFSET_LIMIT.min(reach_down)` parses the same way
    // but reads as `(-6).min(reach_down)`, which at 3 reference repeats gives −6 — a support
    // reaching three copies below an empty tract, silently.
    let lowest = -(reach_down.min(ALLELE_OFFSET_LIMIT));
    (lowest..=ALLELE_OFFSET_LIMIT)
        .map(WholeRepeatOffset)
        .collect()
}

// ---------------------------------------------------------------------
// What a fit returns, and what a person reads instead of it.
//
// This path runs **one fit per (read group, stratum)** — several hundred a
// sample, against the SNP/indel path's four in total. That is the whole
// reason the output has two halves: a per-stratum record, because a fit
// that looks wrong has to be traceable; and a summary per read group,
// because several hundred records are a file nobody opens, and a flag
// nobody reads is how a badly-fitted parameter reaches a caller
// (`spec/parameter_prepass_ssr.md` §4.4).
// ---------------------------------------------------------------------

/// One allele pair the fit placed mass on, and how much: a genotype and its frequency.
///
/// The two alleles are **unordered** — a locus carrying a reference-length allele and one a
/// repeat short is the same genotype either way round — and the convention that keeps one
/// spelling per genotype is that `shorter <= longer`.
///
/// **Named fields rather than the tuple the architecture sketches**, for the reason
/// `fitting/`'s own `WeightedCell` replaced a three-member tuple: nothing in
/// `(WholeRepeatOffset, WholeRepeatOffset, f64)` says which member is which, and two of the
/// three have the same type.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct AllelePairFrequency {
    /// The shorter of the two alleles, as an offset from the reference tract length.
    shorter: WholeRepeatOffset,
    /// The longer of the two. Equal to `shorter` at a homozygous genotype.
    longer: WholeRepeatOffset,
    /// The share of the stratum's loci the fit gives this pair.
    frequency: f64,
}

impl AllelePairFrequency {
    /// The only constructor, and it **establishes the unordered convention rather than
    /// trusting it**: the two alleles are sorted, so one genotype has one spelling however
    /// the caller happened to hold it. Without that, a fit emitting `(-1, 0)` and one
    /// emitting `(0, -1)` would describe the same genotype twice and the frequencies would
    /// not sum to one.
    #[must_use]
    pub fn new(one: WholeRepeatOffset, other: WholeRepeatOffset, frequency: f64) -> Self {
        let (shorter, longer) = if one <= other {
            (one, other)
        } else {
            (other, one)
        };
        Self {
            shorter,
            longer,
            frequency,
        }
    }

    /// The shorter of the two alleles.
    #[inline]
    #[must_use]
    pub fn shorter(self) -> WholeRepeatOffset {
        self.shorter
    }

    /// The longer of the two — equal to [`Self::shorter`] at a homozygous genotype.
    #[inline]
    #[must_use]
    pub fn longer(self) -> WholeRepeatOffset {
        self.longer
    }

    /// The share of the stratum's loci the fit gives this pair.
    #[inline]
    #[must_use]
    pub fn frequency(self) -> f64 {
        self.frequency
    }

    /// Whether both alleles are the same length — a homozygous genotype.
    #[inline]
    #[must_use]
    pub fn is_homozygous(self) -> bool {
        self.shorter == self.longer
    }
}

/// Where one of the search's starting points began, where it ended, and what it scored.
///
/// **Emitted for every start, not just the winner**, and that is the diagnostic the
/// multi-start search exists to produce: an answer with no spread beside it is
/// indistinguishable from a search that never looked. The generic path's inbreeding fit
/// produced exactly that failure — a confident zero on a genome 26% covered by runs — from
/// starts that agreed because they shared one guess at a nuisance axis
/// (`spec/parameter_prepass_ssr.md` §4.2).
///
/// **The reverse trap is just as real and is why the score travels too:** a deterministic
/// optimiser returns the same point from every start wherever the objective is flat, so four
/// starts agreeing to four decimal places is also what a search that never moved produces.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SlippageStart {
    /// Where this start began.
    ///
    /// Named `from` rather than `start`, so that a reader of `start.start` is never asked to
    /// work out which of the two is the noun — the sibling `StartOutcome` in `fitting/` uses
    /// the same word.
    pub from: SlippageModel,
    /// Where the climb from `from` ended.
    pub reached: SlippageModel,
    /// What `reached` scored, as a natural logarithm.
    pub log_likelihood: f64,
}

/// What one stratum's fit returned, for one read group.
///
/// **Two provenance lists and not one**, which is this type's least obvious feature and the
/// one a consumer must not collapse. The level and the two shares are measured from different
/// populations, and at the bottom of the repeat range those populations differ by four orders
/// of magnitude: a stratum of 100,000 loci at five reads each holds **half a million reads**,
/// of which — at the 0.091% level measured below four repeats — about **455** slipped, about
/// 77 of those gained a repeat, and about **5** of those gained two. So the level is a
/// proportion over 500,000 reads and the fall-off's gaining arm stands on 5, which is why the
/// same stratum measures the first to about 5% of itself and the second to about 45%. The
/// honest answer is to keep the level it measured and borrow the two shares it did not
/// (`spec/parameter_prepass_ssr.md` §4.5).
#[derive(Clone, PartialEq, Debug)]
pub struct StratumFit {
    /// Which stratum this is — the period and the reference repeat count.
    pub stratum: Stratum,
    /// The three slippage parameters, with where they came from and how much data stood
    /// behind them.
    pub slippage: Estimate<SlippageModel>,
    /// The per-base substitution rate inside this stratum's tracts. **Fitted separately from
    /// the generic path's rate and never tied to it**: each is the error parameter of its own
    /// model, absorbing what that model cannot otherwise explain, so forcing one number to
    /// carry both would make each wrong in a way neither could report. What separate fits buy
    /// is that the two can be *compared* — and where a stratum barely slips they must agree,
    /// because there the two noise models describe the same thing
    /// (`spec/parameter_prepass_ssr.md` §1.1, §4.5).
    pub substitution: Estimate<ErrorRate>,
    /// The allele-length distribution the fit weighed the genotype against, one entry per
    /// unordered pair. Emitted because it is what a reader needs when a slippage rate looks
    /// wrong, and because the cohort gather has a use for it.
    pub genotypes: Vec<AllelePairFrequency>,
    /// Of the reads that differed from the **reference tract length**, the share that did so
    /// by something other than a whole number of copies. Above [`GUARD_SHARE_LIMIT`] this
    /// model does not describe the stratum, and its fitted slippage is mostly mis-modelled
    /// ordinary indel however many loci stood behind it.
    ///
    /// **The denominator is the reference and the model means the allele**, which the
    /// accumulator cannot know. A real non-reference allele contributes whole-repeat
    /// differences to this denominator and none to its numerator, so the reported share is
    /// diluted relative to the model's, never inflated — a stratum that crosses the limit on
    /// this number has crossed it on the model's too.
    pub not_whole_repeat_share: f64,
    /// The share of this stratum's loci whose shape the **fitted** model calls very unlikely,
    /// per read (`spec/parameter_prepass_ssr.md` §4.6).
    ///
    /// **A different question from the guard share above it**: that one asks whether this
    /// model is right for these tracts, this one asks whether these are all the tracts the
    /// model was told they were. A duplication the reference does not carry moves its reads by
    /// whole copies, so the guard share never sees it, and one locus in a thousand behaving
    /// strangely is invisible in any stratum-wide ratio.
    ///
    /// **Reported, never acted on.** Dropping the loci that score badly is threshold-then-count
    /// — the bias this whole step exists to remove — and it would take real long alleles before
    /// it took artefacts, since both sit where the model has least mass.
    pub unexplained_locus_share: f64,
    /// Every starting point the search tried, best-scoring first.
    pub starts_tried: SmallVec<[SlippageStart; 4]>,
    /// Which strata **the level's** loci came from: this one where it was fitted in place, its
    /// neighbours where it borrowed, both where a merge fired.
    pub fitted_over: SmallVec<[Stratum; 2]>,
    /// Which strata **the direction share and the fall-off** came from, which is not always
    /// the same answer. A stratum of 100,000 loci clears [`MIN_LOCI_TO_FIT`] a hundred times
    /// over and still puts about five reads behind the fall-off's gaining arm, at the level
    /// measured below four repeats.
    pub shares_fitted_over: SmallVec<[Stratum; 2]>,
    /// Reads that showed a length other than the reference tract's. **The count that decides
    /// whether the two shares are measurable at all**, and the one a consumer needs to tell a
    /// level of 0.0003 standing on 4 slipped reads from one standing on 4,000 — which nothing
    /// downstream could otherwise do.
    pub slipped_reads: u64,
}

/// Below this per-read log-likelihood under the fitted model, a locus counts towards
/// [`StratumFit::unexplained_locus_share`].
///
/// **Per read, not per locus, or the statistic measures depth**: a shape's likelihood is a
/// product over its reads, so a locus at twelve reads scores about four times lower than one
/// at three whatever the model thinks of either, and a floor on the total would flag the
/// deepest loci of every stratum by arithmetic.
///
/// **Soft, and to be set by measurement rather than by argument.** No run has yet asked where
/// real loci sit on this scale; the value here flags a tail rather than a body, and Milestone
/// H's contamination experiment is what fixes it, by showing where the two populations
/// separate. A diagnostic that fires on every stratum is as useless as one that fires on none.
pub const UNEXPLAINED_SHAPE_LN_LIMIT: f64 = -3.0;

/// Fewest loci a stratum needs before its slippage is fitted rather than borrowed from a
/// neighbouring repeat count at the same period.
///
/// **Soft, and — unlike its three companions below — carrying no derivation, here or in the
/// design.** The other floors are computed from the precision wanted; this one is a round
/// number, and what it trades is a fit on thin evidence against a borrow that costs 15 to 25%
/// of the level per repeat count borrowed across. Both sides of that trade are recorded and
/// neither has been measured against the other, so the value is a starting point rather than a
/// conclusion.
pub const MIN_LOCI_TO_FIT: u64 = 1_000;

/// Fewest **slipped** reads a stratum needs before its direction share and its fall-off are
/// its own rather than a neighbour's.
///
/// **The level has no such floor**, and the asymmetry is the point: the level is a proportion
/// over every read, so a stratum of 100,000 loci at five reads each measures it to about 5% of
/// itself, while the fall-off — measured only by the reads that actually moved twice — is
/// measured to about 45%.
///
/// **Derived from the precision wanted rather than chosen.** At the values the design
/// measures — a gain share of 0.17 and a fall-off of 0.065 — holding the share to 6% of itself
/// takes about 1,400 slipped reads and holding the fall-off to the same takes about 4,000, so
/// the fall-off binds. *Soft*, and **expected to be missed by every stratum at the bottom of
/// the repeat range**: at a level of 0.091% it takes about 880,000 loci at five reads each,
/// against tomato's 1.73 million STR loci in total. That is the rule working — the alternative
/// is a share fitted on five reads and reported as though it were measured.
pub const MIN_SLIPPED_READS_TO_FIT_SHARES: u64 = 4_000;

/// How far two starting points may land apart, as a ratio in the slippage level, before the
/// fit is reported as not having found an answer.
///
/// **One quarter-Phred — 6% — borrowed from the generic path's error-rate ladder spacing**,
/// which that design argues is the finest difference a caller can feel. This path searches
/// rather than scanning, so it has no rungs of its own to measure a disagreement against.
pub const START_AGREEMENT_LIMIT: f64 = 1.06;

/// One read group's fits, summarised — **the part a person reads**.
///
/// Several hundred per-stratum records are a file nobody opens, so every field here answers a
/// question a reader would otherwise have to grep those records for
/// (`spec/parameter_prepass_ssr.md` §4.4).
#[derive(Clone, PartialEq, Debug, Default)]
pub struct StratumFitSummary {
    /// Strata whose slippage was fitted from their own loci.
    pub strata_fitted_here: u32,
    /// Strata that took the whole slippage model from a neighbour, being below
    /// [`MIN_LOCI_TO_FIT`].
    pub strata_borrowed: u32,
    /// Strata that kept the level they measured and borrowed only the direction share and the
    /// fall-off, being below [`MIN_SLIPPED_READS_TO_FIT_SHARES`].
    ///
    /// **Counted apart from [`Self::strata_borrowed`] because it is a different claim**: these
    /// report a level of their own.
    pub strata_with_borrowed_shares: u32,
    /// The merged sets, named. A merge is a claim about two strata at once, so the summary
    /// says which — it changes both estimates, by about a quarter of the level at a 1.5-fold
    /// difference between neighbours and up to 141% at four-fold.
    pub strata_merged: Vec<SmallVec<[Stratum; 2]>>,
    /// The substitution rate this read group's least-slippery strata returned, and how many
    /// loci stood behind it.
    ///
    /// **Emitted so that step 4's own surface can put it beside the generic path's rate for
    /// the same read group**, which is where the two must agree to a quarter-Phred: where a
    /// stratum barely slips, this path's noise model *is* the generic one. This unit never
    /// sees the generic half, so it emits the operand and makes no comparison.
    ///
    /// **An `Estimate` rather than a rate and a count**, because the comparison is only
    /// meaningful against a rate that was *fitted*: a borrowed one describes a neighbouring
    /// stratum's tracts, and comparing that with the generic path's would test nothing. The
    /// provenance is what says which happened.
    pub low_slippage_substitution: Option<Estimate<ErrorRate>>,
    /// Fits whose starting points landed further apart than [`START_AGREEMENT_LIMIT`] in the
    /// level — the diagnostic the several starts exist to produce.
    pub strata_with_disagreeing_starts: u32,
    /// The worst of those, and by what ratio.
    pub worst_start_disagreement: Option<(Stratum, f64)>,
    /// Strata above [`GUARD_SHARE_LIMIT`] — the ones this noise model does not describe.
    pub strata_above_guard_limit: u32,
    /// The worst of those, and its share.
    pub worst_guard_share: Option<(Stratum, f64)>,
    /// The stratum holding the largest share of loci its own fitted model cannot explain, and
    /// that share (`spec/parameter_prepass_ssr.md` §4.6).
    pub worst_unexplained_locus_share: Option<(Stratum, f64)>,
    /// How many loci stood behind the **thinnest** fit this read group produced.
    ///
    /// Two named fields rather than the `(u64, u64)` the architecture sketches, for the reason
    /// [`AllelePairFrequency`] is not a tuple either: nothing in a pair of `u64`s says which
    /// end is which, and a transposition would report the thinnest fit as the thickest with
    /// every type still matching.
    pub loci_behind_thinnest_fit: u64,
    /// And behind the **thickest**.
    pub loci_behind_thickest_fit: u64,
}

/// Everything the STR path estimates for one sample.
///
/// **No `Default`, deliberately, and the sibling `GenericSampleParameters` has none either.**
/// A default-constructed value is an empty set of fits, which is indistinguishable from a
/// sample whose every stratum was fitted and found nothing — so `.unwrap_or_default()` would
/// turn a failure into a report. A sample that could not be fitted raises
/// [`SsrEstimationError`].
#[derive(Clone, PartialEq, Debug)]
pub struct SsrSampleParameters {
    /// One fit per read group and stratum. Traceable rather than read: a fit that looks wrong
    /// has to be findable, and nothing downstream is expected to walk this.
    pub by_stratum: BTreeMap<(ReadGroupId, Stratum), StratumFit>,
    /// **What a person reads instead.**
    pub summary: BTreeMap<ReadGroupId, StratumFitSummary>,
}

/// Everything the accumulator did to a locus other than enter it as it arrived.
///
/// **Every field is a plain sum, so shards merge; and every field is reported**, because the
/// alternative is a set of parameters that quietly describes a different population of reads
/// from the one the caller will see.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct SsrAccumulationCounts {
    /// Loci this unit subsampled down to [`MAX_LOCUS_READS`]. A run where most loci appear
    /// here is one where the cap, and not the data, is setting the depth.
    pub loci_subsampled_to_cap: u64,
    /// Loci whose reads the **generator's** own reservoir had already subsampled. Entered
    /// anyway, at the depth observed — that reservoir is a random subsample, so a locus it
    /// fired on is a locus observed at a lower depth rather than a locus to skip.
    pub loci_with_upstream_subsample: u64,
    /// Reads that covered a tract and witnessed nothing. Not part of any depth.
    pub reads_without_observation: u64,
    /// Reads whose witness was partial, so their length is a **lower bound**. Excluded, and
    /// deliberately: scoring one as a length reads as a read that lost repeats, which is a
    /// direct bias in the one parameter the design exists to protect. A large share here means
    /// the reads are short against these tracts.
    pub reads_with_partial_witness: u64,
    /// Loci whose reference tract is not a whole number of motif copies, so no stratum holds
    /// them. **Should read near zero** — the classification that admits a locus delimits on
    /// whole copies — and a large count is a bug report against it rather than something this
    /// unit absorbs.
    pub loci_without_whole_repeat_reference: u64,
}

impl SsrAccumulationCounts {
    /// Add another shard's counts to these. Associative and exact, so a region-sharded walk
    /// and a single-threaded one report the same adjustments.
    ///
    /// Destructured exhaustively, so a counter added later is a compile error here rather
    /// than a field that silently stops merging.
    pub fn merge(&mut self, other: &Self) {
        let Self {
            loci_subsampled_to_cap,
            loci_with_upstream_subsample,
            reads_without_observation,
            reads_with_partial_witness,
            loci_without_whole_repeat_reference,
        } = other;
        self.loci_subsampled_to_cap += loci_subsampled_to_cap;
        self.loci_with_upstream_subsample += loci_with_upstream_subsample;
        self.reads_without_observation += reads_without_observation;
        self.reads_with_partial_witness += reads_with_partial_witness;
        self.loci_without_whole_repeat_reference += loci_without_whole_repeat_reference;
    }
}

/// What went wrong while estimating a sample's STR parameters.
///
/// **This path's own enum rather than variants added to the generic path's**, because the two
/// units fail differently: a stratum too thin to fit has neighbours to borrow from, and a
/// sample's heterozygosity does not. Step 4's own surface wraps both for a caller that drove
/// the whole step.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SsrEstimationError {
    /// Every stratum at this period was below [`MIN_LOCI_TO_FIT`], so there is no neighbour to
    /// borrow from.
    ///
    /// **Deliberately has no default to fall back on.** A slippage level spans twenty-two-fold
    /// across repeat counts within one dataset, so any constant would be wrong for most strata
    /// — and wrong in the direction that reads as a measurement.
    ///
    /// **It names how far short the period fell**, because that is what separates the two
    /// remedies: a period whose thickest stratum held 812 loci wants a run over more of the
    /// genome, and one that held 3 wants dropping.
    #[error(
        "sample {sample}, read group {read_group}: period {period} has no stratum with \
         {MIN_LOCI_TO_FIT} loci to fit or borrow from — its thickest holds {thickest_loci}; \
         drop the period or widen the run"
    )]
    NoFittableStratumAtPeriod {
        sample: String,
        /// **Which library**, because the fits are per read group and a sample can carry
        /// several: naming only the sample points at four fits and says nothing about which.
        read_group: ReadGroupId,
        period: SsrPeriod,
        /// How many loci the period's thickest stratum held.
        thickest_loci: u64,
    },

    /// The search reached materially different answers from different starting points, so what
    /// it returned is where it stopped rather than what the data says.
    ///
    /// **Not "too little data"** — that is a borrow, and it is recorded as one.
    #[error(
        "sample {sample}, read group {read_group}: stratum ({stratum}) reached slippage levels \
         spanning {spread:.1}x across {starts} starting points, more than the \
         {START_AGREEMENT_LIMIT} that leaves the level identified — this is a search that did \
         not settle, not a slippage level"
    )]
    SlippageNotIdentified {
        sample: String,
        /// Which library's fit did not settle.
        read_group: ReadGroupId,
        stratum: Stratum,
        /// The highest level any start reached, divided by the lowest.
        spread: f64,
        /// How many starts were tried.
        starts: usize,
    },

    /// A constrained scalar rejected its value while a named fit was running.
    ///
    /// **Not `#[from]`, and not transparent**, for the reason the generic path's twin gives:
    /// the inner error names the quantity and the offending value, and what it cannot know is
    /// whose data and which fit produced it — on a cohort run of hundreds of samples, that is
    /// the half that locates the fault. A `#[from]` conversion would let `?` mint this variant
    /// at every constructor, so the sample and the fit would have to be *remembered* rather
    /// than demanded.
    #[error(
        "sample {sample}, read group {read_group}, {stratum}: {fit} rejected a value — {source}"
    )]
    Domain {
        sample: String,
        /// Which library's fit raised it. Several hundred fits run per read group, so the
        /// sample alone locates nothing.
        read_group: ReadGroupId,
        /// And which stratum inside it.
        stratum: Stratum,
        /// Which fit was running, in the words the emitted summary uses — "the slippage
        /// search", "the substitution rate".
        fit: &'static str,
        source: DomainError,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn period(bases: u8) -> SsrPeriod {
        SsrPeriod::try_new(usize::from(bases)).expect("a period inside the STR scope")
    }

    fn stratum(bases: u8, repeats: u32) -> Stratum {
        Stratum::new(period(bases), RepeatCount(repeats))
    }

    /// The monotonicity rule walks a period's strata in repeat-count order and compares each
    /// fit with the one before it, so this ordering is what puts neighbours next to each
    /// other. Sorting by repeat count first would interleave the periods and compare a
    /// dinucleotide's fit against a hexamer's.
    #[test]
    fn strata_sort_by_period_and_then_by_repeat_count() {
        let mut strata = vec![
            stratum(2, 12),
            stratum(1, 20),
            stratum(2, 6),
            stratum(1, 8),
            stratum(6, 4),
        ];
        strata.sort();

        assert_eq!(
            strata,
            vec![
                stratum(1, 8),
                stratum(1, 20),
                stratum(2, 6),
                stratum(2, 12),
                stratum(6, 4),
            ]
        );
    }

    /// A stratum names itself in the words the summary and the error messages use, so
    /// neither has to spell the pair out again — and a reader of a log sees the pair rather
    /// than a tuple.
    #[test]
    fn a_stratum_names_its_period_and_its_repeat_count() {
        assert_eq!(stratum(2, 6).to_string(), "period 2, 6 repeats");
    }

    /// **A zero period is unrepresentable here rather than rejected here**: `Stratum::new`
    /// takes an already-checked [`SsrPeriod`], and the rejection itself is
    /// `ng::types::tests::ssr_period_accepts_exactly_the_str_scope`. What this module owns is
    /// that a stratum carries back the pair it was built from — a swap between the two would
    /// file every locus under a stratum no model was fitted at.
    #[test]
    fn a_stratum_carries_back_the_period_and_repeat_count_it_was_built_from() {
        let built = stratum(2, 6);
        assert_eq!(built.period.get(), 2);
        assert_eq!(built.repeats.get(), 6);
    }

    /// An offset is read for its **direction** first — a report column of them is how one
    /// tells a stratum losing repeats from one gaining them — so the sign is explicit on
    /// gains and at the origin, not only on losses.
    #[test]
    fn an_offset_renders_with_its_sign_in_both_directions() {
        assert_eq!(WholeRepeatOffset(2).to_string(), "+2");
        assert_eq!(WholeRepeatOffset(-2).to_string(), "-2");
        assert_eq!(WholeRepeatOffset(0).to_string(), "+0");
    }

    /// **The support is every offset from the clip to the limit, with nothing missing and
    /// nothing repeated.** Asserting its length and its two ends is not enough: a support
    /// holding a duplicate and a hole has the right cardinality, so the `A(A+1)/2` genotype
    /// count still comes out right while the fit never considers one allele length and
    /// double-counts another.
    ///
    /// The five strata are the clip's whole range: none at 0 copies, biting at 3 and at 5,
    /// exactly released at 6, and far above it at 20.
    #[test]
    fn the_allele_support_is_every_offset_from_the_clip_to_the_limit() {
        for (repeats, lowest) in [(0u32, 0i8), (3, -3), (5, -5), (6, -6), (20, -6)] {
            assert_eq!(
                allele_support(RepeatCount(repeats)),
                (lowest..=ALLELE_OFFSET_LIMIT)
                    .map(WholeRepeatOffset)
                    .collect::<Vec<_>>(),
                "at {repeats} reference repeats"
            );
        }
    }

    /// The lengths behind the genotype count the design reasons about: 10, 13 and 13 allele
    /// lengths give 55, 91 and 91 unordered pairs.
    #[test]
    fn the_allele_support_clips_below_the_reference_but_not_above() {
        let short = allele_support(RepeatCount(3));
        assert_eq!(short.len(), 10);
        assert_eq!(short.first(), Some(&WholeRepeatOffset(-3)));
        assert_eq!(short.last(), Some(&WholeRepeatOffset(ALLELE_OFFSET_LIMIT)));

        for repeats in [6, 20] {
            let full = allele_support(RepeatCount(repeats));
            assert_eq!(
                full.len(),
                13,
                "at {repeats} repeats the support is the full ±{ALLELE_OFFSET_LIMIT}"
            );
            assert_eq!(full.first(), Some(&WholeRepeatOffset(-ALLELE_OFFSET_LIMIT)));
        }
    }

    /// A tract at the copy floors this caller routes on never sees this, but the arithmetic
    /// that clips the support must not panic on a repeat count that does not fit in the
    /// offset type — a satellite of 20,000 copies is one `u32` value away from any tract.
    #[test]
    fn the_support_survives_a_repeat_count_far_larger_than_an_offset() {
        let huge = allele_support(RepeatCount(u32::MAX));
        assert_eq!(huge.len(), 13);
        assert_eq!(huge.first(), Some(&WholeRepeatOffset(-ALLELE_OFFSET_LIMIT)));
    }

    /// A zero-copy tract is not something region typing emits, but the support must still be
    /// a set the fit can walk: at zero the only allele length that exists is the reference's
    /// own, and every longer one.
    #[test]
    fn a_tract_of_no_repeats_can_only_hold_alleles_at_or_above_the_reference() {
        let none = allele_support(RepeatCount(0));
        assert_eq!(none.first(), Some(&WholeRepeatOffset(0)));
        assert_eq!(
            none.len(),
            usize::try_from(ALLELE_OFFSET_LIMIT).unwrap() + 1
        );
    }

    fn a_slippage(level: f64) -> SlippageModel {
        SlippageModel::try_new(level, 0.17, 0.065).expect("three probabilities")
    }

    fn fitted_here<T>(value: T, observations: u64) -> Estimate<T> {
        Estimate {
            value,
            provenance: crate::ng::parameter_estimation::Provenance::FittedHere,
            observations,
        }
    }

    /// **The two provenance lists are what a consumer reads to tell three different claims
    /// apart**, and a fit with one list could not: a stratum that measured everything itself,
    /// one that measured nothing, and one that measured its level well and its two shares not
    /// at all. The third is the common case at the bottom of the repeat range, where reaching
    /// 4,000 slipped reads would take about 880,000 loci.
    ///
    /// **A shape test, and it says so rather than pretending otherwise.** Nothing produces a
    /// `StratumFit` until Milestone E, so this asserts over a value the test itself built: what
    /// it can fail on is the *type* losing its ability to carry two different answers — a
    /// merged field, a shared list — and what it cannot fail on is a producer filling both
    /// lists alike. That assertion lands with the producer.
    #[test]
    fn a_fit_can_say_separately_where_its_level_and_its_two_shares_came_from() {
        let here = stratum(2, 5);
        let neighbour = stratum(2, 6);

        let level_kept_shares_borrowed = StratumFit {
            stratum: here,
            slippage: fitted_here(a_slippage(0.00091), 500_000),
            substitution: fitted_here(ErrorRate::try_new(0.003).unwrap(), 500_000),
            genotypes: vec![AllelePairFrequency::new(
                WholeRepeatOffset(0),
                WholeRepeatOffset(0),
                1.0,
            )],
            not_whole_repeat_share: 0.02,
            unexplained_locus_share: 0.0,
            starts_tried: SmallVec::new(),
            fitted_over: SmallVec::from_slice(&[here]),
            shares_fitted_over: SmallVec::from_slice(&[neighbour]),
            slipped_reads: 455,
        };

        assert_eq!(
            level_kept_shares_borrowed.fitted_over.as_slice(),
            &[here],
            "the level is this stratum's own"
        );
        assert_eq!(
            level_kept_shares_borrowed.shares_fitted_over.as_slice(),
            &[neighbour],
            "the two shares came from the neighbour that had the slipped reads"
        );
        assert_ne!(
            level_kept_shares_borrowed.fitted_over, level_kept_shares_borrowed.shares_fitted_over,
            "a consumer must be able to see that the two halves have different warrants"
        );
        assert_eq!(
            level_kept_shares_borrowed.slipped_reads, 455,
            "and the count that decides whether the shares were measurable travels with them"
        );
    }

    /// The summary has a **separate counter** for a stratum that borrowed everything and one
    /// that kept its own level, because they are different claims and a reader who cannot tell
    /// them apart cannot act on either.
    ///
    /// **A shape test.** No rule sums these until Milestone E, so what this can fail on is the
    /// two counters being merged into one; whether a producer counts a level-keeping stratum
    /// into `strata_borrowed` as well is a question for the producer's own test.
    #[test]
    fn the_summary_counts_a_whole_borrow_apart_from_a_borrowed_pair_of_shares() {
        let summary = StratumFitSummary {
            strata_fitted_here: 40,
            strata_borrowed: 3,
            strata_with_borrowed_shares: 12,
            ..StratumFitSummary::default()
        };

        assert_ne!(summary.strata_borrowed, summary.strata_with_borrowed_shares);
    }

    /// **A default summary reports nothing rather than something false**, which matters
    /// because `loci_behind_thickest_fit` is the kind of field a fold seeds and then maxes: a
    /// seed of `u64::MAX` would survive every real stratum and report the thinnest fit as the
    /// largest one. Zero on both ends, `None` on all four "worst" fields.
    #[test]
    fn a_default_summary_claims_no_fits_and_names_no_worst_stratum() {
        let empty = StratumFitSummary::default();

        assert_eq!(empty.strata_fitted_here, 0);
        assert_eq!(empty.strata_borrowed, 0);
        assert_eq!(empty.strata_with_borrowed_shares, 0);
        assert_eq!(empty.loci_behind_thinnest_fit, 0);
        assert_eq!(empty.loci_behind_thickest_fit, 0);
        assert!(empty.worst_start_disagreement.is_none());
        assert!(empty.worst_guard_share.is_none());
        assert!(empty.worst_unexplained_locus_share.is_none());
        assert!(empty.low_slippage_substitution.is_none());
        assert!(empty.strata_merged.is_empty());
    }

    /// **One genotype has one spelling**, whichever way round the fit held its two alleles —
    /// otherwise the same pair appears twice in a stratum's frequencies and they no longer sum
    /// to one.
    #[test]
    fn an_allele_pair_orders_its_two_alleles_however_it_was_given_them() {
        let one_way = AllelePairFrequency::new(WholeRepeatOffset(-1), WholeRepeatOffset(2), 0.25);
        let other_way = AllelePairFrequency::new(WholeRepeatOffset(2), WholeRepeatOffset(-1), 0.25);

        assert_eq!(one_way, other_way);
        assert_eq!(one_way.shorter(), WholeRepeatOffset(-1));
        assert_eq!(one_way.longer(), WholeRepeatOffset(2));
        assert_eq!(one_way.frequency(), 0.25);
        assert!(!one_way.is_homozygous());

        let homozygous = AllelePairFrequency::new(WholeRepeatOffset(3), WholeRepeatOffset(3), 0.5);
        assert!(homozygous.is_homozygous());
        assert_eq!(homozygous.shorter(), homozygous.longer());
    }

    /// **The guard share and the unexplained-locus share answer different questions**, so the
    /// summary names the worst stratum on each: one asks whether this model is right for these
    /// tracts, the other whether these are all the tracts the model was told they were. A
    /// stratum can be the worst on one and unremarkable on the other.
    #[test]
    fn the_summary_names_the_worst_stratum_on_each_of_its_two_diagnostics() {
        let summary = StratumFitSummary {
            worst_guard_share: Some((stratum(4, 3), 0.58)),
            worst_unexplained_locus_share: Some((stratum(1, 20), 0.004)),
            ..StratumFitSummary::default()
        };

        assert_ne!(
            summary.worst_guard_share.map(|(s, _)| s),
            summary.worst_unexplained_locus_share.map(|(s, _)| s)
        );
    }

    /// A sharded walk must report the same adjustments as a single-threaded one, so the
    /// counters are a plain sum in every field. **Every field**: the fixture gives all five
    /// different values, because a merge that dropped one would otherwise pass on a fixture
    /// that left it zero.
    #[test]
    fn accumulation_counts_merge_field_by_field() {
        let mut left = SsrAccumulationCounts {
            loci_subsampled_to_cap: 1,
            loci_with_upstream_subsample: 2,
            reads_without_observation: 3,
            reads_with_partial_witness: 4,
            loci_without_whole_repeat_reference: 5,
        };
        let right = SsrAccumulationCounts {
            loci_subsampled_to_cap: 10,
            loci_with_upstream_subsample: 20,
            reads_without_observation: 30,
            reads_with_partial_witness: 40,
            loci_without_whole_repeat_reference: 50,
        };

        left.merge(&right);

        assert_eq!(
            left,
            SsrAccumulationCounts {
                loci_subsampled_to_cap: 11,
                loci_with_upstream_subsample: 22,
                reads_without_observation: 33,
                reads_with_partial_witness: 44,
                loci_without_whole_repeat_reference: 55,
            }
        );
    }

    /// Merging is order-independent, which is what makes a region-sharded walk report the same
    /// adjustments however the genome was cut.
    #[test]
    fn accumulation_counts_merge_in_any_order_to_the_same_totals() {
        let shards = [
            SsrAccumulationCounts {
                loci_subsampled_to_cap: 7,
                ..SsrAccumulationCounts::default()
            },
            SsrAccumulationCounts {
                reads_with_partial_witness: 11,
                ..SsrAccumulationCounts::default()
            },
            SsrAccumulationCounts {
                loci_without_whole_repeat_reference: 2,
                ..SsrAccumulationCounts::default()
            },
        ];

        let mut forwards = SsrAccumulationCounts::default();
        for shard in &shards {
            forwards.merge(shard);
        }
        let mut backwards = SsrAccumulationCounts::default();
        for shard in shards.iter().rev() {
            backwards.merge(shard);
        }

        assert_eq!(forwards, backwards);
        assert_eq!(forwards.loci_subsampled_to_cap, 7);
        assert_eq!(forwards.reads_with_partial_witness, 11);
        assert_eq!(forwards.loci_without_whole_repeat_reference, 2);
    }

    /// **Every message names the sample and the number that was too small**, beside the floor
    /// it fell short of. A step-4 failure is read out of a log on a cohort run of hundreds of
    /// samples, so a message that omits which sample sends the reader back to the data.
    #[test]
    fn each_failure_names_the_sample_the_read_group_and_the_quantity_that_was_wrong() {
        let no_stratum = SsrEstimationError::NoFittableStratumAtPeriod {
            sample: "SL_landrace_07".to_string(),
            read_group: ReadGroupId(3),
            period: period(4),
            thickest_loci: 812,
        };
        let message = no_stratum.to_string();
        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(
            message.contains("read group 3"),
            "which library's fit, since the sample has several: {message}"
        );
        assert!(message.contains("period 4"), "{message}");
        assert!(
            message.contains("812"),
            "how far short it fell — 812 loci wants a wider run, 3 wants dropping: {message}"
        );
        assert!(
            message.contains("1000"),
            "the floor it fell short of: {message}"
        );
    }

    /// **The floors render as the values they are set to**, which is what the messages above
    /// leave unpinned: interpolating a constant means the message follows it silently, so a
    /// floor raised tenfold changes what the caller is told with no test failing. These are
    /// the values the design derives, and moving one is a decision rather than an edit.
    #[test]
    fn the_three_floors_hold_the_values_their_derivations_give() {
        assert_eq!(MIN_LOCI_TO_FIT, 1_000);
        assert_eq!(
            MIN_SLIPPED_READS_TO_FIT_SHARES, 4_000,
            "the fall-off binds: holding it to 6% of itself takes about 4,000 slipped reads, \
             against about 1,400 for the direction share"
        );
        assert_eq!(
            START_AGREEMENT_LIMIT, 1.06,
            "one quarter-Phred, borrowed from the generic path's ladder spacing"
        );
    }

    /// **A search that did not settle is not a thin stratum**, and the message has to say
    /// which happened: a thin stratum borrows and records that it did, while this one returned
    /// the place it stopped. The spread and the limit both render, so a reader sees how far
    /// apart the starts landed without going to find the constant.
    #[test]
    fn the_unsettled_search_message_says_it_is_not_too_little_data() {
        let unsettled = SsrEstimationError::SlippageNotIdentified {
            sample: "SL_landrace_07".to_string(),
            read_group: ReadGroupId(3),
            stratum: stratum(2, 6),
            spread: 333.0,
            starts: 4,
        };
        let message = unsettled.to_string();

        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(message.contains("read group 3"), "{message}");
        assert!(message.contains("period 2, 6 repeats"), "{message}");
        assert!(message.contains("333.0x"), "how far apart: {message}");
        assert!(
            message.contains("across 4 starting points"),
            "how many starts stood behind that spread — two starts landing 333x apart is a \
             different claim from four doing so: {message}"
        );
        assert!(
            message.contains("1.06"),
            "and what that was measured against: {message}"
        );
        assert!(
            message.contains("not a slippage level"),
            "the reader must not take this for a fitted number: {message}"
        );
    }

    /// A domain violation carries three things a reader needs and the inner error has only
    /// one: the quantity comes from the newtype, the sample and the fit have to be attached
    /// here, because on a cohort run those are what say where to look.
    #[test]
    fn a_domain_violation_names_the_sample_and_the_fit_as_well_as_the_quantity() {
        let rejected = SsrEstimationError::Domain {
            sample: "SL_landrace_07".to_string(),
            read_group: ReadGroupId(3),
            stratum: stratum(2, 6),
            fit: "the slippage search",
            source: DomainError::SlipGainShare(1.5),
        };
        let message = rejected.to_string();

        assert!(message.contains("SL_landrace_07"), "the sample: {message}");
        assert!(
            message.contains("read group 3") && message.contains("period 2, 6 repeats"),
            "which of several hundred fits raised it: {message}"
        );
        assert!(
            message.contains("the slippage search"),
            "the fit: {message}"
        );
        assert!(message.contains("1.5"), "the offending value: {message}");
        assert!(
            message.contains("gain share of slipped reads"),
            "the quantity, in the newtype's own words: {message}"
        );
    }
}
