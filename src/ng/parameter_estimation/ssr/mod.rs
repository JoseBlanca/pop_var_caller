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
use std::sync::Arc;

use smallvec::SmallVec;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::parameter_estimation::fitting::FitTermination;
use crate::ng::parameter_estimation::fitting::multistart::{SearchPrecision, fit_by_multistart};
use crate::ng::parameter_estimation::generic::accumulators::PloidyMap;
use crate::ng::parameter_estimation::ssr::locus_offsets::{
    LocusStratum, base_comparison_of, shape_of, stratum_of, tally_of,
};
use crate::ng::parameter_estimation::ssr::slippage::{
    SlippageModel, SsrNoiseModel, slippage_starts,
};
use crate::ng::parameter_estimation::ssr::stratum_table::{
    StratumCell, StratumEntry, StratumTable,
};
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::types::{DomainError, ErrorRate, Ploidy, ReadGroupId, SsrPeriod};

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
    /// **bases** stood behind it.
    ///
    /// Bases and not loci — the architecture's sketch of this field says loci and is behind the
    /// code — because the number exists to be read beside the SNP/indel path's rate for the same
    /// library, whose own warrant counts reads times the sites they covered. Two warrants on
    /// different scales cannot be compared, and a locus count says the wrong thing anyway: what
    /// stands behind a per-base rate is bases.
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

/// What one table of evidence — and the one fit built from it — is *about*: a library, a
/// stratum, and how many genome copies its loci sit on.
///
/// **The ploidy is in the key and not merely carried alongside it**, and that is the one part
/// of this type worth arguing. A table's entries are the same objects whatever the ploidy — a
/// shape is a count of reads at each offset, and nothing about it knows how many chromosomes
/// produced them — so pooling a haploid locus with a diploid one is invisible while the loci
/// are being counted. It becomes wrong at the fit: the fit scores each entry against the
/// genotypes of *one* ploidy, and a pooled table has no ploidy that is true of all of it. The
/// SNP/indel path keys all three of its tables the same way for the same reason
/// (`generic/accumulators.rs`), and there a mixed key was the failure that had to be designed
/// out rather than one that had been observed.
///
/// **It costs nothing on today's runs.** Every genome region is currently declared to have
/// the same ploidy ([`ConstantPloidy`](super::generic::accumulators::ConstantPloidy) is
/// production's only map), so every key carries the same value and the tables are exactly
/// those a two-part key would have built. What it buys is that the first sex chromosome or
/// mixed-ploidy genome to arrive splits the tables instead of silently merging them.
///
/// The field order is the order the fits walk it in: read group, then stratum, then ploidy.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct StratumKey {
    /// Which library. Slippage is a property of the chemistry, so each library is fitted
    /// separately.
    pub read_group: ReadGroupId,
    /// The motif period and the reference repeat count.
    pub stratum: Stratum,
    /// How many genome copies these loci sit on — the set of genotypes the fit will score
    /// each of this table's entries against.
    pub ploidy: Ploidy,
}

impl fmt::Display for StratumKey {
    /// For the error messages and the summary, which name a fit by all three of these.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "read group {}, {}, ploidy {}",
            self.read_group, self.stratum, self.ploidy
        )
    }
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
    /// One fit per [`StratumKey`]. Traceable rather than read: a fit that looks wrong has to be
    /// findable, and nothing downstream is expected to walk this.
    pub by_stratum: BTreeMap<StratumKey, StratumFit>,
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
        "sample {sample}, read group {read_group}: stratum ({stratum}) at ploidy {ploidy} reached \
         slippage levels spanning {spread:.1}x across {starts} starting points, more than the \
         {START_AGREEMENT_LIMIT} that leaves the level identified — this is a search that did \
         not settle, not a slippage level"
    )]
    SlippageNotIdentified {
        sample: String,
        /// Which library's fit did not settle.
        read_group: ReadGroupId,
        stratum: Stratum,
        /// And at which ploidy: one stratum is fitted once per ploidy its loci sat on, so the
        /// three above do not name one fit on a genome that carries more than one.
        ploidy: Ploidy,
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

// ---------------------------------------------------------------------
// The accumulator: one table per StratumKey, filled a locus at a time.
// ---------------------------------------------------------------------

/// Step 4's STR front door: one [`StratumTable`] per [`StratumKey`], and the tally of
/// everything this unit did to a locus other than enter it as it arrived.
///
/// **One per region shard**, merged when the shards are done. Merging is exact — every table's
/// entries are integer counts and every counter is a plain sum — so a genome cut any way and
/// merged in any order gives the tables one walk would have built. That is what the read cap's
/// position-seeded draw is for: without it the two would differ by a few reads at each deep
/// locus.
///
/// **No trait over it**, unlike the fitting seam: nothing generic drives an accumulator, and the
/// walk knows which object it is filling (`arch/parameter_prepass_ssr.md` §5).
pub struct SsrAccumulators {
    by_stratum: BTreeMap<StratumKey, StratumTable>,
    /// How many genome copies a locus sits on. **Read for every locus filed**, because it is
    /// part of the key: an entry looks the same whatever the ploidy, so a table that pooled
    /// two of them would be wrong only at the fit, where there is no ploidy true of all of it
    /// ([`StratumKey`]).
    ploidy: Arc<dyn PloidyMap>,
    counts: SsrAccumulationCounts,
}

impl SsrAccumulators {
    /// One accumulator, for one region shard.
    #[must_use]
    pub fn new(ploidy: Arc<dyn PloidyMap>) -> Self {
        Self {
            by_stratum: BTreeMap::new(),
            ploidy,
            counts: SsrAccumulationCounts::default(),
        }
    }

    /// Add one locus.
    ///
    /// **Borrows it and passes it on untouched** — the caller keeps the locus and hands it to
    /// whatever else reads the stream — and **tallies rather than repairs**: a locus this unit
    /// cannot file is counted and skipped, never rounded into a stratum it does not belong to.
    ///
    /// A locus that is not one repeat tract is passed over in silence. One whose reference tract
    /// is not a whole number of motif copies is counted in
    /// [`loci_without_whole_repeat_reference`](SsrAccumulationCounts::loci_without_whole_repeat_reference),
    /// which should read near zero and is a bug report against the classification that admitted
    /// it if it does not.
    ///
    /// **A locus covered by two read groups makes one entry in each group's table, and that is
    /// sound**: the genotype is drawn once for the locus and enters both through the same
    /// mixture, so the product over them is a composite likelihood — consistent, and losing
    /// precision rather than correctness. What must not be split is a locus's reads *within* one
    /// read group, which is what keying an entry by locus prevents.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations) {
        let stratum = match stratum_of(locus) {
            LocusStratum::Stratified(stratum) => stratum,
            LocusStratum::NotOneRepeatTract => return,
            LocusStratum::WithoutWholeRepeatReference { .. } => {
                self.counts.loci_without_whole_repeat_reference += 1;
                return;
            }
        };

        let ploidy = self.ploidy.ploidy_at(locus.region);

        // Counted once for the locus rather than once per read group: it is a property of the
        // locus, and the generator does not attribute those reads to a library.
        self.counts.reads_without_observation += u64::from(locus.reads_without_observation);
        if locus.reads_discarded_by_cap > 0 {
            self.counts.loci_with_upstream_subsample += 1;
        }

        // The read groups this locus actually witnessed, in a fixed order, so that two shards
        // walking the same locus build the same tables in the same sequence. Taken from the
        // observations rather than from a declared list, because a list is a second place for a
        // read group to be missing from.
        let mut read_groups: SmallVec<[ReadGroupId; 2]> = locus
            .observations
            .iter()
            .map(|observation| observation.read_group)
            .collect();
        read_groups.sort_unstable();
        read_groups.dedup();

        let mut cap_fired = false;
        for read_group in read_groups {
            let tally = tally_of(locus, read_group);
            self.counts.reads_with_partial_witness += u64::from(tally.reads_with_partial_witness());

            let Some(entered) = shape_of(tally, locus.region) else {
                // This library witnessed no length here. The next one may still have.
                continue;
            };
            cap_fired |= entered.subsampled_from().is_some();

            self.by_stratum
                .entry(StratumKey {
                    read_group,
                    stratum,
                    ploidy,
                })
                .or_default()
                .add_locus(entered.shape(), base_comparison_of(locus, read_group));
        }
        // **Once for the locus, however many of its libraries were thinned**, because that is
        // what the counter is named after and what a reader compares against the locus count. The
        // cap itself fires per library — each library's reads are drawn from separately — so
        // counting each firing would let a field called `loci_…` exceed the loci walked.
        if cap_fired {
            self.counts.loci_subsampled_to_cap += 1;
        }
    }

    /// Combine a shard's tables and counters into these. Associative and exact.
    ///
    /// # Panics
    ///
    /// If the two shards were not built on the same ploidy map. The ploidy is part of the key, so
    /// the damage is not a pooled table but a **split** one: the same stratum would arrive
    /// under two keys because the two shards disagreed about the genome rather than because
    /// the genome has two ploidies there, and each half would then be fitted on part of its
    /// evidence. The same guard the SNP/indel path's merge carries, and for the same reason.
    pub fn merge(&mut self, other: Self) {
        // Pointer identity rather than equality, as the sibling path's merge does: a shard's
        // accumulator is built from one shared map, so two shards that disagree here were driven
        // by different configurations rather than by two equal maps.
        assert!(
            Arc::ptr_eq(&self.ploidy, &other.ploidy),
            "these two shards were not built on the same ploidy map, so a stratum they both saw \
             could arrive here under two keys and each half be fitted on part of its evidence"
        );

        for (key, table) in other.by_stratum {
            match self.by_stratum.get_mut(&key) {
                Some(mine) => mine.merge(&table),
                None => {
                    self.by_stratum.insert(key, table);
                }
            }
        }
        self.counts.merge(&other.counts);
    }

    /// Every stratum's evidence, in [`StratumKey`] order — read group, then stratum, then
    /// ploidy — which is the order the fits walk it in, and a property of the contents rather
    /// than of when a locus arrived.
    ///
    /// An iterator rather than the map itself, so that how the tables are stored stays this
    /// unit's business: the same reason `StratumTable` hands out its entries as a list.
    pub fn strata(&self) -> impl Iterator<Item = (StratumKey, &StratumTable)> {
        self.by_stratum.iter().map(|(&key, table)| (key, table))
    }

    /// One key's evidence, or `None` where it holds no loci.
    #[must_use]
    pub fn table_for(&self, key: StratumKey) -> Option<&StratumTable> {
        self.by_stratum.get(&key)
    }

    /// How many keys hold any loci at all.
    #[must_use]
    pub fn stratum_count(&self) -> usize {
        self.by_stratum.len()
    }

    /// Everything this unit did to a locus other than enter it as it arrived.
    #[must_use]
    pub fn adjustments(&self) -> &SsrAccumulationCounts {
        &self.counts
    }
}

// ---------------------------------------------------------------------
// The first of the four fits: the substitution rate, which is a division.
// ---------------------------------------------------------------------

/// One table's substitution rate with its warrant: mismatched bases over compared bases, how
/// many bases that was, and that it was measured here rather than taken from somewhere else.
///
/// **The count beside the rate is bases and not loci**, and the deciding argument is the
/// comparison this rate exists to be put into. A stratum's rate is meant to be read against the
/// SNP/indel path's rate for the same library, where the two must agree wherever a tract barely
/// slips (`spec/parameter_prepass_ssr.md` §4.5); that path counts an error rate's observations
/// as *reads times the sites they covered*, which is one base per read per site — this quantity
/// exactly. A locus count would put the two warrants on different scales, and would say the
/// wrong thing on its own besides: a hundred loci at three reads each stand behind less of this
/// rate than one locus at three hundred.
///
/// **`None` where nothing was compared, never a zero** — the rule
/// [`StratumTable::substitution_rate`] states and the reason it returns an `Option`: a stratum
/// whose reads all matched has a measured rate of zero, and a stratum nothing was compared in
/// has no rate at all. A zero here would be borrowed from, and later compared against the
/// SNP/indel path's rate, as though it had been measured.
///
/// **Always [`Provenance::FittedHere`], and there is no rung below it for this parameter.** The
/// borrowing the design defines is over the slippage model — the level and its two shares — not
/// over the substitution rate, and [`StratumFit::substitution`] is not an `Option`, so a stratum
/// with loci and no compared bases is a case the design has not ruled on. Nothing observed can
/// reach it: a read reaches a table only through a complete witness, which compares its bases.
#[must_use]
pub fn substitution_rate_of(table: &StratumTable) -> Option<Estimate<ErrorRate>> {
    table.substitution_rate().map(|value| Estimate {
        value,
        provenance: Provenance::FittedHere,
        observations: table.bases_compared(),
    })
}

/// Every stratum's substitution rate, one division per [`StratumKey`] — **the first of the four
/// fits, and the only one that is not a search** (`spec/parameter_prepass_ssr.md` §4.2).
///
/// It goes first because it needs none of the other three: a read's mismatch count is binomial
/// at this rate whatever tract length the read showed, so the channel scoring *which length a
/// read showed* and the channel scoring *how its bases read* factorise exactly, and this one
/// closes in a division (`spec/parameter_prepass_ssr.md` §4.1).
///
/// **Nothing is pooled across keys.** Each table's own bases give each key's own rate — libraries
/// differ in chemistry, which is why they are fitted apart, and tracts differ in how much of
/// their length is repeat. A rate pooled over a read group's strata would be a plausible number
/// with no stratum it describes.
///
/// A key whose table compared no bases is **absent from the result** rather than present with a
/// zero, so a consumer that has to distinguish "measured at zero" from "not measured" still can.
#[must_use]
pub fn substitution_rates(
    accumulators: &SsrAccumulators,
) -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
    accumulators
        .strata()
        .filter_map(|(key, table)| substitution_rate_of(table).map(|rate| (key, rate)))
        .collect()
}

// ---------------------------------------------------------------------
// The second fit: the three slippage parameters, searched per stratum.
// ---------------------------------------------------------------------

/// What one stratum's slippage search returned: the three parameters, the genotype frequencies
/// they were fitted against, and enough beside them to tell a fitted answer from a stopped
/// search.
///
/// **Not [`StratumFit`] yet**, and the difference is provenance. Everything here was measured on
/// this stratum's own loci; whether that is what the stratum ends up *reporting* is decided
/// afterwards, by the two floors that may borrow the model from a neighbour and by the
/// monotonicity walk that may merge two strata and refit. This type is what those steps read.
#[derive(Clone, PartialEq, Debug)]
pub struct StratumSlippageFit {
    /// The best-scoring of the starts: how often a read slips, which way, and how far.
    pub model: SlippageModel,
    /// How often each genotype occurred, in the order
    /// [`SsrNoiseModel::for_each_genotype`](slippage::SsrNoiseModel::for_each_genotype) visits
    /// them — which is the walk that gives them meaning, and the reason that method is public.
    ///
    /// **Kept as bare frequencies here** rather than as allele pairs: the walk is over `P`-tuples
    /// at ploidy `P`, and a pair is only the diploid case of one.
    pub genotype_frequencies: Vec<f64>,
    /// The allele lengths those frequencies are over, ascending, as offsets from the reference
    /// tract length.
    ///
    /// **Carried rather than left to be rebuilt from the stratum**, because after a merge the two
    /// are not the same answer: the monotonicity walk pools two neighbouring repeat counts and
    /// refits with the model of the **lower** one, so a consumer rebuilding the support from the
    /// key would index the frequencies by a support two lengths wider than the one they were
    /// fitted over.
    pub allele_support: SmallVec<[WholeRepeatOffset; 13]>,
    /// What [`Self::model`] scored, as a natural logarithm.
    pub log_likelihood: f64,
    /// Every starting point, best-scoring first — the diagnostic several starts exist to produce.
    pub starts_tried: SmallVec<[SlippageStart; 4]>,
    /// The highest slippage level any start reached, divided by the lowest.
    ///
    /// Compared against [`START_AGREEMENT_LIMIT`] by whoever is in a position to name the sample
    /// and the library, which this function is not.
    pub start_spread: f64,
    /// Whether **every** inner climb over the genotype frequencies reached stillness, or some ran
    /// out of passes.
    ///
    /// **False is not an error, and this is where that is decided.** The surface the climb walks
    /// is concave, so a climb that ran out did not find a wrong summit — it ran out of time on the
    /// right one. What differs from the sibling path, where the same condition is treated as a
    /// bug, is how far there is to go: that path's climb has three genotypes and its pass cap was
    /// measured on them, while a stratum here has between 66 and 91, most heading to a frequency
    /// of zero. Making this fatal would end most real runs; ignoring it would hide the one thing
    /// that makes a candidate score below its own summit, so a neighbouring candidate can win on
    /// that alone. So it is reported and never silently dropped, and the summary counts it.
    pub every_climb_settled: bool,
    /// How the **outer** search over the three parameters ended — settled, or out of sweeps — and
    /// how many sweeps that took.
    ///
    /// **A different question from [`Self::every_climb_settled`]**, which is about the climb over
    /// the genotype frequencies inside each candidate. This one has no concavity proof behind it
    /// at all, which is why the search is capped and the best-scoring iterate kept
    /// (`arch/parameter_prepass_ssr.md` §3).
    pub termination: FitTermination,
    /// Reads that showed a length **a whole number of motif copies** away from the reference
    /// tract's — the count that decides whether the direction share and the fall-off are
    /// measurable at all.
    ///
    /// **Whole-repeat movement only, which is narrower than "a length other than the
    /// reference's".** The reads that moved by something else — an ordinary indel inside the
    /// tract, an interruption — are counted apart, in [`Self::guard_reads`], because they carry
    /// nothing about *which way* a slip went or *how far*, and those two are exactly what this
    /// count gates ([`MIN_SLIPPED_READS_TO_FIT_SHARES`]). Pooling them would overstate the
    /// evidence behind the two shares by the guard share of the stratum, which spec §5 measures
    /// at up to 58.5% — a stratum reporting the floor's 4,000 would then have put about 1,660
    /// reads behind its fall-off.
    pub slipped_reads: u64,
    /// Reads that differed from the reference tract length by something **other** than a whole
    /// number of copies. Reported beside the count above so the two are never confused.
    pub guard_reads: u64,
    /// Reads that sat on the whole-repeat grid, and so were scored by the model — the denominator
    /// the slippage level is a share of.
    pub scored_reads: u64,
    /// How many loci stood behind the fit.
    pub loci: u64,
}

/// **Fit one stratum's three slippage parameters from its own loci**, searching from four
/// starting points spread over all three and climbing the genotype frequencies at every trial
/// (`spec/parameter_prepass_ssr.md` §4.2).
///
/// The genotype is summed over rather than guessed at: at each candidate the frequencies of the
/// unordered allele tuples are climbed to their best — the half of the problem with a concavity
/// proof behind it — and only then is the candidate scored. **The frequencies are fitted freely**
/// rather than tied through one allele frequency, matching the SNP/indel path and for the same
/// reason: a Hardy–Weinberg tie presumes the inbreeding coefficient is zero, and that is a
/// quantity this run measures rather than assumes.
///
/// **Four starts and not one, spread over all three parameters.** Starts that disagree about the
/// headline number while sharing one guess at a nuisance axis are how the SNP/indel path's
/// inbreeding fit returned a confident zero on a genome 29% covered by runs of homozygosity.
/// Where the four land is reported whatever they say.
///
/// The four starts are placed around the share of this stratum's reads that moved a whole number
/// of copies, which over-estimates the slippage level because it counts real non-reference alleles
/// too — which is why the four multipliers run below one as well as above.
///
/// **`None` where no read of the stratum sat on the whole-repeat grid**, and that is a guard
/// rather than a tidiness: the model scores only whole-repeat movement, so every genotype and
/// every candidate score alike on such a stratum, and a search over a flat surface returns
/// wherever its steps happened to stop. Measured on 500 loci whose every read carried an ordinary
/// indel inside the tract, the search came back at a level of **0.5976** — the top of its range —
/// with the four starts agreeing to 1.00, so the one diagnostic that guards this fit could not
/// see it. Such a stratum has as many loci as any other, so the borrowing step's floor cannot see
/// it either; refusing the fit here is what makes it borrow.
///
/// # Panics
///
/// If the table holds no entries, or if the search returns frequencies at some ploidy other than
/// the cells'. Neither can arise from an accumulator: a table exists because a locus made it, and
/// every cell of one fit carries the ploidy of the key its table is stored under.
#[must_use]
pub fn fit_slippage(
    table: &StratumTable,
    stratum: Stratum,
    ploidy: Ploidy,
    precision: SearchPrecision,
) -> Option<StratumSlippageFit> {
    let entries = table.entries();
    assert!(
        !entries.is_empty(),
        "a stratum with no entries has nothing to fit: {stratum}"
    );

    let reads = read_counts(&entries);
    if reads.scored == 0 {
        return None;
    }

    let cells: Vec<StratumCell> = entries
        .into_iter()
        .map(|entry| StratumCell::new(entry, ploidy))
        .collect();

    let model = SsrNoiseModel::for_stratum(stratum.repeats);
    // A number, because `scored` is not zero by the guard above; `slippage_starts` refuses a
    // `NaN` in any case.
    let starts = slippage_starts(reads.slipped as f64 / reads.scored as f64);
    let fitted = fit_by_multistart(&model, &cells, &starts, precision);

    let genotype_frequencies = fitted
        .genotype_frequencies
        .get(&ploidy)
        .expect("the fit climbed frequencies at the one ploidy every cell carries")
        .clone();

    Some(StratumSlippageFit {
        model: fitted.best,
        genotype_frequencies,
        allele_support: SmallVec::from_slice(model.allele_support()),
        log_likelihood: fitted.log_likelihood.get(),
        starts_tried: fitted
            .starts
            .into_iter()
            .map(|outcome| SlippageStart {
                from: outcome.from,
                reached: outcome.reached,
                log_likelihood: outcome.log_likelihood.get(),
            })
            .collect(),
        start_spread: fitted.headline_spread,
        every_climb_settled: fitted.every_climb_settled,
        termination: fitted.termination,
        slipped_reads: reads.slipped,
        guard_reads: reads.guard,
        scored_reads: reads.scored,
        loci: table.loci(),
    })
}

/// How a stratum's reads divide, once each entry is weighted by how many loci showed it.
///
/// **Three counts and not two, because "moved" and "moved by a whole number of copies" are
/// different populations** and the difference reaches 58.5% of the moved reads on the strata spec
/// §5 measures. The model explains only the second, so only the second belongs in the level's
/// numerator or in the count that gates the two shares.
struct StratumReads {
    /// Reads a whole number of copies away from the reference tract length.
    slipped: u64,
    /// Reads a length away from it that is not a whole number of copies.
    guard: u64,
    /// Reads on the whole-repeat grid, at the reference length or away from it — every read the
    /// model scores.
    scored: u64,
}

/// Count [`StratumReads`] over a table's entries.
///
/// No sum can overflow: each is bounded by the locus count times [`MAX_LOCUS_READS`], which the
/// table's own per-entry and entry-count limits hold at about 3.3 × 10¹⁶.
fn read_counts(entries: &[StratumEntry]) -> StratumReads {
    let mut reads = StratumReads {
        slipped: 0,
        guard: 0,
        scored: 0,
    };
    for entry in entries {
        let loci = entry.loci;
        let guard = u64::from(entry.shape.reads_not_whole_repeat());
        reads.slipped += (u64::from(entry.shape.reads_off_reference()) - guard) * loci;
        reads.guard += guard * loci;
        reads.scored += u64::from(entry.shape.whole_repeat_depth()) * loci;
    }
    reads
}

/// **Every stratum's slippage, fitted from its own loci** — the second of the four fits, and the
/// one that is a search.
///
/// The strata are walked in [`StratumKey`] order, which is a property of what the tables hold
/// rather than of when a locus arrived, so two runs over the same evidence fit in the same
/// sequence.
///
/// # Errors
///
/// [`SsrEstimationError::SlippageNotIdentified`] where a stratum's four starts reached levels
/// spanning more than [`START_AGREEMENT_LIMIT`]. **That is not "too little data"** — too little
/// data is a borrow, and it is recorded as one — it is a search that did not settle, and the
/// number it would otherwise report is where it stopped rather than what the loci say.
///
/// A stratum **no read of which sat on the whole-repeat grid** is absent from the result rather
/// than fitted, for the reason [`fit_slippage`] gives: there is nothing there for the model to
/// explain, and what a search returns on a flat surface is where its steps stopped. Filling that
/// gap is the borrowing step's business.
pub fn fit_slippage_by_stratum(
    accumulators: &SsrAccumulators,
    sample: &str,
    precision: SearchPrecision,
) -> Result<BTreeMap<StratumKey, StratumSlippageFit>, SsrEstimationError> {
    let mut fits = BTreeMap::new();
    for (key, table) in accumulators.strata() {
        let Some(fit) = fit_slippage(table, key.stratum, key.ploidy, precision) else {
            continue;
        };
        starts_must_agree(&fit, sample, key)?;
        fits.insert(key, fit);
    }
    Ok(fits)
}

/// Refuse a fit whose starting points did not converge on one answer.
///
/// **Separate from the search because only the caller can name whose fit it was.** The search
/// returns the spread and makes no judgement — how far apart two answers may sit before a fit is
/// disowned is this path's call, and a message naming neither the sample nor the library locates
/// nothing on a cohort run of hundreds of samples.
///
/// # Errors
///
/// [`SsrEstimationError::SlippageNotIdentified`], carrying the spread and how many starts produced
/// it.
pub fn starts_must_agree(
    fit: &StratumSlippageFit,
    sample: &str,
    key: StratumKey,
) -> Result<(), SsrEstimationError> {
    if fit.start_spread > START_AGREEMENT_LIMIT {
        return Err(SsrEstimationError::SlippageNotIdentified {
            sample: sample.to_owned(),
            read_group: key.read_group,
            stratum: key.stratum,
            ploidy: key.ploidy,
            spread: fit.start_spread,
            starts: fit.starts_tried.len(),
        });
    }
    Ok(())
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

    // -----------------------------------------------------------------
    // The accumulator.
    // -----------------------------------------------------------------

    use crate::ng::locus_generation::{
        LocusKind, LocusLen, ReadWitness, SequenceObservation, SsrDetail,
    };
    use crate::ng::parameter_estimation::generic::accumulators::ConstantPloidy;
    use crate::ng::types::{ContigId, GenomeRegion, Motif, Ploidy, Position};

    fn diploid() -> Arc<dyn PloidyMap> {
        Arc::new(ConstantPloidy(
            Ploidy::try_new(2).expect("two genome copies"),
        ))
    }

    fn observation(bases: &[u8], group: u32, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(group),
            num_obs,
            num_fwd: num_obs / 2,
            q_sum: -10.0,
            mapq_sum: 60 * num_obs,
            mapq_sum_sq: 3_600 * u64::from(num_obs),
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// A tract at `start`, tiled by `motif`, carrying whatever reads a test wants.
    fn tract(
        start: u64,
        reference_bases: &[u8],
        motif: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(start + reference_bases.len().max(1) as u64 - 1),
            },
            reference_bases: Box::from(reference_bases),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(motif).expect("a motif inside the STR period range"),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        }
    }

    fn dinucleotide(start: u64, reads_at_reference: u32) -> SampleLocusObservations {
        let reference = b"ATATATATATATATATATAT";
        tract(
            start,
            reference,
            b"AT",
            vec![observation(reference, 0, reads_at_reference)],
        )
    }

    /// One stratum's table at the ploidy [`diploid`] gives every locus, for a test that knows
    /// it is there.
    fn table(
        accumulators: &SsrAccumulators,
        group: u32,
        period_bases: u8,
        repeats: u32,
    ) -> &StratumTable {
        accumulators
            .table_for(StratumKey {
                read_group: ReadGroupId(group),
                stratum: stratum(period_bases, repeats),
                ploidy: Ploidy::try_new(2).expect("two genome copies"),
            })
            .expect("a table for that read group, stratum and ploidy")
    }

    /// Two genome copies below `haploid_from`, one at or above it — the only `PloidyMap` in the
    /// tests that returns more than one answer, and the only way to reach the mixed-ploidy case
    /// the key exists for while `ConstantPloidy` is production's only map.
    struct PloidyChangesAt {
        haploid_from: u64,
    }

    impl PloidyMap for PloidyChangesAt {
        fn ploidy_at(&self, region: GenomeRegion) -> Ploidy {
            let copies = if region.start.get() >= self.haploid_from {
                1
            } else {
                2
            };
            Ploidy::try_new(copies).expect("one or two genome copies")
        }
    }

    /// **Loci reduce to entries, filed by read group and stratum**, and two loci of the same
    /// period and reference repeat count land in one table however far apart they sit.
    #[test]
    fn loci_of_the_same_stratum_land_in_one_table() {
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&dinucleotide(1_000, 5));
        accumulators.add_locus(&dinucleotide(9_000, 5));
        accumulators.add_locus(&tract(
            2_000,
            b"AAAAAAAA",
            b"A",
            vec![observation(b"AAAAAAAA", 0, 4)],
        ));

        assert_eq!(
            accumulators.stratum_count(),
            2,
            "two strata, not two loci and not one"
        );
        assert_eq!(table(&accumulators, 0, 2, 10).loci(), 2);
        assert_eq!(table(&accumulators, 0, 1, 8).loci(), 1);
    }

    /// **Two loci of the same period and repeat count sitting on different numbers of genome
    /// copies land in two tables**, not one.
    ///
    /// The failure this makes unreachable is silent everywhere else: the entries a haploid
    /// locus and a diploid one produce are the same kind of object — reads at offsets — so a
    /// pooled table looks healthy and is only wrong at the fit, which scores every entry
    /// against the genotypes of a single ploidy. Nothing counts it, nothing errors, and the
    /// slippage level comes back describing a population that does not exist.
    ///
    /// The two loci are otherwise identical — same motif, same reference length, same read
    /// group, same depth — so the only thing that can separate them is the ploidy map's
    /// answer. With the ploidy left out of the key both would be one table of two loci.
    #[test]
    fn two_ploidies_of_one_stratum_are_two_tables() {
        let mut accumulators = SsrAccumulators::new(Arc::new(PloidyChangesAt {
            haploid_from: 5_000,
        }));
        accumulators.add_locus(&dinucleotide(1_000, 5));
        accumulators.add_locus(&dinucleotide(9_000, 5));

        assert_eq!(
            accumulators.stratum_count(),
            2,
            "one stratum at two ploidies is two tables, not one table of two loci"
        );
        for (copies, start) in [(2u8, 1_000u64), (1, 9_000)] {
            let key = StratumKey {
                read_group: ReadGroupId(0),
                stratum: stratum(2, 10),
                ploidy: Ploidy::try_new(copies).expect("one or two genome copies"),
            };
            let table = accumulators
                .table_for(key)
                .unwrap_or_else(|| panic!("a table at {key}"));
            assert_eq!(
                table.loci(),
                1,
                "only the locus at {start} sits on {copies} genome copies"
            );
        }
    }

    /// **A locus covered by two libraries makes one entry in each**, because slippage is a
    /// property of the chemistry: pooling them would fit one stutter model to two.
    #[test]
    fn a_locus_two_libraries_witnessed_makes_one_entry_in_each_of_their_tables() {
        let reference = b"ATATATATATATATATATAT";
        let locus = tract(
            1_000,
            reference,
            b"AT",
            vec![
                observation(reference, 0, 5),
                observation(b"ATATATATATATATATAT", 1, 4),
            ],
        );

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&locus);

        assert_eq!(accumulators.stratum_count(), 2);
        assert_eq!(table(&accumulators, 0, 2, 10).loci(), 1);
        assert_eq!(table(&accumulators, 1, 2, 10).loci(), 1);
        assert_ne!(
            table(&accumulators, 0, 2, 10).entries(),
            table(&accumulators, 1, 2, 10).entries(),
            "the two libraries saw different lengths, so their entries differ"
        );
    }

    /// **A locus this unit cannot file is counted and skipped, never rounded.** A locus of
    /// another kind is passed over in silence; a tract whose reference length is not a whole
    /// number of copies is a delimiting fault upstream and is reported.
    #[test]
    fn a_locus_with_no_stratum_is_counted_or_passed_over_but_never_entered() {
        let mut accumulators = SsrAccumulators::new(diploid());

        let mut generic = dinucleotide(1_000, 5);
        generic.kind = LocusKind::Generic;
        accumulators.add_locus(&generic);

        accumulators.add_locus(&tract(
            2_000,
            b"CAGCAGCAGCAGC",
            b"CAG",
            vec![observation(b"CAGCAGCAGCAGC", 0, 5)],
        ));

        assert_eq!(accumulators.stratum_count(), 0);
        assert_eq!(
            accumulators
                .adjustments()
                .loci_without_whole_repeat_reference,
            1,
            "the fractional tract is reported and the generic locus is not"
        );
    }

    /// The adjustments are what a run reads to tell a shallow sample from a capped one, and each
    /// counts a different thing: loci the cap thinned, loci the generator had already thinned,
    /// reads that witnessed nothing, and reads whose witness was partial.
    #[test]
    fn the_adjustments_count_what_was_done_to_each_locus() {
        let reference = b"ATATATATATATATATATAT";
        let mut deep = tract(
            1_000,
            reference,
            b"AT",
            vec![
                observation(reference, 0, 200),
                SequenceObservation {
                    read_witness: ReadWitness::from_left(8, LocusLen::from_positions(20))
                        .expect("a run of eight positions"),
                    ..observation(b"ATATATAT", 0, 7)
                },
            ],
        );
        deep.reads_without_observation = 3;
        deep.reads_discarded_by_cap = 11;
        // Both thinning counters read one here, so this fixture cannot tell them apart; the
        // fixture that can is `the_locus_counters_count_loci_however_many_libraries_covered_them`.

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&deep);
        accumulators.add_locus(&dinucleotide(2_000, 4));

        assert_eq!(
            *accumulators.adjustments(),
            SsrAccumulationCounts {
                loci_subsampled_to_cap: 1,
                loci_with_upstream_subsample: 1,
                reads_without_observation: 3,
                reads_with_partial_witness: 7,
                loci_without_whole_repeat_reference: 0,
            }
        );
    }

    /// **A genome cut any way and merged in any order gives the tables one walk would have
    /// built** — entry for entry and counter for counter. This is the property the whole
    /// sharded design rests on, and it is an equality rather than a tolerance because every
    /// number in a table is an integer count and the read cap's draw is seeded from the locus's
    /// own position rather than from the shard's.
    #[test]
    fn a_walk_cut_into_shards_and_merged_equals_the_uncut_walk() {
        let reference = b"ATATATATATATATATATAT";
        let loci: Vec<SampleLocusObservations> = (0..9u64)
            .map(|at| {
                let start = 1_000 + at * 500;
                // Deep enough that the read cap fires, so the draw is under test rather than
                // bypassed: below the cap every implementation agrees.
                tract(
                    start,
                    reference,
                    b"AT",
                    vec![
                        observation(reference, 0, 40 + u32::try_from(at).expect("small")),
                        observation(b"ATATATATATATATATAT", u32::from(at % 2 == 0), 30),
                    ],
                )
            })
            .collect();

        // A mononucleotide stratum at one end only, so a cut can leave it wholly inside one
        // shard — which is what makes the merge's adopt-a-new-key branch run at all.
        let mut loci = loci;
        loci.push(tract(
            9_000,
            b"AAAAAAAAAAAAAAAA",
            b"A",
            vec![observation(b"AAAAAAAAAAAAAAAA", 0, 40)],
        ));

        let mut whole = SsrAccumulators::new(diploid());
        for locus in &loci {
            whole.add_locus(locus);
        }

        let ploidy = diploid();
        for cut in [1usize, 4, 8] {
            let (left, right) = loci.split_at(cut);
            let mut first = SsrAccumulators::new(Arc::clone(&ploidy));
            for locus in left {
                first.add_locus(locus);
            }
            let mut second = SsrAccumulators::new(Arc::clone(&ploidy));
            for locus in right {
                second.add_locus(locus);
            }
            first.merge(second);

            assert_eq!(
                first.stratum_count(),
                whole.stratum_count(),
                "cut after {cut} loci"
            );
            for (key, table) in whole.strata() {
                assert_eq!(
                    first.table_for(key),
                    Some(table),
                    "cut after {cut} loci, at {key}"
                );
            }
            assert_eq!(
                first.adjustments(),
                whole.adjustments(),
                "cut after {cut} loci"
            );
        }
    }

    /// **A stratum only one shard saw survives the merge.** The receiving accumulator holds no
    /// table under that key, so the merge has to adopt the shard's rather than pass it over —
    /// the branch a fixture whose every shard sees every stratum never reaches. Strata are
    /// unevenly spread along a genome, so what this protects is exactly the rarest ones: the
    /// stratum a single shard holds is the one nearest `MIN_LOCI_TO_FIT`.
    #[test]
    fn a_stratum_only_one_shard_saw_survives_the_merge() {
        let ploidy = diploid();
        let mut first = SsrAccumulators::new(Arc::clone(&ploidy));
        first.add_locus(&dinucleotide(1_000, 5));

        let mut second = SsrAccumulators::new(Arc::clone(&ploidy));
        second.add_locus(&tract(
            2_000,
            b"AAAAAAAA",
            b"A",
            vec![observation(b"AAAAAAAA", 0, 4)],
        ));

        first.merge(second);

        assert_eq!(first.stratum_count(), 2);
        assert_eq!(table(&first, 0, 2, 10).loci(), 1);
        assert_eq!(table(&first, 0, 1, 8).loci(), 1);
    }

    /// **Each table holds its own read group's bases**, which is the one thing the entries
    /// cannot show: a shape carries no base counts, so a table handed another library's
    /// comparison looks identical until its substitution rate is read. Under a mix-up every
    /// library's rate becomes the first library's; under a dropped comparison every stratum in
    /// the genome reports no rate at all.
    #[test]
    fn each_read_groups_table_holds_that_read_groups_bases() {
        let reference = b"ATATATATATATATATATAT";
        let locus = tract(
            1_000,
            reference,
            b"AT",
            vec![
                observation(reference, 0, 4),
                observation(b"ACACACACACACACACACAC", 1, 6),
            ],
        );

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&locus);

        let clean = table(&accumulators, 0, 2, 10)
            .substitution_rate()
            .expect("bases were compared");
        let noisy = table(&accumulators, 1, 2, 10)
            .substitution_rate()
            .expect("bases were compared");

        assert!(clean.get().abs() < 1e-12, "{clean:?}");
        assert!(
            (noisy.get() - 0.5).abs() < 1e-12,
            "every second base of that library's reads differs: {noisy:?}"
        );
    }

    /// **Which counters count loci, and which count a locus once per library.** Reads that
    /// witnessed nothing, and a subsample the generator had already made, are properties of the
    /// locus and are counted once however many libraries covered it. The read cap fires per
    /// library — each library's reads are drawn from separately — but the counter is named after
    /// loci and counts them, or a field called `loci_…` could exceed the loci walked.
    ///
    /// The two thinning counters are given **different** values on purpose: they mean opposite
    /// things — this unit's cap against the generator's own reservoir — and a fixture where both
    /// read one passes with them exchanged.
    #[test]
    fn the_locus_counters_count_loci_however_many_libraries_covered_them() {
        let reference = b"ATATATATATATATATATAT";
        let mut two_libraries = tract(
            1_000,
            reference,
            b"AT",
            vec![
                observation(reference, 0, 40),
                observation(b"ATATATATATATATATAT", 0, 31),
                observation(reference, 1, 50),
                observation(b"ATATATATATATATATAT", 1, 23),
            ],
        );
        two_libraries.reads_without_observation = 3;

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&two_libraries);
        // A second locus the generator had already thinned, and this unit's cap did not.
        let mut upstream_thinned = dinucleotide(2_000, 4);
        upstream_thinned.reads_discarded_by_cap = 11;
        accumulators.add_locus(&upstream_thinned);
        accumulators.add_locus(&dinucleotide(3_000, 5));

        assert_eq!(
            *accumulators.adjustments(),
            SsrAccumulationCounts {
                loci_subsampled_to_cap: 1,
                loci_with_upstream_subsample: 1,
                reads_without_observation: 3,
                reads_with_partial_witness: 0,
                loci_without_whole_repeat_reference: 0,
            },
            "one locus thinned by this unit, a different one thinned by the generator, and \
             three reads that witnessed nothing — each counted once for the locus"
        );
    }

    /// A library that witnessed no length at a locus does not take the other libraries' evidence
    /// with it: the loop passes over that library and files the rest.
    #[test]
    fn a_library_that_witnessed_nothing_does_not_cost_the_others_their_entry() {
        let reference = b"ATATATATATATATATATAT";
        let locus = tract(
            1_000,
            reference,
            b"AT",
            vec![
                SequenceObservation {
                    read_witness: ReadWitness::from_left(8, LocusLen::from_positions(20))
                        .expect("a run of eight positions"),
                    ..observation(b"ATATATAT", 0, 5)
                },
                observation(reference, 1, 6),
            ],
        );

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&locus);

        assert_eq!(
            accumulators.stratum_count(),
            1,
            "only the second library files"
        );
        assert_eq!(table(&accumulators, 1, 2, 10).loci(), 1);
        assert_eq!(accumulators.adjustments().reads_with_partial_witness, 5);
    }

    /// **Two shards built on separate ploidy maps cannot be merged, even where the two maps
    /// agree** — which is what this fixture hands it, two allocations of the same constant map.
    /// The guard is pointer identity, deliberately: every shard of a real run is handed one
    /// shared map, so two that differ here were driven by different configurations rather than
    /// by two equal objects, and equality cannot be asked of a trait object anyway.
    ///
    /// What the guard protects against is in `merge`'s own doc: maps that *disagree* file a
    /// stratum both shards saw under two keys, and each half is then fitted on part of its
    /// evidence.
    #[test]
    #[should_panic(expected = "not built on the same ploidy map")]
    fn merging_two_shards_built_on_separate_ploidy_maps_is_refused() {
        let mut first = SsrAccumulators::new(diploid());
        first.add_locus(&dinucleotide(1_000, 5));

        let mut second = SsrAccumulators::new(diploid());
        second.add_locus(&dinucleotide(2_000, 5));

        first.merge(second);
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
            ploidy: Ploidy::try_new(4).expect("four genome copies"),
            spread: 333.0,
            starts: 4,
        };
        let message = unsettled.to_string();

        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(message.contains("read group 3"), "{message}");
        assert!(message.contains("period 2, 6 repeats"), "{message}");
        assert!(
            message.contains("ploidy 4"),
            "one stratum is fitted once per ploidy, so the message names which: {message}"
        );
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

    // -----------------------------------------------------------------
    // The substitution rate, fitted per stratum.
    // -----------------------------------------------------------------

    /// A tract whose reads carry a **known** number of mismatched bases: `clean` reads showing
    /// the reference tract exactly, and `mismatching` reads showing it with the base at
    /// `changed_at` replaced, so each of those contributes exactly one mismatch and no change of
    /// length. The compared bases are `(clean + mismatching) × the tract's length`, and the
    /// mismatched bases are `mismatching`.
    ///
    /// **Counted against the motif tiled from the read's first base, not against `reference`**,
    /// because that is what the accumulator compares against: a read is scored on how far it
    /// departs from a perfect tract, so a reference that is not itself a perfect tiling would
    /// charge its own imperfections to every read. The assertion below is over the tiling, so a
    /// fixture that broke that assumption would fail here rather than quietly report a rate
    /// neither test intended.
    fn tract_at_a_known_rate(
        start: u64,
        reference: &[u8],
        motif: &[u8],
        group: u32,
        clean: u32,
        mismatching: u32,
        changed_at: usize,
    ) -> SampleLocusObservations {
        let mut changed = reference.to_vec();
        // A base the motif cannot show at that position, whatever the motif is, so the read
        // mismatches there and nowhere else.
        changed[changed_at] = if reference[changed_at] == b'G' {
            b'C'
        } else {
            b'G'
        };
        let tiled = |bases: &[u8]| {
            bases
                .iter()
                .zip(motif.iter().cycle())
                .filter(|(shown, tiled)| shown != tiled)
                .count()
        };
        assert_eq!(
            tiled(reference),
            0,
            "the clean reads must show a whole tract"
        );
        assert_eq!(
            tiled(&changed),
            1,
            "each mismatching read must differ from the tiled motif in exactly one base"
        );

        tract(
            start,
            reference,
            motif,
            vec![
                observation(reference, group, clean),
                observation(&changed, group, mismatching),
            ],
        )
    }

    /// **The rate a stratum returns is the rate its reads carried**: 30 mismatched bases in
    /// 10,000 compared comes back as 0.0030, exactly, with the bases as its warrant and marked
    /// as measured here.
    ///
    /// Exactly, and not within a tolerance, because there is nothing here to be approximate
    /// about: the two counters are integers and the rate is their quotient, which is why this
    /// fit is a division rather than a search. The maximum-likelihood property of that division
    /// is proven separately, against a grid of 100,000 rates, in `stratum_table.rs`.
    ///
    /// **The 500 reads are also the check that the read cap does not reach the bases.** The cap
    /// keeps twelve reads a locus for the *shape*; if it also thinned the base counters this
    /// stratum would report its rate over 240 bases — twelve reads of twenty — rather than
    /// 10,000, so the warrant beside every rate in a run would be 41.7 times too small.
    #[test]
    fn a_stratums_substitution_rate_is_the_rate_its_reads_carried() {
        let reference = b"ATATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract_at_a_known_rate(
            1_000, reference, b"AT", 0, 470, 30, 1,
        ));

        let rates = substitution_rates(&accumulators);
        let key = StratumKey {
            read_group: ReadGroupId(0),
            stratum: stratum(2, 10),
            ploidy: Ploidy::try_new(2).expect("two genome copies"),
        };
        let fitted = rates.get(&key).expect("that stratum was fitted");

        assert_eq!(
            fitted.value.get(),
            0.0030,
            "30 mismatched bases in 10,000 compared"
        );
        assert_eq!(
            fitted.observations, 10_000,
            "500 reads of 20 bases, none of them lost to the twelve-read cap on shapes"
        );
        assert_eq!(fitted.provenance, Provenance::FittedHere);
    }

    /// **Each stratum's rate comes from its own bases, and never from the read group's pool.**
    /// The two strata here are ten-fold apart — 0.0030 and 0.0300 — and a fit that pooled a read
    /// group's bases would report 0.0107 for both, which is neither and is a plausible enough
    /// number that nothing downstream would question it.
    ///
    /// Tracts differ in how much of their length is repeat and how well reads align inside them,
    /// so this is not a hypothetical difference: it is the reason the rate is fitted per stratum
    /// rather than per library.
    #[test]
    fn two_strata_of_one_library_keep_their_own_rates() {
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract_at_a_known_rate(
            1_000,
            b"ATATATATATATATATATAT",
            b"AT",
            0,
            470,
            30,
            1,
        ));
        accumulators.add_locus(&tract_at_a_known_rate(
            2_000,
            b"AAAAAAAA",
            b"A",
            0,
            380,
            120,
            3,
        ));

        let rates = substitution_rates(&accumulators);
        let rate_at = |period_bases: u8, repeats: u32| {
            rates
                .get(&StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: stratum(period_bases, repeats),
                    ploidy: Ploidy::try_new(2).expect("two genome copies"),
                })
                .expect("that stratum was fitted")
                .clone()
        };

        assert_eq!(rate_at(2, 10).value.get(), 0.0030);
        assert_eq!(rate_at(2, 10).observations, 10_000);
        assert_eq!(rate_at(1, 8).value.get(), 0.0300);
        assert_eq!(rate_at(1, 8).observations, 4_000);
    }

    /// **And each library's rate comes from its own reads**, at one and the same stratum. Two
    /// tracts of the same period and repeat count, each witnessed by one library, one of the two
    /// ten times noisier; pooling them would hand both the base-weighted mean and hide a library
    /// worth re-running.
    #[test]
    fn two_libraries_at_one_stratum_keep_their_own_rates() {
        let reference = b"ATATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract_at_a_known_rate(
            1_000, reference, b"AT", 0, 470, 30, 1,
        ));
        accumulators.add_locus(&tract_at_a_known_rate(
            2_000, reference, b"AT", 1, 200, 300, 1,
        ));

        let rates = substitution_rates(&accumulators);
        let rate_of = |group: u32| {
            rates
                .get(&StratumKey {
                    read_group: ReadGroupId(group),
                    stratum: stratum(2, 10),
                    ploidy: Ploidy::try_new(2).expect("two genome copies"),
                })
                .expect("that library was fitted")
                .value
                .get()
        };

        assert_eq!(rate_of(0), 0.0030);
        assert_eq!(rate_of(1), 0.0300, "300 mismatches in 10,000 bases");
    }

    /// **A stratum nothing was compared in has no rate, and is absent rather than zero.** The
    /// two are different claims — a stratum whose every read matched measures zero — and a zero
    /// invented here would be borrowed from by the strata around it and, at the end of the run,
    /// compared against the SNP/indel path's rate as though it had been measured.
    #[test]
    fn a_table_with_nothing_compared_yields_no_rate_at_all() {
        assert_eq!(substitution_rate_of(&StratumTable::default()), None);

        let empty = SsrAccumulators::new(diploid());
        assert!(substitution_rates(&empty).is_empty());
    }

    /// **A stratum that holds loci and compared no bases is absent, and that is the case worth
    /// building**, because an empty accumulator is absent under every implementation ever
    /// written: it has no strata to walk. Here the table exists and holds a locus, so an
    /// implementation that answered zero rather than nothing — or that decided on the locus count
    /// instead of the base count — is visible.
    ///
    /// The locus is a tract every read shows as **entirely deleted**: the reads witness the
    /// tract completely, so they are entered and the locus files a shape, and each of them shows
    /// no bases, so there is nothing to compare against the motif.
    #[test]
    fn a_stratum_with_loci_and_no_compared_bases_is_absent_rather_than_zero() {
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract(
            1_000,
            b"ATATATATATATATATATAT",
            b"AT",
            vec![observation(b"", 0, 6)],
        ));

        let key = StratumKey {
            read_group: ReadGroupId(0),
            stratum: stratum(2, 10),
            ploidy: Ploidy::try_new(2).expect("two genome copies"),
        };
        let table = accumulators.table_for(key).expect("the locus was filed");
        assert_eq!(table.loci(), 1, "the stratum holds a locus");
        assert_eq!(table.bases_compared(), 0, "and compared nothing");

        assert_eq!(substitution_rate_of(table), None);
        assert!(
            !substitution_rates(&accumulators).contains_key(&key),
            "a stratum with no compared bases has no rate, not a rate of zero"
        );
    }

    /// **A stratum whose every read matched reports a rate of zero, and reports it as measured.**
    /// This is the other half of the distinction above, and the half a floor or a
    /// drop-the-zeroes filter would quietly remove: zero is what a clean stratum measures, and a
    /// caller that never sees it cannot tell one from a stratum that was never fitted.
    ///
    /// **It also pins that thinness is not this fit's business.** 100 compared bases is a
    /// hundredth of what the fixtures above carry, and the rate still comes back, marked as
    /// fitted here. Deciding that a stratum has too little evidence — and saying where its number
    /// then came from — belongs to the borrowing step, which marks what it does; a gate here
    /// would drop the stratum with nothing recording that it had.
    #[test]
    fn a_stratum_whose_reads_all_matched_reports_a_measured_zero() {
        let reference = b"ATATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract(
            1_000,
            reference,
            b"AT",
            vec![observation(reference, 0, 5)],
        ));

        let fitted = substitution_rates(&accumulators)
            .get(&StratumKey {
                read_group: ReadGroupId(0),
                stratum: stratum(2, 10),
                ploidy: Ploidy::try_new(2).expect("two genome copies"),
            })
            .expect("a clean stratum is still a fitted stratum")
            .clone();

        assert_eq!(fitted.value.get(), 0.0, "every base matched");
        assert_eq!(fitted.observations, 100, "five reads of twenty bases");
        assert_eq!(
            fitted.provenance,
            Provenance::FittedHere,
            "thin evidence is still this stratum's own evidence"
        );
    }

    /// **Two ploidies of one stratum keep their own rates**, which is what putting the ploidy in
    /// the key is for at the fitting end rather than the counting end. The two loci here are the
    /// same period and repeat count and the same library, ten-fold apart in rate, and differ only
    /// in how many genome copies they sit on; under a key that dropped the ploidy one of the two
    /// rates would simply not be reported.
    #[test]
    fn two_ploidies_of_one_stratum_keep_their_own_rates() {
        let reference = b"ATATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(Arc::new(PloidyChangesAt {
            haploid_from: 5_000,
        }));
        accumulators.add_locus(&tract_at_a_known_rate(
            1_000, reference, b"AT", 0, 470, 30, 1,
        ));
        accumulators.add_locus(&tract_at_a_known_rate(
            9_000, reference, b"AT", 0, 200, 300, 1,
        ));

        let rates = substitution_rates(&accumulators);
        assert_eq!(rates.len(), 2, "one stratum at two ploidies is two fits");

        let rate_at = |copies: u8| {
            rates
                .get(&StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: stratum(2, 10),
                    ploidy: Ploidy::try_new(copies).expect("one or two genome copies"),
                })
                .expect("that ploidy was fitted")
                .value
                .get()
        };

        assert_eq!(rate_at(2), 0.0030, "the locus on two genome copies");
        assert_eq!(rate_at(1), 0.0300, "the locus on one");
    }

    // -----------------------------------------------------------------
    // The slippage search, per stratum.
    // -----------------------------------------------------------------

    /// `loci` tracts of ten reads each, of which `losing_reads` show one motif copy fewer than
    /// the reference. Every locus is identical, so the stratum is one entry standing for all of
    /// them — which is what a real stratum mostly is, and what makes these fits cheap enough to
    /// run in a unit test.
    fn stratum_losing_repeats(loci: u64, losing_reads: u32) -> SsrAccumulators {
        let reference = b"ATATATATATATATATATAT";
        let short = b"ATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(diploid());
        for locus in 0..loci {
            let mut observations = vec![observation(reference, 0, 10 - losing_reads)];
            if losing_reads > 0 {
                observations.push(observation(short, 0, losing_reads));
            }
            accumulators.add_locus(&tract(1_000 + locus * 100, reference, b"AT", observations));
        }
        accumulators
    }

    fn only_table(accumulators: &SsrAccumulators) -> (StratumKey, &StratumTable) {
        let mut strata = accumulators.strata();
        let one = strata.next().expect("a stratum");
        assert!(strata.next().is_none(), "exactly one stratum");
        one
    }

    /// **A stratum where every read sat at the reference length comes back barely slipping at
    /// all.** There is nothing for the level to explain, so the search runs it down to its own
    /// floor rather than settling on some middling value — which is what makes this a check on
    /// the search and not on the fixture.
    ///
    /// **The bound is 1e-3, and what it is set against is the other end of the range.** A search
    /// that read nothing from its cells does not hand back a starting point — every axis is
    /// line-searched over its whole range at every sweep — it walks to wherever ties resolve,
    /// which for this golden section is the **top** of the level's range. Measured on a stratum
    /// the model cannot score at all, that is 0.5976. So the failure this bound catches is a
    /// level six hundred times its own size, and the four starting levels (3e-4, 1e-4, 3.3e-5,
    /// 3e-5) are not what it is guarding against.
    #[test]
    fn a_stratum_where_nothing_slipped_is_fitted_at_almost_no_slippage() {
        let accumulators = stratum_losing_repeats(200, 0);
        let (key, table) = only_table(&accumulators);

        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("a stratum whose reads sat on the whole-repeat grid");

        assert!(
            fit.model.slip_rate.get() < 1e-3,
            "nothing slipped, so the level should run down to its floor: {}",
            fit.model
        );
        assert_eq!(fit.slipped_reads, 0);
        assert_eq!(fit.loci, 200);
    }

    /// **A stratum where one read in ten lost a repeat is fitted at about one in ten, losing.**
    /// Both halves matter and they fail differently: a level near 0.1 says the search reads its
    /// cells, and a gain share below one half says it reads which *way* they moved — the one
    /// property the whole design exists to protect, and the one the estimator this replaces gets
    /// backwards.
    ///
    /// **The truth here is exactly 0.1** and needs no tolerance argument: with the genotype
    /// frequencies free, the cheapest explanation of every locus showing nine reads at the
    /// reference and one a copy short is that every locus is homozygous reference and one read in
    /// ten slipped down. **The bound is ±2% of the level**, which is twice the 1% resolution the
    /// coarse search setting resolves a rate to — so it is the search's own step that sets it and
    /// not a guess. Measured, the fit returns 0.09986, 0.14% below the truth.
    #[test]
    fn a_stratum_that_loses_a_repeat_in_one_read_of_ten_is_fitted_at_about_that() {
        let accumulators = stratum_losing_repeats(500, 1);
        let (key, table) = only_table(&accumulators);

        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("a stratum whose reads sat on the whole-repeat grid");

        let level = fit.model.slip_rate.get();
        assert!(
            (0.098..=0.102).contains(&level),
            "one read in ten sat a copy short: {}",
            fit.model
        );
        assert!(
            fit.model.gain_share.get() < 0.1,
            "every read that moved lost a copy, so gains cannot be the commoner direction: {}",
            fit.model
        );
        assert_eq!(
            fit.slipped_reads, 500,
            "one of ten reads at each of 500 loci"
        );
        assert_eq!(
            fit.genotype_frequencies.len(),
            91,
            "13 allele lengths at ploidy 2"
        );
    }

    /// **Every start is reported with where it began, where it ended and what it scored, best
    /// first** — and the spread across them is reported whatever it says. An answer with no
    /// spread beside it cannot be told apart from a search that never looked.
    #[test]
    fn the_search_reports_all_four_starts_best_scoring_first() {
        let accumulators = stratum_losing_repeats(500, 1);
        let (key, table) = only_table(&accumulators);

        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("a stratum whose reads sat on the whole-repeat grid");

        assert_eq!(fit.starts_tried.len(), 4, "four starts, all reported");
        for pair in fit.starts_tried.windows(2) {
            assert!(
                pair[0].log_likelihood >= pair[1].log_likelihood,
                "the starts are ordered best-scoring first: {:?}",
                fit.starts_tried
            );
        }
        // **The four starting points disagree about all three parameters, not only the level**,
        // which is the trap the SNP/indel path's inbreeding fit fell into: starts that share one
        // guess at a nuisance axis agree at the end for that reason alone.
        for axis in [
            |model: SlippageModel| model.slip_rate.get(),
            |model: SlippageModel| model.gain_share.get(),
            |model: SlippageModel| model.step_decay.get(),
        ] {
            let mut values: Vec<f64> = fit.starts_tried.iter().map(|s| axis(s.from)).collect();
            values.sort_by(f64::total_cmp);
            values.dedup();
            assert_eq!(values.len(), 4, "four distinct starting values: {values:?}");
        }
        let levels: Vec<f64> = fit
            .starts_tried
            .iter()
            .map(|start| start.from.slip_rate.get())
            .collect();
        let widest = levels.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            / levels.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            (widest - 10.0).abs() < 1e-9,
            "the highest starting level is ten times the lowest — 3x and 0.3x the stratum's own \
             share of moved reads: {levels:?}"
        );
        // **And that share is this stratum's own, not a constant.** The four starts are placed at
        // 3, 1, 1/3 and 0.3 times it, so the highest is three times the share — 0.3 here, where
        // one read in ten moved. A search started from a fixed level would still show four
        // distinct values ten times apart and would be reading nothing from the table.
        let share = fit.slipped_reads as f64 / fit.scored_reads as f64;
        assert!(
            (levels.iter().copied().fold(f64::NEG_INFINITY, f64::max) - 3.0 * share).abs() < 1e-12,
            "the starts are placed around this stratum's share of moved reads, {share}: {levels:?}"
        );
        assert!(
            fit.start_spread.is_finite() && fit.start_spread >= 1.0,
            "the spread is a ratio of the highest level reached to the lowest: {}",
            fit.start_spread
        );
        assert_eq!(
            fit.log_likelihood, fit.starts_tried[0].log_likelihood,
            "the fit's score is the best start's"
        );
    }

    /// **Starts that landed far apart are refused rather than reported**, and the refusal names
    /// whose fit it was: on a cohort of hundreds of samples, several hundred fits a sample, a
    /// message naming only the quantity locates nothing.
    ///
    /// Built by hand rather than searched for, because what is under test is the judgement and
    /// not the search: this is the one place the limit is applied, and a fixture that had to find
    /// a genuinely unidentified stratum would be testing how hard that is to construct.
    #[test]
    fn a_fit_whose_starts_landed_far_apart_is_refused_and_says_whose_it_was() {
        let far_apart = StratumSlippageFit {
            start_spread: 1.5,
            ..fitted_at(SlippageModel::try_new(0.02, 0.17, 0.065).expect("a model"))
        };
        let key = StratumKey {
            read_group: ReadGroupId(3),
            stratum: stratum(2, 6),
            ploidy: Ploidy::try_new(2).expect("two genome copies"),
        };

        let refused = starts_must_agree(&far_apart, "SL_landrace_07", key)
            .expect_err("1.5 is beyond the limit");
        let message = refused.to_string();

        assert!(message.contains("SL_landrace_07"), "the sample: {message}");
        assert!(message.contains("read group 3"), "the library: {message}");
        assert!(
            message.contains("period 2, 6 repeats"),
            "which stratum: {message}"
        );
        assert!(message.contains("1.5"), "how far apart: {message}");
        assert!(
            message.contains("did not settle"),
            "and that this is not too little data: {message}"
        );

        let close_enough = StratumSlippageFit {
            start_spread: START_AGREEMENT_LIMIT,
            ..far_apart.clone()
        };
        assert!(
            starts_must_agree(&close_enough, "SL_landrace_07", key).is_ok(),
            "the limit itself is agreement, not disagreement"
        );

        // **And a fit whose inner climbs ran out of passes is not refused**, which is the other
        // half of the decision this function makes. The surface those climbs walk is concave, so
        // one that ran out did not find a wrong summit; it ran out of time on the right one. A
        // stratum here has 66 to 91 genotypes against the sibling path's three, so most real
        // strata are expected to arrive with this false, and refusing them would end the run.
        let ran_out = StratumSlippageFit {
            every_climb_settled: false,
            start_spread: 1.0,
            ..far_apart
        };
        assert!(
            starts_must_agree(&ran_out, "SL_landrace_07", key).is_ok(),
            "a climb that ran out of passes is reported, not refused"
        );
    }

    /// A fit with everything but the field a test is about.
    fn fitted_at(model: SlippageModel) -> StratumSlippageFit {
        StratumSlippageFit {
            model,
            genotype_frequencies: vec![1.0],
            allele_support: SmallVec::from_slice(&[WholeRepeatOffset(0)]),
            log_likelihood: -12.0,
            starts_tried: SmallVec::from_vec(vec![SlippageStart {
                from: model,
                reached: model,
                log_likelihood: -12.0,
            }]),
            start_spread: 1.0,
            every_climb_settled: true,
            termination: FitTermination {
                iterations: 3,
                converged: true,
            },
            slipped_reads: 0,
            guard_reads: 0,
            scored_reads: 10,
            loci: 1,
        }
    }

    /// **The walk fits every stratum the accumulator holds, under the accumulator's own keys.**
    /// Two strata of one library here, one of them slipping and one not, so a walk that fitted
    /// one stratum and copied its answer to the rest would be visible.
    #[test]
    fn the_walk_fits_every_stratum_under_its_own_key() {
        let mut accumulators = stratum_losing_repeats(300, 1);
        // A second stratum — a mononucleotide run — whose reads never left the reference.
        for locus in 0..300u64 {
            accumulators.add_locus(&tract(
                500_000 + locus * 100,
                b"AAAAAAAA",
                b"A",
                vec![observation(b"AAAAAAAA", 0, 10)],
            ));
        }

        let fits =
            fit_slippage_by_stratum(&accumulators, "SL_landrace_07", SearchPrecision::fast())
                .expect("both strata settled");

        let keys: Vec<StratumKey> = fits.keys().copied().collect();
        let from_accumulator: Vec<StratumKey> = accumulators.strata().map(|(key, _)| key).collect();
        assert_eq!(keys, from_accumulator, "one fit per stratum, same keys");

        let level_at = |period_bases: u8, repeats: u32| {
            fits[&StratumKey {
                read_group: ReadGroupId(0),
                stratum: stratum(period_bases, repeats),
                ploidy: Ploidy::try_new(2).expect("two genome copies"),
            }]
                .model
                .slip_rate
                .get()
        };

        assert!(
            level_at(2, 10) > 0.05,
            "the dinucleotide stratum lost a repeat in one read of ten"
        );
        assert!(
            level_at(1, 8) < 1e-3,
            "the mononucleotide stratum never left the reference"
        );
    }

    /// **A stratum that says nothing is refused by the walk, and the search's own diagnostics are
    /// what refuse it.** One locus holding one read fifteen copies longer than the reference is
    /// consistent with almost any slippage: the four starts reach levels 67 times apart, far past
    /// the [`START_AGREEMENT_LIMIT`] of 1.06, so what the search returned is where it stopped.
    ///
    /// **This is also the only fixture here where the inner climbs run out of passes**, and it is
    /// asserted, because that flag decides nothing on its own and would otherwise be a field
    /// nobody reads: the fit is returned and the walk does not refuse it for that.
    ///
    /// The search is run at a deliberately coarse setting. Not to change the answer — the spread
    /// is 193 at the ordinary setting — but because a fit whose climbs never settle costs 10,000
    /// passes per candidate, and this fixture takes 22 seconds at the ordinary setting against 4
    /// at this one.
    #[test]
    fn a_stratum_whose_starts_land_far_apart_is_refused_by_the_walk() {
        let reference = b"ATATATATATATATATATAT";
        let far_longer = b"ATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let coarse = SearchPrecision {
            tolerance: 0.1,
            max_axis_steps: 4,
            max_sweeps: 1,
        };

        let mut accumulators = SsrAccumulators::new(diploid());
        accumulators.add_locus(&tract(
            1_000,
            reference,
            b"AT",
            vec![observation(far_longer, 0, 1)],
        ));

        let (key, table) = only_table(&accumulators);
        let fit = fit_slippage(table, key.stratum, key.ploidy, coarse).expect("a scoreable read");
        assert!(
            fit.start_spread > START_AGREEMENT_LIMIT,
            "one read says almost nothing about the level: {}",
            fit.start_spread
        );
        assert!(
            !fit.every_climb_settled,
            "the climb over 91 genotypes on one locus runs out of passes"
        );

        // **This is also the only fixture where the four starts score differently**, so it is the
        // only place "best-scoring first" and "the answer is the best start's" have any content.
        // Where every start reaches the same point — which is what the settled fixtures do, and
        // what the search working looks like — both assertions compare a number with itself.
        let scores: Vec<f64> = fit
            .starts_tried
            .iter()
            .map(|start| start.log_likelihood)
            .collect();
        assert!(
            scores.first() > scores.last(),
            "the starts really do disagree here, and are ordered best first: {scores:?}"
        );
        assert_eq!(
            fit.model, fit.starts_tried[0].reached,
            "the answer is where the best-scoring start ended"
        );
        assert_eq!(
            fit.log_likelihood, scores[0],
            "and its score is that start's score"
        );
        let reached: Vec<f64> = fit
            .starts_tried
            .iter()
            .map(|start| start.reached.slip_rate.get())
            .collect();
        let widest = reached.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            / reached.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            (fit.start_spread - widest).abs() < 1e-9,
            "the spread is the highest level reached over the lowest: {} against {widest}",
            fit.start_spread
        );
        assert!(
            starts_must_agree(&fit, "SL_landrace_07", key).is_err(),
            "and a fit that did not settle is refused"
        );

        let refused = fit_slippage_by_stratum(&accumulators, "SL_landrace_07", coarse)
            .expect_err("the walk refuses it too, rather than reporting where it stopped");
        let message = refused.to_string();
        assert!(message.contains("SL_landrace_07"), "{message}");
        assert!(message.contains("did not settle"), "{message}");
    }

    /// **A tract of four copies cannot carry an allele six copies shorter**, so its support is
    /// clipped at the low end: eleven allele lengths and 66 genotypes, against thirteen and 91
    /// for a tract of six copies or more. The support is emitted beside the frequencies because
    /// after a merge it cannot be rebuilt from the stratum — the merged model is built from the
    /// **lower** of the two repeat counts.
    #[test]
    fn a_short_tract_is_fitted_over_a_clipped_support() {
        let reference = b"ATATATAT";
        let mut accumulators = SsrAccumulators::new(diploid());
        for locus in 0..50u64 {
            accumulators.add_locus(&tract(
                1_000 + locus * 100,
                reference,
                b"AT",
                vec![observation(reference, 0, 6)],
            ));
        }

        let (key, table) = only_table(&accumulators);
        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("a stratum whose reads sat on the whole-repeat grid");

        assert_eq!(
            fit.allele_support.as_slice(),
            &[
                WholeRepeatOffset(-4),
                WholeRepeatOffset(-3),
                WholeRepeatOffset(-2),
                WholeRepeatOffset(-1),
                WholeRepeatOffset(0),
                WholeRepeatOffset(1),
                WholeRepeatOffset(2),
                WholeRepeatOffset(3),
                WholeRepeatOffset(4),
                WholeRepeatOffset(5),
                WholeRepeatOffset(6),
            ],
            "four copies below the reference and six above it"
        );
        assert_eq!(
            fit.genotype_frequencies.len(),
            66,
            "eleven allele lengths make 66 unordered pairs"
        );
    }

    /// **A stratum on one genome copy is fitted over single alleles, not pairs** — thirteen
    /// frequencies against a diploid stratum's 91. The ploidy reaches the fit through the key its
    /// table is stored under, so a fit that assumed two copies would score a haploid locus
    /// against genotypes it cannot have.
    #[test]
    fn a_haploid_stratum_is_fitted_over_single_alleles() {
        let reference = b"ATATATATATATATATATAT";
        let mut accumulators = SsrAccumulators::new(Arc::new(PloidyChangesAt {
            haploid_from: 5_000,
        }));
        for locus in 0..50u64 {
            accumulators.add_locus(&tract(
                9_000 + locus * 100,
                reference,
                b"AT",
                vec![observation(reference, 0, 6)],
            ));
        }

        let (key, table) = only_table(&accumulators);
        assert_eq!(key.ploidy.get(), 1, "the map puts these loci on one copy");

        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("a stratum whose reads sat on the whole-repeat grid");

        assert_eq!(
            fit.genotype_frequencies.len(),
            13,
            "thirteen allele lengths, one copy each"
        );
    }

    /// **Reads that moved by a whole number of copies, reads that moved by something else, and
    /// reads the model can score are three different counts**, and this fixture separates all
    /// three: 800 loci where nothing moved, and 200 where one read of ten sits a copy short and
    /// another carries a single-base deletion inside the tract.
    ///
    /// Why it matters: the count that gates the direction share and the fall-off is the first of
    /// the three, and those two are estimated only from whole-copy movement. Pooling the second
    /// into it overstates the evidence behind them by the stratum's guard share, which spec §5
    /// measures at up to 58.5% of the moved reads.
    ///
    /// **The two entries carry different locus counts on purpose.** Each entry's counts are per
    /// locus, so a walk that forgot to weigh them by how many loci showed that shape would report
    /// 1, 1 and 19 here instead of 200, 200 and 9,800.
    #[test]
    fn the_reads_that_moved_are_counted_apart_from_the_reads_that_moved_off_the_grid() {
        let reference = b"ATATATATATATATATATAT";
        let short_by_a_copy = b"ATATATATATATATATAT";
        let short_by_a_base = b"ATATATATATATATATATA";

        let mut accumulators = SsrAccumulators::new(diploid());
        for locus in 0..800u64 {
            accumulators.add_locus(&tract(
                1_000 + locus * 100,
                reference,
                b"AT",
                vec![observation(reference, 0, 10)],
            ));
        }
        for locus in 0..200u64 {
            accumulators.add_locus(&tract(
                500_000 + locus * 100,
                reference,
                b"AT",
                vec![
                    observation(reference, 0, 8),
                    observation(short_by_a_copy, 0, 1),
                    observation(short_by_a_base, 0, 1),
                ],
            ));
        }

        let (key, table) = only_table(&accumulators);
        assert_eq!(table.entry_count(), 2, "two shapes, unevenly common");

        let fit = fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast())
            .expect("most reads sat on the whole-repeat grid");

        assert_eq!(fit.slipped_reads, 200, "one whole copy short, at 200 loci");
        assert_eq!(fit.guard_reads, 200, "one base short, at the same 200 loci");
        assert_eq!(
            fit.scored_reads, 9_800,
            "8,000 reads at 800 loci plus nine of ten at 200 more"
        );
        assert_eq!(fit.loci, 1_000);
    }

    /// **A stratum no read of which sits on the whole-repeat grid is not fitted at all**, and
    /// that is a guard rather than a tidiness. The model scores only whole-repeat movement, so
    /// every genotype and every candidate score alike on such a stratum; measured, the search
    /// then returns a slippage level of **0.5976** — the top of its range — with the four starts
    /// agreeing to 1.00, so the diagnostic that guards this fit cannot see it. It has as many
    /// loci as any other stratum, so the borrowing step's floor cannot see it either.
    #[test]
    fn a_stratum_no_read_of_which_sits_on_the_grid_is_not_fitted() {
        let reference = b"ATATATATATATATATATAT";
        let short_by_a_base = b"ATATATATATATATATATA";

        let mut accumulators = SsrAccumulators::new(diploid());
        for locus in 0..500u64 {
            accumulators.add_locus(&tract(
                1_000 + locus * 100,
                reference,
                b"AT",
                vec![observation(short_by_a_base, 0, 10)],
            ));
        }

        let (key, table) = only_table(&accumulators);
        assert_eq!(
            table.loci(),
            500,
            "the stratum is not thin — it is uninformative"
        );

        assert_eq!(
            fit_slippage(table, key.stratum, key.ploidy, SearchPrecision::fast()),
            None,
            "nothing here is on the grid the model scores"
        );
        assert!(
            fit_slippage_by_stratum(&accumulators, "SL_landrace_07", SearchPrecision::fast())
                .expect("no fit means nothing to disagree about")
                .is_empty(),
            "and the walk leaves it out rather than reporting a level for it"
        );
    }
}
