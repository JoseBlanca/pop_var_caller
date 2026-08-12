//! What a read at a repeat tract does to the length it shows, and the search that fits it.
//!
//! **Three numbers, and they answer different questions** (`spec/parameter_prepass_ssr.md`
//! §1.1, §3):
//!
//! - **how often a read slips at all** — the level, and the quantity strata are compared on:
//!   9 reads in 10,000 below four repeats against 2 in 100 at six or more;
//! - **which way it slips**, which is strongly asymmetric — a read at a tomato dinucleotide
//!   is 4.9 times as likely to have lost a repeat as to have gained one, and the imbalance
//!   grows with the motif period;
//! - **how far it slips when it does**, which decays — of the reads that slipped by one
//!   repeat, about 7 in 100 slipped by two instead in tomato, about 10 in 100 in human —
//!   and decays the **same way in both directions**.
//!
//! **The direction is asymmetric and the distance is not**, and both halves of that were
//! decided against measurement rather than assumed. The asymmetry is large enough that a
//! model without it cannot describe the data; the gaining arm's decay rests on 3 to 13 reads
//! above dinucleotides, so a free parameter there would fit counting noise rather than a
//! difference — the four rows measured differ by 1.5, 0.9, 1.6 and 0.5 standard errors, and
//! all four point the same way, which pools to about 2. The finding is "no difference we can
//! afford to fit", not "no difference".
//!
//! **The fourth number this path fits is not slippage at all.** A read can also be misread
//! at fixed length, and that per-base substitution rate is a division — mismatched bases over
//! bases compared — not an axis of any search (§4.1).
//!
//! A4 landed the three parameters; D1 lands the kernel they describe —
//! [`SlippageModel::probability_of_slipping_by`], which turns them into how often a read shows
//! exactly `d` whole motif copies more than the allele it came off.

use std::fmt;

use smallvec::SmallVec;

use crate::genetics::lgamma;
use crate::ng::parameter_estimation::fitting::multistart::SearchableNoise;
use crate::ng::parameter_estimation::fitting::{NoiseModel, WeightedCell};
use crate::ng::types::{DomainError, Ploidy, checked_probability};

use super::stratum_table::StratumCell;
use super::{
    ALLELE_OFFSET_LIMIT, OFFSET_BUCKETS, RepeatCount, WholeRepeatOffset, allele_support, bucket_of,
};

/// How often a read shows a length other than its allele's — **the level**, and the number a
/// stratum is chosen by: it spans twenty-two-fold across repeat counts within one dataset,
/// which is why strata exist at all (`spec/parameter_prepass_ssr.md` §4).
///
/// A probability in `[0, 1]`. Zero is a real value here rather than a degenerate one — the
/// bottom of the repeat range sits at 0.00091 — and it is also the boundary a finite
/// stratum's estimate piles up against, which is why the count of reads that actually slipped
/// travels beside every fitted level (§4.5).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct SlipRate(f64);

impl SlipRate {
    /// The only constructor. A rate that is not a probability in `[0, 1]` is rejected rather
    /// than coerced.
    pub fn try_new(rate: f64) -> Result<Self, DomainError> {
        checked_probability(rate, DomainError::SlipRate).map(Self)
    }

    #[inline]
    #[must_use]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// Of the reads that slipped, the share that **gained** repeats rather than losing them.
///
/// A probability in `[0, 1]`, and **far from a half**: 0.17 at tomato dinucleotides, where a
/// read is 4.9 times as likely to have lost a repeat as to have gained one
/// (`spec/parameter_prepass_ssr.md` §3). Half would mean the tract slips symmetrically, which
/// no dataset measured here does.
///
/// **It is also the parameter that collapses when the estimator is wrong**, which is what
/// makes it a diagnostic as well as a parameter. Production's estimator, which pools reads
/// from loci that passed a confident-genotype gate, goes past collapse to inversion — it
/// reports gains as marginally *more* common than losses. Centring each locus on its own
/// modal observed length and scoring it as though the origin were fixed costs about the same
/// **size**: the share comes back at 0.48 against a truth of 0.17, so losses lead by 1.1-fold
/// where they truly lead by 4.9-fold (§4.1). One arrives from thresholding and the other from
/// a keying choice; neither leaves the asymmetry the model exists to carry.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct SlipGainShare(f64);

impl SlipGainShare {
    /// The only constructor. A share that is not a probability in `[0, 1]` is rejected rather
    /// than coerced.
    pub fn try_new(share: f64) -> Result<Self, DomainError> {
        checked_probability(share, DomainError::SlipGainShare).map(Self)
    }

    #[inline]
    #[must_use]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// Given that a read slipped by one repeat, how often it slipped by another — the geometric
/// fall-off, **one number shared by both directions**.
///
/// A probability in `[0, 1]`. Read `0.065` as *of every 100 reads that slipped by one repeat,
/// about 7 slipped by two instead*: 5,072 reads one repeat short against 329 two repeats
/// short, at tomato homopolymers (`spec/parameter_prepass_ssr.md` §3). **Those counts are
/// measured from each unit's own modal observed length, not from the reference** — the
/// origin §4.1 rejects for the accumulator — which is fine for a ratio between two distances
/// and would not be for the level.
///
/// **The value does not transfer between datasets and the structure does** — about 10 in 100
/// in human against about 7 in tomato — so it is fitted per stratum rather than assumed.
/// It is also the parameter that starves first: holding it to 6% of itself takes about 4,000
/// slipped reads, against about 1,400 for the direction share, which is why a stratum can
/// keep its own level and borrow these two (§4.5).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct SlipStepDecay(f64);

impl SlipStepDecay {
    /// The only constructor. A decay that is not a probability in `[0, 1]` is rejected rather
    /// than coerced.
    pub fn try_new(decay: f64) -> Result<Self, DomainError> {
        checked_probability(decay, DomainError::SlipStepDecay).map(Self)
    }

    #[inline]
    #[must_use]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// How a read's repeat count moves away from its allele's: how often, which way, and how far.
///
/// **Three types and not one shared probability**, though all three are fractions in
/// `[0, 1]`: one type would let the gain share be handed to something expecting the level and
/// compile (`arch/parameter_prepass_ssr.md` §2.4). **And the two are not reliably far apart,
/// which is the argument rather than against it** — the gain share is 0.17 everywhere, while
/// the level runs from 0.00091 at the bottom of the repeat range to 0.150 at tomato
/// dinucleotides of 12 to 15 repeats. At the bottom a transposition is a 190-fold error and
/// obvious; at the top the two numbers sit within 1.1-fold of each other, and nothing about
/// the answer would look wrong.
///
/// **The fourth number a stratum emits is not in here.** The per-base substitution rate
/// belongs to the composition channel, which factorises out of this one exactly — a read's
/// mismatch count is binomial whatever length it showed — so it is fitted by a division and
/// never enters this search (`spec/parameter_prepass_ssr.md` §4.1).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SlippageModel {
    /// How often a read shows a length other than its allele's.
    pub slip_rate: SlipRate,
    /// Of the reads that slipped, the share that gained repeats.
    pub gain_share: SlipGainShare,
    /// Of the reads that slipped, the chance of a second step given a first.
    pub step_decay: SlipStepDecay,
}

impl SlippageModel {
    /// Build a slippage model from three already-checked rates.
    #[must_use]
    pub fn new(slip_rate: SlipRate, gain_share: SlipGainShare, step_decay: SlipStepDecay) -> Self {
        Self {
            slip_rate,
            gain_share,
            step_decay,
        }
    }

    /// The same, from three plain fractions — the door a search or a starting point comes
    /// through, where all three arrive together.
    ///
    /// # Errors
    ///
    /// The first of the three that is not a probability in `[0, 1]`, named as the quantity it
    /// was offered for.
    pub fn try_new(slip_rate: f64, gain_share: f64, step_decay: f64) -> Result<Self, DomainError> {
        Ok(Self::new(
            SlipRate::try_new(slip_rate)?,
            SlipGainShare::try_new(gain_share)?,
            SlipStepDecay::try_new(step_decay)?,
        ))
    }

    /// **How often a read shows exactly `whole_repeats` motif copies more than the allele it
    /// came off** — positive for a read that gained copies, negative for one that lost them,
    /// zero for a read that did not slip. The whole noise model in one function.
    ///
    /// Read it as three questions asked in order, which is what the three parameters are
    /// (`spec/parameter_prepass_ssr.md` §3). Did the read slip at all — no, with probability
    /// `1 − slip_rate`. If it did, which way — up with probability `gain_share`, down
    /// otherwise, and that split is far from even: at tomato dinucleotides a read is 4.9 times
    /// as likely to have lost a copy as gained one. And how far — one copy usually, two with
    /// probability `step_decay`, three with `step_decay` again, a geometric fall-off **shared
    /// by both directions**.
    ///
    /// **The fall-off is truncated at [`MAX_SLIP_STEP`] copies and renormalised over what is
    /// left**, so the kernel is a proper distribution: it sums to one over
    /// `−MAX_SLIP_STEP ..= MAX_SLIP_STEP` at every parameter setting, which is the first of the
    /// three algebraic gates the design puts before any fitting (§10.1). Truncating without
    /// renormalising would quietly lose mass instead — at the fall-offs measured on real data,
    /// 7 to 12 reads in 100 taking a second step, the mass beyond eight copies is under 6e-8 of
    /// the slipped reads, so the renormalisation changes no measured number and is what keeps
    /// the gate an identity rather than an approximation.
    ///
    /// **The renormaliser is written as the geometric's partial sum, `1 + f + … + f^{S−1}`, and
    /// not as the closed form `(1 − f)/(1 − f^S)` the reference implementation uses**
    /// ([`examples/shared/stutter_model.rs`](../../../../../examples/shared/stutter_model.rs)).
    /// The two are the same algebra and they are not the same arithmetic. The closed form
    /// subtracts two nearly equal numbers as the fall-off approaches one, so its relative error
    /// runs at about `1e-16 / (S·(1 − f))`: at `1 − f = 1e-9` the kernel sums to 0.9999999965,
    /// which is 3,500 times the tolerance the sums-to-one gate is asserted at, and the gate
    /// stops being an identity over a band of legal fall-offs. The partial sum has no
    /// cancellation in it, and it also has the `f = 1` case built in rather than branching to
    /// it — the sum is then `S`, so each distance gets `1/S`.
    ///
    /// **A fall-off of exactly one is a real input**: [`SlipStepDecay`] accepts 1.0, and the
    /// closed form is `0/0` there. The limit is a *uniform* distance — every distance from one
    /// to [`MAX_SLIP_STEP`] equally likely — and that is what this returns. Getting a `NaN`
    /// instead would not fail loudly: it would reach the likelihood, and the searches in this
    /// crate pick their maximum with `total_cmp`
    /// ([`coupled_fit.rs:1492`](../generic/coupled_fit.rs)), which ranks `NaN` above every
    /// finite score — so a fall-off of one would be *selected* rather than skipped, and
    /// reported with a `NaN` likelihood beside it.
    #[must_use]
    pub fn probability_of_slipping_by(&self, whole_repeats: i32) -> f64 {
        if whole_repeats == 0 {
            return 1.0 - self.slip_rate.get();
        }
        // Range-checked **before** the absolute value rather than after it, because
        // `i32::MIN.abs()` overflows: it panics in a debug build and wraps to `i32::MIN` in a
        // release one, where a negative `copies` slips past a `copies > MAX_SLIP_STEP` guard and
        // the function charges a slip of two billion copies whatever a one-copy slip costs.
        if !(-MAX_SLIP_STEP..=MAX_SLIP_STEP).contains(&whole_repeats) {
            return 0.0;
        }
        let copies = whole_repeats.abs();

        let direction = if whole_repeats > 0 {
            self.gain_share.get()
        } else {
            1.0 - self.gain_share.get()
        };

        let decay = self.step_decay.get();
        // The geometric over `1..=MAX_SLIP_STEP`, renormalised so the truncation loses no mass.
        // `1 + f + … + f^{S−1}`, for the cancellation reason on this function's doc comment.
        let normaliser: f64 = (0..MAX_SLIP_STEP).map(|step| decay.powi(step)).sum();
        let distance = decay.powi(copies - 1) / normaliser;

        self.slip_rate.get() * direction * distance
    }
}

/// The furthest a read is allowed to slip, in whole motif copies. The distance kernel is
/// renormalised over `1..=this`, so it stays a proper distribution
/// ([`SlippageModel::probability_of_slipping_by`]).
///
/// **Eight, and what it costs is below anything that could be measured.** At the fall-offs real
/// data shows — 7 reads in 100 taking a second step in tomato, 10 to 12 in 100 in human — the
/// mass beyond eight copies runs from 3.2e-10 at the bottom of that range to 5.2e-8 at the top
/// (`spec/parameter_prepass_ssr.md` §3's largest measured fall-off, 0.123). At the highest
/// slippage level any stratum reaches — 15.0%, at tomato dinucleotides of 12 to 15 repeats —
/// that is 7.9e-9 of all reads. The renormalisation hands that mass back to the distances
/// inside the range rather than dropping it.
///
/// **It is not the same width as anything else in this module and must not be confused with
/// them.** [`OFFSET_HALF_RANGE`](super::OFFSET_HALF_RANGE) is how far from the reference an
/// entry *records* a read (4); [`ALLELE_OFFSET_LIMIT`](super::ALLELE_OFFSET_LIMIT) is how far
/// from the reference the fit may place an *allele* (6); this is how far a read may slip from
/// **its own allele**, which is a distance between two different pairs of things.
pub const MAX_SLIP_STEP: i32 = 8;

/// How many allele lengths a stratum's fit can ever place mass on:
/// `±ALLELE_OFFSET_LIMIT` around the reference, so 13. A stratum's own support is this or
/// shorter — it clips at the low end, because an allele cannot be shorter than nothing — so
/// this is the ceiling that lets one locus be scored without allocating
/// ([`allele_support`](super::allele_support)).
const MAX_ALLELE_LENGTHS: usize = 2 * ALLELE_OFFSET_LIMIT as usize + 1;

/// How many distances a read can slip, counting both directions and standing still:
/// `−MAX_SLIP_STEP ..= MAX_SLIP_STEP`, so 17.
const SLIP_DISTANCES: usize = 2 * MAX_SLIP_STEP as usize + 1;

/// [`SlippageModel::probability_of_slipping_by`] at every distance at once, indexed by
/// `slip + MAX_SLIP_STEP`.
///
/// **Evaluated once per cell rather than once per allele.** The kernel renormalises its truncation
/// over eight powers of the fall-off on every call, and a cell's score walks it once for each of
/// up to thirteen allele lengths — so computing it here turns 1,872 `powi` calls per cell into
/// 144. Bit-identical: the same distances, the same values, accumulated into the same buckets in
/// the same order.
fn slip_kernel(noise: &SlippageModel) -> [f64; SLIP_DISTANCES] {
    let mut kernel = [0.0; SLIP_DISTANCES];
    for (index, share) in kernel.iter_mut().enumerate() {
        let slip = i32::try_from(index).expect("seventeen distances") - MAX_SLIP_STEP;
        *share = noise.probability_of_slipping_by(slip);
    }
    kernel
}

/// **How a repeat tract's reads land, given the alleles the locus carries** — the STR path's
/// side of the one seam step 4 has (`arch/parameter_prepass_ssr.md` §3).
///
/// It holds one thing: **the allele lengths this stratum's fit may place mass on**, as
/// whole-repeat offsets from the reference tract length. That belongs to the model rather than to
/// a cell because it is a property of the stratum — a tract of four reference copies cannot carry
/// an allele six copies shorter, so its support is 11 lengths where a tract of six or more has
/// the full 13, and the genotype count follows from that.
///
/// **One model per stratum, therefore, and not one per sample**, which is the difference from the
/// sibling path's stateless `SubstitutionNoiseModel`: there a genotype is a count of alternative
/// copies and the same three exist at every site, while here a genotype is a tuple of *lengths*
/// and which lengths exist is what a stratum is.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SsrNoiseModel {
    /// The allele lengths this stratum's genotypes are drawn from, ascending. Never empty.
    allele_support: SmallVec<[WholeRepeatOffset; MAX_ALLELE_LENGTHS]>,
}

impl SsrNoiseModel {
    /// The model for a stratum of tracts holding this many reference copies.
    ///
    /// **The only constructor, and it takes the stratum rather than a list of lengths**, so the
    /// low-end clip cannot be forgotten at a call site: a support reaching below an empty tract
    /// would let the fit explain a locus with an allele of negative length.
    ///
    /// **When strata are merged, the model is built from the lowest repeat count in the merge**
    /// (§4.3's monotonicity walk pools two neighbouring repeat counts and refits). That is the
    /// intersection of their supports rather than the union, and it is the right one for the same
    /// reason the clip exists: the shorter tract cannot carry the longer one's shortest alleles.
    #[must_use]
    pub fn for_stratum(reference_repeats: RepeatCount) -> Self {
        Self {
            allele_support: SmallVec::from_vec(allele_support(reference_repeats)),
        }
    }

    /// The allele lengths this stratum's genotypes are drawn from, ascending — as offsets from
    /// the reference tract length.
    #[inline]
    #[must_use]
    pub fn allele_support(&self) -> &[WholeRepeatOffset] {
        &self.allele_support
    }

    /// **How a single read coming off one allele lands across the buckets**, and the one place
    /// the end-bucket rule lives.
    ///
    /// The read slips by `d` with [`SlippageModel::probability_of_slipping_by`] and shows a tract
    /// `allele + d` copies from the reference length, which the entry records in the bucket that
    /// holds that offset. **The two end buckets absorb everything past them, so their probability
    /// is the sum over every offset they absorb** — never the probability of sitting exactly on
    /// the edge (`spec/parameter_prepass_ssr.md` §4.1). Walking the slip distances and letting
    /// `bucket_of` clamp is what makes that true by construction rather than by a special case at
    /// the two ends.
    ///
    /// **Why it is not the tempting shortcut, measured twice over.** Scoring an end bucket as
    /// though the read sat on the edge fails the sums-to-one gate outright — the buckets come to
    /// 0.9488 at a recorded range of ±1 — and on a stratum whose alleles reach three copies either
    /// side it returns the slippage level **52% low**. Rescaling that plug-in so the buckets do
    /// sum to one repairs the gate and not the bias: it then runs **33% high** where 30 in 100
    /// slipped reads take a second step, which is the regime long tracts sit in. Neither shows
    /// from outside. The marginal rule is exactly unbiased at every recorded range tried
    /// (research note §6.4).
    ///
    /// **This is also what makes the recorded range the cheap width and the allele support the
    /// expensive one.** Against the reference origin an end bucket absorbs whole *alleles* and not
    /// only far slips, and the marginal rule is what lets the fit attribute that mass to a distant
    /// allele instead of to a distant slip — so the range can be narrow and the support cannot.
    /// **The kernel arrives already evaluated** ([`slip_kernel`]), because a cell's score calls
    /// this once per allele length and the kernel does not depend on the allele.
    fn read_bucket_probabilities(
        &self,
        kernel: &[f64; SLIP_DISTANCES],
        allele: WholeRepeatOffset,
    ) -> [f64; OFFSET_BUCKETS] {
        let mut probabilities = [0.0; OFFSET_BUCKETS];
        for slip in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
            // Saturating rather than checked, and it loses nothing: `bucket_of` clamps to
            // ±OFFSET_HALF_RANGE anyway, so an offset past an `i8`'s ends is already inside the
            // end bucket that absorbs it. (It cannot happen at today's widths — an allele reaches
            // ±6 and a slip ±8 — but the arithmetic says so rather than a comment.)
            let landed = i32::from(allele.get()) + slip;
            let shown = WholeRepeatOffset(i8::try_from(landed).unwrap_or(if landed < 0 {
                i8::MIN
            } else {
                i8::MAX
            }));
            let at = usize::try_from(slip + MAX_SLIP_STEP).expect("a distance inside the range");
            probabilities[bucket_of(shown).index()] += kernel[at];
        }
        probabilities
    }
}

/// What the fit needs of a cell besides how to score it: how many genome copies its loci sit on,
/// and how many loci it stands for.
///
/// The impl sits here rather than beside [`StratumCell`] because it is this path's plug into
/// `fitting/`, and the other one — [`NoiseModel`] — is below it.
impl WeightedCell for StratumCell {
    fn ploidy(&self) -> Ploidy {
        StratumCell::ploidy(*self)
    }

    /// **Loci, not reads.** An entry is one locus's whole shape, and what the fit weighs it by is
    /// how many loci showed that shape — which is the design's central choice rather than a
    /// bookkeeping one (`spec/parameter_prepass_ssr.md` §4.1).
    fn sites(&self) -> u64 {
        StratumCell::loci(*self)
    }
}

impl NoiseModel for SsrNoiseModel {
    type Cell = StratumCell;
    type NoiseParams = SlippageModel;

    /// **Unordered tuples of allele lengths** — `C(A + P − 1, P)` of them for `A` lengths and `P`
    /// genome copies, so 91 for a diploid stratum with the full 13 lengths and 66 for one at four
    /// reference copies, whose support clips at −4.
    ///
    /// **This is the number the sibling path's `ploidy + 1` would have got wrong**, which is why
    /// the trait asks the model rather than deriving it: there a genotype is how many copies carry
    /// the alternative allele, and a diploid has three of them however many alleles are in play.
    ///
    /// # Panics
    ///
    /// If the count passes a `usize`, which needs a ploidy in the hundreds against the full
    /// support — `Ploidy` admits one, and a genotype set that large is a run nobody wants
    /// silently wrapped to a small number.
    fn genotypes(&self, ploidy: Ploidy) -> usize {
        genotype_count(self.allele_support.len(), usize::from(ploidy.get()))
    }

    /// `ln L(shape | genotype)` for every genotype of this stratum, appended in the model's own
    /// order: allele-index tuples in non-decreasing order, ascending lexicographically, so a
    /// diploid stratum over lengths `a, b, c` gives `aa, ab, ac, bb, bc, cc`.
    ///
    /// **One multinomial over the buckets, at the bucket probabilities the genotype implies.**
    /// Each read picks one of the locus's `P` copies with equal chance and then slips, so a
    /// bucket's probability is the average of the copies' own bucket distributions
    /// ([`SsrNoiseModel::read_bucket_probabilities`], which is where the end-bucket rule lives),
    /// and the shape's probability is
    ///
    /// ```text
    ///                        n!
    /// L(shape | g)  =  ──────────────  ·  Π  q_b(g)^{n_b}
    ///                    Π_b  n_b!         b
    /// ```
    ///
    /// **The guard reads are not in `n`.** A read that showed a length differing from the
    /// reference by something other than a whole number of copies is modelled as an independent
    /// per-read outcome, so the likelihood splits exactly into *how many reads did that* times
    /// *how the rest fell across the buckets* — nothing about slippage is estimated from the guard
    /// and nothing about the guard disturbs it (`spec/parameter_prepass_ssr.md` §4.1). So `n` here
    /// is the shape's whole-repeat depth, and a locus whose every read landed in the guard is an
    /// empty product — a likelihood of **one**, `ln L = 0`, for every genotype. That is a locus
    /// with nothing to say about which alleles it carries, not a locus that refutes them all.
    ///
    /// **The multinomial coefficient is computed rather than dropped**, though it is the same for
    /// every genotype and so cancels out of the mixture and cannot move a fit. What it buys is
    /// that the scoring rule **sums to one over the shapes at a given whole-repeat depth** — an
    /// identity, and
    /// the first of the three algebraic gates the design puts before any fitting (spec §10). That
    /// identity is what separates this rule from the plug-in it replaced, which fails it at 0.9488
    /// (`arch/parameter_prepass_ssr.md` §3), so it is worth an `ln Γ` per cell to keep it
    /// checkable. The sibling path keeps its own coefficient for the same reason.
    ///
    /// `−∞` is a legal entry and says this genotype cannot have produced this shape: at a slippage
    /// level of zero, a locus with a read anywhere but on its own alleles is exactly that.
    ///
    /// # Panics
    ///
    /// If `ploidy` disagrees with the cell's own. The two arrive separately through the trait, and
    /// scoring a shape against a genotype set built for a different number of genome copies is a
    /// wrong answer that nothing downstream could question.
    fn append_genotype_likelihoods(
        &self,
        cell: &StratumCell,
        noise: &SlippageModel,
        ploidy: Ploidy,
        out: &mut Vec<f64>,
    ) {
        assert_eq!(
            ploidy,
            cell.ploidy(),
            "a cell of ploidy {} scored against the genotypes of ploidy {ploidy}",
            cell.ploidy()
        );

        let shape = cell.shape();
        let reads_by_bucket = shape.reads_by_bucket();
        // Common to every genotype, so it cannot move a fit — kept because it is what makes the
        // sums-to-one gate an identity. See the doc comment.
        let ln_arrangements = ln_factorial(shape.whole_repeat_depth())
            - reads_by_bucket
                .iter()
                .map(|&reads| ln_factorial(u32::from(reads)))
                .sum::<f64>();

        // One row per allele length, built once per cell rather than once per genotype: a diploid
        // stratum has 13 lengths and 91 genotypes, so this is 13 walks of the slip kernel instead
        // of 182. On the stack, because the support cannot exceed `MAX_ALLELE_LENGTHS`.
        let kernel = slip_kernel(noise);
        let mut bucket_probabilities_by_allele = [[0.0; OFFSET_BUCKETS]; MAX_ALLELE_LENGTHS];
        for (row, &allele) in bucket_probabilities_by_allele
            .iter_mut()
            .zip(&self.allele_support)
        {
            *row = self.read_bucket_probabilities(&kernel, allele);
        }

        let copies = f64::from(ploidy.get());
        out.reserve(self.genotypes(ploidy));
        self.for_each_genotype(ploidy, |genotype| {
            let mut ln_likelihood = ln_arrangements;
            for (bucket, &reads) in reads_by_bucket.iter().enumerate() {
                if reads == 0 {
                    // `0 · ln 0` is zero by the multinomial's own convention: an empty bucket
                    // against an impossible outcome contributes nothing, where the arithmetic
                    // would give `NaN`.
                    continue;
                }
                let share: f64 = genotype
                    .iter()
                    .map(|&allele| bucket_probabilities_by_allele[allele][bucket])
                    .sum::<f64>()
                    / copies;
                ln_likelihood += f64::from(reads) * share.ln();
            }
            out.push(ln_likelihood);
        });
    }
}

impl SsrNoiseModel {
    /// Every genotype of this stratum in the model's own order, as indices into
    /// [`allele_support`](Self::allele_support): non-decreasing tuples, ascending
    /// lexicographically.
    ///
    /// **Non-decreasing is what makes a genotype unordered**: a locus carrying the reference
    /// length and one copy short is the same genotype either way round, so exactly one spelling of
    /// it exists. Emitting both would double-count it and leave the fitted frequencies summing to
    /// more than one.
    ///
    /// A callback rather than an iterator, because a genotype is a *slice* whose length is the
    /// ploidy: an iterator would have to hand out an allocation per genotype, and this is called
    /// once per cell per candidate.
    ///
    /// **Public because the order is part of the answer, not an implementation detail.** A fit
    /// returns one frequency per genotype in exactly this sequence, so whatever reports those
    /// frequencies has to walk the same one — and re-deriving the walk in another module is how
    /// two orders that disagree get written.
    pub fn for_each_genotype(&self, ploidy: Ploidy, mut visit: impl FnMut(&[usize])) {
        let alleles = self.allele_support.len();
        let copies = usize::from(ploidy.get());
        // `Ploidy` rejects zero and the support is never empty, so neither guard fires today.
        // What they prevent is a width disagreement rather than a crash: measured, without them
        // the loop emits one bogus genotype and returns, while `genotypes()` counts zero — and the
        // scan sizes its row-major table from `genotypes()`.
        if alleles == 0 || copies == 0 {
            return;
        }

        let mut genotype: SmallVec<[usize; 4]> = SmallVec::from_elem(0, copies);
        loop {
            visit(&genotype);
            // Advance to the next non-decreasing tuple: raise the rightmost position that can
            // still rise, and reset everything to its right to the same value, which is what
            // keeps the tuple non-decreasing.
            let mut at = copies;
            loop {
                if at == 0 {
                    return;
                }
                at -= 1;
                if genotype[at] + 1 < alleles {
                    let raised = genotype[at] + 1;
                    genotype[at..].fill(raised);
                    break;
                }
            }
        }
    }
}

/// `C(alleles + copies − 1, copies)` — how many unordered tuples of `copies` alleles can be drawn
/// from `alleles` lengths, repetition allowed.
///
/// Built one factor at a time so that every partial product is itself a binomial coefficient and
/// so exact in integers: `C(a + i − 1, i)` after step `i`. Computed in `u128` because a caller
/// asking for a high ploidy over a wide support is asking for a big number, and a silently wrapped
/// genotype count would size the fit's table wrongly rather than fail.
fn genotype_count(alleles: usize, copies: usize) -> usize {
    let mut count: u128 = 1;
    for step in 1..=copies {
        count = count * (alleles + step - 1) as u128 / step as u128;
    }
    usize::try_from(count).unwrap_or_else(|_| {
        panic!("{alleles} allele lengths at ploidy {copies} make {count} genotypes")
    })
}

/// `ln n!`, through `ln Γ(n + 1)` — the same helper and the same route the sibling path's noise
/// model takes.
fn ln_factorial(n: u32) -> f64 {
    lgamma(f64::from(n) + 1.0)
}

/// The logit and its inverse — the scale the two share-shaped parameters are searched on, so a
/// search cannot walk out of `(0, 1)` whatever step it takes.
fn logit(share: f64) -> f64 {
    (share / (1.0 - share)).ln()
}

fn expit(scaled: f64) -> f64 {
    1.0 / (1.0 + (-scaled).exp())
}

/// **How wide the search may range on each of the three axes**, on that axis's own scale.
///
/// The level runs from 1e-5 to 0.6 — the measured range is 0.00091 to 0.150, so both ends have
/// room and neither is where a real answer sits. The two shares stop just inside 0 and 1, because
/// the logit of either endpoint is infinite and the golden section needs a finite bracket; a fit
/// that genuinely wants a share of zero rails at 0.005 or 0.002 and is visible as having done so.
const LEVEL_RANGE: (f64, f64) = (1e-5, 0.6);
const GAIN_SHARE_RANGE: (f64, f64) = (0.005, 0.995);
const STEP_DECAY_RANGE: (f64, f64) = (0.002, 0.95);

/// The lowest level [`slippage_starts`] will build a start around — above the search's own lower
/// bound, so that the four multipliers stay four distinct points. See that function.
const STARTING_LEVEL_FLOOR: f64 = 1e-4;

/// **The three slippage parameters as a search walks them**, each on the scale that makes a step
/// mean the same thing everywhere along it.
///
/// The level is searched on a **log** scale, because it spans twenty-two-fold across the repeat
/// counts of one dataset and a step that is 1% of 0.15 is 165% of 0.00091. The two shares are
/// searched on a **logit** scale, because they are fractions in `(0, 1)` and a search on that
/// scale cannot walk out of the range whatever step it takes.
///
/// **The level is axis 0**, which is what makes it the parameter
/// [`MultistartResult::headline_spread`](crate::ng::parameter_estimation::fitting::multistart::MultistartResult::headline_spread)
/// reports the disagreement in — and the level is the number a stratum is chosen by, so it is the
/// right one to judge a fit on.
impl SearchableNoise for SlippageModel {
    fn axes() -> usize {
        3
    }

    fn axis_bounds(axis: usize) -> (f64, f64) {
        match axis {
            0 => (LEVEL_RANGE.0.ln(), LEVEL_RANGE.1.ln()),
            1 => (logit(GAIN_SHARE_RANGE.0), logit(GAIN_SHARE_RANGE.1)),
            2 => (logit(STEP_DECAY_RANGE.0), logit(STEP_DECAY_RANGE.1)),
            other => panic!("a slippage model has three axes and was asked for axis {other}"),
        }
    }

    fn coordinate(&self, axis: usize) -> f64 {
        match axis {
            0 => self.slip_rate.get().ln(),
            1 => logit(self.gain_share.get()),
            2 => logit(self.step_decay.get()),
            other => panic!("a slippage model has three axes and was asked for axis {other}"),
        }
    }

    fn headline(&self) -> f64 {
        self.slip_rate.get()
    }

    fn with_coordinate(&self, axis: usize, value: f64) -> Self {
        // Clamped back into the range each constrained type accepts, because the exponential and
        // the logistic can each return a value one ulp outside `[0, 1]` at the ends of their
        // brackets — and `try_new` would then reject a point the search legitimately reached.
        let moved = match axis {
            0 => Self::try_new(
                value.exp().clamp(0.0, 1.0),
                self.gain_share.get(),
                self.step_decay.get(),
            ),
            1 => Self::try_new(
                self.slip_rate.get(),
                expit(value).clamp(0.0, 1.0),
                self.step_decay.get(),
            ),
            2 => Self::try_new(
                self.slip_rate.get(),
                self.gain_share.get(),
                expit(value).clamp(0.0, 1.0),
            ),
            other => panic!("a slippage model has three axes and was asked for axis {other}"),
        };
        moved.expect("every axis is written back inside its own type's range")
    }
}

/// **Four starting points, each disagreeing about the level, the direction and the fall-off at
/// once** (`spec/parameter_prepass_ssr.md` §4.2).
///
/// **Starts that disagree only about the headline parameter are not a spread**, and that is the
/// trap this set exists to avoid: on the sibling path's inbreeding fit, five starts that disagreed
/// about the headline number while sharing one guess at a nuisance axis returned a confident zero
/// on a genome 29% covered by runs.
///
/// **The level is a multiplier on a starting estimate rather than an absolute value**, because a
/// stratum's level spans twenty-two-fold across repeat counts, so a fixed ladder of absolute rates
/// would begin every stratum in the wrong place. `level_estimate` is meant to be the share of
/// reads sitting off the reference length — an over-estimate, since it counts real alleles too,
/// which is why the four multipliers run below one as well as above it.
///
/// **The floor is 1e-4 and not the search's own lower bound of 1e-5**, which is the reference
/// implementation's value and is what keeps the four starts distinct: `1e-4 × {3, ⅓, 1, 0.3}` is
/// `3e-4, 3.3e-5, 1e-4, 3e-5`, all four clear of the bound. Floored at 1e-5 instead, a stratum
/// where no read left the reference length gives an estimate of zero and **three of the four
/// starts collapse onto the same point** — the set would still be four starts and would no longer
/// be a spread.
///
/// # Panics
///
/// If `level_estimate` is not a number. It is a share of reads, so a caller reaches this by
/// dividing by a zero read count — a stratum with no whole-repeat reads at all — and a `NaN`
/// clamped and handed to the search would make every candidate score alike.
#[must_use]
pub fn slippage_starts(level_estimate: f64) -> SmallVec<[SlippageModel; 4]> {
    assert!(
        !level_estimate.is_nan(),
        "the starting points were asked for at a slippage estimate that is not a number, which \
         is what dividing by a stratum with no whole-repeat reads gives"
    );
    let level = level_estimate.clamp(STARTING_LEVEL_FLOOR, 0.3);
    // (multiplier on the estimate, gain share, step decay)
    [
        (3.0, 0.20, 0.03),
        (1.0 / 3.0, 0.80, 0.40),
        (1.0, 0.50, 0.15),
        (0.3, 0.35, 0.08),
    ]
    .into_iter()
    .map(|(multiplier, gain_share, step_decay)| {
        let scaled = (level * multiplier).clamp(LEVEL_RANGE.0, 0.5);
        SlippageModel::try_new(scaled, gain_share, step_decay)
            .expect("three probabilities, each built inside its own range")
    })
    .collect()
}

impl fmt::Display for SlippageModel {
    /// `0.02010 of reads slipping, 0.170 of those gaining, 0.065 of those taking a further
    /// step` — the shape a summary line over several hundred strata wants.
    ///
    /// **Each number names its own denominator**, because the three are over different
    /// populations and the differences are large: the gain share is over the reads that
    /// slipped rather than over all reads, which at a level of 0.02 is a fiftyfold difference
    /// in what "0.17" would mean.
    ///
    /// **Five decimals on the level and not four**, because four cannot tell a stratum that
    /// barely slips from one that does not slip at all: the bottom of the measured range,
    /// 0.00091, renders as `0.0009` at four and anything under 0.00005 renders as `0.0000` —
    /// the same text as a genuine zero, which [`SlipRate`] documents as a real answer.
    ///
    /// Destructured, so a fourth parameter added to the model is a compile error here rather
    /// than a line that silently describes three quarters of it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self {
            slip_rate,
            gain_share,
            step_decay,
        } = self;
        write!(
            f,
            "{:.5} of reads slipping, {:.3} of those gaining, {:.3} of those taking a further step",
            slip_rate.get(),
            gain_share.get(),
            step_decay.get()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::parameter_estimation::fitting::multistart::{
        SearchPrecision, fit_by_multistart,
    };

    use super::super::OFFSET_HALF_RANGE;
    use super::super::stratum_table::{LocusShape, StratumEntry};

    /// The slippage models the kernel tests sweep, chosen so that each of them **breaks
    /// something** rather than to look like real data. Three groups of fall-offs, and each group
    /// is the only one that reaches its own failure.
    ///
    /// **0.065 to 0.123 — what real data shows.** They already catch a missing renormaliser, but
    /// only just: the mass the truncation moves is `level · f⁸`, which at a level of 0.00091 and a
    /// fall-off of 0.12 leaves the total at 0.9999999999609. That is 39 times the tolerance below,
    /// so it fails — while reading like rounding to whoever has to diagnose it.
    ///
    /// **0.9 and 0.95 — the same failure made unmissable.** At 0.9 the truncation moves 43% of the
    /// slipped mass and at 0.95 it moves 66%, so a kernel that forgot to hand that mass back sums
    /// to 0.57 and 0.34 at a level of one rather than to something in the twelfth decimal.
    ///
    /// **`1 − 1e-8` and `1 − 1e-9` — the band where the closed form loses the identity.** These
    /// reject the reference implementation's `(1 − f)/(1 − f^S)`, which subtracts two nearly equal
    /// numbers there: at a level of one it sums to 1.000000001082 and 0.999999996500, which is
    /// 1,100 and 3,500 times the tolerance below. Nothing else in the sweep sees it — the rows
    /// either side, 0.95 and exactly 1.0, are both fine under the closed form — and the margin is
    /// not uniform across the sweep, because the deviation scales with the level: the thinnest
    /// failing row, at a level of 0.00091, misses by only 3.2 times the tolerance.
    ///
    /// The endpoints 0.0 and 1.0 on every parameter are here because [`SlipRate`],
    /// [`SlipGainShare`] and [`SlipStepDecay`] all accept them, and a fall-off of exactly one is
    /// where the closed form divides zero by zero.
    fn kernels_that_break_things() -> Vec<SlippageModel> {
        let mut models = Vec::new();
        for &rate in &[0.0, 0.00091, 0.0201, 0.15, 1.0] {
            for &gain in &[0.0, 0.17, 0.5, 0.83, 1.0] {
                for &decay in &[
                    0.0,
                    0.065,
                    0.123,
                    0.4,
                    0.9,
                    0.95,
                    1.0 - 1e-8,
                    1.0 - 1e-9,
                    1.0,
                ] {
                    models.push(
                        SlippageModel::try_new(rate, gain, decay).expect("three probabilities"),
                    );
                }
            }
        }
        models
    }

    /// **The first of the three algebraic gates the design puts before any fitting** (§10.1): a
    /// rule whose probabilities do not sum to one is not the likelihood of anything, and no
    /// consistency result covers it. Here it is asked of the kernel alone, over every distance a
    /// read may slip.
    ///
    /// The sweep is deliberately hostile — see [`kernels_that_break_things`], which says what
    /// each group of fall-offs is the only one to catch. The tolerance is 1e-12 against a worst
    /// residual of about 1e-15 on correct code, so it has three orders of headroom and is not
    /// flaky.
    #[test]
    fn the_slip_kernel_is_a_distribution_at_every_parameter_setting() {
        for model in kernels_that_break_things() {
            let total: f64 = (-MAX_SLIP_STEP..=MAX_SLIP_STEP)
                .map(|copies| model.probability_of_slipping_by(copies))
                .sum();

            assert!(
                (total - 1.0).abs() < 1e-12,
                "{model} sums to {total} over its support"
            );
            for copies in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
                let probability = model.probability_of_slipping_by(copies);
                assert!(
                    (0.0..=1.0).contains(&probability),
                    "{model} charges {probability} to a slip of {copies} copies"
                );
            }
        }
    }

    /// **With nothing slipping, every read shows its own allele's length.** The third algebraic
    /// gate, at the kernel: a kernel putting mass anywhere but zero at a level of zero is
    /// describing movement that cannot happen, whatever the other two parameters say.
    #[test]
    fn a_stratum_where_nothing_slips_puts_every_read_on_its_own_allele() {
        for gain in [0.0, 0.17, 1.0] {
            for decay in [0.0, 0.065, 1.0] {
                let silent = SlippageModel::try_new(0.0, gain, decay).expect("three probabilities");

                assert_eq!(silent.probability_of_slipping_by(0), 1.0, "{silent}");
                for copies in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
                    if copies != 0 {
                        assert_eq!(
                            silent.probability_of_slipping_by(copies),
                            0.0,
                            "{silent} moves a read by {copies} copies"
                        );
                    }
                }
            }
        }
    }

    /// **Each further copy costs the fall-off, and it is the same number in both directions** —
    /// which is §3's decision, taken because the gaining arm rests on 3 to 13 reads above
    /// dinucleotides and a free parameter there would fit counting noise.
    ///
    /// **Every rung, not just the second copy against the first.** The sums-to-one gate is
    /// invariant under any *permutation* of the eight distances, and the direction test only
    /// compares the two one-copy arms, so a kernel that exchanged the masses of seven and eight
    /// copies would pass both — and would put a read that slipped eight copies at `1/f` times its
    /// true probability, 2.5-fold at a fall-off of 0.4. Only walking the ladder sees it.
    #[test]
    fn each_further_copy_costs_the_fall_off_whichever_way_the_read_slipped() {
        for decay in [0.065, 0.123, 0.4] {
            let model = SlippageModel::try_new(0.0201, 0.17, decay).expect("three probabilities");

            for copies in 2..=MAX_SLIP_STEP {
                let losing = model.probability_of_slipping_by(-copies)
                    / model.probability_of_slipping_by(-(copies - 1));
                let gaining = model.probability_of_slipping_by(copies)
                    / model.probability_of_slipping_by(copies - 1);

                assert!(
                    (losing - decay).abs() < 1e-12,
                    "losing at {copies} copies: {losing} vs {decay}"
                );
                assert!(
                    (gaining - decay).abs() < 1e-12,
                    "gaining at {copies} copies: {gaining} vs {decay}"
                );
            }
        }
    }

    /// **A read is far likelier to have lost copies than gained them, and the kernel must have
    /// that the way round the data does.**
    ///
    /// This is the assertion the three gates above cannot make. Swap the two arms — charge a
    /// gain to `1 − gain_share` — and the kernel still sums to one, still falls silent at a
    /// level of zero, and still decays by the fall-off in both directions: all three pass, and
    /// the fitted direction split comes back inverted. **Inversion is precisely the failure the
    /// whole step exists to remove** — production's estimator reports gains as marginally more
    /// common than losses (`spec/parameter_prepass.md` §2.2) — so it is worth its own test.
    ///
    /// At the measured tomato dinucleotide split of 0.17, losses lead gains 4.9-fold.
    #[test]
    fn a_read_loses_copies_far_more_often_than_it_gains_them() {
        let model = SlippageModel::try_new(0.0201, 0.17, 0.065).expect("three probabilities");

        let gaining = model.probability_of_slipping_by(1);
        let losing = model.probability_of_slipping_by(-1);

        assert!(
            losing > gaining,
            "a read at a tomato dinucleotide loses copies more often than it gains them: \
             {losing} against {gaining}"
        );
        let imbalance = losing / gaining;
        assert!(
            (imbalance - 0.83 / 0.17).abs() < 1e-12,
            "losses lead gains {imbalance}-fold at a gain share of 0.17"
        );
        // Of all reads, and not only of the slipped ones — so a kernel that lost the level
        // somewhere in the direction split would show here. The two arms share the level and
        // the one-copy term of the renormalised geometric, `1/(1 + f + … + f⁷)`, and split only
        // on the direction, so together they come to the whole of it.
        let one_copy = 1.0
            / (0..MAX_SLIP_STEP)
                .map(|step| 0.065_f64.powi(step))
                .sum::<f64>();
        assert!(
            (gaining + losing - 0.0201 * one_copy).abs() < 1e-15,
            "{gaining} + {losing} against {}",
            0.0201 * one_copy
        );
    }

    /// **Beyond the truncation a read shows nothing**, and the boundary is the copy the
    /// truncation keeps rather than the one after it — an off-by-one either way is a kernel
    /// that either drops a distance it renormalised for or charges one it did not.
    ///
    /// **`i32::MIN` is asserted rather than `i32::MIN + 1`**, and the difference is the whole
    /// point: `i32::MIN.abs()` overflows, so a kernel that takes the absolute value before
    /// checking the range panics in a debug build and, in a release one, wraps to a negative
    /// count that walks straight past the range guard. Measured on that version, a model of
    /// `(0.15, 0.5, 1.0)` charged 0.009375 to a slip of −2,147,483,648 copies — exactly what it
    /// charges to a slip of one.
    ///
    /// **The constant's value is pinned here too**, because every other assertion in the module
    /// spells it symbolically: without this line, setting [`MAX_SLIP_STEP`] to 2 leaves all
    /// fourteen tests green while the kernel returns exactly zero for a read that slipped three
    /// copies, and setting it to 4 leaves it green while silently making it the same width as
    /// [`OFFSET_HALF_RANGE`](super::OFFSET_HALF_RANGE) — the confusion the constant's own doc
    /// comment warns about in prose.
    #[test]
    fn no_read_slips_further_than_the_truncation_and_the_last_copy_inside_it_is_kept() {
        assert_eq!(MAX_SLIP_STEP, 8);

        let model = SlippageModel::try_new(0.15, 0.5, 0.4).expect("three probabilities");

        assert!(model.probability_of_slipping_by(MAX_SLIP_STEP) > 0.0);
        assert!(model.probability_of_slipping_by(-MAX_SLIP_STEP) > 0.0);
        assert_eq!(model.probability_of_slipping_by(MAX_SLIP_STEP + 1), 0.0);
        assert_eq!(model.probability_of_slipping_by(-MAX_SLIP_STEP - 1), 0.0);
        assert_eq!(model.probability_of_slipping_by(i32::MAX), 0.0);
        assert_eq!(model.probability_of_slipping_by(i32::MIN), 0.0);
    }

    /// **A fall-off of one is a legal parameter and must not produce a `NaN`.** The closed form
    /// the reference implementation uses, `(1 − f)·f^{s−1} / (1 − f^S)`, is `0/0` at `f = 1`; the
    /// limit is a slip distance drawn uniformly from one to [`MAX_SLIP_STEP`] copies, and the
    /// partial-sum renormaliser returns it without a special case.
    ///
    /// A `NaN` here would not fail loudly. It reaches the likelihood, and the searches in this
    /// crate pick their maximum with `total_cmp`, which ranks `NaN` above every finite score — so
    /// a fall-off of one would be *selected* rather than skipped, and reported with a `NaN`
    /// likelihood beside it.
    #[test]
    fn a_fall_off_of_one_spreads_the_distance_evenly_instead_of_returning_a_nan() {
        let model = SlippageModel::try_new(0.2, 0.17, 1.0).expect("three probabilities");

        for copies in 1..=MAX_SLIP_STEP {
            let gaining = model.probability_of_slipping_by(copies);
            let losing = model.probability_of_slipping_by(-copies);
            assert!(gaining.is_finite() && losing.is_finite(), "{copies}");
            assert!(
                (gaining - 0.2 * 0.17 / f64::from(MAX_SLIP_STEP)).abs() < 1e-12,
                "at {copies} copies: {gaining}"
            );
            assert!((losing - 0.2 * 0.83 / f64::from(MAX_SLIP_STEP)).abs() < 1e-12);
        }
    }

    /// Both endpoints are real answers rather than degenerate ones, and each has a meaning a
    /// fit can reach: a level of exactly zero is a stratum where nothing slipped, a gain share
    /// of zero is one where every slipped read lost repeats, and a decay of zero is one where
    /// no read took a second step.
    #[test]
    fn every_slippage_rate_accepts_both_endpoints_and_round_trips() {
        for value in [0.0, 1.0, 0.0201, 0.17, 0.065] {
            assert_eq!(SlipRate::try_new(value).unwrap().get(), value);
            assert_eq!(SlipGainShare::try_new(value).unwrap().get(), value);
            assert_eq!(SlipStepDecay::try_new(value).unwrap().get(), value);
        }
    }

    /// A value outside `[0, 1]` is rejected **as the quantity it was offered for**, so a log
    /// line names which of the three parameters of one fit was wrong rather than saying that
    /// some fraction was.
    ///
    /// **Both bounds on all three**, which is the standard `types.rs` sets for its own
    /// constrained rates and states the reason for: a test that only ever crosses one bound
    /// leaves the other free to be widened. Here the three types are structurally identical,
    /// so a widening in one and not its siblings is exactly the drift that goes unseen.
    #[test]
    fn each_slippage_rate_rejects_both_bounds_under_its_own_name() {
        for below in [-0.01, -1.0] {
            assert!(matches!(
                SlipRate::try_new(below),
                Err(DomainError::SlipRate(_))
            ));
            assert!(matches!(
                SlipGainShare::try_new(below),
                Err(DomainError::SlipGainShare(_))
            ));
            assert!(matches!(
                SlipStepDecay::try_new(below),
                Err(DomainError::SlipStepDecay(_))
            ));
        }
        for above in [1.01, 2.0] {
            assert!(matches!(
                SlipRate::try_new(above),
                Err(DomainError::SlipRate(_))
            ));
            assert!(matches!(
                SlipGainShare::try_new(above),
                Err(DomainError::SlipGainShare(_))
            ));
            assert!(matches!(
                SlipStepDecay::try_new(above),
                Err(DomainError::SlipStepDecay(_))
            ));
        }

        let messages = [
            SlipRate::try_new(-0.01).unwrap_err().to_string(),
            SlipGainShare::try_new(1.01).unwrap_err().to_string(),
            SlipStepDecay::try_new(2.0).unwrap_err().to_string(),
        ];
        assert!(messages[0].contains("slippage rate"), "{}", messages[0]);
        assert!(
            messages[1].contains("gain share of slipped reads"),
            "the message says what the share is of: {}",
            messages[1]
        );
        assert!(messages[2].contains("step decay"), "{}", messages[2]);
    }

    /// **The three non-values a search produces when it goes wrong**, and none of them may
    /// become a parameter: a division by zero gives an infinity, `0.0 / 0.0` gives `NaN`, and
    /// a `NaN` that reaches a likelihood makes every candidate score alike, which reads from
    /// the outside as a search that found a flat surface.
    #[test]
    fn no_slippage_rate_admits_a_nan_or_an_infinity() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(SlipRate::try_new(bad).is_err(), "{bad}");
            assert!(SlipGainShare::try_new(bad).is_err(), "{bad}");
            assert!(SlipStepDecay::try_new(bad).is_err(), "{bad}");
        }
    }

    /// The model keeps the three numbers in the roles they were given. A transposition here
    /// would be invisible to the type system — all three are fractions — and would put the
    /// direction share on the axis the level belongs to, which is the failure the three
    /// separate types exist to make impossible at every *other* call site.
    #[test]
    fn a_slippage_model_carries_its_three_numbers_in_their_own_roles() {
        let model = SlippageModel::try_new(0.0201, 0.17, 0.065).expect("three probabilities");

        assert_eq!(model.slip_rate.get(), 0.0201);
        assert_eq!(model.gain_share.get(), 0.17);
        assert_eq!(model.step_decay.get(), 0.065);
    }

    /// The three-fraction door names the bad value as the parameter it was offered for, so a
    /// caller building a starting point from a table of numbers learns which column was
    /// wrong.
    ///
    /// **All three columns, and the first one is the reason this test says so.** With only
    /// the second and third checked, replacing the first check with an unchecked construction
    /// left every test green while `try_new(NaN, 0.17, 0.065)` returned `Ok` carrying a `NaN`
    /// level — a parameter that makes every candidate in a search score alike.
    #[test]
    fn a_slippage_model_refuses_whichever_parameter_is_not_a_probability() {
        assert!(matches!(
            SlippageModel::try_new(f64::NAN, 0.17, 0.065),
            Err(DomainError::SlipRate(_))
        ));
        assert!(matches!(
            SlippageModel::try_new(1.5, 0.17, 0.065),
            Err(DomainError::SlipRate(_))
        ));
        assert!(matches!(
            SlippageModel::try_new(0.02, 1.5, 0.065),
            Err(DomainError::SlipGainShare(_))
        ));
        assert!(matches!(
            SlippageModel::try_new(0.02, 0.17, -0.1),
            Err(DomainError::SlipStepDecay(_))
        ));
    }

    /// **When more than one column is wrong, the leftmost is the one reported** — which is
    /// what the `# Errors` clause promises and what nothing else here would notice, since a
    /// fixture with a single bad column gives the same answer under every ordering of the
    /// three checks.
    #[test]
    fn a_slippage_model_reports_the_leftmost_bad_parameter() {
        assert!(matches!(
            SlippageModel::try_new(2.0, 1.5, -0.1),
            Err(DomainError::SlipRate(_))
        ));
        assert!(matches!(
            SlippageModel::try_new(0.02, 1.5, -0.1),
            Err(DomainError::SlipGainShare(_))
        ));
    }

    /// A fitted model renders in the words that make its three numbers readable — a summary
    /// line over several hundred strata is the only place most of them are ever seen — and
    /// each number names the population it is a share of.
    #[test]
    fn a_slippage_model_renders_each_number_with_what_it_measures() {
        let rendered = SlippageModel::try_new(0.0201, 0.17, 0.065)
            .unwrap()
            .to_string();

        assert_eq!(
            rendered,
            "0.02010 of reads slipping, 0.170 of those gaining, 0.065 of those taking a further step"
        );
    }

    /// **A stratum that barely slips must not read as one that does not slip at all.** The
    /// bottom of the measured range — 0.00091 below four repeats — has to survive the
    /// rendering, and at four decimals it does not: it would print `0.0009`, and a level an
    /// order of magnitude smaller would print `0.0000`, which is the text a genuine zero
    /// gets. Half the loci this path sees sit in strata at that level.
    #[test]
    fn a_barely_slipping_stratum_renders_differently_from_one_that_never_slips() {
        let barely = SlippageModel::try_new(0.00091, 0.17, 0.065).unwrap();
        let fainter = SlippageModel::try_new(0.00003, 0.17, 0.065).unwrap();
        let never = SlippageModel::try_new(0.0, 0.17, 0.065).unwrap();

        assert!(barely.to_string().starts_with("0.00091"), "{barely}");
        assert_ne!(barely.to_string(), never.to_string());
        assert_ne!(
            fainter.to_string(),
            never.to_string(),
            "a level below the measured range still has to be distinguishable from zero"
        );
    }

    // -----------------------------------------------------------------
    // D2 — the noise model, and the three algebraic checks the design
    // puts before anything is fitted (`spec/parameter_prepass_ssr.md` §10).
    // -----------------------------------------------------------------

    /// The slippage the measurements are taken at: 2.0% of reads slipping, 17 in 100 of those
    /// gaining, 9 in 100 of those taking a second step — `spec/parameter_prepass_ssr.md` §3's
    /// tomato dinucleotide row, and the truth every table in §4.1 is generated from.
    fn measured_slippage() -> SlippageModel {
        SlippageModel::try_new(0.0201, 0.17, 0.09).expect("three probabilities")
    }

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("at least one genome copy")
    }

    /// Every way `depth` reads can fall across the nine offset buckets — the **entry space** at
    /// that depth, which is what "sums to one" is a sum over. 165 shapes at three reads, 495 at
    /// four.
    ///
    /// No guard reads, because the guard is not in the length half of the likelihood at all; the
    /// test that says so is [`the_guard_reads_do_not_enter_the_length_likelihood`].
    fn shapes_at_whole_repeat_depth(depth: u32) -> Vec<LocusShape> {
        fn walk(
            remaining: u32,
            bucket: usize,
            counts: &mut [u32; OFFSET_BUCKETS],
            out: &mut Vec<LocusShape>,
        ) {
            if bucket == OFFSET_BUCKETS - 1 {
                counts[bucket] = remaining;
                out.push(LocusShape::try_new(*counts, 0).expect("a shape at the depth asked for"));
                return;
            }
            for here in 0..=remaining {
                counts[bucket] = here;
                walk(remaining - here, bucket + 1, counts, out);
            }
            counts[bucket] = 0;
        }

        let mut out = Vec::new();
        let mut counts = [0; OFFSET_BUCKETS];
        walk(depth, 0, &mut counts, &mut out);
        out
    }

    fn cell_of(shape: LocusShape, ploidy: Ploidy) -> StratumCell {
        StratumCell::new(StratumEntry { shape, loci: 1 }, ploidy)
    }

    /// What each genotype's likelihood adds up to over the whole entry space at one depth. Every
    /// entry must come to one: a genotype produces *some* shape.
    fn mass_over_the_entry_space(
        model: &SsrNoiseModel,
        noise: &SlippageModel,
        ploidy: Ploidy,
        depth: u32,
    ) -> Vec<f64> {
        let genotypes = model.genotypes(ploidy);
        let mut mass = vec![0.0; genotypes];
        let mut row = Vec::with_capacity(genotypes);
        for shape in shapes_at_whole_repeat_depth(depth) {
            row.clear();
            model.append_genotype_likelihoods(&cell_of(shape, ploidy), noise, ploidy, &mut row);
            assert_eq!(row.len(), genotypes, "one likelihood per genotype");
            for (total, &ln_likelihood) in mass.iter_mut().zip(&row) {
                *total += ln_likelihood.exp();
            }
        }
        mass
    }

    /// **Gate one: does the scoring rule sum to one over the entry space?** A rule that does not
    /// is not the likelihood of anything, and no consistency result covers it — so this rejects a
    /// broken rule without fitting anything, which is the point of asking it first
    /// (`spec/parameter_prepass_ssr.md` §10).
    ///
    /// Asked per genotype, because a shape's probability is conditional on one, and over four
    /// depths because the multinomial's arrangements term changes with the depth and a rule that
    /// dropped it would still sum to one at a depth of one.
    #[test]
    fn the_scoring_rule_sums_to_one_over_the_entry_space() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(2));

        for noise in [
            measured_slippage(),
            // **Not padding.** A plug-in rescaled to sum to one is proper, so this gate is the
            // wrong instrument for it — except here: at a level of zero the rescaling divides a
            // row that is one at the allele's own bucket and zero elsewhere, and the two rules
            // then disagree about which bucket that is. Measured, dropping this row lets the
            // rescaled plug-in pass the gate.
            SlippageModel::try_new(0.0, 0.5, 0.1).expect("three probabilities"),
            SlippageModel::try_new(0.30, 0.83, 0.40).expect("three probabilities"),
            SlippageModel::try_new(0.15, 1.0, 0.0).expect("three probabilities"),
        ] {
            for depth in 1..=4 {
                for (genotype, total) in mass_over_the_entry_space(&model, &noise, ploidy(2), depth)
                    .into_iter()
                    .enumerate()
                {
                    assert!(
                        (total - 1.0).abs() < 1e-12,
                        "{noise}: genotype {genotype} at {depth} reads sums to {total}"
                    );
                }
            }
        }
    }

    /// The same gate at **ploidy four**, where a genotype is a tuple of four lengths rather than a
    /// pair and the stratum has 495 of them rather than 45.
    ///
    /// **A separate test because the diploid one cannot see the generalisation**: the average over
    /// the copies, the genotype enumeration and the genotype count all specialise to something
    /// simpler at two, and this path is specified to run at ploidy 4 as well
    /// (`spec/parameter_prepass_ssr.md` §4.2).
    #[test]
    fn the_scoring_rule_sums_to_one_at_a_tetraploid_stratum_too() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(2));
        let noise = measured_slippage();

        assert_eq!(model.genotypes(ploidy(4)), 495);
        for depth in 1..=3 {
            for (genotype, total) in mass_over_the_entry_space(&model, &noise, ploidy(4), depth)
                .into_iter()
                .enumerate()
            {
                assert!(
                    (total - 1.0).abs() < 1e-12,
                    "genotype {genotype} at {depth} reads sums to {total}"
                );
            }
        }
    }

    /// **Gate one is what rejects the tempting shortcut, and this is the shortcut.** Scoring an
    /// end bucket as though every read in it sat exactly on the edge — rather than summing over
    /// everything the bucket absorbs — is improper: the bucket probabilities do not come to one,
    /// so neither does the likelihood over the entry space.
    ///
    /// **The plug-in is built here rather than shipped**, because a rule the library never uses
    /// has no business being in the library; what the library owes is the gate that rejects it.
    ///
    /// The size of the failure depends entirely on how far the genotype's alleles sit from the
    /// recorded range, which is the design's whole point about the two widths. At the reference
    /// length the plug-in is nearly proper — the recorded range of ±4 is much wider than a read
    /// usually slips — and at an allele six copies long, which the fit is allowed to place and the
    /// entry cannot record, it collapses: the marginal puts that allele's reads in the end bucket
    /// that absorbs them, while the plug-in charges each bucket the probability of that bucket's
    /// own offset, which from an allele at +6 means losing between two and ten copies.
    #[test]
    fn scoring_an_end_bucket_at_its_edge_fails_the_sums_to_one_gate() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let noise = measured_slippage();
        let alleles = model.allele_support();
        let longest = *alleles.last().expect("a non-empty support");
        assert_eq!(longest, WholeRepeatOffset(ALLELE_OFFSET_LIMIT));

        // The rule the design adopted, at the same allele: proper, because every offset the end
        // bucket absorbs is counted into it.
        let marginal: f64 = model
            .read_bucket_probabilities(&slip_kernel(&noise), longest)
            .iter()
            .sum();
        assert!(
            (marginal - 1.0).abs() < 1e-12,
            "the marginal rule: {marginal}"
        );

        // The shortcut: bucket `b` is charged the probability of landing exactly on `b`.
        let mut plug_in = [0.0; OFFSET_BUCKETS];
        for (bucket, share) in plug_in.iter_mut().enumerate() {
            let offset =
                i32::try_from(bucket).expect("nine buckets") - i32::from(OFFSET_HALF_RANGE);
            *share = noise.probability_of_slipping_by(offset - i32::from(longest.get()));
        }
        let plugged: f64 = plug_in.iter().sum();

        // What survives the plug-in is exactly the reads that lost between two and eight copies —
        // the only slips that can land a `+6` allele inside a `−4 ..= +4` range when the bucket is
        // read as an exact offset. That is 15 reads in 10,000, so the plug-in accounts for 0.15%
        // of this allele's reads and loses the rest. It does not call them impossible: 98 reads in
        // 100 land at offset +6, and it charges the bucket holding them 670 times too little.
        // Asserted as the mechanism and not as a bound, because a bound is satisfied by any small
        // number however it arose.
        let two_or_more_copies_lost: f64 = (2..=MAX_SLIP_STEP)
            .map(|copies| noise.probability_of_slipping_by(-copies))
            .sum();
        assert!(
            (plugged - two_or_more_copies_lost).abs() < 1e-15,
            "the plug-in keeps {plugged}, and the reads that lost two or more copies are \
             {two_or_more_copies_lost}"
        );
        assert!(plugged < 0.002, "{plugged}");
    }

    /// **Gate two: is any bucket charged a negative number of reads?**
    ///
    /// **It cannot be charged one directly** — a [`LocusShape`] holds its counts in `u8` and its
    /// constructor refuses a shape whose reads exceed the cap, so a negative count is not a value
    /// this rule can be handed. What is left to check is the arithmetic side of the same
    /// statement: every genotype's bucket distribution is a distribution, so no bucket is given a
    /// negative or greater-than-one share of the locus's reads.
    #[test]
    fn no_bucket_is_charged_a_negative_share_of_a_locus_reads() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));

        for noise in [
            measured_slippage(),
            SlippageModel::try_new(1.0, 0.0, 1.0).expect("three probabilities"),
            SlippageModel::try_new(0.0, 1.0, 0.0).expect("three probabilities"),
        ] {
            for &allele in model.allele_support() {
                let shares = model.read_bucket_probabilities(&slip_kernel(&noise), allele);
                let total: f64 = shares.iter().sum();
                assert!((total - 1.0).abs() < 1e-12, "{noise} at {allele}: {total}");
                for (bucket, &share) in shares.iter().enumerate() {
                    assert!(
                        (0.0..=1.0).contains(&share),
                        "{noise} charges bucket {bucket} a share of {share} from allele {allele}"
                    );
                }
            }
        }
    }

    /// **Gate three: with the slippage level at zero, every locus's reads land on its own
    /// alleles.** A rule that puts mass anywhere else is describing movement that cannot happen.
    ///
    /// Stated as sharply as it can be: a homozygous genotype gives the shape with all its reads on
    /// that allele a likelihood of exactly one, and gives **every other shape at that depth**
    /// exactly zero — `−∞` in logs, which the trait documents as "this genotype cannot have
    /// produced this cell".
    #[test]
    fn a_silent_kernel_puts_every_read_on_its_own_allele() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let silent = SlippageModel::try_new(0.0, 0.17, 0.09).expect("three probabilities");
        let diploid = ploidy(2);
        let depth = 3;

        // The homozygote at each allele length, and the shape that has all three reads there.
        for (index, &allele) in model.allele_support().iter().enumerate() {
            let mut counts = [0; OFFSET_BUCKETS];
            counts[bucket_of(allele).index()] = depth;
            let own = LocusShape::try_new(counts, 0).expect("three reads in one bucket");

            let mut row = Vec::new();
            model.append_genotype_likelihoods(&cell_of(own, diploid), &silent, diploid, &mut row);
            // The homozygote's own index in the model's order: tuples starting at `index` begin
            // after every tuple starting lower, and the first of them is `(index, index)`.
            let homozygote =
                model.genotypes(diploid) - genotype_count(model.allele_support().len() - index, 2);
            assert!(
                (row[homozygote] - 0.0).abs() < 1e-12,
                "allele {allele}: ln L = {} where it should be ln 1",
                row[homozygote]
            );

            for shape in shapes_at_whole_repeat_depth(depth) {
                if shape == own {
                    continue;
                }
                row.clear();
                model.append_genotype_likelihoods(
                    &cell_of(shape, diploid),
                    &silent,
                    diploid,
                    &mut row,
                );
                assert_eq!(
                    row[homozygote],
                    f64::NEG_INFINITY,
                    "allele {allele}: a silent kernel made a shape it cannot produce"
                );
            }
        }
    }

    /// **The end-bucket rule agrees with the harness's, to floating point.** The harness
    /// ([`examples/ng_str_stutter_harness.rs`](../../../../../examples/ng_str_stutter_harness.rs))
    /// is what measured the rule exactly unbiased at every recorded range tried, so the library's
    /// answer has to be the harness's answer or that measurement does not transfer.
    ///
    /// Its expression is transcribed here rather than called — a `#[cfg(test)]` module cannot
    /// reach an example — and the transcription is deliberately literal: walk the slip distances,
    /// let the bucket mapping clamp, average the two copies at the end. Reading it beside
    /// [`SsrNoiseModel::read_bucket_probabilities`] is the check.
    #[test]
    fn the_bucket_probabilities_agree_with_the_harness() {
        /// `Slip::p` of `examples/shared/stutter_model.rs`, transcribed with its own closed-form
        /// renormaliser. **Not the library's kernel**: calling that on both sides would leave this
        /// test pinning the clamp loop and nothing else.
        ///
        /// The closed form is what D1 replaced, for a cancellation that only bites as the fall-off
        /// approaches one, so the two agree to about an ulp at every fall-off used here — which is
        /// why the comparison below is to `1e-15` and not to exact equality.
        fn harness_slip_probability(noise: &SlippageModel, step: i32) -> f64 {
            if step == 0 {
                return 1.0 - noise.slip_rate.get();
            }
            let copies = step.abs();
            if copies > MAX_SLIP_STEP {
                return 0.0;
            }
            let direction = if step > 0 {
                noise.gain_share.get()
            } else {
                1.0 - noise.gain_share.get()
            };
            let falloff = noise.step_decay.get();
            let tail = 1.0 - falloff.powi(MAX_SLIP_STEP);
            noise.slip_rate.get() * direction * (1.0 - falloff) * falloff.powi(copies - 1) / tail
        }

        /// `read_bucket_probs` of the same file, at `EdgeScoring::Marginal`.
        fn harness_read_bucket_probs(noise: &SlippageModel, allele: i32) -> Vec<f64> {
            let half_range = i32::from(OFFSET_HALF_RANGE);
            let mut out = vec![0.0; OFFSET_BUCKETS];
            for step in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
                let bucket = (allele + step).clamp(-half_range, half_range) + half_range;
                out[usize::try_from(bucket).expect("clamped into the range")] +=
                    harness_slip_probability(noise, step);
            }
            out
        }

        /// `genotype_bucket_probs` of the same file — the oracle the plan names: each read picks
        /// one of the two copies with equal chance, then slips.
        fn harness_genotype_bucket_probs(
            noise: &SlippageModel,
            first: i32,
            second: i32,
        ) -> Vec<f64> {
            let one = harness_read_bucket_probs(noise, first);
            let other = harness_read_bucket_probs(noise, second);
            one.iter()
                .zip(&other)
                .map(|(a, b)| 0.5 * a + 0.5 * b)
                .collect()
        }

        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let diploid = ploidy(2);
        for noise in [
            measured_slippage(),
            SlippageModel::try_new(0.30, 0.83, 0.40).expect("three probabilities"),
        ] {
            for &allele in model.allele_support() {
                let mine = model.read_bucket_probabilities(&slip_kernel(&noise), allele);
                let theirs = harness_read_bucket_probs(&noise, i32::from(allele.get()));
                for (bucket, (&ours, &reference)) in mine.iter().zip(&theirs).enumerate() {
                    assert!(
                        (ours - reference).abs() < 1e-15,
                        "{noise}, allele {allele}, bucket {bucket}: {ours} against {reference}"
                    );
                }
            }

            // And the genotype's own distribution, read off the likelihood the model returns for a
            // one-read shape: at a depth of one the arrangements term is zero, so `exp(ln L)` for
            // a shape with its single read in bucket `b` *is* that genotype's probability for `b`.
            let alleles: Vec<i32> = model
                .allele_support()
                .iter()
                .map(|offset| i32::from(offset.get()))
                .collect();
            let mut genotypes = Vec::new();
            model.for_each_genotype(diploid, |genotype| genotypes.push(genotype.to_vec()));

            for bucket in 0..OFFSET_BUCKETS {
                let mut counts = [0; OFFSET_BUCKETS];
                counts[bucket] = 1;
                let shape = LocusShape::try_new(counts, 0).expect("one read");
                let mut row = Vec::new();
                model.append_genotype_likelihoods(
                    &cell_of(shape, diploid),
                    &noise,
                    diploid,
                    &mut row,
                );

                for (genotype, &ln_likelihood) in genotypes.iter().zip(&row) {
                    let reference = harness_genotype_bucket_probs(
                        &noise,
                        alleles[genotype[0]],
                        alleles[genotype[1]],
                    )[bucket];
                    assert!(
                        (ln_likelihood.exp() - reference).abs() < 1e-15,
                        "{noise}, genotype {genotype:?}, bucket {bucket}: {} against {reference}",
                        ln_likelihood.exp()
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // D4 — the search, and the control it has to read zero on.
    // -----------------------------------------------------------------

    /// **A stratum's whole entry space, weighted by each shape's exact probability under a known
    /// truth** — the method the whole of Milestone D rests on, and the reason its answers have no
    /// tolerance to argue about.
    ///
    /// Replace each entry's observed locus count with `N ×` its probability under the truth and
    /// the sum the estimator maximises becomes the objective it would climb with an infinite
    /// genome: what comes back is the value a consistent estimator converges to, with **no
    /// sampling noise in it at all**, so "unbiased" is decided rather than estimated. Same method
    /// as [`examples/ng_str_stutter_harness.rs`](../../../../../examples/ng_str_stutter_harness.rs),
    /// against the same truth.
    ///
    /// **`N` is a billion loci and a shape is dropped once its share of them rounds below one.**
    /// A `WeightedCell` counts sites in whole numbers, so the exactness is limited by that
    /// rounding, and the limit is what the *dropped shapes add up to* rather than what each one
    /// carries. Measured on the depth-4 fixture below: 273 of the 495 shapes round to nothing and
    /// together hold **1.1e-8** of the mass, while the rounding over all 495 comes to **6.4e-8**.
    /// That is three orders below the search's own resolution and nearly six below the 0.005 the
    /// shares are asserted to, so it cannot be what any assertion here is reading.
    fn exact_entry_space(
        model: &SsrNoiseModel,
        truth: &SlippageModel,
        allele_frequencies: &[(WholeRepeatOffset, f64)],
        ploidy: Ploidy,
        depth: u32,
    ) -> Vec<StratumCell> {
        const LOCI: f64 = 1e9;

        // The genotype frequencies the truth implies, under Hardy-Weinberg over the alleles
        // named: `p²` for a homozygote and `2pq` for a heterozygote, in the model's own genotype
        // order so they line up with the likelihoods it returns.
        let support = model.allele_support();
        let share_of = |index: usize| {
            allele_frequencies
                .iter()
                .find(|(offset, _)| *offset == support[index])
                .map_or(0.0, |(_, share)| *share)
        };
        let mut genotype_frequencies = Vec::new();
        model.for_each_genotype(ploidy, |genotype| {
            let product: f64 = genotype.iter().map(|&allele| share_of(allele)).product();
            // How many orderings of this multiset there are — 1 for a homozygote, 2 for a
            // diploid heterozygote.
            let orderings = if genotype.windows(2).all(|pair| pair[0] == pair[1]) {
                1.0
            } else {
                2.0
            };
            genotype_frequencies.push(orderings * product);
        });
        let total: f64 = genotype_frequencies.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "the truth's genotype frequencies sum to {total}"
        );

        let mut cells = Vec::new();
        let mut row = Vec::new();
        for shape in shapes_at_whole_repeat_depth(depth) {
            row.clear();
            model.append_genotype_likelihoods(&cell_of(shape, ploidy), truth, ploidy, &mut row);
            let probability: f64 = row
                .iter()
                .zip(&genotype_frequencies)
                .map(|(ln_likelihood, frequency)| frequency * ln_likelihood.exp())
                .sum();
            let loci = (probability * LOCI).round() as u64;
            if loci > 0 {
                cells.push(StratumCell::new(StratumEntry { shape, loci }, ploidy));
            }
        }
        cells
    }

    /// **The control the whole method rests on: generate under the truth, fit, and get the truth
    /// back.** With the entry space weighted by its exact probabilities there is no sampling noise
    /// for a residual to hide in, so anything other than recovery is this code's fault.
    ///
    /// **And the spread is not read alone**, which is the trap that cost a published finding on
    /// 2026-08-06: a deterministic search returns the same point from every start wherever the
    /// objective is flat, so four starts agreeing is also what a search that never looked
    /// produces. It is paired here with the score at the truth — one extra evaluation, and a
    /// correctly specified model cannot be beaten at its own truth, so **a fitted point scoring
    /// above it would be a defect in this test rather than a finding about the estimator**.
    ///
    /// **Ignored by default, and run by hand — measured at 12 minutes in a debug build**
    /// (719.6 s of test time on an 8-core Apple `container` VM):
    ///
    /// ```text
    /// ./scripts/dev.sh cargo test --lib the_search_recovers_a_known_truth -- --ignored
    /// ```
    ///
    /// Every candidate the search tries runs a whole climb over this stratum's 36 genotypes and
    /// the 222 cells its entry space leaves after the rounding below, and the search tries
    /// several hundred: four starts, three axes, and about 29 golden-section steps an axis —
    /// which is where the bracket reaches the tolerance, not the cap of forty — over up to eight
    /// sweeps. **Shortening the search is what a control may least afford**: at a coarser
    /// resolution a genuine 1% bias and the search's own step are the same size, so it would be
    /// reporting its own grid rather than the estimator. Reported at Checkpoint D rather than
    /// quietly weakened.
    #[test]
    #[ignore = "nine minutes in a debug build; the sharp control, run by hand"]
    fn the_search_recovers_a_known_truth_and_no_start_beats_it() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(1));
        let truth = SlippageModel::try_new(0.0201, 0.17, 0.09).expect("three probabilities");
        let diploid = ploidy(2);
        // Three allele lengths at a 30% non-reference **allele** frequency, which under the
        // Hardy-Weinberg weighting below is **51 loci in 100** carrying something other than the
        // reference — enough that the fit has to tell a long allele from a slipped read rather
        // than having only one story available.
        let alleles = [
            (WholeRepeatOffset(-1), 0.15),
            (WholeRepeatOffset(0), 0.70),
            (WholeRepeatOffset(1), 0.15),
        ];
        let cells = exact_entry_space(&model, &truth, &alleles, diploid, 4);
        assert!(
            cells.len() > 20,
            "the entry space collapsed to {} shapes, so this control is not exercising one",
            cells.len()
        );

        // The starts are handed the share of reads sitting off the reference length, which is
        // what a real driver has before it fits anything — an over-estimate of the level, since
        // it counts real alleles too.
        let off_reference: u64 = cells
            .iter()
            .map(|cell| u64::from(cell.shape().reads_off_reference()) * cell.loci())
            .sum();
        let reads: u64 = cells
            .iter()
            .map(|cell| u64::from(cell.shape().whole_repeat_depth()) * cell.loci())
            .sum();
        let starts = slippage_starts(off_reference as f64 / reads as f64);

        let fitted = fit_by_multistart(&model, &cells, &starts, SearchPrecision::fine());

        assert_eq!(fitted.starts.len(), 4, "every start is reported");
        assert!(
            fitted.termination.converged,
            "a start ran out of sweeps after {}",
            fitted.termination.iterations
        );

        let level = fitted.best.slip_rate.get();
        let gain = fitted.best.gain_share.get();
        let decay = fitted.best.step_decay.get();
        assert!(
            (level / truth.slip_rate.get() - 1.0).abs() < 0.01,
            "the level came back at {level} against a truth of {}",
            truth.slip_rate.get()
        );
        assert!(
            (gain - truth.gain_share.get()).abs() < 0.005,
            "the direction split came back at {gain} against a truth of 0.17"
        );
        assert!(
            (decay - truth.step_decay.get()).abs() < 0.005,
            "the fall-off came back at {decay} against a truth of 0.09"
        );

        // **Paired with the score at the truth, and this is the half that is not optional.** The
        // spread says the four starts agreed; the score says they agreed on the right thing.
        let at_truth = score_of(&model, &cells, &truth);
        assert!(
            fitted.log_likelihood.get() <= at_truth + 1e-6,
            "the fitted point scores {} against {at_truth} at the truth, so this test is wrong \
             rather than the estimator",
            fitted.log_likelihood.get()
        );
        assert!(
            fitted.headline_spread < 1.06,
            "the four starts reached levels spanning {}-fold",
            fitted.headline_spread
        );
    }

    /// **What the search reports beside its answer, which is the half a consumer judges it on.**
    /// Cheap enough for the suite — one sweep over a two-read entry space — because what it pins
    /// is the reporting and not the recovery: every start recorded, best-scoring first, the spread
    /// measured on the level, and how the search ended.
    ///
    /// **The sharp recovery control is
    /// [`the_search_recovers_a_known_truth_and_no_start_beats_it`], which is ignored by default**
    /// and takes nine minutes. This test would pass on a search that stopped after one sweep in
    /// the wrong place; that one would not.
    #[test]
    fn every_start_is_reported_with_its_score_and_the_spread_across_them() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(1));
        let truth = SlippageModel::try_new(0.0201, 0.17, 0.09).expect("three probabilities");
        let diploid = ploidy(2);
        let alleles = [(WholeRepeatOffset(0), 1.0)];
        let cells = exact_entry_space(&model, &truth, &alleles, diploid, 2);

        let starts = slippage_starts(0.03);
        assert_eq!(starts.len(), 4);
        // **Each start disagrees with the next about all three parameters at once**, which is the
        // property the set exists for: five starts that disagreed about the headline number while
        // sharing one guess at a nuisance axis returned a confident zero on the sibling path
        // (`spec/parameter_prepass_ssr.md` §4.2).
        for (earlier, later) in starts.iter().zip(starts.iter().skip(1)) {
            assert_ne!(earlier.slip_rate.get(), later.slip_rate.get());
            assert_ne!(earlier.gain_share.get(), later.gain_share.get());
            assert_ne!(earlier.step_decay.get(), later.step_decay.get());
        }

        let fitted = fit_by_multistart(
            &model,
            &cells,
            &starts,
            // **Two golden-section steps an axis and one sweep**, which is all this test needs:
            // it asks what the search *reports*, not where it lands. Each step is a whole climb
            // over the stratum's 36 genotypes, so the cost is set here and nowhere else.
            SearchPrecision {
                tolerance: 0.05,
                max_axis_steps: 2,
                max_sweeps: 1,
            },
        );

        assert_eq!(
            fitted.starts.len(),
            4,
            "every start is reported, not just the winner"
        );
        for (better, worse) in fitted.starts.iter().zip(fitted.starts.iter().skip(1)) {
            assert!(
                better.log_likelihood >= worse.log_likelihood,
                "the starts are not ordered best-scoring first"
            );
        }
        assert_eq!(
            fitted.best, fitted.starts[0].reached,
            "the winner is the best-scoring start"
        );
        assert!(
            (fitted.log_likelihood.get() - fitted.starts[0].log_likelihood.get()).abs() < 1e-9,
            "the reported score is not the winner's"
        );

        // Every start began somewhere different, and the report says where.
        let began: Vec<f64> = fitted
            .starts
            .iter()
            .map(|start| start.from.slip_rate.get())
            .collect();
        for (earlier, later) in began.iter().zip(began.iter().skip(1)) {
            assert_ne!(earlier, later);
        }

        // The spread is the highest-to-lowest ratio of the **level** across what the starts
        // reached — a ratio, because the level spans orders of magnitude across a real ladder of
        // strata.
        let reached: Vec<f64> = fitted
            .starts
            .iter()
            .map(|start| start.reached.slip_rate.get())
            .collect();
        let highest = reached.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let lowest = reached.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            (fitted.headline_spread - highest / lowest).abs() < 1e-9,
            "the spread reads {} where the levels reached span {highest} to {lowest}",
            fitted.headline_spread
        );

        // One sweep is not enough to settle, and the search says so rather than presenting where
        // it stopped as where it converged.
        assert_eq!(fitted.termination.iterations, 1);
    }

    /// What one set of slippage parameters scores over the whole table, with the genotype
    /// frequencies climbed at it — a search of one start and no sweeps, which is the cheapest way
    /// to ask the question through exactly the code the search asks it through.
    fn score_of(model: &SsrNoiseModel, cells: &[StratumCell], noise: &SlippageModel) -> f64 {
        fit_by_multistart(
            model,
            cells,
            std::slice::from_ref(noise),
            SearchPrecision {
                tolerance: 1e-5,
                max_axis_steps: 0,
                max_sweeps: 0,
            },
        )
        .log_likelihood
        .get()
    }

    /// **What the fit weighs a cell by is how many loci showed that shape — not how many reads
    /// they held.** Nothing inside this module reads it, so only a test through the trait does;
    /// measured, replacing it with the shape's read depth left all 567 tests of this path green,
    /// and it is the substitution that turns a table of loci back into a tally of reads, which is
    /// the keying the design rejects at a **333-fold** spread in the fitted level
    /// (`spec/parameter_prepass_ssr.md` §4.1).
    ///
    /// Three distinct numbers — 2 loci, 5 whole-repeat reads, 3 in the guard — so every wrong
    /// answer differs from the right one.
    #[test]
    fn a_cell_counts_the_loci_that_showed_its_shape_and_not_their_reads() {
        let mut counts = [0; OFFSET_BUCKETS];
        counts[bucket_of(WholeRepeatOffset(0)).index()] = 4;
        counts[bucket_of(WholeRepeatOffset(-1)).index()] = 1;
        let shape = LocusShape::try_new(counts, 3).expect("eight reads");
        assert_eq!(shape.whole_repeat_depth(), 5);
        assert_eq!(shape.depth(), 8);

        let cell = StratumCell::new(StratumEntry { shape, loci: 2 }, ploidy(2));

        assert_eq!(WeightedCell::sites(&cell), 2);
        assert_eq!(WeightedCell::ploidy(&cell), ploidy(2));
    }

    /// **A stratum's genotypes are its own allele lengths' unordered tuples**, and the count
    /// follows from the support rather than from the ploidy: 45 at nine lengths, 91 at the full
    /// thirteen, and 66 at a stratum of four reference copies, whose support clips at −4 because
    /// an allele cannot be shorter than nothing.
    ///
    /// **These are the numbers the sibling path's `ploidy + 1` would have got wrong**, which is
    /// why the trait asks the model. Three for a diploid is the right answer there and the wrong
    /// one here by a factor of thirty at the widest stratum and fifteen at the narrowest.
    #[test]
    fn a_stratum_has_one_genotype_per_unordered_tuple_of_its_own_allele_lengths() {
        for (reference_repeats, lengths, diploid_genotypes) in
            [(2, 9, 45), (4, 11, 66), (6, 13, 91)]
        {
            let model = SsrNoiseModel::for_stratum(RepeatCount(reference_repeats));
            assert_eq!(model.allele_support().len(), lengths);
            assert_eq!(model.genotypes(ploidy(2)), diploid_genotypes);
        }

        // A haploid stratum has one genotype per length, which is the identity the genotype count
        // has to reduce to and the one a `C(A+P-1, P)` written the wrong way round loses.
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        assert_eq!(model.genotypes(ploidy(1)), 13);

        // The widest support any stratum can have is the width the per-cell scratch array is
        // sized at. Undersizing it panics on an index; oversizing it is silent and merely wastes
        // stack, so the equality is asserted rather than the inequality.
        assert_eq!(MAX_ALLELE_LENGTHS, 13);
        assert_eq!(model.allele_support().len(), MAX_ALLELE_LENGTHS);
    }

    /// **Each genotype appears exactly once, spelled shortest-allele-first.** A genotype is an
    /// unordered tuple — a locus carrying the reference length and one copy short is one genotype,
    /// not two — so emitting both spellings would double-count it and leave the fitted frequencies
    /// summing to more than one.
    #[test]
    fn every_genotype_is_emitted_once_in_non_decreasing_order() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(2));

        for copies in [1u8, 2, 3] {
            let mut seen = Vec::new();
            model.for_each_genotype(ploidy(copies), |genotype| {
                assert!(
                    genotype.windows(2).all(|pair| pair[0] <= pair[1]),
                    "{genotype:?} is not in non-decreasing order"
                );
                assert_eq!(genotype.len(), usize::from(copies));
                seen.push(genotype.to_vec());
            });

            assert_eq!(seen.len(), model.genotypes(ploidy(copies)));
            let mut unique = seen.clone();
            unique.sort();
            unique.dedup();
            assert_eq!(unique.len(), seen.len(), "a genotype was emitted twice");
            assert!(
                seen.windows(2).all(|pair| pair[0] < pair[1]),
                "not ascending"
            );
        }
    }

    /// **The guard reads are not in the length likelihood**, which is what lets the guard be a
    /// diagnostic rather than a parameter: the likelihood splits exactly into how many reads
    /// showed a non-whole-repeat length times how the rest fell across the buckets, so nothing
    /// about slippage is estimated from the guard and nothing about the guard disturbs it
    /// (`spec/parameter_prepass_ssr.md` §4.1).
    ///
    /// Two shapes with the same bucket counts and different guard counts must therefore score
    /// identically — and a locus whose every read landed in the guard scores an empty product, so
    /// it says nothing about any genotype rather than saying the fit is wrong.
    #[test]
    fn the_guard_reads_do_not_enter_the_length_likelihood() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let noise = measured_slippage();
        let diploid = ploidy(2);

        let mut counts = [0; OFFSET_BUCKETS];
        counts[bucket_of(WholeRepeatOffset(0)).index()] = 3;
        counts[bucket_of(WholeRepeatOffset(-1)).index()] = 1;

        let without = LocusShape::try_new(counts, 0).expect("four reads");
        let with = LocusShape::try_new(counts, 5).expect("four reads and five in the guard");

        let mut plain = Vec::new();
        let mut guarded = Vec::new();
        model.append_genotype_likelihoods(&cell_of(without, diploid), &noise, diploid, &mut plain);
        model.append_genotype_likelihoods(&cell_of(with, diploid), &noise, diploid, &mut guarded);
        assert_eq!(plain, guarded);

        // Every read in the guard: an empty product, which is a likelihood of one for every
        // genotype. Not zero — a locus with no length evidence is not a locus that refutes
        // something.
        let all_guard = LocusShape::try_new([0; OFFSET_BUCKETS], 4).expect("four guard reads");
        let mut silent = Vec::new();
        model.append_genotype_likelihoods(
            &cell_of(all_guard, diploid),
            &noise,
            diploid,
            &mut silent,
        );
        assert!(silent.iter().all(|&ln_likelihood| ln_likelihood == 0.0));
    }

    /// **A cell scored against the genotypes of another ploidy is refused rather than answered.**
    /// The two arrive separately through the trait, so nothing but this assertion stands between a
    /// diploid locus and the 1,820 genotypes this stratum's thirteen allele lengths make at ploidy
    /// four — and the answer would be a plausible number with no way to question it.
    #[test]
    #[should_panic(expected = "a cell of ploidy 2 scored against the genotypes of ploidy 4")]
    fn a_cell_cannot_be_scored_against_another_ploidy_genotypes() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let mut counts = [0; OFFSET_BUCKETS];
        counts[bucket_of(WholeRepeatOffset(0)).index()] = 2;
        let shape = LocusShape::try_new(counts, 0).expect("two reads");

        model.append_genotype_likelihoods(
            &cell_of(shape, ploidy(2)),
            &measured_slippage(),
            ploidy(4),
            &mut Vec::new(),
        );
    }

    /// **A heterozygote is the average of its two copies, not the product and not the sum.** Each
    /// read picks one copy with equal chance and then slips, so the bucket probabilities of the
    /// genotype `(a, b)` sit exactly halfway between those of `(a, a)` and `(b, b)`.
    ///
    /// Worth its own test because the three arithmetics agree at a homozygote — where both copies
    /// are the same allele — so every check that uses only homozygous genotypes is blind to which
    /// one was written.
    #[test]
    fn a_heterozygote_is_the_average_of_its_two_copies() {
        let model = SsrNoiseModel::for_stratum(RepeatCount(6));
        let noise = measured_slippage();
        let diploid = ploidy(2);

        // One read in the bucket one copy short of the reference.
        let mut counts = [0; OFFSET_BUCKETS];
        counts[bucket_of(WholeRepeatOffset(-1)).index()] = 1;
        let shape = LocusShape::try_new(counts, 0).expect("one read");

        let mut row = Vec::new();
        model.append_genotype_likelihoods(&cell_of(shape, diploid), &noise, diploid, &mut row);

        // The support runs −6..=6, so the reference length is index 6 and one copy short is 5.
        let alleles = model.allele_support();
        assert_eq!(alleles[6], WholeRepeatOffset(0));
        assert_eq!(alleles[5], WholeRepeatOffset(-1));

        let index_of = |first: usize, second: usize| {
            let mut found = None;
            let mut at = 0;
            model.for_each_genotype(diploid, |genotype| {
                if genotype == [first, second] {
                    found = Some(at);
                }
                at += 1;
            });
            found.expect("a genotype of this stratum")
        };
        let short_short = index_of(5, 5);
        let short_reference = index_of(5, 6);
        let reference_reference = index_of(6, 6);

        let average = 0.5 * (row[short_short].exp() + row[reference_reference].exp());
        assert!(
            (row[short_reference].exp() - average).abs() < 1e-15,
            "the heterozygote scores {} where the average of its copies is {average}",
            row[short_reference].exp()
        );

        // **The line above pins linearity and not the half**, which is the trap: at one read the
        // arrangements term is zero, so `L(a,b) = ½(L(a,a) + L(b,b))` holds for any constant
        // multiple of the truth — a rule that summed the copies instead of averaging them scores
        // every genotype twice over and still satisfies it. Anchoring one homozygote against the
        // single-copy distribution it is built from is what fixes the scale.
        let reference_bucket = bucket_of(WholeRepeatOffset(-1)).index();
        let from_one_copy = model
            .read_bucket_probabilities(&slip_kernel(&noise), WholeRepeatOffset(0))
            [reference_bucket];
        assert!(
            (row[reference_reference].exp() - from_one_copy).abs() < 1e-15,
            "a homozygote's one read scores {} where one copy gives it {from_one_copy}",
            row[reference_reference].exp()
        );
    }
}
