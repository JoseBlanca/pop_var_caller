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

use crate::ng::types::{DomainError, checked_probability};

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

    /// The slippage models the kernel tests sweep, chosen so that each of them **breaks
    /// something** rather than to look like real data.
    ///
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
}
