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
//! Empty of mathematics until Milestone D; A4 lands the parameters it will search over.

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
