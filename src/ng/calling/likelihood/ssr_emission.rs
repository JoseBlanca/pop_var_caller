//! The STR emission seam — the one surface of this step that swaps.
//!
//! *How probable is one observed sequence, given one candidate allele?* At a repeat tract
//! that question has more than one defensible answer, and the caller must be able to run a
//! second one against the default without touching anything around it. This module holds
//! the **seam**: what an emission is handed, what it returns, and nothing about how it
//! decides (`doc/devel/ng/spec/read_likelihoods.md` §2.4, §4.1).
//!
//! Everything the model reads arrives **per call**, in [`SsrScoringContext`]. Nothing is read
//! from global state and nothing is fitted here — the whole point of the seam is that the
//! EM loop can re-estimate the slippage numbers between iterations with no change on this
//! side (spec §6.1).
//!
//! # The two questions a model answers, and why they are separate methods
//!
//! A read that spanned the whole tract pins what the sample carries there. A read that
//! entered the tract and **ran off its own end** proves only that the tract is *at least* as
//! long as what it showed — the statistician's word is **censored**. Those are different
//! questions about the same candidate, so they are [`SsrEmissionModel::emission`] and
//! [`SsrEmissionModel::censored_emission`]. Routing between them is the row's job, from the
//! observation's witness; a model never inspects one.
//!
//! # What a context is built per, and what must not be hoisted
//!
//! **Per `(read group, candidate)`** — not per locus. A read's chance of slipping is a
//! property of the tract it was copied from, and that is the **candidate** allele: a
//! candidate of 6 repeats and one of 12 at the same locus are drawn from different strata and
//! slip at measurably different rates, about 1.3-fold per repeat count over the measured
//! range (spec §4.4). So the stutter parameters cannot be hoisted out of the candidate loop.
//! **The lookup can** — it is a small table indexed by period and repeat count — and that is
//! the distinction a coder has to keep.

use std::num::{NonZeroU8, NonZeroU32};

use crate::ng::alignment::StutterModel;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::types::{ErrorRate, Motif};

/// One candidate allele, as the emission sees it.
///
/// Two fields, and the second is not derivable from the first: **the repeat count keys the
/// stratum lookup**, and it is the candidate's own count rather than the reference's
/// (spec §4.4). Counting it from `bases` would mean re-measuring a tract the locus generator
/// has already measured, which is the duplication spec §7 puts on the alignment module's side
/// of the boundary — this type *consumes* a measurement and never makes one.
#[derive(Debug, Clone, Copy)]
pub struct SsrCandidate<'a> {
    /// The whole locus as a sample carrying this allele has it — flanks included, not the
    /// tract alone, and not the reference with the tract swapped out.
    pub bases: &'a [u8],
    /// How many repeats this candidate's tract holds. **Non-zero**: a candidate whose tract
    /// holds no repeats is not a candidate, which is the same contract
    /// [`StutterModel::unreachable_mass`] states from the distribution's side.
    pub repeat_count: NonZeroU32,
}

/// Everything a model is handed for one `(read group, candidate)`, and the only channel it
/// has.
///
/// **The tier-two seam** (spec §6.1): every number here may be re-estimated between the
/// caller's iterations, and the emission never asks where any of them came from. That is what
/// makes the EM loop's re-fitting a change of nobody's code but its own.
#[derive(Debug, Clone, Copy)]
pub struct SsrScoringContext<'a> {
    /// The repeating unit, whose length is the period the stutter distribution is indexed by.
    pub motif: &'a Motif,
    /// How likely each length change is, for **this candidate's** stratum — built per
    /// `(read group, candidate)` and never hoisted out of the candidate loop (spec §4.4).
    pub stutter: &'a StutterModel,
    /// The per-base substitution rate for this read group and stratum.
    ///
    /// **Never the SNP/indel path's ε, and never a read's own summed quality** — spec §4.3's
    /// closed question Q6, and the reason is a unit mismatch rather than a preference. A
    /// read's error probability is a per-*read* number, the chance it is wrong somewhere; the
    /// substitution term needs a per-*base* rate, applied once for each of the tract's twenty
    /// or forty bases. Using the first as the second overcharges by the tract's length.
    ///
    /// The two rates are separate fitted parameters that are never tied: each absorbs what
    /// its own model cannot otherwise explain, and forcing one number to carry both would
    /// make each model wrong in a way neither could report.
    pub substitution_rate: ErrorRate,
    /// The mass the stutter distribution cannot place for this candidate
    /// ([`StutterModel::unreachable_mass`]) — computed and carried, never assumed negligible.
    ///
    /// **It travels because the row compares candidates.** A model that loses mass on some
    /// candidates and not others is comparing them on different scales, and at period 1 the
    /// loss is 2 in 100 rather than the 1-in-10¹³ a cutoff tail costs (spec §4.2).
    pub unreachable_mass: f64,
    /// The weakest warrant behind any parameter in this context.
    ///
    /// **The model never branches on it; it propagates** (spec §4.4). A stratum whose numbers
    /// were borrowed is used exactly as a fitted one, with no down-weighting — a borrowed
    /// value is the best estimate available and discounting it would mean inventing a
    /// penalty. But the fact travels, so a call resting on borrowed parameters is
    /// distinguishable in the run's output from one resting on a fit, without re-running
    /// anything.
    pub weakest_provenance: Provenance,
}

impl<'a> SsrScoringContext<'a> {
    /// Build a context for one candidate, taking the unreachable mass from the distribution
    /// rather than from the caller.
    ///
    /// **The mass and the model must agree**, and letting a caller pass both separately is
    /// two chances to disagree. The only inputs that are genuinely the caller's are the
    /// parameters' warrants, which come from the fits the numbers were read out of —
    /// combined here with [`Provenance::weaker_of`], because a context resting on one fitted
    /// number and one borrowed number is a borrowed context.
    #[must_use]
    pub fn new(
        motif: &'a Motif,
        stutter: &'a StutterModel,
        candidate: &SsrCandidate<'_>,
        substitution_rate: ErrorRate,
        parameter_warrants: impl IntoIterator<Item = Provenance>,
    ) -> Self {
        let weakest_provenance = parameter_warrants
            .into_iter()
            .fold(Provenance::FittedHere, Provenance::weaker_of);
        Self {
            motif,
            stutter,
            substitution_rate,
            unreachable_mass: stutter.unreachable_mass(period_of(motif), candidate.repeat_count),
            weakest_provenance,
        }
    }
}

/// A motif's period as the stutter distribution wants it.
///
/// [`Motif::new`] rejects an empty motif and one longer than six bases, so the period is in
/// `1..=6` by construction and neither conversion can fail — which is what
/// [`StutterModel::probability`]'s doc means when it says "the conversion at a real call site
/// cannot fail". Done once here rather than at every consumer.
#[inline]
fn period_of(motif: &Motif) -> NonZeroU8 {
    let period = u8::try_from(motif.period()).expect("a motif is at most six bases");
    NonZeroU8::new(period).expect("a motif is at least one base")
}

/// `Lr(observation | one candidate allele)` — **the only part that differs between models**.
///
/// Everything around it is shared: the copy-weighted mixture over a genotype's alleles, the
/// outlier term, the caching, the logarithm. A second model swaps this trait and nothing
/// else, which is what makes the comparison behind spec §4.1 a fair one.
///
/// # Probability space, not log space
///
/// Both methods return a **linear probability**, floored. The row takes one logarithm per
/// observation per genotype, after the mixture — putting the log inside would mean taking it
/// per allele instead, and spec §2.1's junk term needs a logarithm around a *sum* over
/// alleles, which a per-allele logarithm cannot express.
///
/// # The scratch is the model's own
///
/// An implementation that needs working memory declares its shape; the row owns one and
/// hands it back on every call, so nothing allocates per observation per candidate. A model
/// that needs none says `type Scratch = ()`.
pub trait SsrEmissionModel {
    /// Working memory this model reuses across calls. `()` for a model that needs none.
    type Scratch: Default;

    /// The probability that one copy of `candidate` produced this whole observed sequence.
    ///
    /// `observation` is the read's own bases over the locus — **the whole locus as the read
    /// saw it**, on the same footing as [`SsrCandidate::bases`], so the two are comparable
    /// without either being re-measured.
    fn emission(
        &self,
        observation: &[u8],
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64;

    /// The probability that one copy of `candidate` produced a read that showed **at least**
    /// this much and then ran out — `P(length ≥ what was witnessed | candidate)` times the
    /// letter match over what was witnessed (spec §5.2).
    ///
    /// **Not a shorter complete observation.** A read that ran out inside the tract has not
    /// shown a shorter allele; it has shown a lower bound, and scoring it as though it pinned
    /// a length would let a truncated read out-discriminate a whole one.
    fn censored_emission(
        &self,
        witnessed_prefix: &[u8],
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;

    fn a_motif(bases: &[u8]) -> Motif {
        Motif::new(bases).expect("a valid test motif")
    }

    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a test candidate always holds a repeat")
    }

    fn a_model() -> StutterModel {
        stutter_model_for(&Slippage {
            level: 0.02,
            shorter_share: 0.83,
            fall_off: 0.35,
        })
    }

    /// **The context takes its unreachable mass from the distribution, not from the caller.**
    /// Two ways in would be two chances to disagree, and the row compares candidates on the
    /// strength of this number.
    #[test]
    fn the_context_reads_the_unreachable_mass_off_the_distribution() {
        let motif = a_motif(b"CAG");
        let model = a_model();
        let candidate = SsrCandidate {
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(4),
        };
        let context = SsrScoringContext::new(
            &motif,
            &model,
            &candidate,
            ErrorRate::try_new(0.001).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        assert_eq!(
            context.unreachable_mass,
            model.unreachable_mass(period_of(&motif), repeats(4))
        );
    }

    /// **The mass differs between candidates at one locus**, which is the whole reason a
    /// context is per candidate rather than per locus. A four-repeat tract cannot lose as
    /// many repeats as a thirty-repeat one, so it leaves more unplaced.
    #[test]
    fn two_candidates_at_one_locus_get_different_contexts() {
        let motif = a_motif(b"CAG");
        let model = a_model();
        let rate = ErrorRate::try_new(0.001).expect("a valid rate");

        let short = SsrCandidate {
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(4),
        };
        let long = SsrCandidate {
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(30),
        };
        let short_context =
            SsrScoringContext::new(&motif, &model, &short, rate, [Provenance::FittedHere]);
        let long_context =
            SsrScoringContext::new(&motif, &model, &long, rate, [Provenance::FittedHere]);

        assert!(
            short_context.unreachable_mass > long_context.unreachable_mass,
            "four repeats left {} unplaced, thirty left {}",
            short_context.unreachable_mass,
            long_context.unreachable_mass
        );
    }

    /// **The weakest warrant wins, and one borrowed parameter is enough.** A context resting
    /// on a fitted rate and a borrowed slippage row is a borrowed context; stamping it
    /// `FittedHere` would launder the weaker of the two.
    #[test]
    fn the_context_carries_the_weakest_warrant_that_entered_it() {
        let motif = a_motif(b"CA");
        let model = a_model();
        let candidate = SsrCandidate {
            bases: b"AACACACATT",
            repeat_count: repeats(3),
        };
        let rate = ErrorRate::try_new(0.001).expect("a valid rate");

        for (warrants, expected) in [
            (
                vec![Provenance::FittedHere, Provenance::FittedHere],
                Provenance::FittedHere,
            ),
            (
                vec![Provenance::FittedHere, Provenance::Borrowed],
                Provenance::Borrowed,
            ),
            (
                vec![Provenance::Borrowed, Provenance::Supplied],
                Provenance::Supplied,
            ),
            (
                vec![Provenance::Supplied, Provenance::Defaulted],
                Provenance::Defaulted,
            ),
            (
                vec![Provenance::Defaulted, Provenance::FittedHere],
                Provenance::Defaulted,
            ),
        ] {
            let context =
                SsrScoringContext::new(&motif, &model, &candidate, rate, warrants.clone());
            assert_eq!(
                context.weakest_provenance, expected,
                "warrants {warrants:?} gave {:?}",
                context.weakest_provenance
            );
        }
    }

    /// **No warrants at all is `FittedHere`**, which is the identity of the fold rather than
    /// an opinion: a context nothing weakened is as well founded as its inputs. Pinned
    /// because the alternative — defaulting to `Defaulted` — would mark every context of a
    /// fully-fitted run as a guess.
    #[test]
    fn a_context_with_no_weakening_warrant_is_fitted() {
        let motif = a_motif(b"CA");
        let model = a_model();
        let candidate = SsrCandidate {
            bases: b"AACACACATT",
            repeat_count: repeats(3),
        };
        let context = SsrScoringContext::new(
            &motif,
            &model,
            &candidate,
            ErrorRate::try_new(0.001).expect("a valid rate"),
            [],
        );
        assert_eq!(context.weakest_provenance, Provenance::FittedHere);
    }
}
