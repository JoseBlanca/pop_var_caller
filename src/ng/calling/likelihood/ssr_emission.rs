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
use crate::ng::alignment::emission::{Emission, FlatEmission};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::types::{BaseQual, ErrorRate, Motif};

/// One candidate allele, as the emission sees it.
///
/// Two fields, and the second is not derivable from the first: **the repeat count keys the
/// stratum lookup**, and it is the candidate's own count rather than the reference's
/// (spec §4.4). Counting it from `bases` would mean re-measuring a tract the locus generator
/// has already measured, which is the duplication spec §7 puts on the alignment module's side
/// of the boundary — this type *consumes* a measurement and never makes one. **Nor is it
/// derivable**: an interrupted tract's byte length divided by the period is not its repeat
/// count.
#[derive(Debug, Clone, Copy)]
pub struct SsrCandidate<'a> {
    /// **The tract, and only the tract** — the repeat run as a carrier of this allele has
    /// it, without the flanks.
    ///
    /// **This differs from the generic path deliberately, and the difference is the locus's
    /// own shape.** There an allele is the whole locus as a carrier has it, because the
    /// varying regions sit inside surrounding context that alleles disagree about. Here the
    /// locus *is* the tract: the STR generator slices exactly the tract into an observation's
    /// bases ([`locus_generation/ssr.rs`](../../../../src/ng/locus_generation/ssr.rs)), the
    /// flanks are shared by every candidate by construction, and comparing them would be
    /// comparing two copies of the same context.
    ///
    /// So an observation and a candidate are on the same footing here, which is what lets the
    /// substitution term compare them once the stutter factor has equalised their lengths.
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
    /// `observation` is the read's own bases over the **tract**, which is what the STR
    /// generator slices — the same footing as [`SsrCandidate::bases`], so the two are
    /// comparable without either being re-measured.
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

/// **Model A** — the default STR emission: how likely a length change is, times how likely
/// the letters are.
///
/// # Two factors, each unable to answer the other's question
///
/// The **stutter factor** decides *how long* the read is: it stretches or trims the candidate
/// to the read's own length, using the distribution
/// ([`StutterModel`]) the pre-pass fitted for this candidate's stratum. The **substitution
/// factor** then decides *which letters*, over two sequences that are the same length by
/// construction — so there is nothing left for it to insert or delete.
///
/// **Keeping them apart is what makes their two parameters separately estimable.** If the
/// substitution factor could also delete bases inside the tract, a read one repeat short
/// would have two explanations — a polymerase slip, or a sequencing deletion — and nothing in
/// the data would choose between them: raising the slippage rate and lowering the indel rate
/// would describe the reads exactly as well as the reverse (spec §4.3).
///
/// # The substitution factor is composed, never re-implemented
///
/// It is [`FlatEmission`] under the context's fitted per-stratum rate — `1 − ε` per matching
/// base and `ε/3` per mismatching one, the alignment module's own comparison
/// (spec §4.3, §7, §9). This model resizes and calls; it does not score letters itself.
///
/// # Where a slip can land
///
/// In a pure tract, adding a repeat anywhere gives the same sequence and there is nothing to
/// choose. In an **interrupted** tract the placements give genuinely different sequences, so
/// the model enumerates them and averages with equal weight.
///
/// **ng enumerates for whole-repeat slips only, and resizes a part-repeat change at the
/// tract's end in a single placement.** That is production's split, stated here rather than
/// left for a coder to guess, and it is a simplification of the same class as the fixed
/// part-repeat share (spec §4.2).
///
/// **A slip that cannot be reached at all scores zero** — contracting away more repeats than
/// a run holds. The mass that costs is not assumed away: it is what
/// [`SsrScoringContext::unreachable_mass`] reports.
#[derive(Debug, Default, Clone, Copy)]
pub struct StutterSubstitutionEmission;

/// Working memory Model A reuses across calls: the placements a whole-repeat slip could land
/// in, and the buffer a part-repeat resize is rendered into.
///
/// Held by the row and handed back on every call, so nothing allocates per observation per
/// candidate.
#[derive(Debug, Default)]
pub struct StutterSubstitutionScratch {
    /// One entry per placement-distinct realisation of the slip.
    placements: Vec<Vec<u8>>,
    /// The candidate resized to the read's length, for the single-placement branch.
    resized: Vec<u8>,
}

impl SsrEmissionModel for StutterSubstitutionEmission {
    type Scratch = StutterSubstitutionScratch;

    fn emission(
        &self,
        observation: &[u8],
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64 {
        let period = period_of(context.motif);
        let bp_diff = observation.len() as i64 - candidate.bases.len() as i64;

        // How likely a length change of this size is at all. Zero past the cutoffs, and a
        // read that lands there falls to the row's outlier term rather than being explained
        // away as an implausibly large slip.
        let length_probability = context.stutter.probability(bp_diff, period);
        if length_probability <= 0.0 {
            return 0.0;
        }

        let motif = context.motif.as_bytes();
        let period_bases = i64::from(period.get());
        let letters = if bp_diff % period_bases == 0 {
            // A whole-repeat change: every run the slip could have landed in gives its own
            // sequence, and they are averaged with equal weight.
            enumerate_placements(
                candidate.bases,
                motif,
                bp_diff / period_bases,
                &mut scratch.placements,
            );
            if scratch.placements.is_empty() {
                // Unreachable from this candidate — more repeats contracted than any run
                // holds. `SsrScoringContext::unreachable_mass` is what accounts for it.
                return 0.0;
            }
            let each = 1.0 / scratch.placements.len() as f64;
            scratch
                .placements
                .iter()
                .map(|placement| each * substitution_probability(observation, placement, context))
                .sum()
        } else {
            // A part-repeat change: one placement, resized at the tract's end.
            resize_at_the_end(candidate.bases, motif, bp_diff, &mut scratch.resized);
            substitution_probability(observation, &scratch.resized, context)
        };

        length_probability * letters
    }

    fn censored_emission(
        &self,
        _witnessed_prefix: &[u8],
        _candidate: &SsrCandidate<'_>,
        _context: &SsrScoringContext<'_>,
        _scratch: &mut Self::Scratch,
    ) -> f64 {
        // Milestone G's, and deliberately not a placeholder that scores: a read that ran out
        // is not a shorter complete observation, and treating it as one would let a truncated
        // read out-discriminate a whole one (spec §5.2).
        unimplemented!("the censored term is plan step G1")
    }
}

/// `P(the letters differ this way)` for two sequences the stutter factor has already made the
/// same length — [`FlatEmission`] under the context's fitted rate, in probability space.
///
/// **Unequal lengths cannot arrive here**: every branch above resizes first, which is the
/// separation spec §4.3 requires. A mismatch in length would mean the resize was wrong, so it
/// is an assertion rather than a silent zero.
fn substitution_probability(
    observation: &[u8],
    resized_candidate: &[u8],
    context: &SsrScoringContext<'_>,
) -> f64 {
    debug_assert_eq!(
        observation.len(),
        resized_candidate.len(),
        "the stutter factor equalises the lengths before the letters are compared"
    );
    if observation.len() != resized_candidate.len() {
        return 0.0;
    }

    let flat = FlatEmission::try_new(context.substitution_rate.get())
        .expect("an ErrorRate is a probability by construction");
    // A flat model's two scores do not vary with quality, so they are resolved once and the
    // loop is one comparison per base.
    // A flat model ignores the quality entirely, so any value resolves the same two
    // scores; zero is the one that says "this is not read from the base".
    let scores = flat.scores_for(BaseQual(0));
    let total_ln: f64 = observation
        .iter()
        .zip(resized_candidate)
        .map(|(&read_base, &candidate_base)| scores.pick(read_base, candidate_base))
        .sum();
    total_ln.exp()
}

/// Every placement-distinct realisation of `candidate` with `repeats` whole repeats added to
/// one of its runs, written into `out` (cleared first).
///
/// A pure tract is one run, so it yields one sequence. An interrupted tract yields one per run
/// the slip could land in, **deduplicated** — two runs of the same length in a symmetric tract
/// give the same bytes, and counting that twice would weight it double.
///
/// A run that cannot absorb the change contributes nothing, and a candidate where **no** run
/// can yields an empty list, which the caller reads as unreachable. Two limits make a run
/// unable to: it cannot go below zero repeats, and **the tract as a whole must keep at least
/// one** — the same rule [`StutterModel::unreachable_mass`] applies from the distribution's
/// side, so the mass the report calls unplaced is the mass the scoring actually declines to
/// place.
///
/// **On a pure tract the two agree exactly. On an interrupted one they can differ**, and in a
/// stated direction: a slip lands in *one* run, so a tract of two runs of two repeats cannot
/// lose three even though it holds four — while `unreachable_mass` sees only the total and
/// counts that contraction as reachable. **The report therefore understates the loss on
/// interrupted tracts.** Closing that needs the run structure to reach the distribution,
/// which takes only a period and a repeat count today; it is recorded as an open question
/// rather than papered over.
fn enumerate_placements(candidate: &[u8], motif: &[u8], repeats: i64, out: &mut Vec<Vec<u8>>) {
    out.clear();
    if repeats == 0 {
        out.push(candidate.to_vec());
        return;
    }

    let segments = segment_tract(candidate, motif);
    let run_count = segments
        .iter()
        .filter(|segment| matches!(segment, TractSegment::Run(_)))
        .count();

    let held_in_total: usize = segments
        .iter()
        .filter_map(|segment| match segment {
            TractSegment::Run(repeats) => Some(*repeats),
            TractSegment::Interruption(_) => None,
        })
        .sum();

    for target in 0..run_count {
        let held = segments
            .iter()
            .filter_map(|segment| match segment {
                TractSegment::Run(repeats) => Some(*repeats),
                TractSegment::Interruption(_) => None,
            })
            .nth(target)
            .expect("the run count bounds the index");
        let Some(resized) = (held as i64).checked_add(repeats).filter(|held| *held >= 0) else {
            continue;
        };
        // **The tract must keep a repeat**, which is the same rule
        // [`StutterModel::unreachable_mass`] applies from the distribution's side: a tract
        // holding none is not this locus any more. Without it the two halves disagree by one
        // step on a pure tract — the report would call a total contraction unreachable while
        // the scoring placed mass on it.
        if held_in_total + resized as usize - held == 0 {
            continue;
        }
        let rendered = render_tract(&segments, motif, target, resized as usize);
        if !out.contains(&rendered) {
            out.push(rendered);
        }
    }
}

/// One stretch of a tract: a run of whole motif copies, or the bases between two runs.
#[derive(Debug, Clone, PartialEq, Eq)]
enum TractSegment {
    /// How many whole copies of the motif this run holds.
    Run(usize),
    /// Bases that are not a whole copy — an interruption, or a partial copy at either end.
    Interruption(Vec<u8>),
}

/// Split a tract into runs of the motif and the interruptions between them, matching the
/// motif greedily from the left.
///
/// **Left-to-right and greedy, which is a choice about interrupted tracts.** A different
/// phase would split some tracts differently; this is production's rule, and the two must
/// agree while the models are meant to.
fn segment_tract(candidate: &[u8], motif: &[u8]) -> Vec<TractSegment> {
    let period = motif.len();
    let mut segments = Vec::new();
    let mut run = 0usize;
    let mut interruption: Vec<u8> = Vec::new();
    let mut at = 0usize;

    while at < candidate.len() {
        if at + period <= candidate.len() && &candidate[at..at + period] == motif {
            if !interruption.is_empty() {
                segments.push(TractSegment::Interruption(std::mem::take(
                    &mut interruption,
                )));
            }
            run += 1;
            at += period;
        } else {
            if run > 0 {
                segments.push(TractSegment::Run(run));
                run = 0;
            }
            interruption.push(candidate[at]);
            at += 1;
        }
    }
    if run > 0 {
        segments.push(TractSegment::Run(run));
    }
    if !interruption.is_empty() {
        segments.push(TractSegment::Interruption(interruption));
    }
    segments
}

/// Render the segments back to bytes, with the run at `target` holding `resized` repeats
/// instead of its own.
fn render_tract(segments: &[TractSegment], motif: &[u8], target: usize, resized: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let mut seen = 0usize;
    for segment in segments {
        match segment {
            TractSegment::Run(held) => {
                let repeats = if seen == target { resized } else { *held };
                for _ in 0..repeats {
                    out.extend_from_slice(motif);
                }
                seen += 1;
            }
            TractSegment::Interruption(bases) => out.extend_from_slice(bases),
        }
    }
    out
}

/// Resize a candidate by `bp_diff` bases **at the tract's end**, in one placement — the
/// part-repeat branch.
///
/// Lengthening continues the motif's tiling from the phase the tract ended on, so the added
/// bases are the ones a slipped polymerase would have written. Shortening drops bases from the
/// end. Production's rule, and a simplification: a part-repeat change could have landed
/// anywhere, and enumerating those placements is the follow-up spec §4.2 records.
fn resize_at_the_end(candidate: &[u8], motif: &[u8], bp_diff: i64, out: &mut Vec<u8>) {
    out.clear();
    let period = motif.len();
    if bp_diff >= 0 {
        out.extend_from_slice(candidate);
        for step in 0..bp_diff as usize {
            out.push(motif[(candidate.len() + step) % period]);
        }
    } else {
        let keep = candidate
            .len()
            .saturating_sub(bp_diff.unsigned_abs() as usize);
        out.extend_from_slice(&candidate[..keep]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::MarginalAligner;
    use crate::ng::alignment::ssr_marginal_sequence::{
        SequenceMarginalScratch, SsrSequenceMarginal,
    };
    use crate::ng::alignment::stutter::MAX_WHOLE_REPEAT_SLIP;
    use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;

    /// **Model B — the comparator, and an oracle rather than an alternative.**
    ///
    /// `Σ_n S(n) · avg_v align(observation | candidate ⊕ n repeats)`: the whole-repeat
    /// stutter mass marginalised **outside** a sequence-versus-sequence align that absorbs
    /// whatever length is left over as slop at the tract's ends.
    ///
    /// # How it differs from Model A, which is the whole point
    ///
    /// **Model A picks one length change** — the one the observation actually shows — scores
    /// its probability, and compares letters over sequences it has already made equal. **Model
    /// B sums over every whole-repeat change**, and lets the aligner explain the residual.
    /// Where A has an explicit part-repeat branch with its own two parameters, B has none: a
    /// change that is not a whole number of repeats reaches B only as end-gap slop on a
    /// neighbouring term.
    ///
    /// So the two explain a read's length by genuinely different routes, and that is what
    /// makes agreement between them evidence. **What they share is the placement enumeration**
    /// — production shares it too — so the independence is in how length is explained, not in
    /// where a slip may land.
    ///
    /// **Test-only, deliberately.** It is worth more as a second opinion than as a model
    /// anyone runs; production keeps its own Model B the same way (spec §9).
    #[derive(Debug, Default, Clone, Copy)]
    struct ClassicEmissionOracle;

    #[derive(Debug, Default)]
    struct ClassicOracleScratch {
        placements: Vec<Vec<u8>>,
        aligner: SequenceMarginalScratch,
    }

    impl SsrEmissionModel for ClassicEmissionOracle {
        type Scratch = ClassicOracleScratch;

        fn emission(
            &self,
            observation: &[u8],
            candidate: &SsrCandidate<'_>,
            context: &SsrScoringContext<'_>,
            scratch: &mut Self::Scratch,
        ) -> f64 {
            let period = period_of(context.motif);
            let period_bases = i64::from(period.get());
            let motif = context.motif.as_bytes();
            let aligner = SsrSequenceMarginal::try_new(context.substitution_rate.get())
                .expect("an ErrorRate is a probability by construction");

            let cutoff = i64::from(MAX_WHOLE_REPEAT_SLIP);
            let mut total = 0.0;
            for repeats in -cutoff..=cutoff {
                let length_probability =
                    context.stutter.probability(repeats * period_bases, period);
                if length_probability <= 0.0 {
                    continue;
                }
                enumerate_placements(candidate.bases, motif, repeats, &mut scratch.placements);
                if scratch.placements.is_empty() {
                    continue;
                }
                let each = 1.0 / scratch.placements.len() as f64;
                let letters: f64 = scratch
                    .placements
                    .iter()
                    .map(|placement| {
                        each * aligner
                            .marginal_probability(observation, placement, (), &mut scratch.aligner)
                            .get()
                            .exp()
                    })
                    .sum();
                total += length_probability * letters;
            }
            total
        }

        fn censored_emission(
            &self,
            _witnessed_prefix: &[u8],
            _candidate: &SsrCandidate<'_>,
            _context: &SsrScoringContext<'_>,
            _scratch: &mut Self::Scratch,
        ) -> f64 {
            unimplemented!("the oracle scores complete observations only")
        }
    }

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
    /// Score one observation against one candidate under Model A, at a stated substitution
    /// rate — the shape every test below uses.
    fn score(
        observation: &[u8],
        candidate_bases: &[u8],
        repeat_count: u32,
        motif: &[u8],
        model: &StutterModel,
        substitution_rate: f64,
    ) -> f64 {
        let motif = a_motif(motif);
        let candidate = SsrCandidate {
            bases: candidate_bases,
            repeat_count: repeats(repeat_count),
        };
        let context = SsrScoringContext::new(
            &motif,
            model,
            &candidate,
            ErrorRate::try_new(substitution_rate).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        let mut scratch = StutterSubstitutionScratch::default();
        StutterSubstitutionEmission.emission(observation, &candidate, &context, &mut scratch)
    }

    /// Rank the candidates by one model's emission, best first.
    fn ranking<M: SsrEmissionModel>(
        model: &M,
        observation: &[u8],
        candidates: &[(Vec<u8>, u32)],
        motif: &Motif,
        stutter: &StutterModel,
        epsilon: f64,
        scratch: &mut M::Scratch,
    ) -> (Vec<usize>, Vec<f64>) {
        let scores: Vec<f64> = candidates
            .iter()
            .map(|(bases, count)| {
                let candidate = SsrCandidate {
                    bases,
                    repeat_count: repeats(*count),
                };
                let context = SsrScoringContext::new(
                    motif,
                    stutter,
                    &candidate,
                    ErrorRate::try_new(epsilon).expect("a valid rate"),
                    [Provenance::FittedHere],
                );
                model.emission(observation, &candidate, &context, scratch)
            })
            .collect();
        let mut order: Vec<usize> = (0..scores.len()).collect();
        order.sort_by(|a, b| scores[*b].partial_cmp(&scores[*a]).expect("no NaN scores"));
        (order, scores)
    }

    /// **The independent-implementation check: the two models rank candidates alike wherever
    /// the observation's length differs from the candidate's by a whole number of repeats.**
    ///
    /// That is the case Model B is built for, and the two reach it by genuinely different
    /// routes — A picks the single length change the read shows and compares letters over
    /// equal-length sequences; B sums over *every* whole-repeat change and lets the aligner
    /// absorb the residual. Agreement across that grid is evidence about both.
    ///
    /// Measured here: the full ranking is identical for every observation, and the winning
    /// score agrees to better than one part in ten thousand.
    #[test]
    fn the_two_models_rank_candidates_alike_on_whole_repeat_differences() {
        let stutter = a_model();
        let motif = a_motif(b"CAG");
        let epsilon = 0.001;
        let candidates: Vec<(Vec<u8>, u32)> = (3..=6)
            .map(|count| (b"CAG".repeat(count), count as u32))
            .collect();

        // Whole numbers of repeats, plus one carrying a substitution — every one a whole-repeat
        // difference from every candidate.
        let observations: Vec<Vec<u8>> = vec![
            b"CAGCAGCAG".to_vec(),
            b"CAGCAGCAGCAG".to_vec(),
            b"CAGCAGCAGCAGCAG".to_vec(),
            b"CAGCTGCAGCAG".to_vec(),
        ];

        for observation in &observations {
            let mut a_scratch = StutterSubstitutionScratch::default();
            let mut b_scratch = ClassicOracleScratch::default();
            let (a_order, a_scores) = ranking(
                &StutterSubstitutionEmission,
                observation,
                &candidates,
                &motif,
                &stutter,
                epsilon,
                &mut a_scratch,
            );
            let (b_order, b_scores) = ranking(
                &ClassicEmissionOracle,
                observation,
                &candidates,
                &motif,
                &stutter,
                epsilon,
                &mut b_scratch,
            );

            assert_eq!(
                a_order,
                b_order,
                "the two models ranked {} differently: {a_scores:?} against {b_scores:?}",
                String::from_utf8_lossy(observation)
            );

            let winner = a_order[0];
            let relative = (a_scores[winner] - b_scores[winner]).abs() / a_scores[winner];
            assert!(
                relative < 1e-4,
                "the winning scores differ by {relative} on {}",
                String::from_utf8_lossy(observation)
            );
        }
    }

    /// **Where the two part company, and why that is the choice rather than a defect.**
    ///
    /// A read whose length is *not* a whole number of repeats from the candidate is the one
    /// case the two models explain differently by construction. **Model A charges the fitted
    /// part-repeat share**; **Model B has no such branch at all** and absorbs the odd base as
    /// slop at the tract's end, at the flat per-base rate.
    ///
    /// At ε = 0.001 the slop route is about nine times cheaper than the part-repeat event, so
    /// **B keeps the shorter candidate and A prefers the one needing the smaller stutter
    /// event**. Measured on a ten-base read against a three-repeat and a four-repeat
    /// candidate: A scores 1.094e-4 and 1.869e-4 — the four-repeat candidate wins — while B
    /// scores 9.706e-4 and 1.167e-5, and the three-repeat candidate wins.
    ///
    /// **Model A is the one to trust here**, and this divergence is the reason: an explicit
    /// part-repeat branch with its own two fitted parameters is what the comparison behind
    /// spec §4.1 chose it for. The test exists so the difference stays a measured, deliberate
    /// thing rather than a surprise the first time a ranking disagrees.
    #[test]
    fn the_two_models_part_company_on_a_part_repeat_length() {
        let stutter = a_model();
        let motif = a_motif(b"CAG");
        let epsilon = 0.001;
        let candidates: Vec<(Vec<u8>, u32)> =
            vec![(b"CAGCAGCAG".to_vec(), 3), (b"CAGCAGCAGCAG".to_vec(), 4)];
        // Three repeats and one base over — a part-repeat difference from both candidates.
        let observation = b"CAGCAGCAGC";

        let mut a_scratch = StutterSubstitutionScratch::default();
        let mut b_scratch = ClassicOracleScratch::default();
        let (a_order, a_scores) = ranking(
            &StutterSubstitutionEmission,
            observation,
            &candidates,
            &motif,
            &stutter,
            epsilon,
            &mut a_scratch,
        );
        let (b_order, b_scores) = ranking(
            &ClassicEmissionOracle,
            observation,
            &candidates,
            &motif,
            &stutter,
            epsilon,
            &mut b_scratch,
        );

        assert_eq!(
            a_order[0], 1,
            "Model A should prefer the four-repeat candidate"
        );
        assert_eq!(
            b_order[0], 0,
            "Model B should prefer the three-repeat candidate"
        );

        // The sizes, so the mechanism is checkable and not merely asserted.
        assert!((a_scores[0] - 1.094e-4).abs() < 1e-7, "{:?}", a_scores);
        assert!((a_scores[1] - 1.869e-4).abs() < 1e-7, "{:?}", a_scores);
        assert!((b_scores[0] - 9.706e-4).abs() < 1e-6, "{:?}", b_scores);
        assert!((b_scores[1] - 1.167e-5).abs() < 1e-7, "{:?}", b_scores);

        // B's route is the flat per-base rate; A's is the fitted part-repeat product.
        let slop_route = b_scores[0] / stutter.same_length_share();
        let part_repeat_route =
            stutter.part_repeat_longer_share() * stutter.part_repeat_one_step_share();
        assert!(
            slop_route > part_repeat_route * 8.0,
            "slop {slop_route} against a part-repeat event {part_repeat_route}"
        );
    }

    /// **Spec §12's first test: a read matching its allele is dominated by the same-length
    /// term.** A read identical to the candidate scores at least
    /// `same_length_share × (1 − ε)^length` and within 5% of it — the two factors and nothing
    /// else, with no placement multiplicity to dilute it.
    #[test]
    fn a_read_identical_to_its_candidate_is_dominated_by_the_same_length_term() {
        let model = a_model();
        let tract = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;
        let scored = score(tract, tract, 6, b"CAG", &model, epsilon);

        let floor = model.same_length_share() * (1.0 - epsilon).powi(tract.len() as i32);
        assert!(scored >= floor, "{scored} fell below {floor}");
        assert!(
            scored <= floor * 1.05,
            "{scored} exceeded {floor} by more than a twentieth"
        );
    }

    /// **Spec §12's second test: direction and size are ordered as fitted.** With a
    /// contraction-biased split, a read one repeat short outscores one a repeat long; and one
    /// repeat outscores two. Both are properties of the fitted numbers, so a model that
    /// dropped the direction split or the size decay would fail here rather than silently
    /// score a symmetric distribution.
    #[test]
    fn a_shorter_read_outscores_a_longer_one_and_one_repeat_outscores_two() {
        let model = a_model();
        let candidate = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;

        let one_short = score(b"CAGCAGCAGCAGCAG", candidate, 6, b"CAG", &model, epsilon);
        let one_long = score(
            b"CAGCAGCAGCAGCAGCAGCAG",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        assert!(
            one_short > one_long,
            "a repeat short scored {one_short}, a repeat long {one_long}"
        );

        let two_short = score(b"CAGCAGCAGCAG", candidate, 6, b"CAG", &model, epsilon);
        assert!(
            one_short > two_short,
            "one repeat short scored {one_short}, two short {two_short}"
        );
    }

    /// **Spec §12's third test: a whole repeat beats a single stray base** — and under the
    /// **corrected** condition, which is a comparison of *products* rather than of the
    /// direction shares alone.
    ///
    /// "Any part-repeat share below the whole-repeat one" suffices only while the two one-step
    /// shares are tied, and spec §10 schedules untying them, so the test states the condition
    /// it actually depends on and asserts it holds for this fixture before asserting the
    /// consequence. The specification's own counter-example: a whole-repeat share of 0.02 at a
    /// one-step share of 0.1 gives 0.002, against a part-repeat share of 0.019 at 0.95 giving
    /// 0.018 — nine times higher, both inside the clamps.
    #[test]
    fn a_whole_repeat_beats_a_stray_base_when_the_products_say_so() {
        let model = a_model();
        let candidate = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;

        let whole_product = model.whole_repeat_longer_share() * model.whole_repeat_one_step_share();
        let part_product = model.part_repeat_longer_share() * model.part_repeat_one_step_share();
        assert!(
            part_product < whole_product,
            "this fixture does not meet the condition: {part_product} against {whole_product}"
        );

        let whole_repeat_longer = score(
            b"CAGCAGCAGCAGCAGCAGCAG",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        let one_base_longer = score(
            b"CAGCAGCAGCAGCAGCAGC",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        assert!(
            whole_repeat_longer > one_base_longer,
            "a whole repeat scored {whole_repeat_longer}, a stray base {one_base_longer}"
        );
    }

    /// **An interrupted tract has more than one place a slip could land, and they are averaged
    /// with equal weight.** Two runs either side of an interruption give two distinct
    /// sequences for a one-repeat expansion; a read matching one of them is scored at half the
    /// letters' probability, because the model does not know which run slipped.
    #[test]
    fn an_interrupted_tract_averages_over_the_placements_a_slip_could_land_in() {
        let mut placements = Vec::new();
        enumerate_placements(b"CACACATTCACA", b"CA", 1, &mut placements);
        assert_eq!(placements.len(), 2, "{placements:?}");
        assert!(placements.contains(&b"CACACACATTCACA".to_vec()));
        assert!(placements.contains(&b"CACACATTCACACA".to_vec()));

        // A pure tract has one placement, so nothing is diluted there.
        let mut pure = Vec::new();
        enumerate_placements(b"CACACACACA", b"CA", 1, &mut pure);
        assert_eq!(pure.len(), 1);
        assert_eq!(pure[0], b"CACACACACACA".to_vec());
    }

    /// **A contraction no single run can absorb is unreachable, and scores exactly zero** —
    /// not a small number. The mass that costs is what `SsrScoringContext::unreachable_mass`
    /// reports, which is why it is reported rather than assumed away.
    ///
    /// **It takes an interrupted tract to reach this branch through `emission` at all**, and
    /// that is worth stating: a slip lands in *one* run, so a tract of two runs of two
    /// repeats cannot lose three even though it holds four. On a **pure** tract the branch is
    /// unreachable from `emission`, because the deepest contraction an observation can ask
    /// for is the one that leaves it empty — and a single run can always give up everything
    /// it holds.
    #[test]
    fn a_contraction_no_single_run_can_absorb_scores_zero() {
        let model = a_model();
        // Two runs of two repeats, four in total: no run can give up three.
        let interrupted = b"CACATTCACA";
        let scored = score(b"CACA", interrupted, 4, b"CA", &model, 0.001);
        assert_eq!(scored, 0.0);

        let mut placements = Vec::new();
        enumerate_placements(interrupted, b"CA", -3, &mut placements);
        assert!(placements.is_empty(), "{placements:?}");

        // Two repeats it can give up, from either run — and emptying one run is not the same
        // sequence as emptying the other, so both placements stand: `TTCACA` and `CACATT`.
        enumerate_placements(interrupted, b"CA", -2, &mut placements);
        assert_eq!(placements.len(), 2, "{placements:?}");
        assert!(placements.contains(&b"TTCACA".to_vec()));
        assert!(placements.contains(&b"CACATT".to_vec()));

        // **A pure tract may not vanish entirely**, which is the rule `unreachable_mass`
        // states from the distribution's side. Four repeats can give up three, not four.
        enumerate_placements(b"CACACACA", b"CA", -3, &mut placements);
        assert_eq!(placements.len(), 1, "{placements:?}");
        enumerate_placements(b"CACACACA", b"CA", -4, &mut placements);
        assert!(placements.is_empty(), "{placements:?}");
    }

    /// **A part-repeat change is resized at the tract's end, in one placement** — production's
    /// split, and a simplification the specification records rather than hides. Lengthening
    /// continues the motif's tiling from the phase the tract ended on.
    #[test]
    fn a_part_repeat_change_is_resized_at_the_tract_end() {
        let mut resized = Vec::new();

        resize_at_the_end(b"CAGCAGCAG", b"CAG", 1, &mut resized);
        assert_eq!(resized, b"CAGCAGCAGC".to_vec());

        resize_at_the_end(b"CAGCAGCAG", b"CAG", 2, &mut resized);
        assert_eq!(resized, b"CAGCAGCAGCA".to_vec());

        resize_at_the_end(b"CAGCAGCAG", b"CAG", -1, &mut resized);
        assert_eq!(resized, b"CAGCAGCA".to_vec());
    }

    /// **The two factors are separable, and the substitution one is the alignment module's.**
    /// At a candidate the read matches exactly, the score is the same-length share times the
    /// flat model's own product over the tract — so a change to either factor moves the score
    /// by exactly that factor and nothing else.
    #[test]
    fn the_score_is_the_length_factor_times_the_letter_factor() {
        let model = a_model();
        let tract = b"CAGCAGCAGCAGCAGCAG";

        for epsilon in [1e-4, 1e-3, 1e-2] {
            let scored = score(tract, tract, 6, b"CAG", &model, epsilon);
            let letters = (1.0 - epsilon).powi(tract.len() as i32);
            let expected = model.same_length_share() * letters;
            assert!(
                (scored - expected).abs() <= expected * 1e-12,
                "at ε={epsilon}: {scored} against {expected}"
            );
        }

        // One mismatching base costs exactly one `ε/3` in place of one `1 − ε`.
        let epsilon = 1e-3;
        let mismatched = b"TAGCAGCAGCAGCAGCAG";
        let scored = score(mismatched, tract, 6, b"CAG", &model, epsilon);
        let expected = model.same_length_share()
            * (1.0 - epsilon).powi(tract.len() as i32 - 1)
            * (epsilon / 3.0);
        assert!(
            (scored - expected).abs() <= expected * 1e-12,
            "{scored} against {expected}"
        );
    }
}
