//! Maximising the noise parameters when a ladder is the wrong instrument, and reporting
//! enough beside the answer to tell a fit from a stopped search.
//!
//! **The second of `fitting/`'s two search drivers, and it exists because the first one
//! does not fit the STR path.** [`ladder_scan::fit_by_profile_scan`](super::ladder_scan)
//! steps a prepared list of candidates end to end, which is right for one error rate on a
//! quarter-Phred ladder. Here there are three parameters at once: a flat scan over them is
//! 4.2 million scores **per (read group × stratum)** with several hundred strata a sample,
//! and a quarter-Phred spacing is the wrong ruler for two of the three, which are shares in
//! `(0, 1)` rather than rates spanning orders of magnitude
//! (`spec/parameter_prepass_ssr.md` §4.2).
//!
//! So this searches continuously from several starting points and keeps the best-scoring —
//! **with every start's outcome beside it**, which is not decoration:
//!
//! - **Starts that disagree are what says the answer was found rather than assumed.** The
//!   sibling path's inbreeding fit returned a confident zero on a genome 29% covered by
//!   runs, from five starts that disagreed about the headline number while sharing one
//!   guess at a nuisance axis. Every start here must differ on **every** axis.
//! - **And starts that agree are not evidence either**, which is the trap that cost a
//!   published finding on 2026-08-06: a deterministic search returns the same point from
//!   every start wherever the objective is flat, so four starts agreeing to four decimal
//!   places is exactly what a search that never looked also produces. The spread is
//!   reported to be paired with a score at a known point, never read alone.
//!
//! **Convergence failure is a data condition here, unlike the inner climb.** The climb over
//! the genotype frequencies is provably concave, so a stall there is a bug; this search has
//! no such proof, so it is capped, the best-scoring point is kept, and how it ended is
//! reported (`arch/parameter_prepass_ssr.md` §3).

use std::collections::BTreeMap;

use smallvec::SmallVec;

use super::ladder_scan::{PloidyPlan, ploidy_plans, score_one_candidate};
use super::mixture_weights::climb_mixture_weights;
use super::{FitTermination, NoiseModel};
use crate::ng::types::{LogProb, Ploidy};

/// **A set of noise parameters a coordinate search can walk**, and the scale each of its
/// axes is walked on.
///
/// **The scale is the point of this trait, not the reading and writing.** A slippage level
/// runs from about 0.0009 to 0.15 and is searched on a log scale, where a step is a fixed
/// *fraction* of the rate; a direction share is a fraction in `(0, 1)` and is searched on a
/// logit scale, where the search cannot walk out of the range whatever step it takes. A
/// search that stepped both linearly would resolve the level to nothing at the bottom of its
/// range and would have to be told to stop at 0 and 1.
///
/// **Axis 0 is the headline parameter** — the one [`MultistartResult::headline_spread`]
/// reports across the starts, and the one a caller judges "did this fit find something" on.
pub trait SearchableNoise: Sized {
    /// How many axes the search sweeps.
    fn axes() -> usize;

    /// The lowest and highest value of one axis, **on that axis's own scale**.
    fn axis_bounds(axis: usize) -> (f64, f64);

    /// One axis's current value, on that axis's own scale.
    fn coordinate(&self, axis: usize) -> f64;

    /// The same parameters with one axis moved to `value`, on that axis's own scale.
    fn with_coordinate(&self, axis: usize, value: f64) -> Self;

    /// **The headline parameter on its natural scale** — the quantity
    /// [`MultistartResult::headline_spread`] is a ratio of.
    ///
    /// **Asked for rather than derived from axis 0's coordinate**, and that is not tidiness: a
    /// coordinate is on whatever scale its axis is searched on, so taking `exp` of a difference
    /// of coordinates is the ratio only if that axis happens to be logarithmic. On a logit axis
    /// the same arithmetic returns the *odds* ratio — measured on this driver's own fixture,
    /// 12.17 where the ratio is 8.64, a 41% error reported as a spread. Since the spread is the
    /// one number a caller gates a fit on, the scale is asked for rather than assumed.
    fn headline(&self) -> f64;
}

/// How finely each axis is resolved and how many sweeps the search takes.
///
/// **These are resolution knobs, not accuracy ones, and they are set by what a consumer can
/// feel.** A tolerance of 0.01 on a rate's log scale is 1% of the rate, and
/// `spec/parameter_prepass.md` §3 argues a few percent between neighbouring values is
/// already finer than a genotype likelihood registers — shifting a prior by that much moves
/// a call far less than one more read would. An earlier version of this search resolved to
/// 1e-7, seven digits on a quantity known to two, and spent about four times the work.
#[derive(Copy, Clone, Debug)]
pub struct SearchPrecision {
    /// How finely an axis is resolved, on its own transformed scale. On a log scale this is
    /// a relative resolution: 0.01 is 1% of the value.
    pub tolerance: f64,
    /// How many golden-section steps one axis gets in one sweep.
    pub max_axis_steps: usize,
    /// How many times the search cycles through the axes.
    pub max_sweeps: usize,
}

impl SearchPrecision {
    /// **For a real table — resolution matched to what a consumer can feel.** 1% of the
    /// headline rate, for the reason on the type.
    #[must_use]
    pub fn fast() -> Self {
        Self {
            tolerance: 0.01,
            max_axis_steps: 16,
            max_sweeps: 5,
        }
    }

    /// **For a control, where the point is sharpness rather than speed.** A control's job is
    /// to say the estimator is unbiased, and it can only say so to the resolution the search
    /// itself has — at [`Self::fast`] a genuine 1% bias and the search's own step are the
    /// same size, so the control would be reporting its own grid.
    #[must_use]
    pub fn fine() -> Self {
        Self {
            tolerance: 1e-5,
            max_axis_steps: 40,
            max_sweeps: 8,
        }
    }
}

/// Where one starting point began, where it ended, and what it scored.
///
/// **Emitted for every start and not only the winner.** An answer with no spread beside it
/// is indistinguishable from a search that never looked.
#[derive(Clone, PartialEq, Debug)]
pub struct StartOutcome<P> {
    /// Where this start began.
    ///
    /// Named `from` rather than `start`, so that a reader of `start.start` is never asked to
    /// work out which of the two is the noun.
    pub from: P,
    /// Where the search from `from` ended.
    pub reached: P,
    /// What `reached` scored.
    pub log_likelihood: LogProb,
}

/// What a multi-start search returned: the winner, and enough beside it to tell a fitted
/// number from a stopped search.
#[derive(Clone, PartialEq, Debug)]
pub struct MultistartResult<P> {
    /// The best-scoring parameters any start reached.
    pub best: P,
    /// The genotype frequencies climbed at [`Self::best`], one set per ploidy the cells
    /// cover.
    pub genotype_frequencies: BTreeMap<Ploidy, Vec<f64>>,
    /// What [`Self::best`] scored.
    pub log_likelihood: LogProb,
    /// Every start's outcome, **best-scoring first**. Never empty.
    pub starts: SmallVec<[StartOutcome<P>; 4]>,
    /// The ratio between the highest and lowest value any start reached on **axis 0**, the
    /// headline parameter.
    ///
    /// **A ratio and not a difference**, because the quantity it describes spans orders of
    /// magnitude: a slippage level runs from 0.0009 to 0.15 across the repeat counts of one
    /// dataset, so a fixed difference would be meaningless at one end of that and enormous
    /// at the other.
    ///
    /// **The caller compares it against its own limit rather than this driver doing so**,
    /// because how far apart two answers may sit before the fit is disowned is a
    /// path-specific judgement — the STR path's `START_AGREEMENT_LIMIT` is one quarter-Phred,
    /// borrowed from the sibling path's ladder spacing.
    ///
    /// Infinite where some start reached zero or below, which is a real answer for a slippage
    /// level and not a degenerate one.
    ///
    /// **Read it as a floor on the disagreement, not as four independent searches agreeing.**
    /// Each axis is line-searched over its *whole* range at every sweep, so a start's own value
    /// on the headline axis is overwritten before it is ever read — measured on this driver's
    /// control, three of four starts reach a bit-identical point and the spread is exactly 1.00.
    /// That is the search working, and it is also precisely the reading this module's own
    /// preamble warns against, so **a spread near one is evidence of nothing on its own** and has
    /// to be paired with a score at a point whose answer is known.
    pub headline_spread: f64,
    /// How the **worst** of the starts ended — capped or settled — so that a caller is never
    /// handed a number that ran out of sweeps with nothing to say so.
    ///
    /// **The outer search has no concavity proof**, unlike the climb inside it, so running
    /// out is a data condition rather than a bug and is reported instead of asserted.
    pub termination: FitTermination,
    /// Whether **every** inner climb over the genotype frequencies settled, or some ran out of
    /// passes.
    ///
    /// **Reported rather than asserted, and that is a departure from the sibling scan**, which
    /// treats a stalled climb as a debug-time bug. The surface is concave either way, so a climb
    /// that ran out did not find a wrong summit — it ran out of time on the right one. What
    /// differs is how far there is to go: the sibling path's climb has three genotypes and the
    /// pass cap was measured on it, while a stratum here has between 66 and 91, most of them
    /// heading to a frequency of zero, and reaching stillness on those takes far longer. A
    /// candidate scored short of its own summit is scored low, so a neighbour can win on that
    /// alone — which is why this travels with the answer rather than being swallowed.
    pub every_climb_settled: bool,
}

/// Everything one search carries between candidates: the buffers it reuses and the two facts
/// it accumulates.
///
/// **One value rather than three arguments threaded through every frame**, because all three are
/// the same thing — state that outlives a single candidate — and a function taking nine arguments
/// is one whose next reader has to check the call sites to see which is which.
struct SearchState {
    /// The row-major likelihood table, cleared and refilled once per (candidate, ploidy) rather
    /// than allocated there: a search runs several hundred candidates and this is the largest
    /// thing it builds.
    scratch: Vec<f64>,
    /// Whether every inner climb so far settled rather than running out of passes.
    every_climb_settled: bool,
}

/// **Maximise the noise parameters from several starting points**, climbing the genotype
/// frequencies at every trial, and return the best-scoring with every start's outcome beside
/// it.
///
/// Each start is followed by a coordinate-wise golden-section search on the scale each axis
/// lives on ([`SearchableNoise`]), so the search cannot walk out of a share's range and
/// resolves a rate to a fixed fraction of itself rather than to a fixed amount.
///
/// **Continuous rather than on a ladder**, which is the whole difference from
/// [`fit_by_profile_scan`](super::ladder_scan::fit_by_profile_scan): what is being measured
/// here is a property of the objective, and a ladder would report its own spacing.
///
/// # Panics
///
/// If `cells` or `starts` is empty, if a ploidy present carries no sites, if the model
/// declares no genotypes at a ploidy present, or if it writes something that is not a
/// log-likelihood. The message names the start and the sweep, which is the part no later
/// frame can recover.
#[must_use]
pub fn fit_by_multistart<M>(
    model: &M,
    cells: &[M::Cell],
    starts: &[M::NoiseParams],
    precision: SearchPrecision,
) -> MultistartResult<M::NoiseParams>
where
    M: NoiseModel,
    M::NoiseParams: SearchableNoise + Clone,
{
    assert!(
        !starts.is_empty(),
        "a multi-start search needs at least one starting point"
    );
    let plans = ploidy_plans(model, cells);
    let mut state = SearchState {
        scratch: Vec::new(),
        every_climb_settled: true,
    };

    let mut outcomes: SmallVec<[StartOutcome<M::NoiseParams>; 4]> = SmallVec::new();
    let mut sweeps_used = 0;
    let mut every_start_settled = true;

    for (start, from) in starts.iter().enumerate() {
        let (reached, log_likelihood, settled, sweeps) =
            climb_from(model, cells, &plans, &mut state, from, precision, start);
        sweeps_used = sweeps_used.max(sweeps);
        every_start_settled &= settled;
        outcomes.push(StartOutcome {
            from: from.clone(),
            reached,
            log_likelihood: LogProb(log_likelihood),
        });
    }

    let headline_spread = spread_across(&outcomes);

    // Best-scoring first. A total order rather than `partial_cmp().unwrap()`: a score is
    // never `NaN` — the scoring pass refuses one — but `total_cmp` says so without a panic
    // path, and it ranks `NaN` above every finite score, so a sort that silently promoted
    // one would be visible as an absurd winner rather than as a crash somewhere else.
    outcomes.sort_by(|left, right| {
        right
            .log_likelihood
            .get()
            .total_cmp(&left.log_likelihood.get())
    });

    let best = outcomes[0].reached.clone();
    let (log_likelihood, genotype_frequencies) =
        score_and_climb(model, cells, &plans, &mut state, &best, "the winner");

    MultistartResult {
        best,
        genotype_frequencies,
        log_likelihood: LogProb(log_likelihood),
        starts: outcomes,
        headline_spread,
        termination: FitTermination {
            iterations: u32::try_from(sweeps_used).unwrap_or(u32::MAX),
            converged: every_start_settled,
        },
        every_climb_settled: state.every_climb_settled,
    }
}

/// The highest-to-lowest ratio of the headline parameter across the starts' answers.
///
/// **Infinite where the lowest is zero or below**, which is a real answer for a slippage level
/// and not a degenerate one — and the guard is what stops that case returning `NaN`, which
/// compares `false` against any limit and so would be read as agreement.
fn spread_across<P: SearchableNoise>(outcomes: &[StartOutcome<P>]) -> f64 {
    let mut highest = f64::NEG_INFINITY;
    let mut lowest = f64::INFINITY;
    for outcome in outcomes {
        let reached = outcome.reached.headline();
        highest = highest.max(reached);
        lowest = lowest.min(reached);
    }
    if lowest > 0.0 {
        highest / lowest
    } else {
        f64::INFINITY
    }
}

/// One start's climb: coordinate-wise golden-section, sweeping the axes until nothing moves
/// or the sweeps run out. Returns where it stopped, what that scored, whether it settled,
/// and how many sweeps it took.
fn climb_from<M>(
    model: &M,
    cells: &[M::Cell],
    plans: &[PloidyPlan],
    state: &mut SearchState,
    from: &M::NoiseParams,
    precision: SearchPrecision,
    start: usize,
) -> (M::NoiseParams, f64, bool, usize)
where
    M: NoiseModel,
    M::NoiseParams: SearchableNoise + Clone,
{
    // `1/φ`. Golden-section keeps one of the two interior evaluations at every step, so a
    // bracket costs one score per step rather than two.
    let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
    let mut current = from.clone();
    let mut settled = false;
    let mut sweeps = 0;
    let mut best = from.clone();
    let mut best_score = f64::NEG_INFINITY;

    for sweep in 0..precision.max_sweeps {
        sweeps = sweep + 1;
        let mut furthest_move: f64 = 0.0;

        for axis in 0..M::NoiseParams::axes() {
            let (mut low, mut high) = M::NoiseParams::axis_bounds(axis);
            let mut left = high - inverse_phi * (high - low);
            let mut right = low + inverse_phi * (high - low);
            let score_at = |value: f64, base: &M::NoiseParams, state: &mut SearchState| {
                let trial = base.with_coordinate(axis, value);
                score_and_climb(
                    model,
                    cells,
                    plans,
                    state,
                    &trial,
                    &format!("start {start}, sweep {sweep}, axis {axis}"),
                )
                .0
            };
            let mut at_left = score_at(left, &current, state);
            let mut at_right = score_at(right, &current, state);

            for _ in 0..precision.max_axis_steps {
                if at_left > at_right {
                    high = right;
                    right = left;
                    at_right = at_left;
                    left = high - inverse_phi * (high - low);
                    at_left = score_at(left, &current, state);
                } else {
                    low = left;
                    left = right;
                    at_left = at_right;
                    right = low + inverse_phi * (high - low);
                    at_right = score_at(right, &current, state);
                }
                if (high - low) < precision.tolerance {
                    break;
                }
            }

            let before = current.coordinate(axis);
            let after = 0.5 * (low + high);
            furthest_move = furthest_move.max((after - before).abs());
            current = current.with_coordinate(axis, after);

            // **The best point evaluated is kept, not merely the last one reached**, which is
            // what `arch/parameter_prepass_ssr.md` §3 asks of a capped search. The bracket's
            // midpoint is where the line search stopped; on an axis that is not unimodal a sweep
            // can move it downhill, and a search that reported that would be handing back a point
            // it had already seen beaten.
            let here = score_at(after, &current, state);
            if here > best_score {
                best_score = here;
                best = current.clone();
            }
        }

        if furthest_move < precision.tolerance {
            settled = true;
            break;
        }
    }

    if best_score == f64::NEG_INFINITY {
        // No sweep ran, so nothing has been scored yet: the start is the best point there is.
        best_score = score_and_climb(
            model,
            cells,
            plans,
            state,
            from,
            &format!("start {start}, unsearched"),
        )
        .0;
    }
    (best, best_score, settled, sweeps)
}

/// One candidate's score, with the genotype frequencies climbed at it.
///
/// **The climb is what makes this a marginal likelihood rather than a conditional one**: the
/// genotype is summed over, so the frequencies have to be at their own maximum before two
/// candidates can be compared. That is the alternation `spec/parameter_prepass.md` §3.1
/// settles, and the half of it with a concavity proof behind it.
fn score_and_climb<M>(
    model: &M,
    cells: &[M::Cell],
    plans: &[PloidyPlan],
    state: &mut SearchState,
    noise: &M::NoiseParams,
    at: &str,
) -> (f64, BTreeMap<Ploidy, Vec<f64>>)
where
    M: NoiseModel,
{
    let SearchState {
        scratch,
        every_climb_settled,
    } = state;
    score_one_candidate(
        model,
        cells,
        noise,
        plans,
        scratch,
        at,
        |plan, table, climbed_by_ploidy: &mut BTreeMap<Ploidy, Vec<f64>>| {
            // **Every candidate climbs from the uniform point, and warm-starting from the
            // previous candidate's answer was tried and removed.** It is the reference
            // implementation's speed device and it cannot work through this seam: the climb
            // requires a strictly interior start, and after a few hundred passes a stratum's
            // unused allele lengths — 30 of 36 on the control's own fixture — have underflowed to
            // exactly zero, so the previous answer is no longer a legal start. Measured, the
            // filter that would have caught that rejected **668 of 668** climbs. What is left of
            // the idea needs a floor inside the climb rather than a filter outside it.
            let climbed = climb_mixture_weights(table, &plan.cell_weights, &plan.uniform_start);
            *every_climb_settled &= climbed.converged;
            climbed_by_ploidy.insert(plan.ploidy, climbed.genotype_frequencies);
            climbed.log_likelihood.get()
        },
    )
}
