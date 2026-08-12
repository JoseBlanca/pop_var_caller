//! The scan: step through a ladder of noise parameters, score the cells at every rung,
//! and keep the rung that scores highest.
//!
//! **Every rung is scored and no early exit is taken**, because nobody has shown the
//! curve has a single hump. That is the whole reason for a scan rather than a
//! one-dimensional optimiser, and it is why the cost is stated rather than avoided: the
//! generic ladder is 161 rungs, once per read group.
//!
//! **Two scans, differing only in what a rung's score means.**
//!
//! - [`fit_by_profile_scan`] climbs to the best genotype frequencies at every rung and
//!   scores the rung there. In the standard vocabulary that is a **profile likelihood**
//!   over the noise parameters — the frequencies are maximised out at every value of the
//!   parameter being scanned, leaving a curve in that parameter alone. Splitting the
//!   search this way puts the effort where the difficulty is: the frequencies are
//!   provably concave and never needed searching, while about the noise parameters there
//!   is no proof either way (`spec/parameter_prepass.md` §3.1).
//! - [`fit_by_fixed_frequency_scan`] scores every rung at genotype frequencies the caller
//!   hands it and climbs nothing. That is **not** a profile likelihood, which is why it is
//!   a sibling rather than a mode of the first (owner's call, 2026-08-07): a function
//!   whose name says *profile* and which sometimes is not one is a name that has stopped
//!   being true. It is the `ε` half of Milestone E's coupled alternation — the half the
//!   research harness measured, whose error-rate step holds the frequencies where the
//!   previous iteration left them (`spec/parameter_prepass_generic.md` §5.1,
//!   `examples/ng_multilib_key_harness.rs`, `fit_eps_on_read_group`).
//!
//! Everything else about a rung is the same in both, and shared in one place
//! ([`scan_ladder`]): the per-ploidy plan, the check that the model appended the width it
//! declared, the value checks that name the rung, the tie rule, and the rail flag.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4.2.

use std::collections::BTreeMap;

use crate::ng::parameter_estimation::fitting::mixture_weights::{
    GenotypeLikelihoodTable, climb_mixture_weights, weighted_log_likelihood,
};
use crate::ng::parameter_estimation::fitting::{NoiseModel, WeightedCell};
use crate::ng::types::{LogProb, Ploidy};

/// How far the genotype frequencies handed to [`fit_by_fixed_frequency_scan`] may miss
/// summing to one before they are refused.
///
/// They arrive from a climb over the simplex, which normalises by a division, so exact
/// equality is not available; this is loose enough for that and tight enough that a set
/// built from the wrong denominator — or one genotype short — cannot pass.
const FREQUENCY_SUM_TOLERANCE: f64 = 1e-9;

/// What one scan over a ladder of noise parameters returned.
///
/// Generic over the noise parameters, because the two paths scan different things: one
/// error rate here, three stutter parameters on the STR path.
#[derive(Clone, PartialEq, Debug)]
pub struct ScanResult<P> {
    /// The winning rung.
    pub noise: P,
    /// The genotype frequencies climbed to at that rung, **one set per ploidy the cells
    /// covered**. On the error-rate scan these are a means rather than an output — the
    /// scan is run for `noise` and they are discarded — while the sample's own rates
    /// come from a scan run for these.
    ///
    /// **Keyed by ploidy rather than the single vector the architecture sketches**
    /// (§5.2). A haploid region has two genotype classes and a diploid three,
    /// so they cannot share a weight vector and the scan climbs once per ploidy (§4.2);
    /// a single vector would mean picking one of them to report and dropping the rest,
    /// silently, in the module whose whole difficulty is that its wrong numbers have no
    /// symptom.
    ///
    /// The same quantity `MixtureWeightsFit::genotype_frequencies` carries and under the
    /// same name, because it is the same numbers: the climb returns them and the scan
    /// files them by ploidy.
    pub genotype_frequencies: BTreeMap<Ploidy, Vec<f64>>,
    /// What makes "the best-scoring iterate" a defined comparison in an alternating
    /// fit. A [`LogProb`] rather than a bare `f64` because comparing is the only thing
    /// it is for, and `LogProb` carries `ln(0)` as `-∞` — the score of a rung nothing
    /// could have produced — where a linear probability would reach `0` and be
    /// indistinguishable from a value that merely got too small to represent.
    pub log_likelihood: LogProb,
    /// **Whether the answer sat on the ladder's edge, and it is not decoration.** A
    /// read group whose true rate lies outside the ladder — a bad run, heavy
    /// contamination, or any of the arithmetic failures a scan can suffer — has its
    /// answer silently clamped to an endpoint and emitted as though it were fitted,
    /// with a large observation count behind it. This one bit is what stands between a
    /// railed fit and a plausible-looking number.
    pub argmax_at_ladder_end: bool,
}

/// Step through `ladder`; at each rung climb to the genotype frequencies that best
/// explain `cells`, and return the best-scoring rung with the frequencies found there.
///
/// **Takes a slice of cells, not a histogram and not an iterator.** Not a histogram,
/// because the generic path's concrete table is not a shape the STR accumulator has —
/// that would make the "shared with the STR path" claim false. Not an iterator, because
/// the scan re-walks the cells once per rung, and a slice is re-walkable without asking
/// the caller for `Clone` on a map traversal or re-deriving the keys 161 times. The
/// caller materialises once and lends.
///
/// **Ploidy travels with each cell rather than being one argument**, because one noise
/// parameter is fitted across every ploidy its reads covered — a haploid sex chromosome
/// and the diploid autosomes were prepared by the same chemistry. Each cell is scored
/// against its own genotype set, and the frequencies are climbed **once per ploidy** on
/// that ploidy's cells, because a haploid cell has two genotypes to mix and a diploid
/// three. A rung's score is the sum over ploidies.
///
/// **Ties resolve to the last of the tied rungs**, which is the only rule a scan can
/// state without knowing which way its ladder runs. The generic ladder ascends in Phred
/// and so **descends** in error rate (`generic::error_rate_ladder`), so the last of a
/// tied run is the lowest error rate — which is what
/// `impl_plan/parameter_prepass_generic.md` D3 asks for, and
/// `a_tie_resolves_to_the_lower_error_rate_on_the_generic_ladder` is what binds the two
/// together, so that reversing the ladder cannot quietly reverse the rule.
///
/// **It reports whether its answer sat on the ladder's edge, and that is not
/// decoration.** A read group whose true rate lies outside the ladder — a bad run, heavy
/// contamination, or any of the arithmetic failures a scan can suffer — has its answer
/// silently clamped to an endpoint and emitted as though it were fitted, with a large
/// observation count behind it. That one bit is the only thing standing between a railed
/// fit and a plausible-looking number.
///
/// # Panics
///
/// If `cells` or `ladder` is empty, if any ploidy present carries no sites at all, if the
/// model declares no genotypes at a ploidy present, if it appends a number of entries
/// other than the [`NoiseModel::genotypes`] it declares, or if what it wrote is not a
/// log-likelihood — a `NaN`, a `+∞`, or a whole row of `−∞` on a cell that carries sites.
/// **The last three name the rung**, which is the part no later frame can recover: by the
/// time the climb refuses the table, the noise parameters are gone.
#[must_use]
pub fn fit_by_profile_scan<M>(
    model: &M,
    cells: &[M::Cell],
    ladder: &[M::NoiseParams],
) -> ScanResult<M::NoiseParams>
where
    M: NoiseModel,
    M::NoiseParams: Clone,
{
    let plans = ploidy_plans(model, cells);
    let winner = scan_ladder(
        model,
        cells,
        ladder,
        &plans,
        |rung, plan, table, climbed_by_ploidy: &mut BTreeMap<Ploidy, Vec<f64>>| {
            let climbed = climb_mixture_weights(table, &plan.cell_weights, &plan.uniform_start);
            // **A rung whose climb ran out of passes is scored below its own summit**,
            // so the rung beside it can win on that alone and the argmax stops being
            // the argmax. There is no channel to report it — `FitTermination` covers
            // the outer alternation of Milestone E and not this — and `MAX_CLIMB_PASSES`
            // is measured on four-cell fixtures where the slowest truth took 1,234
            // passes, against a real table of 583 cells climbed 161 times. Debug only,
            // because a slow climb is a data condition rather than a bug
            // (`spec/parameter_prepass.md` §3.1) and the release path must not abort on
            // one.
            debug_assert!(
                climbed.converged,
                "rung {rung}, ploidy {}: the climb used all {} passes, so this rung is \
                 scored short of its summit and the rung beside it may win on that alone",
                plan.ploidy, climbed.passes
            );
            climbed_by_ploidy.insert(plan.ploidy, climbed.genotype_frequencies);
            climbed.log_likelihood.get()
        },
    );

    ScanResult {
        noise: winner.noise,
        genotype_frequencies: winner.scoring_output,
        log_likelihood: winner.log_likelihood,
        argmax_at_ladder_end: winner.argmax_at_ladder_end,
    }
}

/// What [`fit_by_fixed_frequency_scan`] returned: the winning rung and nothing that was
/// merely handed in.
///
/// **No genotype frequencies, deliberately.** The sibling [`ScanResult`] carries the ones
/// its climb arrived at, which are a fitted quantity. This scan climbs nothing, so the
/// only frequencies it could report are the caller's own argument echoed back — and an
/// echo filed under the same name as a fitted value is how a held-fixed number gets read
/// as a result (owner's call, 2026-08-07).
///
/// Generic over the noise parameters for the same reason [`ScanResult`] is: the two paths
/// scan different things.
#[derive(Clone, PartialEq, Debug)]
pub struct FixedFrequencyScanResult<P> {
    /// Where the winner sat on the ladder, counting from zero.
    ///
    /// **Carried because the coupled fit's stopping rule is stated in rungs** — it stops
    /// when every read group's winning rung is the one it had last iteration
    /// (`arch/parameter_prepass_generic.md` §5.2). Comparing the noise parameters instead
    /// would work only as long as they carry an exact equality, and would say "the same
    /// rung" in a vocabulary the design does not use.
    pub rung: usize,
    /// The winning rung's noise parameters.
    pub noise: P,
    /// The weighted log-likelihood of the cells at that rung, **at the frequencies the
    /// caller handed in**. Comparable across rungs of one scan; comparable across two
    /// scans only if they were handed the same frequencies.
    pub log_likelihood: LogProb,
    /// Whether the answer sat on the ladder's edge — the same bit [`ScanResult`] carries
    /// and for the same reason: a read group whose true rate lies outside the ladder has
    /// its answer clamped to an endpoint and emitted as though it were fitted.
    pub argmax_at_ladder_end: bool,
}

/// Step through `ladder` scoring `cells` at the genotype frequencies handed in, and
/// return the best-scoring rung. **Nothing is climbed.**
///
/// **The `ε` half of the coupled alternation** (`spec/parameter_prepass_generic.md` §5.1):
/// each read group's error rate is fitted from its own table with the sample's genotype
/// frequencies held where the previous iteration left them. That is the procedure the
/// research harness measured — `fit_eps_on_read_group(space, freqs)` in
/// `examples/ng_multilib_key_harness.rs` scores each candidate rate at the frequencies it
/// is handed and never re-climbs them — and from a start at three times the true rates and
/// half the true frequencies its fixed point is the truth in all 25 worlds tried, to 0.000
/// rungs and 0.000% (research note §2.6).
///
/// `genotype_frequencies` carries **one set per ploidy the cells cover**, each as many
/// entries as [`NoiseModel::genotypes`] declares at that ploidy, in the model's own
/// genotype order. One shared set across read groups, not one per group: the frequencies
/// are a property of the individual and the error rates are a property of the chemistry.
///
/// Everything else — every rung scored, the tie rule, the rail flag, the value checks that
/// name the rung — is [`fit_by_profile_scan`]'s, shared rather than restated.
///
/// # Panics
///
/// Everything [`fit_by_profile_scan`] panics on, and additionally if
/// `genotype_frequencies` has no entry for a ploidy the cells cover, if an entry is not as
/// wide as the model declares that ploidy's genotype set, or if an entry is not a
/// probability vector — a negative or non-finite weight, or a set that does not sum to
/// one. Each of those is a caller that built the frequencies from a different table from
/// the cells, and each would otherwise score every rung against a mixture that is not one.
///
/// And if **every** rung scores `−∞`, which this scan can reach where the profile scan
/// cannot. [`scan_ladder`] refuses a weighted cell no genotype of the *model* can have
/// produced; the mixture scored here is the model's likelihoods weighted by the caller's
/// frequencies, and a genotype held at zero — which
/// [`check_genotype_frequencies`] deliberately allows, since a climb may legally leave one
/// there — can make a cell impossible that the model says is fine. Left alone, the tie
/// rule then hands back whichever rung came last with `argmax_at_ladder_end` set: the exact
/// shape of a read group whose true rate lies past the ladder, which is the one thing that
/// flag exists to tell apart.
#[must_use]
pub fn fit_by_fixed_frequency_scan<M>(
    model: &M,
    cells: &[M::Cell],
    ladder: &[M::NoiseParams],
    genotype_frequencies: &BTreeMap<Ploidy, Vec<f64>>,
) -> FixedFrequencyScanResult<M::NoiseParams>
where
    M: NoiseModel,
    M::NoiseParams: Clone,
{
    // **Every check on the frequencies happens here, before the first rung**, which is the
    // rule [`PloidyPlan`] exists to keep: nothing that does not move along the ladder is
    // recomputed at each of its 161 rungs. It also puts the panic where a reader can act on
    // it — "ploidy 2 was handed 2 frequencies" arriving at rung 0 of a scan already under
    // way reads as a fault in the scan rather than in the argument.
    //
    // The ploidies **the cells actually cover**, and not every entry in the map: a caller
    // holding a set for a ploidy this read group never saw is not wrong about anything.
    let plans = ploidy_plans(model, cells);
    for plan in &plans {
        let frequencies = genotype_frequencies
            .get(&plan.ploidy)
            .unwrap_or_else(|| panic!("no genotype frequencies for ploidy {}", plan.ploidy));
        assert_eq!(
            frequencies.len(),
            plan.genotypes,
            "ploidy {} was handed {} genotype frequencies where the model scores it \
             against {} genotypes",
            plan.ploidy,
            frequencies.len(),
            plan.genotypes
        );
        check_genotype_frequencies(plan.ploidy, frequencies);
    }

    // Scratch for the one log-sum-exp per cell, resized when the ploidy changes rather
    // than allocated per cell: the scan walks every cell at each of 161 rungs.
    let mut ln_joint: Vec<f64> = Vec::new();

    let winner = scan_ladder(
        model,
        cells,
        ladder,
        &plans,
        |_rung, plan, table, _scoring_output: &mut ()| {
            let frequencies = &genotype_frequencies[&plan.ploidy];
            ln_joint.clear();
            ln_joint.resize(plan.genotypes, 0.0);
            weighted_log_likelihood(table, &plan.cell_weights, frequencies, &mut ln_joint)
        },
    );

    // **Checked on the winner and not per rung, which is the difference between refusing
    // a degenerate scan and refusing a legal one.** A single rung at `−∞` is an ordinary
    // answer — that rung simply loses, and a model whose parameters make one cell
    // impossible at one rung of 161 has said something true. What cannot be an answer is
    // *every* rung at `−∞`: no rung then outscores another, so `winner` is the last one
    // the loop saw and its rail flag is an accident of the tie rule.
    //
    // `>` and not `!= NEG_INFINITY`, so that a `NaN` is refused here too: a `NaN` score
    // fails every comparison, so the first rung would keep the win against any number the
    // other 160 produced.
    assert!(
        winner.log_likelihood.get() > f64::NEG_INFINITY,
        "no rung of the ladder can have produced these cells at the genotype frequencies \
         handed in — every one of the {} scored {}, so the winner is whichever rung came \
         last rather than an argmax. A genotype held at frequency zero is the usual cause",
        ladder.len(),
        winner.log_likelihood.get()
    );

    FixedFrequencyScanResult {
        rung: winner.rung,
        noise: winner.noise,
        log_likelihood: winner.log_likelihood,
        argmax_at_ladder_end: winner.argmax_at_ladder_end,
    }
}

/// Refuse a set of genotype frequencies that is not a probability vector.
///
/// A zero weight is legal — a genotype the climb drove out of the mixture — so this is
/// weaker than [`mixture_weights`]' interior-point check on a climb's start.
///
/// **`is_finite` cannot change whether a set is accepted, and is kept for what it
/// says.** `NaN >= 0.0` and `−∞ >= 0.0` are both false, so the sign check already refuses
/// those; the one non-finite value that passes it is `+∞`, and a set containing one totals
/// `+∞`, which the sum check always refuses. What the predicate buys is the diagnosis:
/// "genotype 2: inf is not a share of the sample's sites" names the entry to go and look
/// at, where "the frequencies sum to inf rather than one" names neither the genotype nor
/// the fault.
///
/// [`mixture_weights`]: super::mixture_weights
fn check_genotype_frequencies(ploidy: Ploidy, frequencies: &[f64]) {
    let mut total = 0.0;
    for (genotype, &frequency) in frequencies.iter().enumerate() {
        assert!(
            frequency >= 0.0 && frequency.is_finite(),
            "ploidy {ploidy}, genotype {genotype}: {frequency} is not a share of the \
             sample's sites"
        );
        total += frequency;
    }
    assert!(
        (total - 1.0).abs() < FREQUENCY_SUM_TOLERANCE,
        "ploidy {ploidy}: the genotype frequencies sum to {total} rather than one"
    );
}

/// The rung both scans keep, with whatever the scoring step produced beside the score.
///
/// Private, and one type rather than each scan's own: the tie rule and the rail flag are
/// decided here, so a second copy of this shape is a second place they could be decided
/// differently.
struct WinningRung<P, S> {
    rung: usize,
    noise: P,
    /// What the scoring step produced at this rung besides the score itself — the
    /// frequencies it climbed to, per ploidy, or `()` when it climbed nothing.
    scoring_output: S,
    log_likelihood: LogProb,
    argmax_at_ladder_end: bool,
}

/// Everything about the cells that does not change along the ladder, one entry per ploidy
/// they cover, in ascending ploidy order.
///
/// **Built by the two scans rather than inside [`scan_ladder`]**, because a scan that has
/// something to check about its own arguments — as
/// [`fit_by_fixed_frequency_scan`] does about the frequencies it is handed — needs the
/// per-ploidy genotype count to check it with, and checking it inside the rung loop would
/// repeat the check at every one of 161 rungs.
///
/// `BTreeMap` rather than `HashMap` because the per-ploidy scores are added into one total
/// and floating-point addition is not associative, so the order has to be the same on every
/// run.
///
/// # Panics
///
/// If `cells` is empty, if the model declares no genotypes at a ploidy present, or if a
/// ploidy present carries no sites at all.
fn ploidy_plans<M>(model: &M, cells: &[M::Cell]) -> Vec<PloidyPlan>
where
    M: NoiseModel,
{
    assert!(!cells.is_empty(), "a scan needs at least one cell to score");

    let mut cells_of_ploidy: BTreeMap<Ploidy, Vec<usize>> = BTreeMap::new();
    for (position, cell) in cells.iter().enumerate() {
        cells_of_ploidy
            .entry(cell.ploidy())
            .or_default()
            .push(position);
    }

    let mut plans: Vec<PloidyPlan> = Vec::with_capacity(cells_of_ploidy.len());
    for (ploidy, positions) in cells_of_ploidy {
        let genotypes = model.genotypes(ploidy);
        assert!(
            genotypes > 0,
            "the model scores a cell of ploidy {ploidy} against no genotypes at all"
        );
        let cell_weights: Vec<f64> = positions
            .iter()
            .map(|&position| cells[position].sites() as f64)
            .collect();
        let total_sites: f64 = cell_weights.iter().sum();
        assert!(
            total_sites > 0.0,
            "every cell of ploidy {ploidy} holds zero sites, so there is nothing to fit"
        );
        plans.push(PloidyPlan {
            ploidy,
            genotypes,
            positions,
            cell_weights,
            uniform_start: vec![1.0 / genotypes as f64; genotypes],
        });
    }
    plans
}

/// The rung loop both scans run, with the scoring step left to the caller.
///
/// `score_one_ploidy` is handed the rung's index, that ploidy's plan and its row-major
/// likelihood table, and a slot to file whatever it produces beside the score; it returns
/// that ploidy's contribution to the rung's score. A rung's score is the sum over ploidies,
/// and `scoring_output` is discarded for every rung but the winner.
///
/// `plans` comes from [`ploidy_plans`] over the same `model` and `cells`.
///
/// # Panics
///
/// If `ladder` is empty, if the model appends a number of entries other than the
/// [`NoiseModel::genotypes`] it declares, or if what it wrote is not a log-likelihood.
fn scan_ladder<M, S, F>(
    model: &M,
    cells: &[M::Cell],
    ladder: &[M::NoiseParams],
    plans: &[PloidyPlan],
    mut score_one_ploidy: F,
) -> WinningRung<M::NoiseParams, S>
where
    M: NoiseModel,
    M::NoiseParams: Clone,
    S: Default,
    F: FnMut(usize, &PloidyPlan, GenotypeLikelihoodTable<'_>, &mut S) -> f64,
{
    assert!(!ladder.is_empty(), "a scan needs at least one rung to try");

    // Scratch, cleared and refilled once per (rung, ploidy) rather than allocated there:
    // the generic path's cell table has room for 583 cells and walks whatever it holds at
    // each of 161 rungs — the ladder's capacity is the bound, not the count, since only
    // non-empty cells are materialised and the attributed arm adds one entry per listing.
    // The model appends
    // straight into this, so what comes out is the row-major table the climb borrows —
    // no per-cell row and no copy.
    let mut ln_likelihood_row_major: Vec<f64> = Vec::new();

    let mut best: Option<WinningRung<M::NoiseParams, S>> = None;

    for (rung, noise) in ladder.iter().enumerate() {
        let mut rung_log_likelihood = 0.0;
        let mut scoring_output = S::default();

        for plan in plans {
            ln_likelihood_row_major.clear();
            for &position in &plan.positions {
                let before = ln_likelihood_row_major.len();
                model.append_genotype_likelihoods(
                    &cells[position],
                    noise,
                    plan.ploidy,
                    &mut ln_likelihood_row_major,
                );
                assert_eq!(
                    ln_likelihood_row_major.len() - before,
                    plan.genotypes,
                    "the model declares {} genotypes at ploidy {} and appended {} for \
                     cell {position}",
                    plan.genotypes,
                    plan.ploidy,
                    ln_likelihood_row_major.len() - before
                );
            }

            // **What the model wrote is checked here, where the rung is still in
            // scope.** The same faults are refused two frames down, by
            // `GenotypeLikelihoodTable::from_natural_logs` and by the climb — but
            // those messages name a row of a buffer and nothing else, and the rung is
            // the one part no later frame can recover, because the noise parameters
            // are gone by then. A model that goes wrong at one rung of 161 is exactly
            // the case worth localising.
            for ((row, &position), &cell_weight) in ln_likelihood_row_major
                .chunks_exact(plan.genotypes)
                .zip(&plan.positions)
                .zip(&plan.cell_weights)
            {
                for (genotype, &entry) in row.iter().enumerate() {
                    assert!(
                        entry.is_finite() || entry == f64::NEG_INFINITY,
                        "rung {rung}, ploidy {}: cell {position} scored {entry} under \
                         genotype {genotype}, which is not a log-likelihood",
                        plan.ploidy
                    );
                }
                // A cell that carries no sites may legally say no genotype produced
                // it — the climb never looks at it — so the weight is what makes this
                // a fault rather than a shape.
                assert!(
                    cell_weight == 0.0 || row.iter().any(|&entry| entry > f64::NEG_INFINITY),
                    "rung {rung}, ploidy {}: no genotype can have produced cell \
                     {position}, which holds {cell_weight} sites",
                    plan.ploidy
                );
            }

            let table = GenotypeLikelihoodTable::from_natural_logs(
                &ln_likelihood_row_major,
                plan.genotypes,
            );
            rung_log_likelihood += score_one_ploidy(rung, plan, table, &mut scoring_output);
        }

        // The rung's answer is built here rather than at the end from a stored index,
        // so that the rail flag is decided where `rung` is in scope and there is no
        // second, private copy of this shape to keep in step with it.
        let scanned = WinningRung {
            rung,
            noise: noise.clone(),
            scoring_output,
            log_likelihood: LogProb(rung_log_likelihood),
            argmax_at_ladder_end: rung == 0 || rung == ladder.len() - 1,
        };

        // `>=` and not `>`: a tie keeps the **later** rung. See the tie rule above — it
        // is what makes the generic ladder, which descends in error rate, resolve a tie
        // to the lower rate.
        let scores_at_least_as_well = best
            .as_ref()
            .is_none_or(|winner| scanned.log_likelihood >= winner.log_likelihood);
        if scores_at_least_as_well {
            best = Some(scanned);
        }
    }

    best.expect("the ladder is not empty, so some rung won")
}

/// Everything about one ploidy's cells that does not change along the ladder.
///
/// Built once before the rungs, because the cells, their weights and the width of their
/// genotype set are the same at every rung — only the noise parameters move.
struct PloidyPlan {
    ploidy: Ploidy,
    /// How many genotypes this ploidy's cells are scored against, as the model declares
    /// it. Not `ploidy + 1`: that is the SNP/indel path's count and not every path's.
    genotypes: usize,
    /// Where this ploidy's cells sit in the caller's slice.
    positions: Vec<usize>,
    /// One site count per entry of `positions`, in the same order.
    ///
    /// **Parallel to `positions` rather than zipped into it**, because
    /// [`climb_mixture_weights`] takes the weights as one contiguous `&[f64]`; a
    /// `Vec<(usize, f64)>` would have to be unzipped into a fresh buffer at every one of
    /// the ladder's 161 rungs. Nothing can be transposed between the two — one is a
    /// `usize` index and the other an `f64` weight — so what has to hold is only that
    /// they stay the same length and the same order, and both are built in one pass
    /// below.
    cell_weights: Vec<f64>,
    /// The interior point the climb starts from. Where it starts does not matter — the
    /// surface is concave — so it is built once rather than at every rung.
    uniform_start: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use super::*;

    /// A cell of a toy table: which genotype set it belongs to, how many sites landed in
    /// it, and — for the fixture models below — which column of a hand-written matrix
    /// says how likely each genotype makes it.
    #[derive(Clone, Debug)]
    struct ToyCell {
        ploidy: Ploidy,
        sites: u64,
        row: usize,
    }

    impl WeightedCell for ToyCell {
        fn ploidy(&self) -> Ploidy {
            self.ploidy
        }
        fn sites(&self) -> u64 {
            self.sites
        }
    }

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    /// The rail flag is the field with teeth, so it is stated rather than left to a
    /// reader to notice: a scan that railed reports the same shape as one that did not.
    #[test]
    fn a_scan_result_reports_whether_its_answer_sat_on_the_ladders_edge() {
        let diploid = ploidy(2);
        let railed = ScanResult {
            noise: 0.1_f64,
            genotype_frequencies: BTreeMap::from([(diploid, [0.98, 0.015, 0.005].to_vec())]),
            log_likelihood: LogProb(-1.2e9),
            argmax_at_ladder_end: true,
        };

        assert!(railed.argmax_at_ladder_end);
        assert_eq!(railed.genotype_frequencies[&diploid].len(), 3);
    }

    /// A model whose likelihoods come from a table indexed by rung and cell, so a test
    /// can state the whole likelihood surface and know what the answer must be.
    ///
    /// It also **records which rungs it was asked about**, which is how
    /// `every_rung_is_scored_and_none_is_skipped` checks the contract that no early exit
    /// is taken.
    struct TableModel {
        /// `[rung][row][genotype]` — the linear probability, not its log.
        likelihood: Vec<Vec<Vec<f64>>>,
        genotypes: usize,
        asked: RefCell<Vec<usize>>,
    }

    impl NoiseModel for TableModel {
        type Cell = ToyCell;
        type NoiseParams = usize;

        fn genotypes(&self, _ploidy: Ploidy) -> usize {
            self.genotypes
        }

        fn append_genotype_likelihoods(
            &self,
            cell: &ToyCell,
            noise: &usize,
            _ploidy: Ploidy,
            out: &mut Vec<f64>,
        ) {
            self.asked.borrow_mut().push(*noise);
            for likelihood in &self.likelihood[*noise][cell.row] {
                out.push(likelihood.ln());
            }
        }
    }

    /// Three rungs over four cells and three genotypes. Rung 1 is the truth: its columns
    /// are the ones the cell weights below are generated from.
    fn three_rung_table() -> Vec<Vec<Vec<f64>>> {
        vec![
            // rung 0 — a poor explanation of the table
            vec![
                vec![0.40, 0.30, 0.30],
                vec![0.30, 0.30, 0.20],
                vec![0.20, 0.20, 0.30],
                vec![0.10, 0.20, 0.20],
            ],
            // rung 1 — the truth
            vec![
                vec![0.70, 0.10, 0.02],
                vec![0.20, 0.40, 0.08],
                vec![0.07, 0.35, 0.30],
                vec![0.03, 0.15, 0.60],
            ],
            // rung 2 — poor in the other direction
            vec![
                vec![0.10, 0.20, 0.60],
                vec![0.20, 0.30, 0.30],
                vec![0.30, 0.30, 0.08],
                vec![0.40, 0.20, 0.02],
            ],
        ]
    }

    /// The exact cell weights an infinite genome at `truth` would produce under a set of
    /// `columns` — the same device the research harnesses use, so the answer is what the
    /// estimator converges to with no sampling noise in it.
    ///
    /// Takes the columns rather than a table and a rung index, so that the rail test can
    /// ask it for the weights under a truth that is on **no** rung of its ladder.
    fn weights_under(columns: &[Vec<f64>], truth: &[f64], sites: f64) -> Vec<u64> {
        columns
            .iter()
            .map(|row| {
                let probability: f64 = row
                    .iter()
                    .zip(truth)
                    .map(|(&likelihood, &frequency)| frequency * likelihood)
                    .sum();
                (sites * probability).round() as u64
            })
            .collect()
    }

    fn diploid_cells(sites: &[u64]) -> Vec<ToyCell> {
        sites
            .iter()
            .enumerate()
            .map(|(row, &sites)| ToyCell {
                ploidy: ploidy(2),
                sites,
                row,
            })
            .collect()
    }

    /// **The claim the step exists for: a table generated at a known rung recovers that
    /// rung.** The cell weights are the infinite-genome ones under rung 1, so rung 1 is
    /// the truth exactly, and the two rungs either side of it are wrong in opposite
    /// directions — a scan that stopped at the first improvement would take rung 0.
    #[test]
    fn a_table_generated_at_a_known_rung_recovers_it() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000_000.0));

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);

        assert_eq!(scan.noise, 1, "the middle rung is the truth");
        assert!(
            !scan.argmax_at_ladder_end,
            "the winner is interior, so nothing railed"
        );
        let frequencies = &scan.genotype_frequencies[&ploidy(2)];
        assert_eq!(frequencies.len(), 3);
        for (genotype, (&got, &want)) in frequencies.iter().zip(&[0.90, 0.07, 0.03]).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "genotype {genotype}: fitted {got}, truth {want}"
            );
        }
    }

    /// **Every rung is scored — no early exit.** The model records which rungs it was
    /// asked about, so this is checked directly rather than inferred from the answer.
    /// Nobody has shown the profile curve has a single hump, and a scan that stopped
    /// climbing would be sound only if it did.
    #[test]
    fn every_rung_is_scored_and_none_is_skipped() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));

        let _ = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);

        let mut asked: Vec<usize> = model.asked.borrow().clone();
        asked.sort_unstable();
        asked.dedup();
        assert_eq!(asked, vec![0, 1, 2], "a rung went unscored");
    }

    /// Four rungs whose profile curve has two humps. The truth is rung 3; rung 0 is made a
    /// local summit by putting a dip at rungs 1 and 2 between them, so a scan that climbed
    /// and stopped at the first non-improving rung takes rung 0.
    ///
    /// Shared by both scans, so that neither can be given a unimodal fixture and still
    /// claim to prove no early exit is taken.
    fn two_humped_table() -> Vec<Vec<Vec<f64>>> {
        let truth = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let near = vec![
            vec![0.55, 0.15, 0.05],
            vec![0.25, 0.35, 0.10],
            vec![0.12, 0.30, 0.28],
            vec![0.08, 0.20, 0.57],
        ];
        let dip = vec![
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
        ];
        vec![near, dip.clone(), dip, truth]
    }

    /// A curve with two humps, which is the case the scan exists for and the case an
    /// optimiser that climbed would get wrong. Rung 0 is a local best — better than rung
    /// 1 beside it — and rung 3 is the global one.
    #[test]
    fn a_curve_with_two_humps_returns_the_higher_one() {
        let table = two_humped_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[3], &[0.90, 0.07, 0.03], 1_000_000.0));

        // First prove the local summit is really there: over the first two rungs alone,
        // rung 0 wins. Without this the test below could pass on a table with one hump.
        let sub = fit_by_profile_scan(&model, &cells, &[0, 1]);
        assert_eq!(
            sub.noise, 0,
            "rung 0 is not a local summit, so there is one hump"
        );

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2, 3]);

        assert_eq!(scan.noise, 3, "the scan settled for the near rung");
        assert!(scan.argmax_at_ladder_end, "rung 3 is the ladder's last");
    }

    /// A ladder of five rungs whose columns march towards `truth` in equal steps, so the
    /// last rung is the closest the ladder gets and the truth itself is past its end.
    /// `reversed` marches the other way, so the closest rung is the first.
    fn ladder_marching_towards(truth: &[Vec<f64>], reversed: bool) -> Vec<Vec<Vec<f64>>> {
        const RUNGS: usize = 5;
        // **One over the number of *cells*, not the number of genotypes.** A column is one
        // genotype's distribution over the cells, so the flat column that sums to one is
        // `1/cells`. At `1/genotypes` the blend below stops preserving the sum and the
        // renormalisation is doing work it should not have to.
        let flat = 1.0 / truth.len() as f64;
        (0..RUNGS)
            .map(|rung| {
                let step = if reversed { RUNGS - 1 - rung } else { rung };
                // The last rung sits four fifths of the way from flat to the truth, so no
                // rung *is* the truth and the closest one is at the end — which is what
                // makes this a railed fit rather than a recovered one.
                let towards = 0.8 * step as f64 / (RUNGS - 1) as f64;
                let mut blended: Vec<Vec<f64>> = truth
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|&likelihood| flat + towards * (likelihood - flat))
                            .collect()
                    })
                    .collect();
                // Both `flat` and every column of `truth` sum to one, so a convex
                // combination of them does too and this loop is a no-op — kept because it
                // is what makes that a checked property of the fixture rather than an
                // assumption about `flat`.
                //
                // **The first version of this fixture used `1/genotypes` and skipped
                // it, and what went wrong is worth recording, because it is not what a
                // reader would guess.** Unnormalised columns did *not* hand the win to the
                // flattest rung: measured, that rung scored **worst of the five**, about
                // 160,000 nats below the winner. What the extra mass did was move the
                // argmax **one rung inward**, from 4 to 3 — a symptom far easier to read
                // as rounding than as a broken fixture.
                let genotypes = blended[0].len();
                for genotype in 0..genotypes {
                    let column: f64 = blended.iter().map(|row| row[genotype]).sum();
                    for row in &mut blended {
                        row[genotype] /= column;
                    }
                }
                blended
            })
            .collect()
    }

    /// **The rail flag, which is the only thing between a railed fit and a
    /// plausible-looking number.** A table whose truth lies past the ladder's end has
    /// its answer clamped to that end, and the answer looks exactly like a fitted one —
    /// same shape, same score, a large observation count behind it.
    ///
    /// Both ends are checked, because a flag that only watched one of them would pass a
    /// test written at the other. The contrast that gives this test teeth is
    /// `a_table_generated_at_a_known_rung_recovers_it`, where the winner is interior and
    /// the flag is **false** — without that, a flag hard-coded to `true` would pass here.
    #[test]
    fn a_table_whose_truth_lies_past_the_ladder_sets_the_rail_flag() {
        let truth_columns = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let truth = [0.90, 0.07, 0.03];

        for reversed in [false, true] {
            let table = ladder_marching_towards(&truth_columns, reversed);
            let model = TableModel {
                likelihood: table.clone(),
                genotypes: 3,
                asked: RefCell::new(Vec::new()),
            };
            // The weights come from the truth's own columns, which are on no rung.
            let sites = weights_under(&truth_columns, &truth, 1_000_000.0);

            let scan = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2, 3, 4]);

            let expected_end = if reversed { 0 } else { 4 };
            assert_eq!(
                scan.noise, expected_end,
                "reversed = {reversed}: the ladder's closest rung should win"
            );
            assert!(
                scan.argmax_at_ladder_end,
                "reversed = {reversed}: the answer sat on the ladder's edge and said \
                 nothing"
            );
        }

        // And the degenerate ladder is all edge, which is why a one-rung ladder would
        // flag every read group — the failure `ERROR_RATE_LADDER_RUNGS` guards against.
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));
        assert!(fit_by_profile_scan(&model, &cells, &[0]).argmax_at_ladder_end);
    }

    /// **One noise parameter across every ploidy its reads covered, and one set of
    /// frequencies per ploidy.** A haploid cell has two genotypes to mix and a diploid
    /// three, so they cannot share a weight vector — the scan climbs once per ploidy and
    /// adds the two scores.
    #[test]
    fn cells_of_two_ploidies_are_scored_against_their_own_genotype_sets() {
        /// A model whose genotype count really does depend on the ploidy, so a scan that
        /// used one width for both would build a table of the wrong shape.
        struct DosageModel;
        impl NoiseModel for DosageModel {
            type Cell = ToyCell;
            type NoiseParams = f64;

            fn genotypes(&self, ploidy: Ploidy) -> usize {
                usize::from(ploidy.get()) + 1
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &f64,
                ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                // A crude but honest per-genotype likelihood: `row` stands for how many
                // reads showed the alternative allele out of four, and `noise` for the
                // per-read error rate.
                for alt_copies in 0..=ploidy.get() {
                    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
                    let p = carried * (1.0 - noise / 3.0) + (1.0 - carried) * noise;
                    let alt = cell.row as u32;
                    out.push(f64::from(alt).mul_add(p.ln(), f64::from(4 - alt) * (1.0 - p).ln()));
                }
            }
        }

        let cells: Vec<ToyCell> = (0..=4)
            .map(|row| ToyCell {
                ploidy: ploidy(2),
                sites: if row == 0 { 10_000 } else { 20 },
                row,
            })
            .chain((0..=4).map(|row| ToyCell {
                ploidy: ploidy(1),
                sites: if row == 0 { 5_000 } else { 10 },
                row,
            }))
            .collect();

        let scan = fit_by_profile_scan(&DosageModel, &cells, &[0.01, 0.001, 0.0001]);

        assert_eq!(
            scan.genotype_frequencies.len(),
            2,
            "one set of frequencies per ploidy"
        );
        assert_eq!(
            scan.genotype_frequencies[&ploidy(1)].len(),
            2,
            "a haploid has two genotype classes"
        );
        assert_eq!(
            scan.genotype_frequencies[&ploidy(2)].len(),
            3,
            "a diploid has three"
        );
        for (ploidy, weights) in &scan.genotype_frequencies {
            let total: f64 = weights.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-9,
                "ploidy {ploidy}: the frequencies sum to {total}"
            );
        }
    }

    /// **A tie resolves to the later rung, and on the generic ladder that is the lower
    /// error rate.** The two halves are asserted together on purpose: the scan cannot
    /// see which way its ladder runs, so the positional rule and the ladder's direction
    /// are one fact split across two files, and reversing either without the other would
    /// quietly reverse the answer.
    #[test]
    fn a_tie_resolves_to_the_lower_error_rate_on_the_generic_ladder() {
        use crate::ng::parameter_estimation::generic::error_rate_ladder;

        let ladder = error_rate_ladder();
        assert!(
            ladder[0].get() > ladder[ladder.len() - 1].get(),
            "the generic ladder descends in error rate, which is what makes the \
             positional tie rule the lower-rate one"
        );

        // Three rungs whose columns are identical, so every rung scores the same and the
        // whole ladder is one tie.
        let flat = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let table = vec![flat.clone(), flat.clone(), flat];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[0], &[0.90, 0.07, 0.03], 1_000.0));

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);
        assert_eq!(scan.noise, 2, "a tie keeps the last of the tied rungs");
    }

    /// The score reported is the sum over ploidies at the winning rung, which is what
    /// makes "best-scoring iterate" a defined comparison in the coupled fit of E2.
    #[test]
    fn the_reported_score_is_the_sum_over_the_ploidies_at_the_winning_rung() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let sites = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0);

        let diploid_only = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2]);

        // The same cells again under a second ploidy label doubles the score, because
        // the two groups are scored independently and added.
        let mut doubled = diploid_cells(&sites);
        doubled.extend(sites.iter().enumerate().map(|(row, &sites)| ToyCell {
            ploidy: ploidy(4),
            sites,
            row,
        }));
        let both = fit_by_profile_scan(&model, &doubled, &[0, 1, 2]);

        assert_eq!(both.noise, diploid_only.noise);
        assert!(
            (both.log_likelihood.get() - 2.0 * diploid_only.log_likelihood.get()).abs() < 1e-6,
            "two identical groups scored {} against one group's {}",
            both.log_likelihood.get(),
            diploid_only.log_likelihood.get()
        );
    }

    /// A model that appends a different number of entries from the count it declares
    /// would reshape the table without changing its length, and the climb would run on
    /// transposed rows. Named rather than absorbed.
    #[test]
    #[should_panic(expected = "declares 3 genotypes at ploidy 2 and appended 2")]
    fn a_model_that_appends_a_width_it_did_not_declare_is_refused() {
        struct Liar;
        impl NoiseModel for Liar {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                _cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                out.extend([0.5_f64.ln(), 0.5_f64.ln()]);
            }
        }

        let _ = fit_by_profile_scan(&Liar, &diploid_cells(&[10, 20]), &[0]);
    }

    /// **The zero-sites guard is per ploidy, and this is what makes its scope
    /// testable.** With one ploidy in the scan, "this ploidy's cells" and "every cell in
    /// the slice" are the same set, so a guard summing the whole slice would pass every
    /// single-ploidy test. A read group covering a haploid contig that contributed no
    /// sites is the ordinary case — one error rate is fitted across every ploidy the
    /// group covered — and the guard is the only thing that names *which* ploidy was
    /// empty. Two frames down the message is "every cell carries zero weight", which
    /// names neither.
    #[test]
    #[should_panic(expected = "every cell of ploidy 1 holds zero sites")]
    fn a_ploidy_holding_no_sites_beside_one_that_does_names_that_ploidy() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let mut cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));
        cells.extend((0..4).map(|row| ToyCell {
            ploidy: ploidy(1),
            sites: 0,
            row,
        }));
        let _ = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);
    }

    /// A model that scores a cell as `NaN` at one rung of three names **that rung**, not
    /// a row of a buffer nobody outside the loop can see. The rung is what no later frame
    /// can recover: by the time the climb refuses the table, the noise parameters are
    /// gone.
    #[test]
    #[should_panic(
        expected = "rung 2, ploidy 2: cell 3 scored NaN under genotype 0, which is not a"
    )]
    fn a_model_that_scores_a_cell_as_nan_names_the_rung_and_the_cell() {
        struct NanAtTheLastRung;
        impl NoiseModel for NanAtTheLastRung {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if *noise == 2 && cell.row == 3 {
                    out.extend([f64::NAN, 0.5_f64.ln(), 0.5_f64.ln()]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let _ = fit_by_profile_scan(
            &NanAtTheLastRung,
            &diploid_cells(&[10, 20, 30, 40]),
            &[0, 1, 2],
        );
    }

    /// A rung at which some weighted cell becomes impossible under every genotype is a
    /// fault in the model, and it is the rung that has to be named — the same cell is
    /// perfectly scorable at the other 160.
    #[test]
    #[should_panic(
        expected = "rung 2, ploidy 2: no genotype can have produced cell 3, which holds 40 sites"
    )]
    fn a_rung_at_which_a_weighted_cell_is_impossible_names_that_rung() {
        struct ImpossibleAtTheLastRung;
        impl NoiseModel for ImpossibleAtTheLastRung {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if *noise == 2 && cell.row == 3 {
                    out.extend([f64::NEG_INFINITY; 3]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let _ = fit_by_profile_scan(
            &ImpossibleAtTheLastRung,
            &diploid_cells(&[10, 20, 30, 40]),
            &[0, 1, 2],
        );
    }

    /// The companion to the two above, and what keeps the impossible-cell check from
    /// being tightened into one that refuses a legal table: a cell carrying **no** sites
    /// may say no genotype produced it, because the climb never looks at it.
    #[test]
    fn a_cell_holding_no_sites_may_say_no_genotype_produced_it() {
        struct ImpossibleWhereEmpty;
        impl NoiseModel for ImpossibleWhereEmpty {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if cell.sites == 0 {
                    out.extend([f64::NEG_INFINITY; 3]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let cells = diploid_cells(&[10, 0, 30, 40]);
        let scan = fit_by_profile_scan(&ImpossibleWhereEmpty, &cells, &[0, 1, 2]);
        assert!(scan.log_likelihood.get().is_finite());
    }

    /// A cell holding no sites is legal on its own — only the ploidy's total is
    /// checked — and it must not move the answer, because a cell nothing landed in is no
    /// evidence.
    #[test]
    fn a_cell_holding_no_sites_reaches_the_climb_without_moving_the_fit() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let sites = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000_000.0);

        let without = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2]);

        let mut with_empty = diploid_cells(&sites);
        with_empty.push(ToyCell {
            ploidy: ploidy(2),
            sites: 0,
            row: 0,
        });
        let with = fit_by_profile_scan(&model, &with_empty, &[0, 1, 2]);

        assert_eq!(with.noise, without.noise);
        assert_eq!(
            with.log_likelihood.get().to_bits(),
            without.log_likelihood.get().to_bits(),
            "an empty cell moved the score"
        );
    }

    /// A model that scores some ploidy against no genotypes at all is named where the
    /// ploidy is still in scope. Two frames down the message knows only that some table
    /// was empty.
    #[test]
    #[should_panic(expected = "scores a cell of ploidy 2 against no genotypes at all")]
    fn a_model_that_scores_a_ploidy_against_no_genotypes_at_all_is_refused() {
        struct NoGenotypes;
        impl NoiseModel for NoGenotypes {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                0
            }

            fn append_genotype_likelihoods(
                &self,
                _cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                _out: &mut Vec<f64>,
            ) {
            }
        }

        let _ = fit_by_profile_scan(&NoGenotypes, &diploid_cells(&[10, 20]), &[0]);
    }

    /// **The per-ploidy scores are added in a fixed order whatever order the cells
    /// arrive in**, which is what the `BTreeMap` is for: floating-point addition is not
    /// associative, so two runs that summed the ploidies in different orders would return
    /// different bits. Asserted on the bits and not against a tolerance, because a
    /// tolerance is exactly what this cannot afford — the scan compares rungs on this
    /// number and adjacent rungs can be close.
    ///
    /// **Six ploidies, and their site counts spread over four orders of magnitude.** Both
    /// are what make the test able to fail. At two groups a `HashMap` happens to iterate
    /// them in the same order often enough to pass; at six the chance is one in 720. And
    /// with equal group scores the sum is order-independent whatever the container does,
    /// so the counts are scaled ×1 to ×10⁵ to force the partial sums apart in the last
    /// bits.
    #[test]
    fn the_per_ploidy_scores_are_summed_in_a_fixed_order_whatever_the_slice_order() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let base = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0);

        let groups: Vec<Vec<ToyCell>> = (1..=6_u8)
            .map(|copies| {
                let scale = 10_u64.pow(u32::from(copies) - 1);
                base.iter()
                    .enumerate()
                    .map(|(row, &sites)| ToyCell {
                        ploidy: ploidy(copies),
                        sites: sites * scale + 1,
                        row,
                    })
                    .collect()
            })
            .collect();

        let grouped: Vec<ToyCell> = groups.iter().flatten().cloned().collect();
        // The same cells, one from each ploidy in turn, so a container that iterated in
        // insertion order would add the six partial sums the other way round.
        let mut interleaved: Vec<ToyCell> = Vec::with_capacity(grouped.len());
        for row in 0..base.len() {
            for group in groups.iter().rev() {
                interleaved.push(group[row].clone());
            }
        }

        let first = fit_by_profile_scan(&model, &grouped, &[0, 1, 2]);
        let second = fit_by_profile_scan(&model, &interleaved, &[0, 1, 2]);

        assert_eq!(first.genotype_frequencies.len(), 6);
        assert_eq!(
            first.log_likelihood.get().to_bits(),
            second.log_likelihood.get().to_bits(),
            "the ploidies were summed in a different order: {} against {}",
            first.log_likelihood.get(),
            second.log_likelihood.get()
        );
    }

    /// A ladder carrying the same rung twice is the other way a tie arises, and the one
    /// an alternating fit re-scanning a narrowed ladder would produce. The rule is the
    /// same: the later of the two wins.
    #[test]
    fn a_ladder_with_a_repeated_rung_returns_the_later_of_the_two() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));

        // Rung 1 is the truth and appears at positions 1 and 2 of the ladder.
        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 1]);

        assert_eq!(scan.noise, 1);
        assert!(
            scan.argmax_at_ladder_end,
            "the later of the two tied rungs is the ladder's last"
        );
    }

    #[test]
    #[should_panic(expected = "at least one cell")]
    fn a_scan_over_no_cells_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &[], &[0]);
    }

    #[test]
    #[should_panic(expected = "at least one rung")]
    fn a_scan_over_no_rungs_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &diploid_cells(&[10, 20]), &[]);
    }

    #[test]
    #[should_panic(expected = "every cell of ploidy 2 holds zero sites")]
    fn a_ploidy_whose_cells_hold_no_sites_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &diploid_cells(&[0, 0]), &[0, 1, 2]);
    }

    // -----------------------------------------------------------------
    // The sibling: a scan at frequencies handed in, climbing nothing.
    // -----------------------------------------------------------------

    /// The frequency map the sibling takes, for one ploidy.
    fn frequencies_at(copies: u8, frequencies: &[f64]) -> BTreeMap<Ploidy, Vec<f64>> {
        BTreeMap::from([(ploidy(copies), frequencies.to_vec())])
    }

    /// **The claim E1 exists for: a table generated at a known rung recovers that rung
    /// when the scan is handed the frequencies it was generated at.** Same fixture as
    /// `a_table_generated_at_a_known_rung_recovers_it`, so the two answers can be read
    /// against each other: the profile scan finds this rung by climbing to those
    /// frequencies, and this one finds it by being told them.
    #[test]
    fn a_table_generated_at_a_known_rung_recovers_it_from_its_own_frequencies() {
        let table = three_rung_table();
        let truth = [0.90, 0.07, 0.03];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &truth, 1_000_000.0));

        let scan =
            fit_by_fixed_frequency_scan(&model, &cells, &[0, 1, 2], &frequencies_at(2, &truth));

        assert_eq!(scan.noise, 1, "the middle rung is the truth");
        assert_eq!(scan.rung, 1, "and it is the middle rung of the ladder");
        assert!(
            !scan.argmax_at_ladder_end,
            "the winner is interior, so nothing railed"
        );
    }

    /// **The frequencies handed in decide the answer, and this is the test that says
    /// so.** Two genotypes, two cells, three rungs: genotype 0 is explained best by rung
    /// 0 and genotype 1 by rung 2, so the winner is whichever genotype the caller says
    /// the sample is made of.
    ///
    /// Two mutants die here — not only here, but this is where the mechanism is legible.
    /// A scan that **ignored** its frequency argument returns one rung for both halves. A
    /// scan that **climbed** returns rung 2 for both, because the two end rungs' scores
    /// maximised over π tie at `ln 0.9` and the tie rule keeps the later — which is also
    /// why the ladder has a middle rung the answer never lands on: without it a mutant
    /// picking "the last rung" would be right half the time by construction.
    #[test]
    fn the_frequencies_handed_in_decide_which_rung_wins() {
        let leaning = |towards_second: f64| {
            vec![
                vec![1.0 - towards_second, towards_second],
                vec![1.0 - towards_second, towards_second],
            ]
        };
        let table = vec![leaning(0.1), leaning(0.5), leaning(0.9)];
        let model = TableModel {
            likelihood: table,
            genotypes: 2,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&[600, 400]);

        let all_first = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[0, 1, 2],
            &frequencies_at(2, &[1.0, 0.0]),
        );
        let all_second = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[0, 1, 2],
            &frequencies_at(2, &[0.0, 1.0]),
        );

        assert_eq!(
            all_first.rung, 0,
            "a sample made of genotype 0 is explained best by the rung that makes \
             genotype 0 likely"
        );
        assert_eq!(
            all_second.rung, 2,
            "and one made of genotype 1 by the rung at the other end"
        );
    }

    /// **Nothing is climbed**, asserted as a score rather than inferred from the absence
    /// of a frequency field. One rung, so both scans score the same table; the profile
    /// scan climbs to the frequencies that maximise it and this one is handed a set far
    /// from them, so its score must be **strictly lower**. A sibling that quietly climbed
    /// would report the same number to the last bit.
    #[test]
    fn a_scan_at_fixed_frequencies_scores_below_the_climb_it_declines_to_do() {
        let table = three_rung_table();
        let truth = [0.90, 0.07, 0.03];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &truth, 1_000.0));

        let climbed = fit_by_profile_scan(&model, &cells, &[1]);
        let told = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[1],
            &frequencies_at(2, &[0.20, 0.30, 0.50]),
        );

        assert!(
            told.log_likelihood.get() < climbed.log_likelihood.get() - 1.0,
            "handed {} against the climb's {}",
            told.log_likelihood.get(),
            climbed.log_likelihood.get()
        );

        // And handed the frequencies the climb arrived at, the two agree — which is what
        // makes the inequality above a statement about the climb and not about this
        // fixture being hard to score.
        let at_the_summit = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[1],
            &BTreeMap::from([(ploidy(2), climbed.genotype_frequencies[&ploidy(2)].clone())]),
        );
        assert!(
            (at_the_summit.log_likelihood.get() - climbed.log_likelihood.get()).abs() < 1e-9,
            "at the climb's own answer this scan scored {} against {}",
            at_the_summit.log_likelihood.get(),
            climbed.log_likelihood.get()
        );
    }

    /// **Every rung is asked about**, checked by asking the model which rungs it was
    /// shown.
    ///
    /// That is narrower than "no early exit" and deliberately so: an exit that breaks out
    /// of the loop *after* scoring a non-improving rung leaves this model having seen every
    /// rung, and the answer still wrong. What refuses that is
    /// `a_two_humped_curve_at_fixed_frequencies_returns_the_higher_hump` below.
    #[test]
    fn every_rung_is_scored_at_fixed_frequencies_too() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));

        let _ = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[0, 1, 2],
            &frequencies_at(2, &[0.90, 0.07, 0.03]),
        );

        let mut asked = model.asked.borrow().clone();
        asked.sort_unstable();
        asked.dedup();
        assert_eq!(asked, vec![0, 1, 2], "every rung was scored");
    }

    /// **No early exit, at fixed frequencies too** — the claim the rung-recording test
    /// above cannot make. Rung 0 is a local summit and rung 3 the global one, so a scan
    /// that stopped at the first non-improving rung returns rung 0 while still having
    /// asked the model about every rung it scored.
    ///
    /// The profile scan has the same test over the same fixture. Both are kept even though
    /// the rung loop is shared, because the day either scan grows a loop of its own is the
    /// day the shared test stops covering it.
    #[test]
    fn a_two_humped_curve_at_fixed_frequencies_returns_the_higher_hump() {
        let table = two_humped_table();
        let truth = [0.90, 0.07, 0.03];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[3], &truth, 1_000_000.0));
        let frequencies = frequencies_at(2, &truth);

        // The local summit is proven to be there first, exactly as the profile scan's
        // version does it: over the first two rungs alone, rung 0 wins.
        let sub = fit_by_fixed_frequency_scan(&model, &cells, &[0, 1], &frequencies);
        assert_eq!(
            sub.rung, 0,
            "rung 0 is not a local summit, so there is one hump"
        );

        let scan = fit_by_fixed_frequency_scan(&model, &cells, &[0, 1, 2, 3], &frequencies);

        assert_eq!(scan.rung, 3, "the scan settled for the near rung");
    }

    /// **The rail flag, on the sibling.** Its own test rather than a note that the loop
    /// is shared: the flag is the one bit between a railed fit and a plausible-looking
    /// number, and E1 emits it per read group.
    ///
    /// Both ends, for the reason the profile scan's version gives — a flag watching one
    /// end passes a test written at the other — and the interior winner of
    /// `a_table_generated_at_a_known_rung_recovers_it_from_its_own_frequencies` is what
    /// stops a flag hard-coded to `true` from passing.
    #[test]
    fn a_fixed_frequency_scan_whose_truth_lies_past_the_ladder_sets_the_rail_flag() {
        let truth_columns = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let truth = [0.90, 0.07, 0.03];

        for reversed in [false, true] {
            let table = ladder_marching_towards(&truth_columns, reversed);
            let model = TableModel {
                likelihood: table,
                genotypes: 3,
                asked: RefCell::new(Vec::new()),
            };
            let sites = weights_under(&truth_columns, &truth, 1_000_000.0);

            let scan = fit_by_fixed_frequency_scan(
                &model,
                &diploid_cells(&sites),
                &[0, 1, 2, 3, 4],
                &frequencies_at(2, &truth),
            );

            let expected_end = if reversed { 0 } else { 4 };
            assert_eq!(scan.rung, expected_end, "reversed = {reversed}");
            assert!(scan.argmax_at_ladder_end, "reversed = {reversed}");
        }
    }

    /// A tie keeps the later rung here too, which on the generic ladder is the lower
    /// error rate. Stated on the sibling as well as on the profile scan because E1 and
    /// E2 read this one, and a rule that held in only one of the two would be a silent
    /// disagreement between the alternation's two halves.
    #[test]
    fn a_tie_in_a_fixed_frequency_scan_keeps_the_later_rung() {
        let flat = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let table = vec![flat.clone(), flat.clone(), flat];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let truth = [0.90, 0.07, 0.03];
        let cells = diploid_cells(&weights_under(&table[0], &truth, 1_000.0));

        let scan =
            fit_by_fixed_frequency_scan(&model, &cells, &[0, 1, 2], &frequencies_at(2, &truth));

        assert_eq!(scan.rung, 2, "a tie keeps the last of the tied rungs");
        assert_eq!(scan.noise, 2, "and `noise` is that rung's parameters");
    }

    /// **One set of frequencies per ploidy, each against its own genotype set.** The two
    /// ploidies here have two and three genotypes, so the sets cannot be interchanged
    /// without a width fault — and changing only the haploid set moves the answer, which
    /// is what says the haploid cells were scored against the haploid frequencies rather
    /// than against whichever set the map happened to yield first.
    #[test]
    fn each_ploidy_is_scored_against_its_own_frequency_set() {
        struct DosageModel;
        impl NoiseModel for DosageModel {
            type Cell = ToyCell;
            type NoiseParams = f64;

            fn genotypes(&self, ploidy: Ploidy) -> usize {
                usize::from(ploidy.get()) + 1
            }

            /// A cell's row index is how many alternative reads it showed out of two, and
            /// `noise` is the per-read error rate: the ordinary binomial, with the
            /// heterozygote at a half.
            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &f64,
                ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                let alt_reads = cell.row as u32;
                for alt_copies in 0..=ploidy.get() {
                    let p = f64::from(alt_copies) / f64::from(ploidy.get());
                    let p = p * (1.0 - noise) + (1.0 - p) * noise;
                    let reads = 2;
                    let ways: f64 = if alt_reads == 1 { 2.0 } else { 1.0 };
                    out.push(
                        ways.ln()
                            + f64::from(alt_reads) * p.ln()
                            + f64::from(reads - alt_reads) * (1.0 - p).ln(),
                    );
                }
            }
        }

        // Three cells per ploidy — zero, one and two alternative reads. **The diploid arm
        // is a hundredth the size of the haploid one on purpose, and it is what makes this
        // test able to discriminate at all**: only the haploid frequencies differ between
        // the two scans below, so the diploid cells score the same at every rung in both
        // arms, and at `[500, 300, 200]` their curvature decides the argmax on its own —
        // measured, both scans then return rate 0.05 and the test fails on correct code.
        //
        // **Interchanging the two sets outright is caught earlier and elsewhere**, by the
        // declared-width check, since the sets are two and three wide. What survives that
        // check, and what this test is the killer of, is a haploid cell that never sees
        // the haploid set.
        let cells: Vec<ToyCell> = [(1u8, [900u64, 60, 40]), (2, [5, 3, 2])]
            .into_iter()
            .flat_map(|(copies, sites)| {
                sites
                    .into_iter()
                    .enumerate()
                    .map(move |(row, sites)| ToyCell {
                        ploidy: ploidy(copies),
                        sites,
                        row,
                    })
            })
            .collect();
        let ladder = [0.02_f64, 0.05, 0.10, 0.20];
        let diploid = [0.50, 0.30, 0.20].to_vec();

        let mostly_reference = BTreeMap::from([
            (ploidy(1), [0.99, 0.01].to_vec()),
            (ploidy(2), diploid.clone()),
        ]);
        let evenly_split =
            BTreeMap::from([(ploidy(1), [0.50, 0.50].to_vec()), (ploidy(2), diploid)]);

        let strict = fit_by_fixed_frequency_scan(&DosageModel, &cells, &ladder, &mostly_reference);
        let loose = fit_by_fixed_frequency_scan(&DosageModel, &cells, &ladder, &evenly_split);

        assert!(
            strict.rung > loose.rung,
            "a haploid sample said to be almost all reference must explain its \
             alternative reads as error and take a higher rate: {} against {}",
            strict.noise,
            loose.noise
        );
    }

    /// The rung reported is a **position in the ladder handed in**, which is the identity
    /// the coupled fit's stopping rule rests on: it compares rungs rather than rates.
    ///
    /// **The ladder is `[0, 2, 1]` and not `[0, 1, 2]`**, so that a position and the
    /// parameters at it are different numbers. On an ascending ladder every wrong answer
    /// of the form "report the parameters as the position" is right by coincidence; here
    /// the truth sits at position 2 carrying parameters `1`, and a scan that confused the
    /// two reports 1 against an expected 2.
    #[test]
    fn the_rung_reported_indexes_the_ladder() {
        let table = three_rung_table();
        let truth = [0.90, 0.07, 0.03];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &truth, 1_000_000.0));
        let ladder = [0, 2, 1];

        let scan = fit_by_fixed_frequency_scan(&model, &cells, &ladder, &frequencies_at(2, &truth));

        assert_eq!(scan.noise, 1, "rung 1 of the fixture is the truth");
        assert_eq!(
            scan.rung, 2,
            "and the truth sits at position 2 of this ladder"
        );
        assert_eq!(ladder[scan.rung], scan.noise);
    }

    #[test]
    #[should_panic(expected = "no genotype frequencies for ploidy 2")]
    fn a_ploidy_the_caller_gave_no_frequencies_for_is_named() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&[10, 20]),
            &[0, 1, 2],
            &frequencies_at(4, &[0.90, 0.07, 0.03]),
        );
    }

    /// A frequency set one genotype short would be scored against the first two of three
    /// genotypes and leave the third unweighted, which is a plausible wrong answer rather
    /// than a crash — so it is named at the ploidy.
    #[test]
    #[should_panic(expected = "ploidy 2 was handed 2 genotype frequencies")]
    fn a_frequency_set_that_is_not_as_wide_as_the_genotype_set_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&[10, 20]),
            &[0, 1, 2],
            &frequencies_at(2, &[0.93, 0.07]),
        );
    }

    #[test]
    #[should_panic(expected = "the genotype frequencies sum to 0.8 rather than one")]
    fn frequencies_that_are_not_a_probability_vector_are_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&[10, 20]),
            &[0, 1, 2],
            &frequencies_at(2, &[0.70, 0.07, 0.03]),
        );
    }

    /// A genotype the climb drove out of the mixture arrives here as a zero, and that is
    /// legal — it is the interior-point rule of a climb's *start* that this must not
    /// inherit, since nothing here starts anywhere.
    #[test]
    fn a_genotype_at_zero_frequency_is_accepted() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.10, 0.0], 1_000.0));

        let scan = fit_by_fixed_frequency_scan(
            &model,
            &cells,
            &[0, 1, 2],
            &frequencies_at(2, &[0.90, 0.10, 0.0]),
        );

        assert!(scan.log_likelihood.get().is_finite());
    }

    #[test]
    #[should_panic(expected = "genotype 1: -0.1 is not a share of the sample's sites")]
    fn a_negative_frequency_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&[10, 20]),
            &[0, 1, 2],
            &frequencies_at(2, &[1.05, -0.1, 0.05]),
        );
    }

    /// **An infinite frequency is named at the genotype that carries it.** Dropping
    /// `is_finite` from the check leaves the set refused either way — `+∞` makes the total
    /// `+∞`, which the sum check catches — so nothing else in the suite notices its
    /// removal. What it changes is the diagnosis, and that is what this pins.
    ///
    /// The `+∞` sits at genotype 2 rather than genotype 0, so that a check which merely
    /// looked at the first entry could not pass.
    #[test]
    #[should_panic(expected = "genotype 2: inf is not a share of the sample's sites")]
    fn an_infinite_frequency_is_named_at_its_genotype() {
        check_genotype_frequencies(ploidy(2), &[0.5, 0.5, f64::INFINITY]);
    }

    /// A table with a fifth cell that only genotype 0 can have produced, and a truth that
    /// holds no genotype 0 at all. `weights_under` gives that cell zero sites for exactly
    /// that reason, so it is legal; `sites[4]` is what the two tests below disagree about.
    fn table_with_a_cell_only_the_first_genotype_explains() -> (Vec<Vec<Vec<f64>>>, [f64; 3]) {
        let mut table = three_rung_table();
        for rung in &mut table {
            rung.push(vec![1.0, 0.0, 0.0]);
        }
        (table, [0.0, 0.70, 0.30])
    }

    /// **A cell no site landed in must not poison the score at fixed frequencies
    /// either.** The profile scan states this in
    /// `a_cell_holding_no_sites_may_say_no_genotype_produced_it`; this scan reaches
    /// [`weighted_log_likelihood`] by a different door and had no fixture that got there,
    /// so removing the zero-weight skip left every one of its tests green.
    ///
    /// **What makes the failure reachable is a genotype at frequency zero**, which is the
    /// door only this scan has: the climb starts at an interior point, so the profile scan
    /// meets an all-`−∞` row only if the *model* writes one, while here the caller's own
    /// frequency set turns a perfectly finite row into one.
    ///
    /// `0 · −∞` is `NaN`, and a `NaN` does not lose loudly — it fails `>=` against every
    /// rung, so the scan would keep rung 0 whatever the others scored. Both halves are
    /// asserted: a score-only test would pass a scan that returned the wrong rung with a
    /// finite score, and a rung-only test would pass one that returned `NaN` at rung 1.
    #[test]
    fn a_cell_no_genotype_in_the_handed_in_mixture_could_explain_does_not_poison_the_score() {
        let (table, truth) = table_with_a_cell_only_the_first_genotype_explains();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let sites = weights_under(&table[1], &truth, 1_000_000.0);
        assert_eq!(
            sites[4], 0,
            "the impossible cell must be the empty one, or the scan is right to refuse it"
        );

        let scan = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&sites),
            &[0, 1, 2],
            &frequencies_at(2, &truth),
        );

        assert!(
            scan.log_likelihood.get().is_finite(),
            "an empty cell no genotype could explain scored {}",
            scan.log_likelihood.get()
        );
        assert_eq!(scan.rung, 1, "the middle rung is still the truth");
    }

    /// **A frequency set that explains none of the cells is a caller fault, not a railed
    /// fit** — the same fixture as above with 100 sites put in the cell that no genotype
    /// the caller weighted can have produced.
    ///
    /// Every rung then scores `−∞`, the tie rule hands back whichever came last, and
    /// `argmax_at_ladder_end` is `true`: the result is shaped exactly like a read group
    /// whose true rate lies past the ladder, which is the one thing that flag is there to
    /// tell apart. The site count is the whole difference between this test and the one
    /// above, and it is what turns a legal empty cell into an unexplainable one.
    #[test]
    #[should_panic(expected = "no rung of the ladder can have produced these cells")]
    fn a_frequency_set_that_explains_none_of_the_weighted_cells_is_refused() {
        let (table, truth) = table_with_a_cell_only_the_first_genotype_explains();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let mut sites = weights_under(&table[1], &truth, 1_000_000.0);
        sites[4] = 100;

        let _ = fit_by_fixed_frequency_scan(
            &model,
            &diploid_cells(&sites),
            &[0, 1, 2],
            &frequencies_at(2, &truth),
        );
    }
}
