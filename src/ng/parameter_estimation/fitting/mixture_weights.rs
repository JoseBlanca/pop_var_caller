//! The concave climb: given how likely each genotype makes each cell, find the
//! genotype frequencies — the **mixture weights** this file is named for — that best
//! explain the whole table.
//!
//! A sample's sites are a mixture: some fraction of them are homozygous for the
//! reference, some heterozygous, some homozygous for something else. Each genotype
//! makes a given cell more or less likely, and those per-genotype likelihoods are
//! fixed while the climb runs; what is being fitted is only how much of the genome
//! each genotype accounts for.
//!
//! Shared by two of the four fits — the error-rate scan, which climbs to the best
//! frequencies at every rung of its ladder, and the sample's own rates, climbed once
//! on the whole-sample table. The surface is concave, so the climb cannot get **stuck**
//! (`spec/parameter_prepass.md` §3.1); what concavity does not promise is that it
//! arrives quickly, and [`MAX_CLIMB_PASSES`] carries the measurement of how slow it
//! actually gets.
//!
//! **Not the runs model.** Its two states share one inbreeding coefficient, which is a
//! two-dimensional surface inside the simplex rather than a free point on it; a
//! free-simplex maximiser cannot impose that tie, and the concavity below is a
//! statement about the simplex that does not transfer to a curve inside it
//! (`arch/parameter_prepass_generic.md` §4.1).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4.1.

use smallvec::SmallVec;

use crate::ng::types::LogProb;

/// How many passes the climb is allowed before it stops without having settled.
///
/// **Concavity rules out a false summit, not a slow one.** `spec/parameter_prepass.md`
/// §3.1 proves the surface has no local maximum that is not also global, so a climb
/// cannot get stuck; the same section says the *rate* of approach is only linear and is
/// set by how much the mixture components overlap (Redner & Walker, 1984) — worst at low
/// coverage, where a heterozygous site and a homozygous-reference one look alike, which
/// is this estimator's own regime. So a climb that exhausts this cap has not found a
/// second summit. It has run out of time on the only one, and that is a data condition
/// rather than a bug.
///
/// **Measured, not guessed, and the first value was too small.** On the four-cell
/// fixture in this file's tests, reaching [`CLIMB_STILLNESS`] takes 257 passes at a truth
/// of `[0.60, 0.35, 0.05]`, 449 at `[0.90, 0.07, 0.03]` and **1,234** at
/// `[0.80, 0.02, 0.18]`. D1 shipped this at 1,000, which cut the third one off: it
/// returned `converged = false` 3.6 × 10⁻¹⁰ from its answer and the test that read it
/// could not tell, because it went through [`fit_mixture_weights`], which discards the
/// flag. Ten thousand is eight times the slowest measured. The harness this file's
/// arithmetic follows caps at 400 and breaks out silently
/// (`examples/ng_multilib_key_harness.rs`, `climb_frequencies`).
pub const MAX_CLIMB_PASSES: u32 = 10_000;

/// How still the weights have to be before the climb is called settled: the largest move
/// **any one** genotype's weight made in a pass.
///
/// The largest over every genotype and not the last one's: a genotype that no cell can
/// have produced reaches weight zero on the first pass and never moves again, so a
/// measure watching only that coordinate calls the climb finished on pass two with the
/// others still far out — 0.7698 against a truth of 0.95 on the table
/// `a_genotype_impossible_everywhere_does_not_end_the_climb_early` builds.
///
/// Absolute rather than relative, and it has to be read against what the weights are.
/// The two that matter here — heterozygosity and the homozygous-non-reference rate —
/// sit near 0.001, so this is a relative stillness of about one part in 10¹⁰, against
/// the 1% the recovery tests of Milestone F assert. Taken from the same harness
/// (`examples/ng_multilib_key_harness.rs`, `climb_frequencies`).
pub const CLIMB_STILLNESS: f64 = 1e-13;

/// How likely each genotype makes each cell — the table the climb is fitted to, and the
/// only thing about the data the climb sees.
///
/// **One buffer, row-major, with the row width named once.** The entries are laid out
/// one cell after another, each cell contributing `genotypes()` entries in genotype
/// order. Borrowed rather than owned, so the caller keeps the buffer and refills it: the
/// error-rate scan rebuilds this table at every one of its 161 rungs, and a shape made
/// of one slice per cell would force it to rebuild an index of row pointers each time as
/// well — a `Vec<&[f64]>` of a few hundred entries per rung, which the borrow checker
/// will not let it hoist out of the loop because refilling the buffer invalidates it.
///
/// **Natural logarithms**, which is why the constructor says so in its name, and a
/// caller handing this linear probabilities would otherwise get a plausible wrong
/// answer. The corners of a deep table are what force it: a site at the depth cap of 124
/// whose every read showed the alternative allele is `124 · ln ε` under the
/// homozygous-reference genotype, which is `−857` at an error rate of 0.001 and `−1428`
/// at the ladder's floor — both zero in linear space, where the ordinary cell beside them
/// (no alternative read at all) is 0.883. At the ladder's floor a depth of 65 already
/// underflows. `−∞` is a legal entry and means this genotype cannot have produced this
/// cell; it is the only non-finite value accepted.
#[derive(Copy, Clone, Debug)]
pub struct GenotypeLikelihoodTable<'a> {
    ln_likelihood_row_major: &'a [f64],
    genotypes: usize,
}

impl<'a> GenotypeLikelihoodTable<'a> {
    /// Borrow a row-major buffer of **natural-log** likelihoods as a table `genotypes`
    /// wide.
    ///
    /// The width is the only number this takes, and the cell count is derived from it,
    /// so there is no pair of `usize`s here to hand over transposed. A width that does
    /// not divide the buffer is what a raggedly-built table becomes: it cannot be
    /// represented, only refused.
    ///
    /// # Panics
    ///
    /// If `genotypes` is zero, if the buffer is empty, if its length is not a multiple
    /// of `genotypes`, or if any entry is `NaN` or `+∞` — each of those would otherwise
    /// leave a `NaN` to travel through the fit as a plausible number.
    #[must_use]
    pub fn from_natural_logs(ln_likelihood_row_major: &'a [f64], genotypes: usize) -> Self {
        assert!(
            genotypes > 0,
            "a mixture cannot be fitted over an empty genotype set"
        );
        assert!(
            !ln_likelihood_row_major.is_empty(),
            "a mixture cannot be fitted from an empty table"
        );
        assert_eq!(
            ln_likelihood_row_major.len() % genotypes,
            0,
            "a table {genotypes} genotypes wide cannot hold {} entries",
            ln_likelihood_row_major.len()
        );
        for (entry, &ln_likelihood) in ln_likelihood_row_major.iter().enumerate() {
            assert!(
                ln_likelihood.is_finite() || ln_likelihood == f64::NEG_INFINITY,
                "cell {}, genotype {}: {ln_likelihood} is not a log-likelihood",
                entry / genotypes,
                entry % genotypes
            );
        }
        Self {
            ln_likelihood_row_major,
            genotypes,
        }
    }

    /// How many genotypes each cell is scored against — the width of a row, and so the
    /// length of the answer: three for a diploid, five for a tetraploid.
    #[must_use]
    pub fn genotypes(&self) -> usize {
        self.genotypes
    }

    /// How many cells the table holds.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.ln_likelihood_row_major.len() / self.genotypes
    }

    /// One row per cell, each `genotypes()` long.
    fn rows(&self) -> std::slice::ChunksExact<'a, f64> {
        self.ln_likelihood_row_major.chunks_exact(self.genotypes)
    }
}

/// What one climb over the genotype frequencies returned.
///
/// **Crate-internal, deliberately.** The public [`fit_mixture_weights`] returns the
/// weights alone, because `converged` is a flag no consumer would read: on a concave
/// surface a stalled climb is a bug rather than a data condition, and a flag that
/// callers routinely ignore is how a badly-fitted parameter reaches a caller
/// (`spec/parameter_prepass.md` §3.1). It is kept here for the two readers that are
/// not consumers — the profile scan, which needs `log_likelihood` to compare rungs and
/// would otherwise re-score the table, and the tests, which assert the convergence the
/// design claims instead of trusting it.
#[derive(Clone, PartialEq, Debug)]
pub(crate) struct MixtureWeightsFit {
    /// What fraction of the sample's sites each genotype accounts for: one entry per
    /// genotype, non-negative and summing to one.
    ///
    /// **Named for what the numbers are and not for the role they play in the
    /// arithmetic.** They are the mixture's weights, which is what this file's procedure
    /// is named for, but to a reader of the fit they are genotype frequencies — and the
    /// climb carries a *cell* weight per row of the table as well, which is never the
    /// same quantity. `ScanResult::genotype_frequencies` is this vector filed by ploidy,
    /// under the same name.
    pub genotype_frequencies: SmallVec<[f64; 3]>,
    /// The weighted log-likelihood **at `genotype_frequencies`** —
    /// `Σ_cells w · ln Σ_genotypes π·L`.
    /// Computed in a final pass after the climb stopped, so it belongs to the
    /// frequencies returned beside it rather than to the pass before them.
    pub log_likelihood: LogProb,
    /// How many passes ran. One pass is one expectation step and one maximization step.
    pub passes: u32,
    /// Whether the weights stopped moving, as against running out of passes. The
    /// implication runs one way only: `false` means `passes == MAX_CLIMB_PASSES`, while
    /// `true` says nothing about the count, since a climb may settle exactly on its last
    /// allowed pass.
    pub converged: bool,
}

/// Fit the mixing weights of a finite mixture whose component likelihoods are already
/// known and fixed. Returns one weight per genotype, non-negative and summing to one.
///
/// Two consumers, both fitting a free point on the simplex: the error-rate scan, which
/// calls this at every rung of its ladder, and the sample's own rates, climbed once on
/// the whole-sample table (`arch/parameter_prepass_generic.md` §5.1, §5.2).
///
/// `ln_likelihood_by_cell_and_genotype` says how likely each genotype makes each cell;
/// its width is the number of genotypes, so it is also the length of the answer.
///
/// `cell_weights` is how much each cell counts for — its site count, or its site count
/// times a posterior. One entry per cell of the table.
///
/// The climb starts from the uniform point. Where it starts does not matter: with the
/// component likelihoods fixed the surface is concave over the simplex, so there is no
/// local summit that is not also the global one (`spec/parameter_prepass.md` §3.1).
///
/// It runs until no genotype's weight moves by more than [`CLIMB_STILLNESS`] in a pass,
/// or for [`MAX_CLIMB_PASSES`] passes, **whichever comes first** — and a climb that
/// stopped on the second returns its weights with nothing to say so. That is deliberate
/// and the reasoning is on [`MAX_CLIMB_PASSES`]: exhausting the cap means the climb ran
/// out of time on the one summit, not that it found a wrong one.
///
/// # Examples
///
/// Two cells and three genotypes, with the cells' likelihoods laid out one cell after
/// another. The second cell is three times as common as the first, and the third
/// genotype explains it best, so that is where the weight goes:
///
/// ```
/// use pop_var_caller::ng::parameter_estimation::fitting::mixture_weights::{
///     fit_mixture_weights, GenotypeLikelihoodTable,
/// };
///
/// let ln_likelihood: Vec<f64> = [
///     0.8_f64, 0.5, 0.1, // cell 0, under each of the three genotypes
///     0.1, 0.4, 0.9, //     cell 1, likewise
/// ]
/// .iter()
/// .map(|p| p.ln())
/// .collect();
///
/// let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 3);
/// let fitted = fit_mixture_weights(table, &[250.0, 750.0]);
///
/// assert_eq!(fitted.len(), 3);
/// assert!((fitted.iter().sum::<f64>() - 1.0).abs() < 1e-12);
/// assert!(fitted[2] > fitted[0], "{fitted:?}");
/// ```
///
/// # Panics
///
/// If the table and `cell_weights` disagree in length, if a weight is negative or not
/// finite, if every weight is zero, or if a cell that carries weight is one no genotype
/// could have produced. Each of those is a caller error that would otherwise leave a
/// `NaN` to propagate through the fit as a plausible number.
#[must_use]
pub fn fit_mixture_weights(
    ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>,
    cell_weights: &[f64],
) -> SmallVec<[f64; 3]> {
    let genotypes = ln_likelihood_by_cell_and_genotype.genotypes();
    let uniform = vec![1.0 / genotypes as f64; genotypes];
    climb_mixture_weights(ln_likelihood_by_cell_and_genotype, cell_weights, &uniform)
        .genotype_frequencies
}

/// [`fit_mixture_weights`], with the start named and the termination reported.
///
/// The scan calls this rather than the public wrapper so that it can take the score
/// without walking the table again; the tests call it to assert that the climb settled
/// and to start it somewhere other than the uniform point.
///
/// `start` must be a point in the **interior** of the simplex — every weight strictly
/// positive, summing to one. A start on a face pins that genotype's weight at zero for
/// the whole climb, since the expectation step multiplies by it.
///
/// # Panics
///
/// As [`fit_mixture_weights`], and additionally if `start` is not an interior point of
/// the simplex of the right width.
#[must_use]
pub(super) fn climb_mixture_weights(
    ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>,
    cell_weights: &[f64],
    start: &[f64],
) -> MixtureWeightsFit {
    climb_with_cap(
        ln_likelihood_by_cell_and_genotype,
        cell_weights,
        start,
        MAX_CLIMB_PASSES,
    )
}

/// [`climb_mixture_weights`] with the pass cap named rather than taken from the
/// constant.
///
/// **The cap is a parameter here for one reason: so that what it costs is a test rather
/// than a recompile.** Two questions can then be asked of a real table — does the answer
/// move when the cap is raised, and what does a climb that ran out of passes report —
/// and both are asked of the fixtures below. Every caller outside the tests goes through
/// [`climb_mixture_weights`] and gets [`MAX_CLIMB_PASSES`].
#[must_use]
fn climb_with_cap(
    ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>,
    cell_weights: &[f64],
    start: &[f64],
    max_passes: u32,
) -> MixtureWeightsFit {
    check_cell_weights(ln_likelihood_by_cell_and_genotype, cell_weights);
    let genotypes = ln_likelihood_by_cell_and_genotype.genotypes();
    check_start(start, genotypes);

    let total_cell_weight: f64 = cell_weights.iter().sum();

    let mut genotype_frequencies: SmallVec<[f64; 3]> = SmallVec::from_slice(start);
    let mut next_genotype_frequencies: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, genotypes);
    // Scratch, loaded and cleared once per cell rather than allocated per cell. The
    // loop below walks every cell of the table — up to 583 on the generic path, which is
    // the cell ladder's capacity rather than a count
    // (`arch/parameter_prepass_generic.md` §4.1) — once per pass, and the profile scan
    // runs one whole climb at each of its 161 rungs.
    let mut ln_joint: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, genotypes);

    let mut passes = 0;
    let mut converged = false;
    let mut previous_score = f64::NEG_INFINITY;

    while passes < max_passes {
        next_genotype_frequencies.fill(0.0);
        let mut score = 0.0;

        for (cell, (row, &cell_weight)) in ln_likelihood_by_cell_and_genotype
            .rows()
            .zip(cell_weights)
            .enumerate()
        {
            if cell_weight == 0.0 {
                continue;
            }
            for (slot, (&genotype_frequency, &ln_likelihood)) in ln_joint
                .iter_mut()
                .zip(genotype_frequencies.iter().zip(row))
            {
                *slot = genotype_frequency.ln() + ln_likelihood;
            }
            let ln_total = ln_sum_exp(&ln_joint);
            assert!(
                ln_total.is_finite(),
                "cell {cell} carries weight {cell_weight} and no genotype can have \
                 produced it"
            );
            score += cell_weight * ln_total;
            // The expectation step's responsibilities, accumulated straight into the
            // maximization step's numerator rather than materialised per cell.
            for (slot, &ln_j) in next_genotype_frequencies.iter_mut().zip(ln_joint.iter()) {
                *slot += cell_weight * (ln_j - ln_total).exp();
            }
        }

        // Monotone ascent is the one thing expectation-maximization is actually proved
        // to give (Dempster, Laird & Rubin 1977), so a pass that lost ground is a bug
        // in the two steps above rather than a hard table. Skipped on the first pass,
        // where there is nothing to compare against — stated rather than left to
        // `-∞ - ∞ = -∞` to make the check vacuous by arithmetic. The slack is relative
        // because the score scales with the total cell weight: −1,001 over a thousand
        // sites and −1,001,017 over a million.
        debug_assert!(
            passes == 0 || score >= previous_score - 1e-9 * previous_score.abs(),
            "a pass lost ground: {previous_score} → {score}"
        );
        previous_score = score;

        let mut largest_move: f64 = 0.0;
        for (slot, &genotype_frequency) in next_genotype_frequencies
            .iter_mut()
            .zip(genotype_frequencies.iter())
        {
            *slot /= total_cell_weight;
            largest_move = largest_move.max((*slot - genotype_frequency).abs());
        }
        genotype_frequencies.copy_from_slice(&next_genotype_frequencies);
        passes += 1;

        if largest_move < CLIMB_STILLNESS {
            converged = true;
            break;
        }
    }

    let log_likelihood = weighted_log_likelihood(
        ln_likelihood_by_cell_and_genotype,
        cell_weights,
        &genotype_frequencies,
        &mut ln_joint,
    );

    MixtureWeightsFit {
        genotype_frequencies,
        log_likelihood: LogProb(log_likelihood),
        passes,
        converged,
    }
}

/// `Σ_cells w · ln Σ_genotypes π·L`, the quantity the climb maximises.
///
/// Kept separate from the climb because the climb's own running score belongs to the
/// weights it started the pass with, and what a caller compares rungs on is the score
/// at the weights it is handed.
fn weighted_log_likelihood(
    ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>,
    cell_weights: &[f64],
    genotype_frequencies: &[f64],
    ln_joint: &mut [f64],
) -> f64 {
    let mut score = 0.0;
    for (row, &cell_weight) in ln_likelihood_by_cell_and_genotype.rows().zip(cell_weights) {
        // Not an optimisation. A cell that carries no weight may legally be one no
        // genotype could have produced, and `0 · −∞` is `NaN` — which does not lose
        // loudly, it simply never wins, so the profile scan would return whichever rung
        // happened to be scored on a table with no empty cells. A posterior-weighted
        // call, which is what `cell_weights` documents, produces zero-weight cells as a
        // matter of course.
        if cell_weight == 0.0 {
            continue;
        }
        for (slot, (&genotype_frequency, &ln_likelihood)) in ln_joint
            .iter_mut()
            .zip(genotype_frequencies.iter().zip(row))
        {
            *slot = genotype_frequency.ln() + ln_likelihood;
        }
        score += cell_weight * ln_sum_exp(ln_joint);
    }
    score
}

/// `ln Σ exp(terms)`, shifted by the largest term so that a table of likelihoods no
/// `f64` can hold in linear space still sums exactly.
///
/// `−∞` terms contribute nothing and an all-`−∞` slice returns `−∞`, which is the
/// caller's signal that no genotype could have produced the cell.
fn ln_sum_exp(terms: &[f64]) -> f64 {
    let largest_term = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if largest_term == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = terms.iter().map(|&term| (term - largest_term).exp()).sum();
    largest_term + sum.ln()
}

/// Check there is one usable weight per cell of the table.
///
/// The table's own shape was settled when it was borrowed
/// ([`GenotypeLikelihoodTable::from_natural_logs`]); what is left is whether the weights
/// beside it line up with it and are counts of sites.
fn check_cell_weights(
    ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>,
    cell_weights: &[f64],
) {
    assert_eq!(
        ln_likelihood_by_cell_and_genotype.cells(),
        cell_weights.len(),
        "one weight per cell: {} cells against {} weights",
        ln_likelihood_by_cell_and_genotype.cells(),
        cell_weights.len()
    );

    let mut total_cell_weight = 0.0;
    for (cell, &cell_weight) in cell_weights.iter().enumerate() {
        assert!(
            cell_weight.is_finite() && cell_weight >= 0.0,
            "cell {cell} carries weight {cell_weight}, which is not a count of sites"
        );
        total_cell_weight += cell_weight;
    }
    assert!(
        total_cell_weight > 0.0,
        "every cell carries zero weight, so there is nothing to fit"
    );
}

/// Check a start is an interior point of the simplex of the right width.
fn check_start(start: &[f64], genotypes: usize) {
    assert_eq!(
        start.len(),
        genotypes,
        "the start lists {} weights where the table has {genotypes} genotypes",
        start.len()
    );
    let mut total = 0.0;
    for (genotype, &weight) in start.iter().enumerate() {
        assert!(
            weight > 0.0 && weight.is_finite(),
            "genotype {genotype} starts at {weight}, which is not inside the simplex"
        );
        total += weight;
    }
    assert!(
        (total - 1.0).abs() < 1e-9,
        "the start sums to {total} rather than one"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A table whose maximiser is known exactly, and known **without a fit**.
    ///
    /// Each column is one genotype's distribution over the four cells, so it sums to
    /// one; weighting each cell by its probability under a chosen truth,
    /// `w_cell = Σ_j π_j · L(cell | j)`, makes that truth the maximiser exactly — this
    /// is the infinite-genome table, with no sampling noise in it to argue about. It is
    /// the same device the research harnesses use
    /// (`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §1, applied
    /// throughout §2).
    ///
    /// **Four cells and three genotypes, deliberately unequal**, so that a fixture built
    /// with the two dimensions swapped is a buffer the width does not divide rather than
    /// a wrong number. The cell weights it produces run 0.6376 down to 0.0555 — an
    /// 11.5-to-1 spread, which is what makes a climb that ignored them fail rather than
    /// merely drift.
    const CELL_GIVEN_GENOTYPE: [[f64; 3]; 4] = [
        // homozygous reference, heterozygous, homozygous non-reference
        [0.70, 0.10, 0.02],
        [0.20, 0.40, 0.08],
        [0.07, 0.35, 0.30],
        [0.03, 0.15, 0.60],
    ];

    /// [`CELL_GIVEN_GENOTYPE`] in logs, row-major, as the climb takes it.
    fn ln_table() -> Vec<f64> {
        natural_logs_of(&CELL_GIVEN_GENOTYPE)
    }

    /// A table of fixed-width rows of linear probabilities, flattened into the row-major
    /// buffer of natural logs [`GenotypeLikelihoodTable`] borrows.
    fn natural_logs_of<const GENOTYPES: usize>(table: &[[f64; GENOTYPES]]) -> Vec<f64> {
        table
            .iter()
            .flat_map(|row| row.iter().map(|likelihood| likelihood.ln()))
            .collect()
    }

    /// The diploid fixture, borrowed as the three-genotype table it is.
    fn diploid_table(ln_likelihood_row_major: &[f64]) -> GenotypeLikelihoodTable<'_> {
        GenotypeLikelihoodTable::from_natural_logs(ln_likelihood_row_major, 3)
    }

    /// The exact cell weights an infinite genome at `truth` would produce, times a
    /// site count so they read as counts rather than probabilities.
    fn weights_under(truth: &[f64; 3], sites: f64) -> Vec<f64> {
        CELL_GIVEN_GENOTYPE
            .iter()
            .map(|row| {
                sites
                    * row
                        .iter()
                        .zip(truth)
                        .map(|(&likelihood, &frequency)| frequency * likelihood)
                        .sum::<f64>()
            })
            .collect()
    }

    /// A sample's genotype frequencies as they actually look: nearly everything
    /// homozygous for the reference, and the two rare classes an order of magnitude
    /// apart from each other.
    const TRUTH: [f64; 3] = [0.90, 0.07, 0.03];

    /// The claim the whole climb rests on: it finds the frequencies that generated the
    /// table, and it finds them exactly, because the table is the infinite-genome one.
    #[test]
    fn the_climb_recovers_the_frequencies_that_generated_the_table() {
        let ln_likelihood = ln_table();
        let fitted = fit_mixture_weights(
            diploid_table(&ln_likelihood),
            &weights_under(&TRUTH, 1_000_000.0),
        );

        assert_eq!(fitted.len(), 3);
        for (genotype, (&got, &want)) in fitted.iter().zip(&TRUTH).enumerate() {
            assert!(
                (got - want).abs() < 1e-9,
                "genotype {genotype}: fitted {got}, truth {want}"
            );
        }
    }

    /// Concavity is the property that makes the start irrelevant, and it is asserted
    /// rather than assumed: five starts — one in each of the three corners'
    /// neighbourhoods, the uniform point, and one ordinary interior point — reach the
    /// same answer. The heterozygous corner is in the list because it is the coordinate
    /// the estimator exists to measure.
    #[test]
    fn every_interior_start_reaches_the_same_summit() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);
        let cell_weights = weights_under(&TRUTH, 1_000_000.0);

        let starts = [
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [0.98, 0.01, 0.01],
            [0.01, 0.98, 0.01],
            [0.01, 0.01, 0.98],
            [0.20, 0.30, 0.50],
        ];
        for start in starts {
            let fit = climb_mixture_weights(table, &cell_weights, &start);
            assert!(fit.converged, "start {start:?} ran out of passes");
            assert!(
                fit.passes < MAX_CLIMB_PASSES,
                "start {start:?} used every pass"
            );
            for (genotype, (&got, &want)) in fit.genotype_frequencies.iter().zip(&TRUTH).enumerate()
            {
                assert!(
                    (got - want).abs() < 1e-9,
                    "start {start:?}, genotype {genotype}: fitted {got}, truth {want}"
                );
            }
        }
    }

    /// The cell weights are what the table's shape is carried in, so a climb that
    /// ignored them would still return something plausible. Two truths generate two
    /// different weight vectors over the *same* likelihood table, and the climb has to
    /// tell them apart.
    #[test]
    fn two_truths_over_one_likelihood_table_give_two_answers() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);

        let inbred = [0.80, 0.02, 0.18];
        let outcrossing = [0.60, 0.35, 0.05];

        // Through the reporting entry point rather than the public one, so that a climb
        // that ran out of passes fails here rather than passing on the margin. It is
        // this fixture that measured the cap: at `MAX_CLIMB_PASSES = 1_000` the inbred
        // truth stopped 3.6 × 10⁻¹⁰ short with `converged = false`, inside the 10⁻⁹ the
        // loops below assert.
        let uniform = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let first = climb_mixture_weights(table, &weights_under(&inbred, 500_000.0), &uniform);
        let second =
            climb_mixture_weights(table, &weights_under(&outcrossing, 500_000.0), &uniform);

        assert!(
            first.converged,
            "the inbred truth ran out of passes at {}",
            first.passes
        );
        assert!(
            second.converged,
            "the outcrossing truth ran out of passes at {}",
            second.passes
        );
        for (genotype, (&got, &want)) in first.genotype_frequencies.iter().zip(&inbred).enumerate()
        {
            assert!((got - want).abs() < 1e-9, "genotype {genotype}: {got}");
        }
        for (genotype, (&got, &want)) in second
            .genotype_frequencies
            .iter()
            .zip(&outcrossing)
            .enumerate()
        {
            assert!((got - want).abs() < 1e-9, "genotype {genotype}: {got}");
        }
    }

    /// **What a climb that ran out of passes reports, which nothing else reaches.** Two
    /// claims live here and both were unheld: that `converged` is ever `false`, and that
    /// `log_likelihood` is the score at the weights returned beside it rather than at the
    /// weights the last pass started from. On a settled climb those two scores agree to
    /// thirteen digits, so only a stopped-early climb can tell them apart — here the gap
    /// is large enough that the longhand check below has teeth.
    #[test]
    fn a_climb_stopped_early_says_so_and_scores_the_weights_it_returns() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);
        let cell_weights = weights_under(&TRUTH, 1_000_000.0);

        let stopped = climb_with_cap(table, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], 5);

        assert!(!stopped.converged, "five passes was expected to be too few");
        assert_eq!(stopped.passes, 5);

        let mut longhand = 0.0;
        for (row, &cell_weight) in CELL_GIVEN_GENOTYPE.iter().zip(&cell_weights) {
            let mixed: f64 = row
                .iter()
                .zip(&stopped.genotype_frequencies)
                .map(|(&likelihood, &genotype_frequency)| genotype_frequency * likelihood)
                .sum();
            longhand += cell_weight * mixed.ln();
        }
        assert!(
            (stopped.log_likelihood.get() - longhand).abs() < 1e-6,
            "reported {}, score at the weights returned {longhand}",
            stopped.log_likelihood.get()
        );
    }

    /// The cap is a stopping rule, not a shaping one: raise it tenfold and the answer is
    /// the same bits. This is the experiment the pass-cap decision rests on, and it is a
    /// test rather than a recompile because `climb_with_cap` takes the cap.
    #[test]
    fn raising_the_pass_cap_does_not_move_the_answer() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);
        let inbred = [0.80, 0.02, 0.18];
        let cell_weights = weights_under(&inbred, 500_000.0);
        let uniform = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

        let standard = climb_with_cap(table, &cell_weights, &uniform, MAX_CLIMB_PASSES);
        let generous = climb_with_cap(table, &cell_weights, &uniform, 10 * MAX_CLIMB_PASSES);

        assert!(standard.converged && generous.converged);
        assert_eq!(standard.genotype_frequencies, generous.genotype_frequencies);
        assert_eq!(standard.passes, generous.passes);
    }

    /// **The stillness test has to watch every genotype, not the last one.** A genotype
    /// no cell can have produced reaches weight zero on the first pass and never moves
    /// again; a measure that looked only at it would call this climb finished on pass
    /// two, with the first weight at 0.7698 against a truth of 0.95. It is also the only
    /// table here where a weight reaches exactly zero *during* the climb, which is the
    /// state `genotype_frequency.ln()` has to survive.
    #[test]
    fn a_genotype_impossible_everywhere_does_not_end_the_climb_early() {
        let cell_given_genotype: [[f64; 3]; 3] =
            [[0.70, 0.10, 0.0], [0.20, 0.40, 0.0], [0.10, 0.50, 0.0]];
        let truth = [0.95, 0.05, 0.0];

        let ln_likelihood = natural_logs_of(&cell_given_genotype);
        let table = diploid_table(&ln_likelihood);
        let cell_weights: Vec<f64> = cell_given_genotype
            .iter()
            .map(|row| {
                1_000_000.0
                    * row
                        .iter()
                        .zip(&truth)
                        .map(|(&likelihood, &frequency)| frequency * likelihood)
                        .sum::<f64>()
            })
            .collect();

        let fit = climb_mixture_weights(table, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);

        assert!(fit.converged);
        assert!(
            fit.passes > 100,
            "the climb finished in {} passes, which is too few to have crawled there",
            fit.passes
        );
        for (genotype, (&got, &want)) in fit.genotype_frequencies.iter().zip(&truth).enumerate() {
            assert!(
                (got - want).abs() < 1e-9,
                "genotype {genotype}: fitted {got}, truth {want}"
            );
        }
    }

    /// The answer really is the summit, checked against the objective rather than
    /// against the construction that produced the fixture: each of the eighteen steps
    /// away from it along an edge of the simplex scores lower. It survives a change of
    /// fixture, which the recovery tests do not.
    ///
    /// **What it does not check is the scorer.** Both points are scored by the same
    /// function, so a `weighted_log_likelihood` that was wrong in a way common to every
    /// point passes this untouched. That half is held by
    /// `the_reported_score_belongs_to_the_weights_returned_with_it`, which recomputes the
    /// score longhand in linear space.
    #[test]
    fn no_move_away_from_the_fitted_weights_scores_higher() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);
        let cell_weights = weights_under(&TRUTH, 1_000_000.0);

        let fit = climb_mixture_weights(table, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        let mut ln_joint = vec![0.0; 3];

        for from in 0..3 {
            for to in 0..3 {
                if from == to {
                    continue;
                }
                for step in [1e-4, 1e-3, 1e-2] {
                    let mut nudged: Vec<f64> = fit.genotype_frequencies.to_vec();
                    if nudged[from] <= step {
                        continue;
                    }
                    nudged[from] -= step;
                    nudged[to] += step;
                    let score =
                        weighted_log_likelihood(table, &cell_weights, &nudged, &mut ln_joint);
                    assert!(
                        score <= fit.log_likelihood.get(),
                        "moving {step} from genotype {from} to {to} scored \
                         {score} against the summit's {}",
                        fit.log_likelihood.get()
                    );
                }
            }
        }
    }

    /// The score reported beside the weights is the score **at** those weights — the
    /// number the profile scan compares rungs on. Checked against the formula written
    /// out longhand, not against the climb's own running total.
    #[test]
    fn the_reported_score_belongs_to_the_weights_returned_with_it() {
        let ln_likelihood = ln_table();
        let table = diploid_table(&ln_likelihood);
        let cell_weights = weights_under(&TRUTH, 1_000.0);

        let fit = climb_mixture_weights(table, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);

        let mut longhand = 0.0;
        for (row, &cell_weight) in CELL_GIVEN_GENOTYPE.iter().zip(&cell_weights) {
            let mixed: f64 = row
                .iter()
                .zip(&fit.genotype_frequencies)
                .map(|(&likelihood, &genotype_frequency)| genotype_frequency * likelihood)
                .sum();
            longhand += cell_weight * mixed.ln();
        }

        assert!(
            (fit.log_likelihood.get() - longhand).abs() < 1e-6,
            "reported {}, longhand {longhand}",
            fit.log_likelihood.get()
        );
    }

    /// The weights are a point on the simplex whatever the table said, which is what
    /// every consumer reads them as. A maximization step that forgot to divide by the
    /// total weight would return counts here and nothing else in the fit would notice.
    #[test]
    fn the_fitted_weights_are_a_point_on_the_simplex() {
        let ln_likelihood = ln_table();
        let cell_weights = weights_under(&TRUTH, 7.0);
        let fitted = fit_mixture_weights(diploid_table(&ln_likelihood), &cell_weights);

        let total: f64 = fitted.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "the weights sum to {total}");
        assert!(fitted.iter().all(|&weight| (0.0..=1.0).contains(&weight)));
    }

    /// Nothing here is diploid. A tetraploid has five genotypes, which is past the
    /// three a `SmallVec<[f64; 3]>` holds inline, so this exercises the spill as well
    /// as the loop bound.
    #[test]
    fn a_tetraploid_table_is_fitted_over_its_own_five_genotypes() {
        // Six cells across, five dosages of the alternative allele down: entry
        // `[cell][dosage]` is how likely that dosage makes that cell, and each dosage's
        // column sums to one. The cells are ordered by how many reads showed the
        // alternative allele, so each dosage concentrates on its own.
        let cell_given_dosage: [[f64; 5]; 6] = [
            [0.90, 0.05, 0.01, 0.00, 0.00],
            [0.08, 0.80, 0.08, 0.01, 0.00],
            [0.02, 0.10, 0.80, 0.04, 0.00],
            [0.00, 0.04, 0.08, 0.80, 0.02],
            [0.00, 0.01, 0.02, 0.10, 0.08],
            [0.00, 0.00, 0.01, 0.05, 0.90],
        ];
        let truth = [0.70, 0.15, 0.08, 0.05, 0.02];

        let ln_likelihood = natural_logs_of(&cell_given_dosage);
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 5);
        let cell_weights: Vec<f64> = cell_given_dosage
            .iter()
            .map(|row| {
                1_000_000.0
                    * row
                        .iter()
                        .zip(&truth)
                        .map(|(&likelihood, &frequency)| frequency * likelihood)
                        .sum::<f64>()
            })
            .collect();

        let uniform = vec![0.2; 5];
        let fit = climb_mixture_weights(table, &cell_weights, &uniform);
        assert!(fit.converged, "the climb ran out of passes");
        let fitted = fit.genotype_frequencies;

        assert_eq!(fitted.len(), 5);
        assert_eq!(
            table.cells(),
            6,
            "six cells, five genotypes, not the reverse"
        );
        for (dosage, (&got, &want)) in fitted.iter().zip(&truth).enumerate() {
            assert!(
                (got - want).abs() < 1e-8,
                "dosage {dosage}: fitted {got}, truth {want}"
            );
        }
    }

    /// A genotype that cannot have produced a cell is `−∞` there, and that is a legal
    /// entry rather than a rejected one: it is what a homozygous-reference genotype
    /// says about a cell where every read showed something else at a zero error rate.
    #[test]
    fn a_genotype_that_cannot_have_produced_a_cell_is_allowed_to_say_so() {
        let ln_likelihood = natural_logs_of(&[[0.0, 0.5, 0.9], [1.0, 0.5, 0.1]]);
        let fitted = fit_mixture_weights(diploid_table(&ln_likelihood), &[100.0, 900.0]);

        assert!(fitted.iter().all(|weight| weight.is_finite()));
        let total: f64 = fitted.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "the weights sum to {total}");
    }

    /// A cell no genotype can have produced is a `NaN` waiting to happen — the mixture
    /// over it is zero, and every responsibility is `0/0`. It is a fault in whatever
    /// built the table, so it is named rather than propagated.
    #[test]
    #[should_panic(expected = "cell 1 carries weight 900 and no genotype can have produced it")]
    fn a_cell_no_genotype_can_produce_is_refused_rather_than_scored() {
        let ln_likelihood = [
            0.5_f64.ln(),
            0.5_f64.ln(),
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let _ = fit_mixture_weights(table, &[100.0, 900.0]);
    }

    /// A cell that carries no weight is skipped, so it may say anything at all —
    /// including that no genotype produced it.
    ///
    /// **And the score has to come back finite**, which is the half a fit that returns
    /// only the weights cannot see. `0 · −∞` is `NaN`, and a `NaN` score does not lose
    /// loudly in the profile scan — it simply never wins. That is why this goes through
    /// the reporting entry point.
    #[test]
    fn a_cell_carrying_no_weight_is_skipped_whatever_it_says() {
        let ln_likelihood = [
            0.5_f64.ln(),
            0.5_f64.ln(),
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let fit = climb_mixture_weights(table, &[100.0, 0.0], &[0.5, 0.5]);

        assert!(
            fit.log_likelihood.get().is_finite(),
            "the score came back {}",
            fit.log_likelihood.get()
        );
        let total: f64 = fit.genotype_frequencies.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "the weights sum to {total}");
    }

    /// Two cells and two genotypes as a flat buffer, for the refusals below.
    fn two_by_two() -> [f64; 4] {
        [0.5_f64.ln(), 0.5_f64.ln(), 0.5_f64.ln(), 0.5_f64.ln()]
    }

    #[test]
    #[should_panic(expected = "one weight per cell: 2 cells against 3 weights")]
    fn a_weight_without_a_cell_is_refused() {
        let ln_likelihood = two_by_two();
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let _ = fit_mixture_weights(table, &[1.0, 1.0, 1.0]);
    }

    /// What a raggedly-built table becomes once the rows live in one buffer: a length
    /// the width does not divide. It cannot reach the climb at all.
    #[test]
    #[should_panic(expected = "a table 2 genotypes wide cannot hold 5 entries")]
    fn a_buffer_the_genotype_count_does_not_divide_is_refused() {
        let ln_likelihood = [0.5_f64.ln(); 5];
        let _ = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
    }

    #[test]
    #[should_panic(expected = "cell 1, genotype 0: NaN is not a log-likelihood")]
    fn a_nan_likelihood_is_refused() {
        let ln_likelihood = [0.5_f64.ln(), 0.5_f64.ln(), f64::NAN, 0.5_f64.ln()];
        let _ = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
    }

    /// `+∞` is not a log-likelihood either, and it is the one that would survive
    /// `ln_sum_exp` as a `NaN` rather than as an obvious wrong answer.
    #[test]
    #[should_panic(expected = "cell 0, genotype 1: inf is not a log-likelihood")]
    fn a_positive_infinity_likelihood_is_refused() {
        let ln_likelihood = [0.5_f64.ln(), f64::INFINITY, 0.5_f64.ln(), 0.5_f64.ln()];
        let _ = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
    }

    #[test]
    #[should_panic(expected = "cell 0 carries weight -1, which is not a count of sites")]
    fn a_negative_cell_weight_is_refused() {
        let ln_likelihood = two_by_two();
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let _ = fit_mixture_weights(table, &[-1.0, 2.0]);
    }

    /// The companion to the negative weight above, and the one that pins the
    /// `is_finite()` half of the check: `+∞ >= 0.0` is true, so a lower bound alone
    /// would let an infinite site count through and the total weight the maximization
    /// step divides by would be `∞`.
    #[test]
    #[should_panic(expected = "cell 1 carries weight inf, which is not a count of sites")]
    fn an_infinite_cell_weight_is_refused() {
        let ln_likelihood = two_by_two();
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let _ = fit_mixture_weights(table, &[1.0, f64::INFINITY]);
    }

    #[test]
    #[should_panic(expected = "every cell carries zero weight")]
    fn a_table_of_empty_cells_is_refused() {
        let ln_likelihood = two_by_two();
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 2);
        let _ = fit_mixture_weights(table, &[0.0, 0.0]);
    }

    #[test]
    #[should_panic(expected = "empty table")]
    fn an_empty_table_is_refused() {
        let _ = GenotypeLikelihoodTable::from_natural_logs(&[], 3);
    }

    /// The panic the committed version documented and no test reached: a table whose
    /// cells are scored against no genotypes at all. It used to be a zero-length uniform
    /// start reaching the climb; now the table cannot be borrowed.
    #[test]
    #[should_panic(expected = "empty genotype set")]
    fn a_table_over_no_genotypes_at_all_is_refused() {
        let _ = GenotypeLikelihoodTable::from_natural_logs(&[0.5_f64.ln()], 0);
    }

    /// A start on a face of the simplex pins that genotype at zero for the whole
    /// climb — the expectation step multiplies by it — so it is refused rather than
    /// silently returning a fit over fewer genotypes than it was asked for.
    #[test]
    #[should_panic(expected = "genotype 1 starts at 0, which is not inside the simplex")]
    fn a_start_on_a_face_of_the_simplex_is_refused() {
        let ln_likelihood = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(
            diploid_table(&ln_likelihood),
            &cell_weights,
            &[0.5, 0.0, 0.5],
        );
    }

    #[test]
    #[should_panic(expected = "the start sums to 1.5 rather than one")]
    fn a_start_that_is_not_a_distribution_is_refused() {
        let ln_likelihood = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(
            diploid_table(&ln_likelihood),
            &cell_weights,
            &[0.5, 0.5, 0.5],
        );
    }

    /// A start of the wrong width is a caller who thinks the table is a different shape
    /// than it is. Refused where the mismatch is legible, rather than left to the
    /// `zip`s in the climb, which would silently fit over the shorter of the two.
    #[test]
    #[should_panic(expected = "the start lists 2 weights where the table has 3 genotypes")]
    fn a_start_of_the_wrong_width_is_refused() {
        let ln_likelihood = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(diploid_table(&ln_likelihood), &cell_weights, &[0.5, 0.5]);
    }

    /// `+∞` in the start is the mirror of `+∞` in the cell weights: it clears the
    /// `> 0.0` bound and would make the sum check pass or fail by accident.
    #[test]
    #[should_panic(expected = "genotype 0 starts at inf, which is not inside the simplex")]
    fn an_infinite_start_weight_is_refused() {
        let ln_likelihood = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(
            diploid_table(&ln_likelihood),
            &cell_weights,
            &[f64::INFINITY, 0.5, 0.5],
        );
    }

    /// **The shape the profile scan will use, compiled.** One buffer, allocated before
    /// the ladder and refilled at every rung; the table borrows it inside the loop and
    /// gives it back at the end of the iteration, so no row index is built and nothing
    /// is allocated per rung. The committed `&[&[f64]]` cannot do this — a
    /// `Vec<&[f64]>` collected before the loop makes the refill a borrow error
    /// (E0502), so it has to be rebuilt inside it, once per rung.
    #[test]
    fn a_rung_loop_refills_one_buffer_and_allocates_nothing_per_rung() {
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let mut ln_likelihood = vec![0.0; CELL_GIVEN_GENOTYPE.len() * TRUTH.len()];

        let mut fitted_per_rung = Vec::new();
        for rung in [1.0, 0.5, 0.25] {
            for (slot, &likelihood) in ln_likelihood
                .iter_mut()
                .zip(CELL_GIVEN_GENOTYPE.as_flattened())
            {
                *slot = (likelihood * rung).ln();
            }
            let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood, 3);
            fitted_per_rung.push(fit_mixture_weights(table, &cell_weights));
        }

        // Scaling every likelihood by the same factor cannot move the maximiser, so all
        // three rungs land on the truth — which is also what says the refill took.
        assert_eq!(fitted_per_rung.len(), 3);
        for fitted in &fitted_per_rung {
            for (genotype, (&got, &want)) in fitted.iter().zip(&TRUTH).enumerate() {
                assert!((got - want).abs() < 1e-9, "genotype {genotype}: {got}");
            }
        }
    }

    /// `ln_sum_exp` is the one place the table's dynamic range is handled, and a table
    /// of cells at depth 124 is far outside what linear space holds. Its answer is
    /// checked against the shift-free arithmetic where that arithmetic still works, and
    /// its behaviour on the range where it does not is stated.
    #[test]
    fn the_log_sum_is_exact_where_the_linear_sum_still_works_and_survives_where_it_does_not() {
        let small = [0.2_f64.ln(), 0.3_f64.ln(), 0.5_f64.ln()];
        assert!((ln_sum_exp(&small) - 1.0_f64.ln()).abs() < 1e-15);

        // Three terms a linear sum would flush to zero, differing by known amounts:
        // the answer is the largest plus `ln(1 + e^-10 + e^-20)`.
        let tiny = [-800.0, -810.0, -820.0];
        let want = -800.0 + (1.0 + (-10.0_f64).exp() + (-20.0_f64).exp()).ln();
        assert!((ln_sum_exp(&tiny) - want).abs() < 1e-12);

        assert_eq!(
            ln_sum_exp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]),
            f64::NEG_INFINITY
        );
        assert_eq!(ln_sum_exp(&[f64::NEG_INFINITY, 0.0]), 0.0);
    }
}
