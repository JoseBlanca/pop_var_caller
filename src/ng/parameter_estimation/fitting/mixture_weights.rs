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
//! on the whole-sample table. The surface is concave, so the climb cannot get stuck
//! on a false summit and a failure to converge is a bug rather than a data condition
//! (`spec/parameter_prepass.md` §3.1).
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

/// How many passes the climb is allowed before it gives up.
///
/// **Generous rather than tight, because exhausting it is not treated as an error.**
/// Expectation-maximization converges linearly at a rate set by how much the mixture
/// components overlap (Redner & Walker, 1984), and they overlap most at low coverage,
/// where a heterozygous site and a homozygous-reference one look alike — which is the
/// regime this estimator was built for. A cap that fired on a slow-but-correct climb
/// would turn shallow data into a crash; the tests assert convergence on tables whose
/// answer is known instead (`spec/parameter_prepass.md` §3.1).
const MAX_CLIMB_PASSES: u32 = 1_000;

/// How still the weights have to be before the climb is called finished: the largest
/// move any one weight made in a pass.
///
/// Absolute rather than relative, and it has to be read against what the weights are.
/// The two that matter here — heterozygosity and the homozygous-non-reference rate —
/// sit near 0.001, so this is a relative stillness of about one part in 10¹⁰, against
/// the 1% the recovery tests of Milestone F assert.
const CLIMB_STILLNESS: f64 = 1e-13;

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
    /// One weight per genotype, non-negative and summing to one.
    pub weights: SmallVec<[f64; 3]>,
    /// The weighted log-likelihood **at `weights`** — `Σ_cells w · ln Σ_genotypes π·L`.
    /// Computed in a final pass after the climb stopped, so it belongs to the weights
    /// returned beside it rather than to the pass before them.
    pub log_likelihood: LogProb,
    /// How many passes ran. One pass is one expectation step and one maximization step.
    pub passes: u32,
    /// Whether the weights stopped moving, as against running out of passes.
    pub converged: bool,
}

/// Fit the mixing weights of a finite mixture whose component likelihoods are already
/// known and fixed. Returns one weight per genotype, non-negative and summing to one.
///
/// Two consumers, both fitting a free point on the simplex: the error-rate scan, which
/// calls this at every rung of its ladder, and the sample's own rates, climbed once on
/// the whole-sample table (`arch/parameter_prepass_generic.md` §5.1, §5.2).
///
/// `ln_likelihood_by_cell_and_genotype` is **one row per cell, one entry per genotype
/// within a row** — natural logarithms, because a cell at depth 124 has a likelihood no
/// `f64` can hold in linear space. Every row must be the same width, and that width is
/// the number of genotypes, so it is the length of the answer: three for a diploid,
/// five for a tetraploid. `−∞` is a legal entry and means this genotype cannot have
/// produced this cell; it is the only non-finite value accepted.
///
/// `cell_weights` is how much each cell counts for — its site count, or its site count
/// times a posterior. One entry per row above.
///
/// The climb starts from the uniform point. Where it starts does not matter: with the
/// component likelihoods fixed the surface is concave over the simplex, so there is no
/// local summit that is not also the global one (`spec/parameter_prepass.md` §3.1).
///
/// # Panics
///
/// If the two arguments disagree in length, if the table is empty or ragged, if a
/// likelihood is `NaN` or `+∞`, if a weight is negative or not finite, if every weight
/// is zero, or if a cell that carries weight is one no genotype could have produced.
/// Each of those is a caller error that would otherwise leave a `NaN` to propagate
/// through the fit as a plausible number.
#[must_use]
pub fn fit_mixture_weights(
    ln_likelihood_by_cell_and_genotype: &[&[f64]],
    cell_weights: &[f64],
) -> SmallVec<[f64; 3]> {
    // The width is read off the table and the checking is left to the climb, so that
    // an empty or ragged table is reported by the one function that knows what is
    // wrong with it rather than twice, differently.
    let genotypes = ln_likelihood_by_cell_and_genotype
        .first()
        .map_or(0, |row| row.len());
    let uniform = vec![1.0 / genotypes as f64; genotypes];
    climb_mixture_weights(ln_likelihood_by_cell_and_genotype, cell_weights, &uniform).weights
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
pub(crate) fn climb_mixture_weights(
    ln_likelihood_by_cell_and_genotype: &[&[f64]],
    cell_weights: &[f64],
    start: &[f64],
) -> MixtureWeightsFit {
    let genotypes = check_table(ln_likelihood_by_cell_and_genotype, cell_weights);
    check_start(start, genotypes);

    let total_weight: f64 = cell_weights.iter().sum();

    let mut weights: SmallVec<[f64; 3]> = SmallVec::from_slice(start);
    let mut next: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, genotypes);
    // Scratch, loaded and cleared once per cell rather than allocated per cell: the
    // inner loop runs over every cell of the table, 161 times per fit.
    let mut ln_joint: SmallVec<[f64; 3]> = SmallVec::from_elem(0.0, genotypes);

    let mut passes = 0;
    let mut converged = false;
    let mut previous_score = f64::NEG_INFINITY;

    while passes < MAX_CLIMB_PASSES {
        next.fill(0.0);
        let mut score = 0.0;

        for (cell, (&row, &cell_weight)) in ln_likelihood_by_cell_and_genotype
            .iter()
            .zip(cell_weights)
            .enumerate()
        {
            if cell_weight == 0.0 {
                continue;
            }
            for (slot, (&weight, &ln_likelihood)) in
                ln_joint.iter_mut().zip(weights.iter().zip(row))
            {
                *slot = weight.ln() + ln_likelihood;
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
            for (slot, &ln_j) in next.iter_mut().zip(ln_joint.iter()) {
                *slot += cell_weight * (ln_j - ln_total).exp();
            }
        }

        // Monotone ascent is the one thing expectation-maximization is actually proved
        // to give (Dempster, Laird & Rubin 1977), so a pass that lost ground is a bug
        // in the two steps above rather than a hard table. Slack for the summation's
        // own rounding, which is a sum over hundreds of cells.
        debug_assert!(
            score >= previous_score - 1e-9 * previous_score.abs(),
            "a pass lost ground: {previous_score} → {score}"
        );
        previous_score = score;

        let mut moved: f64 = 0.0;
        for (slot, &weight) in next.iter_mut().zip(weights.iter()) {
            *slot /= total_weight;
            moved = moved.max((*slot - weight).abs());
        }
        weights.copy_from_slice(&next);
        passes += 1;

        if moved < CLIMB_STILLNESS {
            converged = true;
            break;
        }
    }

    let log_likelihood = weighted_log_likelihood(
        ln_likelihood_by_cell_and_genotype,
        cell_weights,
        &weights,
        &mut ln_joint,
    );

    MixtureWeightsFit {
        weights,
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
    ln_likelihood_by_cell_and_genotype: &[&[f64]],
    cell_weights: &[f64],
    weights: &[f64],
    ln_joint: &mut [f64],
) -> f64 {
    let mut score = 0.0;
    for (&row, &cell_weight) in ln_likelihood_by_cell_and_genotype.iter().zip(cell_weights) {
        if cell_weight == 0.0 {
            continue;
        }
        for (slot, (&weight, &ln_likelihood)) in ln_joint.iter_mut().zip(weights.iter().zip(row)) {
            *slot = weight.ln() + ln_likelihood;
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
    let largest = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if largest == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = terms.iter().map(|&term| (term - largest).exp()).sum();
    largest + sum.ln()
}

/// Check the table and return the number of genotypes — the width every row shares.
///
/// **The width is derived from the table rather than taken as an argument**, which is
/// what stops a cell count and a genotype count, two `usize`s in the same call, from
/// being swapped without a compiler or a test noticing.
fn check_table(ln_likelihood_by_cell_and_genotype: &[&[f64]], cell_weights: &[f64]) -> usize {
    assert_eq!(
        ln_likelihood_by_cell_and_genotype.len(),
        cell_weights.len(),
        "one weight per cell: {} cells against {} weights",
        ln_likelihood_by_cell_and_genotype.len(),
        cell_weights.len()
    );
    let first = ln_likelihood_by_cell_and_genotype
        .first()
        .expect("a mixture cannot be fitted from an empty table");
    let genotypes = first.len();
    assert!(
        genotypes > 0,
        "a mixture cannot be fitted over an empty genotype set"
    );

    for (cell, &row) in ln_likelihood_by_cell_and_genotype.iter().enumerate() {
        assert_eq!(
            row.len(),
            genotypes,
            "cell {cell} lists {} genotypes where cell 0 lists {genotypes}",
            row.len()
        );
        for (genotype, &ln_likelihood) in row.iter().enumerate() {
            assert!(
                ln_likelihood.is_finite() || ln_likelihood == f64::NEG_INFINITY,
                "cell {cell}, genotype {genotype}: {ln_likelihood} is not a log-likelihood"
            );
        }
    }

    let mut total_weight = 0.0;
    for (cell, &cell_weight) in cell_weights.iter().enumerate() {
        assert!(
            cell_weight.is_finite() && cell_weight >= 0.0,
            "cell {cell} carries weight {cell_weight}, which is not a count of sites"
        );
        total_weight += cell_weight;
    }
    assert!(
        total_weight > 0.0,
        "every cell carries zero weight, so there is nothing to fit"
    );

    genotypes
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
    /// (`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §2).
    ///
    /// **Four cells and three genotypes, deliberately unequal.** The two dimensions of
    /// this table are the pair a transposed index would swap, and at 4 × 3 a
    /// transposition is a length mismatch rather than a wrong number.
    const CELL_GIVEN_GENOTYPE: [[f64; 3]; 4] = [
        // homozygous reference, heterozygous, homozygous non-reference
        [0.70, 0.10, 0.02],
        [0.20, 0.40, 0.08],
        [0.07, 0.35, 0.30],
        [0.03, 0.15, 0.60],
    ];

    /// The rows of [`CELL_GIVEN_GENOTYPE`] in logs, as the climb takes them.
    fn ln_table() -> Vec<[f64; 3]> {
        CELL_GIVEN_GENOTYPE
            .iter()
            .map(|row| [row[0].ln(), row[1].ln(), row[2].ln()])
            .collect()
    }

    /// Borrow a table of fixed-width rows as the `&[&[f64]]` the climb takes.
    fn rows<const GENOTYPES: usize>(table: &[[f64; GENOTYPES]]) -> Vec<&[f64]> {
        table.iter().map(|row| row.as_slice()).collect()
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
        let table = ln_table();
        let fitted = fit_mixture_weights(&rows(&table), &weights_under(&TRUTH, 1_000_000.0));

        assert_eq!(fitted.len(), 3);
        for (genotype, (&got, &want)) in fitted.iter().zip(&TRUTH).enumerate() {
            assert!(
                (got - want).abs() < 1e-9,
                "genotype {genotype}: fitted {got}, truth {want}"
            );
        }
    }

    /// Concavity is the property that makes the start irrelevant, and it is asserted
    /// rather than assumed: three starts spanning the simplex — one at each corner's
    /// neighbourhood and one uniform — reach the same point.
    #[test]
    fn every_interior_start_reaches_the_same_summit() {
        let table = ln_table();
        let borrowed = rows(&table);
        let cell_weights = weights_under(&TRUTH, 1_000_000.0);

        let starts = [
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [0.98, 0.01, 0.01],
            [0.01, 0.01, 0.98],
            [0.20, 0.30, 0.50],
        ];
        for start in starts {
            let fit = climb_mixture_weights(&borrowed, &cell_weights, &start);
            assert!(fit.converged, "start {start:?} ran out of passes");
            assert!(
                fit.passes < MAX_CLIMB_PASSES,
                "start {start:?} used every pass"
            );
            for (genotype, (&got, &want)) in fit.weights.iter().zip(&TRUTH).enumerate() {
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
        let table = ln_table();
        let borrowed = rows(&table);

        let inbred = [0.80, 0.02, 0.18];
        let outcrossing = [0.60, 0.35, 0.05];

        let first = fit_mixture_weights(&borrowed, &weights_under(&inbred, 500_000.0));
        let second = fit_mixture_weights(&borrowed, &weights_under(&outcrossing, 500_000.0));

        for (genotype, (&got, &want)) in first.iter().zip(&inbred).enumerate() {
            assert!((got - want).abs() < 1e-9, "genotype {genotype}: {got}");
        }
        for (genotype, (&got, &want)) in second.iter().zip(&outcrossing).enumerate() {
            assert!((got - want).abs() < 1e-9, "genotype {genotype}: {got}");
        }
    }

    /// The answer really is the summit, checked against the objective rather than
    /// against the construction that produced the fixture: every step away from it
    /// along an edge of the simplex scores lower. This is the assertion that survives a
    /// change of fixture.
    #[test]
    fn no_move_away_from_the_fitted_weights_scores_higher() {
        let table = ln_table();
        let borrowed = rows(&table);
        let cell_weights = weights_under(&TRUTH, 1_000_000.0);

        let fit =
            climb_mixture_weights(&borrowed, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        let mut scratch = vec![0.0; 3];

        for from in 0..3 {
            for to in 0..3 {
                if from == to {
                    continue;
                }
                for step in [1e-4, 1e-3, 1e-2] {
                    let mut moved: Vec<f64> = fit.weights.to_vec();
                    if moved[from] <= step {
                        continue;
                    }
                    moved[from] -= step;
                    moved[to] += step;
                    let score =
                        weighted_log_likelihood(&borrowed, &cell_weights, &moved, &mut scratch);
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
        let table = ln_table();
        let borrowed = rows(&table);
        let cell_weights = weights_under(&TRUTH, 1_000.0);

        let fit =
            climb_mixture_weights(&borrowed, &cell_weights, &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);

        let mut longhand = 0.0;
        for (row, &cell_weight) in CELL_GIVEN_GENOTYPE.iter().zip(&cell_weights) {
            let mixed: f64 = row
                .iter()
                .zip(&fit.weights)
                .map(|(&likelihood, &weight)| weight * likelihood)
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
        let table = ln_table();
        let cell_weights = weights_under(&TRUTH, 7.0);
        let fitted = fit_mixture_weights(&rows(&table), &cell_weights);

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

        let table: Vec<[f64; 5]> = cell_given_dosage
            .iter()
            .map(|row| {
                let mut logs = [0.0; 5];
                for (slot, &likelihood) in logs.iter_mut().zip(row) {
                    *slot = likelihood.ln();
                }
                logs
            })
            .collect();
        let borrowed: Vec<&[f64]> = table.iter().map(|row| row.as_slice()).collect();
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
        let fit = climb_mixture_weights(&borrowed, &cell_weights, &uniform);
        assert!(fit.converged, "the climb ran out of passes");
        let fitted = fit.weights;

        assert_eq!(fitted.len(), 5);
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
        let table = [
            [0.0_f64.ln(), 0.5_f64.ln(), 0.9_f64.ln()],
            [1.0_f64.ln(), 0.5_f64.ln(), 0.1_f64.ln()],
        ];
        let borrowed = rows(&table);
        let fitted = fit_mixture_weights(&borrowed, &[100.0, 900.0]);

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
        let table = [
            [0.5_f64.ln(), 0.5_f64.ln()],
            [f64::NEG_INFINITY, f64::NEG_INFINITY],
        ];
        let _ = fit_mixture_weights(&rows(&table), &[100.0, 900.0]);
    }

    /// A cell that carries no weight is skipped, so it may say anything at all —
    /// including that no genotype produced it.
    #[test]
    fn a_cell_carrying_no_weight_is_skipped_whatever_it_says() {
        let table = [
            [0.5_f64.ln(), 0.5_f64.ln()],
            [f64::NEG_INFINITY, f64::NEG_INFINITY],
        ];
        let fitted = fit_mixture_weights(&rows(&table), &[100.0, 0.0]);

        let total: f64 = fitted.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "the weights sum to {total}");
    }

    #[test]
    #[should_panic(expected = "one weight per cell: 2 cells against 3 weights")]
    fn a_weight_without_a_cell_is_refused() {
        let table = [[0.5_f64.ln(), 0.5_f64.ln()], [0.5_f64.ln(), 0.5_f64.ln()]];
        let _ = fit_mixture_weights(&rows(&table), &[1.0, 1.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "cell 1 lists 3 genotypes where cell 0 lists 2")]
    fn a_ragged_table_is_refused() {
        let first = [0.5_f64.ln(), 0.5_f64.ln()];
        let second = [0.3_f64.ln(), 0.3_f64.ln(), 0.4_f64.ln()];
        let borrowed: Vec<&[f64]> = vec![first.as_slice(), second.as_slice()];
        let _ = fit_mixture_weights(&borrowed, &[1.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "is not a log-likelihood")]
    fn a_nan_likelihood_is_refused() {
        let table = [[f64::NAN, 0.5_f64.ln()], [0.5_f64.ln(), 0.5_f64.ln()]];
        let _ = fit_mixture_weights(&rows(&table), &[1.0, 1.0]);
    }

    /// `+∞` is not a log-likelihood either, and it is the one that would survive
    /// `ln_sum_exp` as a `NaN` rather than as an obvious wrong answer.
    #[test]
    #[should_panic(expected = "is not a log-likelihood")]
    fn a_positive_infinity_likelihood_is_refused() {
        let table = [[f64::INFINITY, 0.5_f64.ln()], [0.5_f64.ln(), 0.5_f64.ln()]];
        let _ = fit_mixture_weights(&rows(&table), &[1.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "cell 0 carries weight -1, which is not a count of sites")]
    fn a_negative_cell_weight_is_refused() {
        let table = [[0.5_f64.ln(), 0.5_f64.ln()], [0.5_f64.ln(), 0.5_f64.ln()]];
        let _ = fit_mixture_weights(&rows(&table), &[-1.0, 2.0]);
    }

    #[test]
    #[should_panic(expected = "every cell carries zero weight")]
    fn a_table_of_empty_cells_is_refused() {
        let table = [[0.5_f64.ln(), 0.5_f64.ln()], [0.5_f64.ln(), 0.5_f64.ln()]];
        let _ = fit_mixture_weights(&rows(&table), &[0.0, 0.0]);
    }

    #[test]
    #[should_panic(expected = "empty table")]
    fn an_empty_table_is_refused() {
        let _ = fit_mixture_weights(&[], &[]);
    }

    /// A start on a face of the simplex pins that genotype at zero for the whole
    /// climb — the expectation step multiplies by it — so it is refused rather than
    /// silently returning a fit over fewer genotypes than it was asked for.
    #[test]
    #[should_panic(expected = "genotype 1 starts at 0, which is not inside the simplex")]
    fn a_start_on_a_face_of_the_simplex_is_refused() {
        let table = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(&rows(&table), &cell_weights, &[0.5, 0.0, 0.5]);
    }

    #[test]
    #[should_panic(expected = "the start sums to 1.5 rather than one")]
    fn a_start_that_is_not_a_distribution_is_refused() {
        let table = ln_table();
        let cell_weights = weights_under(&TRUTH, 1_000.0);
        let _ = climb_mixture_weights(&rows(&table), &cell_weights, &[0.5, 0.5, 0.5]);
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
