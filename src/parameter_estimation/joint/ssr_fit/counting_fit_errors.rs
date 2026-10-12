//! **Standard errors for the counting fit** ([`composite_fit`](super::composite_fit)), two ways:
//! resampling the class's tracts, and the sandwich formula.
//!
//! **Why the plain curvature is not enough.** The counting fit scores each sample at a tract as if
//! its alleles were drawn from the class's allele-length variation independently of the other
//! samples there. They are not: samples at one tract share that tract's allele frequencies. The
//! answer is still right on average, but the objective's curvature counts every sample as fresh
//! evidence, so an error read off it alone comes out too small — by more the more samples a tract
//! has and the more tracts differ from one another. Both methods here treat the **tract** as the
//! unit that repeats, which is what is independent.
//!
//! - **The bootstrap** ([`bootstrap_errors`]) draws the class's tracts again with replacement, as
//!   many as it has, refits each draw from the full fit's answer, and reports the spread of the
//!   refitted numbers. It costs one refit a replicate.
//! - **The sandwich** ([`sandwich_errors`]) is the formula a composite likelihood takes (Godambe's):
//!   the inverse of the objective's curvature, times how much the objective's slope varies from
//!   tract to tract, times the inverse again. Where the objective were the true likelihood the
//!   middle term would equal the curvature and the formula would reduce to the plain inverse
//!   curvature, which is also returned, as `naive`, to show the gap.
//!
//! **How the slopes are computed.** Each tract's slope comes from the expected counts the fit's own
//! rounds gather (Fisher's identity: the slope of the log of a sum over hidden genotypes and slips
//! is the expected slope of the log-likelihood that knew them), in one pass over the reads. The
//! slope with respect to the pairwise likelihood's concentration is taken numerically, tract by
//! tract, because the pairwise likelihood has no hidden counts to read it from. The curvature is
//! the change in the summed slope over small steps of each number in turn, two passes a number.
//!
//! **The concentration is estimated by a second objective** when it comes from pairs of samples,
//! so the sandwich is the general one for a set of estimating equations: the matrix in the middle
//! of the curvature is not symmetric, and the formula is `A⁻¹ B A⁻ᵀ` with `A` the slopes' change
//! and `B` the tract-to-tract spread of the slopes.

use super::composite_fit::{
    ConcentrationSource, CountedClass, Counts, PairRow, expected_counts_by_tract,
    pair_rows_by_tract, refit_class_by_counting_from, scaled_ln_pairs_of_tract,
};
use super::*;

/// **The standard errors of one class's counting fit**, each on its number's own scale.
#[derive(Debug, Clone, PartialEq)]
pub struct CountedClassErrors {
    pub level: StratumError,
    pub shorter_share: StratumError,
    pub fall_off: StratumError,
    /// Per allele class, indexed as [`CountedClass::allele_length_variation`] is.
    pub allele_length_variation: Vec<StratumError>,
    /// The Dirichlet concentration's error.
    pub concentration: StratumError,
    /// The between-locus share `F = 1 / (1 + concentration)`'s error.
    pub between_locus_share: StratumError,
}

/// What [`sandwich_errors`] returns: the sandwich, and the plain inverse curvature for comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvatureErrors {
    /// Errors that allow for samples at one tract sharing its alleles.
    pub sandwich: CountedClassErrors,
    /// **Errors as if the objective were the true likelihood**: the inverse of its curvature alone.
    /// Too small wherever samples at a tract share alleles; returned to show by how much. With the
    /// concentration from pairs, the stutter and shares come from the counting objective's
    /// curvature and the concentration from the pairwise objective's.
    pub naive: CountedClassErrors,
}

/// **How a bootstrap is run**: how many replicates, the seed of the draws, and the refit's own
/// settings — the same rounds, tolerance and alternations the full fit was run with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BootstrapSettings {
    pub replicates: usize,
    pub seed: u64,
    pub max_rounds: usize,
    pub tolerance: f64,
    pub alternations: usize,
}

/// A share below this is treated as absent: held at its value, not fitted, and given no error.
const SMALLEST_SHARE: f64 = 1e-10;

/// The step, on the log scale, the pairwise likelihood's slope in the concentration is taken over.
const LN_CONCENTRATION_STEP: f64 = 1e-3;

/// **The numbers whose errors are computed, in the order the matrices hold them**: the three
/// stutter numbers, every share but one (the largest, which is one minus the rest), and the log of
/// the concentration when it is identified.
struct Layout {
    classes: usize,
    /// The share written as one minus the others.
    pivot: usize,
    /// The shares fitted, in order: every class but the pivot and the absent ones.
    free_shares: Vec<usize>,
    with_concentration: bool,
    /// Which of the three stutter numbers sit at an end of their range: held there, with no slope
    /// and no error, since the formula needs room either side.
    held_stutter: [bool; 3],
}

impl Layout {
    fn width(&self) -> usize {
        3 + self.free_shares.len() + usize::from(self.with_concentration)
    }

    fn concentration_index(&self) -> usize {
        3 + self.free_shares.len()
    }
}

/// One point the objective is evaluated at.
#[derive(Debug, Clone)]
struct Point {
    slippage: Slippage,
    shares: Vec<f64>,
    concentration: f64,
}

impl Point {
    /// The point moved by `step` along coordinate `which` of `layout`.
    fn moved(&self, layout: &Layout, which: usize, step: f64) -> Point {
        let mut moved = self.clone();
        match which {
            0 => moved.slippage.level += step,
            1 => moved.slippage.shorter_share += step,
            2 => moved.slippage.fall_off += step,
            _ if layout.with_concentration && which == layout.concentration_index() => {
                moved.concentration *= float::exp(step);
            }
            share => {
                let class = layout.free_shares[share - 3];
                moved.shares[class] += step;
                moved.shares[layout.pivot] -= step;
            }
        }
        moved
    }

    /// How far coordinate `which` may be stepped either side and stay well inside its range.
    fn step_for(&self, layout: &Layout, which: usize) -> f64 {
        let inside = |value: f64| 1e-4 * value.min(1.0 - value).max(1e-9);
        match which {
            0 => inside(self.slippage.level),
            1 => inside(self.slippage.shorter_share),
            2 => inside(self.slippage.fall_off),
            _ if layout.with_concentration && which == layout.concentration_index() => {
                LN_CONCENTRATION_STEP
            }
            share => {
                let class = layout.free_shares[share - 3];
                1e-4 * self.shares[class].min(self.shares[layout.pivot])
            }
        }
    }
}

/// **The counting objective's slope at `point`, from its expected counts there** (Fisher's
/// identity), laid out as `layout` says. The concentration's entry is filled only when it comes
/// from each sample's two copies; from pairs it is the pairwise likelihood's and is filled apart.
///
/// - slip rate `L`: each read reported its allele (chance `1 − L`) or slipped (`L`):
///   `slipped / L − exact / (1 − L)`;
/// - shorter share `q`: `shorter / q − longer / (1 − q)`;
/// - fall-off `f`: a slip of `k` steps has chance `(1 − f) f^(k − 1)`:
///   `(steps − stepped) / f − stepped / (1 − f)`;
/// - share `s_k`, the pivot `s_p` being one minus the rest: `draws_k / s_k − draws_p / s_p`;
/// - log concentration `ln c`: of a sample's two copies not identical by descent, a share
///   `1 / (1 + c)` drew one chromosome lineage of the tract: `pairs / (1 + c) − shared`.
fn counting_slope(
    counts: &Counts,
    point: &Point,
    layout: &Layout,
    source: ConcentrationSource,
) -> Vec<f64> {
    let mut slope = vec![0.0; layout.width()];
    let Slippage {
        level,
        shorter_share,
        fall_off,
    } = point.slippage;
    slope[0] = counts.slipped / level - counts.exact / (1.0 - level);
    slope[1] = counts.shorter / shorter_share - counts.longer / (1.0 - shorter_share);
    slope[2] = (counts.steps - counts.stepped) / fall_off - counts.stepped / (1.0 - fall_off);
    for (index, held) in layout.held_stutter.iter().enumerate() {
        if *held {
            slope[index] = 0.0;
        }
    }
    let pivot = counts.draws[layout.pivot] / point.shares[layout.pivot];
    for (index, class) in layout.free_shares.iter().enumerate() {
        slope[3 + index] = counts.draws[*class] / point.shares[*class] - pivot;
    }
    if layout.with_concentration && source == ConcentrationSource::EachSamplesTwoCopies {
        slope[layout.concentration_index()] =
            counts.pairs_not_by_descent / (1.0 + point.concentration) - counts.shared_in_tract;
    }
    slope
}

/// **Each tract's pairwise slope in the log concentration**, by a central difference.
fn pairwise_slope_by_tract(rows: &[Vec<PairRow>], point: &Point, classes: usize) -> Vec<f64> {
    let up = point.concentration * float::exp(LN_CONCENTRATION_STEP);
    let down = point.concentration * float::exp(-LN_CONCENTRATION_STEP);
    rows.par_iter()
        .map(|samples| {
            (scaled_ln_pairs_of_tract(samples, &point.shares, up, classes)
                - scaled_ln_pairs_of_tract(samples, &point.shares, down, classes))
                / (2.0 * LN_CONCENTRATION_STEP)
        })
        .collect()
}

/// Every tract's slope at `point`, one row a tract laid out as `layout` says.
fn slopes_by_tract(
    evidence: &StratumEvidence,
    point: &Point,
    layout: &Layout,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    rows: Option<&[Vec<PairRow>]>,
) -> Vec<Vec<f64>> {
    let counts = expected_counts_by_tract(
        evidence,
        point.slippage,
        &point.shares,
        point.concentration,
        homozygote_excess,
        config,
    );
    let mut slopes: Vec<Vec<f64>> = counts
        .iter()
        .map(|tract| counting_slope(tract, point, layout, source))
        .collect();
    if layout.with_concentration && source == ConcentrationSource::PairsOfSamples {
        let rows = rows.expect("the pairs' rows are built when the concentration comes from pairs");
        let at = layout.concentration_index();
        for (slope, pairwise) in
            slopes
                .iter_mut()
                .zip(pairwise_slope_by_tract(rows, point, layout.classes))
        {
            slope[at] = pairwise;
        }
    }
    slopes
}

/// The slope summed over every tract.
fn summed(slopes: &[Vec<f64>], width: usize) -> Vec<f64> {
    let mut total = vec![0.0; width];
    for slope in slopes {
        for (into, value) in total.iter_mut().zip(slope) {
            *into += value;
        }
    }
    total
}

/// **The inverse of a square matrix** by Gauss–Jordan elimination with partial pivoting, row-major,
/// or `None` when a pivot falls below `1e-12` of the largest entry (singular to working precision).
fn inverse(matrix: &[f64], size: usize) -> Option<Vec<f64>> {
    let scale = matrix
        .iter()
        .fold(0.0_f64, |most, value| most.max(value.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return None;
    }
    let mut left = matrix.to_vec();
    let mut right = vec![0.0; size * size];
    for index in 0..size {
        right[index * size + index] = 1.0;
    }
    for column in 0..size {
        let pivot_row = (column..size)
            .max_by(|a, b| {
                left[a * size + column]
                    .abs()
                    .total_cmp(&left[b * size + column].abs())
            })
            .expect("the range is not empty");
        if left[pivot_row * size + column].abs() < 1e-12 * scale {
            return None;
        }
        if pivot_row != column {
            for position in 0..size {
                left.swap(pivot_row * size + position, column * size + position);
                right.swap(pivot_row * size + position, column * size + position);
            }
        }
        let pivot = left[column * size + column];
        for position in 0..size {
            left[column * size + position] /= pivot;
            right[column * size + position] /= pivot;
        }
        for row in 0..size {
            if row == column {
                continue;
            }
            let factor = left[row * size + column];
            if factor == 0.0 {
                continue;
            }
            for position in 0..size {
                left[row * size + position] -= factor * left[column * size + position];
                right[row * size + position] -= factor * right[column * size + position];
            }
        }
    }
    Some(right)
}

/// `first · second`, both `size × size`, row-major.
fn product(first: &[f64], second: &[f64], size: usize) -> Vec<f64> {
    let mut out = vec![0.0; size * size];
    for row in 0..size {
        for middle in 0..size {
            let left = first[row * size + middle];
            if left == 0.0 {
                continue;
            }
            for column in 0..size {
                out[row * size + column] += left * second[middle * size + column];
            }
        }
    }
    out
}

/// The transpose of a `size × size` matrix.
fn transpose(matrix: &[f64], size: usize) -> Vec<f64> {
    let mut out = vec![0.0; size * size];
    for row in 0..size {
        for column in 0..size {
            out[column * size + row] = matrix[row * size + column];
        }
    }
    out
}

/// One number's error from its variance: none when the variance is not a positive number.
fn error_of(variance: f64) -> StratumError {
    if variance.is_finite() && variance > 0.0 {
        StratumError::Estimated(variance.sqrt())
    } else {
        StratumError::NotIdentified
    }
}

/// **A covariance over `layout`'s coordinates turned into each reported number's error.** The pivot
/// share's variance is that of one minus the other shares; the concentration's and the
/// between-locus share's are carried from the log concentration's.
fn errors_from_covariance(
    covariance: Option<&[f64]>,
    layout: &Layout,
    fitted: &CountedClass,
) -> CountedClassErrors {
    let size = layout.width();
    let variance = |row: usize, column: usize| -> f64 {
        covariance.map_or(f64::NAN, |matrix| matrix[row * size + column])
    };
    // A stutter number held at an end of its range is not placed by the formula.
    let stutter = |index: usize| {
        if layout.held_stutter[index] {
            StratumError::NotPlaced
        } else {
            error_of(variance(index, index))
        }
    };
    let mut shares = vec![StratumError::NoShare; layout.classes];
    for (index, class) in layout.free_shares.iter().enumerate() {
        shares[*class] = error_of(variance(3 + index, 3 + index));
    }
    let free = layout.free_shares.len();
    let pivot_variance: f64 = (0..free)
        .flat_map(|row| (0..free).map(move |column| (row, column)))
        .map(|(row, column)| variance(3 + row, 3 + column))
        .sum();
    shares[layout.pivot] = if free == 0 {
        StratumError::NotIdentified
    } else {
        error_of(pivot_variance)
    };
    let (concentration, between_locus_share) = if layout.with_concentration {
        let at = layout.concentration_index();
        match error_of(variance(at, at)) {
            StratumError::Estimated(ln_error) if ln_error > NOT_PLACED => {
                (StratumError::NotPlaced, StratumError::NotPlaced)
            }
            StratumError::Estimated(ln_error) => {
                let c = fitted.concentration;
                let share = 1.0 / (1.0 + c);
                (
                    StratumError::Estimated(c * ln_error),
                    StratumError::Estimated(share * (1.0 - share) * ln_error),
                )
            }
            other => (other, other),
        }
    } else {
        (StratumError::NotIdentified, StratumError::NotIdentified)
    };
    CountedClassErrors {
        level: stutter(0),
        shorter_share: stutter(1),
        fall_off: stutter(2),
        allele_length_variation: shares,
        concentration,
        between_locus_share,
    }
}

/// **The sandwich errors of a counting fit, and the plain inverse curvature beside them.**
///
/// `fitted` is the answer of the estimator `source` names, on `evidence`, at `homozygote_excess`
/// and `config` — the same ones it was fitted with. The concentration has an error only where its
/// estimator has something to read: pairs need a tract with two samples, a sample's own two copies
/// need a homozygote excess below one. Otherwise it is [`StratumError::NotIdentified`], and the
/// other numbers' errors are those with the concentration held at its value.
///
/// **A stutter number at an end of its range** — a shorter share of exactly one when no read slipped
/// longer, say — is held there and reported as [`StratumError::NotPlaced`]: the formula needs room
/// either side of the answer. The other numbers' errors are those with it held.
///
/// Costs two passes over the reads for each number fitted, plus one, in the counting fit's own
/// per-round work; with the concentration from pairs, also the pairs' rows at seven stutters.
#[must_use]
pub fn sandwich_errors(
    evidence: &StratumEvidence,
    fitted: &CountedClass,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> CurvatureErrors {
    let classes = allele_classes_of(config);
    let pivot = fitted
        .allele_length_variation
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1))
        .map_or(0, |(index, _)| index);
    let free_shares: Vec<usize> = (0..classes)
        .filter(|class| *class != pivot && fitted.allele_length_variation[*class] > SMALLEST_SHARE)
        .collect();
    let concentration_identified =
        the_concentration_is_identified(evidence, source, homozygote_excess);
    let at_an_end = |value: f64, highest: f64| value <= 1e-6 || value >= highest - 1e-6;
    let held_stutter = [
        at_an_end(fitted.slippage.level, 1.0),
        at_an_end(fitted.slippage.shorter_share, 1.0),
        // The fit holds the fall-off at or below 0.99.
        at_an_end(fitted.slippage.fall_off, 0.99),
    ];
    let layout = Layout {
        classes,
        pivot,
        free_shares,
        with_concentration: concentration_identified,
        held_stutter,
    };
    let size = layout.width();
    let at = Point {
        slippage: fitted.slippage,
        shares: fitted.allele_length_variation.clone(),
        concentration: fitted.concentration,
    };
    let from_pairs = layout.with_concentration && source == ConcentrationSource::PairsOfSamples;
    let rows_at = |slippage: Slippage| {
        from_pairs.then(|| pair_rows_by_tract(evidence, slippage, homozygote_excess, config))
    };
    let rows = rows_at(at.slippage);
    let by_tract = slopes_by_tract(
        evidence,
        &at,
        &layout,
        source,
        homozygote_excess,
        config,
        rows.as_deref(),
    );

    // The slopes' change, column by column: minus the derivative of the summed slope.
    let mut change = vec![0.0; size * size];
    for column in 0..size {
        if column < 3 && layout.held_stutter[column] {
            // A held number is its own block: the identity, so it moves nothing else.
            change[column * size + column] = 1.0;
            continue;
        }
        let step = at.step_for(&layout, column);
        let mut sides = [0.0, 0.0].map(|_| Vec::new());
        for (side, sign) in [1.0, -1.0].into_iter().enumerate() {
            let moved = at.moved(&layout, column, sign * step);
            // The pairs' rows depend on the stutter alone.
            let moved_rows = if column < 3 {
                rows_at(moved.slippage)
            } else {
                None
            };
            let rows_here = moved_rows.as_deref().or(rows.as_deref());
            sides[side] = summed(
                &slopes_by_tract(
                    evidence,
                    &moved,
                    &layout,
                    source,
                    homozygote_excess,
                    config,
                    rows_here,
                ),
                size,
            );
        }
        for row in 0..size {
            if row < 3 && layout.held_stutter[row] {
                continue;
            }
            change[row * size + column] = -(sides[0][row] - sides[1][row]) / (2.0 * step);
        }
    }

    // How much the slope varies from tract to tract.
    let mut spread = vec![0.0; size * size];
    for slope in &by_tract {
        for row in 0..size {
            for column in 0..size {
                spread[row * size + column] += slope[row] * slope[column];
            }
        }
    }

    let sandwich = inverse(&change, size).map(|inverse_change| {
        product(
            &product(&inverse_change, &spread, size),
            &transpose(&inverse_change, size),
            size,
        )
    });

    // The plain inverse curvature: the counting objective's own block, symmetrised, and with the
    // concentration from pairs, the pairwise objective's curvature in it on its own.
    let naive = if from_pairs {
        let counted = size - 1;
        let mut block = vec![0.0; counted * counted];
        for row in 0..counted {
            for column in 0..counted {
                block[row * counted + column] =
                    0.5 * (change[row * size + column] + change[column * size + row]);
            }
        }
        inverse(&block, counted).map(|inverse_block| {
            let mut full = vec![0.0; size * size];
            for row in 0..counted {
                for column in 0..counted {
                    full[row * size + column] = inverse_block[row * counted + column];
                }
            }
            let at = size - 1;
            full[at * size + at] = 1.0 / change[at * size + at];
            full
        })
    } else {
        inverse(&transpose_mean(&change, size), size)
    };

    CurvatureErrors {
        sandwich: errors_from_covariance(sandwich.as_deref(), &layout, fitted),
        naive: errors_from_covariance(naive.as_deref(), &layout, fitted),
    }
}

/// `(M + Mᵀ) / 2`: a curvature measured column by column, made symmetric.
fn transpose_mean(matrix: &[f64], size: usize) -> Vec<f64> {
    let mut out = vec![0.0; size * size];
    for row in 0..size {
        for column in 0..size {
            out[row * size + column] =
                0.5 * (matrix[row * size + column] + matrix[column * size + row]);
        }
    }
    out
}

/// **Whether `source` has anything to read the concentration from**: pairs need a tract with two
/// samples, a sample's own two copies need a homozygote excess below one.
fn the_concentration_is_identified(
    evidence: &StratumEvidence,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
) -> bool {
    match source {
        ConcentrationSource::PairsOfSamples => {
            evidence.tracts.iter().any(|tract| tract.samples.len() >= 2)
        }
        ConcentrationSource::EachSamplesTwoCopies => evidence.tracts.iter().any(|tract| {
            tract
                .samples
                .iter()
                .any(|sample| homozygote_excess[sample.sample as usize] < 1.0)
        }),
    }
}

/// A reproducible stream for drawing tracts, one a replicate.
struct Draws(u64);

impl Draws {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A whole number below `bound`.
    fn below(&mut self, bound: usize) -> usize {
        ((u128::from(self.next_u64()) * bound as u128) >> 64) as usize
    }
}

/// **The class's tracts drawn again with replacement**, as many as it has, by the stream seeded at
/// `seed`. Everything but the tracts is the class's own.
fn resampled(evidence: &StratumEvidence, seed: u64) -> StratumEvidence {
    let mut draws = Draws(seed);
    let count = evidence.tracts.len();
    StratumEvidence {
        tracts: (0..count)
            .map(|_| evidence.tracts[draws.below(count)].clone())
            .collect(),
        ..evidence.clone_without_tracts()
    }
}

impl StratumEvidence {
    /// Everything of the class but its tracts, which are left empty.
    fn clone_without_tracts(&self) -> StratumEvidence {
        StratumEvidence {
            stratum: self.stratum,
            tracts: Vec::new(),
            read_span: self.read_span,
            groups: self.groups,
            libraries: Arc::clone(&self.libraries),
            tracts_over_guard_threshold: self.tracts_over_guard_threshold,
            reads_reaching_not_crossing: self.reads_reaching_not_crossing,
            guard_reads: self.guard_reads,
            bases_compared: self.bases_compared,
            mismatching_bases: self.mismatching_bases,
        }
    }
}

/// The sample standard deviation of `values`, or `None` with fewer than two.
fn standard_deviation(values: &[f64]) -> Option<f64> {
    if values.len() < 2 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let squares: f64 = values
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum();
    Some((squares / (values.len() - 1) as f64).sqrt())
}

/// **The bootstrap errors of a counting fit**: the class's tracts drawn again with replacement
/// `settings.replicates` times, each draw refitted from `fitted` by the estimator `source` names,
/// and each number's error the standard deviation of its refitted values.
///
/// **The concentration has no error where its estimator has nothing to read** (as in
/// [`sandwich_errors`]), whatever the replicates returned.
///
/// **The concentration's error is [`StratumError::NotPlaced`]** when the refitted log
/// concentrations spread by more than [`NOT_PLACED`]: an estimate that runs to either end of the
/// search between replicates has no error worth quoting on its own scale.
///
/// Costs one refit a replicate. A refit starts at the full fit's answer, so it takes fewer rounds
/// than a fit from flat shares; how many fewer depends on the class.
#[must_use]
pub fn bootstrap_errors(
    evidence: &StratumEvidence,
    fitted: &CountedClass,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    settings: BootstrapSettings,
) -> CountedClassErrors {
    let replicates: Vec<CountedClass> = (0..settings.replicates)
        .into_par_iter()
        .map(|replicate| {
            let drawn = resampled(
                evidence,
                settings
                    .seed
                    .wrapping_mul(0x2545_F491_4F6C_DD1D)
                    .wrapping_add(replicate as u64),
            );
            refit_class_by_counting_from(
                &drawn,
                fitted,
                source,
                homozygote_excess,
                config,
                settings.max_rounds,
                settings.tolerance,
                settings.alternations,
            )
        })
        .collect();
    let mut errors = errors_from_replicates(&replicates, fitted);
    if !the_concentration_is_identified(evidence, source, homozygote_excess) {
        errors.concentration = StratumError::NotIdentified;
        errors.between_locus_share = StratumError::NotIdentified;
    }
    errors
}

/// Each number's error as the spread of its value over `replicates`.
fn errors_from_replicates(
    replicates: &[CountedClass],
    fitted: &CountedClass,
) -> CountedClassErrors {
    let spread_of = |value: &dyn Fn(&CountedClass) -> f64| -> StratumError {
        let values: Vec<f64> = replicates.iter().map(value).collect();
        standard_deviation(&values).map_or(StratumError::NotIdentified, error_of_spread)
    };
    let ln_spread = spread_of(&|fit| float::ln(fit.concentration));
    let not_placed = matches!(ln_spread, StratumError::Estimated(error) if error > NOT_PLACED);
    let (concentration, between_locus_share) = if not_placed {
        (StratumError::NotPlaced, StratumError::NotPlaced)
    } else {
        (
            spread_of(&|fit| fit.concentration),
            spread_of(&|fit| 1.0 / (1.0 + fit.concentration)),
        )
    };
    CountedClassErrors {
        level: spread_of(&|fit| fit.slippage.level),
        shorter_share: spread_of(&|fit| fit.slippage.shorter_share),
        fall_off: spread_of(&|fit| fit.slippage.fall_off),
        allele_length_variation: (0..fitted.allele_length_variation.len())
            .map(|class| spread_of(&|fit| fit.allele_length_variation[class]))
            .collect(),
        concentration,
        between_locus_share,
    }
}

/// A spread of exactly zero says the replicates never moved the number, not that it is known
/// exactly.
fn error_of_spread(spread: f64) -> StratumError {
    if spread > 0.0 && spread.is_finite() {
        StratumError::Estimated(spread)
    } else {
        StratumError::NotIdentified
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};
    use crate::parameter_estimation::joint::ssr_fit::composite_fit::{
        fit_class_by_counting, fit_class_by_counting_with_pairs,
    };

    const TRUTH: Slippage = Slippage {
        level: 0.08,
        shorter_share: 0.8,
        fall_off: 0.25,
    };

    /// **The slope read off the expected counts is the composite log-likelihood's slope**, taken
    /// numerically, at a point away from the answer, for every coordinate.
    #[test]
    fn the_slope_from_the_expected_counts_is_the_objectives_slope() {
        let config = SsrFitConfig {
            allele_span: 2,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(TRUTH, &spectrum_of(5), 0.5, 0.4, 60, 6, 5, 2, 3);
        let excess = [0.4; 6];
        let at = Point {
            slippage: Slippage {
                level: 0.12,
                shorter_share: 0.6,
                fall_off: 0.35,
            },
            shares: vec![0.1, 0.2, 0.4, 0.2, 0.1],
            concentration: 0.8,
        };
        let layout = Layout {
            classes: 5,
            pivot: 2,
            free_shares: vec![0, 1, 3, 4],
            with_concentration: true,
            held_stutter: [false; 3],
        };
        let source = ConcentrationSource::EachSamplesTwoCopies;
        let objective = |point: &Point| -> f64 {
            expected_counts_by_tract(
                &evidence,
                point.slippage,
                &point.shares,
                point.concentration,
                &excess,
                &config,
            )
            .iter()
            .map(|counts| counts.log_likelihood)
            .sum()
        };
        let slopes = slopes_by_tract(&evidence, &at, &layout, source, &excess, &config, None);
        let analytic = summed(&slopes, layout.width());
        for (which, value) in analytic.iter().enumerate() {
            let step = 1e-5;
            let numeric = (objective(&at.moved(&layout, which, step))
                - objective(&at.moved(&layout, which, -step)))
                / (2.0 * step);
            assert!(
                (value - numeric).abs() <= 1e-4 * numeric.abs().max(1.0),
                "coordinate {which}: from the counts {value}, numerically {numeric}"
            );
        }
    }

    /// **Where samples at a tract share its alleles, the sandwich's share errors exceed the plain
    /// curvature's**, and the bootstrap's sit near the sandwich's. Twenty samples a tract at a
    /// between-locus share of two thirds: each tract carries far less about the class's shares than
    /// twenty independent samples would.
    #[test]
    fn the_sandwich_widens_the_share_errors_where_tracts_differ() {
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(TRUTH, &spectrum_of(3), 0.5, 0.4, 300, 20, 3, 1, 11);
        let excess = [0.4; 20];
        let start = Slippage {
            level: 0.05,
            shorter_share: 0.5,
            fall_off: 0.3,
        };
        let fitted =
            fit_class_by_counting_with_pairs(&evidence, start, &excess, &config, 2000, 1e-7, 10);
        let source = ConcentrationSource::PairsOfSamples;
        let curvature = sandwich_errors(&evidence, &fitted, source, &excess, &config);
        let bootstrap = bootstrap_errors(
            &evidence,
            &fitted,
            source,
            &excess,
            &config,
            BootstrapSettings {
                replicates: 60,
                seed: 5,
                max_rounds: 2000,
                tolerance: 1e-7,
                alternations: 10,
            },
        );
        let middle = |errors: &CountedClassErrors| {
            errors.allele_length_variation[1]
                .value()
                .expect("the middle share has an error")
        };
        let (sandwich, naive, resampled) = (
            middle(&curvature.sandwich),
            middle(&curvature.naive),
            middle(&bootstrap),
        );
        assert!(sandwich > 2.0 * naive, "sandwich {sandwich}, naive {naive}");
        assert!(
            (0.7..1.4).contains(&(resampled / sandwich)),
            "bootstrap {resampled}, sandwich {sandwich}"
        );
        assert!(curvature.sandwich.concentration.value().is_some());
    }

    /// **One sample has no pair, so the concentration from pairs has no error**, and the other
    /// numbers still get theirs.
    #[test]
    fn one_sample_leaves_the_pairs_concentration_without_an_error() {
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(TRUTH, &spectrum_of(3), 0.5, 0.4, 300, 1, 30, 1, 2);
        let excess = [0.4];
        let fitted = fit_class_by_counting(&evidence, TRUTH, &excess, &config, 2000, 1e-7);
        let errors = sandwich_errors(
            &evidence,
            &fitted,
            ConcentrationSource::PairsOfSamples,
            &excess,
            &config,
        );
        assert_eq!(
            errors.sandwich.concentration,
            StratumError::NotIdentified,
            "{errors:?}"
        );
        assert!(errors.sandwich.level.value().is_some(), "{errors:?}");
    }
}
