//! The integral over a position's allele frequency: a Gauss–Jacobi rule for a Beta density, and
//! the special functions the shapes' step reads.

use crate::float;
use crate::parameter_estimation::joint::linear_algebra;

use super::expectation::{genotype_frequencies, genotype_frequencies_slope_in_frequency};

/// A Gauss–Jacobi rule for `∫₀¹ Beta(f; a, b) · h(f) df`, normalised so the weights sum to one.
///
/// **Gauss–Jacobi rather than a grid**, because a Beta with a shape below one — which is what
/// a population with many rare alleles has — is unbounded at zero, and a uniform grid puts no
/// node where nearly all the mass is. The rule absorbs the singularity into its weight
/// function, so the nodes are placed by the density itself.
pub(super) struct BetaQuadrature {
    pub(super) nodes: Vec<f64>,
    /// The rule's weights, summing to one. **The pass reads `ln_weights`**; these are what the
    /// tests integrate a known Beta with, which is the only check that the rule is the right
    /// rule rather than merely a consistent one.
    #[cfg_attr(
        not(test),
        allow(dead_code, reason = "read by this module's own tests")
    )]
    pub(super) weights: Vec<f64>,
    pub(super) ln_weights: Vec<f64>,
    /// `ln f` at each node, for the test that measures the digamma form ([`Statistics`]).
    #[cfg(test)]
    pub(super) ln_nodes: Vec<f64>,
    /// `[node][sample][genotype]` — how common each genotype is at that node's frequency for
    /// that sample's own inbreeding. **Computed once for the whole pass**, because it depends
    /// on nothing that varies from position to position, and it sits in the innermost loop.
    pub(super) priors: Vec<f64>,
    /// `[node][sample][genotype]` — how that genotype frequency moves with the node's frequency
    /// (`genotype_frequencies_slope_in_frequency`), laid out as `priors`: what the likelihood's
    /// slope in the density's shapes reads as the rule's nodes move.
    pub(super) prior_slopes: Vec<f64>,
}

impl BetaQuadrature {
    /// The rule, plus the genotype priors every position will read.
    pub(super) fn with_genotype_priors(a: f64, b: f64, count: usize, hom_excess: &[f64]) -> Self {
        let mut rule = Self::new(a, b, count);
        rule.priors = Vec::with_capacity(rule.nodes.len() * hom_excess.len() * 3);
        rule.prior_slopes = Vec::with_capacity(rule.nodes.len() * hom_excess.len() * 3);
        for &f in &rule.nodes {
            for &excess in hom_excess {
                rule.priors
                    .extend_from_slice(&genotype_frequencies(f, excess));
                rule.prior_slopes
                    .extend_from_slice(&genotype_frequencies_slope_in_frequency(f, excess));
            }
        }
        rule
    }

    pub(super) fn new(a: f64, b: f64, count: usize) -> Self {
        // Mapping `f = (1 + x)/2` turns the Beta's weight into the Jacobi weight
        // `(1 − x)^(b−1) (1 + x)^(a−1)`.
        let (alpha, beta) = (b - 1.0, a - 1.0);
        let (x, w) = gauss_jacobi(alpha, beta, count);
        let nodes: Vec<f64> = x
            .iter()
            .map(|x| (0.5 * (1.0 + x)).clamp(1e-12, 1.0 - 1e-12))
            .collect();
        let total: f64 = w.iter().sum();
        let weights: Vec<f64> = w.iter().map(|w| w / total).collect();
        #[cfg(test)]
        let ln_nodes = nodes.iter().map(|f| float::ln(*f)).collect();
        let ln_weights = weights
            .iter()
            .map(|w| float::ln(w.max(f64::MIN_POSITIVE)))
            .collect();
        Self {
            nodes,
            weights,
            ln_weights,
            #[cfg(test)]
            ln_nodes,
            priors: Vec::new(),
            prior_slopes: Vec::new(),
        }
    }
}

/// Nodes and weights of the `count`-point Gauss–Jacobi rule on `[-1, 1]` with weight
/// `(1 − x)^α (1 + x)^β`, by the Golub–Welsch construction: the rule's nodes are the
/// eigenvalues of the recurrence's tridiagonal matrix and its weights come from the first
/// component of each eigenvector.
fn gauss_jacobi(alpha: f64, beta: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
    let n = count.max(2);
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        let k = i as f64;
        let s = 2.0 * k + alpha + beta;
        let diagonal = if s.abs() < 1e-12 || (s + 2.0).abs() < 1e-12 {
            (beta - alpha) / (alpha + beta + 2.0)
        } else {
            (beta * beta - alpha * alpha) / (s * (s + 2.0))
        };
        matrix[i * n + i] = diagonal;
        if i + 1 < n {
            let k = (i + 1) as f64;
            let s = 2.0 * k + alpha + beta;
            let numerator = 4.0 * k * (k + alpha) * (k + beta) * (k + alpha + beta);
            let denominator = s * s * (s + 1.0) * (s - 1.0);
            let off = if denominator.abs() < 1e-12 {
                0.0
            } else {
                (numerator / denominator).max(0.0).sqrt()
            };
            matrix[i * n + i + 1] = off;
            matrix[(i + 1) * n + i] = off;
        }
    }
    let (values, vectors) = symmetric_eigen(&matrix, n);
    // The rule's weight is the **first component** of each eigenvector, squared — and
    // `vectors` is laid out eigenvector-major, so that component is at `k * n` and not at `k`.
    let weights: Vec<f64> = (0..n).map(|k| vectors[k * n] * vectors[k * n]).collect();
    (values, weights)
}

/// Eigenvalues and eigenvectors of a small symmetric matrix, by cyclic Jacobi rotations
/// ([`linear_algebra::jacobi_rotations`]), in ascending order of eigenvalue.
///
/// `vectors[k * n + s]` is component `s` of the eigenvector for eigenvalue `k`; the rule above
/// reads only the first component of each, which is why they are returned that way round.
fn symmetric_eigen(matrix: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let (diagonal, v) = linear_algebra::jacobi_rotations(matrix, n, 1e-22, 1e-18);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        diagonal[i]
            .partial_cmp(&diagonal[j])
            .expect("no NaN in a Jacobi recurrence")
    });
    let values = order.iter().map(|&i| diagonal[i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (k, &i) in order.iter().enumerate() {
        for s in 0..n {
            vectors[k * n + s] = v[s * n + i];
        }
    }
    (values, vectors)
}

/// `ψ(x)`, the derivative of `ln Γ(x)`, by recurrence up to eight and then an asymptotic
/// series. Only the tests read it: the shapes' update steps on the rule's own slope instead of
/// solving the digamma form ([`step_beta_shapes`](super::maximisation::step_beta_shapes)).
#[cfg(test)]
pub(super) fn digamma(x: f64) -> f64 {
    let mut x = x;
    let mut result = 0.0;
    while x < 8.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    result + float::ln(x)
        - 0.5 * inverse
        - square * (1.0 / 12.0 - square * (1.0 / 120.0 - square / 252.0))
}

/// `ψ'(x)`, by the same shape.
pub(super) fn trigamma(x: f64) -> f64 {
    let mut x = x;
    let mut result = 0.0;
    while x < 8.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    result + inverse * (1.0 + 0.5 * inverse + square * (1.0 / 6.0 - square / 30.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_quadrature_integrates_a_beta_to_one_and_gets_its_mean_right() {
        for (a, b) in [(0.3, 2.0), (1.0, 1.0), (2.5, 5.0), (0.1, 0.9)] {
            let rule = BetaQuadrature::new(a, b, 24);
            let total: f64 = rule.weights.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-10,
                "Beta({a}, {b}) weights sum to {total}"
            );
            let mean: f64 = rule
                .nodes
                .iter()
                .zip(rule.weights.iter())
                .map(|(f, w)| f * w)
                .sum();
            let truth = a / (a + b);
            assert!(
                (mean - truth).abs() < 1e-8,
                "Beta({a}, {b}) mean {mean} against {truth}"
            );
        }
    }

    /// **The shape a neutral population has is the one a uniform grid cannot integrate**, so
    /// it is the one to check: `Beta(0.2, 1)` puts nearly all its mass below `f = 0.01`.
    #[test]
    fn the_quadrature_handles_the_rare_allele_pile_up() {
        let rule = BetaQuadrature::new(0.2, 1.0, 24);
        // E[2 f (1 − f)] under Beta(a, b) is 2ab / ((a+b)(a+b+1)).
        let heterozygosity: f64 = rule
            .nodes
            .iter()
            .zip(rule.weights.iter())
            .map(|(f, w)| w * 2.0 * f * (1.0 - f))
            .sum();
        let truth = 2.0 * 0.2 * 1.0 / (1.2 * 2.2);
        assert!(
            (heterozygosity - truth).abs() < 1e-8,
            "{heterozygosity} against {truth}"
        );
    }
}
