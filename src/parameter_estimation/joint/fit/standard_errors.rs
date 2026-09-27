//! Each SNP/indel parameter's **standard error**, from the information a pass summed in blocks.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §3.2. Build order:
//! `doc/devel/implementation_plans/fit_precision.md`, step A3.
//!
//! # From information to errors
//!
//! A parameter's standard error is how far its estimate would typically move if the same kind of
//! data were drawn again. It is the square root of the parameter's entry on the diagonal of the
//! **inverse** of the information matrix — not one over its own diagonal entry, because a parameter
//! whose effect another parameter can mimic is less well determined than its own curvature says.
//!
//! The pass keeps the information in blocks ([`InformationSums`]): the cohort's eight parameters
//! with each other (`C`), each sample's three with each other (`A_s`), and each sample's three with
//! the cohort's (`B_s`). Two samples' parameters are never paired, so the matrix is an **arrow**:
//! the cohort's row and column run along one edge and every sample's block sits on the diagonal,
//! alone. An arrow is inverted exactly, block by block, without ever forming it:
//!
//! - the cohort's errors come from `C − Σ_s B_sᵀ A_s⁻¹ B_s` inverted — the cohort's information,
//!   less what each sample's own parameters could explain of it;
//! - a sample's errors come from `A_s⁻¹ + A_s⁻¹ B_s V B_sᵀ A_s⁻¹`, where `V` is the cohort's inverse
//!   just computed — the sample's own uncertainty, plus the cohort's carried through the parameters
//!   the two share.
//!
//! What the arrow leaves out is the direct pairing of two samples' scores, which the pass does not
//! sum; the tests measure what that costs against the full matrix on small cohorts.
//!
//! # When a parameter has no error
//!
//! A parameter is left out of the inversion, and says why ([`StandardError`]), when:
//!
//! - **no position's likelihood depends on it** — its diagonal entry is zero: a sample with no
//!   reads, or the duplicated class when the run does not fit it;
//! - **the fit does not move it** — a sample's homozygote excess in a cohort of one sample
//!   ([`fits_homozygote_excess`](super::fits_homozygote_excess)). Its information is not zero
//!   there, but the fit holds it at its start, so no error describes where it ended;
//! - **the data cannot tell it apart from the parameters before it** — other parameters can mimic
//!   its effect on every position, so everything its curvature says is already said by theirs.
//!   Measured on drawn cohorts with the duplicated class on: at one sample the four
//!   frequency-density parameters cannot all be told apart, since three genotypes give three
//!   frequencies to fit four numbers; at two samples the duplicated class's share and carrier
//!   shapes cannot. **Only the parameter that fails is dropped**, and the others are inverted
//!   without it: their errors are then those the data give with it held where it was fitted.
//!
//! - **the data do not place it anywhere in the interval the fit keeps it in** — its error came out
//!   wider than that whole interval. Measured at two samples with the duplicated class on: once the
//!   carrier Beta's second shape is dropped, the invariant share's error came out at 3.6 though it
//!   lives in [0, 1], and the density's shapes' at 720 and 1,148 though they live in [0.02, 50];
//!   the full matrix, with the products of the two samples' scores in, does not invert over those
//!   parameters at all.
//!
//! A parameter counts as not told apart when the curvature left to it once the parameters before
//! it are accounted for is below [`IDENTIFIED_SHARE`] of its own curvature — which parameter of a
//! mimicking set that is depends on their order, the cohort's first. Without the threshold, a
//! remainder that is only rounding would be inverted into an error of millions (measured: 2 × 10⁶
//! on a carrier shape, 13.5 on a share that lives in [0, 1], at two samples and 300,000 positions).

#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "only this module's tests compute the errors until the fit does, at plan step A5"
    )
)]

use super::information::{COHORT_PARAMETERS, InformationSums, SAMPLE_PARAMETERS, sample};
use super::{
    BETA_SHAPE_BOUNDS, CLEAN_ERROR_BOUNDS, DUPLICATED_SHARE_BOUNDS, HOM_EXCESS_BOUNDS,
    NOISY_ERROR_BOUNDS, NOISY_SHARE_BOUNDS, P_FIXED_ALT_BOUNDS, P_INVARIANT_BOUNDS,
    fits_homozygote_excess,
};

/// One parameter's standard error, or why it has none (module doc).
#[derive(Copy, Clone, Debug, PartialEq)]
pub(super) enum StandardError {
    /// The error, on the parameter's own scale.
    Estimated(f64),
    /// No position's likelihood depends on the parameter.
    NoInformation,
    /// The fit holds the parameter where it starts: a sample's homozygote excess at one sample.
    HeldFixed,
    /// Other parameters can mimic the parameter's effect, so the data cannot tell it apart from
    /// them.
    NotIdentified,
    /// The error came out wider than the whole interval the fit keeps the parameter in: the data
    /// do not place it anywhere in that interval. The width it came out at is kept for reporting;
    /// it is not an error.
    WiderThanItsRange(f64),
}

impl StandardError {
    /// The error, where there is one.
    pub(super) fn value(self) -> Option<f64> {
        match self {
            Self::Estimated(error) => Some(error),
            Self::NoInformation
            | Self::HeldFixed
            | Self::NotIdentified
            | Self::WiderThanItsRange(_) => None,
        }
    }
}

/// Every fitted parameter's standard error, or why it has none.
#[derive(Clone, Debug)]
pub(super) struct StandardErrors {
    /// The cohort's eight, indexed by [`cohort`](super::information::cohort).
    pub cohort: [StandardError; COHORT_PARAMETERS],
    /// Per sample, in the order the fit iterates samples, indexed by [`sample`].
    pub samples: Vec<[StandardError; SAMPLE_PARAMETERS]>,
}

impl StandardErrors {
    /// The errors of every parameter, from the information one pass summed at the parameters
    /// being reported.
    pub(super) fn of(sums: &InformationSums) -> Self {
        let samples = sums.sample_blocks.len();
        let excess_is_fitted = fits_homozygote_excess(samples);
        let mut errors = Self {
            cohort: [StandardError::NoInformation; COHORT_PARAMETERS],
            samples: vec![[StandardError::NoInformation; SAMPLE_PARAMETERS]; samples],
        };

        // Which cohort parameters carry information at all.
        let cohort_informed: Vec<usize> = (0..COHORT_PARAMETERS)
            .filter(|&i| has_information(sums.cohort[i * COHORT_PARAMETERS + i]))
            .collect();
        let c = cohort_informed.len();

        // Each sample's own block, inverted over the parameters it identifies, and `A_s⁻¹ B_s`;
        // the sum `Σ B_sᵀ A_s⁻¹ B_s` over the informed cohort parameters.
        let mut explained = vec![0.0; c * c];
        let mut per_sample: Vec<Option<SampleBlockInverse>> = Vec::with_capacity(samples);
        for s in 0..samples {
            let own_diagonal = |j: usize| sums.sample_blocks[s][j * SAMPLE_PARAMETERS + j];
            let mut informed = Vec::with_capacity(SAMPLE_PARAMETERS);
            for j in 0..SAMPLE_PARAMETERS {
                if !has_information(own_diagonal(j)) {
                    continue;
                }
                if j == sample::HOMOZYGOTE_EXCESS && !excess_is_fitted {
                    errors.samples[s][j] = StandardError::HeldFixed;
                    continue;
                }
                informed.push(j);
            }
            let own: Vec<f64> = square_of(&informed, |row, column| {
                sums.sample_blocks[s][row * SAMPLE_PARAMETERS + column]
            });
            let reference: Vec<f64> = informed.iter().map(|&j| own_diagonal(j)).collect();
            let identified = invert_identified(&own, informed.len(), &reference);
            for &position in &identified.dropped {
                errors.samples[s][informed[position]] = StandardError::NotIdentified;
            }
            let kept: Vec<usize> = identified.kept.iter().map(|&p| informed[p]).collect();
            if kept.is_empty() {
                per_sample.push(None);
                continue;
            }
            let k = kept.len();
            // `B_s` restricted: the sample's kept parameters (rows) with the informed cohort ones.
            let with_cohort: Vec<f64> = kept
                .iter()
                .flat_map(|&row| {
                    cohort_informed.iter().map(move |&column| {
                        sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column]
                    })
                })
                .collect();
            let own_inverse_times_cross = multiply(&identified.inverse, &with_cohort, k, k, c);
            // Σ_s B_sᵀ (A_s⁻¹ B_s), added sample by sample in the fit's order.
            for p in 0..c {
                for q in 0..c {
                    let mut entry = 0.0;
                    for j in 0..k {
                        entry += with_cohort[j * c + p] * own_inverse_times_cross[j * c + q];
                    }
                    explained[p * c + q] += entry;
                }
            }
            per_sample.push(Some(SampleBlockInverse {
                kept,
                own_inverse: identified.inverse,
                own_inverse_times_cross,
            }));
        }

        // The cohort's information less what the samples explain, inverted over the parameters it
        // identifies; each judged against its own curvature before the samples' share is taken.
        let reduced: Vec<f64> = (0..c)
            .flat_map(|p| {
                let explained = &explained;
                let cohort_informed = &cohort_informed;
                (0..c).map(move |q| {
                    sums.cohort[cohort_informed[p] * COHORT_PARAMETERS + cohort_informed[q]]
                        - explained[p * c + q]
                })
            })
            .collect();
        let reference: Vec<f64> = cohort_informed
            .iter()
            .map(|&i| sums.cohort[i * COHORT_PARAMETERS + i])
            .collect();
        let identified = invert_identified(&reduced, c, &reference);
        for &position in &identified.dropped {
            errors.cohort[cohort_informed[position]] = StandardError::NotIdentified;
        }
        let kept = &identified.kept;
        let v = kept.len();
        let covariance = &identified.inverse;
        for (p, &position) in kept.iter().enumerate() {
            let slot = cohort_informed[position];
            errors.cohort[slot] = error_from_variance(covariance[p * v + p], COHORT_BOUNDS[slot]);
        }

        // Each sample's own inverse, plus the kept cohort parameters' uncertainty carried through
        // `A_s⁻¹ B_s`. A dropped cohort parameter is held where it was fitted: it carries none.
        for (s, inverse) in per_sample.iter().enumerate() {
            let Some(inverse) = inverse else { continue };
            let k = inverse.kept.len();
            for (j, &which) in inverse.kept.iter().enumerate() {
                let mut carried = 0.0;
                for (p, &cp) in kept.iter().enumerate() {
                    let mut row_times_covariance = 0.0;
                    for (q, &cq) in kept.iter().enumerate() {
                        row_times_covariance +=
                            covariance[p * v + q] * inverse.own_inverse_times_cross[j * c + cq];
                    }
                    carried += inverse.own_inverse_times_cross[j * c + cp] * row_times_covariance;
                }
                errors.samples[s][which] = error_from_variance(
                    inverse.own_inverse[j * k + j] + carried,
                    SAMPLE_BOUNDS[which],
                );
            }
        }
        errors
    }
}

/// One sample's own block inverted: which of its parameters were kept, the inverse over them, and
/// that inverse times the sample's block with the informed cohort parameters (`A_s⁻¹ B_s`, `k × c`).
struct SampleBlockInverse {
    kept: Vec<usize>,
    own_inverse: Vec<f64>,
    own_inverse_times_cross: Vec<f64>,
}

/// Whether a diagonal entry of the information says the parameter was informed at all.
fn has_information(diagonal: f64) -> bool {
    diagonal > 0.0 && diagonal.is_finite()
}

/// The standard error a variance gives, for a parameter the fit keeps within `bounds`. A variance
/// left not positive by rounding means the parameter was not told apart after all; an error wider
/// than the whole interval says nothing about where in it the parameter lies.
fn error_from_variance(variance: f64, bounds: (f64, f64)) -> StandardError {
    if !(variance > 0.0 && variance.is_finite()) {
        return StandardError::NotIdentified;
    }
    let error = variance.sqrt();
    if error > bounds.1 - bounds.0 {
        StandardError::WiderThanItsRange(error)
    } else {
        StandardError::Estimated(error)
    }
}

/// The interval the fit keeps each cohort parameter in, by slot ([`cohort`](super::information::cohort)).
const COHORT_BOUNDS: [(f64, f64); COHORT_PARAMETERS] = [
    NOISY_SHARE_BOUNDS,
    P_INVARIANT_BOUNDS,
    P_FIXED_ALT_BOUNDS,
    BETA_SHAPE_BOUNDS,
    BETA_SHAPE_BOUNDS,
    DUPLICATED_SHARE_BOUNDS,
    BETA_SHAPE_BOUNDS,
    BETA_SHAPE_BOUNDS,
];

/// The interval the fit keeps each of a sample's parameters in, by slot ([`sample`]).
const SAMPLE_BOUNDS: [(f64, f64); SAMPLE_PARAMETERS] =
    [CLEAN_ERROR_BOUNDS, NOISY_ERROR_BOUNDS, HOM_EXCESS_BOUNDS];

/// The square matrix `entry(row, column)` over `indices`, row-major.
fn square_of(indices: &[usize], entry: impl Fn(usize, usize) -> f64) -> Vec<f64> {
    indices
        .iter()
        .flat_map(|&row| indices.iter().map(move |&column| (row, column)))
        .map(|(row, column)| entry(row, column))
        .collect()
}

/// `left` (`rows × inner`) times `right` (`inner × columns`), both row-major.
fn multiply(left: &[f64], right: &[f64], rows: usize, inner: usize, columns: usize) -> Vec<f64> {
    let mut product = vec![0.0; rows * columns];
    for row in 0..rows {
        for column in 0..columns {
            let mut entry = 0.0;
            for j in 0..inner {
                entry += left[row * inner + j] * right[j * columns + column];
            }
            product[row * columns + column] = entry;
        }
    }
    product
}

/// The share of a parameter's own curvature that must be left once the parameters before it are
/// accounted for, for the data to tell it apart from them.
///
/// The share left is `1 − R²` of the parameter on those before it, so rescaling a parameter does
/// not change it. A remainder that is only rounding sits near 10⁻¹² to 10⁻¹⁶ (measured: 3.2 × 10⁻¹²
/// for a carrier shape two samples cannot determine, at 300,000 positions); a parameter the data
/// determine even poorly keeps far more than 10⁻⁸ — its error would otherwise be ten thousand times
/// what its own curvature alone gives.
pub(super) const IDENTIFIED_SHARE: f64 = 1e-8;

/// A symmetric matrix inverted over the parameters it identifies: which were kept and which dropped
/// (positions in the matrix given), and the inverse over the kept ones, in their order.
pub(super) struct Identified {
    pub kept: Vec<usize>,
    pub dropped: Vec<usize>,
    pub inverse: Vec<f64>,
}

/// Invert a symmetric `n × n` matrix (row-major) over the parameters it identifies: a Cholesky
/// factor is built in order, and a parameter whose pivot keeps less than [`IDENTIFIED_SHARE`] of
/// `reference[i]` — its curvature before anything is taken from it — is dropped and the factor
/// rebuilt without it.
///
/// Only `+ − × ÷` and `sqrt`, each rounded exactly, in a fixed order: the same bits on every
/// platform.
pub(super) fn invert_identified(matrix: &[f64], n: usize, reference: &[f64]) -> Identified {
    let mut kept: Vec<usize> = (0..n).collect();
    let mut dropped = Vec::new();
    loop {
        let m = kept.len();
        let sub = square_of(&kept, |row, column| matrix[row * n + column]);
        let floor: Vec<f64> = kept
            .iter()
            .map(|&i| IDENTIFIED_SHARE * reference[i])
            .collect();
        match inverse_by_cholesky(&sub, m, &floor) {
            Ok(inverse) => {
                return Identified {
                    kept,
                    dropped,
                    inverse,
                };
            }
            Err(position) => dropped.push(kept.remove(position)),
        }
    }
}

/// The inverse of a symmetric positive-definite `n × n` matrix (row-major), by its Cholesky factor;
/// `None` when any pivot is not above zero — the matrix has no inverse a variance could be read
/// from. The full-matrix comparison in the tests uses it.
pub(super) fn inverse_of_positive_definite(matrix: &[f64], n: usize) -> Option<Vec<f64>> {
    inverse_by_cholesky(matrix, n, &vec![0.0; n]).ok()
}

/// The inverse by the Cholesky factor `matrix = L Lᵀ`, or the position of the first pivot not above
/// its `floor` (or not finite).
fn inverse_by_cholesky(matrix: &[f64], n: usize, floor: &[f64]) -> Result<Vec<f64>, usize> {
    let mut lower = vec![0.0; n * n];
    for row in 0..n {
        for column in 0..=row {
            let mut entry = matrix[row * n + column];
            for j in 0..column {
                entry -= lower[row * n + j] * lower[column * n + j];
            }
            if row == column {
                if !(entry > floor[row] && entry.is_finite()) {
                    return Err(row);
                }
                lower[row * n + row] = entry.sqrt();
            } else {
                lower[row * n + column] = entry / lower[column * n + column];
            }
        }
    }
    // `L⁻¹`, lower triangular, by forward substitution on each column of the identity.
    let mut lower_inverse = vec![0.0; n * n];
    for column in 0..n {
        for row in column..n {
            let mut entry = if row == column { 1.0 } else { 0.0 };
            for j in column..row {
                entry -= lower[row * n + j] * lower_inverse[j * n + column];
            }
            lower_inverse[row * n + column] = entry / lower[row * n + row];
        }
    }
    // `matrix⁻¹ = L⁻ᵀ L⁻¹`.
    let mut inverse = vec![0.0; n * n];
    for row in 0..n {
        for column in 0..n {
            let mut entry = 0.0;
            for j in row.max(column)..n {
                entry += lower_inverse[j * n + row] * lower_inverse[j * n + column];
            }
            inverse[row * n + column] = entry;
        }
    }
    Ok(inverse)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::fit::information::cohort;

    /// A random symmetric positive-definite matrix, `G Gᵀ + n·I` with `G` drawn from a fixed
    /// sequence — well conditioned, so its inverse is known to many digits.
    fn positive_definite(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        let mut draw = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let g: Vec<f64> = (0..n * n).map(|_| draw()).collect();
        let mut matrix = multiply(
            &g,
            &(0..n * n)
                .map(|i| g[(i % n) * n + i / n])
                .collect::<Vec<_>>(),
            n,
            n,
            n,
        );
        for i in 0..n {
            matrix[i * n + i] += n as f64;
        }
        matrix
    }

    /// **The Cholesky inverse is an inverse**: the matrix times it is the identity to 10⁻¹².
    #[test]
    fn the_inverse_times_the_matrix_is_the_identity() {
        for n in [1, 3, 8, 20] {
            let matrix = positive_definite(n, 7 + n as u64);
            let inverse = inverse_of_positive_definite(&matrix, n).expect("positive definite");
            let product = multiply(&matrix, &inverse, n, n, n);
            for row in 0..n {
                for column in 0..n {
                    let expected = if row == column { 1.0 } else { 0.0 };
                    assert!(
                        (product[row * n + column] - expected).abs() < 1e-12,
                        "n {n}: ({row}, {column}) is {}",
                        product[row * n + column]
                    );
                }
            }
        }
    }

    /// **A matrix that is not positive definite has no inverse here**: a zero, a negative and a
    /// NaN pivot each give `None`.
    #[test]
    fn a_matrix_without_a_positive_pivot_has_no_inverse() {
        assert!(inverse_of_positive_definite(&[1.0, 1.0, 1.0, 1.0], 2).is_none());
        assert!(inverse_of_positive_definite(&[1.0, 2.0, 2.0, 1.0], 2).is_none());
        assert!(inverse_of_positive_definite(&[f64::NAN], 1).is_none());
    }

    /// Blocks of a random arrow and the whole arrow beside them: a random symmetric matrix of the
    /// full size with every entry between two different samples set to zero, and each diagonal entry
    /// then raised above the sum of its row's other entries' sizes — which makes it positive
    /// definite.
    fn an_arrow(samples: usize, seed: u64) -> (InformationSums, Vec<f64>) {
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let mut full = positive_definite(n, seed);
        let owner = |i: usize| {
            (i >= COHORT_PARAMETERS).then(|| (i - COHORT_PARAMETERS) / SAMPLE_PARAMETERS)
        };
        for row in 0..n {
            for column in 0..n {
                if let (Some(a), Some(b)) = (owner(row), owner(column))
                    && a != b
                {
                    full[row * n + column] = 0.0;
                }
            }
        }
        for row in 0..n {
            let off_diagonal: f64 = (0..n)
                .filter(|&column| column != row)
                .map(|column| full[row * n + column].abs())
                .sum();
            full[row * n + row] = off_diagonal + 1.0;
        }
        // Information of the size a few thousand positions give, so every error is about 10⁻³ —
        // inside every parameter's range.
        for entry in &mut full {
            *entry *= 1e6;
        }
        let mut sums = InformationSums::new(samples);
        for row in 0..COHORT_PARAMETERS {
            for column in 0..COHORT_PARAMETERS {
                sums.cohort[row * COHORT_PARAMETERS + column] = full[row * n + column];
            }
        }
        for s in 0..samples {
            let first = COHORT_PARAMETERS + SAMPLE_PARAMETERS * s;
            for row in 0..SAMPLE_PARAMETERS {
                for column in 0..SAMPLE_PARAMETERS {
                    sums.sample_blocks[s][row * SAMPLE_PARAMETERS + column] =
                        full[(first + row) * n + first + column];
                }
                for column in 0..COHORT_PARAMETERS {
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column] =
                        full[(first + row) * n + column];
                }
            }
        }
        (sums, full)
    }

    /// **The blocks give exactly the errors of the whole arrow inverted at once**: on an arrow of
    /// five samples (23 parameters), every error agrees with the square root of the dense inverse's
    /// diagonal to 10⁻¹² relative. A cohort term with its sign flipped, or a sample's own inverse
    /// without the cohort's uncertainty carried into it, fails.
    /// The error of layout parameter `i` — the cohort's eight, then each sample's three.
    fn error_at(errors: &StandardErrors, i: usize) -> StandardError {
        if i < COHORT_PARAMETERS {
            errors.cohort[i]
        } else {
            let at = i - COHORT_PARAMETERS;
            errors.samples[at / SAMPLE_PARAMETERS][at % SAMPLE_PARAMETERS]
        }
    }

    /// Every parameter's error from the whole arrow inverted densely, with the parameters `removed`
    /// says left out; `None` for those.
    fn dense_errors_without(
        full: &[f64],
        n: usize,
        removed: impl Fn(usize) -> bool,
    ) -> Vec<Option<f64>> {
        let kept: Vec<usize> = (0..n).filter(|&i| !removed(i)).collect();
        let inverse = inverse_of_positive_definite(
            &square_of(&kept, |row, column| full[row * n + column]),
            kept.len(),
        )
        .expect("positive definite");
        let mut errors = vec![None; n];
        for (at, &i) in kept.iter().enumerate() {
            errors[i] = Some(inverse[at * kept.len() + at].sqrt());
        }
        errors
    }

    /// Every error the blocks give, against the dense inverse's, to 10⁻¹² relative; the dropped
    /// parameters must be the ones the dense inverse left out.
    fn assert_errors_match(errors: &StandardErrors, dense: &[Option<f64>]) {
        for (i, expected) in dense.iter().enumerate() {
            match (error_at(errors, i).value(), expected) {
                (Some(got), Some(expected)) => assert!(
                    (got - expected).abs() < 1e-12 * expected,
                    "parameter {i}: {got} against {expected}"
                ),
                (None, None) => {}
                (got, expected) => panic!("parameter {i}: {got:?} against {expected:?}"),
            }
        }
    }

    /// **The blocks give exactly the errors of the whole arrow inverted at once**: on an arrow of
    /// five samples (23 parameters), every error agrees with the square root of the dense inverse's
    /// diagonal to 10⁻¹² relative. A cohort term with its sign flipped, or a sample's own inverse
    /// without the cohort's uncertainty carried into it, fails.
    #[test]
    fn the_blocks_give_the_errors_of_the_whole_arrow() {
        let samples = 5;
        let (sums, full) = an_arrow(samples, 11);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let errors = StandardErrors::of(&sums);
        assert_errors_match(&errors, &dense_errors_without(&full, n, |_| false));
    }

    /// **A parameter with no information has no error, and does not disturb the others'**: the
    /// duplicated class's three rows zeroed, and sample 1's whole block zeroed. Those say
    /// `NoInformation`; every other error equals the one the arrow gives with them removed.
    #[test]
    fn a_parameter_with_no_information_has_no_error() {
        let samples = 3;
        let (mut sums, full) = an_arrow(samples, 23);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let removed = |i: usize| {
            (cohort::DUPLICATED_SHARE..COHORT_PARAMETERS).contains(&i)
                || (COHORT_PARAMETERS + SAMPLE_PARAMETERS
                    ..COHORT_PARAMETERS + 2 * SAMPLE_PARAMETERS)
                    .contains(&i)
        };
        for i in cohort::DUPLICATED_SHARE..COHORT_PARAMETERS {
            for j in 0..COHORT_PARAMETERS {
                sums.cohort[i * COHORT_PARAMETERS + j] = 0.0;
                sums.cohort[j * COHORT_PARAMETERS + i] = 0.0;
            }
            for s in 0..samples {
                for row in 0..SAMPLE_PARAMETERS {
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + i] = 0.0;
                }
            }
        }
        sums.sample_blocks[1] = [0.0; SAMPLE_PARAMETERS * SAMPLE_PARAMETERS];
        sums.sample_cohort_blocks[1] = [0.0; SAMPLE_PARAMETERS * COHORT_PARAMETERS];
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.cohort[cohort::DUPLICATED_SHARE..],
            [StandardError::NoInformation; 3]
        );
        assert_eq!(errors.samples[1], [StandardError::NoInformation; 3]);
        assert_errors_match(&errors, &dense_errors_without(&full, n, removed));
    }

    /// **At one sample the homozygote excess has no error — it is held fixed — though its
    /// information is not zero.** The sample's two error rates and the cohort's parameters still get
    /// theirs: those of the arrow with the excess's row and column removed.
    #[test]
    fn at_one_sample_the_homozygote_excess_is_held_fixed() {
        let (sums, full) = an_arrow(1, 31);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS;
        let excess = COHORT_PARAMETERS + sample::HOMOZYGOTE_EXCESS;
        assert!(full[excess * n + excess] > 0.0, "the excess is informed");
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.samples[0][sample::HOMOZYGOTE_EXCESS],
            StandardError::HeldFixed
        );
        assert_errors_match(&errors, &dense_errors_without(&full, n, |i| i == excess));
    }

    /// **A parameter the others can mimic exactly is dropped, alone, and the rest keep their
    /// errors**: cohort parameter 1 is made a copy of parameter 0 in every block, so the data cannot
    /// tell them apart. Parameter 1 — the later of the two — says `NotIdentified`, and every other
    /// error equals the arrow's with it removed. Refusing every error instead, as a failed inversion
    /// of the whole cohort block would, or inverting the rounding the copy leaves, fails.
    #[test]
    fn a_parameter_the_others_mimic_is_dropped_alone() {
        let samples = 2;
        let (mut sums, mut full) = an_arrow(samples, 41);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        for j in 0..n {
            full[n + j] = full[j];
        }
        for j in 0..n {
            full[j * n + 1] = full[j * n];
        }
        for j in 0..COHORT_PARAMETERS {
            sums.cohort[COHORT_PARAMETERS + j] = full[n + j];
            sums.cohort[j * COHORT_PARAMETERS + 1] = full[j * n + 1];
        }
        for s in 0..samples {
            for row in 0..SAMPLE_PARAMETERS {
                sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + 1] =
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS];
            }
        }
        let errors = StandardErrors::of(&sums);
        assert_eq!(errors.cohort[1], StandardError::NotIdentified);
        assert_errors_match(&errors, &dense_errors_without(&full, n, |i| i == 1));
    }

    /// **An error wider than the interval the fit keeps its parameter in is not reported as an
    /// error**: the duplicated share's information cut ten-thousandfold puts its error a hundred times
    /// higher, past the width of [10⁻⁹, 0.05]; it says so and keeps the width, and the others keep theirs.
    #[test]
    fn an_error_wider_than_its_range_is_not_an_error() {
        let samples = 2;
        let (mut sums, _) = an_arrow(samples, 61);
        let at = cohort::DUPLICATED_SHARE;
        let shrink = 1e-2;
        for j in 0..COHORT_PARAMETERS {
            sums.cohort[at * COHORT_PARAMETERS + j] *= shrink;
            sums.cohort[j * COHORT_PARAMETERS + at] *= shrink;
        }
        for s in 0..samples {
            for row in 0..SAMPLE_PARAMETERS {
                sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + at] *= shrink;
            }
        }
        let errors = StandardErrors::of(&sums);
        match errors.cohort[at] {
            StandardError::WiderThanItsRange(width) => {
                assert!(width > DUPLICATED_SHARE_BOUNDS.1, "{width}");
            }
            other => panic!("the duplicated share's error is {other:?}"),
        }
        for i in (0..COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples).filter(|&i| i != at) {
            assert!(error_at(&errors, i).value().is_some(), "parameter {i}");
        }
    }

    /// **A sample's parameter it cannot tell apart from its others is dropped alone, and the sample
    /// keeps the rest**: sample 0's mismapped rate made a copy of its clean rate in every block. That
    /// rate says `NotIdentified`; the sample's other two parameters, the other sample's and the
    /// cohort's equal the arrow's with it removed. Dropping the whole sample instead — which leaves
    /// the cohort's errors too small by what that sample explains — fails.
    #[test]
    fn a_sample_parameter_it_cannot_tell_apart_is_dropped_alone() {
        let samples = 2;
        let (mut sums, mut full) = an_arrow(samples, 67);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let (clean, noisy) = (COHORT_PARAMETERS, COHORT_PARAMETERS + 1);
        for j in 0..n {
            full[noisy * n + j] = full[clean * n + j];
        }
        for j in 0..n {
            full[j * n + noisy] = full[j * n + clean];
        }
        for row in 0..SAMPLE_PARAMETERS {
            for column in 0..SAMPLE_PARAMETERS {
                sums.sample_blocks[0][row * SAMPLE_PARAMETERS + column] =
                    full[(clean + row) * n + clean + column];
            }
            for column in 0..COHORT_PARAMETERS {
                sums.sample_cohort_blocks[0][row * COHORT_PARAMETERS + column] =
                    full[(clean + row) * n + column];
            }
        }
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.samples[0][sample::NOISY_ERROR_RATE],
            StandardError::NotIdentified
        );
        assert_errors_match(&errors, &dense_errors_without(&full, n, |i| i == noisy));
    }

    /// **The threshold keeps a parameter that keeps a part in a million of its curvature, and drops
    /// one that keeps a part in 10¹⁰** — two parameters correlated so that `1 − R²` is each of those.
    /// A threshold much stricter than [`IDENTIFIED_SHARE`] would drop the first, a looser one keep
    /// the second.
    #[test]
    fn the_threshold_sits_between_poor_information_and_rounding() {
        for (left_share, identified) in [(1e-6, true), (1e-10, false)] {
            let correlation = (1.0_f64 - left_share).sqrt();
            let matrix = [1.0, correlation, correlation, 1.0];
            let result = invert_identified(&matrix, 2, &[1.0, 1.0]);
            assert_eq!(
                result.kept.len() == 2,
                identified,
                "1 − R² of {left_share}: kept {:?}",
                result.kept
            );
        }
    }

    /// **With no cohort information every sample still gets its own block's errors**: the cohort's
    /// blocks zeroed, so each sample's error is the square root of its own block's inverse's diagonal,
    /// and the cohort's say `NoInformation`.
    #[test]
    fn without_cohort_information_the_samples_keep_their_own_errors() {
        let samples = 2;
        let (mut sums, _) = an_arrow(samples, 71);
        sums.cohort = [0.0; COHORT_PARAMETERS * COHORT_PARAMETERS];
        sums.sample_cohort_blocks
            .fill([0.0; SAMPLE_PARAMETERS * COHORT_PARAMETERS]);
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.cohort,
            [StandardError::NoInformation; COHORT_PARAMETERS]
        );
        for (s, own) in sums.sample_blocks.iter().enumerate() {
            let inverse = inverse_of_positive_definite(own, SAMPLE_PARAMETERS).expect("inverts");
            for j in 0..SAMPLE_PARAMETERS {
                let expected = inverse[j * SAMPLE_PARAMETERS + j].sqrt();
                let got = errors.samples[s][j].value().expect("an error");
                assert!((got - expected).abs() < 1e-12 * expected, "sample {s}, {j}");
            }
        }
    }

    /// **A parameter with no information in the middle of a block leaves every other error in its
    /// own slot**: the invariant share (the cohort's slot 1) and sample 0's mismapped rate (its slot
    /// 1, as for a read group with no mismapped positions) zeroed. Both say `NoInformation`, and every
    /// other error equals the arrow's without them — an error written back to its position among the
    /// kept parameters rather than to its own slot lands one slot early and fails.
    #[test]
    fn a_parameter_without_information_mid_block_leaves_the_others_in_place() {
        let samples = 2;
        let (mut sums, full) = an_arrow(samples, 73);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let (share, rate) = (
            cohort::P_INVARIANT,
            COHORT_PARAMETERS + sample::NOISY_ERROR_RATE,
        );
        for j in 0..COHORT_PARAMETERS {
            sums.cohort[share * COHORT_PARAMETERS + j] = 0.0;
            sums.cohort[j * COHORT_PARAMETERS + share] = 0.0;
        }
        for s in 0..samples {
            for row in 0..SAMPLE_PARAMETERS {
                sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + share] = 0.0;
            }
        }
        let which = sample::NOISY_ERROR_RATE;
        for j in 0..SAMPLE_PARAMETERS {
            sums.sample_blocks[0][which * SAMPLE_PARAMETERS + j] = 0.0;
            sums.sample_blocks[0][j * SAMPLE_PARAMETERS + which] = 0.0;
        }
        for column in 0..COHORT_PARAMETERS {
            sums.sample_cohort_blocks[0][which * COHORT_PARAMETERS + column] = 0.0;
        }
        let errors = StandardErrors::of(&sums);
        assert_eq!(errors.cohort[share], StandardError::NoInformation);
        assert_eq!(errors.samples[0][which], StandardError::NoInformation);
        assert_errors_match(
            &errors,
            &dense_errors_without(&full, n, |i| i == share || i == rate),
        );
    }

    /// **A diagonal entry that is not finite counts as no information**: an infinite and a NaN
    /// diagonal say `NoInformation`, and every other error equals the arrow's without them.
    #[test]
    fn a_diagonal_that_is_not_finite_is_no_information() {
        let samples = 2;
        let (mut sums, full) = an_arrow(samples, 79);
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let at = cohort::DENSITY_B;
        sums.cohort[at * COHORT_PARAMETERS + at] = f64::INFINITY;
        let excess = sample::HOMOZYGOTE_EXCESS;
        sums.sample_blocks[1][excess * SAMPLE_PARAMETERS + excess] = f64::NAN;
        let errors = StandardErrors::of(&sums);
        assert_eq!(errors.cohort[at], StandardError::NoInformation);
        assert_eq!(errors.samples[1][excess], StandardError::NoInformation);
        let excess_of_one = COHORT_PARAMETERS + SAMPLE_PARAMETERS + excess;
        assert_errors_match(
            &errors,
            &dense_errors_without(&full, n, |i| i == at || i == excess_of_one),
        );
    }

    /// **The threshold does not depend on a parameter's units**: one cohort parameter's scores
    /// scaled by 10⁻³ — its row and column by that, its own entry by its square — leave it
    /// identified, with an error 10³ times larger, and every other error unchanged.
    #[test]
    fn a_parameter_in_small_units_is_still_identified() {
        let samples = 2;
        let (sums, _) = an_arrow(samples, 53);
        let scale = 1e-3;
        let mut scaled = sums.clone();
        let at = cohort::DENSITY_A;
        for j in 0..COHORT_PARAMETERS {
            scaled.cohort[at * COHORT_PARAMETERS + j] *= scale;
            scaled.cohort[j * COHORT_PARAMETERS + at] *= scale;
        }
        for s in 0..samples {
            for row in 0..SAMPLE_PARAMETERS {
                scaled.sample_cohort_blocks[s][row * COHORT_PARAMETERS + at] *= scale;
            }
        }
        let (plain, small) = (StandardErrors::of(&sums), StandardErrors::of(&scaled));
        for i in 0..COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples {
            let expected = error_at(&plain, i).value().expect("identified")
                / if i == at { scale } else { 1.0 };
            let got = error_at(&small, i).value().expect("still identified");
            assert!(
                (got - expected).abs() < 1e-9 * expected,
                "parameter {i}: {got} against {expected}"
            );
        }
    }
}
