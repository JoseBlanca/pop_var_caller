//! **Small matrix routines for the fits**, pure linear algebra with no model in it:
//!
//! - shared by the fits: a Cholesky factor, and cyclic Jacobi rotations for the eigenvectors of a
//!   symmetric matrix;
//! - used by the repeat-tract (STR) fit alone: its inverse and product of dense matrices, for the
//!   standard errors, and a banded symmetric matrix with its Cholesky factor, for Newton's method
//!   along each curve.
//!
//! Each fit keeps the parts of these that set its numbers: its own thresholds, its own ordering of
//! the eigenpairs, and its own order of summation when it turns a factor into an inverse. Apart
//! from the banded factor's log-determinant, which takes [`float::ln`], only `+ − × ÷` and `sqrt`
//! are used, each rounded exactly by IEEE 754, in a fixed order, so every routine gives the same
//! bits on every platform.
//!
//! A dense matrix is `n × n`, row-major, in one slice; a banded one keeps only its lower band
//! ([`SymmetricBand`]).

use crate::float;

/// **The lower-triangular Cholesky factor `L` of a symmetric matrix**, `matrix = L Lᵀ`, built row
/// by row; or the position of the first pivot not above its `floor` (or not finite), where the
/// factorisation stops.
///
/// A `floor` of zero everywhere asks only that the matrix be positive definite; the SNP/indel fit's
/// standard errors pass a share of each parameter's own curvature instead, to drop a parameter the
/// others mimic.
pub(super) fn cholesky_factor(matrix: &[f64], n: usize, floor: &[f64]) -> Result<Vec<f64>, usize> {
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
    Ok(lower)
}

/// The inverse of a symmetric positive-definite matrix by its Cholesky factor
/// ([`cholesky_factor`]); `None` when it is not positive definite.
///
/// **The factor is shared with the SNP/indel fit's standard errors; the inverse is not.** This one
/// solves for each column of the identity by forward then backward substitution, and the SNP/indel
/// side forms `L⁻ᵀ L⁻¹`. The two sum in different orders, so moving this one to the other's would
/// change the last digits of the STR standard errors.
pub(super) fn dense_inverse(matrix: &[f64], size: usize) -> Option<Vec<f64>> {
    let lower = cholesky_factor(matrix, size, &vec![0.0; size]).ok()?;
    let mut inverse = vec![0.0; size * size];
    let mut column = vec![0.0; size];
    for target in 0..size {
        column.fill(0.0);
        column[target] = 1.0;
        for i in 0..size {
            let known: f64 = (0..i).map(|k| lower[i * size + k] * column[k]).sum();
            column[i] = (column[i] - known) / lower[i * size + i];
        }
        for i in (0..size).rev() {
            let known: f64 = ((i + 1)..size)
                .map(|k| lower[k * size + i] * column[k])
                .sum();
            column[i] = (column[i] - known) / lower[i * size + i];
        }
        for i in 0..size {
            inverse[i * size + target] = column[i];
        }
    }
    Some(inverse)
}

pub(super) fn dense_product(left: &[f64], right: &[f64], size: usize) -> Vec<f64> {
    let mut out = vec![0.0; size * size];
    for i in 0..size {
        for k in 0..size {
            let a = left[i * size + k];
            if a == 0.0 {
                continue;
            }
            for j in 0..size {
                out[i * size + j] += a * right[k * size + j];
            }
        }
    }
    out
}

/// **Cyclic Jacobi rotations of a symmetric matrix until it is diagonal**: its diagonal (the
/// eigenvalues, in no particular order) and the accumulated rotation `v`, whose column `i`
/// (`v[s * n + i]` for each `s`) is the eigenvector of `diagonal[i]`.
///
/// A sweep visits every pair above the diagonal and skips a pair whose entry is below
/// `skip_below` in size; the rotations stop once the sum of the squared entries above the diagonal
/// falls below `stop_below`, or after 100 sweeps. The two fits pass different thresholds (the
/// quadrature rule 10⁻²² and 10⁻¹⁸, the contamination fit's ancestry axes 10⁻¹⁸ and 10⁻¹⁵), so
/// they are arguments rather than constants here.
pub(super) fn jacobi_rotations(
    matrix: &[f64],
    n: usize,
    stop_below: f64,
    skip_below: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut a = matrix.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum();
        if off < stop_below {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p * n + q].abs() < skip_below {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let diagonal = (0..n).map(|i| a[i * n + i]).collect();
    (diagonal, v)
}

/// A symmetric matrix kept as its lower band: entry `(i, i − d)` for `d ≤ width`.
#[derive(Clone)]
pub(super) struct SymmetricBand {
    pub(super) size: usize,
    width: usize,
    lower: Vec<f64>,
}

impl SymmetricBand {
    pub(super) fn zeros(size: usize, width: usize) -> Self {
        Self {
            size,
            width,
            lower: vec![0.0; size * (width + 1)],
        }
    }

    /// Adds `value` at `(row, column)` (and so at `(column, row)`).
    pub(super) fn add(&mut self, row: usize, column: usize, value: f64) {
        let (high, low) = if row >= column {
            (row, column)
        } else {
            (column, row)
        };
        let distance = high - low;
        debug_assert!(distance <= self.width, "outside the band");
        self.lower[high * (self.width + 1) + distance] += value;
    }

    /// The lower-triangular factor `L` with `L Lᵀ` this matrix, or `None` when it is not positive
    /// definite.
    pub(super) fn cholesky(mut self) -> Option<BandCholesky> {
        let w = self.width;
        let at = |lower: &[f64], i: usize, j: usize| lower[i * (w + 1) + (i - j)];
        for i in 0..self.size {
            let first = i.saturating_sub(w);
            for j in first..=i {
                let mut sum = at(&self.lower, i, j);
                for k in first.max(j.saturating_sub(w))..j {
                    sum -= at(&self.lower, i, k) * at(&self.lower, j, k);
                }
                let value = if i == j {
                    if sum <= 0.0 || !sum.is_finite() {
                        return None;
                    }
                    sum.sqrt()
                } else {
                    sum / at(&self.lower, j, j)
                };
                self.lower[i * (w + 1) + (i - j)] = value;
            }
        }
        Some(BandCholesky {
            size: self.size,
            width: w,
            lower: self.lower,
        })
    }
}

pub(super) struct BandCholesky {
    size: usize,
    width: usize,
    lower: Vec<f64>,
}

impl BandCholesky {
    fn at(&self, i: usize, j: usize) -> f64 {
        self.lower[i * (self.width + 1) + (i - j)]
    }

    pub(super) fn ln_determinant(&self) -> f64 {
        2.0 * (0..self.size)
            .map(|i| float::ln(self.at(i, i)))
            .sum::<f64>()
    }

    /// Overwrites `rhs` with the solution of `L Lᵀ x = rhs`.
    pub(super) fn solve(&self, rhs: &mut [f64]) {
        for i in 0..self.size {
            let known: f64 = (i.saturating_sub(self.width)..i)
                .map(|k| self.at(i, k) * rhs[k])
                .sum();
            rhs[i] = (rhs[i] - known) / self.at(i, i);
        }
        for i in (0..self.size).rev() {
            let known: f64 = ((i + 1)..(i + self.width + 1).min(self.size))
                .map(|k| self.at(k, i) * rhs[k])
                .sum();
            rhs[i] = (rhs[i] - known) / self.at(i, i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_factor_times_its_transpose_is_the_matrix() {
        let matrix = [4.0, 2.0, 0.4, 2.0, 5.0, 1.0, 0.4, 1.0, 3.0];
        let lower = cholesky_factor(&matrix, 3, &[0.0; 3]).expect("positive definite");
        for row in 0..3 {
            for column in 0..3 {
                let product: f64 = (0..3)
                    .map(|j| lower[row * 3 + j] * lower[column * 3 + j])
                    .sum();
                assert!((product - matrix[row * 3 + column]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn the_factor_stops_at_the_first_pivot_not_above_its_floor() {
        // The second row repeats the first, so nothing is left of its pivot.
        let matrix = [1.0, 1.0, 1.0, 1.0];
        assert_eq!(cholesky_factor(&matrix, 2, &[0.0; 2]), Err(1));
        assert_eq!(cholesky_factor(&[2.0], 1, &[3.0]), Err(0));
        assert_eq!(cholesky_factor(&[f64::NAN], 1, &[0.0]), Err(0));
    }

    #[test]
    fn the_rotations_find_the_eigenpairs_of_a_symmetric_matrix() {
        // [[2, 1], [1, 2]]: eigenvalues 1 and 3, along (1, −1) and (1, 1).
        let (diagonal, v) = jacobi_rotations(&[2.0, 1.0, 1.0, 2.0], 2, 1e-22, 1e-18);
        let mut values = diagonal.clone();
        values.sort_by(f64::total_cmp);
        assert!((values[0] - 1.0).abs() < 1e-12 && (values[1] - 3.0).abs() < 1e-12);
        for (i, value) in diagonal.iter().enumerate() {
            let vector = [v[i], v[2 + i]];
            let image = [2.0 * vector[0] + vector[1], vector[0] + 2.0 * vector[1]];
            assert!((image[0] - value * vector[0]).abs() < 1e-12);
            assert!((image[1] - value * vector[1]).abs() < 1e-12);
        }
    }

    /// **The banded Cholesky solves what a dense solve does**, on a random positive-definite
    /// matrix of bandwidth four.
    #[test]
    #[allow(
        clippy::needless_range_loop,
        reason = "fills dense[i][j] and dense[j][i] together"
    )]
    fn the_banded_solve_is_the_dense_one() {
        let (size, width) = (12, 4);
        let mut band = SymmetricBand::zeros(size, width);
        let mut dense = vec![vec![0.0; size]; size];
        let mut state = 7_u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        for i in 0..size {
            for j in i.saturating_sub(width)..=i {
                let value = if i == j {
                    10.0 + uniform()
                } else {
                    uniform() - 0.5
                };
                band.add(i, j, value);
                dense[i][j] = value;
                dense[j][i] = value;
            }
        }
        let rhs: Vec<f64> = (0..size).map(|i| i as f64 - 3.0).collect();
        let mut solved = rhs.clone();
        band.cholesky()
            .expect("positive definite")
            .solve(&mut solved);
        for (row, expected) in dense.iter().zip(&rhs) {
            let product: f64 = row.iter().zip(&solved).map(|(a, x)| a * x).sum();
            assert!(
                (product - expected).abs() < 1e-10,
                "{product} against {expected}"
            );
        }
    }
}
