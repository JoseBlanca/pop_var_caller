//! **Small dense matrix routines the fits share**: a Cholesky factor, and cyclic Jacobi rotations
//! for the eigenvectors of a symmetric matrix. Pure linear algebra, with no model in it.
//!
//! Each fit keeps the parts of these that set its numbers: its own thresholds, its own ordering of
//! the eigenpairs, and its own order of summation when it turns a factor into an inverse. Only `+ −
//! × ÷` and `sqrt` are used, each rounded exactly by IEEE 754, in a fixed order, so every routine
//! gives the same bits on every platform.
//!
//! Matrices are `n × n`, row-major, in one slice.

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
            let image = [
                2.0 * vector[0] + vector[1],
                vector[0] + 2.0 * vector[1],
            ];
            assert!((image[0] - value * vector[0]).abs() < 1e-12);
            assert!((image[1] - value * vector[1]).abs() < 1e-12);
        }
    }
}
