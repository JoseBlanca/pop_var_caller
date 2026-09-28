//! Each SNP/indel parameter's **standard error**, from the information a pass summed in blocks.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §3.2–3.3. Build order:
//! `doc/devel/implementation_plans/fit_precision.md`, steps A3 and A5.
//!
//! Every fit computes them once, at the parameters it returns, from the information its final pass
//! sums ([`fit_jointly`](super::fit_jointly)); it prints them as one line of the run's log
//! ([`StandardErrors::described`]) and, when the per-pass trace is on, one row a parameter
//! ([`StandardErrors::named`]).
//!
//! # From information to errors
//!
//! A parameter's standard error is how far its estimate would typically move if the same kind of
//! data were drawn again. It is the square root of the parameter's entry on the diagonal of the
//! **inverse** of the information matrix — not one over its own diagonal entry, because a parameter
//! whose effect another parameter can mimic is less well determined than its own curvature says.
//!
//! The pass keeps the information in blocks ([`InformationSums`]): the cohort's eight parameters
//! with each other (`C`), each sample's own with each other (`A_s`) — three for a sample of one
//! library, `1 + 2k` for one of `k` — and each sample's own with the cohort's (`B_s`). Two samples'
//! parameters are never paired, so the matrix is an **arrow**: the cohort's row and column run
//! along one edge and every sample's block sits on the diagonal, alone. An arrow is inverted
//! exactly, block by block, without ever forming it:
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

use super::information::{COHORT_PARAMETERS, InformationSums, cohort, own_parameters, sample};
use super::{
    BETA_SHAPE_BOUNDS, CLEAN_ERROR_BOUNDS, DUPLICATED_SHARE_BOUNDS, HOM_EXCESS_BOUNDS,
    NOISY_ERROR_BOUNDS, NOISY_SHARE_BOUNDS, P_FIXED_ALT_BOUNDS, P_INVARIANT_BOUNDS, Parameters,
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

    /// Why there is no error, in the words the run's log uses; `None` when there is one.
    fn reason(self) -> Option<&'static str> {
        match self {
            Self::Estimated(_) => None,
            Self::NoInformation => Some(NO_INFORMATION),
            Self::HeldFixed => Some(HELD_FIXED),
            Self::NotIdentified => Some(NOT_IDENTIFIED),
            Self::WiderThanItsRange(_) => Some(WIDER_THAN_ITS_RANGE),
        }
    }

    /// A parameter as the run's log prints it: its `value`, and its standard error or why it has
    /// none.
    fn printed_with(self, value: f64) -> String {
        let why = match self {
            Self::Estimated(error) => return format!("{value:.4e} ± {error:.2e}"),
            Self::WiderThanItsRange(width) => {
                return format!(
                    "{value:.4e}, no standard error (it came out at {width:.2e}, \
                     {WIDER_THAN_ITS_RANGE})"
                );
            }
            Self::NoInformation => NO_INFORMATION,
            Self::HeldFixed => HELD_FIXED,
            Self::NotIdentified => NOT_IDENTIFIED,
        };
        format!("{value:.4e}, no standard error ({why})")
    }
}

/// Why a parameter has no standard error, in the words the run's log uses.
const NO_INFORMATION: &str = "no information";
const HELD_FIXED: &str = "held fixed";
const NOT_IDENTIFIED: &str = "not identified";
const WIDER_THAN_ITS_RANGE: &str = "wider than the parameter's whole range";

/// The reasons in the order the run's log counts them.
const ABSENT_REASONS: [&str; 4] = [
    NO_INFORMATION,
    HELD_FIXED,
    NOT_IDENTIFIED,
    WIDER_THAN_ITS_RANGE,
];

/// What the run's log calls each cohort-level parameter, by slot
/// ([`cohort`](super::information::cohort)).
const COHORT_PARAMETER_NAMES: [&str; COHORT_PARAMETERS] = [
    "mismapped share",
    "invariant share",
    "fixed non-reference share",
    "allele-frequency shape a",
    "allele-frequency shape b",
    "duplicated share",
    "carrier-frequency shape a",
    "carrier-frequency shape b",
];

/// What the run's log calls each library's two rates, by noise class.
const RATE_NAMES: [&str; 2] = [
    "error rates at ordinary positions",
    "error rates at mismapped positions",
];

/// Every fitted parameter's standard error, or why it has none.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct StandardErrors {
    /// The cohort's eight, indexed by [`cohort`](super::information::cohort).
    pub cohort: [StandardError; COHORT_PARAMETERS],
    /// Per sample, in the order the fit iterates samples, in [`sample`]'s layout: its first
    /// library's two rates, its homozygote excess, then each further library's two rates.
    pub samples: Vec<Vec<StandardError>>,
}

impl StandardErrors {
    /// The errors of every parameter, from the information one pass summed at the parameters
    /// being reported.
    pub(super) fn of(sums: &InformationSums) -> Self {
        let samples = sums.sample_blocks.len();
        let excess_is_fitted = fits_homozygote_excess(samples);
        let mut errors = Self {
            cohort: [StandardError::NoInformation; COHORT_PARAMETERS],
            samples: (0..samples)
                .map(|s| vec![StandardError::NoInformation; sums.own_parameters_of(s)])
                .collect(),
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
            let n = sums.own_parameters_of(s);
            let own_diagonal = |j: usize| sums.sample_blocks[s][j * n + j];
            let mut informed = Vec::with_capacity(n);
            for j in 0..n {
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
                sums.sample_blocks[s][row * n + column]
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
                    bounds_of_own(which),
                );
            }
        }
        errors
    }

    /// **Every parameter's error, under the name the per-pass trace gives its value**
    /// ([`Parameters::named`], prefixed `standard_error:`, in its order), and NaN where it has none
    /// — so the trace can be read in units of each parameter's error. `group_index` says which read
    /// groups each sample's sections are, as the pass reads them: each read group carries its own
    /// two rates' errors, and a read group no sample holds none.
    pub(super) fn named(
        &self,
        parameters: &Parameters,
        group_index: &[Vec<usize>],
    ) -> Vec<(String, f64)> {
        debug_assert!(
            self.laid_out_for(group_index),
            "rows laid out for another cohort"
        );
        let value = |error: StandardError| error.value().unwrap_or(f64::NAN);
        let mut errors: Vec<f64> = self.cohort[..cohort_parameters_fitted(parameters)]
            .iter()
            .map(|&error| value(error))
            .collect();
        // Which sample holds each read group, and as which of its libraries.
        let mut held_by = vec![None; parameters.clean.len()];
        for (s, own) in group_index.iter().enumerate() {
            for (section, &group) in own.iter().enumerate() {
                held_by[group] = Some((s, section));
            }
        }
        for holder in held_by {
            let [clean, noisy] = match holder {
                Some((s, section)) => [
                    self.samples[s][sample::rate(section, 0)],
                    self.samples[s][sample::rate(section, 1)],
                ],
                None => [StandardError::NoInformation; 2],
            };
            errors.extend([value(clean), value(noisy)]);
        }
        errors.extend(
            self.samples
                .iter()
                .map(|row| value(row[sample::HOMOZYGOTE_EXCESS])),
        );
        let names = parameters.named();
        assert_eq!(
            names.len(),
            errors.len(),
            "one error for every parameter the trace names"
        );
        names
            .into_iter()
            .zip(errors)
            .map(|((name, _), error)| (format!("standard_error:{name}"), error))
            .collect()
    }

    /// **The errors as one line of the run's log**: each cohort-level parameter's value and
    /// standard error, or why it has none; for each library's two rates and each sample's
    /// homozygote excess, the median and the largest standard error and how many have none, by
    /// reason. `group_index` says which read groups each sample's sections are. The duplicated
    /// class's three are left out when the run does not fit it. The median of an even count is the
    /// lower of the middle two.
    pub(super) fn described(&self, parameters: &Parameters, group_index: &[Vec<usize>]) -> String {
        debug_assert!(
            self.laid_out_for(group_index),
            "rows laid out for another cohort"
        );
        let values = parameters.named();
        let mut parts = vec![
            (0..cohort_parameters_fitted(parameters))
                .map(|i| {
                    format!(
                        "{} {}",
                        COHORT_PARAMETER_NAMES[i],
                        self.cohort[i].printed_with(values[i].1)
                    )
                })
                .collect::<Vec<_>>()
                .join(", "),
        ];
        for (class, name) in RATE_NAMES.iter().enumerate() {
            let rates: Vec<StandardError> = self
                .samples
                .iter()
                .zip(group_index)
                .flat_map(|(row, own)| (0..own.len()).map(move |i| row[sample::rate(i, class)]))
                .collect();
            parts.push(format!(
                "{name} ({} read groups): {}",
                rates.len(),
                summary_of(&rates)
            ));
        }
        parts.push(format!(
            "homozygote excesses ({} samples): {}",
            self.samples.len(),
            summary_of(
                &self
                    .samples
                    .iter()
                    .map(|row| row[sample::HOMOZYGOTE_EXCESS])
                    .collect::<Vec<_>>()
            )
        ));
        parts.join("; ")
    }
}

impl StandardErrors {
    /// Whether each sample's row has the length its libraries in `group_index` give it.
    fn laid_out_for(&self, group_index: &[Vec<usize>]) -> bool {
        self.samples.len() == group_index.len()
            && self
                .samples
                .iter()
                .zip(group_index)
                .all(|(row, own)| row.len() == own_parameters(own.len()))
    }
}

/// How many of the cohort's slots are parameters of this fit: all eight with the duplicated class,
/// the first five without it — the order [`Parameters::named`] lists them in.
fn cohort_parameters_fitted(parameters: &Parameters) -> usize {
    if parameters.duplicated.is_some() {
        COHORT_PARAMETERS
    } else {
        cohort::DUPLICATED_SHARE
    }
}

/// The median and the largest of `errors`' standard errors, and how many have none, by reason.
fn summary_of(errors: &[StandardError]) -> String {
    let mut values: Vec<f64> = errors.iter().filter_map(|error| error.value()).collect();
    values.sort_by(f64::total_cmp);
    let reasons: Vec<String> = ABSENT_REASONS
        .iter()
        .map(|&reason| {
            let count = errors
                .iter()
                .filter(|error| error.reason() == Some(reason))
                .count();
            (count, reason)
        })
        .filter(|&(count, _)| count > 0)
        .map(|(count, reason)| format!("{count} {reason}"))
        .collect();
    let reasons = reasons.join(", ");
    match values.as_slice() {
        [] => format!("none has a standard error ({reasons})"),
        some => {
            let spread = format!(
                "median standard error {:.2e}, largest {:.2e}",
                some[(some.len() - 1) / 2],
                some[some.len() - 1]
            );
            match errors.len() - some.len() {
                0 => format!("{spread}, none missing"),
                missing => format!("{spread}, {missing} missing ({reasons})"),
            }
        }
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

/// The interval the fit keeps a sample's own parameter in, by its slot ([`sample`]).
fn bounds_of_own(slot: usize) -> (f64, f64) {
    match sample::class_of(slot) {
        None => HOM_EXCESS_BOUNDS,
        Some(0) => CLEAN_ERROR_BOUNDS,
        Some(_) => NOISY_ERROR_BOUNDS,
    }
}

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
#[cfg(test)]
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
    use crate::parameter_estimation::joint::fit::information::{
        ONE_LIBRARY_SAMPLE_PARAMETERS, cohort, own_parameters,
    };

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

    /// Blocks of a random arrow over `samples` samples of one library each, and the whole arrow
    /// beside them ([`an_arrow_of`]).
    fn an_arrow(samples: usize, seed: u64) -> (InformationSums, Vec<f64>) {
        an_arrow_of(&vec![1; samples], seed)
    }

    /// Which read groups each sample reads, for samples of `libraries` libraries each, numbered in
    /// sample order.
    fn group_index_of(libraries: &[usize]) -> Vec<Vec<usize>> {
        let mut next = 0;
        libraries
            .iter()
            .map(|&count| {
                let own: Vec<usize> = (next..next + count).collect();
                next += count;
                own
            })
            .collect()
    }

    /// Blocks of a random arrow over samples of `libraries` libraries each, and the whole arrow
    /// beside them — the cohort's eight, then each sample's own parameters in turn: a random
    /// symmetric matrix of the full size with every entry between two different samples set to
    /// zero, and each diagonal entry then raised above the sum of its row's other entries' sizes —
    /// which makes it positive definite.
    fn an_arrow_of(libraries: &[usize], seed: u64) -> (InformationSums, Vec<f64>) {
        let sizes: Vec<usize> = libraries.iter().map(|&k| own_parameters(k)).collect();
        let n = COHORT_PARAMETERS + sizes.iter().sum::<usize>();
        let mut full = positive_definite(n, seed);
        let firsts: Vec<usize> = sizes
            .iter()
            .scan(COHORT_PARAMETERS, |at, &size| {
                let first = *at;
                *at += size;
                Some(first)
            })
            .collect();
        let owner = |i: usize| {
            (i >= COHORT_PARAMETERS)
                .then(|| firsts.iter().rposition(|&first| first <= i).expect("owned"))
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
        let mut sums = InformationSums::new(&group_index_of(libraries));
        for row in 0..COHORT_PARAMETERS {
            for column in 0..COHORT_PARAMETERS {
                sums.cohort[row * COHORT_PARAMETERS + column] = full[row * n + column];
            }
        }
        for (s, (&first, &size)) in firsts.iter().zip(&sizes).enumerate() {
            for row in 0..size {
                for column in 0..size {
                    sums.sample_blocks[s][row * size + column] =
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

    /// The error of layout parameter `i` — the cohort's eight, then each sample's own in turn.
    fn error_at(errors: &StandardErrors, i: usize) -> StandardError {
        if i < COHORT_PARAMETERS {
            return errors.cohort[i];
        }
        let mut at = i - COHORT_PARAMETERS;
        for row in &errors.samples {
            if at < row.len() {
                return row[at];
            }
            at -= row.len();
        }
        panic!("parameter {i} is past the layout")
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let errors = StandardErrors::of(&sums);
        assert_errors_match(&errors, &dense_errors_without(&full, n, |_| false));
    }

    /// **The same holds where samples read different numbers of libraries**: samples of one, two,
    /// three and one library — blocks of 3, 5, 7 and 3 own parameters, 26 with the cohort's — every
    /// error against the dense inverse to 10⁻¹² relative. A sample's block read at another sample's
    /// size, or a later library's error written to the first library's slot, fails. **And each
    /// slot is judged against its own parameter's interval**: a second library's rate whose error
    /// comes out at 0.3 is an error as a mismapped rate and wider than the whole range as a clean
    /// one.
    #[test]
    fn the_blocks_give_the_errors_of_an_arrow_of_mixed_libraries() {
        let libraries = [1, 2, 3, 1];
        let (sums, full) = an_arrow_of(&libraries, 13);
        let n = COHORT_PARAMETERS + libraries.iter().map(|&k| own_parameters(k)).sum::<usize>();
        assert_eq!(n, COHORT_PARAMETERS + 18);
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.samples.iter().map(Vec::len).collect::<Vec<_>>(),
            [3, 5, 7, 3]
        );
        assert_errors_match(&errors, &dense_errors_without(&full, n, |_| false));

        // **Each slot is judged against its own parameter's interval.** Sample 1's second library's
        // rates rescaled so each error comes out at 0.3: inside the mismapped rate's interval
        // (0.0001 to 0.45, width 0.4499) and outside the clean rate's (0.000001 to 0.2).
        for (class, expect_estimated) in [(1, true), (0, false)] {
            let slot = sample::rate(1, class);
            let size = sums.own_parameters_of(1);
            let factor = 0.3 / errors.samples[1][slot].value().expect("an error");
            let mut rescaled = sums.clone();
            for j in 0..size {
                rescaled.sample_blocks[1][slot * size + j] /= factor;
                rescaled.sample_blocks[1][j * size + slot] /= factor;
            }
            for column in 0..COHORT_PARAMETERS {
                rescaled.sample_cohort_blocks[1][slot * COHORT_PARAMETERS + column] /= factor;
            }
            let wider = StandardErrors::of(&rescaled).samples[1][slot];
            match (wider, expect_estimated) {
                (StandardError::Estimated(error), true)
                | (StandardError::WiderThanItsRange(error), false) => {
                    assert!((error - 0.3).abs() < 1e-9, "rate {class}: {error}");
                }
                (other, _) => panic!("rate {class} of a second library at 0.3: {other:?}"),
            }
        }
    }

    /// **A parameter with no information has no error, and does not disturb the others'**: the
    /// duplicated class's three rows zeroed, and sample 1's whole block zeroed. Those say
    /// `NoInformation`; every other error equals the one the arrow gives with them removed.
    #[test]
    fn a_parameter_with_no_information_has_no_error() {
        let samples = 3;
        let (mut sums, full) = an_arrow(samples, 23);
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let removed = |i: usize| {
            (cohort::DUPLICATED_SHARE..COHORT_PARAMETERS).contains(&i)
                || (COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS
                    ..COHORT_PARAMETERS + 2 * ONE_LIBRARY_SAMPLE_PARAMETERS)
                    .contains(&i)
        };
        for i in cohort::DUPLICATED_SHARE..COHORT_PARAMETERS {
            for j in 0..COHORT_PARAMETERS {
                sums.cohort[i * COHORT_PARAMETERS + j] = 0.0;
                sums.cohort[j * COHORT_PARAMETERS + i] = 0.0;
            }
            for s in 0..samples {
                for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + i] = 0.0;
                }
            }
        }
        sums.sample_blocks[1].fill(0.0);
        sums.sample_cohort_blocks[1].fill(0.0);
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS;
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
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
            for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
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
            for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
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
        for i in
            (0..COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples).filter(|&i| i != at)
        {
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let (clean, noisy) = (COHORT_PARAMETERS, COHORT_PARAMETERS + 1);
        for j in 0..n {
            full[noisy * n + j] = full[clean * n + j];
        }
        for j in 0..n {
            full[j * n + noisy] = full[j * n + clean];
        }
        for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
            for column in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                sums.sample_blocks[0][row * ONE_LIBRARY_SAMPLE_PARAMETERS + column] =
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
        for block in &mut sums.sample_cohort_blocks {
            block.fill(0.0);
        }
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.cohort,
            [StandardError::NoInformation; COHORT_PARAMETERS]
        );
        for (s, own) in sums.sample_blocks.iter().enumerate() {
            let inverse =
                inverse_of_positive_definite(own, ONE_LIBRARY_SAMPLE_PARAMETERS).expect("inverts");
            for j in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                let expected = inverse[j * ONE_LIBRARY_SAMPLE_PARAMETERS + j].sqrt();
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let (share, rate) = (
            cohort::P_INVARIANT,
            COHORT_PARAMETERS + sample::NOISY_ERROR_RATE,
        );
        for j in 0..COHORT_PARAMETERS {
            sums.cohort[share * COHORT_PARAMETERS + j] = 0.0;
            sums.cohort[j * COHORT_PARAMETERS + share] = 0.0;
        }
        for s in 0..samples {
            for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + share] = 0.0;
            }
        }
        let which = sample::NOISY_ERROR_RATE;
        for j in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
            sums.sample_blocks[0][which * ONE_LIBRARY_SAMPLE_PARAMETERS + j] = 0.0;
            sums.sample_blocks[0][j * ONE_LIBRARY_SAMPLE_PARAMETERS + which] = 0.0;
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let at = cohort::DENSITY_B;
        sums.cohort[at * COHORT_PARAMETERS + at] = f64::INFINITY;
        let excess = sample::HOMOZYGOTE_EXCESS;
        sums.sample_blocks[1][excess * ONE_LIBRARY_SAMPLE_PARAMETERS + excess] = f64::NAN;
        let errors = StandardErrors::of(&sums);
        assert_eq!(errors.cohort[at], StandardError::NoInformation);
        assert_eq!(errors.samples[1][excess], StandardError::NoInformation);
        let excess_of_one = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS + excess;
        assert_errors_match(
            &errors,
            &dense_errors_without(&full, n, |i| i == at || i == excess_of_one),
        );
    }

    /// A fit's parameters over `groups` read groups and `samples` samples, with or without the
    /// duplicated class; the cohort's values all differ.
    fn parameters_for_the_log(duplicated: bool, groups: usize, samples: usize) -> Parameters {
        use crate::parameter_estimation::joint::fit::{DuplicatedPositions, FrequencyDensity};
        Parameters {
            clean: vec![0.001; groups],
            noisy: vec![0.05; groups],
            noisy_share: 0.0312,
            density: FrequencyDensity {
                p_invariant: 0.9,
                p_fixed_alt: 0.0125,
                a: 0.5,
                b: 61.0,
            },
            hom_excess: vec![0.1; samples],
            duplicated: duplicated.then_some(DuplicatedPositions {
                share: 0.01,
                carrier_a: 1.5,
                carrier_b: 5.0,
            }),
        }
    }

    /// **The run's log line gives every cohort parameter's value and standard error by name, and
    /// summarises each kind of a sample's own**: the median (the lower middle of an even count) and
    /// the largest standard error, and how many are missing, by reason — with values chosen so a
    /// median taken as the mean, the upper middle, or the largest taken as the last listed, prints a
    /// different line. The duplicated class's three appear only when the run fits it. **The rates
    /// are counted a read group each**: with sample 0 read from two libraries, its second library's
    /// two rates join the counts, from their own slots.
    #[test]
    fn the_logged_line_summarises_each_kind() {
        use StandardError::{
            Estimated, HeldFixed, NoInformation, NotIdentified, WiderThanItsRange,
        };
        let cohort = [
            Estimated(2e-4),
            Estimated(3.25e-3),
            NotIdentified,
            WiderThanItsRange(61.0),
            Estimated(0.5),
            Estimated(1e-3),
            NoInformation,
            NotIdentified,
        ];
        let errors = StandardErrors {
            cohort,
            samples: vec![
                vec![Estimated(4e-4), Estimated(0.01), HeldFixed],
                vec![Estimated(1e-4), NoInformation, HeldFixed],
                vec![Estimated(3e-4), NotIdentified, HeldFixed],
                vec![Estimated(2e-4), Estimated(0.02), HeldFixed],
            ],
        };
        let one_group_each = vec![vec![0], vec![1], vec![2], vec![3]];
        assert_eq!(
            errors.described(&parameters_for_the_log(false, 4, 4), &one_group_each),
            "mismapped share 3.1200e-2 ± 2.00e-4, invariant share 9.0000e-1 ± 3.25e-3, fixed \
             non-reference share 1.2500e-2, no standard error (not identified), allele-frequency \
             shape a 5.0000e-1, no standard error (it came out at 6.10e1, wider than the \
             parameter's whole range), allele-frequency shape b 6.1000e1 ± 5.00e-1; error rates at \
             ordinary positions (4 read groups): median standard error 2.00e-4, largest 4.00e-4, \
             none missing; error rates at mismapped positions (4 read groups): median standard \
             error 1.00e-2, largest 2.00e-2, 2 missing (1 no information, 1 not identified); \
             homozygote excesses (4 samples): none has a standard error (4 held fixed)"
        );
        let with_the_class = errors.described(&parameters_for_the_log(true, 4, 4), &one_group_each);
        assert!(
            with_the_class.contains(
                "allele-frequency shape b 6.1000e1 ± 5.00e-1, duplicated share 1.0000e-2 ± 1.00e-3, \
                 carrier-frequency shape a 1.5000e0, no standard error (no information), \
                 carrier-frequency shape b 5.0000e0, no standard error (not identified); error \
                 rates at ordinary positions"
            ),
            "{with_the_class}"
        );
        let mut two_libraries = errors.clone();
        two_libraries.samples[0].extend([Estimated(5e-4), NoInformation]);
        let second_group_of_sample_0 = vec![vec![0, 2], vec![1], vec![3], vec![4]];
        let line = two_libraries.described(
            &parameters_for_the_log(false, 5, 4),
            &second_group_of_sample_0,
        );
        assert!(
            line.ends_with(
                "error rates at ordinary positions (5 read groups): median standard error \
                 3.00e-4, largest 5.00e-4, none missing; error rates at mismapped positions (5 read \
                 groups): median standard error 1.00e-2, largest 2.00e-2, 3 missing (2 no \
                 information, 1 not identified); homozygote excesses (4 samples): none has a \
                 standard error (4 held fixed)"
            ),
            "{line}"
        );
    }

    /// **The trace names each error after the value it belongs to**: the same names as
    /// [`Parameters::named`](crate::parameter_estimation::joint::fit::Parameters), in the same order,
    /// with and without the duplicated class. Two samples over four read groups, the first sample
    /// holding groups 0 and 2 and no sample group 3: group 0 carries sample 0's first library's
    /// errors, group 2 its second library's (from the row's fourth and fifth slots), group 1 sample
    /// 1's, and group 3 NaN.
    #[test]
    fn the_trace_names_each_error_after_its_value() {
        use StandardError::{Estimated, HeldFixed, NoInformation};
        let errors = StandardErrors {
            cohort: [Estimated(1.0); COHORT_PARAMETERS],
            samples: vec![
                vec![
                    Estimated(0.1),
                    Estimated(0.2),
                    Estimated(0.3),
                    Estimated(0.5),
                    Estimated(0.6),
                ],
                vec![Estimated(0.4), NoInformation, HeldFixed],
            ],
        };
        let group_index = vec![vec![0, 2], vec![1]];
        for duplicated in [false, true] {
            let parameters = parameters_for_the_log(duplicated, 4, 2);
            let named = errors.named(&parameters, &group_index);
            let expected: Vec<String> = parameters
                .named()
                .into_iter()
                .map(|(name, _)| format!("standard_error:{name}"))
                .collect();
            let got: Vec<String> = named.iter().map(|(name, _)| name.clone()).collect();
            assert_eq!(got, expected);
            let value_of = |name: &str| {
                named
                    .iter()
                    .find(|(n, _)| n == &format!("standard_error:{name}"))
                    .map(|(_, value)| *value)
                    .expect("named")
            };
            assert_eq!(value_of("clean_error_0"), 0.1);
            assert_eq!(value_of("noisy_error_0"), 0.2);
            assert_eq!(value_of("clean_error_2"), 0.5);
            assert_eq!(value_of("noisy_error_2"), 0.6);
            assert!(value_of("clean_error_3").is_nan());
            assert!(value_of("noisy_error_3").is_nan());
            assert_eq!(value_of("clean_error_1"), 0.4);
            assert!(value_of("noisy_error_1").is_nan());
            assert_eq!(value_of("hom_excess_0"), 0.3);
            assert!(value_of("hom_excess_1").is_nan());
        }
    }

    /// **Each slot of a sample's row is the parameter the layout says, with that parameter's
    /// bounds**: for samples of one to four libraries, every library's two rates land in distinct
    /// slots that, with the homozygote excess, fill the row exactly; a rate's slot names its own noise
    /// class and carries that class's bounds, and the excess's slot names none and carries the
    /// excess's. A further library's mismapped rate given the clean rate's bounds fails it.
    #[test]
    fn each_slot_of_a_samples_row_is_its_own_parameter_with_its_own_bounds() {
        for libraries in 1..=4 {
            let n = own_parameters(libraries);
            let mut filled = vec![false; n];
            filled[sample::HOMOZYGOTE_EXCESS] = true;
            assert_eq!(sample::class_of(sample::HOMOZYGOTE_EXCESS), None);
            assert_eq!(bounds_of_own(sample::HOMOZYGOTE_EXCESS), HOM_EXCESS_BOUNDS);
            for library in 0..libraries {
                for (class, bounds) in [(0, CLEAN_ERROR_BOUNDS), (1, NOISY_ERROR_BOUNDS)] {
                    let slot = sample::rate(library, class);
                    assert!(
                        slot < n && !filled[slot],
                        "{libraries} libraries: library {library}'s rate {class} at slot {slot}"
                    );
                    filled[slot] = true;
                    assert_eq!(sample::class_of(slot), Some(class), "slot {slot}");
                    assert_eq!(bounds_of_own(slot), bounds, "slot {slot}");
                }
            }
            assert!(
                filled.iter().all(|&slot| slot),
                "{libraries} libraries: {filled:?}"
            );
        }
        assert_eq!(own_parameters(0), ONE_LIBRARY_SAMPLE_PARAMETERS);
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
            for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                scaled.sample_cohort_blocks[s][row * COHORT_PARAMETERS + at] *= scale;
            }
        }
        let (plain, small) = (StandardErrors::of(&sums), StandardErrors::of(&scaled));
        for i in 0..COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples {
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
