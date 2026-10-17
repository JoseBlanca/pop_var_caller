//! Each SNP/indel parameter's **standard error**, from the information a pass summed — in blocks,
//! and as the whole matrix for a small cohort.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §3.2–3.3. Build order:
//! `doc/devel/implementation_plans/fit_precision.md`, steps A3, A5, A8 and B1.
//!
//! Every fit computes them at the parameters it returns, from the information its final pass
//! sums ([`fit_jointly`](super::fit_jointly)); it prints them as one line of the run's log
//! ([`StandardErrors::described`]) and, when the per-pass trace is on, one row a parameter
//! ([`StandardErrors::named`]). While it runs, once its log-likelihood has nearly stopped moving, it
//! computes them every cycle beside each parameter's distance to the maximum ([`newton_step`], the
//! same matrix solved against the summed scores), and stops when every distance is below a tenth of
//! its error ([`settled`](super::settled), [`StandardErrors::by_coordinate`]).
//!
//! # From information to errors
//!
//! A parameter's standard error is how far its estimate would typically move if the same kind of
//! data were drawn again. It is the square root of the parameter's entry on the diagonal of the
//! **inverse** of the information matrix — not one over its own diagonal entry, because a parameter
//! whose effect another parameter can mimic is less well determined than its own curvature says.
//!
//! **Two ways, chosen by the cohort's size** ([`StandardErrors::of`]):
//!
//! - **the whole matrix**, for a cohort of at most
//!   [`FULL_MATRIX_SAMPLES`](super::information::FULL_MATRIX_SAMPLES) samples and
//!   [`FULL_MATRIX_PARAMETERS`](super::information::FULL_MATRIX_PARAMETERS) parameters: every
//!   parameter's score paired with every other's, two samples' included ([`FullInformation`]),
//!   inverted at once;
//! - **the blocks**, above that ([`InformationSums`]): the cohort's eight parameters with each other
//!   (`C`), each sample's own with each other (`A_s`) — three for a sample of one library, `1 + 2k`
//!   for one of `k` — and each sample's own with the cohort's (`B_s`). Two samples' parameters are
//!   not paired, so the matrix is an **arrow**: the cohort's row and column run along one edge and
//!   every sample's block sits on the diagonal, alone. An arrow is inverted exactly, block by block,
//!   without ever forming it:
//!   - the cohort's errors come from `C − Σ_s B_sᵀ A_s⁻¹ B_s` inverted — the cohort's information,
//!     less what each sample's own parameters could explain of it;
//!   - a sample's errors come from `A_s⁻¹ + A_s⁻¹ B_s V B_sᵀ A_s⁻¹`, where `V` is the cohort's
//!     inverse just computed — the sample's own uncertainty, plus the cohort's carried through the
//!     parameters the two share.
//!
//! What the arrow leaves out is the pairing of two samples' scores, and **at small cohorts that
//! matters**: at 4 samples and 3 reads a position the error rates and the mismapped share scattered
//! 1.23 to 1.67 times the blocks' errors and 0.89 to 1.02 times the whole matrix's (plan step A4).
//! At 20 samples the two came within 11% of each other (1.067 against 0.965, the mismapped rates at
//! 3 reads), and the whole matrix grows as the square of its parameters, so above 20 samples, or
//! 188 parameters, the blocks are kept.
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
//!   lives in [0, 1], and the density's shapes' at 720 and 1,148 though they live in [0.02, 50].
//!   That was on the blocks; a two-sample cohort now takes the whole matrix, which, with the
//!   products of the two samples' scores in, cannot tell those parameters apart at all and drops
//!   them as not identified.
//!
//! A parameter counts as not told apart when the curvature left to it once the parameters before
//! it are accounted for is below [`IDENTIFIED_SHARE`] of its own curvature — which parameter of a
//! mimicking set that is depends on their order: each sample's own first, the samples in turn,
//! then the cohort's, on both ways. On an arrow the two ways drop the same parameters. On the whole
//! matrix a sample's parameter is also judged against earlier samples', so where a parameter of one
//! sample mimics one of another, the later sample's is the one dropped. Without the threshold, a
//! remainder that is only rounding would be inverted into an error of millions (measured: 2 × 10⁶
//! on a carrier shape, 13.5 on a share that lives in [0, 1], at two samples and 300,000 positions).

use crate::float;

use super::determined::{DensityReport, DeterminedShare, DuplicatedReport, STATED_FREQUENCY};
use super::information::{
    COHORT_PARAMETERS, FullInformation, InformationSums, cohort, own_parameters, sample,
};
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
pub(super) const COHORT_PARAMETER_NAMES: [&str; COHORT_PARAMETERS] = [
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

/// Whose parameter a row of the whole matrix is: the cohort's, or one sample's.
#[derive(Copy, Clone, Debug)]
enum Owner {
    Cohort,
    Sample(usize),
}

/// A parameter the whole matrix is inverted over: whose it is, its slot among its owner's
/// parameters (a [`cohort`] slot or a [`sample`] slot), and its row in the matrix.
#[derive(Copy, Clone, Debug)]
struct InvertedParameter {
    owner: Owner,
    slot: usize,
    row: usize,
}

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
    /// being reported: from the whole matrix when the pass kept it — a cohort of at most
    /// [`FULL_MATRIX_SAMPLES`](super::information::FULL_MATRIX_SAMPLES) samples and
    /// [`FULL_MATRIX_PARAMETERS`](super::information::FULL_MATRIX_PARAMETERS) parameters — and from
    /// the blocks otherwise.
    pub(super) fn of(sums: &InformationSums) -> Self {
        match &sums.full {
            Some(full) => Self::of_the_whole_matrix(full),
            None => Self::of_the_blocks(sums),
        }
    }

    /// **The errors from the whole matrix**, two samples' parameters paired included: the matrix
    /// inverted at once over every parameter that has information and that the fit moves.
    ///
    /// The parameters are taken in the order the blocks take them — each sample's own, the samples
    /// in turn, then the cohort's — so a parameter the others mimic is dropped by the same rule
    /// ([`invert_identified`]): a sample's own is judged against the parameters before it, the
    /// cohort's against everything the samples explain. On an arrow, where no two samples' scores
    /// are paired, the two ways drop the same parameters and give the same errors.
    fn of_the_whole_matrix(full: &FullInformation) -> Self {
        let samples = full.samples();
        let excess_is_fitted = fits_homozygote_excess(samples_informing_their_excess_in(full));
        let mut errors = Self {
            cohort: [StandardError::NoInformation; COHORT_PARAMETERS],
            samples: (0..samples)
                .map(|s| vec![StandardError::NoInformation; full.rows_of_sample(s).len()])
                .collect(),
        };
        let mut inverted: Vec<InvertedParameter> = Vec::with_capacity(full.side());
        for s in 0..samples {
            for (slot, row) in full.rows_of_sample(s).enumerate() {
                if !has_information(full.entry(row, row)) {
                    continue;
                }
                if slot == sample::HOMOZYGOTE_EXCESS && !excess_is_fitted {
                    errors.samples[s][slot] = StandardError::HeldFixed;
                    continue;
                }
                inverted.push(InvertedParameter {
                    owner: Owner::Sample(s),
                    slot,
                    row,
                });
            }
        }
        for slot in 0..COHORT_PARAMETERS {
            if has_information(full.entry(slot, slot)) {
                inverted.push(InvertedParameter {
                    owner: Owner::Cohort,
                    slot,
                    row: slot,
                });
            }
        }
        let rows: Vec<usize> = inverted.iter().map(|parameter| parameter.row).collect();
        let matrix = square_of(&rows, |row, column| full.entry(row, column));
        let reference: Vec<f64> = rows.iter().map(|&row| full.entry(row, row)).collect();
        let identified = invert_identified(&matrix, rows.len(), &reference);
        let mut record = |parameter: InvertedParameter, error: StandardError| match parameter.owner
        {
            Owner::Sample(s) => errors.samples[s][parameter.slot] = error,
            Owner::Cohort => errors.cohort[parameter.slot] = error,
        };
        for &position in &identified.dropped {
            record(inverted[position], StandardError::NotIdentified);
        }
        let kept = identified.kept.len();
        for (p, &position) in identified.kept.iter().enumerate() {
            let parameter = inverted[position];
            let bounds = match parameter.owner {
                Owner::Sample(_) => bounds_of_own(parameter.slot),
                Owner::Cohort => COHORT_BOUNDS[parameter.slot],
            };
            record(
                parameter,
                error_from_variance(identified.inverse[p * kept + p], bounds),
            );
        }
        errors
    }

    /// **The errors from the blocks**: the arrow inverted block by block (module doc).
    pub(super) fn of_the_blocks(sums: &InformationSums) -> Self {
        let samples = sums.sample_blocks.len();
        let excess_is_fitted = fits_homozygote_excess(samples_informing_their_excess(sums));
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
        let errors = self.by_coordinate(parameters, group_index);
        let names = parameters.named();
        assert_eq!(
            names.len(),
            errors.len(),
            "one error for every parameter the trace names"
        );
        names
            .into_iter()
            .zip(errors)
            .map(|((name, _), error)| {
                (
                    format!("standard_error:{name}"),
                    error.value().unwrap_or(f64::NAN),
                )
            })
            .collect()
    }

    /// **Every parameter's error in the order the fit's own vector lists them**
    /// ([`Parameters::coordinates`], which [`Parameters::named`] follows): the cohort's fitted
    /// parameters, each read group's two rates, then each sample's homozygote excess. `group_index`
    /// says which read groups each sample's sections are, as the pass reads them; a read group no
    /// sample holds has no information.
    pub(super) fn by_coordinate(
        &self,
        parameters: &Parameters,
        group_index: &[Vec<usize>],
    ) -> Vec<StandardError> {
        debug_assert!(
            self.laid_out_for(group_index),
            "rows laid out for another cohort"
        );
        in_coordinate_order(
            &self.cohort,
            &self.samples,
            parameters,
            group_index,
            StandardError::NoInformation,
        )
    }

    /// **The errors as one line of the run's log**: each cohort-level parameter's value and
    /// standard error, or why it has none; for each library's two rates and each sample's
    /// homozygote excess, the median and the largest standard error and how many have none, by
    /// reason. `group_index` says which read groups each sample's sections are. The duplicated
    /// class's three are left out when the run does not fit it. The median of an even count is the
    /// lower of the middle two.
    ///
    /// **What the data determine** (`parameter_prepass_joint_fit.md` §2.1.3), from `density` and,
    /// when the run fits the class, `duplicated`: where a pair's two numbers are not separately
    /// determined, each is printed with that verdict in place of its error, the first with the
    /// share the data do fix and what trading the pair costs; where they are, that share and the
    /// cost follow the pair's last number.
    pub(super) fn described(
        &self,
        parameters: &Parameters,
        group_index: &[Vec<usize>],
        density: &DensityReport,
        duplicated: Option<&DuplicatedReport>,
    ) -> String {
        debug_assert!(
            self.laid_out_for(group_index),
            "rows laid out for another cohort"
        );
        let values = parameters.named();
        let fitted = cohort_parameters_fitted(parameters);
        let mut in_place: [Option<String>; COHORT_PARAMETERS] = Default::default();
        let mut after: [Option<String>; COHORT_PARAMETERS] = Default::default();
        let pairs = [
            Some((&DENSITY_PAIR, density.determined)),
            duplicated
                .filter(|_| fitted == COHORT_PARAMETERS)
                .map(|report| (&DUPLICATED_PAIR, report.determined)),
        ];
        for (pair, determined) in pairs.into_iter().flatten() {
            if determined.separately_determined {
                after[pair.last_of_group] = Some(pair.determined_share_entry(&determined));
            } else {
                in_place[pair.first] = Some(pair.verdict_on_the_first(&determined));
                in_place[pair.second] = Some(pair.verdict_on_the_second());
            }
        }
        let mut cohort_entries = Vec::with_capacity(fitted + 2);
        for i in 0..fitted {
            let name = COHORT_PARAMETER_NAMES[i];
            let value = values[i].1;
            cohort_entries.push(match &in_place[i] {
                Some(verdict) => format!("{name} {value:.4e} ({verdict})"),
                None => format!("{name} {}", self.cohort[i].printed_with(value)),
            });
            cohort_entries.extend(after[i].take());
        }
        let mut parts = vec![cohort_entries.join(", ")];
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

/// **How the run's log speaks of one pair of numbers the data may not separate**
/// (`parameter_prepass_joint_fit.md` §2.1.3): which two cohort slots the pair is, after which slot
/// a separately determined pair's share is printed, and the words for the share the trade moves
/// and the share the data see.
struct PairWording {
    first: usize,
    second: usize,
    last_of_group: usize,
    /// The share the trade moves, `σ`.
    traded_share: &'static str,
    /// The share the data see, as a noun phrase ending where the stated frequency is printed.
    share_seen: &'static str,
    /// What follows the stated frequency in [`Self::share_seen`].
    share_seen_ends: &'static str,
    /// What the data do not fix.
    not_seen: &'static str,
}

/// The density's pair: the invariant share against the first shape.
const DENSITY_PAIR: PairWording = PairWording {
    first: cohort::P_INVARIANT,
    second: cohort::DENSITY_A,
    last_of_group: cohort::DENSITY_B,
    traded_share: "the segregating share",
    share_seen: "the share segregating above",
    share_seen_ends: "",
    not_seen: "not how much lies below it",
};

/// The duplicated class's pair: its share against the carrier Beta's first shape.
const DUPLICATED_PAIR: PairWording = PairWording {
    first: cohort::DUPLICATED_SHARE,
    second: cohort::CARRIER_A,
    last_of_group: cohort::CARRIER_B,
    traded_share: "the duplicated share",
    share_seen: "the share of positions duplicated in more than",
    share_seen_ends: " of the panel",
    not_seen: "not how many are duplicated in fewer",
};

impl PairWording {
    /// The share the data see, named with the stated frequency: "the share segregating above 1%".
    fn share_seen(&self) -> String {
        format!(
            "{} {}%{}",
            self.share_seen,
            STATED_FREQUENCY * 100.0,
            self.share_seen_ends
        )
    }

    /// **The first number's verdict when the pair is not separately determined**: which number it
    /// trades against, the share the data fix with its error, and what trading costs.
    fn verdict_on_the_first(&self, determined: &DeterminedShare) -> String {
        let cost = match determined.trade_cost {
            Some(cost) => format!(
                "halving or doubling {} costs {}",
                self.traded_share,
                a_cost(cost)
            ),
            None => format!(
                "what halving or doubling {} costs cannot be computed, as the data do not \
                 determine every number it is computed from",
                self.traded_share
            ),
        };
        format!(
            "not separately determined from {}: the data fix {}, {}, {}; {cost}",
            COHORT_PARAMETER_NAMES[self.second],
            self.share_seen(),
            printed_share(determined),
            self.not_seen
        )
    }

    /// The second number's verdict when the pair is not separately determined.
    fn verdict_on_the_second(&self) -> String {
        format!(
            "not separately determined from the {}",
            COHORT_PARAMETER_NAMES[self.first]
        )
    }

    /// **The share the data see, as its own entry, when the pair is separately determined**, with
    /// what trading the traded share against it would cost.
    fn determined_share_entry(&self, determined: &DeterminedShare) -> String {
        let cost = match determined.trade_cost {
            Some(cost) if cost.is_infinite() => {
                format!("{} cannot change without changing it", self.traded_share)
            }
            Some(cost) => format!(
                "halving or doubling {} with it held costs {}",
                self.traded_share,
                a_cost(cost)
            ),
            None => "what trading them costs cannot be computed".to_owned(),
        };
        let share_seen = self.share_seen();
        format!(
            "{} {} (the {} and {} are separately determined: {cost})",
            share_seen.trim_start_matches("the "),
            printed_share(determined),
            COHORT_PARAMETER_NAMES[self.first],
            COHORT_PARAMETER_NAMES[self.second],
        )
    }
}

/// The share the data see and its error, as the run's log prints them: `8.45e-3 ± 2.0e-4`.
fn printed_share(determined: &DeterminedShare) -> String {
    let share = determined.share_above;
    match share.standard_error {
        Some(error) => format!("{:.2e} ± {error:.1e}", share.value),
        None => format!("{:.2e} with no standard error", share.value),
    }
}

/// **A log-likelihood cost as the run's log prints it**: in three significant figures, and above a
/// million as just that. A cost that large is rounding: where the two shares move together, as on
/// the 4 oracle samples, `V_trade` comes out a unit in the last place above zero and the cost at
/// 10¹⁶ or so, and printing that figure would invite a reader to make something of it.
fn a_cost(cost: f64) -> String {
    if cost > 1e6 {
        "over a million log-likelihood units".to_owned()
    } else {
        format!("{} log-likelihood units", in_three_figures(cost))
    }
}

/// A log-likelihood cost in three significant figures — `1.35`, `0.00549`, `324` — and in
/// scientific notation from a thousand up, and for anything not above zero.
fn in_three_figures(cost: f64) -> String {
    if cost.is_nan() || cost <= 0.0 || cost >= 1_000.0 {
        return format!("{cost:.2e}");
    }
    let decimals = (2.0 - float::log10(cost).floor()).max(0.0) as usize;
    format!("{cost:.decimals$}")
}

/// **Per-parameter values in the order the fit's own vector lists them**
/// ([`Parameters::coordinates`]), from the cohort's eight slots and each sample's own row in the
/// [`sample`] layout: the cohort's fitted parameters, each read group's two rates, then each
/// sample's homozygote excess. A read group no sample holds gets `absent`.
fn in_coordinate_order<T: Copy>(
    cohort: &[T; COHORT_PARAMETERS],
    samples: &[Vec<T>],
    parameters: &Parameters,
    group_index: &[Vec<usize>],
    absent: T,
) -> Vec<T> {
    let mut ordered: Vec<T> = cohort[..cohort_parameters_fitted(parameters)].to_vec();
    // Which sample holds each read group, and as which of its libraries.
    let mut held_by = vec![None; parameters.clean.len()];
    for (s, own) in group_index.iter().enumerate() {
        for (section, &group) in own.iter().enumerate() {
            held_by[group] = Some((s, section));
        }
    }
    for holder in held_by {
        ordered.extend(match holder {
            Some((s, section)) => [
                samples[s][sample::rate(section, 0)],
                samples[s][sample::rate(section, 1)],
            ],
            None => [absent; 2],
        });
    }
    ordered.extend(samples.iter().map(|row| row[sample::HOMOZYGOTE_EXCESS]));
    ordered
}

/// **How far each parameter is from the likelihood's maximum, as a Newton step estimates it**, in
/// the order the fit's vector lists them ([`Parameters::coordinates`]), on each parameter's own
/// scale; `None` where the parameter is not solved for — no information, held fixed by the fit, or
/// not told apart from the others.
///
/// **The step is the one to the maximum of the quadratic model within the parameters' bounds**
/// ([`bounded_step`]): the information `I` and the summed scores `g` of one pass describe the
/// log-likelihood near the parameters as `g·d − ½ dᵀ I d`, and a parameter whose maximum lies past
/// one end of the interval the fit keeps it in stops at that end, its distance the way there. Without
/// the bounds the step is `I⁻¹ g`, which carries a parameter the fit is walking towards a bound far
/// past it and, through the parameters it is correlated with, the others' steps too.
///
/// The same matrix the errors come from ([`StandardErrors::of`]): the whole one where the pass kept
/// it, the arrow of blocks otherwise, and the same parameters left out of it.
///
/// **It also says what the step is worth** ([`NewtonStep::decrement`]), from the same solve and
/// `errors`, the standard errors the same pass's sums give ([`StandardErrors::of`]): the stopping
/// rule reads that, and the distances are what the step moves.
pub(super) fn newton_step(
    sums: &InformationSums,
    errors: &StandardErrors,
    parameters: &Parameters,
    group_index: &[Vec<usize>],
) -> NewtonStep {
    let layout = Layout::of(sums);
    let values = layout.values_of(parameters, group_index);
    let bounds = layout.bounds();
    let step = bounded_step(sums, &layout, &values, &bounds);
    let has_error: Vec<bool> = errors
        .cohort
        .iter()
        .chain(errors.samples.iter().flatten())
        .map(|error| error.value().is_some())
        .collect();
    assert_eq!(
        has_error.len(),
        layout.side,
        "the errors are laid out as the sums they came from"
    );
    let decrement = if step.finished {
        determined_decrement(sums, &layout, &step.moves, &has_error)
            .unwrap_or_else(|| decrement_of(sums, &layout, &step.moves))
    } else {
        f64::NAN
    };
    let (cohort, samples) = layout.split(&step.moves);
    NewtonStep {
        distances: in_coordinate_order(&cohort, &samples, parameters, group_index, None),
        decrement,
    }
}

/// **One judging pass's Newton step**: how far each parameter is from the maximum, and how much
/// log-likelihood the quadratic model says the whole step would gain.
pub(super) struct NewtonStep {
    /// Each parameter's distance to the maximum of the quadratic model within the bounds, in the
    /// order the fit's vector lists them ([`Parameters::coordinates`]); `None` where the parameter is
    /// not solved for.
    pub(super) distances: Vec<Option<f64>>,
    /// **The Newton decrement, over the directions the data determine** (`fit_precision.md` §2,
    /// amendment of 2026-10-17). When every parameter the step moves has a standard error, it is the
    /// log-likelihood the step would gain under the quadratic model, `g·d − ½ dᵀ I d` for the slope
    /// `g`, the information `I` and the step `d` with each parameter not solved for moving 0 —
    /// `½ gᵀ I⁻¹ g` when no bound binds. When some have none, it is the step's length in standard
    /// errors with those parameters profiled out ([`determined_decrement`]). Never below zero
    /// either way; **not a number when the box solve reached its round limit**, so that no start
    /// stops on a step that is not the box's maximum.
    pub(super) decrement: f64,
}

/// The cohort's eight parameters' covariance, by slot ([`cohort`](super::information::cohort)), row
/// by row; `None` for an entry whose row or column was not solved for.
pub(super) type CohortCovariance = [[Option<f64>; COHORT_PARAMETERS]; COHORT_PARAMETERS];

/// **The covariance of the cohort's eight parameters**: the inverse information restricted to them,
/// every sample's own parameters accounted for (`parameter_prepass_joint_fit.md` §2.1.3). Its
/// diagonal is the square of the cohort's standard errors ([`StandardErrors::of`]); its off-diagonal
/// entries are what a derived number such as the share segregating above 1% needs, since its error
/// depends on how the parameters it is built from move together.
///
/// Column `k` is the solve against a slope of one on parameter `k` and zero elsewhere, by the same
/// solve [`newton_step`] uses — the whole matrix where the pass kept it ([`whole_matrix_solve`]),
/// the arrow otherwise ([`arrow_solve`]) — with nothing held. So the parameters left out are the
/// ones the errors leave out: no information, or not told apart from the others; their rows and
/// columns are `None`. A parameter whose error came out wider than its range still has its entries.
/// Eight solves, each costing what one Newton step's first solve does: linear in the samples.
pub(super) fn cohort_covariance(sums: &InformationSums) -> CohortCovariance {
    let layout = Layout::of(sums);
    let held = vec![false; layout.side];
    let mut covariance = [[None; COHORT_PARAMETERS]; COHORT_PARAMETERS];
    let mut unit_slope = vec![0.0; layout.side];
    for k in 0..COHORT_PARAMETERS {
        unit_slope.fill(0.0);
        unit_slope[k] = 1.0;
        let column = match &sums.full {
            Some(full) => whole_matrix_solve(full, &unit_slope, &held),
            None => arrow_solve(sums, &layout, &unit_slope, &held),
        };
        // A parameter not solved for has no column: the solve ignored its slope. Its row is `None`
        // in every other column, since every solve leaves out the same parameters.
        if column[k].is_none() {
            continue;
        }
        for (row, entry) in covariance.iter_mut().zip(&column) {
            row[k] = *entry;
        }
    }
    covariance
}

/// **The box solve's answer**: each parameter's move in the flat layout, `None` where it was never
/// solved for, and whether the solve reached the box's maximum before its round limit.
struct BoundedStep {
    moves: Vec<Option<f64>>,
    finished: bool,
}

/// How far a fixed parameter's slope must point back into its interval, as a share of the two terms
/// it is the difference of, for the box solve to free it: below this the sign is rounding, and
/// freeing on it could free and re-fix the same parameter until the round limit.
const SLOPE_POINTING_IN: f64 = 1e-9;

/// **What the bounded step `step` is worth under the quadratic model**, `g·d − ½ dᵀ I d`, in the
/// flat layout, a parameter not solved for moving 0.
fn decrement_of(sums: &InformationSums, layout: &Layout, step: &[Option<f64>]) -> f64 {
    let moved: Vec<f64> = step.iter().map(|moved| moved.unwrap_or(0.0)).collect();
    let slope = layout.slope_of(sums);
    let curved = information_times(sums, layout, &moved);
    let gained: f64 = slope.iter().zip(&moved).map(|(g, d)| g * d).sum();
    let lost: f64 = moved.iter().zip(&curved).map(|(d, c)| d * c).sum();
    gained - 0.5 * lost
}

/// **What the bounded step `step` is worth along the directions the data determine**
/// (`fit_precision.md` §2, amendment of 2026-10-17), in the flat layout; `None` when every
/// parameter the step moves has a standard error (`has_error`), where the decrement is
/// [`decrement_of`]'s.
///
/// `½ dᵀ I d` is the step's length in standard errors: no quantity built from the parameters moves
/// by more than `√(dᵀ I d)` of its own error. Split the parameters into those the step moves that
/// have an error, `K`, and those that the data leave without one, `U` — a parameter the step moves
/// with no error, or one with an error the step does not move. This is **the smallest that length
/// can be when the parameters in `K` move as the step says and those in `U` may move anywhere**:
///
/// ```text
/// λ_K = min over e_U of ½ (d_K, e_U)ᵀ I (d_K, e_U) = ½ d_Kᵀ (I_KK − I_KU I_UU⁻¹ I_UK) d_K
/// ```
///
/// the largest `½ (Δh / SE_h)²` over every quantity `h` built from `K` alone. The minimum is at
/// `e_U = −I_UU⁻¹ I_UK d_K`: one solve over `U` of what moving `K` pulls on it, by the same block
/// elimination the step uses on the arrow ([`arrow_solve`]) or the whole matrix
/// ([`whole_matrix_solve`]), so nothing is inverted densely and the cost grows with the samples as
/// the step's does. A parameter in neither set — no information, held fixed, or not told apart from
/// the others — stays put, as it does in the step; so does one of `U` the solve over `U` cannot
/// tell apart from the rest of `U`.
fn determined_decrement(
    sums: &InformationSums,
    layout: &Layout,
    step: &[Option<f64>],
    has_error: &[bool],
) -> Option<f64> {
    let profiled: Vec<bool> = step
        .iter()
        .zip(has_error)
        .map(|(moved, has_error)| moved.is_some() != *has_error)
        .collect();
    if !profiled.contains(&true) {
        return None;
    }
    // The step's determined part, `(d_K, 0)`, and what it pulls on every parameter, `I (d_K, 0)`.
    let determined: Vec<f64> = step
        .iter()
        .zip(has_error)
        .map(|(moved, has_error)| {
            if *has_error {
                moved.unwrap_or(0.0)
            } else {
                0.0
            }
        })
        .collect();
    let pulled = information_times(sums, layout, &determined);
    // `I_UU⁻¹ I_UK d_K`, solved over the profiled parameters alone.
    let held: Vec<bool> = profiled.iter().map(|profiled| !profiled).collect();
    let followed = match &sums.full {
        Some(full) => whole_matrix_solve(full, &pulled, &held),
        None => arrow_solve(sums, layout, &pulled, &held),
    };
    let moves: Vec<f64> = determined
        .iter()
        .zip(&followed)
        .map(|(determined, followed)| determined - followed.unwrap_or(0.0))
        .collect();
    let curved = information_times(sums, layout, &moves);
    Some(0.5 * moves.iter().zip(&curved).map(|(m, c)| m * c).sum::<f64>())
}

/// **The Newton step to the maximum of the quadratic model within the bounds**, by the primal
/// active-set method for a box (`fit_precision.md` §2, amendment of 2026-10-17; Nocedal & Wright,
/// *Numerical Optimization*, §16.5). From staying put, which is inside the box:
///
/// 1. solve the free parameters' Newton step from the current point, the fixed ones held where they
///    are — the slope there is the summed scores less the information times the moves so far;
/// 2. move along it until the first free parameter meets an end of its interval, and fix it there,
///    with any other that meets an end at the same point; then solve again;
/// 3. after a whole step, free the one fixed parameter whose slope points furthest back into its
///    interval ([`SLOPE_POINTING_IN`]), and solve again; when none does, the point is the box's
///    maximum and the solve has finished.
///
/// Every point it passes through is inside the box, and the model's value `g·d − ½ dᵀ I d` never
/// falls along the way, so the step is worth at least nothing. When no bound binds, the first solve
/// is the whole step. It gives up after ten rounds a parameter, plus ten, with the point it has
/// reached, still inside the box, and says so.
fn bounded_step(
    sums: &InformationSums,
    layout: &Layout,
    values: &[f64],
    bounds: &[(f64, f64)],
) -> BoundedStep {
    let slope = layout.slope_of(sums);
    let mut moved = vec![0.0; layout.side];
    // Where each fixed parameter is fixed: the end of its interval, or `None` while it is free.
    let mut fixed_at: Vec<Option<f64>> = vec![None; layout.side];
    // Whether the solve has given the parameter a move at least once.
    let mut solved = vec![false; layout.side];
    let mut finished = false;
    for _ in 0..10 * layout.side + 10 {
        // The slope where the moves so far have taken the parameters, and the free parameters'
        // Newton step from there.
        let curved = information_times(sums, layout, &moved);
        let slope_here: Vec<f64> = slope.iter().zip(&curved).map(|(g, c)| g - c).collect();
        let held: Vec<bool> = fixed_at.iter().map(Option::is_some).collect();
        let towards = match &sums.full {
            Some(full) => whole_matrix_solve(full, &slope_here, &held),
            None => arrow_solve(sums, layout, &slope_here, &held),
        };
        // How far along it each free parameter can go before it meets an end, and the shortest.
        let room: Vec<Option<(f64, f64)>> = towards
            .iter()
            .enumerate()
            .map(|(j, toward)| {
                let toward = (*toward)?;
                if toward == 0.0 {
                    return None;
                }
                let end = if toward > 0.0 {
                    bounds[j].1
                } else {
                    bounds[j].0
                };
                Some((((end - values[j] - moved[j]) / toward).max(0.0), end))
            })
            .collect();
        let reach = room
            .iter()
            .flatten()
            .map(|&(fraction, _)| fraction)
            .fold(1.0_f64, f64::min);
        for (j, toward) in towards.iter().enumerate() {
            if let Some(toward) = toward {
                solved[j] = true;
                moved[j] += reach * toward;
            }
        }
        if reach < 1.0 {
            for (j, room) in room.iter().enumerate() {
                if let Some((fraction, end)) = *room
                    && fraction <= reach
                {
                    moved[j] = end - values[j];
                    fixed_at[j] = Some(end);
                }
            }
            continue;
        }
        // A whole step: the maximum with these parameters fixed. Free the one whose slope points
        // furthest back into its interval, if any does.
        let curved = information_times(sums, layout, &moved);
        let mut freed: Option<(usize, f64)> = None;
        for (j, end) in fixed_at.iter().enumerate() {
            let Some(end) = end else { continue };
            let slope_there = slope[j] - curved[j];
            let inwards = if *end == bounds[j].0 {
                slope_there
            } else {
                -slope_there
            };
            let rounding = SLOPE_POINTING_IN * (slope[j].abs() + curved[j].abs());
            if inwards > rounding && freed.is_none_or(|(_, most)| inwards > most) {
                freed = Some((j, inwards));
            }
        }
        match freed {
            Some((j, _)) => fixed_at[j] = None,
            None => {
                finished = true;
                break;
            }
        }
    }
    BoundedStep {
        moves: moved
            .iter()
            .zip(&solved)
            .map(|(moved, solved)| solved.then_some(*moved))
            .collect(),
        finished,
    }
}

/// **Where each parameter sits in one flat list**: the cohort's eight, then each sample's own in turn
/// — the whole matrix's layout ([`FullInformation`]), used for the blocks too.
struct Layout {
    /// Where each sample's own parameters start.
    starts: Vec<usize>,
    /// How many parameters in all.
    side: usize,
}

impl Layout {
    fn of(sums: &InformationSums) -> Self {
        let mut starts = Vec::with_capacity(sums.sample_blocks.len());
        let mut side = COHORT_PARAMETERS;
        for s in 0..sums.sample_blocks.len() {
            starts.push(side);
            side += sums.own_parameters_of(s);
        }
        Self { starts, side }
    }

    /// Sample `s`'s own parameters' positions.
    fn of_sample(&self, s: usize) -> std::ops::Range<usize> {
        let end = self.starts.get(s + 1).copied().unwrap_or(self.side);
        self.starts[s]..end
    }

    /// The summed scores, flat.
    fn slope_of(&self, sums: &InformationSums) -> Vec<f64> {
        let mut slope = sums.slope.cohort.to_vec();
        for row in &sums.slope.samples {
            slope.extend_from_slice(row);
        }
        debug_assert_eq!(slope.len(), self.side, "the slope fills the layout");
        slope
    }

    /// The fit's parameters, flat: the duplicated class's three at zero when it is not fitted.
    fn values_of(&self, parameters: &Parameters, group_index: &[Vec<usize>]) -> Vec<f64> {
        let mut values = vec![0.0; self.side];
        values[cohort::NOISY_SHARE] = parameters.noisy_share;
        values[cohort::P_INVARIANT] = parameters.density.p_invariant;
        values[cohort::P_FIXED_ALT] = parameters.density.p_fixed_alt;
        values[cohort::DENSITY_A] = parameters.density.a;
        values[cohort::DENSITY_B] = parameters.density.b;
        if let Some(duplicated) = parameters.duplicated {
            values[cohort::DUPLICATED_SHARE] = duplicated.share;
            values[cohort::CARRIER_A] = duplicated.carrier_a;
            values[cohort::CARRIER_B] = duplicated.carrier_b;
        }
        for (s, own) in group_index.iter().enumerate() {
            let start = self.starts[s];
            values[start + sample::HOMOZYGOTE_EXCESS] = parameters.hom_excess[s];
            for (section, &group) in own.iter().enumerate() {
                values[start + sample::rate(section, 0)] = parameters.clean[group];
                values[start + sample::rate(section, 1)] = parameters.noisy[group];
            }
        }
        values
    }

    /// The interval the fit keeps each parameter in, flat.
    fn bounds(&self) -> Vec<(f64, f64)> {
        let mut bounds = COHORT_BOUNDS.to_vec();
        for s in 0..self.starts.len() {
            bounds.extend(
                self.of_sample(s)
                    .map(|at| bounds_of_own(at - self.starts[s])),
            );
        }
        bounds
    }

    /// A flat list cut back into the cohort's eight and each sample's own.
    fn split<T: Copy>(&self, flat: &[T]) -> ([T; COHORT_PARAMETERS], Vec<Vec<T>>) {
        let cohort = std::array::from_fn(|slot| flat[slot]);
        let samples = (0..self.starts.len())
            .map(|s| flat[self.of_sample(s)].to_vec())
            .collect();
        (cohort, samples)
    }
}

/// The information times a flat vector: on the whole matrix every pairing, on the arrow the cohort's
/// block, each sample's own, and each sample's with the cohort's both ways.
fn information_times(sums: &InformationSums, layout: &Layout, vector: &[f64]) -> Vec<f64> {
    if let Some(full) = &sums.full {
        return (0..layout.side)
            .map(|row| {
                (0..layout.side)
                    .filter(|&column| vector[column] != 0.0)
                    .map(|column| full.entry(row, column) * vector[column])
                    .sum()
            })
            .collect();
    }
    let mut product = vec![0.0; layout.side];
    for (row, into) in product[..COHORT_PARAMETERS].iter_mut().enumerate() {
        let entries = &sums.cohort[row * COHORT_PARAMETERS..(row + 1) * COHORT_PARAMETERS];
        for (entry, value) in entries.iter().zip(&vector[..COHORT_PARAMETERS]) {
            *into += entry * value;
        }
    }
    for s in 0..layout.starts.len() {
        let rows = layout.of_sample(s);
        let (start, n) = (rows.start, rows.len());
        for j in 0..n {
            let mut own = 0.0;
            for k in 0..n {
                own += sums.sample_blocks[s][j * n + k] * vector[start + k];
            }
            for c in 0..COHORT_PARAMETERS {
                let cross = sums.sample_cohort_blocks[s][j * COHORT_PARAMETERS + c];
                own += cross * vector[c];
                product[c] += cross * vector[start + j];
            }
            product[start + j] += own;
        }
    }
    product
}

/// `I⁻¹ g` on the whole matrix, over the parameters [`StandardErrors::of_the_whole_matrix`] inverts
/// less those `held`, from the flat `slope`; flat, `None` for the rest.
fn whole_matrix_solve(full: &FullInformation, slope: &[f64], held: &[bool]) -> Vec<Option<f64>> {
    let excess_is_fitted = fits_homozygote_excess(samples_informing_their_excess_in(full));
    let mut solved: Vec<usize> = Vec::with_capacity(full.side());
    for s in 0..full.samples() {
        for (slot, row) in full.rows_of_sample(s).enumerate() {
            if has_information(full.entry(row, row))
                && (slot != sample::HOMOZYGOTE_EXCESS || excess_is_fitted)
                && !held[row]
            {
                solved.push(row);
            }
        }
    }
    for (slot, &is_held) in held.iter().enumerate().take(COHORT_PARAMETERS) {
        if has_information(full.entry(slot, slot)) && !is_held {
            solved.push(slot);
        }
    }
    let matrix = square_of(&solved, |row, column| full.entry(row, column));
    let reference: Vec<f64> = solved.iter().map(|&row| full.entry(row, row)).collect();
    let identified = invert_identified(&matrix, solved.len(), &reference);
    let kept = identified.kept.len();
    let mut step = vec![None; full.side()];
    for (p, &position) in identified.kept.iter().enumerate() {
        let moved: f64 = identified
            .kept
            .iter()
            .enumerate()
            .map(|(q, &other)| identified.inverse[p * kept + q] * slope[solved[other]])
            .sum();
        step[solved[position]] = Some(moved);
    }
    step
}

/// One sample's part of the arrow's solve: its parameters kept, `A_s⁻¹ g_s` and `A_s⁻¹ B_s`.
type SampleSolve = (Vec<usize>, Vec<f64>, Vec<f64>);

/// `I⁻¹ g` on the arrow of blocks, by the same elimination the errors use
/// ([`StandardErrors::of_the_blocks`]), over the parameters it inverts less those `held`: each
/// sample's own parameters solved given the cohort's, the cohort's from what the samples leave of its
/// slope, and each sample's own step then corrected by the cohort's. Flat, `None` for the rest.
fn arrow_solve(
    sums: &InformationSums,
    layout: &Layout,
    slope: &[f64],
    held: &[bool],
) -> Vec<Option<f64>> {
    let samples = sums.sample_blocks.len();
    let excess_is_fitted = fits_homozygote_excess(samples_informing_their_excess(sums));
    let mut step = vec![None; layout.side];
    let cohort_informed: Vec<usize> = (0..COHORT_PARAMETERS)
        .filter(|&i| has_information(sums.cohort[i * COHORT_PARAMETERS + i]) && !held[i])
        .collect();
    let c = cohort_informed.len();
    let mut explained = vec![0.0; c * c];
    // The cohort's slope less what each sample's own parameters take up of it: `g_c − Σ B_sᵀ A_s⁻¹ g_s`.
    let mut cohort_slope_left: Vec<f64> = cohort_informed.iter().map(|&i| slope[i]).collect();
    // Per sample: its kept parameters, `A_s⁻¹ g_s` and `A_s⁻¹ B_s`.
    let mut per_sample: Vec<Option<SampleSolve>> = Vec::with_capacity(samples);
    for s in 0..samples {
        let start = layout.starts[s];
        let n = sums.own_parameters_of(s);
        let own_diagonal = |j: usize| sums.sample_blocks[s][j * n + j];
        let informed: Vec<usize> = (0..n)
            .filter(|&j| {
                has_information(own_diagonal(j))
                    && (j != sample::HOMOZYGOTE_EXCESS || excess_is_fitted)
                    && !held[start + j]
            })
            .collect();
        let own = square_of(&informed, |row, column| {
            sums.sample_blocks[s][row * n + column]
        });
        let reference: Vec<f64> = informed.iter().map(|&j| own_diagonal(j)).collect();
        let identified = invert_identified(&own, informed.len(), &reference);
        let kept: Vec<usize> = identified.kept.iter().map(|&p| informed[p]).collect();
        if kept.is_empty() {
            per_sample.push(None);
            continue;
        }
        let k = kept.len();
        let with_cohort: Vec<f64> = kept
            .iter()
            .flat_map(|&row| {
                cohort_informed.iter().map(move |&column| {
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column]
                })
            })
            .collect();
        let own_inverse_times_cross = multiply(&identified.inverse, &with_cohort, k, k, c);
        let own_slope: Vec<f64> = kept.iter().map(|&j| slope[start + j]).collect();
        let own_inverse_times_slope = multiply(&identified.inverse, &own_slope, k, k, 1);
        for p in 0..c {
            for q in 0..c {
                let mut entry = 0.0;
                for j in 0..k {
                    entry += with_cohort[j * c + p] * own_inverse_times_cross[j * c + q];
                }
                explained[p * c + q] += entry;
            }
            let mut taken = 0.0;
            for j in 0..k {
                taken += with_cohort[j * c + p] * own_inverse_times_slope[j];
            }
            cohort_slope_left[p] -= taken;
        }
        per_sample.push(Some((
            kept,
            own_inverse_times_slope,
            own_inverse_times_cross,
        )));
    }
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
    let v = identified.kept.len();
    // The cohort's step over the informed slots, zero for one dropped as not identified.
    let mut cohort_solved = vec![0.0; c];
    for (p, &position) in identified.kept.iter().enumerate() {
        let moved: f64 = identified
            .kept
            .iter()
            .enumerate()
            .map(|(q, &other)| identified.inverse[p * v + q] * cohort_slope_left[other])
            .sum();
        cohort_solved[position] = moved;
        step[cohort_informed[position]] = Some(moved);
    }
    for (s, solved) in per_sample.into_iter().enumerate() {
        let Some((kept, own_inverse_times_slope, own_inverse_times_cross)) = solved else {
            continue;
        };
        for (j, &which) in kept.iter().enumerate() {
            let carried: f64 = (0..c)
                .map(|p| own_inverse_times_cross[j * c + p] * cohort_solved[p])
                .sum();
            step[layout.starts[s] + which] = Some(own_inverse_times_slope[j] - carried);
        }
    }
    step
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

/// **How many samples carry information on their own homozygote excess**, from the blocks — the
/// count the excess is fitted across ([`fits_homozygote_excess`]). A sample with no read at an
/// ordinary position carries none, and does not make a one-sample cohort a two-sample one.
pub(super) fn samples_informing_their_excess(sums: &InformationSums) -> usize {
    (0..sums.sample_blocks.len())
        .filter(|&s| {
            let n = sums.own_parameters_of(s);
            has_information(
                sums.sample_blocks[s][sample::HOMOZYGOTE_EXCESS * n + sample::HOMOZYGOTE_EXCESS],
            )
        })
        .count()
}

/// The same, from the whole matrix.
fn samples_informing_their_excess_in(full: &FullInformation) -> usize {
    (0..full.samples())
        .filter(|&s| {
            let row = full.rows_of_sample(s).start + sample::HOMOZYGOTE_EXCESS;
            has_information(full.entry(row, row))
        })
        .count()
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
pub(in crate::parameter_estimation::joint) struct Identified {
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
pub(in crate::parameter_estimation::joint) fn invert_identified(
    matrix: &[f64],
    n: usize,
    reference: &[f64],
) -> Identified {
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
    use crate::parameter_estimation::joint::fit::ReportedNumber;
    use crate::parameter_estimation::joint::fit::information::{
        ONE_LIBRARY_SAMPLE_PARAMETERS, PositionScores, cohort, own_parameters,
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
        (blocks_of_dense(libraries, &full), full)
    }

    /// The blocks of a dense matrix in the layout — the cohort's eight, then each sample's own in
    /// turn — for samples of `libraries` libraries each; what pairs two samples is left out.
    fn blocks_of_dense(libraries: &[usize], dense: &[f64]) -> InformationSums {
        let sizes: Vec<usize> = libraries.iter().map(|&k| own_parameters(k)).collect();
        let n = COHORT_PARAMETERS + sizes.iter().sum::<usize>();
        let mut sums = InformationSums::new(&group_index_of(libraries));
        for row in 0..COHORT_PARAMETERS {
            for column in 0..COHORT_PARAMETERS {
                sums.cohort[row * COHORT_PARAMETERS + column] = dense[row * n + column];
            }
        }
        let mut first = COHORT_PARAMETERS;
        for (s, &size) in sizes.iter().enumerate() {
            for row in 0..size {
                for column in 0..size {
                    sums.sample_blocks[s][row * size + column] =
                        dense[(first + row) * n + first + column];
                }
                for column in 0..COHORT_PARAMETERS {
                    sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column] =
                        dense[(first + row) * n + column];
                }
            }
            first += size;
        }
        sums
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

    /// A random dense information matrix in the layout — every parameter paired with every other,
    /// two samples' included — for samples of `libraries` libraries each, and the sums holding it as
    /// the whole matrix. Scaled as [`an_arrow_of`] is, so every error is about 10⁻³.
    fn a_whole_matrix_of(libraries: &[usize], seed: u64) -> (InformationSums, Vec<f64>) {
        let n = COHORT_PARAMETERS + libraries.iter().map(|&k| own_parameters(k)).sum::<usize>();
        let dense: Vec<f64> = positive_definite(n, seed)
            .into_iter()
            .map(|entry| entry * 1e6)
            .collect();
        (with_the_whole_matrix(libraries, &dense), dense)
    }

    /// Sums holding `dense` as the whole matrix, for samples of `libraries` libraries each.
    fn with_the_whole_matrix(libraries: &[usize], dense: &[f64]) -> InformationSums {
        let group_index = group_index_of(libraries);
        let mut sums = InformationSums::new(&group_index);
        sums.full = Some(FullInformation::of_dense(&group_index, dense));
        sums
    }

    /// **The whole matrix gives the errors of its own dense inverse** (plan step A8): a random
    /// matrix with every pair of parameters informed — two samples' included, which the blocks
    /// never see — over samples of one, two, three and one library, 26 parameters; every error
    /// agrees with the square root of the dense inverse's diagonal to 10⁻¹² relative. An entry read
    /// from the wrong triangle, a sample's slots read at another's offset, or the blocks used
    /// instead, fails.
    #[test]
    fn the_whole_matrix_gives_the_errors_of_its_dense_inverse() {
        let libraries = [1, 2, 3, 1];
        let (sums, dense) = a_whole_matrix_of(&libraries, 17);
        let n = COHORT_PARAMETERS + 18;
        let errors = StandardErrors::of(&sums);
        assert_eq!(
            errors.samples.iter().map(Vec::len).collect::<Vec<_>>(),
            [3, 5, 7, 3]
        );
        assert_errors_match(&errors, &dense_errors_without(&dense, n, |_| false));
    }

    /// **On an arrow the whole matrix and the blocks give the same errors**: where no two samples'
    /// scores are paired the two ways invert the same matrix, and every error agrees to 10⁻¹²
    /// relative — the whole matrix only adds what the arrow leaves out.
    #[test]
    fn on_an_arrow_the_whole_matrix_and_the_blocks_agree() {
        let libraries = [1, 2, 3, 1];
        let (mut sums, dense) = an_arrow_of(&libraries, 19);
        let n = COHORT_PARAMETERS + 18;
        let blocks = StandardErrors::of_the_blocks(&sums);
        sums.full = Some(FullInformation::of_dense(
            &group_index_of(&libraries),
            &dense,
        ));
        let whole = StandardErrors::of(&sums);
        for i in 0..n {
            match (error_at(&blocks, i), error_at(&whole, i)) {
                (StandardError::Estimated(from_blocks), StandardError::Estimated(from_whole)) => {
                    assert!(
                        (from_blocks - from_whole).abs() < 1e-12 * from_blocks,
                        "parameter {i}: {from_blocks} from the blocks, {from_whole} from the whole"
                    );
                }
                (from_blocks, from_whole) => {
                    panic!(
                        "parameter {i}: {from_blocks:?} from the blocks, {from_whole:?} from the whole"
                    )
                }
            }
        }
    }

    /// **On an arrow the whole matrix and the blocks drop the same parameters**, which rests on the
    /// order the whole matrix takes them in — each sample's own first, the samples in turn, then the
    /// cohort's — since the blocks judge the cohort's only after the samples' share is taken out.
    /// On an arrow over samples of one, two and one library, the density's first shape is made to
    /// copy sample 0's clean rate (twice its row), and sample 1's second library's mismapped rate
    /// to copy that library's clean rate (three times its row). Both ways drop the density's shape,
    /// keep sample 0's rate, drop the second library's mismapped rate, and agree on every other
    /// error to 10⁻⁹ relative. Taking the cohort's parameters first would keep the shape and drop
    /// the sample's rate instead.
    #[test]
    fn on_an_arrow_the_whole_matrix_and_the_blocks_drop_the_same_parameters() {
        let libraries = [1, 2, 1];
        let (_, mut dense) = an_arrow_of(&libraries, 23);
        let n = COHORT_PARAMETERS + 11;
        let sample_rate = COHORT_PARAMETERS + sample::CLEAN_ERROR_RATE;
        let second_library = COHORT_PARAMETERS + 3;
        let (second_clean, second_mismapped) = (
            second_library + sample::rate(1, 0),
            second_library + sample::rate(1, 1),
        );
        for (copy, of, factor) in [
            (cohort::DENSITY_A, sample_rate, 2.0),
            (second_mismapped, second_clean, 3.0),
        ] {
            for k in 0..n {
                dense[copy * n + k] = factor * dense[of * n + k];
            }
            for k in 0..n {
                dense[k * n + copy] = factor * dense[k * n + of];
            }
        }
        let mut sums = blocks_of_dense(&libraries, &dense);
        let blocks = StandardErrors::of_the_blocks(&sums);
        sums.full = Some(FullInformation::of_dense(
            &group_index_of(&libraries),
            &dense,
        ));
        let whole = StandardErrors::of(&sums);
        for errors in [&blocks, &whole] {
            assert_eq!(
                errors.cohort[cohort::DENSITY_A],
                StandardError::NotIdentified
            );
            assert!(matches!(
                error_at(errors, sample_rate),
                StandardError::Estimated(_)
            ));
            assert_eq!(
                error_at(errors, second_mismapped),
                StandardError::NotIdentified
            );
        }
        for i in 0..n {
            match (error_at(&blocks, i), error_at(&whole, i)) {
                (StandardError::Estimated(from_blocks), StandardError::Estimated(from_whole)) => {
                    assert!(
                        (from_blocks - from_whole).abs() < 1e-9 * from_blocks,
                        "parameter {i}: {from_blocks} from the blocks, {from_whole} from the whole"
                    );
                }
                (from_blocks, from_whole) => assert_eq!(from_blocks, from_whole, "parameter {i}"),
            }
        }
    }

    /// **The whole matrix judges each slot against its own parameter's interval**, as the blocks do:
    /// a second library's mismapped rate rescaled so its error comes out at 0.3 has that error —
    /// inside the mismapped rate's interval (0.0001 to 0.45), outside the clean rate's (0.000001 to
    /// 0.2), so judging every slot by the clean rate's would report it wider than its range.
    #[test]
    fn the_whole_matrix_judges_each_slot_against_its_own_interval() {
        let libraries = [1, 2];
        let (sums, dense) = a_whole_matrix_of(&libraries, 47);
        let n = COHORT_PARAMETERS + 8;
        let slot = sample::rate(1, 1);
        let row = COHORT_PARAMETERS + 3 + slot;
        let factor = 0.3
            / StandardErrors::of(&sums).samples[1][slot]
                .value()
                .expect("an error");
        let mut rescaled = dense;
        for j in 0..n {
            rescaled[row * n + j] /= factor;
            rescaled[j * n + row] /= factor;
        }
        let rescaled_error =
            StandardErrors::of(&with_the_whole_matrix(&libraries, &rescaled)).samples[1][slot];
        assert!(
            matches!(rescaled_error, StandardError::Estimated(error) if (error - 0.3).abs() < 1e-9),
            "{rescaled_error:?}"
        );
    }

    /// **The whole matrix says why a parameter has no error, by the blocks' rules**: a sample whose
    /// rows are empty has no information; the duplicated share made twice the mismapped share's
    /// row cannot be told apart from it and is dropped alone, the mismapped share keeping its
    /// error; at one sample the homozygote excess is held fixed. Every other error equals the dense
    /// inverse's without those parameters.
    #[test]
    fn the_whole_matrix_says_why_a_parameter_has_no_error() {
        let libraries = [1, 2, 1];
        let (_, mut dense) = a_whole_matrix_of(&libraries, 29);
        let n = COHORT_PARAMETERS + 11;
        let silent: Vec<usize> = (COHORT_PARAMETERS + 3..COHORT_PARAMETERS + 8).collect();
        for &i in &silent {
            for k in 0..n {
                dense[i * n + k] = 0.0;
                dense[k * n + i] = 0.0;
            }
        }
        let (copy, of) = (cohort::DUPLICATED_SHARE, cohort::NOISY_SHARE);
        for k in 0..n {
            dense[copy * n + k] = 2.0 * dense[of * n + k];
        }
        for k in 0..n {
            dense[k * n + copy] = 2.0 * dense[k * n + of];
        }
        let errors = StandardErrors::of(&with_the_whole_matrix(&libraries, &dense));
        assert_eq!(errors.samples[1], [StandardError::NoInformation; 5]);
        assert_eq!(errors.cohort[copy], StandardError::NotIdentified);
        assert!(matches!(errors.cohort[of], StandardError::Estimated(_)));
        assert_errors_match(
            &errors,
            &dense_errors_without(&dense, n, |i| silent.contains(&i) || i == copy),
        );

        let (one, dense_of_one) = a_whole_matrix_of(&[1], 31);
        let excess = COHORT_PARAMETERS + sample::HOMOZYGOTE_EXCESS;
        let errors = StandardErrors::of(&one);
        assert_eq!(
            errors.samples[0][sample::HOMOZYGOTE_EXCESS],
            StandardError::HeldFixed
        );
        assert_errors_match(
            &errors,
            &dense_errors_without(&dense_of_one, COHORT_PARAMETERS + 3, |i| i == excess),
        );
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
    /// diagonal say `NoInformation`, and every other error equals the arrow's without them. Three
    /// samples, so that with the second's excess uninformed two still inform theirs and the excess
    /// stays fitted ([`fits_homozygote_excess`]).
    #[test]
    fn a_diagonal_that_is_not_finite_is_no_information() {
        let samples = 3;
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

    /// **Two samples, one of which informs nothing about its excess, hold the other's excess**: the
    /// excess is fitted across the samples that inform theirs, and one is not enough. On both ways
    /// of inverting, the informed sample's excess is held fixed.
    #[test]
    fn one_sample_informing_its_excess_holds_it() {
        let excess = sample::HOMOZYGOTE_EXCESS;
        let (mut sums, _) = an_arrow(2, 79);
        sums.sample_blocks[1][excess * ONE_LIBRARY_SAMPLE_PARAMETERS + excess] = 0.0;
        assert_eq!(samples_informing_their_excess(&sums), 1);
        let errors = StandardErrors::of(&sums);
        assert_eq!(errors.samples[0][excess], StandardError::HeldFixed);
        assert_eq!(errors.samples[1][excess], StandardError::NoInformation);
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

    /// What the data determine about one pair, by hand: the share above 1% and its error, and the
    /// trade's cost, from which the verdict follows as `determined` gives it.
    fn determined(share: f64, error: Option<f64>, trade_cost: Option<f64>) -> DeterminedShare {
        DeterminedShare {
            share_above: ReportedNumber {
                value: share,
                standard_error: error,
            },
            trade_cost,
            separately_determined: trade_cost.is_some_and(|cost| cost >= 2.0),
        }
    }

    /// The two reports `described` reads for `parameters`, with these verdicts; the numbers'
    /// errors, which the log takes from the standard errors and not from the reports, left out.
    fn reports_for_the_log(
        parameters: &Parameters,
        density: DeterminedShare,
        duplicated: Option<DeterminedShare>,
    ) -> (DensityReport, Option<DuplicatedReport>) {
        let bare = |value: f64| ReportedNumber {
            value,
            standard_error: None,
        };
        let values = parameters.density;
        (
            DensityReport {
                invariant_share: bare(values.p_invariant),
                fixed_non_reference_share: bare(values.p_fixed_alt),
                segregating_shape_a: bare(values.a),
                segregating_shape_b: bare(values.b),
                expected_heterozygosity: bare(values.expected_heterozygosity()),
                expected_alternative_frequency: bare(values.expected_alternative_frequency()),
                determined: density,
            },
            duplicated.map(|determined| {
                let values = parameters.duplicated.expect("the class is fitted");
                DuplicatedReport {
                    duplicated_share: bare(values.share),
                    carrier_shape_a: bare(values.carrier_a),
                    carrier_shape_b: bare(values.carrier_b),
                    determined,
                }
            }),
        )
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
        let no_trade = determined(0.06, Some(2e-3), Some(f64::INFINITY));
        let without_the_class = parameters_for_the_log(false, 4, 4);
        let (density, _) = reports_for_the_log(&without_the_class, no_trade, None);
        assert_eq!(
            errors.described(&without_the_class, &one_group_each, &density, None),
            "mismapped share 3.1200e-2 ± 2.00e-4, invariant share 9.0000e-1 ± 3.25e-3, fixed \
             non-reference share 1.2500e-2, no standard error (not identified), allele-frequency \
             shape a 5.0000e-1, no standard error (it came out at 6.10e1, wider than the \
             parameter's whole range), allele-frequency shape b 6.1000e1 ± 5.00e-1, share \
             segregating above 1% 6.00e-2 ± 2.0e-3 (the invariant share and allele-frequency shape \
             a are separately determined: the segregating share cannot change without changing \
             it); error rates at \
             ordinary positions (4 read groups): median standard error 2.00e-4, largest 4.00e-4, \
             none missing; error rates at mismapped positions (4 read groups): median standard \
             error 1.00e-2, largest 2.00e-2, 2 missing (1 no information, 1 not identified); \
             homozygote excesses (4 samples): none has a standard error (4 held fixed)"
        );
        let parameters = parameters_for_the_log(true, 4, 4);
        let (density, duplicated) = reports_for_the_log(
            &parameters,
            no_trade,
            Some(determined(8e-3, Some(1e-4), Some(5.3))),
        );
        let with_the_class =
            errors.described(&parameters, &one_group_each, &density, duplicated.as_ref());
        assert!(
            with_the_class.contains(
                "allele-frequency shape b 6.1000e1 ± 5.00e-1, share segregating above 1% 6.00e-2 ± \
                 2.0e-3 (the invariant share and allele-frequency shape a are separately determined: \
                 the segregating share cannot change without changing it), duplicated share \
                 1.0000e-2 ± 1.00e-3, carrier-frequency shape a 1.5000e0, no standard error (no \
                 information), carrier-frequency shape b 5.0000e0, no standard error (not \
                 identified), share of positions duplicated in more than 1% of the panel 8.00e-3 ± \
                 1.0e-4 (the duplicated share and carrier-frequency shape a are separately \
                 determined: halving or doubling the duplicated share with it held costs 5.30 \
                 log-likelihood units); error rates at ordinary positions"
            ),
            "{with_the_class}"
        );
        let mut two_libraries = errors.clone();
        two_libraries.samples[0].extend([Estimated(5e-4), NoInformation]);
        let second_group_of_sample_0 = vec![vec![0, 2], vec![1], vec![3], vec![4]];
        let line = two_libraries.described(
            &parameters_for_the_log(false, 5, 4),
            &second_group_of_sample_0,
            &density,
            None,
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

    /// **A pair the data do not separate is printed with the verdict in place of its errors**
    /// (spec §2.1.3, *What is reported*): the first number with the share the data fix, its error,
    /// and what halving or doubling the traded share costs; the second naming the first. With the
    /// 59-sample numbers of the spec's example — the density's cost 1.35, the class's 0.116 — and,
    /// at one sample, a cost that cannot be computed and a share with no error. The number's own
    /// error, which the standard errors still hold, is not printed.
    #[test]
    fn a_pair_not_separately_determined_is_printed_with_the_verdict() {
        use StandardError::{Estimated, HeldFixed, NotIdentified};
        let errors = StandardErrors {
            cohort: [
                Estimated(1.4e-4),
                Estimated(0.068),
                Estimated(4.05e-6),
                Estimated(0.0089),
                Estimated(0.0605),
                Estimated(0.023),
                Estimated(0.031),
                Estimated(0.0466),
            ],
            samples: vec![vec![Estimated(4e-4), Estimated(0.01), HeldFixed]],
        };
        let mut parameters = parameters_for_the_log(true, 1, 1);
        parameters.density.p_invariant = 0.84159;
        parameters.density.a = 0.02;
        if let Some(duplicated) = parameters.duplicated.as_mut() {
            duplicated.share = 0.015427;
            duplicated.carrier_a = 0.02;
        }
        let (density, duplicated) = reports_for_the_log(
            &parameters,
            determined(8.447e-3, Some(2.04e-4), Some(1.3467)),
            Some(determined(1.2122e-3, Some(3.38e-5), Some(0.11571))),
        );
        let line = errors.described(&parameters, &[vec![0]], &density, duplicated.as_ref());
        assert!(
            line.starts_with(
                "mismapped share 3.1200e-2 ± 1.40e-4, invariant share 8.4159e-1 (not separately \
                 determined from allele-frequency shape a: the data fix the share segregating above \
                 1%, 8.45e-3 ± 2.0e-4, not how much lies below it; halving or doubling the \
                 segregating share costs 1.35 log-likelihood units), fixed non-reference share \
                 1.2500e-2 ± 4.05e-6, allele-frequency shape a 2.0000e-2 (not separately determined \
                 from the invariant share), allele-frequency shape b 6.1000e1 ± 6.05e-2, duplicated \
                 share 1.5427e-2 (not separately determined from carrier-frequency shape a: the \
                 data fix the share of positions duplicated in more than 1% of the panel, 1.21e-3 ± \
                 3.4e-5, not how many are duplicated in fewer; halving or doubling the duplicated \
                 share costs 0.116 log-likelihood units), carrier-frequency shape a 2.0000e-2 (not \
                 separately determined from the duplicated share), carrier-frequency shape b \
                 5.0000e0 ± 4.66e-2; error rates"
            ),
            "{line}"
        );

        let mut at_one_sample = errors.clone();
        at_one_sample.cohort[cohort::DENSITY_A] = NotIdentified;
        at_one_sample.cohort[cohort::DENSITY_B] = NotIdentified;
        let (density, _) = reports_for_the_log(&parameters, determined(2.33e-3, None, None), None);
        let without_the_class = Parameters {
            duplicated: None,
            ..parameters
        };
        let line = at_one_sample.described(&without_the_class, &[vec![0]], &density, None);
        assert!(
            line.contains(
                "invariant share 8.4159e-1 (not separately determined from allele-frequency shape \
                 a: the data fix the share segregating above 1%, 2.33e-3 with no standard error, \
                 not how much lies below it; what halving or doubling the segregating share costs \
                 cannot be computed, as the data do not determine every number it is computed \
                 from), fixed non-reference share 1.2500e-2 ± 4.05e-6, allele-frequency shape a \
                 2.0000e-2 (not separately determined from the invariant share), allele-frequency \
                 shape b 6.1000e1, no standard error (not identified); error rates"
            ),
            "{line}"
        );
    }

    /// **The duplicated class's part is absent when the class is off**: no duplicated number and no
    /// verdict on its pair, even when a report for it is handed over; and a separately determined
    /// density's share comes after its second shape, its cost in three figures.
    #[test]
    fn the_duplicated_part_is_absent_when_the_class_is_off() {
        use StandardError::{Estimated, HeldFixed};
        let errors = StandardErrors {
            cohort: [Estimated(1e-3); COHORT_PARAMETERS],
            samples: vec![vec![Estimated(4e-4), Estimated(0.01), HeldFixed]],
        };
        let with_the_class = parameters_for_the_log(true, 1, 1);
        let (density, duplicated) = reports_for_the_log(
            &with_the_class,
            determined(0.0806, Some(3.2e-3), Some(324.25)),
            Some(determined(1e-3, Some(1e-4), Some(0.5))),
        );
        let without_the_class = parameters_for_the_log(false, 1, 1);
        for handed in [None, duplicated.as_ref()] {
            let line = errors.described(&without_the_class, &[vec![0]], &density, handed);
            assert!(
                !line.contains("duplicated") && !line.contains("carrier"),
                "{line}"
            );
            assert!(
                line.contains(
                    "allele-frequency shape b 6.1000e1 ± 1.00e-3, share segregating above 1% \
                     8.06e-2 ± 3.2e-3 (the invariant share and allele-frequency shape a are \
                     separately determined: halving or doubling the segregating share with it held \
                     costs 324 log-likelihood units); error rates"
                ),
                "{line}"
            );
        }
    }

    /// **A cost in three significant figures**, in scientific notation from a thousand up, and
    /// above a million as just that.
    #[test]
    fn a_cost_is_printed_in_three_figures() {
        assert_eq!(in_three_figures(1.346728), "1.35");
        assert_eq!(in_three_figures(0.005494), "0.00549");
        assert_eq!(in_three_figures(0.11571), "0.116");
        assert_eq!(in_three_figures(324.25), "324");
        assert_eq!(in_three_figures(4302.7), "4.30e3");
        assert_eq!(a_cost(1.346728), "1.35 log-likelihood units");
        assert_eq!(a_cost(4.66e16), "over a million log-likelihood units");
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

    /// A fit's parameters for samples of `libraries` libraries each, every one well inside its
    /// interval, with the duplicated class fitted.
    fn parameters_inside(libraries: &[usize]) -> Parameters {
        let mut parameters = parameters_for_the_log(true, libraries.iter().sum(), libraries.len());
        parameters.density.b = 5.0;
        parameters
    }

    /// A slope in the errors' layout for samples of `libraries` libraries: every entry different,
    /// of the size a few thousand positions give, and the same numbers flattened in the layout's order.
    fn a_slope_of(libraries: &[usize], seed: u64) -> (PositionScores, Vec<f64>) {
        let group_index = group_index_of(libraries);
        let mut slope = PositionScores::new(&group_index, 0);
        let mut state = seed;
        let mut draw = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 20.0
        };
        for entry in slope.cohort.iter_mut() {
            *entry = draw();
        }
        for row in &mut slope.samples {
            for entry in row.iter_mut() {
                *entry = draw();
            }
        }
        let flat = slope
            .cohort
            .iter()
            .chain(slope.samples.iter().flatten())
            .copied()
            .collect();
        (slope, flat)
    }

    /// `matrix⁻¹ slope` over the parameters `kept` says to solve for, by the dense inverse, in the
    /// errors' layout; `None` for the rest.
    fn dense_step(matrix: &[f64], slope: &[f64], kept: impl Fn(usize) -> bool) -> Vec<Option<f64>> {
        let n = slope.len();
        let kept: Vec<usize> = (0..n).filter(|&i| kept(i)).collect();
        let inverse = inverse_of_positive_definite(
            &square_of(&kept, |row, column| matrix[row * n + column]),
            kept.len(),
        )
        .expect("positive definite");
        let mut step = vec![None; n];
        for (at, &i) in kept.iter().enumerate() {
            step[i] = Some(
                kept.iter()
                    .enumerate()
                    .map(|(other, &j)| inverse[at * kept.len() + other] * slope[j])
                    .sum(),
            );
        }
        step
    }

    /// A step in the errors' layout reordered as the fit's vector lists its parameters.
    fn in_the_fits_order(
        step: &[Option<f64>],
        libraries: &[usize],
        parameters: &Parameters,
    ) -> Vec<Option<f64>> {
        let cohort: [Option<f64>; COHORT_PARAMETERS] = std::array::from_fn(|slot| step[slot]);
        let mut at = COHORT_PARAMETERS;
        let samples: Vec<Vec<Option<f64>>> = libraries
            .iter()
            .map(|&k| {
                let row = step[at..at + own_parameters(k)].to_vec();
                at += own_parameters(k);
                row
            })
            .collect();
        in_coordinate_order(
            &cohort,
            &samples,
            parameters,
            &group_index_of(libraries),
            None,
        )
    }

    fn assert_steps_match(got: &[Option<f64>], expected: &[Option<f64>], what: &str) {
        assert_eq!(got.len(), expected.len(), "{what}");
        for (i, (got, expected)) in got.iter().zip(expected).enumerate() {
            match (got, expected) {
                (Some(got), Some(expected)) => assert!(
                    (got - expected).abs() <= 1e-9 * expected.abs().max(1e-12),
                    "{what}, parameter {i}: {got} against {expected}"
                ),
                (None, None) => {}
                _ => panic!("{what}, parameter {i}: {got:?} against {expected:?}"),
            }
        }
    }

    /// **On the whole matrix the Newton step is the dense solve**: samples of one, two and three
    /// libraries, every parameter informed, every step inside its interval — each parameter's
    /// distance equals the dense inverse times the slope to 10⁻⁹, reordered as the fit lists them.
    #[test]
    fn on_the_whole_matrix_the_newton_step_is_the_dense_solve() {
        let libraries = [1, 2, 3, 1];
        let (mut sums, dense) = a_whole_matrix_of(&libraries, 91);
        let (slope, flat) = a_slope_of(&libraries, 92);
        sums.slope = slope;
        let parameters = parameters_inside(&libraries);
        let got = newton_step(
            &sums,
            &StandardErrors::of(&sums),
            &parameters,
            &group_index_of(&libraries),
        )
        .distances;
        let expected = in_the_fits_order(
            &dense_step(&dense, &flat, |_| true),
            &libraries,
            &parameters,
        );
        assert!(expected.iter().all(Option::is_some));
        assert_steps_match(&got, &expected, "whole matrix");
    }

    /// **On the arrow the Newton step is the dense solve of the whole arrow**: the blocks solved one
    /// sample at a time, then the cohort's, then each sample's corrected — the same numbers as
    /// inverting the whole arrow at once.
    #[test]
    fn on_the_arrow_the_newton_step_is_the_dense_solve() {
        let libraries = [1, 2, 3, 1];
        let (mut sums, dense) = an_arrow_of(&libraries, 93);
        assert!(sums.full.is_none(), "the blocks only");
        let (slope, flat) = a_slope_of(&libraries, 94);
        sums.slope = slope;
        let parameters = parameters_inside(&libraries);
        let got = newton_step(
            &sums,
            &StandardErrors::of(&sums),
            &parameters,
            &group_index_of(&libraries),
        )
        .distances;
        let expected = in_the_fits_order(
            &dense_step(&dense, &flat, |_| true),
            &libraries,
            &parameters,
        );
        assert_steps_match(&got, &expected, "arrow");
    }

    /// Set the parameter at flat position `row` of the errors' layout — a cohort slot, or the first
    /// sample's homozygote excess — to `value`.
    fn set_at(parameters: &mut Parameters, row: usize, value: f64) {
        match row {
            cohort::P_FIXED_ALT => parameters.density.p_fixed_alt = value,
            cohort::DENSITY_B => parameters.density.b = value,
            _ if row == COHORT_PARAMETERS + sample::HOMOZYGOTE_EXCESS => {
                parameters.hom_excess[0] = value
            }
            _ => unreachable!("the tests move only these"),
        }
    }

    /// **A parameter whose free step would carry it past an end of its interval stops at that end**,
    /// its distance the way there, and every other parameter's distance is the dense solve without it,
    /// its move taken out of their scores — for a cohort share, a density shape and a sample's
    /// homozygote excess, on both matrices, at whichever end the free step points to, with the slope
    /// and its negative so each parameter is tried at both ends; the parameter placed a quarter of its
    /// free step inside the end, and again exactly on it.
    #[test]
    fn a_parameter_carried_past_an_end_stops_there() {
        let libraries = [1, 2, 1];
        let group_index = group_index_of(&libraries);
        let excess = COHORT_PARAMETERS + sample::HOMOZYGOTE_EXCESS;
        let mut ends_reached = [false; 2];
        for row in [cohort::P_FIXED_ALT, cohort::DENSITY_B, excess] {
            let bounds = if row == excess {
                HOM_EXCESS_BOUNDS
            } else {
                COHORT_BOUNDS[row]
            };
            for sign in [1.0, -1.0] {
                for (label, (mut sums, dense)) in [
                    ("whole matrix", a_whole_matrix_of(&libraries, 95)),
                    ("arrow", an_arrow_of(&libraries, 97)),
                ] {
                    let (mut slope, mut flat) = a_slope_of(&libraries, 96);
                    for entry in slope
                        .cohort
                        .iter_mut()
                        .chain(slope.samples.iter_mut().flatten())
                    {
                        *entry *= sign;
                    }
                    flat.iter_mut().for_each(|entry| *entry *= sign);
                    sums.slope = slope;
                    let n = flat.len();
                    let free = dense_step(&dense, &flat, |_| true)[row].expect("solved");
                    let end = if free < 0.0 { bounds.0 } else { bounds.1 };
                    ends_reached[usize::from(free > 0.0)] = true;
                    for inside in [0.25 * free.abs(), 0.0] {
                        let value = if free < 0.0 {
                            end + inside
                        } else {
                            end - inside
                        };
                        let mut parameters = parameters_inside(&libraries);
                        set_at(&mut parameters, row, value);
                        let moved = end - value;
                        let slope_left: Vec<f64> = (0..n)
                            .map(|i| flat[i] - dense[i * n + row] * moved)
                            .collect();
                        let mut expected = dense_step(&dense, &slope_left, |i| i != row);
                        expected[row] = Some(moved);
                        let expected = in_the_fits_order(&expected, &libraries, &parameters);
                        let got = newton_step(
                            &sums,
                            &StandardErrors::of(&sums),
                            &parameters,
                            &group_index,
                        )
                        .distances;
                        assert_steps_match(
                            &got,
                            &expected,
                            &format!("{label}, row {row}, sign {sign}, {inside} inside"),
                        );
                    }
                }
            }
        }
        assert_eq!(ends_reached, [true, true], "both ends reached");
    }

    /// **A parameter on an end whose step points back into its interval is left free**: its distance
    /// and every other one is the unbounded dense solve.
    #[test]
    fn a_parameter_on_an_end_stepping_inwards_is_free() {
        let libraries = [1, 2, 1];
        let group_index = group_index_of(&libraries);
        let row = cohort::P_FIXED_ALT;
        for (label, (mut sums, dense)) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 95)),
            ("arrow", an_arrow_of(&libraries, 97)),
        ] {
            let (slope, flat) = a_slope_of(&libraries, 96);
            sums.slope = slope;
            let free = dense_step(&dense, &flat, |_| true);
            let step = free[row].expect("solved");
            // The end the step points away from.
            let end = if step > 0.0 {
                COHORT_BOUNDS[row].0
            } else {
                COHORT_BOUNDS[row].1
            };
            let mut parameters = parameters_inside(&libraries);
            set_at(&mut parameters, row, end);
            let got =
                newton_step(&sums, &StandardErrors::of(&sums), &parameters, &group_index).distances;
            assert_steps_match(
                &got,
                &in_the_fits_order(&free, &libraries, &parameters),
                label,
            );
        }
    }

    /// **A clean error rate the golden section leaves beside its floor or its ceiling is fixed at
    /// that end**, with a distance of the gap the search left: the search stops within 10⁻⁹ of an end
    /// rather than on it (it left one rate 4.4 × 10⁻¹⁰ above its floor of 10⁻⁶ in review), and a rate
    /// held only when exactly at its end never settled.
    #[test]
    fn a_rate_the_golden_section_leaves_beside_an_end_stops_there() {
        use crate::parameter_estimation::joint::fit::{CLEAN_ERROR_BOUNDS, golden_section};
        let (low, high) = CLEAN_ERROR_BOUNDS;
        let beside = [
            golden_section(&|rate: f64| -rate, low, high),
            golden_section(&|rate: f64| rate, low, high),
        ];
        assert!(
            beside[0] > low && beside[1] < high,
            "the search stops beside the ends: {beside:?}"
        );
        let libraries = [1, 2, 1];
        for sign in [1.0, -1.0] {
            let (mut sums, dense) = a_whole_matrix_of(&libraries, 95);
            let (mut slope, mut flat) = a_slope_of(&libraries, 96);
            for entry in slope
                .cohort
                .iter_mut()
                .chain(slope.samples.iter_mut().flatten())
            {
                *entry *= sign;
            }
            flat.iter_mut().for_each(|entry| *entry *= sign);
            sums.slope = slope;
            let row = COHORT_PARAMETERS + sample::rate(0, 0);
            let free = dense_step(&dense, &flat, |_| true)[row].expect("solved");
            let mut parameters = parameters_inside(&libraries);
            parameters.clean[0] = if free < 0.0 { beside[0] } else { beside[1] };
            let got = newton_step(
                &sums,
                &StandardErrors::of(&sums),
                &parameters,
                &group_index_of(&libraries),
            )
            .distances;
            // In the fit's order: the cohort's eight, then read group 0's clean rate.
            let distance = got[COHORT_PARAMETERS].expect("solved");
            assert!(
                distance.abs() < 1e-9,
                "sign {sign}: a rate at {} moves {distance}",
                parameters.clean[0]
            );
        }
    }

    /// **At one sample the homozygote excess is not solved for** — the fit holds it at its start —
    /// and every other distance is the dense solve without it, on both matrices; the excess set
    /// inside its interval, so no end can mask a mistake.
    #[test]
    fn at_one_sample_the_newton_step_leaves_the_held_excess_out() {
        let libraries = [2];
        let excess = COHORT_PARAMETERS + sample::HOMOZYGOTE_EXCESS;
        for (label, (mut sums, dense)) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 101)),
            ("arrow", an_arrow_of(&libraries, 103)),
        ] {
            let (slope, flat) = a_slope_of(&libraries, 102);
            sums.slope = slope;
            let mut parameters = parameters_inside(&libraries);
            parameters.hom_excess[0] = 0.5;
            let got = newton_step(
                &sums,
                &StandardErrors::of(&sums),
                &parameters,
                &group_index_of(&libraries),
            )
            .distances;
            let expected = in_the_fits_order(
                &dense_step(&dense, &flat, |i| i != excess),
                &libraries,
                &parameters,
            );
            assert_eq!(expected.last(), Some(&None), "the excess is not solved for");
            assert_steps_match(&got, &expected, label);
        }
    }

    /// **A parameter the Newton step cannot solve for has no distance, and the rest are solved
    /// without it**, on both matrices: the second sample's first mismapped rate given no information
    /// (its row, column and score zero, as for a library with no reads), and the fixed non-reference
    /// share made a copy of the invariant share in every entry and in its score (the data cannot tell
    /// them apart). Both come back `None`, as their errors do, and every other distance is the dense
    /// solve with the two left out — so neither leaks a stale step into the others'.
    #[test]
    fn a_parameter_the_newton_step_cannot_solve_has_no_distance() {
        let libraries = [1, 2, 1];
        let n = COHORT_PARAMETERS + libraries.iter().map(|&k| own_parameters(k)).sum::<usize>();
        let uninformed = COHORT_PARAMETERS + own_parameters(1) + sample::rate(0, 1);
        let (copied, copy_of) = (cohort::P_FIXED_ALT, cohort::P_INVARIANT);
        for (label, (_, mut dense), whole) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 111), true),
            ("arrow", an_arrow_of(&libraries, 113), false),
        ] {
            let (mut slope, mut flat) = a_slope_of(&libraries, 112);
            for j in 0..n {
                dense[uninformed * n + j] = 0.0;
                dense[j * n + uninformed] = 0.0;
            }
            flat[uninformed] = 0.0;
            for j in 0..n {
                dense[copied * n + j] = dense[copy_of * n + j];
            }
            for j in 0..n {
                dense[j * n + copied] = dense[j * n + copy_of];
            }
            flat[copied] = flat[copy_of];
            slope.cohort[copied] = slope.cohort[copy_of];
            slope.samples[1][sample::rate(0, 1)] = 0.0;
            let mut sums = if whole {
                with_the_whole_matrix(&libraries, &dense)
            } else {
                blocks_of_dense(&libraries, &dense)
            };
            sums.slope = slope;
            let errors = StandardErrors::of(&sums);
            assert_eq!(
                error_at(&errors, uninformed),
                StandardError::NoInformation,
                "{label}"
            );
            assert_eq!(
                error_at(&errors, copied),
                StandardError::NotIdentified,
                "{label}"
            );
            let parameters = parameters_inside(&libraries);
            let got = newton_step(
                &sums,
                &StandardErrors::of(&sums),
                &parameters,
                &group_index_of(&libraries),
            )
            .distances;
            let in_layout = dense_step(&dense, &flat, |i| i != uninformed && i != copied);
            let expected = in_the_fits_order(&in_layout, &libraries, &parameters);
            assert_eq!(
                expected.iter().filter(|step| step.is_none()).count(),
                2,
                "{label}: the two left out"
            );
            assert_steps_match(&got, &expected, label);
        }
    }

    /// What a step in the errors' layout is worth under the quadratic model of `matrix` and
    /// `slope`, `g·d − ½ dᵀ M d`, densely; a parameter without a step moves 0.
    fn dense_decrement(matrix: &[f64], slope: &[f64], step: &[Option<f64>]) -> f64 {
        let n = slope.len();
        let moved: Vec<f64> = step.iter().map(|moved| moved.unwrap_or(0.0)).collect();
        let gained: f64 = (0..n).map(|i| slope[i] * moved[i]).sum();
        let lost: f64 = (0..n)
            .map(|i| moved[i] * (0..n).map(|j| matrix[i * n + j] * moved[j]).sum::<f64>())
            .sum();
        gained - 0.5 * lost
    }

    /// **With no bound binding, the decrement is `½ gᵀ I⁻¹ g`**, on the whole matrix and on the
    /// arrow: the inputs of the two dense-solve tests above, whose steps no bound stops.
    #[test]
    fn with_no_bound_binding_the_decrement_is_half_the_slope_through_the_inverse() {
        let libraries = [1, 2, 3, 1];
        for (label, (mut sums, dense), slope_seed) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 91), 92),
            ("arrow", an_arrow_of(&libraries, 93), 94),
        ] {
            let (slope, flat) = a_slope_of(&libraries, slope_seed);
            sums.slope = slope;
            let parameters = parameters_inside(&libraries);
            let got = newton_step(
                &sums,
                &StandardErrors::of(&sums),
                &parameters,
                &group_index_of(&libraries),
            )
            .decrement;
            let free = dense_step(&dense, &flat, |_| true);
            let half_slope_through_inverse = 0.5
                * flat
                    .iter()
                    .zip(&free)
                    .map(|(g, d)| g * d.expect("solved"))
                    .sum::<f64>();
            assert!(half_slope_through_inverse > 0.0, "{label}");
            assert!(
                (got - half_slope_through_inverse).abs() <= 1e-9 * half_slope_through_inverse,
                "{label}: {got} against {half_slope_through_inverse}"
            );
        }
    }

    /// **With a bound binding, the decrement is `g·d − ½ dᵀ I d` for the box step `d`, and below the
    /// unbounded `½ gᵀ I⁻¹ g`**, on both matrices: the density's second shape placed a quarter of
    /// its free step inside the end that step points to, so the box stops it there (the case
    /// `a_parameter_carried_past_an_end_stops_there` checks the distances of).
    #[test]
    fn with_a_bound_binding_the_decrement_is_the_box_steps_worth() {
        let libraries = [1, 2, 1];
        let group_index = group_index_of(&libraries);
        let row = cohort::DENSITY_B;
        for (label, (mut sums, dense)) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 95)),
            ("arrow", an_arrow_of(&libraries, 97)),
        ] {
            let (slope, flat) = a_slope_of(&libraries, 96);
            sums.slope = slope;
            let n = flat.len();
            let free_step = dense_step(&dense, &flat, |_| true);
            let free = free_step[row].expect("solved");
            let end = if free < 0.0 {
                COHORT_BOUNDS[row].0
            } else {
                COHORT_BOUNDS[row].1
            };
            let value = end - 0.25 * free;
            let mut parameters = parameters_inside(&libraries);
            set_at(&mut parameters, row, value);
            let moved = end - value;
            let slope_left: Vec<f64> = (0..n)
                .map(|i| flat[i] - dense[i * n + row] * moved)
                .collect();
            let mut box_step = dense_step(&dense, &slope_left, |i| i != row);
            box_step[row] = Some(moved);
            let step = newton_step(&sums, &StandardErrors::of(&sums), &parameters, &group_index);
            assert_steps_match(
                &step.distances,
                &in_the_fits_order(&box_step, &libraries, &parameters),
                label,
            );
            let expected = dense_decrement(&dense, &flat, &box_step);
            let unbounded = dense_decrement(&dense, &flat, &free_step);
            assert!(
                (step.decrement - expected).abs() <= 1e-9 * expected.abs(),
                "{label}: {} against {expected}",
                step.decrement
            );
            assert!(
                expected > 0.0 && expected < unbounded,
                "{label}: the box step's {expected} against the unbounded {unbounded}"
            );
        }
    }

    /// `½ d_Kᵀ (Σ_KK)⁻¹ d_K` densely: `Σ` the inverse of `matrix` over every parameter with an
    /// error or a move, `K` those with both.
    fn dense_determined_decrement(matrix: &[f64], step: &[Option<f64>], has_error: &[bool]) -> f64 {
        let n = step.len();
        let solved: Vec<usize> = (0..n)
            .filter(|&i| step[i].is_some() || has_error[i])
            .collect();
        let covariance = inverse_of_positive_definite(
            &square_of(&solved, |row, column| matrix[row * n + column]),
            solved.len(),
        )
        .expect("positive definite");
        let determined: Vec<usize> = (0..solved.len())
            .filter(|&p| step[solved[p]].is_some() && has_error[solved[p]])
            .collect();
        let k = determined.len();
        let block = square_of(&determined, |p, q| covariance[p * solved.len() + q]);
        let precision = inverse_of_positive_definite(&block, k).expect("positive definite");
        let moved: Vec<f64> = determined
            .iter()
            .map(|&p| step[solved[p]].expect("determined"))
            .collect();
        0.5 * (0..k)
            .map(|a| moved[a] * (0..k).map(|b| precision[a * k + b] * moved[b]).sum::<f64>())
            .sum::<f64>()
    }

    /// **Measured over the determined parameters, the decrement is `½ d_Kᵀ (Σ_KK)⁻¹ d_K`**, the
    /// step's length in the determined parameters' own joint covariance — on the whole matrix, the
    /// arrow and an arrow of samples of one to three libraries, for an arbitrary step and three
    /// choices of which parameters have no error: the density's two shapes; those, a sample's
    /// mismapped rate and another's homozygote excess; and those, with a parameter that has an error
    /// but no move (profiled too) and one with neither (left out). Against a dense inverse, to 10⁻⁹
    /// relative. Holding the profiled parameters still instead — `½ d_Kᵀ I_KK d_K` — fails, and so
    /// does the whole step's length, which the profiled one is checked to stay below.
    #[test]
    fn the_determined_decrement_is_the_step_through_the_determined_covariance() {
        let libraries = [1, 2, 3, 1];
        let n = COHORT_PARAMETERS + libraries.iter().map(|&k| own_parameters(k)).sum::<usize>();
        let second_sample = COHORT_PARAMETERS + own_parameters(1);
        let third_sample = second_sample + own_parameters(2);
        let shapes = [cohort::DENSITY_A, cohort::DENSITY_B];
        let with_samples = [
            cohort::DENSITY_A,
            cohort::DENSITY_B,
            second_sample + sample::rate(1, 1),
            third_sample + sample::HOMOZYGOTE_EXCESS,
        ];
        for (label, (sums, dense)) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 121)),
            ("arrow", an_arrow_of(&libraries, 123)),
            ("arrow of one library a sample", an_arrow(4, 125)),
        ] {
            let n = if label.ends_with("sample") {
                COHORT_PARAMETERS + 4 * ONE_LIBRARY_SAMPLE_PARAMETERS
            } else {
                n
            };
            let layout = Layout::of(&sums);
            assert_eq!(layout.side, n, "{label}");
            let mut draw = uniform_draws(127);
            let step: Vec<Option<f64>> = (0..n).map(|_| Some(draw() - 0.5)).collect();
            for (case, undetermined) in [&shapes[..], &with_samples[..]].into_iter().enumerate() {
                let undetermined: Vec<usize> =
                    undetermined.iter().copied().filter(|&i| i < n).collect();
                let has_error: Vec<bool> = (0..n).map(|i| !undetermined.contains(&i)).collect();
                let got = determined_decrement(&sums, &layout, &step, &has_error)
                    .expect("some parameter is profiled");
                let expected = dense_determined_decrement(&dense, &step, &has_error);
                // `½ dᵀ I d`: the quadratic model's loss at the step, with no slope.
                let whole = -dense_decrement(&dense, &vec![0.0; n], &step);
                assert!(
                    (got - expected).abs() <= 1e-9 * expected,
                    "{label}, case {case}: {got} against {expected}"
                );
                assert!(
                    expected > 0.0 && expected < whole,
                    "{label}, case {case}: profiling lowers the step's length ({expected} \
                     against {whole})"
                );
            }
            // A parameter with an error and no move is profiled; one with neither stays put.
            let mut step_without = step.clone();
            let (error_no_move, neither) = (cohort::P_FIXED_ALT, cohort::NOISY_SHARE);
            step_without[error_no_move] = None;
            step_without[neither] = None;
            let has_error: Vec<bool> = (0..n)
                .map(|i| i != neither && !shapes.contains(&i))
                .collect();
            let got = determined_decrement(&sums, &layout, &step_without, &has_error)
                .expect("some parameter is profiled");
            let expected = dense_determined_decrement(&dense, &step_without, &has_error);
            assert!(
                (got - expected).abs() <= 1e-9 * expected,
                "{label}, with a parameter of each kind left: {got} against {expected}"
            );
        }
    }

    /// **When every parameter the step moves has an error the decrement is the old one, bit for
    /// bit, and when none has one it is zero**: on the whole matrix and the arrow, the decrement
    /// `newton_step` gives equals `g·d − ½ dᵀ I d` of its own box step exactly; with every error
    /// taken away the determined decrement is exactly zero — nothing is left to measure the step
    /// in.
    #[test]
    fn with_every_error_the_decrement_is_the_old_one_and_with_none_it_is_zero() {
        let libraries = [1, 2, 3, 1];
        let group_index = group_index_of(&libraries);
        for (label, (mut sums, _), slope_seed) in [
            ("whole matrix", a_whole_matrix_of(&libraries, 91), 92),
            ("arrow", an_arrow_of(&libraries, 93), 94),
        ] {
            let (slope, _) = a_slope_of(&libraries, slope_seed);
            sums.slope = slope;
            let parameters = parameters_inside(&libraries);
            let errors = StandardErrors::of(&sums);
            assert!(
                errors
                    .cohort
                    .iter()
                    .chain(errors.samples.iter().flatten())
                    .all(|error| error.value().is_some()),
                "{label}: every parameter has an error"
            );
            let layout = Layout::of(&sums);
            let step = bounded_step(
                &sums,
                &layout,
                &layout.values_of(&parameters, &group_index),
                &layout.bounds(),
            );
            let old = decrement_of(&sums, &layout, &step.moves);
            let got = newton_step(&sums, &errors, &parameters, &group_index).decrement;
            assert_eq!(got.to_bits(), old.to_bits(), "{label}: {got} against {old}");
            assert_eq!(
                determined_decrement(&sums, &layout, &step.moves, &vec![false; layout.side]),
                Some(0.0),
                "{label}"
            );
        }
    }

    /// A draw in [0, 1) from a fixed sequence.
    fn uniform_draws(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// **An information matrix summed the way a pass sums one**, densely in the layout: `positions`
    /// outer products of score vectors, each position scoring the cohort's eight parameters and one
    /// sample's own, so two samples are never paired and the matrix is an arrow; then `floor` times
    /// its largest diagonal entry added to the diagonal. With fewer positions than parameters it is
    /// singular but for that floor, so its condition number is about `1 / floor`.
    fn summed_from_positions(
        libraries: &[usize],
        positions: usize,
        floor: f64,
        draw: &mut impl FnMut() -> f64,
    ) -> Vec<f64> {
        let sizes: Vec<usize> = libraries.iter().map(|&k| own_parameters(k)).collect();
        let n = COHORT_PARAMETERS + sizes.iter().sum::<usize>();
        let mut dense = vec![0.0; n * n];
        for _ in 0..positions {
            let s = ((draw() * sizes.len() as f64) as usize).min(sizes.len() - 1);
            let first = COHORT_PARAMETERS + sizes[..s].iter().sum::<usize>();
            let touched: Vec<usize> = (0..COHORT_PARAMETERS)
                .chain(first..first + sizes[s])
                .collect();
            let score: Vec<f64> = touched.iter().map(|_| draw() - 0.5).collect();
            for (a, &row) in touched.iter().enumerate() {
                for (b, &column) in touched.iter().enumerate() {
                    dense[row * n + column] += score[a] * score[b];
                }
            }
        }
        let largest = (0..n).map(|i| dense[i * n + i]).fold(0.0_f64, f64::max);
        for i in 0..n {
            dense[i * n + i] += floor * largest;
        }
        dense.iter().map(|entry| entry * 1e6).collect()
    }

    /// **The box solve stays inside the box, is worth at least nothing, and finishes** — on 400
    /// random cases: one to five samples of one to three libraries, on the whole matrix and on the
    /// arrow, information well conditioned (three times as many positions as parameters) or nearly
    /// singular (half as many, condition numbers about 10⁶ and 10¹²), random slopes, and every
    /// parameter placed on its lower end, on its upper end, just inside one or in the middle of its
    /// interval. Every move lands inside the interval (to the rounding of `value + move`), the step's
    /// worth `g·d − ½ dᵀ I d` is not negative, and the solve says it finished. On the well-conditioned
    /// cases with every parameter solved for, the step is also the box's maximum: the slope left at
    /// the step, `g − I d`, is zero for a parameter strictly inside, and points out of the interval
    /// for one on an end. The solve that preceded the primal method failed this on cases like the
    /// review's two-sample fixtures (`fit_precision.md` §2, amendment of 2026-10-17).
    #[test]
    fn the_box_step_stays_inside_the_box_and_is_worth_at_least_nothing() {
        let mut draw = uniform_draws(0xB0C5_0000_0000_0017);
        let mut maxima_checked = 0;
        let mut bounds_binding = 0;
        for case in 0..400 {
            let samples = 1 + (draw() * 5.0) as usize;
            let libraries: Vec<usize> = (0..samples).map(|_| 1 + (draw() * 3.0) as usize).collect();
            let n = COHORT_PARAMETERS + libraries.iter().map(|&k| own_parameters(k)).sum::<usize>();
            let (positions, floor) = match case % 3 {
                0 => (3 * n, 1e-9),
                1 => (n / 2, 1e-6),
                _ => (n / 2, 1e-12),
            };
            let dense = summed_from_positions(&libraries, positions, floor, &mut draw);
            let whole = case % 2 == 0;
            let mut sums = if whole {
                with_the_whole_matrix(&libraries, &dense)
            } else {
                blocks_of_dense(&libraries, &dense)
            };
            let (slope, flat) = a_slope_of(&libraries, 1_000 + case);
            sums.slope = slope;
            let layout = Layout::of(&sums);
            let bounds = layout.bounds();
            let values: Vec<f64> = bounds
                .iter()
                .map(|&(low, high)| {
                    let inside = crate::float::powf(10.0, -7.0 + 3.0 * draw()) * (high - low);
                    match (draw() * 5.0) as usize {
                        0 => low,
                        1 => high,
                        2 => low + inside,
                        3 => high - inside,
                        _ => low + (0.2 + 0.6 * draw()) * (high - low),
                    }
                })
                .collect();
            let label = format!(
                "case {case}: libraries {libraries:?}, {positions} positions, floor {floor:e}, \
                 {}",
                if whole { "whole matrix" } else { "arrow" }
            );
            let step = bounded_step(&sums, &layout, &values, &bounds);
            assert!(step.finished, "{label}: the solve did not finish");
            let mut on_an_end = vec![None; n];
            for (j, moved) in step.moves.iter().enumerate() {
                let Some(moved) = moved else { continue };
                let (low, high) = bounds[j];
                let landing = values[j] + moved;
                let rounding = 4.0 * f64::EPSILON * values[j].abs().max(high.abs());
                assert!(
                    landing >= low - rounding && landing <= high + rounding,
                    "{label}: parameter {j} at {} moved {moved} to {landing}, outside {:?}",
                    values[j],
                    bounds[j]
                );
                if (landing - low).abs() <= rounding {
                    on_an_end[j] = Some(low);
                } else if (landing - high).abs() <= rounding {
                    on_an_end[j] = Some(high);
                }
            }
            let worth = dense_decrement(&dense, &flat, &step.moves);
            assert!(worth >= 0.0, "{label}: the step is worth {worth}");
            if case % 3 == 0 && step.moves.iter().all(Option::is_some) {
                maxima_checked += 1;
                let moved: Vec<f64> = step.moves.iter().map(|m| m.unwrap_or(0.0)).collect();
                for j in 0..n {
                    let curved: f64 = (0..n).map(|k| dense[j * n + k] * moved[k]).sum();
                    let left = flat[j] - curved;
                    let scale = 1e-7 * (flat[j].abs() + curved.abs()) + 1e-9;
                    match on_an_end[j] {
                        None => assert!(
                            left.abs() <= scale,
                            "{label}: parameter {j} inside its interval with slope {left}"
                        ),
                        Some(end) => {
                            bounds_binding += 1;
                            let outwards = if end == bounds[j].0 { -left } else { left };
                            assert!(
                                outwards >= -scale,
                                "{label}: parameter {j} on {end} with slope {left} pointing in"
                            );
                        }
                    }
                }
            }
        }
        assert!(
            maxima_checked > 50 && bounds_binding > 100,
            "the cases exercise the box: {maxima_checked} maxima checked, {bounds_binding} ends"
        );
    }

    /// The dense inverse of `full` (side `n`) without the parameters `removed` says, restricted to
    /// the cohort's eight: `None` for an entry whose row or column was removed.
    fn dense_cohort_covariance_without(
        full: &[f64],
        n: usize,
        removed: impl Fn(usize) -> bool,
    ) -> CohortCovariance {
        let kept: Vec<usize> = (0..n).filter(|&i| !removed(i)).collect();
        let m = kept.len();
        let inverse = inverse_of_positive_definite(
            &square_of(&kept, |row, column| full[row * n + column]),
            m,
        )
        .expect("positive definite");
        let at = |i: usize| kept.iter().position(|&k| k == i);
        std::array::from_fn(|j| std::array::from_fn(|k| Some(inverse[at(j)? * m + at(k)?])))
    }

    /// Every entry of `got` against `expected`, to 10⁻¹² of the larger diagonal entry of its row
    /// and column; `None` exactly where `expected` has none.
    fn assert_covariances_match(got: &CohortCovariance, expected: &CohortCovariance, what: &str) {
        for j in 0..COHORT_PARAMETERS {
            for k in 0..COHORT_PARAMETERS {
                match (got[j][k], expected[j][k]) {
                    (Some(got), Some(dense)) => {
                        let scale = (expected[j][j].expect("diagonal")
                            * expected[k][k].expect("diagonal"))
                        .sqrt();
                        assert!(
                            (got - dense).abs() <= 1e-12 * scale,
                            "{what}: ({j}, {k}) is {got:e}, the dense inverse {dense:e}"
                        );
                    }
                    (None, None) => {}
                    (got, dense) => panic!("{what}: ({j}, {k}) is {got:?}, dense {dense:?}"),
                }
            }
        }
    }

    /// **The cohort's covariance is the dense inverse restricted to the cohort**, every other
    /// parameter accounted for: on the whole matrix of samples of one, two, three and one library
    /// (26 parameters), on an arrow of five one-library samples (23), and on an arrow of the same
    /// mixed libraries — every entry, off the diagonal too, to 10⁻¹² of its row's and column's
    /// errors. **Its diagonal's square roots are the cohort's standard errors**, to 10⁻¹² relative.
    /// Solving the unit slopes with the samples' own parameters held — the cohort's block inverted
    /// alone — fails (checked).
    #[test]
    fn the_cohort_covariance_is_the_dense_inverse_restricted_to_the_cohort() {
        let mixed = [1, 2, 3, 1];
        let cases = [
            ("the whole matrix", a_whole_matrix_of(&mixed, 71), 26),
            ("an arrow", an_arrow(5, 73), 23),
            ("an arrow of mixed libraries", an_arrow_of(&mixed, 79), 26),
        ];
        for (what, (sums, dense), n) in cases {
            let covariance = cohort_covariance(&sums);
            assert_covariances_match(
                &covariance,
                &dense_cohort_covariance_without(&dense, n, |_| false),
                what,
            );
            let errors = StandardErrors::of(&sums);
            for (slot, (row, error)) in covariance.iter().zip(&errors.cohort).enumerate() {
                let from_covariance = row[slot].expect("solved").sqrt();
                let error = error.value().expect("an error");
                assert!(
                    (from_covariance - error).abs() <= 1e-12 * error,
                    "{what}: slot {slot}, {from_covariance:e} against the error {error:e}"
                );
            }
        }
    }

    /// **A cohort parameter the errors leave out has no row or column**: on the arrow and on the
    /// whole matrix, the duplicated class's three without information and the density's second
    /// shape made a copy of its first (so not told apart) are `None` across, and every other entry is
    /// the dense inverse's without them. **One whose error came out wider than its range keeps its
    /// entries**: the duplicated share's information cut ten-thousandfold.
    #[test]
    fn a_cohort_parameter_without_an_error_has_no_covariance_unless_it_is_only_wide() {
        let libraries = [1, 2, 1];
        let n = COHORT_PARAMETERS + 11;
        let (_, mut dense) = an_arrow_of(&libraries, 83);
        for i in cohort::DUPLICATED_SHARE..COHORT_PARAMETERS {
            for k in 0..n {
                dense[i * n + k] = 0.0;
                dense[k * n + i] = 0.0;
            }
        }
        let (copy, of) = (cohort::DENSITY_B, cohort::DENSITY_A);
        for k in 0..n {
            dense[copy * n + k] = dense[of * n + k];
        }
        for k in 0..n {
            dense[k * n + copy] = dense[k * n + of];
        }
        let removed =
            |i: usize| i == copy || (cohort::DUPLICATED_SHARE..COHORT_PARAMETERS).contains(&i);
        let expected = dense_cohort_covariance_without(&dense, n, removed);
        let arrow = blocks_of_dense(&libraries, &dense);
        let whole = with_the_whole_matrix(&libraries, &dense);
        for (what, sums) in [("the arrow", &arrow), ("the whole matrix", &whole)] {
            let errors = StandardErrors::of(sums);
            assert_eq!(errors.cohort[copy], StandardError::NotIdentified, "{what}");
            assert_eq!(
                errors.cohort[cohort::DUPLICATED_SHARE..],
                [StandardError::NoInformation; 3],
                "{what}"
            );
            let covariance = cohort_covariance(sums);
            for slot in (0..COHORT_PARAMETERS).filter(|&slot| removed(slot)) {
                assert!(
                    covariance[slot].iter().all(Option::is_none),
                    "{what}: row {slot}"
                );
                assert!(
                    covariance.iter().all(|row| row[slot].is_none()),
                    "{what}: column {slot}"
                );
            }
            assert_covariances_match(&covariance, &expected, what);
        }

        let (mut sums, dense) = an_arrow(2, 61);
        let at = cohort::DUPLICATED_SHARE;
        let shrink = 1e-2;
        for j in 0..COHORT_PARAMETERS {
            sums.cohort[at * COHORT_PARAMETERS + j] *= shrink;
            sums.cohort[j * COHORT_PARAMETERS + at] *= shrink;
        }
        for block in &mut sums.sample_cohort_blocks {
            for row in 0..ONE_LIBRARY_SAMPLE_PARAMETERS {
                block[row * COHORT_PARAMETERS + at] *= shrink;
            }
        }
        assert!(matches!(
            StandardErrors::of(&sums).cohort[at],
            StandardError::WiderThanItsRange(_)
        ));
        let n = COHORT_PARAMETERS + 2 * ONE_LIBRARY_SAMPLE_PARAMETERS;
        let mut shrunk = dense;
        for k in 0..n {
            shrunk[at * n + k] *= shrink;
            shrunk[k * n + at] *= shrink;
        }
        assert_covariances_match(
            &cohort_covariance(&sums),
            &dense_cohort_covariance_without(&shrunk, n, |_| false),
            "a wide error",
        );
    }
}
