//! **Every repeat class of one motif length fitted at once, each number a smooth curve across the
//! reference repeat count** — penalised expectation maximisation, a research prototype
//! (`doc/devel/ng/research/str_fit_research_plan.md`, avenue D).
//!
//! The counting fit ([`composite_fit`](super::composite_fit)) fits each class — one motif length,
//! one reference repeat count — on its own, and the curves across repeat count are drawn afterwards
//! through its answers. Here each round gathers every class's expected counts exactly as the
//! counting fit does ([`expected_counts`]), and then chooses all the classes' numbers together: each
//! number is the best fit to the counts **subject to a penalty on how sharply it bends from one
//! repeat count to the next**. A thin class then borrows from its neighbours during the fit, by as
//! much as its own counts leave undecided, instead of being blended with a curve afterwards.
//!
//! **How sharply a curve may bend is measured, not chosen.** Each curve's penalty strength is the
//! one that maximises an approximation to the marginal likelihood of its counts (the Laplace
//! approximation that restricted maximum likelihood rests on, Wood 2011), re-chosen every round.
//!
//! **A strong, well-measured kink survives.** The penalty is quadratic in each bend up to
//! [`PenalisedFitConfig::kink_knee`] standard deviations of what the curve's smoothness predicts,
//! and grows only logarithmically beyond it, so a class whose own reads insist on a step pays a
//! bounded price for it. The weight on such a bend is divided by `(bend / knee)²`, the rule the
//! slippage curve's blend already uses (`str_slippage_level_curve.md` §7.2).

use std::time::Instant;

use super::composite_fit::{ClassCounts, SamplePairs, expected_counts};
use super::*;

/// **How the penalised fit runs.**
#[derive(Debug, Clone, PartialEq)]
pub struct PenalisedFitConfig {
    /// Where every class's stutter starts. The allele length variation starts flat and the
    /// concentration at one, as in the counting fit.
    pub start: Slippage,
    /// At most this many rounds.
    pub max_rounds: usize,
    /// Stop once no number moves by more than this in a round (measured as the counting fit
    /// measures it).
    pub tolerance: f64,
    /// How strongly each curve is smoothed.
    pub smoothing: Smoothing,
    /// **How many standard deviations a bend may reach before its penalty stops growing as its
    /// square**: 2.5, the knee of the slippage curve's blend. `f64::INFINITY` makes every bend
    /// purely quadratic.
    pub kink_knee: f64,
    /// **How many times the between-locus share is read off pairs of samples at a locus**, each time
    /// followed by penalised counting with it held. Zero reads it only from each sample's own two
    /// copies, which says almost nothing in an inbred cohort.
    pub pair_alternations: usize,
}

impl Default for PenalisedFitConfig {
    fn default() -> Self {
        Self {
            start: Slippage {
                level: 0.2,
                shorter_share: 0.5,
                fall_off: 0.5,
            },
            max_rounds: 500,
            tolerance: 1e-6,
            smoothing: Smoothing::MarginalLikelihood,
            kink_knee: 2.5,
            pair_alternations: 0,
        }
    }
}

/// **How a curve's penalty strength is set.**
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Smoothing {
    /// Chosen every round by the Laplace approximation to the marginal likelihood of the curve's
    /// expected counts, searched between 10⁻⁴ and 10¹⁰.
    MarginalLikelihood,
    /// Held at this strength for every curve. A strength near zero fits every class on its own,
    /// which is the counting fit.
    Fixed(f64),
}

/// The penalty strengths searched, as natural logarithms.
const LN_SMOOTHING_RANGE: (f64, f64) = (-9.210_340_371_976_182, 23.025_850_929_940_457);

/// **A ridge on every unbounded number**, so that a number with no counts behind it at any class
/// stays finite: a class with no slipped reads would otherwise send the logit of its slip rate to
/// minus infinity. At 10⁻⁶ it moves a number with even one count behind it by less than one part in
/// a million.
const RIDGE: f64 = 1e-6;

/// Unbounded numbers are held within ±this many units, beyond which every probability here is
/// zero or one to working precision.
const LARGEST_UNBOUNDED: f64 = 40.0;

/// One class's fitted numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct PenalisedClass {
    /// Which class.
    pub stratum: Stratum,
    /// The class's stutter.
    pub slippage: Slippage,
    /// How the class's chromosomes are spread over the allele lengths, `-allele_span ..=
    /// +allele_span`.
    pub allele_length_variation: Vec<f64>,
    /// The Dirichlet concentration, as the counting fit reports it.
    pub concentration: f64,
    /// The same, as the between-locus share `1 / (1 + concentration)`.
    pub between_locus_share: f64,
}

/// **One number's curve across the repeat counts of a motif length**, on its natural scale.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterCurve {
    /// The repeat count of `values[0]`; the curve runs over consecutive repeat counts from here.
    pub first_repeats: u64,
    /// The curve's value at each repeat count, from the lowest class with loci to the highest. A
    /// repeat count with no loci between them takes the curve's value there.
    pub values: Vec<f64>,
    /// The penalty strength the curve was fitted at.
    pub smoothing: f64,
    /// **The repeat counts where the curve was let bend**: the middle repeat count of each bend
    /// whose penalty was eased below half its full weight because the counts insisted on it.
    pub eased_at: Vec<u64>,
}

impl ParameterCurve {
    /// **The curve at `repeats`**, held at the nearest end outside the range it was fitted over.
    #[must_use]
    pub fn at(&self, repeats: u64) -> f64 {
        let last = self.first_repeats + self.values.len() as u64 - 1;
        let index = repeats.clamp(self.first_repeats, last) - self.first_repeats;
        self.values[index as usize]
    }
}

/// **The curves of one motif length.**
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothCurves {
    pub level: ParameterCurve,
    pub shorter_share: ParameterCurve,
    pub fall_off: ParameterCurve,
    pub between_locus_share: ParameterCurve,
    /// One curve per allele class, `-allele_span ..= +allele_span`: that class's share. **The
    /// shares' curves share one penalty strength**, and each offset's bends are eased on their own.
    pub allele_length_variation: Vec<ParameterCurve>,
}

/// What [`fit_period_by_penalised_counting`] returns.
#[derive(Debug, Clone, PartialEq)]
pub struct PenalisedPeriodFit {
    /// One entry per class handed in, in the same order: `None` for a class with no reads at any
    /// locus, which is not fitted. The curves still give such a class a value
    /// ([`ParameterCurve::at`]).
    pub classes: Vec<Option<PenalisedClass>>,
    /// `None` when no class had reads.
    pub curves: Option<SmoothCurves>,
    /// How many rounds ran, and the largest change of any number in the last.
    pub rounds: usize,
    pub last_change: f64,
    /// The counting fit's objective summed over the classes, at the numbers of the last round.
    pub composite_log_likelihood: f64,
    /// Wall time spent gathering expected counts and choosing numbers, over every round.
    pub expectation_seconds: f64,
    pub maximisation_seconds: f64,
}

/// **Fits every class of one motif length at once, by penalised counting.**
///
/// `classes` are the classes of one motif length, any order, each repeat count at most once.
/// `homozygote_excess` is indexed by sample. Every library is scored at its class's stutter
/// (multipliers are not fitted here). The allele window is `config.allele_span`.
///
/// With [`PenalisedFitConfig::pair_alternations`] above zero, the between-locus share is then read
/// off pairs of samples at a locus, as [`fit_class_by_counting_with_pairs`](super::composite_fit::fit_class_by_counting_with_pairs)
/// does, its curve penalised like every other: each round of pairs alternates with penalised
/// counting at the shares held.
///
/// # Panics
///
/// When the classes are of more than one motif length, or two share a repeat count.
#[must_use]
pub fn fit_period_by_penalised_counting(
    classes: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    fit: &PenalisedFitConfig,
) -> PenalisedPeriodFit {
    if let Some(first) = classes.first() {
        assert!(
            classes
                .iter()
                .all(|class| class.stratum.period == first.stratum.period),
            "one motif length at a time"
        );
    }
    let Some(mut period) = Period::new(classes, config, fit) else {
        return PenalisedPeriodFit {
            classes: vec![None; classes.len()],
            curves: None,
            rounds: 0,
            last_change: 0.0,
            composite_log_likelihood: 0.0,
            expectation_seconds: 0.0,
            maximisation_seconds: 0.0,
        };
    };
    period.count(classes, homozygote_excess, config, fit, false);
    for alternation in 0..fit.pair_alternations {
        let started = Instant::now();
        let Some(tables) = period.pair_tables(classes, homozygote_excess, config) else {
            break;
        };
        let spreads = period.between_curve.refit(&tables, fit, alternation == 0);
        let mut moved: f64 = 0.0;
        for (point, class) in period.class_at_point.iter().enumerate() {
            if class.is_some() {
                let concentration = float::exp(-spreads[point]).clamp(1e-3, 1e4);
                moved = moved
                    .max((float::ln(concentration) - float::ln(period.concentration[point])).abs());
                period.concentration[point] = concentration;
            }
        }
        period.maximisation_seconds += started.elapsed().as_secs_f64();
        period.count(classes, homozygote_excess, config, fit, true);
        if moved < 0.01 {
            break;
        }
    }
    period.into_fit(classes)
}

/// The ln-spread grid a class's pairwise likelihood is tabulated on, as the logit of the
/// between-locus share (`−ln spread`): 57 points from a spread of 10⁴ to 10⁻³.
fn between_share_grid() -> Vec<f64> {
    let (low, high) = (-float::ln(1e4), -float::ln(1e-3));
    (0..57)
        .map(|step| low + (high - low) * f64::from(step) / 56.0)
        .collect()
}

/// **One motif length's classes as the fit carries them from round to round**: every number at
/// every repeat count of the grid, and each curve's state.
struct Period {
    /// Per grid point (consecutive repeat counts from `lowest`), the class there, if any.
    class_at_point: Vec<Option<usize>>,
    lowest: u64,
    alleles: usize,
    reference_allele: usize,
    slippage: Vec<Slippage>,
    shares: Vec<Vec<f64>>,
    concentration: Vec<f64>,
    level_curve: CurveState,
    shorter_curve: CurveState,
    fall_off_curve: CurveState,
    between_curve: CurveState,
    shares_curve: CurveState,
    rounds: usize,
    last_change: f64,
    composite_log_likelihood: f64,
    expectation_seconds: f64,
    maximisation_seconds: f64,
}

impl Period {
    /// The grid over the classes with reads, every number at its start; `None` when no class has
    /// reads.
    fn new(
        classes: &[StratumEvidence],
        config: &SsrFitConfig,
        fit: &PenalisedFitConfig,
    ) -> Option<Self> {
        let with_reads: Vec<usize> = (0..classes.len())
            .filter(|&index| classes[index].rows_with_reads() > 0)
            .collect();
        let repeats = |index: &usize| classes[*index].stratum.reference_repeats;
        let lowest = with_reads.iter().map(repeats).min()?;
        let highest = with_reads.iter().map(repeats).max()?;
        let points = (highest - lowest + 1) as usize;
        let mut class_at_point: Vec<Option<usize>> = vec![None; points];
        for &index in &with_reads {
            assert!(
                class_at_point[(repeats(&index) - lowest) as usize]
                    .replace(index)
                    .is_none(),
                "two classes at one repeat count"
            );
        }
        let alleles = allele_classes_of(config);
        Some(Self {
            class_at_point,
            lowest,
            alleles,
            reference_allele: config.allele_span as usize,
            slippage: vec![fit.start; points],
            shares: vec![vec![1.0 / alleles as f64; alleles]; points],
            concentration: vec![1.0; points],
            level_curve: CurveState::new(points, 1),
            shorter_curve: CurveState::new(points, 1),
            fall_off_curve: CurveState::new(points, 1),
            between_curve: CurveState::new(points, 1),
            shares_curve: CurveState::new(points, alleles - 1),
            rounds: 0,
            last_change: f64::INFINITY,
            composite_log_likelihood: f64::NEG_INFINITY,
            expectation_seconds: 0.0,
            maximisation_seconds: 0.0,
        })
    }

    fn points(&self) -> usize {
        self.class_at_point.len()
    }

    /// **Rounds of penalised counting** until no number moves by more than the tolerance or the
    /// rounds run out. With `spread_held`, the concentration stays where it is.
    fn count(
        &mut self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        config: &SsrFitConfig,
        fit: &PenalisedFitConfig,
        spread_held: bool,
    ) {
        let points = self.points();
        let alleles = self.alleles;
        let mut rounds_here = 0;
        while self.rounds < fit.max_rounds {
            self.rounds += 1;
            rounds_here += 1;
            let first_round = rounds_here == 1;
            let started = Instant::now();
            let counts: Vec<Option<ClassCounts>> = self
                .class_at_point
                .iter()
                .enumerate()
                .map(|(point, class)| {
                    class.map(|index| {
                        expected_counts(
                            &classes[index],
                            self.slippage[point],
                            &self.shares[point],
                            self.concentration[point],
                            homozygote_excess,
                            config,
                        )
                    })
                })
                .collect();
            self.expectation_seconds += started.elapsed().as_secs_f64();
            let started = Instant::now();
            self.composite_log_likelihood = counts.iter().flatten().map(|c| c.log_likelihood).sum();
            let binomial = |successes: &dyn Fn(&ClassCounts) -> f64,
                            trials: &dyn Fn(&ClassCounts) -> f64| {
                let (mut ys, mut ns) = (vec![0.0; points], vec![0.0; points]);
                for (point, count) in counts.iter().enumerate() {
                    if let Some(count) = count {
                        let n = trials(count).max(0.0);
                        ns[point] = n;
                        ys[point] = successes(count).clamp(0.0, n);
                    }
                }
                CurveData::Binomial {
                    successes: ys,
                    trials: ns,
                }
            };
            let level = self.level_curve.refit(
                &binomial(&|c| c.slipped, &|c| c.slipped + c.exact),
                fit,
                first_round,
            );
            let shorter = self.shorter_curve.refit(
                &binomial(&|c| c.shorter, &|c| c.shorter + c.longer),
                fit,
                first_round,
            );
            // A slip of k steps has probability (1 − f) f^(k − 1): k − 1 "one step further" against
            // one "stopped".
            let fall_off = self.fall_off_curve.refit(
                &binomial(&|c| c.steps - c.stepped, &|c| c.steps),
                fit,
                first_round,
            );
            let between = if spread_held {
                None
            } else {
                Some(self.between_curve.refit(
                    &binomial(&|c| c.shared_in_tract, &|c| c.pairs_not_by_descent),
                    fit,
                    first_round,
                ))
            };
            let draws: Vec<Vec<f64>> = counts
                .iter()
                .map(|count| count.as_ref().map_or_else(Vec::new, |c| c.draws.clone()))
                .collect();
            let share_logits = self.shares_curve.refit(
                &CurveData::Multinomial {
                    counts: draws,
                    reference: self.reference_allele,
                },
                fit,
                first_round,
            );
            self.maximisation_seconds += started.elapsed().as_secs_f64();

            let mut change: f64 = 0.0;
            for point in 0..points {
                if self.class_at_point[point].is_none() {
                    continue;
                }
                let new_slippage = Slippage {
                    level: logistic(level[point]),
                    shorter_share: logistic(shorter[point]),
                    fall_off: logistic(fall_off[point]).clamp(0.0, 0.99),
                };
                let new_concentration = between.as_ref().map_or(self.concentration[point], |b| {
                    let share = logistic(b[point]).clamp(1.0 / (1.0 + 1e4), 1.0 / (1.0 + 1e-3));
                    (1.0 / share - 1.0).clamp(1e-3, 1e4)
                });
                let new_shares = softmax_with_reference(
                    &share_logits[point * (alleles - 1)..(point + 1) * (alleles - 1)],
                    self.reference_allele,
                );
                let old = self.slippage[point];
                change = new_shares
                    .iter()
                    .zip(&self.shares[point])
                    .map(|(new, old)| (new - old).abs())
                    .chain([
                        (new_slippage.level - old.level).abs(),
                        (new_slippage.shorter_share - old.shorter_share).abs(),
                        (new_slippage.fall_off - old.fall_off).abs(),
                        (float::ln(new_concentration) - float::ln(self.concentration[point])).abs()
                            * 0.01,
                    ])
                    .fold(change, f64::max);
                self.slippage[point] = new_slippage;
                self.shares[point] = new_shares;
                self.concentration[point] = new_concentration;
            }
            self.last_change = change;
            if self.last_change < fit.tolerance {
                break;
            }
        }
    }

    /// **Every class's pairwise log-likelihood tabulated over the between-locus share**, at the
    /// class's current stutter and shares; `None` when no class has a locus with two samples.
    fn pair_tables(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        config: &SsrFitConfig,
    ) -> Option<CurveData> {
        let grid = between_share_grid();
        let values: Vec<Vec<f64>> = self
            .class_at_point
            .iter()
            .enumerate()
            .map(|(point, class)| {
                let Some(index) = class else {
                    return Vec::new();
                };
                let pairs = SamplePairs::of(
                    &classes[*index],
                    self.slippage[point],
                    homozygote_excess,
                    config,
                );
                if pairs.is_empty() {
                    return Vec::new();
                }
                grid.iter()
                    .map(|logit| pairs.log_likelihood(&self.shares[point], -logit))
                    .collect()
            })
            .collect();
        values
            .iter()
            .any(|table| !table.is_empty())
            .then_some(CurveData::Tabulated { grid, values })
    }

    fn into_fit(self, classes: &[StratumEvidence]) -> PenalisedPeriodFit {
        let alleles = self.alleles;
        let reference_allele = self.reference_allele;
        let lowest = self.lowest;
        let point_of = |index: usize| {
            self.class_at_point
                .iter()
                .position(|class| *class == Some(index))
        };
        let fitted_classes = (0..classes.len())
            .map(|index| {
                let point = point_of(index)?;
                Some(PenalisedClass {
                    stratum: classes[index].stratum,
                    slippage: self.slippage[point],
                    allele_length_variation: self.shares[point].clone(),
                    concentration: self.concentration[point],
                    between_locus_share: 1.0 / (1.0 + self.concentration[point]),
                })
            })
            .collect();
        let scalar_curve = |state: &CurveState, transform: &dyn Fn(f64) -> f64| ParameterCurve {
            first_repeats: lowest,
            values: state.values.iter().map(|value| transform(*value)).collect(),
            smoothing: float::exp(state.ln_smoothing),
            eased_at: state.eased_points(0, lowest),
        };
        let share_values: Vec<Vec<f64>> = (0..self.points())
            .map(|point| {
                softmax_with_reference(
                    &self.shares_curve.values[point * (alleles - 1)..(point + 1) * (alleles - 1)],
                    reference_allele,
                )
            })
            .collect();
        let curves = SmoothCurves {
            level: scalar_curve(&self.level_curve, &logistic),
            shorter_share: scalar_curve(&self.shorter_curve, &logistic),
            fall_off: scalar_curve(&self.fall_off_curve, &|value| {
                logistic(value).clamp(0.0, 0.99)
            }),
            between_locus_share: scalar_curve(&self.between_curve, &|value| {
                logistic(value).clamp(1.0 / (1.0 + 1e4), 1.0 / (1.0 + 1e-3))
            }),
            allele_length_variation: (0..alleles)
                .map(|allele| ParameterCurve {
                    first_repeats: lowest,
                    values: share_values.iter().map(|shares| shares[allele]).collect(),
                    smoothing: float::exp(self.shares_curve.ln_smoothing),
                    eased_at: match allele.cmp(&reference_allele) {
                        std::cmp::Ordering::Equal => Vec::new(),
                        std::cmp::Ordering::Less => self.shares_curve.eased_points(allele, lowest),
                        std::cmp::Ordering::Greater => {
                            self.shares_curve.eased_points(allele - 1, lowest)
                        }
                    },
                })
                .collect(),
        };
        PenalisedPeriodFit {
            classes: fitted_classes,
            curves: Some(curves),
            rounds: self.rounds,
            last_change: self.last_change,
            composite_log_likelihood: self.composite_log_likelihood,
            expectation_seconds: self.expectation_seconds,
            maximisation_seconds: self.maximisation_seconds,
        }
    }
}

fn logistic(value: f64) -> f64 {
    1.0 / (1.0 + float::exp(-value))
}

/// The shares from the log-ratios of every allele class but `reference` to it.
fn softmax_with_reference(logits: &[f64], reference: usize) -> Vec<f64> {
    let full: Vec<f64> = (0..=logits.len())
        .map(|allele| match allele.cmp(&reference) {
            std::cmp::Ordering::Less => logits[allele],
            std::cmp::Ordering::Equal => 0.0,
            std::cmp::Ordering::Greater => logits[allele - 1],
        })
        .collect();
    let largest = full.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exponentials: Vec<f64> = full.iter().map(|x| float::exp(x - largest)).collect();
    let total: f64 = exponentials.iter().sum();
    exponentials.iter().map(|x| x / total).collect()
}

// ---------------------------------------------------------------------
// One curve's penalised fit
// ---------------------------------------------------------------------

/// **The counts one curve is fitted to**, per grid point (a repeat count). A point whose counts are
/// all zero carries no evidence and takes the curve's value.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum CurveData {
    /// `successes` out of `trials` at each point; the unbounded number is the logit of the
    /// probability.
    Binomial {
        successes: Vec<f64>,
        trials: Vec<f64>,
    },
    /// Counts in each category at each point (an empty vector for a point with none); the
    /// unbounded numbers are each category's log-ratio to `reference`'s.
    Multinomial {
        counts: Vec<Vec<f64>>,
        reference: usize,
    },
    /// A log-likelihood tabulated at each point over one unbounded number, on `grid` (ascending,
    /// evenly spaced); an empty table for a point with none. Read between the nodes by a cubic
    /// with the table's slopes, and held at the end node beyond them.
    Tabulated {
        grid: Vec<f64>,
        values: Vec<Vec<f64>>,
    },
}

impl CurveData {
    /// **One point's log-likelihood at `values`**, with its gradient and its information (minus
    /// the second derivative), the information `block × block` row by row. Zero, with no
    /// information, at a point with no counts.
    fn evaluate(
        &self,
        point: usize,
        values: &[f64],
        gradient: &mut [f64],
        information: &mut [f64],
    ) -> f64 {
        gradient.fill(0.0);
        information.fill(0.0);
        match self {
            CurveData::Binomial { successes, trials } => {
                let (y, n) = (successes[point], trials[point]);
                if n <= 0.0 {
                    return 0.0;
                }
                let eta = values[0];
                let p = logistic(eta);
                gradient[0] = y - n * p;
                information[0] = n * p * (1.0 - p);
                // y·η − n·ln(1 + e^η), computed stably.
                let ln_one_plus = if eta > 0.0 {
                    eta + float::ln_1p(float::exp(-eta))
                } else {
                    float::ln_1p(float::exp(eta))
                };
                y * eta - n * ln_one_plus
            }
            CurveData::Multinomial { counts, reference } => {
                let counts = &counts[point];
                let total: f64 = counts.iter().sum();
                if total <= 0.0 {
                    return 0.0;
                }
                let block = values.len();
                let shares = softmax_with_reference(values, *reference);
                let others: Vec<usize> = (0..counts.len()).filter(|a| a != reference).collect();
                let mut value = 0.0;
                for (allele, count) in counts.iter().enumerate() {
                    if *count > 0.0 {
                        value += count * float::ln(shares[allele].max(1e-300));
                    }
                }
                for (row, &first) in others.iter().enumerate() {
                    gradient[row] = counts[first] - total * shares[first];
                    for (column, &second) in others.iter().enumerate() {
                        let same = if row == column { shares[first] } else { 0.0 };
                        information[row * block + column] =
                            total * (same - shares[first] * shares[second]);
                    }
                }
                value
            }
            CurveData::Tabulated {
                grid,
                values: tables,
            } => {
                let table = &tables[point];
                if table.is_empty() {
                    return 0.0;
                }
                let spacing = grid[1] - grid[0];
                let last = grid.len() - 1;
                let x = values[0];
                if x <= grid[0] {
                    return table[0];
                }
                if x >= grid[last] {
                    return table[last];
                }
                let node = (((x - grid[0]) / spacing) as usize).min(last - 1);
                let t = (x - grid[node]) / spacing;
                // Each node's slope from its neighbours, one-sided at the ends.
                let slope = |at: usize| {
                    if at == 0 {
                        table[1] - table[0]
                    } else if at == last {
                        table[last] - table[last - 1]
                    } else {
                        0.5 * (table[at + 1] - table[at - 1])
                    }
                };
                let (y0, y1) = (table[node], table[node + 1]);
                let (m0, m1) = (slope(node), slope(node + 1));
                let (t2, t3) = (t * t, t * t * t);
                let value = (2.0 * t3 - 3.0 * t2 + 1.0) * y0
                    + (t3 - 2.0 * t2 + t) * m0
                    + (-2.0 * t3 + 3.0 * t2) * y1
                    + (t3 - t2) * m1;
                let first = (6.0 * t2 - 6.0 * t) * y0
                    + (3.0 * t2 - 4.0 * t + 1.0) * m0
                    + (-6.0 * t2 + 6.0 * t) * y1
                    + (3.0 * t2 - 2.0 * t) * m1;
                let second = (12.0 * t - 6.0) * y0
                    + (6.0 * t - 4.0) * m0
                    + (-12.0 * t + 6.0) * y1
                    + (6.0 * t - 2.0) * m1;
                gradient[0] = first / spacing;
                // Where the table bends the wrong way, a little information keeps the step finite;
                // the step halving keeps the climb honest.
                information[0] = (-second / (spacing * spacing)).max(1e-3);
                value
            }
        }
    }

    fn has_counts(&self, point: usize) -> bool {
        match self {
            CurveData::Binomial { trials, .. } => trials[point] > 0.0,
            CurveData::Multinomial { counts, .. } => counts[point].iter().sum::<f64>() > 0.0,
            CurveData::Tabulated { values, .. } => !values[point].is_empty(),
        }
    }
}

/// **One curve's state carried from round to round**: its unbounded values at every grid point
/// (`block` numbers a point, point-major), its penalty strength, and how much each bend's penalty
/// is eased.
#[derive(Debug, Clone)]
pub(crate) struct CurveState {
    points: usize,
    block: usize,
    values: Vec<f64>,
    ln_smoothing: f64,
    /// Per bend (second difference `i` of component `k`, at `i * block + k`), its weight in the
    /// penalty, in `(0, 1]`.
    bend_weights: Vec<f64>,
}

/// What a fit at one penalty strength reached.
#[derive(Debug, Clone)]
struct CurveFit {
    values: Vec<f64>,
    bend_weights: Vec<f64>,
    /// The Laplace approximation to the log marginal likelihood at this strength, up to a
    /// constant that does not depend on it.
    marginal: f64,
}

/// **The penalty on one bend, in standard deviations `size` of the curve's smoothness**: `size² /
/// 2` up to `knee`, and beyond it growing only as `knee² · ln(size / knee)`. Its derivative over
/// `size` is `size · weight` with the bend's weight `min(1, (knee / size)²)`.
///
/// As a prior on a bend this is a normal distribution with tails falling as a power, `size^−knee²`,
/// so a bend many deviations wide is rare but not ruled out. Its normalising constant does not
/// depend on the strength once the bend is measured in deviations, which is what lets the marginal
/// likelihood compare strengths.
fn bend_penalty(size: f64, knee: f64) -> f64 {
    let size = size.abs();
    if size <= knee {
        0.5 * size * size
    } else {
        0.5 * knee * knee + knee * knee * float::ln(size / knee)
    }
}

impl CurveState {
    pub(crate) fn new(points: usize, block: usize) -> Self {
        Self {
            points,
            block,
            values: vec![0.0; points * block],
            ln_smoothing: 0.0,
            bend_weights: vec![1.0; points.saturating_sub(2) * block],
        }
    }

    fn bends(&self) -> usize {
        self.points.saturating_sub(2)
    }

    /// The second difference of `component` at `bend`.
    fn bend(&self, values: &[f64], bend: usize, component: usize) -> f64 {
        let at = |point: usize| values[point * self.block + component];
        at(bend) - 2.0 * at(bend + 1) + at(bend + 2)
    }

    /// **Refits the curve to `data`** and returns its unbounded values, point-major.
    ///
    /// With [`Smoothing::MarginalLikelihood`] the penalty strength is searched over the whole range
    /// when `search_widely`, and otherwise within two natural-log units of the last round's,
    /// widening when the best lands on an edge.
    pub(crate) fn refit(
        &mut self,
        data: &CurveData,
        config: &PenalisedFitConfig,
        search_widely: bool,
    ) -> Vec<f64> {
        let knee = config.kink_knee;
        let with_counts = (0..self.points).filter(|&p| data.has_counts(p)).count();
        let fitted = match config.smoothing {
            Smoothing::Fixed(strength) => self.fit_robustly_at(
                data,
                float::ln(strength),
                &self.values,
                &self.bend_weights,
                knee,
            ),
            // Two points or fewer: the line through them is the curve, whatever the strength.
            Smoothing::MarginalLikelihood if with_counts <= 2 || self.bends() == 0 => {
                let fitted = self.fit_robustly_at(
                    data,
                    LN_SMOOTHING_RANGE.1,
                    &self.values,
                    &vec![1.0; self.bend_weights.len()],
                    f64::INFINITY,
                );
                self.ln_smoothing = LN_SMOOTHING_RANGE.1;
                fitted
            }
            Smoothing::MarginalLikelihood => {
                let (ln_smoothing, fitted) = self.choose_strength(data, knee, search_widely);
                self.ln_smoothing = ln_smoothing;
                fitted
            }
        };
        if let Smoothing::Fixed(strength) = config.smoothing {
            self.ln_smoothing = float::ln(strength);
        }
        self.values = fitted.values;
        self.bend_weights = fitted.bend_weights;
        self.values.clone()
    }

    /// The middle repeat count of every bend of `component` eased below half its weight.
    fn eased_points(&self, component: usize, first_repeats: u64) -> Vec<u64> {
        (0..self.bends())
            .filter(|bend| self.bend_weights[bend * self.block + component] < 0.5)
            .map(|bend| first_repeats + bend as u64 + 1)
            .collect()
    }

    /// **The penalty strength whose marginal likelihood is highest, and the fit there.**
    ///
    /// A grid from the strongest penalty to the weakest, each fit starting from the last one's
    /// values and eased bends, then golden-section search around the best. **Strong first, so that
    /// a kink is found where the curve is stiff**: at a weak penalty the curve can bend through a
    /// kink without easing anything, and a search started there stays there.
    fn choose_strength(&self, data: &CurveData, knee: f64, search_widely: bool) -> (f64, CurveFit) {
        let (low, high) = if search_widely {
            LN_SMOOTHING_RANGE
        } else {
            (
                (self.ln_smoothing - 2.0).max(LN_SMOOTHING_RANGE.0),
                (self.ln_smoothing + 2.0).min(LN_SMOOTHING_RANGE.1),
            )
        };
        let steps = (((high - low) / 1.15).ceil() as usize).max(2);
        let grid: Vec<f64> = (0..=steps)
            .rev()
            .map(|step| low + (high - low) * step as f64 / steps as f64)
            .collect();
        let mut best: Option<(usize, CurveFit)> = None;
        let mut from = (self.values.clone(), self.bend_weights.clone());
        for (index, &ln_smoothing) in grid.iter().enumerate() {
            let fitted = self.fit_robustly_at(data, ln_smoothing, &from.0, &from.1, knee);
            from = (fitted.values.clone(), fitted.bend_weights.clone());
            if best
                .as_ref()
                .is_none_or(|(_, held)| fitted.marginal > held.marginal)
            {
                best = Some((index, fitted));
            }
        }
        // PANIC-FREE: the grid has at least three strengths.
        let (best_index, best_fit) = best.expect("a strength was scored");
        // On a narrow search, an edge best means the strength moved further than two units: search
        // the whole range instead.
        if !search_widely
            && ((best_index == steps && low > LN_SMOOTHING_RANGE.0)
                || (best_index == 0 && high < LN_SMOOTHING_RANGE.1))
        {
            return self.choose_strength(data, knee, true);
        }
        let (mut left_end, mut right_end) = (
            grid[(best_index + 1).min(steps)],
            grid[best_index.saturating_sub(1)],
        );
        let score = |ln_smoothing: f64| {
            self.fit_robustly_at(
                data,
                ln_smoothing,
                &best_fit.values,
                &best_fit.bend_weights,
                knee,
            )
        };
        let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
        let (mut left, mut right) = (
            right_end - ratio * (right_end - left_end),
            left_end + ratio * (right_end - left_end),
        );
        let (mut left_fit, mut right_fit) = (score(left), score(right));
        for _ in 0..12 {
            if left_fit.marginal > right_fit.marginal {
                right_end = right;
                right = left;
                right_fit = left_fit;
                left = right_end - ratio * (right_end - left_end);
                left_fit = score(left);
            } else {
                left_end = left;
                left = right;
                left_fit = right_fit;
                right = left_end + ratio * (right_end - left_end);
                right_fit = score(right);
            }
        }
        let refined = if left_fit.marginal > right_fit.marginal {
            (left, left_fit)
        } else {
            (right, right_fit)
        };
        if refined.1.marginal >= best_fit.marginal {
            refined
        } else {
            (grid[best_index], best_fit)
        }
    }

    /// **The fit at one strength with every bend's weight settled**: Newton's method at the weights,
    /// then the weights from the bends it reached, until they stop moving (at most eight times).
    /// Each pass lowers the penalised objective's surrogate, so the passes climb the objective.
    fn fit_robustly_at(
        &self,
        data: &CurveData,
        ln_smoothing: f64,
        from: &[f64],
        weights: &[f64],
        knee: f64,
    ) -> CurveFit {
        let strength = float::exp(ln_smoothing);
        let deviation = float::exp(0.5 * ln_smoothing);
        let mut values = from.to_vec();
        let mut weights = if knee.is_finite() {
            weights.to_vec()
        } else {
            vec![1.0; weights.len()]
        };
        let passes = if knee.is_finite() { 8 } else { 1 };
        let mut log_likelihood = 0.0;
        let mut ln_determinant = 0.0;
        for _ in 0..passes {
            (values, log_likelihood, ln_determinant) =
                self.newton(data, strength, &values, &weights);
            if !knee.is_finite() {
                break;
            }
            let mut moved: f64 = 0.0;
            for bend in 0..self.bends() {
                for component in 0..self.block {
                    let size = (self.bend(&values, bend, component) * deviation).abs();
                    let weight = if size > knee {
                        float::powi(knee / size, 2)
                    } else {
                        1.0
                    };
                    let held = &mut weights[bend * self.block + component];
                    moved = moved.max((weight - *held).abs());
                    *held = weight;
                }
            }
            if moved < 1e-3 {
                break;
            }
        }
        let mut prior = 0.0;
        for bend in 0..self.bends() {
            for component in 0..self.block {
                prior -= bend_penalty(self.bend(&values, bend, component) * deviation, knee);
            }
        }
        let rank = (self.bends() * self.block) as f64;
        let ridge = 0.5 * RIDGE * values.iter().map(|v| v * v).sum::<f64>();
        let marginal =
            log_likelihood + prior - ridge + 0.5 * rank * ln_smoothing - 0.5 * ln_determinant;
        CurveFit {
            values,
            bend_weights: weights,
            marginal,
        }
    }

    /// The quadratic penalty `½ λ Σ weight · bend²` at `values`, without `λ`.
    fn roughness(&self, values: &[f64], weights: &[f64]) -> f64 {
        let mut total = 0.0;
        for bend in 0..self.bends() {
            for component in 0..self.block {
                let size = self.bend(values, bend, component);
                total += weights[bend * self.block + component] * size * size;
            }
        }
        0.5 * total
    }

    /// **The penalised fit at one strength and one set of bend weights**, by Newton's method with
    /// step halving from `from`. Returns the values, their log-likelihood, and ln |information +
    /// penalty + ridge| there.
    fn newton(
        &self,
        data: &CurveData,
        strength: f64,
        from: &[f64],
        weights: &[f64],
    ) -> (Vec<f64>, f64, f64) {
        let size = self.points * self.block;
        let mut values = from.to_vec();
        let mut gradient = vec![0.0; self.block];
        let mut information = vec![0.0; self.block * self.block];
        let log_likelihood_at = |values: &[f64], gradient: &mut [f64], information: &mut [f64]| {
            (0..self.points)
                .map(|point| {
                    data.evaluate(
                        point,
                        &values[point * self.block..(point + 1) * self.block],
                        gradient,
                        information,
                    )
                })
                .sum::<f64>()
        };
        let objective = |log_likelihood: f64, values: &[f64]| {
            log_likelihood
                - strength * self.roughness(values, weights)
                - 0.5 * RIDGE * values.iter().map(|v| v * v).sum::<f64>()
        };
        let mut log_likelihood = log_likelihood_at(&values, &mut gradient, &mut information);
        let mut current = objective(log_likelihood, &values);
        let mut ln_determinant = 0.0;
        for iteration in 0..100 {
            // The system: (information + λ·penalty + ridge) · step = gradient of the objective.
            let mut matrix = self.penalty_matrix(strength, weights);
            let mut right_hand = vec![0.0; size];
            for point in 0..self.points {
                let block = point * self.block..(point + 1) * self.block;
                data.evaluate(
                    point,
                    &values[block.clone()],
                    &mut gradient,
                    &mut information,
                );
                for row in 0..self.block {
                    right_hand[point * self.block + row] = gradient[row];
                    for column in 0..=row {
                        matrix.add(
                            point * self.block + row,
                            point * self.block + column,
                            information[row * self.block + column],
                        );
                    }
                }
            }
            for (index, value) in values.iter().enumerate() {
                right_hand[index] -= RIDGE * value;
            }
            let penalty_gradient = self.penalty_times(&values, weights);
            for (into, penalty) in right_hand.iter_mut().zip(&penalty_gradient) {
                *into -= strength * penalty;
            }
            let Some(factor) = matrix.cholesky() else {
                break;
            };
            ln_determinant = factor.ln_determinant();
            if iteration > 0 && right_hand.iter().all(|g| g.abs() < 1e-9) {
                break;
            }
            let mut step = right_hand;
            factor.solve(&mut step);
            let mut length = 1.0;
            let mut accepted = false;
            for _ in 0..40 {
                let trial: Vec<f64> = values
                    .iter()
                    .zip(&step)
                    .map(|(v, s)| (v + length * s).clamp(-LARGEST_UNBOUNDED, LARGEST_UNBOUNDED))
                    .collect();
                let trial_likelihood = log_likelihood_at(&trial, &mut gradient, &mut information);
                let score = objective(trial_likelihood, &trial);
                if score >= current - 1e-12 * current.abs() {
                    let moved = trial
                        .iter()
                        .zip(&values)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0, f64::max);
                    values = trial;
                    current = score;
                    log_likelihood = trial_likelihood;
                    accepted = moved > 1e-10;
                    break;
                }
                length *= 0.5;
            }
            if !accepted {
                break;
            }
        }
        (values, log_likelihood, ln_determinant)
    }

    /// `λ · (Dᵀ W D)` per component, plus the ridge, as a banded matrix: `D` takes second
    /// differences along the points, `W` the bends' weights.
    fn penalty_matrix(&self, strength: f64, weights: &[f64]) -> SymmetricBand {
        let size = self.points * self.block;
        let mut matrix = SymmetricBand::zeros(size, 2 * self.block);
        for index in 0..size {
            matrix.add(index, index, RIDGE);
        }
        let stencil = [1.0, -2.0, 1.0];
        for bend in 0..self.bends() {
            for component in 0..self.block {
                let weight = strength * weights[bend * self.block + component];
                for (i, a) in stencil.iter().enumerate() {
                    for (j, b) in stencil.iter().enumerate().take(i + 1) {
                        matrix.add(
                            (bend + i) * self.block + component,
                            (bend + j) * self.block + component,
                            weight * a * b,
                        );
                    }
                }
            }
        }
        matrix
    }

    /// `Dᵀ W D · values`, without `λ`.
    fn penalty_times(&self, values: &[f64], weights: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; values.len()];
        let stencil = [1.0, -2.0, 1.0];
        for bend in 0..self.bends() {
            for component in 0..self.block {
                let weighted =
                    weights[bend * self.block + component] * self.bend(values, bend, component);
                for (offset, coefficient) in stencil.iter().enumerate() {
                    out[(bend + offset) * self.block + component] += coefficient * weighted;
                }
            }
        }
        out
    }
}

// ---------------------------------------------------------------------
// A symmetric banded matrix and its Cholesky factor
// ---------------------------------------------------------------------

/// A symmetric matrix kept as its lower band: entry `(i, i − d)` for `d ≤ width`.
struct SymmetricBand {
    size: usize,
    width: usize,
    lower: Vec<f64>,
}

impl SymmetricBand {
    fn zeros(size: usize, width: usize) -> Self {
        Self {
            size,
            width,
            lower: vec![0.0; size * (width + 1)],
        }
    }

    /// Adds `value` at `(row, column)` (and so at `(column, row)`).
    fn add(&mut self, row: usize, column: usize, value: f64) {
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
    fn cholesky(mut self) -> Option<BandCholesky> {
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

struct BandCholesky {
    size: usize,
    width: usize,
    lower: Vec<f64>,
}

impl BandCholesky {
    fn at(&self, i: usize, j: usize) -> f64 {
        self.lower[i * (self.width + 1) + (i - j)]
    }

    fn ln_determinant(&self) -> f64 {
        2.0 * (0..self.size)
            .map(|i| float::ln(self.at(i, i)))
            .sum::<f64>()
    }

    /// Overwrites `rhs` with the solution of `L Lᵀ x = rhs`.
    fn solve(&self, rhs: &mut [f64]) {
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
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};
    use crate::parameter_estimation::joint::ssr_fit::composite_fit::fit_class_by_counting;

    /// **The banded Cholesky solves what a dense solve does**, on a random positive-definite
    /// matrix of bandwidth four.
    #[test]
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

    /// Binomial successes out of `trials` at each probability, with sampling noise from a fixed
    /// stream (the normal approximation, which is close at these counts).
    fn noisy_successes(probabilities: &[f64], trials: &[f64], seed: u64) -> Vec<f64> {
        let mut state = seed;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
        };
        probabilities
            .iter()
            .zip(trials)
            .map(|(p, n)| {
                let normal = (-2.0 * float::ln(uniform())).sqrt()
                    * float::cos(std::f64::consts::TAU * uniform());
                (n * p + normal * (n * p * (1.0 - p)).sqrt()).clamp(0.0, *n)
            })
            .collect()
    }

    /// **A well-measured kink survives the penalty, and a thin point is pulled onto the curve.**
    /// Binomial counts along a slip rate that rises and bends smoothly from 8 to 20 repeats, with
    /// the point at 9 repeats drawn at 0.53 of the curve on 400,000 reads (about 3,000 slipped, near
    /// HG002's 9-repeat homopolymers), every other point on 3,000 reads, one at 15 repeats on 40
    /// reads that came out with no slipped read, and sampling noise everywhere. Five draws.
    #[test]
    fn a_well_measured_kink_survives_and_a_thin_point_borrows() {
        let repeats: Vec<u64> = (8..=20).collect();
        let truth: Vec<f64> = repeats
            .iter()
            .map(|n| {
                let x = *n as f64 - 8.0;
                logistic(-5.3 + 0.45 * x - 0.012 * x * x)
            })
            .collect();
        let mut drawn = truth.clone();
        drawn[1] *= 0.53;
        let trials: Vec<f64> = repeats
            .iter()
            .map(|n| match n {
                9 => 400_000.0,
                15 => 40.0,
                _ => 3_000.0,
            })
            .collect();
        for seed in 1..=5 {
            let mut successes = noisy_successes(&drawn, &trials, seed);
            successes[7] = 0.0;
            let data = CurveData::Binomial {
                successes,
                trials: trials.clone(),
            };
            let config = PenalisedFitConfig::default();
            let mut state = CurveState::new(repeats.len(), 1);
            let fitted = state.refit(&data, &config, true);
            let at_nine = logistic(fitted[1]);
            let at_fifteen = logistic(fitted[7]);

            assert!(
                (at_nine / drawn[1] - 1.0).abs() < 0.05,
                "the kink at 9 repeats: {at_nine} against {}",
                drawn[1]
            );
            assert!(
                (at_fifteen / truth[7] - 1.0).abs() < 0.25,
                "the thin point: {at_fifteen} against {}",
                truth[7]
            );
        }
    }

    /// **With the penalty near zero the fit is the counting fit, class by class**: three classes of
    /// one motif length, each fitted on its own by counting and all together with a penalty of
    /// 10⁻⁸.
    #[test]
    fn with_no_penalty_the_fit_is_the_counting_fit() {
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let spectrum = spectrum_of(3);
        let classes: Vec<StratumEvidence> = [(8_u64, 0.03), (9, 0.05), (10, 0.08)]
            .iter()
            .enumerate()
            .map(|(index, (repeats, level))| {
                let truth = Slippage {
                    level: *level,
                    shorter_share: 0.8,
                    fall_off: 0.25,
                };
                let mut evidence =
                    draw_stratum(truth, &spectrum, 0.5, 0.4, 300, 10, 8, 1, 31 + index as u64);
                evidence.stratum = Stratum {
                    period: 1,
                    reference_repeats: *repeats,
                };
                evidence
            })
            .collect();
        let excess = [0.4; 10];
        let fit = PenalisedFitConfig {
            smoothing: Smoothing::Fixed(1e-8),
            kink_knee: f64::INFINITY,
            tolerance: 1e-10,
            max_rounds: 2000,
            ..PenalisedFitConfig::default()
        };
        let together = fit_period_by_penalised_counting(&classes, &excess, &config, &fit);
        for (evidence, fitted) in classes.iter().zip(&together.classes) {
            let fitted = fitted.as_ref().expect("every class has reads");
            let alone = fit_class_by_counting(evidence, fit.start, &excess, &config, 2000, 1e-10);
            assert!(
                (fitted.slippage.level - alone.slippage.level).abs() < 1e-5,
                "{fitted:?} against {alone:?}"
            );
            assert!(
                (fitted.slippage.shorter_share - alone.slippage.shorter_share).abs() < 1e-4,
                "{fitted:?} against {alone:?}"
            );
            assert!(
                (fitted.slippage.fall_off - alone.slippage.fall_off).abs() < 1e-4,
                "{fitted:?} against {alone:?}"
            );
            for (a, b) in fitted
                .allele_length_variation
                .iter()
                .zip(&alone.allele_length_variation)
            {
                assert!((a - b).abs() < 1e-4, "{fitted:?} against {alone:?}");
            }
            assert!(
                (float::ln(fitted.concentration) - float::ln(alone.concentration)).abs() < 1e-3,
                "{fitted:?} against {alone:?}"
            );
        }
    }
}
