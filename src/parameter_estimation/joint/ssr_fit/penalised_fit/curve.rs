//! **One curve's penalised fit**: the counts a curve is fitted to ([`CurveData`]), the prior on
//! how sharply it bends ([`BendPrior`]), and Newton's method on the penalised objective at each
//! penalty strength, with the search for the strength the reads support ([`CurveState::refit`]).
//!
//! Each round of the repeat-tract fit hands every curve its counts here and reads back its values;
//! nothing here knows about classes, loci or rounds.

use super::*;
use crate::parameter_estimation::joint::linear_algebra::{BandCholesky, SymmetricBand};

/// **The slip rate's curve data with read groups apart** (spec `str_library_stutter.md` §3.3): at
/// each repeat count with counts, each read group's slipped and exactly-read reads at its
/// multiplier, and the stutter prior, `prior_reads` read at the rate itself with `centre` of it
/// slipped. A repeat count `from_the_curve` holds no data and takes the curve's value (spec
/// `str_penalised_fit.md` §2.11).
pub(super) fn multiplied_binomial(
    counts: &[Option<ClassCounts>],
    from_the_curve: &[bool],
    multipliers: &[f64],
    prior_reads: f64,
    centre: f64,
) -> CurveData {
    let points = counts.len();
    let (mut read_groups, mut successes, mut trials) = (
        vec![Vec::new(); points],
        vec![0.0; points],
        vec![0.0; points],
    );
    for (point, count) in counts.iter().enumerate() {
        let Some(count) = count else {
            continue;
        };
        if from_the_curve[point] {
            continue;
        }
        successes[point] = prior_reads * centre;
        trials[point] = prior_reads;
        read_groups[point] = count
            .read_groups
            .iter()
            .zip(multipliers)
            .filter(|(reads, _)| reads.exact + reads.slipped > 0.0)
            .map(|(reads, multiplier)| MultipliedReads {
                multiplier: *multiplier,
                slipped: reads.slipped.max(0.0),
                exact: reads.exact.max(0.0),
            })
            .collect();
    }
    CurveData::MultipliedBinomial {
        read_groups,
        successes,
        trials,
    }
}

/// **The counting objective's slope in the slip rate's logit `eta`, read groups apart**: each read
/// group's slope at the class's rate times its multiplier, summed ([`multiplied_at_logit`]). At
/// every multiplier one it is the binomial's, `slipped − p · (slipped + exact)`.
pub(super) fn multiplied_level_slope(counts: &ClassCounts, multipliers: &[f64], eta: f64) -> f64 {
    counts
        .read_groups
        .iter()
        .zip(multipliers)
        .map(|(reads, multiplier)| {
            multiplied_at_logit(
                MultipliedReads {
                    multiplier: *multiplier,
                    slipped: reads.slipped,
                    exact: reads.exact,
                },
                eta,
            )
            .1
        })
        .sum()
}

/// **The counts one curve is fitted to**, per grid point (a repeat count). A point whose counts are
/// all zero carries no evidence and takes the curve's value.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum CurveData {
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
    /// **A slip rate whose reads come from read groups at different multipliers** (spec
    /// `str_library_stutter.md` §3.3): at each point, each read group's multiplier with its slipped
    /// and exactly-read reads, a read group's reads slipping at the point's rate times its
    /// multiplier; and `successes` out of `trials` read at the rate itself, which carry the stutter
    /// prior. The unbounded number is the logit of the rate.
    MultipliedBinomial {
        read_groups: Vec<Vec<MultipliedReads>>,
        successes: Vec<f64>,
        trials: Vec<f64>,
    },
    /// A log-likelihood tabulated at each point over one unbounded number, on `grid` (ascending,
    /// evenly spaced); an empty table for a point with none. Read between the nodes by a cubic
    /// with the table's slopes, and held at the end node beyond them.
    Tabulated {
        grid: Vec<f64>,
        values: Vec<Vec<f64>>,
    },
}

/// **One read group's reads at one point of a [`CurveData::MultipliedBinomial`]**: its multiplier,
/// and its expected slipped and exactly-read reads there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MultipliedReads {
    pub multiplier: f64,
    pub slipped: f64,
    pub exact: f64,
}

/// **The log-likelihood of `slipped` out of `slipped + exact` reads at the logit `eta`**, with its
/// slope and information in `eta`: `y·η − n·ln(1 + e^η)`, computed stably.
fn binomial_at_logit(slipped: f64, trials: f64, eta: f64) -> (f64, f64, f64) {
    let p = logistic(eta);
    let ln_one_plus = if eta > 0.0 {
        eta + float::ln_1p(float::exp(-eta))
    } else {
        float::ln_1p(float::exp(eta))
    };
    (
        slipped * eta - trials * ln_one_plus,
        slipped - trials * p,
        trials * p * (1.0 - p),
    )
}

/// **One read group's reads at a rate `p = logistic(eta)` times `multiplier`**, with the slope and
/// information in `eta`: `S·ln(m p) + N·ln(1 − m p)`. Where `m p` reaches the ceiling the read
/// group's rate is held there ([`read_group_level`]), so it adds its value and no slope.
///
/// With `λ = m p`: the slope is `S(1 − p) − N m p (1 − p) / (1 − λ)`, and the information
/// `S p(1 − p) + N m p (1 − p) (1 − 2p + m p²) / (1 − λ)²`, which at `m = 1` is the binomial's
/// `n p (1 − p)`.
fn multiplied_at_logit(reads: MultipliedReads, eta: f64) -> (f64, f64, f64) {
    let MultipliedReads {
        multiplier,
        slipped,
        exact,
    } = reads;
    let p = logistic(eta);
    let level = read_group_level(p, multiplier);
    if level < multiplier * p {
        // Held at the ceiling.
        return (
            slipped * float::ln(level) + exact * float::ln_1p(-level),
            0.0,
            0.0,
        );
    }
    let ln_p = if eta > 0.0 {
        -float::ln_1p(float::exp(-eta))
    } else {
        eta - float::ln_1p(float::exp(eta))
    };
    let value = slipped * (float::ln(multiplier) + ln_p) + exact * float::ln_1p(-level);
    let spread = p * (1.0 - p);
    let slope = slipped * (1.0 - p) - exact * multiplier * spread / (1.0 - level);
    let information = slipped * spread
        + exact * multiplier * spread * (1.0 - 2.0 * p + multiplier * p * p)
            / ((1.0 - level) * (1.0 - level));
    (value, slope, information)
}

impl CurveData {
    /// **One point's log-likelihood at `values`**, with its gradient and its information (minus
    /// the second derivative), the information `block × block` row by row. Zero, with no
    /// information, at a point with no counts.
    pub(super) fn evaluate(
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
                let (value, slope, info) = binomial_at_logit(y, n, values[0]);
                gradient[0] = slope;
                information[0] = info;
                value
            }
            CurveData::MultipliedBinomial {
                read_groups,
                successes,
                trials,
            } => {
                let eta = values[0];
                let (mut value, mut slope, mut info) = if trials[point] > 0.0 {
                    binomial_at_logit(successes[point], trials[point], eta)
                } else {
                    (0.0, 0.0, 0.0)
                };
                for reads in &read_groups[point] {
                    let (v, g, i) = multiplied_at_logit(*reads, eta);
                    value += v;
                    slope += g;
                    info += i;
                }
                gradient[0] = slope;
                // Where read groups above one near the ceiling bend it the wrong way, a little
                // information keeps the step finite, as for a table; the step halving keeps the
                // climb honest.
                information[0] = if read_groups[point].is_empty() {
                    info
                } else {
                    info.max(1e-3)
                };
                value
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
            CurveData::MultipliedBinomial {
                read_groups,
                trials,
                ..
            } => {
                trials[point] > 0.0
                    || read_groups[point]
                        .iter()
                        .any(|reads| reads.slipped + reads.exact > 0.0)
            }
            CurveData::Multinomial { counts, .. } => counts[point].iter().sum::<f64>() > 0.0,
            CurveData::Tabulated { values, .. } => !values[point].is_empty(),
        }
    }
}

/// **One curve's state carried from round to round**: its unbounded values at every grid point
/// (`block` numbers a point, point-major), its penalty strength, and how much each bend's penalty
/// is eased.
#[derive(Debug, Clone)]
pub(super) struct CurveState {
    points: usize,
    block: usize,
    pub(super) values: Vec<f64>,
    pub(super) ln_smoothing: f64,
    /// Per bend (second difference `i` of component `k`, at `i * block + k`), its weight in the
    /// penalty, in `(0, 1]`.
    pub(super) bend_weights: Vec<f64>,
    /// Whether the strength has stopped moving ([`PenalisedFitConfig::strength_steady_within`]),
    /// after which a search near it keeps it.
    steady: bool,
    /// How many strengths tried in searches so far could not be scored, because the fit there
    /// could not be factorised even with a larger ridge, and so were never chosen.
    pub(super) strengths_rejected: usize,
    /// How many searches so far found no strength that could be scored, and so kept the strength
    /// held before.
    pub(super) searches_without_a_strength: usize,
    /// How many fits at a strength no search chose could not be factorised, so the curve kept the
    /// values Newton's method had reached when it stopped.
    pub(super) unfactorised_fits_kept: usize,
}

/// What a search for a curve's penalty strength found.
pub(super) struct StrengthChoice {
    /// The strength chosen, as a natural logarithm, and the fit there; `None` when not one
    /// strength tried could be scored.
    pub(super) chosen: Option<(f64, CurveFit)>,
    /// How many strengths tried could not be scored, and so were never chosen.
    pub(super) rejected: usize,
}

/// What a fit at one penalty strength reached.
#[derive(Debug, Clone)]
pub(super) struct CurveFit {
    values: Vec<f64>,
    bend_weights: Vec<f64>,
    /// The Laplace approximation to the log marginal likelihood at this strength, up to a
    /// constant that does not depend on it.
    pub(super) marginal: f64,
}

/// **What a bend costs, in standard deviations `size` of the curve's smoothness**: `size² / 2` up to
/// `knee`, and beyond it growing only as `tail · ln(size / knee)`. Its derivative over `size` is
/// `size · weight`, the bend's weight one up to the knee and `tail / size²` beyond.
///
/// As a prior on a bend this is a normal distribution with tails falling as a power, `size^−tail`,
/// so a bend many deviations wide is rare but not ruled out; `tail` must exceed one for the prior to
/// be a distribution. Its normalising constant does not depend on the strength once the bend is
/// measured in deviations, which is what lets the marginal likelihood compare strengths.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct BendPrior {
    knee: f64,
    tail: f64,
}

impl BendPrior {
    /// Every bend costs its square: no knee.
    pub(super) const QUADRATIC: Self = Self {
        knee: f64::INFINITY,
        tail: f64::INFINITY,
    };

    fn eases(self) -> bool {
        self.knee.is_finite()
    }

    fn penalty(self, size: f64) -> f64 {
        let size = size.abs();
        if size <= self.knee {
            0.5 * size * size
        } else {
            0.5 * self.knee * self.knee + self.tail * float::ln(size / self.knee)
        }
    }

    fn weight(self, size: f64) -> f64 {
        if size > self.knee {
            self.tail / (size * size)
        } else {
            1.0
        }
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
            steady: false,
            strengths_rejected: 0,
            searches_without_a_strength: 0,
            unfactorised_fits_kept: 0,
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
    /// With [`Smoothing::MarginalLikelihood`] the penalty strength is chosen as `search` says.
    pub(crate) fn refit(
        &mut self,
        data: &CurveData,
        config: &PenalisedFitConfig,
        search: StrengthSearch,
    ) -> Vec<f64> {
        let knee = BendPrior {
            knee: config.kink_knee,
            tail: config.kink_tail,
        };
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
                    BendPrior::QUADRATIC,
                );
                self.ln_smoothing = LN_SMOOTHING_RANGE.1;
                fitted
            }
            Smoothing::MarginalLikelihood
                if search == StrengthSearch::Held
                    || (search == StrengthSearch::Near && self.steady) =>
            {
                self.fit_robustly_at(
                    data,
                    self.ln_smoothing,
                    &self.values,
                    &self.bend_weights,
                    knee,
                )
            }
            Smoothing::MarginalLikelihood => {
                let choice = self.choose_strength(data, knee, search == StrengthSearch::Wide);
                self.strengths_rejected += choice.rejected;
                if let Some((ln_smoothing, fitted)) = choice.chosen {
                    if search == StrengthSearch::Near
                        && config.strength_steady_within > 0.0
                        && (ln_smoothing - self.ln_smoothing).abs() < config.strength_steady_within
                    {
                        self.steady = true;
                    }
                    self.ln_smoothing = ln_smoothing;
                    fitted
                } else {
                    // No strength could be scored: the strength held is kept, and the curve is
                    // fitted at it as a held strength would be.
                    self.searches_without_a_strength += 1;
                    self.fit_robustly_at(
                        data,
                        self.ln_smoothing,
                        &self.values,
                        &self.bend_weights,
                        knee,
                    )
                }
            }
        };
        // A search never chooses a fit that could not be factorised, so this counts only fits at
        // a strength no search chose. Their values are kept all the same: there is no better fit
        // at that strength to put in their place.
        if fitted.marginal == f64::NEG_INFINITY {
            self.unfactorised_fits_kept += 1;
        }
        if let Smoothing::Fixed(strength) = config.smoothing {
            self.ln_smoothing = float::ln(strength);
        }
        self.values = fitted.values;
        self.bend_weights = fitted.bend_weights;
        self.values.clone()
    }

    /// The middle repeat count of every bend of `component` eased below half its weight.
    pub(super) fn eased_points(&self, component: usize, first_repeats: u64) -> Vec<u64> {
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
    ///
    /// **A strength whose fit could not be factorised is never chosen** (its marginal likelihood
    /// is minus infinity, [`Self::fit_robustly_at`]); each one is counted in
    /// [`StrengthChoice::rejected`]. When no strength tried can be scored, a narrow search
    /// searches the whole range instead, and nothing is chosen only when that scores none either.
    pub(super) fn choose_strength(
        &self,
        data: &CurveData,
        knee: BendPrior,
        search_widely: bool,
    ) -> StrengthChoice {
        let rejected = std::cell::Cell::new(0_usize);
        let fit_at = |ln_smoothing: f64, from: &[f64], weights: &[f64]| {
            let fitted = self.fit_robustly_at(data, ln_smoothing, from, weights, knee);
            if fitted.marginal == f64::NEG_INFINITY {
                rejected.set(rejected.get() + 1);
            }
            fitted
        };
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
            let fitted = fit_at(ln_smoothing, &from.0, &from.1);
            from = (fitted.values.clone(), fitted.bend_weights.clone());
            if best
                .as_ref()
                .is_none_or(|(_, held)| fitted.marginal > held.marginal)
            {
                best = Some((index, fitted));
            }
        }
        // PANIC-FREE: the grid has at least three strengths.
        let (best_index, best_fit) = best.expect("a strength was tried");
        if best_fit.marginal == f64::NEG_INFINITY {
            // Not one strength on the grid could be scored. Searching between two of them would
            // compare minus infinity with minus infinity. A narrow search tries the whole range
            // first, since the strengths it skipped may score; only when they do not either is
            // nothing chosen.
            if !search_widely && (low > LN_SMOOTHING_RANGE.0 || high < LN_SMOOTHING_RANGE.1) {
                let mut wider = self.choose_strength(data, knee, true);
                wider.rejected += rejected.get();
                return wider;
            }
            return StrengthChoice {
                chosen: None,
                rejected: rejected.get(),
            };
        }
        // On a narrow search, an edge best means the strength moved further than two units: search
        // the whole range instead.
        if !search_widely
            && ((best_index == steps && low > LN_SMOOTHING_RANGE.0)
                || (best_index == 0 && high < LN_SMOOTHING_RANGE.1))
        {
            let mut wider = self.choose_strength(data, knee, true);
            wider.rejected += rejected.get();
            return wider;
        }
        let (mut left_end, mut right_end) = (
            grid[(best_index + 1).min(steps)],
            grid[best_index.saturating_sub(1)],
        );
        let score =
            |ln_smoothing: f64| fit_at(ln_smoothing, &best_fit.values, &best_fit.bend_weights);
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
        let chosen = if refined.1.marginal >= best_fit.marginal {
            refined
        } else {
            (grid[best_index], best_fit)
        };
        StrengthChoice {
            chosen: Some(chosen),
            rejected: rejected.get(),
        }
    }

    /// **The fit at one strength with every bend's weight settled**: Newton's method at the weights,
    /// then the weights from the bends it reached, until they stop moving (at most eight times).
    /// Each pass lowers the penalised objective's surrogate, so the passes climb the objective.
    ///
    /// Its marginal likelihood is minus infinity when the last pass's Newton system could not be
    /// factorised ([`Self::newton`]).
    pub(super) fn fit_robustly_at(
        &self,
        data: &CurveData,
        ln_smoothing: f64,
        from: &[f64],
        weights: &[f64],
        knee: BendPrior,
    ) -> CurveFit {
        let strength = float::exp(ln_smoothing);
        let deviation = float::exp(0.5 * ln_smoothing);
        let mut values = from.to_vec();
        let mut weights = if knee.eases() {
            weights.to_vec()
        } else {
            vec![1.0; weights.len()]
        };
        let passes = if knee.eases() { 8 } else { 1 };
        let mut log_likelihood = 0.0;
        let mut ln_determinant = None;
        for _ in 0..passes {
            (values, log_likelihood, ln_determinant) =
                self.newton(data, strength, &values, &weights);
            if !knee.eases() {
                break;
            }
            let mut moved: f64 = 0.0;
            for bend in 0..self.bends() {
                for component in 0..self.block {
                    let size = (self.bend(&values, bend, component) * deviation).abs();
                    let weight = knee.weight(size);
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
                prior -= knee.penalty(self.bend(&values, bend, component) * deviation);
            }
        }
        let rank = (self.bends() * self.block) as f64;
        let ridge = 0.5 * RIDGE * values.iter().map(|v| v * v).sum::<f64>();
        // A fit whose system could not be factorised has no determinant, so no marginal
        // likelihood: it is scored below every strength that has one and is never chosen.
        let marginal = ln_determinant.map_or(f64::NEG_INFINITY, |ln_determinant| {
            log_likelihood + prior - ridge + 0.5 * rank * ln_smoothing - 0.5 * ln_determinant
        });
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
    ///
    /// **The log-determinant is `None` when a step's system could not be factorised**, even with
    /// the larger ridges of [`factorise_with_retries`]; the values are then those of the last step
    /// taken, `from` if none was. A strength whose fit has no log-determinant cannot be scored,
    /// and [`Self::fit_robustly_at`] scores it as never to be chosen.
    fn newton(
        &self,
        data: &CurveData,
        strength: f64,
        from: &[f64],
        weights: &[f64],
    ) -> (Vec<f64>, f64, Option<f64>) {
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
        let mut ln_determinant = None;
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
            let factor = factorise_with_retries(&matrix);
            // A test can make the factorisation fail at chosen strengths, which real counts
            // reach only through rounding the retry ridges almost always rescue.
            #[cfg(test)]
            let factor = factor.filter(|_| !tests::factorisation_fails_at(strength));
            let Some(factor) = factor else {
                ln_determinant = None;
                break;
            };
            ln_determinant = Some(factor.ln_determinant());
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
        // The band stores one triangle: the entries on or below the diagonal.
        self.for_each_bend_entry(strength, weights, |row, column, component, value| {
            if column <= row {
                matrix.add(
                    row * self.block + component,
                    column * self.block + component,
                    value,
                );
            }
        });
        matrix
    }

    /// **Every entry of `λ · (Dᵀ W D)`, the penalty on the curve's bends**, handed to
    /// `add(row point, column point, component, value)`: for each bend — the second difference
    /// `1, −2, 1` over three neighbouring points — and each component, the stencil times itself,
    /// times the strength `λ` and the bend's weight. Bend by bend, component by component, row by
    /// row, so a caller adding into a matrix adds in one fixed order.
    pub(super) fn for_each_bend_entry(
        &self,
        strength: f64,
        weights: &[f64],
        mut add: impl FnMut(usize, usize, usize, f64),
    ) {
        let stencil = [1.0, -2.0, 1.0];
        for bend in 0..self.bends() {
            for component in 0..self.block {
                let weight = strength * weights[bend * self.block + component];
                for (i, a) in stencil.iter().enumerate() {
                    for (j, b) in stencil.iter().enumerate() {
                        add(bend + i, bend + j, component, weight * a * b);
                    }
                }
            }
        }
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

/// **The ridges a Newton system is retried with when it cannot be factorised at the fit's own
/// [`RIDGE`]**, each ten times the last.
///
/// At a strength of 10¹⁰ the rounding in the band factorisation's last pivots is about
/// 10¹⁰ × 2.2×10⁻¹⁶ × the stencil's diagonal of 6, about 1.3×10⁻⁵, so a pivot that the reads
/// barely inform can come out zero or negative at a ridge of 10⁻⁶. A ridge in the system alone only shortens
/// the step: the objective climbed is unchanged, and so is the point Newton's method converges
/// to. It does enter the log-determinant, by at most ln 1000 for each number no read informs,
/// which lowers that strength's score; a retried strength is therefore never favoured by its
/// retry.
const NEWTON_RETRY_RIDGES: [f64; 3] = [1e-5, 1e-4, 1e-3];

/// **The Cholesky factor of a Newton system**, at the ridge it was built with and then, if that
/// fails, at each of [`NEWTON_RETRY_RIDGES`] in turn; `None` when every one fails.
fn factorise_with_retries(matrix: &SymmetricBand) -> Option<BandCholesky> {
    if let Some(factor) = matrix.clone().cholesky() {
        return Some(factor);
    }
    NEWTON_RETRY_RIDGES.iter().find_map(|&ridge| {
        let mut retried = matrix.clone();
        for index in 0..retried.size {
            retried.add(index, index, ridge - RIDGE);
        }
        retried.cholesky()
    })
}
