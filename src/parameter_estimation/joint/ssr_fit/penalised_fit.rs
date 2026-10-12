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
//! approximation that restricted maximum likelihood rests on, Wood 2011), re-chosen for the first
//! rounds and then held.
//!
//! **A strong, well-measured kink survives.** The penalty is quadratic in each bend up to
//! [`PenalisedFitConfig::kink_knee`] standard deviations of what the curve's smoothness predicts,
//! and grows only logarithmically beyond it ([`PenalisedFitConfig::kink_tail`]), so a class whose own
//! reads insist on a step pays a bounded price for it. The knee is the slippage curve's blend's
//! (`str_slippage_level_curve.md` §7.2).

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
    /// **How much a bend beyond the knee costs per natural-log unit of its size: 2**, so that a bend's
    /// prior falls as `size^−2` past the knee. It must exceed one. Measured on drawn motif lengths
    /// (track 2's report, §3.3): at the knee's square, 6.25, which keeps the penalty's slope
    /// continuous, a kink at 9 repeats was cheaper to absorb by smoothing the whole slip-rate curve
    /// 60 to 100 times less, and its thin classes lost what they had borrowed; at 2 the kink is eased
    /// and they keep it, and where no bend passes the knee the two give the same answer.
    pub kink_tail: f64,
    /// **Allele draws added to every allele length of every class before its shares are chosen:
    /// one draw a class, spread evenly over the lengths** (1/19 of a draw each at ±9), after the
    /// draws are weighed as clustered by locus. Without it a length the reads stop supporting in one
    /// round gets a share near zero, so the next round's E-step gives it no draws and it can never
    /// return: EM has an absorbing boundary at zero, which costs a held-out locus carrying that
    /// length hundreds of log-likelihood units.
    ///
    /// **Why not half a draw a length (the Jeffreys prior).** At one sample and 100 loci a class
    /// holds about 120 draws once clustering is allowed for, and 19 half-draws spread evenly moved
    /// enough mass into the tails to bias the reference length's share by −0.13 (one-base runs) and
    /// −0.20 (two-base runs) on the fitting bench's tomato-like truth.
    pub allele_pseudo_draws: f64,
    /// **Whether each class's allele draws are divided by the design effect of draws clustered by
    /// locus** before its allele-length curve is fitted (`Period::allele_draw_inflation`). **Off by
    /// default**, because measured both ways on drawn motif lengths (track 2's report, §3.5) it helps
    /// the allele length variation and can wreck the stutter: at 50 samples × 30 reads it halved
    /// the thin classes' allele-length error (0.090 to 0.050) with the slip rate unchanged, but at
    /// 2,000 samples × 3 reads the honest information left a straight allele-length curve, slightly
    /// wrong at rich classes (0.018 against 0.009), and with millions of reads the stutter absorbed
    /// the misfit: the slip rate at rich classes came out 12% wrong against 1.7% without it. The
    /// strength is chosen from the allele counts alone, which cannot see what the reads say about
    /// the spread through the stutter.
    pub weigh_allele_draws_by_locus: bool,
    /// **How many times the between-locus share is read off pairs of samples at a locus**, each time
    /// followed by penalised counting with it held. Zero reads it only from each sample's own two
    /// copies, which says almost nothing in an inbred cohort.
    pub pair_alternations: usize,
    /// **Each curve's penalty strength is chosen afresh every round until the fit has nearly
    /// settled** — no number moving by more than `settled_change_for_strength` (10⁻⁵) — **or for
    /// at most `rounds_choosing_strength` rounds (200)**, then held, so that the remaining rounds
    /// are penalised counting at fixed strengths, which climbs its objective every round.
    ///
    /// Both limits are needed. Choosing afresh every round can leave a strength alternating
    /// between two values and the fit never settling: one draw in five at 50 samples × 30 reads ran
    /// 500 rounds that way. Holding after a fixed 30 rounds freezes strengths chosen on counts still
    /// shaped by the flat starting shares: at 2,000 samples × 3 reads, where the fit took 303
    /// rounds, the allele-length curve was held straight and the slip rate at rich classes came out
    /// 12% wrong against the counting fit's 1.7%.
    pub rounds_choosing_strength: usize,
    pub settled_change_for_strength: f64,
    /// **A curve whose strength moved by less than this many natural-log units when last chosen
    /// keeps it from then on: 0.1**, and only the fit's first round searches the whole range (later
    /// starts of penalised counting, after each reading of F off pairs, search near the strength
    /// held). Zero chooses every curve's strength every round until the limits above, and searches
    /// the whole range at every start: the first version's schedule, kept as the measuring arm.
    /// On the fitting bench's cells the steady schedule cut the strength search's time 4.5 to 16-fold
    /// and, at 60 samples × 10 reads on two-base runs, the rounds from 534 to 288 (track 2's
    /// report, follow-up section).
    pub strength_steady_within: f64,
    /// **A class whose numbers moved by less than the tolerance in the last round keeps its last
    /// expected counts** instead of gathering them again, until the curves move it by more. Each
    /// round still chooses every class's numbers. Without it every class is counted every round
    /// until the slowest settles: at 2,000 samples × 3 reads on one-base runs that was 4,298
    /// class-rounds where the counting fit, which stops each class when it settles, needed 2,741.
    pub reuse_settled_counts: bool,
    /// **Whether standard errors are computed** for every number the fit reports (see
    /// [`PenalisedClassErrors`]).
    pub standard_errors: bool,
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
            kink_tail: 2.0,
            allele_pseudo_draws: 1.0 / 19.0,
            weigh_allele_draws_by_locus: false,
            pair_alternations: 0,
            rounds_choosing_strength: 200,
            settled_change_for_strength: 1e-5,
            strength_steady_within: 0.1,
            reuse_settled_counts: true,
            standard_errors: false,
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

/// **How a refit chooses its curve's penalty strength**, under [`Smoothing::MarginalLikelihood`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StrengthSearch {
    /// Over the whole range searched.
    Wide,
    /// Within two natural-log units of the strength held, widening when the best lands on an edge.
    Near,
    /// The strength held, unchanged.
    Held,
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
    /// The Dirichlet concentration, as the counting fit reports it. **`None` when no locus of the
    /// motif length has two samples with reads**: at one sample how much loci differ cannot be told
    /// from that sample's own inbreeding, and the project reports such a number as absent rather
    /// than as a fitted value. The fit still carries a working value internally, read off the
    /// sample's own two copies, because the genotype prior needs one; it is not reported.
    pub concentration: Option<f64>,
    /// The same, as the between-locus share `1 / (1 + concentration)`.
    pub between_locus_share: Option<f64>,
    /// **The concentration the fit worked at**, reported or not: what the genotype prior used. Not
    /// an estimate of how much loci differ when [`concentration`](Self::concentration) is `None`;
    /// there for a scorer that needs every number to compute a likelihood.
    pub working_concentration: f64,
    /// Standard errors, when [`PenalisedFitConfig::standard_errors`] asked for them.
    pub errors: Option<PenalisedClassErrors>,
}

/// **One class's standard errors**, each on its number's own scale; `None` for a number the fit
/// does not report or could not place.
///
/// **What they are.** The fit's numbers are the peak of the counting fit's objective summed over the
/// classes, minus each curve's bend penalty. Its curvature there, `H` (the objective's curvature
/// `A` plus the penalty's `P`), would give the Bayesian covariance a penalised fit reports, `H⁻¹`, if
/// the objective were a true likelihood. It is not: it treats the samples at one locus as
/// independent, so `A` overstates what the reads say. The spread of the objective's slope from locus
/// to locus, `B`, measures it honestly, and the covariance reported is `H⁻¹ (B + P) H⁻¹`: `H⁻¹` with
/// `B` in place of `A`, the same substitution the sandwich formula makes for an unpenalised fit. The
/// penalty's part keeps the allowance a Bayesian interval makes for the curve's own bias.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PenalisedClassErrors {
    pub level: Option<f64>,
    pub shorter_share: Option<f64>,
    pub fall_off: Option<f64>,
    /// Of the concentration `c`, and of the between-locus share `1 / (1 + c)`.
    pub concentration: Option<f64>,
    pub between_locus_share: Option<f64>,
    /// One per allele class, `-allele_span ..= +allele_span`.
    pub allele_length_variation: Vec<f64>,
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
    /// `None` when no locus has two samples with reads, as [`PenalisedClass::concentration`].
    pub between_locus_share: Option<ParameterCurve>,
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
    /// Of the M-step's time, the part spent in rounds that chose strengths, and how many such
    /// rounds there were.
    pub strength_search_seconds: f64,
    pub rounds_choosing_strength: usize,
    /// Time reading F off pairs of samples (building the tables and fitting F's curve), and how
    /// many times it was read.
    pub pair_seconds: f64,
    pub pair_alternations: usize,
    /// Rounds of each run of penalised counting: the first, then one after each reading of F.
    pub rounds_per_counting: Vec<usize>,
    /// Class-rounds whose expected counts were gathered (the rest reused a settled class's).
    pub class_rounds_counted: usize,
    /// Time computing standard errors.
    pub error_seconds: f64,
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
            strength_search_seconds: 0.0,
            rounds_choosing_strength: 0,
            pair_seconds: 0.0,
            pair_alternations: 0,
            rounds_per_counting: Vec::new(),
            class_rounds_counted: 0,
            error_seconds: 0.0,
        };
    };
    period.count(classes, homozygote_excess, config, fit, false);
    for alternation in 0..fit.pair_alternations {
        let started = Instant::now();
        let Some(tables) = period.pair_tables(classes, homozygote_excess, config) else {
            break;
        };
        let spreads = period.between_curve.refit(
            &tables,
            fit,
            if alternation == 0 {
                StrengthSearch::Wide
            } else {
                StrengthSearch::Near
            },
        );
        period.pair_alternations += 1;
        let mut moved: f64 = 0.0;
        for (point, class) in period.class_at_point.iter().enumerate() {
            if class.is_some() {
                let concentration = float::exp(-spreads[point]).clamp(1e-3, 1e4);
                moved = moved
                    .max((float::ln(concentration) - float::ln(period.concentration[point])).abs());
                period.concentration[point] = concentration;
            }
        }
        period.pair_seconds += started.elapsed().as_secs_f64();
        period.count(classes, homozygote_excess, config, fit, true);
        if moved < 0.01 {
            break;
        }
    }
    let started = Instant::now();
    let errors = if fit.standard_errors {
        period.standard_errors(classes, homozygote_excess, config)
    } else {
        vec![None; period.points()]
    };
    period.error_seconds = started.elapsed().as_secs_f64();
    period.into_fit(classes, errors)
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
    /// Per grid point, the allele draws a locus of its class holds on average: two per sample with
    /// reads there. Zero where there is no class.
    allele_draws_per_locus: Vec<f64>,
    /// Whether any locus of the motif length has two samples with reads, without which the
    /// between-locus share is not reported.
    between_locus_share_measurable: bool,
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
    strength_search_seconds: f64,
    rounds_choosing_strength: usize,
    pair_seconds: f64,
    pair_alternations: usize,
    rounds_per_counting: Vec<usize>,
    class_rounds_counted: usize,
    /// Per grid point, how far the class's numbers moved in the last round, and the expected
    /// counts last gathered for it.
    class_change: Vec<f64>,
    last_counts: Vec<Option<ClassCounts>>,
    error_seconds: f64,
    /// Whether F was read off pairs of samples, so is not a number of the counting objective.
    spread_from_pairs: bool,
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
            between_locus_share_measurable: with_reads.iter().any(|&index| {
                classes[index]
                    .tracts
                    .iter()
                    .any(|tract| tract.samples.len() >= 2)
            }),
            allele_draws_per_locus: class_at_point
                .iter()
                .map(|class| {
                    class.map_or(0.0, |index| {
                        let evidence = &classes[index];
                        2.0 * evidence.rows_with_reads() as f64
                            / evidence.tracts_with_reads().max(1) as f64
                    })
                })
                .collect(),
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
            strength_search_seconds: 0.0,
            rounds_choosing_strength: 0,
            pair_seconds: 0.0,
            pair_alternations: 0,
            rounds_per_counting: Vec::new(),
            class_rounds_counted: 0,
            class_change: vec![f64::INFINITY; points],
            last_counts: vec![None; points],
            error_seconds: 0.0,
            spread_from_pairs: false,
        })
    }

    fn points(&self) -> usize {
        self.class_at_point.len()
    }

    /// **How much more the allele draws of one locus vary than independent draws would**:
    /// `1 + (m − 1) · F`, with `m` the draws a locus holds and `F` the class's between-locus
    /// share. It is the variance of a Dirichlet-multinomial count over a multinomial one, the
    /// design effect of draws clustered by locus.
    ///
    /// **Why the allele counts are divided by it.** The E-step counts every sample's alleles as
    /// independent draws from the class, which is what the counting fit's objective assumes; but
    /// samples at one locus share its frequencies. At 2,000 samples a 5-locus class then holds
    /// 20,000 draws and looks measured to the curve fit, which chose almost no smoothing for the
    /// allele length variation and left its thin classes no better than the counting fit's
    /// (track 2's report, §3). Dividing every count of a class by the same number leaves the
    /// class's own shares where they were and gives its curve the information it actually has.
    fn allele_draw_inflation(&self, point: usize) -> f64 {
        let between_locus_share = 1.0 / (1.0 + self.concentration[point]);
        1.0 + (self.allele_draws_per_locus[point] - 1.0).max(0.0) * between_locus_share
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
            let search =
                if self.rounds == 1 || (rounds_here == 1 && fit.strength_steady_within <= 0.0) {
                    StrengthSearch::Wide
                } else if self.rounds <= fit.rounds_choosing_strength
                    && self.last_change >= fit.settled_change_for_strength
                {
                    StrengthSearch::Near
                } else {
                    StrengthSearch::Held
                };
            let started = Instant::now();
            // Classes side by side, each also spreading its loci over the pool: at one sample a
            // class's loci alone leave most of the machine idle.
            let reuse = |point: usize| {
                fit.reuse_settled_counts
                    && rounds_here > 1
                    && self.class_change[point] < fit.tolerance
                    && self.last_counts[point].is_some()
            };
            self.class_rounds_counted += self
                .class_at_point
                .iter()
                .enumerate()
                .filter(|(point, class)| class.is_some() && !reuse(*point))
                .count();
            let counts: Vec<Option<ClassCounts>> = self
                .class_at_point
                .par_iter()
                .enumerate()
                .map(|(point, class)| {
                    if reuse(point) {
                        return self.last_counts[point].clone();
                    }
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
                search,
            );
            let shorter = self.shorter_curve.refit(
                &binomial(&|c| c.shorter, &|c| c.shorter + c.longer),
                fit,
                search,
            );
            // A slip of k steps has probability (1 − f) f^(k − 1): k − 1 "one step further" against
            // one "stopped".
            let fall_off = self.fall_off_curve.refit(
                &binomial(&|c| c.steps - c.stepped, &|c| c.steps),
                fit,
                search,
            );
            let between = if spread_held {
                None
            } else {
                Some(self.between_curve.refit(
                    &binomial(&|c| c.shared_in_tract, &|c| c.pairs_not_by_descent),
                    fit,
                    search,
                ))
            };
            // Samples at one locus share its frequencies, so a class's allele draws carry less than
            // independent draws would: divide them by the clustering's variance inflation.
            let draws: Vec<Vec<f64>> = counts
                .iter()
                .enumerate()
                .map(|(point, count)| {
                    count.as_ref().map_or_else(Vec::new, |c| {
                        let inflation = if fit.weigh_allele_draws_by_locus {
                            self.allele_draw_inflation(point)
                        } else {
                            1.0
                        };
                        c.draws
                            .iter()
                            .map(|d| d / inflation + fit.allele_pseudo_draws)
                            .collect()
                    })
                })
                .collect();
            let share_logits = self.shares_curve.refit(
                &CurveData::Multinomial {
                    counts: draws,
                    reference: self.reference_allele,
                },
                fit,
                search,
            );
            let elapsed = started.elapsed().as_secs_f64();
            self.maximisation_seconds += elapsed;
            if search != StrengthSearch::Held {
                self.strength_search_seconds += elapsed;
                self.rounds_choosing_strength += 1;
            }

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
                let class_change = new_shares
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
                    .fold(0.0, f64::max);
                self.class_change[point] = class_change;
                change = change.max(class_change);
                self.slippage[point] = new_slippage;
                self.shares[point] = new_shares;
                self.concentration[point] = new_concentration;
            }
            self.last_counts = counts;
            self.last_change = change;
            if self.last_change < fit.tolerance {
                break;
            }
        }
        self.rounds_per_counting.push(rounds_here);
        if spread_held {
            self.spread_from_pairs = true;
        }
    }

    /// **The standard errors of every number at every repeat count with a class**
    /// ([`PenalisedClassErrors`]): `H⁻¹ (B + P) H⁻¹` over the curves' unbounded values, carried to
    /// each number's own scale.
    ///
    /// - `A`, the counting objective's curvature, class by class: the change in the class's slope
    ///   over a small step of each of its numbers either way, the slope from the expected counts
    ///   (Fisher's identity, as track 1 computes the counting fit's errors);
    /// - `B`, the sum over loci of each locus's slope times itself: the locus is the unit that
    ///   repeats independently;
    /// - `P`, each curve's bend penalty at the strength and bend weights it was fitted at.
    ///
    /// F read off pairs of samples is the answer to a second objective, so it is handled on its own:
    /// the same formula over F's curve, its slope and curvature taken numerically locus by locus
    /// from the pairwise likelihood, and the other numbers' errors taken with F held. The half
    /// allele draw a class adds is left out of `H` and `P`.
    fn standard_errors(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        config: &SsrFitConfig,
    ) -> Vec<Option<PenalisedClassErrors>> {
        let points = self.points();
        let shares_block = self.alleles - 1;
        let f_in_counting = !self.spread_from_pairs && self.between_locus_share_measurable;
        let share_offset = 3 + usize::from(f_in_counting);
        let width = share_offset + shares_block;
        let size = points * width;
        let reference = self.reference_allele;
        let theta_at = |point: usize| -> Vec<f64> {
            let mut theta = vec![
                self.level_curve.values[point],
                self.shorter_curve.values[point],
                self.fall_off_curve.values[point],
            ];
            if f_in_counting {
                theta.push(self.between_curve.values[point]);
            }
            theta.extend_from_slice(
                &self.shares_curve.values[point * shares_block..(point + 1) * shares_block],
            );
            theta
        };
        let numbers_of = |point: usize, theta: &[f64]| {
            let slippage = Slippage {
                level: logistic(theta[0]),
                shorter_share: logistic(theta[1]),
                fall_off: logistic(theta[2]).clamp(0.0, 0.99),
            };
            let concentration = if f_in_counting {
                (1.0 / logistic(theta[3]) - 1.0).clamp(1e-3, 1e4)
            } else {
                self.concentration[point]
            };
            let shares = softmax_with_reference(&theta[share_offset..], reference);
            (slippage, shares, concentration)
        };
        // The slope of the counting objective in the unbounded numbers, from the expected counts.
        let slope_of = |counts: &ClassCounts, theta: &[f64]| -> Vec<f64> {
            let mut slope = Vec::with_capacity(width);
            let level = logistic(theta[0]);
            slope.push(counts.slipped - level * (counts.slipped + counts.exact));
            let shorter = logistic(theta[1]);
            slope.push(counts.shorter - shorter * (counts.shorter + counts.longer));
            let fall_off = logistic(theta[2]);
            slope.push((counts.steps - counts.stepped) - fall_off * counts.steps);
            if f_in_counting {
                let share = logistic(theta[3]);
                slope.push(counts.shared_in_tract - share * counts.pairs_not_by_descent);
            }
            let shares = softmax_with_reference(&theta[share_offset..], reference);
            let total: f64 = counts.draws.iter().sum();
            for (allele, draws) in counts.draws.iter().enumerate() {
                if allele != reference {
                    slope.push(draws - total * shares[allele]);
                }
            }
            slope
        };
        let class_slope = |evidence: &StratumEvidence, point: usize, theta: &[f64]| {
            let (slippage, shares, concentration) = numbers_of(point, theta);
            let counts = expected_counts(
                evidence,
                slippage,
                &shares,
                concentration,
                homozygote_excess,
                config,
            );
            slope_of(&counts, theta)
        };
        let step = 1e-3;
        // Each class's curvature and spread blocks, classes side by side.
        let blocks: Vec<Option<(Vec<f64>, Vec<f64>)>> = self
            .class_at_point
            .par_iter()
            .enumerate()
            .map(|(point, class)| {
                let evidence = &classes[(*class)?];
                let theta = theta_at(point);
                let mut curvature = vec![0.0; width * width];
                for column in 0..width {
                    let mut up = theta.clone();
                    up[column] += step;
                    let mut down = theta.clone();
                    down[column] -= step;
                    let (above, below) = (
                        class_slope(evidence, point, &up),
                        class_slope(evidence, point, &down),
                    );
                    for row in 0..width {
                        let value = -(above[row] - below[row]) / (2.0 * step);
                        // Symmetrised: the objective's curvature is symmetric, its numerical
                        // estimate nearly so.
                        curvature[row * width + column] += 0.5 * value;
                        curvature[column * width + row] += 0.5 * value;
                    }
                }
                let by_locus: Vec<Vec<f64>> = evidence
                    .tracts
                    .par_iter()
                    .filter(|tract| !tract.samples.is_empty())
                    .map(|tract| {
                        let one = StratumEvidence {
                            stratum: evidence.stratum,
                            tracts: vec![tract.clone()],
                            read_span: evidence.read_span,
                            groups: evidence.groups,
                            libraries: Arc::clone(&evidence.libraries),
                            tracts_over_guard_threshold: 0,
                            reads_reaching_not_crossing: 0,
                            guard_reads: 0,
                            bases_compared: 0,
                            mismatching_bases: 0,
                        };
                        class_slope(&one, point, &theta)
                    })
                    .collect();
                let mut spread = vec![0.0; width * width];
                for slope in &by_locus {
                    for row in 0..width {
                        for column in 0..width {
                            spread[row * width + column] += slope[row] * slope[column];
                        }
                    }
                }
                Some((curvature, spread))
            })
            .collect();
        let mut curvature = vec![0.0; size * size];
        let mut spread = vec![0.0; size * size];
        for (point, block) in blocks.iter().enumerate() {
            let Some((block_curvature, block_spread)) = block else {
                continue;
            };
            let base = point * width;
            for row in 0..width {
                for column in 0..width {
                    curvature[(base + row) * size + base + column] =
                        block_curvature[row * width + column];
                    spread[(base + row) * size + base + column] =
                        block_spread[row * width + column];
                }
            }
        }
        // The bend penalties, at each curve's strength and bend weights.
        let mut penalty = vec![0.0; size * size];
        let mut add_curve = |curve: &CurveState, first_component: usize| {
            let strength = float::exp(curve.ln_smoothing);
            let stencil = [1.0, -2.0, 1.0];
            for bend in 0..curve.bends() {
                for component in 0..curve.block {
                    let weight = strength * curve.bend_weights[bend * curve.block + component];
                    for (i, a) in stencil.iter().enumerate() {
                        for (j, b) in stencil.iter().enumerate() {
                            let row = (bend + i) * width + first_component + component;
                            let column = (bend + j) * width + first_component + component;
                            penalty[row * size + column] += weight * a * b;
                        }
                    }
                }
            }
        };
        add_curve(&self.level_curve, 0);
        add_curve(&self.shorter_curve, 1);
        add_curve(&self.fall_off_curve, 2);
        if f_in_counting {
            add_curve(&self.between_curve, 3);
        }
        add_curve(&self.shares_curve, share_offset);
        let Some(covariance) = penalised_sandwich(&curvature, &spread, &penalty, size) else {
            return vec![None; points];
        };
        // F off pairs of samples: the same formula over F's curve alone.
        let pairs_variance: Option<Vec<f64>> = (self.spread_from_pairs
            && self.between_locus_share_measurable)
            .then(|| {
                let mut curvature = vec![0.0; points * points];
                let mut spread = vec![0.0; points * points];
                let step = 0.02;
                for (point, class) in self.class_at_point.iter().enumerate() {
                    let Some(index) = class else { continue };
                    let pairs = SamplePairs::of(
                        &classes[*index],
                        self.slippage[point],
                        homozygote_excess,
                        config,
                    );
                    if pairs.is_empty() {
                        continue;
                    }
                    let theta = self.between_curve.values[point];
                    let at =
                        |value: f64| pairs.log_likelihood_by_tract(&self.shares[point], -value);
                    let (above, middle, below) = (at(theta + step), at(theta), at(theta - step));
                    let mut bend = 0.0;
                    let mut slopes = 0.0;
                    for ((up, mid), down) in above.iter().zip(&middle).zip(&below) {
                        bend += (up - 2.0 * mid + down) / (step * step);
                        let slope = (up - down) / (2.0 * step);
                        slopes += slope * slope;
                    }
                    curvature[point * points + point] = -bend;
                    spread[point * points + point] = slopes;
                }
                let mut penalty = vec![0.0; points * points];
                let strength = float::exp(self.between_curve.ln_smoothing);
                let stencil = [1.0, -2.0, 1.0];
                for bend in 0..self.between_curve.bends() {
                    let weight = strength * self.between_curve.bend_weights[bend];
                    for (i, a) in stencil.iter().enumerate() {
                        for (j, b) in stencil.iter().enumerate() {
                            penalty[(bend + i) * points + bend + j] += weight * a * b;
                        }
                    }
                }
                penalised_sandwich(&curvature, &spread, &penalty, points)
                    .map(|covariance| (0..points).map(|p| covariance[p * points + p]).collect())
            })
            .flatten();
        (0..points)
            .map(|point| {
                self.class_at_point[point]?;
                let variance = |component: usize| {
                    let index = point * width + component;
                    covariance[index * size + index].max(0.0)
                };
                let probability_error = |component: usize, value: f64| {
                    Some(value * (1.0 - value) * variance(component).sqrt())
                };
                let slippage = self.slippage[point];
                let concentration = self.concentration[point];
                let share = 1.0 / (1.0 + concentration);
                let f_variance = if !self.between_locus_share_measurable {
                    None
                } else if f_in_counting {
                    Some(variance(3))
                } else {
                    pairs_variance.as_ref().map(|v| v[point].max(0.0))
                };
                // The shares through the softmax: d share_a / d ratio_b = share_a (δ_ab − share_b).
                let shares = &self.shares[point];
                let others: Vec<usize> = (0..self.alleles).filter(|a| *a != reference).collect();
                let allele_length_variation = (0..self.alleles)
                    .map(|allele| {
                        let gradient: Vec<f64> = others
                            .iter()
                            .map(|&other| {
                                let same = if other == allele { 1.0 } else { 0.0 };
                                shares[allele] * (same - shares[other])
                            })
                            .collect();
                        let mut total = 0.0;
                        for (b, gb) in gradient.iter().enumerate() {
                            for (c, gc) in gradient.iter().enumerate() {
                                let row = point * width + share_offset + b;
                                let column = point * width + share_offset + c;
                                total += gb * covariance[row * size + column] * gc;
                            }
                        }
                        total.max(0.0).sqrt()
                    })
                    .collect();
                Some(PenalisedClassErrors {
                    level: probability_error(0, slippage.level),
                    shorter_share: probability_error(1, slippage.shorter_share),
                    fall_off: probability_error(2, slippage.fall_off),
                    // c = e^(−logit F), so dc / dlogit F = −c.
                    concentration: f_variance.map(|v| concentration * v.sqrt()),
                    between_locus_share: f_variance.map(|v| share * (1.0 - share) * v.sqrt()),
                    allele_length_variation,
                })
            })
            .collect()
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

    fn into_fit(
        self,
        classes: &[StratumEvidence],
        mut errors: Vec<Option<PenalisedClassErrors>>,
    ) -> PenalisedPeriodFit {
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
                    concentration: self
                        .between_locus_share_measurable
                        .then_some(self.concentration[point]),
                    between_locus_share: self
                        .between_locus_share_measurable
                        .then(|| 1.0 / (1.0 + self.concentration[point])),
                    working_concentration: self.concentration[point],
                    errors: errors[point].take(),
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
            between_locus_share: self.between_locus_share_measurable.then(|| {
                scalar_curve(&self.between_curve, &|value| {
                    logistic(value).clamp(1.0 / (1.0 + 1e4), 1.0 / (1.0 + 1e-3))
                })
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
            strength_search_seconds: self.strength_search_seconds,
            rounds_choosing_strength: self.rounds_choosing_strength,
            pair_seconds: self.pair_seconds,
            pair_alternations: self.pair_alternations,
            rounds_per_counting: self.rounds_per_counting,
            class_rounds_counted: self.class_rounds_counted,
            error_seconds: self.error_seconds,
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
    /// Whether the strength has stopped moving ([`PenalisedFitConfig::strength_steady_within`]),
    /// after which a search near it keeps it.
    steady: bool,
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

/// **What a bend costs, in standard deviations `size` of the curve's smoothness**: `size² / 2` up to
/// `knee`, and beyond it growing only as `tail · ln(size / knee)`. Its derivative over `size` is
/// `size · weight`, the bend's weight one up to the knee and `tail / size²` beyond.
///
/// As a prior on a bend this is a normal distribution with tails falling as a power, `size^−tail`,
/// so a bend many deviations wide is rare but not ruled out; `tail` must exceed one for the prior to
/// be a distribution. Its normalising constant does not depend on the strength once the bend is
/// measured in deviations, which is what lets the marginal likelihood compare strengths.
#[derive(Debug, Clone, Copy, PartialEq)]
struct BendPrior {
    knee: f64,
    tail: f64,
}

impl BendPrior {
    /// Every bend costs its square: no knee.
    const QUADRATIC: Self = Self {
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
                let (ln_smoothing, fitted) =
                    self.choose_strength(data, knee, search == StrengthSearch::Wide);
                if search == StrengthSearch::Near
                    && config.strength_steady_within > 0.0
                    && (ln_smoothing - self.ln_smoothing).abs() < config.strength_steady_within
                {
                    self.steady = true;
                }
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
    fn choose_strength(
        &self,
        data: &CurveData,
        knee: BendPrior,
        search_widely: bool,
    ) -> (f64, CurveFit) {
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
        let mut ln_determinant = 0.0;
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
// The penalised sandwich, in dense matrices
// ---------------------------------------------------------------------

/// **`H⁻¹ (B + P) H⁻¹` with `H = A + P` plus the ridge**, every matrix `size × size` row by row.
/// `None` when `H` cannot be made positive definite by a ridge up to a thousandth of its mean
/// diagonal.
fn penalised_sandwich(
    curvature: &[f64],
    spread: &[f64],
    penalty: &[f64],
    size: usize,
) -> Option<Vec<f64>> {
    let mean_diagonal = (0..size)
        .map(|i| curvature[i * size + i] + penalty[i * size + i])
        .sum::<f64>()
        / size as f64;
    let mut inverse = None;
    // The fit's own ridge first; then larger ones, as shares of the mean diagonal.
    let scale = mean_diagonal.abs().max(1.0);
    for ridge in [
        RIDGE,
        1e-9 * scale,
        1e-7 * scale,
        1e-5 * scale,
        1e-3 * scale,
    ] {
        let mut matrix: Vec<f64> = curvature.iter().zip(penalty).map(|(a, p)| a + p).collect();
        for i in 0..size {
            matrix[i * size + i] += ridge;
        }
        if let Some(found) = dense_inverse(&matrix, size) {
            inverse = Some(found);
            break;
        }
    }
    let inverse = inverse?;
    let middle: Vec<f64> = spread.iter().zip(penalty).map(|(b, p)| b + p).collect();
    let left = dense_product(&inverse, &middle, size);
    Some(dense_product(&left, &inverse, size))
}

/// The inverse of a symmetric positive-definite matrix by its Cholesky factor; `None` when it is
/// not positive definite.
fn dense_inverse(matrix: &[f64], size: usize) -> Option<Vec<f64>> {
    let mut lower = vec![0.0; size * size];
    for i in 0..size {
        for j in 0..=i {
            let mut sum = matrix[i * size + j];
            for k in 0..j {
                sum -= lower[i * size + k] * lower[j * size + k];
            }
            if i == j {
                if sum <= 0.0 || !sum.is_finite() {
                    return None;
                }
                lower[i * size + i] = sum.sqrt();
            } else {
                lower[i * size + j] = sum / lower[j * size + j];
            }
        }
    }
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

fn dense_product(left: &[f64], right: &[f64], size: usize) -> Vec<f64> {
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
            let fitted = state.refit(&data, &config, StrengthSearch::Wide);
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
            allele_pseudo_draws: 0.0,
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
                (float::ln(fitted.concentration.expect("ten samples"))
                    - float::ln(alone.concentration))
                .abs()
                    < 1e-3,
                "{fitted:?} against {alone:?}"
            );
        }
    }
}
