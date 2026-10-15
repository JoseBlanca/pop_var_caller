//! **Every repeat class of one motif length fitted at once, each number a smooth curve across the
//! reference repeat count** — penalised expectation maximisation, a research prototype
//! (`doc/devel/ng/research/str_fit_research_plan.md`, avenue D), being promoted to the caller's
//! repeat-tract fit (spec `doc/devel/ng/spec/str_penalised_fit.md`): [`fit_repeat_tracts`] is the
//! entry point that fits a whole cohort and returns what the rest of the caller reads.
//!
//! The counting fit ([`expected_counts`](mod@super::expected_counts)) fits each class — one motif
//! length, one reference repeat count — on its own, and the curves across repeat count are drawn
//! afterwards through its answers. Here each round gathers every class's expected counts exactly as
//! the counting fit does ([`expected_counts`](fn@expected_counts)), and then chooses all the
//! classes' numbers together: each number is the best fit to the counts **subject to a penalty on
//! how sharply it bends from one repeat count to the next**. A thin class then borrows from its neighbours during the fit, by as
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

use super::expected_counts::{
    CONCENTRATION_SEARCH, ClassCounts, DistinctRows, SamplePairs,
    concentration_stopped_inside_its_search, expected_counts, expected_counts_by_distinct_rows,
};
use super::*;
use crate::parameter_estimation::joint::slippage_curve::PenalisedPoint;

/// **How the penalised fit runs.** [`production`](Self::production) is what `estimate-parameters`
/// runs (spec `str_penalised_fit.md` §2.7); [`prototype`](Self::prototype) keeps the research
/// prototype's own settings, which the research examples start from and which differ from
/// production's. There is deliberately no `Default`, so that no caller gets either by accident.
#[derive(Debug, Clone, PartialEq)]
pub struct PenalisedFitConfig {
    /// **How far from the reference length, in repeats, a class's alleles may sit**: the fit
    /// places every class's chromosomes over `-allele_span ..= +allele_span`, the outermost two
    /// meaning "this far or further". Production's is [`ALLELE_SPAN`], 9.
    ///
    /// **Coupled to [`allele_pseudo_draws`](Self::allele_pseudo_draws)**, which holds one draw a
    /// class spread over this window's lengths. Change the two together, through
    /// [`with_allele_span`](Self::with_allele_span).
    pub allele_span: i32,
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
    /// one draw a class, spread evenly over the lengths** (1/19 of a draw each at ±9). Without it
    /// a length the reads stop supporting in one round gets a share near zero, so the next round's
    /// E-step gives it no draws and it can never return: EM has an absorbing boundary at zero,
    /// which costs a held-out locus carrying that length hundreds of log-likelihood units.
    ///
    /// **The value is per length, so it is coupled to [`allele_span`](Self::allele_span)**: the
    /// fit adds it to each of the window's lengths and does not rescale it. Setting the window
    /// alone changes the floor a class gets (at ±4, 1/19 a length is 9/19 of a draw a class).
    /// Change the two together, through [`with_allele_span`](Self::with_allele_span).
    ///
    /// **Why not half a draw a length (the Jeffreys prior).** At one sample and 100 loci a class
    /// holds about 120 draws once clustering is allowed for, and 19 half-draws spread evenly moved
    /// enough mass into the tails to bias the reference length's share by −0.13 (one-base runs) and
    /// −0.20 (two-base runs) on the fitting bench's tomato-like truth.
    pub allele_pseudo_draws: f64,
    /// **A second floor under the allele draws, sized as the counting fit sizes its own: this many
    /// loci's worth of a class's allele draws, spread evenly over the lengths** (the counting fit's
    /// `SHARE_PRIOR_TRACTS` is one). Production's is [`ALLELE_PRIOR_LOCI`], one. The floor of
    /// [`allele_pseudo_draws`](Self::allele_pseudo_draws) is one draw a class whatever the cohort,
    /// so from tens of samples up it no longer holds off zero a length that no locus of the class
    /// carries; this one grows with the cohort and does (spec §2.4). Zero, the prototype's default,
    /// leaves only the first floor.
    pub allele_prior_tracts: f64,
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

impl PenalisedFitConfig {
    /// **What `estimate-parameters` runs** (spec `str_penalised_fit.md` §2.7): the settings the
    /// fitting bench measured as the lead candidate.
    ///
    /// | setting | value |
    /// |---|---|
    /// | allele window | ±[`ALLELE_SPAN`], 9 |
    /// | start | [`STARTING_STUTTER`]: slip rate 0.05, shorter share 0.5, fall-off 0.3 |
    /// | rounds at most / stop when no number moves more than | [`MAX_ROUNDS`] 2,000 / [`TOLERANCE`] 10⁻⁷ |
    /// | floor under the allele lengths | [`ALLELE_PRIOR_LOCI`], one locus's worth, plus [`ALLELE_PSEUDO_DRAWS_A_CLASS`], one draw a class |
    /// | F read off pairs of samples, at most | [`PAIR_ALTERNATIONS`], 10 times |
    /// | kink knee / tail | [`KINK_KNEE`] 2.5 / [`KINK_TAIL`] 2 |
    /// | a curve's strength held once it moves less than | [`STRENGTH_STEADY_WITHIN`], 0.1 natural-log units |
    /// | a class whose numbers stopped moving reuses its counts | on |
    /// | standard errors | on |
    ///
    /// Two settings §2.7 does not list keep the prototype's values: strengths are chosen for at
    /// most 200 rounds, or until no number moves by more than 10⁻⁵. The rule that withholds errors
    /// on thin classes ([`LOCI_FOR_AN_ERROR`], [`LOCI_FOR_A_SHARE_ERROR`]) is not a setting of the
    /// fit: it is applied when the fit's numbers are handed to the rest of the caller.
    #[must_use]
    pub fn production() -> Self {
        Self {
            allele_span: ALLELE_SPAN,
            start: STARTING_STUTTER,
            max_rounds: MAX_ROUNDS,
            tolerance: TOLERANCE,
            smoothing: Smoothing::MarginalLikelihood,
            kink_knee: KINK_KNEE,
            kink_tail: KINK_TAIL,
            allele_pseudo_draws: pseudo_draws_a_length(ALLELE_SPAN),
            allele_prior_tracts: ALLELE_PRIOR_LOCI,
            pair_alternations: PAIR_ALTERNATIONS,
            rounds_choosing_strength: 200,
            settled_change_for_strength: 1e-5,
            strength_steady_within: STRENGTH_STEADY_WITHIN,
            reuse_settled_counts: true,
            standard_errors: true,
        }
    }

    /// **The research prototype's own settings, not production's** (production's are
    /// [`production`](Self::production)), kept for the research examples and a few tests: a looser
    /// start, at most 500 rounds, stop below 10⁻⁶, no second floor under the allele lengths, F not
    /// read off pairs, and no standard errors.
    #[must_use]
    pub fn prototype() -> Self {
        Self {
            allele_span: ALLELE_SPAN,
            start: Slippage {
                level: 0.2,
                shorter_share: 0.5,
                fall_off: 0.5,
            },
            max_rounds: 500,
            tolerance: 1e-6,
            smoothing: Smoothing::MarginalLikelihood,
            kink_knee: KINK_KNEE,
            kink_tail: KINK_TAIL,
            allele_pseudo_draws: pseudo_draws_a_length(ALLELE_SPAN),
            allele_prior_tracts: 0.0,
            pair_alternations: 0,
            rounds_choosing_strength: 200,
            settled_change_for_strength: 1e-5,
            strength_steady_within: STRENGTH_STEADY_WITHIN,
            reuse_settled_counts: true,
            standard_errors: false,
        }
    }

    /// **These settings with the allele window set to `allele_span` and the first floor rescaled to
    /// match**: [`allele_pseudo_draws`](Self::allele_pseudo_draws) becomes
    /// [`ALLELE_PSEUDO_DRAWS_A_CLASS`] spread over the new window's lengths, whatever it held
    /// before. At ±9 that is 1/19, the value [`production`](Self::production) holds, to the bit.
    #[must_use]
    pub fn with_allele_span(self, allele_span: i32) -> Self {
        Self {
            allele_span,
            allele_pseudo_draws: pseudo_draws_a_length(allele_span),
            ..self
        }
    }
}

/// **The first floor's draws a length in a window of `allele_span`**: one draw a class
/// ([`ALLELE_PSEUDO_DRAWS_A_CLASS`]) divided by the number of lengths, `2 × allele_span + 1`.
fn pseudo_draws_a_length(allele_span: i32) -> f64 {
    ALLELE_PSEUDO_DRAWS_A_CLASS / allele_classes_in(allele_span) as f64
}

/// **Where every class's stutter starts in production**: slip rate 0.05, shorter share 0.5,
/// fall-off 0.3 (spec §2.7). The bench's values; soft.
pub const STARTING_STUTTER: Slippage = Slippage {
    level: 0.05,
    shorter_share: 0.5,
    fall_off: 0.3,
};
/// **At most this many rounds of the loop**, over every run of penalised counting together
/// (spec §2.7). The bench's value; soft.
pub const MAX_ROUNDS: usize = 2000;
/// **The loop stops once no number moves by more than this in a round** (spec §2.7). The bench's
/// value; soft.
pub const TOLERANCE: f64 = 1e-7;
/// **One draw a class, spread evenly over the allele lengths**: the first floor under the allele
/// draws ([`PenalisedFitConfig::allele_pseudo_draws`] holds it per length). Spec §2.4, measured.
pub const ALLELE_PSEUDO_DRAWS_A_CLASS: f64 = 1.0;
/// **One locus's worth of allele draws a class, spread evenly over the lengths**: the second
/// floor ([`PenalisedFitConfig::allele_prior_tracts`]). Spec §2.4, measured.
pub const ALLELE_PRIOR_LOCI: f64 = 1.0;
/// Readings of F off pairs of samples, at most (spec §2.5). The bench's value; soft.
pub const PAIR_ALTERNATIONS: usize = 10;
/// **How many standard deviations a bend may reach before its penalty stops growing as its
/// square** (spec §2.3): the knee of the slippage curve's blend.
pub const KINK_KNEE: f64 = 2.5;
/// **What a bend beyond the knee costs per natural-log unit of its size** (spec §2.3); see
/// [`PenalisedFitConfig::kink_tail`] for how it was measured.
pub const KINK_TAIL: f64 = 2.0;
/// **A curve whose strength moved by less than this many natural-log units keeps it** (spec
/// §2.7). Measured.
pub const STRENGTH_STEADY_WITHIN: f64 = 0.1;
/// A share's error is withheld when the class's loci with reads times the share is below this
/// (spec §7, owner 2026-10-14). Rests on two settings; soft.
pub const LOCI_FOR_A_SHARE_ERROR: f64 = 3.0;
/// **Every error of a class is withheld below this many loci with reads** (spec §2.6; confirmed
/// by the owner, 2026-10-14). Errors on thinner classes were never measured; soft.
pub const LOCI_FOR_AN_ERROR: usize = 8;

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
    /// **How many of the class's reads slipped, by the fit's last expected counts**: each read's
    /// chance of having slipped, given its sample's genotype probabilities, summed over the class's
    /// reads. These are the counts the returned numbers were fitted to; they were gathered at the
    /// numbers one round before, which differ from the returned ones by less than the fit's
    /// tolerance when it converged.
    ///
    /// **Not the level times the reads.** The level is a point on the penalised curve, so for a
    /// thin class it is mostly the neighbours'; this count is what the class's own reads said.
    /// `None` only when no round ran (`max_rounds` of 0).
    pub expected_slipped_reads: Option<f64>,
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
/// (multipliers are not fitted here). The allele window is `fit.allele_span`.
///
/// With [`PenalisedFitConfig::pair_alternations`] above zero, the between-locus share is then read
/// off pairs of samples at a locus, as [`fit_class_by_counting_with_pairs`](super::expected_counts::fit_class_by_counting_with_pairs)
/// does, its curve penalised like every other: each round of pairs alternates with penalised
/// counting at the shares held.
///
/// **Memory: one class's evidence table at a time.** Every round, and the standard errors, visit
/// the classes one after another and use the threads inside each class, over its loci and distinct
/// rows. So at most one class's table of distinct rows is alive, whatever the thread count (in the
/// standard errors, also one locus's small table a busy thread, all from the same class). Between
/// classes the fit keeps only what does not grow with samples or reads: each class's numbers and
/// expected counts, a few dozen floats a class, and in the standard errors its curvature and spread
/// blocks, `width × width` floats each: about 900 floats a class at the default allele window of ±9
/// (`width` is 21, or 22 when the between-locus share is counted).
///
/// # Panics
///
/// When the classes are of more than one motif length, or two share a repeat count.
#[must_use]
pub fn fit_period_by_penalised_counting(
    classes: &[StratumEvidence],
    homozygote_excess: &[f64],
    fit: &PenalisedFitConfig,
) -> PenalisedPeriodFit {
    let classes: Vec<&StratumEvidence> = classes.iter().collect();
    fit_one_period(&classes, homozygote_excess, fit)
}

/// [`fit_period_by_penalised_counting`] over borrowed classes, so that [`fit_repeat_tracts`] can
/// hand one motif length's classes over from wherever they sit in its input without copying
/// their reads.
fn fit_one_period(
    classes: &[&StratumEvidence],
    homozygote_excess: &[f64],
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
    let Some(mut period) = Period::new(classes, fit) else {
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
    period.count(classes, homozygote_excess, fit, false);
    for alternation in 0..fit.pair_alternations {
        let started = Instant::now();
        let Some(tables) = period.pair_tables(classes, homozygote_excess, fit.allele_span) else {
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
                let concentration = float::exp(-spreads[point])
                    .clamp(CONCENTRATION_SEARCH.0, CONCENTRATION_SEARCH.1);
                moved = moved
                    .max((float::ln(concentration) - float::ln(period.concentration[point])).abs());
                period.concentration[point] = concentration;
            }
        }
        period.pair_seconds += started.elapsed().as_secs_f64();
        period.count(classes, homozygote_excess, fit, true);
        if moved < 0.01 {
            break;
        }
    }
    let started = Instant::now();
    let errors = if fit.standard_errors {
        period.standard_errors(classes, homozygote_excess, fit.allele_span)
    } else {
        vec![None; period.points()]
    };
    period.error_seconds = started.elapsed().as_secs_f64();
    period.into_fit(classes, errors)
}

/// **Fit every repeat-tract class of a cohort by penalised EM, one motif length at a time**
/// (spec `str_penalised_fit.md` §2), and return one [`StratumOutcome`] per class handed in, in
/// the same order. This is the repeat-tract fit `estimate-parameters` is to run.
///
/// The classes are grouped by motif length, and each group is fitted together by
/// [`fit_period_by_penalised_counting`], since every curve runs across a motif length's repeat
/// counts. What each class gets back:
///
/// - **A class with reads at even one locus is [`StratumOutcome::Fitted`].** Its numbers are
///   points on its motif length's penalised curves: the level and both shares carry
///   [`LevelSource::PenalisedCurve`] / [`ShareSource::PenalisedCurve`], inside the curves' range,
///   with no parametric curve beside them. It was read from every sample, so
///   [`StratumFit::samples_fitted_on`] is `None`, and [`StratumFit::climb`] is `None`. There is no
///   floor of loci below which a class is refused.
/// - **A class with no reads at any locus is [`StratumOutcome::Refused`]** for
///   [`StratumRefusal::NoSpanningReads`], and takes no value off the curves.
/// - **No class is [`StratumOutcome::Derived`].**
///
/// A class's concentration is `None` when no locus of its motif length holds two samples with
/// reads, or when the estimate stopped within 1% of an end of its search
/// ([`concentration_stopped_inside_its_search`], spec §2.5). The concentration the fit worked at
/// is never reported in its place.
///
/// Each number carries its standard error or the reason it has none, by the table of arch
/// `str_penalised_fit.md` §1.4: [`StratumError::Estimated`] only for a positive, finite error that
/// places the number; [`StratumError::TooFewLoci`] for every number of a class with fewer than
/// [`LOCI_FOR_AN_ERROR`] loci with reads, and for a share worth fewer than
/// [`LOCI_FOR_A_SHARE_ERROR`] loci; [`StratumError::NotPlaced`] for the rest, an absent
/// concentration's error included.
///
/// **Memory: one motif length at a time, and one class at a time within it.** Each motif
/// length's classes are handed to its fit by reference, never copied, and everything that fit
/// built beyond the outcomes — its curves, its timings — is dropped before the next motif length
/// starts. Within a motif length, see [`fit_period_by_penalised_counting`].
///
/// `homozygote_excess` is indexed by sample, as the evidence's samples are.
///
/// # Panics
///
/// - When a class's evidence declares more than one slippage group. The fit scores every library
///   at one stutter, and a class's [`StratumFit::slippage`] holds that one; stutter per slippage
///   group is not built (`str_library_stutter.md`).
/// - When a class's libraries do not fit together ([`StratumEvidence::check_libraries`]), or two
///   classes with reads share a motif length and a repeat count. Two classes at one repeat count
///   where one has no reads do not panic: each comes back, the one without reads refused.
#[must_use]
pub fn fit_repeat_tracts(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &PenalisedFitConfig,
) -> Vec<StratumOutcome> {
    for evidence in strata {
        assert_eq!(
            evidence.groups, 1,
            "period {} at {} repeats: the evidence declares {} slippage groups, and the \
             penalised fit fits one stutter for every library; stutter per slippage group is not \
             built (str_library_stutter.md)",
            evidence.stratum.period, evidence.stratum.reference_repeats, evidence.groups
        );
        evidence.check_libraries();
    }
    let mut classes_of_each_period: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for (index, evidence) in strata.iter().enumerate() {
        classes_of_each_period
            .entry(evidence.stratum.period)
            .or_default()
            .push(index);
    }
    let mut outcomes: Vec<Option<StratumOutcome>> = vec![None; strata.len()];
    for indices in classes_of_each_period.values() {
        let classes: Vec<&StratumEvidence> = indices.iter().map(|&index| &strata[index]).collect();
        let PenalisedPeriodFit {
            classes: fitted,
            curves,
            ..
        } = fit_one_period(&classes, homozygote_excess, config);
        for ((&index, evidence), fitted) in indices.iter().zip(&classes).zip(fitted) {
            outcomes[index] = Some(outcome_of_a_class(
                evidence,
                fitted,
                curves.as_ref(),
                config.standard_errors,
            ));
        }
    }
    outcomes
        .into_iter()
        .map(|outcome| outcome.expect("every class is in its motif length's group"))
        .collect()
}

/// **One class's [`StratumOutcome`], from what its motif length's fit returned for it.**
///
/// `fitted` is `None` for a class with no reads at any locus, which is refused. `curves` are the
/// motif length's, `None` only when no class of it had reads. `standard_errors` says whether the
/// fit was asked for errors: when it was not, the outcome carries none
/// ([`StratumFit::standard_errors`] is `None`); when it was, every number has an error or a reason
/// for having none ([`errors_of_a_class`]).
fn outcome_of_a_class(
    evidence: &StratumEvidence,
    fitted: Option<PenalisedClass>,
    curves: Option<&SmoothCurves>,
    standard_errors: bool,
) -> StratumOutcome {
    let Some(fitted) = fitted else {
        return StratumOutcome::Refused {
            stratum: evidence.stratum,
            tracts: evidence.tracts_with_reads(),
            reason: StratumRefusal::NoSpanningReads,
        };
    };
    let curves = curves.expect(
        "a fitted class has curves: its motif length's fit returns no curves only when no class of \
         it had reads, and then it fits no class",
    );
    let repeats = evidence.stratum.reference_repeats;
    let point_on = |curve: &ParameterCurve| PenalisedPoint {
        strength: curve.smoothing,
        eased_here: curve.eased_at.contains(&repeats),
    };
    let share_on = |curve: &ParameterCurve| ShareProvenance {
        source: ShareSource::PenalisedCurve(point_on(curve)),
        curve: None,
        reach: Some(CurveReach::Inside),
    };
    let concentration = fitted
        .concentration
        .filter(|estimate| concentration_stopped_inside_its_search(*estimate));
    let loci_with_reads = evidence.tracts_with_reads();
    let reads_crossing = evidence.spanning_reads();
    // The class's own reads that slipped, from the fit's expected counts, not the curve's level
    // times the reads: for a thin class that level is mostly the neighbours'.
    let slipped_reads = fitted.expected_slipped_reads;
    let standard_errors =
        standard_errors.then(|| errors_of_a_class(&fitted, concentration, loci_with_reads));
    StratumOutcome::Fitted(Box::new(StratumFit {
        stratum: evidence.stratum,
        slippage: vec![Some(fitted.slippage)],
        allele_length_variation: fitted.allele_length_variation,
        concentration,
        climb: None,
        samples_fitted_on: None,
        tracts_of_its_own: loci_with_reads,
        reads_crossing,
        level_provenance: vec![Some(LevelProvenance {
            source: LevelSource::PenalisedCurve(point_on(&curves.level)),
            curve: None,
            reach: Some(CurveReach::Inside),
            slipped_reads,
        })],
        shares_provenance: vec![Some(SharesProvenance {
            slipped_reads,
            shorter_share: share_on(&curves.shorter_share),
            fall_off: share_on(&curves.fall_off),
        })],
        standard_errors,
    }))
}

/// **Every standard error of one fitted class as the rest of the caller reads it**, or why a
/// number has none (arch `str_penalised_fit.md` §1.4). Each number is judged by the first rule
/// that applies:
///
/// | case | error |
/// |---|---|
/// | the concentration, where `concentration` is `None` | [`StratumError::NotPlaced`] |
/// | any number of a class with fewer than [`LOCI_FOR_AN_ERROR`] loci with reads | [`StratumError::TooFewLoci`] |
/// | a share whose expected count of loci, share × loci with reads, is below [`LOCI_FOR_A_SHARE_ERROR`] | [`StratumError::TooFewLoci`] |
/// | the fit gave no error, or one that is zero or not finite | [`StratumError::NotPlaced`] |
/// | an error above [`NOT_PLACED`] on the scale the number is fitted on | [`StratumError::NotPlaced`] |
/// | otherwise | [`StratumError::Estimated`] |
///
/// **The scale a number is fitted on** is the logit for the three stutter numbers and the
/// shares, so the error is divided by `value × (1 − value)`, and the log for the concentration, so
/// it is divided by the concentration. **No error is ever `Estimated(0)`**: the fit sets a
/// negative variance to zero, and a zero error is not placed.
///
/// `concentration` is the class's reported concentration, already `None` where it stopped at an
/// end of its search; the fit's own error beside an absent concentration is ignored. A class with
/// no errors from the fit (the motif length's covariance could not be inverted) gets every number
/// `NotPlaced`, the withheld ones aside.
fn errors_of_a_class(
    fitted: &PenalisedClass,
    concentration: Option<f64>,
    loci_with_reads: usize,
) -> StratumErrors {
    let none = PenalisedClassErrors::default();
    let errors = fitted.errors.as_ref().unwrap_or(&none);
    let too_few_loci = loci_with_reads < LOCI_FOR_AN_ERROR;
    let of_a_probability = |error: Option<f64>, value: f64| {
        if too_few_loci {
            StratumError::TooFewLoci
        } else {
            an_error_of(error, value * (1.0 - value))
        }
    };
    let slippage = fitted.slippage;
    StratumErrors {
        slippage: vec![Some(SlippageErrors {
            level: of_a_probability(errors.level, slippage.level),
            shorter_share: of_a_probability(errors.shorter_share, slippage.shorter_share),
            fall_off: of_a_probability(errors.fall_off, slippage.fall_off),
        })],
        allele_length_variation: fitted
            .allele_length_variation
            .iter()
            .enumerate()
            .map(|(allele, &share)| {
                if share * (loci_with_reads as f64) < LOCI_FOR_A_SHARE_ERROR {
                    StratumError::TooFewLoci
                } else {
                    of_a_probability(errors.allele_length_variation.get(allele).copied(), share)
                }
            })
            .collect(),
        concentration: match concentration {
            None => StratumError::NotPlaced,
            Some(_) if too_few_loci => StratumError::TooFewLoci,
            Some(concentration) => an_error_of(errors.concentration, concentration),
        },
    }
}

/// **One number's error, `error` on the number's own scale, judged on the scale it is fitted
/// on**: `error / slope` is its size there, `slope` being how fast the number moves with that
/// scale. [`StratumError::Estimated`] only for a positive, finite error whose size on the fitted
/// scale is at most [`NOT_PLACED`]; [`StratumError::NotPlaced`] for anything else, an absent error
/// included.
fn an_error_of(error: Option<f64>, slope: f64) -> StratumError {
    match error {
        Some(error) if error.is_finite() && error > 0.0 && error / slope <= NOT_PLACED => {
            StratumError::Estimated(error)
        }
        _ => StratumError::NotPlaced,
    }
}

/// The ln-spread grid a class's pairwise likelihood is tabulated on, as the logit of the
/// between-locus share (`−ln spread`): 57 points from a spread of 10⁴ to 10⁻³.
fn between_share_grid() -> Vec<f64> {
    let (low, high) = (
        -float::ln(CONCENTRATION_SEARCH.1),
        -float::ln(CONCENTRATION_SEARCH.0),
    );
    (0..57)
        .map(|step| low + (high - low) * f64::from(step) / 56.0)
        .collect()
}

/// **One motif length's classes as the fit carries them from round to round**: every number at
/// every repeat count of the grid, and each curve's state.
struct Period {
    /// Per grid point (consecutive repeat counts from `lowest`), the class there, if any.
    class_at_point: Vec<Option<usize>>,
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
    fn new(classes: &[&StratumEvidence], fit: &PenalisedFitConfig) -> Option<Self> {
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
        let alleles = allele_classes_in(fit.allele_span);
        Some(Self {
            between_locus_share_measurable: with_reads.iter().any(|&index| {
                classes[index]
                    .tracts
                    .iter()
                    .any(|tract| tract.samples.len() >= 2)
            }),
            class_at_point,
            lowest,
            alleles,
            reference_allele: fit.allele_span as usize,
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

    /// **Rounds of penalised counting** until no number moves by more than the tolerance or the
    /// rounds run out. With `spread_held`, the concentration stays where it is.
    fn count(
        &mut self,
        classes: &[&StratumEvidence],
        homozygote_excess: &[f64],
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
            // **One class at a time**, the threads working inside it (`expected_counts` spreads the
            // class's distinct rows over the pool). A class's table of distinct rows, with each
            // row's genotype columns scored at the round's numbers, is the round's largest memory,
            // and is dropped before the next class builds its own, so one table is alive at any
            // thread count (spec §2.2, §3). A parallel loop over classes
            // would build one table a busy thread, and more when a thread waiting on a class's
            // inner work picks up another class. What stays from class to class is its counts: a
            // few numbers and one per allele length.
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
                .iter()
                .enumerate()
                .map(|(point, class)| {
                    if reuse(point) {
                        return self.last_counts[point].clone();
                    }
                    class.map(|index| {
                        expected_counts(
                            classes[index],
                            self.slippage[point],
                            &self.shares[point],
                            self.concentration[point],
                            homozygote_excess,
                            fit.allele_span,
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
            // Each class's allele draws with the two floors under them (spec §2.4): one draw a
            // class, and `allele_prior_tracts` loci's worth of its draws, each spread evenly over
            // the lengths.
            let draws: Vec<Vec<f64>> = counts
                .iter()
                .enumerate()
                .map(|(point, count)| {
                    count.as_ref().map_or_else(Vec::new, |c| {
                        let prior_total = if fit.allele_prior_tracts > 0.0 {
                            let tracts = self.class_at_point[point]
                                .map_or(1, |index| classes[index].tracts_with_reads().max(1))
                                as f64;
                            let total: f64 = c.draws.iter().sum();
                            fit.allele_prior_tracts * total / tracts
                        } else {
                            0.0
                        };
                        let prior = prior_total / c.draws.len().max(1) as f64;
                        c.draws
                            .iter()
                            .map(|d| d + fit.allele_pseudo_draws + prior)
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
                    let share = logistic(b[point]).clamp(
                        1.0 / (1.0 + CONCENTRATION_SEARCH.1),
                        1.0 / (1.0 + CONCENTRATION_SEARCH.0),
                    );
                    (1.0 / share - 1.0).clamp(CONCENTRATION_SEARCH.0, CONCENTRATION_SEARCH.1)
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
        classes: &[&StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
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
                (1.0 / logistic(theta[3]) - 1.0)
                    .clamp(CONCENTRATION_SEARCH.0, CONCENTRATION_SEARCH.1)
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
        // `distinct` is `evidence`'s table of distinct rows, built once and reused for every
        // step: the table depends on the reads and the homozygote excesses, neither of which a
        // step moves.
        let class_slope =
            |evidence: &StratumEvidence, distinct: &DistinctRows, point: usize, theta: &[f64]| {
                let (slippage, shares, concentration) = numbers_of(point, theta);
                let counts = expected_counts_by_distinct_rows(
                    evidence,
                    distinct,
                    slippage,
                    &shares,
                    concentration,
                    allele_span,
                );
                slope_of(&counts, theta)
            };
        let step = 1e-3;
        // Each class's curvature and spread blocks, **one class at a time** as in `count`, so one
        // class's table of distinct rows is alive at once; the threads work inside the class,
        // over its loci and its distinct rows.
        let blocks: Vec<Option<(Vec<f64>, Vec<f64>)>> = self
            .class_at_point
            .iter()
            .enumerate()
            .map(|(point, class)| {
                let evidence = classes[(*class)?];
                let theta = theta_at(point);
                let distinct = DistinctRows::of(evidence, homozygote_excess);
                let mut curvature = vec![0.0; width * width];
                for column in 0..width {
                    let mut up = theta.clone();
                    up[column] += step;
                    let mut down = theta.clone();
                    down[column] -= step;
                    let (above, below) = (
                        class_slope(evidence, &distinct, point, &up),
                        class_slope(evidence, &distinct, point, &down),
                    );
                    for row in 0..width {
                        let value = -(above[row] - below[row]) / (2.0 * step);
                        // Symmetrised: the objective's curvature is symmetric, its numerical
                        // estimate nearly so.
                        curvature[row * width + column] += 0.5 * value;
                        curvature[column * width + row] += 0.5 * value;
                    }
                }
                // Gone before the loci build their own small tables.
                drop(distinct);
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
                        class_slope(
                            &one,
                            &DistinctRows::of(&one, homozygote_excess),
                            point,
                            &theta,
                        )
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
                        classes[*index],
                        self.slippage[point],
                        &self.shares[point],
                        homozygote_excess,
                        allele_span,
                    );
                    if pairs.is_empty() {
                        continue;
                    }
                    let theta = self.between_curve.values[point];
                    let at = |value: f64| pairs.log_likelihood_by_tract(-value);
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
        classes: &[&StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
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
                    classes[*index],
                    self.slippage[point],
                    &self.shares[point],
                    homozygote_excess,
                    allele_span,
                );
                if pairs.is_empty() {
                    return Vec::new();
                }
                grid.iter()
                    .map(|logit| pairs.log_likelihood(-logit))
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
        classes: &[&StratumEvidence],
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
                    expected_slipped_reads: self.last_counts[point]
                        .as_ref()
                        .map(|counts| counts.slipped),
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
                    logistic(value).clamp(
                        1.0 / (1.0 + CONCENTRATION_SEARCH.1),
                        1.0 / (1.0 + CONCENTRATION_SEARCH.0),
                    )
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
    use crate::parameter_estimation::joint::ssr_fit::expected_counts::fit_class_by_counting_at_share_prior;

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
            let config = PenalisedFitConfig::prototype();
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
            allele_span: config.allele_span,
            smoothing: Smoothing::Fixed(1e-8),
            allele_pseudo_draws: 0.0,
            kink_knee: f64::INFINITY,
            tolerance: 1e-10,
            max_rounds: 2000,
            ..PenalisedFitConfig::prototype()
        };
        let together = fit_period_by_penalised_counting(&classes, &excess, &fit);
        for (evidence, fitted) in classes.iter().zip(&together.classes) {
            let fitted = fitted.as_ref().expect("every class has reads");
            // The counting fit without its shares' prior: with no floor here either, the two fits
            // maximise the same objective.
            let alone = fit_class_by_counting_at_share_prior(
                evidence, fit.start, &excess, &config, 2000, 1e-10, 0.0,
            );
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
                    - float::ln(alone.concentration_counted_at))
                .abs()
                    < 1e-3,
                "{fitted:?} against {alone:?}"
            );
        }
    }

    /// **The same cohort gives the same bits on every run and at every thread count** (spec §3):
    /// four classes of one-base runs, 12 samples at 5 reads over 60 loci each, fitted with
    /// production's settings (F read off pairs, standard errors on) twice on 18 threads and once
    /// on one. Every number and every count must agree; only the times may differ.
    #[test]
    fn the_fit_is_the_same_on_every_run_and_thread_count() {
        let spectrum = spectrum_of(3);
        let classes: Vec<StratumEvidence> = [(8_u64, 0.03), (9, 0.05), (10, 0.08), (11, 0.12)]
            .iter()
            .enumerate()
            .map(|(index, (repeats, level))| {
                let truth = Slippage {
                    level: *level,
                    shorter_share: 0.8,
                    fall_off: 0.25,
                };
                let mut evidence =
                    draw_stratum(truth, &spectrum, 2.0, 0.3, 60, 12, 5, 1, 71 + index as u64);
                evidence.stratum = Stratum {
                    period: 1,
                    reference_repeats: *repeats,
                };
                evidence
            })
            .collect();
        let excess = [0.3; 12];
        let fit = PenalisedFitConfig::production().with_allele_span(2);
        let on_threads = |threads: usize| {
            let fitted = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("a thread pool builds")
                .install(|| fit_period_by_penalised_counting(&classes, &excess, &fit));
            PenalisedPeriodFit {
                expectation_seconds: 0.0,
                maximisation_seconds: 0.0,
                strength_search_seconds: 0.0,
                pair_seconds: 0.0,
                error_seconds: 0.0,
                ..fitted
            }
        };
        let first = on_threads(18);
        assert!(
            first.pair_alternations > 0,
            "F is read off pairs, so the pairs' sums are exercised"
        );
        // Debug prints every number in the shortest form that reads back to the same bits.
        let first = format!("{first:?}");
        assert_eq!(
            first,
            format!("{:?}", on_threads(18)),
            "two runs on 18 threads"
        );
        assert_eq!(
            first,
            format!("{:?}", on_threads(1)),
            "18 threads against one"
        );
    }

    /// **Production runs the settings of spec §2.7**, each from its named constant, and the
    /// prototype's two settings §2.7 does not list.
    #[test]
    fn production_runs_the_settings_of_the_spec() {
        let production = PenalisedFitConfig::production();
        let spec = PenalisedFitConfig {
            allele_span: 9,
            start: Slippage {
                level: 0.05,
                shorter_share: 0.5,
                fall_off: 0.3,
            },
            max_rounds: 2000,
            tolerance: 1e-7,
            smoothing: Smoothing::MarginalLikelihood,
            kink_knee: 2.5,
            kink_tail: 2.0,
            allele_pseudo_draws: 1.0 / 19.0,
            allele_prior_tracts: 1.0,
            pair_alternations: 10,
            rounds_choosing_strength: 200,
            settled_change_for_strength: 1e-5,
            strength_steady_within: 0.1,
            reuse_settled_counts: true,
            standard_errors: true,
        };
        assert_eq!(production, spec);
    }

    /// **A class drawn at a known stutter**: one-base or two-base runs, `loci` loci of `samples`
    /// samples at `depth` reads each, the three-length spectrum every drawn test here uses, and a
    /// slip rate rising with the repeat count. With `loci` of 0 the class has no reads at all.
    fn a_drawn_class(
        period: u8,
        repeats: u64,
        loci: usize,
        samples: usize,
        depth: u32,
        seed: u64,
    ) -> StratumEvidence {
        let truth = Slippage {
            level: 0.01 * (repeats as f64 - 5.0),
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let mut evidence = draw_stratum(
            truth,
            &spectrum_of(3),
            2.0,
            0.3,
            loci,
            samples,
            depth,
            1,
            seed,
        );
        evidence.stratum = Stratum {
            period,
            reference_repeats: repeats,
        };
        evidence
    }

    /// The fit `fit_repeat_tracts` runs in these tests: production's settings at a window of ±2,
    /// so a test takes seconds.
    fn a_narrow_production_fit() -> PenalisedFitConfig {
        PenalisedFitConfig::production().with_allele_span(2)
    }

    /// Every error a class carries, slippage, shares and concentration together.
    fn every_error_of(errors: &StratumErrors) -> Vec<StratumError> {
        let mut every: Vec<StratumError> = errors
            .slippage
            .iter()
            .flatten()
            .flat_map(|group| [group.level, group.shorter_share, group.fall_off])
            .collect();
        every.extend(&errors.allele_length_variation);
        every.push(errors.concentration);
        every
    }

    /// **Every class comes back in the order it was handed in**, fitted when it has reads and
    /// refused when it has none, never derived, and each fitted class's numbers are the points of
    /// its own motif length's fit.
    ///
    /// Six classes of three motif lengths, interleaved: one-base runs at 8, 10 (no reads) and 9
    /// repeats; two-base runs at 9 and 10; three-base runs at 7, with no reads, the only class of
    /// its motif length. Twelve samples at 5 reads, 30 loci a class.
    #[test]
    fn every_class_comes_back_in_order_fitted_with_reads_and_refused_without() {
        let strata = vec![
            a_drawn_class(2, 9, 30, 12, 5, 11),
            a_drawn_class(1, 8, 30, 12, 5, 12),
            a_drawn_class(1, 10, 0, 12, 5, 13),
            a_drawn_class(3, 7, 0, 12, 5, 14),
            a_drawn_class(2, 10, 30, 12, 5, 15),
            a_drawn_class(1, 9, 30, 12, 5, 16),
        ];
        let excess = [0.3; 12];
        let config = a_narrow_production_fit();
        let outcomes = fit_repeat_tracts(&strata, &excess, &config);

        assert_eq!(outcomes.len(), strata.len());
        for (outcome, evidence) in outcomes.iter().zip(&strata) {
            assert_eq!(outcome.stratum(), evidence.stratum, "in input order");
            match outcome {
                StratumOutcome::Refused {
                    tracts,
                    reason: StratumRefusal::NoSpanningReads,
                    ..
                } => {
                    assert_eq!(evidence.rows_with_reads(), 0, "{:?}", evidence.stratum);
                    assert_eq!(*tracts, 0);
                }
                StratumOutcome::Fitted(fit) => {
                    assert!(evidence.rows_with_reads() > 0);
                    assert_eq!(fit.climb, None);
                    assert_eq!(fit.samples_fitted_on, None);
                    assert_eq!(fit.slippage.len(), 1, "one slippage group");
                    assert_eq!(fit.tracts_of_its_own, 30);
                    assert_eq!(fit.reads_crossing, evidence.spanning_reads());
                    let level = fit.level_provenance[0].expect("group 0 has numbers");
                    assert!(matches!(level.source, LevelSource::PenalisedCurve(_)));
                    assert_eq!(level.curve, None);
                    assert_eq!(level.reach, Some(CurveReach::Inside));
                    let slipped = level
                        .slipped_reads
                        .expect("a fitted class counted its reads");
                    assert!(
                        (0.0..=evidence.spanning_reads() as f64).contains(&slipped),
                        "{slipped}"
                    );
                    let shares = fit.shares_provenance[0].expect("group 0 has numbers");
                    assert_eq!(shares.slipped_reads, Some(slipped));
                    for share in [shares.shorter_share, shares.fall_off] {
                        assert!(matches!(share.source, ShareSource::PenalisedCurve(_)));
                        assert_eq!(share.curve, None);
                        assert_eq!(share.reach, Some(CurveReach::Inside));
                    }
                    let errors = fit
                        .standard_errors
                        .as_ref()
                        .expect("production asks for errors");
                    for error in every_error_of(errors) {
                        if let StratumError::Estimated(value) = error {
                            assert!(value.is_finite() && value > 0.0, "{error:?}");
                        }
                    }
                }
                StratumOutcome::Derived(_) => panic!("the penalised fit derives no class"),
                StratumOutcome::Refused { reason, .. } => {
                    panic!("refused only for no spanning reads, not {reason:?}")
                }
            }
        }
        assert!(matches!(outcomes[1], StratumOutcome::Fitted(_)));
        assert!(matches!(outcomes[2], StratumOutcome::Refused { .. }));
        assert!(matches!(outcomes[3], StratumOutcome::Refused { .. }));

        // The one-base classes, fitted on their own in input order, give the same numbers.
        let one_base = [strata[1].clone(), strata[2].clone(), strata[5].clone()];
        let alone = fit_period_by_penalised_counting(&one_base, &excess, &config);
        for (outcome, fitted) in [&outcomes[1], &outcomes[2], &outcomes[5]]
            .into_iter()
            .zip(&alone.classes)
        {
            match (outcome, fitted) {
                (StratumOutcome::Fitted(fit), Some(fitted)) => {
                    assert_eq!(fit.slippage[0], Some(fitted.slippage));
                    assert_eq!(fit.allele_length_variation, fitted.allele_length_variation);
                    let level = fit.level_provenance[0].expect("group 0 has numbers");
                    assert_eq!(level.slipped_reads, fitted.expected_slipped_reads);
                }
                (StratumOutcome::Refused { .. }, None) => {}
                other => panic!("the motif length's own fit disagrees: {other:?}"),
            }
        }
    }

    /// **At one sample no class has a concentration**, and its error says it was not placed: no
    /// locus holds two samples, so how much loci differ cannot be told apart from the sample's
    /// own inbreeding (spec §2.5). Three one-base classes of 40 loci at 10 reads.
    #[test]
    fn at_one_sample_no_class_has_a_concentration() {
        let strata: Vec<StratumEvidence> = (8..=10)
            .map(|repeats| a_drawn_class(1, repeats, 40, 1, 10, 20 + repeats))
            .collect();
        let outcomes = fit_repeat_tracts(&strata, &[0.3], &a_narrow_production_fit());
        for outcome in &outcomes {
            let StratumOutcome::Fitted(fit) = outcome else {
                panic!("every class has reads: {outcome:?}");
            };
            assert_eq!(fit.concentration, None);
            let errors = fit
                .standard_errors
                .as_ref()
                .expect("production asks for errors");
            assert_eq!(errors.concentration, StratumError::NotPlaced);
            for error in every_error_of(errors) {
                assert_ne!(error, StratumError::Estimated(0.0));
            }
        }
    }

    /// **A class with fewer than 8 loci with reads has every error withheld, and a rare share has
    /// its own withheld** (spec §2.6). Two-base runs at 8, 9 and 10 repeats with 40, 5 and 40
    /// loci, 12 samples at 8 reads. The thin class's every number is `TooFewLoci` (its
    /// concentration `NotPlaced` instead if the fit could not place it); on the rich classes a
    /// share below 3 loci's worth, 3/40, is `TooFewLoci` and a larger one is not, and some number
    /// carries an error.
    #[test]
    fn a_thin_class_and_a_rare_share_have_their_errors_withheld() {
        let strata = vec![
            a_drawn_class(2, 8, 40, 12, 8, 31),
            a_drawn_class(2, 9, 5, 12, 8, 32),
            a_drawn_class(2, 10, 40, 12, 8, 33),
        ];
        let outcomes = fit_repeat_tracts(&strata, &[0.3; 12], &a_narrow_production_fit());
        let fits: Vec<&StratumFit> = outcomes
            .iter()
            .map(|outcome| match outcome {
                StratumOutcome::Fitted(fit) => fit.as_ref(),
                other => panic!("every class has reads: {other:?}"),
            })
            .collect();

        let thin = fits[1];
        assert_eq!(thin.tracts_of_its_own, 5);
        let errors = thin
            .standard_errors
            .as_ref()
            .expect("production asks for errors");
        let group = errors.slippage[0].expect("group 0 has numbers");
        for error in [group.level, group.shorter_share, group.fall_off] {
            assert_eq!(error, StratumError::TooFewLoci);
        }
        assert!(
            errors
                .allele_length_variation
                .iter()
                .all(|error| *error == StratumError::TooFewLoci)
        );
        let expected_concentration_error = if thin.concentration.is_some() {
            StratumError::TooFewLoci
        } else {
            StratumError::NotPlaced
        };
        assert_eq!(errors.concentration, expected_concentration_error);

        let mut rare_shares = 0;
        for rich in [fits[0], fits[2]] {
            let errors = rich
                .standard_errors
                .as_ref()
                .expect("production asks for errors");
            assert!(errors.any(), "a rich class has some error: {errors:?}");
            for (share, error) in rich
                .allele_length_variation
                .iter()
                .zip(&errors.allele_length_variation)
            {
                if share * 40.0 < LOCI_FOR_A_SHARE_ERROR {
                    rare_shares += 1;
                    assert_eq!(*error, StratumError::TooFewLoci, "share {share}");
                } else {
                    assert_ne!(*error, StratumError::TooFewLoci, "share {share}");
                }
            }
        }
        assert!(rare_shares > 0, "the window's outer lengths are rare");
    }

    /// **A class's slipped reads are the fit's expected count, not the level times the reads.**
    /// One-base runs at 8 to 12 repeats, 40 loci each, drawn at a slip rate of 0.01 per repeat
    /// above 5, except the class at 10 repeats: 3 loci, drawn at 0.08 where its neighbours' line
    /// passes through 0.05. Twelve samples at 8 reads.
    ///
    /// - On a well-measured class (9 repeats) the field is the expected count at the returned
    ///   numbers, to within the one round's lag the fit's tolerance allows.
    /// - On the thin class the penalised level is pulled towards its neighbours, while the
    ///   expected count stays with the class's own reads, so the two numbers differ.
    #[test]
    fn slipped_reads_are_the_fits_expected_count() {
        let mut strata: Vec<StratumEvidence> = (8..=12)
            .map(|repeats| a_drawn_class(1, repeats, 40, 12, 8, 40 + repeats))
            .collect();
        let mut thin = draw_stratum(
            Slippage {
                level: 0.08,
                shorter_share: 0.8,
                fall_off: 0.25,
            },
            &spectrum_of(3),
            2.0,
            0.3,
            3,
            12,
            8,
            1,
            61,
        );
        thin.stratum = Stratum {
            period: 1,
            reference_repeats: 10,
        };
        strata[2] = thin;
        let excess = [0.3; 12];
        let config = a_narrow_production_fit();
        let outcomes = fit_repeat_tracts(&strata, &excess, &config);
        let fits: Vec<&StratumFit> = outcomes
            .iter()
            .map(|outcome| match outcome {
                StratumOutcome::Fitted(fit) => fit.as_ref(),
                other => panic!("every class has reads: {other:?}"),
            })
            .collect();
        let slipped_reads_of = |fit: &StratumFit| {
            let level = fit.level_provenance[0].expect("group 0 has numbers");
            let shares = fit.shares_provenance[0].expect("group 0 has numbers");
            assert_eq!(level.slipped_reads, shares.slipped_reads);
            level
                .slipped_reads
                .expect("a fitted class counted its reads")
        };

        let rich = fits[1];
        let counts = expected_counts(
            &strata[1],
            rich.slippage[0].expect("group 0 has numbers"),
            &rich.allele_length_variation,
            rich.concentration
                .expect("twelve samples place the concentration"),
            &excess,
            config.allele_span,
        );
        let field = slipped_reads_of(rich);
        eprintln!(
            "rich class: expected slipped reads {field:.4}, recounted at the returned numbers {:.4}",
            counts.slipped
        );
        assert!(
            (field - counts.slipped).abs() <= 1e-3 * counts.slipped,
            "field {field}, expected count at the returned numbers {}",
            counts.slipped
        );

        let thin = fits[2];
        let field = slipped_reads_of(thin);
        let level = thin.slippage[0].expect("group 0 has numbers").level;
        let level_times_reads = level * thin.reads_crossing as f64;
        eprintln!(
            "thin class: {} reads crossing, level {level:.4}; expected slipped reads {field:.2}, \
             level × reads crossing {level_times_reads:.2}",
            thin.reads_crossing
        );
        assert!(
            field > 1.2 * level_times_reads,
            "the class's own reads slipped more than the curve's level says: {field} against \
             {level_times_reads}"
        );
    }

    /// A class as the fit returns it, for the conversion's own tests: 5% slip rate, every error
    /// given, five allele lengths.
    fn a_fitted_class(concentration: Option<f64>) -> PenalisedClass {
        PenalisedClass {
            stratum: Stratum {
                period: 2,
                reference_repeats: 9,
            },
            slippage: Slippage {
                level: 0.05,
                shorter_share: 0.8,
                fall_off: 0.25,
            },
            allele_length_variation: vec![0.01, 0.2, 0.5, 0.25, 0.04],
            concentration,
            between_locus_share: concentration.map(|c| 1.0 / (1.0 + c)),
            working_concentration: concentration.unwrap_or(1.0),
            expected_slipped_reads: Some(50.0),
            errors: Some(PenalisedClassErrors {
                level: Some(0.01),
                shorter_share: Some(0.05),
                fall_off: Some(0.03),
                concentration: Some(0.5),
                between_locus_share: Some(0.01),
                allele_length_variation: vec![0.002, 0.03, 0.04, 0.03, 0.01],
            }),
        }
    }

    /// **The error mapping follows its table, rule by rule** (arch `str_penalised_fit.md` §1.4),
    /// on a class built by hand.
    #[test]
    fn the_error_mapping_follows_its_table() {
        use StratumError::{Estimated, NotPlaced, TooFewLoci};
        // 100 loci: every error given is placed except the rarest share's (1% of 100 loci is 1
        // locus's worth, under 3). The last share, 4% of 100 = 4 loci's worth, keeps its error.
        let class = a_fitted_class(Some(2.0));
        let errors = errors_of_a_class(&class, Some(2.0), 100);
        assert_eq!(
            errors.slippage,
            vec![Some(SlippageErrors {
                level: Estimated(0.01),
                shorter_share: Estimated(0.05),
                fall_off: Estimated(0.03),
            })]
        );
        assert_eq!(
            errors.allele_length_variation,
            vec![
                TooFewLoci,
                Estimated(0.03),
                Estimated(0.04),
                Estimated(0.03),
                Estimated(0.01)
            ]
        );
        assert_eq!(errors.concentration, Estimated(0.5));

        // Seven loci: every number withheld, the concentration too.
        let errors = errors_of_a_class(&class, Some(2.0), 7);
        assert!(every_error_of(&errors).iter().all(|e| *e == TooFewLoci));
        // An absent concentration is not placed, before any other rule.
        let errors = errors_of_a_class(&class, None, 7);
        assert_eq!(errors.concentration, NotPlaced);
        let errors = errors_of_a_class(&class, None, 100);
        assert_eq!(errors.concentration, NotPlaced);

        // Zero, not a number, infinite and absent errors are not placed; so is an error wider
        // than 3 on the fitted scale: 0.2 on a slip rate of 5% is 0.2 / 0.0475 ≈ 4.2 in logit,
        // and 7 on a concentration of 2 is 3.5 in log.
        let mut bad = a_fitted_class(Some(2.0));
        bad.errors = Some(PenalisedClassErrors {
            level: Some(0.2),
            shorter_share: Some(0.0),
            fall_off: Some(f64::NAN),
            concentration: Some(7.0),
            between_locus_share: None,
            allele_length_variation: vec![0.002, f64::INFINITY, 0.0],
        });
        let errors = errors_of_a_class(&bad, Some(2.0), 100);
        assert!(
            every_error_of(&errors)
                .iter()
                .all(|e| matches!(e, NotPlaced | TooFewLoci)),
            "{errors:?}"
        );
        assert_eq!(errors.concentration, NotPlaced);
        // Just inside: 0.14 on a slip rate of 5% is about 2.9 in logit.
        bad.errors = Some(PenalisedClassErrors {
            level: Some(0.14),
            ..PenalisedClassErrors::default()
        });
        let errors = errors_of_a_class(&bad, Some(2.0), 100);
        assert_eq!(
            errors.slippage[0].map(|group| group.level),
            Some(Estimated(0.14))
        );

        // A class the fit gave no errors at all: not placed, the withheld aside.
        let mut none = a_fitted_class(Some(2.0));
        none.errors = None;
        let errors = errors_of_a_class(&none, Some(2.0), 100);
        assert_eq!(errors.allele_length_variation[0], TooFewLoci);
        assert!(
            every_error_of(&errors)
                .iter()
                .all(|e| matches!(e, NotPlaced | TooFewLoci))
        );
    }

    /// One curve over 8 to 10 repeats, at `strength`, eased at `eased_at`.
    fn a_curve(strength: f64, eased_at: Vec<u64>) -> ParameterCurve {
        ParameterCurve {
            first_repeats: 8,
            values: vec![0.05; 3],
            smoothing: strength,
            eased_at,
        }
    }

    /// **A concentration at an end of its search is reported absent, with its error not
    /// placed**, and one inside is kept; and each number's point records its own curve's strength
    /// and whether that curve was eased at the class's repeat count (9).
    #[test]
    fn a_concentration_at_an_end_of_its_search_is_absent() {
        let curves = SmoothCurves {
            level: a_curve(3.0, vec![9]),
            shorter_share: a_curve(40.0, Vec::new()),
            fall_off: a_curve(500.0, vec![9]),
            between_locus_share: Some(a_curve(7.0, Vec::new())),
            allele_length_variation: Vec::new(),
        };
        let evidence = a_drawn_class(2, 9, 20, 4, 5, 41);
        for (estimate, kept) in [
            (1e-3, false),
            (1.01e-3, false),
            (1e4, false),
            (9_901.0, false),
            (2.0, true),
        ] {
            let outcome = outcome_of_a_class(
                &evidence,
                Some(a_fitted_class(Some(estimate))),
                Some(&curves),
                true,
            );
            let StratumOutcome::Fitted(fit) = outcome else {
                panic!("a fitted class is fitted");
            };
            let errors = fit.standard_errors.expect("errors were asked for");
            if kept {
                assert_eq!(fit.concentration, Some(estimate));
                assert_eq!(errors.concentration, StratumError::Estimated(0.5));
            } else {
                assert_eq!(fit.concentration, None, "{estimate}");
                assert_eq!(errors.concentration, StratumError::NotPlaced, "{estimate}");
            }
            let level = fit.level_provenance[0].expect("one group");
            assert_eq!(
                level.source,
                LevelSource::PenalisedCurve(PenalisedPoint {
                    strength: 3.0,
                    eased_here: true
                })
            );
            let shares = fit.shares_provenance[0].expect("one group");
            assert_eq!(
                shares.shorter_share.source,
                ShareSource::PenalisedCurve(PenalisedPoint {
                    strength: 40.0,
                    eased_here: false
                })
            );
            assert_eq!(
                shares.fall_off.source,
                ShareSource::PenalisedCurve(PenalisedPoint {
                    strength: 500.0,
                    eased_here: true
                })
            );
        }
        // A fit asked for no errors carries none.
        let StratumOutcome::Fitted(fit) = outcome_of_a_class(
            &evidence,
            Some(a_fitted_class(Some(2.0))),
            Some(&curves),
            false,
        ) else {
            panic!("a fitted class is fitted");
        };
        assert_eq!(fit.standard_errors, None);
    }

    /// **A drawn cohort whose loci all carry the same spread of lengths puts F at the bottom of
    /// its search**, and the class reports no concentration: 40 samples at 20 reads over 40 loci
    /// drawn at a concentration of 10⁶, so the pairs of samples at a locus share nothing beyond
    /// the motif length's spectrum.
    #[test]
    fn a_drawn_cohort_with_f_at_the_end_of_its_search_reports_none() {
        let strata: Vec<StratumEvidence> = (8..=9)
            .map(|repeats| {
                let truth = Slippage {
                    level: 0.03,
                    shorter_share: 0.8,
                    fall_off: 0.25,
                };
                let mut evidence = draw_stratum(
                    truth,
                    &spectrum_of(3),
                    1e6,
                    0.0,
                    40,
                    40,
                    20,
                    1,
                    50 + repeats,
                );
                evidence.stratum = Stratum {
                    period: 1,
                    reference_repeats: repeats,
                };
                evidence
            })
            .collect();
        let config = a_narrow_production_fit();
        // The fit itself stopped at the end: the class's estimate is the search's top.
        let raw = fit_period_by_penalised_counting(&strata, &[0.0; 40], &config);
        for fitted in &raw.classes {
            let estimate = fitted
                .as_ref()
                .and_then(|class| class.concentration)
                .expect("every locus holds two samples, so the fit reports an estimate");
            assert!(estimate >= 1e4 / 1.01, "{estimate}");
        }
        let outcomes = fit_repeat_tracts(&strata, &[0.0; 40], &config);
        for outcome in &outcomes {
            let StratumOutcome::Fitted(fit) = outcome else {
                panic!("every class has reads");
            };
            assert_eq!(fit.concentration, None);
            let errors = fit
                .standard_errors
                .as_ref()
                .expect("production asks for errors");
            assert_eq!(errors.concentration, StratumError::NotPlaced);
        }
    }

    /// **Evidence that declares two slippage groups is refused**, by name: the fit has one
    /// stutter for every library.
    #[test]
    #[should_panic(expected = "declares 2 slippage groups")]
    fn evidence_with_two_slippage_groups_is_refused() {
        let mut evidence = a_drawn_class(1, 8, 5, 2, 3, 61);
        evidence.groups = 2;
        let _ = fit_repeat_tracts(&[evidence], &[0.3; 2], &a_narrow_production_fit());
    }

    /// **Setting the window rescales the first floor to one draw a class**, and at production's
    /// own window it changes nothing, to the bit.
    #[test]
    fn a_new_window_keeps_one_draw_a_class() {
        let production = PenalisedFitConfig::production();
        assert_eq!(production.clone().with_allele_span(ALLELE_SPAN), production);
        let narrow = PenalisedFitConfig::production().with_allele_span(4);
        assert_eq!(narrow.allele_span, 4);
        assert_eq!(narrow.allele_pseudo_draws, 1.0 / 9.0);
        assert_eq!(
            PenalisedFitConfig {
                allele_pseudo_draws: 0.0,
                ..PenalisedFitConfig::prototype()
            }
            .with_allele_span(9)
            .allele_pseudo_draws,
            1.0 / 19.0
        );
    }
}
