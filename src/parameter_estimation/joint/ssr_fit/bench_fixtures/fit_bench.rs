//! **The library half of the repeat-tract fitting bench**: the interface every fitting method
//! implements and a held-out score that does not depend on the method
//! (`doc/devel/ng/research/str_fit_research_plan.md` §2, §5).
//!
//! The bench itself is `examples/ng_str_fit_bench.rs`: it draws a cohort with
//! `cohort_draw`, runs each method through this interface, and writes one row
//! per method, class and number.
//!
//! # Adding a method
//!
//! A method implements [`PeriodFit`]: it is handed **every class of one motif length at once** —
//! each class's evidence, as the census builds it, and each sample's homozygote excess — and returns
//! one [`ClassEstimate`] per class, in the same order. Taking a whole motif length lets a method fit
//! all its classes together (penalised EM fits every number as a curve across repeat count); a
//! method that fits class by class just loops.
//!
//! A [`ClassEstimate`] holds the three stutter numbers, the allele-length shares over
//! `-allele_span ..= +allele_span` (the outermost two meaning "this far or further"), the
//! concentration `c` of the Dirichlet each locus's frequencies are drawn from, and a standard error
//! for any of them the method can give. **Absent is a legitimate answer**: a number the method
//! cannot estimate is `None`, never a fitted zero. Each method's answers are reported under its
//! own name.
//!
//! Two methods are wired in: [`ThePenalisedFit`] (`D`, and `D_tract_prior`, the fit
//! `estimate-parameters` runs) and [`TheCountingFit`] (`B`, expectation maximisation per class
//! with the spread read from pairs of samples, the baseline the held-out comparisons are made
//! against). **Two methods went with the climb on 2026-10-16** (`arch/str_penalised_fit.md`
//! §4.1): `A`, the climb itself, and the `_curves` variant of any method that did not draw its own
//! curves, which applied the climb's curves to its answers (`B_curves`); the counting fit's
//! answers, which were reported as `B_raw` beside it, are now reported as `B`. Their results stay
//! in the research reports.
//!
//! # The held-out score
//!
//! Each class's loci are split into folds, **whole loci**, because samples at one locus share its
//! alleles. A method is fitted on all folds but one and scored on the loci of the one left out, by
//! the log of the probability of every sample's reads at the locus under the class's fitted
//! numbers ([`score_locus`]):
//!
//! `P(reads at the locus) = E over the locus's frequencies p ~ Dirichlet(c × shares) of
//! Π over samples Σ over genotypes P(genotype | p, homozygote excess) P(the sample's reads | genotype)`.
//!
//! The expectation over `p` has no closed form beyond a couple of samples, so:
//!
//! - **up to [`HeldOutSettings::exact_up_to_samples`] samples with reads (2 by default) it is
//!   summed exactly**: each sample's genotype probability is a polynomial in `p`, and the
//!   Dirichlet's moments are known;
//! - **beyond that it is a Monte Carlo average with importance sampling**: frequencies are drawn
//!   from Dirichlets placed where the locus's reads put them — the class's Dirichlet with the
//!   locus's expected allele counts added, at full, half and quarter weight — mixed with one draw in
//!   ten from the class's own Dirichlet, so that no draw can carry an unbounded weight. Draws are
//!   added until the Monte Carlo standard error of the log-likelihood is below
//!   [`HeldOutSettings::target_standard_error`] (0.05 by default) or the draw limit is reached; the
//!   error is reported either way.
//!
//! **The climb's 256-point average was not used**: at 13 allele lengths it is several
//! log-likelihood units a locus wrong (`fit_precision_quadrature_average_2026-10-01.md`).
//!
//! Every method's score at a locus uses the same random numbers, so differences between methods
//! carry less Monte Carlo noise than either score alone.

use std::time::Instant;

use super::cohort_draw::ClassTruth;
use super::*;
use crate::parameter_estimation::joint::ssr_fit::expected_counts::{
    DistinctRows, ScoredDistinctRows,
};
use crate::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, PenalisedPeriodFit, fit_period_by_penalised_counting,
};

// ---------------------------------------------------------------------
// What a method returns
// ---------------------------------------------------------------------

/// **One class's numbers as a fitting method returns them.** Every number may be absent.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassEstimate {
    pub stratum: Stratum,
    /// The class's stutter at multiplier one.
    pub slippage: Option<Slippage>,
    /// The allele-length shares, `-allele_span ..= +allele_span`.
    pub allele_length_variation: Option<Vec<f64>>,
    /// The concentration `c`; the between-locus share is `1 / (1 + c)`.
    pub concentration: Option<f64>,
    pub errors: ClassErrors,
}

/// **The standard errors a method gives**, each on its number's own scale; `None` where it gives
/// none.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ClassErrors {
    pub level: Option<f64>,
    pub shorter_share: Option<f64>,
    pub fall_off: Option<f64>,
    pub concentration: Option<f64>,
    /// Indexed as the shares are; empty when the method gives none.
    pub allele_length_variation: Vec<Option<f64>>,
}

/// The numbers the held-out score needs, all present.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassNumbers {
    pub slippage: Slippage,
    pub allele_length_variation: Vec<f64>,
    pub concentration: f64,
}

impl ClassEstimate {
    /// A class the method gave no answer for.
    #[must_use]
    pub fn absent(evidence: &StratumEvidence) -> Self {
        Self {
            stratum: evidence.stratum,
            slippage: None,
            allele_length_variation: None,
            concentration: None,
            errors: ClassErrors::default(),
        }
    }

    /// Every number the held-out score needs, or `None` when one is absent.
    #[must_use]
    pub fn numbers(&self) -> Option<ClassNumbers> {
        Some(ClassNumbers {
            slippage: self.slippage?,
            allele_length_variation: self.allele_length_variation.clone()?,
            concentration: self.concentration?,
        })
    }
}

/// One class's numbers as the held-out score uses them, and whether they came from the caller's
/// fallback rung.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoredNumbers {
    pub numbers: ClassNumbers,
    /// **The class had no allele-length shares or no concentration of its own**, so it is scored
    /// at the stated concentration — at a flat shape by default, or at its own shares under
    /// [`ScoringFallback::keep_own_shares`].
    pub on_the_fallback: bool,
}

/// **The numbers each class is scored at, with a fallback where a class lacks its own.**
///
/// Calling seeds a tract's genotype prior from the tract ladder
/// ([`StratumFits::allele_length_variation_at`](crate::parameter_estimation::joint::stratum_fits::StratumFits::allele_length_variation_at)):
/// the class's own shares and concentration where it has both; the class's own shares at the
/// stated concentration where it has shares and no concentration; and otherwise a flat shape at
/// the stated concentration. The stated concentration is the median of the concentrations the
/// run's classes fitted, or
/// [`STATED_FLAT_CONCENTRATION`](crate::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION)
/// where none fitted one. (The ladder's pooled rung, a period's pooled tracts, is reached only when
/// a run asks for it, and the bench does not.)
///
/// **The default scoring still models the ladder as it was before the second rung was added**
/// (`doc/devel/ng/spec/str_fit_to_calling.md` §3.2): a class with shares but no concentration —
/// the counting fit at one sample — is scored at a flat shape and the stated concentration, not
/// skipped, and marked. That keeps the bench's numbers comparable with the research results
/// already recorded. [`ScoringFallback::keep_own_shares`] scores such a class the way calling now
/// seeds it.
///
/// The median is over `estimates`, the classes of one motif length fitted together: the bench's
/// stand-in for a run. A class with no stutter (no answer and no curve) is still not scored.
///
/// `fallback` can depart from calling, for measurement only ([`ScoringFallback`]).
#[must_use]
pub fn numbers_for_scoring(
    estimates: &[ClassEstimate],
    fallback: ScoringFallback,
) -> Vec<Option<ScoredNumbers>> {
    let mut fitted: Vec<f64> = estimates
        .iter()
        .filter(|estimate| estimate.allele_length_variation.is_some())
        .filter_map(|estimate| estimate.concentration)
        .collect();
    fitted.sort_by(f64::total_cmp);
    let median = match fitted.len() {
        0 => crate::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION,
        count if count % 2 == 1 => fitted[count / 2],
        count => 0.5 * (fitted[count / 2 - 1] + fitted[count / 2]),
    };
    let stated = fallback.stated_concentration.unwrap_or(median);
    estimates
        .iter()
        .map(|estimate| {
            let slippage = estimate.slippage?;
            Some(
                match (&estimate.allele_length_variation, estimate.concentration) {
                    (Some(shares), Some(concentration)) => ScoredNumbers {
                        numbers: ClassNumbers {
                            slippage,
                            allele_length_variation: shares.clone(),
                            concentration,
                        },
                        on_the_fallback: false,
                    },
                    (Some(shares), None) if fallback.keep_own_shares => ScoredNumbers {
                        numbers: ClassNumbers {
                            slippage,
                            allele_length_variation: shares.clone(),
                            concentration: stated,
                        },
                        on_the_fallback: true,
                    },
                    (shares, _) => {
                        let lengths = shares
                            .as_ref()
                            .map_or(2 * ALLELE_SPAN as usize + 1, Vec::len);
                        ScoredNumbers {
                            numbers: ClassNumbers {
                                slippage,
                                allele_length_variation: vec![1.0 / lengths as f64; lengths],
                                concentration: stated,
                            },
                            on_the_fallback: true,
                        }
                    }
                },
            )
        })
        .collect()
}

/// **How [`numbers_for_scoring`] fills a class that lacks its own numbers.** The default is what
/// calling did before its ladder gained the rung that keeps a class's own shares; the two
/// departures exist to measure what a fallback choice costs, and calling is not changed by them.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ScoringFallback {
    /// For a class with its own shares but no concentration, keep its shares and take only the
    /// stated concentration. This is what calling now does; the default replaces both by a flat
    /// shape, which is what calling did when the recorded research results were measured.
    pub keep_own_shares: bool,
    /// The stated concentration in place of the period's median (or
    /// [`STATED_FLAT_CONCENTRATION`](crate::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION)).
    pub stated_concentration: Option<f64>,
}

impl ClassNumbers {
    /// The truth as a fit at `allele_span` sees it.
    #[must_use]
    pub fn of_truth(truth: &ClassTruth, allele_span: i32) -> Self {
        Self {
            slippage: truth.slippage,
            allele_length_variation: truth.shares_in_window(allele_span),
            concentration: truth.concentration,
        }
    }
}

/// **A way of fitting every class of one motif length.** See the module documentation.
pub trait PeriodFit: Sync {
    /// The method's name in the bench's output, e.g. `B`.
    fn name(&self) -> String;

    /// Fit `classes`, every class of one motif length, each sample's homozygote excess held at
    /// `homozygote_excess`, the shares placed over `-allele_span ..= +allele_span`. One estimate per
    /// class, in the order of `classes`.
    fn fit_period(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
    ) -> Vec<ClassEstimate>;
}

// ---------------------------------------------------------------------
// The two methods wired in
// ---------------------------------------------------------------------

/// **Method B: the counting fit**, class by class
/// (`fit_class_by_counting_with_pairs`
/// with `pairs`, the within-sample version without), as `examples/ng_str_counting_fit.rs` runs it.
/// A class with fewer loci carrying reads than
/// `COUNTING_FIT_FLOOR` (8, the climb's
/// refusal floor) is left absent, as that example leaves it. Its concentration is absent where the
/// fit cannot place it.
///
/// **Its standard errors** come from `counting_fit_errors`, by the method `errors` names.
#[derive(Debug, Clone, Copy)]
pub struct TheCountingFit {
    /// Read the spread off pairs of samples at a locus (what an inbred cohort needs).
    pub pairs: bool,
    pub max_rounds: usize,
    pub tolerance: f64,
    /// Rounds of spread-then-counting, at most, when `pairs`.
    pub alternations: usize,
    pub errors: CountingFitErrorMethod,
    /// The shares' prior weight, in tracts
    /// (`SHARE_PRIOR_TRACTS` by default; zero is
    /// the fit before the prior). Read only with `pairs`.
    pub share_prior_tracts: f64,
}

/// Which standard errors [`TheCountingFit`] reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountingFitErrorMethod {
    None,
    /// The composite likelihood's sandwich formula.
    Sandwich,
    /// The tracts resampled this many times, each draw refitted.
    Bootstrap(usize),
}

/// **Method D: penalised counting**, every class of the motif length at once, each number a
/// penalised curve across repeat count ([`fit_period_by_penalised_counting`]), the spread read off
/// pairs of samples as method B reads it. It fits every
/// class with reads (no refusal floor: a thin class takes what its counts and its neighbours say),
/// and gives standard errors for every number it reports (`PenalisedClassErrors`).
///
/// **The concentration is absent when no locus has two samples with reads**, as the project
/// reports a number it cannot estimate. With `report_working_concentration` it reports the value
/// the fit worked at instead, so that the held-out score, which needs every number, can be computed
/// at one sample; that variant is `D_working_F` and its concentration is not an estimate.
#[derive(Debug, Clone)]
pub struct ThePenalisedFit {
    pub config: PenalisedFitConfig,
    pub report_working_concentration: bool,
}

/// **Method `D`: production's settings without the second floor under the allele lengths**
/// (`allele_prior_tracts` zero), as the fit was first measured. [`ThePenalisedFit::production`] is
/// `D_tract_prior`.
impl Default for ThePenalisedFit {
    fn default() -> Self {
        Self {
            config: PenalisedFitConfig {
                allele_prior_tracts: 0.0,
                ..PenalisedFitConfig::production()
            },
            report_working_concentration: false,
        }
    }
}

impl ThePenalisedFit {
    /// **Method `D_tract_prior`: the fit as `estimate-parameters` runs it**
    /// ([`PenalisedFitConfig::production`]).
    #[must_use]
    pub fn production() -> Self {
        Self {
            config: PenalisedFitConfig::production(),
            report_working_concentration: false,
        }
    }

    /// **The whole fit of one motif length**: every class's numbers and errors, the curves with
    /// their strengths, and how many rounds each part of the loop ran. [`PeriodFit::fit_period`]
    /// keeps only the classes' numbers of this same fit; the bench's reference fixture records all
    /// of it, so a later change to the fit can be checked against every number and every round
    /// count.
    #[must_use]
    pub fn fit_period_in_full(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
    ) -> PenalisedPeriodFit {
        // The bench chooses the window for every method alike, with the floor that goes with it;
        // the rest is this method's own.
        let config = self.config.clone().with_allele_span(allele_span);
        fit_period_by_penalised_counting(classes, homozygote_excess, &config)
    }
}

impl PeriodFit for ThePenalisedFit {
    fn name(&self) -> String {
        if self.report_working_concentration {
            "D_working_F".to_owned()
        } else if self.config.allele_prior_tracts > 0.0 {
            "D_tract_prior".to_owned()
        } else {
            "D".to_owned()
        }
    }

    fn fit_period(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
    ) -> Vec<ClassEstimate> {
        let fitted = self.fit_period_in_full(classes, homozygote_excess, allele_span);
        classes
            .iter()
            .zip(fitted.classes)
            .map(|(evidence, class)| match class {
                None => ClassEstimate::absent(evidence),
                Some(class) => ClassEstimate {
                    stratum: evidence.stratum,
                    slippage: Some(class.slippage),
                    allele_length_variation: Some(class.allele_length_variation),
                    concentration: if self.report_working_concentration {
                        Some(class.working_concentration)
                    } else {
                        class.concentration
                    },
                    errors: class
                        .errors
                        .map_or_else(ClassErrors::default, |errors| ClassErrors {
                            level: errors.level,
                            shorter_share: errors.shorter_share,
                            fall_off: errors.fall_off,
                            concentration: errors.concentration,
                            allele_length_variation: errors
                                .allele_length_variation
                                .into_iter()
                                .map(Some)
                                .collect(),
                        }),
                },
            })
            .collect()
    }
}

// ---------------------------------------------------------------------
// Folds
// ---------------------------------------------------------------------

/// How the held-out score is taken.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeldOutSettings {
    /// How many folds each class's loci are split into; 0 switches the held-out score off.
    pub folds: usize,
    /// At most this many held-out loci with reads are scored per class, spread evenly over the
    /// folds. Scoring a locus of 2,000 samples costs about a tenth of a second.
    pub loci_scored_per_class: usize,
    /// Monte Carlo draws are added until the log-likelihood's standard error is below this.
    pub target_standard_error: f64,
    /// ... or this many draws have been taken.
    pub max_draws: usize,
    /// A locus with at most this many samples with reads is summed exactly.
    pub exact_up_to_samples: usize,
    /// Which sampler a locus with more samples gets.
    pub sampler: Sampler,
    pub seed: u64,
}

/// Which Monte Carlo sampler integrates a locus's frequencies out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Sampler {
    /// The Dirichlet importance sampler, then the sequential one where the first falls short of
    /// the error target ([`sampled_integral`]).
    #[default]
    Automatic,
    /// The Dirichlet importance sampler alone.
    DirichletOnly,
    /// The sequential sampler alone.
    SequentialOnly,
}

impl Default for HeldOutSettings {
    fn default() -> Self {
        Self {
            folds: 5,
            loci_scored_per_class: 100,
            target_standard_error: 0.05,
            max_draws: 65_536,
            exact_up_to_samples: 2,
            sampler: Sampler::Automatic,
            seed: 0x5EED,
        }
    }
}

/// **Which fold each locus of a class is in, and whether it is scored**: the loci in a seeded
/// random order, fold = position modulo `folds`, scored when the position is below the cap among
/// loci with reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FoldPlan {
    /// Per locus, its fold.
    pub fold_of_locus: Vec<usize>,
    /// Per locus, whether its held-out score is taken.
    pub scored: Vec<bool>,
}

impl FoldPlan {
    /// The plan for `evidence`, `class_seed` distinguishing classes.
    #[must_use]
    pub fn of(evidence: &StratumEvidence, settings: &HeldOutSettings, class_seed: u64) -> Self {
        let loci = evidence.tracts.len();
        let folds = settings.folds.max(1);
        let mut order: Vec<usize> = (0..loci).collect();
        let mut rng = Rng(settings.seed ^ class_seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        // Fisher–Yates.
        for index in (1..loci).rev() {
            let other = (rng.next_u64() % (index as u64 + 1)) as usize;
            order.swap(index, other);
        }
        let mut fold_of_locus = vec![0; loci];
        let mut scored = vec![false; loci];
        let mut taken = 0;
        for (position, locus) in order.into_iter().enumerate() {
            fold_of_locus[locus] = position % folds;
            if taken < settings.loci_scored_per_class && !evidence.tracts[locus].samples.is_empty()
            {
                scored[locus] = true;
                taken += 1;
            }
        }
        Self {
            fold_of_locus,
            scored,
        }
    }
}

/// **A class's seed in a cohort read from psps**, where one repeat count recurs at several motif
/// lengths: the motif length is mixed in, so that the one-base and two-base classes at ten repeats
/// are not shuffled alike. It seeds the class's folds and its loci's Monte Carlo draws.
///
/// The drawn bench, one motif length a run, seeds a class by its repeat count alone, and keeps
/// doing so: changing that would move every recorded drawn result.
#[must_use]
pub fn class_seed_across_motif_lengths(stratum: Stratum) -> u64 {
    (u64::from(stratum.period) << 16) ^ stratum.reference_repeats
}

/// **Every class's folds, for a cohort of several motif lengths** — the bench's psp mode.
///
/// Each class is folded on its own by [`FoldPlan::of`], **whole loci**: a locus's every sample
/// lands in one fold, since samples at one locus share its alleles, and every locus is in exactly
/// one fold, so over the folds each locus is held out once. The order is fixed by
/// `settings.seed` and [`class_seed_across_motif_lengths`]; the same seed folds the same way.
#[must_use]
pub fn fold_plans_across_motif_lengths(
    classes: &[StratumEvidence],
    settings: &HeldOutSettings,
) -> Vec<FoldPlan> {
    classes
        .iter()
        .map(|evidence| {
            FoldPlan::of(
                evidence,
                settings,
                class_seed_across_motif_lengths(evidence.stratum),
            )
        })
        .collect()
}

/// **The class with fold `fold`'s loci taken out**, for fitting.
#[must_use]
pub fn without_fold(evidence: &StratumEvidence, plan: &FoldPlan, fold: usize) -> StratumEvidence {
    StratumEvidence {
        stratum: evidence.stratum,
        tracts: evidence
            .tracts
            .iter()
            .zip(&plan.fold_of_locus)
            .filter(|(_, in_fold)| **in_fold != fold)
            .map(|(tract, _)| tract.clone())
            .collect(),
        read_span: evidence.read_span,
        groups: evidence.groups,
        libraries: Arc::clone(&evidence.libraries),
        // **These five are the whole class's totals, held-out loci included**, copied as they
        // are. No fit and not the score reads them today, so nothing leaks. A fit that starts
        // reading one (the mismatch rate from `bases_compared`, say) must recount it over the
        // kept tracts here, or the held-out loci inform their own score.
        tracts_over_guard_threshold: evidence.tracts_over_guard_threshold,
        reads_reaching_not_crossing: evidence.reads_reaching_not_crossing,
        guard_reads: evidence.guard_reads,
        bases_compared: evidence.bases_compared,
        mismatching_bases: evidence.mismatching_bases,
    }
}

// ---------------------------------------------------------------------
// The held-out run of one method
// ---------------------------------------------------------------------

/// A locus's Monte Carlo seed, the same for every method so their scores share their draws.
#[must_use]
pub fn locus_seed(settings: &HeldOutSettings, class_seed: u64, locus: usize) -> u64 {
    settings.seed ^ (class_seed << 40) ^ (locus as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

/// **One motif length's classes, folded, with everything the held-out score reads.**
pub struct HeldOutClasses<'a> {
    classes: &'a [StratumEvidence],
    /// Per sample, indexed as the classes' evidence numbers samples.
    homozygote_excess: &'a [f64],
    plans: Vec<FoldPlan>,
    /// Per class, the seed its loci's Monte Carlo draws are taken from ([`locus_seed`]).
    class_seeds: Vec<u64>,
    scored: Vec<FoldedLocus>,
}

/// **One locus whose held-out score is taken**: its class, its index in the class, its fold, and
/// the truth's score there — `NaN` until [`HeldOutClasses::with_truth_scores`] gives it, and for
/// good without a truth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FoldedLocus {
    pub class: usize,
    pub locus: usize,
    pub fold: usize,
    pub truth_ln: f64,
}

impl<'a> HeldOutClasses<'a> {
    /// `classes` folded by `plans` (one a class), each class's loci drawn from `class_seeds`.
    ///
    /// # Panics
    ///
    /// When `plans` or `class_seeds` do not hold one entry a class.
    #[must_use]
    pub fn new(
        classes: &'a [StratumEvidence],
        homozygote_excess: &'a [f64],
        plans: Vec<FoldPlan>,
        class_seeds: Vec<u64>,
    ) -> Self {
        assert_eq!(plans.len(), classes.len(), "one fold plan a class");
        assert_eq!(class_seeds.len(), classes.len(), "one seed a class");
        let scored = plans
            .iter()
            .enumerate()
            .flat_map(|(class, plan)| {
                plan.scored
                    .iter()
                    .enumerate()
                    .filter(|(_, scored)| **scored)
                    .map(move |(locus, _)| FoldedLocus {
                        class,
                        locus,
                        fold: plan.fold_of_locus[locus],
                        truth_ln: f64::NAN,
                    })
            })
            .collect();
        Self {
            classes,
            homozygote_excess,
            plans,
            class_seeds,
            scored,
        }
    }

    /// Every locus whose held-out score is taken, class by class.
    #[must_use]
    pub fn scored(&self) -> &[FoldedLocus] {
        &self.scored
    }

    /// The truth's score at every scored locus, in the order of [`Self::scored`].
    ///
    /// # Panics
    ///
    /// When `truth_scores` does not hold one score a scored locus.
    #[must_use]
    pub fn with_truth_scores(mut self, truth_scores: &[f64]) -> Self {
        assert_eq!(truth_scores.len(), self.scored.len());
        for (scored, truth) in self.scored.iter_mut().zip(truth_scores) {
            scored.truth_ln = *truth;
        }
        self
    }

    /// One locus's held-out score at `numbers`.
    #[must_use]
    pub fn score(
        &self,
        class: usize,
        locus: usize,
        numbers: &ClassNumbers,
        settings: &HeldOutSettings,
    ) -> LocusScore {
        let evidence = &self.classes[class];
        score_locus(
            &evidence.tracts[locus],
            numbers,
            &evidence.libraries,
            self.homozygote_excess,
            evidence.read_span,
            settings,
            locus_seed(settings, self.class_seeds[class], locus),
        )
    }
}

/// Held-out sums for one method variant over one class.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HeldOutSums {
    pub ln: f64,
    pub variance: f64,
    pub truth_ln: f64,
    pub reads: u64,
    pub loci: usize,
    /// Of `loci`, those whose class was scored at calling's fallback rung (no concentration or no
    /// shares of its own).
    pub fallback_loci: usize,
}

/// **One held-out locus's score under one variant**, for the bench's `--locus-out` file.
#[derive(Debug, Clone, PartialEq)]
pub struct HeldOutLocusScore {
    /// Index into [`HeldOutRun::variants`].
    pub variant: usize,
    pub class: usize,
    pub locus: usize,
    /// **The fold the fit that scored this locus was made without** — the locus's own fold, which
    /// is what makes the score held out.
    pub fit_without_fold: usize,
    pub ln_likelihood: f64,
    pub truth_ln: f64,
    pub reads: u64,
    /// The read offset most of the locus's reads report ([`most_reported_offset`]).
    pub most_reported_offset: i32,
    /// The share the fit gave that offset's allele class.
    pub fitted_share: f64,
}

/// What one method's held-out run over one motif length produced: per fallback (each a variant,
/// named in `variants`), per class, the summed scores; the time its fits and its scoring took; and,
/// when asked for, every scored locus.
#[derive(Debug, Clone)]
pub struct HeldOutRun {
    pub variants: Vec<String>,
    pub sums: Vec<Vec<HeldOutSums>>,
    pub fit_seconds: f64,
    pub score_seconds: f64,
    /// Empty unless the run was asked to keep them; fold by fold, then variant by variant.
    pub loci: Vec<HeldOutLocusScore>,
}

/// One held-out locus's score under a variant: its class, its log-likelihood and that score's
/// Monte Carlo variance, the truth's log-likelihood, and its reads.
type HeldOutLocus = (usize, f64, f64, f64, u64);

/// One fallback's per-class numbers and the scores of the fold's loci at them, kept so that the
/// next fallback reuses a score where a class's numbers did not change.
type ScoredAtNumbers = (Vec<Option<ClassNumbers>>, Vec<Option<HeldOutLocus>>);

/// **The held-out score of one method over one motif length**: for each fold, the method fitted on
/// the other folds' loci, then every scored locus of the fold scored at the fitted numbers of its
/// class, under each fallback (`fallbacks` pairs each with the suffix the method's name gets to
/// name its variant). A class without its own shares or concentration is scored at the fallback rung and
/// counted, not skipped. With `keep_loci`, every scored locus is kept in [`HeldOutRun::loci`].
#[must_use]
pub fn held_out_run(
    method: &dyn PeriodFit,
    folded: &HeldOutClasses,
    settings: &HeldOutSettings,
    fallbacks: &[(&str, ScoringFallback)],
    keep_loci: bool,
) -> HeldOutRun {
    let name = method.name();
    let variants: Vec<String> = fallbacks
        .iter()
        .map(|(suffix, _)| format!("{name}{suffix}"))
        .collect();
    let classes = folded.classes;
    let mut sums = vec![vec![HeldOutSums::default(); classes.len()]; variants.len()];
    let mut loci = Vec::new();
    let mut fit_seconds = 0.0;
    let mut score_seconds = 0.0;
    for fold in 0..settings.folds {
        let training: Vec<StratumEvidence> = classes
            .iter()
            .zip(&folded.plans)
            .map(|(evidence, plan)| without_fold(evidence, plan, fold))
            .collect();
        let started = Instant::now();
        let estimates = method.fit_period(&training, folded.homozygote_excess, ALLELE_SPAN);
        fit_seconds += started.elapsed().as_secs_f64();
        let started = Instant::now();
        // **Only the loci of the fold the fit was made without** are scored with it.
        let in_fold: Vec<&FoldedLocus> = folded
            .scored
            .iter()
            .filter(|scored| scored.fold == fold)
            .collect();
        // The numbers and scores of the previous fallback: a class scored at the same numbers
        // again reuses its scores.
        let mut previous: Option<ScoredAtNumbers> = None;
        for (output, (_, fallback)) in fallbacks.iter().enumerate() {
            let scored_numbers = numbers_for_scoring(&estimates, *fallback);
            let numbers: Vec<Option<ClassNumbers>> = scored_numbers
                .iter()
                .map(|scored| scored.as_ref().map(|s| s.numbers.clone()))
                .collect();
            let on_the_fallback: Vec<bool> = scored_numbers
                .iter()
                .map(|scored| scored.as_ref().is_some_and(|s| s.on_the_fallback))
                .collect();
            let scores: Vec<Option<HeldOutLocus>> = in_fold
                .par_iter()
                .enumerate()
                .map(|(index, scored)| {
                    let numbers_here = numbers[scored.class].as_ref()?;
                    if let Some((earlier_numbers, earlier_scores)) = &previous
                        && earlier_numbers[scored.class].as_ref() == Some(numbers_here)
                    {
                        return earlier_scores[index];
                    }
                    let score = folded.score(scored.class, scored.locus, numbers_here, settings);
                    Some((
                        scored.class,
                        score.ln_likelihood,
                        score.variance,
                        scored.truth_ln,
                        score.reads,
                    ))
                })
                .collect();
            if keep_loci {
                for (scored, score) in in_fold.iter().zip(&scores) {
                    let Some((_, ln, _, truth, reads)) = score else {
                        continue;
                    };
                    let evidence = &classes[scored.class];
                    let (offset, share) = most_reported_offset(
                        &evidence.tracts[scored.locus],
                        evidence.read_span,
                        numbers[scored.class].as_ref(),
                    );
                    loci.push(HeldOutLocusScore {
                        variant: output,
                        class: scored.class,
                        locus: scored.locus,
                        fit_without_fold: fold,
                        ln_likelihood: *ln,
                        truth_ln: *truth,
                        reads: *reads,
                        most_reported_offset: offset,
                        fitted_share: share,
                    });
                }
            }
            for (class, ln, variance, truth, reads) in scores.iter().flatten().copied() {
                let into = &mut sums[output][class];
                into.ln += ln;
                into.variance += variance;
                into.truth_ln += truth;
                into.reads += reads;
                into.loci += 1;
                into.fallback_loci += usize::from(on_the_fallback[class]);
            }
            previous = Some((numbers, scores));
        }
        score_seconds += started.elapsed().as_secs_f64();
    }
    HeldOutRun {
        variants,
        sums,
        fit_seconds,
        score_seconds,
        loci,
    }
}

/// **The read offset most of a locus's reads report**, summed over its samples, and the share the
/// fit gave that offset's allele class; at an end bucket, the smaller of the two classes it can
/// hold (`±8` and `±9`), since either could be the allele. The share is `NaN` without numbers.
#[must_use]
pub fn most_reported_offset(
    tract: &TractReads,
    read_span: i32,
    numbers: Option<&ClassNumbers>,
) -> (i32, f64) {
    let mut totals = vec![0_u64; (2 * read_span + 1) as usize];
    for sample in &tract.samples {
        for (_, counts) in &sample.by_library {
            for (into, count) in totals.iter_mut().zip(counts) {
                *into += u64::from(*count);
            }
        }
    }
    let bucket = totals
        .iter()
        .enumerate()
        .max_by_key(|(_, count)| **count)
        .map_or(0, |(index, _)| index);
    let offset = bucket as i32 - read_span;
    let share = numbers.map_or(f64::NAN, |numbers| {
        let span = (numbers.allele_length_variation.len() as i32 - 1) / 2;
        let at = |offset: i32| numbers.allele_length_variation[(offset + span) as usize];
        if offset.abs() == read_span && span > read_span {
            at(offset).min(at(offset + offset.signum()))
        } else {
            at(offset)
        }
    });
    (offset, share)
}

// ---------------------------------------------------------------------
// The held-out score of one locus
// ---------------------------------------------------------------------

/// One locus's held-out score.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocusScore {
    /// `ln P(every sample's reads at the locus | the class's numbers)`.
    pub ln_likelihood: f64,
    /// The Monte Carlo variance of `ln_likelihood`; zero when it was summed exactly.
    pub variance: f64,
    /// Monte Carlo draws taken; zero when summed exactly.
    pub draws: usize,
    /// Spanning reads at the locus.
    pub reads: u64,
}

/// One sample's reads at a locus, in the terms the score needs: its genotype likelihoods (scaled
/// so the largest is one) as a symmetric matrix over the allele classes kept, its diagonal, and its
/// homozygote excess.
struct SampleTerms {
    matrix: Vec<f64>,
    diagonal: Vec<f64>,
    excess: f64,
}

/// **The held-out score of one locus**: `ln P(every sample's reads | numbers)`, the locus's
/// frequencies integrated over the class's Dirichlet — exactly for few samples, by importance
/// sampling otherwise (module documentation). `seed` fixes the Monte Carlo draws, so two methods
/// scored with the same seed share them.
///
/// # Panics
///
/// When the shares are not an odd number of classes, or a sample has no homozygote excess.
#[must_use]
pub fn score_locus(
    tract: &TractReads,
    numbers: &ClassNumbers,
    libraries: &[Library],
    homozygote_excess: &[f64],
    read_span: i32,
    settings: &HeldOutSettings,
    seed: u64,
) -> LocusScore {
    let classes = numbers.allele_length_variation.len();
    assert!(classes % 2 == 1, "the shares run from -span to +span");
    let allele_span = (classes / 2) as i32;
    let reads: u64 = tract
        .samples
        .iter()
        .flat_map(|sample| sample.by_library.iter())
        .flat_map(|(_, counts)| counts.iter())
        .map(|count| u64::from(*count))
        .sum();
    if tract.samples.is_empty() {
        return LocusScore {
            ln_likelihood: 0.0,
            variance: 0.0,
            draws: 0,
            reads,
        };
    }
    let genotypes = genotype_pairs(classes);
    let scored_libraries: Vec<Library> = libraries
        .iter()
        .map(|library| Library {
            slippage_group: 0,
            ..*library
        })
        .collect();
    let distributions = ReadDistributionsByLevel::of(
        &[numbers.slippage],
        &scored_libraries,
        allele_span,
        read_span,
    );
    // The classes the Dirichlet can put frequency on.
    let kept: Vec<usize> = (0..classes)
        .filter(|class| numbers.allele_length_variation[*class] > 0.0)
        .collect();
    let weights: Vec<f64> = kept
        .iter()
        .map(|class| numbers.concentration * numbers.allele_length_variation[*class])
        .collect();
    // **Only the genotypes the shares allow are scored**, so each sample's likelihoods are scaled
    // to the best of those. Scaled to the best of all genotypes, a sample of 300 reads whose own
    // allele has no share would see every allowed genotype underflow to zero and score −∞, where
    // the true log-likelihood is finite, if very low.
    let allowed: Vec<(usize, usize)> = genotypes
        .iter()
        .filter(|(first, second)| {
            numbers.allele_length_variation[*first] > 0.0
                && numbers.allele_length_variation[*second] > 0.0
        })
        .copied()
        .collect();
    let likelihoods = TractLikelihoods::of(tract, &distributions, &allowed, homozygote_excess);
    let mut position = vec![0_usize; classes * classes];
    for (index, (first, second)) in allowed.iter().enumerate() {
        position[first * classes + second] = index;
        position[second * classes + first] = index;
    }
    let width = kept.len();
    let samples: Vec<SampleTerms> = likelihoods
        .scaled
        .chunks_exact(likelihoods.width)
        .zip(&likelihoods.homozygote_excess)
        .map(|(row, excess)| {
            let mut matrix = vec![0.0; width * width];
            for (i, first) in kept.iter().enumerate() {
                for (j, second) in kept.iter().enumerate() {
                    matrix[i * width + j] = row[position[first * classes + second]];
                }
            }
            let diagonal = (0..width).map(|i| matrix[i * width + i]).collect();
            SampleTerms {
                matrix,
                diagonal,
                excess: *excess,
            }
        })
        .collect();

    let (ln_integral, variance, draws) =
        if !numbers.concentration.is_finite() || numbers.concentration <= 0.0 {
            (f64::NAN, 0.0, 0)
        } else if samples.len() <= settings.exact_up_to_samples {
            (float::ln(exact_integral(&samples, &weights)), 0.0, 0)
        } else {
            match settings.sampler {
                Sampler::Automatic => sampled_integral(&samples, &weights, settings, seed),
                Sampler::DirichletOnly => {
                    importance_sampled_integral(&samples, &weights, settings, seed)
                }
                Sampler::SequentialOnly => sequential_integral(&samples, &weights, settings, seed),
            }
        };
    LocusScore {
        ln_likelihood: ln_integral + likelihoods.ln_offset,
        variance,
        draws,
        reads,
    }
}

/// **The integral summed exactly**: each sample's genotype probability expanded as a polynomial in
/// the frequencies — `h p_a` for one draw, `(1 − h) p_a p_b` for two — and every product of terms
/// replaced by its Dirichlet moment, `Π_a α_a (α_a + 1) … / A (A + 1) …` with `A = Σ α`.
/// `(classes + classes²)^samples` products: for few samples only.
fn exact_integral(samples: &[SampleTerms], weights: &[f64]) -> f64 {
    let width = weights.len();
    let total_weight: f64 = weights.iter().sum();
    // Each sample's terms: a coefficient and the one or two classes it multiplies.
    let terms: Vec<Vec<(f64, usize, Option<usize>)>> = samples
        .iter()
        .map(|sample| {
            let mut terms = Vec::new();
            for a in 0..width {
                let one = sample.excess * sample.diagonal[a];
                if one > 0.0 {
                    terms.push((one, a, None));
                }
                for b in 0..width {
                    let two = (1.0 - sample.excess) * sample.matrix[a * width + b];
                    if two > 0.0 {
                        terms.push((two, a, Some(b)));
                    }
                }
            }
            terms
        })
        .collect();
    let mut counts = vec![0_u32; width];
    let mut sum = 0.0;
    fn walk(
        terms: &[Vec<(f64, usize, Option<usize>)>],
        sample: usize,
        coefficient: f64,
        counts: &mut [u32],
        weights: &[f64],
        total_weight: f64,
        sum: &mut f64,
    ) {
        if sample == terms.len() {
            let mut moment = 1.0;
            let mut drawn = 0_u32;
            for (class, count) in counts.iter().enumerate() {
                for step in 0..*count {
                    moment *= weights[class] + f64::from(step);
                }
                drawn += count;
            }
            for step in 0..drawn {
                moment /= total_weight + f64::from(step);
            }
            *sum += coefficient * moment;
            return;
        }
        for (term, a, b) in &terms[sample] {
            counts[*a] += 1;
            if let Some(b) = b {
                counts[*b] += 1;
            }
            walk(
                terms,
                sample + 1,
                coefficient * term,
                counts,
                weights,
                total_weight,
                sum,
            );
            counts[*a] -= 1;
            if let Some(b) = b {
                counts[*b] -= 1;
            }
        }
    }
    walk(&terms, 0, 1.0, &mut counts, weights, total_weight, &mut sum);
    sum
}

/// A drawn frequency below this is left out of [`ln_reads_given_frequencies`]'s sums: every term it
/// enters is at most the frequency itself (the scaled likelihoods are at most one), so a sample's
/// sum moves by at most `classes × 10⁻¹⁵` — below any probability a draw that matters gives it.
const NEGLIGIBLE_FREQUENCY: f64 = 1e-15;

/// `Σ over samples of ln(h p·d + (1 − h) pᵀ M p)`: the log of the reads' probability at frequencies
/// `p`, every sample's genotype summed over.
///
/// **Only the classes `p` gives a frequency above [`NEGLIGIBLE_FREQUENCY`] are summed**, listed once
/// for all samples in `active`. Under a class whose loci are mostly fixed at one length nearly every
/// drawn frequency but one or two is far below it, so this sums a few classes rather than nineteen,
/// and the triangle of the symmetric table rather than all of it.
fn ln_reads_given_frequencies(
    samples: &[SampleTerms],
    frequencies: &[f64],
    active: &mut Vec<usize>,
) -> f64 {
    let width = frequencies.len();
    active.clear();
    active.extend((0..width).filter(|class| frequencies[*class] > NEGLIGIBLE_FREQUENCY));
    let mut total = 0.0;
    for sample in samples {
        let mut one = 0.0;
        let mut two = 0.0;
        for (position, a) in active.iter().enumerate() {
            let row = &sample.matrix[a * width..(a + 1) * width];
            let pa = frequencies[*a];
            let mut off_diagonal = 0.0;
            for b in &active[position + 1..] {
                off_diagonal += row[*b] * frequencies[*b];
            }
            two += pa * (row[*a] * pa + 2.0 * off_diagonal);
            one += pa * sample.diagonal[*a];
        }
        total += float::ln((sample.excess * one + (1.0 - sample.excess) * two).max(1e-300));
    }
    total
}

/// A Dirichlet's log-density constant, `ln Γ(Σβ) − Σ ln Γ(β)`.
fn ln_dirichlet_constant(weights: &[f64]) -> f64 {
    ln_gamma(weights.iter().sum()) - weights.iter().map(|w| ln_gamma(*w)).sum::<f64>()
}

/// A Dirichlet's log-density at log-frequencies `ln_p`.
fn ln_dirichlet_density(weights: &[f64], constant: f64, ln_p: &[f64]) -> f64 {
    constant
        + weights
            .iter()
            .zip(ln_p)
            .map(|(weight, ln)| (weight - 1.0) * ln)
            .sum::<f64>()
}

/// One Dirichlet draw, **as log-frequencies**, so a class with a tiny weight is drawn at
/// `e^-10⁸` rather than rounded to zero and its density stays finite.
fn draw_ln_dirichlet(rng: &mut Rng, weights: &[f64], ln_p: &mut [f64]) {
    for (weight, out) in weights.iter().zip(ln_p.iter_mut()) {
        *out = if *weight >= 1.0 {
            float::ln(rng.gamma(*weight).max(1e-300))
        } else {
            // G(β) = G(β + 1) · U^(1/β), taken in logs.
            float::ln(rng.gamma(weight + 1.0).max(1e-300))
                + float::ln(rng.uniform().max(1e-300)) / weight
        };
    }
    let total = ln_sum_exp(ln_p);
    for out in ln_p.iter_mut() {
        *out -= total;
    }
}

/// **The proposal's components**: how much of the locus's expected allele counts each adds to the
/// class's Dirichlet weights, and the share of draws it gets. The reads leave each sample's genotype
/// uncertain, so the frequencies are less sharply placed than the expected counts would place them
/// were they observed: the half and quarter counts cover that. The last component, no counts, is
/// the class's own Dirichlet, which bounds every draw's weight at ten times the reads' probability.
const PROPOSAL: [(f64, f64); 4] = [(1.0, 0.45), (0.5, 0.30), (0.25, 0.15), (0.0, 0.10)];

/// Draws added at a time, until the error target is met.
const DRAWS_A_BATCH: usize = 1024;

/// **The integral by importance sampling**: `(ln of the estimate, its variance, draws)`.
///
/// The frequencies are drawn from a mixture of Dirichlets ([`PROPOSAL`]) built from the locus's
/// expected allele counts, which a few rounds of counting find; each draw is weighted by the
/// class's Dirichlet density times the reads' probability, over the mixture's density.
///
/// **Not a Dirichlet fitted to the weighted draws' moments**: a class with a tiny share has a
/// Dirichlet weight far below one, and matching moments sets it from the main class's spread —
/// 10⁻¹⁰ where the class's own is 10⁻², which made some draws' weights e^(10⁸) times others' and
/// left half of the loci of a 10-sample cohort short of the error target after 65,536 draws.
fn importance_sampled_integral(
    samples: &[SampleTerms],
    prior: &[f64],
    settings: &HeldOutSettings,
    seed: u64,
) -> (f64, f64, usize) {
    let width = prior.len();
    let mut rng = Rng(seed);
    let prior_total: f64 = prior.iter().sum();

    // The locus's expected allele counts: rounds of counting from the class's shares, each
    // sample's genotype weighed at the current frequencies.
    let mut frequencies: Vec<f64> = prior.iter().map(|w| w / prior_total).collect();
    let mut expected = vec![0.0; width];
    for _ in 0..30 {
        expected.iter_mut().for_each(|count| *count = 0.0);
        for sample in samples {
            let mut total = 0.0;
            for a in 0..width {
                total += sample.excess * frequencies[a] * sample.diagonal[a];
                for b in 0..width {
                    total += (1.0 - sample.excess)
                        * frequencies[a]
                        * frequencies[b]
                        * sample.matrix[a * width + b];
                }
            }
            if total <= 0.0 {
                continue;
            }
            for a in 0..width {
                expected[a] += sample.excess * frequencies[a] * sample.diagonal[a] / total;
                for b in 0..width {
                    // A pair of draws: one count for each copy, `a` here and `b` when the loop
                    // reaches it as the first.
                    expected[a] += 2.0
                        * (1.0 - sample.excess)
                        * frequencies[a]
                        * frequencies[b]
                        * sample.matrix[a * width + b]
                        / total;
                }
            }
        }
        let drawn: f64 = expected.iter().sum();
        for a in 0..width {
            frequencies[a] = (prior[a] + expected[a]) / (prior_total + drawn);
        }
    }
    let components: Vec<(Vec<f64>, f64, f64)> = PROPOSAL
        .iter()
        .map(|(counts, share)| {
            let weights: Vec<f64> = prior
                .iter()
                .zip(&expected)
                .map(|(w, n)| w + counts * n)
                .collect();
            let constant = ln_dirichlet_constant(&weights);
            (weights, constant, float::ln(*share))
        })
        .collect();
    let prior_constant = ln_dirichlet_constant(prior);

    let mut ln_p = vec![0.0; width];
    let mut p = vec![0.0; width];
    let mut active: Vec<usize> = Vec::with_capacity(width);
    let mut component_densities = vec![0.0; components.len()];
    let mut ln_weights: Vec<f64> = Vec::new();
    let mut estimate = (f64::NAN, f64::INFINITY);
    while ln_weights.len() < settings.max_draws.max(DRAWS_A_BATCH) {
        for _ in 0..DRAWS_A_BATCH {
            let mut pick = rng.uniform();
            let mut chosen = components.len() - 1;
            for (index, (_, share)) in PROPOSAL.iter().enumerate() {
                pick -= share;
                if pick <= 0.0 {
                    chosen = index;
                    break;
                }
            }
            draw_ln_dirichlet(&mut rng, &components[chosen].0, &mut ln_p);
            for (out, ln) in p.iter_mut().zip(ln_p.iter()) {
                *out = float::exp(*ln);
            }
            for ((weights, constant, ln_share), density) in
                components.iter().zip(component_densities.iter_mut())
            {
                *density = ln_share + ln_dirichlet_density(weights, *constant, &ln_p);
            }
            let ln_proposal = ln_sum_exp(&component_densities);
            let ln_prior = ln_dirichlet_density(prior, prior_constant, &ln_p);
            ln_weights.push(
                ln_prior + ln_reads_given_frequencies(samples, &p, &mut active) - ln_proposal,
            );
        }
        estimate = ln_mean_and_variance(&ln_weights);
        if estimate.1.sqrt() < settings.target_standard_error {
            break;
        }
    }
    (estimate.0, estimate.1, ln_weights.len())
}

/// Dirichlet draws tried before the sequential sampler takes over.
const FIRST_STAGE_DRAWS: usize = 8192;

/// **The integral by sampling**: the Dirichlet importance sampler first, and where it falls short
/// of the error target within [`FIRST_STAGE_DRAWS`], the sequential sampler
/// ([`sequential_integral`]); whichever ends with the smaller error is kept, its draws added to the
/// first's.
///
/// The Dirichlet sampler is cheap and right where the reads put the locus's frequencies in one
/// place. It fails where the reads allow two or more placements about equally — a locus whose
/// reads all sit in the −8 bucket, which an allele at −8 explains as well as one at −9 — because
/// under a class with a large between-locus share the locus is then "fixed at −8" or "fixed at −9",
/// two peaks a single Dirichlet cannot cover. On ten samples at three reads this was 91 loci in 700.
fn sampled_integral(
    samples: &[SampleTerms],
    prior: &[f64],
    settings: &HeldOutSettings,
    seed: u64,
) -> (f64, f64, usize) {
    let first_settings = HeldOutSettings {
        max_draws: settings.max_draws.min(FIRST_STAGE_DRAWS),
        ..*settings
    };
    let first = importance_sampled_integral(samples, prior, &first_settings, seed);
    if first.1.sqrt() < settings.target_standard_error {
        return first;
    }
    let second = sequential_integral(samples, prior, settings, seed ^ 0xA5A5_5A5A_F00D_BEEF);
    let draws = first.2 + second.2;
    if second.1 < first.1 {
        (second.0, second.1, draws)
    } else {
        (first.0, first.1, draws)
    }
}

/// Particles in one run of the sequential sampler.
const PARTICLES: usize = 256;
/// Independent runs added at a time, until the error target is met.
const RUNS_A_BATCH: usize = 8;

/// **The integral by sequential Monte Carlo over the samples**: `(ln of the estimate, its
/// variance, particles × runs)`.
///
/// The locus's frequencies are integrated out exactly by drawing the samples' alleles one after
/// another from a Pólya urn: given the alleles drawn so far, the next is allele `a` with chance
/// `(α_a + earlier draws of a) / (A + earlier draws)`. Each particle takes the samples in turn,
/// draws each sample's genotype in proportion to its urn chance times its reads' likelihood, and is
/// weighted by the sum it drew from — the probability of that sample's reads given the alleles
/// before it. The product of those averages over the particles is an unbiased estimate of the
/// integral; particles are resampled when a few carry most of the weight.
///
/// The error comes from repeating the whole run independently: the estimate is the mean of the
/// runs' estimates, and its error their spread. Runs are added [`RUNS_A_BATCH`] at a time until the
/// error is below the target or `max_draws` particles have been spent.
fn sequential_integral(
    samples: &[SampleTerms],
    prior: &[f64],
    settings: &HeldOutSettings,
    seed: u64,
) -> (f64, f64, usize) {
    let mut rng = Rng(seed);
    let mut runs: Vec<f64> = Vec::new();
    let mut estimate = (f64::NAN, f64::INFINITY);
    let most_runs = (settings.max_draws / PARTICLES).max(RUNS_A_BATCH);
    while runs.len() < most_runs {
        for _ in 0..RUNS_A_BATCH {
            runs.push(one_sequential_run(&mut rng, samples, prior));
        }
        estimate = ln_mean_and_variance(&runs);
        if estimate.1.sqrt() < settings.target_standard_error {
            break;
        }
    }
    (estimate.0, estimate.1, runs.len() * PARTICLES)
}

/// One run of [`sequential_integral`]: the log of its estimate.
fn one_sequential_run(rng: &mut Rng, samples: &[SampleTerms], prior: &[f64]) -> f64 {
    let width = prior.len();
    let prior_total: f64 = prior.iter().sum();
    let mut counts = vec![0.0_f64; PARTICLES * width];
    let mut drawn = vec![0.0_f64; PARTICLES];
    let mut ln_weights = vec![0.0_f64; PARTICLES];
    let mut ln_estimate = 0.0;
    // Per particle, the terms of its sum: `width` for one draw, then `width²` for two.
    let mut terms = vec![0.0_f64; width + width * width];
    for sample in samples {
        let largest = ln_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let normalised: Vec<f64> = ln_weights.iter().map(|w| float::exp(w - largest)).collect();
        let normaliser: f64 = normalised.iter().sum();
        let mut step = 0.0;
        for particle in 0..PARTICLES {
            let held = &mut counts[particle * width..(particle + 1) * width];
            let first_total = prior_total + drawn[particle];
            let mut total = 0.0;
            for a in 0..width {
                let first = (prior[a] + held[a]) / first_total;
                let one = sample.excess * first * sample.diagonal[a];
                terms[a] = one;
                total += one;
                let row = &sample.matrix[a * width..(a + 1) * width];
                for b in 0..width {
                    let earlier = held[b] + if a == b { 1.0 } else { 0.0 };
                    let second = (prior[b] + earlier) / (first_total + 1.0);
                    let two = (1.0 - sample.excess) * first * second * row[b];
                    terms[width + a * width + b] = two;
                    total += two;
                }
            }
            step += normalised[particle] / normaliser * total;
            if total > 0.0 {
                let mut pick = rng.uniform() * total;
                let mut chosen = terms.len() - 1;
                for (index, term) in terms.iter().enumerate() {
                    pick -= term;
                    if pick <= 0.0 {
                        chosen = index;
                        break;
                    }
                }
                if chosen < width {
                    held[chosen] += 1.0;
                    drawn[particle] += 1.0;
                } else {
                    let pair = chosen - width;
                    held[pair / width] += 1.0;
                    held[pair % width] += 1.0;
                    drawn[particle] += 2.0;
                }
            }
            ln_weights[particle] += float::ln(total.max(1e-300));
        }
        ln_estimate += float::ln(step.max(1e-300));
        // Resample, systematically, when fewer than half the particles carry the weight.
        let largest = ln_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let weights: Vec<f64> = ln_weights.iter().map(|w| float::exp(w - largest)).collect();
        let sum: f64 = weights.iter().sum();
        let effective = sum * sum / weights.iter().map(|w| w * w).sum::<f64>();
        if effective < PARTICLES as f64 / 2.0 {
            let start = rng.uniform() / PARTICLES as f64;
            let mut cumulative = 0.0;
            let mut source = 0;
            let mut new_counts = vec![0.0_f64; PARTICLES * width];
            let mut new_drawn = vec![0.0_f64; PARTICLES];
            for target in 0..PARTICLES {
                let point = start + target as f64 / PARTICLES as f64;
                while source < PARTICLES - 1 && cumulative + weights[source] / sum < point {
                    cumulative += weights[source] / sum;
                    source += 1;
                }
                new_counts[target * width..(target + 1) * width]
                    .copy_from_slice(&counts[source * width..(source + 1) * width]);
                new_drawn[target] = drawn[source];
            }
            counts = new_counts;
            drawn = new_drawn;
            ln_weights.iter_mut().for_each(|w| *w = 0.0);
        }
    }
    ln_estimate
}

/// The log of the mean of `exp(ln_weights)`, and that log's variance by the delta method,
/// `Var(mean) / mean²`.
fn ln_mean_and_variance(ln_weights: &[f64]) -> (f64, f64) {
    let largest = ln_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() {
        return (f64::NEG_INFINITY, f64::INFINITY);
    }
    let count = ln_weights.len() as f64;
    let weights: Vec<f64> = ln_weights.iter().map(|w| float::exp(w - largest)).collect();
    let mean = weights.iter().sum::<f64>() / count;
    let variance = weights
        .iter()
        .map(|w| float::powi(w - mean, 2))
        .sum::<f64>()
        / (count - 1.0);
    (float::ln(mean) + largest, variance / count / (mean * mean))
}

// ---------------------------------------------------------------------
// Genotype calls under two sets of numbers (plan step D6)
// ---------------------------------------------------------------------

/// **A class's samples, ready to be called under any set of numbers**: each (locus, sample) with
/// reads given the genotype its weights put highest ([`Self::under`]). A proxy for calling, for the
/// question of spec `str_penalised_fit.md` §5.1 — whether a slip rate a little off changes calls —
/// on cohorts too large for the project to simulate reads for.
///
/// **The weights are the fit's own**, the expectation step's: a genotype's prior from the class's
/// shares and spread and the sample's homozygote excess, times how well it explains the sample's
/// reads. Calling also borrows from the other samples at the locus; this does not.
pub struct ClassCaller<'a> {
    evidence: &'a StratumEvidence,
    distinct: DistinctRows,
}

impl<'a> ClassCaller<'a> {
    /// `evidence`'s distinct rows of reads, gathered once for every set of numbers.
    #[must_use]
    pub fn of(evidence: &'a StratumEvidence, homozygote_excess: &[f64]) -> Self {
        Self {
            evidence,
            distinct: DistinctRows::of(evidence, homozygote_excess),
        }
    }

    /// **Each (locus, sample) with reads called under `numbers`**, in the order of the class's
    /// loci then each locus's samples: the genotype as its two allele classes over the shares'
    /// window (`0` is `-allele_span`), the first not above the second, or `None` where no genotype
    /// explains the reads.
    #[must_use]
    pub fn under(&self, numbers: &ClassNumbers) -> Vec<Option<(usize, usize)>> {
        let allele_span = (numbers.allele_length_variation.len() / 2) as i32;
        ScoredDistinctRows::of(self.evidence, &self.distinct, numbers.slippage, allele_span)
            .genotypes_weighted_highest(&numbers.allele_length_variation, numbers.concentration)
    }
}

/// **The genotypes a class's samples drew, as [`ClassCaller::under`] names genotypes**: per locus,
/// per sample listed at it (`DrawnPeriod::genotypes` of one class), each copy's offset folded into
/// `-allele_span ..= +allele_span` (the outermost classes meaning "this far or further") and the
/// pair put in order, flattened in the order of the class's loci then each locus's samples.
#[must_use]
pub fn drawn_genotypes_in_window(
    genotypes: &[Vec<(i8, i8)>],
    allele_span: i32,
) -> Vec<(usize, usize)> {
    let class =
        |offset: i8| (i32::from(offset).clamp(-allele_span, allele_span) + allele_span) as usize;
    genotypes
        .iter()
        .flatten()
        .map(|(first, second)| {
            let (first, second) = (class(*first), class(*second));
            (first.min(second), first.max(second))
        })
        .collect()
}

/// **How the calls under one set of numbers differ from the calls under another**, the
/// *reference*, each held against the genotype drawn. A missing call is a wrong one.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CallChanges {
    /// (locus, sample) pairs with reads.
    pub pairs: usize,
    /// Pairs whose call under the reference is not the drawn genotype: the error already there.
    pub wrong_under_reference: usize,
    /// Pairs whose call differs between the two.
    pub changed: usize,
    /// Changed pairs the reference called right.
    pub right_to_wrong: usize,
    /// Changed pairs the other numbers called right.
    pub wrong_to_right: usize,
}

impl CallChanges {
    /// Changed pairs that both called wrong, each differently.
    #[must_use]
    pub fn wrong_to_wrong(&self) -> usize {
        self.changed - self.right_to_wrong - self.wrong_to_right
    }

    /// Adds `other`'s pairs to these.
    pub fn add(&mut self, other: &Self) {
        self.pairs += other.pairs;
        self.wrong_under_reference += other.wrong_under_reference;
        self.changed += other.changed;
        self.right_to_wrong += other.right_to_wrong;
        self.wrong_to_right += other.wrong_to_right;
    }
}

/// **`other`'s calls against `reference`'s**, pair by pair, with `drawn` the genotypes drawn, all
/// three in the same order ([`ClassCaller::under`], [`drawn_genotypes_in_window`]).
///
/// # Panics
///
/// When the three are not the same length.
#[must_use]
pub fn call_changes(
    reference: &[Option<(usize, usize)>],
    other: &[Option<(usize, usize)>],
    drawn: &[(usize, usize)],
) -> CallChanges {
    assert!(
        reference.len() == other.len() && other.len() == drawn.len(),
        "one call under each set of numbers and one drawn genotype a pair"
    );
    let mut changes = CallChanges {
        pairs: drawn.len(),
        ..CallChanges::default()
    };
    for ((reference, other), drawn) in reference.iter().zip(other).zip(drawn) {
        let reference_right = *reference == Some(*drawn);
        if !reference_right {
            changes.wrong_under_reference += 1;
        }
        if reference != other {
            changes.changed += 1;
            if reference_right {
                changes.right_to_wrong += 1;
            } else if *other == Some(*drawn) {
                changes.wrong_to_right += 1;
            }
        }
    }
    changes
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The brute-force sum the exact path must equal**: every sample's identity-by-descent coin
    /// and ordered alleles enumerated, the frequencies integrated as a Pólya urn — each draw takes
    /// class `a` with chance `(α_a + earlier draws of a) / (A + earlier draws)` — and the reads'
    /// probability computed straight from [`Slippage::read_probabilities`].
    fn brute_force(
        tract: &TractReads,
        numbers: &ClassNumbers,
        excess: &[f64],
        read_span: i32,
    ) -> f64 {
        let classes = numbers.allele_length_variation.len();
        let span = (classes / 2) as i32;
        let weights: Vec<f64> = numbers
            .allele_length_variation
            .iter()
            .map(|share| numbers.concentration * share)
            .collect();
        let per_allele: Vec<Vec<f64>> = (0..classes)
            .map(|class| {
                numbers
                    .slippage
                    .read_probabilities(class as i32 - span, read_span)
            })
            .collect();
        let reads_given = |sample: &SampleTractReads, a: usize, b: usize| -> f64 {
            let counts = &sample.by_library[0].1;
            counts
                .iter()
                .enumerate()
                .map(|(bucket, count)| {
                    float::powi(
                        0.5 * (per_allele[a][bucket] + per_allele[b][bucket]),
                        *count as i32,
                    )
                })
                .product()
        };
        // Each sample: one draw (excess) or two; each case enumerated.
        fn recurse(
            tract: &TractReads,
            sample: usize,
            urn: &mut Vec<f64>,
            weight: f64,
            excess: &[f64],
            reads_given: &dyn Fn(&SampleTractReads, usize, usize) -> f64,
            total: &mut f64,
        ) {
            if sample == tract.samples.len() {
                *total += weight;
                return;
            }
            let reads = &tract.samples[sample];
            let h = excess[reads.sample as usize];
            let classes = urn.len();
            let urn_total: f64 = urn.iter().sum();
            for a in 0..classes {
                let pa = urn[a] / urn_total;
                // One draw.
                urn[a] += 1.0;
                recurse(
                    tract,
                    sample + 1,
                    urn,
                    weight * h * pa * reads_given(reads, a, a),
                    excess,
                    reads_given,
                    total,
                );
                // Two draws.
                let after_a: f64 = urn.iter().sum();
                for b in 0..classes {
                    let pb = urn[b] / after_a;
                    urn[b] += 1.0;
                    recurse(
                        tract,
                        sample + 1,
                        urn,
                        weight * (1.0 - h) * pa * pb * reads_given(reads, a, b),
                        excess,
                        reads_given,
                        total,
                    );
                    urn[b] -= 1.0;
                }
                urn[a] -= 1.0;
            }
        }
        let mut urn = weights;
        let mut total = 0.0;
        recurse(tract, 0, &mut urn, 1.0, excess, &reads_given, &mut total);
        float::ln(total)
    }

    fn tiny_locus(samples: &[[u32; 3]]) -> TractReads {
        TractReads {
            samples: samples
                .iter()
                .enumerate()
                .map(|(sample, counts)| SampleTractReads {
                    sample: sample as u32,
                    by_library: vec![(LibraryIndex(0), counts.to_vec())],
                })
                .collect(),
        }
    }

    fn tiny_numbers() -> ClassNumbers {
        ClassNumbers {
            slippage: Slippage {
                level: 0.15,
                shorter_share: 0.7,
                fall_off: 0.3,
            },
            allele_length_variation: vec![0.2, 0.7, 0.1],
            concentration: 1.5,
        }
    }

    fn one_library() -> Vec<Library> {
        vec![Library::in_group(ReadGroupId(0), 0)]
    }

    /// **The exact score is the brute-force sum**, on two samples and three allele classes.
    #[test]
    fn the_exact_score_is_the_brute_force_sum() {
        let tract = tiny_locus(&[[1, 3, 0], [0, 2, 2]]);
        let numbers = tiny_numbers();
        let excess = [0.3, 0.8];
        let settings = HeldOutSettings::default();
        let scored = score_locus(&tract, &numbers, &one_library(), &excess, 1, &settings, 1);
        let expected = brute_force(&tract, &numbers, &excess, 1);
        assert_eq!(scored.draws, 0, "two samples are summed exactly");
        assert!(
            (scored.ln_likelihood - expected).abs() < 1e-10,
            "{} against {expected}",
            scored.ln_likelihood
        );
        assert_eq!(scored.reads, 8);
    }

    /// **A deep sample whose own allele has no share scores a finite log-likelihood**: 300 reads at
    /// the one length the shares leave out, so every allowed genotype is hundreds of log units less
    /// likely than the excluded one.
    #[test]
    fn a_deep_sample_at_an_allele_without_a_share_scores_finitely() {
        let numbers = ClassNumbers {
            allele_length_variation: vec![0.5, 0.0, 0.5],
            ..tiny_numbers()
        };
        let tract = tiny_locus(&[[2, 300, 1]]);
        let scored = score_locus(
            &tract,
            &numbers,
            &one_library(),
            &[0.5],
            1,
            &HeldOutSettings::default(),
            1,
        );
        assert!(
            scored.ln_likelihood.is_finite() && scored.ln_likelihood < -500.0,
            "{}",
            scored.ln_likelihood
        );
    }

    /// **Both Monte Carlo samplers agree with the brute-force sum within their own error**, and
    /// their error is below the target: on the two-sample locus forced through them, and on a
    /// five-sample one where the brute force is the only exact answer.
    #[test]
    fn the_monte_carlo_score_agrees_with_the_brute_force_sum() {
        let numbers = tiny_numbers();
        for (sampler, (counts, excess)) in [Sampler::DirichletOnly, Sampler::SequentialOnly]
            .into_iter()
            .flat_map(|sampler| {
                [
                    (vec![[1, 3, 0], [0, 2, 2]], vec![0.3, 0.8]),
                    (
                        vec![[0, 4, 0], [2, 1, 0], [0, 0, 3], [1, 1, 1], [0, 5, 0]],
                        vec![0.0, 0.5, 0.9, 0.2, 0.95],
                    ),
                ]
                .into_iter()
                .map(move |case| (sampler, case))
            })
        {
            let forced = HeldOutSettings {
                exact_up_to_samples: 0,
                sampler,
                ..HeldOutSettings::default()
            };
            let tract = tiny_locus(&counts);
            let expected = brute_force(&tract, &numbers, &excess, 1);
            let scored = score_locus(&tract, &numbers, &one_library(), &excess, 1, &forced, 9);
            let error = scored.variance.sqrt();
            assert!(scored.draws > 0);
            assert!(error < forced.target_standard_error, "error {error}");
            assert!(
                (scored.ln_likelihood - expected).abs() < 4.0 * error.max(1e-3),
                "{sampler:?}: {} against {expected}, error {error}",
                scored.ln_likelihood
            );
        }
    }

    /// **Where two allele classes explain the reads about equally, under a class whose loci are
    /// mostly fixed at one length**, the frequencies' posterior has two peaks ("fixed at −2" and
    /// "fixed at −1": reads are recorded to ±1, so both land in the −1 bucket). The automatic choice
    /// and the sequential sampler both reach the target and agree with the brute-force sum.
    #[test]
    fn the_samplers_agree_where_two_alleles_explain_the_reads_equally() {
        let numbers = ClassNumbers {
            slippage: Slippage {
                level: 0.1,
                shorter_share: 0.6,
                fall_off: 0.3,
            },
            allele_length_variation: vec![0.1, 0.1, 0.6, 0.1, 0.1],
            concentration: 0.3,
        };
        let tract = tiny_locus(&[[3, 0, 0], [2, 1, 0], [4, 0, 0], [3, 0, 0], [1, 0, 0]]);
        let excess = [0.9, 0.95, 0.9, 0.8, 0.99];
        let expected = brute_force(&tract, &numbers, &excess, 1);
        for sampler in [Sampler::Automatic, Sampler::SequentialOnly] {
            let settings = HeldOutSettings {
                exact_up_to_samples: 0,
                sampler,
                ..HeldOutSettings::default()
            };
            let scored = score_locus(&tract, &numbers, &one_library(), &excess, 1, &settings, 5);
            let error = scored.variance.sqrt();
            assert!(
                error < settings.target_standard_error,
                "{sampler:?}: error {error}"
            );
            assert!(
                (scored.ln_likelihood - expected).abs() < 4.0 * error.max(1e-3),
                "{sampler:?}: {} against {expected}, error {error}",
                scored.ln_likelihood
            );
        }
    }

    /// **The sampler reaches its error target where samples are many and the reads pin the
    /// frequencies down** — 400 samples at 30 reads, where a draw from the class's own Dirichlet
    /// almost never lands where the reads put the frequencies — and two seeds agree within their
    /// errors.
    #[test]
    fn the_sampler_reaches_its_target_on_a_deep_wide_locus() {
        let numbers = tiny_numbers();
        let samples: Vec<[u32; 3]> = (0..400)
            .map(|index| match index % 4 {
                0 => [3, 25, 2],
                1 => [1, 27, 2],
                2 => [14, 14, 2],
                _ => [0, 28, 2],
            })
            .collect();
        let tract = tiny_locus(&samples);
        let excess = vec![0.5; 400];
        let settings = HeldOutSettings::default();
        let scored = score_locus(&tract, &numbers, &one_library(), &excess, 1, &settings, 3);
        assert!(scored.ln_likelihood.is_finite());
        assert!(
            scored.variance.sqrt() < settings.target_standard_error,
            "error {} after {} draws",
            scored.variance.sqrt(),
            scored.draws
        );
        // Two independent seeds land within their errors of each other.
        let again = score_locus(&tract, &numbers, &one_library(), &excess, 1, &settings, 4);
        let joint = (scored.variance + again.variance).sqrt();
        assert!((scored.ln_likelihood - again.ln_likelihood).abs() < 4.0 * joint);
    }

    #[test]
    fn folds_cover_every_locus_once_and_cap_the_scored_loci() {
        let evidence = bench_fixtures::draw_stratum(
            Slippage {
                level: 0.05,
                shorter_share: 0.6,
                fall_off: 0.3,
            },
            &spectrum_of(3),
            1.0,
            0.5,
            53,
            3,
            2,
            1,
            5,
        );
        let settings = HeldOutSettings {
            loci_scored_per_class: 20,
            ..HeldOutSettings::default()
        };
        let plan = FoldPlan::of(&evidence, &settings, 7);
        let sizes: Vec<usize> = (0..5)
            .map(|fold| plan.fold_of_locus.iter().filter(|f| **f == fold).count())
            .collect();
        assert_eq!(sizes.iter().sum::<usize>(), 53);
        assert!(sizes.iter().all(|size| (10..=11).contains(size)));
        assert_eq!(plan.scored.iter().filter(|s| **s).count(), 20);
        let training = without_fold(&evidence, &plan, 2);
        assert_eq!(training.tracts.len(), 53 - sizes[2]);
    }

    /// **The psp mode's folds: whole loci, each held out exactly once, fixed by the seed.** Two
    /// classes at the same repeat count and different motif lengths, as a real cohort has them.
    /// For every fold, the loci fitted are the class's loci outside the fold with every sample's
    /// reads intact, and the loci held out are the rest; over the folds every locus is held out
    /// once. The same seed folds the same way; another seed, and the other motif length at the
    /// same repeat count, fold differently.
    #[test]
    fn the_psp_modes_folds_hold_each_whole_locus_out_once_and_follow_the_seed() {
        let drawn = |period: u8, seed: u64| {
            let mut evidence = bench_fixtures::draw_stratum(
                Slippage {
                    level: 0.05,
                    shorter_share: 0.6,
                    fall_off: 0.3,
                },
                &spectrum_of(3),
                1.0,
                0.5,
                40,
                4,
                3,
                1,
                seed,
            );
            evidence.stratum = Stratum {
                period,
                reference_repeats: 10,
            };
            evidence
        };
        let classes = vec![drawn(1, 5), drawn(2, 6)];
        let settings = HeldOutSettings {
            folds: 5,
            ..HeldOutSettings::default()
        };
        let plans = fold_plans_across_motif_lengths(&classes, &settings);
        assert_eq!(plans.len(), 2);
        for (evidence, plan) in classes.iter().zip(&plans) {
            let loci = evidence.tracts.len();
            assert_eq!(plan.fold_of_locus.len(), loci);
            let mut times_held_out = vec![0_usize; loci];
            for fold in 0..settings.folds {
                let held_out: Vec<usize> = (0..loci)
                    .filter(|&locus| plan.fold_of_locus[locus] == fold)
                    .collect();
                assert!(!held_out.is_empty(), "fold {fold} holds no locus");
                for &locus in &held_out {
                    times_held_out[locus] += 1;
                }
                // Whole loci: the fitted set is the class's other tracts, each with all its
                // samples, in their own order.
                let expected: Vec<TractReads> = (0..loci)
                    .filter(|&locus| plan.fold_of_locus[locus] != fold)
                    .map(|locus| evidence.tracts[locus].clone())
                    .collect();
                assert_eq!(without_fold(evidence, plan, fold).tracts, expected);
            }
            assert!(times_held_out.iter().all(|&times| times == 1));
        }
        // Fixed by the seed.
        assert_eq!(fold_plans_across_motif_lengths(&classes, &settings), plans);
        let reseeded = fold_plans_across_motif_lengths(
            &classes,
            &HeldOutSettings {
                seed: settings.seed + 1,
                ..settings
            },
        );
        assert_ne!(reseeded[0].fold_of_locus, plans[0].fold_of_locus);
        // The two motif lengths at ten repeats are not shuffled alike.
        assert_ne!(
            FoldPlan::of(
                &classes[1],
                &settings,
                class_seed_across_motif_lengths(classes[1].stratum)
            )
            .fold_of_locus,
            FoldPlan::of(
                &classes[1],
                &settings,
                class_seed_across_motif_lengths(classes[0].stratum)
            )
            .fold_of_locus,
        );
    }

    /// **A stand-in method that records what each of its fits was given**, and answers with numbers
    /// that tell its fits apart: fit `k` (counting from 0) returns concentration `1 + k`.
    struct RecordingFit {
        /// Per fit, in the order they were made: per class, the tracts it was handed.
        handed: std::sync::Mutex<Vec<Vec<Vec<TractReads>>>>,
    }

    impl RecordingFit {
        fn numbers_of_fit(fit: usize) -> ClassNumbers {
            ClassNumbers {
                slippage: Slippage {
                    level: 0.1,
                    shorter_share: 0.6,
                    fall_off: 0.3,
                },
                allele_length_variation: vec![0.25, 0.5, 0.25],
                concentration: 1.0 + fit as f64,
            }
        }
    }

    impl PeriodFit for RecordingFit {
        fn name(&self) -> String {
            "recording".to_owned()
        }

        fn fit_period(
            &self,
            classes: &[StratumEvidence],
            _homozygote_excess: &[f64],
            _allele_span: i32,
        ) -> Vec<ClassEstimate> {
            let mut handed = self.handed.lock().expect("no fit panicked");
            let fit = handed.len();
            handed.push(
                classes
                    .iter()
                    .map(|evidence| evidence.tracts.clone())
                    .collect(),
            );
            let numbers = Self::numbers_of_fit(fit);
            classes
                .iter()
                .map(|evidence| ClassEstimate {
                    slippage: Some(numbers.slippage),
                    allele_length_variation: Some(numbers.allele_length_variation.clone()),
                    concentration: Some(numbers.concentration),
                    ..ClassEstimate::absent(evidence)
                })
                .collect()
        }
    }

    /// **Every held-out locus is scored by the fit made without it.** A stand-in method records
    /// the loci each fold's fit was handed and answers with numbers unique to that fit. Then, for
    /// every locus scored: the fit that scored it was made without the locus's fold, was not
    /// handed the locus, and its numbers give exactly the recorded score. Every locus picked for
    /// scoring is scored once, and each fit was handed every locus outside its fold.
    #[test]
    fn every_held_out_locus_is_scored_by_the_fit_made_without_it() {
        // Two classes of 12 loci, two samples each, every locus's reads distinct, so a locus can
        // be recognised among the tracts a fit was handed.
        let class = |period: u8, first: u32| {
            let mut evidence = bench_fixtures::draw_stratum(
                Slippage {
                    level: 0.05,
                    shorter_share: 0.6,
                    fall_off: 0.3,
                },
                &spectrum_of(3),
                1.0,
                0.5,
                12,
                2,
                3,
                1,
                u64::from(first),
            );
            evidence.stratum = Stratum {
                period,
                reference_repeats: 10,
            };
            evidence.tracts = (0..12)
                .map(|locus| tiny_locus(&[[first + locus, 2, 0], [0, 3, locus + 1]]))
                .collect();
            evidence
        };
        let classes = vec![class(1, 1), class(2, 100)];
        let excess = [0.3, 0.8];
        let settings = HeldOutSettings {
            folds: 4,
            ..HeldOutSettings::default()
        };
        let folded = HeldOutClasses::new(
            &classes,
            &excess,
            fold_plans_across_motif_lengths(&classes, &settings),
            classes
                .iter()
                .map(|evidence| class_seed_across_motif_lengths(evidence.stratum))
                .collect(),
        );
        assert_eq!(folded.scored().len(), 24, "every locus has reads");
        let method = RecordingFit {
            handed: std::sync::Mutex::new(Vec::new()),
        };
        let run = held_out_run(
            &method,
            &folded,
            &settings,
            &[("", ScoringFallback::default())],
            true,
        );
        let handed = method.handed.into_inner().expect("no fit panicked");
        assert_eq!(handed.len(), settings.folds, "one fit a fold");

        let plan_of = fold_plans_across_motif_lengths(&classes, &settings);
        // Each fit was handed exactly the loci outside its fold.
        for (fit, by_class) in handed.iter().enumerate() {
            for (class, tracts) in by_class.iter().enumerate() {
                let expected: Vec<TractReads> = (0..12)
                    .filter(|&locus| plan_of[class].fold_of_locus[locus] != fit)
                    .map(|locus| classes[class].tracts[locus].clone())
                    .collect();
                assert_eq!(tracts, &expected, "fit {fit}, class {class}");
            }
        }
        // Each scored locus: scored by the fit made without its fold, which never saw it, at that
        // fit's numbers.
        let mut times_scored = vec![vec![0_usize; 12]; 2];
        for scored in &run.loci {
            let fit = scored.fit_without_fold;
            assert_eq!(
                fit, plan_of[scored.class].fold_of_locus[scored.locus],
                "locus {} of class {} scored by the fit without fold {fit}",
                scored.locus, scored.class
            );
            assert!(
                !handed[fit][scored.class].contains(&classes[scored.class].tracts[scored.locus]),
                "the fit that scored locus {} of class {} was handed it",
                scored.locus,
                scored.class
            );
            let expected = folded.score(
                scored.class,
                scored.locus,
                &RecordingFit::numbers_of_fit(fit),
                &settings,
            );
            assert_eq!(
                scored.ln_likelihood.to_bits(),
                expected.ln_likelihood.to_bits(),
                "locus {} of class {}: not the score at fit {fit}'s numbers",
                scored.locus,
                scored.class
            );
            times_scored[scored.class][scored.locus] += 1;
        }
        assert!(times_scored.iter().flatten().all(|&times| times == 1));
        assert_eq!(run.variants, ["recording"]);
        assert_eq!(
            run.sums[0].iter().map(|sums| sums.loci).sum::<usize>(),
            24,
            "every scored locus is in the sums"
        );
    }

    /// **A class without a concentration is scored at calling's bottom rung, not skipped**: a flat
    /// shape at the median of the other classes' concentrations, marked; a class with no stutter is
    /// still not scored; with no concentration anywhere the stated constant is used.
    #[test]
    fn a_class_without_a_concentration_is_scored_at_the_fallback() {
        let slippage = Slippage {
            level: 0.05,
            shorter_share: 0.7,
            fall_off: 0.2,
        };
        let stratum = |reference_repeats| Stratum {
            period: 1,
            reference_repeats,
        };
        let estimate = |repeats, concentration: Option<f64>, with_stutter: bool| ClassEstimate {
            stratum: stratum(repeats),
            slippage: with_stutter.then_some(slippage),
            allele_length_variation: Some(vec![0.2, 0.6, 0.2]),
            concentration,
            errors: ClassErrors::default(),
        };
        let estimates = [
            estimate(8, Some(0.2), true),
            estimate(9, Some(0.4), true),
            estimate(10, Some(3.0), true),
            estimate(11, None, true),
            estimate(12, None, false),
        ];
        let scored = numbers_for_scoring(&estimates, ScoringFallback::default());
        let own = scored[1]
            .as_ref()
            .expect("a class with its numbers is scored");
        assert!(!own.on_the_fallback);
        assert_eq!(own.numbers.concentration, 0.4);
        let fallback = scored[3]
            .as_ref()
            .expect("a class without a concentration is scored");
        assert!(fallback.on_the_fallback);
        assert_eq!(
            fallback.numbers.concentration, 0.4,
            "the median of 0.2, 0.4 and 3.0"
        );
        assert_eq!(fallback.numbers.allele_length_variation, vec![1.0 / 3.0; 3]);
        assert!(scored[4].is_none());

        let none = numbers_for_scoring(&[estimate(8, None, true)], ScoringFallback::default());
        assert_eq!(
            none[0].as_ref().map(|s| s.numbers.concentration),
            Some(crate::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION)
        );

        // Kept shares, and a stated concentration in place of the median.
        let kept = numbers_for_scoring(
            &estimates,
            ScoringFallback {
                keep_own_shares: true,
                stated_concentration: Some(2.0),
            },
        );
        let fallback = kept[3].as_ref().expect("scored");
        assert!(fallback.on_the_fallback);
        assert_eq!(
            fallback.numbers.allele_length_variation,
            vec![0.2, 0.6, 0.2]
        );
        assert_eq!(fallback.numbers.concentration, 2.0);
        assert_eq!(kept[1].as_ref().map(|s| s.numbers.concentration), Some(0.4));
    }

    /// **The same stutter changes no call; a grossly wrong one changes some.** Two drawn two-base
    /// classes, 40 samples at 3 reads: called under the truth twice, from two separate scorings,
    /// no call moves; called with the slip rate five times the truth's, some do, and more of the
    /// moved calls go from right to wrong than back.
    ///
    /// **The pairing of calls with drawn genotypes.** Most samples here are homozygous for the
    /// reference, so a call held against another sample's genotype is still often right; "most
    /// calls are right" alone would pass with the pairing broken (measured: 1,743 of 2,293 right
    /// with every locus's genotypes shifted by one sample). So the test holds the same calls
    /// against the drawn genotypes shifted that way too, and asks the shifted pairing for more
    /// than twice the wrong calls (measured: 550 against 142).
    #[test]
    fn the_same_stutter_changes_no_call_and_a_wrong_one_changes_some() {
        use super::super::cohort_draw::{CohortModel, LociPerClass, StutterRegime, draw_period};
        let mut model = CohortModel::tomato_like(2, StutterRegime::Pcr);
        model.samples = 40;
        model.depth = 3.0;
        model.loci = LociPerClass::Each(30);
        model.repeat_counts = Some(vec![10, 16]);
        let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
        let mut same = CallChanges::default();
        let mut wrong = CallChanges::default();
        let mut shifted = CallChanges::default();
        for ((evidence, truth), drawn_genotypes) in
            drawn.strata.iter().zip(&drawn.truth).zip(&drawn.genotypes)
        {
            let caller = ClassCaller::of(evidence, &drawn.homozygote_excess);
            let numbers = ClassNumbers::of_truth(truth, ALLELE_SPAN);
            let genotypes = drawn_genotypes_in_window(drawn_genotypes, ALLELE_SPAN);
            let under_truth = caller.under(&numbers);
            // Each locus's genotypes moved one sample along: every call paired with its
            // neighbour's genotype.
            let shifted_genotypes: Vec<Vec<(i8, i8)>> = drawn_genotypes
                .iter()
                .map(|locus| {
                    let mut locus = locus.clone();
                    if !locus.is_empty() {
                        locus.rotate_left(1);
                    }
                    locus
                })
                .collect();
            shifted.add(&call_changes(
                &under_truth,
                &under_truth,
                &drawn_genotypes_in_window(&shifted_genotypes, ALLELE_SPAN),
            ));
            same.add(&call_changes(
                &under_truth,
                &caller.under(&numbers.clone()),
                &genotypes,
            ));
            let too_much_stutter = ClassNumbers {
                slippage: Slippage {
                    level: 5.0 * numbers.slippage.level,
                    ..numbers.slippage
                },
                ..numbers.clone()
            };
            wrong.add(&call_changes(
                &under_truth,
                &caller.under(&too_much_stutter),
                &genotypes,
            ));
        }
        assert!(same.pairs > 1_000, "{same:?}");
        assert_eq!(same.changed, 0, "{same:?}");
        assert!(
            2 * same.wrong_under_reference < shifted.wrong_under_reference,
            "calls paired with the drawn genotypes are wrong far less often than with each \
             locus's genotypes shifted by one sample: {same:?} against {shifted:?}"
        );
        assert!(wrong.changed > 0, "{wrong:?}");
        assert!(wrong.right_to_wrong > wrong.wrong_to_right, "{wrong:?}");
        assert_eq!(
            wrong.right_to_wrong + wrong.wrong_to_right + wrong.wrong_to_wrong(),
            wrong.changed
        );
    }

    /// **A drawn copy past the window counts at the window's end, and a pair is put in order**,
    /// as [`ClassCaller::under`] names genotypes: at span 9, -12 is class 0 and 10 is class 18,
    /// and (5, -1) becomes (8, 14). Loci are flattened in order, an empty locus adding nothing.
    #[test]
    fn drawn_genotypes_in_window_folds_far_copies_to_the_end_classes_and_orders_the_pair() {
        let drawn = vec![vec![(-12, 3), (5, -1), (10, 10)], vec![], vec![(0, 0)]];
        assert_eq!(
            drawn_genotypes_in_window(&drawn, 9),
            vec![(0, 12), (8, 14), (18, 18), (9, 9)]
        );
    }
}
