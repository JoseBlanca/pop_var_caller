//! **The library half of the repeat-tract fitting bench**: the interface every fitting method
//! implements, the existing curves applied to any method's answers, and a held-out score that does
//! not depend on the method (`doc/devel/ng/research/str_fit_research_plan.md` §2, §5).
//!
//! The bench itself is `examples/ng_str_fit_bench.rs`: it draws a cohort with
//! [`cohort_draw`](super::cohort_draw), runs each method through this interface, and writes one row
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
//! cannot estimate is `None`, never a fitted zero. A method that draws its own curves across
//! repeat count says so with [`PeriodFit::draws_its_own_curves`]; otherwise the bench scores it
//! twice, as it comes (`_raw`) and after the curves the caller draws today (`_curves`,
//! [`with_the_curves`]).
//!
//! Two methods are wired in: [`TheClimb`] (`A`, the fit the caller ships) and [`TheCountingFit`]
//! (`B`, expectation maximisation per class with the spread read from pairs of samples).
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
//! **The climb's own 256-point average is not used**: at 13 allele lengths it is several
//! log-likelihood units a locus wrong (`fit_precision_quadrature_average_2026-10-01.md`).
//!
//! Every method's score at a locus uses the same random numbers, so differences between methods
//! carry less Monte Carlo noise than either score alone.

use super::cohort_draw::ClassTruth;
use super::*;

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
    /// **The reads the class's answer rests on**, which the curves weigh its own stutter by: the
    /// class's spanning reads, unless the method read only some of them.
    pub reads_crossing: u64,
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
            reads_crossing: evidence.spanning_reads(),
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
    /// The method's name in the bench's output, e.g. `A`.
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

    /// Whether the method draws its own curves across repeat count, so the caller's curves are not
    /// applied after it.
    fn draws_its_own_curves(&self) -> bool {
        false
    }
}

// ---------------------------------------------------------------------
// The two methods wired in
// ---------------------------------------------------------------------

/// **Method A: the climb**, as the caller runs it ([`fit_strata_on_sample_subsets`], so a cohort
/// of more than 256 samples is fitted on a subset grown until each class's level is measured to
/// 2%), with the curves left for the bench to apply.
#[derive(Debug, Clone)]
pub struct TheClimb {
    pub config: SsrFitConfig,
}

impl PeriodFit for TheClimb {
    fn name(&self) -> String {
        "A".to_owned()
    }

    fn fit_period(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
    ) -> Vec<ClassEstimate> {
        let mut config = self.config.clone();
        config.allele_span = allele_span;
        config.curve.draw_curves = false;
        // The drawn samples are exchangeable, so their own order is as good as the run's.
        let order: Vec<usize> = (0..homozygote_excess.len()).collect();
        let outcomes = fit_strata_on_sample_subsets(classes, homozygote_excess, &order, &config);
        classes
            .iter()
            .zip(outcomes)
            .map(|(evidence, outcome)| match outcome {
                StratumOutcome::Fitted(fit) => {
                    let value = |error: StratumError| error.value();
                    let errors =
                        fit.standard_errors
                            .as_ref()
                            .map_or_else(ClassErrors::default, |errors| {
                                let slippage = errors.slippage.first().copied().flatten();
                                ClassErrors {
                                    level: slippage.and_then(|group| value(group.level)),
                                    shorter_share: slippage
                                        .and_then(|group| value(group.shorter_share)),
                                    fall_off: slippage.and_then(|group| value(group.fall_off)),
                                    concentration: value(errors.concentration),
                                    allele_length_variation: errors
                                        .allele_length_variation
                                        .iter()
                                        .map(|error| value(*error))
                                        .collect(),
                                }
                            });
                    ClassEstimate {
                        stratum: evidence.stratum,
                        slippage: fit.slippage.first().copied().flatten(),
                        allele_length_variation: Some(fit.allele_length_variation.clone()),
                        concentration: Some(fit.concentration),
                        errors,
                        reads_crossing: fit.reads_crossing,
                    }
                }
                StratumOutcome::Derived(_) | StratumOutcome::Refused { .. } => {
                    ClassEstimate::absent(evidence)
                }
            })
            .collect()
    }
}

/// **Method B: the counting fit**, class by class
/// ([`fit_class_by_counting_with_pairs`](super::super::composite_fit::fit_class_by_counting_with_pairs)
/// with `pairs`, the within-sample version without), as `examples/ng_str_counting_fit.rs` runs it.
/// A class with fewer loci carrying reads than the climb's refusal floor is left absent, as that
/// example leaves it. It gives no standard errors.
#[derive(Debug, Clone, Copy)]
pub struct TheCountingFit {
    /// Read the spread off pairs of samples at a locus (what an inbred cohort needs).
    pub pairs: bool,
    pub max_rounds: usize,
    pub tolerance: f64,
    /// Rounds of spread-then-counting, at most, when `pairs`.
    pub alternations: usize,
}

impl Default for TheCountingFit {
    fn default() -> Self {
        Self {
            pairs: true,
            max_rounds: 2000,
            tolerance: 1e-7,
            alternations: 10,
        }
    }
}

impl PeriodFit for TheCountingFit {
    fn name(&self) -> String {
        "B".to_owned()
    }

    fn fit_period(
        &self,
        classes: &[StratumEvidence],
        homozygote_excess: &[f64],
        allele_span: i32,
    ) -> Vec<ClassEstimate> {
        let config = SsrFitConfig {
            allele_span,
            ..SsrFitConfig::default()
        };
        let start = Slippage {
            level: 0.05,
            shorter_share: 0.5,
            fall_off: 0.3,
        };
        // The classes are independent, so they are fitted side by side; each also spreads its
        // loci over the pool.
        classes
            .par_iter()
            .map(|evidence| {
                if evidence.tracts_with_reads() < DEFAULT_REFUSAL_FLOOR {
                    return ClassEstimate::absent(evidence);
                }
                let fit = if self.pairs {
                    super::super::composite_fit::fit_class_by_counting_with_pairs(
                        evidence,
                        start,
                        homozygote_excess,
                        &config,
                        self.max_rounds,
                        self.tolerance,
                        self.alternations,
                    )
                } else {
                    super::super::composite_fit::fit_class_by_counting(
                        evidence,
                        start,
                        homozygote_excess,
                        &config,
                        self.max_rounds,
                        self.tolerance,
                    )
                };
                ClassEstimate {
                    stratum: evidence.stratum,
                    slippage: Some(fit.slippage),
                    allele_length_variation: Some(fit.allele_length_variation),
                    concentration: Some(fit.concentration),
                    errors: ClassErrors::default(),
                    reads_crossing: evidence.spanning_reads(),
                }
            })
            .collect()
    }
}

/// **Method D: penalised counting**, every class of the motif length at once, each number a
/// penalised curve across repeat count
/// ([`fit_period_by_penalised_counting`](super::super::penalised_fit::fit_period_by_penalised_counting)),
/// the spread read off pairs of samples as method B reads it. It draws its own curves, fits every
/// class with reads (no refusal floor: a thin class takes what its counts and its neighbours say),
/// and gives no standard errors.
///
/// **The concentration is absent when no locus has two samples with reads**, as the project
/// reports a number it cannot estimate. With `report_working_concentration` it reports the value
/// the fit worked at instead, so that the held-out score, which needs every number, can be computed
/// at one sample; that variant is `D_working_F` and its concentration is not an estimate.
#[derive(Debug, Clone)]
pub struct ThePenalisedFit {
    pub config: super::super::penalised_fit::PenalisedFitConfig,
    pub report_working_concentration: bool,
}

impl Default for ThePenalisedFit {
    fn default() -> Self {
        Self {
            config: super::super::penalised_fit::PenalisedFitConfig {
                start: Slippage {
                    level: 0.05,
                    shorter_share: 0.5,
                    fall_off: 0.3,
                },
                max_rounds: 2000,
                tolerance: 1e-7,
                pair_alternations: 10,
                ..super::super::penalised_fit::PenalisedFitConfig::default()
            },
            report_working_concentration: false,
        }
    }
}

impl PeriodFit for ThePenalisedFit {
    fn name(&self) -> String {
        if self.report_working_concentration {
            "D_working_F".to_owned()
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
        let config = SsrFitConfig {
            allele_span,
            ..SsrFitConfig::default()
        };
        let fitted = super::super::penalised_fit::fit_period_by_penalised_counting(
            classes,
            homozygote_excess,
            &config,
            &self.config,
        );
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
                    errors: ClassErrors::default(),
                    reads_crossing: evidence.spanning_reads(),
                },
            })
            .collect()
    }

    fn draws_its_own_curves(&self) -> bool {
        true
    }
}

// ---------------------------------------------------------------------
// The curves, after any method
// ---------------------------------------------------------------------

/// **Any method's answers with the caller's curves applied**, exactly as [`fit_strata`] applies
/// them to the climb's ([`draw_and_apply_curves`]): each class's stutter blended with its period's
/// curves by how many of its reads slipped, and a class with no answer of its own given the curves'
/// stutter (its shares and concentration stay absent). The shares and concentration of a class
/// that has them are not touched.
///
/// **A number the curves moved loses its standard error**: a blended number has none
/// (`str_slippage_level_curve.md` §5.2).
///
/// # Panics
///
/// When `raw` and `classes` are not the same classes in the same order.
#[must_use]
pub fn with_the_curves(
    raw: &[ClassEstimate],
    classes: &[StratumEvidence],
    config: &SsrFitConfig,
) -> Vec<ClassEstimate> {
    assert_eq!(raw.len(), classes.len(), "one estimate per class");
    let mut outcomes: Vec<StratumOutcome> = raw
        .iter()
        .zip(classes)
        .map(|(estimate, evidence)| {
            assert_eq!(
                estimate.stratum, evidence.stratum,
                "estimates in class order"
            );
            let tracts = evidence.tracts_with_reads();
            let Some(slippage) = estimate.slippage else {
                return StratumOutcome::Refused {
                    stratum: evidence.stratum,
                    tracts,
                    reason: if tracts == 0 {
                        StratumRefusal::NoSpanningReads
                    } else {
                        StratumRefusal::BelowTheFloor {
                            tracts,
                            floor: config.refusal_floor,
                        }
                    },
                };
            };
            let slipped_reads = slippage.level * estimate.reads_crossing as f64;
            StratumOutcome::Fitted(Box::new(StratumFit {
                stratum: evidence.stratum,
                slippage: vec![Some(slippage)],
                allele_length_variation: estimate
                    .allele_length_variation
                    .clone()
                    .unwrap_or_default(),
                concentration: estimate.concentration.unwrap_or(f64::NAN),
                log_likelihood_a_tract: f64::NAN,
                tracts_fitted: tracts,
                borrowed: Vec::new(),
                // Not read by the curves; a method that does not climb has no walk to report.
                ending: ClimbEnding::Settled,
                walks: Vec::new(),
                samples_fitted_on: None,
                tracts_of_its_own: tracts,
                reads_crossing: estimate.reads_crossing,
                level_provenance: vec![Some(LevelProvenance {
                    source: LevelSource::Cell,
                    curve: None,
                    reach: None,
                    slipped_reads: Some(slipped_reads),
                })],
                shares_provenance: vec![Some(SharesProvenance::own(slipped_reads))],
                standard_errors: None,
            }))
        })
        .collect();
    draw_and_apply_curves(&mut outcomes, classes, config);
    raw.iter()
        .zip(outcomes)
        .map(|(estimate, outcome)| {
            let smoothed = outcome.slippage().first().copied().flatten();
            let kept = |own: Option<f64>, after: Option<f64>, error: Option<f64>| {
                (own.is_some() && own.map(f64::to_bits) == after.map(f64::to_bits))
                    .then_some(error)
                    .flatten()
            };
            let own = estimate.slippage;
            ClassEstimate {
                slippage: smoothed,
                errors: ClassErrors {
                    level: kept(
                        own.map(|s| s.level),
                        smoothed.map(|s| s.level),
                        estimate.errors.level,
                    ),
                    shorter_share: kept(
                        own.map(|s| s.shorter_share),
                        smoothed.map(|s| s.shorter_share),
                        estimate.errors.shorter_share,
                    ),
                    fall_off: kept(
                        own.map(|s| s.fall_off),
                        smoothed.map(|s| s.fall_off),
                        estimate.errors.fall_off,
                    ),
                    ..estimate.errors.clone()
                },
                ..estimate.clone()
            }
        })
        .collect()
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
        tracts_over_guard_threshold: evidence.tracts_over_guard_threshold,
        reads_reaching_not_crossing: evidence.reads_reaching_not_crossing,
        guard_reads: evidence.guard_reads,
        bases_compared: evidence.bases_compared,
        mismatching_bases: evidence.mismatching_bases,
    }
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

    /// **The curves the bench applies to the climb's raw answers are the ones [`fit_strata`] applies
    /// itself**, so `A_curves` is the climb as the caller runs it: five classes of one period, one
    /// of them below the refusal floor so that it is furnished from the curves.
    #[test]
    fn the_curves_after_the_climb_are_the_climbs_own() {
        let strata: Vec<StratumEvidence> = (0..5_u64)
            .map(|index| {
                let mut evidence = bench_fixtures::draw_stratum(
                    Slippage {
                        level: 0.01 * (index + 1) as f64,
                        shorter_share: 0.6,
                        fall_off: 0.3,
                    },
                    &spectrum_of(3),
                    1.0,
                    0.5,
                    if index == 2 { 5 } else { 15 },
                    6,
                    3,
                    1,
                    100 + index,
                );
                evidence.stratum = Stratum {
                    period: 2,
                    reference_repeats: 6 + index,
                };
                evidence
            })
            .collect();
        let excess = vec![0.5; 6];
        // What is tested is the curves' plumbing, not how well the climb fits, so a short climb
        // from one start does.
        let config = SsrFitConfig {
            allele_span: 1,
            starting_points: vec![StartingPoint::spanning_the_monomorphic_range()[1]],
            max_rounds: 5,
            ..SsrFitConfig::default()
        };
        let by_the_fit = fit_strata(&strata, &excess, &config);
        let raw = TheClimb {
            config: config.clone(),
        }
        .fit_period(&strata, &excess, 1);
        let by_the_bench = with_the_curves(&raw, &strata, &config);
        assert!(
            raw[2].slippage.is_none(),
            "the thin class has no answer of its own"
        );
        assert!(
            by_the_bench[2].slippage.is_some(),
            "and the curves furnish one"
        );
        for (outcome, estimate) in by_the_fit.iter().zip(&by_the_bench) {
            assert_eq!(
                outcome.slippage().first().copied().flatten(),
                estimate.slippage
            );
        }
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
}
