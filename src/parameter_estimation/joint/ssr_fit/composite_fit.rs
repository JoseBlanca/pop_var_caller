//! **A repeat class's stutter and allele-length shares by expectation maximisation, without the
//! average over tract frequencies** — the owner's routes 1 and 2 together (2026-10-11), a prototype.
//!
//! The fit proper averages each tract's likelihood over how its alleles could be shared among its
//! chromosomes, and that average, over 256 points, is several log-likelihood units a tract wrong at
//! thirteen allele lengths (`fit_precision_quadrature_average_2026-10-01.md`). It also makes the fit
//! slow. This estimates the same stutter and the class's allele-length shares from a simpler
//! objective: a sample's alleles at a random tract of the class are drawn from the class's shares,
//! whatever that tract's own frequencies, so
//!
//! `Σ over (tract, sample) of ln Σ_g prior(g | shares, F) L(g | stutter)`
//!
//! is maximised by the class's numbers. It ignores that samples at one tract share that tract's
//! frequencies — which costs precision, not consistency — and so says nothing about how much tracts
//! differ (the concentration).
//!
//! Each round shares every read between "read exactly" and "slipped this far, this way", given each
//! sample's genotype probabilities, and every sample's alleles among the lengths; the slip rate, the
//! shorter share, the fall-off and the shares then follow from the expected counts.

use super::*;

/// What [`fit_class_by_counting`] returns.
#[derive(Debug, Clone, PartialEq)]
pub struct CountedClass {
    /// The class's stutter.
    pub slippage: Slippage,
    /// How the class's chromosomes are spread over the allele lengths, `-allele_span ..= +allele_span`.
    pub allele_length_variation: Vec<f64>,
    /// How much the class's tracts differ from one another: the Dirichlet concentration, read off
    /// how much more often a sample's two copies match than its homozygote excess alone explains.
    pub concentration: f64,
    /// How many rounds ran, and the largest change of any number in the last.
    pub rounds: usize,
    pub last_change: f64,
    /// The objective at the answer.
    pub composite_log_likelihood: f64,
}

/// The expected counts one round gathers over a class, or over one tract of it.
///
/// **They are also the objective's slope.** The slope of the log of a sum over hidden genotypes
/// and slips is the expected slope of the log-likelihood that knew them (Fisher's identity), and
/// that slope is a function of these counts alone (`counting_fit_errors`).
#[derive(Debug, Clone, Default)]
pub(super) struct Counts {
    /// Per allele class, how many independent allele draws.
    pub(super) draws: Vec<f64>,
    /// Of each sample's two copies not identical by descent, how many pairs, and how many of those
    /// were one draw from the tract's frequencies.
    pub(super) pairs_not_by_descent: f64,
    pub(super) shared_in_tract: f64,
    /// Reads that reported their allele's length, and reads that slipped.
    pub(super) exact: f64,
    pub(super) slipped: f64,
    /// Slipped reads by direction; `steps` is the repeats they slipped summed, over the `stepped`
    /// reads whose direction and length were placed.
    pub(super) shorter: f64,
    pub(super) longer: f64,
    pub(super) steps: f64,
    pub(super) stepped: f64,
    /// The objective summed over what was counted.
    pub(super) log_likelihood: f64,
}

impl Counts {
    pub(super) fn new(classes: usize) -> Self {
        Self {
            draws: vec![0.0; classes],
            ..Self::default()
        }
    }

    pub(super) fn add(&mut self, other: &Counts) {
        for (into, from) in self.draws.iter_mut().zip(&other.draws) {
            *into += from;
        }
        self.exact += other.exact;
        self.pairs_not_by_descent += other.pairs_not_by_descent;
        self.shared_in_tract += other.shared_in_tract;
        self.slipped += other.slipped;
        self.shorter += other.shorter;
        self.longer += other.longer;
        self.steps += other.steps;
        self.stepped += other.stepped;
        self.log_likelihood += other.log_likelihood;
    }
}

/// The slipped reads an allele sends into one read bucket, by direction: each direction's
/// probability among all reads, and `steps`, the sum over both of probability times the repeats
/// slipped.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct SlipsIntoBucket {
    shorter: f64,
    longer: f64,
    steps: f64,
}

/// **Where an allele's slips land, bucket by bucket**, in the terms of
/// [`Slippage::read_probabilities`]: the slipped part of each bucket's probability there,
/// split by direction, with the expected steps.
///
/// **Measured from the allele's own offset, not from its bucket.** An allele past the recorded
/// range shares its end bucket with the allele at the edge, and that bucket also holds its
/// shorter slips and its one-step longer slip. A slip into an end bucket went at least that far;
/// a geometric fall-off has no memory, so its expected length is the least step that reaches
/// the bucket plus `f / (1 − f)`.
fn slips_into_each_bucket(
    slippage: Slippage,
    allele_offset: i32,
    read_span: i32,
) -> Vec<SlipsIntoBucket> {
    let mut out = vec![SlipsIntoBucket::default(); (2 * read_span + 1) as usize];
    let slot = |offset: i32| (offset.clamp(-read_span, read_span) + read_span) as usize;
    let fall_off = slippage.fall_off;
    let beyond_the_least_step = fall_off / (1.0 - fall_off);
    for (direction, share) in [
        (-1_i32, slippage.shorter_share),
        (1, 1.0 - slippage.shorter_share),
    ] {
        let mut add = |bucket: usize, probability: f64, steps: f64| {
            let into = &mut out[bucket];
            if direction < 0 {
                into.shorter += probability;
            } else {
                into.longer += probability;
            }
            into.steps += probability * steps;
        };
        let inside_steps = if direction < 0 {
            (allele_offset + read_span).max(0)
        } else {
            (read_span - allele_offset).max(0)
        };
        for step in 1..=inside_steps {
            let weight = (1.0 - fall_off) * float::powi(fall_off, step - 1);
            add(
                slot(allele_offset + direction * step),
                slippage.level * share * weight,
                f64::from(step),
            );
        }
        let tail = float::powi(fall_off, inside_steps);
        add(
            slot(allele_offset + direction * (inside_steps + 1)),
            slippage.level * share * tail,
            f64::from(inside_steps + 1) + beyond_the_least_step,
        );
    }
    out
}

/// **Fits one class's stutter and shares by counting**, from `start`'s stutter and flat shares,
/// until no number moves by more than `tolerance` or `max_rounds` pass. Every library is scored at
/// the class's stutter (multipliers are not fitted here).
#[must_use]
pub fn fit_class_by_counting(
    evidence: &StratumEvidence,
    start: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
) -> CountedClass {
    count_class(
        evidence,
        start,
        None,
        None,
        1.0,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
    )
}

/// The rounds of counting, from `start` and `start_shares` (flat when `None`). With `held_spread`
/// the spread is held there rather than counted from each sample's two copies, which otherwise
/// starts at `start_spread`.
#[allow(clippy::too_many_arguments, reason = "a prototype's internal entry")]
fn count_class(
    evidence: &StratumEvidence,
    start: Slippage,
    start_shares: Option<&[f64]>,
    held_spread: Option<f64>,
    start_spread: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
) -> CountedClass {
    let classes = allele_classes_of(config);
    let mut slippage = start;
    let mut shares =
        start_shares.map_or_else(|| vec![1.0 / classes as f64; classes], <[f64]>::to_vec);
    let mut concentration = held_spread.unwrap_or(start_spread);
    let mut rounds = 0;
    let mut last_change = f64::INFINITY;
    let mut log_likelihood = f64::NEG_INFINITY;
    while rounds < max_rounds {
        rounds += 1;
        let per_tract = expected_counts_by_tract(
            evidence,
            slippage,
            &shares,
            concentration,
            homozygote_excess,
            config,
        );
        let mut counts = Counts::new(classes);
        for tract in &per_tract {
            counts.add(tract);
        }
        log_likelihood = counts.log_likelihood;
        let draws: f64 = counts.draws.iter().sum();
        let new_shares: Vec<f64> = counts.draws.iter().map(|d| d / draws).collect();
        let new_concentration = if let Some(held) = held_spread {
            held
        } else if counts.shared_in_tract > 0.0 {
            (counts.pairs_not_by_descent / counts.shared_in_tract - 1.0).clamp(1e-3, 1e4)
        } else {
            1e4
        };
        let new_slippage = Slippage {
            level: counts.slipped / (counts.slipped + counts.exact),
            shorter_share: if counts.shorter + counts.longer > 0.0 {
                counts.shorter / (counts.shorter + counts.longer)
            } else {
                slippage.shorter_share
            },
            // A step k has probability (1 − f) f^(k − 1), whose mean is 1 / (1 − f).
            fall_off: if counts.steps > 0.0 {
                (1.0 - counts.stepped / counts.steps).clamp(0.0, 0.99)
            } else {
                slippage.fall_off
            },
        };
        last_change = new_shares
            .iter()
            .zip(&shares)
            .map(|(new, old)| (new - old).abs())
            .chain([
                (new_slippage.level - slippage.level).abs(),
                (new_slippage.shorter_share - slippage.shorter_share).abs(),
                (new_slippage.fall_off - slippage.fall_off).abs(),
                (float::ln(new_concentration) - float::ln(concentration)).abs() * 0.01,
            ])
            .fold(0.0, f64::max);
        shares = new_shares;
        slippage = new_slippage;
        concentration = new_concentration;
        if last_change < tolerance {
            break;
        }
    }
    CountedClass {
        slippage,
        allele_length_variation: shares,
        concentration,
        rounds,
        last_change,
        composite_log_likelihood: log_likelihood,
    }
}

/// **One round's expected counts, tract by tract**, in the order of `evidence.tracts`, at the
/// stutter `slippage`, the allele-length shares `shares` and the spread `concentration`.
///
/// Weighs every genotype each sample could carry by its prior times how well it explains the
/// sample's reads, then shares each read between "reported its allele's length" and "slipped this
/// far, this way", and each sample's alleles among the lengths.
pub(super) fn expected_counts_by_tract(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    concentration: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<Counts> {
    let classes = allele_classes_of(config);
    let genotypes = genotype_pairs(classes);
    let span = config.allele_span;
    let read_span = evidence.read_span;
    let last_bucket = (2 * read_span) as usize;
    let slot =
        |class: usize| ((class as i32 - span).clamp(-read_span, read_span) + read_span) as usize;
    let libraries = one_stutter_for_every_library(evidence);
    let distributions = ReadDistributionsByLevel::of(&[slippage], &libraries, span, read_span);
    let per_allele = &distributions.read_distribution_by_allele[0];
    let stay = 1.0 - slippage.level;
    let slips_by_bucket: Vec<Vec<SlipsIntoBucket>> = (0..classes)
        .map(|class| slips_into_each_bucket(slippage, class as i32 - span, read_span))
        .collect();
    // Each genotype's prior, given a sample's homozygote excess, split into the part where its two
    // copies are one draw (identical by descent) and the part where they are two.
    let one_draw: Vec<f64> = genotypes
        .iter()
        .map(|(a, b)| if a == b { shares[*a] } else { 0.0 })
        .collect();
    let two_draws: Vec<f64> = genotypes
        .iter()
        .map(|(a, b)| {
            if a == b {
                shares[*a] * shares[*a]
            } else {
                2.0 * shares[*a] * shares[*b]
            }
        })
        .collect();
    evidence
        .tracts
        .par_iter()
        .map(|tract| {
            let mut counts = Counts::new(classes);
            for sample in &tract.samples {
                let excess = homozygote_excess[sample.sample as usize];
                let reads = pooled_reads(sample, last_bucket + 1);
                // Each genotype's log-likelihood, then its posterior.
                let log_likelihoods: Vec<f64> = genotypes
                    .iter()
                    .map(|(a, b)| {
                        let mut total = 0.0;
                        add_ln_reads(&mut total, per_allele, &reads, *a, *b);
                        total
                    })
                    .collect();
                let largest = log_likelihoods
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                // A sample's two copies are one draw if identical by descent, or, failing that, if they
                // drew the same chromosome lineage from the tract's own frequencies: chance
                // 1 / (1 + concentration) under a Dirichlet.
                let from_the_tract = (1.0 - excess) / (1.0 + concentration);
                let one_draw_chance = excess + from_the_tract;
                let two_draw_chance = (1.0 - excess) * concentration / (1.0 + concentration);
                let joint_one: Vec<f64> = log_likelihoods
                    .iter()
                    .zip(&one_draw)
                    .map(|(ll, prior)| one_draw_chance * prior * float::exp(ll - largest))
                    .collect();
                let joint_two: Vec<f64> = log_likelihoods
                    .iter()
                    .zip(&two_draws)
                    .map(|(ll, prior)| two_draw_chance * prior * float::exp(ll - largest))
                    .collect();
                let total: f64 = joint_one.iter().sum::<f64>() + joint_two.iter().sum::<f64>();
                if total <= 0.0 {
                    continue;
                }
                counts.log_likelihood += float::ln(total) + largest;
                for (index, (a, b)) in genotypes.iter().enumerate() {
                    let one = joint_one[index] / total;
                    let two = joint_two[index] / total;
                    let weight = one + two;
                    if weight < 1e-12 {
                        continue;
                    }
                    counts.draws[*a] += one + two;
                    counts.draws[*b] += two;
                    let shared = one * from_the_tract / one_draw_chance.max(1e-300);
                    counts.shared_in_tract += shared;
                    counts.pairs_not_by_descent += shared + two;
                    for (bucket, count) in reads.iter().enumerate() {
                        if *count == 0 {
                            continue;
                        }
                        let reported = 0.5 * (per_allele[*a][bucket] + per_allele[*b][bucket]);
                        if reported <= 0.0 {
                            continue;
                        }
                        for allele in [*a, *b] {
                            let from_allele = 0.5 * per_allele[allele][bucket] / reported;
                            let share = weight * f64::from(*count) * from_allele;
                            let unslipped = if slot(allele) == bucket {
                                stay / per_allele[allele][bucket]
                            } else {
                                0.0
                            };
                            let exact = share * unslipped.min(1.0);
                            let slipped = share - exact;
                            counts.exact += exact;
                            counts.slipped += slipped;
                            // Which way the slipped reads went and how far, split by the
                            // model's own terms: an end bucket can hold slips of both
                            // directions from an allele at or past the recorded range.
                            let into = slips_by_bucket[allele][bucket];
                            let landed = into.shorter + into.longer;
                            if slipped > 0.0 && landed > 0.0 {
                                counts.shorter += slipped * into.shorter / landed;
                                counts.longer += slipped * into.longer / landed;
                                counts.steps += slipped * into.steps / landed;
                                counts.stepped += slipped;
                            }
                        }
                    }
                }
            }
            counts
        })
        .collect()
}

/// Every library of `evidence` scored at the class's one stutter: one slippage group, every
/// multiplier one. The counting fit does not fit libraries apart.
fn one_stutter_for_every_library(evidence: &StratumEvidence) -> Vec<Library> {
    evidence
        .libraries
        .iter()
        .map(|library| Library {
            slippage_group: 0,
            stutter_multiplier: 1.0,
            ..*library
        })
        .collect()
}

/// One sample's reads at a tract, its libraries summed bucket by bucket.
fn pooled_reads(sample: &SampleTractReads, buckets: usize) -> Vec<u32> {
    let mut reads = vec![0_u32; buckets];
    for (_, library_reads) in &sample.by_library {
        for (into, count) in reads.iter_mut().zip(library_reads) {
            *into += count;
        }
    }
    reads
}

/// **The class fitted by counting, its spread read off pairs of samples at a tract** — what an
/// inbred cohort needs, where a sample's two copies are nearly always one draw and so say almost
/// nothing about how much tracts differ.
///
/// Alternates: the stutter and shares counted with the spread held, then the spread that maximises
/// the pairwise likelihood ([`fit_spread_by_pairs`]) with the rest held, `alternations` times or
/// until the spread moves by less than 1% — starting from the within-sample counting fit.
#[must_use]
pub fn fit_class_by_counting_with_pairs(
    evidence: &StratumEvidence,
    start: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    alternations: usize,
) -> CountedClass {
    let fitted = fit_class_by_counting(
        evidence,
        start,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
    );
    alternate_with_pairs(
        evidence,
        fitted,
        1.0,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
        alternations,
    )
}

/// **Where a class's concentration is read from**: which of the counting fit's two estimators made
/// a [`CountedClass`], so that its errors and its refits use the same one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConcentrationSource {
    /// From how often each sample's own two copies match ([`fit_class_by_counting`]). Says almost
    /// nothing in an inbred cohort, where the two copies are nearly always one draw.
    EachSamplesTwoCopies,
    /// From pairs of samples at a tract ([`fit_class_by_counting_with_pairs`]). Says nothing at one
    /// sample, where there is no pair.
    PairsOfSamples,
}

/// **The class refitted from an earlier answer `from`** rather than from flat shares: what a
/// bootstrap replicate needs, whose answer lies near the full fit's. Uses the estimator `source`
/// names, with the same rounds, tolerance and alternations as a fit from scratch.
#[allow(
    clippy::too_many_arguments,
    reason = "the counting fit's own settings, passed through"
)]
#[must_use]
pub fn refit_class_by_counting_from(
    evidence: &StratumEvidence,
    from: &CountedClass,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    alternations: usize,
) -> CountedClass {
    match source {
        ConcentrationSource::EachSamplesTwoCopies => count_class(
            evidence,
            from.slippage,
            Some(&from.allele_length_variation),
            None,
            from.concentration,
            homozygote_excess,
            config,
            max_rounds,
            tolerance,
        ),
        ConcentrationSource::PairsOfSamples => alternate_with_pairs(
            evidence,
            from.clone(),
            from.concentration,
            homozygote_excess,
            config,
            max_rounds,
            tolerance,
            alternations,
        ),
    }
}

/// Alternates, from `fitted` and the spread `spread`: the spread that maximises the pairwise
/// likelihood with the rest held, then the stutter and shares counted with that spread held,
/// `alternations` times or until the spread moves by less than 1%.
#[allow(clippy::too_many_arguments, reason = "a prototype's internal entry")]
fn alternate_with_pairs(
    evidence: &StratumEvidence,
    mut fitted: CountedClass,
    mut spread: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    alternations: usize,
) -> CountedClass {
    for _ in 0..alternations {
        let next = fit_spread_by_pairs(
            evidence,
            fitted.slippage,
            &fitted.allele_length_variation,
            homozygote_excess,
            config,
        );
        fitted = count_class(
            evidence,
            fitted.slippage,
            Some(&fitted.allele_length_variation),
            Some(next),
            1.0,
            homozygote_excess,
            config,
            max_rounds,
            tolerance,
        );
        let settled = (float::ln(next) - float::ln(spread)).abs() < 0.01;
        spread = next;
        if settled {
            break;
        }
    }
    fitted
}

/// **The spread that maximises the pairwise likelihood of the class**, stutter and shares held.
///
/// At a tract, every sample with reads is paired with the next (and the last with the first). A
/// pair's likelihood sums both samples' reads over their alleles, each sample's copies one draw with
/// chance its homozygote excess and two otherwise, the draws from the tract's frequencies weighted
/// as a Pólya urn: a draw takes length `a` with chance `(spread · share_a + earlier draws of a) /
/// (spread + earlier draws)`. Searched on the log scale between 10⁻³ and 10⁴.
#[must_use]
pub fn fit_spread_by_pairs(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> f64 {
    let classes = allele_classes_of(config);
    let per_tract = pair_rows_by_tract(evidence, slippage, homozygote_excess, config);
    let objective = |ln_spread: f64| -> f64 {
        let spread = float::exp(ln_spread);
        per_tract
            .par_iter()
            .map(|samples| scaled_ln_pairs_of_tract(samples, shares, spread, classes))
            .sum()
    };
    // A coarse grid, then golden-section search around its best point.
    let grid: Vec<f64> = (0..=28)
        .map(|step| float::ln(1e-3) + f64::from(step) * (float::ln(1e4) - float::ln(1e-3)) / 28.0)
        .collect();
    let scores: Vec<f64> = grid.iter().map(|x| objective(*x)).collect();
    let best = scores
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1))
        .map_or(0, |(index, _)| index);
    let (mut low, mut high) = (
        grid[best.saturating_sub(1)],
        grid[(best + 1).min(grid.len() - 1)],
    );
    let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut left, mut right) = (high - ratio * (high - low), low + ratio * (high - low));
    let (mut left_score, mut right_score) = (objective(left), objective(right));
    for _ in 0..30 {
        if left_score > right_score {
            high = right;
            right = left;
            right_score = left_score;
            left = high - ratio * (high - low);
            left_score = objective(left);
        } else {
            low = left;
            left = right;
            left_score = right_score;
            right = low + ratio * (high - low);
            right_score = objective(right);
        }
    }
    float::exp(0.5 * (low + high))
}

/// One sample's reads at a tract, in the form a pair's likelihood is summed from.
pub(super) struct PairRow {
    /// The likelihood of the sample's reads under every **ordered** allele pair, `first ·
    /// classes + second`, divided by the largest of them.
    pub(super) ordered: Vec<f64>,
    /// The sample's homozygote excess.
    pub(super) excess: f64,
}

/// **Every tract's samples as [`PairRow`]s at the stutter `slippage`**, in the order of
/// `evidence.tracts`; a tract with fewer than two samples has none, since it has no pair.
pub(super) fn pair_rows_by_tract(
    evidence: &StratumEvidence,
    slippage: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<Vec<PairRow>> {
    let classes = allele_classes_of(config);
    let genotypes = genotype_pairs(classes);
    let buckets = (2 * evidence.read_span + 1) as usize;
    let libraries = one_stutter_for_every_library(evidence);
    let distributions = ReadDistributionsByLevel::of(
        &[slippage],
        &libraries,
        config.allele_span,
        evidence.read_span,
    );
    let per_allele = &distributions.read_distribution_by_allele[0];
    evidence
        .tracts
        .par_iter()
        .map(|tract| {
            if tract.samples.len() < 2 {
                return Vec::new();
            }
            tract
                .samples
                .iter()
                .map(|sample| {
                    let reads = pooled_reads(sample, buckets);
                    let log_likelihoods: Vec<f64> = genotypes
                        .iter()
                        .map(|(first, second)| {
                            let mut total = 0.0;
                            add_ln_reads(&mut total, per_allele, &reads, *first, *second);
                            total
                        })
                        .collect();
                    let largest = log_likelihoods
                        .iter()
                        .copied()
                        .fold(f64::NEG_INFINITY, f64::max);
                    let mut ordered = vec![0.0; classes * classes];
                    for ((first, second), value) in genotypes.iter().zip(&log_likelihoods) {
                        let scaled = float::exp(value - largest);
                        ordered[first * classes + second] = scaled;
                        ordered[second * classes + first] = scaled;
                    }
                    PairRow {
                        ordered,
                        excess: *homozygote_excess
                            .get(sample.sample as usize)
                            .expect("every sample in a tract has a homozygote excess"),
                    }
                })
                .collect()
        })
        .collect()
}

/// **One tract's pairwise log-likelihood, less its rows' offsets**: every sample paired with the
/// next and the last with the first (one pair when there are two), at the allele-length `shares`
/// and the spread `spread`. The rows' offsets are left out: they do not depend on the shares or
/// the spread.
pub(super) fn scaled_ln_pairs_of_tract(
    samples: &[PairRow],
    shares: &[f64],
    spread: f64,
    classes: usize,
) -> f64 {
    let count = samples.len();
    if count < 2 {
        return 0.0;
    }
    let weights: Vec<f64> = shares.iter().map(|share| spread * share).collect();
    let pairs = if count == 2 { 1 } else { count };
    let mut total = 0.0;
    for index in 0..pairs {
        let one = &samples[index];
        let other = &samples[(index + 1) % count];
        let pair = pair_likelihood(
            &one.ordered,
            one.excess,
            &other.ordered,
            other.excess,
            &weights,
            spread,
            classes,
        );
        total += float::ln(pair.max(1e-300));
    }
    total
}

/// One pair of samples' likelihood at a tract: both samples' reads summed over their alleles, each
/// sample one draw (chance its excess) or two, the draws weighted as a Pólya urn.
///
/// **Summed in closed form, a pair in `classes²` steps rather than `classes⁴`.** With `w = spread ·
/// share`, the urn's weight for a later draw is `w_x` plus how many earlier draws were `x`, so every
/// sum over the later sample's alleles reduces to a few totals of its row: `Q = Σ L(c,d) w_c w_d`,
/// `R = Σ L(c,c) w_c`, `u_x = Σ_c L(c,x) w_c` and the diagonal `D_x = L(x,x)`. Checked against the
/// four nested sums in `the_closed_form_pair_likelihood_is_the_four_nested_sums`.
fn pair_likelihood(
    one: &[f64],
    one_excess: f64,
    other: &[f64],
    other_excess: f64,
    weights: &[f64],
    spread: f64,
    classes: usize,
) -> f64 {
    let totals = |row: &[f64]| {
        let diagonal: Vec<f64> = (0..classes).map(|x| row[x * classes + x]).collect();
        let u: Vec<f64> = (0..classes)
            .map(|x| {
                (0..classes)
                    .map(|c| row[c * classes + x] * weights[c])
                    .sum()
            })
            .collect();
        let q: f64 = (0..classes).map(|d| u[d] * weights[d]).sum();
        let r: f64 = (0..classes).map(|c| diagonal[c] * weights[c]).sum();
        (diagonal, u, q, r)
    };
    let (d1, u1, q1, r1) = totals(one);
    let (d2, u2, q2, r2) = totals(other);
    let alpha = spread;
    let two = alpha * (alpha + 1.0);
    let three = two * (alpha + 2.0);
    let four = three * (alpha + 3.0);

    // Each sample one draw: Σ_a Σ_c h1_a h2_c w_a (w_c + [c = a]).
    let one_one = ((0..classes).map(|a| d1[a] * weights[a]).sum::<f64>()
        * (0..classes).map(|c| d2[c] * weights[c]).sum::<f64>()
        + (0..classes)
            .map(|a| d1[a] * d2[a] * weights[a])
            .sum::<f64>())
        / two;
    // The first sample one draw `a`, then the other's two: Σ_a h1_a w_a (Q2 + R2 + 2 u2_a + 2 D2_a).
    let one_two = (0..classes)
        .map(|a| d1[a] * weights[a] * (q2 + r2 + 2.0 * u2[a] + 2.0 * d2[a]))
        .sum::<f64>()
        / three;
    // The first sample two draws, then the other's one `a`: Σ_a h2_a w_a (Q1 + R1 + 2 u1_a + 2 D1_a).
    let two_one = (0..classes)
        .map(|a| d2[a] * weights[a] * (q1 + r1 + 2.0 * u1[a] + 2.0 * d1[a]))
        .sum::<f64>()
        / three;
    // Each sample two draws: Σ_ab L1(a,b) w_a (w_b + [b = a]) ·
    //   (Q2 + R2 + 2 (u2_a + u2_b) + 2 (D2_a + D2_b) + 2 L2(a,b)).
    let mut two_two = 0.0;
    for a in 0..classes {
        for b in 0..classes {
            let first = one[a * classes + b];
            if first == 0.0 {
                continue;
            }
            let urn = weights[a] * (weights[b] + f64::from(u8::from(a == b)));
            let rest = q2
                + r2
                + 2.0 * (u2[a] + u2[b])
                + 2.0 * (d2[a] + d2[b])
                + 2.0 * other[a * classes + b];
            two_two += first * urn * rest;
        }
    }
    two_two /= four;
    one_excess * other_excess * one_one
        + one_excess * (1.0 - other_excess) * one_two
        + (1.0 - one_excess) * other_excess * two_one
        + (1.0 - one_excess) * (1.0 - other_excess) * two_two
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};

    /// One pair of samples' likelihood at a tract: both samples' reads summed over their alleles, each
    /// sample one draw (chance its excess) or two, the draws weighted as a Pólya urn.
    fn pair_likelihood_by_nested_sums(
        one: &[f64],
        one_excess: f64,
        other: &[f64],
        other_excess: f64,
        weights: &[f64],
        spread: f64,
        classes: usize,
    ) -> f64 {
        // The chance of each ordered sequence of draws, built draw by draw.
        let urn = |drawn: &[usize]| -> f64 {
            let mut chance = 1.0;
            for (position, allele) in drawn.iter().enumerate() {
                let earlier = drawn[..position]
                    .iter()
                    .filter(|previous| *previous == allele)
                    .count() as f64;
                chance *= (weights[*allele] + earlier) / (spread + position as f64);
            }
            chance
        };
        let homozygous = |row: &[f64], allele: usize| row[allele * classes + allele];
        let (mut one_one, mut one_two, mut two_one, mut two_two) = (0.0, 0.0, 0.0, 0.0);
        for a in 0..classes {
            let one_a = homozygous(one, a);
            for c in 0..classes {
                // Each sample one draw.
                one_one += one_a * homozygous(other, c) * urn(&[a, c]);
                for d in 0..classes {
                    // The first one draw, the second two; and the other way round.
                    one_two += one_a * other[c * classes + d] * urn(&[a, c, d]);
                    two_one += one[c * classes + d] * homozygous(other, a) * urn(&[c, d, a]);
                }
            }
        }
        if one_excess < 1.0 && other_excess < 1.0 {
            for a in 0..classes {
                for b in 0..classes {
                    let one_ab = one[a * classes + b];
                    if one_ab < 1e-12 {
                        continue;
                    }
                    for c in 0..classes {
                        for d in 0..classes {
                            let other_cd = other[c * classes + d];
                            if other_cd < 1e-12 {
                                continue;
                            }
                            two_two += one_ab * other_cd * urn(&[a, b, c, d]);
                        }
                    }
                }
            }
        }
        one_excess * other_excess * one_one
            + one_excess * (1.0 - other_excess) * one_two
            + (1.0 - one_excess) * other_excess * two_one
            + (1.0 - one_excess) * (1.0 - other_excess) * two_two
    }

    /// **The closed-form pair likelihood is the four nested sums over both samples' alleles**, to
    /// within rounding, on random rows at thirteen classes, every excess and several spreads.
    #[test]
    fn the_closed_form_pair_likelihood_is_the_four_nested_sums() {
        let classes = 13;
        let mut state = 0x1234_5678_u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut row = || {
            let mut values = vec![0.0; classes * classes];
            for a in 0..classes {
                for b in a..classes {
                    let value = uniform();
                    values[a * classes + b] = value;
                    values[b * classes + a] = value;
                }
            }
            values
        };
        let (one, other) = (row(), row());
        let shares: Vec<f64> = {
            let raw: Vec<f64> = (0..classes).map(|_| uniform() + 0.05).collect();
            let total: f64 = raw.iter().sum();
            raw.iter().map(|value| value / total).collect()
        };
        for spread in [0.01, 0.5, 3.0, 200.0] {
            let weights: Vec<f64> = shares.iter().map(|share| spread * share).collect();
            for (one_excess, other_excess) in [(0.0, 0.0), (0.3, 0.9), (0.98, 0.5), (1.0, 1.0)] {
                let fast = pair_likelihood(
                    &one,
                    one_excess,
                    &other,
                    other_excess,
                    &weights,
                    spread,
                    classes,
                );
                let slow = pair_likelihood_by_nested_sums(
                    &one,
                    one_excess,
                    &other,
                    other_excess,
                    &weights,
                    spread,
                    classes,
                );
                assert!(
                    (fast - slow).abs() <= 1e-10 * slow.abs(),
                    "spread {spread}, excesses {one_excess} {other_excess}: {fast} against {slow}"
                );
            }
        }
    }

    /// **On a drawn class the counting fit returns the stutter and the shares it was drawn at**,
    /// within sampling error, from a start far from them. At a concentration of 0.5 each tract is nearly fixed
    /// at one length, so 2,000 tracts pin a share to about ±0.011.
    #[test]
    fn counting_recovers_a_drawn_class() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let spectrum = spectrum_of(3);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 2000, 20, 10, 1, 13);
        let start = Slippage {
            level: 0.2,
            shorter_share: 0.5,
            fall_off: 0.5,
        };
        let fitted = fit_class_by_counting(&evidence, start, &[0.4; 20], &config, 500, 1e-9);
        assert!((fitted.slippage.level - 0.08).abs() < 0.01, "{fitted:?}");
        assert!(
            (fitted.slippage.shorter_share - 0.8).abs() < 0.05,
            "{fitted:?}"
        );
        for (fitted, truth) in fitted.allele_length_variation.iter().zip(&spectrum) {
            assert!((fitted - truth).abs() < 0.035, "{fitted} against {truth}");
        }
        assert!(
            (0.3..0.8).contains(&fitted.concentration),
            "concentration {} against 0.5",
            fitted.concentration
        );
    }

    /// **A class whose alleles reach 8 repeats short needs a window past the recorded range**
    /// (`str_allele_window.md` §6.1). A fifth of the chromosomes sit 7 or 8 repeats short and
    /// reads are recorded to ±8. With a window of ±6 those alleles have no class of their own,
    /// the −6 class takes them, and their reads can only be explained as slips, so the slip rate
    /// comes out high. At the shipped window of ±9 the slip rate is the one drawn and the far
    /// alleles keep their own shares.
    #[test]
    fn alleles_past_a_narrow_window_are_counted_as_stutter_and_the_shipped_window_keeps_them() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let recorded = RECORDED_OFFSET_RANGE;
        // Over −8 ..= +8, with a small floor so that every class can be drawn.
        let mut spectrum = vec![0.002; (2 * recorded + 1) as usize];
        let at = |offset: i32| (offset + recorded) as usize;
        spectrum[at(-8)] = 0.15;
        spectrum[at(-7)] = 0.05;
        spectrum[at(-1)] = 0.10;
        spectrum[at(0)] = 0.60;
        spectrum[at(1)] = 0.10;
        normalise(&mut spectrum);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 2000, 20, 10, recorded, 19);
        let start = Slippage {
            level: 0.2,
            shorter_share: 0.5,
            fall_off: 0.5,
        };
        let fit_at = |span: i32| {
            let config = SsrFitConfig {
                allele_span: span,
                ..SsrFitConfig::default()
            };
            fit_class_by_counting(&evidence, start, &[0.4; 20], &config, 500, 1e-9)
        };

        let narrow = fit_at(6);
        assert!(narrow.slippage.level > 0.1, "{narrow:?}");
        assert!(narrow.allele_length_variation[0] > 0.15, "{narrow:?}");

        assert_eq!(SsrFitConfig::default().allele_span, recorded + 1);
        let shipped = fit_at(SsrFitConfig::default().allele_span);
        assert!((shipped.slippage.level - 0.08).abs() < 0.01, "{shipped:?}");
        // Which way and how far the far alleles' reads slipped is read from the model, not from
        // their end bucket, which they share with the edge allele's own reads.
        assert!(
            (shipped.slippage.shorter_share - 0.8).abs() < 0.03,
            "{shipped:?}"
        );
        assert!(
            (shipped.slippage.fall_off - 0.25).abs() < 0.04,
            "{shipped:?}"
        );
        // Index 0 is −9 and index 1 is −8. They are held together, not apart: an allele at −8 and
        // one at −9 differ only in reads that slipped inwards, so the two trade mass (about 0.03 of
        // the 0.15 drawn at −8 sits at −9 on this draw).
        let far = shipped.allele_length_variation[0] + shipped.allele_length_variation[1];
        assert!((far - spectrum[at(-8)]).abs() < 0.035, "{shipped:?}");
    }

    /// **In an inbred cohort the spread comes from pairs of samples, not from each sample's two
    /// copies**: at a homozygote excess of 0.95 and three reads a sample, the pairs return the spread a
    /// class was drawn at, where a sample's own two copies carry too little to.
    #[test]
    fn an_inbred_classs_spread_comes_from_pairs_of_samples() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let spectrum = spectrum_of(3);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.95, 2000, 20, 3, 1, 17);
        let start = Slippage {
            level: 0.2,
            shorter_share: 0.5,
            fall_off: 0.5,
        };
        let fitted =
            fit_class_by_counting_with_pairs(&evidence, start, &[0.95; 20], &config, 500, 1e-9, 10);
        assert!(
            (0.35..0.7).contains(&fitted.concentration),
            "spread {} against 0.5: {fitted:?}",
            fitted.concentration
        );
        assert!((fitted.slippage.level - 0.08).abs() < 0.012, "{fitted:?}");
    }
}
