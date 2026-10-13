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
    /// **Never exactly zero**: each length carries a small prior share ([`SHARE_PRIOR_TRACTS`]).
    pub allele_length_variation: Vec<f64>,
    /// How much the class's tracts differ from one another: the Dirichlet concentration, the
    /// between-locus share being `1 / (1 + concentration)`.
    ///
    /// **`None` where the data cannot place it** ([`concentration_is_placed`]): the estimator has
    /// nothing to read — no sample whose two copies can differ, for the own-copies estimator — or the
    /// estimate stopped at an end of its search, `10⁻³` or `10⁴`. A class at one sample and an excess
    /// of one gets `None`, never the search's bound.
    pub concentration: Option<f64>,
    /// The concentration the stutter and shares were counted at, placed or not. Only for refitting
    /// and for the errors, which must evaluate the objective where the fit did.
    pub concentration_counted_at: f64,
    /// Which estimator the concentration came from.
    pub concentration_from: ConcentrationSource,
    /// The shares' prior weight the fit used, in tracts ([`SHARE_PRIOR_TRACTS`] unless a
    /// measurement asked for another); a refit uses the same.
    pub share_prior_tracts: f64,
    /// How many rounds ran, and the largest change of any number in the last.
    pub rounds: usize,
    pub last_change: f64,
    /// Rounds of counting over the whole fit: with the spread read off pairs of samples, every
    /// counting after each reading, where [`rounds`](Self::rounds) is the last one's.
    pub rounds_in_all: usize,
    /// How many times the spread was read off pairs of samples (each reading builds the pairs and
    /// searches the spread); zero for the within-sample estimator. For measuring the fit's cost.
    pub spread_readings: usize,
    /// The objective at the answer.
    pub composite_log_likelihood: f64,
}

impl CountedClass {
    /// The between-locus share `F = 1 / (1 + concentration)`, or `None` with the concentration.
    #[must_use]
    pub fn between_locus_share(&self) -> Option<f64> {
        self.concentration.map(|c| 1.0 / (1.0 + c))
    }
}

/// **How much prior weight each allele length's share gets, in tracts: one tract's worth of allele
/// draws spread evenly over the lengths.**
///
/// Without it a length with no draw gets a share of exactly zero, and once zero it can never come
/// back: its prior is zero, so it draws nothing in the next round. A held-out locus fixed at that
/// length is then impossible, which cost the fit thousands of log-likelihood units a locus on the
/// fitting bench (`research/str_fit/track1_counting_fit_errors_2026-10-12.md` §6). Counting the
/// weight in tracts, not in allele draws, keeps it the same size against the evidence whatever the
/// cohort: at 2,000 samples one tract holds about 4,000 draws, and a prior of a few draws would be
/// lost among them. Each share moves by at most about `1 / tracts` — 0.005 at 200 tracts, a fifth
/// of a standard error at the bench's base setting.
pub const SHARE_PRIOR_TRACTS: f64 = 1.0;

/// The ends of the concentration's search, where an estimate says the data did not place it.
const CONCENTRATION_SEARCH: (f64, f64) = (1e-3, 1e4);

/// **Whether a concentration `estimate` is placed**: its estimator `source` had something to read
/// in `evidence` ([`concentration_is_identified`]), and the estimate stopped inside its search, not
/// within 1% of either end.
#[must_use]
pub fn concentration_is_placed(
    evidence: &StratumEvidence,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
    estimate: f64,
) -> bool {
    let (lowest, highest) = CONCENTRATION_SEARCH;
    concentration_is_identified(evidence, source, homozygote_excess)
        && estimate > lowest * 1.01
        && estimate < highest / 1.01
}

/// **Whether `source` has anything to read the concentration from**: pairs need a tract with two
/// samples; a sample's own two copies need a sample at some tract with a homozygote excess below
/// one, since a sample whose copies are always one draw says nothing about how tracts differ.
#[must_use]
pub fn concentration_is_identified(
    evidence: &StratumEvidence,
    source: ConcentrationSource,
    homozygote_excess: &[f64],
) -> bool {
    match source {
        ConcentrationSource::PairsOfSamples => {
            evidence.tracts.iter().any(|tract| tract.samples.len() >= 2)
        }
        ConcentrationSource::EachSamplesTwoCopies => evidence.tracts.iter().any(|tract| {
            tract
                .samples
                .iter()
                .any(|sample| homozygote_excess[sample.sample as usize] < 1.0)
        }),
    }
}

/// **The expected counts one round gathers over a class, or over one tract of it**
/// ([`expected_counts`], [`expected_counts_by_tract`]): every read shared between "read exactly"
/// and "slipped", every sample's alleles shared among the lengths, given the class's current
/// numbers. Each of the class's numbers is a closed-form function of them.
///
/// **They are also the objective's slope.** The slope of the log of a sum over hidden genotypes
/// and slips is the expected slope of the log-likelihood that knew them (Fisher's identity), and
/// that slope is a function of these counts alone (`counting_fit_errors`).
#[derive(Debug, Clone, Default)]
pub struct ClassCounts {
    /// Per allele class, how many independent allele draws.
    pub draws: Vec<f64>,
    /// Of each sample's two copies not identical by descent, how many pairs, and how many of those
    /// were one draw from the tract's frequencies.
    pub pairs_not_by_descent: f64,
    pub shared_in_tract: f64,
    /// Reads that reported their allele's length, and reads that slipped.
    pub exact: f64,
    pub slipped: f64,
    /// Slipped reads by direction; `steps` is the repeats they slipped summed, over the `stepped`
    /// reads whose direction and length were placed.
    pub shorter: f64,
    pub longer: f64,
    pub steps: f64,
    pub stepped: f64,
    /// The objective summed over what was counted.
    pub log_likelihood: f64,
}

/// The counts one round gathers, by the name the counting fit's own code uses.
pub(super) type Counts = ClassCounts;

impl ClassCounts {
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
        SHARE_PRIOR_TRACTS,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
    )
}

/// [`fit_class_by_counting`] with the shares' prior weight `share_prior_tracts` in place of
/// [`SHARE_PRIOR_TRACTS`] — zero is the fit before the prior, the plain maximum of its objective.
/// **For measuring and for checks**, not for a run.
#[must_use]
pub fn fit_class_by_counting_at_share_prior(
    evidence: &StratumEvidence,
    start: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    share_prior_tracts: f64,
) -> CountedClass {
    count_class(
        evidence,
        start,
        None,
        None,
        1.0,
        share_prior_tracts,
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
    share_prior_tracts: f64,
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
    let distinct = DistinctRows::of(evidence, homozygote_excess);
    while rounds < max_rounds {
        rounds += 1;
        let counts = expected_counts_by_distinct_rows(
            evidence,
            &distinct,
            slippage,
            &shares,
            concentration,
            config,
        );
        log_likelihood = counts.log_likelihood;
        // Each length's prior weight: one tract's worth of draws, spread evenly.
        let draws: f64 = counts.draws.iter().sum();
        let tracts = evidence.tracts_with_reads().max(1) as f64;
        let prior = share_prior_tracts * draws / tracts / classes as f64;
        let new_shares: Vec<f64> = counts
            .draws
            .iter()
            .map(|d| (d + prior) / (draws + prior * classes as f64))
            .collect();
        let new_concentration = if let Some(held) = held_spread {
            held
        } else if counts.shared_in_tract > 0.0 {
            (counts.pairs_not_by_descent / counts.shared_in_tract - 1.0)
                .clamp(CONCENTRATION_SEARCH.0, CONCENTRATION_SEARCH.1)
        } else {
            CONCENTRATION_SEARCH.1
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
    let concentration_from = if held_spread.is_some() {
        ConcentrationSource::PairsOfSamples
    } else {
        ConcentrationSource::EachSamplesTwoCopies
    };
    CountedClass {
        slippage,
        allele_length_variation: shares,
        concentration: concentration_is_placed(
            evidence,
            concentration_from,
            homozygote_excess,
            concentration,
        )
        .then_some(concentration),
        concentration_counted_at: concentration,
        concentration_from,
        share_prior_tracts,
        rounds,
        last_change,
        rounds_in_all: rounds,
        spread_readings: 0,
        composite_log_likelihood: log_likelihood,
    }
}

/// **A class's rows of reads, each distinct one held once**: a row is one sample's reads at one
/// tract, the count in each read-length bucket with its libraries pooled. Two samples whose reads
/// fall in the same buckets in the same numbers hold the same row. On tomato2's 59 samples,
/// 754,601 rows hold 8,780 distinct ones, class by class.
///
/// **Why the round can use it.** A sample's reads enter a round of counting only through its row
/// and its homozygote excess. The costly part, each genotype's likelihood of the reads, depends on
/// the row alone; the excess only mixes two sums of it. So the genotypes are scored once for each
/// distinct row, and each sample adds three numbers to its row's totals ([`RowWeights`]).
pub(super) struct DistinctRows {
    /// Each distinct row once, in the order first met.
    reads: Vec<Vec<u32>>,
    /// Per (tract, sample) with reads, in the order of `evidence.tracts` then each tract's
    /// samples: which distinct row it holds, and the sample's homozygote excess.
    rows: Vec<(u32, f64)>,
    /// Where each tract's rows start in `rows`, and one past the last.
    tract_starts: Vec<usize>,
}

/// One block of tracts' distinct rows of reads, and each of its (tract, sample)'s place among them
/// with the sample's homozygote excess.
type RowsOfABlock = (Vec<Vec<u32>>, Vec<(u32, f64)>);

impl DistinctRows {
    /// **Built in blocks of tracts side by side**, each block finding its own distinct rows, then
    /// the blocks' rows merged in order: the same rows, numbered in the same order, as one pass
    /// over the tracts would give.
    pub(super) fn of(evidence: &StratumEvidence, homozygote_excess: &[f64]) -> Self {
        let buckets = (2 * evidence.read_span + 1) as usize;
        // Per block: its distinct rows in the order first met, and each (tract, sample)'s place
        // among them with the sample's excess.
        let blocks: Vec<RowsOfABlock> = evidence
            .tracts
            .par_chunks(64)
            .map(|tracts| {
                let mut index_of: std::collections::HashMap<Vec<u32>, u32> =
                    std::collections::HashMap::new();
                let mut reads: Vec<Vec<u32>> = Vec::new();
                let mut rows = Vec::new();
                for sample in tracts.iter().flat_map(|tract| &tract.samples) {
                    let pooled = pooled_reads(sample, buckets);
                    let next = reads.len() as u32;
                    let index = *index_of.entry(pooled).or_insert_with_key(|pooled| {
                        reads.push(pooled.clone());
                        next
                    });
                    rows.push((index, homozygote_excess[sample.sample as usize]));
                }
                (reads, rows)
            })
            .collect();
        let mut index_of: std::collections::HashMap<Vec<u32>, u32> =
            std::collections::HashMap::new();
        let mut reads: Vec<Vec<u32>> = Vec::new();
        let renumbered: Vec<Vec<u32>> = blocks
            .iter()
            .map(|(block_reads, _)| {
                block_reads
                    .iter()
                    .map(|pooled| {
                        let next = reads.len() as u32;
                        *index_of.entry(pooled.clone()).or_insert_with_key(|pooled| {
                            reads.push(pooled.clone());
                            next
                        })
                    })
                    .collect()
            })
            .collect();
        let rows: Vec<(u32, f64)> = blocks
            .par_iter()
            .zip(&renumbered)
            .flat_map_iter(|((_, rows), global)| {
                rows.iter()
                    .map(move |(local, excess)| (global[*local as usize], *excess))
            })
            .collect();
        let mut tract_starts = Vec::with_capacity(evidence.tracts.len() + 1);
        let mut start = 0;
        for tract in &evidence.tracts {
            tract_starts.push(start);
            start += tract.samples.len();
        }
        tract_starts.push(start);
        Self {
            reads,
            rows,
            tract_starts,
        }
    }
}

/// One distinct row of reads scored at a round's numbers: the largest genotype log-likelihood,
/// and each genotype's prior times its likelihood (scaled by the largest), for the sample's copies
/// as one draw and as two, with their sums.
struct ScoredRow {
    largest: f64,
    one: Vec<f64>,
    two: Vec<f64>,
    one_total: f64,
    two_total: f64,
}

/// **What a set of samples holding one distinct row add up to**: with `c₁` and `c₂` a sample's
/// chances that its copies are one draw and two, and `T = c₁ Σ one + c₂ Σ two` the sample's
/// likelihood, the sums over the samples of `c₁ / T` and `c₂ / T` turn the row's genotype columns
/// into the samples' summed posterior, and the sum of `(from the tract) / T` gives the one-draw
/// pairs that are two copies drawn alike from the tract.
#[derive(Clone, Copy, Default)]
struct RowWeights {
    as_one: f64,
    as_two: f64,
    alike: f64,
}

/// **One round's fixed inputs**: the genotypes, the read distributions at the round's stutter and
/// each genotype's prior at its shares — what turns a distinct row's weights into counts.
struct RoundModel {
    classes: usize,
    span: i32,
    read_span: i32,
    genotypes: Vec<(usize, usize)>,
    per_allele: Vec<Vec<f64>>,
    stay: f64,
    slips_by_bucket: Vec<Vec<SlipsIntoBucket>>,
    one_draw: Vec<f64>,
    two_draws: Vec<f64>,
    concentration: f64,
}

impl RoundModel {
    fn of(
        evidence: &StratumEvidence,
        slippage: Slippage,
        shares: &[f64],
        concentration: f64,
        config: &SsrFitConfig,
    ) -> Self {
        let classes = allele_classes_of(config);
        let genotypes = genotype_pairs(classes);
        let span = config.allele_span;
        let read_span = evidence.read_span;
        let libraries = one_stutter_for_every_library(evidence);
        let mut distributions =
            ReadDistributionsByLevel::of(&[slippage], &libraries, span, read_span);
        let per_allele = distributions.read_distribution_by_allele.swap_remove(0);
        let slips_by_bucket = (0..classes)
            .map(|class| slips_into_each_bucket(slippage, class as i32 - span, read_span))
            .collect();
        // Each genotype's prior, given a sample's homozygote excess, split into the part where
        // its two copies are one draw (identical by descent) and the part where they are two.
        let one_draw = genotypes
            .iter()
            .map(|(a, b)| if a == b { shares[*a] } else { 0.0 })
            .collect();
        let two_draws = genotypes
            .iter()
            .map(|(a, b)| {
                if a == b {
                    shares[*a] * shares[*a]
                } else {
                    2.0 * shares[*a] * shares[*b]
                }
            })
            .collect();
        Self {
            classes,
            span,
            read_span,
            genotypes,
            per_allele,
            stay: 1.0 - slippage.level,
            slips_by_bucket,
            one_draw,
            two_draws,
            concentration,
        }
    }

    fn score(&self, reads: &[u32]) -> ScoredRow {
        let log_likelihoods: Vec<f64> = self
            .genotypes
            .iter()
            .map(|(a, b)| {
                let mut total = 0.0;
                add_ln_reads(&mut total, &self.per_allele, reads, *a, *b);
                total
            })
            .collect();
        let largest = log_likelihoods
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let scaled: Vec<f64> = log_likelihoods
            .iter()
            .map(|ll| float::exp(ll - largest))
            .collect();
        let one: Vec<f64> = scaled
            .iter()
            .zip(&self.one_draw)
            .map(|(l, p)| l * p)
            .collect();
        let two: Vec<f64> = scaled
            .iter()
            .zip(&self.two_draws)
            .map(|(l, p)| l * p)
            .collect();
        ScoredRow {
            largest,
            one_total: one.iter().sum(),
            two_total: two.iter().sum(),
            one,
            two,
        }
    }

    /// Adds one sample holding `row` at homozygote excess `excess` to `into`, and returns its
    /// log-likelihood; `None` for a sample whose reads no genotype explains.
    fn add_sample(&self, row: &ScoredRow, excess: f64, into: &mut RowWeights) -> Option<f64> {
        // A sample's two copies are one draw if identical by descent, or, failing that, if they
        // drew the same chromosome lineage from the tract's own frequencies: chance
        // 1 / (1 + concentration) under a Dirichlet.
        let from_the_tract = (1.0 - excess) / (1.0 + self.concentration);
        let one_draw_chance = excess + from_the_tract;
        let two_draw_chance = (1.0 - excess) * self.concentration / (1.0 + self.concentration);
        let total = one_draw_chance * row.one_total + two_draw_chance * row.two_total;
        if total <= 0.0 {
            return None;
        }
        into.as_one += one_draw_chance / total;
        into.as_two += two_draw_chance / total;
        into.alike += from_the_tract / total;
        Some(float::ln(total) + row.largest)
    }

    /// **The counts of every sample that holds `reads`**, from the row's genotype columns and the
    /// samples' summed [`RowWeights`]: each read shared between "reported its allele's length" and
    /// "slipped this far, this way", and each sample's alleles among the lengths.
    fn counts_of(&self, reads: &[u32], row: &ScoredRow, weights: RowWeights, into: &mut Counts) {
        let per_allele = &self.per_allele;
        let slot = |class: usize| {
            ((class as i32 - self.span).clamp(-self.read_span, self.read_span) + self.read_span)
                as usize
        };
        for (index, (a, b)) in self.genotypes.iter().enumerate() {
            let one = row.one[index] * weights.as_one;
            let two = row.two[index] * weights.as_two;
            let weight = one + two;
            if weight < 1e-12 {
                continue;
            }
            into.draws[*a] += one + two;
            into.draws[*b] += two;
            let shared = row.one[index] * weights.alike;
            into.shared_in_tract += shared;
            into.pairs_not_by_descent += shared + two;
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
                        self.stay / per_allele[allele][bucket]
                    } else {
                        0.0
                    };
                    let exact = share * unslipped.min(1.0);
                    let slipped = share - exact;
                    into.exact += exact;
                    into.slipped += slipped;
                    // Which way the slipped reads went and how far, split by the model's own
                    // terms: an end bucket can hold slips of both directions from an allele at
                    // or past the recorded range.
                    let landed_in = self.slips_by_bucket[allele][bucket];
                    let landed = landed_in.shorter + landed_in.longer;
                    if slipped > 0.0 && landed > 0.0 {
                        into.shorter += slipped * landed_in.shorter / landed;
                        into.longer += slipped * landed_in.longer / landed;
                        into.steps += slipped * landed_in.steps / landed;
                        into.stepped += slipped;
                    }
                }
            }
        }
    }
}

/// **One round's expected counts over the whole class**, at the stutter `slippage`, the
/// allele-length shares `shares` and the spread `concentration`: the counts of
/// [`expected_counts_by_tract`] summed, with every distinct row of reads scored once however many
/// samples hold it ([`DistinctRows`]).
///
/// **What is not exactly the same as summing the tracts.** The sums are taken in another order,
/// and a genotype is left out of a row's counts when its weight summed over the samples holding the
/// row is below 10⁻¹², where the tract-by-tract count leaves it out of one sample's counts when its
/// posterior is. Both move a count by about the rounding of a double.
pub(super) fn expected_counts_by_distinct_rows(
    evidence: &StratumEvidence,
    distinct: &DistinctRows,
    slippage: Slippage,
    shares: &[f64],
    concentration: f64,
    config: &SsrFitConfig,
) -> Counts {
    let model = RoundModel::of(evidence, slippage, shares, concentration, config);
    let scored: Vec<ScoredRow> = distinct
        .reads
        .par_iter()
        .map(|reads| model.score(reads))
        .collect();
    // In fixed chunks summed in order, so that the answer does not depend on the thread pool.
    let width = distinct.reads.len();
    let chunk = (distinct.rows.len() / 256).max(4096);
    let partials: Vec<(Vec<RowWeights>, f64)> = distinct
        .rows
        .par_chunks(chunk)
        .map(|rows| {
            let mut weights = vec![RowWeights::default(); width];
            let mut log_likelihood = 0.0;
            for (index, excess) in rows {
                let index = *index as usize;
                if let Some(ll) = model.add_sample(&scored[index], *excess, &mut weights[index]) {
                    log_likelihood += ll;
                }
            }
            (weights, log_likelihood)
        })
        .collect();
    let mut weights = vec![RowWeights::default(); width];
    let mut log_likelihood = 0.0;
    for (part, part_ll) in &partials {
        for (into, from) in weights.iter_mut().zip(part) {
            into.as_one += from.as_one;
            into.as_two += from.as_two;
            into.alike += from.alike;
        }
        log_likelihood += part_ll;
    }
    let per_row: Vec<Counts> = distinct
        .reads
        .par_iter()
        .zip(&scored)
        .zip(&weights)
        .map(|((reads, row), weights)| {
            let mut counts = Counts::new(model.classes);
            model.counts_of(reads, row, *weights, &mut counts);
            counts
        })
        .collect();
    let mut counts = Counts::new(model.classes);
    for row in &per_row {
        counts.add(row);
    }
    counts.log_likelihood = log_likelihood;
    counts
}

/// **One round's expected counts, tract by tract**, in the order of `evidence.tracts`, at the
/// stutter `slippage`, the allele-length shares `shares` and the spread `concentration`.
///
/// Weighs every genotype each sample could carry by its prior times how well it explains the
/// sample's reads, then shares each read between "reported its allele's length" and "slipped this
/// far, this way", and each sample's alleles among the lengths. Every distinct row of reads in the
/// class is scored once ([`DistinctRows`]), and the samples of a tract that hold the same row are
/// counted together.
pub(super) fn expected_counts_by_tract(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    concentration: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<Counts> {
    let distinct = DistinctRows::of(evidence, homozygote_excess);
    let model = RoundModel::of(evidence, slippage, shares, concentration, config);
    let scored: Vec<ScoredRow> = distinct
        .reads
        .par_iter()
        .map(|reads| model.score(reads))
        .collect();
    distinct
        .tract_starts
        .par_windows(2)
        .map(|bounds| {
            let mut counts = Counts::new(model.classes);
            // The tract's distinct rows in the order first met, each with its samples' weights.
            let mut held: Vec<(usize, RowWeights)> = Vec::new();
            for (index, excess) in &distinct.rows[bounds[0]..bounds[1]] {
                let index = *index as usize;
                let at = held
                    .iter()
                    .position(|(row, _)| *row == index)
                    .unwrap_or_else(|| {
                        held.push((index, RowWeights::default()));
                        held.len() - 1
                    });
                if let Some(ll) = model.add_sample(&scored[index], *excess, &mut held[at].1) {
                    counts.log_likelihood += ll;
                }
            }
            for (index, weights) in held {
                model.counts_of(&distinct.reads[index], &scored[index], weights, &mut counts);
            }
            counts
        })
        .collect()
}

/// **[`expected_counts_by_tract`] as it was, sample by sample**: the check that counting each
/// distinct row once changes nothing but rounding.
///
/// One round's expected counts, tract by tract, in the order of `evidence.tracts`, at the
/// stutter `slippage`, the allele-length shares `shares` and the spread `concentration`.
///
/// Weighs every genotype each sample could carry by its prior times how well it explains the
/// sample's reads, then shares each read between "reported its allele's length" and "slipped this
/// far, this way", and each sample's alleles among the lengths.
#[cfg(test)]
fn expected_counts_by_tract_row_by_row(
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

/// **One round's expected counts over a class, at the class's current numbers** — the counting
/// fit's E-step summed over the tracts, so that a fit across several classes can gather each
/// class's counts and choose the numbers itself.
///
/// `slippage`, `shares` (one per allele class, `-allele_span ..= +allele_span`) and
/// `concentration` are the class's numbers this round; `homozygote_excess` is indexed by sample.
/// Every library is scored at `slippage` (multipliers are not fitted here).
#[must_use]
pub fn expected_counts(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    concentration: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> ClassCounts {
    expected_counts_by_distinct_rows(
        evidence,
        &DistinctRows::of(evidence, homozygote_excess),
        slippage,
        shares,
        concentration,
        config,
    )
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
    fit_class_by_counting_with_pairs_at_share_prior(
        evidence,
        start,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
        alternations,
        SHARE_PRIOR_TRACTS,
    )
}

/// [`fit_class_by_counting_with_pairs`] with the shares' prior weight `share_prior_tracts` in
/// place of [`SHARE_PRIOR_TRACTS`] — zero is the fit before the prior. **For measuring what the
/// prior does**, not for a run.
#[allow(
    clippy::too_many_arguments,
    reason = "the counting fit's own settings, passed through"
)]
#[must_use]
pub fn fit_class_by_counting_with_pairs_at_share_prior(
    evidence: &StratumEvidence,
    start: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    alternations: usize,
    share_prior_tracts: f64,
) -> CountedClass {
    let fitted = count_class(
        evidence,
        start,
        None,
        None,
        1.0,
        share_prior_tracts,
        homozygote_excess,
        config,
        max_rounds,
        tolerance,
    );
    // With no tract of two samples there is no pair to read the spread from. The search would
    // return its upper end, and the stutter counted at that spread is wrong — at one sample and
    // three reads a slip rate a third of the truth (`track1_counting_fit_errors_2026-10-12.md`
    // §7) — so the spread stays where each sample's own two copies put it.
    if !concentration_is_identified(
        evidence,
        ConcentrationSource::PairsOfSamples,
        homozygote_excess,
    ) {
        return fitted;
    }
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
    /// sample, where there is no pair; there [`fit_class_by_counting_with_pairs`] uses the other.
    PairsOfSamples,
}

/// **The class refitted from an earlier answer `from`** rather than from flat shares: what a
/// bootstrap replicate needs, whose answer lies near the full fit's. Uses the estimator `from` was
/// made by, with the same rounds, tolerance and alternations as a fit from scratch.
#[must_use]
pub fn refit_class_by_counting_from(
    evidence: &StratumEvidence,
    from: &CountedClass,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
    alternations: usize,
) -> CountedClass {
    match from.concentration_from {
        ConcentrationSource::EachSamplesTwoCopies => count_class(
            evidence,
            from.slippage,
            Some(&from.allele_length_variation),
            None,
            from.concentration_counted_at,
            from.share_prior_tracts,
            homozygote_excess,
            config,
            max_rounds,
            tolerance,
        ),
        ConcentrationSource::PairsOfSamples => alternate_with_pairs(
            evidence,
            from.clone(),
            from.concentration_counted_at,
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
    let mut rounds_in_all = fitted.rounds_in_all;
    let mut spread_readings = fitted.spread_readings;
    for _ in 0..alternations {
        spread_readings += 1;
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
            fitted.share_prior_tracts,
            homozygote_excess,
            config,
            max_rounds,
            tolerance,
        );
        rounds_in_all += fitted.rounds;
        let settled = (float::ln(next) - float::ln(spread)).abs() < 0.01;
        spread = next;
        if settled {
            break;
        }
    }
    fitted.rounds_in_all = rounds_in_all;
    fitted.spread_readings = spread_readings;
    fitted
}

/// **The spread that maximises the pairwise likelihood of the class**, stutter and shares held.
///
/// At a tract, every sample with reads is paired with the next (and the last with the first). A
/// pair's likelihood sums both samples' reads over their alleles, each sample's copies one draw with
/// chance its homozygote excess and two otherwise, the draws from the tract's frequencies weighted
/// as a Pólya urn: a draw takes length `a` with chance `(spread · share_a + earlier draws of a) /
/// (spread + earlier draws)`. Searched on the log scale between 10⁻³ and 10⁴.
///
/// **Each pair is reduced to four numbers before the search** ([`PairCubic`]): with the shares
/// held, a pair's likelihood is a ratio of two cubics in the spread, so the search's 61 evaluations
/// cost four multiply-adds a pair rather than a sweep of both samples' allele-pair tables, and a
/// table is held once for each distinct row of reads ([`pair_cubics_by_tract`]), not for each
/// sample.
#[must_use]
pub fn fit_spread_by_pairs(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> f64 {
    let per_tract = pair_cubics_by_tract(evidence, slippage, shares, homozygote_excess, config);
    let objective = |ln_spread: f64| -> f64 {
        let spread = float::exp(ln_spread);
        per_tract
            .par_iter()
            .map(|cubics| scaled_ln_pairs_of_tract(cubics, spread))
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

/// **One pair of samples' likelihood at a tract as a function of the spread `α` alone**, the
/// shares held: `N(α) / ((α + 1)(α + 2)(α + 3))`, with `N(α) = c₀ + c₁α + c₂α² + c₃α³` and these
/// the four coefficients, lowest power first. Every coefficient is a sum of products of
/// likelihoods, shares and excesses, so none is negative and `N` is evaluated without cancellation.
///
/// **Why a cubic.** Write `w = α · s` for the urn's weights. Every term of the pair's likelihood
/// ([`pair_likelihood`]) is a sum over alleles of a row's likelihoods times a product of weights,
/// over the urn's denominator `α(α + 1)…`: each sample one draw over `α(α + 1)`, one and two over
/// `α(α + 1)(α + 2)`, two and two over `α(α + 1)(α + 2)(α + 3)`. Pulling `α` out of each weight
/// leaves sums that depend on the shares only, a polynomial in `α` over each denominator, and
/// putting the four terms over the largest denominator gives a cubic over `(α + 1)(α + 2)(α + 3)`.
/// The sums, with `L` a sample's row, `D_x = L(x, x)`, `U_x = Σ_c L(c, x) s_c`, `Q = Σ_x U_x s_x`
/// and `R = Σ_x D_x s_x`:
///
/// - `S = Σ_a D¹_a D²_a s_a`, `X₁₂ = Σ_a D¹_a s_a U²_a`, `X₂₁ = Σ_a D²_a s_a U¹_a`,
///   `Y = Σ_a s_a U¹_a U²_a`, `Z = Σ_ab L¹(a, b) L²(a, b) s_a s_b`;
/// - each one draw: `(α R¹R² + S) / (α + 1)`;
/// - the first one draw, the other two: `(α² R¹Q² + α (R¹R² + 2 X₁₂) + 2 S) / ((α + 1)(α + 2))`,
///   and the other way round with 1 and 2 exchanged;
/// - each two draws: `(α³ Q¹Q² + α² (Q¹R² + R¹Q² + 4 Y) + α (R¹R² + 4 X₁₂ + 4 X₂₁ + 2 Z) + 6 S) /
///   ((α + 1)(α + 2)(α + 3))`.
///
/// Checked against [`pair_likelihood`] in `a_pairs_cubic_is_its_closed_form_likelihood`.
pub(super) type PairCubic = [f64; 4];

/// A sample's row of ordered allele-pair likelihoods summarised at the shares: its diagonal `D`,
/// `U_x = Σ_c L(c, x) s_c`, `Q = Σ_x U_x s_x` and `R = Σ_x D_x s_x` ([`PairCubic`]).
struct RowAtShares {
    diagonal: Vec<f64>,
    u: Vec<f64>,
    q: f64,
    r: f64,
}

impl RowAtShares {
    fn of(row: &[f64], shares: &[f64], classes: usize) -> Self {
        let diagonal: Vec<f64> = (0..classes).map(|x| row[x * classes + x]).collect();
        let u: Vec<f64> = (0..classes)
            .map(|x| (0..classes).map(|c| row[c * classes + x] * shares[c]).sum())
            .collect();
        let q = (0..classes).map(|x| u[x] * shares[x]).sum();
        let r = (0..classes).map(|x| diagonal[x] * shares[x]).sum();
        Self { diagonal, u, q, r }
    }
}

/// **The five sums of two rows that a pair's [`PairCubic`] needs beyond each row's own**: `S`,
/// `X₁₂`, `X₂₁`, `Y` and `Z` of [`PairCubic`], at `shares`. They depend on the two rows and the
/// shares only, not on the samples' excesses, so two pairs of samples holding the same two rows of
/// reads share them.
#[derive(Debug, Clone, Copy)]
struct PairSums {
    s: f64,
    x12: f64,
    x21: f64,
    y: f64,
    z: f64,
}

impl PairSums {
    fn of(
        (one, one_at): (&[f64], &RowAtShares),
        (other, other_at): (&[f64], &RowAtShares),
        shares: &[f64],
        classes: usize,
    ) -> Self {
        let (d1, u1) = (&one_at.diagonal, &one_at.u);
        let (d2, u2) = (&other_at.diagonal, &other_at.u);
        let mut s = 0.0;
        let mut x12 = 0.0;
        let mut x21 = 0.0;
        let mut y = 0.0;
        for a in 0..classes {
            s += d1[a] * d2[a] * shares[a];
            x12 += d1[a] * shares[a] * u2[a];
            x21 += d2[a] * shares[a] * u1[a];
            y += shares[a] * u1[a] * u2[a];
        }
        let mut z = 0.0;
        for a in 0..classes {
            let mut inner = 0.0;
            for b in 0..classes {
                inner += one[a * classes + b] * other[a * classes + b] * shares[b];
            }
            z += inner * shares[a];
        }
        Self { s, x12, x21, y, z }
    }
}

/// **One pair's [`PairCubic`]** from its two rows' [`PairSums`], each row's own sums and the two
/// samples' excesses.
fn pair_cubic_of_sums(
    sums: PairSums,
    (one_at, one_excess): (&RowAtShares, f64),
    (other_at, other_excess): (&RowAtShares, f64),
) -> PairCubic {
    let PairSums { s, x12, x21, y, z } = sums;
    let (q1, r1) = (one_at.q, one_at.r);
    let (q2, r2) = (other_at.q, other_at.r);
    let rr = r1 * r2;
    // Each one draw, over (α + 1): times (α + 2)(α + 3) = α² + 5α + 6.
    let both_one = one_excess * other_excess;
    let (e1, e0) = (both_one * rr, both_one * s);
    // One draw and two, over (α + 1)(α + 2): times (α + 3).
    let first_one = one_excess * (1.0 - other_excess);
    let second_one = (1.0 - one_excess) * other_excess;
    let m2 = first_one * r1 * q2 + second_one * r2 * q1;
    let m1 = first_one * (rr + 2.0 * x12) + second_one * (rr + 2.0 * x21);
    let m0 = (first_one + second_one) * 2.0 * s;
    // Each two draws, over (α + 1)(α + 2)(α + 3).
    let both_two = (1.0 - one_excess) * (1.0 - other_excess);
    let t3 = both_two * q1 * q2;
    let t2 = both_two * (q1 * r2 + r1 * q2 + 4.0 * y);
    let t1 = both_two * (rr + 4.0 * x12 + 4.0 * x21 + 2.0 * z);
    let t0 = both_two * 6.0 * s;
    [
        6.0 * e0 + 3.0 * m0 + t0,
        6.0 * e1 + 5.0 * e0 + 3.0 * m1 + m0 + t1,
        5.0 * e1 + e0 + 3.0 * m2 + m1 + t2,
        e1 + m2 + t3,
    ]
}

/// **One pair's [`PairCubic`]** from both samples' rows and excesses at `shares`.
#[cfg(test)]
fn pair_cubic(
    (one, one_at, one_excess): (&[f64], &RowAtShares, f64),
    (other, other_at, other_excess): (&[f64], &RowAtShares, f64),
    shares: &[f64],
    classes: usize,
) -> PairCubic {
    pair_cubic_of_sums(
        PairSums::of((one, one_at), (other, other_at), shares, classes),
        (one_at, one_excess),
        (other_at, other_excess),
    )
}

/// **One tract's pairwise log-likelihood at the spread `spread`, less its rows' offsets**, from its
/// pairs' [`PairCubic`]s. The rows' offsets are left out: they depend on neither the shares nor the
/// spread.
pub(super) fn scaled_ln_pairs_of_tract(cubics: &[PairCubic], spread: f64) -> f64 {
    if cubics.is_empty() {
        return 0.0;
    }
    let below = (spread + 1.0) * (spread + 2.0) * (spread + 3.0);
    cubics
        .iter()
        .map(|[c0, c1, c2, c3]| {
            let above = ((c3 * spread + c2) * spread + c1) * spread + c0;
            float::ln((above / below).max(1e-300))
        })
        .sum()
}

/// **A row of reads' likelihood under every ordered allele pair**, `first · classes + second`,
/// divided by the largest of them.
fn ordered_row(
    reads: &[u32],
    per_allele: &[Vec<f64>],
    genotypes: &[(usize, usize)],
    classes: usize,
) -> Vec<f64> {
    let log_likelihoods: Vec<f64> = genotypes
        .iter()
        .map(|(first, second)| {
            let mut total = 0.0;
            add_ln_reads(&mut total, per_allele, reads, *first, *second);
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
    ordered
}

/// **Every tract's pairs as [`PairCubic`]s at the stutter `slippage` and the shares `shares`**, in
/// the order of `evidence.tracts`; a tract with fewer than two samples has none, since it has no
/// pair. At a tract every sample with reads is paired with the next, and the last with the first
/// (one pair when there are two).
///
/// **Each distinct row of reads is scored once** ([`DistinctRows`]), and each distinct pair of
/// rows has its [`PairSums`] taken once, so the work grows with the distinct rows and pairs of rows
/// rather than with the samples; 32 bytes a pair are kept, and `8 · classes²` (2.9 kB at ±9) a
/// distinct row.
pub(super) fn pair_cubics_by_tract(
    evidence: &StratumEvidence,
    slippage: Slippage,
    shares: &[f64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<Vec<PairCubic>> {
    let classes = allele_classes_of(config);
    let genotypes = genotype_pairs(classes);
    let libraries = one_stutter_for_every_library(evidence);
    let distributions = ReadDistributionsByLevel::of(
        &[slippage],
        &libraries,
        config.allele_span,
        evidence.read_span,
    );
    let per_allele = &distributions.read_distribution_by_allele[0];
    let distinct = DistinctRows::of(evidence, homozygote_excess);
    let ordered: Vec<Vec<f64>> = distinct
        .reads
        .par_iter()
        .map(|reads| ordered_row(reads, per_allele, &genotypes, classes))
        .collect();
    let at: Vec<RowAtShares> = ordered
        .par_iter()
        .map(|row| RowAtShares::of(row, shares, classes))
        .collect();
    // Each tract's pairs, as positions in `distinct.rows`.
    let pairs_of = |bounds: &[usize]| -> Vec<(usize, usize)> {
        let (start, count) = (bounds[0], bounds[1] - bounds[0]);
        match count {
            0 | 1 => Vec::new(),
            2 => vec![(start, start + 1)],
            _ => (0..count)
                .map(|index| (start + index, start + (index + 1) % count))
                .collect(),
        }
    };
    let row_of = |position: usize| distinct.rows[position].0;
    let mut row_pairs: Vec<(u32, u32)> = distinct
        .tract_starts
        .windows(2)
        .flat_map(&pairs_of)
        .map(|(one, other)| (row_of(one), row_of(other)))
        .collect();
    row_pairs.par_sort_unstable();
    row_pairs.dedup();
    let sums: Vec<PairSums> = row_pairs
        .par_iter()
        .map(|(one, other)| {
            let (one, other) = (*one as usize, *other as usize);
            PairSums::of(
                (&ordered[one], &at[one]),
                (&ordered[other], &at[other]),
                shares,
                classes,
            )
        })
        .collect();
    distinct
        .tract_starts
        .par_windows(2)
        .map(|bounds| {
            pairs_of(bounds)
                .into_iter()
                .map(|(one, other)| {
                    let (one_row, one_excess) = distinct.rows[one];
                    let (other_row, other_excess) = distinct.rows[other];
                    // PANIC-FREE: every pair of rows a tract makes was put in the list above.
                    let index = row_pairs
                        .binary_search(&(one_row, other_row))
                        .expect("every pair of rows is listed");
                    pair_cubic_of_sums(
                        sums[index],
                        (&at[one_row as usize], one_excess),
                        (&at[other_row as usize], other_excess),
                    )
                })
                .collect()
        })
        .collect()
}

/// **Every pair of samples a class's spread is read from**, at the stutter and shares it was built
/// at: at a tract, every sample with reads paired with the next (and the last with the first), each
/// pair held as its [`PairCubic`]. Built once, then [`log_likelihood`](Self::log_likelihood) scores
/// any spread.
pub struct SamplePairs {
    /// Per tract with two samples or more, its pairs.
    per_tract: Vec<Vec<PairCubic>>,
}

impl SamplePairs {
    /// The pairs of `evidence`'s tracts, each sample's reads scored at `slippage`, the pairs reduced
    /// at the allele-length shares `shares`.
    #[must_use]
    pub fn of(
        evidence: &StratumEvidence,
        slippage: Slippage,
        shares: &[f64],
        homozygote_excess: &[f64],
        config: &SsrFitConfig,
    ) -> Self {
        let per_tract = pair_cubics_by_tract(evidence, slippage, shares, homozygote_excess, config)
            .into_iter()
            .filter(|pairs| !pairs.is_empty())
            .collect();
        Self { per_tract }
    }

    /// **The pairwise log-likelihood of each tract with two samples or more**, in the order
    /// [`log_likelihood`](Self::log_likelihood) sums them: what a standard error that treats the
    /// tract as the independent unit needs.
    #[must_use]
    pub fn log_likelihood_by_tract(&self, ln_spread: f64) -> Vec<f64> {
        let spread = float::exp(ln_spread);
        self.per_tract
            .par_iter()
            .map(|cubics| scaled_ln_pairs_of_tract(cubics, spread))
            .collect()
    }

    /// Whether any tract has two samples with reads.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.per_tract.is_empty()
    }

    /// **The pairwise log-likelihood at spread `exp(ln_spread)`**, at the shares the pairs were built
    /// at, up to a constant that depends on neither.
    #[must_use]
    pub fn log_likelihood(&self, ln_spread: f64) -> f64 {
        let spread = float::exp(ln_spread);
        self.per_tract
            .par_iter()
            .map(|cubics| scaled_ln_pairs_of_tract(cubics, spread))
            .sum()
    }
}

/// One pair of samples' likelihood at a tract: both samples' reads summed over their alleles, each
/// sample one draw (chance its excess) or two, the draws weighted as a Pólya urn.
///
/// **Summed in closed form, a pair in `classes²` steps rather than `classes⁴`.** With `w = spread ·
/// share`, the urn's weight for a later draw is `w_x` plus how many earlier draws were `x`, so every
/// sum over the later sample's alleles reduces to a few totals of its row: `Q = Σ L(c,d) w_c w_d`,
/// `R = Σ L(c,c) w_c`, `u_x = Σ_c L(c,x) w_c` and the diagonal `D_x = L(x,x)`. Checked against the
/// four nested sums in `the_closed_form_pair_likelihood_is_the_four_nested_sums`.
#[cfg(test)]
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

    /// **Counting each distinct row of reads once gives the counts sample by sample did**, tract
    /// by tract and over the class, to within rounding: on a drawn class whose samples each carry
    /// their own homozygote excess, at a point away from the truth.
    #[test]
    fn distinct_rows_count_what_sample_by_sample_counted() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 2,
            ..SsrFitConfig::default()
        };
        let spectrum = spectrum_of(5);
        let samples = 20;
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.9, 300, samples, 3, 2, 23);
        let excess: Vec<f64> = (0..samples).map(|s| 0.98 * s as f64 / 19.0).collect();
        let at = Slippage {
            level: 0.12,
            shorter_share: 0.6,
            fall_off: 0.4,
        };
        let shares = [0.1, 0.15, 0.4, 0.2, 0.15];
        let close = |left: f64, right: f64| (left - right).abs() <= 1e-11 * right.abs().max(1e-300);
        let same = |left: &Counts, right: &Counts| {
            left.draws
                .iter()
                .zip(&right.draws)
                .all(|(l, r)| close(*l, *r))
                && close(left.pairs_not_by_descent, right.pairs_not_by_descent)
                && close(left.shared_in_tract, right.shared_in_tract)
                && close(left.exact, right.exact)
                && close(left.slipped, right.slipped)
                && close(left.shorter, right.shorter)
                && close(left.longer, right.longer)
                && close(left.steps, right.steps)
                && close(left.stepped, right.stepped)
                && close(left.log_likelihood, right.log_likelihood)
        };
        let by_row =
            expected_counts_by_tract_row_by_row(&evidence, at, &shares, 0.7, &excess, &config);
        let by_tract = expected_counts_by_tract(&evidence, at, &shares, 0.7, &excess, &config);
        assert_eq!(by_row.len(), by_tract.len());
        for (index, (row, tract)) in by_row.iter().zip(&by_tract).enumerate() {
            assert!(same(tract, row), "tract {index}: {tract:?} against {row:?}");
        }
        let mut summed = Counts::new(5);
        for tract in &by_row {
            summed.add(tract);
        }
        let distinct = DistinctRows::of(&evidence, &excess);
        assert!(distinct.reads.len() < evidence.rows_with_reads());
        let class =
            expected_counts_by_distinct_rows(&evidence, &distinct, at, &shares, 0.7, &config);
        assert!(same(&class, &summed), "{class:?} against {summed:?}");
    }

    /// **A pair's cubic in the spread is its closed-form likelihood**, to within rounding, on random
    /// rows at nineteen classes (the shipped ±9), every excess, spreads across the search's range,
    /// and a row with exact zeros, as a sample whose reads rule alleles out has.
    #[test]
    fn a_pairs_cubic_is_its_closed_form_likelihood() {
        let classes = 19;
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut row = |zeros: bool| {
            let mut values = vec![0.0; classes * classes];
            for a in 0..classes {
                for b in a..classes {
                    let value = if zeros && (a + b) % 3 == 0 {
                        0.0
                    } else {
                        uniform()
                    };
                    values[a * classes + b] = value;
                    values[b * classes + a] = value;
                }
            }
            values
        };
        let (one, other) = (row(false), row(true));
        let shares: Vec<f64> = {
            let raw: Vec<f64> = (0..classes).map(|_| uniform().powi(3) + 1e-4).collect();
            let total: f64 = raw.iter().sum();
            raw.iter().map(|value| value / total).collect()
        };
        let (one_at, other_at) = (
            RowAtShares::of(&one, &shares, classes),
            RowAtShares::of(&other, &shares, classes),
        );
        for spread in [1e-3, 0.01, 0.5, 3.0, 200.0, 1e4] {
            let weights: Vec<f64> = shares.iter().map(|share| spread * share).collect();
            for (one_excess, other_excess) in [(0.0, 0.0), (0.3, 0.9), (0.98, 0.5), (1.0, 1.0)] {
                let closed = pair_likelihood(
                    &one,
                    one_excess,
                    &other,
                    other_excess,
                    &weights,
                    spread,
                    classes,
                );
                let [c0, c1, c2, c3] = pair_cubic(
                    (&one, &one_at, one_excess),
                    (&other, &other_at, other_excess),
                    &shares,
                    classes,
                );
                let cubic = (((c3 * spread + c2) * spread + c1) * spread + c0)
                    / ((spread + 1.0) * (spread + 2.0) * (spread + 3.0));
                assert!(
                    (cubic - closed).abs() <= 1e-12 * closed.abs(),
                    "spread {spread}, excesses {one_excess} {other_excess}: {cubic} against {closed}"
                );
            }
        }
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
            fitted
                .concentration
                .is_some_and(|c| (0.3..0.8).contains(&c)),
            "concentration {} against 0.5",
            fitted.concentration_counted_at
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
            fitted
                .concentration
                .is_some_and(|c| (0.35..0.7).contains(&c)),
            "spread {} against 0.5: {fitted:?}",
            fitted.concentration_counted_at
        );
        assert!((fitted.slippage.level - 0.08).abs() < 0.012, "{fitted:?}");
    }

    /// **A length no tract carries keeps a small share, never exactly zero**: drawn at three
    /// lengths and fitted at five, the two outer lengths get about the prior's weight — one tract's
    /// worth of draws spread over five lengths, so about `1 / (5 · 300)` — and the drawn shares are
    /// still recovered.
    #[test]
    fn a_length_no_tract_carries_keeps_a_small_share() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let spectrum = [0.0, 0.25, 0.5, 0.25, 0.0];
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 300, 20, 10, 2, 23);
        let config = SsrFitConfig {
            allele_span: 2,
            ..SsrFitConfig::default()
        };
        let fitted =
            fit_class_by_counting_with_pairs(&evidence, truth, &[0.4; 20], &config, 2000, 1e-7, 10);
        for outer in [0, 4] {
            let share = fitted.allele_length_variation[outer];
            assert!(
                (1e-4..2e-3).contains(&share),
                "length {outer}: share {share}: {fitted:?}"
            );
        }
        for (fitted, truth) in fitted.allele_length_variation.iter().zip(&spectrum) {
            assert!((fitted - truth).abs() < 0.06, "{fitted} against {truth}");
        }
    }

    /// **At one sample the spread comes from the sample's own two copies, and is absent when they
    /// cannot differ.** The pairs fit has no pair at one sample: it uses the own-copies estimator,
    /// which places the spread when the homozygote excess is below one and returns `None` at one.
    #[test]
    fn at_one_sample_the_spread_is_read_from_its_own_copies_or_absent() {
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
        let outbred = draw_stratum(truth, &spectrum, 0.5, 0.4, 2000, 1, 30, 1, 29);
        let fitted =
            fit_class_by_counting_with_pairs(&outbred, truth, &[0.4], &config, 2000, 1e-7, 10);
        assert_eq!(
            fitted.concentration_from,
            ConcentrationSource::EachSamplesTwoCopies
        );
        assert!(
            fitted
                .concentration
                .is_some_and(|c| (0.3..0.8).contains(&c)),
            "{fitted:?}"
        );
        assert!((fitted.slippage.level - 0.08).abs() < 0.01, "{fitted:?}");

        let inbred = draw_stratum(truth, &spectrum, 0.5, 1.0, 2000, 1, 30, 1, 31);
        let fitted =
            fit_class_by_counting_with_pairs(&inbred, truth, &[1.0], &config, 2000, 1e-7, 10);
        assert_eq!(fitted.concentration, None, "{fitted:?}");
        assert_eq!(fitted.between_locus_share(), None);
        assert!((fitted.slippage.level - 0.08).abs() < 0.01, "{fitted:?}");
    }
}
