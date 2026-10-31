//! **A repeat class's stutter and allele-length shares by expectation maximisation, without the
//! average over tract frequencies** — the owner's routes 1 and 2 together (2026-10-11), a prototype.
//!
//! The climb, the fit this replaced (deleted 2026-10-16), averaged each tract's likelihood over how
//! its alleles could be shared among its chromosomes, and that average, over 256 points, is several
//! log-likelihood units a tract wrong at thirteen allele lengths
//! (`fit_precision_quadrature_average_2026-10-01.md`). It also made the fit slow. This estimates the same stutter and the class's allele-length shares from a simpler
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
//!
//! **Why the file is named for the expected counts.** The penalised fit ([`penalised_fit`])
//! becomes the production fit, and what it takes from here is that first step — `expected_counts` and [`ClassCounts`] — and the
//! likelihood of pairs of samples ([`SamplePairs`]). The counting fit's own drivers
//! (`fit_class_by_counting*`) stay here as the zero-penalty answer the penalised fit's tests
//! compare against.

use ahash::AHashMap;

use super::*;

/// **The ends of the concentration's search**, where an estimate says the data did not place it.
/// The penalised fit clamps its concentration to these same ends, so an estimate pinned at its
/// clamp is the one [`concentration_stopped_inside_its_search`] reports as absent.
pub(super) const CONCENTRATION_SEARCH: (f64, f64) = (1e-3, 1e4);

/// **Whether a concentration `estimate` stopped inside its search**, not within 1% of either end
/// of `10⁻³ ..= 10⁴`. An estimate at an end is where the search ran out of room, not where the
/// data put it, so it is reported as absent (spec `str_penalised_fit.md` §2.5). Not a number is
/// not inside.
#[must_use]
pub fn concentration_stopped_inside_its_search(estimate: f64) -> bool {
    let (lowest, highest) = CONCENTRATION_SEARCH;
    estimate > lowest * 1.01 && estimate < highest / 1.01
}

/// **The expected counts one round gathers over a class, or over one tract of it**
/// (`expected_counts`, `expected_counts_by_tract`): every read shared
/// between "read exactly" and "slipped", every sample's alleles shared among the lengths, given
/// the class's current numbers. Each of the class's numbers is a closed-form function of them.
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

/// About how many (tract, sample) rows one block of tracts holds while [`DistinctRows::of`] builds
/// a class's table.
const ROWS_PER_BLOCK: usize = 8192;

/// One block of tracts' distinct rows of reads, and each of its (tract, sample)'s place among them
/// with the sample's homozygote excess.
type RowsOfABlock = (Vec<Vec<u32>>, Vec<(u32, f64)>);

impl DistinctRows {
    /// **Built in blocks of tracts side by side**, each block finding its own distinct rows, then
    /// the blocks' rows merged in order: the same rows, numbered in the same order, as one pass
    /// over the tracts would give, whatever the size of a block.
    ///
    /// A block holds about [`ROWS_PER_BLOCK`] rows, and from 1 to 64 tracts: at 2,000 samples a
    /// block of 64 tracts left a class of 100 loci two blocks, so two threads built its table
    /// while the rest waited. That mattered once the penalised fit visits classes one at a time.
    pub(super) fn of(evidence: &StratumEvidence, homozygote_excess: &[f64]) -> Self {
        let rows_in_class: usize = evidence.tracts.iter().map(|t| t.samples.len()).sum();
        let tracts_per_block = (ROWS_PER_BLOCK * evidence.tracts.len())
            .checked_div(rows_in_class)
            .unwrap_or(64)
            .clamp(1, 64);
        Self::in_blocks_of(evidence, homozygote_excess, tracts_per_block)
    }

    /// **How many distinct rows the class holds**: the rows whose genotypes a round scores, one
    /// for every set of (tract, sample) pairs with the same reads.
    pub(super) fn distinct_rows(&self) -> usize {
        self.reads.len()
    }

    /// [`of`](Self::of) with `tracts_per_block` tracts to a block, at least one.
    fn in_blocks_of(
        evidence: &StratumEvidence,
        homozygote_excess: &[f64],
        tracts_per_block: usize,
    ) -> Self {
        let buckets = (2 * evidence.read_span + 1) as usize;
        // Per block: its distinct rows in the order first met, and each (tract, sample)'s place
        // among them with the sample's excess. Every sample's reads are pooled into one scratch
        // row and looked up by value; only a row not met before is copied.
        let blocks: Vec<RowsOfABlock> = evidence
            .tracts
            .par_chunks(tracts_per_block)
            .map(|tracts| {
                let mut index_of: AHashMap<Vec<u32>, u32> = AHashMap::new();
                let mut rows = Vec::new();
                let mut pooled = vec![0_u32; buckets];
                for sample in tracts.iter().flat_map(|tract| &tract.samples) {
                    pool_reads_into(sample, &mut pooled);
                    let index = match index_of.get(pooled.as_slice()) {
                        Some(index) => *index,
                        None => {
                            let next = index_of.len() as u32;
                            index_of.insert(pooled.clone(), next);
                            next
                        }
                    };
                    rows.push((index, homozygote_excess[sample.sample as usize]));
                }
                // Each row put at its own number, so the map's order is never read.
                let mut reads = vec![Vec::new(); index_of.len()];
                for (row, index) in index_of {
                    reads[index as usize] = row;
                }
                (reads, rows)
            })
            .collect();
        // Keyed by the blocks' own rows, so a row met in an earlier block is not copied again.
        let mut index_of: AHashMap<&[u32], u32> = AHashMap::new();
        let mut reads: Vec<Vec<u32>> = Vec::new();
        let renumbered: Vec<Vec<u32>> = blocks
            .iter()
            .map(|(block_reads, _)| {
                block_reads
                    .iter()
                    .map(|pooled| {
                        *index_of.entry(pooled.as_slice()).or_insert_with(|| {
                            reads.push(pooled.clone());
                            reads.len() as u32 - 1
                        })
                    })
                    .collect()
            })
            .collect();
        // Freed before the rows' table is built.
        drop(index_of);
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

/// **The log chance that one read reports each bucket's length, for each genotype**: one row per
/// genotype in `genotypes` order, one column per bucket, each entry
/// `ln max(½ (per_allele[a][bucket] + per_allele[b][bucket]), 10⁻³⁰⁰)` for the genotype `(a, b)`.
///
/// **Taken once for a round, not once for each row of reads.** The log depends on the genotype and
/// the bucket alone, yet scoring every distinct row took it afresh, and `ln` was a quarter of the
/// research bench's CPU samples at 2,000 samples of 30 reads. [`ln_reads`] reads the table and adds
/// the same terms in the same order, so a row's log-likelihood is the same to the bit.
fn ln_read_chances(per_allele: &[Vec<f64>], genotypes: &[(usize, usize)]) -> Vec<Vec<f64>> {
    genotypes
        .iter()
        .map(|(a, b)| {
            per_allele[*a]
                .iter()
                .zip(&per_allele[*b])
                .map(|(from_a, from_b)| float::ln((0.5 * (from_a + from_b)).max(1e-300)))
                .collect()
        })
        .collect()
}

/// **The log-likelihood of a row of reads under one genotype**, from that genotype's row of
/// [`ln_read_chances`]: each non-empty bucket's read count times its log chance, summed bucket by
/// bucket from zero.
#[inline]
fn ln_reads(ln_chances: &[f64], reads: &[u32]) -> f64 {
    debug_assert_eq!(ln_chances.len(), reads.len(), "one log chance a bucket");
    let mut total = 0.0;
    for (reads, ln_chance) in reads.iter().zip(ln_chances) {
        if *reads == 0 {
            continue;
        }
        total += f64::from(*reads) * ln_chance;
    }
    total
}

/// **One distinct row of reads scored at one stutter**: the largest genotype log-likelihood, and
/// each genotype's likelihood divided by the largest, in `genotypes` order.
///
/// It holds no prior. A genotype's weight in a round is this likelihood times the genotype's prior
/// ([`GenotypePrior`]), multiplied where it is read, so a step that moves only the shares or the
/// spread reuses the row as it is.
struct ScoredRow {
    largest: f64,
    scaled: Vec<f64>,
}

/// **A scored row's likelihoods summed under the prior**: over the genotypes, each scaled
/// likelihood times the genotype's prior with the sample's copies one draw, and times its prior
/// with them two draws. A sample's likelihood mixes the two by its chances of each.
#[derive(Clone, Copy)]
struct PriorWeightedTotals {
    one_draw: f64,
    two_draws: f64,
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

/// **What scores a row of reads at one stutter**: the genotypes, each allele's read-length
/// distribution, the log chance of a read in each bucket under each genotype, and where the
/// slipped reads land. Everything in it depends on the stutter alone, so the rows scored with it
/// serve every set of shares and spread ([`ScoredDistinctRows`]).
struct ReadModel {
    classes: usize,
    span: i32,
    read_span: i32,
    genotypes: Vec<(usize, usize)>,
    per_allele: Vec<Vec<f64>>,
    /// [`ln_read_chances`] of `per_allele`: a row per genotype in `genotypes` order.
    ln_read_chances: Vec<Vec<f64>>,
    stay: f64,
    slips_by_bucket: Vec<Vec<SlipsIntoBucket>>,
}

impl ReadModel {
    fn of(evidence: &StratumEvidence, slippage: Slippage, allele_span: i32) -> Self {
        let classes = allele_classes_in(allele_span);
        let genotypes = genotype_pairs(classes);
        let span = allele_span;
        let read_span = evidence.read_span;
        let libraries = one_stutter_for_every_library(evidence);
        let mut distributions =
            ReadDistributionsByLevel::of(&[slippage], &libraries, span, read_span);
        let per_allele = distributions.read_distribution_by_allele.swap_remove(0);
        let ln_read_chances = ln_read_chances(&per_allele, &genotypes);
        let slips_by_bucket = (0..classes)
            .map(|class| slips_into_each_bucket(slippage, class as i32 - span, read_span))
            .collect();
        Self {
            classes,
            span,
            read_span,
            genotypes,
            per_allele,
            ln_read_chances,
            stay: 1.0 - slippage.level,
            slips_by_bucket,
        }
    }

    /// `reads` scored under every genotype, each log-likelihood read off [`Self::ln_read_chances`].
    fn score(&self, reads: &[u32]) -> ScoredRow {
        // Each genotype's log-likelihood, then its likelihood scaled by the largest, in place.
        let mut scaled: Vec<f64> = self
            .ln_read_chances
            .iter()
            .map(|ln_chances| ln_reads(ln_chances, reads))
            .collect();
        let largest = scaled.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for value in &mut scaled {
            *value = float::exp(*value - largest);
        }
        ScoredRow { largest, scaled }
    }

    /// **The counts of every sample that holds `reads`**, from the row's scaled likelihoods, the
    /// genotypes' `prior` and the samples' summed [`RowWeights`]: each read shared between
    /// "reported its allele's length" and "slipped this far, this way", and each sample's alleles
    /// among the lengths.
    fn counts_of(
        &self,
        prior: &GenotypePrior,
        reads: &[u32],
        row: &ScoredRow,
        weights: RowWeights,
        into: &mut Counts,
    ) {
        let per_allele = &self.per_allele;
        let slot = |class: usize| {
            ((class as i32 - self.span).clamp(-self.read_span, self.read_span) + self.read_span)
                as usize
        };
        for (index, (a, b)) in self.genotypes.iter().enumerate() {
            let as_one_draw = row.scaled[index] * prior.one_draw[index];
            let one = as_one_draw * weights.as_one;
            let two = (row.scaled[index] * prior.two_draws[index]) * weights.as_two;
            let weight = one + two;
            if weight < 1e-12 {
                continue;
            }
            into.draws[*a] += one + two;
            into.draws[*b] += two;
            let shared = as_one_draw * weights.alike;
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

/// **Each genotype's prior at a round's allele-length shares and spread**, split into the part
/// where a sample's two copies are one draw (identical by descent, or drawn alike from the tract)
/// and the part where they are two. A sample's homozygote excess, which weighs the two parts,
/// comes in with the sample ([`Self::add_sample`]).
struct GenotypePrior {
    /// The allele's share for a homozygote, zero for a heterozygote; in `genotypes` order.
    one_draw: Vec<f64>,
    /// The genotype's Hardy–Weinberg chance from the shares; in `genotypes` order.
    two_draws: Vec<f64>,
    concentration: f64,
}

impl GenotypePrior {
    fn of(genotypes: &[(usize, usize)], shares: &[f64], concentration: f64) -> Self {
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
            one_draw,
            two_draws,
            concentration,
        }
    }

    /// `row`'s scaled likelihoods summed under this prior, genotype by genotype in order.
    fn totals(&self, row: &ScoredRow) -> PriorWeightedTotals {
        PriorWeightedTotals {
            one_draw: row
                .scaled
                .iter()
                .zip(&self.one_draw)
                .map(|(l, p)| l * p)
                .sum(),
            two_draws: row
                .scaled
                .iter()
                .zip(&self.two_draws)
                .map(|(l, p)| l * p)
                .sum(),
        }
    }

    /// Adds one sample holding a row with largest log-likelihood `largest` and `totals`, at
    /// homozygote excess `excess`, to `into`, and returns its log-likelihood; `None` for a sample
    /// whose reads no genotype explains.
    fn add_sample(
        &self,
        largest: f64,
        totals: PriorWeightedTotals,
        excess: f64,
        into: &mut RowWeights,
    ) -> Option<f64> {
        // A sample's two copies are one draw if identical by descent, or, failing that, if they
        // drew the same chromosome lineage from the tract's own frequencies: chance
        // 1 / (1 + concentration) under a Dirichlet.
        let from_the_tract = (1.0 - excess) / (1.0 + self.concentration);
        let one_draw_chance = excess + from_the_tract;
        let two_draw_chance = (1.0 - excess) * self.concentration / (1.0 + self.concentration);
        let total = one_draw_chance * totals.one_draw + two_draw_chance * totals.two_draws;
        if total <= 0.0 {
            return None;
        }
        into.as_one += one_draw_chance / total;
        into.as_two += two_draw_chance / total;
        into.alike += from_the_tract / total;
        Some(float::ln(total) + largest)
    }
}

/// **A class's distinct rows of reads, each scored at one stutter**: the [`ReadModel`] and a
/// [`ScoredRow`] for every row of `distinct`, in its order.
///
/// **Why it is kept apart from the shares and the spread.** Scoring the rows is the costly part of
/// a round and depends on the stutter alone. The standard errors take steps in the shares and in F
/// at fixed stutter, and these all read one scoring of the class instead of scoring it again for
/// each step ([`Self::expected_counts`]).
pub(super) struct ScoredDistinctRows<'a> {
    distinct: &'a DistinctRows,
    model: ReadModel,
    scored: Vec<ScoredRow>,
}

impl<'a> ScoredDistinctRows<'a> {
    /// Every row of `distinct`, a table of `evidence`, scored at `slippage`.
    pub(super) fn of(
        evidence: &StratumEvidence,
        distinct: &'a DistinctRows,
        slippage: Slippage,
        allele_span: i32,
    ) -> Self {
        let model = ReadModel::of(evidence, slippage, allele_span);
        // [`ln_reads`] reads a row against its log table bucket by bucket and would drop a bucket
        // the table lacks; every row of a table is built `2 · read_span + 1` buckets long.
        let buckets = (2 * evidence.read_span + 1) as usize;
        assert!(
            distinct.reads.iter().all(|reads| reads.len() == buckets),
            "every distinct row holds one count a bucket"
        );
        let scored = distinct
            .reads
            .par_iter()
            .map(|reads| model.score(reads))
            .collect();
        Self {
            distinct,
            model,
            scored,
        }
    }

    /// **One round's expected counts over the whole class** at the scored stutter, the
    /// allele-length shares `shares` and the spread `concentration`
    /// ([`expected_counts_by_distinct_rows`]).
    pub(super) fn expected_counts(&self, shares: &[f64], concentration: f64) -> Counts {
        let (distinct, model) = (self.distinct, &self.model);
        let prior = GenotypePrior::of(&model.genotypes, shares, concentration);
        let totals: Vec<PriorWeightedTotals> = self
            .scored
            .par_iter()
            .map(|row| prior.totals(row))
            .collect();
        // In fixed chunks summed in order, so that the answer does not depend on the thread pool.
        //
        // A chunk keeps only the distinct rows its samples hold, not one entry for every distinct
        // row: with about 256 chunks, a dense entry each made the partials about 256 times the
        // distinct rows. The research bench on one class of 1,000 loci, 2,000 samples of 30
        // reads, peaked at 6.3 GB with them and 2.9 GB without. The sums are the same to the bit.
        // Each row's chunk sums are still added in chunk order; what is left out is the `+0.0`
        // an untouched chunk added. Adding `+0.0` changes nothing but `-0.0`, and a
        // sum begun at `+0.0` is never `-0.0`: a floating-point sum is `-0.0` only when both its
        // terms are.
        let chunk = (distinct.rows.len() / 256).max(4096);
        let partials: Vec<(Vec<(u32, RowWeights)>, f64)> = distinct
            .rows
            .par_chunks(chunk)
            .map(|rows| {
                // The chunk's distinct rows in the order first met, each with its samples' sums.
                let mut held: Vec<(u32, RowWeights)> = Vec::new();
                let mut place_of: AHashMap<u32, usize> = AHashMap::new();
                let mut log_likelihood = 0.0;
                for (index, excess) in rows {
                    let at = *place_of.entry(*index).or_insert_with(|| {
                        held.push((*index, RowWeights::default()));
                        held.len() - 1
                    });
                    let index = *index as usize;
                    let largest = self.scored[index].largest;
                    if let Some(ll) =
                        prior.add_sample(largest, totals[index], *excess, &mut held[at].1)
                    {
                        log_likelihood += ll;
                    }
                }
                (held, log_likelihood)
            })
            .collect();
        let mut weights = vec![RowWeights::default(); distinct.reads.len()];
        let mut log_likelihood = 0.0;
        for (part, part_ll) in &partials {
            for (index, from) in part {
                let into = &mut weights[*index as usize];
                into.as_one += from.as_one;
                into.as_two += from.as_two;
                into.alike += from.alike;
            }
            log_likelihood += part_ll;
        }
        // Freed before the rows' counts are built.
        drop(partials);
        let per_row: Vec<Counts> = distinct
            .reads
            .par_iter()
            .zip(&self.scored)
            .zip(&weights)
            .map(|((reads, row), weights)| {
                let mut counts = Counts::new(model.classes);
                model.counts_of(&prior, reads, row, *weights, &mut counts);
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
}

/// For the fitting bench's proxy for calling (plan step D6), not for a run.
#[cfg(any(test, feature = "bench-fixtures"))]
impl ScoredDistinctRows<'_> {
    /// **The genotype each sample's weights put highest**: for every (tract, sample) with reads,
    /// in the order of `evidence.tracts` then each tract's samples, the genotype with the largest
    /// weight at the scored stutter, the shares `shares` and the spread `concentration`.
    ///
    /// The weights are the expectation step's own ([`ReadModel::counts_of`]): the row's
    /// likelihood times the genotype's prior, its one-draw and two-draw parts mixed by the
    /// sample's homozygote excess ([`GenotypePrior::add_sample`]). A genotype is its two allele
    /// classes, the first not above the second; a tie goes to the genotype first in
    /// [`genotype_pairs`] order. `None` for a sample whose reads no genotype explains.
    pub(super) fn genotypes_weighted_highest(
        &self,
        shares: &[f64],
        concentration: f64,
    ) -> Vec<Option<(usize, usize)>> {
        let model = &self.model;
        let prior = GenotypePrior::of(&model.genotypes, shares, concentration);
        self.distinct
            .rows
            .par_iter()
            .map(|(index, excess)| {
                let row = &self.scored[*index as usize];
                let mut weights = RowWeights::default();
                prior.add_sample(row.largest, prior.totals(row), *excess, &mut weights)?;
                let mut highest: Option<(usize, f64)> = None;
                for (genotype, scaled) in row.scaled.iter().enumerate() {
                    let weight = scaled
                        * (prior.one_draw[genotype] * weights.as_one
                            + prior.two_draws[genotype] * weights.as_two);
                    if highest.is_none_or(|(_, top)| weight > top) {
                        highest = Some((genotype, weight));
                    }
                }
                highest.map(|(genotype, _)| model.genotypes[genotype])
            })
            .collect()
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
    allele_span: i32,
) -> Counts {
    ScoredDistinctRows::of(evidence, distinct, slippage, allele_span)
        .expected_counts(shares, concentration)
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

/// One sample's reads at a tract, its libraries summed bucket by bucket, written over `reads`,
/// which holds one count a bucket.
fn pool_reads_into(sample: &SampleTractReads, reads: &mut [u32]) {
    reads.fill(0);
    for (_, library_reads) in &sample.by_library {
        for (into, count) in reads.iter_mut().zip(library_reads) {
            *into += count;
        }
    }
}

/// [`pool_reads_into`] into a new row of `buckets` counts.
#[cfg(test)]
fn pooled_reads(sample: &SampleTractReads, buckets: usize) -> Vec<u32> {
    let mut reads = vec![0_u32; buckets];
    pool_reads_into(sample, &mut reads);
    reads
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
/// divided by the largest of them; `ln_read_chances` is [`ln_read_chances`] for `genotypes`.
fn ordered_row(
    reads: &[u32],
    ln_read_chances: &[Vec<f64>],
    genotypes: &[(usize, usize)],
    classes: usize,
) -> Vec<f64> {
    let log_likelihoods: Vec<f64> = ln_read_chances
        .iter()
        .map(|ln_chances| ln_reads(ln_chances, reads))
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
    allele_span: i32,
) -> Vec<Vec<PairCubic>> {
    let classes = allele_classes_in(allele_span);
    let genotypes = genotype_pairs(classes);
    let libraries = one_stutter_for_every_library(evidence);
    let distributions =
        ReadDistributionsByLevel::of(&[slippage], &libraries, allele_span, evidence.read_span);
    let ln_read_chances =
        ln_read_chances(&distributions.read_distribution_by_allele[0], &genotypes);
    let distinct = DistinctRows::of(evidence, homozygote_excess);
    let ordered: Vec<Vec<f64>> = distinct
        .reads
        .par_iter()
        .map(|reads| ordered_row(reads, &ln_read_chances, &genotypes, classes))
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
        allele_span: i32,
    ) -> Self {
        let per_tract =
            pair_cubics_by_tract(evidence, slippage, shares, homozygote_excess, allele_span)
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
        ln_pairs_by_tract(&self.per_tract, float::exp(ln_spread))
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
        ln_pairs_summed_in_tract_order(&self.per_tract, float::exp(ln_spread))
    }
}

/// **Every tract's pairwise log-likelihood at `spread`, summed in the order of the tracts.** The
/// tracts are scored side by side, then added one after another, so the total has the same bits
/// whatever the number of threads. Summed by the thread pool instead, the order of the additions
/// followed how the work was split, and the total moved in its last places from one run to the
/// next: enough to move F's curve, and through it every other number of the fit (spec
/// `str_penalised_fit.md` §3).
fn ln_pairs_summed_in_tract_order(per_tract: &[Vec<PairCubic>], spread: f64) -> f64 {
    ln_pairs_by_tract(per_tract, spread).iter().sum()
}

/// **Each tract's pairwise log-likelihood at `spread`, in the order of the tracts.** The tracts are
/// scored side by side and collected in order. Both the total
/// ([`ln_pairs_summed_in_tract_order`]) and the per-tract values a standard error needs
/// ([`SamplePairs::log_likelihood_by_tract`]) come from here, so the two hold the same numbers in
/// the same order.
fn ln_pairs_by_tract(per_tract: &[Vec<PairCubic>], spread: f64) -> Vec<f64> {
    per_tract
        .par_iter()
        .map(|cubics| scaled_ln_pairs_of_tract(cubics, spread))
        .collect()
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

    /// **A sample's call is the genotype its posterior puts highest**, the posterior computed here
    /// sample by sample as the row-by-row expected counts compute it: each genotype's read
    /// likelihood times its prior, the prior mixing "the two copies are one draw" and "two draws"
    /// by the sample's homozygote excess and the spread. The 20 samples carry excesses from 0 to
    /// 0.98, so a call that dropped either part of the prior, or swapped them, would put some
    /// samples' weight on another genotype. A tie may go either way: the call's posterior must be
    /// the largest to within rounding.
    #[test]
    fn genotypes_weighted_highest_matches_the_argmax_of_a_brute_force_posterior() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let allele_span = 2;
        let samples = 20;
        let evidence = draw_stratum(truth, &spectrum_of(5), 0.5, 0.9, 300, samples, 3, 2, 23);
        let excess: Vec<f64> = (0..samples).map(|s| 0.98 * s as f64 / 19.0).collect();
        let at = Slippage {
            level: 0.12,
            shorter_share: 0.6,
            fall_off: 0.4,
        };
        let shares = [0.1, 0.15, 0.4, 0.2, 0.15];
        let concentration = 0.7;

        let distinct = DistinctRows::of(&evidence, &excess);
        let calls = ScoredDistinctRows::of(&evidence, &distinct, at, allele_span)
            .genotypes_weighted_highest(&shares, concentration);

        let classes = allele_classes_in(allele_span);
        let genotypes = genotype_pairs(classes);
        let read_span = evidence.read_span;
        let libraries = one_stutter_for_every_library(&evidence);
        let distributions = ReadDistributionsByLevel::of(&[at], &libraries, allele_span, read_span);
        let per_allele = &distributions.read_distribution_by_allele[0];
        let mut posteriors: Vec<Vec<f64>> = Vec::new();
        for tract in &evidence.tracts {
            for sample in &tract.samples {
                let excess = excess[sample.sample as usize];
                let reads = pooled_reads(sample, (2 * read_span + 1) as usize);
                let one_draw_chance = excess + (1.0 - excess) / (1.0 + concentration);
                let two_draw_chance = (1.0 - excess) * concentration / (1.0 + concentration);
                let posterior: Vec<f64> = genotypes
                    .iter()
                    .map(|(a, b)| {
                        let mut ln_likelihood = 0.0;
                        add_ln_reads(&mut ln_likelihood, per_allele, &reads, *a, *b);
                        let prior = if a == b {
                            one_draw_chance * shares[*a] + two_draw_chance * shares[*a] * shares[*a]
                        } else {
                            two_draw_chance * 2.0 * shares[*a] * shares[*b]
                        };
                        prior * float::exp(ln_likelihood)
                    })
                    .collect();
                posteriors.push(posterior);
            }
        }

        assert_eq!(calls.len(), posteriors.len());
        let mut heterozygous = 0;
        for (index, (call, posterior)) in calls.iter().zip(&posteriors).enumerate() {
            let highest = posterior.iter().copied().fold(0.0, f64::max);
            assert!(
                highest > 0.0,
                "sample row {index}: no genotype explains the reads"
            );
            let call = call.expect("every sample's reads are explained");
            let called = genotypes
                .iter()
                .position(|genotype| *genotype == call)
                .expect("the call is a genotype");
            assert!(
                posterior[called] >= highest * (1.0 - 1e-9),
                "sample row {index}: called {call:?} at {}, the highest is {highest}",
                posterior[called]
            );
            if call.0 != call.1 {
                heterozygous += 1;
            }
        }
        // Both kinds of call are made, so both parts of the prior were weighed.
        assert!(
            heterozygous > 0 && heterozygous < calls.len(),
            "{heterozygous} of {} calls heterozygous",
            calls.len()
        );
    }

    /// **The table of distinct rows does not depend on the size of a block**: built in blocks of
    /// one tract, of seven and of all the tracts at once, the rows are the same and numbered in the
    /// same order, so the expected counts summed over them are the same bits.
    #[test]
    fn the_distinct_rows_are_the_same_in_blocks_of_any_size() {
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let evidence = draw_stratum(truth, &spectrum_of(5), 2.0, 0.3, 40, 12, 5, 1, 43);
        let excess = [0.3; 12];
        let whole = DistinctRows::in_blocks_of(&evidence, &excess, evidence.tracts.len());
        // Some row must repeat across tracts, so that at one tract a block it is met in two
        // blocks: only then does a merge that skipped de-duplication between blocks fail below.
        let tracts_of_row = |row: u32| {
            whole
                .tract_starts
                .windows(2)
                .filter(|range| {
                    whole.rows[range[0]..range[1]]
                        .iter()
                        .any(|(index, _)| *index == row)
                })
                .count()
        };
        assert!(
            (0..whole.reads.len() as u32).any(|row| tracts_of_row(row) > 1),
            "no row repeats across tracts"
        );
        for tracts_per_block in [1, 7] {
            let blocks = DistinctRows::in_blocks_of(&evidence, &excess, tracts_per_block);
            assert_eq!(
                blocks.reads, whole.reads,
                "{tracts_per_block} tracts a block"
            );
            assert_eq!(blocks.rows, whole.rows, "{tracts_per_block} tracts a block");
            assert_eq!(blocks.tract_starts, whole.tract_starts);
        }
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
            let raw: Vec<f64> = (0..classes)
                .map(|_| float::powi(uniform(), 3) + 1e-4)
                .collect();
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
}
