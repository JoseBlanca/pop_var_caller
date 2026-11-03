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
use smallvec::SmallVec;

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
    /// **Per read group, its share of `exact` and `slipped`**, indexed as
    /// [`StratumEvidence::read_groups`] lists them — what a read group's stutter multiplier is
    /// fitted from (`str_library_stutter.md` §3.4). **Empty when the rows were pooled**
    /// ([`RowKeying::Pooled`]), where a sample's read groups are no longer told apart.
    pub read_groups: Vec<ReadGroupSlips>,
}

/// **One read group's expected reads in a class**: those that reported their allele's length and
/// those that slipped, each read shared between the two at its own read group's slip rate.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ReadGroupSlips {
    pub exact: f64,
    pub slipped: f64,
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
        if self.read_groups.len() < other.read_groups.len() {
            self.read_groups
                .resize(other.read_groups.len(), ReadGroupSlips::default());
        }
        for (into, from) in self.read_groups.iter_mut().zip(&other.read_groups) {
            into.exact += from.exact;
            into.slipped += from.slipped;
        }
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
/// tract, the count in each read-length bucket with its read groups pooled. Two samples whose reads
/// fall in the same buckets in the same numbers hold the same row. On tomato2's 59 samples,
/// 754,601 rows hold 8,780 distinct ones, class by class.
///
/// **Why the round can use it.** A sample's reads enter a round of counting only through its row
/// and its homozygote excess. The costly part, each genotype's likelihood of the reads, depends on
/// the row alone; the excess only mixes two sums of it. So the genotypes are scored once for each
/// distinct row, and each sample adds three numbers to its row's totals ([`RowWeights`]).
///
/// **Keyed two ways** ([`RowKeying`]). Pooled, a row is one count a bucket, as above. Keyed by read
/// group, a row is each of the sample's read groups with its own counts, so two samples hold the same
/// row only if their reads came through the same read groups in the same numbers: a fit scoring
/// each read group at its own slip rate must not share a row's score between read groups. On
/// tomato2 that is 61,474 distinct rows where pooling gives 9,199 (`str_library_stutter.md` §3.2).
pub(super) struct DistinctRows {
    /// How each row is keyed, and so how [`Self::reads`] is laid out.
    keying: RowKeying,
    /// How many read-length buckets a row counts, `2 × read_span + 1`.
    buckets: usize,
    /// **The multiplier each read group's reads are scored at**, indexed as
    /// [`StratumEvidence::read_groups`]; `None` to score at the evidence's own. Read only by a
    /// keyed table: a pooled one scores every read group at the class's stutter.
    multipliers: Option<Vec<f64>>,
    /// Each distinct row once, in the order first met. **Pooled**, one count a bucket. **Keyed by
    /// read group**, one part for each of the sample's read groups in increasing order, each part
    /// the read group's index followed by its count in each bucket ([`Self::parts`]).
    reads: Vec<Vec<u32>>,
    /// Per (tract, sample) with reads, in the order of `evidence.tracts` then each tract's
    /// samples: which distinct row it holds, and the sample's homozygote excess.
    rows: Vec<(u32, f64)>,
    /// Where each tract's rows start in `rows`, and one past the last.
    tract_starts: Vec<usize>,
}

/// **How a class's rows of reads are told apart** ([`DistinctRows`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum RowKeying {
    /// A sample's read groups summed bucket by bucket: every read group is scored at the class's
    /// one stutter.
    Pooled,
    /// Each read group of the sample apart, so that each can be scored at its own slip rate.
    ByReadGroup,
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
        Self::keyed(evidence, homozygote_excess, RowKeying::Pooled)
    }

    /// [`of`](Self::of), with each row keyed as `keying` says.
    pub(super) fn keyed(
        evidence: &StratumEvidence,
        homozygote_excess: &[f64],
        keying: RowKeying,
    ) -> Self {
        let rows_in_class: usize = evidence.tracts.iter().map(|t| t.samples.len()).sum();
        let tracts_per_block = (ROWS_PER_BLOCK * evidence.tracts.len())
            .checked_div(rows_in_class)
            .unwrap_or(64)
            .clamp(1, 64);
        Self::in_blocks_of(evidence, homozygote_excess, tracts_per_block, keying)
    }

    /// **A table keyed by read group whose rows are scored at `multipliers`**, one a read group of
    /// `evidence`, rather than at the evidence's own: what the fit's rounds score at while they
    /// move the multipliers.
    pub(super) fn keyed_at(
        evidence: &StratumEvidence,
        homozygote_excess: &[f64],
        multipliers: &[f64],
    ) -> Self {
        assert_eq!(
            multipliers.len(),
            evidence.read_groups.len(),
            "one multiplier a read group"
        );
        let mut table = Self::keyed(evidence, homozygote_excess, RowKeying::ByReadGroup);
        table.multipliers = Some(multipliers.to_vec());
        table
    }

    /// **The parts of a distinct row**: each with the read group it came through, `None` for a
    /// pooled row, and its count in each bucket.
    fn parts<'r>(
        &self,
        row: &'r [u32],
    ) -> impl Iterator<Item = (Option<ReadGroupIndex>, &'r [u32])> {
        let pooled = self.keying == RowKeying::Pooled;
        let stride = if pooled {
            self.buckets
        } else {
            self.buckets + 1
        };
        row.chunks_exact(stride).map(move |part| {
            if pooled {
                (None, part)
            } else {
                (Some(ReadGroupIndex(part[0])), &part[1..])
            }
        })
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
        keying: RowKeying,
    ) -> Self {
        let buckets = (2 * evidence.read_span + 1) as usize;
        // Per block: its distinct rows in the order first met, and each (tract, sample)'s place
        // among them with the sample's excess. Every sample's reads are written into one scratch
        // row, pooled or by read group, and looked up by value; only a row not met before is
        // copied.
        let blocks: Vec<RowsOfABlock> = evidence
            .tracts
            .par_chunks(tracts_per_block)
            .map(|tracts| {
                let mut index_of: AHashMap<Vec<u32>, u32> = AHashMap::new();
                let mut rows = Vec::new();
                let mut pooled = vec![0_u32; buckets];
                for sample in tracts.iter().flat_map(|tract| &tract.samples) {
                    match keying {
                        RowKeying::Pooled => pool_reads_into(sample, &mut pooled),
                        RowKeying::ByReadGroup => read_groups_apart_into(sample, &mut pooled),
                    }
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
            keying,
            buckets,
            multipliers: None,
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
/// distribution at every slip rate a read group is scored at, the log chance of a read in each
/// bucket under each genotype, and where the slipped reads land. Everything in it depends on the
/// stutter and the read groups' multipliers alone, so the rows scored with it serve every set of
/// shares and spread ([`ScoredDistinctRows`]).
///
/// **One slip rate for a pooled table, one for each distinct multiplier for a keyed one.** Read
/// groups at the same multiplier share their rate's distributions, so a keyed table at every
/// multiplier one builds one, as a pooled table does.
struct ReadModel {
    classes: usize,
    span: i32,
    read_span: i32,
    genotypes: Vec<(usize, usize)>,
    /// One for each distinct slip rate; a pooled table's is the class's.
    rates: Vec<ReadModelAtRate>,
    /// Per read group of the evidence, which of `rates` its reads are scored at. Empty for a
    /// pooled table, whose rows are scored at `rates[0]`.
    rate_of_read_group: Vec<usize>,
    /// Where an allele's slipped reads land, by direction. A slip rate scales every entry alike and
    /// only their ratios are read, so it serves every rate.
    slips_by_bucket: Vec<Vec<SlipsIntoBucket>>,
}

/// **The part of a [`ReadModel`] that depends on the slip rate**: each allele's read-length
/// distribution, its log chances under each genotype, and the chance a read reports its allele.
struct ReadModelAtRate {
    per_allele: Vec<Vec<f64>>,
    /// [`ln_read_chances`] of `per_allele`: a row per genotype in `genotypes` order.
    ln_read_chances: Vec<Vec<f64>>,
    stay: f64,
}

impl ReadModelAtRate {
    fn of(slippage: Slippage, per_allele: Vec<Vec<f64>>, genotypes: &[(usize, usize)]) -> Self {
        let ln_read_chances = ln_read_chances(&per_allele, genotypes);
        Self {
            per_allele,
            ln_read_chances,
            stay: 1.0 - slippage.level,
        }
    }
}

/// A row's parts, each with the slip rate it is scored at: almost every row has one or two.
type PartsAtRates<'m, 'r> = SmallVec<[(&'m ReadModelAtRate, Option<ReadGroupIndex>, &'r [u32]); 2]>;

impl ReadModel {
    fn of(
        evidence: &StratumEvidence,
        slippage: Slippage,
        allele_span: i32,
        keying: RowKeying,
        multipliers: Option<&[f64]>,
    ) -> Self {
        let classes = allele_classes_in(allele_span);
        let genotypes = genotype_pairs(classes);
        let span = allele_span;
        let read_span = evidence.read_span;
        let (rates, rate_of_read_group) = match keying {
            RowKeying::Pooled => {
                let read_groups = one_stutter_for_every_read_group(evidence);
                let mut distributions =
                    ReadDistributionsByLevel::of(&[slippage], &read_groups, span, read_span);
                let per_allele = distributions.read_distribution_by_allele.swap_remove(0);
                (
                    vec![ReadModelAtRate::of(slippage, per_allele, &genotypes)],
                    Vec::new(),
                )
            }
            RowKeying::ByReadGroup => {
                debug_assert!(
                    evidence
                        .read_groups
                        .iter()
                        .all(|read_group| read_group.slippage_group == 0),
                    "the repeat-tract fit scores one slippage group"
                );
                // Each read group's multiplier: the table's, or the evidence's own.
                let of_read_group: Vec<f64> = multipliers.map_or_else(
                    || {
                        evidence
                            .read_groups
                            .iter()
                            .map(|read_group| read_group.stutter_multiplier)
                            .collect()
                    },
                    <[f64]>::to_vec,
                );
                // Each distinct multiplier once, by its bits, in increasing order.
                let mut multipliers: Vec<u64> = of_read_group
                    .iter()
                    .map(|multiplier| multiplier.to_bits())
                    .collect();
                multipliers
                    .sort_unstable_by(|a, b| f64::from_bits(*a).total_cmp(&f64::from_bits(*b)));
                multipliers.dedup();
                let rates = multipliers
                    .iter()
                    .map(|multiplier| {
                        let scored = Slippage {
                            level: read_group_level(slippage.level, f64::from_bits(*multiplier)),
                            ..slippage
                        };
                        let per_allele = (-span..=span)
                            .map(|allele| scored.read_probabilities(allele, read_span))
                            .collect();
                        ReadModelAtRate::of(scored, per_allele, &genotypes)
                    })
                    .collect();
                let rate_of_read_group = of_read_group
                    .iter()
                    .map(|multiplier| {
                        let bits = multiplier.to_bits();
                        // PANIC-FREE: every read group's multiplier was put in the list above.
                        multipliers
                            .binary_search_by(|held| {
                                f64::from_bits(*held).total_cmp(&f64::from_bits(bits))
                            })
                            .expect("every read group's multiplier is listed")
                    })
                    .collect();
                (rates, rate_of_read_group)
            }
        };
        let slips_by_bucket = (0..classes)
            .map(|class| slips_into_each_bucket(slippage, class as i32 - span, read_span))
            .collect();
        Self {
            classes,
            span,
            read_span,
            genotypes,
            rates,
            rate_of_read_group,
            slips_by_bucket,
        }
    }

    /// How many read groups a keyed table's counts are kept for; zero for a pooled one.
    fn read_groups(&self) -> usize {
        self.rate_of_read_group.len()
    }

    /// `row`'s parts, a row of `distinct`, each with the slip rate its read group is scored at.
    fn parts_at_rates<'m, 'r>(
        &'m self,
        distinct: &DistinctRows,
        row: &'r [u32],
    ) -> PartsAtRates<'m, 'r> {
        distinct
            .parts(row)
            .map(|(read_group, counts)| {
                let rate =
                    read_group.map_or(0, |read_group| self.rate_of_read_group[read_group.get()]);
                (&self.rates[rate], read_group, counts)
            })
            .collect()
    }

    /// **A row scored under every genotype**: each genotype's log-likelihood is the sum over the
    /// row's parts of that part's reads read off its own rate's [`ReadModelAtRate::ln_read_chances`].
    /// A row of one part, every row of a pooled table, takes that part's sum as it is, so a pooled
    /// row scores to the bit as before read groups were told apart.
    fn score(&self, parts: &PartsAtRates<'_, '_>) -> ScoredRow {
        // Each genotype's log-likelihood, then its likelihood scaled by the largest, in place.
        let mut scaled = self.log_likelihoods(parts);
        let largest = scaled.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for value in &mut scaled {
            *value = float::exp(*value - largest);
        }
        ScoredRow { largest, scaled }
    }

    /// **Each genotype's log-likelihood of a row**, in `genotypes` order: the sum over the row's
    /// parts of the part's reads read off its own rate's log chances, a row of one part taking that
    /// part's sum as it is.
    fn log_likelihoods(&self, parts: &PartsAtRates<'_, '_>) -> Vec<f64> {
        (0..self.genotypes.len())
            .map(|genotype| {
                let mut each = parts
                    .iter()
                    .map(|(rate, _, counts)| ln_reads(&rate.ln_read_chances[genotype], counts));
                let first = each.next().unwrap_or(0.0);
                each.fold(first, |total, part| total + part)
            })
            .collect()
    }

    /// **The counts of every sample that holds a row**, from the row's `parts`, its scaled
    /// likelihoods, the genotypes' `prior` and the samples' summed [`RowWeights`]: each read shared
    /// between "reported its allele's length" and "slipped this far, this way" at its own read
    /// group's slip rate, and each sample's alleles among the lengths.
    ///
    /// **Each part's exactly-read and slipped reads are also added to `by_part`**, one entry a part
    /// in the order of `parts`, for a keyed table's per-read-group counts.
    fn counts_of(
        &self,
        prior: &GenotypePrior,
        parts: &PartsAtRates<'_, '_>,
        row: &ScoredRow,
        weights: RowWeights,
        into: &mut Counts,
        by_part: &mut [ReadGroupSlips],
    ) {
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
            for ((rate, _, reads), part) in parts.iter().zip(by_part.iter_mut()) {
                let per_allele = &rate.per_allele;
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
                            rate.stay / per_allele[allele][bucket]
                        } else {
                            0.0
                        };
                        let exact = share * unslipped.min(1.0);
                        let slipped = share - exact;
                        into.exact += exact;
                        into.slipped += slipped;
                        part.exact += exact;
                        part.slipped += slipped;
                        // Which way the slipped reads went and how far, split by the model's own
                        // terms: an end bucket can hold slips of both directions from an allele
                        // at or past the recorded range.
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
        let model = ReadModel::of(
            evidence,
            slippage,
            allele_span,
            distinct.keying,
            distinct.multipliers.as_deref(),
        );
        // [`ln_reads`] reads a row against its log table bucket by bucket and would drop a bucket
        // the table lacks; every part of a row is built `2 · read_span + 1` buckets long.
        let buckets = (2 * evidence.read_span + 1) as usize;
        let stride = match distinct.keying {
            RowKeying::Pooled => buckets,
            RowKeying::ByReadGroup => buckets + 1,
        };
        assert!(
            distinct.buckets == buckets
                && distinct
                    .reads
                    .iter()
                    .all(|reads| !reads.is_empty() && reads.len() % stride == 0),
            "every part of a distinct row holds one count a bucket"
        );
        let scored = distinct
            .reads
            .par_iter()
            .map(|reads| model.score(&model.parts_at_rates(distinct, reads)))
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
        // Each row's counts, and each of its parts' exactly-read and slipped reads with the read
        // group it came through; a pooled row's part names none.
        type PartSlips = SmallVec<[(Option<ReadGroupIndex>, ReadGroupSlips); 2]>;
        let per_row: Vec<(Counts, PartSlips)> = distinct
            .reads
            .par_iter()
            .zip(&self.scored)
            .zip(&weights)
            .map(|((reads, row), weights)| {
                let parts = model.parts_at_rates(distinct, reads);
                let mut by_part: SmallVec<[ReadGroupSlips; 2]> =
                    smallvec::smallvec![ReadGroupSlips::default(); parts.len()];
                let mut counts = Counts::new(model.classes);
                model.counts_of(&prior, &parts, row, *weights, &mut counts, &mut by_part);
                let slips = parts
                    .iter()
                    .zip(by_part)
                    .map(|((_, read_group, _), slips)| (*read_group, slips))
                    .collect();
                (counts, slips)
            })
            .collect();
        let mut counts = Counts::new(model.classes);
        counts.read_groups = vec![ReadGroupSlips::default(); model.read_groups()];
        for (row, parts) in &per_row {
            counts.add(row);
            for (read_group, slips) in parts {
                if let Some(read_group) = read_group {
                    let into = &mut counts.read_groups[read_group.get()];
                    into.exact += slips.exact;
                    into.slipped += slips.slipped;
                }
            }
        }
        counts.log_likelihood = log_likelihood;
        counts
    }
}

/// For the fitting bench's proxy for calling (plan step D6), not for a run.
#[cfg(test)]
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

/// The stutter-tail rule's genotype (`parameter_prepass_ssr.md` §4.6.1).
impl ScoredDistinctRows<'_> {
    /// **The genotype each sample's reads favour alone**: for every (tract, sample) with reads, in
    /// the order of `evidence.tracts` then each tract's samples, the genotype whose likelihood of
    /// the reads is highest at the scored stutter, with no prior. A tie goes to the genotype first
    /// in [`genotype_pairs`] order.
    ///
    /// **What the stutter-tail rule counts off-length reads against** (`parameter_prepass_ssr.md`
    /// §4.6.1). With no prior, a sample's reads at two lengths are always two alleles, so a tract
    /// polymorphic between samples is never counted as stutter, whatever its class's shares say.
    pub(super) fn genotypes_the_reads_favour(&self) -> Vec<Option<(usize, usize)>> {
        let model = &self.model;
        self.distinct
            .rows
            .par_iter()
            .map(|(index, _)| {
                let row = &self.scored[*index as usize];
                let mut highest: Option<(usize, f64)> = None;
                for (genotype, scaled) in row.scaled.iter().enumerate() {
                    if *scaled > 0.0 && highest.is_none_or(|(_, top)| *scaled > top) {
                        highest = Some((genotype, *scaled));
                    }
                }
                highest.map(|(genotype, _)| model.genotypes[genotype])
            })
            .collect()
    }
}

/// **Each tract's reads, and how many of them sat off both alleles of the genotype their
/// sample's reads favour**, summed over the tract's samples: one `(crossing reads, off-length reads)` a tract, in
/// the order of `evidence.tracts`. What the stutter-tail rule tests (`parameter_prepass_ssr.md`
/// §4.6.1).
///
/// A sample's genotype is the one its reads favour alone, with no prior
/// ([`ScoredDistinctRows::genotypes_the_reads_favour`]). **Not the one the fit weighs highest**: that
/// one's prior is the class's allele shares, not the tract's own allele frequencies, so at a few
/// reads a sample carrying a length rare in its class was given the class's common genotype and its
/// reads counted as stutter. On tomato2's 59 inbred accessions that left out 438 tracts, 86 to 90 in
/// 100 of them polymorphic between accessions; with the reads alone, 3
/// (`ng_stutter_fit_outliers_2026-11-02.md` §5). **A read is off-length when its bucket is neither
/// allele's**, and the end buckets saturate, so a read and an allele both at or past the recorded
/// range share a bucket and the read is not counted. A sample whose reads no genotype explains adds
/// its reads and none of them as off-length.
pub(super) fn reads_off_the_genotype_the_reads_favour(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    slippage: Slippage,
    allele_span: i32,
    multipliers: Option<&[f64]>,
) -> Vec<TractOffLengthReads> {
    // The table pairs each row with its sample's excess; with no prior, nothing reads it.
    let distinct = match multipliers {
        None => DistinctRows::of(evidence, homozygote_excess),
        Some(multipliers) => DistinctRows::keyed_at(evidence, homozygote_excess, multipliers),
    };
    let scored = ScoredDistinctRows::of(evidence, &distinct, slippage, allele_span);
    let genotypes = scored.genotypes_the_reads_favour();
    let read_span = evidence.read_span;
    let slot = |class: usize| {
        ((class as i32 - allele_span).clamp(-read_span, read_span) + read_span) as usize
    };
    (0..evidence.tracts.len())
        .map(|tract| {
            let rows = distinct.tract_starts[tract]..distinct.tract_starts[tract + 1];
            let mut crossing = 0_u64;
            let mut off_length = 0_u64;
            let mut expected_slipped = multipliers.map(|_| 0.0);
            for row in rows {
                let reads = &distinct.reads[distinct.rows[row].0 as usize];
                for (read_group, counts) in distinct.parts(reads) {
                    let part_reads = counts.iter().map(|count| u64::from(*count)).sum::<u64>();
                    crossing += part_reads;
                    if let (Some(expected), Some(read_group), Some(multipliers)) =
                        (expected_slipped.as_mut(), read_group, multipliers)
                    {
                        *expected += part_reads as f64
                            * read_group_level(slippage.level, multipliers[read_group.get()]);
                    }
                    if let Some((a, b)) = genotypes[row] {
                        let (a, b) = (slot(a), slot(b));
                        off_length += counts
                            .iter()
                            .enumerate()
                            .filter(|(bucket, _)| *bucket != a && *bucket != b)
                            .map(|(_, count)| u64::from(*count))
                            .sum::<u64>();
                    }
                }
            }
            TractOffLengthReads {
                crossing,
                off_length,
                expected_slipped,
            }
        })
        .collect()
}

/// **One tract's reads for the stutter-tail rule**, summed over its samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct TractOffLengthReads {
    /// Reads that crossed the tract.
    pub crossing: u64,
    /// Of those, reads at neither allele of the genotype their sample's reads favour.
    pub off_length: u64,
    /// **The slipped reads the read groups' own rates predict**, `Σᵣ nᵣ λᵣ`, when the rows were
    /// scored at read groups' multipliers; `None` when every read was scored at the class's rate.
    pub expected_slipped: Option<f64>,
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

/// Every read group of `evidence` scored at the class's one stutter: one slippage group, every
/// multiplier one. The counting fit does not fit read groups apart.
fn one_stutter_for_every_read_group(evidence: &StratumEvidence) -> Vec<ReadGroupStutter> {
    evidence
        .read_groups
        .iter()
        .map(|read_group| ReadGroupStutter {
            slippage_group: 0,
            stutter_multiplier: 1.0,
            ..*read_group
        })
        .collect()
}

/// One sample's reads at a tract, its read groups summed bucket by bucket, written over `reads`,
/// which holds one count a bucket.
fn pool_reads_into(sample: &SampleTractReads, reads: &mut [u32]) {
    reads.fill(0);
    for (_, read_group_reads) in &sample.by_read_group {
        for (into, count) in reads.iter_mut().zip(read_group_reads) {
            *into += count;
        }
    }
}

/// **One sample's reads at a tract, each read group apart**, written over `row`: for each read group
/// in increasing order, its index and then its count in each bucket, the layout of a keyed
/// [`DistinctRows`] row.
fn read_groups_apart_into(sample: &SampleTractReads, row: &mut Vec<u32>) {
    row.clear();
    for (read_group, counts) in &sample.by_read_group {
        row.push(read_group.0);
        row.extend_from_slice(counts);
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
fn ordered_row(log_likelihoods: &[f64], genotypes: &[(usize, usize)], classes: usize) -> Vec<f64> {
    let largest = log_likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut ordered = vec![0.0; classes * classes];
    for ((first, second), value) in genotypes.iter().zip(log_likelihoods) {
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
    multipliers: Option<&[f64]>,
) -> Vec<Vec<PairCubic>> {
    let classes = allele_classes_in(allele_span);
    let distinct = match multipliers {
        None => DistinctRows::of(evidence, homozygote_excess),
        Some(multipliers) => DistinctRows::keyed_at(evidence, homozygote_excess, multipliers),
    };
    let model = ReadModel::of(
        evidence,
        slippage,
        allele_span,
        distinct.keying,
        multipliers,
    );
    let genotypes = &model.genotypes;
    let ordered: Vec<Vec<f64>> = distinct
        .reads
        .par_iter()
        .map(|reads| {
            let log_likelihoods = model.log_likelihoods(&model.parts_at_rates(&distinct, reads));
            ordered_row(&log_likelihoods, genotypes, classes)
        })
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
    /// at the allele-length shares `shares`. With `multipliers`, one a read group, each read
    /// group's reads are scored at the class's rate times its multiplier.
    #[must_use]
    pub fn of(
        evidence: &StratumEvidence,
        slippage: Slippage,
        shares: &[f64],
        homozygote_excess: &[f64],
        allele_span: i32,
        multipliers: Option<&[f64]>,
    ) -> Self {
        let per_tract = pair_cubics_by_tract(
            evidence,
            slippage,
            shares,
            homozygote_excess,
            allele_span,
            multipliers,
        )
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
        let read_groups = one_stutter_for_every_read_group(&evidence);
        let distributions =
            ReadDistributionsByLevel::of(&[at], &read_groups, allele_span, read_span);
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

    /// **The table of distinct rows does not depend on the size of a block**, pooled or keyed by
    /// read group: built in blocks of one tract, of seven and of all the tracts at once, the rows
    /// are the same and numbered in the same order, so the expected counts summed over them are
    /// the same bits.
    #[test]
    fn the_distinct_rows_are_the_same_in_blocks_of_any_size() {
        for keying in [RowKeying::Pooled, RowKeying::ByReadGroup] {
            distinct_rows_are_the_same_in_blocks_of_any_size(keying);
        }
    }

    fn distinct_rows_are_the_same_in_blocks_of_any_size(keying: RowKeying) {
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        // Two read groups a sample, so that keyed rows differ from pooled ones.
        let evidence = two_read_groups_a_sample(draw_stratum(
            truth,
            &spectrum_of(5),
            2.0,
            0.3,
            40,
            12,
            5,
            1,
            43,
        ));
        let excess = [0.3; 12];
        let whole = DistinctRows::in_blocks_of(&evidence, &excess, evidence.tracts.len(), keying);
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
            let blocks = DistinctRows::in_blocks_of(&evidence, &excess, tracts_per_block, keying);
            assert_eq!(
                blocks.reads, whole.reads,
                "{tracts_per_block} tracts a block"
            );
            assert_eq!(blocks.rows, whole.rows, "{tracts_per_block} tracts a block");
            assert_eq!(blocks.tract_starts, whole.tract_starts);
        }
    }

    /// **Each sample's reads split between two read groups**, `2s` and `2s + 1`, every bucket's
    /// reads half to each (the odd one to the first), every read group at a multiplier of one in
    /// slippage group 0. The pooled reads are the evidence's own.
    fn two_read_groups_a_sample(mut evidence: StratumEvidence) -> StratumEvidence {
        let samples = evidence.read_groups.len();
        for row in evidence
            .tracts
            .iter_mut()
            .flat_map(|tract| &mut tract.samples)
        {
            let buckets = (2 * evidence.read_span + 1) as usize;
            let pooled = pooled_reads(row, buckets);
            let first: Vec<u32> = pooled.iter().map(|count| count.div_ceil(2)).collect();
            let second: Vec<u32> = pooled
                .iter()
                .zip(&first)
                .map(|(all, one)| all - one)
                .collect();
            let sample = row.sample;
            row.by_read_group = [(2 * sample, first), (2 * sample + 1, second)]
                .into_iter()
                .filter(|(_, counts)| counts.iter().any(|count| *count > 0))
                .map(|(read_group, counts)| (ReadGroupIndex(read_group), counts))
                .collect();
        }
        evidence.read_groups = (0..2 * samples as u32)
            .map(|read_group| ReadGroupStutter::in_group(ReadGroupId(read_group), 0))
            .collect();
        evidence.check_read_groups();
        evidence
    }

    /// **At every multiplier one, a table keyed by read group counts what the pooled table
    /// counts**, to within rounding — the read groups' reads are scored at one rate, and a sum of
    /// their log chances is the log chance of their sum — and its read groups' exactly-read and
    /// slipped reads add up to the class's. A pooled table keeps no read group's counts.
    #[test]
    fn a_keyed_table_at_every_multiplier_one_counts_as_the_pooled_one() {
        let truth = Slippage {
            level: 0.06,
            shorter_share: 0.7,
            fall_off: 0.3,
        };
        let allele_span = 2;
        let samples = 15;
        let evidence = two_read_groups_a_sample(draw_stratum(
            truth,
            &spectrum_of(5),
            0.8,
            0.5,
            200,
            samples,
            4,
            2,
            7,
        ));
        let excess: Vec<f64> = (0..samples).map(|s| 0.9 * s as f64 / 14.0).collect();
        let at = Slippage {
            level: 0.1,
            shorter_share: 0.6,
            fall_off: 0.35,
        };
        let shares = [0.1, 0.2, 0.4, 0.2, 0.1];
        let counts_of = |keying| {
            let distinct = DistinctRows::keyed(&evidence, &excess, keying);
            expected_counts_by_distinct_rows(&evidence, &distinct, at, &shares, 0.6, allele_span)
        };
        let pooled = counts_of(RowKeying::Pooled);
        let keyed = counts_of(RowKeying::ByReadGroup);
        assert!(
            DistinctRows::keyed(&evidence, &excess, RowKeying::ByReadGroup).distinct_rows()
                > DistinctRows::of(&evidence, &excess).distinct_rows(),
            "the split must make keyed rows that pooling merges"
        );
        let close = |what: &str, a: f64, b: f64| {
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1.0),
                "{what}: pooled {a}, keyed {b}"
            );
        };
        close(
            "log-likelihood",
            pooled.log_likelihood,
            keyed.log_likelihood,
        );
        close("exact", pooled.exact, keyed.exact);
        close("slipped", pooled.slipped, keyed.slipped);
        close("shorter", pooled.shorter, keyed.shorter);
        close("longer", pooled.longer, keyed.longer);
        close("steps", pooled.steps, keyed.steps);
        close("stepped", pooled.stepped, keyed.stepped);
        close(
            "shared in tract",
            pooled.shared_in_tract,
            keyed.shared_in_tract,
        );
        close(
            "pairs",
            pooled.pairs_not_by_descent,
            keyed.pairs_not_by_descent,
        );
        for (p, k) in pooled.draws.iter().zip(&keyed.draws) {
            close("draws", *p, *k);
        }
        assert!(pooled.read_groups.is_empty());
        assert_eq!(keyed.read_groups.len(), 2 * samples);
        let summed = |pick: fn(&ReadGroupSlips) -> f64| keyed.read_groups.iter().map(pick).sum();
        close("read groups' exact", keyed.exact, summed(|r| r.exact));
        close("read groups' slipped", keyed.slipped, summed(|r| r.slipped));
    }

    /// **Each read group's reads are scored at its own slip rate**: one sample through two read
    /// groups at multipliers one and three, its copies one draw (an excess of one), so its
    /// likelihood is `Σₐ sₐ Πᵣ Π_bucket pᵣ(bucket | a)^reads`, each read group's `p` at the class
    /// rate times its multiplier. The log-likelihood and each read group's slipped reads equal that
    /// sum and its posterior, worked out here directly; a fit scoring both read groups at one rate
    /// fails both.
    #[test]
    fn each_read_group_is_scored_at_its_own_slip_rate() {
        let allele_span = 2;
        let read_span = 2;
        let buckets = (2 * read_span + 1) as usize;
        let class = Slippage {
            level: 0.05,
            shorter_share: 0.75,
            fall_off: 0.3,
        };
        let multipliers = [1.0, 3.0];
        // Read group 0 reads the reference length and one repeat short; read group 1 mostly the
        // reference length, one short and two long.
        let reads = [vec![0, 1, 6, 0, 0], vec![0, 2, 5, 0, 1]];
        let mut evidence = draw_stratum(class, &spectrum_of(5), 1.0, 1.0, 1, 1, 3, 2, 1);
        evidence.tracts = vec![TractReads {
            samples: vec![SampleTractReads {
                sample: 0,
                by_read_group: vec![
                    (ReadGroupIndex(0), reads[0].clone()),
                    (ReadGroupIndex(1), reads[1].clone()),
                ],
            }],
        }];
        evidence.read_groups = multipliers
            .iter()
            .enumerate()
            .map(|(index, multiplier)| ReadGroupStutter {
                stutter_multiplier: *multiplier,
                ..ReadGroupStutter::in_group(ReadGroupId(index as u32), 0)
            })
            .collect();
        evidence.check_read_groups();
        let shares = [0.1, 0.2, 0.4, 0.2, 0.1];
        let distinct = DistinctRows::keyed(&evidence, &[1.0], RowKeying::ByReadGroup);
        let counts = expected_counts_by_distinct_rows(
            &evidence,
            &distinct,
            class,
            &shares,
            0.6,
            allele_span,
        );

        // Directly: each allele's chance of the reads, each read group at its own rate.
        let per_allele = |multiplier: f64, allele: usize| {
            Slippage {
                level: read_group_level(class.level, multiplier),
                ..class
            }
            .read_probabilities(allele as i32 - allele_span, read_span)
        };
        let likelihood = |allele: usize| -> f64 {
            (0..2)
                .map(|group| {
                    let p = per_allele(multipliers[group], allele);
                    (0..buckets)
                        .map(|bucket| float::powi(p[bucket], reads[group][bucket] as i32))
                        .product::<f64>()
                })
                .product()
        };
        let total: f64 = (0..5)
            .map(|allele| shares[allele] * likelihood(allele))
            .sum();
        assert!(
            (counts.log_likelihood - float::ln(total)).abs() < 1e-10,
            "log-likelihood {} against {}",
            counts.log_likelihood,
            float::ln(total)
        );
        for group in 0..2 {
            let stay = 1.0 - read_group_level(class.level, multipliers[group]);
            let slipped: f64 = (0..5)
                .map(|allele| {
                    let posterior = shares[allele] * likelihood(allele) / total;
                    let p = per_allele(multipliers[group], allele);
                    let own = (allele as i32 - allele_span + read_span) as usize;
                    posterior
                        * (0..buckets)
                            .map(|bucket| {
                                let reads = f64::from(reads[group][bucket]);
                                if bucket == own {
                                    reads * (1.0 - (stay / p[bucket]).min(1.0))
                                } else {
                                    reads
                                }
                            })
                            .sum::<f64>()
                })
                .sum();
            let got = counts.read_groups[group].slipped;
            assert!(
                (got - slipped).abs() < 1e-10,
                "read group {group}: {got} slipped reads against {slipped}"
            );
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
