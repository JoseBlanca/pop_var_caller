//! **Tracts the stutter fit leaves out because their reads come from more copies than the
//! reference holds** (spec `parameter_prepass_ssr.md` §4.6.1).
//!
//! Copies of different lengths put one tract's reads at several lengths, and the fit has one place
//! for that: slippage. On HG002 chr20 at 30× a satellite III array beside the centromere made a
//! 5-base motif of 9 copies fit at 23 slips in 100 reads, where tracts known to be homozygous show
//! none. Two rules leave such tracts out:
//!
//! - **the depth rule**: a tract whose crossing reads, summed over every sample, far exceed what
//!   one copy of it gives ([`leave_out_deep_tracts`]). It runs once, before the fit, and reads no
//!   fitted number.
//! - **the stutter-tail rule**: a tract where, against the genotype each sample's reads favour, far more
//!   reads sat off both alleles than the class's slip rate allows ([`tracts_in_the_stutter_tail`]).
//!   It runs between fits of a motif length, against the numbers just fitted, until it flags none.
//!
//! A tract left out is still called. It only stops teaching the fit what slippage looks like.

use super::expected_counts::reads_off_the_genotype_the_reads_favour;
use super::*;
use crate::genetics::lgamma;

/// **Which of the two rules run, and at what settings.** [`production`](Self::production) is what
/// `estimate-parameters` runs; [`none`](Self::none) fits every tract the guard kept, as the fit did
/// before 2026-11-02.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TractExclusion {
    pub depth: Option<DepthRule>,
    pub stutter_tail: Option<StutterTailRule>,
}

impl TractExclusion {
    /// Both rules at the settings of spec §4.6.1.
    #[must_use]
    pub fn production() -> Self {
        Self {
            depth: Some(DepthRule::production()),
            stutter_tail: Some(StutterTailRule::production()),
        }
    }

    /// Neither rule.
    #[must_use]
    pub fn none() -> Self {
        Self {
            depth: None,
            stutter_tail: None,
        }
    }
}

/// **The depth rule's settings.** A tract is left out when its crossing reads exceed the larger
/// of `times_expected` times what one copy gives and the count one copy exceeds by chance only
/// `chance` of the time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DepthRule {
    /// 2: the satellite array sits at 2.8 times chr20's typical depth.
    pub times_expected: f64,
    /// 0.001: at a few reads a tract, twice the expected count is reached by chance too often —
    /// at 5× it left out 5 tracts in 100 — so the cutoff rises to the Poisson 0.999 quantile,
    /// 2.6 times the expected count at 5 reads and 3.3 times at 3.
    pub chance: f64,
}

impl DepthRule {
    #[must_use]
    pub fn production() -> Self {
        Self {
            times_expected: 2.0,
            chance: 0.001,
        }
    }
}

/// **The stutter-tail rule's settings.** A tract is left out when its off-length reads are at least
/// `times_rate` times what the class's slip rate predicts **and** that many or more would occur by
/// slippage less than `chance` of the time; tested in rounds, at most `max_rounds` of them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StutterTailRule {
    /// 2. Without it, a tract summed over thousands of samples is flagged for slipping only a
    /// little more than its class: the binomial test grows stricter with every read.
    pub times_rate: f64,
    /// 0.001.
    pub chance: f64,
    /// 10. A tract left out stays out, so the rounds end; a motif length still flagging tracts
    /// after this many is reported as not settled.
    pub max_rounds: usize,
}

impl StutterTailRule {
    #[must_use]
    pub fn production() -> Self {
        Self {
            times_rate: 2.0,
            chance: 0.001,
            max_rounds: 10,
        }
    }
}

/// **What the stutter-tail rule did to one motif length.**
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StutterTailOutcome {
    /// Tracts it left out, over the motif length's classes.
    pub tracts: u64,
    /// Rounds that left a tract out.
    pub rounds: usize,
    /// Whether the last fit's test flagged none. False when the rule's rounds ran out first: the
    /// numbers returned are then from a fit some of whose tracts the rule would still flag.
    pub settled: bool,
}

/// **What the two rules left out of a run's fit, and what they were set to** — what the
/// parameters file's `[repeat_tracts_left_out]` records (spec `parameter_prepass_ssr.md` §4.6.1).
#[derive(Debug, Clone, PartialEq)]
pub struct TractsLeftOut {
    /// The rules as the fit ran them.
    pub exclusion: TractExclusion,
    /// The depth rule's typical depth; `None` where it did not run or no tract was crossed.
    pub typical_depth: Option<f64>,
    /// Per motif length the tail rule ran on, in increasing order, what it did.
    pub stutter_tail_by_period: Vec<(u8, StutterTailOutcome)>,
    /// One entry a class, in the order the fit was handed them.
    pub by_stratum: Vec<StratumTractsLeftOut>,
}

/// **One class's tracts: how many the fit used, and how many each rule left out.**
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StratumTractsLeftOut {
    pub stratum: Stratum,
    /// Tracts the last fit read.
    pub fitted: u64,
    /// Left out before the fit by the guard: more than one read in ten of those differing from the
    /// reference by a non-whole number of repeats.
    pub over_guard_threshold: u64,
    pub too_deep: u64,
    pub in_stutter_tail: u64,
}

impl StratumTractsLeftOut {
    /// What `evidence` says, once both rules have run on it.
    #[must_use]
    pub fn of(evidence: &StratumEvidence) -> Self {
        Self {
            stratum: evidence.stratum,
            fitted: evidence.tracts.len() as u64,
            over_guard_threshold: evidence.tracts_over_guard_threshold,
            too_deep: evidence.tracts_too_deep,
            in_stutter_tail: evidence.tracts_in_stutter_tail,
        }
    }
}

/// What the depth rule did over a run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeepTractsLeftOut {
    /// **How many reads the typical tract would have if every read reaching it crossed it**: the
    /// median, over every tract of a stratum some read crossed, of its crossing reads divided by
    /// its stratum's crossing share. `None` when no tract had a stratum with a crossing share.
    pub typical_depth: Option<f64>,
    /// Tracts left out, over every stratum.
    pub tracts: u64,
}

/// **Of the reads that reached a stratum's tracts, the share that crossed one**, over every tract
/// the stratum holds, those the guard left out included; `None` when no read reached any.
///
/// Every tract of a stratum is the same length, so one share serves all of them, and a
/// duplicated tract adds reads to both counts, so it barely moves the share. Reads the guard
/// caught reached the tract and did not report a whole-repeat length, so they count as reaching.
fn crossing_share(evidence: &StratumEvidence) -> Option<f64> {
    let crossing = evidence
        .tracts
        .iter()
        .map(TractReads::crossing_reads)
        .sum::<u64>()
        + evidence.reads_crossing_tracts_over_guard_threshold;
    let reaching = crossing + evidence.guard_reads + evidence.reads_reaching_not_crossing;
    (crossing > 0).then(|| crossing as f64 / reaching as f64)
}

/// **The count a Poisson variable of mean `mean` exceeds with chance below `chance`**: the
/// smallest `c` with `P(X > c) < chance`.
///
/// Summed term by term in logs, so a mean of tens of thousands of reads — a cohort's sum — does
/// not underflow. Called once a stratum, never once a tract.
fn poisson_count_exceeded_by_chance(mean: f64, chance: f64) -> u64 {
    if mean <= 0.0 {
        return 0;
    }
    let ln_mean = float::ln(mean);
    let mut at_most = 0.0_f64;
    let mut count = 0_u64;
    loop {
        at_most += float::exp(count as f64 * ln_mean - mean - lgamma(count as f64 + 1.0));
        if 1.0 - at_most < chance {
            return count;
        }
        count += 1;
    }
}

/// **Leave out every tract whose crossing reads far exceed what one copy gives** (spec
/// `parameter_prepass_ssr.md` §4.6.1, rule 1), counting each stratum's in
/// [`StratumEvidence::tracts_too_deep`].
///
/// A tract's expected count, if it is one copy, is its stratum's [crossing share](crossing_share)
/// times the run's typical depth. It is left out when its crossing reads exceed the larger of
/// `rule.times_expected` times that and the count one copy exceeds by chance only `rule.chance`
/// of the time. The crossing share is what puts strata of different lengths on one scale: a long
/// tract is crossed by fewer reads than a short one at the same depth.
pub fn leave_out_deep_tracts(
    strata: &mut [StratumEvidence],
    rule: &DepthRule,
) -> DeepTractsLeftOut {
    let shares: Vec<Option<f64>> = strata.iter().map(crossing_share).collect();
    let mut depths: Vec<f64> = strata
        .iter()
        .zip(&shares)
        .filter_map(|(evidence, share)| share.map(|share| (evidence, share)))
        .flat_map(|(evidence, share)| {
            evidence
                .tracts
                .iter()
                .map(move |tract| tract.crossing_reads() as f64 / share)
        })
        .collect();
    if depths.is_empty() {
        return DeepTractsLeftOut {
            typical_depth: None,
            tracts: 0,
        };
    }
    depths.sort_unstable_by(f64::total_cmp);
    let middle = depths.len() / 2;
    let typical_depth = if depths.len() % 2 == 1 {
        depths[middle]
    } else {
        0.5 * (depths[middle - 1] + depths[middle])
    };
    let mut left_out = 0;
    for (evidence, share) in strata.iter_mut().zip(&shares) {
        let Some(share) = share else { continue };
        let expected = share * typical_depth;
        let cutoff = (rule.times_expected * expected)
            .max(poisson_count_exceeded_by_chance(expected, rule.chance) as f64);
        let before = evidence.tracts.len();
        evidence
            .tracts
            .retain(|tract| tract.crossing_reads() as f64 <= cutoff);
        let removed = (before - evidence.tracts.len()) as u64;
        evidence.tracts_too_deep += removed;
        left_out += removed;
    }
    DeepTractsLeftOut {
        typical_depth: Some(typical_depth),
        tracts: left_out,
    }
}

/// **`P(X ≥ at_least)` for `X ~ Binomial(trials, chance)`**, summed upward from `at_least` in logs.
///
/// Called only where `at_least` is above the mean, where the terms fall from the first, so the
/// sum stops once a term no longer moves it.
fn binomial_upper_tail(trials: u64, chance: f64, at_least: u64) -> f64 {
    if at_least == 0 {
        return 1.0;
    }
    if at_least > trials {
        return 0.0;
    }
    if chance >= 1.0 {
        return 1.0;
    }
    let n = trials as f64;
    let ln_choose = |k: f64| lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
    let (ln_p, ln_q) = (float::ln(chance), float::ln(1.0 - chance));
    let mut total = 0.0;
    for k in at_least..=trials {
        let k = k as f64;
        let term = float::exp(ln_choose(k) + k * ln_p + (n - k) * ln_q);
        total += term;
        if term < 1e-17 * total {
            break;
        }
    }
    total.min(1.0)
}

/// **Which of a class's tracts the stutter-tail rule leaves out at the numbers just fitted**:
/// indices into `evidence.tracts`, in order (spec `parameter_prepass_ssr.md` §4.6.1, rule 2).
///
/// Per tract, over its samples summed, `N` crossing reads and `K` off-length reads against the
/// genotype each sample's reads favour alone ([`reads_off_the_genotype_the_reads_favour`]). Left out
/// when `K ≥ rule.times_rate × level × N` and `P(Binomial(N, level) ≥ K) < rule.chance`, with
/// `level` the class's slip rate.
///
/// **Why the reads alone, with no prior.** Two lengths in one sample are then always two alleles,
/// so a tract polymorphic between samples is never counted as stutter. A prior from the class's
/// allele shares did count it, in a cohort at a few reads a sample: on tomato2's 59 accessions it
/// left out 438 tracts, nearly all polymorphic, against 3 with the reads alone, with identical calls
/// on HG002 and the stuttering simulation (`ng_stutter_fit_outliers_2026-11-02.md` §5, §7). The cost:
/// a duplication that puts exactly two lengths in a sample looks heterozygous, so at a few reads a
/// sample only a third length shows; the depth rule is what catches the rest.
///
/// **The test leans towards keeping a tract.** A read that slips from one allele onto the other's
/// length is not counted, so on a tract the model describes `K` falls short of `level × N` on
/// average.
///
/// **With `multipliers`, each read group's reads at its own rate** (spec `str_library_stutter.md`
/// §5): the genotype is read off the reads scored at each read group's rate, and `level` is the
/// read-weighted mean of the read groups' rates at the tract, `Σᵣ nᵣ λᵣ / N`. A sum of reads
/// slipping at different rates is less spread than a binomial at their mean (Hoeffding 1956), so
/// the binomial at the mean overstates the chance of many off-length reads, and the test still
/// leans towards keeping. Without them, every read is at the class's rate.
pub(super) fn tracts_in_the_stutter_tail(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    slippage: Slippage,
    allele_span: i32,
    rule: &StutterTailRule,
    multipliers: Option<&[f64]>,
) -> Vec<usize> {
    let level = slippage.level;
    if !(level > 0.0) {
        return Vec::new();
    }
    reads_off_the_genotype_the_reads_favour(
        evidence,
        homozygote_excess,
        slippage,
        allele_span,
        multipliers,
    )
    .into_iter()
    .enumerate()
    .filter(|(_, reads)| {
        let (crossing, off_length) = (reads.crossing, reads.off_length);
        match reads.expected_slipped {
            None => {
                off_length as f64 >= rule.times_rate * level * crossing as f64
                    && off_length > 0
                    && binomial_upper_tail(crossing, level, off_length) < rule.chance
            }
            Some(expected) => {
                let mean_rate = expected / crossing as f64;
                off_length as f64 >= rule.times_rate * expected
                    && off_length > 0
                    && mean_rate > 0.0
                    && binomial_upper_tail(crossing, mean_rate, off_length) < rule.chance
            }
        }
    })
    .map(|(tract, _)| tract)
    .collect()
}

/// PROBE ONLY (tail-rule-probe): one line a tract of `evidence` for this round's test, appended to
/// the file `TAIL_DUMP` names.
#[allow(clippy::too_many_arguments)]
pub(super) fn probe_dump_round(
    round: usize,
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    slippage: Slippage,
    allele_span: i32,
    rule: &StutterTailRule,
    multipliers: Option<&[f64]>,
    flagged: &[usize],
) {
    use std::io::Write;
    let Ok(path) = std::env::var("TAIL_DUMP") else {
        return;
    };
    let reads = reads_off_the_genotype_the_reads_favour(
        evidence,
        homozygote_excess,
        slippage,
        allele_span,
        multipliers,
    );
    let origin = super::PROBE_TRACT_ORIGIN.lock().unwrap();
    let mut text = String::new();
    for (index, (tract, counted)) in evidence.tracts.iter().zip(&reads).enumerate() {
        let locus = origin
            .as_ref()
            .and_then(|map| map.get(&(tract.samples.as_ptr() as usize)))
            .map_or(-1, |(_, _, locus)| *locus as i64);
        // Summed over samples: reads in each length bucket; samples with reads; samples whose reads
        // sit at three or more lengths; samples whose second-commonest length holds at least 3
        // reads and 30 in 100 of the sample's.
        let mut histogram: Vec<u64> = Vec::new();
        let (mut samples, mut three_lengths, mut two_strong) = (0_u64, 0_u64, 0_u64);
        for sample in &tract.samples {
            let mut summed: Vec<u64> = Vec::new();
            for (_, counts) in &sample.by_read_group {
                if summed.len() < counts.len() {
                    summed.resize(counts.len(), 0);
                }
                for (bucket, count) in counts.iter().enumerate() {
                    summed[bucket] += u64::from(*count);
                }
            }
            if histogram.len() < summed.len() {
                histogram.resize(summed.len(), 0);
            }
            for (bucket, count) in summed.iter().enumerate() {
                histogram[bucket] += count;
            }
            let total: u64 = summed.iter().sum();
            if total == 0 {
                continue;
            }
            samples += 1;
            if summed.iter().filter(|count| **count > 0).count() >= 3 {
                three_lengths += 1;
            }
            let mut sorted = summed.clone();
            sorted.sort_unstable_by(|a, b| b.cmp(a));
            if sorted.len() > 1 && sorted[1] >= 3 && sorted[1] as f64 >= 0.3 * total as f64 {
                two_strong += 1;
            }
        }
        let mean_rate = counted
            .expected_slipped
            .map_or(slippage.level, |expected| expected / counted.crossing.max(1) as f64);
        let p = binomial_upper_tail(counted.crossing, mean_rate.min(1.0), counted.off_length);
        let histogram: Vec<String> = histogram.iter().map(u64::to_string).collect();
        text.push_str(&format!(
            "{round}\t{}\t{}\t{index}\t{locus}\t{}\t{}\t{}\t{}\t{}\t{p:e}\t{}\t{samples}\t{three_lengths}\t{two_strong}\t{}\n",
            evidence.stratum.period,
            evidence.stratum.reference_repeats,
            slippage.level,
            mean_rate,
            counted.crossing,
            counted.off_length,
            rule.times_rate * mean_rate * counted.crossing as f64,
            u8::from(flagged.binary_search(&index).is_ok()),
            histogram.join(","),
        ));
    }
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .expect("TAIL_DUMP opens");
    file.write_all(text.as_bytes()).expect("TAIL_DUMP writes");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_poisson_cutoff_matches_the_spec_table() {
        // Spec §4.6.1's table: the count one copy exceeds with chance below 1 in 1,000.
        for (mean, cutoff) in [(3.0, 10), (5.0, 13), (10.0, 21), (30.0, 48), (300.0, 355)] {
            assert_eq!(
                poisson_count_exceeded_by_chance(mean, 0.001),
                cutoff,
                "mean {mean}"
            );
        }
    }

    #[test]
    fn the_binomial_tail_is_the_sum_of_its_terms() {
        // 3 reads at 0.01: P(X ≥ 2) = 3·0.01²·0.99 + 0.01³.
        let expected = 3.0 * 0.0001 * 0.99 + 0.000_001;
        assert!((binomial_upper_tail(3, 0.01, 2) - expected).abs() < 1e-15);
        assert_eq!(binomial_upper_tail(3, 0.01, 0), 1.0);
        assert_eq!(binomial_upper_tail(3, 0.01, 4), 0.0);
        // A deep tract: 10,000 reads at 0.02, mean 200; 300 or more is far out.
        let far = binomial_upper_tail(10_000, 0.02, 300);
        assert!(far > 0.0 && far < 1e-9, "{far}");
    }

    // ---- what the tests below are built from ----

    /// The recorded buckets, `-4 … +4`.
    const BUCKETS: usize = (2 * RECORDED_OFFSET_RANGE + 1) as usize;

    /// One sample's reads, `(offset, reads)` pairs, as a bucket row.
    fn row(reads: &[(i32, u32)]) -> Vec<u32> {
        let mut row = vec![0; BUCKETS];
        for (offset, count) in reads {
            row[(offset + RECORDED_OFFSET_RANGE) as usize] += count;
        }
        row
    }

    /// A tract whose sample `s` read `rows[s]`, each sample its own read group.
    fn tract(rows: &[Vec<u32>]) -> TractReads {
        TractReads {
            samples: rows
                .iter()
                .enumerate()
                .map(|(sample, reads)| SampleTractReads {
                    sample: sample as u32,
                    by_read_group: vec![(ReadGroupIndex(sample as u32), reads.clone())],
                })
                .collect(),
        }
    }

    /// A class at `repeats` repeats of a 2-base motif, holding `tracts` read by `samples` samples.
    fn class(repeats: u64, tracts: Vec<TractReads>, samples: usize) -> StratumEvidence {
        StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: repeats,
            },
            tracts,
            read_span: RECORDED_OFFSET_RANGE,
            groups: 1,
            read_groups: bench_fixtures::one_read_group_a_sample(samples),
            tracts_over_guard_threshold: 0,
            reads_crossing_tracts_over_guard_threshold: 0,
            tracts_too_deep: 0,
            tracts_in_stutter_tail: 0,
            reads_reaching_not_crossing: 0,
            guard_reads: 0,
            bases_compared: 0,
            mismatching_bases: 0,
        }
    }

    /// `count` one-sample tracts with every read at the reference length.
    fn reference_tracts(count: usize, reads: u32) -> Vec<TractReads> {
        (0..count).map(|_| tract(&[row(&[(0, reads)])])).collect()
    }

    // ---- the depth rule ----

    #[test]
    fn a_long_tract_twice_as_deep_as_its_stratum_is_left_out_through_the_crossing_share() {
        // Short tracts: every read reaching one crosses it, 30 reads each. Long tracts: a third
        // of the reads reaching one cross it, so one copy shows 10, and a duplicated one 25.
        // Compared with the run's 30 crossing reads the duplicated long tract looks ordinary;
        // compared with its own stratum's 10 it is two and a half times one copy.
        let short = class(6, reference_tracts(99, 30), 1);
        let mut long_tracts = reference_tracts(20, 10);
        long_tracts.push(tract(&[row(&[(0, 25)])]));
        let mut long = class(20, long_tracts, 1);
        long.reads_reaching_not_crossing = 2 * (20 * 10 + 25);
        let mut strata = vec![short, long];

        let deep = leave_out_deep_tracts(&mut strata, &DepthRule::production());

        assert_eq!(deep.typical_depth, Some(30.0));
        assert_eq!(strata[0].tracts_too_deep, 0);
        assert_eq!(strata[1].tracts_too_deep, 1, "the duplicated long tract");
        assert_eq!(strata[1].tracts.len(), 20);
    }

    #[test]
    fn at_three_reads_the_cutoff_rises_to_what_one_copy_reaches_by_chance() {
        // One copy gives 3 reads. Twice that is 6, which one copy exceeds about 1 time in 30;
        // the cutoff is the 0.999 quantile, 10. So 7 reads stay and 11 go.
        let mut tracts = reference_tracts(100, 3);
        tracts.push(tract(&[row(&[(0, 7)])]));
        tracts.push(tract(&[row(&[(0, 11)])]));
        let mut strata = vec![class(6, tracts, 1)];

        let deep = leave_out_deep_tracts(&mut strata, &DepthRule::production());

        assert_eq!(deep.tracts, 1);
        assert!(
            strata[0]
                .tracts
                .iter()
                .any(|tract| tract.crossing_reads() == 7),
            "the 7-read tract is within chance of one copy"
        );
    }

    // ---- the stutter-tail rule ----

    /// A class's stutter in these tests: one read in 100 slips.
    fn stutter(level: f64) -> Slippage {
        Slippage {
            level,
            shorter_share: 0.5,
            fall_off: 0.3,
        }
    }

    fn flagged(evidence: &StratumEvidence, level: f64) -> Vec<usize> {
        let samples = evidence.read_groups.len();
        tracts_in_the_stutter_tail(
            evidence,
            &vec![0.0; samples],
            stutter(level),
            2,
            &StutterTailRule::production(),
            None,
        )
    }

    /// **A tract read by a read group that slips three times as often as typical is tested at that
    /// read group's rate** (spec `str_library_stutter.md` §5): 1,000 reads, 60 of them one repeat
    /// short. At the class's rate of 2 in 100 that is three times what it predicts and far past
    /// chance, so a test at the class's rate leaves the tract out; at the read group's own 6 in 100
    /// it is what slippage gives, and the tract stays. Sample 0, the typical read group, reads
    /// another tract.
    #[test]
    fn a_tract_read_by_a_stuttering_read_group_is_tested_at_its_rate() {
        let evidence = class(
            10,
            vec![
                tract(&[row(&[(0, 100)])]),
                TractReads {
                    samples: vec![SampleTractReads {
                        sample: 1,
                        by_read_group: vec![(ReadGroupIndex(1), row(&[(0, 940), (-1, 60)]))],
                    }],
                },
            ],
            2,
        );
        let test = |multipliers: Option<&[f64]>| {
            tracts_in_the_stutter_tail(
                &evidence,
                &[0.0, 0.0],
                stutter(0.02),
                2,
                &StutterTailRule::production(),
                multipliers,
            )
        };
        assert_eq!(test(None), vec![1], "at the class's rate");
        assert!(
            test(Some(&[1.0, 3.0])).is_empty(),
            "at the read group's own"
        );
    }

    #[test]
    fn a_heterozygote_is_not_stutter_and_a_third_length_is() {
        // Half the reads one repeat long is a heterozygote, and the genotype the reads favour says
        // so: none of its reads is off-length. A third length, a third of the reads two repeats
        // short, is off both alleles of any genotype.
        let evidence = class(
            10,
            vec![
                tract(&[row(&[(0, 15), (1, 15)])]),
                tract(&[row(&[(0, 10), (1, 10), (-2, 10)])]),
            ],
            1,
        );
        assert_eq!(flagged(&evidence, 0.01), vec![1]);
    }

    #[test]
    fn a_duplication_three_reads_a_sample_shows_only_summed_over_samples() {
        // Every sample reads three lengths, one read each: a third length no genotype of two
        // alleles carries, so one read in three is off-length. One sample's 1 in 3 happens by
        // slippage about 3 times in 100; fifty samples' 50 in 150 never does.
        let a_sample = row(&[(0, 1), (2, 1), (-2, 1)]);
        let one = class(10, vec![tract(std::slice::from_ref(&a_sample))], 1);
        assert!(flagged(&one, 0.01).is_empty(), "one sample's three reads");
        let fifty = class(10, vec![tract(&vec![a_sample; 50])], 50);
        assert_eq!(
            flagged(&fifty, 0.01),
            vec![0],
            "fifty samples' reads summed"
        );
    }

    #[test]
    fn a_tract_polymorphic_between_inbred_samples_is_not_stutter() {
        // Fifty inbred samples: 25 read the reference length three times, and 25 carry a length
        // two repeats long and show it in their one read. Each sample's reads sit at one length,
        // so none is off-length. A genotype drawn from a class whose shares put that length at 1 in
        // 1,000 gave each one-read sample the reference homozygote, and its read counted as
        // stutter: 25 of 100 reads (tomato2, report §5).
        let mut samples = vec![row(&[(0, 3)]); 25];
        samples.extend(vec![row(&[(2, 1)]); 25]);
        let evidence = class(10, vec![tract(&samples)], 50);
        assert!(flagged(&evidence, 0.01).is_empty());
    }

    /// **A drawn class with two tiers of tracts no genotype explains**: 300 ordinary tracts
    /// slipping 1 read in 100, 30 tracts with reads split evenly over three lengths, and 30 with 5
    /// of their 30 reads a repeat short. The second tier is improbable only at the true rate, not
    /// at the rate the first tier inflates. One sample at 30 reads.
    fn two_tiers_of_tracts_no_genotype_explains() -> Vec<StratumEvidence> {
        let truth = Slippage {
            level: 0.01,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let spectrum = [0.002, 0.008, 0.98, 0.008, 0.002];
        let mut evidence =
            bench_fixtures::draw_stratum(truth, &spectrum, 50.0, 0.3, 300, 1, 30, 2, 17);
        let one_sample = |reads: &[(i32, u32)]| {
            let mut row = vec![0; 5];
            for (offset, count) in reads {
                row[(offset + 2) as usize] += count;
            }
            TractReads {
                samples: vec![SampleTractReads {
                    sample: 0,
                    by_read_group: vec![(ReadGroupIndex(0), row)],
                }],
            }
        };
        for _ in 0..30 {
            evidence
                .tracts
                .push(one_sample(&[(-2, 10), (0, 10), (2, 10)]));
            evidence.tracts.push(one_sample(&[(0, 25), (-1, 5)]));
        }
        vec![evidence]
    }

    /// The class's fitted slip rate, and the tracts the stutter-tail rule left out, when the rule
    /// may run `max_rounds` rounds; `None` fits every tract.
    fn level_and_tracts_left_out(max_rounds: Option<usize>) -> (f64, u64) {
        use super::super::penalised_fit::{PenalisedFitConfig, fit_repeat_tracts};
        let mut strata = two_tiers_of_tracts_no_genotype_explains();
        let exclusion = TractExclusion {
            depth: None,
            stutter_tail: max_rounds.map(|max_rounds| StutterTailRule {
                max_rounds,
                ..StutterTailRule::production()
            }),
        };
        let outcomes = fit_repeat_tracts(
            &mut strata,
            &[0.3],
            &PenalisedFitConfig::production().with_allele_span(2),
            &exclusion,
        )
        .outcomes;
        let level = outcomes[0].slippage()[0].expect("fitted").level;
        (level, strata[0].tracts_in_stutter_tail)
    }

    #[test]
    fn the_rounds_reach_tracts_the_first_rate_hid() {
        let (every_tract, _) = level_and_tracts_left_out(None);
        let (one_round, after_one) = level_and_tracts_left_out(Some(1));
        let (settled, after_all) = level_and_tracts_left_out(Some(10));
        eprintln!(
            "every tract {every_tract:.4}; one round {one_round:.4} ({after_one} out); \
             settled {settled:.4} ({after_all} out)"
        );
        assert!(after_one >= 30 && after_one < 60, "{after_one}");
        assert!(after_all >= 60, "{after_all}");
        assert!(settled < 0.015, "{settled}");
        assert!(one_round > 1.2 * settled, "{one_round} against {settled}");
        assert!(every_tract > one_round, "{every_tract} against {one_round}");
    }

    #[test]
    fn a_deep_tract_slipping_one_and_a_half_times_its_class_is_kept() {
        // 10,000 reads at a slip rate of 2 in 100 predict 200 off-length. 300 is far beyond chance
        // — a binomial test alone would flag it — but only 1.5 times the rate; 500 is 2.5 times.
        let deep = |off: u32| tract(&[row(&[(0, 10_000 - off), (-2, off)])]);
        let evidence = class(10, vec![deep(300), deep(500)], 1);
        assert_eq!(flagged(&evidence, 0.02), vec![1]);
    }
}
