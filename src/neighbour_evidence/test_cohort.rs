//! **A small drawn cohort for the pass's tests**: biallelic records along one contig, samples whose
//! copies are stretches of a few founder haplotypes, reads drawn at a chosen depth, and the `PL`,
//! `GP` and `AF` calling would write from them. Every draw comes from one seed, so a test sees the
//! same cohort on every run.

use super::records::{SampleColumns, SiteRecord};
use crate::float;

/// **What to draw.**
#[derive(Clone, Copy, Debug)]
pub(crate) struct CohortShape {
    pub samples: usize,
    pub founders: usize,
    pub records: usize,
    /// Bases between neighbouring records.
    pub spacing: u64,
    /// Mean stretch one copy copies from one founder, in bases.
    pub segment_bases: f64,
    /// Whether each sample's two copies are one copy counted twice everywhere.
    pub inbred: bool,
    /// Mean reads a sample a record.
    pub depth: f64,
    /// Per-read chance of showing the other allele.
    pub read_error: f64,
    pub seed: u64,
    /// A sample that gets no reads from a record on: `(sample, record)`.
    pub silent_from: Option<(usize, usize)>,
}

/// **The drawn cohort**: the records as the pass reads them, the true genotype index of every
/// sample at every record, and the genotype calling would write — the argmax of its `GP`.
pub(crate) struct TestCohort {
    pub records: Vec<SiteRecord>,
    pub truth: Vec<Vec<usize>>,
    pub called: Vec<Vec<usize>>,
}

struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut mixed = self.0;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        mixed ^ (mixed >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }

    fn poisson(&mut self, mean: f64) -> usize {
        let limit = float::exp(-mean);
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit {
            product *= self.uniform();
            count += 1;
        }
        count
    }
}

/// **Draws a cohort of `shape`.**
pub(crate) fn test_cohort(shape: CohortShape) -> TestCohort {
    let mut draw = Draw(shape.seed);
    let frequencies: Vec<f64> = (0..shape.records)
        .map(|_| 0.05 + 0.45 * draw.uniform())
        .collect();
    let founder_alleles: Vec<Vec<usize>> = (0..shape.founders)
        .map(|_| {
            frequencies
                .iter()
                .map(|&frequency| usize::from(draw.uniform() < frequency))
                .collect()
        })
        .collect();
    let inbreeding = if shape.inbred { 0.95 } else { 0.0 };
    let switch = 1.0 - float::exp(-(shape.spacing as f64) / shape.segment_bases);

    let mut copies: Vec<[usize; 2]> = (0..shape.samples)
        .map(|_| [draw.below(shape.founders), draw.below(shape.founders)])
        .collect();
    let mut records = Vec::with_capacity(shape.records);
    let mut truth = Vec::with_capacity(shape.records);
    let mut called = Vec::with_capacity(shape.records);
    for record in 0..shape.records {
        let mut phred = Vec::with_capacity(shape.samples);
        let mut posteriors = Vec::with_capacity(shape.samples);
        let mut true_genotypes = Vec::with_capacity(shape.samples);
        let mut called_genotypes = Vec::with_capacity(shape.samples);
        for (sample, founders) in copies.iter_mut().enumerate() {
            for copy in founders.iter_mut() {
                if record > 0 && draw.uniform() < switch {
                    *copy = draw.below(shape.founders);
                }
            }
            if shape.inbred {
                founders[1] = founders[0];
            }
            let alleles = [
                founder_alleles[founders[0]][record],
                founder_alleles[founders[1]][record],
            ];
            true_genotypes.push(alleles[0] + alleles[1]);
            let mut reads = draw.poisson(shape.depth);
            if shape
                .silent_from
                .is_some_and(|(silent, from)| silent == sample && record >= from)
            {
                reads = 0;
            }
            let mut alternative_reads = 0;
            for _ in 0..reads {
                let mut allele = alleles[draw.below(2)];
                if draw.uniform() < shape.read_error {
                    allele = 1 - allele;
                }
                alternative_reads += allele;
            }
            let log_likelihoods: Vec<f64> = [0.0, 0.5, 1.0]
                .iter()
                .map(|&alternative_share| {
                    let alternative = alternative_share * (1.0 - shape.read_error)
                        + (1.0 - alternative_share) * shape.read_error;
                    alternative_reads as f64 * float::ln(alternative)
                        + (reads - alternative_reads) as f64 * float::ln(1.0 - alternative)
                })
                .collect();
            let best = log_likelihoods.iter().copied().fold(f64::MIN, f64::max);
            let sample_phred: Vec<u32> = log_likelihoods
                .iter()
                .map(|&log_likelihood| {
                    (-10.0 * (log_likelihood - best) / std::f64::consts::LN_10).round() as u32
                })
                .collect();
            let frequency = frequencies[record];
            let prior = [
                (1.0 - inbreeding) * (1.0 - frequency) * (1.0 - frequency)
                    + inbreeding * (1.0 - frequency),
                (1.0 - inbreeding) * 2.0 * frequency * (1.0 - frequency),
                (1.0 - inbreeding) * frequency * frequency + inbreeding * frequency,
            ];
            let unscaled: Vec<f64> = sample_phred
                .iter()
                .zip(prior)
                .map(|(&phred, prior)| float::powf(10.0, -f64::from(phred) / 10.0) * prior)
                .collect();
            let total: f64 = unscaled.iter().sum();
            let sample_posteriors: Vec<f64> = unscaled
                .iter()
                .map(|value| (value / total * 1000.0).round() / 1000.0)
                .collect();
            let call = (0..3).fold(0, |best, genotype| {
                if sample_posteriors[genotype] > sample_posteriors[best] {
                    genotype
                } else {
                    best
                }
            });
            called_genotypes.push(call);
            phred.push(sample_phred);
            posteriors.push(sample_posteriors);
        }
        let columns: Vec<Option<SampleColumns<'_>>> = phred
            .iter()
            .zip(&posteriors)
            .map(|(phred, posteriors)| {
                Some(SampleColumns {
                    phred_likelihoods: phred,
                    genotype_posteriors: posteriors,
                })
            })
            .collect();
        records.push(SiteRecord::new(
            1 + record as u64 * shape.spacing,
            &[frequencies[record]],
            &columns,
        ));
        truth.push(true_genotypes);
        called.push(called_genotypes);
    }
    TestCohort {
        records,
        truth,
        called,
    }
}
