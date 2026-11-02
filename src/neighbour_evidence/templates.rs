//! **Which other samples a sample copies from in a window** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §3.3).
//!
//! A sample's templates are the others that disagree with it least, judged only where both have
//! reads, with more data preferred between two that disagree about as little (owner, 2026-11-02;
//! spec §3.3). At each record where both are informative the two samples' disagreement is half
//! their summed allele-share difference, `½ Σ_a |q_self(a) − q_other(a)|`, from 0 to 1. Over the
//! `n` shared records the summed disagreement `d` gives a rate, and the rate is scored
//! pessimistically: the mean of its Beta(`d + 1`, `n − d + 1`) posterior plus two of its standard
//! deviations, [`disagreement_bound`]. The `K` with the lowest scores are the templates.
//!
//! Why only where both have reads: a sample without reads at a record carries the cohort's
//! frequency as its genotype, and counting that as a disagreement would rank a shallow close
//! relative below a deep distant one. Why the pessimistic score: a low-coverage relative's rate rests
//! on fewer records and the score charges for that, so a deep sample wins a near tie — and it is
//! the deep sample whose own data fills a site the target has none at.
//!
//! The comparison runs over at most [`AGREEMENT_RECORDS`] of the sample's informative records,
//! spread evenly by record order, so the choice costs the same in a dense window as in a sparse one
//! and is the same on every run. The background template is not chosen here; it joins every walk
//! (§3.3).

use super::model::GenotypePairs;
use super::records::{SampleRow, SiteRecord};

/// The most templates a walk copies from, the background aside (`K`, §3.3) — what the prototype
/// ran; how accuracy and cost move with it is a measurement (§8).
pub const MAX_TEMPLATES: usize = 16;

/// The most of a sample's informative records its templates are chosen over (§3.3).
pub const AGREEMENT_RECORDS: usize = 500;

/// **The templates sample `sample` copies from in the window `records`**, best first:
/// `min(max_templates, samples − 1)` other samples, ties broken by sample order. `None` when the
/// sample has no informative record in the window — then there is nothing to choose by, and the
/// sample is not re-called there (§3.3).
#[must_use]
pub fn choose_templates(
    records: &[SiteRecord],
    sample: usize,
    max_templates: usize,
) -> Option<Vec<usize>> {
    let compared: Vec<&SiteRecord> = compared_records(records, sample)
        .into_iter()
        .map(|index| &records[index])
        .collect();
    if compared.is_empty() {
        return None;
    }
    let samples = compared[0].sample_count();
    let mut differences: Vec<(f64, usize)> = (0..samples)
        .filter(|&other| other != sample)
        .map(|other| {
            let (mut shared, mut disagreement) = (0_usize, 0.0);
            for record in compared
                .iter()
                .filter(|record| record.row(other) == SampleRow::Informative)
            {
                shared += 1;
                disagreement += 0.5
                    * record
                        .template_shares(sample)
                        .zip(record.template_shares(other))
                        .map(|(own, theirs)| (own - theirs).abs())
                        .sum::<f64>();
            }
            (disagreement_bound(disagreement, shared), other)
        })
        .collect();
    differences.sort_by(|left, right| left.0.total_cmp(&right.0).then(left.1.cmp(&right.1)));
    Some(
        differences
            .into_iter()
            .take(max_templates)
            .map(|(_, other)| other)
            .collect(),
    )
}

/// **The records a sample's templates are chosen over**, as indices into `records`: its
/// informative sites, at most [`AGREEMENT_RECORDS`] of them spread evenly by record order (§3.3).
/// Empty when it has none.
#[must_use]
pub fn compared_records(records: &[SiteRecord], sample: usize) -> Vec<usize> {
    let informative: Vec<usize> = (0..records.len())
        .filter(|&index| {
            records[index].is_site() && records[index].row(sample) == SampleRow::Informative
        })
        .collect();
    if informative.len() <= AGREEMENT_RECORDS {
        informative
    } else {
        (0..AGREEMENT_RECORDS)
            .map(|draw| informative[draw * informative.len() / AGREEMENT_RECORDS])
            .collect()
    }
}

/// **The haplotypes sample `sample` copies from in a phasing round** (phasing spec §2.3), best
/// first, as `(sample, copy)`: at most `max_templates`, at most one of each other sample.
///
/// A haplotype's score is the rate at which it carries an allele the sample's genotype lacks — at
/// each record `Σ_a h(a) · P(genotype without a)`, `h` the haplotype and the genotype probabilities
/// the sample's posterior from the previous round, `previous(record)` — summed over the `compared`
/// records where the haplotype's sample is informative too, and scored by
/// [`disagreement_bound`]. Each other sample offers its better haplotype (copy 0 on a tie); the
/// offers are ranked by score, then sample, then copy. `haplotype(other, copy, record)` reads the
/// previous round's haplotypes.
///
/// One of each sample because without it the slots fill with copies of the commonest haplotype:
/// 871 wrong genotypes against 689 on the prototype's first cohort (spec §2.3).
#[must_use]
pub fn choose_haplotypes<'a>(
    records: &[SiteRecord],
    sample: usize,
    compared: &[usize],
    previous: &dyn Fn(usize) -> &'a [f64],
    haplotype: &dyn Fn(usize, usize, usize) -> &'a [f64],
    max_templates: usize,
) -> Vec<(usize, usize)> {
    let Some(&first) = compared.first() else {
        return Vec::new();
    };
    let samples = records[first].sample_count();
    let mut pairs = GenotypePairs::default();
    // P(genotype without a) at each compared record, for every allele.
    let lacking: Vec<Vec<f64>> = compared
        .iter()
        .map(|&index| {
            let alleles = records[index].allele_count;
            let posterior = previous(index);
            (0..alleles)
                .map(|allele| {
                    pairs
                        .of(alleles)
                        .iter()
                        .zip(posterior)
                        .filter(|&(&(a, b), _)| a != allele && b != allele)
                        .map(|(_, &probability)| probability)
                        .sum()
                })
                .collect()
        })
        .collect();
    let mut offers: Vec<(f64, usize, usize)> = (0..samples)
        .filter(|&other| other != sample)
        .map(|other| {
            let mut shared = 0_usize;
            let mut disagreement = [0.0_f64; 2];
            for (position, &index) in compared.iter().enumerate() {
                if records[index].row(other) != SampleRow::Informative {
                    continue;
                }
                shared += 1;
                for (copy, total) in disagreement.iter_mut().enumerate() {
                    *total += haplotype(other, copy, index)
                        .iter()
                        .zip(&lacking[position])
                        .map(|(carried, lacked)| carried * lacked)
                        .sum::<f64>();
                }
            }
            let scores = disagreement.map(|total| disagreement_bound(total, shared));
            if scores[1] < scores[0] {
                (scores[1], other, 1)
            } else {
                (scores[0], other, 0)
            }
        })
        .collect();
    offers.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then(left.1.cmp(&right.1))
            .then(left.2.cmp(&right.2))
    });
    offers
        .into_iter()
        .take(max_templates)
        .map(|(_, other, copy)| (other, copy))
        .collect()
}

/// **A pessimistic disagreement rate** from `disagreement` summed over `shared` records: the mean
/// of the rate's Beta(`disagreement + 1`, `shared − disagreement + 1`) posterior, from a flat prior,
/// plus two of its standard deviations. With no shared record it is ½ + 2·√(1/12), above any rate
/// measured on records.
#[must_use]
pub fn disagreement_bound(disagreement: f64, shared: usize) -> f64 {
    let shared = shared as f64;
    let mean = (disagreement + 1.0) / (shared + 2.0);
    mean + 2.0 * (mean * (1.0 - mean) / (shared + 3.0)).sqrt()
}

#[cfg(test)]
mod tests;
