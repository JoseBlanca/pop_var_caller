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
    let informative: Vec<&SiteRecord> = records
        .iter()
        .filter(|record| record.is_site() && record.row(sample) == SampleRow::Informative)
        .collect();
    if informative.is_empty() {
        return None;
    }
    let compared: Vec<&SiteRecord> = if informative.len() <= AGREEMENT_RECORDS {
        informative
    } else {
        (0..AGREEMENT_RECORDS)
            .map(|draw| informative[draw * informative.len() / AGREEMENT_RECORDS])
            .collect()
    };
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
