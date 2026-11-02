//! **Which other samples a sample copies from in a window** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §3.3).
//!
//! A sample's templates are the others whose allele shares look most like its own over the window:
//! for each other sample, the summed difference `Σ_records Σ_a |q_self(a) − q_other(a)|` over the
//! sample's own informative records, and the `K` with the smallest sums. The comparison runs over
//! at most [`AGREEMENT_RECORDS`] of those records, spread evenly by record order, so the choice
//! costs the same in a dense window as in a sparse one and is the same on every run. The background
//! template is not chosen here; it joins every walk (§3.3).

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
            let difference = compared
                .iter()
                .map(|record| {
                    record
                        .template_shares(sample)
                        .zip(record.template_shares(other))
                        .map(|(own, theirs)| (own - theirs).abs())
                        .sum::<f64>()
                })
                .sum();
            (difference, other)
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

#[cfg(test)]
mod tests;
