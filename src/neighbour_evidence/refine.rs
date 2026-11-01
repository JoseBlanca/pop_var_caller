//! **Re-calling every sample in one window** (`doc/devel/ng/spec/neighbour_evidence_pass.md` §2.2,
//! §3, §4.2).
//!
//! A window is a run of consecutive records: a *core*, whose genotypes this solve gives, and
//! *margins* either side read only for context. For each sample the solve chooses its templates
//! over the whole window, fits or reads its walk parameters, walks every site of the window, and
//! keeps the posteriors at the core's records. Samples are solved independently of one another,
//! so they are shared among threads and the result does not depend on how many there are.
//!
//! **A sample whose reads said nothing at a record is called there only with evidence nearby**
//! (§4.2): an informative record of its own within one segment length `S` on each side, or on one
//! side where the contig ends within `S` on the other. Without that its posterior is the cohort's
//! frequency reached through its templates — a confident `0/0` with nothing behind it — and it
//! stays `./.`, with the posterior as its `GP`.

use std::ops::Range;

use rayon::prelude::*;

use super::fit::SampleFit;
use super::model::{WalkSite, walk};
use super::records::{SampleRow, SiteRecord};
use super::templates::choose_templates;

/// The chance a template's allele is not the one the copy carries (`ε`, §3.4) — the prototype's.
pub const MISMATCH: f64 = 0.01;

/// **What the pass did to one sample at one record.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Not re-called: the sample is absent at the record, the record is no site, or the sample has
    /// no informative record in the window. Its column is written as it came.
    Unchanged,
    /// Re-called: written with the argmax of the new posterior.
    Called,
    /// Re-called without evidence of its own nearby (§4.2): written `./.`, the new posterior as
    /// its `GP`.
    Withheld,
}

/// **The pass's result at one core record**, every sample.
#[derive(Clone, Debug, PartialEq)]
pub struct RefinedRecord {
    /// What happened to each sample.
    pub outcomes: Vec<Outcome>,
    /// Each sample's new genotype posterior, `genotype_count` numbers a sample in the genotype
    /// table's order; meaningful only where the outcome is not [`Outcome::Unchanged`].
    pub posteriors: Vec<f64>,
}

/// One sample's posteriors at the core's sites.
struct SampleSolution {
    /// Index into the window's records, the outcome, and where the posterior starts in
    /// `posteriors`.
    sites: Vec<(usize, Outcome, usize)>,
    posteriors: Vec<f64>,
}

/// **Solves one window**: `records` is the whole window, `core` the indices of the records it
/// writes, `contig_length` the contig's length if the header gives it (for §4.2's rule at the
/// contig's end), `fits` one per sample, updated while the samples are still fitting. Returns one
/// [`RefinedRecord`] per core record, in order.
///
/// # Panics
///
/// If `fits` does not hold one entry per sample, or the records' positions descend.
#[must_use]
pub fn refine_window(
    records: &[SiteRecord],
    core: Range<usize>,
    contig_length: Option<u64>,
    fits: &mut [SampleFit],
    max_templates: usize,
) -> Vec<RefinedRecord> {
    assert!(
        records
            .windows(2)
            .all(|pair| pair[0].position <= pair[1].position),
        "a window's records must not descend in position"
    );
    let samples = records.first().map_or(0, SiteRecord::sample_count);
    assert_eq!(fits.len(), samples, "one fit a sample");
    let solutions: Vec<Option<SampleSolution>> = fits
        .par_iter_mut()
        .enumerate()
        .map(|(sample, fit)| {
            refine_sample(
                records,
                core.clone(),
                contig_length,
                sample,
                fit,
                max_templates,
            )
        })
        .collect();

    let mut refined: Vec<RefinedRecord> = records[core.clone()]
        .iter()
        .map(|record| RefinedRecord {
            outcomes: vec![Outcome::Unchanged; samples],
            posteriors: vec![0.0; samples * record.genotype_count()],
        })
        .collect();
    for (sample, solution) in solutions.iter().enumerate() {
        let Some(solution) = solution else {
            continue;
        };
        for &(record, outcome, start) in &solution.sites {
            let genotypes = records[record].genotype_count();
            let out = &mut refined[record - core.start];
            out.outcomes[sample] = outcome;
            out.posteriors[sample * genotypes..(sample + 1) * genotypes]
                .copy_from_slice(&solution.posteriors[start..start + genotypes]);
        }
    }
    refined
}

/// One sample's walk over the window, or `None` if it has no informative record there.
fn refine_sample(
    records: &[SiteRecord],
    core: Range<usize>,
    contig_length: Option<u64>,
    sample: usize,
    fit: &mut SampleFit,
    max_templates: usize,
) -> Option<SampleSolution> {
    let templates = choose_templates(records, sample, max_templates)?;
    let walked: Vec<usize> = (0..records.len())
        .filter(|&index| records[index].is_site())
        .collect();

    // Every site's likelihoods and template shares, the background last, in two flat buffers the
    // walk's sites point into.
    let mut likelihoods: Vec<f64> = Vec::new();
    let mut shares: Vec<f64> = Vec::new();
    let mut spans = Vec::with_capacity(walked.len());
    for &index in &walked {
        let record = &records[index];
        let likelihood_start = likelihoods.len();
        likelihoods.extend(record.likelihoods(sample));
        let share_start = shares.len();
        for &template in &templates {
            shares.extend(record.template_shares(template));
        }
        shares.extend(record.background_shares());
        spans.push((
            likelihood_start..likelihoods.len(),
            share_start..shares.len(),
        ));
    }
    let sites: Vec<WalkSite<'_>> = walked
        .iter()
        .zip(&spans)
        .map(|(&index, (likelihood_span, share_span))| WalkSite {
            position: records[index].position,
            allele_count: records[index].allele_count,
            likelihoods: &likelihoods[likelihood_span.clone()],
            template_shares: &shares[share_span.clone()],
        })
        .collect();
    let walk_templates = templates.len() + 1;

    let informative_positions: Vec<u64> = walked
        .iter()
        .filter(|&&index| records[index].row(sample) == SampleRow::Informative)
        .map(|&index| records[index].position)
        .collect();
    let parameters = fit.fit_window(
        &sites,
        walk_templates,
        informative_positions.len(),
        MISMATCH,
    );
    let result = walk(&sites, walk_templates, parameters);
    let reach = parameters.segment_bases as u64;

    let mut solution = SampleSolution {
        sites: Vec::new(),
        posteriors: Vec::new(),
    };
    for (site, &index) in walked.iter().enumerate() {
        if !core.contains(&index) {
            continue;
        }
        let record = &records[index];
        let outcome = match record.row(sample) {
            SampleRow::Absent => continue,
            SampleRow::Informative => Outcome::Called,
            SampleRow::Flat => {
                if has_evidence_nearby(
                    &informative_positions,
                    record.position,
                    reach,
                    contig_length,
                ) {
                    Outcome::Called
                } else {
                    Outcome::Withheld
                }
            }
        };
        solution
            .sites
            .push((index, outcome, solution.posteriors.len()));
        solution
            .posteriors
            .extend_from_slice(&result.posteriors[site]);
    }
    Some(solution)
}

/// **§4.2's rule for a sample whose reads said nothing at `position`**: an informative record of
/// its own within `reach` bases on each side, or on one side where the contig ends within `reach`
/// on the other. `informative_positions` ascend. A record at the same position counts on neither
/// side: it shares the site's reads, so it does not bracket the site.
#[must_use]
pub fn has_evidence_nearby(
    informative_positions: &[u64],
    position: u64,
    reach: u64,
    contig_length: Option<u64>,
) -> bool {
    let before = informative_positions.partition_point(|&other| other < position);
    let left = before > 0 && position - informative_positions[before - 1] <= reach;
    let after = informative_positions.partition_point(|&other| other <= position);
    let right =
        after < informative_positions.len() && informative_positions[after] - position <= reach;
    // Positions are 1-based: the contig starts `position − 1` bases to the left.
    let start_near = position.saturating_sub(1) < reach;
    let end_near = contig_length.is_some_and(|length| length.saturating_sub(position) < reach);
    (left && right) || (left && end_near) || (right && start_near)
}

#[cfg(test)]
mod tests;
