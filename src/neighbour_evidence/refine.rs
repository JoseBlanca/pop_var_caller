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
//! side where the contig ends within `S` on the other — and only if the call's `GQ` is at least
//! [`NO_READ_MIN_GENOTYPE_QUALITY`] and it is not probably copying the background
//! ([`NO_READ_MAX_BACKGROUND_COPY`]). Without that its posterior is the cohort's frequency reached
//! through its templates — a confident `0/0` with nothing behind it — and it stays `./.`, with the
//! posterior as its `GP`.
//!
//! **Each sample's fills are checked on its own data** ([`FillCheck`]): at the core's records
//! where its reads settle its genotype, the genotype its neighbours alone predict — what it would
//! be written with had it no reads there — is compared with the reads'.

use std::ops::Range;

use rayon::prelude::*;

use super::fit::SampleFit;
use super::model::{GenotypePairs, WalkSite, walk};
use super::records::{SampleRow, SiteRecord};
use super::templates::choose_templates;
use crate::calling::quality::score_best_genotype;

/// The chance a template's allele is not the one the copy carries (`ε`, §3.4) — the prototype's.
pub const MISMATCH: f64 = 0.01;

/// **The least `GQ` a sample whose reads said nothing is written with a genotype at** (§4.2, owner
/// 2026-11-02): below it, a one-in-ten chance of being wrong, the sample stays `./.` with the pass's
/// posterior as `GP`. On the founder simulator at 1× the calls below it were wrong 28 times in 100,
/// most of them a common alternative called homozygous where the sample carries the reference
/// (`doc/devel/reports/implementations/neighbour_pass_phase4_2026-11-02.md`, decision 1).
pub const NO_READ_MIN_GENOTYPE_QUALITY: f32 = 10.0;

/// **When a sample's own evidence settles its genotype**, for [`FillCheck`]: the genotype's
/// posterior from its reads and the prior calling uses — the cohort's allele frequencies and the
/// sample's inbreeding — is at least this, `GQ` 30. The prior has to be in it: at 3× the reads alone
/// never tell a homozygote from a heterozygote with confidence (ten reads of one allele are needed
/// for a ratio of a thousand), and a check that needs them never runs on the shallow samples, the
/// ones that are filled.
pub const SETTLED_POSTERIOR: f64 = 0.999;

/// **The most probable a sample without reads may be copying the background template** — copying no
/// sampled haplotype (§3.3) — and still be written with a genotype (§4.2, owner 2026-11-02). Where
/// no other sample carries its haplotype the call is the cohort's allele frequency, not a copy. On
/// the founder simulator's inbred cohorts, the fills with the probability between 0.1 and 0.5 were
/// wrong 13 times in 165, against 1 to 2 in 100 for the fills kept; at 10 samples and 1×, the limit
/// keeps 101 of 294 fills with 2 wrong, from 47
/// (`doc/devel/reports/implementations/neighbour_pass_phase4_2026-11-02.md`).
pub const NO_READ_MAX_BACKGROUND_COPY: f64 = 0.1;

/// **When a sample without reads at a record is written with a genotype** (§4.2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NoReadRule {
    /// The least `GQ` the pass's call must reach.
    pub min_genotype_quality: f32,
    /// The most probable the sample may be copying the background template at the record.
    /// `None` applies no such limit.
    pub max_background_copy: Option<f64>,
}

impl Default for NoReadRule {
    fn default() -> Self {
        Self {
            min_genotype_quality: NO_READ_MIN_GENOTYPE_QUALITY,
            max_background_copy: Some(NO_READ_MAX_BACKGROUND_COPY),
        }
    }
}

/// **One sample's fills, and how often a fill would be wrong, measured on its own data** — the
/// test of hiding genotypes the reads settled and asking the pass for them, done at every such
/// record at no cost: the walk's prediction without the sample's own reads at a record is what a
/// no-read call there would be written with.
///
/// `checked` counts the core's records where the sample's own evidence settles its genotype
/// ([`SETTLED_POSTERIOR`]) and the [`NoReadRule`] would write the prediction as a call —
/// evidence on both sides, `GQ` at the floor, the background below its limit; `disagreed` those
/// where the prediction's genotype is not the reads'. `disagreed / checked` estimates the sample's
/// fill error rate where it had reads to spare — a best case, since a site no read reached can be
/// harder. `filled` counts its no-read calls written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FillCheck {
    pub filled: u32,
    pub checked: u32,
    pub disagreed: u32,
}

impl FillCheck {
    /// **The least the sample's fill error rate plausibly is**: the mean of the rate's
    /// Beta(`disagreed + 1`, `checked − disagreed + 1`) posterior, from a flat prior, less two of
    /// its standard deviations. Above a limit, the excess is not the chance of a few checked sites:
    /// one wrong in 30 is a share of 3.3% and a bound of 0.
    #[must_use]
    pub fn wrong_share_at_least(&self) -> f64 {
        let checked = f64::from(self.checked);
        let mean = (f64::from(self.disagreed) + 1.0) / (checked + 2.0);
        mean - 2.0 * (mean * (1.0 - mean) / (checked + 3.0)).sqrt()
    }
}

impl std::ops::AddAssign for FillCheck {
    fn add_assign(&mut self, other: Self) {
        self.filled += other.filled;
        self.checked += other.checked;
        self.disagreed += other.disagreed;
    }
}

/// **The pass's result over one window's core**: each record's, and each sample's fill check.
#[derive(Clone, Debug, PartialEq)]
pub struct WindowRefinement {
    /// One per core record, in order.
    pub records: Vec<RefinedRecord>,
    /// One per sample, over the core; zero for a sample not re-called in the window.
    pub checks: Vec<FillCheck>,
}

/// **What the pass did to one sample at one record.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Not re-called: the sample is absent at the record, the record is no site, or the sample has
    /// no informative record in the window. Its column is written as it came.
    Unchanged,
    /// Re-called: written with the argmax of the new posterior.
    Called,
    /// Re-called without reads, and without evidence of its own nearby, below
    /// [`NO_READ_MIN_GENOTYPE_QUALITY`] or above [`NO_READ_MAX_BACKGROUND_COPY`] (§4.2): written
    /// `./.`, the new posterior as its `GP`.
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
    /// Each sample's probability of copying the background template here — of having no sampled
    /// haplotype to copy (§3.3); meaningful only where the outcome is not [`Outcome::Unchanged`].
    pub background_copy: Vec<f64>,
}

/// One sample's posteriors at the core's sites.
struct SampleSolution {
    /// Index into the window's records, the outcome, where the posterior starts in `posteriors`,
    /// and the probability of copying the background there.
    sites: Vec<(usize, Outcome, usize, f64)>,
    posteriors: Vec<f64>,
    check: FillCheck,
}

/// **Solves one window**: `records` is the whole window, `core` the indices of the records it
/// writes, `contig_length` the contig's length if the header gives it (for §4.2's rule at the
/// contig's end), `fits` one per sample, updated while the samples are still fitting, `no_read`
/// when a sample without reads is written with a genotype. Returns one [`RefinedRecord`] per core
/// record, in order, and each sample's [`FillCheck`] over the core.
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
    no_read: NoReadRule,
) -> WindowRefinement {
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
                no_read,
            )
        })
        .collect();

    let mut refined: Vec<RefinedRecord> = records[core.clone()]
        .iter()
        .map(|record| RefinedRecord {
            outcomes: vec![Outcome::Unchanged; samples],
            posteriors: vec![0.0; samples * record.genotype_count()],
            background_copy: vec![0.0; samples],
        })
        .collect();
    let mut checks = vec![FillCheck::default(); samples];
    for (sample, solution) in solutions.iter().enumerate() {
        let Some(solution) = solution else {
            continue;
        };
        checks[sample] = solution.check;
        for &(record, outcome, start, background) in &solution.sites {
            let genotypes = records[record].genotype_count();
            let out = &mut refined[record - core.start];
            out.outcomes[sample] = outcome;
            out.background_copy[sample] = background;
            out.posteriors[sample * genotypes..(sample + 1) * genotypes]
                .copy_from_slice(&solution.posteriors[start..start + genotypes]);
        }
    }
    WindowRefinement {
        records: refined,
        checks,
    }
}

/// One sample's walk over the window, or `None` if it has no informative record there.
fn refine_sample(
    records: &[SiteRecord],
    core: Range<usize>,
    contig_length: Option<u64>,
    sample: usize,
    fit: &mut SampleFit,
    max_templates: usize,
    no_read: NoReadRule,
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
    let mut pairs = GenotypePairs::default();
    let reach = parameters.segment_bases as u64;

    // Whether a no-read call at `site`, with genotype posterior `posterior`, is written. For a
    // sample with reads there, `posterior` is the neighbours' prediction and the rule asks whether
    // it would have been; the record's own position counts on neither side of it.
    let would_write = |site: usize, posterior: &[f64]| {
        has_evidence_nearby(
            &informative_positions,
            sites[site].position,
            reach,
            contig_length,
        ) && score_best_genotype(posterior).1.get() >= no_read.min_genotype_quality
            && no_read
                .max_background_copy
                .is_none_or(|most| result.last_template_copy[site] <= most)
    };

    let mut solution = SampleSolution {
        sites: Vec::new(),
        posteriors: Vec::new(),
        check: FillCheck::default(),
    };
    for (site, &index) in walked.iter().enumerate() {
        if !core.contains(&index) {
            continue;
        }
        let record = &records[index];
        let outcome = match record.row(sample) {
            SampleRow::Absent => continue,
            SampleRow::Informative => {
                let prediction = &result.neighbour_predictions[site];
                let walk_site = &sites[site];
                let frequencies =
                    &walk_site.template_shares[templates.len() * walk_site.allele_count..];
                if let Some(settled) = settled_genotype(
                    walk_site.likelihoods,
                    frequencies,
                    parameters.autozygous_share,
                    pairs.of(walk_site.allele_count),
                ) && would_write(site, prediction)
                {
                    solution.check.checked += 1;
                    solution.check.disagreed += u32::from(argmax(prediction) != settled);
                }
                Outcome::Called
            }
            SampleRow::Flat => {
                if would_write(site, &result.posteriors[site]) {
                    solution.check.filled += 1;
                    Outcome::Called
                } else {
                    Outcome::Withheld
                }
            }
        };
        solution.sites.push((
            index,
            outcome,
            solution.posteriors.len(),
            result.last_template_copy[site],
        ));
        solution
            .posteriors
            .extend_from_slice(&result.posteriors[site]);
    }
    Some(solution)
}

/// The genotype a sample's own evidence settles ([`SETTLED_POSTERIOR`]), or `None`: its reads'
/// `likelihoods` times the prior of the allele `frequencies` at inbreeding `inbreeding` — a
/// homozygote `F·p + (1 − F)·p²`, a heterozygote `2·(1 − F)·p·q` — with `pairs` each genotype's
/// two alleles.
fn settled_genotype(
    likelihoods: &[f64],
    frequencies: &[f64],
    inbreeding: f64,
    pairs: &[(usize, usize)],
) -> Option<usize> {
    let joint: Vec<f64> = pairs
        .iter()
        .zip(likelihoods)
        .map(|(&(first, second), &likelihood)| {
            let prior = if first == second {
                inbreeding * frequencies[first]
                    + (1.0 - inbreeding) * frequencies[first] * frequencies[first]
            } else {
                2.0 * (1.0 - inbreeding) * frequencies[first] * frequencies[second]
            };
            prior * likelihood
        })
        .collect();
    let total: f64 = joint.iter().sum();
    let best = argmax(&joint);
    (total > 0.0 && joint[best] >= SETTLED_POSTERIOR * total).then_some(best)
}

/// The index of the largest value, the first among equals.
fn argmax(values: &[f64]) -> usize {
    (0..values.len()).fold(0, |best, index| {
        if values[index] > values[best] {
            index
        } else {
            best
        }
    })
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
