//! **Re-calling every sample in one window** (`doc/devel/ng/spec/neighbour_evidence_pass.md` §2.2,
//! §3, §4.2).
//!
//! A window is a run of consecutive records: a *core*, whose genotypes this solve gives, and
//! *margins* either side read only for context. For each sample the solve chooses its templates
//! over the whole window, fits or reads its walk parameters, walks every site of the window, and
//! keeps the posteriors at the core's records. Samples are solved independently of one another,
//! so they are shared among threads and the result does not depend on how many there are.
//!
//! **With phasing the window is solved in rounds** (`neighbour_evidence_phasing.md` §2):
//! round 1 is the solve above; in each later round every sample copies other samples' haplotypes
//! from the round before, draws a path and re-estimates its own two haplotypes along it, and the
//! written result is the mean of rounds 2 to [`PHASING_ROUNDS`] ([`refine_window_in_rounds`]).
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
use super::model::{GenotypePairs, HaplotypePair, WalkParameters, WalkSite, walk, walk_and_draw};
use super::random::{StreamPurpose, Uniforms};
use super::records::{SampleRow, SiteRecord};
use super::templates::{choose_haplotypes, choose_templates, compared_records};
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

/// **How many rounds a window is solved in, and what keys its random draws** (phasing spec §2.1,
/// §5). One round is the pass without phasing; more are phasing, its genotypes averaged over rounds
/// 2 to `rounds`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Phasing {
    /// Rounds; 1 writes round 1, the pass without phasing.
    pub rounds: usize,
    /// The run-wide number of the window's first core site, which keys the window's draws.
    pub window_key: u64,
}

impl Phasing {
    /// No phasing: one round, the pass as it was before phasing.
    #[must_use]
    pub fn none() -> Self {
        Self {
            rounds: 1,
            window_key: 0,
        }
    }

    /// The phased pass's rounds, for the window whose first core site is run-wide site `window_key`.
    #[must_use]
    pub fn rounds_for(window_key: u64) -> Self {
        Self {
            rounds: PHASING_ROUNDS,
            window_key,
        }
    }
}

/// **The rounds a phased window is solved in** (phasing spec §2.2: a starting value, measured on one
/// cohort; 3 and 8 to be measured).
pub const PHASING_ROUNDS: usize = 5;

/// **How many core records of each sample are hidden from the rounds for the fill check**
/// (phasing spec §4.1), at most: a run of one window must be able to reach the check's 30 sites.
pub const HIDDEN_RECORDS_PER_SAMPLE: usize = 50;

/// **The largest share of a sample's informative core records hidden**, so a short window does not
/// lose much of a sample's evidence to the check.
pub const MAX_HIDDEN_SHARE: f64 = 0.05;

/// **Solves one window**: `records` is the whole window, `core` the indices of the records it
/// writes, `contig_length` the contig's length if the header gives it (for §4.2's rule at the
/// contig's end), `fits` one per sample, updated while the samples are still fitting, `no_read`
/// when a sample without reads is written with a genotype. One round: the pass without phasing
/// (see [`refine_window_in_rounds`]). Returns one [`RefinedRecord`] per core record, in order, and
/// each sample's [`FillCheck`] over the core.
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
    refine_window_in_rounds(
        records,
        core,
        contig_length,
        fits,
        max_templates,
        no_read,
        Phasing::none(),
    )
}

/// One sample's state across a window's rounds.
struct SampleRounds {
    /// Round 1's template samples; `None` when the sample has no informative record in the window,
    /// so it is not walked and carries its `GP` shares as both haplotypes (phasing spec §2.3).
    templates: Option<Vec<usize>>,
    /// The records its templates are chosen over.
    compared: Vec<usize>,
    /// Its previous round's posterior at each compared record.
    previous: Vec<Vec<f64>>,
    /// Core records hidden from the rounds for the fill check, ascending; empty with one round.
    hidden: Vec<usize>,
    /// Positions of its informative records, for §4.2's evidence-nearby rule.
    informative_positions: Vec<u64>,
    /// The parameters round 1 fitted or read.
    parameters: Option<WalkParameters>,
    /// Per core walked site, in order: the summed posterior, prediction and background-copy
    /// probability over the rounds averaged.
    sums: Vec<(Vec<f64>, Vec<f64>, f64)>,
}

/// **Solves one window in `phasing.rounds` rounds** (phasing spec §2.1). Round 1 is the pass without
/// phasing; from round 2 each sample's templates are other samples' haplotypes from the round before
/// (`choose_haplotypes`), every sample walks and draws a path, and its two haplotypes are
/// re-estimated along it; all samples' haplotypes are replaced together at the end of a round. The
/// result is the mean of rounds 2 to `rounds` (round 1's alone with one round), to which the no-read
/// rules and the fill check are applied once. With more than one round, [`HIDDEN_RECORDS_PER_SAMPLE`]
/// core records of each sample where its reads are informative are hidden from every round, and the
/// fill check runs at them only (§4.1).
///
/// # Panics
///
/// As [`refine_window`]; and if `phasing.rounds` is 0.
#[must_use]
pub fn refine_window_in_rounds(
    records: &[SiteRecord],
    core: Range<usize>,
    contig_length: Option<u64>,
    fits: &mut [SampleFit],
    max_templates: usize,
    no_read: NoReadRule,
    phasing: Phasing,
) -> WindowRefinement {
    assert!(
        records
            .windows(2)
            .all(|pair| pair[0].position <= pair[1].position),
        "a window's records must not descend in position"
    );
    assert!(
        phasing.rounds > 0,
        "a window is solved in at least one round"
    );
    let samples = records.first().map_or(0, SiteRecord::sample_count);
    assert_eq!(fits.len(), samples, "one fit a sample");
    let rounds = phasing.rounds;

    let walked: Vec<usize> = (0..records.len())
        .filter(|&index| records[index].is_site())
        .collect();
    // Where each walked site's two haplotypes start in a sample's store: `2 · A` numbers a site.
    let mut offsets = Vec::with_capacity(walked.len() + 1);
    offsets.push(0);
    for &index in &walked {
        offsets.push(offsets.last().copied().unwrap_or(0) + 2 * records[index].allele_count);
    }
    let core_walked: Vec<usize> = (0..walked.len())
        .filter(|&site| core.contains(&walked[site]))
        .collect();

    let mut states: Vec<SampleRounds> = (0..samples)
        .into_par_iter()
        .map(|sample| {
            let templates = choose_templates(records, sample, max_templates);
            let informative_positions = walked
                .iter()
                .filter(|&&index| records[index].row(sample) == SampleRow::Informative)
                .map(|&index| records[index].position)
                .collect();
            let hidden = if rounds > 1 && templates.is_some() {
                hidden_records(records, &core_walked, &walked, sample, phasing.window_key)
            } else {
                Vec::new()
            };
            SampleRounds {
                templates,
                compared: compared_records(records, sample),
                previous: Vec::new(),
                hidden,
                informative_positions,
                parameters: None,
                sums: Vec::new(),
            }
        })
        .collect();

    // Every sample's two haplotypes, starting from its `GP` shares for both copies.
    let mut store: Vec<Vec<f64>> = (0..samples)
        .map(|sample| {
            let mut haplotypes = Vec::with_capacity(offsets[walked.len()]);
            for &index in &walked {
                let shares: Vec<f64> = records[index].template_shares(sample).collect();
                haplotypes.extend_from_slice(&shares);
                haplotypes.extend_from_slice(&shares);
            }
            haplotypes
        })
        .collect();

    for round in 1..=rounds {
        let solved: Vec<Option<RoundSolve>> = states
            .par_iter_mut()
            .zip(fits.par_iter_mut())
            .enumerate()
            .map(|(sample, (state, fit))| {
                solve_round(
                    records,
                    &walked,
                    &offsets,
                    &core_walked,
                    &store,
                    sample,
                    state,
                    fit,
                    max_templates,
                    round,
                    phasing,
                )
            })
            .collect();
        for (sample, solve) in solved.into_iter().enumerate() {
            let Some(solve) = solve else {
                continue;
            };
            let state = &mut states[sample];
            state.previous = solve.compared_posteriors;
            if rounds == 1 || round >= 2 {
                if state.sums.is_empty() {
                    state.sums = solve.core;
                } else {
                    for (sum, (posterior, prediction, background)) in
                        state.sums.iter_mut().zip(solve.core)
                    {
                        add_to(&mut sum.0, &posterior);
                        add_to(&mut sum.1, &prediction);
                        sum.2 += background;
                    }
                }
            }
            if let Some(haplotypes) = solve.haplotypes {
                let own = &mut store[sample];
                for (site, pair) in haplotypes.iter().enumerate() {
                    let alleles = records[walked[site]].allele_count;
                    if records[walked[site]].row(sample) == SampleRow::Absent {
                        // The store holds 1/A at absent rows (phasing spec §5).
                        own[offsets[site]..offsets[site + 1]].fill(1.0 / alleles as f64);
                    } else {
                        own[offsets[site]..offsets[site] + alleles].copy_from_slice(&pair[0]);
                        own[offsets[site] + alleles..offsets[site + 1]].copy_from_slice(&pair[1]);
                    }
                }
            }
        }
    }

    let averaged = if rounds == 1 { 1 } else { rounds - 1 } as f64;
    let solutions: Vec<Option<SampleSolution>> = states
        .par_iter()
        .enumerate()
        .map(|(sample, state)| {
            final_solution(
                records,
                &walked,
                &core_walked,
                sample,
                state,
                averaged,
                rounds,
                contig_length,
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

fn add_to(sum: &mut [f64], values: &[f64]) {
    for (total, value) in sum.iter_mut().zip(values) {
        *total += value;
    }
}

/// **The core records of `sample` hidden from the rounds** (phasing spec §4.1): up to
/// [`HIDDEN_RECORDS_PER_SAMPLE`] of the core walked sites where its reads are informative, and at
/// most [`MAX_HIDDEN_SHARE`] of them, drawn
/// without replacement — a partial Fisher–Yates shuffle of them in record order, one uniform a draw
/// — and returned ascending.
fn hidden_records(
    records: &[SiteRecord],
    core_walked: &[usize],
    walked: &[usize],
    sample: usize,
    window_key: u64,
) -> Vec<usize> {
    let mut candidates: Vec<usize> = core_walked
        .iter()
        .map(|&site| walked[site])
        .filter(|&index| records[index].row(sample) == SampleRow::Informative)
        .collect();
    let mut uniforms = Uniforms::new(StreamPurpose::HiddenRecords, window_key, sample, 0);
    let count = HIDDEN_RECORDS_PER_SAMPLE
        .min((MAX_HIDDEN_SHARE * candidates.len() as f64) as usize)
        .min(candidates.len());
    for draw in 0..count {
        let remaining = candidates.len() - draw;
        let picked =
            draw + ((uniforms.next_uniform() * remaining as f64) as usize).min(remaining - 1);
        candidates.swap(draw, picked);
    }
    let mut hidden = candidates[..count].to_vec();
    hidden.sort_unstable();
    hidden
}

/// What one sample's walk in one round gives back.
struct RoundSolve {
    /// Its two haplotypes at every walked site; `None` with one round (no path is drawn).
    haplotypes: Option<Vec<HaplotypePair>>,
    /// Its posterior at each compared record, for the next round's haplotype choice.
    compared_posteriors: Vec<Vec<f64>>,
    /// At each core walked site, in order: the posterior, the prediction and the background-copy
    /// probability.
    core: Vec<(Vec<f64>, Vec<f64>, f64)>,
}

/// **One sample's walk in one round**, or `None` when it is not walked in the window.
#[allow(clippy::too_many_arguments)]
fn solve_round(
    records: &[SiteRecord],
    walked: &[usize],
    offsets: &[usize],
    core_walked: &[usize],
    store: &[Vec<f64>],
    sample: usize,
    state: &mut SampleRounds,
    fit: &mut SampleFit,
    max_templates: usize,
    round: usize,
    phasing: Phasing,
) -> Option<RoundSolve> {
    let round_one = state.templates.as_ref()?;
    let haplotype = |other: usize, copy: usize, index: usize| -> &[f64] {
        let site = walked
            .binary_search(&index)
            .expect("a compared record is walked");
        let alleles = records[index].allele_count;
        let start = offsets[site] + copy * alleles;
        &store[other][start..start + alleles]
    };
    // Round 1: template samples' GP shares; later rounds: other samples' haplotypes.
    let haplotype_templates: Vec<(usize, usize)> = if round == 1 {
        Vec::new()
    } else {
        let previous = |index: usize| -> &[f64] {
            let position = state
                .compared
                .binary_search(&index)
                .expect("the choice reads compared records");
            &state.previous[position]
        };
        choose_haplotypes(
            records,
            sample,
            &state.compared,
            &previous,
            &haplotype,
            max_templates,
        )
    };
    let template_count = if round == 1 {
        round_one.len()
    } else {
        haplotype_templates.len()
    };

    let mut likelihoods: Vec<f64> = Vec::new();
    let mut shares: Vec<f64> = Vec::new();
    let mut spans = Vec::with_capacity(walked.len());
    for &index in walked {
        let record = &records[index];
        let likelihood_start = likelihoods.len();
        if state.hidden.binary_search(&index).is_ok() {
            likelihoods.extend(std::iter::repeat_n(1.0, record.genotype_count()));
        } else {
            likelihoods.extend(record.likelihoods(sample));
        }
        let share_start = shares.len();
        if round == 1 {
            for &template in round_one {
                shares.extend(record.template_shares(template));
            }
        } else {
            for &(other, copy) in &haplotype_templates {
                shares.extend_from_slice(haplotype(other, copy, index));
            }
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
    let walk_templates = template_count + 1;
    let parameters = match state.parameters {
        Some(parameters) => parameters,
        None => {
            // Round 1 only: fitting again would add the window's log-likelihoods twice.
            let parameters = fit.fit_window(
                &sites,
                walk_templates,
                state.informative_positions.len(),
                MISMATCH,
            );
            state.parameters = Some(parameters);
            parameters
        }
    };
    let (result, haplotypes) = if phasing.rounds == 1 {
        (walk(&sites, walk_templates, parameters), None)
    } else {
        let mut uniforms = Uniforms::new(StreamPurpose::Path, phasing.window_key, sample, round);
        let (result, haplotypes) = walk_and_draw(&sites, walk_templates, parameters, &mut uniforms);
        (result, Some(haplotypes))
    };
    let compared_posteriors = state
        .compared
        .iter()
        .map(|index| {
            let site = walked
                .binary_search(index)
                .expect("a compared record is walked");
            result.posteriors[site].clone()
        })
        .collect();
    let core = core_walked
        .iter()
        .map(|&site| {
            (
                result.posteriors[site].clone(),
                result.neighbour_predictions[site].clone(),
                result.last_template_copy[site],
            )
        })
        .collect();
    Some(RoundSolve {
        haplotypes,
        compared_posteriors,
        core,
    })
}

/// **One sample's written result from its averaged rounds**: the no-read rules and the fill check
/// (pass §4.2, §4.3; phasing spec §4, §4.1). `None` when the sample was not walked.
#[allow(clippy::too_many_arguments)]
fn final_solution(
    records: &[SiteRecord],
    walked: &[usize],
    core_walked: &[usize],
    sample: usize,
    state: &SampleRounds,
    averaged: f64,
    rounds: usize,
    contig_length: Option<u64>,
    no_read: NoReadRule,
) -> Option<SampleSolution> {
    let parameters = state.parameters?;
    let reach = parameters.segment_bases as u64;
    let mut pairs = GenotypePairs::default();
    let would_write = |position: u64, posterior: &[f64], background: f64| {
        has_evidence_nearby(&state.informative_positions, position, reach, contig_length)
            && score_best_genotype(posterior).1.get() >= no_read.min_genotype_quality
            && no_read
                .max_background_copy
                .is_none_or(|most| background <= most)
    };
    let mut solution = SampleSolution {
        sites: Vec::new(),
        posteriors: Vec::new(),
        check: FillCheck::default(),
    };
    for (&site, (posterior_sum, prediction_sum, background_sum)) in
        core_walked.iter().zip(&state.sums)
    {
        let index = walked[site];
        let record = &records[index];
        let mean = |sum: &[f64]| -> Vec<f64> { sum.iter().map(|value| value / averaged).collect() };
        let prediction = mean(prediction_sum);
        let background = background_sum / averaged;
        let mut posterior = mean(posterior_sum);
        let outcome = match record.row(sample) {
            SampleRow::Absent => continue,
            SampleRow::Informative => {
                let hidden = state.hidden.binary_search(&index).is_ok();
                // Without phasing every informative record is checked, its prediction free of its
                // own reads; with phasing only hidden records are (phasing spec §4.1).
                if rounds == 1 || hidden {
                    let likelihoods: Vec<f64> = record.likelihoods(sample).collect();
                    let frequencies: Vec<f64> = record.background_shares().collect();
                    if let Some(settled) = settled_genotype(
                        &likelihoods,
                        &frequencies,
                        parameters.autozygous_share,
                        pairs.of(record.allele_count),
                    ) && would_write(record.position, &prediction, background)
                    {
                        solution.check.checked += 1;
                        solution.check.disagreed += u32::from(argmax(&prediction) != settled);
                    }
                    if hidden {
                        // The clean prediction times the sample's own reads, renormalised: its
                        // reads enter once.
                        posterior = prediction
                            .iter()
                            .zip(&likelihoods)
                            .map(|(weight, likelihood)| weight * likelihood)
                            .collect();
                        let total: f64 = posterior.iter().sum();
                        if total > 0.0 {
                            for value in &mut posterior {
                                *value /= total;
                            }
                        }
                    }
                }
                Outcome::Called
            }
            SampleRow::Flat => {
                if would_write(record.position, &posterior, background) {
                    solution.check.filled += 1;
                    Outcome::Called
                } else {
                    Outcome::Withheld
                }
            }
        };
        solution
            .sites
            .push((index, outcome, solution.posteriors.len(), background));
        solution.posteriors.extend_from_slice(&posterior);
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
