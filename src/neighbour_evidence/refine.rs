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
//! **Each sample's fills are checked on its own data** ([`FillCheck`]): at core records where its
//! reads settle its genotype, the genotype its neighbours alone predict — what it would be written
//! with had it no reads there — is compared with the reads'. With phasing those are only the
//! records hidden from the rounds, the only ones where the prediction holds none of the sample's own
//! reads (phasing spec §4.1).

use std::ops::Range;

use rayon::prelude::*;

use super::fit::SampleFit;
use super::model::{
    GenotypePairs, HaplotypePair, WalkParameters, WalkResult, WalkSite, walk, walk_and_draw,
};
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
/// test of hiding genotypes the reads settled and asking the pass for them: the walk's prediction
/// without the sample's own reads at a record is what a no-read call there would be written with.
/// Without phasing every informative record serves, at no cost; with phasing only the records
/// hidden from the rounds do ([`HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW`] a window at most), because
/// elsewhere the sample's reads come back to it through the haplotypes it lent other samples
/// (phasing spec §4.1, §7).
///
/// `checked` counts those core records where the sample's own evidence settles its genotype
/// ([`SETTLED_POSTERIOR`]) and the [`NoReadRule`] would write the prediction as a call —
/// evidence on both sides, `GQ` at the floor, the background below its limit; `disagreed` those
/// where the prediction's genotype is not the reads'. `disagreed / checked` estimates the sample's
/// fill error rate where it had reads to spare — a best case, since a site no read reached can be
/// harder. `filled` counts its no-read calls written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FillCheck {
    /// No-read calls written.
    pub filled: u32,
    /// Records where the check ran.
    pub checked: u32,
    /// Of those, where the neighbours' prediction was not the settled genotype.
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

/// **How a window is solved: in one round, or in phasing's rounds** (phasing spec §2.1, §5). One
/// round is the pass without phasing. With phasing the genotypes are averaged over rounds 2 to the
/// last, and the window's draws are keyed by `window_key`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Phasing {
    rounds: usize,
    window_key: u64,
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

    /// The phased pass's [`PHASING_ROUNDS`], for the window whose draws `window_key` keys
    /// ([`window_key`]).
    #[must_use]
    pub fn rounds_for(window_key: u64) -> Self {
        Self::new(PHASING_ROUNDS, window_key)
    }

    /// `rounds` rounds keyed by `window_key`; one round is [`Phasing::none`].
    ///
    /// # Panics
    ///
    /// If `rounds` is 0.
    #[must_use]
    pub fn new(rounds: usize, window_key: u64) -> Self {
        assert!(rounds > 0, "a window is solved in at least one round");
        Self { rounds, window_key }
    }

    fn is_on(self) -> bool {
        self.rounds > 1
    }

    /// Whether `round`'s results enter the written mean: round 1 alone without phasing, rounds 2
    /// on with it — round 1 is the pass without phasing.
    fn averages(self, round: usize) -> bool {
        !self.is_on() || round >= 2
    }

    /// How many rounds the written mean is over.
    fn averaged_rounds(self) -> f64 {
        if self.is_on() {
            (self.rounds - 1) as f64
        } else {
            1.0
        }
    }
}

/// **The key of a window's draws** (phasing spec §5): its contig and the position of its first
/// core site. It depends on the window alone — not on which other contigs or regions the run
/// called — so a region's genotypes are the same whether it is called alone or with others.
#[must_use]
pub fn window_key(contig: u32, position: u64) -> u64 {
    (u64::from(contig) << 40) ^ position
}

/// **The rounds a phased window is solved in** (phasing spec §2.2: a starting value, measured on one
/// cohort; 3 and 8 to be measured).
pub const PHASING_ROUNDS: usize = 5;

/// **How many core records of each sample are hidden from the rounds in each window, at most**, for
/// the fill check (phasing spec §4.1): a run of one window must be able to reach the check's 30
/// sites.
pub const HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW: usize = 50;

/// **The largest share of a sample's informative core records hidden**, so a short window does not
/// lose much of a sample's evidence to the check. Set at the first build with the 50 above; on the
/// 2-Mb simulated cohorts it moved wrong genotypes within the seed spread (phasing spec §4.1).
pub const MAX_HIDDEN_SHARE: f64 = 0.05;

/// **Solves one window**: `records` is the whole window, `core` the indices of the records it
/// writes, `contig_length` the contig's length if the header gives it (for §4.2's rule at the
/// contig's end), `fits` one per sample, updated while the samples are still fitting, `no_read`
/// when a sample without reads is written with a genotype. One round: the pass without phasing,
/// which the run no longer uses — it phases ([`refine_window_in_rounds`]); this is the reference
/// round 1 is checked against. Returns one [`RefinedRecord`] per core record, in order, and each
/// sample's [`FillCheck`] over the core.
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

/// **Where a window's records sit in its walks.** Three index spaces meet here, and every
/// conversion goes through this type: a *record index* into the window's records, a *site* — a
/// walked record's position among the walked ones — and a *core site*, a site in the core.
struct WindowLayout {
    /// The record index of each site, ascending.
    record_of_site: Vec<usize>,
    /// Each record's site, or `None` for a record no walk uses.
    site_of_record: Vec<Option<usize>>,
    /// The core's sites, ascending.
    core_sites: Vec<usize>,
}

impl WindowLayout {
    fn new(records: &[SiteRecord], core: &Range<usize>) -> Self {
        let record_of_site: Vec<usize> = (0..records.len())
            .filter(|&index| records[index].is_site())
            .collect();
        let mut site_of_record = vec![None; records.len()];
        for (site, &index) in record_of_site.iter().enumerate() {
            site_of_record[index] = Some(site);
        }
        let core_sites = (0..record_of_site.len())
            .filter(|&site| core.contains(&record_of_site[site]))
            .collect();
        Self {
            record_of_site,
            site_of_record,
            core_sites,
        }
    }

    fn sites(&self) -> usize {
        self.record_of_site.len()
    }

    /// The site of a record a walk uses.
    fn site(&self, record: usize) -> usize {
        // PANIC-FREE: callers pass only records `templates::compared_records` returned, which are
        // sites, or records of `record_of_site`.
        self.site_of_record[record].expect("a record read as a site is one")
    }
}

/// **Every sample's two haplotypes at every site of a window** (phasing spec §5), double-buffered:
/// a round reads [`HaplotypeStore::current`] and each sample writes its own row of `next`, which
/// replaces `current` at the round's end. A row holds, site by site, copy 0's `A` allele
/// probabilities then copy 1's; at a sample's absent rows it holds `1/A`.
struct HaplotypeStore {
    /// Where each site's two haplotypes start in a row; one past the last at the end.
    offsets: Vec<usize>,
    current: Vec<Vec<f64>>,
    next: Vec<Vec<f64>>,
}

impl HaplotypeStore {
    /// Every sample's `GP` allele shares as both its haplotypes (phasing spec §2.1).
    fn new(records: &[SiteRecord], layout: &WindowLayout, samples: usize) -> Self {
        let mut offsets = Vec::with_capacity(layout.sites() + 1);
        offsets.push(0);
        for &index in &layout.record_of_site {
            let last = offsets.last().copied().unwrap_or(0);
            offsets.push(last + 2 * records[index].allele_count);
        }
        let rows: Vec<Vec<f64>> = (0..samples)
            .map(|sample| {
                let mut row = Vec::with_capacity(offsets[layout.sites()]);
                for &index in &layout.record_of_site {
                    let start = row.len();
                    row.extend(records[index].template_shares(sample));
                    row.extend_from_within(start..);
                }
                row
            })
            .collect();
        Self {
            offsets,
            next: rows.clone(),
            current: rows,
        }
    }

    /// Copy `copy` of `sample`'s haplotype at `site`, from the round being read.
    fn haplotype(&self, sample: usize, copy: usize, site: usize) -> &[f64] {
        let alleles = (self.offsets[site + 1] - self.offsets[site]) / 2;
        let start = self.offsets[site] + copy * alleles;
        &self.current[sample][start..start + alleles]
    }

    /// Writes one sample's drawn haplotypes into its row of the round being written; `absent(site)`
    /// says where the sample has no column, which keeps `1/A`.
    fn write(
        offsets: &[usize],
        row: &mut [f64],
        haplotypes: &[HaplotypePair],
        absent: impl Fn(usize) -> bool,
    ) {
        assert_eq!(haplotypes.len() + 1, offsets.len(), "one drawn pair a site");
        for (site, [first, second]) in haplotypes.iter().enumerate() {
            let alleles = first.len();
            let slot = &mut row[offsets[site]..offsets[site + 1]];
            if absent(site) {
                slot.fill(1.0 / alleles as f64);
            } else {
                slot[..alleles].copy_from_slice(first);
                slot[alleles..].copy_from_slice(second);
            }
        }
    }

    /// The round's end: the rows written become the rows read.
    fn finish_round(&mut self) {
        std::mem::swap(&mut self.current, &mut self.next);
    }
}

/// **One core site's quantities from one round, or their sums over the averaged rounds.**
struct CoreSite {
    posterior: Vec<f64>,
    prediction: Vec<f64>,
    background_copy: f64,
}

impl CoreSite {
    fn add(&mut self, other: CoreSite) {
        let CoreSite {
            posterior,
            prediction,
            background_copy,
        } = other;
        add_to(&mut self.posterior, &posterior);
        add_to(&mut self.prediction, &prediction);
        self.background_copy += background_copy;
    }
}

/// One sample's state across a window's rounds.
struct SampleRounds {
    /// Round 1's template samples; `None` when the sample has no informative record in the window,
    /// so it is not walked and carries its `GP` shares as both haplotypes (phasing spec §2.3).
    templates: Option<Vec<usize>>,
    /// The record indices its templates are chosen over.
    compared: Vec<usize>,
    /// Its previous round's posterior at each compared record, in `compared`'s order.
    previous_posteriors: Vec<Vec<f64>>,
    /// Record indices hidden from the rounds for the fill check, ascending; empty without phasing.
    hidden: Vec<usize>,
    /// Positions of its informative records, for §4.2's evidence-nearby rule.
    informative_positions: Vec<u64>,
    /// The parameters round 1 fitted or read.
    parameters: Option<WalkParameters>,
    /// One per core site, in order: the sums over the averaged rounds; empty until the first.
    sums: Vec<CoreSite>,
}

/// **Solves one window in `phasing`'s rounds** (phasing spec §2.1). Round 1 is the pass without
/// phasing; from round 2 each sample's templates are other samples' haplotypes from the round before
/// (`choose_haplotypes`), every sample walks and draws a path, and its two haplotypes are
/// re-estimated along it; all samples' haplotypes are replaced together at the end of a round. The
/// result is the mean of rounds 2 to the last (round 1's alone without phasing), to which the no-read
/// rules and the fill check are applied once. With phasing, up to
/// [`HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW`] core records of each sample where its reads are
/// informative are hidden from every round, and the fill check runs at them only (§4.1).
///
/// # Panics
///
/// As [`refine_window`].
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
    let samples = records.first().map_or(0, SiteRecord::sample_count);
    assert_eq!(fits.len(), samples, "one fit a sample");
    let layout = WindowLayout::new(records, &core);

    let mut states: Vec<SampleRounds> = (0..samples)
        .into_par_iter()
        .map(|sample| {
            let templates = choose_templates(records, sample, max_templates);
            let informative_positions = layout
                .record_of_site
                .iter()
                .filter(|&&index| records[index].row(sample) == SampleRow::Informative)
                .map(|&index| records[index].position)
                .collect();
            let hidden = if phasing.is_on() && templates.is_some() {
                draw_hidden_records(records, &layout, sample, phasing.window_key)
            } else {
                Vec::new()
            };
            SampleRounds {
                templates,
                compared: compared_records(records, sample),
                previous_posteriors: Vec::new(),
                hidden,
                informative_positions,
                parameters: None,
                sums: Vec::new(),
            }
        })
        .collect();

    let mut store = HaplotypeStore::new(records, &layout, samples);
    for round in 1..=phasing.rounds {
        let mut next = std::mem::take(&mut store.next);
        let store_view = &store;
        // Each sample solves its round and writes its own row of the round being written; it reads
        // only the round before, so the order samples are solved in changes nothing.
        states
            .par_iter_mut()
            .zip(fits.par_iter_mut())
            .zip(next.par_iter_mut())
            .enumerate()
            .for_each(|(sample, ((state, fit), row))| {
                solve_round(
                    records,
                    &layout,
                    store_view,
                    row,
                    sample,
                    state,
                    fit,
                    max_templates,
                    round,
                    phasing,
                );
            });
        store.next = next;
        store.finish_round();
    }

    let averaged = phasing.averaged_rounds();
    let solutions: Vec<Option<SampleSolution>> = states
        .par_iter()
        .enumerate()
        .map(|(sample, state)| {
            final_solution(
                records,
                &layout,
                sample,
                state,
                averaged,
                phasing,
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

/// **Draws the core records of `sample` hidden from the rounds** (phasing spec §4.1): up to
/// [`HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW`] of the core sites where its reads are informative, and
/// at most [`MAX_HIDDEN_SHARE`] of them, drawn without replacement — a partial Fisher–Yates shuffle
/// of them in record order, one uniform a draw from the window's [`StreamPurpose::HiddenRecords`]
/// stream — and returned as record indices, ascending.
fn draw_hidden_records(
    records: &[SiteRecord],
    layout: &WindowLayout,
    sample: usize,
    window_key: u64,
) -> Vec<usize> {
    let mut candidates: Vec<usize> = layout
        .core_sites
        .iter()
        .map(|&site| layout.record_of_site[site])
        .filter(|&index| records[index].row(sample) == SampleRow::Informative)
        .collect();
    let mut uniforms = Uniforms::new(StreamPurpose::HiddenRecords, window_key, sample, 0);
    let count = HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW
        .min((MAX_HIDDEN_SHARE * candidates.len() as f64) as usize)
        .min(candidates.len());
    for draw in 0..count {
        let remaining = candidates.len() - draw;
        // The largest uniform is 1 − 2⁻⁵³, so the product stays below `remaining`; the clamp only
        // states it.
        let picked =
            draw + ((uniforms.next_uniform() * remaining as f64) as usize).min(remaining - 1);
        candidates.swap(draw, picked);
    }
    let mut hidden = candidates[..count].to_vec();
    hidden.sort_unstable();
    hidden
}

/// **One sample's walk in one round**: its posteriors at its compared records and, in averaged
/// rounds, its core sites' quantities go into `state`; with phasing its drawn haplotypes go into
/// `row`, its row of the round being written. Nothing happens for a sample not walked in the
/// window, whose row keeps its `GP` shares.
#[allow(clippy::too_many_arguments)] // the window's data, one sample's state, and the round
fn solve_round(
    records: &[SiteRecord],
    layout: &WindowLayout,
    store: &HaplotypeStore,
    row: &mut [f64],
    sample: usize,
    state: &mut SampleRounds,
    fit: &mut SampleFit,
    max_templates: usize,
    round: usize,
    phasing: Phasing,
) {
    let Some(round_one) = state.templates.as_ref() else {
        return;
    };
    // Round 1: template samples' GP shares; later rounds: other samples' haplotypes.
    let haplotype_templates: Vec<(usize, usize)> = if round == 1 {
        Vec::new()
    } else {
        let haplotype = |other: usize, copy: usize, record: usize| -> &[f64] {
            store.haplotype(other, copy, layout.site(record))
        };
        choose_haplotypes(
            records,
            sample,
            &state.compared,
            &state.previous_posteriors,
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
    let mut spans = Vec::with_capacity(layout.sites());
    for (site, &index) in layout.record_of_site.iter().enumerate() {
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
                shares.extend_from_slice(store.haplotype(other, copy, site));
            }
        }
        shares.extend(record.background_shares());
        spans.push((
            likelihood_start..likelihoods.len(),
            share_start..shares.len(),
        ));
    }
    let sites: Vec<WalkSite<'_>> = layout
        .record_of_site
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
            // Round 1 only: fitting again would add the window's log-likelihoods twice. The hidden
            // records are not informative in the walk, so they do not count towards freezing.
            let parameters = fit.fit_window(
                &sites,
                walk_templates,
                state.informative_positions.len() - state.hidden.len(),
                MISMATCH,
            );
            state.parameters = Some(parameters);
            parameters
        }
    };
    let result = if phasing.is_on() {
        let mut uniforms = Uniforms::new(StreamPurpose::Path, phasing.window_key, sample, round);
        let (result, haplotypes) = walk_and_draw(&sites, walk_templates, parameters, &mut uniforms);
        HaplotypeStore::write(&store.offsets, row, &haplotypes, |site| {
            records[layout.record_of_site[site]].row(sample) == SampleRow::Absent
        });
        result
    } else {
        walk(&sites, walk_templates, parameters)
    };
    let WalkResult {
        mut posteriors,
        mut neighbour_predictions,
        last_template_copy,
        ..
    } = result;
    state.previous_posteriors = state
        .compared
        .iter()
        .map(|&record| posteriors[layout.site(record)].clone())
        .collect();
    if phasing.averages(round) {
        let core = layout.core_sites.iter().map(|&site| CoreSite {
            posterior: std::mem::take(&mut posteriors[site]),
            prediction: std::mem::take(&mut neighbour_predictions[site]),
            background_copy: last_template_copy[site],
        });
        if state.sums.is_empty() {
            state.sums = core.collect();
        } else {
            assert_eq!(
                state.sums.len(),
                layout.core_sites.len(),
                "one sum a core site"
            );
            for (sum, site) in state.sums.iter_mut().zip(core) {
                sum.add(site);
            }
        }
    }
}

/// **One sample's written result from its averaged rounds**: the no-read rules and the fill check
/// (pass §4.2, §4.3; phasing spec §4, §4.1). `None` when the sample was not walked.
#[allow(clippy::too_many_arguments)] // the window's data, one sample's state, and the rules
fn final_solution(
    records: &[SiteRecord],
    layout: &WindowLayout,
    sample: usize,
    state: &SampleRounds,
    averaged: f64,
    phasing: Phasing,
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
    assert_eq!(
        state.sums.len(),
        layout.core_sites.len(),
        "one sum a core site"
    );
    for (&site, sum) in layout.core_sites.iter().zip(&state.sums) {
        let index = layout.record_of_site[site];
        let record = &records[index];
        let mean =
            |values: &[f64]| -> Vec<f64> { values.iter().map(|value| value / averaged).collect() };
        let prediction = mean(&sum.prediction);
        let background = sum.background_copy / averaged;
        let mut posterior = mean(&sum.posterior);
        let outcome = match record.row(sample) {
            SampleRow::Absent => continue,
            SampleRow::Informative => {
                let hidden = state.hidden.binary_search(&index).is_ok();
                // Without phasing every informative record is checked, its prediction free of its
                // own reads; with phasing only hidden records are (phasing spec §4.1).
                if !phasing.is_on() || hidden {
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
                        // reads enter once. The total is positive: the likelihoods' largest is 1
                        // and every prediction weight is, through the mismatch blur.
                        posterior = prediction
                            .iter()
                            .zip(&likelihoods)
                            .map(|(weight, likelihood)| weight * likelihood)
                            .collect();
                        let total: f64 = posterior.iter().sum();
                        debug_assert!(total > 0.0, "a hidden record's posterior has weight");
                        for value in &mut posterior {
                            *value /= total;
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
