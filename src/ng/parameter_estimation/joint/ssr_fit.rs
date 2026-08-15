//! The repeat-tract half of the joint fit: how often a polymerase slips, per stratum.
//!
//! While a polymerase copies a repeat tract it sometimes **slips**, adding or dropping a whole
//! repeat unit, so a read reports a tract one unit longer or shorter than the DNA it came
//! from. This module estimates how often that happens, from every sample's reads at the same
//! tracts at once. It is not the variant caller: it runs once over the cohort, before anything
//! is called, and hands the caller the numbers it will assume.
//!
//! Design: [`spec/parameter_prepass_joint_fit.md`] §4, §4.1 and §4.2. The
//! [`generic`](super::fit) half runs first and this one depends on it in **one direction
//! only** — it takes each sample's homozygote excess and gives nothing back — so a run may
//! drop the ordinary-position records before reading a single tract.
//!
//! # Vocabulary
//!
//! - **Tract** — one repeat region of the reference; what the records call an STR locus.
//! - **Stratum** — every tract sharing a motif length and a reference repeat count. Slippage
//!   depends on repeat count more than on anything else, so a stratum is the unit that is
//!   fitted. Tomato holds 462,701 kept tracts in 141 strata.
//! - **Offset** — a read's tract length minus the reference tract's, in whole repeat units.
//!   The records store `-4 … +4` with the ends saturating
//!   ([`RECORDED_OFFSET_RANGE`](super::census::RECORDED_OFFSET_RANGE)).
//! - **Length spectrum** — how the stratum's chromosomes are spread over the tract lengths.
//! - **Concentration** — how monomorphic the stratum's tracts are. Small means most tracts are
//!   fixed at one length while the stratum as a whole spans many.
//!
//! # What is fitted
//!
//! Per (read group × stratum), three numbers describing slippage: **how often** a read slips,
//! **which way** — at tomato dinucleotides a slipped read shows a shorter tract 4.9 times as
//! often as a longer one — and **how fast** two-unit slips fall off against one-unit slips.
//! Per stratum, the length spectrum and the concentration. **None of them is per tract.**
//!
//! # Two things the design settles, and this module obeys
//!
//! - **A tract's own length frequencies are a latent vector**, drawn from a fitted per-stratum
//!   Dirichlet and integrated away — never a parameter of the tract. Fitting them per tract
//!   directly moves the slippage level 333-fold depending only on where the search starts
//!   (spec §4.1).
//! - **The integral is a fixed 256-point numerical integration** over that Dirichlet. The
//!   earlier design enumerated the cases where a tract is fixed at one length or segregates
//!   exactly two; it cannot represent a tract carrying three, and the fitted slippage absorbs
//!   the difference — +23.7% where 18% of tracts carry three or more, +722% where nearly all
//!   do. It is withdrawn (spec §4.2).
//!
//! # Where the estimator came from
//!
//! It is lifted from `examples/ng_joint_str_harness.rs`, the program the design was measured
//! with, rather than written again — a second implementation is two things to keep agreeing.
//! The harness's `library` mode fits the same draw both ways and prints the two side by side.
//!
//! Three things are generalised in the lift, and each is the specification catching up with
//! the harness rather than a new idea:
//!
//! 1. **Alleles reach further than the read buckets.** The records store `±4`; the lengths the
//!    fit may place allele mass on reach `±6` (`parameter_prepass_joint_records.md` §3.2),
//!    which is what lets an end bucket be attributed to a distant allele rather than to a far
//!    slip. With the two spans equal the arithmetic is the harness's exactly.
//! 2. **The homozygote excess is per sample**, as it arrives from the ordinary-position half,
//!    where the harness supplied one number for the whole panel.
//! 3. **Slippage is per read group**, as spec §4 has it, where the harness had one set of
//!    slippage numbers. Read groups are named in **slippage groups** so a run may pool them;
//!    one group per read group is the specified default and the widest one.
//!
//! [`spec/parameter_prepass_joint_fit.md`]: ../../../../doc/devel/ng/spec/parameter_prepass_joint_fit.md

use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::ng::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, RECORDED_OFFSET_RANGE, SsrEvidence, SsrLocusState,
};
use crate::ng::parameter_estimation::joint::loci::CensusLoci;
use crate::ng::types::{ContigId, ReadGroupId};

// ---------------------------------------------------------------------
// What a read does
// ---------------------------------------------------------------------

/// The three numbers describing how a polymerase slips, for one read group in one stratum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slippage {
    /// How often a read reports a tract length other than its allele's.
    pub level: f64,
    /// Of the reads that slip, the share showing a **shorter** tract. Tomato dinucleotides sit
    /// near 0.83 — 2,438 shorter against 501 longer.
    pub shorter_share: f64,
    /// How fast two-repeat slips fall off against one-repeat slips.
    pub fall_off: f64,
}

impl Slippage {
    /// `P(a read reports each bucket | the allele sits at `allele_offset`)`.
    ///
    /// Both the allele and the buckets are counted in whole repeat units from the reference
    /// tract length. The buckets run `-read_span … +read_span`; the allele may sit outside
    /// them, which is the point of the two spans being different.
    ///
    /// **The end buckets get their marginal**, never the probability of sitting exactly on the
    /// edge: the outermost bucket is *at least this many repeats short*, so every step that
    /// would land at or beyond it — including the tail of the geometric fall-off, and the
    /// unslipped read of an allele that is itself outside the recorded range — is summed into
    /// it. Measured at a recorded range of `±1` on a stratum whose alleles reach three repeats
    /// either side, the marginal rule returns the slippage level to within 0.05% where
    /// plugging in the edge costs 33% of it (`parameter_prepass_joint_records.md` §3.2).
    pub fn read_probabilities(&self, allele_offset: i32, read_span: i32) -> Vec<f64> {
        let buckets = (2 * read_span + 1) as usize;
        let mut out = vec![0.0; buckets];
        let slot = |offset: i32| (offset.clamp(-read_span, read_span) + read_span) as usize;

        out[slot(allele_offset)] += 1.0 - self.level;

        // Shorter, then longer. Each direction gives its own weight to every step that lands
        // inside the buckets, and the whole remaining tail to the end bucket it saturates
        // into. `steps` is how many steps are still inside; a step past that cannot be told
        // from any later one.
        for (direction, share) in [(-1_i32, self.shorter_share), (1, 1.0 - self.shorter_share)] {
            let inside_steps = if direction < 0 {
                (allele_offset + read_span).max(0)
            } else {
                (read_span - allele_offset).max(0)
            };
            for step in 1..=inside_steps {
                let weight = (1.0 - self.fall_off) * self.fall_off.powi(step - 1);
                out[slot(allele_offset + direction * step)] += self.level * share * weight;
            }
            // Everything at least `inside_steps + 1` steps away, whose weight telescopes to
            // `fall_off^inside_steps`, lands in the end bucket.
            let tail = self.fall_off.powi(inside_steps);
            out[slot(allele_offset + direction * (inside_steps + 1))] += self.level * share * tail;
        }
        out
    }
}

// ---------------------------------------------------------------------
// The evidence one stratum brings
// ---------------------------------------------------------------------

/// Which stratum a tract is in: its motif length and the reference's repeat count.
///
/// **It lives in the census** ([`census::Stratum`](super::census::Stratum)), because the
/// census is keyed by it — two of a tract record's counts are held per stratum rather than per
/// locus — and re-exported here, where the fit reads it, so a consumer of either module names
/// one type.
pub use super::census::Stratum;

/// One sample's spanning reads at one tract, split by the slippage group that produced them.
///
/// **Only the samples that put a read on the tract are here.** A sample with no read
/// contributes a likelihood of exactly one whatever the parameters are, so leaving it out is
/// not an approximation — and at three reads a site it is most of the panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleTractReads {
    pub sample: u32,
    /// Per slippage group that put a read here, its reads in each bucket.
    pub by_group: Vec<(u32, Vec<u32>)>,
}

/// Every sample's reads at one tract.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TractReads {
    pub samples: Vec<SampleTractReads>,
}

/// One stratum's tracts, and what the fit needs to know about the shape of them.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumEvidence {
    pub stratum: Stratum,
    pub tracts: Vec<TractReads>,
    /// How many buckets a read's offset was recorded in, `2 × span + 1`.
    pub read_span: i32,
    /// How many slippage groups the cohort declares. Every group is fitted whether or not it
    /// put a read in this stratum; one that did not is returned as not fitted.
    pub groups: usize,
    /// Tracts left out because more than one read in ten of those differing from the reference
    /// length differed by a **non-whole** number of repeat units — the guard's threshold. Such
    /// a tract is not something this noise model describes, and the fit says so rather than
    /// fitting it (`parameter_prepass_joint_records.md` §3.3).
    pub tracts_over_guard_threshold: u64,
    /// Reads that reached a tract and crossed no whole copy of it, so reported no length.
    /// **Not evidence about slippage and not dropped silently**: a tract longer than a read is
    /// never crossed in any sample at any depth, so this count runs along the repeat-count
    /// axis and a stratum unreadable at this read length must not look like one that was
    /// merely unlucky with coverage.
    pub reads_reaching_not_crossing: u64,
    /// Reads whose tract differed from the reference by a non-whole number of repeat units, at
    /// tracts that stayed in. A diagnostic; nothing about slippage is estimated from them.
    pub guard_reads: u64,
}

impl StratumEvidence {
    /// Tracts carrying at least one spanning read — the count the per-stratum floor is
    /// measured in.
    pub fn tracts_with_reads(&self) -> usize {
        self.tracts.iter().filter(|t| !t.samples.is_empty()).count()
    }

    /// Spanning reads, over every tract and sample.
    pub fn spanning_reads(&self) -> u64 {
        self.tracts
            .iter()
            .flat_map(|tract| tract.samples.iter())
            .flat_map(|sample| sample.by_group.iter())
            .flat_map(|(_, counts)| counts.iter())
            .map(|reads| u64::from(*reads))
            .sum()
    }

    /// Which slippage groups put a read in this stratum.
    fn groups_with_reads(&self) -> Vec<bool> {
        let mut seen = vec![false; self.groups];
        for tract in &self.tracts {
            for sample in &tract.samples {
                for (group, counts) in &sample.by_group {
                    if counts.iter().any(|reads| *reads > 0) {
                        seen[*group as usize] = true;
                    }
                }
            }
        }
        seen
    }

    /// The two strata's tracts in one, for borrowing. **The stratum key of the result is the
    /// receiving one's**; what the pooled set is made of travels in [`StratumFit::borrowed`].
    fn pooled_with(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.tracts.extend(other.tracts.iter().cloned());
        out.tracts_over_guard_threshold += other.tracts_over_guard_threshold;
        out.reads_reaching_not_crossing += other.reads_reaching_not_crossing;
        out.guard_reads += other.guard_reads;
        out
    }
}

// ---------------------------------------------------------------------
// What comes back
// ---------------------------------------------------------------------

/// Everything one stratum's fit produces.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumFit {
    pub stratum: Stratum,
    /// Per slippage group, its three slippage numbers — `None` where that group put no read in
    /// this stratum. **An absent group is not a fitted zero**: a group with no reads here has
    /// no slippage estimate, and saying so is the difference between missing and quiet.
    pub slippage: Vec<Option<Slippage>>,
    /// How the stratum's chromosomes are spread over the allele lengths, indexed from
    /// `-allele_span` to `+allele_span` in whole repeat units.
    pub length_spectrum: Vec<f64>,
    /// How monomorphic the stratum's tracts are. Small means most tracts carry one length.
    pub concentration: f64,
    /// The mean log-likelihood a tract, at the returned parameters.
    pub log_likelihood_a_tract: f64,
    /// Tracts the fit actually read — its own if it stood alone, its own plus its neighbours'
    /// if it borrowed.
    pub tracts_fitted: usize,
    /// The neighbouring repeat counts this stratum borrowed tracts from, empty when it stood
    /// on its own. Same period throughout: slippage is not comparable across motif lengths.
    pub borrowed: Vec<u64>,
    /// Whether the climb settled or ran out of rounds. **Running out is never reported as
    /// convergence.**
    pub converged: bool,
}

/// Why a stratum produced no fit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StratumRefusal {
    /// Not one tract in the stratum carried a spanning read.
    NoSpanningReads,
    /// The stratum, with everything it could borrow, still holds fewer tracts than the floor
    /// says can carry an answer.
    BelowTheFloor { tracts: usize, floor: usize },
}

/// One stratum's answer: what was fitted, or why nothing was.
#[derive(Debug, Clone, PartialEq)]
pub enum StratumOutcome {
    Fitted(Box<StratumFit>),
    Refused {
        stratum: Stratum,
        tracts: usize,
        reason: StratumRefusal,
    },
}

// ---------------------------------------------------------------------
// What the run was asked for
// ---------------------------------------------------------------------

/// How many points the integral over a tract's length frequencies is taken on.
///
/// **Fixed, whatever the number of length classes.** A tensor grid needs
/// `nodes^(classes − 1)` points — 576 at three classes and 1.1 × 10¹¹ at thirteen — where this
/// needs the same 256 at every class count. Measured against the grid at three classes it
/// returns the same answers to within 0.3 percentage points on the concentration and 0.2 on
/// every slippage number, in half the time (spec §4.2).
pub const QUADRATURE_POINTS: usize = 256;

/// How far either side of the reference length the fit may place allele mass.
///
/// **Wider than the recorded offsets' `±4`**, and that is what lets an end bucket be
/// attributed to a distant allele rather than to a far slip
/// (`parameter_prepass_joint_records.md` §3.2).
pub const ALLELE_SPAN: i32 = 6;

/// Where one run of the climb starts.
///
/// **The starting points are part of the estimator, not a tuning detail.** They are spread
/// over the slippage level and over how monomorphic tracts are assumed to be, which are the
/// two the objective is least likely to have a single peak in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StartingPoint {
    pub slippage_level: f64,
    pub concentration: f64,
}

impl StartingPoint {
    /// The three the harness climbs from, and the ones every measurement on this path was made
    /// with.
    pub fn spanning_the_monomorphic_range() -> Vec<Self> {
        vec![
            Self {
                slippage_level: 0.02,
                concentration: 0.3,
            },
            Self {
                slippage_level: 0.10,
                concentration: 3.0,
            },
            Self {
                slippage_level: 0.30,
                concentration: 30.0,
            },
        ]
    }
}

#[derive(Debug, Clone)]
pub struct SsrFitConfig {
    /// How far either side of the reference length allele mass may sit.
    pub allele_span: i32,
    pub quadrature_points: usize,
    pub starting_points: Vec<StartingPoint>,
    /// How many rounds of coordinate ascent one start gets.
    pub max_rounds: u32,
    /// A round stops the climb when it improved the mean log-likelihood a tract by less than
    /// this.
    pub stillness: f64,
    /// How many tracts a stratum needs before it is fitted on its own. Below it the stratum
    /// borrows from its neighbouring repeat counts; see [`fit_strata`].
    pub borrowing_floor: usize,
    /// How many tracts a stratum must reach, its borrowings included, before anything is
    /// fitted at all.
    pub refusal_floor: usize,
}

impl Default for SsrFitConfig {
    fn default() -> Self {
        Self {
            allele_span: ALLELE_SPAN,
            quadrature_points: QUADRATURE_POINTS,
            starting_points: StartingPoint::spanning_the_monomorphic_range(),
            max_rounds: 5,
            stillness: 1e-6,
            borrowing_floor: DEFAULT_BORROWING_FLOOR,
            refusal_floor: DEFAULT_REFUSAL_FLOOR,
        }
    }
}

/// How few tracts a stratum may hold and still be fitted on its own.
///
/// **1,000, and it is a decision about which number breaks first rather than about accuracy in
/// general.** Measured on drawn strata at three reads a site
/// (`doc/devel/ng/reports/str_stratum_size_sweep_2026-08-13.md`), between 1,000 tracts and 250
/// the fall-off's spread between draws goes from 5.5% to 13.2% and the concentration's from
/// 2.8% to 7.9%, while the slippage level's only goes from 2.5% to 3.9%. Above 1,000 every one
/// of the five fitted numbers is within a few percent both on average and between draws.
///
/// **Why not the 5,000-tract floor the same sweep names.** That floor is where a stratum
/// carries its own numbers to within a couple of percent; borrowing is not free, because
/// neighbouring repeat counts have genuinely different slippage, so the choice is between a
/// stratum's own noisy answer and its neighbours' steadier but displaced one. Where those two
/// costs cross is measured in `doc/devel/ng/reports/str_fit_on_real_records_2026-08-13.md`.
pub const DEFAULT_BORROWING_FLOOR: usize = 1_000;

/// How few tracts leave nothing worth returning even after borrowing.
///
/// **50, and it is deliberately far below the borrowing floor.** A stratum this thin is one a
/// consumer must be told has no answer rather than handed a number whose spread between draws
/// is a third of its own size — at 50 tracts and three reads a site the twelve draws put the
/// fall-off anywhere from 0.139 to 0.338 against a truth of 0.250.
pub const DEFAULT_REFUSAL_FLOOR: usize = 50;

// ---------------------------------------------------------------------
// Fitting one stratum
// ---------------------------------------------------------------------

/// What the climb carries between evaluations.
#[derive(Debug, Clone, PartialEq)]
struct Parameters {
    /// Per slippage group.
    slippage: Vec<Slippage>,
    /// Over the allele classes, `-allele_span … +allele_span`.
    length_spectrum: Vec<f64>,
    concentration: f64,
}

impl Parameters {
    fn start(start: StartingPoint, groups: usize, classes: usize) -> Self {
        Self {
            slippage: vec![
                Slippage {
                    level: start.slippage_level,
                    shorter_share: 0.5,
                    fall_off: 0.3,
                };
                groups
            ],
            length_spectrum: vec![1.0 / classes as f64; classes],
            concentration: start.concentration,
        }
    }
}

/// Fit one stratum, on its own tracts alone.
///
/// `homozygote_excess` is one number a sample, indexed the way [`SampleTractReads::sample`] is:
/// how far short of the heterozygote proportions the population's allele frequencies predict
/// that sample falls. It arrives from the ordinary-position half and is **held, never fitted
/// here**.
///
/// # Panics
///
/// When a tract names a sample the homozygote-excess list does not reach — the two are one
/// cohort described twice, and a mismatch is a wiring error rather than data.
pub fn fit_stratum(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Option<StratumFit> {
    fit_pooled(evidence, &[], homozygote_excess, config)
}

/// Fit `evidence`, whose tracts may already include borrowed ones, recording where from.
fn fit_pooled(
    evidence: &StratumEvidence,
    borrowed: &[u64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Option<StratumFit> {
    let classes = (2 * config.allele_span + 1) as usize;
    let genotypes = genotype_pairs(classes);
    let live_groups = evidence.groups_with_reads();
    if !live_groups.iter().any(|live| *live) {
        return None;
    }

    let mut best: Option<(Parameters, f64, bool)> = None;
    for start in &config.starting_points {
        let mut parameters = Parameters::start(*start, evidence.groups, classes);
        let mut scorer = Scorer::new(evidence, homozygote_excess, &genotypes, config);
        let mut score = scorer.score(&parameters);
        let mut converged = false;
        for _ in 0..config.max_rounds {
            let before = score;
            climb_one_round(&mut parameters, &mut scorer, &live_groups, classes);
            score = scorer.score(&parameters);
            if score - before < config.stillness {
                converged = true;
                break;
            }
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_score, _)| score > *best_score)
        {
            best = Some((parameters, score, converged));
        }
    }

    let (parameters, score, converged) = best.expect("at least one starting point");
    Some(StratumFit {
        stratum: evidence.stratum,
        slippage: live_groups
            .iter()
            .enumerate()
            .map(|(group, live)| live.then_some(parameters.slippage[group]))
            .collect(),
        length_spectrum: parameters.length_spectrum,
        concentration: parameters.concentration,
        log_likelihood_a_tract: score,
        tracts_fitted: evidence.tracts_with_reads(),
        borrowed: borrowed.to_vec(),
        converged,
    })
}

/// One pass of coordinate ascent over everything the stratum fits.
///
/// **The slippage numbers move first and the frequencies after**, because moving slippage is
/// what invalidates the cached read likelihoods; doing it the other way round would rebuild
/// them on every spectrum coordinate.
fn climb_one_round(
    parameters: &mut Parameters,
    scorer: &mut Scorer<'_>,
    live_groups: &[bool],
    classes: usize,
) {
    for (group, live) in live_groups.iter().enumerate() {
        if !live {
            continue;
        }
        for which in 0..3 {
            let current = read_slippage(&parameters.slippage[group], which);
            let moved = climb_scalar(
                |x| {
                    let mut trial = parameters.clone();
                    write_slippage(&mut trial.slippage[group], which, expit(x));
                    scorer.score(&trial)
                },
                logit(current),
                3.0,
            );
            write_slippage(&mut parameters.slippage[group], which, expit(moved));
        }
    }

    // The spectrum, one class at a time on a log scale, renormalised each time.
    for class in 0..classes {
        let current = parameters.length_spectrum[class].max(1e-9).ln();
        let moved = climb_scalar(
            |x| {
                let mut trial = parameters.clone();
                trial.length_spectrum[class] = x.exp();
                normalise(&mut trial.length_spectrum);
                scorer.score(&trial)
            },
            current,
            2.0,
        );
        parameters.length_spectrum[class] = moved.exp();
        normalise(&mut parameters.length_spectrum);
    }

    let moved = climb_scalar(
        |x| {
            let mut trial = parameters.clone();
            trial.concentration = x.exp();
            scorer.score(&trial)
        },
        parameters.concentration.ln(),
        2.5,
    );
    parameters.concentration = moved.exp();
}

fn read_slippage(slippage: &Slippage, which: usize) -> f64 {
    match which {
        0 => slippage.level,
        1 => slippage.shorter_share,
        _ => slippage.fall_off,
    }
}

fn write_slippage(slippage: &mut Slippage, which: usize, value: f64) {
    match which {
        0 => slippage.level = value,
        1 => slippage.shorter_share = value,
        _ => slippage.fall_off = value,
    }
}

fn normalise(weights: &mut [f64]) {
    let total: f64 = weights.iter().sum();
    for weight in weights.iter_mut() {
        *weight /= total;
    }
}

// ---------------------------------------------------------------------
// Thin strata: borrowing along the repeat-count axis
// ---------------------------------------------------------------------

/// Fit every stratum, letting the thin ones borrow tracts from their neighbouring repeat
/// counts.
///
/// **A stratum holding at least [`SsrFitConfig::borrowing_floor`] tracts is fitted on its
/// own.** A thinner one takes in its nearest neighbours of the same motif length, nearest
/// repeat count first and alternating either side, until the pooled set reaches the floor or
/// there is nothing left to take. Motif length is never crossed: a dinucleotide tract and a
/// trinucleotide tract of the same repeat count slip at different rates and pooling them would
/// mix two answers.
///
/// **What borrowing costs is bias and what it buys is spread**, and the two run opposite ways:
/// a neighbouring repeat count genuinely slips at a different rate, so a borrowed answer is
/// steady but displaced, where a thin stratum's own answer is unbiased but moves from draw to
/// draw. The floor is where those two costs cross.
pub fn fit_strata(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<StratumOutcome> {
    // **Two strata that borrow the identical set are one fit, not two.** Where a whole motif
    // length is thinner than the floor, every stratum in it pools every other, so the answers
    // are the same object computed as many times as there are strata; on tomato that is 71
    // fits over 4,164 tracts collapsing to a handful.
    let mut done: BTreeMap<Vec<Stratum>, Option<StratumFit>> = BTreeMap::new();
    let mut outcomes = Vec::with_capacity(strata.len());
    for evidence in strata {
        if evidence.tracts_with_reads() == 0 {
            outcomes.push(StratumOutcome::Refused {
                stratum: evidence.stratum,
                tracts: 0,
                reason: StratumRefusal::NoSpanningReads,
            });
            continue;
        }

        let (pooled, borrowed) = borrow_up_to_the_floor(evidence, strata, config);
        let reached = pooled.tracts_with_reads();
        if reached < config.refusal_floor {
            outcomes.push(StratumOutcome::Refused {
                stratum: evidence.stratum,
                tracts: reached,
                reason: StratumRefusal::BelowTheFloor {
                    tracts: reached,
                    floor: config.refusal_floor,
                },
            });
            continue;
        }

        let mut key: Vec<Stratum> = borrowed
            .iter()
            .map(|repeats| Stratum {
                period: evidence.stratum.period,
                reference_repeats: *repeats,
            })
            .chain(std::iter::once(evidence.stratum))
            .collect();
        key.sort_unstable();
        let shared = done
            .entry(key)
            .or_insert_with(|| fit_pooled(&pooled, &borrowed, homozygote_excess, config));

        outcomes.push(match shared {
            Some(fit) => {
                let mut mine = fit.clone();
                mine.stratum = evidence.stratum;
                mine.borrowed = borrowed;
                StratumOutcome::Fitted(Box::new(mine))
            }
            None => StratumOutcome::Refused {
                stratum: evidence.stratum,
                tracts: reached,
                reason: StratumRefusal::NoSpanningReads,
            },
        });
    }
    outcomes
}

/// Take in neighbours of the same motif length, nearest repeat count first, until the floor is
/// reached or the neighbours run out.
fn borrow_up_to_the_floor(
    evidence: &StratumEvidence,
    strata: &[StratumEvidence],
    config: &SsrFitConfig,
) -> (StratumEvidence, Vec<u64>) {
    let mut pooled = evidence.clone();
    let mut borrowed = Vec::new();
    if pooled.tracts_with_reads() >= config.borrowing_floor {
        return (pooled, borrowed);
    }

    // Same period, ordered by how far their repeat count is from this one.
    let here = evidence.stratum;
    let mut neighbours: Vec<&StratumEvidence> = strata
        .iter()
        .filter(|other| {
            other.stratum.period == here.period
                && other.stratum.reference_repeats != here.reference_repeats
                && other.tracts_with_reads() > 0
        })
        .collect();
    neighbours.sort_by_key(|other| {
        let distance = other
            .stratum
            .reference_repeats
            .abs_diff(here.reference_repeats);
        (distance, other.stratum.reference_repeats)
    });

    // **Both sides of a repeat count are taken together, never one and then a test.** Slippage
    // rises with repeat count — roughly 1.3-fold a count (`parameter_prepass_ssr.md` §4.3) — so
    // a rule that stops as soon as the floor is reached would keep whichever neighbour it
    // happened to reach first and carry that neighbour's whole displacement. Taking the ring
    // leaves the two displacements pointing opposite ways.
    let mut taken = 0;
    while taken < neighbours.len() {
        if pooled.tracts_with_reads() >= config.borrowing_floor {
            break;
        }
        let distance = neighbours[taken]
            .stratum
            .reference_repeats
            .abs_diff(here.reference_repeats);
        while taken < neighbours.len()
            && neighbours[taken]
                .stratum
                .reference_repeats
                .abs_diff(here.reference_repeats)
                == distance
        {
            pooled = pooled.pooled_with(neighbours[taken]);
            borrowed.push(neighbours[taken].stratum.reference_repeats);
            taken += 1;
        }
    }
    (pooled, borrowed)
}

// ---------------------------------------------------------------------
// The likelihood
// ---------------------------------------------------------------------

/// The genotypes of a diploid over `classes` allele lengths, unordered pairs.
fn genotype_pairs(classes: usize) -> Vec<(usize, usize)> {
    (0..classes)
        .flat_map(|first| (first..classes).map(move |second| (first, second)))
        .collect()
}

/// One tract's read likelihoods, in the form the integral sweeps over cheaply.
///
/// **Rescaled out of log space, once.** The inner loop is *for each point of the integral, for
/// each sample, sum over genotypes*, so doing it in logs would cost one `exp` per genotype per
/// sample **per point** — which at 256 points is where the whole run's time would go.
/// Subtracting each sample's largest log-likelihood makes the inner sum a plain dot product,
/// and the offsets are added back at the end.
struct TractLikelihoods {
    /// `scaled[i]` is the `i`-th sample-with-reads' row over the genotypes.
    scaled: Vec<Vec<f64>>,
    /// That sample's homozygote excess, in the same order.
    homozygote_excess: Vec<f64>,
    /// Σ over those samples of the log-likelihood each row was divided by.
    ln_offset: f64,
}

impl TractLikelihoods {
    fn of(
        tract: &TractReads,
        per_group_allele: &[Vec<Vec<f64>>],
        genotypes: &[(usize, usize)],
        homozygote_excess: &[f64],
    ) -> Self {
        let mut scaled = Vec::with_capacity(tract.samples.len());
        let mut excess = Vec::with_capacity(tract.samples.len());
        let mut ln_offset = 0.0;
        for sample in &tract.samples {
            let logs: Vec<f64> = genotypes
                .iter()
                .map(|(first, second)| {
                    let mut total = 0.0;
                    for (group, counts) in &sample.by_group {
                        let per_allele = &per_group_allele[*group as usize];
                        for (bucket, reads) in counts.iter().enumerate() {
                            if *reads == 0 {
                                continue;
                            }
                            let probability =
                                0.5 * (per_allele[*first][bucket] + per_allele[*second][bucket]);
                            total += f64::from(*reads) * probability.max(1e-300).ln();
                        }
                    }
                    total
                })
                .collect();
            let largest = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ln_offset += largest;
            scaled.push(logs.into_iter().map(|log| (log - largest).exp()).collect());
            excess.push(
                *homozygote_excess
                    .get(sample.sample as usize)
                    .expect("every sample in a tract has a homozygote excess"),
            );
        }
        Self {
            scaled,
            homozygote_excess: excess,
            ln_offset,
        }
    }
}

/// One point of the simplex a tract's length frequencies are integrated over, with the weight
/// it stands for.
struct Quadrature {
    frequencies: Vec<Vec<f64>>,
    ln_weight: f64,
}

/// `ln P(this tract's panel | parameters)`.
fn ln_tract(
    likelihoods: &TractLikelihoods,
    quadrature: &Quadrature,
    genotypes: &[(usize, usize)],
) -> f64 {
    let mut terms = Vec::with_capacity(quadrature.frequencies.len());
    // The genotype prior splits into the part that comes from the two copies being identical
    // by descent and the part that comes from two independent draws, so a sample's own
    // homozygote excess weights two dot products rather than rebuilding the prior per sample.
    let mut identical = vec![0.0_f64; genotypes.len()];
    let mut independent = vec![0.0_f64; genotypes.len()];
    for point in &quadrature.frequencies {
        for (slot, (first, second)) in genotypes.iter().enumerate() {
            if first == second {
                identical[slot] = point[*first];
                independent[slot] = point[*first] * point[*first];
            } else {
                identical[slot] = 0.0;
                independent[slot] = 2.0 * point[*first] * point[*second];
            }
        }
        let mut total = quadrature.ln_weight;
        for (row, excess) in likelihoods
            .scaled
            .iter()
            .zip(&likelihoods.homozygote_excess)
        {
            let mut by_descent = 0.0;
            let mut at_random = 0.0;
            for slot in 0..genotypes.len() {
                by_descent += identical[slot] * row[slot];
                at_random += independent[slot] * row[slot];
            }
            let sum = excess * by_descent + (1.0 - excess) * at_random;
            if sum <= 0.0 {
                total = f64::NEG_INFINITY;
                break;
            }
            total += sum.ln();
        }
        terms.push(total);
    }
    ln_sum_exp(&terms) + likelihoods.ln_offset
}

/// The evidence with the read likelihoods already computed for one set of slippage numbers.
///
/// **Held between evaluations**, because the climb moves the spectrum and the concentration
/// far more often than it moves slippage, and neither changes a read likelihood. Recomputing
/// them inside the search is what once made this program too slow to run.
struct Prepared {
    slippage: Vec<Slippage>,
    tracts: Vec<TractLikelihoods>,
}

/// The evidence plus whichever slippage the last question was about.
struct Scorer<'a> {
    evidence: &'a StratumEvidence,
    homozygote_excess: &'a [f64],
    genotypes: &'a [(usize, usize)],
    allele_span: i32,
    quadrature_points: usize,
    prepared: Option<Prepared>,
    /// The last integral built, with the spectrum and concentration it was built for.
    ///
    /// **The slippage climb asks about dozens of parameter sets that leave the tract's length
    /// frequencies alone**, and building the integral costs 19 ms at thirteen allele classes —
    /// which on a thin stratum is more than the likelihood it feeds.
    held_quadrature: Option<(Vec<f64>, f64, Quadrature)>,
}

impl<'a> Scorer<'a> {
    fn new(
        evidence: &'a StratumEvidence,
        homozygote_excess: &'a [f64],
        genotypes: &'a [(usize, usize)],
        config: &SsrFitConfig,
    ) -> Self {
        Self {
            evidence,
            homozygote_excess,
            genotypes,
            allele_span: config.allele_span,
            quadrature_points: config.quadrature_points,
            prepared: None,
            held_quadrature: None,
        }
    }

    /// The mean log-likelihood a tract at these parameters.
    fn score(&mut self, parameters: &Parameters) -> f64 {
        self.refresh(&parameters.slippage);
        let stale = !self
            .held_quadrature
            .as_ref()
            .is_some_and(|(spectrum, held, _)| {
                *held == parameters.concentration && *spectrum == parameters.length_spectrum
            });
        if stale {
            self.held_quadrature = Some((
                parameters.length_spectrum.clone(),
                parameters.concentration,
                dirichlet_points(
                    &parameters.length_spectrum,
                    parameters.concentration,
                    self.quadrature_points,
                ),
            ));
        }
        let (_, _, quadrature) = self.held_quadrature.as_ref().expect("built above");
        let prepared = self.prepared.as_ref().expect("refreshed above");
        // Across tracts in parallel, but summed back in tract order: a parallel float sum
        // reorders the additions run to run, and a fit's whole output is a difference between
        // one set of parameters and another.
        let per_tract: Vec<f64> = prepared
            .tracts
            .par_iter()
            .map(|tract| ln_tract(tract, quadrature, self.genotypes))
            .collect();
        per_tract.iter().sum::<f64>() / per_tract.len().max(1) as f64
    }

    fn refresh(&mut self, slippage: &[Slippage]) {
        if self
            .prepared
            .as_ref()
            .is_some_and(|held| held.slippage == slippage)
        {
            return;
        }
        let per_group_allele: Vec<Vec<Vec<f64>>> = slippage
            .iter()
            .map(|group| {
                (-self.allele_span..=self.allele_span)
                    .map(|allele| group.read_probabilities(allele, self.evidence.read_span))
                    .collect()
            })
            .collect();
        let tracts = self
            .evidence
            .tracts
            .par_iter()
            .map(|tract| {
                TractLikelihoods::of(
                    tract,
                    &per_group_allele,
                    self.genotypes,
                    self.homozygote_excess,
                )
            })
            .collect();
        self.prepared = Some(Prepared {
            slippage: slippage.to_vec(),
            tracts,
        });
    }
}

// ---------------------------------------------------------------------
// The integral over a tract's length frequencies
// ---------------------------------------------------------------------

/// The first primes, one a stick-breaking dimension. Thirteen allele classes need twelve.
const HALTON_BASES: [usize; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];

/// A Dirichlet with this mean and concentration, over a **fixed** low-discrepancy point set
/// pushed through the stick-breaking Beta quantiles.
///
/// **The points are fixed and the quantile map is continuous in the concentration**, so the
/// objective is a smooth function of it rather than a jittery one — the same uniforms are
/// reused at every value the climb tries, which is what makes this quadrature rather than
/// Monte Carlo and what stops the search chasing sampling noise.
fn dirichlet_points(spectrum: &[f64], concentration: f64, points: usize) -> Quadrature {
    let classes = spectrum.len();
    let alpha: Vec<f64> = spectrum
        .iter()
        .map(|weight| (concentration * weight).max(1e-3))
        .collect();
    let frequencies: Vec<Vec<f64>> = (0..points)
        .into_par_iter()
        .map(|point| {
            let mut remaining = 1.0;
            let mut piece = Vec::with_capacity(classes);
            for stick in 0..classes - 1 {
                let a = alpha[stick];
                let b: f64 = alpha[stick + 1..].iter().sum();
                let uniform =
                    van_der_corput(point + 1, HALTON_BASES[stick.min(HALTON_BASES.len() - 1)]);
                let share = beta_quantile(uniform, a, b);
                piece.push(remaining * share);
                remaining *= 1.0 - share;
            }
            piece.push(remaining);
            piece
        })
        .collect();
    Quadrature {
        frequencies,
        ln_weight: -(points as f64).ln(),
    }
}

/// The `index`-th term of the van der Corput sequence in `base` — one coordinate of a Halton
/// point.
fn van_der_corput(mut index: usize, base: usize) -> f64 {
    let (mut fraction, mut result) = (1.0 / base as f64, 0.0);
    while index > 0 {
        result += (index % base) as f64 * fraction;
        index /= base;
        fraction /= base as f64;
    }
    result
}

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() {
        return largest;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9 — the same coefficients the harnesses carry.
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * std::f64::consts::TAU.ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// `I_x(a, b)`, by its continued fraction — enough for a Beta quantile by bisection.
fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (a * x.ln() + b * (1.0 - x).ln() + ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)).exp() / a;
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularised_incomplete_beta(1.0 - x, b, a);
    }
    let (mut f, mut c, mut d) = (1.0_f64, 1.0_f64, 0.0_f64);
    for index in 0..=200 {
        let m = index / 2;
        let numerator = if index == 0 {
            1.0
        } else if index % 2 == 0 {
            let m = m as f64;
            (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        } else {
            let m = m as f64;
            -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        };
        d = 1.0 + numerator * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;
        c = 1.0 + numerator / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        let step = c * d;
        f *= step;
        if (1.0 - step).abs() < 1e-12 {
            break;
        }
    }
    front * (f - 1.0)
}

/// The `p`-th quantile of `Beta(a, b)`, by bisection on the cumulative distribution.
fn beta_quantile(p: f64, a: f64, b: f64) -> f64 {
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        if regularised_incomplete_beta(middle, a, b) < p {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Golden-section on one coordinate, over a bracket of `span` either side of `start`.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> f64 {
    const GOLDEN: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (start - span, start + span);
    let (mut left, mut right) = (high - GOLDEN * (high - low), low + GOLDEN * (high - low));
    let (mut at_left, mut at_right) = (score(left), score(right));
    for _ in 0..16 {
        if at_left > at_right {
            high = right;
            right = left;
            at_right = at_left;
            left = high - GOLDEN * (high - low);
            at_left = score(left);
        } else {
            low = left;
            left = right;
            at_left = at_right;
            right = low + GOLDEN * (high - low);
            at_right = score(right);
        }
    }
    if at_left > at_right { left } else { right }
}

// ---------------------------------------------------------------------
// Reading the records
// ---------------------------------------------------------------------

/// Which stratum each kept STR locus is in, in the order the records index them.
///
/// **The record entry carries an index and nothing else** — no coordinates, no stratum — so
/// the order has to be rebuilt from the same kept-loci object the writer was given. That is
/// the loci of every stratum flattened together and sorted by contig, start and end.
///
/// A locus whose contig `contig_of` does not resolve is dropped, exactly as the writer drops
/// it, so the two lists stay the same length.
pub fn strata_of_kept_loci(
    loci: &CensusLoci,
    contig_of: &dyn Fn(&str) -> Option<ContigId>,
) -> Vec<Stratum> {
    let mut with_position: Vec<((u32, u64, u64), Stratum)> = loci
        .ssr()
        .iter_sorted()
        .into_iter()
        .flat_map(|((period, reference_repeats), segments)| {
            segments.iter().filter_map(move |segment| {
                contig_of(segment.chrom()).map(|contig| {
                    (
                        (contig.get(), segment.start(), segment.end()),
                        Stratum {
                            period,
                            reference_repeats,
                        },
                    )
                })
            })
        })
        .collect();
    with_position.sort_unstable_by_key(|(position, _)| *position);
    with_position
        .into_iter()
        .map(|(_, stratum)| stratum)
        .collect()
}

/// Gather one stratum's evidence a locus at a time, from a whole cohort's records.
///
/// `slippage_group_of` names, for each read group, which set of slippage numbers its reads are
/// drawn under. **One group per read group is the specified grain**; a run that knows several
/// read groups ran on one machine may pool them, and one that pools everything is saying it
/// cannot tell them apart.
///
/// Every sample must hold evidence for the same STR loci in the same order, which the
/// recording-terms check the cohort makes at its door has already refused to let fail silently.
///
/// **`strata` is one entry per kept tract, in genome order** — the stratum each tract is in.
/// The census stores a tract under an index within its own stratum, so this list is also what
/// says how many tracts each stratum holds, and a section of a different length means the loci
/// and the evidence were built from different selections.
///
/// # Panics
///
/// When a sample's section for a stratum is not as long as that stratum's share of `strata`.
pub fn gather_strata(
    cohort: &mut CohortCensusEvidence,
    strata: &[Stratum],
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
) -> Result<Vec<StratumEvidence>, CensusError> {
    let groups = slippage_group_of
        .values()
        .map(|group| *group as usize + 1)
        .max()
        .unwrap_or(1);
    let buckets = (2 * RECORDED_OFFSET_RANGE + 1) as usize;

    // How many tracts each stratum holds, which is the length its sections are built at.
    let mut tracts_in: BTreeMap<Stratum, usize> = BTreeMap::new();
    for stratum in strata {
        *tracts_in.entry(*stratum).or_insert(0) += 1;
    }
    let band: Vec<Stratum> = tracts_in.keys().copied().collect();
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();

    // **The whole band at once, and that is this step's shape rather than the design's.** The
    // slippage fit borrows a thin stratum from its neighbours across the whole set
    // (`fit_strata`), so it is handed every stratum together; how many may be resident at once
    // is a measurement the fit specification owns (§11, questions 8 and 10) and not something
    // this call can decide.
    cohort.with_strata(&band, |lent| {
        // Each sample's tract sections, gathered by stratum. **One row a sample and not one
        // value**, because a stratum is fitted from every sample with reads in it at once.
        let by_sample: Vec<BTreeMap<Stratum, Vec<(ReadGroupId, &SsrEvidence)>>> = lent
            .iter()
            .zip(&names)
            .map(|(sections, name)| {
                let mut by_stratum: BTreeMap<Stratum, Vec<(ReadGroupId, &SsrEvidence)>> =
                    BTreeMap::new();
                for (group, stratum, records) in sections {
                    assert_eq!(
                        records.len(),
                        tracts_in.get(stratum).copied().unwrap_or(0),
                        "sample {} holds {} tracts at period {} and {} repeats where the \
                         selection has {}",
                        name,
                        records.len(),
                        stratum.period,
                        stratum.reference_repeats,
                        tracts_in.get(stratum).copied().unwrap_or(0)
                    );
                    by_stratum
                        .entry(*stratum)
                        .or_default()
                        .push((*group, records));
                }
                by_stratum
            })
            .collect();

        tracts_in
            .iter()
            .map(|(stratum, tracts)| {
                let mut evidence = StratumEvidence {
                    stratum: *stratum,
                    tracts: Vec::new(),
                    read_span: RECORDED_OFFSET_RANGE,
                    groups,
                    tracts_over_guard_threshold: 0,
                    reads_reaching_not_crossing: 0,
                    guard_reads: 0,
                };
                // **The reads that reached a tract and crossed none of it are counted per
                // stratum by the writer**, so they are read once here rather than accumulated a
                // locus at a time (spec §3). Every sample's every read group contributes its own
                // total.
                for sample in &by_sample {
                    for (_, records) in sample.get(stratum).into_iter().flatten() {
                        evidence.reads_reaching_not_crossing += records.covering_not_crossing();
                    }
                }

                for locus in 0..*tracts {
                    let mut reads = TractReads::default();
                    let mut over_guard = false;
                    for (sample_index, sample) in by_sample.iter().enumerate() {
                        let mut by_group: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
                        for (read_group, records) in sample.get(stratum).into_iter().flatten() {
                            if records.guard_is_over_threshold(locus) {
                                over_guard = true;
                            }
                            evidence.guard_reads += records
                                .guard()
                                .iter()
                                .filter(|entry| entry.locus as usize == locus)
                                .map(|entry| u64::from(entry.reads))
                                .sum::<u64>();
                            if records.state(locus) != SsrLocusState::Crossed {
                                continue;
                            }
                            let group = *slippage_group_of.get(read_group).unwrap_or(&0);
                            let counts = by_group.entry(group).or_insert_with(|| vec![0; buckets]);
                            let offsets = records.offsets(locus);
                            for (bucket, count) in counts.iter_mut().enumerate() {
                                *count +=
                                    u32::from(offsets.at(bucket as i32 - RECORDED_OFFSET_RANGE));
                            }
                        }
                        if !by_group.is_empty() {
                            reads.samples.push(SampleTractReads {
                                sample: sample_index as u32,
                                by_group: by_group.into_iter().collect(),
                            });
                        }
                    }
                    if over_guard {
                        evidence.tracts_over_guard_threshold += 1;
                        continue;
                    }
                    evidence.tracts.push(reads);
                }
                evidence
            })
            .collect()
    })
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// A reproducible stream, the same one the harnesses use.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        fn uniform(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
        }

        fn gamma(&mut self, shape: f64) -> f64 {
            if shape < 1.0 {
                let u = self.uniform().max(1e-300);
                return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let x = self.normal();
                let v = (1.0 + c * x).powi(3);
                if v <= 0.0 {
                    continue;
                }
                let u = self.uniform().max(1e-300);
                if u.ln() < 0.5 * x * x + d - d * v + d * (v.ln()) {
                    return d * v;
                }
            }
        }

        fn normal(&mut self) -> f64 {
            let u1 = self.uniform().max(1e-300);
            let u2 = self.uniform();
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        }

        fn dirichlet(&mut self, alpha: &[f64]) -> Vec<f64> {
            let draws: Vec<f64> = alpha.iter().map(|a| self.gamma(*a).max(1e-300)).collect();
            let total: f64 = draws.iter().sum();
            draws.into_iter().map(|d| d / total).collect()
        }

        fn categorical(&mut self, weights: &[f64]) -> usize {
            let mut u = self.uniform();
            for (index, weight) in weights.iter().enumerate() {
                u -= weight;
                if u <= 0.0 {
                    return index;
                }
            }
            weights.len() - 1
        }
    }

    /// Draw one stratum: `tracts` tracts, `samples` samples, `depth` reads a sample a tract.
    fn draw_stratum(
        slippage: Slippage,
        spectrum: &[f64],
        concentration: f64,
        homozygote_excess: f64,
        tracts: usize,
        samples: usize,
        depth: u32,
        span: i32,
        seed: u64,
    ) -> StratumEvidence {
        let classes = spectrum.len();
        let buckets = (2 * span + 1) as usize;
        let per_allele: Vec<Vec<f64>> = (0..classes)
            .map(|class| slippage.read_probabilities(class as i32 - span, span))
            .collect();
        let mut rng = Rng(seed);
        let mut drawn = Vec::with_capacity(tracts);
        for _ in 0..tracts {
            let alpha: Vec<f64> = spectrum.iter().map(|q| concentration * q).collect();
            let frequencies = rng.dirichlet(&alpha);
            let mut reads = TractReads::default();
            for sample in 0..samples {
                let first = rng.categorical(&frequencies);
                let second = if rng.uniform() < homozygote_excess {
                    first
                } else {
                    rng.categorical(&frequencies)
                };
                let mut counts = vec![0_u32; buckets];
                for _ in 0..depth {
                    let allele = if rng.uniform() < 0.5 { first } else { second };
                    counts[rng.categorical(&per_allele[allele])] += 1;
                }
                reads.samples.push(SampleTractReads {
                    sample: sample as u32,
                    by_group: vec![(0, counts)],
                });
            }
            drawn.push(reads);
        }
        StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 10,
            },
            tracts: drawn,
            read_span: span,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 0,
            guard_reads: 0,
        }
    }

    /// The three-class spectrum every measurement on this path was made against: most
    /// chromosomes at the reference length, one repeat either side carrying most of the rest.
    fn spectrum_of(classes: usize) -> Vec<f64> {
        let middle = classes / 2;
        let mut spectrum: Vec<f64> = (0..classes)
            .map(|class| 0.55_f64.powi((class as i32 - middle as i32).abs()))
            .collect();
        normalise(&mut spectrum);
        spectrum
    }

    /// Every read distribution is a distribution: it sums to one, whatever the allele's
    /// distance from the recorded range.
    #[test]
    fn a_read_distribution_sums_to_one_from_every_allele() {
        let slippage = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        for allele in -6..=6 {
            let total: f64 = slippage.read_probabilities(allele, 4).iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-12,
                "allele {allele} gave {total}, not one"
            );
        }
    }

    /// **The end bucket carries the marginal**: an allele sitting three repeats outside the
    /// recorded range puts nearly everything in the end bucket rather than nothing there.
    #[test]
    fn an_allele_outside_the_recorded_range_lands_in_the_end_bucket() {
        let slippage = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let reads = slippage.read_probabilities(-6, 4);
        assert!(
            reads[0] > 0.9,
            "the shortest bucket took {}, not nearly everything",
            reads[0]
        );
        assert!(reads[8] < 0.01, "the far end took {}", reads[8]);
    }

    /// **The positive control.** A stratum drawn with a known truth comes back with it: this
    /// is the run that says the estimator has power, and without it a clean-looking answer
    /// cannot be told from one with no information in it.
    #[test]
    fn a_drawn_stratum_returns_the_numbers_it_was_drawn_with() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 1_500, 20, 6, 1, 11);
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let fitted = fit_stratum(&evidence, &vec![0.4; 20], &config).expect("reads were drawn");
        let slippage = fitted.slippage[0].expect("the one group has reads");

        assert!(
            (slippage.level - truth.level).abs() / truth.level < 0.10,
            "slippage level {} against a truth of {}",
            slippage.level,
            truth.level
        );
        assert!(
            (slippage.shorter_share - truth.shorter_share).abs() < 0.05,
            "shorter-share {} against a truth of {}",
            slippage.shorter_share,
            truth.shorter_share
        );
        assert!(
            (fitted.concentration - 0.5).abs() / 0.5 < 0.25,
            "concentration {} against a truth of 0.5",
            fitted.concentration
        );
    }

    /// A stratum with no reads at all is refused rather than fitted, and it is refused by name.
    #[test]
    fn a_stratum_with_no_reads_is_refused_and_not_fitted() {
        let evidence = StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 9,
            },
            tracts: vec![TractReads::default(); 20],
            read_span: 4,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 40,
            guard_reads: 0,
        };
        let outcomes = fit_strata(&[evidence], &[0.4], &SsrFitConfig::default());
        assert!(matches!(
            &outcomes[0],
            StratumOutcome::Refused {
                reason: StratumRefusal::NoSpanningReads,
                ..
            }
        ));
    }

    /// **A thin stratum borrows from its nearest repeat counts and never across motif
    /// lengths.** The borrowing is what a run has to be able to read off the answer, so it is
    /// recorded rather than only acted on.
    #[test]
    fn a_thin_stratum_borrows_nearest_first_within_its_own_motif_length() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let mut strata = Vec::new();
        for (period, repeats, tracts, seed) in [
            (2_u8, 10_u64, 60_usize, 3_u64),
            (2, 11, 400, 5),
            (2, 13, 400, 7),
            (3, 10, 400, 9),
        ] {
            let mut evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, tracts, 8, 6, 1, seed);
            evidence.stratum = Stratum {
                period,
                reference_repeats: repeats,
            };
            strata.push(evidence);
        }
        let config = SsrFitConfig {
            allele_span: 1,
            borrowing_floor: 300,
            max_rounds: 1,
            starting_points: vec![StartingPoint {
                slippage_level: 0.10,
                concentration: 3.0,
            }],
            ..SsrFitConfig::default()
        };
        let outcomes = fit_strata(&strata, &vec![0.4; 8], &config);
        let StratumOutcome::Fitted(thin) = &outcomes[0] else {
            panic!("the thin stratum should be fitted from borrowed tracts");
        };
        assert_eq!(
            thin.borrowed,
            vec![11],
            "the nearest repeat count of the same motif length, and only as far as the floor"
        );
        assert_eq!(thin.tracts_fitted, 460);

        let StratumOutcome::Fitted(fat) = &outcomes[1] else {
            panic!("a stratum above the floor is fitted");
        };
        assert!(
            fat.borrowed.is_empty(),
            "a stratum above the floor borrows nothing"
        );
    }
}
