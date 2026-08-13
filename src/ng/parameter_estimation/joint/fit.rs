//! The estimator: every parameter fitted once, over every sample at the same loci.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_fit.md`. Types:
//! `doc/devel/ng/arch/parameter_prepass_joint_fit.md`. It reads the records of
//! [`records`](super::records) at the loci of [`loci`](super::loci) and nothing else.
//!
//! **This is the generic path.** The repeat-tract half of the spec (§4) is not here.
//!
//! # What having every sample at one locus buys, in one sentence
//!
//! A position's own allele frequency in the population becomes a quantity the fit can weight
//! a genotype against — so a position where every sample looks heterozygous, which is an
//! artefact, is told from one where a quarter do, which is a variant. A route that folds each
//! sample into a histogram has forgotten which position was which and cannot ask.
//!
//! # What is free, and what is not
//!
//! Two numbers per read group (how often a read misreads at an ordinary position, and at a
//! mismapped one), one cohort-wide share of mismapped positions, four numbers describing how
//! the population's allele frequency is distributed, and one number per sample saying how
//! much less heterozygous it is than random mating would predict. **A position holds no
//! parameter of its own**: its frequency is summed over rather than fitted, which is what
//! keeps the parameter count from growing with the data (spec §2.1).
//!
//! # Why this is an EM and not a climb over the likelihood
//!
//! The architecture proposes coordinate climbs — one bounded search per parameter, each
//! costing about twenty evaluations of the whole likelihood. **One evaluation is a pass over
//! two million positions times fifty samples times the quadrature**, so a hundred parameters
//! searched that way is a thousand passes over the data and the program cannot be run. That
//! is the third of the four traps recorded in
//! `doc/devel/ng/reports/joint_route_research_narrative_2026-08-13.md` §7, and it is avoided
//! the same way: **one pass over the data per iteration**, which accumulates the counts every
//! parameter's own maximisation needs, and then every maximisation runs over those counts
//! with the reads untouched. The searches that remain — the error rates, the Beta's two
//! shapes, each sample's homozygote excess — are one- and two-dimensional and see a few dozen
//! accumulated numbers apiece.
//!
//! Each iteration cannot lower the likelihood, which is what expectation-maximisation gives
//! and a coordinate climb over a non-concave surface does not.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rayon::prelude::*;

use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::types::{Ploidy, ReadGroupId};

use super::contamination::{ContaminationEstimate, fit_contamination};
use super::records::{DepthCode, SampleRecords};

// ---------------------------------------------------------------------
// What the route emits
// ---------------------------------------------------------------------

/// A read group's noise, as the classes a position can be drawn from.
///
/// **`clean` and `noisy` are the pair the histogram route already fits.** A clean position
/// disagrees with the reference only when the sequencer misreads a base; a noisy one is
/// mismapped and disagrees far more often. Which of the two a position is, is a property of
/// **the position** and not of the sample — a collapsed duplicate in the reference mismaps in
/// everybody — so the share of noisy positions is one cohort-wide number and lives on
/// [`JointFit`], while the two rates are chemistry and live here.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SiteClassNoise {
    /// How often a read misreads a base at an ordinary position.
    pub clean: f64,
    /// How often a read disagrees with the reference at a mismapped one.
    pub noisy: f64,
}

/// How the population's allele frequency is distributed at an ordinary position: a mass on
/// *the population carries only the reference base*, a mass on *it carries only a
/// non-reference base*, and a Beta over what actually segregates.
///
/// **Four fitted numbers, and the grid the integral is taken on is accuracy rather than
/// freedom** — doubling the nodes costs time and adds no parameter (spec §2.1.2). A free
/// weight per grid point would be an unregularised deconvolution, which collapses onto spikes.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct FrequencyDensity {
    /// The share of positions where the population carries only the reference base.
    pub p_invariant: f64,
    /// The share where it carries only a non-reference base — the reference accession's own
    /// private alleles, which on a crop reference is not a rounding term.
    pub p_fixed_alt: f64,
    /// The Beta's shape over the positions that segregate. Below one it reproduces the
    /// rare-allele pile-up a neutral population has.
    pub a: f64,
    pub b: f64,
}

impl FrequencyDensity {
    /// The share of positions that segregate — one minus the two masses.
    pub fn p_segregating(&self) -> f64 {
        (1.0 - self.p_invariant - self.p_fixed_alt).max(0.0)
    }

    /// `∫ π(f) · 2 f (1 − f) df` — the population's expected heterozygosity, with no
    /// finite-sample correction because there is no panel in it (spec §5.3).
    ///
    /// Only the segregating part contributes: a population carrying one allele has none.
    pub fn expected_heterozygosity(&self) -> f64 {
        let (a, b) = (self.a, self.b);
        // E[2f(1−f)] under Beta(a, b) is 2·a·b / ((a+b)(a+b+1)).
        self.p_segregating() * 2.0 * a * b / ((a + b) * (a + b + 1.0))
    }
}

/// How much less heterozygous an individual is than random mating in the panel would predict.
///
/// **A different quantity from the autozygosity coefficient the caller's genotype prior
/// multiplies**, which is the fraction of a genome lying in runs of homozygosity. A consumer
/// handed the wrong one gets a plausible answer, which is why they are two types and not two
/// values of one (spec §5.1).
///
/// **Bounded below by zero.** An unconstrained fit will go negative under a heterozygote
/// excess, and a negative value books one sample's mismapping as biology, with a plausible
/// number and no error.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct HomozygoteExcess(f64);

impl HomozygoteExcess {
    pub fn try_new(value: f64) -> Option<Self> {
        (0.0..=1.0).contains(&value).then_some(Self(value))
    }

    pub fn get(self) -> f64 {
        self.0
    }
}

/// What a sample's genotypes came out as, once the fit converged.
///
/// **Derived from the posteriors rather than fitted** (spec §3.2), which is a real difference
/// from the per-sample route and the reason a disagreement between the two is not
/// automatically an error in either.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SampleGenotypeRates {
    /// The mean over kept positions of the posterior that the sample is heterozygous.
    pub heterozygous: f64,
    /// …and that every copy is non-reference.
    pub homozygous_alt: f64,
    /// **How many kept positions carried at least one read, and at least two.**
    ///
    /// Reported beside the rates because a position where a sample has no read has a
    /// posterior equal to its prior, so it contributes the answer rather than evidence for
    /// it. At three reads a position about 5 in 100 have no read at all and 15 in 100 have
    /// exactly one (spec §3.2), and the rates must not be read as counts of observed
    /// heterozygotes when a fifth of their support saw nothing.
    pub positions_with_reads: u64,
    pub positions_with_two_reads: u64,
}

/// Every parameter this route produces, for the whole cohort, in one value.
#[derive(Clone, PartialEq, Debug)]
pub struct JointFit {
    /// Chemistry: two rates per read group.
    pub noise: BTreeMap<ReadGroupId, Estimate<SiteClassNoise>>,
    /// The share of positions drawn from the noisy class. **One number for the cohort**, not
    /// one per read group: a position is mismapped or it is not, and it cannot be mismapped
    /// for one library and clean for another at the same position. *This departs from spec
    /// §3.1's table, which lists it per read group; that grain belongs to the histogram
    /// route, where a sample's own libraries are all it ever sees at once.*
    pub noisy_share: f64,
    /// The population's allele-frequency density — this route's own parameter.
    pub density: Estimate<FrequencyDensity>,
    /// Per sample, the departure from the Hardy–Weinberg proportions the density predicts.
    pub hom_excess: BTreeMap<String, Estimate<HomozygoteExcess>>,
    /// Per sample, derived from the converged posteriors rather than fitted.
    pub rates: BTreeMap<String, Estimate<SampleGenotypeRates>>,
    /// Per sample, the fraction of reads from another individual.
    pub contamination: BTreeMap<String, ContaminationEstimate>,
    /// The population's expected heterozygosity, read off the density.
    pub expected_heterozygosity: f64,
    /// How the fit ended. **Running out of passes is never reported as convergence.**
    pub passes: u32,
    pub converged: bool,
    /// The log-likelihood at the returned parameters.
    pub log_likelihood: f64,
}

/// Why a fit could not be produced.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum JointFitError {
    /// Samples that did not keep the same loci. **Refuses rather than averaging**, and runs
    /// before any arithmetic: a run that would fail on the fiftieth sample fails before the
    /// first likelihood evaluation.
    #[error("samples {first} and {second} disagree on {field}; they did not keep the same loci")]
    IdentityMismatch {
        first: String,
        second: String,
        field: &'static str,
    },
    #[error("the joint fit needs at least one sample")]
    NoSamples,
    #[error(
        "this route fits diploids; sample {sample} was given ploidy {ploidy}, and the \
         homozygote excess has no agreed form above two copies"
    )]
    NotDiploid { sample: String, ploidy: u8 },
}

// ---------------------------------------------------------------------
// What the run was asked for, beside the records
// ---------------------------------------------------------------------

/// Where one run of the maximisation starts.
///
/// **The starting points are part of the estimator and not a tuning detail.** A start that
/// puts the ordinary and the mismapped class close together collapses them into one and
/// reports convergence — measured at 46% of the clean error rate, with the population's
/// heterozygosity 10.6% high — while a start that puts them far apart returns the clean rate
/// to within 1.3% on the identical data (spec §11 question 3). So a run enumerates several and
/// keeps the best-scoring fit rather than the last one.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct StartingPoint {
    pub clean: f64,
    pub noisy: f64,
    pub noisy_share: f64,
    pub p_invariant: f64,
    pub p_fixed_alt: f64,
    pub a: f64,
    pub b: f64,
}

impl StartingPoint {
    /// Three starts spanning the separation between the two classes: touching, an order of
    /// magnitude apart, and two orders apart.
    pub fn spanning_the_class_separation() -> Vec<Self> {
        [(0.002, 0.004), (0.002, 0.05), (0.0005, 0.15)]
            .into_iter()
            .map(|(clean, noisy)| Self {
                clean,
                noisy,
                noisy_share: 0.01,
                p_invariant: 0.97,
                p_fixed_alt: 0.005,
                a: 0.5,
                b: 2.0,
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct JointFitConfig {
    pub ploidy: Ploidy,
    /// How many nodes the integral over a position's allele frequency is taken on. **An
    /// accuracy knob and not a parameter count.**
    pub quadrature_nodes: usize,
    pub starting_points: Vec<StartingPoint>,
    pub max_passes: u32,
    /// The pass stops when no parameter moves by more than this, relative to itself.
    pub stillness: f64,
    /// The ladder the records' depth codes index. Two samples binned under different edges
    /// hold codes that mean different depths, which the identity check already refuses.
    pub edges: Arc<DepthBinEdges>,
    /// How many axes of variation each sample's own allele frequency is a straight line in,
    /// which is what makes contamination measurable on a diverged panel
    /// ([`contamination`](super::contamination)). Zero turns it off and leaves every sample
    /// *not identified*.
    pub components: usize,
}

impl Default for JointFitConfig {
    fn default() -> Self {
        Self {
            ploidy: Ploidy::try_new(2).expect("two is a ploidy"),
            quadrature_nodes: 16,
            starting_points: StartingPoint::spanning_the_class_separation(),
            max_passes: 200,
            stillness: 1e-4,
            edges: Arc::new(DepthBinEdges::new()),
            components: crate::ng::parameter_estimation::joint::contamination::DEFAULT_COMPONENTS,
        }
    }
}

// ---------------------------------------------------------------------
// The evidence, laid out the way the pass over it wants
// ---------------------------------------------------------------------

/// One sample's reads at one position, as the likelihood reads them.
///
/// Counts, not codes: `depth` is the ladder's own answer for the stored code, and `on[k]` is
/// how many reads showed base `k`. **The reference base is not here and is not needed** — the
/// records list only what disagreed with it, so reads on the reference are what is left over.
#[derive(Copy, Clone, Default, Debug)]
struct SampleAtPosition {
    depth: f64,
    /// Reads on each of the four bases that are **not** the reference, indexed by allele code;
    /// the reference base's own entry is always zero, and so is any base no read showed.
    on: [f64; 5],
}

impl SampleAtPosition {
    /// Reads that disagreed with the reference in any way.
    fn non_reference(&self) -> f64 {
        self.on.iter().sum()
    }
}

/// The cohort at one position.
struct PositionEvidence {
    /// One entry per sample, in the order the fit iterates samples.
    samples: Vec<SampleAtPosition>,
    /// Which non-reference bases any sample showed a read on. **The candidates the fit sums
    /// over are the three bases that are not the reference**, and the ones nobody showed a
    /// read on all give the identical term, so they are counted rather than enumerated.
    observed_alternatives: Vec<usize>,
}

/// Walks every sample's records once, in position order, handing one position at a time.
///
/// **A cursor per sample rather than a search per position.** Each sample's non-reference
/// observations are sorted by position, so advancing a cursor costs one comparison; a binary
/// search per sample per position would cost twenty-one, two million times over.
struct EvidenceCursor<'a> {
    samples: &'a [SampleRecords],
    edges: &'a DepthBinEdges,
    /// Per sample, per read group, where that group's sparse list has been read to.
    at: Vec<Vec<usize>>,
    positions: usize,
    next: usize,
}

impl<'a> EvidenceCursor<'a> {
    fn position_count(samples: &[SampleRecords]) -> usize {
        samples
            .first()
            .and_then(|s| s.generic.values().next())
            .map_or(0, |g| g.depth().len())
    }

    /// A cursor over `first..end` only — **what lets one pass over the positions be split
    /// across cores.** Each sample's sparse list is binary-searched once to find where the
    /// chunk begins, and walked with a cursor from there, so a chunk costs one search per
    /// sample rather than one per position.
    fn over(
        samples: &'a [SampleRecords],
        edges: &'a DepthBinEdges,
        first: usize,
        end: usize,
    ) -> Self {
        let at = samples
            .iter()
            .map(|sample| {
                sample
                    .generic
                    .values()
                    .map(|records| {
                        records
                            .non_reference()
                            .partition_point(|entry| (entry.index as usize) < first)
                    })
                    .collect()
            })
            .collect();
        Self {
            samples,
            edges,
            at,
            positions: end,
            next: first,
        }
    }

    fn next_position(&mut self, into: &mut PositionEvidence) -> bool {
        if self.next >= self.positions {
            return false;
        }
        let index = self.next;
        self.next += 1;
        into.observed_alternatives.clear();
        let mut seen = [false; 5];
        for (s, sample) in self.samples.iter().enumerate() {
            let slot = &mut into.samples[s];
            *slot = SampleAtPosition::default();
            for (g, records) in sample.generic.values().enumerate() {
                if let DepthCode::Binned(bin) = records.depth().get(index) {
                    slot.depth += self.edges.representative_depth(bin);
                }
                let cursor = &mut self.at[s][g];
                let entries = records.non_reference();
                while *cursor < entries.len() && (entries[*cursor].index as usize) < index {
                    *cursor += 1;
                }
                while *cursor < entries.len() && entries[*cursor].index as usize == index {
                    let entry = entries[*cursor];
                    let code = usize::from(entry.allele.code());
                    slot.on[code] += f64::from(entry.reads);
                    seen[code] = true;
                    *cursor += 1;
                }
            }
            // A binned depth can read below the alternative count it has to carry, at the top
            // of the ladder where a bin is wide. The reads are the harder evidence, so the
            // depth gives way rather than charging a negative count of reference reads.
            slot.depth = slot.depth.max(slot.non_reference());
        }
        for (code, was_seen) in seen.iter().enumerate() {
            // `Other` — an indel or a spanning deletion — is not a base and cannot be the
            // segregating allele of a substitution model. Its reads are held out of both the
            // numerator and the denominator by `SampleAtPosition::reference_reads`.
            if *was_seen && code < 4 {
                into.observed_alternatives.push(code);
            }
        }
        true
    }
}

/// How many non-reference bases a position has, whatever the reference base is: three.
const CANDIDATE_ALTERNATIVES: usize = 3;

// ---------------------------------------------------------------------
// The likelihood of one sample's reads under one genotype
// ---------------------------------------------------------------------

/// `ln P(this sample's reads | j of its copies carry base k, error rate ε)`.
///
/// A read is drawn from one of the individual's copies and may be misread into any of the
/// other three bases, so the four bases form a multinomial with three distinct
/// probabilities: the candidate allele, the reference base, and the two bases that are
/// neither. The multinomial's coefficient is the same however the four categories are
/// labelled, so it is the same for every genotype, every class **and every candidate
/// allele**, and is dropped rather than computed.
///
/// The chance a read shows *anything* other than the reference base is
/// `carried·(1 − ε/3) + (1 − carried)·ε`, which is
/// [`alternative_read_probability`](crate::ng::parameter_estimation::generic::noise_model)
/// exactly — the two routes score a read the same way, which is what makes the comparison
/// between them a comparison.
fn ln_reads_given_genotype(
    sample: &SampleAtPosition,
    alt_copies: u8,
    ploidy: Ploidy,
    alternative: usize,
    error_rate: f64,
) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    let on_candidate = carried * (1.0 - error_rate) + (1.0 - carried) * error_rate / 3.0;
    let on_reference = (1.0 - carried) * (1.0 - error_rate) + carried * error_rate / 3.0;
    let on_neither = error_rate / 3.0;

    let candidate_reads = sample.on[alternative];
    // `Other` is neither the candidate nor the reference; it is held out of the model
    // entirely, so the depth it occupied is removed with it.
    let other_reads = sample.non_reference() - candidate_reads - sample.on[4];
    let reference_reads = (sample.depth - sample.non_reference()).max(0.0);

    count_times_ln(candidate_reads, on_candidate)
        + count_times_ln(other_reads, on_neither)
        + count_times_ln(reference_reads, on_reference)
}

/// `ln P(this sample's reads | every copy is the reference base)` — the term the invariant
/// branch needs, and it does not depend on which base would have been the alternative.
fn ln_reads_given_all_reference(sample: &SampleAtPosition, error_rate: f64) -> f64 {
    let non_reference = sample.non_reference() - sample.on[4];
    let reference_reads = (sample.depth - sample.non_reference()).max(0.0);
    count_times_ln(non_reference, error_rate / 3.0)
        + count_times_ln(reference_reads, 1.0 - error_rate)
}

/// `count · ln p`, with a count of zero contributing nothing — otherwise a category no read
/// fell into, at probability zero, would be `0 · −∞ = NaN`, and a `NaN` does not lose loudly.
fn count_times_ln(count: f64, probability: f64) -> f64 {
    if count == 0.0 {
        0.0
    } else {
        count * probability.max(f64::MIN_POSITIVE).ln()
    }
}

/// How common each genotype is in a diploid drawn from a population at frequency `f`, when
/// the individual is inbred by `excess`.
///
/// The heterozygote is depressed by `1 − excess` and the mass moves to the two homozygotes in
/// proportion to the frequency, which is what an inbreeding coefficient means. **It is this
/// departure from a pair of independent draws that forces the whole route to work with a
/// frequency rather than a count in the panel** (spec §2.1.1): the cancellation the count form
/// needs is exactly the factorisation that this breaks.
fn genotype_frequencies(f: f64, excess: f64) -> [f64; 3] {
    let heterozygous = 2.0 * f * (1.0 - f) * (1.0 - excess);
    let shift = excess * f * (1.0 - f);
    [(1.0 - f) * (1.0 - f) + shift, heterozygous, f * f + shift]
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if largest == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

// ---------------------------------------------------------------------
// The parameters, and one pass over the data
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Parameters {
    /// Per read group, in the order `groups` lists them.
    clean: Vec<f64>,
    noisy: Vec<f64>,
    noisy_share: f64,
    density: FrequencyDensity,
    /// Per sample.
    hom_excess: Vec<f64>,
}

/// What one pass over the data accumulates, and every maximisation then reads instead of the
/// reads.
struct Statistics {
    log_likelihood: f64,
    positions: f64,
    /// Posterior mass on the noisy class, and on each of the three branches.
    noisy: f64,
    invariant: f64,
    fixed_alt: f64,
    segregating: f64,
    /// Σ posterior · ln f and Σ posterior · ln(1 − f) over the quadrature, which is all a
    /// Beta's two shapes need.
    sum_ln_f: f64,
    sum_ln_one_minus_f: f64,
    /// Per read group × class × genotype: expected reads on the candidate allele, on neither
    /// base, and on the reference. **This is what replaces a data pass per candidate error
    /// rate**: the rate's own maximisation is over these nine numbers.
    reads: Vec<[[ReadTally; 3]; 2]>,
    /// Per sample × quadrature node × genotype: expected genotype counts on the segregating
    /// branch, which is all a homozygote excess needs.
    genotypes: Vec<Vec<[f64; 3]>>,
    /// Per sample, the posterior rates, summed over positions.
    heterozygous: Vec<f64>,
    homozygous_alt: Vec<f64>,
    with_reads: Vec<u64>,
    with_two_reads: Vec<u64>,
}

#[derive(Copy, Clone, Default)]
struct ReadTally {
    candidate: f64,
    neither: f64,
    reference: f64,
}

impl Statistics {
    fn new(groups: usize, samples: usize, nodes: usize) -> Self {
        Self {
            log_likelihood: 0.0,
            positions: 0.0,
            noisy: 0.0,
            invariant: 0.0,
            fixed_alt: 0.0,
            segregating: 0.0,
            sum_ln_f: 0.0,
            sum_ln_one_minus_f: 0.0,
            reads: vec![[[ReadTally::default(); 3]; 2]; groups],
            genotypes: vec![vec![[0.0; 3]; nodes]; samples],
            heterozygous: vec![0.0; samples],
            homozygous_alt: vec![0.0; samples],
            with_reads: vec![0; samples],
            with_two_reads: vec![0; samples],
        }
    }
}

impl Statistics {
    /// Add another chunk's counts. **Every field is a sum over positions**, which is what lets
    /// the pass be split across cores at all.
    fn absorb(&mut self, other: &Self) {
        self.log_likelihood += other.log_likelihood;
        self.positions += other.positions;
        self.noisy += other.noisy;
        self.invariant += other.invariant;
        self.fixed_alt += other.fixed_alt;
        self.segregating += other.segregating;
        self.sum_ln_f += other.sum_ln_f;
        self.sum_ln_one_minus_f += other.sum_ln_one_minus_f;
        for (into, from) in self.reads.iter_mut().zip(other.reads.iter()) {
            for (into, from) in into.iter_mut().zip(from.iter()) {
                for (into, from) in into.iter_mut().zip(from.iter()) {
                    into.candidate += from.candidate;
                    into.neither += from.neither;
                    into.reference += from.reference;
                }
            }
        }
        for (into, from) in self.genotypes.iter_mut().zip(other.genotypes.iter()) {
            for (into, from) in into.iter_mut().zip(from.iter()) {
                for (into, from) in into.iter_mut().zip(from.iter()) {
                    *into += from;
                }
            }
        }
        for (into, from) in self.heterozygous.iter_mut().zip(other.heterozygous.iter()) {
            *into += from;
        }
        for (into, from) in self
            .homozygous_alt
            .iter_mut()
            .zip(other.homozygous_alt.iter())
        {
            *into += from;
        }
        for (into, from) in self.with_reads.iter_mut().zip(other.with_reads.iter()) {
            *into += from;
        }
        for (into, from) in self
            .with_two_reads
            .iter_mut()
            .zip(other.with_two_reads.iter())
        {
            *into += from;
        }
    }
}

/// The read counts one sample contributes to a read group's tally, under one candidate allele.
fn tally_of(sample: &SampleAtPosition, alternative: usize) -> ReadTally {
    let candidate = sample.on[alternative];
    ReadTally {
        candidate,
        neither: sample.non_reference() - candidate - sample.on[4],
        reference: (sample.depth - sample.non_reference()).max(0.0),
    }
}

// ---------------------------------------------------------------------
// The entry point
// ---------------------------------------------------------------------

/// Fit every parameter from every sample's records, once.
///
/// **There is no per-sample entry point, because there is no per-sample answer.** Fitting
/// against a position's own allele frequency means a sample's evidence cannot be reduced
/// alone.
///
/// # Errors
///
/// [`JointFitError::IdentityMismatch`] before any arithmetic, when two samples did not keep
/// the same loci — the refusal [`loci`](super::loci) defines and this call enforces.
pub fn fit_jointly(
    samples: &[SampleRecords],
    config: &JointFitConfig,
) -> Result<JointFit, JointFitError> {
    let first = samples.first().ok_or(JointFitError::NoSamples)?;
    for other in &samples[1..] {
        if let Some(field) = first.identity.first_disagreement(&other.identity) {
            return Err(JointFitError::IdentityMismatch {
                first: first.sample.clone(),
                second: other.sample.clone(),
                field,
            });
        }
    }
    if config.ploidy.get() != 2 {
        return Err(JointFitError::NotDiploid {
            sample: first.sample.clone(),
            ploidy: config.ploidy.get(),
        });
    }

    let groups: Vec<ReadGroupId> = samples
        .iter()
        .flat_map(|sample| sample.generic.keys().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    // Which read group each sample's own groups are, in the order the cursor visits them.
    let group_index: Vec<Vec<usize>> = samples
        .iter()
        .map(|sample| {
            sample
                .generic
                .keys()
                .map(|id| {
                    groups
                        .iter()
                        .position(|g| g == id)
                        .expect("every group came from this list")
                })
                .collect()
        })
        .collect();

    let mut best: Option<(f64, Parameters, Statistics, u32, bool)> = None;
    for start in &config.starting_points {
        let (parameters, statistics, passes, converged) =
            maximise(samples, config, &groups, &group_index, start);
        let score = statistics.log_likelihood;
        if best.as_ref().is_none_or(|(current, ..)| score > *current) {
            best = Some((score, parameters, statistics, passes, converged));
        }
    }
    let (score, parameters, statistics, passes, converged) =
        best.expect("a run always has at least one starting point");

    let observations = statistics.positions as u64;
    let noise = groups
        .iter()
        .enumerate()
        .map(|(g, id)| {
            (
                *id,
                Estimate {
                    value: SiteClassNoise {
                        clean: parameters.clean[g],
                        noisy: parameters.noisy[g],
                    },
                    provenance: Provenance::FittedHere,
                    observations,
                },
            )
        })
        .collect();
    let hom_excess = samples
        .iter()
        .enumerate()
        .map(|(s, sample)| {
            (
                sample.sample.clone(),
                Estimate {
                    value: HomozygoteExcess::try_new(parameters.hom_excess[s])
                        .expect("the maximisation is confined to [0, 1]"),
                    provenance: if samples.len() >= 2 {
                        Provenance::FittedHere
                    } else {
                        Provenance::Defaulted
                    },
                    observations: statistics.with_reads[s],
                },
            )
        })
        .collect();
    let rates = samples
        .iter()
        .enumerate()
        .map(|(s, sample)| {
            (
                sample.sample.clone(),
                Estimate {
                    value: SampleGenotypeRates {
                        heterozygous: statistics.heterozygous[s] / statistics.positions,
                        homozygous_alt: statistics.homozygous_alt[s] / statistics.positions,
                        positions_with_reads: statistics.with_reads[s],
                        positions_with_two_reads: statistics.with_two_reads[s],
                    },
                    provenance: Provenance::FittedHere,
                    observations: statistics.with_reads[s],
                },
            )
        })
        .collect();
    // **Contamination is fitted after the alternation, not inside it** (spec §3.4). It reads
    // the converged error rates and the converged homozygote excess, and nothing it produces
    // feeds back into them — a sample's stray reads are a property of the tube it was in
    // rather than of the population, so the density has no business being told about them.
    let per_sample_error: Vec<f64> = samples
        .iter()
        .enumerate()
        .map(|(s, _)| parameters.clean[group_index[s][0]])
        .collect();
    let contamination = samples
        .iter()
        .map(|sample| sample.sample.clone())
        .zip(fit_contamination(
            samples,
            &config.edges,
            &per_sample_error,
            &parameters.hom_excess,
            config.components,
        ))
        .collect();

    Ok(JointFit {
        noise,
        noisy_share: parameters.noisy_share,
        density: Estimate {
            value: parameters.density,
            provenance: Provenance::FittedHere,
            observations,
        },
        hom_excess,
        rates,
        contamination,
        expected_heterozygosity: parameters.density.expected_heterozygosity(),
        passes,
        converged,
        log_likelihood: score,
    })
}

/// One run of the alternation, from one starting point.
fn maximise(
    samples: &[SampleRecords],
    config: &JointFitConfig,
    groups: &[ReadGroupId],
    group_index: &[Vec<usize>],
    start: &StartingPoint,
) -> (Parameters, Statistics, u32, bool) {
    let mut parameters = Parameters {
        clean: vec![start.clean; groups.len()],
        noisy: vec![start.noisy; groups.len()],
        noisy_share: start.noisy_share,
        density: FrequencyDensity {
            p_invariant: start.p_invariant,
            p_fixed_alt: start.p_fixed_alt,
            a: start.a,
            b: start.b,
        },
        hom_excess: vec![0.0; samples.len()],
    };
    let mut statistics;
    let mut converged = false;
    let mut passes = 0;
    let mut previous = f64::NEG_INFINITY;
    for pass in 1..=config.max_passes {
        passes = pass;
        statistics = expectation(samples, config, group_index, &parameters);
        let moved = maximisation(&mut parameters, &statistics, config);
        let gain = statistics.log_likelihood - previous;
        previous = statistics.log_likelihood;
        if moved < config.stillness && gain.abs() < config.stillness * statistics.positions.max(1.0)
        {
            converged = true;
            break;
        }
    }
    // The reported statistics must be the ones the reported parameters produce, so the last
    // maximisation is followed by one more pass rather than by the pass that preceded it.
    statistics = expectation(samples, config, group_index, &parameters);
    (parameters, statistics, passes, converged)
}

/// One pass over every position: the posteriors, and every count the maximisations need.
///
/// **Split across cores by position.** Positions are independent given the parameters, and a
/// chunk's counts add to another chunk's, so the pass is a map and a sum.
fn expectation(
    samples: &[SampleRecords],
    config: &JointFitConfig,
    group_index: &[Vec<usize>],
    parameters: &Parameters,
) -> Statistics {
    let quadrature = BetaQuadrature::with_genotype_priors(
        parameters.density.a,
        parameters.density.b,
        config.quadrature_nodes,
        &parameters.hom_excess,
    );
    let positions = EvidenceCursor::position_count(samples);
    let chunk = POSITIONS_PER_CHUNK.min(positions.div_ceil(rayon::current_num_threads()).max(1));
    let bounds: Vec<(usize, usize)> = (0..positions)
        .step_by(chunk)
        .map(|first| (first, (first + chunk).min(positions)))
        .collect();

    bounds
        .into_par_iter()
        .map(|(first, end)| {
            let mut statistics = Statistics::new(
                parameters.clean.len(),
                samples.len(),
                quadrature.nodes.len(),
            );
            let mut cursor = EvidenceCursor::over(samples, &config.edges, first, end);
            let mut scratch = Scratch::new(samples.len(), quadrature.nodes.len());
            while cursor.next_position(&mut scratch.evidence) {
                one_position(
                    &mut scratch,
                    group_index,
                    config.ploidy,
                    &quadrature,
                    parameters,
                    &mut statistics,
                );
            }
            statistics
        })
        .reduce(
            || {
                Statistics::new(
                    parameters.clean.len(),
                    samples.len(),
                    quadrature.nodes.len(),
                )
            },
            |mut into, from| {
                into.absorb(&from);
                into
            },
        )
}

/// How many positions one core takes at a time. Large enough that the per-chunk binary search
/// into each sample's sparse list disappears against the work, small enough that a cohort of
/// fifty still fills every core.
const POSITIONS_PER_CHUNK: usize = 16_384;

/// At most three non-reference bases can be the segregating one, so the candidate list never
/// exceeds three however many the samples showed reads on.
const MAX_CANDIDATES: usize = CANDIDATE_ALTERNATIVES;

/// Everything one position needs, sized once and reused, so a two-million-position pass
/// allocates nothing.
struct Scratch {
    evidence: PositionEvidence,
    samples: usize,
    nodes: usize,
    /// `[class][candidate][sample][genotype]` — `ln P(reads | genotype)`, and the same thing
    /// exponentiated after its own per-sample maximum is taken out.
    ///
    /// **This is the whole reason the program can be run.** The read likelihoods depend on the
    /// genotype and the error rate but **not on the allele frequency**, so they are computed
    /// once per candidate and reused across every quadrature node — which is the fix the
    /// previous session had to make to its repeat-tract program after it turned out to take
    /// hours (`joint_route_research_narrative_2026-08-13.md` §7).
    ell: Vec<f64>,
    lik: Vec<f64>,
    /// `[class][candidate][sample]` — the maximum taken out of `lik`.
    lik_max: Vec<f64>,
    /// `[class][candidate][node]`, and `[class][candidate]`.
    node_ln: Vec<f64>,
    fixed_ln: Vec<f64>,
    /// `[class]` — the invariant branch, and the three branches combined.
    invariant_ln: Vec<f64>,
    branch_ln: Vec<f64>,
    class_ln: Vec<f64>,
    /// The candidate alleles this position sums over, and how many alleles each stands for.
    candidates: Vec<usize>,
    multiplicity: Vec<f64>,
    /// `[sample][genotype]` — the segregating branch's genotype weight, collapsed over the
    /// nodes and candidates it was spread across.
    genotype_weight: Vec<f64>,
    /// `[candidate][sample][genotype]`, the same before the candidates are collapsed, because
    /// which reads count as the allele depends on which allele it is.
    per_candidate_weight: Vec<f64>,
    shares: Vec<f64>,
}

impl Scratch {
    fn new(samples: usize, nodes: usize) -> Self {
        Self {
            evidence: PositionEvidence {
                samples: vec![SampleAtPosition::default(); samples],
                observed_alternatives: Vec::with_capacity(MAX_CANDIDATES),
            },
            samples,
            nodes,
            ell: vec![0.0; 2 * MAX_CANDIDATES * samples * 3],
            lik: vec![0.0; 2 * MAX_CANDIDATES * samples * 3],
            lik_max: vec![0.0; 2 * MAX_CANDIDATES * samples],
            node_ln: vec![f64::NEG_INFINITY; 2 * MAX_CANDIDATES * nodes],
            fixed_ln: vec![f64::NEG_INFINITY; 2 * MAX_CANDIDATES],
            invariant_ln: vec![0.0; 2],
            branch_ln: vec![0.0; 6],
            class_ln: vec![0.0; 2],
            candidates: Vec::with_capacity(MAX_CANDIDATES),
            multiplicity: Vec::with_capacity(MAX_CANDIDATES),
            genotype_weight: vec![0.0; samples * 3],
            per_candidate_weight: vec![0.0; MAX_CANDIDATES * samples * 3],
            shares: Vec::with_capacity(MAX_CANDIDATES * nodes),
        }
    }

    fn ell_at(&self, class: usize, candidate: usize, sample: usize) -> usize {
        ((class * MAX_CANDIDATES + candidate) * self.samples + sample) * 3
    }

    fn max_at(&self, class: usize, candidate: usize, sample: usize) -> usize {
        (class * MAX_CANDIDATES + candidate) * self.samples + sample
    }

    fn node_at(&self, class: usize, candidate: usize, node: usize) -> usize {
        (class * MAX_CANDIDATES + candidate) * self.nodes + node
    }
}

/// The scale a running product is multiplied back up by when it is about to underflow, and its
/// logarithm. Sixty-three samples' likelihoods multiplied together will underflow a `f64`
/// several times over; rescaling keeps the product in range with **one logarithm per node
/// instead of one per sample**, which is a sixty-fold saving in the innermost loop.
const RESCALE: f64 = 1e150;
const LN_RESCALE: f64 = 345.398_899_014_487; // ln(1e150)

/// One position: its likelihood, and its contribution to every accumulated count.
#[allow(
    clippy::needless_range_loop,
    reason = "the sample index addresses four parallel arrays at once — the evidence, the \
              scratch, the read-group map and the accumulated counts — and zipping them would \
              hide which of the four the loop is really walking"
)]
fn one_position(
    scratch: &mut Scratch,
    group_index: &[Vec<usize>],
    ploidy: Ploidy,
    quadrature: &BetaQuadrature,
    parameters: &Parameters,
    statistics: &mut Statistics,
) {
    let samples = scratch.samples;
    let nodes = scratch.nodes;
    statistics.positions += 1.0;
    for s in 0..samples {
        let depth = scratch.evidence.samples[s].depth;
        if depth >= 1.0 {
            statistics.with_reads[s] += 1;
        }
        if depth >= 2.0 {
            statistics.with_two_reads[s] += 1;
        }
    }

    // ---- which alleles the position sums over ---------------------------------------------
    //
    // The three non-reference bases, with the ones no sample showed a read on folded into one
    // term counted as many times as there are of them: they give the identical likelihood.
    scratch.candidates.clear();
    scratch.multiplicity.clear();
    for &code in &scratch.evidence.observed_alternatives {
        scratch.candidates.push(code);
        scratch.multiplicity.push(1.0);
    }
    let unobserved = CANDIDATE_ALTERNATIVES.saturating_sub(scratch.candidates.len());
    if unobserved > 0 {
        let stand_in = (0..4)
            .find(|code| !scratch.candidates.contains(code))
            .expect("with an unobserved allele left there is a base no read fell on");
        scratch.candidates.push(stand_in);
        scratch.multiplicity.push(unobserved as f64);
    }
    let candidates = scratch.candidates.len();

    // ---- the read likelihoods, once per candidate and reused across every node -------------
    for class in 0..2 {
        let mut invariant = 0.0;
        for s in 0..samples {
            let rate = class_rate(parameters, class, group_index[s][0]);
            invariant += ln_reads_given_all_reference(&scratch.evidence.samples[s], rate);
        }
        scratch.invariant_ln[class] = invariant;

        for candidate in 0..candidates {
            let allele = scratch.candidates[candidate];
            let mut fixed = 0.0;
            for s in 0..samples {
                let rate = class_rate(parameters, class, group_index[s][0]);
                let sample = &scratch.evidence.samples[s];
                let base = scratch.ell_at(class, candidate, s);
                let mut largest = f64::NEG_INFINITY;
                for j in 0..3 {
                    let value = ln_reads_given_genotype(sample, j as u8, ploidy, allele, rate);
                    scratch.ell[base + j] = value;
                    largest = largest.max(value);
                }
                fixed += scratch.ell[base + 2];
                let slot = scratch.max_at(class, candidate, s);
                scratch.lik_max[slot] = largest;
                for j in 0..3 {
                    scratch.lik[base + j] = (scratch.ell[base + j] - largest).exp();
                }
            }
            scratch.fixed_ln[class * MAX_CANDIDATES + candidate] = fixed;
        }
    }

    // ---- the integral over the position's own allele frequency ------------------------------
    for class in 0..2 {
        for candidate in 0..candidates {
            let mut offset = 0.0;
            for s in 0..samples {
                offset += scratch.lik_max[scratch.max_at(class, candidate, s)];
            }
            for node in 0..nodes {
                let mut product = 1.0_f64;
                let mut scale = 0.0_f64;
                for s in 0..samples {
                    let prior = &quadrature.priors[(node * samples + s) * 3..][..3];
                    let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                    let term = prior[0] * lik[0] + prior[1] * lik[1] + prior[2] * lik[2];
                    product *= term;
                    if product < 1.0 / RESCALE {
                        product *= RESCALE;
                        scale -= LN_RESCALE;
                    }
                }
                let slot = scratch.node_at(class, candidate, node);
                scratch.node_ln[slot] = if product <= 0.0 {
                    f64::NEG_INFINITY
                } else {
                    product.ln() + scale + offset + quadrature.ln_weights[node]
                };
            }
        }
    }

    // ---- the three branches, and the two classes --------------------------------------------
    let ln_three = (CANDIDATE_ALTERNATIVES as f64).ln();
    for class in 0..2 {
        scratch.shares.clear();
        for candidate in 0..candidates {
            scratch.shares.push(
                scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                    + scratch.multiplicity[candidate].ln(),
            );
        }
        let fixed_alt = ln_sum_exp(&scratch.shares) - ln_three;
        scratch.shares.clear();
        for candidate in 0..candidates {
            let multiplicity = scratch.multiplicity[candidate].ln();
            for node in 0..nodes {
                scratch
                    .shares
                    .push(scratch.node_ln[scratch.node_at(class, candidate, node)] + multiplicity);
            }
        }
        let segregating = ln_sum_exp(&scratch.shares) - ln_three;
        let density = &parameters.density;
        scratch.branch_ln[class * 3] =
            density.p_invariant.max(f64::MIN_POSITIVE).ln() + scratch.invariant_ln[class];
        scratch.branch_ln[class * 3 + 1] =
            density.p_fixed_alt.max(f64::MIN_POSITIVE).ln() + fixed_alt;
        scratch.branch_ln[class * 3 + 2] =
            density.p_segregating().max(f64::MIN_POSITIVE).ln() + segregating;
        let share = if class == 0 {
            1.0 - parameters.noisy_share
        } else {
            parameters.noisy_share
        };
        scratch.class_ln[class] =
            share.max(f64::MIN_POSITIVE).ln() + ln_sum_exp(&scratch.branch_ln[class * 3..][..3]);
    }
    let position_ln = ln_sum_exp(&scratch.class_ln);
    if !position_ln.is_finite() {
        return;
    }
    statistics.log_likelihood += position_ln;

    // ---- attribute it -------------------------------------------------------------------------
    for class in 0..2 {
        let class_posterior = (scratch.class_ln[class] - position_ln).exp();
        if class_posterior <= 1e-12 {
            continue;
        }
        if class == 1 {
            statistics.noisy += class_posterior;
        }
        let branches = &scratch.branch_ln[class * 3..][..3];
        let within = ln_sum_exp(branches);
        let branch = [
            class_posterior * (branches[0] - within).exp(),
            class_posterior * (branches[1] - within).exp(),
            class_posterior * (branches[2] - within).exp(),
        ];
        statistics.invariant += branch[0];
        statistics.fixed_alt += branch[1];
        statistics.segregating += branch[2];

        // The invariant branch: every sample is homozygous reference, and every read that is
        // not on the reference base is an error, so none of them is "the allele".
        if branch[0] > 1e-12 {
            for s in 0..samples {
                let sample = &scratch.evidence.samples[s];
                let neither = sample.non_reference() - sample.on[4];
                let reference = (sample.depth - sample.non_reference()).max(0.0);
                for &g in &group_index[s] {
                    let tally = &mut statistics.reads[g][class][0];
                    tally.neither += branch[0] * neither;
                    tally.reference += branch[0] * reference;
                }
            }
        }

        // The fixed-alternative branch: every sample carries two copies of the candidate.
        if branch[1] > 1e-12 {
            scratch.shares.clear();
            for candidate in 0..candidates {
                scratch.shares.push(
                    scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                        + scratch.multiplicity[candidate].ln(),
                );
            }
            let total = ln_sum_exp(&scratch.shares);
            for candidate in 0..candidates {
                let share = branch[1] * (scratch.shares[candidate] - total).exp();
                if share <= 1e-12 {
                    continue;
                }
                let allele = scratch.candidates[candidate];
                for s in 0..samples {
                    let counts = tally_of(&scratch.evidence.samples[s], allele);
                    for &g in &group_index[s] {
                        let tally = &mut statistics.reads[g][class][2];
                        tally.candidate += share * counts.candidate;
                        tally.neither += share * counts.neither;
                        tally.reference += share * counts.reference;
                    }
                    statistics.homozygous_alt[s] += share;
                }
            }
        }

        // The segregating branch. The genotype weights are collapsed over the nodes before the
        // read counts are touched: which reads are "the allele" depends on the candidate and
        // not on the node, so the inner loop over read groups runs once per candidate rather
        // than once per candidate per node.
        if branch[2] > 1e-12 {
            scratch.shares.clear();
            for candidate in 0..candidates {
                let multiplicity = scratch.multiplicity[candidate].ln();
                for node in 0..nodes {
                    scratch.shares.push(
                        scratch.node_ln[scratch.node_at(class, candidate, node)] + multiplicity,
                    );
                }
            }
            let total = ln_sum_exp(&scratch.shares);
            scratch.per_candidate_weight[..candidates * samples * 3].fill(0.0);
            scratch.genotype_weight.fill(0.0);
            for candidate in 0..candidates {
                for node in 0..nodes {
                    let share =
                        branch[2] * (scratch.shares[candidate * nodes + node] - total).exp();
                    if share <= 1e-12 {
                        continue;
                    }
                    statistics.sum_ln_f += share * quadrature.ln_nodes[node];
                    statistics.sum_ln_one_minus_f += share * quadrature.ln_one_minus_nodes[node];
                    for s in 0..samples {
                        let prior = &quadrature.priors[(node * samples + s) * 3..][..3];
                        let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                        let joint = [prior[0] * lik[0], prior[1] * lik[1], prior[2] * lik[2]];
                        let total = joint[0] + joint[1] + joint[2];
                        if total <= 0.0 {
                            continue;
                        }
                        let base = (candidate * samples + s) * 3;
                        for j in 0..3 {
                            let weight = share * joint[j] / total;
                            scratch.per_candidate_weight[base + j] += weight;
                            scratch.genotype_weight[s * 3 + j] += weight;
                            statistics.genotypes[s][node][j] += weight;
                        }
                    }
                }
            }
            for candidate in 0..candidates {
                let allele = scratch.candidates[candidate];
                for s in 0..samples {
                    let counts = tally_of(&scratch.evidence.samples[s], allele);
                    let base = (candidate * samples + s) * 3;
                    for j in 0..3 {
                        let weight = scratch.per_candidate_weight[base + j];
                        if weight <= 0.0 {
                            continue;
                        }
                        for &g in &group_index[s] {
                            let tally = &mut statistics.reads[g][class][j];
                            tally.candidate += weight * counts.candidate;
                            tally.neither += weight * counts.neither;
                            tally.reference += weight * counts.reference;
                        }
                    }
                }
            }
            for s in 0..samples {
                statistics.heterozygous[s] += scratch.genotype_weight[s * 3 + 1];
                statistics.homozygous_alt[s] += scratch.genotype_weight[s * 3 + 2];
            }
        }
    }
}

fn class_rate(parameters: &Parameters, class: usize, group: usize) -> f64 {
    if class == 0 {
        parameters.clean[group]
    } else {
        parameters.noisy[group]
    }
}

/// Every parameter's own maximisation, over the counts one pass accumulated. Returns the
/// largest relative move any parameter made.
fn maximisation(
    parameters: &mut Parameters,
    statistics: &Statistics,
    config: &JointFitConfig,
) -> f64 {
    let mut moved: f64 = 0.0;
    let mut note = |before: f64, after: f64| {
        let scale = before.abs().max(1e-6);
        moved = moved.max((after - before).abs() / scale);
    };

    // The share of positions that are mismapped, and the density's two masses: closed form.
    let positions = statistics.positions.max(1.0);
    let before = parameters.noisy_share;
    parameters.noisy_share = (statistics.noisy / positions).clamp(1e-6, 0.5);
    note(before, parameters.noisy_share);

    let branch_total = (statistics.invariant + statistics.fixed_alt + statistics.segregating)
        .max(f64::MIN_POSITIVE);
    let before = parameters.density.p_invariant;
    parameters.density.p_invariant = (statistics.invariant / branch_total).clamp(1e-9, 1.0 - 1e-9);
    note(before, parameters.density.p_invariant);
    let before = parameters.density.p_fixed_alt;
    parameters.density.p_fixed_alt = (statistics.fixed_alt / branch_total).clamp(1e-12, 0.5);
    note(before, parameters.density.p_fixed_alt);

    // The Beta's two shapes, from the mean of `ln f` and `ln (1 − f)` under the posterior.
    if statistics.segregating > 0.0 {
        let mean_ln_f = statistics.sum_ln_f / statistics.segregating;
        let mean_ln_one_minus = statistics.sum_ln_one_minus_f / statistics.segregating;
        let (a, b) = fit_beta_shapes(
            mean_ln_f,
            mean_ln_one_minus,
            parameters.density.a,
            parameters.density.b,
        );
        note(parameters.density.a, a);
        note(parameters.density.b, b);
        parameters.density.a = a;
        parameters.density.b = b;
    }

    // Each read group's two error rates, each over nine accumulated read counts.
    for group in 0..parameters.clean.len() {
        for (class, bounds) in [(0_usize, (1e-6, 0.2)), (1, (1e-4, 0.45))] {
            let tallies = &statistics.reads[group][class];
            let current = if class == 0 {
                parameters.clean[group]
            } else {
                parameters.noisy[group]
            };
            let fitted = maximise_error_rate(tallies, config.ploidy, bounds, current);
            note(current, fitted);
            if class == 0 {
                parameters.clean[group] = fitted;
            } else {
                parameters.noisy[group] = fitted;
            }
        }
    }
    // A class that emptied cannot be told apart from one at the other's rate, and the
    // separation is what the starting points exist to span. Keep them ordered so a swap does
    // not read as a fit.
    for group in 0..parameters.clean.len() {
        if parameters.noisy[group] < parameters.clean[group] {
            parameters.noisy.swap(group, group);
            let (clean, noisy) = (parameters.noisy[group], parameters.clean[group]);
            parameters.clean[group] = clean;
            parameters.noisy[group] = noisy;
        }
    }

    // Each sample's homozygote excess, over its expected genotype counts at each node.
    let quadrature = BetaQuadrature::new(
        parameters.density.a,
        parameters.density.b,
        config.quadrature_nodes,
    );
    // **At one sample the homozygote excess is not identified and is not moved.** Nothing
    // separates "this individual is inbred" from "the population's frequencies are what they
    // are" when there is one individual, so a fit that searched it anyway would wander without
    // converging and hand back a plausible number. `fit_jointly` marks it as not fitted.
    if statistics.genotypes.len() >= 2 {
        for (s, counts) in statistics.genotypes.iter().enumerate() {
            let fitted = maximise_hom_excess(counts, &quadrature, parameters.hom_excess[s]);
            note(parameters.hom_excess[s], fitted);
            parameters.hom_excess[s] = fitted;
        }
    }
    moved
}

/// The error rate that best explains the expected read counts, by golden-section search over
/// a concave one-dimensional objective.
///
/// **This is where the cost of the route is contained.** The objective reads nine numbers, so
/// the twenty evaluations a search costs are twenty passes over nine numbers rather than
/// twenty passes over the reads.
fn maximise_error_rate(
    tallies: &[ReadTally; 3],
    ploidy: Ploidy,
    bounds: (f64, f64),
    current: f64,
) -> f64 {
    let score = |rate: f64| {
        let mut total = 0.0;
        for (j, tally) in tallies.iter().enumerate() {
            let carried = j as f64 / f64::from(ploidy.get());
            let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
            let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
            total += count_times_ln(tally.candidate, on_candidate)
                + count_times_ln(tally.neither, rate / 3.0)
                + count_times_ln(tally.reference, on_reference);
        }
        total
    };
    let observations: f64 = tallies
        .iter()
        .map(|t| t.candidate + t.neither + t.reference)
        .sum();
    if observations <= 0.0 {
        return current;
    }
    golden_section(&score, bounds.0, bounds.1)
}

/// The homozygote excess that best explains one sample's expected genotype counts.
fn maximise_hom_excess(counts: &[[f64; 3]], quadrature: &BetaQuadrature, current: f64) -> f64 {
    let observations: f64 = counts.iter().flatten().sum();
    if observations <= 0.0 {
        return current;
    }
    let score = |excess: f64| {
        let mut total = 0.0;
        for (node, row) in counts.iter().enumerate() {
            let prior = genotype_frequencies(quadrature.nodes[node], excess);
            for (j, count) in row.iter().enumerate() {
                total += count_times_ln(*count, prior[j]);
            }
        }
        total
    };
    golden_section(&score, 0.0, 1.0)
}

/// The maximum of a unimodal function on `[low, high]`, to about six figures.
fn golden_section(score: &dyn Fn(f64) -> f64, low: f64, high: f64) -> f64 {
    const INVERSE_PHI: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (low, high);
    let mut c = high - INVERSE_PHI * (high - low);
    let mut d = low + INVERSE_PHI * (high - low);
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..60 {
        if high - low < 1e-9 * (1.0 + high.abs()) {
            break;
        }
        if fc > fd {
            high = d;
            d = c;
            fd = fc;
            c = high - INVERSE_PHI * (high - low);
            fc = score(c);
        } else {
            low = c;
            c = d;
            fc = fd;
            d = low + INVERSE_PHI * (high - low);
            fd = score(d);
        }
    }
    0.5 * (low + high)
}

/// The Beta whose `E[ln f]` and `E[ln (1 − f)]` match the two accumulated means.
///
/// The two equations are `ψ(a) − ψ(a+b) = mean ln f` and `ψ(b) − ψ(a+b) = mean ln (1−f)`;
/// they are solved by a damped Newton step from the current shapes, which is enough because
/// one expectation-maximisation pass never moves them far.
fn fit_beta_shapes(mean_ln_f: f64, mean_ln_one_minus: f64, a0: f64, b0: f64) -> (f64, f64) {
    let (mut a, mut b) = (a0, b0);
    for _ in 0..50 {
        let ab = a + b;
        let g1 = digamma(a) - digamma(ab) - mean_ln_f;
        let g2 = digamma(b) - digamma(ab) - mean_ln_one_minus;
        let h11 = trigamma(a) - trigamma(ab);
        let h22 = trigamma(b) - trigamma(ab);
        let h12 = -trigamma(ab);
        let determinant = h11 * h22 - h12 * h12;
        if determinant.abs() < 1e-12 {
            break;
        }
        let step_a = (g1 * h22 - g2 * h12) / determinant;
        let step_b = (g2 * h11 - g1 * h12) / determinant;
        let next_a = (a - step_a).clamp(0.02, 50.0);
        let next_b = (b - step_b).clamp(0.02, 50.0);
        let moved = (next_a - a).abs() + (next_b - b).abs();
        a = next_a;
        b = next_b;
        if moved < 1e-10 {
            break;
        }
    }
    (a, b)
}

// ---------------------------------------------------------------------
// The integral over a position's allele frequency
// ---------------------------------------------------------------------

/// A Gauss–Jacobi rule for `∫₀¹ Beta(f; a, b) · h(f) df`, normalised so the weights sum to one.
///
/// **Gauss–Jacobi rather than a grid**, because a Beta with a shape below one — which is what
/// a population with many rare alleles has — is unbounded at zero, and a uniform grid puts no
/// node where nearly all the mass is. The rule absorbs the singularity into its weight
/// function, so the nodes are placed by the density itself.
struct BetaQuadrature {
    nodes: Vec<f64>,
    /// The rule's weights, summing to one. **The pass reads `ln_weights`**; these are what the
    /// tests integrate a known Beta with, which is the only check that the rule is the right
    /// rule rather than merely a consistent one.
    #[cfg_attr(
        not(test),
        allow(dead_code, reason = "read by this module's own tests")
    )]
    weights: Vec<f64>,
    ln_weights: Vec<f64>,
    ln_nodes: Vec<f64>,
    ln_one_minus_nodes: Vec<f64>,
    /// `[node][sample][genotype]` — how common each genotype is at that node's frequency for
    /// that sample's own inbreeding. **Computed once for the whole pass**, because it depends
    /// on nothing that varies from position to position, and it sits in the innermost loop.
    priors: Vec<f64>,
}

impl BetaQuadrature {
    /// The rule, plus the genotype priors every position will read.
    fn with_genotype_priors(a: f64, b: f64, count: usize, hom_excess: &[f64]) -> Self {
        let mut rule = Self::new(a, b, count);
        rule.priors = Vec::with_capacity(rule.nodes.len() * hom_excess.len() * 3);
        for &f in &rule.nodes {
            for &excess in hom_excess {
                rule.priors
                    .extend_from_slice(&genotype_frequencies(f, excess));
            }
        }
        rule
    }

    fn new(a: f64, b: f64, count: usize) -> Self {
        // Mapping `f = (1 + x)/2` turns the Beta's weight into the Jacobi weight
        // `(1 − x)^(b−1) (1 + x)^(a−1)`.
        let (alpha, beta) = (b - 1.0, a - 1.0);
        let (x, w) = gauss_jacobi(alpha, beta, count);
        let nodes: Vec<f64> = x
            .iter()
            .map(|x| (0.5 * (1.0 + x)).clamp(1e-12, 1.0 - 1e-12))
            .collect();
        let total: f64 = w.iter().sum();
        let weights: Vec<f64> = w.iter().map(|w| w / total).collect();
        let ln_nodes = nodes.iter().map(|f| f.ln()).collect();
        let ln_one_minus_nodes = nodes.iter().map(|f| (1.0 - f).ln()).collect();
        let ln_weights = weights
            .iter()
            .map(|w| w.max(f64::MIN_POSITIVE).ln())
            .collect();
        Self {
            nodes,
            weights,
            ln_weights,
            ln_nodes,
            ln_one_minus_nodes,
            priors: Vec::new(),
        }
    }
}

/// Nodes and weights of the `count`-point Gauss–Jacobi rule on `[-1, 1]` with weight
/// `(1 − x)^α (1 + x)^β`, by the Golub–Welsch construction: the rule's nodes are the
/// eigenvalues of the recurrence's tridiagonal matrix and its weights come from the first
/// component of each eigenvector.
fn gauss_jacobi(alpha: f64, beta: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
    let n = count.max(2);
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        let k = i as f64;
        let s = 2.0 * k + alpha + beta;
        let diagonal = if s.abs() < 1e-12 || (s + 2.0).abs() < 1e-12 {
            (beta - alpha) / (alpha + beta + 2.0)
        } else {
            (beta * beta - alpha * alpha) / (s * (s + 2.0))
        };
        matrix[i * n + i] = diagonal;
        if i + 1 < n {
            let k = (i + 1) as f64;
            let s = 2.0 * k + alpha + beta;
            let numerator = 4.0 * k * (k + alpha) * (k + beta) * (k + alpha + beta);
            let denominator = s * s * (s + 1.0) * (s - 1.0);
            let off = if denominator.abs() < 1e-12 {
                0.0
            } else {
                (numerator / denominator).max(0.0).sqrt()
            };
            matrix[i * n + i + 1] = off;
            matrix[(i + 1) * n + i] = off;
        }
    }
    let (values, vectors) = symmetric_eigen(&matrix, n);
    // The rule's weight is the **first component** of each eigenvector, squared — and
    // `vectors` is laid out eigenvector-major, so that component is at `k * n` and not at `k`.
    let weights: Vec<f64> = (0..n).map(|k| vectors[k * n] * vectors[k * n]).collect();
    (values, weights)
}

/// Eigenvalues and eigenvectors of a small symmetric matrix, by cyclic Jacobi rotations.
///
/// `vectors[k * n + s]` is component `s` of the eigenvector for eigenvalue `k`; the rule above
/// reads only the first component of each, which is why they are returned that way round.
fn symmetric_eigen(matrix: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = matrix.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum();
        if off < 1e-22 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p * n + q].abs() < 1e-18 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        a[i * n + i]
            .partial_cmp(&a[j * n + j])
            .expect("no NaN in a Jacobi recurrence")
    });
    let values = order.iter().map(|&i| a[i * n + i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (k, &i) in order.iter().enumerate() {
        for s in 0..n {
            vectors[k * n + s] = v[s * n + i];
        }
    }
    (values, vectors)
}

/// `ψ(x)`, the derivative of `ln Γ(x)`, by recurrence up to eight and then an asymptotic
/// series.
fn digamma(x: f64) -> f64 {
    let mut x = x;
    let mut result = 0.0;
    while x < 8.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    result + x.ln()
        - 0.5 * inverse
        - square * (1.0 / 12.0 - square * (1.0 / 120.0 - square / 252.0))
}

/// `ψ'(x)`, by the same shape.
fn trigamma(x: f64) -> f64 {
    let mut x = x;
    let mut result = 0.0;
    while x < 8.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    result + inverse * (1.0 + 0.5 * inverse + square * (1.0 / 6.0 - square / 30.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- the pieces, each on its own ------------------------------------------

    #[test]
    fn the_quadrature_integrates_a_beta_to_one_and_gets_its_mean_right() {
        for (a, b) in [(0.3, 2.0), (1.0, 1.0), (2.5, 5.0), (0.1, 0.9)] {
            let rule = BetaQuadrature::new(a, b, 24);
            let total: f64 = rule.weights.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-10,
                "Beta({a}, {b}) weights sum to {total}"
            );
            let mean: f64 = rule
                .nodes
                .iter()
                .zip(rule.weights.iter())
                .map(|(f, w)| f * w)
                .sum();
            let truth = a / (a + b);
            assert!(
                (mean - truth).abs() < 1e-8,
                "Beta({a}, {b}) mean {mean} against {truth}"
            );
        }
    }

    /// **The shape a neutral population has is the one a uniform grid cannot integrate**, so
    /// it is the one to check: `Beta(0.2, 1)` puts nearly all its mass below `f = 0.01`.
    #[test]
    fn the_quadrature_handles_the_rare_allele_pile_up() {
        let rule = BetaQuadrature::new(0.2, 1.0, 24);
        // E[2 f (1 − f)] under Beta(a, b) is 2ab / ((a+b)(a+b+1)).
        let heterozygosity: f64 = rule
            .nodes
            .iter()
            .zip(rule.weights.iter())
            .map(|(f, w)| w * 2.0 * f * (1.0 - f))
            .sum();
        let truth = 2.0 * 0.2 * 1.0 / (1.2 * 2.2);
        assert!(
            (heterozygosity - truth).abs() < 1e-8,
            "{heterozygosity} against {truth}"
        );
    }

    #[test]
    fn the_beta_shapes_come_back_from_their_own_log_means() {
        for (a, b) in [(0.4, 1.5), (2.0, 2.0), (0.9, 6.0)] {
            let rule = BetaQuadrature::new(a, b, 48);
            let mean_ln_f: f64 = rule
                .nodes
                .iter()
                .zip(rule.weights.iter())
                .map(|(f, w)| w * f.ln())
                .sum();
            let mean_ln_one_minus: f64 = rule
                .nodes
                .iter()
                .zip(rule.weights.iter())
                .map(|(f, w)| w * (1.0 - f).ln())
                .sum();
            let (fitted_a, fitted_b) = fit_beta_shapes(mean_ln_f, mean_ln_one_minus, 1.0, 1.0);
            assert!(
                (fitted_a - a).abs() < 0.02 && (fitted_b - b).abs() < 0.05,
                "Beta({a}, {b}) came back as Beta({fitted_a}, {fitted_b})"
            );
        }
    }

    #[test]
    fn a_read_scored_the_two_routes_ways_agrees_on_the_non_reference_chance() {
        // The chance a read shows anything but the reference base, summed out of this
        // module's three categories, must equal the histogram route's own `p_j(ε)`.
        let ploidy = Ploidy::try_new(2).expect("two");
        for rate in [0.001, 0.01, 0.1] {
            for j in 0..=2_u8 {
                let carried = f64::from(j) / 2.0;
                let theirs = carried * (1.0 - rate / 3.0) + (1.0 - carried) * rate;
                let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
                let ours = on_candidate + 2.0 * rate / 3.0;
                assert!(
                    (ours - theirs).abs() < 1e-12,
                    "at {j} copies and rate {rate}: {ours} against {theirs}"
                );
                let _ = ploidy;
            }
        }
    }

    #[test]
    fn genotype_frequencies_sum_to_one_and_the_excess_moves_the_heterozygote() {
        for f in [0.01, 0.3, 0.5, 0.9] {
            for excess in [0.0, 0.3, 1.0] {
                let g = genotype_frequencies(f, excess);
                let total: f64 = g.iter().sum();
                assert!((total - 1.0).abs() < 1e-12, "{g:?} sums to {total}");
            }
            let none = genotype_frequencies(f, 0.0);
            let whole = genotype_frequencies(f, 1.0);
            assert!(whole[1] < none[1] || none[1] == 0.0);
            assert!(
                (whole[1] - 0.0).abs() < 1e-12,
                "a fully inbred plant has no heterozygotes"
            );
        }
    }

    #[test]
    fn the_expected_heterozygosity_is_the_densitys_own() {
        let density = FrequencyDensity {
            p_invariant: 0.9,
            p_fixed_alt: 0.01,
            a: 0.5,
            b: 2.0,
        };
        let truth = 0.09 * 2.0 * 0.5 * 2.0 / (2.5 * 3.5);
        assert!((density.expected_heterozygosity() - truth).abs() < 1e-12);
    }

    // ---- the whole fit, against a cohort whose truth is known --------------------
    //
    // **The oracle is a drawn cohort and not a fixture** (spec §12.1): records are filled
    // from parameters chosen here, with no reads and no alignments, and the fit has to return
    // what was drawn.

    use crate::ng::parameter_estimation::joint::loci::{
        CatalogBuildSettings, KeptLociDigester, ReferenceDigest, RegionSetDigest, SelectionIdentity,
    };
    use crate::ng::parameter_estimation::joint::records::{
        AlleleObservation, DepthCap, DepthLadderDigest, GenericRecords, ObservedAllele,
        PackedDepthCodes, ReadCap, RecordIdentity,
    };
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::tandem_repeat::ScanParams;

    fn selection_identity() -> SelectionIdentity {
        SelectionIdentity {
            seed: 42,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: 2_000_000,
            ssr_cap: 1_000,
        }
    }

    /// The stream every drawn number comes from. Deterministic, so a failure is reproducible.
    struct Draw(u64);

    impl Draw {
        fn uniform(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }

        fn pick(&mut self, weights: &[f64]) -> usize {
            let total: f64 = weights.iter().sum();
            let mut cut = self.uniform() * total;
            for (index, weight) in weights.iter().enumerate() {
                cut -= weight;
                if cut <= 0.0 {
                    return index;
                }
            }
            weights.len() - 1
        }

        /// A Beta draw, by the ratio of two Gammas, each by Marsaglia–Tsang.
        fn beta(&mut self, a: f64, b: f64) -> f64 {
            let x = self.gamma(a);
            let y = self.gamma(b);
            (x / (x + y)).clamp(1e-9, 1.0 - 1e-9)
        }

        fn gamma(&mut self, shape: f64) -> f64 {
            if shape < 1.0 {
                let u = self.uniform().max(1e-12);
                return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let z = self.normal();
                let v = (1.0 + c * z).powi(3);
                if v <= 0.0 {
                    continue;
                }
                let u = self.uniform().max(1e-12);
                if u.ln() < 0.5 * z * z + d - d * v + d * (v.ln()) {
                    return d * v;
                }
            }
        }

        fn normal(&mut self) -> f64 {
            let u1 = self.uniform().max(1e-12);
            let u2 = self.uniform();
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }

        fn poisson(&mut self, mean: f64) -> u32 {
            let limit = (-mean).exp();
            let mut product = self.uniform();
            let mut count = 0;
            while product > limit && count < 200 {
                count += 1;
                product *= self.uniform();
            }
            count
        }
    }

    struct DrawnCohort {
        samples: Vec<SampleRecords>,
        clean: f64,
        noisy: f64,
        noisy_share: f64,
        density: FrequencyDensity,
        hom_excess: Vec<f64>,
        heterozygous: Vec<f64>,
    }

    /// Draw a cohort at known parameters and write it into records the fit will read.
    fn draw_cohort(
        samples: usize,
        positions: usize,
        mean_depth: f64,
        truth: (f64, f64, f64),
        density: FrequencyDensity,
        hom_excess: f64,
        seed: u64,
    ) -> DrawnCohort {
        let (clean, noisy, noisy_share) = truth;
        let mut draw = Draw(seed);
        let excess = vec![hom_excess; samples];
        let mut depth: Vec<PackedDepthCodes> = (0..samples)
            .map(|_| PackedDepthCodes::never_walked(positions))
            .collect();
        let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
        let mut heterozygous = vec![0_u64; samples];
        let edges = DepthBinEdges::new();

        for index in 0..positions {
            let rate = if draw.uniform() < noisy_share {
                noisy
            } else {
                clean
            };
            let branch = draw.pick(&[
                density.p_invariant,
                density.p_fixed_alt,
                density.p_segregating(),
            ]);
            // The reference base is code 0 by construction here, so the three candidates are
            // codes 1..=3 and the fit's own sum over them is being asked to find the right one.
            let allele = 1 + (draw.uniform() * 3.0) as usize % 3;
            let f = match branch {
                0 => 0.0,
                1 => 1.0,
                _ => draw.beta(density.a, density.b),
            };
            for sample in 0..samples {
                let reads = draw.poisson(mean_depth);
                let genotype = match branch {
                    0 => 0_usize,
                    1 => 2,
                    _ => draw.pick(&genotype_frequencies(f, excess[sample])),
                };
                if genotype == 1 {
                    heterozygous[sample] += 1;
                }
                let carried = genotype as f64 / 2.0;
                let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
                let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
                let mut counts = [0_u32; 5];
                for _ in 0..reads {
                    let u = draw.uniform();
                    let code = if u < on_candidate {
                        allele
                    } else if u < on_candidate + on_reference {
                        0
                    } else if u < on_candidate + on_reference + rate / 3.0 {
                        (allele % 3) + 1
                    } else {
                        ((allele + 1) % 3) + 1
                    };
                    counts[code] += 1;
                }
                depth[sample].set(index, DepthCode::Binned(edges.bin_for(reads)));
                for (code, count) in counts.iter().enumerate() {
                    if code == 0 || *count == 0 {
                        continue;
                    }
                    sparse[sample].push(AlleleObservation {
                        index: index as u32,
                        allele: match code {
                            1 => ObservedAllele::C,
                            2 => ObservedAllele::G,
                            _ => ObservedAllele::T,
                        },
                        reads: *count,
                    });
                }
            }
        }

        let identity = RecordIdentity {
            selection: selection_identity(),
            kept_loci: KeptLociDigester::new().finish(),
            ssr_stratum_counts: Default::default(),
            read_cap: ReadCap(100),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::new()),
            depth_cap: DepthCap(124),
            coverage_window: None,
        };
        let records = (0..samples)
            .map(|s| SampleRecords {
                sample: format!("s{s}"),
                generic: [(
                    ReadGroupId(0),
                    GenericRecords::from_parts(
                        std::mem::replace(&mut depth[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    ),
                )]
                .into_iter()
                .collect(),
                ssr: BTreeMap::new(),
                coverage: None,
                identity: identity.clone(),
            })
            .collect();
        DrawnCohort {
            samples: records,
            clean,
            noisy,
            noisy_share,
            density,
            hom_excess: excess,
            heterozygous: heterozygous
                .iter()
                .map(|n| *n as f64 / positions as f64)
                .collect(),
        }
    }

    /// **The test that says whether any of this works.** A cohort drawn at known parameters
    /// must come back at them.
    #[test]
    fn a_drawn_cohort_comes_back_at_the_parameters_it_was_drawn_at() {
        let density = FrequencyDensity {
            p_invariant: 0.90,
            p_fixed_alt: 0.01,
            a: 0.7,
            b: 2.5,
        };
        let cohort = draw_cohort(
            10,
            3_000,
            8.0,
            (0.002, 0.06, 0.02),
            density,
            0.2,
            0x9E37_79B9_7F4A_7C15,
        );
        let config = JointFitConfig {
            quadrature_nodes: 12,
            max_passes: 40,
            starting_points: vec![StartingPoint {
                clean: 0.006,
                noisy: 0.12,
                noisy_share: 0.05,
                p_invariant: 0.8,
                p_fixed_alt: 0.02,
                a: 1.0,
                b: 1.0,
            }],
            ..JointFitConfig::default()
        };
        let fit = fit_jointly(&cohort.samples, &config).expect("the cohort pools");

        eprintln!(
            "clean {:.5} (drawn {:.5})  noisy {:.4} (drawn {:.4})  share {:.4} (drawn {:.4})",
            fit.noise[&ReadGroupId(0)].value.clean,
            cohort.clean,
            fit.noise[&ReadGroupId(0)].value.noisy,
            cohort.noisy,
            fit.noisy_share,
            cohort.noisy_share
        );
        eprintln!(
            "p_inv {:.4} (drawn {:.4})  p_fixed {:.5} (drawn {:.5})  a {:.3} (drawn {:.3})  b {:.3} (drawn {:.3})",
            fit.density.value.p_invariant,
            cohort.density.p_invariant,
            fit.density.value.p_fixed_alt,
            cohort.density.p_fixed_alt,
            fit.density.value.a,
            cohort.density.a,
            fit.density.value.b,
            cohort.density.b
        );
        eprintln!(
            "Hexp {:.5}; F {:.3} (drawn {:.3}); Hobs {:.5} (drawn {:.5}); passes {} converged {}",
            fit.expected_heterozygosity,
            fit.hom_excess["s0"].value.get(),
            cohort.hom_excess[0],
            fit.rates["s0"].value.heterozygous,
            cohort.heterozygous[0],
            fit.passes,
            fit.converged
        );
        let clean = fit.noise[&ReadGroupId(0)].value.clean;
        assert!(
            (clean / cohort.clean - 1.0).abs() < 0.20,
            "clean error rate {clean} against {}",
            cohort.clean
        );
        let invariant = fit.density.value.p_invariant;
        assert!(
            (invariant - cohort.density.p_invariant).abs() < 0.02,
            "the invariant share came back {invariant} against {}",
            cohort.density.p_invariant
        );
        let hexp = fit.expected_heterozygosity;
        let drawn: f64 = cohort.heterozygous.iter().sum::<f64>() / cohort.heterozygous.len() as f64;
        // `Hexp` is the population's, and the drawn samples are inbred by 0.2, so the two are
        // related by that factor rather than equal.
        let predicted = hexp * (1.0 - cohort.hom_excess[0]);
        assert!(
            (predicted / drawn - 1.0).abs() < 0.20,
            "predicted heterozygosity {predicted} against the drawn {drawn}"
        );
        let excess = fit.hom_excess["s0"].value.get();
        assert!(
            (excess - cohort.hom_excess[0]).abs() < 0.05,
            "homozygote excess {excess} against the drawn {}",
            cohort.hom_excess[0]
        );
        // **The two classes must stay apart.** A fit that collapsed them would report
        // convergence with one emptied, which is the failure the starting points exist to
        // avoid, and every number above would still look reasonable.
        let noisy = fit.noise[&ReadGroupId(0)].value.noisy;
        assert!(
            noisy > 5.0 * clean && (noisy / cohort.noisy - 1.0).abs() < 0.35,
            "noisy rate {noisy} against the drawn {} and the clean {clean}",
            cohort.noisy
        );
        assert!(
            (fit.noisy_share / cohort.noisy_share - 1.0).abs() < 1.0,
            "noisy share {} against the drawn {}",
            fit.noisy_share,
            cohort.noisy_share
        );
        assert!(fit.converged, "the alternation ran out of passes");
    }

    /// **The control that says the last test measured something.** A cohort drawn with no
    /// inbreeding at all must come back with none — otherwise the agreement above would be
    /// consistent with an estimator that returns whatever it was started at, and a fit that
    /// invents inbreeding books a sample's mismapping as biology.
    #[test]
    fn a_cohort_drawn_without_inbreeding_does_not_invent_it() {
        let density = FrequencyDensity {
            p_invariant: 0.90,
            p_fixed_alt: 0.01,
            a: 0.7,
            b: 2.5,
        };
        let cohort = draw_cohort(
            10,
            3_000,
            8.0,
            (0.002, 0.06, 0.02),
            density,
            0.0,
            0x51_ED_27_09,
        );
        let config = JointFitConfig {
            quadrature_nodes: 12,
            max_passes: 40,
            starting_points: vec![StartingPoint {
                clean: 0.002,
                noisy: 0.06,
                noisy_share: 0.02,
                p_invariant: 0.9,
                p_fixed_alt: 0.01,
                // Started at a *wrong* Beta on purpose: a start at the truth would make this
                // test pass for an estimator that never moved.
                a: 1.5,
                b: 1.0,
            }],
            ..JointFitConfig::default()
        };
        let fit = fit_jointly(&cohort.samples, &config).expect("the cohort pools");
        for sample in cohort.samples.iter() {
            let excess = fit.hom_excess[&sample.sample].value.get();
            assert!(
                excess < 0.08,
                "{} came back inbred by {excess} where the truth is none",
                sample.sample
            );
        }
        // **The shape has to have moved off its start**, or this test would pass for an
        // estimator that never touched it. It is the least well recovered of the fitted
        // numbers — from a start of 1.5 against a truth of 0.7 it comes back near 1.0, so
        // about two thirds of the way — and the assertion says that rather than pretending
        // to a precision the measurement does not have.
        let travelled = (1.5 - fit.density.value.a) / (1.5 - cohort.density.a);
        assert!(
            travelled > 0.5,
            "the Beta's shape moved {:.0}% of the way from its start of 1.5 to the drawn {},              ending at {}",
            100.0 * travelled,
            cohort.density.a,
            fit.density.value.a
        );
    }

    /// **The route runs at one sample and says what it cannot fit there.** Nothing separates a
    /// sample's own inbreeding from the population's frequency density when there is one
    /// sample, so the number that comes back is marked as not fitted rather than reported as a
    /// measurement.
    #[test]
    fn one_sample_fits_the_density_and_marks_what_it_could_not_fit() {
        let density = FrequencyDensity {
            p_invariant: 0.90,
            p_fixed_alt: 0.01,
            a: 0.7,
            b: 2.5,
        };
        let cohort = draw_cohort(1, 3_000, 20.0, (0.002, 0.06, 0.02), density, 0.0, 99);
        let fit = fit_jointly(
            &cohort.samples,
            &JointFitConfig {
                quadrature_nodes: 12,
                max_passes: 30,
                starting_points: vec![StartingPoint::spanning_the_class_separation()[1]],
                ..JointFitConfig::default()
            },
        )
        .expect("one sample is a cohort of one");
        assert_eq!(fit.hom_excess["s0"].provenance, Provenance::Defaulted);
        assert!(
            matches!(
                fit.contamination["s0"],
                ContaminationEstimate::NotIdentified {
                    reason: super::super::contamination::NotIdentifiedReason::NoPanel
                }
            ),
            "with one sample there is no panel to be surprised by, and saying so is the point"
        );
        let clean = fit.noise[&ReadGroupId(0)].value.clean;
        assert!(
            (clean / cohort.clean - 1.0).abs() < 0.35,
            "the error rate is still fitted at one sample: {clean} against {}",
            cohort.clean
        );
    }

    /// **The refusal, before any arithmetic.** Two samples that did not keep the same loci
    /// cannot be pooled, and nothing in the data would look wrong.
    #[test]
    fn samples_that_disagree_on_the_ladder_are_refused_and_the_field_is_named() {
        let density = FrequencyDensity {
            p_invariant: 0.95,
            p_fixed_alt: 0.005,
            a: 0.7,
            b: 2.0,
        };
        let mut cohort = draw_cohort(3, 50, 4.0, (0.002, 0.05, 0.01), density, 0.0, 7);
        cohort.samples[1].identity.depth_cap = DepthCap(60);
        let error = fit_jointly(&cohort.samples, &JointFitConfig::default())
            .expect_err("the samples did not record the same evidence");
        match error {
            JointFitError::IdentityMismatch { field, .. } => {
                assert_eq!(field, "per-position depth cap")
            }
            other => panic!("{other}"),
        }
    }

    #[test]
    fn the_homozygote_excess_refuses_a_heterozygote_excess() {
        assert!(HomozygoteExcess::try_new(-0.1).is_none());
        assert!(HomozygoteExcess::try_new(0.0).is_some());
        assert!(HomozygoteExcess::try_new(1.0).is_some());
        assert!(HomozygoteExcess::try_new(1.2).is_none());
    }
}
