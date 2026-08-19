//! The estimator: every parameter fitted once, over every sample at the same loci.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_fit.md`. Types:
//! `doc/devel/ng/arch/parameter_prepass_joint_fit.md`. It reads the evidence of
//! [`census`](super::census) at the loci of [`loci`](super::loci) and nothing else.
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

use std::collections::BTreeMap;
use std::sync::Arc;

use rayon::prelude::*;

use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::types::{Ploidy, ReadGroupId};

use super::census::{
    CensusError, CohortCensusEvidence, DepthCap, DepthCode, SampleGenericSections,
    TermsDisagreement,
};
use super::contamination::{ContaminationConfig, ContaminationEstimate, fit_contamination_over};

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

/// Positions a **sample** carries more copies of than the reference does.
///
/// Where a plant holds two copies of a stretch the reference holds once, both copies' reads pile
/// up at the same place, and wherever the copies differ from each other about half the reads
/// disagree with the reference. The two-class model has no home for that except *heterozygous*,
/// and heterozygosity is one of the numbers this pass exists to produce: ignoring the class puts
/// observed heterozygosity **50.6% above the truth** on a fifty-sample selfing panel at three
/// reads a position, while expected heterozygosity rises only 10.6%, so the fitted homozygote
/// excess reads **0.4471 where the truth is 0.5942** — a quarter of it gone, at a value nothing
/// refuses (`duplicated_class_identification_2026-08-13.md`).
///
/// **The class is an ordinary variant with one genotype removed.** Each sample is either a
/// carrier, with about half its reads disagreeing, or homozygous reference; what a duplication
/// has no room for is a sample homozygous for the non-reference allele. A real variant at a
/// frequency of a half leaves about a quarter of the panel there and a duplication leaves none,
/// and across a cohort that absence is the whole evidence.
///
/// **It needs samples rather than reads.** Twenty-five reads a position buys the pattern nothing
/// over three; ten samples leaves 39 of every 100 carrier positions behind and fifty leaves 9.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct DuplicatedPositions {
    /// The share of kept positions drawn from the class. **Not a measurement of how much
    /// duplication a cohort carries** — it comes back about twice the truth while sorting the
    /// positions correctly, so what carries the quantity is the per-position posterior.
    pub share: f64,
    /// The Beta the share of the panel carrying a given duplication is integrated over.
    pub carrier_a: f64,
    pub carrier_b: f64,
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
    ///
    /// **A position the sample carries an extra copy of is not counted here**, which is the
    /// whole reason [`DuplicatedPositions`] exists.
    pub heterozygous: f64,
    /// …and that every copy is non-reference.
    pub homozygous_alt: f64,
    /// …and that the sample carries more copies of the position than the reference does. Zero
    /// when the class is not fitted.
    pub duplicated_carrier: f64,
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
    /// Positions a sample carries more copies of than the reference. `None` when the run did
    /// not fit the class.
    pub duplicated: Option<Estimate<DuplicatedPositions>>,
    /// Per sample, the departure from the Hardy–Weinberg proportions the density predicts.
    pub hom_excess: BTreeMap<String, Estimate<HomozygoteExcess>>,
    /// Per sample, derived from the converged posteriors rather than fitted.
    pub rates: BTreeMap<String, Estimate<SampleGenotypeRates>>,
    /// Per sample, the fraction of reads from another individual.
    pub contamination: BTreeMap<String, ContaminationEstimate>,
    /// The population's expected heterozygosity, read off the density.
    pub expected_heterozygosity: f64,
    /// **For every kept position, in position order, the probability that it is mismapped** —
    /// that two stretches of genome the reference holds once are both piling reads up here.
    ///
    /// This is the one quantity that only exists because every sample is present at the same
    /// position: with one sample a few disagreeing reads are indistinguishable from a rare
    /// heterozygote, and with sixty-three a position that reads part non-reference *in
    /// everybody* has nowhere else to go. It is a per-position posterior and not a parameter,
    /// and it is on this value because two consumers need it —
    /// [`contamination`](super::contamination), which must not measure a stray-read fraction
    /// over positions where every sample has stray reads, and the caller downstream.
    ///
    /// Four bytes a position: 8 MB at the two-million-position budget.
    pub noisy_posterior: Vec<f32>,
    /// **For every kept position, in position order, each sample's posterior that it is
    /// heterozygous there, that both its copies are non-reference, and that it carries an
    /// extra copy of the position** — three values a sample, the samples in the order they
    /// were handed in. Empty unless the run asked for it, which is
    /// [`JointFitConfig::genotype_posteriors`].
    ///
    /// The per-sample rates are these values' means, so this is what a comparison against a
    /// benchmark VCF needs in order to ask *which* positions the two disagree at rather than
    /// only by how much.
    pub genotype_posterior: Vec<f32>,
    /// **For every kept position, in position order, the probability that the position is
    /// drawn from the duplicated class** — that some sample carries more copies of it than the
    /// reference does. Empty unless the run asked for it
    /// ([`JointFitConfig::genotype_posteriors`]), and zero at every position when the run did
    /// not fit the class.
    ///
    /// The class's own `share` is this value's mean. Kept separately because a share says how
    /// many positions the class claims and says nothing about *which*, and every check that
    /// the claimed positions are duplications rather than rare variants — their read depth,
    /// their carrier count — is a check on the list.
    ///
    /// Four bytes a position: 8 MB at the two-million-position budget.
    pub duplicated_posterior: Vec<f32>,
    /// One entry per pass of the alternation, when the run asked for it. Empty otherwise.
    pub trace: Vec<PassSummary>,
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
    ///
    /// **Raised where the cohort is built, not here — 2026-08-14.**
    /// [`CohortCensusEvidence::new`](super::census::CohortCensusEvidence::new) makes the check
    /// before a single section is read, and a caller turns its refusal into this with `?`. The
    /// variant keeps its name because
    /// [`parameter_prepass_joint_loci.md`](../../../doc/devel/ng/spec/parameter_prepass_joint_loci.md)
    /// specifies it.
    #[error("samples {first} and {second} disagree on {field}; they did not keep the same loci")]
    IdentityMismatch {
        first: String,
        second: String,
        field: &'static str,
    },
    /// A sample whose census is a file the fit could not read. **The estimator's own failure to
    /// obtain evidence, not a property of the evidence** — see [`CensusError`].
    #[error("a sample's census could not be read")]
    Census(#[from] CensusError),
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
    /// Where the duplicated class starts, when the run fits one. The tomato measurement puts
    /// about one duplicated position in every thousand, carried by a tenth of the panel on
    /// average — `Beta(1.2, 9.5)` — so that is where the search begins.
    pub duplicated_share: f64,
    pub carrier_a: f64,
    pub carrier_b: f64,
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
                duplicated_share: 0.001,
                carrier_a: 1.2,
                carrier_b: 9.5,
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
    /// hold codes that mean different depths, which the recording-terms check already refuses.
    pub edges: Arc<DepthBinEdges>,
    /// How the share of a sample's reads that came from another individual is measured
    /// ([`contamination`](super::contamination)), which the fit runs once it has converged.
    pub contamination: ContaminationConfig,
    /// Whether positions a sample carries more copies of than the reference get a class of
    /// their own ([`DuplicatedPositions`]).
    ///
    /// **On by default, and the cost of turning it off is a quarter of the homozygote excess**
    /// on a selfing panel — at a value nothing refuses. It costs about one further pass's worth
    /// of arithmetic per pass, because the class integrates over a second frequency.
    pub duplicated_positions: bool,
    /// How much more likely each sample's reads around each kept position are if the sample has
    /// two copies of it rather than one — `ln P(coverage | two) − ln P(coverage | one)`, per
    /// sample, indexed by kept position.
    ///
    /// **The second discriminator, and the only one that works at one sample.** The pattern
    /// across a cohort needs about twenty-five samples before the absence of non-reference
    /// homozygotes means anything; a window's read depth says the same thing whoever else is in
    /// the run. Empty — the default — leaves the class identified by the cohort alone.
    pub coverage_odds: Vec<Arc<[f32]>>,
    /// Whether to keep, for each sample at each kept position, the posterior that it is
    /// heterozygous there, that both its copies are non-reference, and that it carries an
    /// extra copy of the position — and, once per position, the probability that the position
    /// belongs to the duplicated class at all.
    ///
    /// **Off by default because of what it weighs**: twelve bytes a position a sample, which
    /// is 16 MB for three samples over the benchmark trio's positions and 1.5 GB for fifty
    /// samples over two million. What it is for is checking the fitted rates against a truth
    /// set position by position rather than as one mean.
    pub genotype_posteriors: bool,
    /// Whether to keep one line per pass of the alternation ([`PassSummary`]). Costs nothing
    /// but a few hundred bytes; off by default because a converged fit has nothing to say
    /// with it.
    pub pass_trace: bool,
    /// Whether a stored depth code is read as the range it stands for, or as the middle of
    /// that range.
    ///
    /// **The range is what ships and the middle is kept only so the two can be measured on the
    /// same reads.** Reading the middle puts a heterozygote's read share away from a half for a
    /// reason that has nothing to do with the sample, and what that costs grows with depth:
    /// on a drawn cohort at thirty reads a position it returns heterozygosity 4 to 7% above
    /// what the cohort was drawn with, where the range returns it within 2%.
    pub depth_as_a_range: bool,
}

impl Default for JointFitConfig {
    fn default() -> Self {
        Self {
            ploidy: Ploidy::try_new(2).expect("two is a ploidy"),
            quadrature_nodes: 16,
            starting_points: StartingPoint::spanning_the_class_separation(),
            max_passes: 200,
            stillness: 1e-4,
            // **The census ladder and not the histogram one.** This fit reads census codes and
            // nothing else, and since 2026-08-16 the two ladders differ above 8 reads a
            // position rather than only above 124 — so the histogram ladder here would not be
            // a coarser reading of the same codes, it would be the wrong depths.
            edges: Arc::new(DepthBinEdges::for_census()),
            contamination: ContaminationConfig::default(),
            duplicated_positions: true,
            coverage_odds: Vec::new(),
            genotype_posteriors: false,
            pass_trace: false,
            depth_as_a_range: true,
        }
    }
}

// ---------------------------------------------------------------------
// The evidence, laid out the way the pass over it wants
// ---------------------------------------------------------------------

/// One sample's reads at one position, as the likelihood reads them.
///
/// Counts, not codes: `on[k]` is how many reads showed base `k`. **The reference base is not
/// here and is not needed** — the records list only what disagreed with it, so reads on the
/// reference are what is left over.
///
/// # The depth is a range, and the reads that disagree are exact
///
/// A position's read count is stored as a code on the census ladder. **To the cap of 124 reads
/// each code is one exact depth**; above it a code stands for a range that widens as it climbs,
/// to about 1,500. So a range is now what a position deeper than the cap has — where its allele
/// counts have been thinned and the record is approximate however the depth is written —
/// and every position below the cap has a depth, not a range.
///
/// The range still has to be summed over rather than read at one depth, because **the count of
/// reads that disagreed with the reference is exact and a range-valued depth is not**, and the
/// reference count is the difference between them: taking one depth out of the range lands a
/// heterozygote's read share away from a half for a reason that has nothing to do with the
/// sample. What that cost while the ladder widened from nine reads upwards was measured twice —
/// a drawn panel with nobody contaminated read a median 2.5% contaminated at ten reads a
/// position from this alone (`contamination_floor_and_duplicated_class_2026-08-13.md` §4), and
/// a genuinely 3%-contaminated sample at 30 reads read 0.0120 against the 0.0263 an exact
/// ladder returns (`census_depth_resolution_2026-08-16.md`).
///
/// So the likelihood **sums over every depth the code could stand for**, giving each equal
/// weight. Below the cap that range is one value and this is the plain multinomial it always
/// was, which is where every cohort this caller has been run on spends nearly all its
/// positions.
#[derive(Copy, Clone, Default, Debug)]
struct SampleAtPosition {
    /// The mean of the depths the code could stand for. **What the read tallies the error rate
    /// is maximised over are attributed at, and not what the likelihood uses** — that sums over
    /// the range.
    depth: f64,
    /// Reads on each of the four bases that are **not** the reference, indexed by allele code;
    /// the reference base's own entry is always zero, and so is any base no read showed.
    on: [f64; 5],
    /// The fewest reference reads the stored code allows, and how many counts the range holds
    /// counting from there. A range that would demand fewer reference reads than zero is cut
    /// at zero: the reads are the harder evidence, so the depth gives way.
    fewest_reference: f64,
    spread: usize,
}

impl SampleAtPosition {
    /// Reads that disagreed with the reference in any way.
    fn non_reference(&self) -> f64 {
        self.on.iter().sum()
    }
}

/// The widest range of depths one stored code can stand for, plus room for a ladder someone
/// later widens. The adopted ladder's widest recorded range is 76–97, which is twenty-two.
const MAX_RECORDED_SPREAD: usize = 32;

/// A cohort refused at the door, as the fit's own error.
///
/// **The check belongs to the census and the name belongs to the fit.**
/// [`CohortCensusEvidence::new`](super::census::CohortCensusEvidence::new) compares the twelve
/// recording terms before a section is decoded; the variant a caller sees is the one
/// [`parameter_prepass_joint_loci.md`](../../../doc/devel/ng/spec/parameter_prepass_joint_loci.md)
/// specifies.
impl From<TermsDisagreement> for JointFitError {
    fn from(refusal: TermsDisagreement) -> Self {
        Self::IdentityMismatch {
            first: refusal.first,
            second: refusal.second,
            field: refusal.field,
        }
    }
}

/// The cohort at one position.
struct PositionEvidence {
    /// One entry per sample, in the order the fit iterates samples.
    samples: Vec<SampleAtPosition>,
    /// Per sample, `MAX_RECORDED_SPREAD` slots holding how much weight each depth in that
    /// sample's range carries, relative to the shallowest — see
    /// [`ln_reference_reads`]. Flat rather than per sample so that one position's evidence is
    /// one allocation reused over two million positions.
    depth_weights: Vec<f64>,
    /// Which non-reference bases any sample showed a read on. **The candidates the fit sums
    /// over are the three bases that are not the reference**, and the ones nobody showed a
    /// read on all give the identical term, so they are counted rather than enumerated.
    observed_alternatives: Vec<usize>,
}

impl PositionEvidence {
    fn weights_of(&self, sample: usize) -> &[f64] {
        &self.depth_weights[sample * MAX_RECORDED_SPREAD..][..MAX_RECORDED_SPREAD]
    }
}

/// Walks every sample's records once, in position order, handing one position at a time.
///
/// **A cursor per sample rather than a search per position.** Each sample's non-reference
/// observations are sorted by position, so advancing a cursor costs one comparison; a binary
/// search per sample per position would cost twenty-one, two million times over.
struct EvidenceCursor<'a> {
    /// Per sample, the ordinary-position sections the census lent for this call — **borrowed
    /// and never held**: the cursor lives inside the scoped call that produced them.
    samples: &'a [SampleGenericSections<'a>],
    edges: &'a DepthBinEdges,
    /// The cap the allele counts were thinned to. **Held beside the ladder because the two
    /// answer different questions**: the ladder says what depth a stored code stands for, and
    /// this puts that depth into the units the counts beside it were taken in.
    cap: DepthCap,
    /// Each sample's own mean read depth over the kept positions — the centre of the Poisson
    /// that [`fill_depth_weights`] weights a stored code's range with.
    coverage: &'a [f64],
    /// False restores the point-read a stored code used to be given: the middle of its range,
    /// one number, with no sum over the depths it stands for.
    as_a_range: bool,
    /// Per sample, per read group, where that group's sparse list has been read to.
    at: Vec<Vec<usize>>,
    positions: usize,
    next: usize,
}

impl<'a> EvidenceCursor<'a> {
    fn position_count(samples: &[SampleGenericSections<'_>]) -> usize {
        samples
            .first()
            .and_then(|sections| sections.first())
            .map_or(0, |(_, records)| records.depth().len())
    }

    /// Each sample's mean read depth over the positions it was walked at.
    ///
    /// **One number per sample for the whole run, and coverage is not one number.** A sample
    /// reads deeper in some parts of a genome than others, and this cannot see that: it is the
    /// centre of the prior a stored code's range is weighted with, so where a position's own
    /// coverage is far from the sample's mean, the range is weighted by the wrong Poisson.
    /// What it costs is bounded by the width of a bin, which is why a single number is enough
    /// to be going on with — the per-window coverage summary the records already carry is
    /// where a better one would come from.
    fn mean_depth(
        samples: &[SampleGenericSections<'_>],
        edges: &DepthBinEdges,
        cap: DepthCap,
    ) -> Vec<f64> {
        samples
            .iter()
            .map(|sections| {
                let positions = sections
                    .first()
                    .map_or(0, |(_, records)| records.depth().len());
                let mut total = 0.0_f64;
                let mut counted = 0_u64;
                for index in 0..positions {
                    // A sample's depth at a position is its read groups' depths added, so the
                    // mean has to be over positions and not over read-group entries.
                    let mut here = 0.0_f64;
                    let mut walked = false;
                    for (_, records) in sections.iter() {
                        if let DepthCode::Binned(bin) = records.depth().get(index) {
                            let range = cap.denominator_for(edges.depth_range(bin));
                            here += 0.5 * f64::from(*range.start() + *range.end());
                            walked = true;
                        }
                    }
                    if walked {
                        total += here;
                        counted += 1;
                    }
                }
                if counted == 0 {
                    1.0
                } else {
                    (total / counted as f64).max(1e-3)
                }
            })
            .collect()
    }

    /// A cursor over `first..end` only — **what lets one pass over the positions be split
    /// across cores.** Each sample's sparse list is binary-searched once to find where the
    /// chunk begins, and walked with a cursor from there, so a chunk costs one search per
    /// sample rather than one per position.
    fn over(
        samples: &'a [SampleGenericSections<'a>],
        edges: &'a DepthBinEdges,
        cap: DepthCap,
        coverage: &'a [f64],
        as_a_range: bool,
        first: usize,
        end: usize,
    ) -> Self {
        let at = samples
            .iter()
            .map(|sections| {
                sections
                    .iter()
                    .map(|(_, records)| {
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
            cap,
            coverage,
            as_a_range,
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
        for (s, sections) in self.samples.iter().enumerate() {
            let slot = &mut into.samples[s];
            *slot = SampleAtPosition::default();
            // The depths this sample's stored codes could stand for, added across its read
            // groups. Two groups' ranges are summed endpoint to endpoint, which is the range
            // the total depth lies in.
            let (mut shallowest, mut deepest) = (0_u32, 0_u32);
            for (g, (_, records)) in sections.iter().enumerate() {
                if let DepthCode::Binned(bin) = records.depth().get(index) {
                    // **The counts' own denominator, not the position's depth.** The reads
                    // subtracted below were thinned to the cap; subtracting them from an
                    // unthinned depth would charge the position reference reads it never
                    // had, and at a few hundred reads a position that is most of them.
                    let range = self.cap.denominator_for(self.edges.depth_range(bin));
                    shallowest += *range.start();
                    deepest += *range.end();
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
            // Reads on an insertion, a deletion or a spanning deletion are not a base and are
            // held out of the model entirely, so the depth they occupied goes with them.
            let held_out = slot.on[4];
            let disagreeing = slot.non_reference() - held_out;
            // How many reference reads each end of the range implies. A range that would
            // demand fewer than zero is cut at zero: the reads are the harder evidence, so the
            // depth gives way rather than charging a negative count of reference reads.
            let (shallowest, deepest) = if self.as_a_range {
                (f64::from(shallowest), f64::from(deepest))
            } else {
                // The point-read this replaced: the middle of the range, given to the
                // likelihood as though it were the depth.
                let middle = 0.5 * f64::from(shallowest + deepest);
                (middle, middle)
            };
            let fewest = (shallowest - held_out - disagreeing).max(0.0);
            let most = (deepest - held_out - disagreeing).max(0.0);
            slot.fewest_reference = fewest;
            slot.spread = (most - fewest) as usize + 1;
            assert!(
                slot.spread <= MAX_RECORDED_SPREAD,
                "a stored depth code stands for {} depths, more than the {MAX_RECORDED_SPREAD} \
                 this fit reserves room for",
                slot.spread
            );
            slot.depth = held_out + disagreeing + 0.5 * (fewest + most);
            let weights = &mut into.depth_weights[s * MAX_RECORDED_SPREAD..];
            fill_depth_weights(&mut weights[..slot.spread], self.coverage[s], fewest);
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
/// **The reference-read term and the non-reference total are handed in rather than taken
/// here**, because neither depends on which base the alternative is. [`ln_reference_reads`]
/// reads the sample, its depth weights and how many copies of the alternative it carries — and
/// it is the one logarithm in this block — so the caller takes it once for the three copy
/// counts and reuses it across every candidate allele and across the invariant branch beside
/// them. Taken inside this function it was ten calls a sample where three do.
fn ln_reads_given_genotype(
    sample: &SampleAtPosition,
    non_reference: f64,
    alt_copies: u8,
    alternative: usize,
    logs: &ReadLogs,
    reference_term: f64,
) -> f64 {
    let copies = alt_copies as usize;
    let candidate_reads = sample.on[alternative];
    // `Other` is neither the candidate nor the reference; it is held out of the model
    // entirely, so the depth it occupied is removed with it.
    let other_reads = non_reference - candidate_reads - sample.on[4];

    // **A count of zero needs no guard now.** `count_times_ln` branched on it only to keep
    // `0 · −∞` out of the sum; the logarithms below are clamped at `MIN_POSITIVE` when they are
    // built, so every one of them is finite and `0 · finite` is zero.
    candidate_reads * logs.ln_candidate[copies] + other_reads * logs.ln_neither + reference_term
}

/// `ln P(this sample's reads | every copy is the reference base)` — the term the invariant
/// branch needs, and it does not depend on which base would have been the alternative. Its
/// `reference_term` is the zero-copy entry of the same triple.
fn ln_reads_given_all_reference(
    sample: &SampleAtPosition,
    non_reference: f64,
    logs: &ReadLogs,
    reference_term: f64,
) -> f64 {
    (non_reference - sample.on[4]) * logs.ln_neither + reference_term
}

/// `ln P(the reference reads | this sample carries 0, 1 or 2 copies of the alternative)`.
///
/// **Candidate-invariant, which is why it is a function of its own.** Nothing in
/// [`ln_reference_reads`] reads *which* base the alternative is — only the sample's reference
/// count, the depths its stored code allows and the per-copy probability from its read group's
/// table — so one triple serves all three candidate alleles and the invariant branch.
fn reference_terms(sample: &SampleAtPosition, depth_weights: &[f64], logs: &ReadLogs) -> [f64; 3] {
    [
        ln_reference_reads(
            sample,
            depth_weights,
            logs.reference[0],
            logs.ln_reference[0],
        ),
        ln_reference_reads(
            sample,
            depth_weights,
            logs.reference[1],
            logs.ln_reference[1],
        ),
        ln_reference_reads(
            sample,
            depth_weights,
            logs.reference[2],
            logs.ln_reference[2],
        ),
    ]
}

/// Every logarithm the per-position kernel needs that depends on a read group's error rate
/// rather than on the position.
///
/// **Built once for the whole pass**, for exactly the reason [`BetaQuadrature`] beside it is: an
/// error rate is fixed for a pass, so `ln P(a read shows …)` is fixed with it. Taken per position
/// it was up to eighteen logarithms a sample — two noise classes × three candidate alleles ×
/// three genotypes — at every kept position of every pass.
struct ReadLogs {
    /// `ln P(a read shows the candidate allele)`, indexed by how many copies of it the sample
    /// carries.
    ln_candidate: [f64; 3],
    /// `P(a read shows the reference base)` by the same index — the power series over a depth
    /// range needs the probability itself and not only its logarithm.
    reference: [f64; 3],
    ln_reference: [f64; 3],
    /// `ln P(a read shows one of the two bases that are neither)`. No genotype moves it.
    ln_neither: f64,
}

impl ReadLogs {
    fn of(error_rate: f64, ploidy: Ploidy) -> Self {
        // The same clamp `count_times_ln` applied before taking a logarithm, so the table holds
        // the values the per-call arithmetic held.
        let ln = |p: f64| p.max(f64::MIN_POSITIVE).ln();
        let mut ln_candidate = [0.0; 3];
        let mut reference = [0.0; 3];
        let mut ln_reference = [0.0; 3];
        for copies in 0..3_usize {
            let carried = copies as f64 / f64::from(ploidy.get());
            let on_candidate = carried * (1.0 - error_rate) + (1.0 - carried) * error_rate / 3.0;
            let on_reference = reference_read_probability(copies as u8, ploidy, error_rate);
            ln_candidate[copies] = ln(on_candidate);
            reference[copies] = on_reference;
            ln_reference[copies] = ln(on_reference);
        }
        Self {
            ln_candidate,
            reference,
            ln_reference,
            ln_neither: ln(error_rate / 3.0),
        }
    }
}

/// How much weight each depth in a sample's range carries, relative to the shallowest.
///
/// # Why the depths in a range are not equally likely
///
/// Summing over the depths a stored code could stand for needs a statement about **which of
/// them the position more probably had**, and *all of them equally* is not it. Two things pull
/// the other way and they nearly cancel, which is why the answer is short: a deeper position
/// has more ways to have produced the reads that were seen — the multinomial's own coefficient
/// — and a deeper position is rarer, because a sample's read count at a position is a Poisson
/// draw around its own coverage. Writing both down, the factorials cancel and what is left is
/// that **the reference reads are themselves Poisson**, at this sample's coverage times the
/// chance a read shows the reference base, cut to the depths the code allows.
///
/// So the weight on the `i`-th depth above the shallowest the code allows is
/// `coverage^i / (fewest + i)!` divided through by the first, which needs no factorial: each
/// step is one multiplication. Dividing through is a constant this sample contributes to every
/// genotype alike, so it never has to be put back.
///
/// **Neither pull may be dropped.** With the coefficient alone the sum collapses onto the
/// deepest depth in the range, and with the Poisson alone onto the shallowest — either way the
/// fix becomes a second point-read at an edge of the bin, which is worse than the midpoint it
/// replaced.
fn fill_depth_weights(weights: &mut [f64], coverage: f64, fewest_reference: f64) {
    let mut running = 1.0;
    for (step, slot) in weights.iter_mut().enumerate() {
        if step > 0 {
            running *= coverage / (fewest_reference + step as f64);
        }
        *slot = running;
    }
}

/// The chance one read shows the reference base, from an individual carrying `alt_copies`
/// copies of the candidate allele.
fn reference_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    (1.0 - carried) * (1.0 - error_rate) + carried * error_rate / 3.0
}

/// **How many reference reads this sample is expected to have had**, given its genotype and
/// the range its stored code allows.
///
/// The counterpart of [`ln_reference_reads`], and it has to exist for the same reason: the
/// error rate is maximised over expected read counts, and an expectation taken at the middle
/// of the range while the likelihood sums over the range is not the same statement twice. Left
/// inconsistent, the two disagree in one direction — the likelihood prefers a deeper position
/// for a homozygous-reference sample that showed a couple of disagreeing reads, the middle of
/// the range books it fewer reference reads than that, and the error rate comes back **24%
/// above the truth on a drawn cohort at eight reads a position** where the consistent pair
/// returns it to within 6%.
fn expected_reference_reads(
    sample: &SampleAtPosition,
    depth_weights: &[f64],
    on_reference: f64,
) -> f64 {
    if sample.spread <= 1 {
        return sample.fewest_reference;
    }
    let probability = on_reference.max(f64::MIN_POSITIVE);
    let (mut total, mut weighted, mut power) = (0.0, 0.0, 1.0);
    for (step, weight) in depth_weights[..sample.spread].iter().enumerate() {
        let term = weight * power;
        total += term;
        weighted += term * (sample.fewest_reference + step as f64);
        power *= probability;
    }
    if total > 0.0 {
        weighted / total
    } else {
        sample.fewest_reference
    }
}

/// `ln P(the reference reads | a read shows the reference base with probability p)`, summed
/// over every depth the sample's stored code could have come from.
///
/// The shallowest depth's term is factored out, which is what keeps the sum in range: what is
/// left runs from one upwards and cannot underflow however small `p` is.
fn ln_reference_reads(
    sample: &SampleAtPosition,
    depth_weights: &[f64],
    on_reference: f64,
    ln_on_reference: f64,
) -> f64 {
    let shallowest = sample.fewest_reference * ln_on_reference;
    if sample.spread <= 1 {
        return shallowest;
    }
    let probability = on_reference.max(f64::MIN_POSITIVE);
    let mut total = 0.0;
    let mut power = 1.0;
    for weight in &depth_weights[..sample.spread] {
        total += weight * power;
        power *= probability;
    }
    shallowest + total.ln()
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
    /// `None` when the run does not fit the class.
    duplicated: Option<DuplicatedPositions>,
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
    /// Posterior mass on the duplicated class, and the same two sums for its carrier Beta.
    duplicated: f64,
    sum_ln_q: f64,
    sum_ln_one_minus_q: f64,
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
    /// Per sample, the posterior that it carries an extra copy, summed over positions.
    carrier: Vec<f64>,
    with_reads: Vec<u64>,
    with_two_reads: Vec<u64>,
    /// **One entry per position visited, and empty unless the pass was asked for it.** Every
    /// other field here is a sum over positions, which is what lets a chunk's counts be added
    /// to another chunk's; this one is a value *per* position, so a chunk's entries are
    /// appended rather than added and the pass that collects them keeps the chunks in
    /// position order. Only the last pass of a run is asked, so the cost is one 8 MB vector
    /// and not one per iteration.
    noisy_posterior: Vec<f32>,
    collect_noisy_posterior: bool,
    /// **Per position visited, then per sample, the posterior that the sample is heterozygous
    /// there, that both its copies are non-reference, and that it carries an extra copy** —
    /// the same position-order list as `noisy_posterior` with three values a sample instead of
    /// one. Empty unless the run asked for it.
    genotype_posterior: Vec<f32>,
    /// **Per position visited, the posterior that the position is drawn from the duplicated
    /// class**, in the same position order. Collected with `genotype_posterior`.
    duplicated_posterior: Vec<f32>,
    collect_genotype_posterior: bool,
}

/// What one pass of the alternation ended at. **Collected only when a run asks for it**, and
/// what it is for is telling a fit that has settled from one that ran out of passes still
/// moving — and, when it is still moving, in which direction.
#[derive(Clone, PartialEq, Debug)]
pub struct PassSummary {
    pub pass: u32,
    pub log_likelihood: f64,
    /// The largest relative move any parameter made in this pass's maximisation.
    pub largest_move: f64,
    /// Per sample, in the order the fit iterates samples.
    pub heterozygous: Vec<f64>,
    pub homozygous_alt: Vec<f64>,
    pub hom_excess: Vec<f64>,
    pub expected_heterozygosity: f64,
    pub noisy_share: f64,
    pub density_a: f64,
    pub density_b: f64,
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
            duplicated: 0.0,
            sum_ln_q: 0.0,
            sum_ln_one_minus_q: 0.0,
            sum_ln_f: 0.0,
            sum_ln_one_minus_f: 0.0,
            reads: vec![[[ReadTally::default(); 3]; 2]; groups],
            genotypes: vec![vec![[0.0; 3]; nodes]; samples],
            heterozygous: vec![0.0; samples],
            homozygous_alt: vec![0.0; samples],
            carrier: vec![0.0; samples],
            with_reads: vec![0; samples],
            with_two_reads: vec![0; samples],
            noisy_posterior: Vec::new(),
            collect_noisy_posterior: false,
            genotype_posterior: Vec::new(),
            duplicated_posterior: Vec::new(),
            collect_genotype_posterior: false,
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
        self.duplicated += other.duplicated;
        self.sum_ln_q += other.sum_ln_q;
        self.sum_ln_one_minus_q += other.sum_ln_one_minus_q;
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
        for (into, from) in self.carrier.iter_mut().zip(other.carrier.iter()) {
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
        self.noisy_posterior
            .extend_from_slice(&other.noisy_posterior);
        self.genotype_posterior
            .extend_from_slice(&other.genotype_posterior);
        self.duplicated_posterior
            .extend_from_slice(&other.duplicated_posterior);
    }
}

/// The read counts one sample contributes to a read group's tally, under one candidate allele
/// and one genotype.
///
/// The counts on the candidate allele and on the two bases that are neither are exact; the
/// reference count is the one that has to be taken in expectation, because the depth it is
/// derived from is a range.
fn tally_of(
    sample: &SampleAtPosition,
    depth_weights: &[f64],
    alternative: usize,
    on_reference: f64,
) -> ReadTally {
    let candidate = sample.on[alternative];
    ReadTally {
        candidate,
        neither: sample.non_reference() - candidate - sample.on[4],
        reference: expected_reference_reads(sample, depth_weights, on_reference),
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
/// [`JointFitError::NoSamples`] on an empty cohort and [`JointFitError::NotDiploid`] on a
/// ploidy this estimator does not model. **The refusal for samples that did not record the same
/// thing has already happened**: building a [`CohortCensusEvidence`] is what makes it, before a
/// section is read.
pub fn fit_jointly(
    cohort: &mut CohortCensusEvidence,
    config: &JointFitConfig,
) -> Result<JointFit, JointFitError> {
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();
    let first = names.first().ok_or(JointFitError::NoSamples)?.clone();
    if config.ploidy.get() != 2 {
        return Err(JointFitError::NotDiploid {
            sample: first,
            ploidy: config.ploidy.get(),
        });
    }
    // **The cap every sample recorded under.** One number, because the cohort refused any panel
    // whose samples disagree on it before a section was read.
    let depth_cap = cohort
        .terms()
        .map_or(DepthCap::MAX, |terms| terms.depth_cap);
    let groups: Vec<ReadGroupId> = cohort.read_groups().to_vec();

    // **Every section the generic half needs, lent for the length of one call.** The estimator
    // reads a position from every sample at once, so what it borrows is one row a sample — and
    // when this call returns, a file-backed census has nothing left decoded.
    let (score, parameters, statistics, passes, converged, trace, contamination) = cohort
        .with_generic(&groups, |lent| {
            // Which read group each sample's own groups are, in the order the cursor visits them.
            let group_index: Vec<Vec<usize>> = lent
                .iter()
                .map(|sections| {
                    sections
                        .iter()
                        .map(|(id, _)| {
                            groups
                                .iter()
                                .position(|g| g == id)
                                .expect("every group came from this list")
                        })
                        .collect()
                })
                .collect();

            // Each sample's own mean depth, read once from the codes: the centre of the prior a
            // stored code's range is weighted with when the likelihood sums over it.
            let coverage = EvidenceCursor::mean_depth(lent, &config.edges, depth_cap);

            let mut best: Option<(f64, Parameters, Statistics, u32, bool, Vec<PassSummary>)> = None;
            for start in &config.starting_points {
                let (parameters, statistics, passes, converged, trace) = maximise(
                    lent,
                    depth_cap,
                    config,
                    &groups,
                    &group_index,
                    &coverage,
                    start,
                );
                let score = statistics.log_likelihood;
                if best.as_ref().is_none_or(|(current, ..)| score > *current) {
                    best = Some((score, parameters, statistics, passes, converged, trace));
                }
            }
            let (score, parameters, statistics, passes, converged, trace) =
                best.expect("a run always has at least one starting point");

            // **Contamination is fitted after the alternation, not inside it** (spec §3.4). It reads
            // the converged error rates and the converged homozygote excess, and nothing it produces
            // feeds back into them — a sample's stray reads are a property of the tube it was in
            // rather than of the population, so the density has no business being told about them.
            // It runs here, inside the same call, so that the sections are lent once rather than
            // twice.
            let per_sample_error: Vec<f64> = (0..lent.len())
                .map(|s| parameters.clean[group_index[s][0]])
                .collect();
            // **The mismapped positions are kept out of it.** A position where two stretches of
            // genome pile up on one place puts a small share of unexpected reads into *every*
            // sample, which is the contamination signature exactly; measured over 63 tomato
            // accessions with those positions left in, the median accession came back 6.5%
            // contaminated.
            let contamination = fit_contamination_over(
                lent,
                depth_cap,
                &config.edges,
                &per_sample_error,
                &parameters.hom_excess,
                &statistics.noisy_posterior,
                &config.contamination,
            );
            (
                score,
                parameters,
                statistics,
                passes,
                converged,
                trace,
                contamination,
            )
        })?;
    let mut statistics = statistics;
    let genotype_posterior = std::mem::take(&mut statistics.genotype_posterior);
    let duplicated_posterior = std::mem::take(&mut statistics.duplicated_posterior);

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
    let hom_excess = names
        .iter()
        .enumerate()
        .map(|(s, name)| {
            (
                name.clone(),
                Estimate {
                    value: HomozygoteExcess::try_new(parameters.hom_excess[s])
                        .expect("the maximisation is confined to [0, 1]"),
                    provenance: if names.len() >= 2 {
                        Provenance::FittedHere
                    } else {
                        Provenance::Defaulted
                    },
                    observations: statistics.with_reads[s],
                },
            )
        })
        .collect();
    let rates = names
        .iter()
        .enumerate()
        .map(|(s, name)| {
            (
                name.clone(),
                Estimate {
                    value: SampleGenotypeRates {
                        heterozygous: statistics.heterozygous[s] / statistics.positions,
                        homozygous_alt: statistics.homozygous_alt[s] / statistics.positions,
                        duplicated_carrier: statistics.carrier[s] / statistics.positions,
                        positions_with_reads: statistics.with_reads[s],
                        positions_with_two_reads: statistics.with_two_reads[s],
                    },
                    provenance: Provenance::FittedHere,
                    observations: statistics.with_reads[s],
                },
            )
        })
        .collect();
    let contamination = names.iter().cloned().zip(contamination).collect();

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
        duplicated: parameters.duplicated.map(|value| Estimate {
            value,
            provenance: Provenance::FittedHere,
            observations,
        }),
        noisy_posterior: statistics.noisy_posterior,
        genotype_posterior,
        duplicated_posterior,
        trace,
        passes,
        converged,
        log_likelihood: score,
    })
}

/// One run of the alternation, from one starting point.
fn maximise(
    samples: &[SampleGenericSections<'_>],
    depth_cap: DepthCap,
    config: &JointFitConfig,
    groups: &[ReadGroupId],
    group_index: &[Vec<usize>],
    coverage: &[f64],
    start: &StartingPoint,
) -> (Parameters, Statistics, u32, bool, Vec<PassSummary>) {
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
        duplicated: config.duplicated_positions.then_some(DuplicatedPositions {
            share: start.duplicated_share,
            carrier_a: start.carrier_a,
            carrier_b: start.carrier_b,
        }),
    };
    let mut statistics;
    let mut converged = false;
    let mut passes = 0;
    let mut previous = f64::NEG_INFINITY;
    let mut trace = Vec::new();
    for pass in 1..=config.max_passes {
        passes = pass;
        statistics = expectation(
            samples,
            depth_cap,
            config,
            group_index,
            coverage,
            &parameters,
        );
        // The parameters as they stood when this pass read the data — the ones its rates
        // belong to. The maximisation below moves them, and the next entry carries where to.
        let entering = config.pass_trace.then(|| {
            (
                parameters.hom_excess.clone(),
                parameters.density,
                parameters.noisy_share,
            )
        });
        let moved = maximisation(&mut parameters, &statistics, config);
        if let Some((hom_excess, density, noisy_share)) = entering {
            let positions = statistics.positions.max(1.0);
            trace.push(PassSummary {
                pass,
                log_likelihood: statistics.log_likelihood,
                largest_move: moved,
                heterozygous: statistics
                    .heterozygous
                    .iter()
                    .map(|h| h / positions)
                    .collect(),
                homozygous_alt: statistics
                    .homozygous_alt
                    .iter()
                    .map(|h| h / positions)
                    .collect(),
                hom_excess,
                expected_heterozygosity: density.expected_heterozygosity(),
                noisy_share,
                density_a: density.a,
                density_b: density.b,
            });
        }
        let gain = statistics.log_likelihood - previous;
        previous = statistics.log_likelihood;
        if moved < config.stillness && gain.abs() < config.stillness * statistics.positions.max(1.0)
        {
            converged = true;
            break;
        }
    }
    // The reported statistics must be the ones the reported parameters produce, so the last
    // maximisation is followed by one more pass rather than by the pass that preceded it —
    // and it is this pass, at the parameters that will be reported, that keeps each
    // position's probability of being mismapped.
    statistics = expectation_pass(
        samples,
        depth_cap,
        config,
        group_index,
        coverage,
        &parameters,
        true,
    );
    (parameters, statistics, passes, converged, trace)
}

/// One pass over every position: the posteriors, and every count the maximisations need.
///
/// **Split across cores by position.** Positions are independent given the parameters, and a
/// chunk's counts add to another chunk's, so the pass is a map and a sum.
fn expectation(
    samples: &[SampleGenericSections<'_>],
    depth_cap: DepthCap,
    config: &JointFitConfig,
    group_index: &[Vec<usize>],
    coverage: &[f64],
    parameters: &Parameters,
) -> Statistics {
    expectation_pass(
        samples,
        depth_cap,
        config,
        group_index,
        coverage,
        parameters,
        false,
    )
}

/// The same pass, told whether to keep each position's probability of being mismapped.
///
/// **Kept only on the last pass of a run.** Keeping it costs one four-byte value a position,
/// and the chunks have to be held until they can be joined in position order rather than
/// summed as they arrive — so the iterating passes use the streaming sum and pay neither.
fn expectation_pass(
    samples: &[SampleGenericSections<'_>],
    depth_cap: DepthCap,
    config: &JointFitConfig,
    group_index: &[Vec<usize>],
    coverage: &[f64],
    parameters: &Parameters,
    collect_noisy_posterior: bool,
) -> Statistics {
    let quadrature = BetaQuadrature::with_genotype_priors(
        parameters.density.a,
        parameters.density.b,
        config.quadrature_nodes,
        &parameters.hom_excess,
    );
    let carrier = parameters.duplicated.map(|duplicated| {
        BetaQuadrature::new(
            duplicated.carrier_a,
            duplicated.carrier_b,
            config.quadrature_nodes,
        )
    });
    // One table a noise class a read group, beside the quadratures and for the same reason.
    let read_logs: Vec<Vec<ReadLogs>> = (0..2)
        .map(|class| {
            (0..parameters.clean.len())
                .map(|group| ReadLogs::of(class_rate(parameters, class, group), config.ploidy))
                .collect()
        })
        .collect();
    let positions = EvidenceCursor::position_count(samples);
    let chunk = POSITIONS_PER_CHUNK.min(positions.div_ceil(rayon::current_num_threads()).max(1));
    let bounds: Vec<(usize, usize)> = (0..positions)
        .step_by(chunk)
        .map(|first| (first, (first + chunk).min(positions)))
        .collect();

    // The per-position lists are kept only on the pass whose parameters will be reported.
    let collect_genotype_posterior = collect_noisy_posterior && config.genotype_posteriors;
    let empty = || {
        let mut statistics = Statistics::new(
            parameters.clean.len(),
            samples.len(),
            quadrature.nodes.len(),
        );
        statistics.collect_noisy_posterior = collect_noisy_posterior;
        statistics.collect_genotype_posterior = collect_genotype_posterior;
        statistics
    };
    let one_chunk = |(first, end): (usize, usize)| {
        let mut statistics = empty();
        if collect_noisy_posterior {
            statistics.noisy_posterior.reserve(end - first);
        }
        if collect_genotype_posterior {
            statistics
                .genotype_posterior
                .reserve((end - first) * samples.len() * 3);
            statistics.duplicated_posterior.reserve(end - first);
        }
        let mut cursor = EvidenceCursor::over(
            samples,
            &config.edges,
            depth_cap,
            coverage,
            config.depth_as_a_range,
            first,
            end,
        );
        let mut scratch = Scratch::new(samples.len(), quadrature.nodes.len());
        let mut odds = vec![
            1.0_f64;
            if config.coverage_odds.is_empty() {
                0
            } else {
                samples.len()
            }
        ];
        let mut index = first;
        while cursor.next_position(&mut scratch.evidence) {
            for (slot, sample) in odds.iter_mut().enumerate() {
                *sample = f64::from(
                    config.coverage_odds[slot]
                        .get(index)
                        .copied()
                        .unwrap_or(0.0),
                )
                .exp();
            }
            one_position(
                &mut scratch,
                group_index,
                config.ploidy,
                &read_logs,
                &quadrature,
                carrier.as_ref(),
                &odds,
                parameters,
                &mut statistics,
            );
            index += 1;
        }
        statistics
    };

    if collect_noisy_posterior {
        // `reduce` may join chunks in any order it likes, which is fine for a sum and wrong
        // for a list in position order, so the chunks are collected first.
        let chunks: Vec<Statistics> = bounds.into_par_iter().map(one_chunk).collect();
        let mut into = empty();
        into.noisy_posterior.reserve(positions);
        if collect_genotype_posterior {
            into.genotype_posterior
                .reserve(positions * samples.len() * 3);
            into.duplicated_posterior.reserve(positions);
        }
        for chunk in &chunks {
            into.absorb(chunk);
        }
        into
    } else {
        bounds
            .into_par_iter()
            .map(one_chunk)
            .reduce(empty, |mut into, from| {
                into.absorb(&from);
                into
            })
    }
}

/// How many positions one core takes at a time. Large enough that the per-chunk binary search
/// into each sample's sparse list disappears against the work, small enough that a cohort of
/// fifty still fills every core.
const POSITIONS_PER_CHUNK: usize = 16_384;

/// At most three non-reference bases can be the segregating one, so the candidate list never
/// exceeds three however many the samples showed reads on.
const MAX_CANDIDATES: usize = CANDIDATE_ALTERNATIVES;

/// What a position can be, within a noise class: the population carries one allele, it is fixed
/// on a non-reference one, it segregates, or it is a stretch some samples carry twice.
const BRANCHES: usize = 4;
const DUPLICATED: usize = 3;

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
    /// The same, over the carrier frequency the duplicated class integrates.
    carrier_node_ln: Vec<f64>,
    fixed_ln: Vec<f64>,
    /// `[class]` — the invariant branch, and the four branches combined.
    invariant_ln: Vec<f64>,
    branch_ln: Vec<f64>,
    class_ln: Vec<f64>,
    /// The candidate alleles this position sums over, and how many alleles each stands for.
    candidates: Vec<usize>,
    multiplicity: Vec<f64>,
    /// `[sample][genotype]` — the segregating branch's genotype weight, collapsed over the
    /// nodes and candidates it was spread across.
    genotype_weight: Vec<f64>,
    /// `[sample]` — the same for the duplicated branch, where a sample has two states rather
    /// than three: it carries the extra copy or it does not.
    carrier_weight: Vec<f64>,
    /// `[sample][heterozygous, homozygous non-reference, carries an extra copy]` — this one
    /// position's posteriors, summed over both classes and every branch, so that a run asking
    /// to keep them has a value to keep. Zeroed per position.
    position_genotype: Vec<f64>,
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
                depth_weights: vec![0.0; samples * MAX_RECORDED_SPREAD],
                observed_alternatives: Vec::with_capacity(MAX_CANDIDATES),
            },
            samples,
            nodes,
            ell: vec![0.0; 2 * MAX_CANDIDATES * samples * 3],
            lik: vec![0.0; 2 * MAX_CANDIDATES * samples * 3],
            lik_max: vec![0.0; 2 * MAX_CANDIDATES * samples],
            node_ln: vec![f64::NEG_INFINITY; 2 * MAX_CANDIDATES * nodes],
            carrier_node_ln: vec![f64::NEG_INFINITY; 2 * MAX_CANDIDATES * nodes],
            fixed_ln: vec![f64::NEG_INFINITY; 2 * MAX_CANDIDATES],
            invariant_ln: vec![0.0; 2],
            branch_ln: vec![0.0; BRANCHES * 2],
            class_ln: vec![0.0; 2],
            candidates: Vec::with_capacity(MAX_CANDIDATES),
            multiplicity: Vec::with_capacity(MAX_CANDIDATES),
            genotype_weight: vec![0.0; samples * 3],
            carrier_weight: vec![0.0; samples],
            position_genotype: vec![0.0; samples * 3],
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
#[allow(
    clippy::too_many_arguments,
    reason = "one position's whole problem: its evidence, the two quadratures, the coverage \
              readings, the parameters and every count they feed"
)]
fn one_position(
    scratch: &mut Scratch,
    group_index: &[Vec<usize>],
    ploidy: Ploidy,
    read_logs: &[Vec<ReadLogs>],
    quadrature: &BetaQuadrature,
    carrier: Option<&BetaQuadrature>,
    coverage_odds: &[f64],
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
    //
    // **The sample loop is on the outside so that the two candidate-invariant quantities are
    // taken once.** [`reference_terms`] holds this block's only logarithms and
    // `non_reference()` is a sum over four counts; neither reads which base the alternative
    // is, so with the candidate loop outermost both were recomputed for every candidate — ten
    // reference terms a sample where three do the same job. The sums are unchanged: `invariant`
    // and each candidate's `fixed` still accumulate over the samples in the same order.
    for class in 0..2 {
        let mut invariant = 0.0;
        let mut fixed = [0.0_f64; MAX_CANDIDATES];
        for s in 0..samples {
            let logs = &read_logs[class][group_index[s][0]];
            let sample = &scratch.evidence.samples[s];
            let weights = scratch.evidence.weights_of(s);
            let non_reference = sample.non_reference();
            let reference_term = reference_terms(sample, weights, logs);
            invariant +=
                ln_reads_given_all_reference(sample, non_reference, logs, reference_term[0]);
            for candidate in 0..candidates {
                let allele = scratch.candidates[candidate];
                let base = scratch.ell_at(class, candidate, s);
                let mut largest = f64::NEG_INFINITY;
                for j in 0..3 {
                    let value = ln_reads_given_genotype(
                        sample,
                        non_reference,
                        j as u8,
                        allele,
                        logs,
                        reference_term[j],
                    );
                    scratch.ell[base + j] = value;
                    largest = largest.max(value);
                }
                fixed[candidate] += scratch.ell[base + 2];
                let slot = scratch.max_at(class, candidate, s);
                scratch.lik_max[slot] = largest;
                for j in 0..3 {
                    scratch.lik[base + j] = (scratch.ell[base + j] - largest).exp();
                }
            }
        }
        scratch.invariant_ln[class] = invariant;
        for candidate in 0..candidates {
            scratch.fixed_ln[class * MAX_CANDIDATES + candidate] = fixed[candidate];
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

    // ---- the integral over how much of the panel carries the extra copy ---------------------
    //
    // The same shape as the branch above with one difference: a sample here has two states
    // rather than three. It carries the duplication, in which case about half its reads
    // disagree — the identical read likelihood a heterozygote has, already computed — or it
    // does not, and is homozygous reference. **There is no third state, and that absence is
    // what tells the class from a real variant across a cohort.**
    if let Some(carrier) = carrier {
        for class in 0..2 {
            for candidate in 0..candidates {
                let mut offset = 0.0;
                for s in 0..samples {
                    offset += scratch.lik_max[scratch.max_at(class, candidate, s)];
                }
                for node in 0..carrier.nodes.len() {
                    let q = carrier.nodes[node];
                    let mut product = 1.0_f64;
                    let mut scale = 0.0_f64;
                    for s in 0..samples {
                        let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                        // Where the run has a coverage summary, this is the only place it
                        // enters: how much more likely this sample's read depth around here is
                        // if it has two copies rather than one. Every other branch has every
                        // sample single-copy, so the one-copy factor is common to all of them
                        // and cancels out of the position entirely.
                        let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                        let term = (1.0 - q) * lik[0] + q * odds * lik[1];
                        product *= term;
                        if product < 1.0 / RESCALE {
                            product *= RESCALE;
                            scale -= LN_RESCALE;
                        }
                    }
                    let slot = scratch.node_at(class, candidate, node);
                    scratch.carrier_node_ln[slot] = if product <= 0.0 {
                        f64::NEG_INFINITY
                    } else {
                        product.ln() + scale + offset + carrier.ln_weights[node]
                    };
                }
            }
        }
    }

    // ---- the four branches, and the two classes ---------------------------------------------
    let ln_three = (CANDIDATE_ALTERNATIVES as f64).ln();
    let duplicated_share = parameters.duplicated.map_or(0.0, |d| d.share);
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
        // The three ordinary branches share what is left when the duplicated class has taken
        // its share, so the four still sum to one.
        let ordinary = (1.0 - duplicated_share).max(f64::MIN_POSITIVE).ln();
        scratch.branch_ln[class * BRANCHES] = ordinary
            + density.p_invariant.max(f64::MIN_POSITIVE).ln()
            + scratch.invariant_ln[class];
        scratch.branch_ln[class * BRANCHES + 1] =
            ordinary + density.p_fixed_alt.max(f64::MIN_POSITIVE).ln() + fixed_alt;
        scratch.branch_ln[class * BRANCHES + 2] =
            ordinary + density.p_segregating().max(f64::MIN_POSITIVE).ln() + segregating;
        scratch.branch_ln[class * BRANCHES + DUPLICATED] = if let Some(carrier) = carrier {
            scratch.shares.clear();
            for candidate in 0..candidates {
                let multiplicity = scratch.multiplicity[candidate].ln();
                for node in 0..carrier.nodes.len() {
                    scratch.shares.push(
                        scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                            + multiplicity,
                    );
                }
            }
            duplicated_share.max(f64::MIN_POSITIVE).ln() + ln_sum_exp(&scratch.shares) - ln_three
        } else {
            f64::NEG_INFINITY
        };
        let share = if class == 0 {
            1.0 - parameters.noisy_share
        } else {
            parameters.noisy_share
        };
        scratch.class_ln[class] = share.max(f64::MIN_POSITIVE).ln()
            + ln_sum_exp(&scratch.branch_ln[class * BRANCHES..][..BRANCHES]);
    }
    let position_ln = ln_sum_exp(&scratch.class_ln);
    if !position_ln.is_finite() {
        // A position whose likelihood underflowed contributes nothing to any parameter, and it
        // must still contribute an entry, or every position after it in the chunk would be
        // attributed to its neighbour.
        if statistics.collect_noisy_posterior {
            statistics.noisy_posterior.push(0.0);
        }
        if statistics.collect_genotype_posterior {
            statistics
                .genotype_posterior
                .extend(std::iter::repeat_n(0.0_f32, samples * 3));
            statistics.duplicated_posterior.push(0.0);
        }
        return;
    }
    scratch.position_genotype.fill(0.0);
    let mut duplicated_here = 0.0_f64;
    statistics.log_likelihood += position_ln;
    if statistics.collect_noisy_posterior {
        statistics
            .noisy_posterior
            .push((scratch.class_ln[1] - position_ln).exp() as f32);
    }

    // ---- attribute it -------------------------------------------------------------------------
    for class in 0..2 {
        let class_posterior = (scratch.class_ln[class] - position_ln).exp();
        if class_posterior <= 1e-12 {
            continue;
        }
        if class == 1 {
            statistics.noisy += class_posterior;
        }
        let branches = &scratch.branch_ln[class * BRANCHES..][..BRANCHES];
        let within = ln_sum_exp(branches);
        let branch = [
            class_posterior * (branches[0] - within).exp(),
            class_posterior * (branches[1] - within).exp(),
            class_posterior * (branches[2] - within).exp(),
            class_posterior * (branches[DUPLICATED] - within).exp(),
        ];
        statistics.invariant += branch[0];
        statistics.fixed_alt += branch[1];
        statistics.segregating += branch[2];
        statistics.duplicated += branch[DUPLICATED];
        duplicated_here += branch[DUPLICATED];

        // The invariant branch: every sample is homozygous reference, and every read that is
        // not on the reference base is an error, so none of them is "the allele".
        if branch[0] > 1e-12 {
            for s in 0..samples {
                let sample = &scratch.evidence.samples[s];
                let rate = class_rate(parameters, class, group_index[s][0]);
                let neither = sample.non_reference() - sample.on[4];
                let reference =
                    expected_reference_reads(sample, scratch.evidence.weights_of(s), 1.0 - rate);
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
                    let rate = class_rate(parameters, class, group_index[s][0]);
                    let counts = tally_of(
                        &scratch.evidence.samples[s],
                        scratch.evidence.weights_of(s),
                        allele,
                        reference_read_probability(2, ploidy, rate),
                    );
                    for &g in &group_index[s] {
                        let tally = &mut statistics.reads[g][class][2];
                        tally.candidate += share * counts.candidate;
                        tally.neither += share * counts.neither;
                        tally.reference += share * counts.reference;
                    }
                    statistics.homozygous_alt[s] += share;
                    scratch.position_genotype[s * 3 + 1] += share;
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
                    let rate = class_rate(parameters, class, group_index[s][0]);
                    let base = (candidate * samples + s) * 3;
                    for j in 0..3 {
                        let weight = scratch.per_candidate_weight[base + j];
                        if weight <= 0.0 {
                            continue;
                        }
                        let counts = tally_of(
                            &scratch.evidence.samples[s],
                            scratch.evidence.weights_of(s),
                            allele,
                            reference_read_probability(j as u8, ploidy, rate),
                        );
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
                scratch.position_genotype[s * 3] += scratch.genotype_weight[s * 3 + 1];
                scratch.position_genotype[s * 3 + 1] += scratch.genotype_weight[s * 3 + 2];
            }
        }

        // The duplicated branch. A carrier's reads are scored exactly as a heterozygote's, so
        // they go into the same tally the error rate is maximised over — what differs is where
        // the position's own weight is booked, and **a carrier is not counted heterozygous**.
        if let Some(carrier) = carrier.filter(|_| branch[DUPLICATED] > 1e-12) {
            {
                let carrier_nodes = carrier.nodes.len();
                scratch.shares.clear();
                for candidate in 0..candidates {
                    let multiplicity = scratch.multiplicity[candidate].ln();
                    for node in 0..carrier_nodes {
                        scratch.shares.push(
                            scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                                + multiplicity,
                        );
                    }
                }
                let total = ln_sum_exp(&scratch.shares);
                scratch.per_candidate_weight[..candidates * samples * 3].fill(0.0);
                scratch.carrier_weight.fill(0.0);
                for candidate in 0..candidates {
                    for node in 0..carrier_nodes {
                        let share = branch[DUPLICATED]
                            * (scratch.shares[candidate * carrier_nodes + node] - total).exp();
                        if share <= 1e-12 {
                            continue;
                        }
                        let q = carrier.nodes[node];
                        statistics.sum_ln_q += share * carrier.ln_nodes[node];
                        statistics.sum_ln_one_minus_q += share * carrier.ln_one_minus_nodes[node];
                        for s in 0..samples {
                            let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                            let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                            let joint = [(1.0 - q) * lik[0], q * odds * lik[1]];
                            let total = joint[0] + joint[1];
                            if total <= 0.0 {
                                continue;
                            }
                            let base = (candidate * samples + s) * 3;
                            scratch.per_candidate_weight[base] += share * joint[0] / total;
                            scratch.per_candidate_weight[base + 1] += share * joint[1] / total;
                            scratch.carrier_weight[s] += share * joint[1] / total;
                        }
                    }
                }
                for candidate in 0..candidates {
                    let allele = scratch.candidates[candidate];
                    for s in 0..samples {
                        let rate = class_rate(parameters, class, group_index[s][0]);
                        let base = (candidate * samples + s) * 3;
                        for j in 0..2 {
                            let weight = scratch.per_candidate_weight[base + j];
                            if weight <= 0.0 {
                                continue;
                            }
                            let counts = tally_of(
                                &scratch.evidence.samples[s],
                                scratch.evidence.weights_of(s),
                                allele,
                                reference_read_probability(j as u8, ploidy, rate),
                            );
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
                    statistics.carrier[s] += scratch.carrier_weight[s];
                    scratch.position_genotype[s * 3 + 2] += scratch.carrier_weight[s];
                }
            }
        }
    }

    if statistics.collect_genotype_posterior {
        statistics
            .genotype_posterior
            .extend(scratch.position_genotype.iter().map(|v| *v as f32));
        statistics.duplicated_posterior.push(duplicated_here as f32);
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
    let moved = std::cell::Cell::new(0.0_f64);
    let note = |before: f64, after: f64| {
        let scale = before.abs().max(1e-6);
        moved.set(moved.get().max((after - before).abs() / scale));
    };
    // **A share is judged against a floor rather than against itself.** A class the data does
    // not want shrinks geometrically — halving on every pass — and halving is a relative move
    // of one half however small the number has become, so a run on a cohort with no such
    // positions would never report convergence. One position in ten thousand is the point below
    // which a share stops being a quantity anyone reads.
    let note_share = |before: f64, after: f64| {
        let scale = before.abs().max(1e-4);
        moved.set(moved.get().max((after - before).abs() / scale));
    };

    // The share of positions that are mismapped, and the density's two masses: closed form.
    let positions = statistics.positions.max(1.0);
    let before = parameters.noisy_share;
    parameters.noisy_share = (statistics.noisy / positions).clamp(1e-6, 0.5);
    note_share(before, parameters.noisy_share);

    let branch_total = (statistics.invariant + statistics.fixed_alt + statistics.segregating)
        .max(f64::MIN_POSITIVE);
    let before = parameters.density.p_invariant;
    parameters.density.p_invariant = (statistics.invariant / branch_total).clamp(1e-9, 1.0 - 1e-9);
    note(before, parameters.density.p_invariant);
    let before = parameters.density.p_fixed_alt;
    parameters.density.p_fixed_alt = (statistics.fixed_alt / branch_total).clamp(1e-12, 0.5);
    note(before, parameters.density.p_fixed_alt);

    // The duplicated class's share of positions, and the Beta its carrier frequency is drawn
    // from. **Bounded well below a half**: a class that grows without bound starts explaining
    // real variants, which is the failure mode measured when duplications are mostly private.
    if let Some(duplicated) = parameters.duplicated.as_mut() {
        let before = duplicated.share;
        duplicated.share = (statistics.duplicated / positions).clamp(1e-9, 0.05);
        note_share(before, duplicated.share);
        if statistics.duplicated > 0.0 {
            let (a, b) = fit_beta_shapes(
                statistics.sum_ln_q / statistics.duplicated,
                statistics.sum_ln_one_minus_q / statistics.duplicated,
                duplicated.carrier_a,
                duplicated.carrier_b,
            );
            note(duplicated.carrier_a, a);
            note(duplicated.carrier_b, b);
            duplicated.carrier_a = a;
            duplicated.carrier_b = b;
        }
    }

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
    moved.get()
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

    // The whole fit, against a cohort whose truth is known, is in `whole_fit_tests` below —
    // the drawn-cohort generator it needs sits between the two modules, because the benchmark
    // reaches it too.
}

/// A cohort drawn at known parameters, for anything that has to fit evidence whose answer is
/// already known.
///
/// **Compiled under `cfg(test)` and under the `bench-fixtures` feature, and nowhere else.** The
/// two callers are this module's own oracle — a cohort drawn at a chosen error rate,
/// heterozygosity and inbreeding must come back at them — and the benchmark, which needs a
/// resident census so that [`fit_jointly`] can be timed with no CRAM and no reference genome
/// between the clock and the estimator. They draw from one generator because a benchmark drawn
/// differently from the oracle would be timing a workload no test has ever checked.
///
/// Nothing here is production code: a release build without the feature compiles none of it.
#[cfg(any(test, feature = "bench-fixtures"))]
pub mod bench_fixtures {
    use super::*;

    use crate::ng::parameter_estimation::joint::census::{
        AlleleObservation, DepthCap, DepthLadderDigest, GenericEvidence, ObservedAllele,
        PackedDepthCodes, ReadCap, RecordingTerms, SampleCensusEvidence, Section, SectionKey,
        SelectionTermsDigest,
    };
    use crate::ng::parameter_estimation::joint::loci::{
        CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
    };
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::tandem_repeat::ScanParams;

    fn selection_terms() -> SelectionTerms {
        SelectionTerms {
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

    /// The drawn samples as the cohort the fit takes. **Cloned**, because a drawn cohort is
    /// refitted several times in one test at different sample counts.
    ///
    /// A benchmark builds the cohort **once**, outside the timed region, and re-fits the same
    /// one: a resident census lends its sections and takes them back, so nothing is consumed
    /// and the clone here is setup a bench must not pay per iteration.
    pub fn as_cohort(samples: &[SampleCensusEvidence]) -> CohortCensusEvidence {
        CohortCensusEvidence::new(samples.to_vec()).expect("a drawn cohort records one way")
    }

    /// A drawn cohort's records, beside the parameters they were drawn at.
    pub struct DrawnCohort {
        pub samples: Vec<SampleCensusEvidence>,
        pub clean: f64,
        pub noisy: f64,
        pub noisy_share: f64,
        pub density: FrequencyDensity,
        pub hom_excess: Vec<f64>,
        pub heterozygous: Vec<f64>,
    }

    /// Draw a cohort at known parameters and write it into records the fit will read.
    pub fn draw_cohort(
        samples: usize,
        positions: usize,
        mean_depth: f64,
        truth: (f64, f64, f64),
        density: FrequencyDensity,
        hom_excess: f64,
        seed: u64,
    ) -> DrawnCohort {
        draw_cohort_with_duplications(
            samples, positions, mean_depth, truth, density, hom_excess, 0.0, seed,
        )
    }

    /// The same, with a share of positions at which some samples carry an extra copy of the
    /// stretch: **twice the reads and about half of them disagreeing, and never a sample
    /// homozygous for the non-reference allele.** Who carries one is drawn from the carrier
    /// frequency the eight-sample tomato counts were fitted to.
    #[allow(
        clippy::too_many_arguments,
        reason = "the drawn cohort's own parameters"
    )]
    pub fn draw_cohort_with_duplications(
        samples: usize,
        positions: usize,
        mean_depth: f64,
        truth: (f64, f64, f64),
        density: FrequencyDensity,
        hom_excess: f64,
        duplicated_share: f64,
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
        let edges = DepthBinEdges::for_census();

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
            let (duplicated, carrier_frequency) = if duplicated_share > 0.0 {
                let duplicated = draw.uniform() < duplicated_share;
                (duplicated, draw.beta(1.19, 9.55).max(0.2))
            } else {
                (false, 0.0)
            };
            for sample in 0..samples {
                let copies = if duplicated && draw.uniform() < carrier_frequency {
                    2.0
                } else {
                    1.0
                };
                let _ = carrier_frequency;
                let reads = draw.poisson(mean_depth * copies);
                let genotype = if duplicated {
                    usize::from(copies > 1.0)
                } else {
                    match branch {
                        0 => 0_usize,
                        1 => 2,
                        _ => draw.pick(&genotype_frequencies(f, excess[sample])),
                    }
                };
                if genotype == 1 && !duplicated {
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
                        reads: u8::try_from(*count)
                            .expect("a drawn count fits the census's one-byte field"),
                    });
                }
            }
        }

        let terms = RecordingTerms {
            selection: SelectionTermsDigest::of(&selection_terms()),
            kept_loci: CensusLociDigester::new().finish(),
            ssr_stratum_counts: Default::default(),
            read_cap: ReadCap(100),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
            depth_cap: DepthCap::new(124),
        };
        let records = (0..samples)
            .map(|s| {
                SampleCensusEvidence::resident(
                    format!("s{s}"),
                    terms.clone(),
                    BTreeMap::from([(
                        SectionKey::Generic(ReadGroupId(0)),
                        Section::Generic(GenericEvidence::from_parts(
                            std::mem::replace(&mut depth[s], PackedDepthCodes::never_walked(0)),
                            std::mem::take(&mut sparse[s]),
                        )),
                    )]),
                )
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
}

/// The whole fit, against cohorts whose truth is known.
///
/// **A second test module and not a second file**: the generator these tests draw from is
/// [`bench_fixtures`], which has to be compiled outside `cfg(test)` so the benchmark can reach
/// it, and a module cannot be half inside a `cfg(test)` one. The tests themselves are
/// unchanged.
#[cfg(test)]
mod whole_fit_tests {
    use super::bench_fixtures::{as_cohort, draw_cohort, draw_cohort_with_duplications};
    use super::*;
    use crate::ng::parameter_estimation::joint::census::DepthCap;

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
            // **This test is about the ordinary-position estimator**, and ten samples is well
            // under the twenty-five the duplicated class needs before the absence of
            // non-reference homozygotes means anything. Left on, it claims 1 position in 217 of
            // a cohort that has none. What the class costs and buys is measured in
            // `examples/ng_joint_duplicated_in_fit.rs`, at the panel sizes it is meant for.
            duplicated_positions: false,
            starting_points: vec![StartingPoint {
                clean: 0.006,
                noisy: 0.12,
                noisy_share: 0.05,
                p_invariant: 0.8,
                p_fixed_alt: 0.02,
                a: 1.0,
                b: 1.0,
                duplicated_share: 0.001,
                carrier_a: 1.2,
                carrier_b: 9.5,
            }],
            ..JointFitConfig::default()
        };
        let fit = fit_jointly(&mut as_cohort(&cohort.samples), &config).expect("the cohort pools");

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

    /// **A stretch some samples carry twice must not be read as heterozygosity.**
    ///
    /// Where a plant holds two copies of a stretch the reference holds once, both copies' reads
    /// land at the same place and about half disagree with the reference — which is what a
    /// heterozygote looks like. A fit with nowhere else to put it books it as one, and the
    /// homozygote excess is what pays: on a fifty-sample selfing panel at three reads it reads
    /// 0.42 where the truth is 0.60. The evidence that separates the two is that **a
    /// duplication leaves nobody homozygous for the non-reference allele** where a real variant
    /// at that frequency leaves a quarter of the panel there, and that is a statement about the
    /// cohort rather than about depth.
    ///
    /// Measured across panel sizes in `examples/ng_joint_duplicated_in_fit.rs`; this asserts the
    /// direction and the size at the one panel a unit test can afford.
    #[test]
    fn a_stretch_some_samples_carry_twice_is_not_read_as_heterozygosity() {
        let density = FrequencyDensity {
            p_invariant: 0.95,
            p_fixed_alt: 0.002,
            a: 0.5,
            b: 2.0,
        };
        let cohort = draw_cohort_with_duplications(
            30,
            4_000,
            3.0,
            (0.002, 0.05, 0.01),
            density,
            0.6,
            0.004,
            0x51ED_2709,
        );
        let drawn: f64 = cohort.heterozygous.iter().sum::<f64>() / cohort.heterozygous.len() as f64;
        let fit_with = |class: bool| {
            let config = JointFitConfig {
                quadrature_nodes: 12,
                max_passes: 120,
                duplicated_positions: class,
                ..JointFitConfig::default()
            };
            let fit =
                fit_jointly(&mut as_cohort(&cohort.samples), &config).expect("the cohort pools");
            let het: f64 = cohort
                .samples
                .iter()
                .map(|sample| fit.rates[&sample.sample].value.heterozygous)
                .sum::<f64>()
                / cohort.samples.len() as f64;
            let excess: f64 = cohort
                .samples
                .iter()
                .map(|sample| fit.hom_excess[&sample.sample].value.get())
                .sum::<f64>()
                / cohort.samples.len() as f64;
            (het, excess)
        };
        let (without, excess_without) = fit_with(false);
        let (with, excess_with) = fit_with(true);
        eprintln!(
            "drawn heterozygosity {drawn:.5}; without the class {without:.5} (excess \
             {excess_without:.3}), with it {with:.5} (excess {excess_with:.3}); drawn excess 0.6"
        );
        assert!(
            without / drawn > 1.15,
            "without the class heterozygosity came back at {without} against a drawn {drawn}, so \
             this cohort does not carry the artefact the test is about"
        );
        assert!(
            (with / drawn - 1.0).abs() < 0.10,
            "with the class heterozygosity came back at {with} against a drawn {drawn}"
        );
        // **Heterozygosity is what the assertion rests on**, because it is what the class
        // protects and the one the panel is large enough to measure. The homozygote excess is
        // printed rather than asserted: at thirty samples and four thousand positions it comes
        // back 0.657 with the class and 0.547 without, either side of the drawn 0.6, and the
        // difference between those two is smaller than the scatter between draws.
        assert!(
            excess_with > 0.45 && excess_with < 0.75,
            "the homozygote excess came back at {excess_with} against a drawn 0.6"
        );
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
                duplicated_share: 0.001,
                carrier_a: 1.2,
                carrier_b: 9.5,
            }],
            ..JointFitConfig::default()
        };
        let fit = fit_jointly(&mut as_cohort(&cohort.samples), &config).expect("the cohort pools");
        for sample in cohort.samples.iter() {
            let excess = fit.hom_excess[&sample.sample].value.get();
            // **Ten samples over three thousand positions is a small cohort and the invented
            // excess is its noise floor, not a bias.** Seven of the ten come back at or near
            // zero and the mean over all ten is 0.022, whether the depth is read as the range
            // it stands for or as the middle of that range; what differs is the worst single
            // sample, 0.070 reading the middle against 0.086 reading the range. The bound is
            // set above that rather than between the two, because a bound that one draw's
            // noisiest sample crosses is a bound on the draw.
            assert!(
                excess < 0.10,
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
            &mut as_cohort(&cohort.samples),
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
        cohort.samples[1].terms.depth_cap = DepthCap::new(60);
        // **The refusal has moved to the door.** Building the cohort is what makes the check,
        // before a section is read, and its refusal becomes the fit's own error unchanged.
        let refusal = CohortCensusEvidence::new(cohort.samples)
            .expect_err("the samples did not record the same evidence");
        match JointFitError::from(refusal) {
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
