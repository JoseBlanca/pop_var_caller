//! How sharply the SNP/indel fit's log-likelihood falls away from where it stands, one position at
//! a time — the raw material of every standard error the fit reports.
//!
//! Design: `doc/devel/ng/spec/fit_precision.md` §3.2. Build order:
//! `doc/devel/implementation_plans/fit_precision.md`, Milestone A.
//!
//! # What a score is, and why one position's is enough
//!
//! A parameter's **score** at a position is the slope of that position's log-likelihood in the
//! parameter: how much more (or less) likely the position's reads become as the parameter is
//! nudged up. The fit's log-likelihood is a sum of one term per position, each computed on its
//! own, so the whole slope is the sum of the positions' slopes — and, more usefully, the positions'
//! slopes multiplied pairwise and summed estimate how sharply the log-likelihood curves near its
//! peak. That curvature is the **information**, and one over its square root is the standard
//! error (spec §3.2's "outer product of gradients").
//!
//! # How each slope is taken
//!
//! **Each slope is the posterior expectation of the slope the parameter would have if the hidden
//! states were known** — which genotype each sample holds, whether the position is mismapped,
//! what its allele frequency is (Fisher's identity). The expectation step already computes every
//! one of those posteriors, so no new pass structure is needed: only new sums over what the pass
//! has in hand.
//!
//! | parameter | its slope at one position, if the hidden states were known |
//! |---|---|
//! | the four shares | `1/share` for the class or branch the position is in, minus `1/(1 − share)` (or the share it is traded against) for the others |
//! | the Beta shapes `a`, `b` | with `ψ` the digamma function (the slope of `ln Γ`), `ln f − ψ(a) + ψ(a + b)` and `ln(1 − f) − ψ(b) + ψ(a + b)` at the position's frequency `f`; the carrier Beta's the same at its carrier frequency — **but see below for how it is integrated** |
//! | a read group's error rate | for each of the sample's reads, the slope of the log-probability of the base it showed, in the rate |
//! | a sample's homozygote excess | the slope of the log-prior of its genotype, in the excess, on the segregating branch |
//!
//! **The expectation is taken over the posterior, never plugged in at the posterior mean.** For
//! a Beta shape the two differ by the gap between `E[ln f]` and `ln E[f]`, which at three reads a
//! position is not small (spec §3.2's first trap).
//!
//! # The Beta shapes' slopes are taken on the rule the likelihood is computed on
//!
//! The pass integrates over a position's frequency on a Gauss–Jacobi rule: sixteen nodes placed
//! by the Beta itself, which is what lets it integrate a rare-allele pile-up that a grid cannot.
//! **The rule integrates the likelihood well and `ln f` badly.** `ln f` is unbounded at zero,
//! exactly where a Beta with a shape below one puts its mass, and no polynomial rule follows it
//! there. Measured on the drawn cohort of this module's tests, at `a = 0.9`: integrating
//! `ln f − ψ(a) + ψ(a + b)` on the rule gives the slope in `a` 2.9% short of the log-likelihood's
//! own slope at the sixteen nodes the fit ships with, 4.7% at 12, 1.4% at 24, 0.4% at 48 and 0.05%
//! at 160, while a finite difference of
//! the log-likelihood reads the same number, −20.9407, at every one of those node counts.
//!
//! So the slope in a shape is taken as **the derivative of the rule's own answer**: as the shape
//! moves, each node moves and each weight changes, and the integrand is read at the moved node.
//! How the nodes and weights move depends on the shapes alone, so it is computed once a pass
//! ([`RuleSlopes`]); how the integrand changes along a node's move is the derivative of each
//! sample's genotype prior in the frequency, one more sum beside the one the pass already takes.
//! It is the same quantity as the table's — the slope of the likelihood in the shape — computed
//! on the likelihood the fit actually maximises.
//!
//! # Which parameters a sample carries
//!
//! **A sample's reads are scored under its first read group's two error rates**, and under no
//! other read group's: [`one_position`](super::one_position) pools a sample's reads over its read
//! groups and reads the rate of `group_index[s][0]`. So the likelihood does not depend on the
//! rates of a sample's second and later read groups at all, and their slope is zero at every
//! position — they carry no information, and no standard error can be computed for them. What a
//! sample contributes to this module is therefore three parameters, whatever its read groups:
//! its first read group's clean and mismapped error rates, and its homozygote excess.
//!
//! # From scores to information
//!
//! The information is the sum, over positions, of each position's scores multiplied pairwise
//! ([`InformationSums`]). **Only the blocks the standard errors need are kept** (spec §3.2's block
//! approximation): the cohort's eight parameters with each other, each sample's three with each
//! other, and each sample's three with the cohort's eight. A product of two different samples'
//! scores is never formed — the full matrix at 2,169 samples would be about 450 MB, these blocks
//! are 264 bytes a sample.

use crate::float;

use super::{
    BRANCHES, BetaQuadrature, CANDIDATE_ALTERNATIVES, MAX_CANDIDATES, Parameters, PassModel,
    Scratch, candidate_read_probability, candidate_read_probability_slope, class_rate,
    expected_reference_reads, genotype_frequencies_slope_in_excess,
    genotype_frequencies_slope_in_frequency, ln_sum_exp, reference_read_probability,
    reference_read_probability_slope,
};

/// How far each shape is moved, relative to itself, when [`RuleSlopes`] differences the rule.
///
/// The rule's nodes are eigenvalues and its weights come from eigenvectors found by iteration,
/// so a step much smaller than this would difference their rounding, and one much larger would
/// difference curvature. At this step the slopes agree with a finite difference of the whole
/// log-likelihood to better than one part in a million (this module's tests).
const RULE_STEP: f64 = 1e-4;

/// How a quadrature rule's nodes and weights move as its Beta's two shapes move.
///
/// **A property of the rule and not of any position**, so it is computed once a pass, from the
/// rule at each shape nudged up and down. It records the shapes it was computed at, so a scorer
/// handed the slopes of another rule can be caught.
#[derive(Clone, Debug)]
pub(super) struct RuleSlopes {
    /// The two shapes the rule was built at.
    shapes: (f64, f64),
    /// `[shape][node]` — how far a node moves per unit of the shape; shape 0 is `a`, 1 is `b`.
    node_slope: [Vec<f64>; 2],
    /// `[shape][node]` — how the logarithm of a node's weight moves per unit of the shape.
    ln_weight_slope: [Vec<f64>; 2],
}

impl RuleSlopes {
    /// The slopes of the `count`-node rule for `Beta(a, b)`, by central differences of the rule.
    pub(super) fn of(a: f64, b: f64, count: usize) -> Self {
        let along = |step_a: f64, step_b: f64| {
            let up = BetaQuadrature::new(a + step_a, b + step_b, count);
            let down = BetaQuadrature::new(a - step_a, b - step_b, count);
            let width = 2.0 * (step_a + step_b);
            let node = up
                .nodes
                .iter()
                .zip(&down.nodes)
                .map(|(up, down)| (up - down) / width)
                .collect();
            let ln_weight = up
                .ln_weights
                .iter()
                .zip(&down.ln_weights)
                .map(|(up, down)| (up - down) / width)
                .collect();
            (node, ln_weight)
        };
        let (node_a, ln_weight_a) = along(RULE_STEP * a, 0.0);
        let (node_b, ln_weight_b) = along(0.0, RULE_STEP * b);
        Self {
            shapes: (a, b),
            node_slope: [node_a, node_b],
            ln_weight_slope: [ln_weight_a, ln_weight_b],
        }
    }

    /// Whether these are the slopes of `rule`, built at `shapes`.
    fn describes(&self, rule: &BetaQuadrature, shapes: (f64, f64)) -> bool {
        self.shapes == shapes && self.node_slope[0].len() == rule.nodes.len()
    }
}

/// **What the scorer needs that depends on the parameters alone**, built once a pass: how the two
/// rules move with their shapes (the frequency density's and, when the run fits the duplicated
/// class, the carrier Beta's), and the log of each branch's share.
#[derive(Debug)]
pub(super) struct ScoringTables {
    density: RuleSlopes,
    carrier: Option<RuleSlopes>,
    branches: BranchShares,
}

impl ScoringTables {
    /// The tables for the rules and parameters `model` carries, at the node count its rules were
    /// built with.
    pub(super) fn of(model: &PassModel<'_>) -> Self {
        let parameters = model.parameters;
        let count = model.quadrature.nodes.len();
        Self {
            density: RuleSlopes::of(parameters.density.a, parameters.density.b, count),
            carrier: parameters
                .duplicated
                .map(|d| RuleSlopes::of(d.carrier_a, d.carrier_b, count)),
            branches: BranchShares::of(model),
        }
    }
}

/// How many cohort-level parameters a score row carries, whether or not the run fits the
/// duplicated class: a run that does not has zeros in the last three.
pub(super) const COHORT_PARAMETERS: usize = 8;

/// Where each cohort-level parameter sits in a row. **Always eight slots**: with the duplicated
/// class fitted this is also the order [`Parameters::coordinates`](super::Parameters) lists the
/// cohort's parameters in; without it the last three slots stay zero here while `coordinates`
/// simply leaves them out, so the two orders then part company after slot 4.
pub(super) mod cohort {
    /// The share of positions that are mismapped.
    pub const NOISY_SHARE: usize = 0;
    /// The share of positions where the population carries only the reference base.
    pub const P_INVARIANT: usize = 1;
    /// The share where it carries only a non-reference base.
    pub const P_FIXED_ALT: usize = 2;
    /// The first shape of the Beta over the frequencies that segregate.
    pub const DENSITY_A: usize = 3;
    /// The second shape of that Beta.
    pub const DENSITY_B: usize = 4;
    /// The share of positions some sample carries more copies of than the reference.
    pub const DUPLICATED_SHARE: usize = 5;
    /// The first shape of the Beta over how much of the panel carries such a copy.
    pub const CARRIER_A: usize = 6;
    /// The second shape of that Beta.
    pub const CARRIER_B: usize = 7;
}

/// How many parameters a sample carries of its own.
pub(super) const SAMPLE_PARAMETERS: usize = 3;

/// Where each of a sample's own parameters sits in its row.
pub(super) mod sample {
    /// Its first read group's error rate at an ordinary position.
    pub const CLEAN_ERROR_RATE: usize = 0;
    /// …and at a mismapped one.
    pub const NOISY_ERROR_RATE: usize = 1;
    /// How much less heterozygous the sample is than random mating predicts.
    pub const HOMOZYGOTE_EXCESS: usize = 2;
}

/// Every parameter's slope at one position: the cohort's eight, and each sample's three.
///
/// **One value reused over every position**: [`score_position`] clears and refills it, so a
/// pass allocates it once per chunk. It also holds the table of each read's slope the error
/// rates' scores are built from, for the same reason.
#[derive(Clone, Debug)]
pub(super) struct PositionScores {
    /// The cohort's parameters, indexed by [`cohort`].
    pub cohort: [f64; COHORT_PARAMETERS],
    /// Per sample, in the order the fit iterates samples, indexed by [`sample`].
    pub samples: Vec<[f64; SAMPLE_PARAMETERS]>,
    /// `[class][candidate][sample][genotype]` — the slope of the log-probability of one sample's
    /// reads in its error rate, if the sample held that genotype with that candidate allele.
    read_slopes: Vec<f64>,
}

impl PositionScores {
    pub(super) fn new(samples: usize) -> Self {
        Self {
            cohort: [0.0; COHORT_PARAMETERS],
            samples: vec![[0.0; SAMPLE_PARAMETERS]; samples],
            read_slopes: vec![0.0; read_slope_table_len(samples)],
        }
    }

    fn clear(&mut self) {
        self.cohort = [0.0; COHORT_PARAMETERS];
        self.samples.fill([0.0; SAMPLE_PARAMETERS]);
    }

    /// The slope of sample `s`'s reads under genotype `genotype` with the position's
    /// `candidate`-th candidate allele, in noise class `class`'s error rate.
    fn read_slope(&self, class: usize, candidate: usize, s: usize, genotype: usize) -> f64 {
        self.read_slopes[read_slope_index(self.samples.len(), class, candidate, s, genotype)]
    }
}

/// How many entries the read-slope table holds for `samples` samples.
const fn read_slope_table_len(samples: usize) -> usize {
    2 * MAX_CANDIDATES * samples * 3
}

/// Where the slope for (class, candidate, sample, genotype) sits in the read-slope table.
const fn read_slope_index(
    samples: usize,
    class: usize,
    candidate: usize,
    s: usize,
    genotype: usize,
) -> usize {
    ((class * MAX_CANDIDATES + candidate) * samples + s) * 3 + genotype
}

/// Every parameter's slope at the position [`one_position`](super::one_position) has just
/// scored, written into `into`.
///
/// **Reads the scratch `one_position` left behind and changes nothing in it**, so calling it or
/// not calling it cannot move a fitted number; `model` must be the one `one_position` was given,
/// and `tables` built from it. A position whose likelihood underflowed carries no slope
/// at all: it contributes nothing to any parameter, exactly as it contributes nothing to the
/// counts.
///
/// **No posterior is dropped for being small.** When the pass credits each branch and node its
/// share of the position's counts, it skips any below one part in 10¹², which costs those counts
/// nothing; a slope divides by a share that can itself be 10⁻¹², so here every posterior is kept.
/// Each is formed as one exponential of a difference of logarithms, which never divides by a
/// share at all.
pub(super) fn score_position(
    scratch: &Scratch,
    model: &PassModel<'_>,
    tables: &ScoringTables,
    coverage_odds: &[f64],
    into: &mut PositionScores,
) {
    let parameters = model.parameters;
    debug_assert!(
        tables.density.describes(
            &model.quadrature,
            (parameters.density.a, parameters.density.b)
        ),
        "the density rule's slopes were built for another rule"
    );
    debug_assert_eq!(
        model.carrier.is_some(),
        tables.carrier.is_some(),
        "the carrier rule and its slopes must come together"
    );
    into.clear();
    if !scratch.position_ln.is_finite() {
        return;
    }
    fill_read_slopes(scratch, model, into);
    let branches = &tables.branches;
    for class in 0..2 {
        let class_terms = NoiseClassTerms::of(scratch, branches, parameters, class);
        score_the_shares(scratch, branches, &class_terms, into);
        score_the_invariant_and_fixed_branches(scratch, branches, &class_terms, into);
        score_the_segregating_branch(
            scratch,
            model,
            &tables.density,
            branches,
            &class_terms,
            into,
        );
        if let (Some(rule), Some(rule_slopes)) = (model.carrier.as_ref(), &tables.carrier) {
            debug_assert!(
                parameters
                    .duplicated
                    .is_some_and(|d| rule_slopes.describes(rule, (d.carrier_a, d.carrier_b))),
                "the carrier rule's slopes were built for another rule"
            );
            score_the_duplicated_branch(
                scratch,
                rule,
                rule_slopes,
                coverage_odds,
                branches,
                &class_terms,
                into,
            );
        }
    }
}

/// The log of each branch's share, the same for both noise classes and every position of a pass.
#[derive(Debug)]
struct BranchShares {
    /// Whether the run fits the duplicated class.
    duplicated_fitted: bool,
    ln_ordinary: f64,
    ln_invariant: f64,
    ln_fixed_alt: f64,
    ln_segregating: f64,
    ln_duplicated: f64,
    ln_three: f64,
}

impl BranchShares {
    fn of(model: &PassModel<'_>) -> Self {
        let ln = |p: f64| float::ln(p.max(f64::MIN_POSITIVE));
        let density = &model.parameters.density;
        let duplicated_share = model.parameters.duplicated.map_or(0.0, |d| d.share);
        Self {
            duplicated_fitted: model.carrier.is_some(),
            ln_ordinary: ln(1.0 - duplicated_share),
            ln_invariant: ln(density.p_invariant),
            ln_fixed_alt: ln(density.p_fixed_alt),
            ln_segregating: ln(density.p_segregating()),
            ln_duplicated: ln(duplicated_share),
            ln_three: float::ln(CANDIDATE_ALTERNATIVES as f64),
        }
    }
}

/// One noise class at this position: its index, the log of its share, and the log of any
/// ordinary branch's posterior with that branch's own share divided back out.
struct NoiseClassTerms {
    class: usize,
    ln_class: f64,
    /// `ln(class share) + ln(1 − duplicated share) − ln(position likelihood)`: added to a branch's
    /// log-likelihood it gives that branch's posterior over its own share — exactly the slope of
    /// the position's log-likelihood in that share.
    ln_posterior_per_share: f64,
}

impl NoiseClassTerms {
    fn of(
        scratch: &Scratch,
        branches: &BranchShares,
        parameters: &Parameters,
        class: usize,
    ) -> Self {
        let share = if class == 0 {
            1.0 - parameters.noisy_share
        } else {
            parameters.noisy_share
        };
        let ln_class = float::ln(share.max(f64::MIN_POSITIVE));
        Self {
            class,
            ln_class,
            ln_posterior_per_share: ln_class + branches.ln_ordinary - scratch.position_ln,
        }
    }
}

/// Each read's slope in its error rate, for every (class, candidate, sample, genotype), into
/// `into.read_slopes`.
///
/// The slope a genotype's reads contribute is
/// `candidate reads · candidate[j] + other reads / ε + reference[j]`: the candidate allele's reads
/// each move by the slope of their log-probability, the reads on the two bases that are neither
/// move by `1/ε` (each is `ε/3`), and the reference reads move by their expected count times their
/// own log-probability's slope — the expected count, because the depth a stored code stands for is
/// a range and the likelihood sums over it. Reads held out of the model (indels, spanning
/// deletions, `N`) do not move at all.
#[allow(
    clippy::needless_range_loop,
    reason = "the sample index addresses the evidence, the read-group map and the table at once"
)]
fn fill_read_slopes(scratch: &Scratch, model: &PassModel<'_>, into: &mut PositionScores) {
    let samples = scratch.samples;
    for s in 0..samples {
        let sample = &scratch.evidence.samples[s];
        let weights = scratch.evidence.weights_of(s);
        for class in 0..2 {
            let rate = class_rate(model.parameters, class, model.group_index[s][0]);
            let per_other_read = 1.0 / rate;
            let mut candidate = [0.0; 3];
            let mut reference = [0.0; 3];
            for j in 0..3_u8 {
                let on_candidate = candidate_read_probability(j, model.ploidy, rate);
                candidate[usize::from(j)] = candidate_read_probability_slope(j, model.ploidy)
                    / on_candidate.max(f64::MIN_POSITIVE);
                let on_reference = reference_read_probability(j, model.ploidy, rate);
                let expected = expected_reference_reads(sample, weights, on_reference);
                reference[usize::from(j)] = expected
                    * reference_read_probability_slope(j, model.ploidy)
                    / on_reference.max(f64::MIN_POSITIVE);
            }
            for (position_candidate, &allele) in scratch.candidates.iter().enumerate() {
                let candidate_reads = sample.on[allele];
                let other_reads = sample.non_reference() - candidate_reads - sample.on[4];
                for j in 0..3 {
                    into.read_slopes[read_slope_index(samples, class, position_candidate, s, j)] =
                        candidate_reads * candidate[j]
                            + other_reads * per_other_read
                            + reference[j];
                }
            }
        }
    }
}

/// The class share's slope, and the three shares the class's branches are weighted by.
fn score_the_shares(
    scratch: &Scratch,
    branches: &BranchShares,
    class_terms: &NoiseClassTerms,
    into: &mut PositionScores,
) {
    let class = class_terms.class;
    let position_ln = scratch.position_ln;
    let within = ln_sum_exp(&scratch.branch_ln[class * BRANCHES..][..BRANCHES]);
    let slope_of_class = float::exp(within - position_ln);
    into.cohort[cohort::NOISY_SHARE] += if class == 1 {
        slope_of_class
    } else {
        -slope_of_class
    };
    let on_segregating_per_share =
        float::exp(class_terms.ln_posterior_per_share + scratch.segregating_ln[class]);
    into.cohort[cohort::P_INVARIANT] +=
        float::exp(class_terms.ln_posterior_per_share + scratch.invariant_ln[class])
            - on_segregating_per_share;
    into.cohort[cohort::P_FIXED_ALT] +=
        float::exp(class_terms.ln_posterior_per_share + scratch.fixed_alt_ln[class])
            - on_segregating_per_share;
    if branches.duplicated_fitted {
        let ordinary = ln_sum_exp(&[
            branches.ln_invariant + scratch.invariant_ln[class],
            branches.ln_fixed_alt + scratch.fixed_alt_ln[class],
            branches.ln_segregating + scratch.segregating_ln[class],
        ]);
        into.cohort[cohort::DUPLICATED_SHARE] +=
            float::exp(class_terms.ln_class + scratch.duplicated_ln[class] - position_ln)
                - float::exp(class_terms.ln_class + ordinary - position_ln);
    }
}

/// The error rates' slopes on the two branches where every sample holds the same genotype: all
/// homozygous for the reference, or all for a candidate allele.
fn score_the_invariant_and_fixed_branches(
    scratch: &Scratch,
    branches: &BranchShares,
    class_terms: &NoiseClassTerms,
    into: &mut PositionScores,
) {
    // A class's error rate is that class's slot in the sample's row: the clean rate is class
    // zero's, the mismapped rate class one's.
    const _: () = assert!(sample::CLEAN_ERROR_RATE == 0 && sample::NOISY_ERROR_RATE == 1);
    let class = class_terms.class;
    let invariant = float::exp(
        class_terms.ln_posterior_per_share + branches.ln_invariant + scratch.invariant_ln[class],
    );
    // At genotype 0 no read's slope depends on which allele is the candidate — the candidate's
    // `(1/3)/(ε/3)` is the other bases' `1/ε` — so the first candidate stands for all (in exact
    // arithmetic; another candidate's table entry can differ in its last bits).
    for s in 0..scratch.samples {
        let slope = into.read_slope(class, 0, s, 0);
        into.samples[s][class] += invariant * slope;
    }
    for candidate in 0..scratch.candidates.len() {
        let fixed = float::exp(
            class_terms.ln_posterior_per_share
                + branches.ln_fixed_alt
                + scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                + float::ln(scratch.multiplicity[candidate])
                - branches.ln_three,
        );
        for s in 0..scratch.samples {
            let slope = into.read_slope(class, candidate, s, 2);
            into.samples[s][class] += fixed * slope;
        }
    }
}

/// The segregating branch, node by node: each sample's error rate and homozygote excess, and the
/// density's two shapes.
fn score_the_segregating_branch(
    scratch: &Scratch,
    model: &PassModel<'_>,
    rule_slopes: &RuleSlopes,
    branches: &BranchShares,
    class_terms: &NoiseClassTerms,
    into: &mut PositionScores,
) {
    let class = class_terms.class;
    let quadrature = &model.quadrature;
    let samples = scratch.samples;
    for candidate in 0..scratch.candidates.len() {
        let multiplicity = float::ln(scratch.multiplicity[candidate]);
        for node in 0..scratch.nodes {
            let at_node = float::exp(
                class_terms.ln_posterior_per_share
                    + branches.ln_segregating
                    + scratch.node_ln[scratch.node_at(class, candidate, node)]
                    + multiplicity
                    - branches.ln_three,
            );
            if at_node == 0.0 {
                continue;
            }
            let f = quadrature.nodes[node];
            let prior_along_excess = genotype_frequencies_slope_in_excess(f);
            // How the node's whole product over samples moves with the frequency — each
            // sample's term's slope in `f` over the term, summed.
            let mut along_frequency = 0.0;
            for s in 0..samples {
                let prior = &quadrature.priors[(node * samples + s) * 3..][..3];
                let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                let joint = [prior[0] * lik[0], prior[1] * lik[1], prior[2] * lik[2]];
                let total = joint[0] + joint[1] + joint[2];
                if total <= 0.0 {
                    continue;
                }
                let excess_slope = (prior_along_excess[0] * lik[0]
                    + prior_along_excess[1] * lik[1]
                    + prior_along_excess[2] * lik[2])
                    / total;
                into.samples[s][sample::HOMOZYGOTE_EXCESS] += at_node * excess_slope;
                let prior_along_f =
                    genotype_frequencies_slope_in_frequency(f, model.parameters.hom_excess[s]);
                along_frequency += (prior_along_f[0] * lik[0]
                    + prior_along_f[1] * lik[1]
                    + prior_along_f[2] * lik[2])
                    / total;
                let mut reads = 0.0;
                for (j, weight) in joint.iter().enumerate() {
                    reads += weight * into.read_slope(class, candidate, s, j);
                }
                into.samples[s][class] += at_node * reads / total;
            }
            for (shape, slot) in [cohort::DENSITY_A, cohort::DENSITY_B]
                .into_iter()
                .enumerate()
            {
                into.cohort[slot] += at_node
                    * (rule_slopes.ln_weight_slope[shape][node]
                        + along_frequency * rule_slopes.node_slope[shape][node]);
            }
        }
    }
}

/// The duplicated branch, node by node over the carrier frequency: each sample's error rate and
/// the carrier Beta's two shapes.
fn score_the_duplicated_branch(
    scratch: &Scratch,
    carrier: &BetaQuadrature,
    rule_slopes: &RuleSlopes,
    coverage_odds: &[f64],
    branches: &BranchShares,
    class_terms: &NoiseClassTerms,
    into: &mut PositionScores,
) {
    let class = class_terms.class;
    for candidate in 0..scratch.candidates.len() {
        let multiplicity = float::ln(scratch.multiplicity[candidate]);
        for node in 0..carrier.nodes.len() {
            let at_node = float::exp(
                class_terms.ln_class + branches.ln_duplicated - scratch.position_ln
                    + scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                    + multiplicity
                    - branches.ln_three,
            );
            if at_node == 0.0 {
                continue;
            }
            let q = carrier.nodes[node];
            let mut along_frequency = 0.0;
            for s in 0..scratch.samples {
                let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                let joint = [(1.0 - q) * lik[0], q * odds * lik[1]];
                let total = joint[0] + joint[1];
                if total <= 0.0 {
                    continue;
                }
                along_frequency += (odds * lik[1] - lik[0]) / total;
                let reads = joint[0] * into.read_slope(class, candidate, s, 0)
                    + joint[1] * into.read_slope(class, candidate, s, 1);
                into.samples[s][class] += at_node * reads / total;
            }
            for (shape, slot) in [cohort::CARRIER_A, cohort::CARRIER_B]
                .into_iter()
                .enumerate()
            {
                into.cohort[slot] += at_node
                    * (rule_slopes.ln_weight_slope[shape][node]
                        + along_frequency * rule_slopes.node_slope[shape][node]);
            }
        }
    }
}

/// **The information, in blocks**: each position's scores multiplied pairwise and summed over
/// positions, keeping only the products spec §3.2's block approximation reads.
///
/// Every block is row-major. **A sum over positions, like every other count of the pass**, so a
/// chunk's sums add to another chunk's, and the pass joins them in its fixed order.
#[derive(Clone, Debug)]
pub(super) struct InformationSums {
    /// The cohort's eight parameters with each other, `[row * 8 + column]`, indexed by [`cohort`].
    pub cohort: [f64; COHORT_PARAMETERS * COHORT_PARAMETERS],
    /// Per sample, its three parameters with each other, `[row * 3 + column]`, indexed by
    /// [`sample`].
    pub sample_blocks: Vec<[f64; SAMPLE_PARAMETERS * SAMPLE_PARAMETERS]>,
    /// Per sample, its three parameters (rows) with the cohort's eight (columns),
    /// `[row * 8 + column]`.
    pub sample_cohort_blocks: Vec<[f64; SAMPLE_PARAMETERS * COHORT_PARAMETERS]>,
}

impl InformationSums {
    pub(super) fn new(samples: usize) -> Self {
        Self {
            cohort: [0.0; COHORT_PARAMETERS * COHORT_PARAMETERS],
            sample_blocks: vec![[0.0; SAMPLE_PARAMETERS * SAMPLE_PARAMETERS]; samples],
            sample_cohort_blocks: vec![[0.0; SAMPLE_PARAMETERS * COHORT_PARAMETERS]; samples],
        }
    }

    /// Add one position's products. A sample whose three slopes are all zero there — no reads,
    /// or a position that underflowed — would add only zeros, so it is skipped: a saving in time,
    /// not in what is added.
    pub(super) fn add_position(&mut self, scores: &PositionScores) {
        let cohort = &scores.cohort;
        for (row, &left) in cohort.iter().enumerate() {
            for (column, &right) in cohort.iter().enumerate() {
                self.cohort[row * COHORT_PARAMETERS + column] += left * right;
            }
        }
        for ((own, with_cohort), row_scores) in self
            .sample_blocks
            .iter_mut()
            .zip(self.sample_cohort_blocks.iter_mut())
            .zip(&scores.samples)
        {
            if row_scores.iter().all(|&score| score == 0.0) {
                continue;
            }
            for (row, &left) in row_scores.iter().enumerate() {
                for (column, &right) in row_scores.iter().enumerate() {
                    own[row * SAMPLE_PARAMETERS + column] += left * right;
                }
                for (column, &right) in cohort.iter().enumerate() {
                    with_cohort[row * COHORT_PARAMETERS + column] += left * right;
                }
            }
        }
    }

    /// Add another chunk's sums.
    pub(super) fn absorb(&mut self, other: &Self) {
        debug_assert_eq!(
            self.sample_blocks.len(),
            other.sample_blocks.len(),
            "two chunks of one pass hold the same samples"
        );
        for (into, from) in self.cohort.iter_mut().zip(&other.cohort) {
            *into += from;
        }
        for (into, from) in self.sample_blocks.iter_mut().zip(&other.sample_blocks) {
            for (into, from) in into.iter_mut().zip(from) {
                *into += from;
            }
        }
        for (into, from) in self
            .sample_cohort_blocks
            .iter_mut()
            .zip(&other.sample_cohort_blocks)
        {
            for (into, from) in into.iter_mut().zip(from) {
                *into += from;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::parameter_estimation::joint::census::{DepthCap, SampleGenericSections};
    use crate::parameter_estimation::joint::fit::bench_fixtures::{
        DrawnCohort, as_cohort, draw_cohort_with_duplications,
    };
    use crate::parameter_estimation::joint::fit::standard_errors::{
        StandardError, StandardErrors, inverse_of_positive_definite,
    };
    use crate::parameter_estimation::joint::fit::{
        DuplicatedPositions, EvidenceCursor, FrequencyDensity, JointFit, JointFitConfig,
        MAX_RECORDED_SPREAD, MIN_POSITIONS_PER_CHUNK, Parameters, PassKeeps, SampleAtPosition,
        Statistics, candidate_read_probability, digamma, expectation_pass, fit_jointly,
        genotype_frequencies, one_position,
    };

    use super::*;

    /// The tolerance every slope is held to, relative to the slope or to 10⁻³ where the slope is
    /// smaller — so a carrier shape's slope of 0.06 is held to 6 × 10⁻⁸ and not to an absolute
    /// 10⁻⁵. The largest disagreement any test here reaches is 1.7 × 10⁻⁷ (the carrier Beta's first
    /// shape); the next, 1.54 × 10⁻⁷, is a homozygote excess.
    const TOLERANCE: f64 = 1e-6;

    fn relative_disagreement(analytic: f64, numeric: f64) -> f64 {
        (analytic - numeric).abs() / numeric.abs().max(1e-3)
    }

    /// Lend a drawn cohort's sections to `f`, with what the pass needs beside them: the cap, which
    /// read group each sample's sections are, and each sample's mean depth.
    fn with_sections<R>(
        cohort: &DrawnCohort,
        cap: Option<DepthCap>,
        f: impl FnOnce(&[SampleGenericSections<'_>], DepthCap, &[Vec<usize>], &[f64]) -> R,
    ) -> R {
        let mut samples = cohort.samples.clone();
        if let Some(cap) = cap {
            for sample in &mut samples {
                sample.terms.depth_cap = cap;
            }
        }
        let mut census = as_cohort(&samples);
        let groups = census.read_groups().to_vec();
        let depth_cap = census
            .terms()
            .map_or(DepthCap::MAX, |terms| terms.depth_cap);
        let edges = JointFitConfig::default().edges;
        census
            .with_generic(&groups, |lent| {
                let group_index: Vec<Vec<usize>> = lent
                    .iter()
                    .map(|sections| {
                        sections
                            .iter()
                            .map(|(id, _)| groups.iter().position(|g| g == id).expect("listed"))
                            .collect()
                    })
                    .collect();
                let coverage = EvidenceCursor::mean_depth(lent, &edges, depth_cap);
                f(lent, depth_cap, &group_index, &coverage)
            })
            .expect("a drawn cohort lends its sections")
    }

    /// One pass over every position, in position order and on one thread: the log-likelihood,
    /// and every position's slopes.
    fn scores_at_every_position(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        config: &JointFitConfig,
        group_index: &[Vec<usize>],
        coverage: &[f64],
        parameters: &Parameters,
    ) -> (f64, Vec<PositionScores>) {
        let model = PassModel::new(parameters, config, group_index);
        let slopes = ScoringTables::of(&model);
        let positions = EvidenceCursor::position_count(lent);
        let mut cursor = EvidenceCursor::over(
            lent,
            &config.edges,
            depth_cap,
            coverage,
            config.depth_as_a_range,
            0,
            positions,
        );
        let nodes = model.quadrature.nodes.len();
        let mut scratch = Scratch::new(lent.len(), nodes);
        let mut statistics = Statistics::new(parameters.clean.len(), lent.len(), nodes);
        let mut here = PositionScores::new(lent.len());
        let mut every = Vec::with_capacity(positions);
        while cursor.next_position(&mut scratch.evidence) {
            one_position(&mut scratch, &model, &[], &mut statistics);
            score_position(&scratch, &model, &slopes, &[], &mut here);
            every.push(here.clone());
        }
        (statistics.log_likelihood, every)
    }

    /// The same pass: the log-likelihood, and every parameter's slope summed over the positions.
    fn summed_scores(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        config: &JointFitConfig,
        group_index: &[Vec<usize>],
        coverage: &[f64],
        parameters: &Parameters,
    ) -> (f64, PositionScores) {
        let (log_likelihood, every) =
            scores_at_every_position(lent, depth_cap, config, group_index, coverage, parameters);
        let mut total = PositionScores::new(lent.len());
        for here in &every {
            for (into, from) in total.cohort.iter_mut().zip(&here.cohort) {
                *into += from;
            }
            for (into, from) in total.samples.iter_mut().zip(&here.samples) {
                for (into, from) in into.iter_mut().zip(from) {
                    *into += from;
                }
            }
        }
        (log_likelihood, total)
    }

    /// The pass's own log-likelihood — the function the slopes are slopes of.
    fn log_likelihood(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        config: &JointFitConfig,
        group_index: &[Vec<usize>],
        coverage: &[f64],
        parameters: &Parameters,
    ) -> f64 {
        expectation_pass(
            lent,
            depth_cap,
            config,
            group_index,
            coverage,
            parameters,
            PassKeeps::SUMS_ONLY,
        )
        .log_likelihood
    }

    /// A drawn cohort of four samples, and a set of parameters **deliberately away from the
    /// maximum**, so that the slopes are far from zero — from 0.35 to about 36,000 log-likelihood
    /// units a unit of the parameter, over every test here, for every one but the carrier shapes'
    /// (−0.20 and 0.063). At the
    /// maximum each is near zero, and a slope of the wrong size agrees with a finite difference of
    /// a flat function.
    ///
    /// Every sample's homozygote excess and every error rate differs from every other, so a slope
    /// credited to the wrong sample or the wrong class is a different number.
    ///
    /// **Above the census's cap the cohort is drawn with no duplicated positions**: a carrier's
    /// reads double, and at a hundred-odd reads a position that puts its stored code past the
    /// widest range the fit reserves room for.
    fn a_cohort_and_parameters_off_the_maximum(mean_depth: f64) -> (DrawnCohort, Parameters) {
        let duplicated_share = if mean_depth > 100.0 { 0.0 } else { 0.01 };
        let cohort = draw_cohort_with_duplications(
            4,
            600,
            mean_depth,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            duplicated_share,
            0x0F17_5C03_E5A1_0001,
        );
        let parameters = Parameters {
            clean: vec![0.0021, 0.0043, 0.0035, 0.0057],
            noisy: vec![0.041, 0.083, 0.062, 0.097],
            noisy_share: 0.047,
            density: FrequencyDensity {
                p_invariant: 0.83,
                p_fixed_alt: 0.023,
                a: 0.9,
                b: 1.6,
            },
            hom_excess: vec![0.12, 0.47, 0.05, 0.71],
            duplicated: Some(DuplicatedPositions {
                share: 0.006,
                carrier_a: 1.7,
                carrier_b: 6.3,
            }),
        };
        (cohort, parameters)
    }

    /// The central difference of the log-likelihood in coordinate `index` of
    /// [`Parameters::coordinates`], at a step of `relative` times the parameter's own value.
    #[allow(
        clippy::too_many_arguments,
        reason = "the pass's own arguments, and which coordinate to move by how much"
    )]
    fn central_difference(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        config: &JointFitConfig,
        group_index: &[Vec<usize>],
        coverage: &[f64],
        parameters: &Parameters,
        index: usize,
        relative: f64,
    ) -> f64 {
        let values: Vec<f64> = parameters.coordinates().iter().map(|c| c.value).collect();
        let step = relative * values[index].abs();
        let moved = |by: f64| {
            let mut changed = values.clone();
            changed[index] += by;
            log_likelihood(
                lent,
                depth_cap,
                config,
                group_index,
                coverage,
                &parameters.with_coordinates(&changed),
            )
        };
        (moved(step) - moved(-step)) / (2.0 * step)
    }

    /// How many of [`Parameters::coordinates`]' leading coordinates are the cohort's: eight with
    /// the duplicated class fitted, five without it.
    fn cohort_coordinates(parameters: &Parameters) -> usize {
        if parameters.duplicated.is_some() {
            COHORT_PARAMETERS
        } else {
            cohort::DUPLICATED_SHARE
        }
    }

    /// The coordinate of [`Parameters::coordinates`] that sample `s`'s own parameter `which`
    /// (one of [`sample`]) is, for the drawn cohort's shape: one read group a sample, sample `s`'s
    /// being group `s`.
    fn sample_coordinate(parameters: &Parameters, s: usize, which: usize) -> usize {
        let groups = parameters.clean.len();
        let first_group_coordinate = cohort_coordinates(parameters);
        if which == sample::HOMOZYGOTE_EXCESS {
            first_group_coordinate + 2 * groups + s
        } else {
            first_group_coordinate + 2 * s + which
        }
    }

    /// Where a coordinate of [`Parameters::coordinates`] sits in a [`PositionScores`] — the
    /// inverse of [`sample_coordinate`], and the identity on the cohort's slots.
    fn score_of(scores: &PositionScores, parameters: &Parameters, index: usize) -> f64 {
        let cohort_count = cohort_coordinates(parameters);
        if index < cohort_count {
            return scores.cohort[index];
        }
        let groups = scores.samples.len();
        let group_part = index - cohort_count;
        if group_part < 2 * groups {
            // `clean_g, noisy_g` in turn, and group `g` is sample `g`'s only one.
            return scores.samples[group_part / 2][group_part % 2];
        }
        scores.samples[group_part - 2 * groups][sample::HOMOZYGOTE_EXCESS]
    }

    /// Every sample's two error rates, named, as coordinates.
    fn every_error_rate(parameters: &Parameters) -> Vec<(usize, String)> {
        (0..parameters.clean.len())
            .flat_map(|s| {
                [
                    (
                        sample_coordinate(parameters, s, sample::CLEAN_ERROR_RATE),
                        format!("sample {s} clean rate"),
                    ),
                    (
                        sample_coordinate(parameters, s, sample::NOISY_ERROR_RATE),
                        format!("sample {s} mismapped rate"),
                    ),
                ]
            })
            .collect()
    }

    /// Check the summed slope of every named coordinate against the central difference, and
    /// return the largest relative disagreement. `fitted_duplication` false runs the fit's
    /// configuration with no duplicated class.
    fn largest_disagreement(
        mean_depth: f64,
        cap: Option<DepthCap>,
        nodes: usize,
        fitted_duplication: bool,
        coordinates: &[(usize, String)],
    ) -> f64 {
        let (cohort, mut parameters) = a_cohort_and_parameters_off_the_maximum(mean_depth);
        if !fitted_duplication {
            parameters.duplicated = None;
        }
        let config = JointFitConfig {
            quadrature_nodes: nodes,
            duplicated_positions: fitted_duplication,
            ..JointFitConfig::default()
        };
        with_sections(&cohort, cap, |lent, depth_cap, group_index, coverage| {
            let (_, scores) =
                summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
            let mut worst = 0.0_f64;
            for (index, name) in coordinates {
                let index = *index;
                let analytic = score_of(&scores, &parameters, index);
                let numeric = central_difference(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                    index,
                    1e-5,
                );
                let disagreement = relative_disagreement(analytic, numeric);
                eprintln!(
                    "{name}: summed slope {analytic:.8e}, central difference {numeric:.8e}, \
                     relative disagreement {disagreement:.2e}"
                );
                worst = worst.max(disagreement);
            }
            worst
        })
    }

    /// **The four shares' slopes are the derivative of the log-likelihood** — the mismapped share,
    /// the density's two point masses and the duplicated class's share, at parameters off the
    /// maximum where the slopes run from −168 to 240 log-likelihood units a unit of share (the
    /// duplicated share's is −13).
    ///
    /// **What this fixture cannot test** is the scorer's rule of keeping the posteriors below
    /// 10⁻¹² that the pass skips when crediting counts: every share here is at least 0.006, where
    /// such a posterior
    /// carries no weight, and near a share's bound a finite difference of the whole cohort's
    /// log-likelihood moves by less than its own rounding. Adopting the pass's skip was measured
    /// to move these slopes by about 10⁻¹¹ relative or not at all.
    #[test]
    fn the_shares_slopes_are_the_derivative_of_the_log_likelihood() {
        let worst = largest_disagreement(
            8.0,
            None,
            12,
            true,
            &[
                (cohort::NOISY_SHARE, "mismapped share".to_owned()),
                (cohort::P_INVARIANT, "invariant share".to_owned()),
                (cohort::P_FIXED_ALT, "fixed non-reference share".to_owned()),
                (cohort::DUPLICATED_SHARE, "duplicated share".to_owned()),
            ],
        );
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );
    }

    /// **Without the duplicated class every slope still holds**, and the class's three slots stay
    /// zero — the configuration a run takes when it switches the class off.
    #[test]
    fn the_slopes_hold_with_the_duplicated_class_off() {
        let (_, mut parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        parameters.duplicated = None;
        let mut coordinates = vec![
            (cohort::NOISY_SHARE, "mismapped share".to_owned()),
            (cohort::P_INVARIANT, "invariant share".to_owned()),
            (cohort::P_FIXED_ALT, "fixed non-reference share".to_owned()),
            (cohort::DENSITY_A, "density shape a".to_owned()),
            (cohort::DENSITY_B, "density shape b".to_owned()),
        ];
        coordinates.extend(every_error_rate(&parameters));
        coordinates.extend((0..4).map(|s| {
            (
                sample_coordinate(&parameters, s, sample::HOMOZYGOTE_EXCESS),
                format!("sample {s} homozygote excess"),
            )
        }));
        let worst = largest_disagreement(8.0, None, 12, false, &coordinates);
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );

        let (cohort, _) = a_cohort_and_parameters_off_the_maximum(8.0);
        let config = JointFitConfig {
            quadrature_nodes: 12,
            duplicated_positions: false,
            ..JointFitConfig::default()
        };
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let (_, scores) =
                summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
            assert_eq!(
                scores.cohort[cohort::DUPLICATED_SHARE..],
                [0.0; 3],
                "a class that is not fitted has no slope"
            );
        });
    }

    /// **The Beta shapes' slopes are the derivative of the log-likelihood**, at the sixteen nodes
    /// the pass ships with — both the frequency density's shapes and the carrier Beta's. The
    /// carrier's slopes are small here (−0.20 and 0.063), which is why the tolerance is relative to
    /// the slope down to 10⁻³.
    #[test]
    fn the_beta_shapes_slopes_are_the_derivative_of_the_log_likelihood() {
        let worst = largest_disagreement(
            8.0,
            None,
            JointFitConfig::default().quadrature_nodes,
            true,
            &[
                (cohort::DENSITY_A, "density shape a".to_owned()),
                (cohort::DENSITY_B, "density shape b".to_owned()),
                (cohort::CARRIER_A, "carrier shape a".to_owned()),
                (cohort::CARRIER_B, "carrier shape b".to_owned()),
            ],
        );
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );
    }

    /// **Why the shapes' slopes are taken on the rule: the digamma form, integrated on it, is
    /// several percent off** — at twelve nodes and at the sixteen the pass ships with.
    ///
    /// The maximisation's own sums are that form — `Σ posterior · ln f` over the rule's nodes, less
    /// the segregating posterior times `ψ(a) − ψ(a + b)`, where `ψ` is the digamma function, the
    /// slope of `ln Γ` — so the pass's statistics give it directly. It lands more than 3% from the
    /// finite difference of the log-likelihood at twelve nodes and more than 2% at sixteen, while
    /// the slope this module takes lands within a part in a million at both. If the rule ever
    /// integrated `ln f` well, the first assertions would fail and the module doc's explanation
    /// would be out of date.
    #[test]
    fn the_digamma_form_integrated_on_the_rule_misses_the_slope_in_a_shape() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        for (nodes, digamma_at_least) in [(12, 0.03), (16, 0.02)] {
            let config = JointFitConfig {
                quadrature_nodes: nodes,
                ..JointFitConfig::default()
            };
            with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                let statistics = expectation_pass(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                    PassKeeps::SUMS_ONLY,
                );
                let (a, b) = (parameters.density.a, parameters.density.b);
                let digamma_form =
                    statistics.sum_ln_f - statistics.segregating * (digamma(a) - digamma(a + b));
                let numeric = central_difference(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                    cohort::DENSITY_A,
                    1e-5,
                );
                let (_, scores) =
                    summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
                let digamma_off = relative_disagreement(digamma_form, numeric);
                let ours_off = relative_disagreement(scores.cohort[cohort::DENSITY_A], numeric);
                eprintln!(
                    "{nodes} nodes, slope in a: central difference {numeric:.6e}; digamma form \
                     {digamma_form:.6e} ({digamma_off:.2e} off); this module's {:.6e} \
                     ({ours_off:.2e} off)",
                    scores.cohort[cohort::DENSITY_A]
                );
                assert!(
                    digamma_off > digamma_at_least,
                    "at {nodes} nodes the digamma form is only {digamma_off:.2e} off"
                );
                assert!(
                    ours_off < TOLERANCE,
                    "at {nodes} nodes this module's slope is {ours_off:.2e} off"
                );
            });
        }
    }

    /// **Each sample's two error rates' slopes are the derivative of the log-likelihood** — both
    /// classes, all four samples, each at its own rate. The per-sample clean-rate slopes are
    /// 2,134, −1,816, −277 and −2,176, so a slope credited to the wrong sample is a different
    /// number.
    #[test]
    fn the_error_rates_slopes_are_the_derivative_of_the_log_likelihood() {
        let (_, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        let worst = largest_disagreement(8.0, None, 12, true, &every_error_rate(&parameters));
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );
    }

    /// **The error rates' slopes are right where a stored depth is a range.**
    ///
    /// Under the census's shipped cap of 124 reads every stored code is one depth, so the reference
    /// reads' expected count is the count itself. Recorded under a higher cap, a code above 124
    /// stands for a range and the likelihood sums over it, and the slope must use the expected
    /// count under the range. The census ladder's first rung above 124 runs to 159, so a cap of 140
    /// makes it a range of sixteen depths — inside the widest the fit reserves room for — and at a
    /// mean of 132 reads 1,729 of the 2,400 sample-positions are ranged. Reading the reference
    /// count at the range's shallowest depth instead passes the test above and fails this one.
    #[test]
    fn the_error_rates_slopes_are_right_where_a_depth_is_a_range() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(132.0);
        let ranged = with_sections(
            &cohort,
            Some(DepthCap::new(140)),
            |lent, depth_cap, _, coverage| {
                let config = JointFitConfig::default();
                let mut cursor = EvidenceCursor::over(
                    lent,
                    &config.edges,
                    depth_cap,
                    coverage,
                    true,
                    0,
                    EvidenceCursor::position_count(lent),
                );
                let mut scratch = Scratch::new(lent.len(), 1);
                let mut ranged = 0;
                while cursor.next_position(&mut scratch.evidence) {
                    ranged += scratch
                        .evidence
                        .samples
                        .iter()
                        .filter(|sample| sample.spread > 1)
                        .count();
                }
                ranged
            },
        );
        assert!(
            ranged > 1_200,
            "the fixture must put more than half its 2,400 sample-positions in a range; {ranged} \
             were"
        );
        let worst = largest_disagreement(
            132.0,
            Some(DepthCap::new(140)),
            12,
            true,
            &every_error_rate(&parameters),
        );
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );
    }

    /// **Each sample's homozygote-excess slope is the derivative of the log-likelihood** — the
    /// four samples at four different excesses, one of them near the top of its range.
    #[test]
    fn the_homozygote_excess_slopes_are_the_derivative_of_the_log_likelihood() {
        let (_, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        let coordinates: Vec<(usize, String)> = (0..4)
            .map(|s| {
                (
                    sample_coordinate(&parameters, s, sample::HOMOZYGOTE_EXCESS),
                    format!("sample {s} homozygote excess"),
                )
            })
            .collect();
        let worst = largest_disagreement(8.0, None, 12, true, &coordinates);
        assert!(
            worst < TOLERANCE,
            "largest relative disagreement {worst:.2e}"
        );
    }

    /// One position's evidence, kept so that it can be altered and scored again.
    struct StoredPosition {
        samples: Vec<SampleAtPosition>,
        depth_weights: Vec<f64>,
        observed_alternatives: Vec<usize>,
    }

    /// The depth-8 fixture's positions, with one reference read of every fourth sample-position
    /// turned into a read held out of the model — an insertion, a deletion or a spanning deletion,
    /// which the drawn cohort never produces.
    fn positions_with_held_out_reads(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        coverage: &[f64],
    ) -> Vec<StoredPosition> {
        let config = JointFitConfig::default();
        let mut cursor = EvidenceCursor::over(
            lent,
            &config.edges,
            depth_cap,
            coverage,
            true,
            0,
            EvidenceCursor::position_count(lent),
        );
        let mut scratch = Scratch::new(lent.len(), 1);
        let (mut stored, mut index) = (Vec::new(), 0_usize);
        while cursor.next_position(&mut scratch.evidence) {
            let mut samples = scratch.evidence.samples.clone();
            for (s, sample) in samples.iter_mut().enumerate() {
                if (index + s) % 4 == 0 && sample.fewest_reference >= 1.0 && sample.spread == 1 {
                    sample.on[4] += 1.0;
                    sample.fewest_reference -= 1.0;
                }
            }
            stored.push(StoredPosition {
                samples,
                depth_weights: scratch.evidence.depth_weights.clone(),
                observed_alternatives: scratch.evidence.observed_alternatives.clone(),
            });
            index += 1;
        }
        stored
    }

    /// The stored positions' summed log-likelihood and, when asked, every slope summed over them,
    /// with each sample's coverage odds as given.
    fn over_stored(
        stored: &[StoredPosition],
        group_index: &[Vec<usize>],
        odds: &[f64],
        parameters: &Parameters,
        mut scores: Option<&mut PositionScores>,
    ) -> f64 {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        let samples = group_index.len();
        let model = PassModel::new(parameters, &config, group_index);
        let slopes = ScoringTables::of(&model);
        let nodes = model.quadrature.nodes.len();
        let mut scratch = Scratch::new(samples, nodes);
        let mut statistics = Statistics::new(parameters.clean.len(), samples, nodes);
        let mut here = PositionScores::new(samples);
        let mut total = 0.0;
        for position in stored {
            scratch.evidence.samples.clone_from(&position.samples);
            scratch
                .evidence
                .depth_weights
                .clone_from(&position.depth_weights);
            scratch
                .evidence
                .observed_alternatives
                .clone_from(&position.observed_alternatives);
            one_position(&mut scratch, &model, odds, &mut statistics);
            total += scratch.position_ln;
            if let Some(into) = scores.as_deref_mut() {
                score_position(&scratch, &model, &slopes, odds, &mut here);
                for (into, from) in into.cohort.iter_mut().zip(&here.cohort) {
                    *into += from;
                }
                for (into, from) in into.samples.iter_mut().zip(&here.samples) {
                    for (into, from) in into.iter_mut().zip(from) {
                        *into += from;
                    }
                }
            }
        }
        total
    }

    /// **The slopes hold where reads are held out of the model and where a sample's read depth
    /// argues for or against it carrying an extra copy** — two inputs the drawn cohort never
    /// produces. One reference read in every fourth sample-position becomes a held-out read, and
    /// each sample's coverage odds differ (a quarter, three, one and eight). Counting held-out
    /// reads as ordinary errors, or dropping the odds from the carrier shapes' slope, fails it.
    #[test]
    fn the_slopes_hold_with_held_out_reads_and_coverage_odds() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, coverage);
            let held_out: usize = stored
                .iter()
                .map(|position| {
                    position
                        .samples
                        .iter()
                        .filter(|sample| sample.on[4] > 0.0)
                        .count()
                })
                .sum();
            assert!(
                held_out > 400,
                "only {held_out} sample-positions carry a held-out read"
            );
            let odds = [0.25, 3.0, 1.0, 8.0];
            let mut scores = PositionScores::new(group_index.len());
            over_stored(&stored, group_index, &odds, &parameters, Some(&mut scores));
            let values: Vec<f64> = parameters.coordinates().iter().map(|c| c.value).collect();
            let mut worst = 0.0_f64;
            let coordinates = [
                cohort::CARRIER_A,
                cohort::CARRIER_B,
                cohort::DUPLICATED_SHARE,
            ]
            .into_iter()
            .chain(
                every_error_rate(&parameters)
                    .into_iter()
                    .map(|(index, _)| index),
            );
            for index in coordinates {
                let step = 1e-5 * values[index].abs();
                let at = |by: f64| {
                    let mut changed = values.clone();
                    changed[index] += by;
                    over_stored(
                        &stored,
                        group_index,
                        &odds,
                        &parameters.with_coordinates(&changed),
                        None,
                    )
                };
                let numeric = (at(step) - at(-step)) / (2.0 * step);
                worst = worst.max(relative_disagreement(
                    score_of(&scores, &parameters, index),
                    numeric,
                ));
            }
            eprintln!("held-out reads at {held_out} sample-positions; worst {worst:.2e}");
            assert!(
                worst < TOLERANCE,
                "largest relative disagreement {worst:.2e}"
            );
        });
    }

    /// **A sample's second read group has no slope, because the likelihood never reads its
    /// rates** — the reason a sample carries three parameters here whatever its read groups.
    ///
    /// Sample 0 is given a second read group at rates unlike its first's. Moving either of that
    /// group's rates leaves the log-likelihood bit-for-bit unchanged, while the sample's clean-rate
    /// slope still matches a finite difference in its first group's rate. A likelihood that
    /// pooled a sample's groups, or read the last, fails the first assertion.
    #[test]
    fn a_second_read_groups_rates_have_no_slope() {
        let (cohort, mut parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        parameters.clean.push(0.0079);
        parameters.noisy.push(0.113);
        let second_group = 4;
        with_sections(&cohort, None, |lent, depth_cap, _, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, coverage);
            let group_index = vec![vec![0, second_group], vec![1], vec![2], vec![3]];
            let mut scores = PositionScores::new(group_index.len());
            over_stored(&stored, &group_index, &[], &parameters, Some(&mut scores));
            let values: Vec<f64> = parameters.coordinates().iter().map(|c| c.value).collect();
            let at = |index: usize, by: f64| {
                let mut changed = values.clone();
                changed[index] += by;
                over_stored(
                    &stored,
                    &group_index,
                    &[],
                    &parameters.with_coordinates(&changed),
                    None,
                )
            };
            // Coordinates list each group's clean and mismapped rate in turn after the cohort's.
            for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                let index = COHORT_PARAMETERS + 2 * second_group + which;
                let step = 1e-3 * values[index];
                assert_eq!(
                    at(index, step),
                    at(index, -step),
                    "the second group's rate {which} moved the log-likelihood"
                );
            }
            let first_clean = COHORT_PARAMETERS + sample::CLEAN_ERROR_RATE;
            let step = 1e-5 * values[first_clean];
            let numeric = (at(first_clean, step) - at(first_clean, -step)) / (2.0 * step);
            let disagreement =
                relative_disagreement(scores.samples[0][sample::CLEAN_ERROR_RATE], numeric);
            assert!(
                disagreement < TOLERANCE,
                "sample 0's clean-rate slope is {disagreement:.2e} off its first group's"
            );
        });
    }

    /// The information blocks formed from every position's scores in position order, each product
    /// written out — the independent account a pass's [`InformationSums`] is checked against.
    fn information_by_hand(every: &[PositionScores]) -> InformationSums {
        let samples = every.first().map_or(0, |scores| scores.samples.len());
        let mut sums = InformationSums::new(samples);
        for here in every {
            for row in 0..COHORT_PARAMETERS {
                for column in 0..COHORT_PARAMETERS {
                    sums.cohort[row * COHORT_PARAMETERS + column] +=
                        here.cohort[row] * here.cohort[column];
                }
            }
            for s in 0..samples {
                for row in 0..SAMPLE_PARAMETERS {
                    for column in 0..SAMPLE_PARAMETERS {
                        sums.sample_blocks[s][row * SAMPLE_PARAMETERS + column] +=
                            here.samples[s][row] * here.samples[s][column];
                    }
                    for column in 0..COHORT_PARAMETERS {
                        sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column] +=
                            here.samples[s][row] * here.cohort[column];
                    }
                }
            }
        }
        sums
    }

    /// Every entry of two sets of blocks, beside the Cauchy–Schwarz bound on it — the square root
    /// of the two diagonal entries its row and column meet — which is the scale its rounding is
    /// measured against when positive and negative products cancel.
    fn entries_with_scale(pass: &InformationSums, hand: &InformationSums) -> Vec<(f64, f64, f64)> {
        let cohort_diagonal = |i: usize| hand.cohort[i * COHORT_PARAMETERS + i];
        let mut entries = Vec::new();
        for row in 0..COHORT_PARAMETERS {
            for column in 0..COHORT_PARAMETERS {
                let at = row * COHORT_PARAMETERS + column;
                let scale = (cohort_diagonal(row) * cohort_diagonal(column)).sqrt();
                entries.push((pass.cohort[at], hand.cohort[at], scale));
            }
        }
        for s in 0..hand.sample_blocks.len() {
            let own_diagonal = |i: usize| hand.sample_blocks[s][i * SAMPLE_PARAMETERS + i];
            for row in 0..SAMPLE_PARAMETERS {
                for column in 0..SAMPLE_PARAMETERS {
                    let at = row * SAMPLE_PARAMETERS + column;
                    let scale = (own_diagonal(row) * own_diagonal(column)).sqrt();
                    entries.push((pass.sample_blocks[s][at], hand.sample_blocks[s][at], scale));
                }
                for column in 0..COHORT_PARAMETERS {
                    let at = row * COHORT_PARAMETERS + column;
                    let scale = (own_diagonal(row) * cohort_diagonal(column)).sqrt();
                    entries.push((
                        pass.sample_cohort_blocks[s][at],
                        hand.sample_cohort_blocks[s][at],
                        scale,
                    ));
                }
            }
        }
        entries
    }

    /// **A pass asked for the information sums each position's scores multiplied pairwise**, in
    /// the three kinds of block the standard errors read, **and moves no count it keeps anyway**
    /// — on both ways a pass joins its chunks: in position order when it also keeps per-position
    /// lists (a run's final pass), and in the fixed halving tree otherwise.
    ///
    /// The 600 positions are cut into three chunks of at most 256, so the pass's sums are joined
    /// across chunks while the account here adds position by position: every entry agrees to 10⁻¹²
    /// of the bound its row and column put on it. A pass that paired a sample's scores with
    /// another sample's, dropped a chunk, or transposed a sample's block with the cohort fails
    /// it; the other two blocks are symmetric, so transposing them changes no bit. **It runs on the
    /// four-sample fixture and on one sample alone**, the smallest cohort the fit takes.
    #[test]
    fn a_pass_sums_the_scores_multiplied_pairwise() {
        let (four, parameters_of_four) = a_cohort_and_parameters_off_the_maximum(8.0);
        let one = draw_cohort_with_duplications(
            1,
            600,
            8.0,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_E5A1_0003,
        );
        let mut parameters_of_one = parameters_of_four.clone();
        parameters_of_one.clean.truncate(1);
        parameters_of_one.noisy.truncate(1);
        parameters_of_one.hom_excess.truncate(1);
        for (cohort, parameters) in [(&four, &parameters_of_four), (&one, &parameters_of_one)] {
            check_a_pass_against_its_scores(cohort, parameters);
        }
    }

    /// The body of [`a_pass_sums_the_scores_multiplied_pairwise`], for one cohort.
    fn check_a_pass_against_its_scores(cohort: &DrawnCohort, parameters: &Parameters) {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        with_sections(cohort, None, |lent, depth_cap, group_index, coverage| {
            let positions = EvidenceCursor::position_count(lent);
            assert!(
                positions > 2 * MIN_POSITIONS_PER_CHUNK,
                "the pass must join more than two chunks; {positions} positions"
            );
            let (_, every) = scores_at_every_position(
                lent,
                depth_cap,
                &config,
                group_index,
                coverage,
                parameters,
            );
            let by_hand = information_by_hand(&every);
            for per_position_posteriors in [false, true] {
                let pass = |information: bool| {
                    expectation_pass(
                        lent,
                        depth_cap,
                        &config,
                        group_index,
                        coverage,
                        parameters,
                        PassKeeps {
                            per_position_posteriors,
                            information,
                        },
                    )
                };
                let (with, without) = (pass(true), pass(false));
                assert!(without.information.is_none(), "nothing asked, nothing kept");
                for (name, asked, not_asked) in [
                    (
                        "log-likelihood",
                        with.log_likelihood,
                        without.log_likelihood,
                    ),
                    ("segregating mass", with.segregating, without.segregating),
                    ("Σ ln f", with.sum_ln_f, without.sum_ln_f),
                    ("duplicated mass", with.duplicated, without.duplicated),
                ] {
                    assert_eq!(
                        asked.to_bits(),
                        not_asked.to_bits(),
                        "asking for the information moved the {name}"
                    );
                }
                let summed = with.information.expect("the pass was asked for it");
                let mut worst = 0.0_f64;
                for (pass, hand, scale) in entries_with_scale(&summed, &by_hand) {
                    if scale == 0.0 {
                        assert_eq!(pass, 0.0, "an entry whose diagonal is empty must be zero");
                        continue;
                    }
                    worst = worst.max((pass - hand).abs() / scale);
                }
                eprintln!(
                    "information (per-position lists kept: {per_position_posteriors}): worst \
                     disagreement {worst:.2e} of each entry's bound"
                );
                assert!(worst < 1e-12, "worst disagreement {worst:.2e}");
                // Every sample has reads here, so every diagonal entry is positive.
                for block in &summed.sample_blocks {
                    for which in 0..SAMPLE_PARAMETERS {
                        assert!(block[which * SAMPLE_PARAMETERS + which] > 0.0);
                    }
                }
            }
        });
    }

    /// **The information is the same bits at one, four and eight threads**, on both ways a pass
    /// joins its chunks. The cohort is long enough for thirteen chunks.
    #[test]
    fn the_information_is_the_same_bits_at_any_pool_width() {
        let positions = 12 * MIN_POSITIONS_PER_CHUNK + 100;
        let cohort = draw_cohort_with_duplications(
            3,
            positions,
            8.0,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_E5A1_0002,
        );
        let (_, mut parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        parameters.clean.truncate(3);
        parameters.noisy.truncate(3);
        parameters.hom_excess.truncate(3);
        let config = JointFitConfig {
            quadrature_nodes: 8,
            ..JointFitConfig::default()
        };
        for per_position_posteriors in [false, true] {
            let at = |threads: usize| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("a pool of the asked-for width");
                with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                    let statistics = pool.install(|| {
                        expectation_pass(
                            lent,
                            depth_cap,
                            &config,
                            group_index,
                            coverage,
                            &parameters,
                            PassKeeps {
                                per_position_posteriors,
                                information: true,
                            },
                        )
                    });
                    let information = statistics
                        .information
                        .expect("the pass was asked for the information");
                    // `Debug` prints every finite `f64` as its shortest round-trip form, so two
                    // equal strings are two equal sets of bits.
                    format!("{information:?}")
                })
            };
            let at_one = at(1);
            for threads in [4, 8] {
                assert!(
                    at_one == at(threads),
                    "the information differs between one thread and {threads} (per-position \
                     lists kept: {per_position_posteriors})"
                );
            }
        }
    }

    /// **A sample with no reads carries no information**: its own block and its block with the
    /// cohort are exactly zero, so its standard errors will be reported as absent rather than
    /// computed from rounding (spec §3.2, "No information"). Sample 2's evidence is emptied at
    /// every position.
    #[test]
    fn a_sample_without_reads_carries_no_information() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        let silent = 2;
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let mut stored = positions_with_held_out_reads(lent, depth_cap, coverage);
            for position in &mut stored {
                position.samples[silent] = SampleAtPosition {
                    depth: 0.0,
                    on: [0.0; 5],
                    fewest_reference: 0.0,
                    spread: 1,
                };
                let weights = &mut position.depth_weights[silent * MAX_RECORDED_SPREAD..]
                    [..MAX_RECORDED_SPREAD];
                weights.fill(0.0);
                weights[0] = 1.0;
            }
            let model = PassModel::new(&parameters, &config, group_index);
            let slopes = ScoringTables::of(&model);
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(group_index.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), group_index.len(), nodes);
            let mut here = PositionScores::new(group_index.len());
            let mut information = InformationSums::new(group_index.len());
            for position in &stored {
                scratch.evidence.samples.clone_from(&position.samples);
                scratch
                    .evidence
                    .depth_weights
                    .clone_from(&position.depth_weights);
                scratch
                    .evidence
                    .observed_alternatives
                    .clone_from(&position.observed_alternatives);
                one_position(&mut scratch, &model, &[], &mut statistics);
                score_position(&scratch, &model, &slopes, &[], &mut here);
                information.add_position(&here);
            }
            assert_eq!(information.sample_blocks[silent], [0.0; 9]);
            assert_eq!(information.sample_cohort_blocks[silent], [0.0; 24]);
            assert!(
                information.sample_blocks[0][0] > 0.0,
                "a sample with reads does carry information"
            );
        });
    }

    /// **A scorer handed tables built for other shapes refuses them**, in the builds that carry
    /// debug checks: the tables here were built at a density shape `a` a hundredth higher than the
    /// pass's own, which would give every shape's slope from the wrong rule without failing.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "built for another rule")]
    fn tables_built_for_another_rule_are_refused() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        let mut other = parameters.clone();
        other.density.a += 0.01;
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let model = PassModel::new(&parameters, &config, group_index);
            let tables = ScoringTables::of(&PassModel::new(&other, &config, group_index));
            let mut cursor = EvidenceCursor::over(
                lent,
                &config.edges,
                depth_cap,
                coverage,
                config.depth_as_a_range,
                0,
                EvidenceCursor::position_count(lent),
            );
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(lent.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), lent.len(), nodes);
            let mut here = PositionScores::new(lent.len());
            assert!(cursor.next_position(&mut scratch.evidence));
            one_position(&mut scratch, &model, &[], &mut statistics);
            score_position(&scratch, &model, &tables, &[], &mut here);
        });
    }

    /// A drawn cohort fitted to its maximum, and the fitted parameters in the pass's own layout.
    fn fitted_parameters(cohort: &DrawnCohort, config: &JointFitConfig) -> Parameters {
        fitted_parameters_and_convergence(cohort, config).0
    }

    /// The same, and whether the fit converged rather than stopping at its pass limit.
    fn fitted_parameters_and_convergence(
        cohort: &DrawnCohort,
        config: &JointFitConfig,
    ) -> (Parameters, bool) {
        let (fit, parameters) = fitted(cohort, config);
        (parameters, fit.converged)
    }

    /// A drawn cohort's whole fit, and its fitted parameters in the pass's own layout.
    fn fitted(cohort: &DrawnCohort, config: &JointFitConfig) -> (JointFit, Parameters) {
        let mut census = as_cohort(&cohort.samples);
        let names: Vec<String> = census.sample_names().map(str::to_string).collect();
        let groups = census.read_groups().to_vec();
        let fit = fit_jointly(&mut census, config).expect("a drawn cohort fits");
        let parameters = Parameters {
            clean: groups.iter().map(|g| fit.noise[g].value.clean).collect(),
            noisy: groups.iter().map(|g| fit.noise[g].value.noisy).collect(),
            noisy_share: fit.noisy_share,
            density: fit.density.value,
            hom_excess: names
                .iter()
                .map(|name| fit.hom_excess[name].value.get())
                .collect(),
            duplicated: fit.duplicated.as_ref().map(|d| d.value),
        };
        (fit, parameters)
    }

    /// The kinds of parameter the comparisons report on.
    #[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    enum Kind {
        Shares,
        DensityShapes,
        CarrierShapes,
        CleanErrorRates,
        MismappedErrorRates,
        HomozygoteExcess,
    }

    /// Which kind each parameter of the score layout is, for reporting.
    fn kind_of(index: usize) -> Kind {
        match index {
            cohort::NOISY_SHARE
            | cohort::P_INVARIANT
            | cohort::P_FIXED_ALT
            | cohort::DUPLICATED_SHARE => Kind::Shares,
            cohort::DENSITY_A | cohort::DENSITY_B => Kind::DensityShapes,
            cohort::CARRIER_A | cohort::CARRIER_B => Kind::CarrierShapes,
            _ => match (index - COHORT_PARAMETERS) % SAMPLE_PARAMETERS {
                sample::HOMOZYGOTE_EXCESS => Kind::HomozygoteExcess,
                sample::CLEAN_ERROR_RATE => Kind::CleanErrorRates,
                _ => Kind::MismappedErrorRates,
            },
        }
    }

    /// Per parameter kind, the smallest, the median and the largest of `estimate / reference − 1`:
    /// negative where the estimate is the smaller error.
    fn by_kind(pairs: &[(usize, f64, f64)]) -> Vec<(Kind, [f64; 3], usize)> {
        let mut kinds: Vec<Kind> = pairs.iter().map(|(i, ..)| kind_of(*i)).collect();
        kinds.sort_unstable();
        kinds.dedup();
        kinds
            .into_iter()
            .map(|kind| {
                let mut gaps: Vec<f64> = pairs
                    .iter()
                    .filter(|(i, ..)| kind_of(*i) == kind)
                    .map(|(_, estimate, reference)| estimate / reference - 1.0)
                    .collect();
                gaps.sort_by(f64::total_cmp);
                let median = gaps[(gaps.len() - 1) / 2];
                (kind, [gaps[0], median, gaps[gaps.len() - 1]], gaps.len())
            })
            .collect()
    }

    /// **The full matrix** over the layout parameters `kept`: every position's scores multiplied
    /// pairwise and summed, the products between two samples' scores included — what the blocks
    /// leave out. Row-major, in `kept`'s order.
    fn full_matrix(every: &[PositionScores], kept: &[usize]) -> Vec<f64> {
        let samples = every.first().map_or(0, |scores| scores.samples.len());
        let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
        let m = kept.len();
        let mut full = vec![0.0; m * m];
        let mut flat = vec![0.0; n];
        for here in every {
            flat[..COHORT_PARAMETERS].copy_from_slice(&here.cohort);
            for (s, row) in here.samples.iter().enumerate() {
                let first = COHORT_PARAMETERS + SAMPLE_PARAMETERS * s;
                flat[first..first + SAMPLE_PARAMETERS].copy_from_slice(row);
            }
            for (p, &i) in kept.iter().enumerate() {
                for (q, &j) in kept.iter().enumerate() {
                    full[p * m + q] += flat[i] * flat[j];
                }
            }
        }
        full
    }

    /// The errors three ways, at a cohort's fitted maximum: from the blocks (what the fit reports),
    /// from the full matrix of every position's scores multiplied pairwise — including the products
    /// between two samples the blocks leave out — and by the spec's first sketch (the cohort's from
    /// its own block inverted, each sample's from `A_s − B_s C⁻¹ B_sᵀ` inverted). Returns, per kind,
    /// the blocks' and the sketch's largest and median relative gap from the full matrix.
    #[allow(
        clippy::type_complexity,
        reason = "two per-kind tables, read only by the test that prints them"
    )]
    fn errors_three_ways(
        cohort: &DrawnCohort,
        config: &JointFitConfig,
    ) -> (Vec<(Kind, [f64; 3], usize)>, Vec<(Kind, [f64; 3], usize)>) {
        let parameters = fitted_parameters(cohort, config);
        with_sections(cohort, None, |lent, depth_cap, group_index, coverage| {
            let (_, every) = scores_at_every_position(
                lent,
                depth_cap,
                config,
                group_index,
                coverage,
                &parameters,
            );
            let samples = lent.len();
            let sums = information_by_hand(&every);
            let blocks = StandardErrors::of(&sums);

            // The parameters the blocks inverted, in the layout cohort then samples: those with an
            // error and those whose error came out wider than their range, which stay in the
            // inversion so their uncertainty still widens the others'. Only the first are compared.
            let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
            let block_error = |i: usize| {
                if i < COHORT_PARAMETERS {
                    blocks.cohort[i]
                } else {
                    let at = i - COHORT_PARAMETERS;
                    blocks.samples[at / SAMPLE_PARAMETERS][at % SAMPLE_PARAMETERS]
                }
            };
            let kept: Vec<usize> = (0..n)
                .filter(|&i| {
                    matches!(
                        block_error(i),
                        StandardError::Estimated(_) | StandardError::WiderThanItsRange(_)
                    )
                })
                .collect();

            let m = kept.len();
            let full = full_matrix(&every, &kept);
            // How strongly two different samples' scores for the same kind of parameter move
            // together over positions: the mean correlation over every pair of samples.
            for (which, name) in [
                (sample::CLEAN_ERROR_RATE, "clean rates"),
                (sample::NOISY_ERROR_RATE, "mismapped rates"),
                (sample::HOMOZYGOTE_EXCESS, "homozygote excesses"),
            ] {
                let at: Vec<usize> = kept
                    .iter()
                    .enumerate()
                    .filter(|&(_, &i)| {
                        i >= COHORT_PARAMETERS
                            && (i - COHORT_PARAMETERS) % SAMPLE_PARAMETERS == which
                    })
                    .map(|(p, _)| p)
                    .collect();
                let (mut total, mut squares, mut pairs) = (0.0, 0.0, 0);
                for (a, &p) in at.iter().enumerate() {
                    for &q in &at[a + 1..] {
                        let correlation =
                            full[p * m + q] / (full[p * m + p] * full[q * m + q]).sqrt();
                        total += correlation;
                        squares += correlation * correlation;
                        pairs += 1;
                    }
                }
                if pairs > 0 {
                    eprintln!(
                        "{samples} samples: two samples' {name} scores correlate at {:+.3} on \
                         average, {:.3} root-mean-square, over {pairs} pairs",
                        total / f64::from(pairs),
                        (squares / f64::from(pairs)).sqrt()
                    );
                }
            }
            let full_inverse = inverse_of_positive_definite(&full, m).expect("full matrix inverts");
            let full_error = |p: usize| full_inverse[p * m + p].sqrt();

            // The spec's first sketch.
            let cohort_kept: Vec<usize> = kept
                .iter()
                .copied()
                .filter(|&i| i < COHORT_PARAMETERS)
                .collect();
            let c = cohort_kept.len();
            let cohort_block: Vec<f64> = cohort_kept
                .iter()
                .flat_map(|&p| {
                    let sums = &sums;
                    cohort_kept
                        .iter()
                        .map(move |&q| sums.cohort[p * COHORT_PARAMETERS + q])
                })
                .collect();
            let cohort_inverse =
                inverse_of_positive_definite(&cohort_block, c).expect("cohort block inverts");
            let mut sketch = vec![0.0; n];
            for (p, &i) in cohort_kept.iter().enumerate() {
                sketch[i] = cohort_inverse[p * c + p].sqrt();
            }
            for s in 0..samples {
                let own: Vec<usize> = kept
                    .iter()
                    .copied()
                    .filter(|&i| {
                        i >= COHORT_PARAMETERS && (i - COHORT_PARAMETERS) / SAMPLE_PARAMETERS == s
                    })
                    .map(|i| (i - COHORT_PARAMETERS) % SAMPLE_PARAMETERS)
                    .collect();
                let k = own.len();
                let mut reduced = vec![0.0; k * k];
                for (a, &ja) in own.iter().enumerate() {
                    for (b, &jb) in own.iter().enumerate() {
                        let mut explained = 0.0;
                        for (p, &ip) in cohort_kept.iter().enumerate() {
                            for (q, &iq) in cohort_kept.iter().enumerate() {
                                explained += sums.sample_cohort_blocks[s]
                                    [ja * COHORT_PARAMETERS + ip]
                                    * cohort_inverse[p * c + q]
                                    * sums.sample_cohort_blocks[s][jb * COHORT_PARAMETERS + iq];
                            }
                        }
                        reduced[a * k + b] =
                            sums.sample_blocks[s][ja * SAMPLE_PARAMETERS + jb] - explained;
                    }
                }
                let inverse = inverse_of_positive_definite(&reduced, k).expect("sample inverts");
                for (a, &ja) in own.iter().enumerate() {
                    sketch[COHORT_PARAMETERS + SAMPLE_PARAMETERS * s + ja] =
                        inverse[a * k + a].sqrt();
                }
            }

            let compared = || {
                kept.iter()
                    .enumerate()
                    .filter_map(|(p, &i)| block_error(i).value().map(|error| (p, i, error)))
            };
            let blocks_against_full: Vec<(usize, f64, f64)> = compared()
                .map(|(p, i, error)| (i, error, full_error(p)))
                .collect();
            let sketch_against_full: Vec<(usize, f64, f64)> = compared()
                .map(|(p, i, _)| (i, sketch[i], full_error(p)))
                .collect();
            (by_kind(&blocks_against_full), by_kind(&sketch_against_full))
        })
    }

    /// **The block errors against the full matrix, on cohorts of 4 and 20 samples** (plan step A3,
    /// spec §3.6 item 2), fitted to their maximum from drawn data at eight reads a position, 3,000
    /// positions, with the duplicated class on.
    ///
    /// The blocks leave out one thing the full matrix has: the products of two samples' scores at
    /// the same position. Measured here, `blocks / full − 1`, over the parameters that have an error
    /// (the carrier Beta's shapes, and at 4 samples the duplicated share, come out wider than their
    /// ranges and are not compared, though both matrices keep them in the inversion):
    ///
    /// | kind | 4 samples | 20 samples |
    /// |---|---|---|
    /// | shares | −0.203 to +0.023 | −0.003 to 0.000 |
    /// | density shapes | +0.011 to +0.021 | −0.001 to +0.003 |
    /// | clean error rates | −0.084 to −0.041 | −0.011 to −0.004 |
    /// | homozygote excess | −0.022 to +0.089 | −0.065 to −0.026 |
    /// | mismapped error rates | −0.244 to −0.043 (median −0.166) | −0.229 to −0.093 (median −0.145) |
    ///
    /// **At 20 samples the cohort's parameters match the full matrix, and the samples' own gaps are
    /// chance.** Two samples' scores for the same kind are uncorrelated on average (between +0.000
    /// and +0.001 over 190 pairs, for each kind), but single pairs correlate by chance, and most
    /// for the kind the fewest positions inform — root-mean-square 0.020 for the clean rates, 0.065
    /// for the excesses and 0.121 for the mismapped rates — and the gaps follow that order. The
    /// review measured them shrink with the positions: the mismapped rates' median gap from −0.145
    /// at 3,000 positions to −0.014 at 30,000.
    ///
    /// **At 4 samples the gaps are not chance.** Measured in review, the mismapped rates' median gap
    /// went from −0.166 to −0.135 to −0.109 at 3,000, 30,000 and 100,000 positions, and every pair's
    /// clean-rate scores correlated at about −0.012; at 3 reads a position the clean rates' blocks sat
    /// 0.19 below the full matrix at 30,000 positions. Two samples share a position's unknown
    /// frequency and class, and with few samples that shared part is a real share of each one's
    /// information, which the blocks leave out. Against the spread of the estimates
    /// (`the_errors_mean_what_they_say`), at 4 samples and 3 reads the estimates scatter 1.23 to 1.67
    /// times the blocks' errors and 0.89 to 1.02 times the full matrix's: the blocks' are too small
    /// and the full matrix's close, on the wide side. At 4 samples and 30 reads both hold; at 20
    /// samples the two are about 10% apart, either side of right.
    ///
    /// The spec's first sketch (the cohort's errors from its own block, each sample's from
    /// `A_s − B_s C⁻¹ B_sᵀ`) is no closer on any kind and further on the cohort's (up to 0.086 at 20
    /// samples), so the blocks are what the fit reports.
    ///
    /// The assertions hold each kind to its measured range with room for another draw. They pin the
    /// comparison, not the algebra — the arrow's inverse is pinned by `standard_errors`' own tests,
    /// which a dropped term fails and this test, at these sizes, may not.
    #[test]
    fn the_block_errors_against_the_full_matrix() {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        for samples in [4, 20] {
            let cohort = draw_cohort_with_duplications(
                samples,
                3_000,
                8.0,
                (0.003, 0.06, 0.03),
                FrequencyDensity {
                    p_invariant: 0.88,
                    p_fixed_alt: 0.01,
                    a: 0.6,
                    b: 2.2,
                },
                0.3,
                0.01,
                0x0F17_5C03_E5A1_0010 + samples as u64,
            );
            let (blocks, sketch) = errors_three_ways(&cohort, &config);
            for (
                (kind, [low, median, high], count),
                (_, [sketch_low, sketch_median, sketch_high], _),
            ) in blocks.iter().zip(&sketch)
            {
                eprintln!(
                    "{samples} samples, {kind:?} ({count}): blocks / full − 1 from {low:+.3} to \
                     {high:+.3}, median {median:+.3}; spec sketch / full − 1 from \
                     {sketch_low:+.3} to {sketch_high:+.3}, median {sketch_median:+.3}"
                );
                // How far each kind may land from the full matrix: the measured range above, with
                // room for another draw.
                let allowed = match (samples, *kind) {
                    (20, Kind::Shares | Kind::DensityShapes | Kind::CarrierShapes) => 0.01,
                    (20, Kind::CleanErrorRates) => 0.03,
                    (20, Kind::HomozygoteExcess) => 0.12,
                    _ => 0.35,
                };
                assert!(
                    low.abs() < allowed && high.abs() < allowed,
                    "{samples} samples, {kind:?}: blocks / full − 1 from {low:+.3} to {high:+.3}, \
                     allowed ±{allowed}"
                );
            }
        }
    }

    /// **At one and two samples every error rate keeps its error, and no parameter gets an error
    /// made of rounding** — the smallest cohorts, fitted with the duplicated class on as runs are,
    /// at eight reads a position over 30,000 positions.
    ///
    /// Some cohort parameters cannot be told apart there (module doc of `standard_errors`); they
    /// must say so and be dropped alone. Before that rule, a failed inversion took every error with
    /// it, and a rounding remainder that passed a looser threshold gave a share an error of 13.5
    /// though a share lives in [0, 1].
    #[test]
    fn the_smallest_cohorts_keep_their_error_rates() {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        for samples in [1, 2] {
            let cohort = draw_cohort_with_duplications(
                samples,
                30_000,
                8.0,
                (0.003, 0.06, 0.03),
                FrequencyDensity {
                    p_invariant: 0.88,
                    p_fixed_alt: 0.01,
                    a: 0.6,
                    b: 2.2,
                },
                0.3,
                0.01,
                0x0F17_5C03_E5A1_0020 + samples as u64,
            );
            let parameters = fitted_parameters(&cohort, &config);
            let errors = with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                let statistics = expectation_pass(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                    PassKeeps {
                        per_position_posteriors: false,
                        information: true,
                    },
                );
                StandardErrors::of(&statistics.information.expect("asked for"))
            });
            eprintln!("{samples} sample(s): cohort {:?}", errors.cohort);
            eprintln!("{samples} sample(s): samples {:?}", errors.samples);
            for row in &errors.samples {
                for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                    assert!(
                        row[which].value().is_some(),
                        "{samples} sample(s): an error rate lost its error: {row:?}"
                    );
                }
            }
            for share in [
                cohort::NOISY_SHARE,
                cohort::P_INVARIANT,
                cohort::P_FIXED_ALT,
                cohort::DUPLICATED_SHARE,
            ] {
                if let Some(error) = errors.cohort[share].value() {
                    assert!(
                        error < 1.0,
                        "{samples} sample(s): share {share} has an error of {error}"
                    );
                }
            }
        }
    }

    /// **A fit carries the errors at the parameters it returns** (plan step A5, spec §3.3): the
    /// errors on the `JointFit` are, bit for bit, those of one more pass at its returned parameters
    /// that sums the information — so they belong to the returned values, and not to another start
    /// or an earlier pass. Four samples with the duplicated class on, as runs are, where three starts
    /// end close together but not at the same bits, and one sample. At four samples every error rate
    /// and every homozygote excess has an error; at one, the excess says it is held fixed and the
    /// two rates keep theirs.
    #[test]
    fn a_fit_carries_the_errors_at_the_parameters_it_returns() {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        for samples in [1, 4] {
            let cohort = draw_cohort_with_duplications(
                samples,
                3_000,
                8.0,
                (0.003, 0.06, 0.03),
                FrequencyDensity {
                    p_invariant: 0.88,
                    p_fixed_alt: 0.01,
                    a: 0.6,
                    b: 2.2,
                },
                0.3,
                0.01,
                0x0F17_5C03_A500_0000 + samples as u64,
            );
            let (fit, parameters) = fitted(&cohort, &config);
            let at_the_returned_parameters =
                with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                    let statistics = expectation_pass(
                        lent,
                        depth_cap,
                        &config,
                        group_index,
                        coverage,
                        &parameters,
                        PassKeeps {
                            per_position_posteriors: true,
                            information: true,
                        },
                    );
                    StandardErrors::of(&statistics.information.expect("asked for"))
                });
            // `Debug` prints every `f64` as its shortest round-trip form: equal strings are equal bits.
            assert_eq!(
                format!("{:?}", fit.standard_errors),
                format!("{at_the_returned_parameters:?}"),
                "{samples} sample(s)"
            );
            for row in &fit.standard_errors.samples {
                for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                    assert!(row[which].value().is_some(), "{samples} sample(s): {row:?}");
                }
                if samples == 1 {
                    assert_eq!(row[sample::HOMOZYGOTE_EXCESS], StandardError::HeldFixed);
                } else {
                    assert!(
                        row[sample::HOMOZYGOTE_EXCESS].value().is_some(),
                        "{samples} samples: {row:?}"
                    );
                }
            }
        }
    }

    /// **The trace files what the fit returns under the winning start, one pass after its last**:
    /// on the four-sample fixture of the test above, whose winner is start 2 of 3 (so a slip to the
    /// first or the last start is a different start), the rows after the winner's last pass are
    /// exactly the returned values and their errors, under the names the value rows use, with the
    /// final pass's log-likelihood; and the winner's last pass left the returned values too (its
    /// last accelerated step was kept, not refused at the pass limit).
    #[test]
    fn the_trace_files_the_returned_fit_after_the_winning_starts_last_pass() {
        use crate::parameter_estimation::joint::fit_trace::captured;
        let config = JointFitConfig {
            quadrature_nodes: 12,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        let cohort = draw_cohort_with_duplications(
            4,
            3_000,
            8.0,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_A500_0004,
        );
        let ((fit, parameters), rows) = captured(|| fitted(&cohort, &config));
        let is_error = |name: &str| name.starts_with("standard_error:");
        let error_rows: Vec<_> = rows.iter().filter(|row| is_error(&row.3)).collect();
        let (start, pass) = (error_rows[0].0, error_rows[0].1);
        assert_eq!(start, 2, "the fixture's winner");
        assert_eq!(pass, fit.passes + 1);
        let values_at = |at: u32| -> Vec<(String, f64)> {
            rows.iter()
                .filter(|row| (row.0, row.1) == (start, at) && !is_error(&row.3))
                .map(|row| (row.3.clone(), row.4))
                .collect()
        };
        // `Debug` prints every `f64` as its shortest round-trip form: equal strings are equal bits.
        let returned = format!("{:?}", parameters.named());
        assert_eq!(format!("{:?}", values_at(pass)), returned);
        assert_eq!(format!("{:?}", values_at(fit.passes)), returned);
        let group_index: Vec<Vec<usize>> = (0..4).map(|s| vec![s]).collect();
        let errors: Vec<(String, f64)> = error_rows
            .iter()
            .map(|row| (row.3.clone(), row.4))
            .collect();
        assert_eq!(
            format!("{errors:?}"),
            format!("{:?}", fit.standard_errors.named(&parameters, &group_index))
        );
        assert!(
            rows.iter()
                .filter(|row| (row.0, row.1) == (start, pass))
                .all(|row| row.2 == fit.log_likelihood),
            "the returned fit's rows carry the final pass's log-likelihood"
        );
    }

    /// **A fit without the duplicated class says its three have no information, and the log leaves
    /// them out** — two samples, the class off in the draw and the fit.
    #[test]
    fn a_fit_without_the_duplicated_class_says_so() {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            duplicated_positions: false,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        let cohort = draw_cohort_with_duplications(
            2,
            3_000,
            8.0,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.0,
            0x0F17_5C03_A500_0102,
        );
        let (fit, parameters) = fitted(&cohort, &config);
        assert_eq!(
            fit.standard_errors.cohort[cohort::DUPLICATED_SHARE..],
            [StandardError::NoInformation; 3]
        );
        let line = with_sections(&cohort, None, |_, _, group_index, _| {
            fit.standard_errors.described(&parameters, group_index)
        });
        assert!(
            !line.contains("duplicated") && !line.contains("carrier"),
            "{line}"
        );
    }

    /// **A fit that runs no pass still reports its errors, at its starting point** — a pass limit
    /// below one cycle's three passes, so every start returns where it began; each sample's two
    /// error rates still have an error.
    #[test]
    fn a_fit_with_no_pass_reports_errors_at_its_start() {
        let config = JointFitConfig {
            quadrature_nodes: 12,
            max_passes: 2,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        let cohort = draw_cohort_with_duplications(
            2,
            3_000,
            8.0,
            (0.003, 0.06, 0.03),
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_A500_0202,
        );
        let (fit, _) = fitted(&cohort, &config);
        assert_eq!(fit.passes, 0);
        for row in &fit.standard_errors.samples {
            for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                assert!(row[which].value().is_some(), "{row:?}");
            }
        }
    }

    /// **The full matrix pairs every parameter with every other**, two samples' included: on two
    /// positions of two samples with scores chosen by hand, a product between the samples' own
    /// parameters, one between a sample's and the cohort's, and one on the diagonal each equal their
    /// sums over the positions, in the order `kept` lists them. A full matrix that left out the
    /// products between two samples — the blocks — gives zero for the first.
    #[test]
    fn the_full_matrix_sums_the_products_between_two_samples() {
        let mut first = PositionScores::new(2);
        let mut second = PositionScores::new(2);
        for (i, score) in first.cohort.iter_mut().enumerate() {
            *score = 1.0 + i as f64;
        }
        first.samples = vec![[2.0, 3.0, 5.0], [7.0, 11.0, 13.0]];
        second.cohort = [0.5; COHORT_PARAMETERS];
        second.samples = vec![[-1.0, 4.0, 0.25], [6.0, -2.0, 3.0]];
        let sample_slot =
            |s: usize, which: usize| COHORT_PARAMETERS + SAMPLE_PARAMETERS * s + which;
        let kept = [
            cohort::DENSITY_A,
            sample_slot(0, sample::NOISY_ERROR_RATE),
            sample_slot(1, sample::CLEAN_ERROR_RATE),
        ];
        let full = full_matrix(&[first, second], &kept);
        let m = kept.len();
        // Sample 0's mismapped rate with sample 1's clean rate: 3·7 + 4·6.
        assert_eq!(full[m + 2], 45.0);
        assert_eq!(full[2 * m + 1], 45.0);
        // The density's first shape with sample 1's clean rate: 4·7 + 0.5·6.
        assert_eq!(full[2], 31.0);
        // Sample 0's mismapped rate with itself: 3² + 4².
        assert_eq!(full[m + 1], 25.0);
    }

    /// The values the coverage cohorts are drawn from. The duplicated class is not drawn: the
    /// generator floors each carrier frequency at 0.2 (every draw below it becomes exactly 0.2), which
    /// is outside the model's family, so no carrier shape would be true.
    struct TrueValues {
        noisy_share: f64,
        density: FrequencyDensity,
        clean_error_rate: f64,
        noisy_error_rate: f64,
        homozygote_excess: f64,
    }

    const TRUTH: TrueValues = TrueValues {
        noisy_share: 0.03,
        density: FrequencyDensity {
            p_invariant: 0.88,
            p_fixed_alt: 0.01,
            a: 0.6,
            b: 2.2,
        },
        clean_error_rate: 0.003,
        noisy_error_rate: 0.06,
        homozygote_excess: 0.3,
    };

    impl TrueValues {
        /// The true values as the fit's parameters for `samples` samples, one read group each.
        fn as_parameters(&self, samples: usize) -> Parameters {
            Parameters {
                clean: vec![self.clean_error_rate; samples],
                noisy: vec![self.noisy_error_rate; samples],
                noisy_share: self.noisy_share,
                density: self.density,
                hom_excess: vec![self.homozygote_excess; samples],
                duplicated: None,
            }
        }
    }

    /// The value at layout `slot` in `parameters` — the cohort's slots by name, a sample's by kind.
    /// A drawn cohort gives each sample one read group, numbered in sample order, so sample `s`'s
    /// rates are read group `s`'s. Every sample is drawn at the same values, so a slip that read one
    /// sample's estimate against another's error would barely move the coverage; drawing each at
    /// its own rates would catch it, and would redraw every cohort (step A4 review, deferred).
    fn value_at(parameters: &Parameters, slot: usize) -> f64 {
        match slot {
            cohort::NOISY_SHARE => parameters.noisy_share,
            cohort::P_INVARIANT => parameters.density.p_invariant,
            cohort::P_FIXED_ALT => parameters.density.p_fixed_alt,
            cohort::DENSITY_A => parameters.density.a,
            cohort::DENSITY_B => parameters.density.b,
            cohort::DUPLICATED_SHARE | cohort::CARRIER_A | cohort::CARRIER_B => {
                unreachable!("the duplicated class is not fitted here")
            }
            _ => {
                let at = slot - COHORT_PARAMETERS;
                let s = at / SAMPLE_PARAMETERS;
                match at % SAMPLE_PARAMETERS {
                    sample::CLEAN_ERROR_RATE => parameters.clean[s],
                    sample::NOISY_ERROR_RATE => parameters.noisy[s],
                    sample::HOMOZYGOTE_EXCESS => parameters.hom_excess[s],
                    _ => unreachable!("a sample carries three parameters"),
                }
            }
        }
    }

    /// How a set of distances, each in units of an error, fell.
    #[derive(Default, Debug, Clone, Copy)]
    struct Tally {
        /// How many distances.
        count: usize,
        /// Of those, how many were at most one, and at most two, in size.
        within_one: usize,
        within_two: usize,
        /// Their sum and their sum of squares.
        sum: f64,
        sum_of_squares: f64,
    }

    impl Tally {
        fn add(&mut self, distance: f64) {
            self.count += 1;
            self.within_one += usize::from(distance.abs() <= 1.0);
            self.within_two += usize::from(distance.abs() <= 2.0);
            self.sum += distance;
            self.sum_of_squares += distance * distance;
        }

        fn share_within_one(&self) -> f64 {
            self.within_one as f64 / self.count.max(1) as f64
        }

        fn share_within_two(&self) -> f64 {
            self.within_two as f64 / self.count.max(1) as f64
        }

        /// How far off centre the distances sit.
        fn mean(&self) -> f64 {
            self.sum / self.count.max(1) as f64
        }

        /// How widely they scatter about their mean: their standard deviation, which for a centred
        /// estimate whose error means what it says is about one.
        fn spread(&self) -> f64 {
            let mean = self.mean();
            (self.sum_of_squares / self.count.max(1) as f64 - mean * mean)
                .max(0.0)
                .sqrt()
        }
    }

    /// How the estimates of one kind of parameter fell against their errors, over many cohorts.
    #[derive(Default, Debug)]
    struct Coverage {
        /// `(estimate − truth)` in units of the blocks' error — the error the fit reports.
        blocks: Tally,
        /// The same in units of the full matrix's error, over the estimates it inverted for.
        full: Tally,
        /// `(estimate + Newton step − truth)` in units of the full matrix's error: where one Newton
        /// step — the full inverse times the log-likelihood's slopes summed over the positions —
        /// puts the likelihood's maximum, against the truth. Over the same estimates as `full`.
        at_flat_point: Tally,
        /// Parameters of this kind without a blocks' error.
        absent: usize,
    }

    impl Coverage {
        /// How far the likelihood's flat point lies from where the fit stopped, in errors, on
        /// average: the Newton step's mean, which is the two tallies' difference.
        fn flat_point(&self) -> f64 {
            self.at_flat_point.mean() - self.full.mean()
        }
    }

    /// Draw `cohorts` cohorts of `samples` samples at `mean_depth` from the true values, fit each by
    /// the fit's own rule, and tally how far each estimate lands from the truth in units of its
    /// error — the blocks' error, and the full matrix's over the same parameters — and how far the
    /// likelihood's maximum lies from it. Returns the tallies by kind, and how many fits converged.
    fn coverage_of(
        samples: usize,
        positions: usize,
        mean_depth: f64,
        cohorts: usize,
        seed: u64,
    ) -> (std::collections::BTreeMap<&'static str, Coverage>, usize) {
        let config = JointFitConfig {
            duplicated_positions: false,
            estimate_contamination: false,
            ..JointFitConfig::default()
        };
        let truth = TRUTH.as_parameters(samples);
        let mut tally: std::collections::BTreeMap<&'static str, Coverage> = Default::default();
        let mut converged = 0;
        for draw in 0..cohorts {
            let cohort = draw_cohort_with_duplications(
                samples,
                positions,
                mean_depth,
                (
                    TRUTH.clean_error_rate,
                    TRUTH.noisy_error_rate,
                    TRUTH.noisy_share,
                ),
                TRUTH.density,
                TRUTH.homozygote_excess,
                0.0,
                seed + draw as u64,
            );
            let (parameters, this_converged) = fitted_parameters_and_convergence(&cohort, &config);
            converged += usize::from(this_converged);
            with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                let (_, every) = scores_at_every_position(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                );
                let blocks = StandardErrors::of(&information_by_hand(&every));
                let n = COHORT_PARAMETERS + SAMPLE_PARAMETERS * samples;
                let block_error = |slot: usize| {
                    if slot < COHORT_PARAMETERS {
                        blocks.cohort[slot]
                    } else {
                        let at = slot - COHORT_PARAMETERS;
                        blocks.samples[at / SAMPLE_PARAMETERS][at % SAMPLE_PARAMETERS]
                    }
                };
                // The duplicated class is off: its three slots are not parameters of this fit.
                let fitted: Vec<usize> = (0..n)
                    .filter(|&slot| !(cohort::DUPLICATED_SHARE..COHORT_PARAMETERS).contains(&slot))
                    .collect();
                let kept: Vec<usize> = fitted
                    .iter()
                    .copied()
                    .filter(|&slot| {
                        matches!(
                            block_error(slot),
                            StandardError::Estimated(_) | StandardError::WiderThanItsRange(_)
                        )
                    })
                    .collect();
                let m = kept.len();
                let full_inverse = inverse_of_positive_definite(&full_matrix(&every, &kept), m);
                // The log-likelihood's slope in each kept parameter, summed over the positions.
                let slope: Vec<f64> = kept
                    .iter()
                    .map(|&slot| {
                        every
                            .iter()
                            .map(|here| {
                                if slot < COHORT_PARAMETERS {
                                    here.cohort[slot]
                                } else {
                                    let at = slot - COHORT_PARAMETERS;
                                    here.samples[at / SAMPLE_PARAMETERS][at % SAMPLE_PARAMETERS]
                                }
                            })
                            .sum()
                    })
                    .collect();
                for &slot in &fitted {
                    let entry = tally.entry(label_of(slot)).or_default();
                    let Some(block) = block_error(slot).value() else {
                        entry.absent += 1;
                        continue;
                    };
                    let off = value_at(&parameters, slot) - value_at(&truth, slot);
                    entry.blocks.add(off / block);
                    if let Some(inverse) = &full_inverse {
                        let p = kept.iter().position(|&k| k == slot).expect("kept");
                        let error = inverse[p * m + p].sqrt();
                        let step: f64 = (0..m).map(|q| inverse[p * m + q] * slope[q]).sum();
                        entry.full.add(off / error);
                        entry.at_flat_point.add((off + step) / error);
                    }
                }
            });
        }
        (tally, converged)
    }

    /// What the coverage tallies a parameter under: each cohort-level parameter by its own name,
    /// each sample-level kind pooled over the samples.
    fn label_of(slot: usize) -> &'static str {
        match slot {
            cohort::NOISY_SHARE => "mismapped share",
            cohort::P_INVARIANT => "invariant share",
            cohort::P_FIXED_ALT => "fixed share",
            cohort::DENSITY_A => "density shape a",
            cohort::DENSITY_B => "density shape b",
            _ => match kind_of(slot) {
                Kind::CleanErrorRates => "clean error rates",
                Kind::MismappedErrorRates => "mismapped error rates",
                Kind::HomozygoteExcess => "homozygote excess",
                _ => unreachable!("the duplicated class is not fitted here"),
            },
        }
    }

    /// Every kind the coverage tallies, as [`label_of`] names them, in the tally's order.
    const COVERAGE_KINDS: [&str; 8] = [
        "clean error rates",
        "density shape a",
        "density shape b",
        "fixed share",
        "homozygote excess",
        "invariant share",
        "mismapped error rates",
        "mismapped share",
    ];

    /// How many cohorts a regime draws unless `NG_FIT_PRECISION_COVERAGE_COHORTS` says otherwise;
    /// the assertions on the measured coverage apply only at this count.
    const COVERAGE_COHORTS: usize = 200;

    /// **What the coverage run must show**, so that a change to what it measures fails it rather
    /// than printing different numbers under a green result.
    ///
    /// At any cohort count: every kind is tallied, and every estimate has an error from both
    /// matrices. At the full count only — the draws are fixed by their seeds, so the numbers move only
    /// when the code does — each regime is held to what it measured, with room to spare:
    ///
    /// - at 4 and 20 samples, at most five fits in a regime stop at the pass limit;
    /// - where the errors cover (20 samples at both depths, 4 samples at 30 reads), the error rates,
    ///   the homozygote excess and the mismapped share — and at 20 samples the fixed share — land
    ///   between 0.58 and 0.78 within one of the blocks' errors and between 0.90 and 0.99 within
    ///   two;
    /// - at 4 samples and 3 reads the full matrix's errors put at least 0.08 more of the error rates
    ///   within one error than the blocks' do;
    /// - at 20 samples the Newton step carries both density shapes to within half an error of the
    ///   truth on average, which a step of the wrong sign or size does not.
    ///
    /// The 2-sample regime is printed and held only to the checks that apply at any count: 34 of its
    /// 200 fits stop at the pass limit, and neither matrix's errors cover there.
    fn check_the_coverage(
        samples: usize,
        mean_depth: f64,
        cohorts: usize,
        converged: usize,
        tally: &std::collections::BTreeMap<&'static str, Coverage>,
    ) {
        let regime = format!("{samples} samples, {mean_depth} reads");
        assert_eq!(
            tally.keys().copied().collect::<Vec<_>>(),
            COVERAGE_KINDS,
            "{regime}: the kinds tallied"
        );
        for (kind, coverage) in tally {
            assert!(
                coverage.blocks.count > 0 && coverage.full.count == coverage.blocks.count,
                "{regime}, {kind}: {} estimates with a blocks' error, {} with a full-matrix error",
                coverage.blocks.count,
                coverage.full.count
            );
        }
        if cohorts != COVERAGE_COHORTS {
            return;
        }
        assert!(
            samples == 2 || converged + 5 >= cohorts,
            "{regime}: {converged} of {cohorts} fits converged"
        );
        let depth = mean_depth as usize;
        if samples == 20 || (samples == 4 && depth == 30) {
            let mut covering = vec![
                "clean error rates",
                "mismapped error rates",
                "homozygote excess",
                "mismapped share",
            ];
            if samples == 20 {
                covering.push("fixed share");
            }
            for kind in covering {
                let blocks = &tally[kind].blocks;
                assert!(
                    (0.58..=0.78).contains(&blocks.share_within_one())
                        && (0.90..=0.99).contains(&blocks.share_within_two()),
                    "{regime}, {kind}: {:.3} within one of the blocks' errors, {:.3} within two",
                    blocks.share_within_one(),
                    blocks.share_within_two()
                );
            }
        }
        if samples == 4 && depth == 3 {
            for kind in ["clean error rates", "mismapped error rates"] {
                let (blocks, full) = (&tally[kind].blocks, &tally[kind].full);
                assert!(
                    full.share_within_one() >= blocks.share_within_one() + 0.08,
                    "{regime}, {kind}: within one error {:.3} (full) against {:.3} (blocks)",
                    full.share_within_one(),
                    blocks.share_within_one()
                );
            }
        }
        if samples == 20 {
            for kind in ["density shape a", "density shape b"] {
                let at_flat_point = tally[kind].at_flat_point.mean();
                assert!(
                    at_flat_point.abs() < 0.5,
                    "{regime}, {kind}: the Newton-corrected estimate sits {at_flat_point:+.3} errors \
                     from the truth"
                );
            }
        }
    }

    /// **The errors mean what they say** (plan step A4, spec §3.6 item 3): 200 cohorts drawn from
    /// known values in each regime, each fitted by the fit's own rule; of every estimate with an
    /// error, the share within one error of the truth and within two, which for an error that means
    /// what it says are about 68 in 100 and 95 in 100, and the spread of the distances, about one.
    /// Measured for the errors the fit reports (the blocks) and for the full matrix's over the same
    /// parameters, on 4 samples over 20,000 positions and 20 over 5,000 at 3 and at 30 reads a
    /// position, and on 2 samples over 30,000 at 3 reads. The duplicated class is off in the draw and
    /// the fit ([`TRUTH`]).
    /// Also printed, per kind: how far the log-likelihood's own flat point lies from where the fit
    /// stopped, and the coverage there — both by one Newton step.
    ///
    /// What it measured, and what it holds each regime to, are in [`check_the_coverage`] and the step
    /// A4 report, which has every number.
    ///
    /// Ignored by default: some twenty minutes in the container. Run with
    /// `scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture`;
    /// `scripts/dev.sh env NG_FIT_PRECISION_COVERAGE_COHORTS=3 cargo test …` shortens it to a minute
    /// or so, with only the checks that hold at any count.
    #[test]
    #[ignore = "a measurement over 1,000 fitted cohorts; run at a checkpoint"]
    fn the_errors_mean_what_they_say() {
        let cohorts = match std::env::var("NG_FIT_PRECISION_COVERAGE_COHORTS") {
            Err(_) => COVERAGE_COHORTS,
            Ok(value) => value
                .trim()
                .parse()
                .ok()
                .filter(|&count: &usize| count > 0)
                .unwrap_or_else(|| {
                    panic!(
                        "NG_FIT_PRECISION_COVERAGE_COHORTS={value:?} is not a positive whole number"
                    )
                }),
        };
        for (samples, positions, mean_depth) in [
            (2, 30_000, 3.0),
            (4, 20_000, 3.0),
            (4, 20_000, 30.0),
            (20, 5_000, 3.0),
            (20, 5_000, 30.0),
        ] {
            let (tally, converged) = coverage_of(
                samples,
                positions,
                mean_depth,
                cohorts,
                0x0F17_5C03_A400_0000 + (samples as u64) * 1_000 + mean_depth as u64,
            );
            eprintln!(
                "COVERAGE {samples} samples, {mean_depth} reads: {converged} of {cohorts} fits \
                     converged"
            );
            for (kind, coverage) in &tally {
                let (blocks, full, at_flat) =
                    (&coverage.blocks, &coverage.full, &coverage.at_flat_point);
                eprintln!(
                    "COVERAGE {samples} samples, {mean_depth} reads, {kind}: {} estimates ({} \
                         with a full-matrix error), {} without an error; within one error {:.3} \
                         (blocks) {:.3} (full), within two {:.3} {:.3}; mean distance {:+.3} \
                         {:+.3}; spread {:.3} {:.3}; the likelihood's flat point lies {:+.3} errors \
                         from the fit, on average; there, within one {:.3}, within two {:.3}, mean \
                         {:+.3}",
                    blocks.count,
                    full.count,
                    coverage.absent,
                    blocks.share_within_one(),
                    full.share_within_one(),
                    blocks.share_within_two(),
                    full.share_within_two(),
                    blocks.mean(),
                    full.mean(),
                    blocks.spread(),
                    full.spread(),
                    coverage.flat_point(),
                    at_flat.share_within_one(),
                    at_flat.share_within_two(),
                    at_flat.mean(),
                );
            }
            check_the_coverage(samples, mean_depth, cohorts, converged, &tally);
        }
    }

    /// **The model's slopes kept beside `fit.rs`'s functions are the derivatives of those
    /// functions** — the genotype prior in the frequency and in the excess, and a read's two
    /// probabilities in the error rate — by a central difference, over a grid of each argument.
    #[test]
    fn the_models_slopes_are_the_derivatives_of_its_functions() {
        let ploidy = JointFitConfig::default().ploidy;
        let step = 1e-6;
        for f in [0.01, 0.2, 0.5, 0.83] {
            for excess in [0.0, 0.3, 0.9] {
                let along_f = genotype_frequencies_slope_in_frequency(f, excess);
                let along_excess = genotype_frequencies_slope_in_excess(f);
                let up_f = genotype_frequencies(f + step, excess);
                let down_f = genotype_frequencies(f - step, excess);
                let up_excess = genotype_frequencies(f, excess + step);
                let down_excess = genotype_frequencies(f, excess - step);
                for j in 0..3 {
                    let numeric_f = (up_f[j] - down_f[j]) / (2.0 * step);
                    let numeric_excess = (up_excess[j] - down_excess[j]) / (2.0 * step);
                    assert!(
                        (along_f[j] - numeric_f).abs() < 1e-7,
                        "genotype {j} at f {f}, excess {excess}: {} against {numeric_f}",
                        along_f[j]
                    );
                    assert!(
                        (along_excess[j] - numeric_excess).abs() < 1e-7,
                        "genotype {j} at f {f}, excess {excess}: {} against {numeric_excess}",
                        along_excess[j]
                    );
                }
            }
        }
        for rate in [1e-4, 0.003, 0.08] {
            for copies in 0..3_u8 {
                let candidate = (candidate_read_probability(copies, ploidy, rate + step)
                    - candidate_read_probability(copies, ploidy, rate - step))
                    / (2.0 * step);
                let reference = (reference_read_probability(copies, ploidy, rate + step)
                    - reference_read_probability(copies, ploidy, rate - step))
                    / (2.0 * step);
                assert!(
                    (candidate_read_probability_slope(copies, ploidy) - candidate).abs() < 1e-7,
                    "the candidate's slope at {copies} copies"
                );
                assert!(
                    (reference_read_probability_slope(copies, ploidy) - reference).abs() < 1e-7,
                    "the reference's slope at {copies} copies"
                );
            }
        }
    }
}
