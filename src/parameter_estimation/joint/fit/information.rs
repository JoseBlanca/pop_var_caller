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
//! **Each library's reads are scored under that library's own two error rates**
//! ([`one_position`](super::one_position)); a sample's libraries share its genotype and nothing
//! else. So a sample read from `k` libraries carries `1 + 2k` parameters of its own: its
//! homozygote excess, and each library's clean and mismapped error rates.
//!
//! A library's rate's slope at a position is its own reads' slope in that rate, weighted by the
//! posterior of the genotype its sample's libraries share. A sample's row lays its parameters out
//! as [`sample`] says: its first library's two rates, its excess, then each further library's two,
//! so a sample of one library has the same three slots whatever the rest of the cohort holds.
//!
//! # From scores to information
//!
//! The information is the sum, over positions, of each position's scores multiplied pairwise
//! ([`InformationSums`]). **Only the blocks the standard errors need are kept** (spec §3.2's block
//! approximation): the cohort's eight parameters with each other, each sample's own with each
//! other, and each sample's own with the cohort's eight. A product of two different samples'
//! scores is never formed — the full matrix at 2,169 samples would be about 450 MB, these blocks
//! are 264 bytes a sample of one library and `8(n² + 8n)` bytes a sample of `n` own parameters
//! (520 at two libraries).

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
    pub(super) node_slope: [Vec<f64>; 2],
    /// `[shape][node]` — how the logarithm of a node's weight moves per unit of the shape.
    pub(super) ln_weight_slope: [Vec<f64>; 2],
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
        Self {
            density: model.density_slopes.clone(),
            carrier: model.carrier.as_ref().map(|carrier| carrier.slopes.clone()),
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

/// How many parameters a sample read from one library carries of its own: that library's two
/// error rates and the sample's homozygote excess.
pub(super) const ONE_LIBRARY_SAMPLE_PARAMETERS: usize = 3;

/// How many parameters a sample read from `libraries` libraries carries of its own: two error rates
/// a library and its homozygote excess. A sample holding no library keeps the one-library row, its
/// two rate slots at zero.
pub(super) const fn own_parameters(libraries: usize) -> usize {
    if libraries <= 1 {
        ONE_LIBRARY_SAMPLE_PARAMETERS
    } else {
        1 + 2 * libraries
    }
}

/// Where each of a sample's own parameters sits in its row: **its first library's two error rates,
/// its homozygote excess, and then each further library's two rates**, in the order the sample's
/// sections are lent. A sample read from one library has exactly the first three.
///
/// The first library comes before the excess so that a one-library sample's row is the same three
/// slots whatever the cohort's other samples hold.
pub(super) mod sample {
    /// Its first library's error rate at an ordinary position.
    pub const CLEAN_ERROR_RATE: usize = 0;
    /// …and at a mismapped one.
    pub const NOISY_ERROR_RATE: usize = 1;
    /// How much less heterozygous the sample is than random mating predicts.
    pub const HOMOZYGOTE_EXCESS: usize = 2;

    /// The slot of the rate in noise class `class` (0 ordinary, 1 mismapped) of the sample's
    /// `library`-th library, counted from 0 in the order its sections are lent.
    pub const fn rate(library: usize, class: usize) -> usize {
        if library == 0 {
            class
        } else {
            1 + 2 * library + class
        }
    }

    /// For a slot holding an error rate, its noise class (0 ordinary, 1 mismapped); `None` for the
    /// homozygote excess.
    pub const fn class_of(slot: usize) -> Option<usize> {
        match slot {
            HOMOZYGOTE_EXCESS => None,
            0 | 1 => Some(slot),
            _ => Some((slot - super::ONE_LIBRARY_SAMPLE_PARAMETERS) % 2),
        }
    }
}

/// Every parameter's slope at one position: the cohort's eight, and each sample's own.
///
/// **One value reused over every position**: [`score_position`] clears and refills it, so a
/// pass allocates it once per chunk. It also holds the table of each read's slope the error
/// rates' scores are built from, for the same reason.
#[derive(Clone, Debug)]
pub(super) struct PositionScores {
    /// The cohort's parameters, indexed by [`cohort`].
    pub cohort: [f64; COHORT_PARAMETERS],
    /// Per sample, in the order the fit iterates samples, its own parameters in [`sample`]'s
    /// layout: [`own_parameters`] of its library count.
    pub samples: Vec<Vec<f64>>,
    /// `[class][candidate][library][genotype]` — the slope of the log-probability of one library's
    /// reads in its error rate, if its sample held that genotype with that candidate allele.
    /// Libraries are indexed as the fit indexes read groups.
    read_slopes: Vec<f64>,
    /// How many read groups the fit has, which the table is laid out over.
    libraries: usize,
}

impl PositionScores {
    /// A row for each sample over its libraries (`group_index`, as the fit lends them), and a
    /// read-slope table over the fit's `libraries` read groups.
    pub(super) fn new(group_index: &[Vec<usize>], libraries: usize) -> Self {
        Self {
            cohort: [0.0; COHORT_PARAMETERS],
            samples: group_index
                .iter()
                .map(|own| vec![0.0; own_parameters(own.len())])
                .collect(),
            read_slopes: vec![0.0; read_slope_table_len(libraries)],
            libraries,
        }
    }

    fn clear(&mut self) {
        self.cohort = [0.0; COHORT_PARAMETERS];
        for row in &mut self.samples {
            row.fill(0.0);
        }
    }

    /// The slope of library `library`'s reads under genotype `genotype` with the position's
    /// `candidate`-th candidate allele, in noise class `class`'s error rate.
    fn read_slope(&self, class: usize, candidate: usize, library: usize, genotype: usize) -> f64 {
        self.read_slopes[read_slope_index(self.libraries, class, candidate, library, genotype)]
    }
}

/// How many entries the read-slope table holds for `libraries` read groups.
const fn read_slope_table_len(libraries: usize) -> usize {
    2 * MAX_CANDIDATES * libraries * 3
}

/// Where the slope for (class, candidate, library, genotype) sits in the read-slope table.
const fn read_slope_index(
    libraries: usize,
    class: usize,
    candidate: usize,
    library: usize,
    genotype: usize,
) -> usize {
    ((class * MAX_CANDIDATES + candidate) * libraries + library) * 3 + genotype
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
        score_the_invariant_and_fixed_branches(
            scratch,
            model.group_index,
            branches,
            &class_terms,
            into,
        );
        score_the_segregating_branch(
            scratch,
            model,
            &tables.density,
            branches,
            &class_terms,
            into,
        );
        if let (Some(rule), Some(rule_slopes)) = (
            model.carrier.as_ref().map(|carrier| &carrier.quadrature),
            &tables.carrier,
        ) {
            debug_assert!(
                parameters
                    .duplicated
                    .is_some_and(|d| rule_slopes.describes(rule, (d.carrier_a, d.carrier_b))),
                "the carrier rule's slopes were built for another rule"
            );
            score_the_duplicated_branch(
                scratch,
                model.group_index,
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
///
/// **One slope a library**, under that library's own rates: a sample's libraries share its
/// genotype, and each library's reads move with its own rates alone.
fn fill_read_slopes(scratch: &Scratch, model: &PassModel<'_>, into: &mut PositionScores) {
    let libraries = into.libraries;
    for &library in model.group_index.iter().flatten() {
        let reads = &scratch.evidence.libraries[library];
        let weights = scratch.evidence.weights_of(library);
        for class in 0..2 {
            let rate = class_rate(model.parameters, class, library);
            let per_other_read = 1.0 / rate;
            let mut candidate = [0.0; 3];
            let mut reference = [0.0; 3];
            for j in 0..3_u8 {
                let on_candidate = candidate_read_probability(j, model.ploidy, rate);
                candidate[usize::from(j)] = candidate_read_probability_slope(j, model.ploidy)
                    / on_candidate.max(f64::MIN_POSITIVE);
                let on_reference = reference_read_probability(j, model.ploidy, rate);
                let expected = expected_reference_reads(reads, weights, on_reference);
                reference[usize::from(j)] = expected
                    * reference_read_probability_slope(j, model.ploidy)
                    / on_reference.max(f64::MIN_POSITIVE);
            }
            for (position_candidate, &allele) in scratch.candidates.iter().enumerate() {
                let candidate_reads = reads.on[allele];
                let other_reads = reads.non_reference() - candidate_reads - reads.on[4];
                for j in 0..3 {
                    into.read_slopes
                        [read_slope_index(libraries, class, position_candidate, library, j)] =
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
/// homozygous for the reference, or all for a candidate allele. Each of a sample's libraries
/// (`group_index`) takes its own reads' slope, into its own slot of the sample's row.
fn score_the_invariant_and_fixed_branches(
    scratch: &Scratch,
    group_index: &[Vec<usize>],
    branches: &BranchShares,
    class_terms: &NoiseClassTerms,
    into: &mut PositionScores,
) {
    // A class's error rate is that class's slot in a library's pair: the clean rate is class
    // zero's, the mismapped rate class one's.
    const _: () = assert!(
        sample::rate(0, 0) == sample::CLEAN_ERROR_RATE
            && sample::rate(0, 1) == sample::NOISY_ERROR_RATE
    );
    let class = class_terms.class;
    let invariant = float::exp(
        class_terms.ln_posterior_per_share + branches.ln_invariant + scratch.invariant_ln[class],
    );
    // At genotype 0 no read's slope depends on which allele is the candidate — the candidate's
    // `(1/3)/(ε/3)` is the other bases' `1/ε` — so the first candidate stands for all (in exact
    // arithmetic; another candidate's table entry can differ in its last bits).
    for (s, own) in group_index.iter().enumerate() {
        for (section, &library) in own.iter().enumerate() {
            let slope = into.read_slope(class, 0, library, 0);
            into.samples[s][sample::rate(section, class)] += invariant * slope;
        }
    }
    for candidate in 0..scratch.candidates.len() {
        let fixed = float::exp(
            class_terms.ln_posterior_per_share
                + branches.ln_fixed_alt
                + scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                + float::ln(scratch.multiplicity[candidate])
                - branches.ln_three,
        );
        for (s, own) in group_index.iter().enumerate() {
            for (section, &library) in own.iter().enumerate() {
                let slope = into.read_slope(class, candidate, library, 2);
                into.samples[s][sample::rate(section, class)] += fixed * slope;
            }
        }
    }
}

/// The segregating branch, node by node: each library's error rate, each sample's homozygote
/// excess, and the density's two shapes.
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
                // Each library's reads move with its own rate, weighted by the genotype posterior
                // its sample's libraries share.
                for (section, &library) in model.group_index[s].iter().enumerate() {
                    let mut reads = 0.0;
                    for (j, weight) in joint.iter().enumerate() {
                        reads += weight * into.read_slope(class, candidate, library, j);
                    }
                    into.samples[s][sample::rate(section, class)] += at_node * reads / total;
                }
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

/// The duplicated branch, node by node over the carrier frequency: each library's error rate and
/// the carrier Beta's two shapes.
#[allow(
    clippy::too_many_arguments,
    reason = "the position, the carrier rule and its slopes, the coverage odds, the library map, \
              the branch's and the class's terms, and the row written"
)]
fn score_the_duplicated_branch(
    scratch: &Scratch,
    group_index: &[Vec<usize>],
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
            for (s, own) in group_index.iter().enumerate() {
                let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                let joint = [(1.0 - q) * lik[0], q * odds * lik[1]];
                let total = joint[0] + joint[1];
                if total <= 0.0 {
                    continue;
                }
                along_frequency += (odds * lik[1] - lik[0]) / total;
                for (section, &library) in own.iter().enumerate() {
                    let reads = joint[0] * into.read_slope(class, candidate, library, 0)
                        + joint[1] * into.read_slope(class, candidate, library, 1);
                    into.samples[s][sample::rate(section, class)] += at_node * reads / total;
                }
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

/// **Up to this many samples the information is also summed as the whole matrix** — every
/// parameter paired with every other, two samples' parameters included — and the standard errors
/// come from it; above, from the blocks alone (checkpoint A's decision 3, plan step A8). At 4
/// samples and 3 reads the blocks' errors were too small — the error rates and the mismapped share
/// scattered 1.23 to 1.67 times them, and 0.89 to 1.02 times the whole matrix's — while at 20
/// samples the two came within 11% of each other. What the whole matrix costs grows as the square
/// of the parameters: at 20 samples of one library each, 68 parameters, a position adds 2,346
/// products where the blocks add 724. So the parameters are limited too
/// ([`FULL_MATRIX_PARAMETERS`]).
pub(super) const FULL_MATRIX_SAMPLES: usize = 20;

/// **And only while the whole matrix has at most this many parameters** (checkpoint A′, owner,
/// 2026-09-29): 20 samples of four libraries each, whose final pass holds 18.3 MB of matrices and
/// took 1.10 to 1.12 times the blocks' pass (plan step A8's review). Samples of more libraries make
/// the matrix grow as the square of their parameters — 20 samples of 16 libraries have 668 and
/// would hold 231 MB — so a cohort above the limit takes the blocks.
pub(super) const FULL_MATRIX_PARAMETERS: usize = 188;

/// **The information, in blocks**: each position's scores multiplied pairwise and summed over
/// positions, keeping only the products spec §3.2's block approximation reads — and, for a cohort
/// of at most [`FULL_MATRIX_SAMPLES`] samples, every product ([`FullInformation`]).
///
/// Every block is row-major. **A sum over positions, like every other count of the pass**, so a
/// chunk's sums add to another chunk's, and the pass joins them in its fixed order.
#[derive(Clone, Debug)]
pub(super) struct InformationSums {
    /// The cohort's eight parameters with each other, `[row * 8 + column]`, indexed by [`cohort`].
    pub cohort: [f64; COHORT_PARAMETERS * COHORT_PARAMETERS],
    /// Per sample, its own parameters with each other, `[row * n + column]` with `n` its row's
    /// length ([`own_parameters`]), indexed by [`sample`].
    pub sample_blocks: Vec<Vec<f64>>,
    /// Per sample, its own parameters (rows) with the cohort's eight (columns),
    /// `[row * 8 + column]`.
    pub sample_cohort_blocks: Vec<Vec<f64>>,
    /// The whole matrix, when the sums were made for a cohort small enough to hold it
    /// ([`InformationSums::for_a_cohort_of`]); the standard errors are then read from it.
    pub full: Option<FullInformation>,
}

impl InformationSums {
    /// Empty sums for samples reading the libraries `group_index` lists, in blocks only.
    pub(super) fn new(group_index: &[Vec<usize>]) -> Self {
        let own = |libraries: &Vec<usize>| own_parameters(libraries.len());
        Self {
            cohort: [0.0; COHORT_PARAMETERS * COHORT_PARAMETERS],
            sample_blocks: group_index
                .iter()
                .map(|libraries| vec![0.0; own(libraries) * own(libraries)])
                .collect(),
            sample_cohort_blocks: group_index
                .iter()
                .map(|libraries| vec![0.0; own(libraries) * COHORT_PARAMETERS])
                .collect(),
            full: None,
        }
    }

    /// **What a pass sums for the standard errors**: the blocks, and the whole matrix too when the
    /// cohort has at most [`FULL_MATRIX_SAMPLES`] samples and the matrix at most
    /// [`FULL_MATRIX_PARAMETERS`] parameters.
    pub(super) fn for_a_cohort_of(group_index: &[Vec<usize>]) -> Self {
        let mut sums = Self::new(group_index);
        let parameters = COHORT_PARAMETERS
            + group_index
                .iter()
                .map(|libraries| own_parameters(libraries.len()))
                .sum::<usize>();
        if group_index.len() <= FULL_MATRIX_SAMPLES && parameters <= FULL_MATRIX_PARAMETERS {
            sums.full = Some(FullInformation::new(group_index));
        }
        sums
    }

    /// How many of its own parameters sample `s` carries.
    pub(super) fn own_parameters_of(&self, s: usize) -> usize {
        self.sample_cohort_blocks[s].len() / COHORT_PARAMETERS
    }

    /// Add one position's products. A sample whose own slopes are all zero there — no reads, or a
    /// position that underflowed — would add only zeros, so it is skipped: a saving in time, not
    /// in what is added.
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
            let n = row_scores.len();
            debug_assert_eq!(
                own.len(),
                n * n,
                "a sample's scores and block disagree on its size"
            );
            if row_scores.iter().all(|&score| score == 0.0) {
                continue;
            }
            for (row, &left) in row_scores.iter().enumerate() {
                for (column, &right) in row_scores.iter().enumerate() {
                    own[row * n + column] += left * right;
                }
                for (column, &right) in cohort.iter().enumerate() {
                    with_cohort[row * COHORT_PARAMETERS + column] += left * right;
                }
            }
        }
        if let Some(full) = &mut self.full {
            full.add_position(scores);
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
            debug_assert_eq!(
                into.len(),
                from.len(),
                "two chunks size a sample's block alike"
            );
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
        match (&mut self.full, &other.full) {
            (Some(into), Some(from)) => into.absorb(from),
            (None, None) => {}
            // A total missing a chunk's pairings would give errors from part of the positions.
            _ => unreachable!("two chunks of one pass both keep the whole matrix or neither does"),
        }
    }
}

/// **The whole information matrix**: every parameter's score paired with every other's and summed
/// over the positions, in the layout the standard errors are reported in — the cohort's eight
/// ([`cohort`]), then each sample's own parameters in turn ([`sample`]) — so two samples' parameters
/// are paired too, which the blocks leave out. Kept only for a cohort of at most
/// [`FULL_MATRIX_SAMPLES`] samples and [`FULL_MATRIX_PARAMETERS`] parameters. It is symmetric, so
/// only its upper triangle is held, row after row. The pass that sums it keeps one for each chunk,
/// at most 128, and one total they are added into, which at 20 samples of one library is
/// 129 × 2,346 numbers, 2.4 MB.
#[derive(Clone)]
pub(super) struct FullInformation {
    /// How many parameters the matrix pairs: its side.
    side: usize,
    /// Where each sample's own parameters start in the layout.
    sample_starts: Vec<usize>,
    /// The upper triangle, row after row, each row from its diagonal entry to the end:
    /// `side · (side + 1) / 2` sums.
    upper: Vec<f64>,
    /// One position's scores in the layout's order. Scratch, refilled at every position and not
    /// part of what the sums are.
    row: Vec<f64>,
}

impl std::fmt::Debug for FullInformation {
    /// The sums only, not the scratch row: two sums that print alike are the same bits.
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("FullInformation")
            .field("side", &self.side)
            .field("sample_starts", &self.sample_starts)
            .field("upper", &self.upper)
            .finish_non_exhaustive()
    }
}

impl FullInformation {
    /// Empty sums for samples reading the libraries `group_index` lists.
    fn new(group_index: &[Vec<usize>]) -> Self {
        let mut sample_starts = Vec::with_capacity(group_index.len());
        let mut side = COHORT_PARAMETERS;
        for libraries in group_index {
            sample_starts.push(side);
            side += own_parameters(libraries.len());
        }
        Self {
            side,
            sample_starts,
            upper: vec![0.0; side * (side + 1) / 2],
            row: Vec::with_capacity(side),
        }
    }

    /// How many parameters the matrix pairs.
    pub(super) fn side(&self) -> usize {
        self.side
    }

    /// The rows of sample `s`'s own parameters in the matrix, in its [`sample`] order.
    pub(super) fn rows_of_sample(&self, s: usize) -> std::ops::Range<usize> {
        let end = self.sample_starts.get(s + 1).copied().unwrap_or(self.side);
        self.sample_starts[s]..end
    }

    /// How many samples the matrix pairs the parameters of.
    pub(super) fn samples(&self) -> usize {
        self.sample_starts.len()
    }

    /// Where row `row`'s part of the upper triangle starts: after the rows before it, each one
    /// entry shorter than the last.
    fn row_start(&self, row: usize) -> usize {
        row * self.side - (row * row - row) / 2
    }

    /// The sum at `(row, column)`, in either order.
    pub(super) fn entry(&self, row: usize, column: usize) -> f64 {
        let (row, column) = (row.min(column), row.max(column));
        self.upper[self.row_start(row) + column - row]
    }

    /// Add one position's products. A score of zero adds only zeros, so its row is skipped.
    fn add_position(&mut self, scores: &PositionScores) {
        let Self {
            side, upper, row, ..
        } = self;
        row.clear();
        row.extend_from_slice(&scores.cohort);
        for own in &scores.samples {
            row.extend_from_slice(own);
        }
        debug_assert_eq!(row.len(), *side, "a position's scores fill the layout");
        let mut start = 0;
        for (at, &left) in row.iter().enumerate() {
            let length = *side - at;
            if left != 0.0 {
                for (into, &right) in upper[start..start + length].iter_mut().zip(&row[at..]) {
                    *into += left * right;
                }
            }
            start += length;
        }
    }

    /// Add another chunk's sums.
    fn absorb(&mut self, other: &Self) {
        assert_eq!(
            self.side, other.side,
            "two chunks of one pass pair the same parameters"
        );
        for (into, from) in self.upper.iter_mut().zip(&other.upper) {
            *into += from;
        }
    }

    /// The whole matrix from a dense one in the same layout, `[row * side + column]`, for the tests
    /// that compare the errors it gives with a dense inverse's.
    #[cfg(test)]
    pub(super) fn of_dense(group_index: &[Vec<usize>], dense: &[f64]) -> Self {
        let mut full = Self::new(group_index);
        let side = full.side;
        assert_eq!(
            dense.len(),
            side * side,
            "a dense matrix of the layout's size"
        );
        let mut at = 0;
        for row in 0..side {
            for column in row..side {
                full.upper[at] = dense[row * side + column];
                at += 1;
            }
        }
        full
    }
}

#[cfg(test)]
mod tests {
    use crate::parameter_estimation::joint::census::{DepthCap, SampleGenericSections};
    use crate::parameter_estimation::joint::fit::bench_fixtures::{
        DrawnCohort, DrawnLibrary, as_cohort, draw_cohort_of_libraries,
        draw_cohort_with_duplications,
    };
    use crate::parameter_estimation::joint::fit::standard_errors::{
        StandardError, StandardErrors, inverse_of_positive_definite,
    };
    use crate::parameter_estimation::joint::fit::{
        DuplicatedPositions, EvidenceCursor, FrequencyDensity, JointFit, JointFitConfig,
        LibraryAtPosition, MAX_RECORDED_SPREAD, MIN_POSITIONS_PER_CHUNK, Parameters, PassKeeps,
        Statistics, candidate_read_probability, digamma, expectation_pass, fill_depth_weights,
        fit_jointly, genotype_frequencies, one_position,
    };

    use super::*;

    /// The tolerance every slope is held to, relative to the slope or to 10⁻³ where the slope is
    /// smaller — so a carrier shape's slope of 0.06 is held to 6 × 10⁻⁸ and not to an absolute
    /// 10⁻⁵. The largest disagreement any test here reaches is 6.8 × 10⁻⁷, a homozygote excess in the
    /// two-library fixture (`a_samples_libraries_are_scored_each_under_its_own_rates`); the next are
    /// 1.7 × 10⁻⁷ (the carrier Beta's first shape) and 1.6 × 10⁻⁷ (a ranged library's excess).
    const TOLERANCE: f64 = 1e-6;

    fn relative_disagreement(analytic: f64, numeric: f64) -> f64 {
        (analytic - numeric).abs() / numeric.abs().max(1e-3)
    }

    /// Lend a drawn cohort's sections to `f`, with what the pass needs beside them: the cap, which
    /// read group each sample's sections are, and each library's mean depth.
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
                let coverage = EvidenceCursor::mean_depth_of_each_library(
                    lent,
                    &group_index,
                    groups.len(),
                    &edges,
                    depth_cap,
                );
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
            group_index,
            &config.edges,
            depth_cap,
            coverage,
            config.depth_as_a_range,
            0,
            positions,
        );
        let nodes = model.quadrature.nodes.len();
        let mut scratch = Scratch::new(lent.len(), parameters.clean.len(), nodes);
        let mut statistics = Statistics::new(parameters.clean.len(), lent.len(), nodes);
        let mut here = PositionScores::new(group_index, parameters.clean.len());
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
        let mut total = PositionScores::new(group_index, parameters.clean.len());
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
    /// That form is `Σ posterior · ln f` over the rule's nodes, less the segregating posterior times
    /// `ψ(a) − ψ(a + b)`, where `ψ` is the digamma function, the slope of `ln Γ` — the form the
    /// maximisation's shapes' update solved until plan step A7, and the sum a test build's pass still
    /// keeps (`Statistics::sum_ln_f`), so the statistics give it directly. It lands more than 3%
    /// from the finite difference of the log-likelihood at twelve nodes and more than 2% at
    /// sixteen, while the slope this module takes lands within a part in a million at both. If the
    /// rule ever integrated `ln f` well, the first assertions would fail and the module doc's
    /// explanation would be out of date.
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
            |lent, depth_cap, group_index, coverage| {
                let config = JointFitConfig::default();
                let mut cursor = EvidenceCursor::over(
                    lent,
                    group_index,
                    &config.edges,
                    depth_cap,
                    coverage,
                    true,
                    0,
                    EvidenceCursor::position_count(lent),
                );
                let mut scratch = Scratch::new(lent.len(), coverage.len(), 1);
                let mut ranged = 0;
                while cursor.next_position(&mut scratch.evidence) {
                    ranged += scratch
                        .evidence
                        .libraries
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
    #[derive(Clone)]
    struct StoredPosition {
        libraries: Vec<LibraryAtPosition>,
        depth_weights: Vec<f64>,
        observed_alternatives: Vec<usize>,
    }

    /// The depth-8 fixture's positions, with one reference read of every fourth library-position
    /// turned into a read held out of the model — an insertion, a deletion or a spanning deletion,
    /// which the drawn cohort never produces. The fixture has one library a sample, library `s`
    /// being sample `s`'s.
    fn positions_with_held_out_reads(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        group_index: &[Vec<usize>],
        coverage: &[f64],
    ) -> Vec<StoredPosition> {
        let config = JointFitConfig::default();
        let mut cursor = EvidenceCursor::over(
            lent,
            group_index,
            &config.edges,
            depth_cap,
            coverage,
            true,
            0,
            EvidenceCursor::position_count(lent),
        );
        let mut scratch = Scratch::new(lent.len(), coverage.len(), 1);
        let (mut stored, mut index) = (Vec::new(), 0_usize);
        while cursor.next_position(&mut scratch.evidence) {
            let mut libraries = scratch.evidence.libraries.clone();
            for (l, library) in libraries.iter_mut().enumerate() {
                if (index + l) % 4 == 0 && library.fewest_reference >= 1.0 && library.spread == 1 {
                    library.on[4] += 1.0;
                    library.fewest_reference -= 1.0;
                }
            }
            stored.push(StoredPosition {
                libraries,
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
        let mut scratch = Scratch::new(samples, parameters.clean.len(), nodes);
        let mut statistics = Statistics::new(parameters.clean.len(), samples, nodes);
        let mut here = PositionScores::new(group_index, parameters.clean.len());
        let mut total = 0.0;
        for position in stored {
            scratch.evidence.libraries.clone_from(&position.libraries);
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
            let stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let held_out: usize = stored
                .iter()
                .map(|position| {
                    position
                        .libraries
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
            let mut scores = PositionScores::new(group_index, parameters.clean.len());
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

    /// The stored fixture with sample 0's reads split into two libraries: the first keeps read
    /// group 0, the second is read group `second`, and the second takes half of every count,
    /// rounded down — reference reads included, reads held out of the model left with the first.
    /// Every position of the depth-8 fixture has a depth that is one value, which the split keeps.
    fn split_sample_zero(stored: &[StoredPosition], second: usize) -> Vec<StoredPosition> {
        stored
            .iter()
            .map(|position| {
                let mut split = position.clone();
                let whole = position.libraries[0];
                assert_eq!(whole.spread, 1, "the fixture's depths are single values");
                let (mut kept, mut moved) = (whole, LibraryAtPosition::default());
                for code in 0..4 {
                    moved.on[code] = (whole.on[code] / 2.0).floor();
                    kept.on[code] -= moved.on[code];
                }
                moved.fewest_reference = (whole.fewest_reference / 2.0).floor();
                kept.fewest_reference -= moved.fewest_reference;
                moved.spread = 1;
                for library in [&mut kept, &mut moved] {
                    library.depth = library.non_reference() + library.fewest_reference;
                }
                split.libraries[0] = kept;
                if split.libraries.len() <= second {
                    split
                        .libraries
                        .resize(second + 1, LibraryAtPosition::default());
                    split
                        .depth_weights
                        .resize((second + 1) * MAX_RECORDED_SPREAD, 0.0);
                }
                split.libraries[second] = moved;
                split.depth_weights[second * MAX_RECORDED_SPREAD] = 1.0;
                split
            })
            .collect()
    }

    /// **Each of a sample's libraries is scored under its own rates** (plan step A6): sample 0's
    /// reads are split into two libraries, the second a read group of its own.
    ///
    /// - **At equal rates the split changes nothing**: a read's probability depends on the rate
    ///   and the sample's genotype, not on which library it came from, so the log-likelihood is
    ///   the pooled one to rounding.
    /// - **At different rates each library's rate moves the log-likelihood**: a nudge to the
    ///   second library's rates moves it, where the likelihood before this step scored every
    ///   library of a sample under its first library's rates and left it bit-for-bit unchanged.
    /// - **Every library's two rates' slopes are the derivative of the log-likelihood**, each read
    ///   from its own slot of its sample's row — both of sample 0's libraries at rates unlike each
    ///   other's and every other sample's — and so is every sample's homozygote excess's.
    #[test]
    fn a_samples_libraries_are_scored_each_under_its_own_rates() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let second = parameters.clean.len();
            let split = split_sample_zero(&stored, second);
            let split_index = vec![vec![0, second], vec![1], vec![2], vec![3]];

            let pooled = over_stored(&stored, group_index, &[], &parameters, None);
            let mut alike = parameters.clone();
            alike.clean.push(parameters.clean[0]);
            alike.noisy.push(parameters.noisy[0]);
            let apart = over_stored(&split, &split_index, &[], &alike, None);
            eprintln!("pooled {pooled:.10e}, two libraries at one rate {apart:.10e}");
            assert!(
                (apart - pooled).abs() < 1e-10 * pooled.abs(),
                "two libraries at one rate give {apart} against the pooled {pooled}"
            );

            let mut own = parameters.clone();
            own.clean.push(0.0079);
            own.noisy.push(0.113);
            // Sample 2's mismapped rate moved off the value where its slope is near zero (−0.038 at
            // 0.062), which a finite difference cannot resolve from the likelihood's rounding.
            own.noisy[2] = 0.03;
            let values: Vec<f64> = own.coordinates().iter().map(|c| c.value).collect();
            let at = |index: usize, by: f64| {
                let mut changed = values.clone();
                changed[index] += by;
                over_stored(
                    &split,
                    &split_index,
                    &[],
                    &own.with_coordinates(&changed),
                    None,
                )
            };
            for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                let index = COHORT_PARAMETERS + 2 * second + which;
                let step = 1e-5 * values[index];
                let slope = (at(index, step) - at(index, -step)) / (2.0 * step);
                eprintln!("the second library's rate {which}: slope {slope:.6e}");
                assert!(
                    slope.abs() > 10.0,
                    "the second library's rate {which} barely moves the log-likelihood: {slope}"
                );
            }

            let mut scores = PositionScores::new(&split_index, own.clean.len());
            over_stored(&split, &split_index, &[], &own, Some(&mut scores));
            let mut worst = 0.0_f64;
            let mut check = |name: String, analytic: f64, index: usize| {
                let step = 1e-5 * values[index];
                let numeric = (at(index, step) - at(index, -step)) / (2.0 * step);
                let disagreement = relative_disagreement(analytic, numeric);
                eprintln!(
                    "{name}: summed slope {analytic:.8e}, central difference {numeric:.8e}, \
                     relative disagreement {disagreement:.2e}"
                );
                worst = worst.max(disagreement);
            };
            for (s, own_libraries) in split_index.iter().enumerate() {
                for (section, &library) in own_libraries.iter().enumerate() {
                    for class in 0..2 {
                        check(
                            format!("sample {s}'s library {library}, rate {class}"),
                            scores.samples[s][sample::rate(section, class)],
                            COHORT_PARAMETERS + 2 * library + class,
                        );
                    }
                }
                check(
                    format!("sample {s}'s homozygote excess"),
                    scores.samples[s][sample::HOMOZYGOTE_EXCESS],
                    sample_coordinate(&own, s, sample::HOMOZYGOTE_EXCESS),
                );
            }
            assert!(
                worst < TOLERANCE,
                "largest relative disagreement {worst:.2e}"
            );
        });
    }

    /// The stored fixture with sample 0 given a second library, read group `second`, that holds no
    /// read at any position.
    fn with_an_empty_second_library(
        stored: &[StoredPosition],
        second: usize,
    ) -> Vec<StoredPosition> {
        stored
            .iter()
            .map(|position| {
                let mut with = position.clone();
                with.libraries
                    .resize(second + 1, LibraryAtPosition::default());
                with.libraries[second].spread = 1;
                with.depth_weights
                    .resize((second + 1) * MAX_RECORDED_SPREAD, 0.0);
                with.depth_weights[second * MAX_RECORDED_SPREAD] = 1.0;
                with
            })
            .collect()
    }

    /// **Each library's reads are scored under its own rates, and under no other library's** —
    /// two exact consequences, where a nudge's size would let a wrong library's rates through.
    ///
    /// - A second library that holds no read leaves sample 0's likelihood exactly what its first
    ///   library alone gives, bit for bit, whatever the second's rates: scoring the first library's
    ///   reads under the second's rates (or under the sample's last library's) moves it.
    /// - Exchanging the two libraries' reads together with their rates leaves the likelihood exactly
    ///   where it was: scoring both libraries under one library's rates moves it.
    #[test]
    fn a_librarys_reads_are_scored_under_its_own_rates_and_no_other() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let second = parameters.clean.len();
            let index = vec![vec![0, second], vec![1], vec![2], vec![3]];
            let alone = over_stored(&stored, group_index, &[], &parameters, None);
            let empty = with_an_empty_second_library(&stored, second);
            for (clean, noisy) in [(0.0079, 0.113), (0.0003, 0.0021)] {
                let mut own = parameters.clone();
                own.clean.push(clean);
                own.noisy.push(noisy);
                let with_empty = over_stored(&empty, &index, &[], &own, None);
                assert_eq!(
                    with_empty.to_bits(),
                    alone.to_bits(),
                    "an empty second library at rates {clean}, {noisy} moved sample 0's \
                     likelihood from {alone} to {with_empty}"
                );
            }

            let split = split_sample_zero(&stored, second);
            let exchanged: Vec<StoredPosition> = split
                .iter()
                .map(|position| {
                    let mut swapped = position.clone();
                    swapped.libraries.swap(0, second);
                    swapped
                })
                .collect();
            let mut own = parameters.clone();
            own.clean.push(0.0079);
            own.noisy.push(0.113);
            let mut exchanged_rates = own.clone();
            exchanged_rates.clean.swap(0, second);
            exchanged_rates.noisy.swap(0, second);
            let as_split = over_stored(&split, &index, &[], &own, None);
            let as_exchanged = over_stored(&exchanged, &index, &[], &exchanged_rates, None);
            eprintln!("split {as_split:.10e}, exchanged {as_exchanged:.10e}");
            assert_eq!(
                as_split.to_bits(),
                as_exchanged.to_bits(),
                "exchanging two libraries' reads with their rates moved the likelihood"
            );
        });
    }

    /// **Every branch credits each library's tallies with that library's own reads, and no other
    /// library's.** Over one position, a library's read tallies summed over both classes, all four
    /// branches and every genotype hold exactly its own reads — the reads on a base that is not the
    /// reference once (candidate plus neither), and its reference reads once — because the posteriors
    /// they are weighted by sum to one. Sample 0's reads are split over two libraries at different
    /// rates, with the duplicated class fitted, so all four branches carry weight. Crediting any one
    /// branch with the sample's pooled reads, or with the first library's, breaks it.
    #[test]
    fn each_librarys_tallies_hold_its_own_reads_in_every_branch() {
        let (cohort, mut parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        assert!(
            parameters.duplicated.is_some(),
            "all four branches must carry weight"
        );
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let second = parameters.clean.len();
            let split = split_sample_zero(&stored, second);
            let index = vec![vec![0, second], vec![1], vec![2], vec![3]];
            parameters.clean.push(0.0079);
            parameters.noisy.push(0.113);
            let config = JointFitConfig::default();
            let model = PassModel::new(&parameters, &config, &index);
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(index.len(), parameters.clean.len(), nodes);
            let (mut worst, mut mass) = (0.0_f64, [0.0_f64; 4]);
            for position in &split {
                let mut statistics = Statistics::new(parameters.clean.len(), index.len(), nodes);
                scratch.evidence.libraries.clone_from(&position.libraries);
                scratch
                    .evidence
                    .depth_weights
                    .clone_from(&position.depth_weights);
                scratch
                    .evidence
                    .observed_alternatives
                    .clone_from(&position.observed_alternatives);
                one_position(&mut scratch, &model, &[], &mut statistics);
                for (slot, branch) in mass.iter_mut().zip([
                    statistics.invariant,
                    statistics.fixed_alt,
                    statistics.segregating,
                    statistics.duplicated,
                ]) {
                    *slot += branch;
                }
                for library in [0, second] {
                    let reads = &position.libraries[library];
                    let own = [reads.non_reference() - reads.on[4], reads.fewest_reference];
                    let tallies = statistics.reads[library].iter().flatten();
                    let credited = tallies.fold([0.0, 0.0], |[off, on], tally| {
                        [off + tally.candidate + tally.neither, on + tally.reference]
                    });
                    for (credited, own) in credited.into_iter().zip(own) {
                        worst = worst.max((credited - own).abs() / own.max(1.0));
                    }
                }
            }
            eprintln!(
                "branch mass over {} positions: invariant {:.2}, fixed {:.2}, segregating {:.2}, \
                 duplicated {:.2}; worst relative gap {worst:.2e}",
                split.len(),
                mass[0],
                mass[1],
                mass[2],
                mass[3]
            );
            assert!(
                mass.iter().all(|&m| m > 1.0),
                "every branch must carry weight: {mass:?}"
            );
            assert!(
                worst < 1e-9,
                "a library's tallies are off its own reads by {worst:.2e}"
            );
        });
    }

    /// **Each library's stored depth range is weighted around that library's own mean depth.** Two
    /// samples, each read from a library at about 40 reads a position and one at about 132, recorded
    /// under a cap of 140 so the deep library's depths above 124 are ranges: each library's mean
    /// depth comes back near its own draw, and every ranged library-position's weights are
    /// [`fill_depth_weights`] at that library's own mean — which the sample's summed mean, or the
    /// other library's, is not.
    #[test]
    fn each_library_is_weighted_around_its_own_depth() {
        let libraries = |shallow: f64, deep: f64| {
            vec![
                DrawnLibrary {
                    clean: 0.003,
                    noisy: 0.06,
                    mean_depth: shallow,
                },
                DrawnLibrary {
                    clean: 0.006,
                    noisy: 0.09,
                    mean_depth: deep,
                },
            ]
        };
        let cohort = draw_cohort_of_libraries(
            &[libraries(40.0, 132.0), libraries(40.0, 132.0)],
            0.03,
            400,
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.0,
            0x0F17_5C03_A600_0003,
        );
        let ranged = with_sections(
            &cohort,
            Some(DepthCap::new(140)),
            |lent, depth_cap, group_index, coverage| {
                for (library, (&mean, drawn)) in
                    coverage.iter().zip([40.0, 132.0, 40.0, 132.0]).enumerate()
                {
                    eprintln!("library {library}: mean depth {mean:.2}, drawn at {drawn}");
                    assert!(
                        (mean / drawn - 1.0).abs() < 0.05,
                        "library {library}'s mean depth {mean} against the {drawn} it was drawn at"
                    );
                }
                let config = JointFitConfig::default();
                let mut cursor = EvidenceCursor::over(
                    lent,
                    group_index,
                    &config.edges,
                    depth_cap,
                    coverage,
                    true,
                    0,
                    EvidenceCursor::position_count(lent),
                );
                let mut scratch = Scratch::new(lent.len(), coverage.len(), 1);
                let mut ranged = 0;
                let mut expected = [0.0; MAX_RECORDED_SPREAD];
                while cursor.next_position(&mut scratch.evidence) {
                    for (library, reads) in scratch.evidence.libraries.iter().enumerate() {
                        if reads.spread <= 1 {
                            continue;
                        }
                        ranged += 1;
                        fill_depth_weights(
                            &mut expected[..reads.spread],
                            coverage[library],
                            reads.fewest_reference,
                        );
                        assert_eq!(
                            scratch.evidence.weights_of(library)[..reads.spread],
                            expected[..reads.spread],
                            "library {library}'s range is weighted around another depth"
                        );
                    }
                }
                ranged
            },
        );
        eprintln!("{ranged} ranged library-positions");
        assert!(ranged > 200, "only {ranged} library-positions were ranges");
    }

    /// **Each library's rate slopes hold where its depths are ranges**: two samples, each read from
    /// a library at about 40 reads a position and one at about 132, recorded under a cap of 140 so
    /// the deep library's depths above 124 are ranges; every library's two slopes and each sample's
    /// excess against a central difference, at parameters away from the maximum and unlike from
    /// library to library.
    #[test]
    fn each_librarys_slopes_hold_where_its_depths_are_ranges() {
        let library = |clean: f64, noisy: f64, mean_depth: f64| DrawnLibrary {
            clean,
            noisy,
            mean_depth,
        };
        let cohort = draw_cohort_of_libraries(
            &[
                vec![library(0.003, 0.06, 40.0), library(0.006, 0.09, 132.0)],
                vec![library(0.004, 0.05, 40.0), library(0.002, 0.08, 132.0)],
            ],
            0.03,
            400,
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.0,
            0x0F17_5C03_A600_0004,
        );
        let parameters = Parameters {
            clean: vec![0.0021, 0.0071, 0.0043, 0.0013],
            noisy: vec![0.041, 0.113, 0.083, 0.062],
            noisy_share: 0.047,
            density: FrequencyDensity {
                p_invariant: 0.83,
                p_fixed_alt: 0.023,
                a: 0.9,
                b: 1.6,
            },
            hom_excess: vec![0.12, 0.47],
            duplicated: None,
        };
        let config = JointFitConfig {
            quadrature_nodes: 12,
            duplicated_positions: false,
            ..JointFitConfig::default()
        };
        with_sections(
            &cohort,
            Some(DepthCap::new(140)),
            |lent, depth_cap, group_index, coverage| {
                let (_, scores) =
                    summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
                let cohort_count = cohort_coordinates(&parameters);
                let mut worst = 0.0_f64;
                for (s, own) in group_index.iter().enumerate() {
                    let mut slots: Vec<(usize, usize)> = own
                        .iter()
                        .enumerate()
                        .flat_map(|(section, &library)| {
                            (0..2).map(move |class| {
                                (
                                    sample::rate(section, class),
                                    cohort_count + 2 * library + class,
                                )
                            })
                        })
                        .collect();
                    slots.push((
                        sample::HOMOZYGOTE_EXCESS,
                        cohort_count + 2 * parameters.clean.len() + s,
                    ));
                    for (slot, index) in slots {
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
                        let disagreement = relative_disagreement(scores.samples[s][slot], numeric);
                        eprintln!(
                            "sample {s}, slot {slot}: summed slope {:.8e}, central difference \
                             {numeric:.8e}, relative disagreement {disagreement:.2e}",
                            scores.samples[s][slot]
                        );
                        worst = worst.max(disagreement);
                    }
                }
                assert!(
                    worst < TOLERANCE,
                    "largest relative disagreement {worst:.2e}"
                );
            },
        );
    }

    /// **A sample of three libraries has every slope right**, through a whole pass: one sample of
    /// one library beside one of three, every library at its own rates and the duplicated class
    /// fitted, so all four branches credit a third library's slot. Each library's two rate slopes
    /// and each excess against a central difference of the log-likelihood.
    #[test]
    fn a_sample_of_three_libraries_has_every_slope_right() {
        let library = |clean: f64, noisy: f64| DrawnLibrary {
            clean,
            noisy,
            mean_depth: 4.0,
        };
        let cohort = draw_cohort_of_libraries(
            &[
                vec![library(0.003, 0.06)],
                vec![
                    library(0.002, 0.05),
                    library(0.006, 0.09),
                    library(0.004, 0.12),
                ],
            ],
            0.03,
            600,
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_A600_0006,
        );
        let parameters = Parameters {
            clean: vec![0.0021, 0.0043, 0.0071, 0.0013],
            noisy: vec![0.041, 0.083, 0.113, 0.062],
            noisy_share: 0.047,
            density: FrequencyDensity {
                p_invariant: 0.83,
                p_fixed_alt: 0.023,
                a: 0.9,
                b: 1.6,
            },
            hom_excess: vec![0.12, 0.47],
            duplicated: Some(DuplicatedPositions {
                share: 0.006,
                carrier_a: 1.7,
                carrier_b: 6.3,
            }),
        };
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            assert_eq!(group_index[1].len(), 3, "sample 1 reads three libraries");
            let (_, scores) =
                summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
            let cohort_count = cohort_coordinates(&parameters);
            let mut worst = 0.0_f64;
            for (s, own) in group_index.iter().enumerate() {
                assert_eq!(scores.samples[s].len(), own_parameters(own.len()));
                let mut slots: Vec<(usize, usize)> = own
                    .iter()
                    .enumerate()
                    .flat_map(|(section, &library)| {
                        (0..2).map(move |class| {
                            (
                                sample::rate(section, class),
                                cohort_count + 2 * library + class,
                            )
                        })
                    })
                    .collect();
                slots.push((
                    sample::HOMOZYGOTE_EXCESS,
                    cohort_count + 2 * parameters.clean.len() + s,
                ));
                for (slot, index) in slots {
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
                    let disagreement = relative_disagreement(scores.samples[s][slot], numeric);
                    eprintln!(
                        "sample {s}, slot {slot}: summed slope {:.8e}, central difference \
                         {numeric:.8e}, relative disagreement {disagreement:.2e}",
                        scores.samples[s][slot]
                    );
                    worst = worst.max(disagreement);
                }
            }
            assert!(
                worst < TOLERANCE,
                "largest relative disagreement {worst:.2e}"
            );
        });
    }

    /// **A sample holding no library carries no parameter information, and changes nothing else**
    /// — a sample whose census holds repeat tracts only, which the census keeps. Sample 2 is given
    /// no library: the log-likelihood is bit-for-bit the one where its library holds no read, its
    /// row is three zeros at every position, its three errors say no information, and the log counts
    /// the three read groups the cohort holds, not four.
    #[test]
    fn a_sample_holding_no_library_carries_no_information() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let mut stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let silent = 2;
            for position in &mut stored {
                position.libraries[silent] = LibraryAtPosition {
                    spread: 1,
                    ..LibraryAtPosition::default()
                };
                position.depth_weights[silent * MAX_RECORDED_SPREAD] = 1.0;
            }
            let without = vec![vec![0], vec![1], vec![], vec![3]];
            let empty_library = over_stored(&stored, group_index, &[], &parameters, None);
            let mut scores = PositionScores::new(&without, parameters.clean.len());
            let no_library = over_stored(&stored, &without, &[], &parameters, Some(&mut scores));
            assert_eq!(no_library.to_bits(), empty_library.to_bits());
            assert_eq!(scores.samples[silent], [0.0; 3]);

            let config = JointFitConfig {
                quadrature_nodes: 12,
                ..JointFitConfig::default()
            };
            let model = PassModel::new(&parameters, &config, &without);
            let tables = ScoringTables::of(&model);
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(without.len(), parameters.clean.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), without.len(), nodes);
            let mut here = PositionScores::new(&without, parameters.clean.len());
            let mut information = InformationSums::new(&without);
            for position in &stored {
                scratch.evidence.libraries.clone_from(&position.libraries);
                scratch
                    .evidence
                    .depth_weights
                    .clone_from(&position.depth_weights);
                scratch
                    .evidence
                    .observed_alternatives
                    .clone_from(&position.observed_alternatives);
                one_position(&mut scratch, &model, &[], &mut statistics);
                score_position(&scratch, &model, &tables, &[], &mut here);
                information.add_position(&here);
            }
            let errors = StandardErrors::of(&information);
            assert_eq!(errors.samples[silent], [StandardError::NoInformation; 3]);
            let line = errors.described(&parameters, &without);
            assert!(
                line.contains("error rates at ordinary positions (3 read groups)"),
                "{line}"
            );
        });
    }

    /// **The slopes the shapes' update steps on are the likelihood's slopes**: a pass's summed
    /// `density_shape_slopes` and `carrier_shape_slopes` equal the scorer's summed slopes in the
    /// same four shapes — which the finite-difference tests above hold to the log-likelihood — with
    /// the duplicated class fitted and without it, over a census the pass cuts into three chunks and
    /// joins. The pass skips nodes whose posterior is below 10⁻¹² where the scorer keeps them, so the
    /// two differ by rounding and those nodes (measured 1.8 × 10⁻⁹ relative; held to 10⁻⁷).
    /// Dropping the along-the-frequency term, flipping its sign, reading another shape's or another
    /// node's slope, or losing a chunk's slopes when chunks are joined fails it.
    #[test]
    fn the_passs_shape_slopes_are_the_likelihoods() {
        for duplicated in [true, false] {
            let (cohort, mut parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
            if !duplicated {
                parameters.duplicated = None;
            }
            let config = JointFitConfig {
                duplicated_positions: duplicated,
                ..JointFitConfig::default()
            };
            with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                assert!(EvidenceCursor::position_count(lent) > 2 * MIN_POSITIONS_PER_CHUNK);
                let statistics = expectation_pass(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                    PassKeeps::SUMS_ONLY,
                );
                let (_, scores) =
                    summed_scores(lent, depth_cap, &config, group_index, coverage, &parameters);
                let pairs = [
                    (
                        "density a",
                        statistics.density_shape_slopes[0],
                        cohort::DENSITY_A,
                    ),
                    (
                        "density b",
                        statistics.density_shape_slopes[1],
                        cohort::DENSITY_B,
                    ),
                    (
                        "carrier a",
                        statistics.carrier_shape_slopes[0],
                        cohort::CARRIER_A,
                    ),
                    (
                        "carrier b",
                        statistics.carrier_shape_slopes[1],
                        cohort::CARRIER_B,
                    ),
                ];
                for (name, pass, slot) in pairs {
                    let scorer = scores.cohort[slot];
                    let off = relative_disagreement(pass, scorer);
                    eprintln!(
                        "duplicated class {duplicated}, {name}: pass {pass:.10e}, scorer \
                         {scorer:.10e}, {off:.2e} apart"
                    );
                    assert!(
                        off < 1e-7,
                        "duplicated class {duplicated}, {name}: the pass's slope {pass} against \
                         the scorer's {scorer}"
                    );
                }
                if !duplicated {
                    assert_eq!(statistics.carrier_shape_slopes, [0.0; 2]);
                }
            });
        }
    }

    /// **The carrier shapes' slopes follow each sample's coverage odds**: with the odds set (a
    /// quarter, three, one and eight), a pass's carrier slopes equal the scorer's, position by
    /// position summed — the one input the drawn cohort never sets, and dropping the odds from the
    /// pass's slope fails it.
    #[test]
    fn the_passs_carrier_slopes_follow_the_coverage_odds() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            let odds = [0.25, 3.0, 1.0, 8.0];
            let config = JointFitConfig::default();
            let samples = group_index.len();
            let model = PassModel::new(&parameters, &config, group_index);
            let tables = ScoringTables::of(&model);
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(samples, parameters.clean.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), samples, nodes);
            let mut here = PositionScores::new(group_index, parameters.clean.len());
            let mut summed = [0.0_f64; 2];
            for position in &stored {
                scratch.evidence.libraries.clone_from(&position.libraries);
                scratch
                    .evidence
                    .depth_weights
                    .clone_from(&position.depth_weights);
                scratch
                    .evidence
                    .observed_alternatives
                    .clone_from(&position.observed_alternatives);
                one_position(&mut scratch, &model, &odds, &mut statistics);
                score_position(&scratch, &model, &tables, &odds, &mut here);
                summed[0] += here.cohort[cohort::CARRIER_A];
                summed[1] += here.cohort[cohort::CARRIER_B];
            }
            for (shape, (&pass, &scorer)) in statistics
                .carrier_shape_slopes
                .iter()
                .zip(&summed)
                .enumerate()
            {
                let off = relative_disagreement(pass, scorer);
                eprintln!(
                    "carrier shape {shape}: pass {pass:.10e}, scorer {scorer:.10e}, {off:.2e}"
                );
                assert!(
                    off < 1e-7,
                    "carrier shape {shape}: pass {pass}, scorer {scorer}"
                );
            }
        });
    }

    /// The information blocks formed from every position's scores in position order, each product
    /// written out — the independent account a pass's [`InformationSums`] is checked against.
    fn information_by_hand(every: &[PositionScores]) -> InformationSums {
        // Each sample's library count, read back off its row's length.
        let shape: Vec<Vec<usize>> = every.first().map_or_else(Vec::new, |scores| {
            scores
                .samples
                .iter()
                .map(|row| vec![0; (row.len() - 1) / 2])
                .collect()
        });
        let mut sums = InformationSums::new(&shape);
        for here in every {
            for row in 0..COHORT_PARAMETERS {
                for column in 0..COHORT_PARAMETERS {
                    sums.cohort[row * COHORT_PARAMETERS + column] +=
                        here.cohort[row] * here.cohort[column];
                }
            }
            for (s, own) in here.samples.iter().enumerate() {
                let n = own.len();
                for row in 0..n {
                    for column in 0..n {
                        sums.sample_blocks[s][row * n + column] += own[row] * own[column];
                    }
                    for column in 0..COHORT_PARAMETERS {
                        sums.sample_cohort_blocks[s][row * COHORT_PARAMETERS + column] +=
                            own[row] * here.cohort[column];
                    }
                }
            }
        }
        sums
    }

    /// Every position's scores laid out flat — the cohort's eight, then each sample's own in turn
    /// — multiplied pairwise and summed: the whole matrix, dense and row-major, and its side.
    fn whole_matrix_by_hand(every: &[PositionScores]) -> (Vec<f64>, usize) {
        let side = every.first().map_or(COHORT_PARAMETERS, |scores| {
            COHORT_PARAMETERS + scores.samples.iter().map(Vec::len).sum::<usize>()
        });
        let mut dense = vec![0.0; side * side];
        for here in every {
            let flat: Vec<f64> = here
                .cohort
                .iter()
                .chain(here.samples.iter().flatten())
                .copied()
                .collect();
            for (row, &left) in flat.iter().enumerate() {
                for (column, &right) in flat.iter().enumerate() {
                    dense[row * side + column] += left * right;
                }
            }
        }
        (dense, side)
    }

    /// **The pass keeps the whole matrix up to [`FULL_MATRIX_SAMPLES`] samples and
    /// [`FULL_MATRIX_PARAMETERS`] parameters, and not above either**: sums made for 20 samples hold
    /// it, sized for their own parameters — a sample of two libraries counting five — and sums made
    /// for 21 do not; 20 samples of four libraries each (188 parameters) hold it, and one library
    /// more (190) does not.
    #[test]
    fn the_whole_matrix_is_kept_up_to_twenty_samples() {
        let mut twenty: Vec<Vec<usize>> = (0..FULL_MATRIX_SAMPLES).map(|s| vec![s]).collect();
        twenty[3].push(FULL_MATRIX_SAMPLES);
        let sums = InformationSums::for_a_cohort_of(&twenty);
        let full = sums
            .full
            .as_ref()
            .expect("twenty samples keep the whole matrix");
        assert_eq!(full.side(), COHORT_PARAMETERS + 3 * 19 + 5);
        assert_eq!(
            full.rows_of_sample(3),
            COHORT_PARAMETERS + 9..COHORT_PARAMETERS + 14
        );
        let twenty_one: Vec<Vec<usize>> = (0..=FULL_MATRIX_SAMPLES).map(|s| vec![s]).collect();
        assert!(
            InformationSums::for_a_cohort_of(&twenty_one).full.is_none(),
            "twenty-one samples keep the blocks only"
        );
        assert!(InformationSums::new(&twenty).full.is_none());
        let mut four_each: Vec<Vec<usize>> = (0..FULL_MATRIX_SAMPLES)
            .map(|s| (4 * s..4 * s + 4).collect())
            .collect();
        let full = InformationSums::for_a_cohort_of(&four_each).full;
        assert_eq!(
            full.map(|full| full.side()),
            Some(FULL_MATRIX_PARAMETERS),
            "twenty samples of four libraries keep the whole matrix"
        );
        four_each[0].push(4 * FULL_MATRIX_SAMPLES);
        assert!(
            InformationSums::for_a_cohort_of(&four_each).full.is_none(),
            "190 parameters keep the blocks only"
        );
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
            let n = hand.own_parameters_of(s);
            let own_diagonal = |i: usize| hand.sample_blocks[s][i * n + i];
            for row in 0..n {
                for column in 0..n {
                    let at = row * n + column;
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
    /// four-sample fixture, on one sample alone** — the smallest cohort the fit takes — **and on
    /// samples of one, two and three libraries**, whose blocks are 3, 5 and 7 wide.
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
        // And a cohort whose samples read one, two and three libraries, so blocks of 3, 5 and 7.
        let library = |clean: f64, noisy: f64| DrawnLibrary {
            clean,
            noisy,
            mean_depth: 4.0,
        };
        let mixed = draw_cohort_of_libraries(
            &[
                vec![library(0.003, 0.06)],
                vec![library(0.002, 0.05), library(0.006, 0.09)],
                vec![
                    library(0.004, 0.07),
                    library(0.001, 0.04),
                    library(0.008, 0.11),
                ],
            ],
            0.03,
            600,
            FrequencyDensity {
                p_invariant: 0.88,
                p_fixed_alt: 0.01,
                a: 0.6,
                b: 2.2,
            },
            0.3,
            0.01,
            0x0F17_5C03_A600_0005,
        );
        let mut parameters_of_mixed = parameters_of_four.clone();
        parameters_of_mixed.clean = vec![0.0021, 0.0043, 0.0035, 0.0057, 0.0012, 0.0074];
        parameters_of_mixed.noisy = vec![0.041, 0.083, 0.062, 0.097, 0.053, 0.121];
        parameters_of_mixed.hom_excess.truncate(3);
        for (cohort, parameters) in [
            (&four, &parameters_of_four),
            (&one, &parameters_of_one),
            (&mixed, &parameters_of_mixed),
        ] {
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
                // The whole matrix, every pairing including two samples', against the same
                // products summed by hand.
                let full = summed
                    .full
                    .as_ref()
                    .expect("a cohort this small keeps the whole matrix");
                let (by_hand_full, side) = whole_matrix_by_hand(&every);
                assert_eq!(full.side(), side, "the whole matrix's size");
                let mut worst_full = 0.0_f64;
                for row in 0..side {
                    for column in 0..side {
                        let hand = by_hand_full[row * side + column];
                        let scale = (by_hand_full[row * side + row]
                            * by_hand_full[column * side + column])
                            .sqrt();
                        if scale == 0.0 {
                            assert_eq!(full.entry(row, column), 0.0, "({row}, {column})");
                            continue;
                        }
                        worst_full = worst_full.max((full.entry(row, column) - hand).abs() / scale);
                    }
                }
                eprintln!("the whole matrix: worst disagreement {worst_full:.2e}");
                assert!(worst_full < 1e-12, "the whole matrix: {worst_full:.2e}");
                // Every sample has reads here, so every diagonal entry is positive.
                for (s, block) in summed.sample_blocks.iter().enumerate() {
                    let n = summed.own_parameters_of(s);
                    for which in 0..n {
                        assert!(block[which * n + which] > 0.0, "sample {s}, slot {which}");
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

    /// **A sample with no reads has no row in the whole matrix either**: through the sums a pass
    /// makes for a cohort this small ([`InformationSums::for_a_cohort_of`]), sample 2's rows and
    /// columns — its pairings with every other parameter — are exactly zero, every other entry equals
    /// the products summed by hand to 10⁻¹² of its bound, the errors say it has no information, and
    /// the other samples keep their rates' errors.
    #[test]
    fn a_sample_without_reads_has_no_row_in_the_whole_matrix() {
        let (cohort, parameters) = a_cohort_and_parameters_off_the_maximum(8.0);
        let config = JointFitConfig {
            quadrature_nodes: 12,
            ..JointFitConfig::default()
        };
        let silent = 2;
        with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
            let mut stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            for position in &mut stored {
                position.libraries[silent] = LibraryAtPosition {
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
            let mut scratch = Scratch::new(group_index.len(), parameters.clean.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), group_index.len(), nodes);
            let mut here = PositionScores::new(group_index, parameters.clean.len());
            let mut information = InformationSums::for_a_cohort_of(group_index);
            let mut every = Vec::new();
            for position in &stored {
                scratch.evidence.libraries.clone_from(&position.libraries);
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
                every.push(here.clone());
            }
            let full = information.full.as_ref().expect("kept");
            let (dense, side) = whole_matrix_by_hand(&every);
            let silent_rows = full.rows_of_sample(silent);
            let mut worst = 0.0_f64;
            for row in 0..side {
                for column in 0..side {
                    if silent_rows.contains(&row) || silent_rows.contains(&column) {
                        assert_eq!(full.entry(row, column), 0.0, "({row}, {column})");
                    }
                    let scale = (dense[row * side + row] * dense[column * side + column]).sqrt();
                    if scale > 0.0 {
                        worst = worst.max(
                            (full.entry(row, column) - dense[row * side + column]).abs() / scale,
                        );
                    }
                }
            }
            assert!(worst < 1e-12, "{worst:.2e}");
            let errors = StandardErrors::of(&information);
            assert_eq!(errors.samples[silent], [StandardError::NoInformation; 3]);
            for s in [0, 1, 3] {
                for which in [sample::CLEAN_ERROR_RATE, sample::NOISY_ERROR_RATE] {
                    assert!(
                        errors.samples[s][which].value().is_some(),
                        "{s}: {:?}",
                        errors.samples[s]
                    );
                }
            }
        });
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
            let mut stored = positions_with_held_out_reads(lent, depth_cap, group_index, coverage);
            for position in &mut stored {
                position.libraries[silent] = LibraryAtPosition {
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
            let mut scratch = Scratch::new(group_index.len(), parameters.clean.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), group_index.len(), nodes);
            let mut here = PositionScores::new(group_index, parameters.clean.len());
            let mut information = InformationSums::new(group_index);
            for position in &stored {
                scratch.evidence.libraries.clone_from(&position.libraries);
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
                group_index,
                &config.edges,
                depth_cap,
                coverage,
                config.depth_as_a_range,
                0,
                EvidenceCursor::position_count(lent),
            );
            let nodes = model.quadrature.nodes.len();
            let mut scratch = Scratch::new(lent.len(), parameters.clean.len(), nodes);
            let mut statistics = Statistics::new(parameters.clean.len(), lent.len(), nodes);
            let mut here = PositionScores::new(group_index, parameters.clean.len());
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
            _ => match (index - COHORT_PARAMETERS) % ONE_LIBRARY_SAMPLE_PARAMETERS {
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
        let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
        let m = kept.len();
        let mut full = vec![0.0; m * m];
        let mut flat = vec![0.0; n];
        for here in every {
            flat[..COHORT_PARAMETERS].copy_from_slice(&here.cohort);
            for (s, row) in here.samples.iter().enumerate() {
                let first = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * s;
                flat[first..first + ONE_LIBRARY_SAMPLE_PARAMETERS].copy_from_slice(row);
            }
            for (p, &i) in kept.iter().enumerate() {
                for (q, &j) in kept.iter().enumerate() {
                    full[p * m + q] += flat[i] * flat[j];
                }
            }
        }
        full
    }

    /// The errors three ways, at a cohort's fitted maximum: from the blocks (what the fit reports
    /// above [`FULL_MATRIX_SAMPLES`] samples),
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
            let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
            let block_error = |i: usize| {
                if i < COHORT_PARAMETERS {
                    blocks.cohort[i]
                } else {
                    let at = i - COHORT_PARAMETERS;
                    blocks.samples[at / ONE_LIBRARY_SAMPLE_PARAMETERS]
                        [at % ONE_LIBRARY_SAMPLE_PARAMETERS]
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
                            && (i - COHORT_PARAMETERS) % ONE_LIBRARY_SAMPLE_PARAMETERS == which
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
                        i >= COHORT_PARAMETERS
                            && (i - COHORT_PARAMETERS) / ONE_LIBRARY_SAMPLE_PARAMETERS == s
                    })
                    .map(|i| (i - COHORT_PARAMETERS) % ONE_LIBRARY_SAMPLE_PARAMETERS)
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
                        reduced[a * k + b] = sums.sample_blocks[s]
                            [ja * ONE_LIBRARY_SAMPLE_PARAMETERS + jb]
                            - explained;
                    }
                }
                let inverse = inverse_of_positive_definite(&reduced, k).expect("sample inverts");
                for (a, &ja) in own.iter().enumerate() {
                    sketch[COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * s + ja] =
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
    /// samples), so the blocks are what the fit reports above [`FULL_MATRIX_SAMPLES`] samples, and
    /// the full matrix below (plan step A8).
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
    /// two rates keep theirs. **They are the whole matrix's**, which a cohort this small keeps: at
    /// four samples the blocks' differ.
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
                    let sums = statistics.information.expect("asked for");
                    assert!(
                        sums.full.is_some(),
                        "a cohort this small keeps the whole matrix"
                    );
                    (
                        StandardErrors::of(&sums),
                        StandardErrors::of_the_blocks(&sums),
                    )
                });
            let (at_the_returned_parameters, from_the_blocks) = at_the_returned_parameters;
            // `Debug` prints every `f64` as its shortest round-trip form: equal strings are equal bits.
            assert_eq!(
                format!("{:?}", fit.standard_errors),
                format!("{at_the_returned_parameters:?}"),
                "{samples} sample(s)"
            );
            // **The errors are the whole matrix's** (plan step A8): at four samples the blocks,
            // which leave two samples' pairings out, give others.
            if samples > 1 {
                assert_ne!(
                    format!("{:?}", fit.standard_errors),
                    format!("{from_the_blocks:?}"),
                    "{samples} samples: the fit's errors are the blocks'"
                );
            }
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
    /// on the four-sample fixture of the test above, with three starts, the rows after the
    /// winner's last pass are exactly the returned values and their errors, under the names the
    /// value rows use, with the final pass's log-likelihood; the winner's last pass left the
    /// returned values too (its last accelerated step was kept, not refused at the pass limit);
    /// and **no other start's last pass did**, so rows filed under another start are a different
    /// start's. Which start wins moves with the draw and with the fit's arithmetic, so the test
    /// reads it from the rows rather than naming it.
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
        assert_eq!(pass, fit.passes + 1);
        let values_of = |of: usize, at: u32| -> Vec<(String, f64)> {
            rows.iter()
                .filter(|row| (row.0, row.1) == (of, at) && !is_error(&row.3))
                .map(|row| (row.3.clone(), row.4))
                .collect()
        };
        let values_at = |at: u32| values_of(start, at);
        // `Debug` prints every `f64` as its shortest round-trip form: equal strings are equal bits.
        let returned = format!("{:?}", parameters.named());
        assert_eq!(format!("{:?}", values_at(pass)), returned);
        assert_eq!(format!("{:?}", values_at(fit.passes)), returned);
        let starts: std::collections::BTreeSet<usize> = rows.iter().map(|row| row.0).collect();
        assert_eq!(starts.len(), 3, "the fixture runs three starts");
        let ending_at_the_returned_values: Vec<usize> = starts
            .into_iter()
            .filter(|&of| {
                let last = rows
                    .iter()
                    .filter(|row| row.0 == of && !is_error(&row.3) && (of, row.1) != (start, pass))
                    .map(|row| row.1)
                    .max()
                    .expect("every start writes its passes");
                format!("{:?}", values_of(of, last)) == returned
            })
            .collect();
        assert_eq!(
            ending_at_the_returned_values,
            [start],
            "only the winner's last pass left the returned values"
        );
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
        let two_samples = [vec![0], vec![1]];
        let mut first = PositionScores::new(&two_samples, 2);
        let mut second = PositionScores::new(&two_samples, 2);
        for (i, score) in first.cohort.iter_mut().enumerate() {
            *score = 1.0 + i as f64;
        }
        first.samples = vec![vec![2.0, 3.0, 5.0], vec![7.0, 11.0, 13.0]];
        second.cohort = [0.5; COHORT_PARAMETERS];
        second.samples = vec![vec![-1.0, 4.0, 0.25], vec![6.0, -2.0, 3.0]];
        let sample_slot =
            |s: usize, which: usize| COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * s + which;
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
                let s = at / ONE_LIBRARY_SAMPLE_PARAMETERS;
                match at % ONE_LIBRARY_SAMPLE_PARAMETERS {
                    sample::CLEAN_ERROR_RATE => parameters.clean[s],
                    sample::NOISY_ERROR_RATE => parameters.noisy[s],
                    sample::HOMOZYGOTE_EXCESS => parameters.hom_excess[s],
                    _ => unreachable!("a sample of one library carries three parameters"),
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
        /// `(estimate − truth)` in units of the error the fit reports — the whole matrix's, since
        /// every regime here has at most [`FULL_MATRIX_SAMPLES`] samples (plan step A8).
        reported: Tally,
        /// The same in units of the blocks' error for the same cohort — the way a cohort above
        /// [`FULL_MATRIX_SAMPLES`] samples is given its errors — over the estimates that have both
        /// errors, so the two tallies compare like with like where the blocks give one.
        blocks: Tally,
        /// `(estimate + Newton step − truth)` in units of the reported error: where one Newton step
        /// — the whole matrix's inverse times the log-likelihood's slopes summed over the positions
        /// — puts the likelihood's maximum, against the truth. Over the same estimates as
        /// `reported`.
        at_flat_point: Tally,
        /// Parameters of this kind without a reported error.
        absent: usize,
    }

    impl Coverage {
        /// How far the likelihood's flat point lies from where the fit stopped, in errors, on
        /// average: the Newton step's mean, which is the two tallies' difference.
        fn flat_point(&self) -> f64 {
            self.at_flat_point.mean() - self.reported.mean()
        }
    }

    /// The entry at layout `slot` of a set of errors — the cohort's by slot, a sample's by kind.
    fn error_in(errors: &StandardErrors, slot: usize) -> StandardError {
        if slot < COHORT_PARAMETERS {
            errors.cohort[slot]
        } else {
            let at = slot - COHORT_PARAMETERS;
            errors.samples[at / ONE_LIBRARY_SAMPLE_PARAMETERS][at % ONE_LIBRARY_SAMPLE_PARAMETERS]
        }
    }

    /// Draw `cohorts` cohorts of `samples` samples at `mean_depth` from the true values, fit each by
    /// the fit's own rule, and tally how far each estimate lands from the truth in units of its
    /// error — the error the fit itself reports, and the blocks' — and how far the likelihood's
    /// maximum lies from it. Returns the tallies by kind, and how many fits converged.
    ///
    /// **The fit's reported errors are checked to be the whole matrix's**: the test sums every
    /// position's scores multiplied pairwise itself, inverts that over the parameters the fit gave
    /// an error, and requires every reported error to agree to 10⁻⁶ relative.
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
            let (fit, parameters) = fitted(&cohort, &config);
            converged += usize::from(fit.converged);
            let reported = &fit.standard_errors;
            with_sections(&cohort, None, |lent, depth_cap, group_index, coverage| {
                let (_, every) = scores_at_every_position(
                    lent,
                    depth_cap,
                    &config,
                    group_index,
                    coverage,
                    &parameters,
                );
                let blocks = StandardErrors::of_the_blocks(&information_by_hand(&every));
                let n = COHORT_PARAMETERS + ONE_LIBRARY_SAMPLE_PARAMETERS * samples;
                // The duplicated class is off: its three slots are not parameters of this fit.
                let fitted: Vec<usize> = (0..n)
                    .filter(|&slot| !(cohort::DUPLICATED_SHARE..COHORT_PARAMETERS).contains(&slot))
                    .collect();
                let kept: Vec<usize> = fitted
                    .iter()
                    .copied()
                    .filter(|&slot| {
                        matches!(
                            error_in(reported, slot),
                            StandardError::Estimated(_) | StandardError::WiderThanItsRange(_)
                        )
                    })
                    .collect();
                let m = kept.len();
                let full_inverse = inverse_of_positive_definite(&full_matrix(&every, &kept), m)
                    .expect("the whole matrix inverts over the parameters the fit gave an error");
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
                                    here.samples[at / ONE_LIBRARY_SAMPLE_PARAMETERS]
                                        [at % ONE_LIBRARY_SAMPLE_PARAMETERS]
                                }
                            })
                            .sum()
                    })
                    .collect();
                for &slot in &fitted {
                    let entry = tally.entry(label_of(slot)).or_default();
                    let off = value_at(&parameters, slot) - value_at(&truth, slot);
                    let Some(error) = error_in(reported, slot).value() else {
                        entry.absent += 1;
                        continue;
                    };
                    if let Some(block) = error_in(&blocks, slot).value() {
                        entry.blocks.add(off / block);
                    }
                    let p = kept.iter().position(|&k| k == slot).expect("kept");
                    let whole = full_inverse[p * m + p].sqrt();
                    assert!(
                        (error - whole).abs() < 1e-6 * whole,
                        "{samples} samples, slot {slot}: the fit reports {error}, the whole matrix \
                         gives {whole}"
                    );
                    let step: f64 = (0..m).map(|q| full_inverse[p * m + q] * slope[q]).sum();
                    entry.reported.add(off / error);
                    entry.at_flat_point.add((off + step) / error);
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
    /// At any cohort count: every kind is tallied, and every estimate the fit gives an error has one
    /// from the whole matrix the test sums itself, equal to it ([`coverage_of`]). At the full count
    /// only — the draws are fixed by their seeds, so the numbers move only when the code does — each
    /// regime is held to what it measured (plan step A8's run,
    /// `tmp/fit_precision/a8_coverage_full.log`), with room to spare:
    ///
    /// - at 4 and 20 samples, at most five fits in a regime stop at the pass limit (2 did, at 4
    ///   samples and 3 reads);
    /// - where the reported errors cover — **every kind** at 20 samples at both depths and at 4
    ///   samples and 30 reads, and the error rates and the homozygote excess at 4 samples and 3
    ///   reads — each lands between 0.58 and 0.78 within one error and between 0.90 and 0.99 within
    ///   two (measured 0.650 to 0.741, and 0.925 to 0.980);
    /// - **the fit stops at the likelihood's flat point**: at 20 samples every kind's flat point lies
    ///   within 0.02 errors of where the fit stopped (measured at most 0.002), and at 4 samples the
    ///   density's shapes' and the invariant share's within 0.25 (at most 0.108). With the shapes'
    ///   update that solved the digamma form, those three lay 0.65 to 1.35 errors away at 20 samples
    ///   and 1.8 to 2.3 at 4 (step A4's report, §2);
    /// - at 20 samples both density shapes sit within 0.25 errors of the truth on average (measured
    ///   at most 0.033), where the digamma form left them 1.05 to 1.29 above it;
    /// - at 4 samples and 3 reads the reported errors put at least 0.08 more of the error rates
    ///   within one error than the blocks' do (measured 0.150 and 0.148 more).
    ///
    /// The 2-sample regime is printed and held only to the checks that apply at any count: 54 of its
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
                coverage.reported.count > 0
                    && coverage.at_flat_point.count == coverage.reported.count,
                "{regime}, {kind}: {} estimates with a reported error",
                coverage.reported.count,
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
        let covering: &[&str] = match (samples, depth) {
            (20, _) | (4, 30) => &COVERAGE_KINDS,
            (4, 3) => &[
                "clean error rates",
                "mismapped error rates",
                "homozygote excess",
            ],
            _ => &[],
        };
        for &kind in covering {
            let reported = &tally[kind].reported;
            assert!(
                (0.58..=0.78).contains(&reported.share_within_one())
                    && (0.90..=0.99).contains(&reported.share_within_two()),
                "{regime}, {kind}: {:.3} within one of the reported errors, {:.3} within two",
                reported.share_within_one(),
                reported.share_within_two()
            );
        }
        let flat_point_bound = match samples {
            20 => Some((COVERAGE_KINDS.as_slice(), 0.02)),
            4 => Some((
                ["density shape a", "density shape b", "invariant share"].as_slice(),
                0.25,
            )),
            _ => None,
        };
        if let Some((kinds, bound)) = flat_point_bound {
            for &kind in kinds {
                let flat_point = tally[kind].flat_point();
                assert!(
                    flat_point.abs() < bound,
                    "{regime}, {kind}: the likelihood's flat point lies {flat_point:+.3} errors from \
                     where the fit stopped"
                );
            }
        }
        if samples == 20 {
            for kind in ["density shape a", "density shape b"] {
                let off_centre = tally[kind].reported.mean();
                assert!(
                    off_centre.abs() < 0.25,
                    "{regime}, {kind}: the estimates sit {off_centre:+.3} errors from the truth"
                );
            }
        }
        if samples == 4 && depth == 3 {
            for kind in ["clean error rates", "mismapped error rates"] {
                let (reported, blocks) = (&tally[kind].reported, &tally[kind].blocks);
                assert!(
                    reported.share_within_one() >= blocks.share_within_one() + 0.08,
                    "{regime}, {kind}: within one error {:.3} (reported) against {:.3} (blocks)",
                    reported.share_within_one(),
                    blocks.share_within_one()
                );
            }
        }
    }

    /// **The errors mean what they say** (plan step A4, spec §3.6 item 3): 200 cohorts drawn from
    /// known values in each regime, each fitted by the fit's own rule; of every estimate with an
    /// error, the share within one error of the truth and within two, which for an error that means
    /// what it says are about 68 in 100 and 95 in 100, and the spread of the distances, about one.
    /// Measured for the errors the fit reports — the whole matrix's, since every regime has at most
    /// [`FULL_MATRIX_SAMPLES`] samples (plan step A8) — and for the blocks', which a larger cohort
    /// would be given, on 4 samples over 20,000 positions and 20 over 5,000 at 3 and at 30 reads a
    /// position, and on 2 samples over 30,000 at 3 reads. The duplicated class is off in the draw and
    /// the fit ([`TRUTH`]), so nothing here checks where the carrier Beta's shapes come to rest: the
    /// slopes they step on are held to the scorer's by `the_passs_shape_slopes_are_the_likelihoods`
    /// and `the_passs_carrier_slopes_follow_the_coverage_odds`.
    /// Also printed, per kind: how far the log-likelihood's own flat point lies from where the fit
    /// stopped, and the coverage there — both by one Newton step.
    ///
    /// What it measured, and what it holds each regime to, are in [`check_the_coverage`] and the
    /// reports of plan steps A4, A7 and A8.
    ///
    /// Ignored by default: some thirty-five minutes in the container. Run with
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
                let (reported, blocks, at_flat) = (
                    &coverage.reported,
                    &coverage.blocks,
                    &coverage.at_flat_point,
                );
                eprintln!(
                    "COVERAGE {samples} samples, {mean_depth} reads, {kind}: {} estimates ({} \
                         with a blocks' error), {} without an error; within one error {:.3} \
                         (reported) {:.3} (blocks), within two {:.3} {:.3}; mean distance {:+.3} \
                         {:+.3}; spread {:.3} {:.3}; the likelihood's flat point lies {:+.3} errors \
                         from the fit, on average; there, within one {:.3}, within two {:.3}, mean \
                         {:+.3}",
                    reported.count,
                    blocks.count,
                    coverage.absent,
                    reported.share_within_one(),
                    blocks.share_within_one(),
                    reported.share_within_two(),
                    blocks.share_within_two(),
                    reported.mean(),
                    blocks.mean(),
                    reported.spread(),
                    blocks.spread(),
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
