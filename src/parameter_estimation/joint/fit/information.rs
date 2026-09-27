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

#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "only this module's tests call the scorer until the expectation pass does, at \
                  plan step A2"
    )
)]

use crate::float;
use crate::types::Ploidy;

use super::{
    BRANCHES, BetaQuadrature, CANDIDATE_ALTERNATIVES, MAX_CANDIDATES, Parameters, Scratch,
    candidate_read_probability, candidate_read_probability_slope, class_rate,
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
/// rule at each shape nudged up and down.
#[derive(Clone, Debug)]
pub(super) struct RuleSlopes {
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
            node_slope: [node_a, node_b],
            ln_weight_slope: [ln_weight_a, ln_weight_b],
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
/// pass allocates it once per chunk. It also holds the per-sample working numbers the error
/// rates' slopes are built from, for the same reason.
#[derive(Clone, Debug)]
pub(super) struct PositionScores {
    /// The cohort's parameters, indexed by [`cohort`].
    pub cohort: [f64; COHORT_PARAMETERS],
    /// Per sample, in the order the fit iterates samples, indexed by [`sample`].
    pub samples: Vec<[f64; SAMPLE_PARAMETERS]>,
    /// Per sample and noise class, the parts of a read's log-probability slope in the error rate
    /// that depend on the genotype but not on the candidate allele (see [`ReadSlopes`]).
    slopes: Vec<[ReadSlopes; 2]>,
}

/// The pieces of `d ln P(this sample's reads | genotype) / d ε` that do not depend on which base
/// the alternative is, for one sample, one noise class and each of the three genotypes.
///
/// The slope a whole genotype contributes is
/// `candidate reads · candidate[j] + other reads / ε + reference[j]`: the candidate allele's reads
/// each move by the slope of their log-probability, the reads on the two bases that are neither
/// move by `1/ε` (each is `ε/3`), and the reference reads move by their expected count times
/// their own log-probability's slope — the expected count, because the depth a stored code
/// stands for is a range and the likelihood sums over it.
#[derive(Copy, Clone, Debug, Default)]
struct ReadSlopes {
    /// `d ln P(a read shows the candidate) / d ε`, by how many copies of it the sample carries.
    candidate: [f64; 3],
    /// The expected reference reads times `d ln P(a read shows the reference) / d ε`.
    reference: [f64; 3],
    /// `1/ε` — the slope of `ln(ε/3)`, which every read on a base that is neither shares.
    per_other_read: f64,
}

impl PositionScores {
    pub(super) fn new(samples: usize) -> Self {
        Self {
            cohort: [0.0; COHORT_PARAMETERS],
            samples: vec![[0.0; SAMPLE_PARAMETERS]; samples],
            slopes: vec![[ReadSlopes::default(); 2]; samples],
        }
    }

    fn clear(&mut self) {
        self.cohort = [0.0; COHORT_PARAMETERS];
        self.samples.fill([0.0; SAMPLE_PARAMETERS]);
    }
}

/// Every parameter's slope at the position [`one_position`](super::one_position) has just
/// scored, written into `into`.
///
/// **Reads the scratch `one_position` left behind and changes nothing in it**, so calling it or
/// not calling it cannot move a fitted number. A position whose likelihood underflowed carries
/// no slope at all: it contributes nothing to any parameter, exactly as it contributes nothing to
/// the counts.
///
/// **No posterior is dropped for being small.** When the pass credits each branch and node its
/// share of the position's counts, it skips any below one part in 10¹², which costs those counts
/// nothing; a slope divides by a share that
/// can itself be 10⁻¹², so here every posterior is kept. Each is formed as one exponential of a
/// difference of logarithms, which never divides by a share at all.
#[allow(
    clippy::too_many_arguments,
    reason = "one position's whole problem, as `one_position` takes it"
)]
#[allow(
    clippy::needless_range_loop,
    reason = "the sample index addresses the scratch, the read-group map and the score rows at \
              once, as in `one_position`"
)]
pub(super) fn score_position(
    scratch: &Scratch,
    group_index: &[Vec<usize>],
    ploidy: Ploidy,
    density_rule: SlopedRule<'_>,
    carrier_rule: Option<SlopedRule<'_>>,
    coverage_odds: &[f64],
    parameters: &Parameters,
    into: &mut PositionScores,
) {
    let quadrature = density_rule.rule;
    let carrier = carrier_rule.map(|carrier| carrier.rule);
    into.clear();
    let position_ln = scratch.position_ln;
    if !position_ln.is_finite() {
        return;
    }
    let samples = scratch.samples;
    let nodes = scratch.nodes;
    let candidates = scratch.candidates.len();
    let ln = |p: f64| float::ln(p.max(f64::MIN_POSITIVE));
    let ln_three = float::ln(CANDIDATE_ALTERNATIVES as f64);

    let density = &parameters.density;
    let duplicated_share = parameters.duplicated.map_or(0.0, |d| d.share);
    let ln_ordinary = ln(1.0 - duplicated_share);
    let ln_invariant = ln(density.p_invariant);
    let ln_fixed_alt = ln(density.p_fixed_alt);
    let ln_segregating = ln(density.p_segregating());
    let ln_duplicated = ln(duplicated_share);

    // ---- the per-sample slopes the error rates are built from -------------------------------
    for s in 0..samples {
        let sample = &scratch.evidence.samples[s];
        let weights = scratch.evidence.weights_of(s);
        for class in 0..2 {
            let rate = class_rate(parameters, class, group_index[s][0]);
            let slopes = &mut into.slopes[s][class];
            slopes.per_other_read = 1.0 / rate;
            for j in 0..3_u8 {
                let on_candidate = candidate_read_probability(j, ploidy, rate);
                slopes.candidate[usize::from(j)] = candidate_read_probability_slope(j, ploidy)
                    / on_candidate.max(f64::MIN_POSITIVE);
                let on_reference = reference_read_probability(j, ploidy, rate);
                let expected = expected_reference_reads(sample, weights, on_reference);
                slopes.reference[usize::from(j)] = expected
                    * reference_read_probability_slope(j, ploidy)
                    / on_reference.max(f64::MIN_POSITIVE);
            }
        }
    }
    // The slope of one sample's reads under genotype `j` with candidate allele `allele`.
    let read_slope = |into: &PositionScores, s: usize, class: usize, j: usize, allele: usize| {
        let sample = &scratch.evidence.samples[s];
        let slopes = &into.slopes[s][class];
        let candidate_reads = sample.on[allele];
        let other_reads = sample.non_reference() - candidate_reads - sample.on[4];
        candidate_reads * slopes.candidate[j]
            + other_reads * slopes.per_other_read
            + slopes.reference[j]
    };

    // A class's error rate is that class's slot in the sample's row: the clean rate is class
    // zero's, the mismapped rate class one's. The slopes are accumulated straight into it.
    const _: () = assert!(sample::CLEAN_ERROR_RATE == 0 && sample::NOISY_ERROR_RATE == 1);

    for class in 0..2 {
        let share = if class == 0 {
            1.0 - parameters.noisy_share
        } else {
            parameters.noisy_share
        };
        let ln_class = ln(share);
        // Each branch's posterior with its own share divided back out, which is exactly the
        // slope of the position's log-likelihood in that share.
        let ln_posterior_per_share = ln_class + ln_ordinary - position_ln;

        // ---- the class share, and the density's two masses ------------------------------------
        let within = ln_sum_exp(&scratch.branch_ln[class * BRANCHES..][..BRANCHES]);
        let slope_of_class = float::exp(within - position_ln);
        into.cohort[cohort::NOISY_SHARE] += if class == 1 {
            slope_of_class
        } else {
            -slope_of_class
        };
        let on_segregating_per_share =
            float::exp(ln_posterior_per_share + scratch.segregating_ln[class]);
        into.cohort[cohort::P_INVARIANT] +=
            float::exp(ln_posterior_per_share + scratch.invariant_ln[class])
                - on_segregating_per_share;
        into.cohort[cohort::P_FIXED_ALT] +=
            float::exp(ln_posterior_per_share + scratch.fixed_alt_ln[class])
                - on_segregating_per_share;
        if carrier.is_some() {
            let ordinary = ln_sum_exp(&[
                ln_invariant + scratch.invariant_ln[class],
                ln_fixed_alt + scratch.fixed_alt_ln[class],
                ln_segregating + scratch.segregating_ln[class],
            ]);
            into.cohort[cohort::DUPLICATED_SHARE] +=
                float::exp(ln_class + scratch.duplicated_ln[class] - position_ln)
                    - float::exp(ln_class + ordinary - position_ln);
        }

        // ---- the invariant branch: every sample homozygous reference ---------------------------
        let invariant =
            float::exp(ln_posterior_per_share + ln_invariant + scratch.invariant_ln[class]);
        let any_allele = scratch.candidates[0];
        for s in 0..samples {
            let slope = read_slope(into, s, class, 0, any_allele);
            into.samples[s][class] += invariant * slope;
        }

        // ---- the fixed branch: every sample carries two copies of the candidate ----------------
        for candidate in 0..candidates {
            let allele = scratch.candidates[candidate];
            let fixed = float::exp(
                ln_posterior_per_share
                    + ln_fixed_alt
                    + scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                    + float::ln(scratch.multiplicity[candidate])
                    - ln_three,
            );
            for s in 0..samples {
                let slope = read_slope(into, s, class, 2, allele);
                into.samples[s][class] += fixed * slope;
            }
        }

        // ---- the segregating branch, node by node ----------------------------------------------
        for candidate in 0..candidates {
            let allele = scratch.candidates[candidate];
            let multiplicity = float::ln(scratch.multiplicity[candidate]);
            for node in 0..nodes {
                let at_node = float::exp(
                    ln_posterior_per_share
                        + ln_segregating
                        + scratch.node_ln[scratch.node_at(class, candidate, node)]
                        + multiplicity
                        - ln_three,
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
                        genotype_frequencies_slope_in_frequency(f, parameters.hom_excess[s]);
                    along_frequency += (prior_along_f[0] * lik[0]
                        + prior_along_f[1] * lik[1]
                        + prior_along_f[2] * lik[2])
                        / total;
                    let mut reads = 0.0;
                    for (j, weight) in joint.iter().enumerate() {
                        reads += weight * read_slope(into, s, class, j, allele);
                    }
                    into.samples[s][class] += at_node * reads / total;
                }
                for (shape, slot) in [cohort::DENSITY_A, cohort::DENSITY_B]
                    .into_iter()
                    .enumerate()
                {
                    let slopes = density_rule.slopes;
                    into.cohort[slot] += at_node
                        * (slopes.ln_weight_slope[shape][node]
                            + along_frequency * slopes.node_slope[shape][node]);
                }
            }
        }

        // ---- the duplicated branch, node by node over the carrier frequency --------------------
        if let Some(carrier_rule) = carrier_rule {
            let carrier = carrier_rule.rule;
            for candidate in 0..candidates {
                let allele = scratch.candidates[candidate];
                let multiplicity = float::ln(scratch.multiplicity[candidate]);
                for node in 0..carrier.nodes.len() {
                    let at_node = float::exp(
                        ln_class + ln_duplicated - position_ln
                            + scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                            + multiplicity
                            - ln_three,
                    );
                    if at_node == 0.0 {
                        continue;
                    }
                    let q = carrier.nodes[node];
                    let mut along_frequency = 0.0;
                    for s in 0..samples {
                        let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                        let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                        let joint = [(1.0 - q) * lik[0], q * odds * lik[1]];
                        let total = joint[0] + joint[1];
                        if total <= 0.0 {
                            continue;
                        }
                        along_frequency += (odds * lik[1] - lik[0]) / total;
                        let reads = joint[0] * read_slope(into, s, class, 0, allele)
                            + joint[1] * read_slope(into, s, class, 1, allele);
                        into.samples[s][class] += at_node * reads / total;
                    }
                    for (shape, slot) in [cohort::CARRIER_A, cohort::CARRIER_B]
                        .into_iter()
                        .enumerate()
                    {
                        let slopes = carrier_rule.slopes;
                        into.cohort[slot] += at_node
                            * (slopes.ln_weight_slope[shape][node]
                                + along_frequency * slopes.node_slope[shape][node]);
                    }
                }
            }
        }
    }
}

/// A quadrature rule beside how its nodes and weights move with its two shapes — the pair a
/// shape's slope is read from.
#[derive(Copy, Clone)]
pub(super) struct SlopedRule<'a> {
    pub rule: &'a BetaQuadrature,
    pub slopes: &'a RuleSlopes,
}

#[cfg(test)]
mod tests {
    use crate::parameter_estimation::joint::census::{DepthCap, SampleGenericSections};
    use crate::parameter_estimation::joint::fit::bench_fixtures::{
        DrawnCohort, as_cohort, draw_cohort_with_duplications,
    };
    use crate::parameter_estimation::joint::fit::{
        DuplicatedPositions, EvidenceCursor, FrequencyDensity, JointFitConfig, ReadLogs,
        SampleAtPosition, Statistics, candidate_read_probability, digamma, expectation_pass,
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
    /// and every parameter's slope summed over the positions.
    fn summed_scores(
        lent: &[SampleGenericSections<'_>],
        depth_cap: DepthCap,
        config: &JointFitConfig,
        group_index: &[Vec<usize>],
        coverage: &[f64],
        parameters: &Parameters,
    ) -> (f64, PositionScores) {
        let quadrature = BetaQuadrature::with_genotype_priors(
            parameters.density.a,
            parameters.density.b,
            config.quadrature_nodes,
            &parameters.hom_excess,
        );
        let density_slopes = RuleSlopes::of(
            parameters.density.a,
            parameters.density.b,
            config.quadrature_nodes,
        );
        let carrier = parameters
            .duplicated
            .map(|d| BetaQuadrature::new(d.carrier_a, d.carrier_b, config.quadrature_nodes));
        let carrier_slopes = parameters
            .duplicated
            .map(|d| RuleSlopes::of(d.carrier_a, d.carrier_b, config.quadrature_nodes));
        let read_logs: Vec<Vec<ReadLogs>> = (0..2)
            .map(|class| {
                (0..parameters.clean.len())
                    .map(|group| ReadLogs::of(class_rate(parameters, class, group), config.ploidy))
                    .collect()
            })
            .collect();
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
        let mut scratch = Scratch::new(lent.len(), quadrature.nodes.len());
        let mut statistics =
            Statistics::new(parameters.clean.len(), lent.len(), quadrature.nodes.len());
        let mut here = PositionScores::new(lent.len());
        let mut total = PositionScores::new(lent.len());
        while cursor.next_position(&mut scratch.evidence) {
            one_position(
                &mut scratch,
                group_index,
                config.ploidy,
                &read_logs,
                &quadrature,
                carrier.as_ref(),
                &[],
                parameters,
                &mut statistics,
            );
            score_position(
                &scratch,
                group_index,
                config.ploidy,
                SlopedRule {
                    rule: &quadrature,
                    slopes: &density_slopes,
                },
                carrier
                    .as_ref()
                    .zip(carrier_slopes.as_ref())
                    .map(|(rule, slopes)| SlopedRule { rule, slopes }),
                &[],
                parameters,
                &mut here,
            );
            for (into, from) in total.cohort.iter_mut().zip(&here.cohort) {
                *into += from;
            }
            for (into, from) in total.samples.iter_mut().zip(&here.samples) {
                for (into, from) in into.iter_mut().zip(from) {
                    *into += from;
                }
            }
        }
        (statistics.log_likelihood, total)
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
            false,
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
                    false,
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
        const NODES: usize = 12;
        let ploidy = JointFitConfig::default().ploidy;
        let samples = group_index.len();
        let quadrature = BetaQuadrature::with_genotype_priors(
            parameters.density.a,
            parameters.density.b,
            NODES,
            &parameters.hom_excess,
        );
        let density_slopes = RuleSlopes::of(parameters.density.a, parameters.density.b, NODES);
        let carrier = parameters
            .duplicated
            .map(|d| BetaQuadrature::new(d.carrier_a, d.carrier_b, NODES));
        let carrier_slopes = parameters
            .duplicated
            .map(|d| RuleSlopes::of(d.carrier_a, d.carrier_b, NODES));
        let read_logs: Vec<Vec<ReadLogs>> = (0..2)
            .map(|class| {
                (0..parameters.clean.len())
                    .map(|group| ReadLogs::of(class_rate(parameters, class, group), ploidy))
                    .collect()
            })
            .collect();
        let mut scratch = Scratch::new(samples, NODES);
        let mut statistics = Statistics::new(parameters.clean.len(), samples, NODES);
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
            one_position(
                &mut scratch,
                group_index,
                ploidy,
                &read_logs,
                &quadrature,
                carrier.as_ref(),
                odds,
                parameters,
                &mut statistics,
            );
            total += scratch.position_ln;
            if let Some(into) = scores.as_deref_mut() {
                score_position(
                    &scratch,
                    group_index,
                    ploidy,
                    SlopedRule {
                        rule: &quadrature,
                        slopes: &density_slopes,
                    },
                    carrier
                        .as_ref()
                        .zip(carrier_slopes.as_ref())
                        .map(|(rule, slopes)| SlopedRule { rule, slopes }),
                    odds,
                    parameters,
                    &mut here,
                );
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
