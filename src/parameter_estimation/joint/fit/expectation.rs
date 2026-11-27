//! The expectation step: one pass over every census position at fixed parameters, which computes
//! each position's likelihood and adds its share to every count the maximisations read. Also the
//! read model the pass scores a library's reads with, which the maximisation and the standard
//! errors' scores ([`information`](super::information)) read too.

use crate::float;
use crate::types::Ploidy;

use super::super::census::{DepthCap, SampleGenericSections};
use super::beta_quadrature::BetaQuadrature;
use super::information::{
    InformationSums, PositionScores, RuleSlopes, ScoringTables, score_position,
};
use super::{
    EvidenceCursor, JointFitConfig, LibraryAtPosition, MAX_RECORDED_SPREAD, Parameters,
    PositionEvidence, ReadTally, Statistics,
};

/// How many non-reference bases a position has, whatever the reference base is: three.
const CANDIDATE_ALTERNATIVES: usize = 3;

// ---------------------------------------------------------------------
// The likelihood of one library's reads under one genotype
// ---------------------------------------------------------------------

/// `ln P(this library's reads | j of the sample's copies carry base k, the library's error rate
/// ε)`. A sample's reads under a genotype are the sum of this over its libraries: the genotype is
/// the sample's, the error rate each library's own.
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
/// [`alternative_read_probability`](crate::parameter_estimation::generic::noise_model)
/// exactly — the two routes score a read the same way, which is what makes the comparison
/// between them a comparison.
/// **The reference-read term and the non-reference total are handed in rather than taken
/// here**, because neither depends on which base the alternative is. [`ln_reference_reads`]
/// reads the library, its depth weights and how many copies of the alternative the sample
/// carries — and it is the one logarithm in this block — so the caller takes it once for the
/// three copy counts and reuses it across every candidate allele and across the invariant branch
/// beside them. Taken inside this function it was ten calls a library where three do.
fn ln_reads_given_genotype(
    library: &LibraryAtPosition,
    non_reference: f64,
    alt_copies: u8,
    alternative: usize,
    logs: &ReadLogs,
    reference_term: f64,
) -> f64 {
    let copies = alt_copies as usize;
    let candidate_reads = library.on[alternative];
    // `Other` is neither the candidate nor the reference; it is held out of the model
    // entirely, so the depth it occupied is removed with it.
    let other_reads = non_reference - candidate_reads - library.on[4];

    // **A count of zero needs no guard now.** `count_times_ln` branched on it only to keep
    // `0 · −∞` out of the sum; the logarithms below are clamped at `MIN_POSITIVE` when they are
    // built, so every one of them is finite and `0 · finite` is zero.
    candidate_reads * logs.ln_candidate[copies] + other_reads * logs.ln_neither + reference_term
}

/// `ln P(this library's reads | every copy is the reference base)` — the term the invariant
/// branch needs, and it does not depend on which base would have been the alternative. Its
/// `reference_term` is the zero-copy entry of the same triple.
fn ln_reads_given_all_reference(
    library: &LibraryAtPosition,
    non_reference: f64,
    logs: &ReadLogs,
    reference_term: f64,
) -> f64 {
    (non_reference - library.on[4]) * logs.ln_neither + reference_term
}

/// `ln P(the reference reads | the sample carries 0, 1 or 2 copies of the alternative)`.
///
/// **Candidate-invariant, which is why it is a function of its own.** Nothing in
/// [`ln_reference_reads`] reads *which* base the alternative is — only the library's reference
/// count, the depths its stored code allows and the per-copy probability from its read group's
/// table — so one triple serves all three candidate alleles and the invariant branch.
fn reference_terms(
    library: &LibraryAtPosition,
    depth_weights: &[f64],
    logs: &ReadLogs,
) -> [f64; 3] {
    [
        ln_reference_reads(
            library,
            depth_weights,
            logs.reference[0],
            logs.ln_reference[0],
        ),
        ln_reference_reads(
            library,
            depth_weights,
            logs.reference[1],
            logs.ln_reference[1],
        ),
        ln_reference_reads(
            library,
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
        let ln = |p: f64| float::ln(p.max(f64::MIN_POSITIVE));
        let mut ln_candidate = [0.0; 3];
        let mut reference = [0.0; 3];
        let mut ln_reference = [0.0; 3];
        for copies in 0..3_usize {
            let on_candidate = candidate_read_probability(copies as u8, ploidy, error_rate);
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
pub(super) fn fill_depth_weights(weights: &mut [f64], coverage: f64, fewest_reference: f64) {
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
pub(super) fn reference_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    (1.0 - carried) * (1.0 - error_rate) + carried * error_rate / 3.0
}

/// How [`reference_read_probability`] moves with the error rate. Kept beside it so the standard
/// errors' slopes (`information`) follow any change to the read model.
pub(super) fn reference_read_probability_slope(alt_copies: u8, ploidy: Ploidy) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    -(1.0 - carried) + carried / 3.0
}

/// The chance one read shows the candidate allele, from an individual carrying `alt_copies`
/// copies of it — the probability [`ReadLogs::of`] takes the logarithm of.
pub(super) fn candidate_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    carried * (1.0 - error_rate) + (1.0 - carried) * error_rate / 3.0
}

/// How [`candidate_read_probability`] moves with the error rate.
pub(super) fn candidate_read_probability_slope(alt_copies: u8, ploidy: Ploidy) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    -carried + (1.0 - carried) / 3.0
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
pub(super) fn expected_reference_reads(
    library: &LibraryAtPosition,
    depth_weights: &[f64],
    on_reference: f64,
) -> f64 {
    if library.spread <= 1 {
        return library.fewest_reference;
    }
    let probability = on_reference.max(f64::MIN_POSITIVE);
    let (mut total, mut weighted, mut power) = (0.0, 0.0, 1.0);
    for (step, weight) in depth_weights[..library.spread].iter().enumerate() {
        let term = weight * power;
        total += term;
        weighted += term * (library.fewest_reference + step as f64);
        power *= probability;
    }
    if total > 0.0 {
        weighted / total
    } else {
        library.fewest_reference
    }
}

/// `ln P(the reference reads | a read shows the reference base with probability p)`, summed
/// over every depth the library's stored code could have come from.
///
/// The shallowest depth's term is factored out, which is what keeps the sum in range: what is
/// left runs from one upwards and cannot underflow however small `p` is.
fn ln_reference_reads(
    library: &LibraryAtPosition,
    depth_weights: &[f64],
    on_reference: f64,
    ln_on_reference: f64,
) -> f64 {
    let shallowest = library.fewest_reference * ln_on_reference;
    if library.spread <= 1 {
        return shallowest;
    }
    let probability = on_reference.max(f64::MIN_POSITIVE);
    let mut total = 0.0;
    let mut power = 1.0;
    for weight in &depth_weights[..library.spread] {
        total += weight * power;
        power *= probability;
    }
    shallowest + float::ln(total)
}

/// `count · ln p`, with a count of zero contributing nothing — otherwise a category no read
/// fell into, at probability zero, would be `0 · −∞ = NaN`, and a `NaN` does not lose loudly.
pub(super) fn count_times_ln(count: f64, probability: f64) -> f64 {
    if count == 0.0 {
        0.0
    } else {
        count * float::ln(probability.max(f64::MIN_POSITIVE))
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
pub(super) fn genotype_frequencies(f: f64, excess: f64) -> [f64; 3] {
    let heterozygous = 2.0 * f * (1.0 - f) * (1.0 - excess);
    let shift = excess * f * (1.0 - f);
    [(1.0 - f) * (1.0 - f) + shift, heterozygous, f * f + shift]
}

/// How [`genotype_frequencies`] moves as the frequency `f` moves, genotype by genotype — what the
/// standard errors' slopes in the Beta shapes read (`information`). **Kept beside the function it
/// differentiates**, so a change to one is made to both; a test compares the two by a finite
/// difference.
pub(super) fn genotype_frequencies_slope_in_frequency(f: f64, excess: f64) -> [f64; 3] {
    [
        -2.0 * (1.0 - f) + excess * (1.0 - 2.0 * f),
        2.0 * (1.0 - 2.0 * f) * (1.0 - excess),
        2.0 * f + excess * (1.0 - 2.0 * f),
    ]
}

/// How [`genotype_frequencies`] moves as the homozygote excess moves: the heterozygote loses
/// `2f(1 − f)` and each homozygote gains `f(1 − f)`. Kept beside it for the same reason as its
/// sibling above.
pub(super) fn genotype_frequencies_slope_in_excess(f: f64) -> [f64; 3] {
    let moved = f * (1.0 - f);
    [moved, -2.0 * moved, moved]
}

/// The read counts one library contributes to its read group's tally, under one candidate allele
/// and one genotype.
///
/// The counts on the candidate allele and on the two bases that are neither are exact; the
/// reference count is the one that has to be taken in expectation, because the depth it is
/// derived from is a range.
fn tally_of(
    library: &LibraryAtPosition,
    depth_weights: &[f64],
    alternative: usize,
    on_reference: f64,
) -> ReadTally {
    let candidate = library.on[alternative];
    ReadTally {
        candidate,
        neither: library.non_reference() - candidate - library.on[4],
        reference: expected_reference_reads(library, depth_weights, on_reference),
    }
}

/// **What a pass holds fixed while it visits positions**: the parameters, which read group each
/// sample's reads are scored under, and the tables built from the parameters once a pass — each
/// read group's read log-probabilities and the two quadrature rules.
///
/// One value, handed both to [`one_position`] and to the scorer that reads what it leaves behind
/// ([`information::score_position`](super::information::score_position)), so the two cannot be
/// given different inputs.
pub(super) struct PassModel<'a> {
    /// Per sample, the read groups its sections are, in the order the cursor visits them.
    pub(super) group_index: &'a [Vec<usize>],
    pub(super) ploidy: Ploidy,
    /// `[class][read group]`.
    read_logs: Vec<Vec<ReadLogs>>,
    /// The frequency density's rule, with each sample's genotype priors at its nodes.
    pub(super) quadrature: BetaQuadrature,
    /// The carrier Beta's rule and its slopes, when the run fits the duplicated class.
    pub(super) carrier: Option<CarrierRule>,
    /// How the density's rule's nodes and weights move with its two shapes — what the
    /// likelihood's slope in the shapes is made of, taken every pass for the shapes' update.
    pub(super) density_slopes: RuleSlopes,
    /// The log of each branch's and each noise class's share, and the `ln 3` and `ln k` the
    /// candidates are weighted by.
    pub(super) branches: BranchShares,
    pub(super) parameters: &'a Parameters,
}

/// **The logarithms a position's branches are weighted by, which depend on the parameters
/// alone**: each branch's share, each noise class's share, and how many alleles a candidate
/// stands for. The same at every position of a pass, so they are taken once a pass, beside the
/// read tables. Each share is floored at `f64::MIN_POSITIVE` before its logarithm is taken.
#[derive(Debug)]
pub(super) struct BranchShares {
    /// Whether the run fits the duplicated class.
    pub(super) duplicated_fitted: bool,
    /// `ln(1 − duplicated share)`: what the three ordinary branches share between them.
    pub(super) ln_ordinary: f64,
    /// The log of the invariant, fixed-non-reference, segregating and duplicated shares.
    pub(super) ln_invariant: f64,
    pub(super) ln_fixed_alt: f64,
    pub(super) ln_segregating: f64,
    pub(super) ln_duplicated: f64,
    /// The log of the three non-reference bases a fixed or segregating branch spreads over.
    pub(super) ln_three: f64,
    /// `[class]` — the log of the ordinary (0) and the mismapped (1) class's share.
    pub(super) ln_class: [f64; 2],
    /// `ln k`, indexed by how many alleles `k` a candidate stands for (one to three).
    ln_multiplicity: [f64; CANDIDATE_ALTERNATIVES + 1],
}

impl BranchShares {
    fn of(parameters: &Parameters, duplicated_fitted: bool) -> Self {
        let ln = |p: f64| float::ln(p.max(f64::MIN_POSITIVE));
        let density = &parameters.density;
        let duplicated_share = parameters.duplicated.map_or(0.0, |d| d.share);
        Self {
            duplicated_fitted,
            ln_ordinary: ln(1.0 - duplicated_share),
            ln_invariant: ln(density.p_invariant),
            ln_fixed_alt: ln(density.p_fixed_alt),
            ln_segregating: ln(density.p_segregating()),
            ln_duplicated: ln(duplicated_share),
            ln_three: float::ln(CANDIDATE_ALTERNATIVES as f64),
            ln_class: [ln(1.0 - parameters.noisy_share), ln(parameters.noisy_share)],
            ln_multiplicity: std::array::from_fn(|k| float::ln(k as f64)),
        }
    }
}

/// **The carrier Beta's rule and how it moves with the carrier's two shapes**, held together so
/// that a pass cannot have one without the other: a rule without its slopes would leave the
/// carrier's shapes where they started, and nothing would say so.
pub(super) struct CarrierRule {
    pub(super) quadrature: BetaQuadrature,
    pub(super) slopes: RuleSlopes,
}

impl<'a> PassModel<'a> {
    pub(super) fn new(
        parameters: &'a Parameters,
        config: &JointFitConfig,
        group_index: &'a [Vec<usize>],
    ) -> Self {
        let quadrature = BetaQuadrature::with_genotype_priors(
            parameters.density.a,
            parameters.density.b,
            config.quadrature_nodes,
            &parameters.hom_excess,
        );
        let carrier = parameters.duplicated.map(|duplicated| {
            let quadrature = BetaQuadrature::new(
                duplicated.carrier_a,
                duplicated.carrier_b,
                config.quadrature_nodes,
            );
            let slopes = RuleSlopes::of(
                duplicated.carrier_a,
                duplicated.carrier_b,
                quadrature.nodes.len(),
            );
            CarrierRule { quadrature, slopes }
        });
        // One table a noise class a read group, beside the quadratures and for the same reason.
        let read_logs: Vec<Vec<ReadLogs>> = (0..2)
            .map(|class| {
                (0..parameters.clean.len())
                    .map(|group| ReadLogs::of(class_rate(parameters, class, group), config.ploidy))
                    .collect()
            })
            .collect();
        let density_slopes = RuleSlopes::of(
            parameters.density.a,
            parameters.density.b,
            quadrature.nodes.len(),
        );
        Self {
            group_index,
            ploidy: config.ploidy,
            read_logs,
            branches: BranchShares::of(parameters, carrier.is_some()),
            quadrature,
            carrier,
            density_slopes,
            parameters,
        }
    }
}

/// One pass over every position: the posteriors, and every count the maximisations need.
///
/// **Split across cores by position.** Positions are independent given the parameters, and a
/// chunk's counts add to another chunk's, so the pass is a map and a sum. The fit's own passes
/// call [`expectation_pass`], which may also keep the information; this plain one is the tests'.
#[cfg(test)]
pub(super) fn expectation(
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
        PassKeeps::SUMS_ONLY,
    )
}

/// **What a pass keeps beyond the sums every maximisation reads.**
#[derive(Copy, Clone, Debug)]
pub(super) struct PassKeeps {
    /// Each position's probability of being mismapped and, when the run asks for them, each
    /// sample's genotype posteriors there — one value a position, kept only on the pass whose
    /// parameters are reported. A pass that keeps them joins its chunks in position order.
    pub(super) per_position_posteriors: bool,
    /// The information the standard errors are computed from ([`InformationSums`]): a score per
    /// position and parameter, summed into a few hundred bytes a sample. Kept only when a pass
    /// asks for it.
    pub(super) information: bool,
}

impl PassKeeps {
    /// The sums, and nothing else — what the tests' plain pass keeps.
    #[cfg(test)]
    pub(super) const SUMS_ONLY: Self = Self {
        per_position_posteriors: false,
        information: false,
    };
}

/// The same pass, told what to keep beyond its sums ([`PassKeeps`]).
///
/// Every pass joins its chunks' totals in an order fixed by the number of positions, so the
/// totals do not depend on the pool's width.
#[cfg(test)]
pub(super) fn expectation_pass(
    samples: &[SampleGenericSections<'_>],
    depth_cap: DepthCap,
    config: &JointFitConfig,
    group_index: &[Vec<usize>],
    coverage: &[f64],
    parameters: &Parameters,
    keeps: PassKeeps,
) -> Statistics {
    expectation_pass_every(
        samples,
        depth_cap,
        config,
        group_index,
        coverage,
        parameters,
        keeps,
        1,
    )
}

/// The same pass over only every `stride`-th position ([`EvidenceCursor::every`]); one is every
/// position.
#[allow(
    clippy::too_many_arguments,
    reason = "expectation_pass's arguments and the stride"
)]
pub(super) fn expectation_pass_every(
    samples: &[SampleGenericSections<'_>],
    depth_cap: DepthCap,
    config: &JointFitConfig,
    group_index: &[Vec<usize>],
    coverage: &[f64],
    parameters: &Parameters,
    keeps: PassKeeps,
    stride: usize,
) -> Statistics {
    let collect_noisy_posterior = keeps.per_position_posteriors;
    let model = PassModel::new(parameters, config, group_index);
    let quadrature = &model.quadrature;
    // What the scorer needs that depends on the parameters alone: once a pass.
    let scoring_tables = keeps.information.then(|| ScoringTables::of(&model));
    let positions = EvidenceCursor::position_count(samples);
    // **The chunks do not depend on the pool's width.** Every total below is a floating-point sum,
    // and a sum's last bits depend on where it is split and in what order the parts are added.
    // Chunks sized from the thread count, joined by `reduce` in whatever tree the pool chose, made
    // the same census fit to a different file at one, four and eight threads (measured on the
    // two-sample fixture of `cli::cross_platform_digests`). The test
    // `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` holds it.
    let chunk = positions
        .div_ceil(CHUNKS_PER_CENSUS)
        .max(MIN_POSITIONS_PER_CHUNK);
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
        if keeps.information {
            statistics.information = Some(InformationSums::for_a_cohort_of(group_index));
        }
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
            group_index,
            &config.edges,
            depth_cap,
            coverage,
            config.depth_as_a_range,
            first,
            end,
        )
        .every(stride);
        let mut scratch = Scratch::new(
            samples.len(),
            parameters.clean.len(),
            quadrature.nodes.len(),
        );
        let mut odds = vec![
            1.0_f64;
            if config.coverage_odds.is_empty() {
                0
            } else {
                samples.len()
            }
        ];
        // Whether this pass scores its positions, and the row it scores them into.
        let mut scoring = scoring_tables.as_ref().map(|tables| {
            (
                tables,
                PositionScores::new(group_index, parameters.clean.len()),
            )
        });
        let mut index = first.next_multiple_of(stride);
        while cursor.next_position(&mut scratch.evidence) {
            for (slot, sample) in odds.iter_mut().enumerate() {
                *sample = float::exp(f64::from(
                    config.coverage_odds[slot]
                        .get(index)
                        .copied()
                        .unwrap_or(0.0),
                ));
            }
            one_position(&mut scratch, &model, &odds, &mut statistics);
            if let Some((tables, scores)) = scoring.as_mut() {
                score_position(&scratch, &model, tables, &odds, scores);
                statistics
                    .information
                    .as_mut()
                    .expect("a pass that scores its positions keeps their sums")
                    .add_position(scores);
            }
            index += stride;
        }
        statistics
    };

    if bounds.is_empty() {
        return empty();
    }
    // **The chunks' totals are joined in an order that depends only on how many chunks there are.**
    // `reduce` may join them in any tree it likes, which is wrong for the per-position lists and,
    // for the sums, makes their last bits depend on the pool.
    if collect_noisy_posterior {
        // The reported pass: the chunks are added once each, in position order, into a total sized
        // for the whole census, so each per-position list is copied once.
        //
        // **A window of chunks at a time, not every chunk at once.** Collecting all of them before
        // adding any held every chunk's per-position lists — a second copy of the whole census's,
        // 8 MB at two million positions — and every chunk's sums, whose genotype counts and
        // information grow with the cohort: about a kilobyte a sample a chunk, so 0.2 to 0.3 GB of
        // chunks at 2,000 samples. The window only decides how many chunks are held at once; they
        // are still added in position order, so the totals are the same bits at any window and
        // any pool.
        use rayon::prelude::*;
        let mut into = empty();
        into.noisy_posterior.reserve(positions);
        if collect_genotype_posterior {
            into.genotype_posterior
                .reserve(positions * samples.len() * 3);
            into.duplicated_posterior.reserve(positions);
        }
        let window = 2 * rayon::current_num_threads().max(1);
        for some in bounds.chunks(window) {
            let chunks: Vec<Statistics> = some.par_iter().copied().map(one_chunk).collect();
            for chunk in &chunks {
                into.absorb(chunk);
            }
        }
        return into;
    }
    // The iterating passes, which keep sums only: the chunk list is halved at fixed points, the
    // halves are computed with `rayon::join`, and the right half's totals are added to the left's.
    // The same additions in the same order at any width, and one partial total held per level.
    fn joined<S, F>(bounds: &[(usize, usize)], one_chunk: &F, absorb: fn(&mut S, &S)) -> S
    where
        S: Send,
        F: Fn((usize, usize)) -> S + Sync,
    {
        if let [only] = bounds {
            return one_chunk(*only);
        }
        let (left, right) = bounds.split_at(bounds.len() / 2);
        let (mut into, from) = rayon::join(
            || joined(left, one_chunk, absorb),
            || joined(right, one_chunk, absorb),
        );
        absorb(&mut into, &from);
        into
    }
    joined(&bounds, &one_chunk, Statistics::absorb)
}

/// How many chunks a census is cut into. **The chunk size comes from the census's length and
/// these two constants, never from the pool's width**, so a census is cut at the same positions on
/// every machine.
///
/// **What sets the number: the wait at the end of every pass.** A pass ends when its last chunk
/// does, and the chunks cost the same to within about 2%, so with too few of them most threads
/// finish their last chunk and wait while a few run one more. At 128 chunks on 18 threads (7.1
/// chunks a thread) the workers idled 6.6% of a 59-sample fit; on four samples, 5.8 of the 6.4
/// idle points were that tail. 512 chunks are a quarter the size, so the tail is about a quarter
/// as long, and took the 59-sample fit (18 threads, two runs each, interleaved) from 880.9 and
/// 866.9 s to 827.8 and 829.6 s with the same instructions retired
/// (`doc/devel/reports/reviews/perf_snp-fit_2026-10-18.md`, L3). That leaves about 28 chunks a
/// thread at 18 threads and 4 at 128. Each chunk's setup — a scratch, a set of totals, one binary
/// search into each sample's sparse list — stays small against its work at 512: about 3,900
/// positions a chunk at two million positions; a census under 131,072 positions gets fewer chunks
/// than 512, because no chunk is shorter than `MIN_POSITIONS_PER_CHUNK`. The iterating passes'
/// join tree is two levels deeper than at 128, so each thread holds up to two more partial totals.
/// Far more chunks are too many: 7,813 chunks of 256 positions made the four-sample fit 28% slower
/// than 123 chunks of 16,384 (one run each, not interleaved).
const CHUNKS_PER_CENSUS: usize = 512;

/// The smallest chunk, so a census of a few thousand positions is not cut into chunks whose setup
/// outweighs their work. Measured on the joint-fit bench at four threads: 256-position chunks
/// matched the thread-sized chunks they replaced to within about 2% on every SNP/indel case.
pub(super) const MIN_POSITIONS_PER_CHUNK: usize = 256;

/// At most three non-reference bases can be the segregating one, so the candidate list never
/// exceeds three however many the samples showed reads on.
pub(super) const MAX_CANDIDATES: usize = CANDIDATE_ALTERNATIVES;

/// What a position can be, within a noise class: the population carries one allele, it is fixed
/// on a non-reference one, it segregates, or it is a stretch some samples carry twice.
const BRANCHES: usize = 4;
const DUPLICATED: usize = 3;

/// Everything one position needs, sized once and reused, so a two-million-position pass
/// allocates nothing.
pub(super) struct Scratch {
    pub(super) evidence: PositionEvidence,
    pub(super) samples: usize,
    pub(super) nodes: usize,
    /// `[class][candidate][sample][genotype]` — `ln P(reads | genotype)`, and the same thing
    /// exponentiated after its own per-sample maximum is taken out.
    ///
    /// **This is the whole reason the program can be run.** The read likelihoods depend on the
    /// genotype and the error rate but **not on the allele frequency**, so they are computed
    /// once per candidate and reused across every quadrature node — which is the fix the
    /// previous session had to make to its repeat-tract program after it turned out to take
    /// hours (`joint_route_research_narrative_2026-08-13.md` §7).
    ell: Vec<f64>,
    pub(super) lik: Vec<f64>,
    /// `[class][candidate][sample]` — the maximum taken out of `lik`.
    lik_max: Vec<f64>,
    /// `[class][candidate][node]`, and `[class][candidate]`.
    pub(super) node_ln: Vec<f64>,
    /// The same, over the carrier frequency the duplicated class integrates.
    pub(super) carrier_node_ln: Vec<f64>,
    pub(super) fixed_ln: Vec<f64>,
    /// `[class]` — the invariant branch, and the four branches combined.
    pub(super) invariant_ln: Vec<f64>,
    branch_ln: Vec<f64>,
    class_ln: Vec<f64>,
    /// `[class]` — the log-likelihood of the branch where every sample is fixed for a
    /// non-reference base, before its share multiplies it. **Kept for the scores**
    /// ([`information`](super::information)), which differentiate the position's likelihood in
    /// the shares and need each branch apart from the share it is weighted by; likewise the next
    /// two fields.
    pub(super) fixed_alt_ln: [f64; 2],
    /// `[class]` — the integral over a segregating frequency, before its share multiplies it.
    pub(super) segregating_ln: [f64; 2],
    /// `[class]` — the integral over how much of the panel carries an extra copy, before the
    /// duplicated share multiplies it; `−∞` when the class is not fitted.
    pub(super) duplicated_ln: [f64; 2],
    /// `[class]` — each branch's log-sum over its candidates (and nodes) before `ln 3` is taken
    /// off: **the totals the attribution divides each share by**, kept so it does not sum the
    /// same terms a second time. Likewise the next two fields.
    fixed_sum_ln: [f64; 2],
    segregating_sum_ln: [f64; 2],
    /// `−∞` when the duplicated class is not fitted.
    duplicated_sum_ln: [f64; 2],
    /// `[class]` — the four branches' log-sum, which the attribution and the scorer divide by.
    pub(super) within_ln: [f64; 2],
    /// The position's whole log-likelihood, once [`one_position`] has computed it.
    pub(super) position_ln: f64,
    /// The candidate alleles this position sums over, and the log of how many alleles each stands
    /// for — the only form the multiplicity is ever read in.
    pub(super) candidates: Vec<usize>,
    pub(super) ln_multiplicity: Vec<f64>,
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
    /// Room for `samples` samples holding `libraries` read groups between them, at a rule of
    /// `nodes` nodes.
    pub(super) fn new(samples: usize, libraries: usize, nodes: usize) -> Self {
        Self {
            evidence: PositionEvidence {
                libraries: vec![LibraryAtPosition::default(); libraries],
                depth_weights: vec![0.0; libraries * MAX_RECORDED_SPREAD],
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
            fixed_alt_ln: [f64::NEG_INFINITY; 2],
            segregating_ln: [f64::NEG_INFINITY; 2],
            duplicated_ln: [f64::NEG_INFINITY; 2],
            fixed_sum_ln: [f64::NEG_INFINITY; 2],
            segregating_sum_ln: [f64::NEG_INFINITY; 2],
            duplicated_sum_ln: [f64::NEG_INFINITY; 2],
            within_ln: [f64::NEG_INFINITY; 2],
            position_ln: f64::NEG_INFINITY,
            candidates: Vec::with_capacity(MAX_CANDIDATES),
            ln_multiplicity: Vec::with_capacity(MAX_CANDIDATES),
            genotype_weight: vec![0.0; samples * 3],
            carrier_weight: vec![0.0; samples],
            position_genotype: vec![0.0; samples * 3],
            per_candidate_weight: vec![0.0; MAX_CANDIDATES * samples * 3],
            shares: Vec::with_capacity(MAX_CANDIDATES * nodes),
        }
    }

    pub(super) fn ell_at(&self, class: usize, candidate: usize, sample: usize) -> usize {
        ((class * MAX_CANDIDATES + candidate) * self.samples + sample) * 3
    }

    fn max_at(&self, class: usize, candidate: usize, sample: usize) -> usize {
        (class * MAX_CANDIDATES + candidate) * self.samples + sample
    }

    pub(super) fn node_at(&self, class: usize, candidate: usize, node: usize) -> usize {
        (class * MAX_CANDIDATES + candidate) * self.nodes + node
    }
}

/// The scale a running product is multiplied back up by when it is about to underflow, and its
/// logarithm. Sixty-three samples' likelihoods multiplied together will underflow a `f64`
/// several times over; rescaling keeps the product in range with **one logarithm per node
/// instead of one per sample**, which is a sixty-fold saving in the innermost loop.
const RESCALE: f64 = 1e150;
const LN_RESCALE: f64 = 345.398_899_014_487; // ln(1e150)

/// **`ln Π terms`, the product kept in range by [`RESCALE`]**: multiplied back up whenever it is
/// about to underflow, and its logarithm taken once at the end. `None` where the product is zero
/// or negative — a node no sample's reads allow.
#[inline(always)]
fn ln_of_product(terms: impl Iterator<Item = f64>) -> Option<f64> {
    let mut product = 1.0_f64;
    let mut scale = 0.0_f64;
    for term in terms {
        product *= term;
        if product < 1.0 / RESCALE {
            product *= RESCALE;
            scale -= LN_RESCALE;
        }
    }
    if product <= 0.0 {
        None
    } else {
        Some(float::ln(product) + scale)
    }
}

/// One position: its likelihood, and its contribution to every accumulated count.
#[allow(
    clippy::needless_range_loop,
    reason = "the sample index addresses four parallel arrays at once — the evidence, the \
              scratch, the read-group map and the accumulated counts — and zipping them would \
              hide which of the four the loop is really walking"
)]
pub(super) fn one_position(
    scratch: &mut Scratch,
    model: &PassModel<'_>,
    coverage_odds: &[f64],
    statistics: &mut Statistics,
) {
    let group_index = model.group_index;
    let ploidy = model.ploidy;
    let read_logs = model.read_logs.as_slice();
    let quadrature = &model.quadrature;
    let carrier = model.carrier.as_ref().map(|rule| &rule.quadrature);
    let parameters = model.parameters;
    let shares = &model.branches;
    let samples = scratch.samples;
    let nodes = scratch.nodes;
    statistics.positions += 1.0;
    for s in 0..samples {
        let depth: f64 = group_index[s]
            .iter()
            .map(|&library| scratch.evidence.libraries[library].depth)
            .sum();
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
    scratch.ln_multiplicity.clear();
    for &code in &scratch.evidence.observed_alternatives {
        scratch.candidates.push(code);
        scratch.ln_multiplicity.push(shares.ln_multiplicity[1]);
    }
    let unobserved = CANDIDATE_ALTERNATIVES.saturating_sub(scratch.candidates.len());
    if unobserved > 0 {
        let stand_in = (0..4)
            .find(|code| !scratch.candidates.contains(code))
            .expect("with an unobserved allele left there is a base no read fell on");
        scratch.candidates.push(stand_in);
        scratch
            .ln_multiplicity
            .push(shares.ln_multiplicity[unobserved]);
    }
    let candidates = scratch.candidates.len();

    // ---- the read likelihoods, once per candidate and reused across every node -------------
    //
    // **The library loop is outside the candidate loop so that the two candidate-invariant
    // quantities are taken once.** [`reference_terms`] holds this block's only logarithms and
    // `non_reference()` is a sum over four counts; neither reads which base the alternative
    // is, so with the candidate loop outermost both were recomputed for every candidate — ten
    // reference terms a library where three do the same job. `invariant` and each candidate's
    // `fixed` accumulate over the samples in the samples' order.
    //
    // **A sample's reads under a genotype are its libraries' reads under that genotype, each
    // library scored under its own two error rates**, and added: the genotype is one a sample,
    // the rates are one a library. The first library's terms are taken as they are and each
    // later library's added to them, so a sample of one library is scored with exactly the
    // arithmetic it was scored with before its libraries were kept apart.
    for class in 0..2 {
        let mut invariant = 0.0;
        let mut fixed = [0.0_f64; MAX_CANDIDATES];
        for s in 0..samples {
            let mut sample_invariant = 0.0;
            if group_index[s].is_empty() {
                // A sample with no library has no reads: every genotype explains them alike.
                for candidate in 0..candidates {
                    let base = scratch.ell_at(class, candidate, s);
                    scratch.ell[base..base + 3].fill(0.0);
                }
            }
            for (section, &library) in group_index[s].iter().enumerate() {
                let logs = &read_logs[class][library];
                let reads = &scratch.evidence.libraries[library];
                let weights = scratch.evidence.weights_of(library);
                let non_reference = reads.non_reference();
                let reference_term = reference_terms(reads, weights, logs);
                let all_reference =
                    ln_reads_given_all_reference(reads, non_reference, logs, reference_term[0]);
                sample_invariant = if section == 0 {
                    all_reference
                } else {
                    sample_invariant + all_reference
                };
                for candidate in 0..candidates {
                    let allele = scratch.candidates[candidate];
                    let base = scratch.ell_at(class, candidate, s);
                    for j in 0..3 {
                        let value = ln_reads_given_genotype(
                            reads,
                            non_reference,
                            j as u8,
                            allele,
                            logs,
                            reference_term[j],
                        );
                        scratch.ell[base + j] = if section == 0 {
                            value
                        } else {
                            scratch.ell[base + j] + value
                        };
                    }
                }
            }
            invariant += sample_invariant;
            for candidate in 0..candidates {
                let base = scratch.ell_at(class, candidate, s);
                let largest = scratch.ell[base..base + 3]
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                fixed[candidate] += scratch.ell[base + 2];
                let slot = scratch.max_at(class, candidate, s);
                scratch.lik_max[slot] = largest;
                for j in 0..3 {
                    scratch.lik[base + j] = float::exp(scratch.ell[base + j] - largest);
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
                let product_ln = ln_of_product((0..samples).map(|s| {
                    let prior = &quadrature.priors[(node * samples + s) * 3..][..3];
                    let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                    prior[0] * lik[0] + prior[1] * lik[1] + prior[2] * lik[2]
                }));
                let slot = scratch.node_at(class, candidate, node);
                scratch.node_ln[slot] = match product_ln {
                    Some(product_ln) => product_ln + offset + quadrature.ln_weights[node],
                    None => f64::NEG_INFINITY,
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
                    let product_ln = ln_of_product((0..samples).map(|s| {
                        let lik = &scratch.lik[scratch.ell_at(class, candidate, s)..][..3];
                        // Where the run has a coverage summary, this is the only place it
                        // enters: how much more likely this sample's read depth around here is
                        // if it has two copies rather than one. Every other branch has every
                        // sample single-copy, so the one-copy factor is common to all of them
                        // and cancels out of the position entirely.
                        let odds = coverage_odds.get(s).copied().unwrap_or(1.0);
                        (1.0 - q) * lik[0] + q * odds * lik[1]
                    }));
                    let slot = scratch.node_at(class, candidate, node);
                    scratch.carrier_node_ln[slot] = match product_ln {
                        Some(product_ln) => product_ln + offset + carrier.ln_weights[node],
                        None => f64::NEG_INFINITY,
                    };
                }
            }
        }
    }

    // ---- the four branches, and the two classes ---------------------------------------------
    let ln_three = shares.ln_three;
    for class in 0..2 {
        scratch.shares.clear();
        for candidate in 0..candidates {
            scratch.shares.push(
                scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                    + scratch.ln_multiplicity[candidate],
            );
        }
        let fixed_sum = float::ln_sum_exp(&scratch.shares);
        scratch.fixed_sum_ln[class] = fixed_sum;
        let fixed_alt = fixed_sum - ln_three;
        scratch.shares.clear();
        for candidate in 0..candidates {
            let multiplicity = scratch.ln_multiplicity[candidate];
            for node in 0..nodes {
                scratch
                    .shares
                    .push(scratch.node_ln[scratch.node_at(class, candidate, node)] + multiplicity);
            }
        }
        let segregating_sum = float::ln_sum_exp(&scratch.shares);
        scratch.segregating_sum_ln[class] = segregating_sum;
        let segregating = segregating_sum - ln_three;
        scratch.fixed_alt_ln[class] = fixed_alt;
        scratch.segregating_ln[class] = segregating;
        // The three ordinary branches share what is left when the duplicated class has taken
        // its share, so the four still sum to one.
        let ordinary = shares.ln_ordinary;
        scratch.branch_ln[class * BRANCHES] =
            ordinary + shares.ln_invariant + scratch.invariant_ln[class];
        scratch.branch_ln[class * BRANCHES + 1] = ordinary + shares.ln_fixed_alt + fixed_alt;
        scratch.branch_ln[class * BRANCHES + 2] = ordinary + shares.ln_segregating + segregating;
        scratch.branch_ln[class * BRANCHES + DUPLICATED] = if let Some(carrier) = carrier {
            scratch.shares.clear();
            for candidate in 0..candidates {
                let multiplicity = scratch.ln_multiplicity[candidate];
                for node in 0..carrier.nodes.len() {
                    scratch.shares.push(
                        scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                            + multiplicity,
                    );
                }
            }
            let summed = float::ln_sum_exp(&scratch.shares);
            scratch.duplicated_sum_ln[class] = summed;
            // `duplicated_ln` is kept for the scores and not reused in the branch's own sum:
            // `ln d + summed − ln 3` is added left to right, and taking `summed − ln 3` first
            // would move the branch's last bits.
            scratch.duplicated_ln[class] = summed - ln_three;
            shares.ln_duplicated + summed - ln_three
        } else {
            scratch.duplicated_sum_ln[class] = f64::NEG_INFINITY;
            scratch.duplicated_ln[class] = f64::NEG_INFINITY;
            f64::NEG_INFINITY
        };
        let within = float::ln_sum_exp(&scratch.branch_ln[class * BRANCHES..][..BRANCHES]);
        scratch.within_ln[class] = within;
        scratch.class_ln[class] = shares.ln_class[class] + within;
    }
    let position_ln = float::ln_sum_exp(&scratch.class_ln);
    scratch.position_ln = position_ln;
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
        // **And it counts as a position for the two population moments, carrying nothing.** Both
        // routes to those moments have to see the same census: the stored array records this
        // position as all-zero posteriors, so the running sums must record a position of no
        // alternative copies rather than skip it, or the two would divide by different counts.
        //
        // **⚠ No test enters this branch, and none of the fixtures can.** Every logarithm on the
        // way here is taken through `float::ln(p.max(f64::MIN_POSITIVE))`, which floors a read's
        // term at about −708, and `ln_sum_exp` returns `−∞` only when every input already is one —
        // so reaching `!position_ln.is_finite()` from well-formed evidence would take on the order
        // of 1e305 reads at a position. What the two lines below maintain is tested where it can be:
        // an all-zero position counts and contributes nothing, in `census_moments`'s
        // all-reference fixtures, which build the identical buffer. Replacing this block with a
        // `panic!` leaves the suite green, and that is a statement about the branch's
        // reachability, not about the bookkeeping.
        scratch.position_genotype.fill(0.0);
        statistics
            .census_moments
            .add_position(&scratch.position_genotype);
        return;
    }
    scratch.position_genotype.fill(0.0);
    let mut duplicated_here = 0.0_f64;
    statistics.log_likelihood += position_ln;
    if statistics.collect_noisy_posterior {
        statistics
            .noisy_posterior
            .push(float::exp(scratch.class_ln[1] - position_ln) as f32);
    }

    // ---- attribute it -------------------------------------------------------------------------
    for class in 0..2 {
        let class_posterior = float::exp(scratch.class_ln[class] - position_ln);
        if class_posterior <= 1e-12 {
            continue;
        }
        if class == 1 {
            statistics.noisy += class_posterior;
        }
        let branches = &scratch.branch_ln[class * BRANCHES..][..BRANCHES];
        let within = scratch.within_ln[class];
        let branch = [
            class_posterior * float::exp(branches[0] - within),
            class_posterior * float::exp(branches[1] - within),
            class_posterior * float::exp(branches[2] - within),
            class_posterior * float::exp(branches[DUPLICATED] - within),
        ];
        statistics.invariant += branch[0];
        statistics.fixed_alt += branch[1];
        statistics.segregating += branch[2];
        statistics.duplicated += branch[DUPLICATED];
        duplicated_here += branch[DUPLICATED];

        // The invariant branch: every sample is homozygous reference, and every read that is
        // not on the reference base is an error, so none of them is "the allele".
        //
        // **Here and in every branch below, a library's tally is credited with that library's
        // reads alone, at its own rate** — so each library's rates are fitted from its own
        // evidence, where crediting every library of a sample with the sample's reads gave them
        // one rate between them.
        if branch[0] > 1e-12 {
            for &library in group_index.iter().flatten() {
                let reads = &scratch.evidence.libraries[library];
                let rate = class_rate(parameters, class, library);
                let neither = reads.non_reference() - reads.on[4];
                let reference = expected_reference_reads(
                    reads,
                    scratch.evidence.weights_of(library),
                    1.0 - rate,
                );
                let tally = &mut statistics.reads[library][class][0];
                tally.neither += branch[0] * neither;
                tally.reference += branch[0] * reference;
            }
        }

        // The fixed-alternative branch: every sample carries two copies of the candidate.
        if branch[1] > 1e-12 {
            let branch_sum_ln = scratch.fixed_sum_ln[class];
            for candidate in 0..candidates {
                let share = branch[1]
                    * float::exp(
                        scratch.fixed_ln[class * MAX_CANDIDATES + candidate]
                            + scratch.ln_multiplicity[candidate]
                            - branch_sum_ln,
                    );
                if share <= 1e-12 {
                    continue;
                }
                let allele = scratch.candidates[candidate];
                for s in 0..samples {
                    for &library in &group_index[s] {
                        let rate = class_rate(parameters, class, library);
                        let counts = tally_of(
                            &scratch.evidence.libraries[library],
                            scratch.evidence.weights_of(library),
                            allele,
                            reference_read_probability(2, ploidy, rate),
                        );
                        statistics.reads[library][class][2].credit(share, counts);
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
            let branch_sum_ln = scratch.segregating_sum_ln[class];
            scratch.per_candidate_weight[..candidates * samples * 3].fill(0.0);
            scratch.genotype_weight.fill(0.0);
            for candidate in 0..candidates {
                let multiplicity = scratch.ln_multiplicity[candidate];
                for node in 0..nodes {
                    let share = branch[2]
                        * float::exp(
                            scratch.node_ln[scratch.node_at(class, candidate, node)] + multiplicity
                                - branch_sum_ln,
                        );
                    if share <= 1e-12 {
                        continue;
                    }
                    #[cfg(test)]
                    {
                        statistics.sum_ln_f += share * quadrature.ln_nodes[node];
                    }
                    // How this node's product over samples moves as its frequency moves: each
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
                        let prior_slope = &quadrature.prior_slopes[(node * samples + s) * 3..][..3];
                        along_frequency += (prior_slope[0] * lik[0]
                            + prior_slope[1] * lik[1]
                            + prior_slope[2] * lik[2])
                            / total;
                        let base = (candidate * samples + s) * 3;
                        for j in 0..3 {
                            let weight = share * joint[j] / total;
                            scratch.per_candidate_weight[base + j] += weight;
                            scratch.genotype_weight[s * 3 + j] += weight;
                            statistics.genotypes[s][node][j] += weight;
                        }
                    }
                    // The node's share of the likelihood's slope in each shape: its weight moves
                    // and it moves along the frequency axis ([`information::RuleSlopes`]).
                    let rule = &model.density_slopes;
                    for (shape, slope) in statistics.density_shape_slopes.iter_mut().enumerate() {
                        *slope += share
                            * (rule.ln_weight_slope[shape][node]
                                + along_frequency * rule.node_slope[shape][node]);
                    }
                }
            }
            for candidate in 0..candidates {
                let allele = scratch.candidates[candidate];
                for s in 0..samples {
                    let base = (candidate * samples + s) * 3;
                    for j in 0..3 {
                        let weight = scratch.per_candidate_weight[base + j];
                        if weight <= 0.0 {
                            continue;
                        }
                        for &library in &group_index[s] {
                            let rate = class_rate(parameters, class, library);
                            let counts = tally_of(
                                &scratch.evidence.libraries[library],
                                scratch.evidence.weights_of(library),
                                allele,
                                reference_read_probability(j as u8, ploidy, rate),
                            );
                            statistics.reads[library][class][j].credit(weight, counts);
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
        if let Some(CarrierRule {
            quadrature: carrier,
            slopes: rule,
        }) = model
            .carrier
            .as_ref()
            .filter(|_| branch[DUPLICATED] > 1e-12)
        {
            {
                let carrier_nodes = carrier.nodes.len();
                let branch_sum_ln = scratch.duplicated_sum_ln[class];
                scratch.per_candidate_weight[..candidates * samples * 3].fill(0.0);
                scratch.carrier_weight.fill(0.0);
                for candidate in 0..candidates {
                    let multiplicity = scratch.ln_multiplicity[candidate];
                    for node in 0..carrier_nodes {
                        let share = branch[DUPLICATED]
                            * float::exp(
                                scratch.carrier_node_ln[scratch.node_at(class, candidate, node)]
                                    + multiplicity
                                    - branch_sum_ln,
                            );
                        if share <= 1e-12 {
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
                            let base = (candidate * samples + s) * 3;
                            scratch.per_candidate_weight[base] += share * joint[0] / total;
                            scratch.per_candidate_weight[base + 1] += share * joint[1] / total;
                            scratch.carrier_weight[s] += share * joint[1] / total;
                        }
                        for (shape, slope) in statistics.carrier_shape_slopes.iter_mut().enumerate()
                        {
                            *slope += share
                                * (rule.ln_weight_slope[shape][node]
                                    + along_frequency * rule.node_slope[shape][node]);
                        }
                    }
                }
                for candidate in 0..candidates {
                    let allele = scratch.candidates[candidate];
                    for s in 0..samples {
                        let base = (candidate * samples + s) * 3;
                        for j in 0..2 {
                            let weight = scratch.per_candidate_weight[base + j];
                            if weight <= 0.0 {
                                continue;
                            }
                            for &library in &group_index[s] {
                                let rate = class_rate(parameters, class, library);
                                let counts = tally_of(
                                    &scratch.evidence.libraries[library],
                                    scratch.evidence.weights_of(library),
                                    allele,
                                    reference_read_probability(j as u8, ploidy, rate),
                                );
                                statistics.reads[library][class][j].credit(weight, counts);
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

    // **The two population moments, summed here rather than stored and reduced afterwards.**
    // `scratch.position_genotype` is the buffer this position just filled; nothing is copied and
    // nothing is kept (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5).
    statistics
        .census_moments
        .add_position(&scratch.position_genotype);
    if statistics.collect_genotype_posterior {
        statistics
            .genotype_posterior
            .extend(scratch.position_genotype.iter().map(|v| *v as f32));
        statistics.duplicated_posterior.push(duplicated_here as f32);
    }
}

pub(super) fn class_rate(parameters: &Parameters, class: usize, group: usize) -> f64 {
    if class == 0 {
        parameters.clean[group]
    } else {
        parameters.noisy[group]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
