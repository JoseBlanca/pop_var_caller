//! **One sample's walk along a window** — the copying model with an autozygous state
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §3).
//!
//! At each record the hidden state is either *the two copies copy templates `k₁` and `k₂`* (K²
//! states) or *the two copies are one, copying template `k`* (K states). A forward and a backward
//! walk give each state's probability given the whole window; the genotype posterior is the
//! sample's own likelihood weighted by what the states say the alleles are (§3.6).
//!
//! Diploid only: the pair state is two copies (§5). The arithmetic is the prototype's
//! (`doc/devel/ng/research/neighbour_evidence/hmm.py`, `haplotypes`), checked against it to 1e-9 by
//! this module's tests.

use super::random::{Uniforms, pick, pick_at, pick_in_running};
use crate::calling::genotype_table::{GenotypeIdx, GenotypeTable};
use crate::float;
use crate::types::Ploidy;
use std::collections::HashMap;

/// **The four numbers a walk runs with** (spec §3.4, §3.5).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WalkParameters {
    /// Mean length, in bases, of a stretch one copy copies from one template (`S`).
    pub segment_bases: f64,
    /// Mean length, in bases, of an autozygous run (`R`).
    pub run_bases: f64,
    /// Long-run share of the sample that is autozygous (`F`), strictly between 0 and 1.
    pub autozygous_share: f64,
    /// The chance a template's allele is not the one the copy carries (`ε`).
    pub mismatch: f64,
}

/// **One record as one sample's walk sees it.**
#[derive(Clone, Copy, Debug)]
pub struct WalkSite<'a> {
    /// The record's position on the contig; sites must ascend.
    pub position: u64,
    /// How many alleles the record writes.
    pub allele_count: usize,
    /// The sample's likelihood of each genotype, in the genotype table's order, on any common
    /// scale — read off `PL` (§3.2).
    pub likelihoods: &'a [f64],
    /// Each template's allele shares, template by template: `templates × allele_count` numbers,
    /// each template's summing to one — read off its `GP` (§3.2).
    pub template_shares: &'a [f64],
}

/// **What a walk gives back**, one entry per site.
#[derive(Clone, Debug, PartialEq)]
pub struct WalkResult {
    /// The sample's genotype posterior at each site, in the genotype table's order.
    pub posteriors: Vec<Vec<f64>>,
    /// **What the sample's neighbours alone say at each site**: the genotype distribution the
    /// states give before the sample's own likelihood is multiplied in (§3.6), normalised. At a
    /// site where the sample has no reads it is the posterior; where it has reads, comparing the
    /// two measures how well the neighbours predict this sample.
    pub neighbour_predictions: Vec<Vec<f64>>,
    /// The probability the sample is autozygous at each site.
    pub autozygous: Vec<f64>,
    /// **The probability that a copy copies the last template** at each site — at least one of the
    /// two, or the one an autozygous state copies. The window solve puts the background template
    /// last (§3.3), so this is the chance the sample has no sampled haplotype to copy there.
    pub last_template_copy: Vec<f64>,
    /// The log-likelihood of the sample's evidence over the window under these parameters — what
    /// the parameter grid is chosen by (§3.5).
    pub log_likelihood: f64,
}

/// The two copies' alleles behind each genotype index, read off the shared genotype table — never
/// an index formula of our own (spec §7) — and kept per allele count.
#[derive(Default)]
pub(crate) struct GenotypePairs {
    by_allele_count: HashMap<usize, Vec<(usize, usize)>>,
}

impl GenotypePairs {
    pub(crate) fn of(&mut self, allele_count: usize) -> &[(usize, usize)] {
        self.by_allele_count.entry(allele_count).or_insert_with(|| {
            let diploid = Ploidy::try_new(2).expect("two copies is a ploidy");
            let table = GenotypeTable::build(diploid, allele_count);
            let table = table.view();
            (0..table.genotype_count())
                .map(|genotype| {
                    // PANIC-FREE: the index ranges over the table's own count.
                    let copies = table
                        .allele_counts_of(GenotypeIdx(genotype as u32))
                        .expect("an index below the count names a genotype");
                    let mut alleles = copies
                        .iter()
                        .enumerate()
                        .flat_map(|(allele, &count)| std::iter::repeat_n(allele, count as usize));
                    let first = alleles.next().expect("a diploid genotype has a first copy");
                    let second = alleles
                        .next()
                        .expect("a diploid genotype has a second copy");
                    (first, second)
                })
                .collect()
        })
    }
}

/// **One site's emissions and blurred shares**, kept from the forward walk for the backward one.
struct SiteTerms {
    /// `q̃`: `templates × allele_count`, each template's shares blurred by the mismatch rate.
    blurred: Vec<f64>,
    /// The emission of each autozygous state, `K`.
    autozygous: Vec<f64>,
    /// The emission of each pair state, `K × K`, row `k₁`.
    pairs: Vec<f64>,
}

fn site_terms(
    site: &WalkSite<'_>,
    templates: usize,
    mismatch: f64,
    pairs: &[(usize, usize)],
) -> SiteTerms {
    let alleles = site.allele_count;
    assert_eq!(
        site.likelihoods.len(),
        pairs.len(),
        "a site at {} has {} alleles, so {} diploid genotypes, and {} likelihoods",
        site.position,
        alleles,
        pairs.len(),
        site.likelihoods.len()
    );
    assert_eq!(
        site.template_shares.len(),
        templates * alleles,
        "a site at {} needs {templates} templates × {alleles} allele shares",
        site.position
    );
    let blurred: Vec<f64> = site
        .template_shares
        .iter()
        .map(|share| (1.0 - mismatch) * share + mismatch / alleles as f64)
        .collect();
    // The likelihood as a symmetric allele-by-allele matrix, so the pair emission is two matrix
    // products: Σ_a Σ_b q̃_k₁(a) · L(a/b) · q̃_k₂(b).
    let mut by_alleles = vec![0.0; alleles * alleles];
    for (&(first, second), &likelihood) in pairs.iter().zip(site.likelihoods) {
        by_alleles[first * alleles + second] = likelihood;
        by_alleles[second * alleles + first] = likelihood;
    }
    let mut weighted = vec![0.0; templates * alleles];
    for k in 0..templates {
        for b in 0..alleles {
            weighted[k * alleles + b] = (0..alleles)
                .map(|a| blurred[k * alleles + a] * by_alleles[a * alleles + b])
                .sum();
        }
    }
    let mut pair_emissions = vec![0.0; templates * templates];
    for first in 0..templates {
        for second in 0..templates {
            pair_emissions[first * templates + second] = (0..alleles)
                .map(|b| weighted[first * alleles + b] * blurred[second * alleles + b])
                .sum();
        }
    }
    let autozygous = (0..templates)
        .map(|k| {
            (0..alleles)
                .map(|a| blurred[k * alleles + a] * by_alleles[a * alleles + a])
                .sum()
        })
        .collect();
    SiteTerms {
        blurred,
        autozygous,
        pairs: pair_emissions,
    }
}

/// The three persistence probabilities between two sites `gap` bases apart: a copy keeps its
/// template, the autozygous state persists, the other state persists.
///
/// **The two states are one two-state Markov chain** (spec §3.5): runs of mean `R`, gaps of mean
/// `R′ = R · (1 − F) / F`, so the chance of staying is `F + (1 − F) · e^(−d/τ)` in a run and
/// `(1 − F) + F · e^(−d/τ)` outside one, with `1/τ = 1/R + 1/R′`. Across a gap far longer than
/// both the state returns to its long-run share `F` — not, as two independent `e^(−d/R)` stays
/// would make it, to the opposite state (found by review, 2026-11-01).
fn persistence(gap: f64, parameters: &WalkParameters) -> (f64, f64, f64) {
    let share = parameters.autozygous_share;
    let gap_bases = parameters.run_bases * (1.0 - share) / share;
    let decay = float::exp(-gap * (1.0 / parameters.run_bases + 1.0 / gap_bases));
    (
        float::exp(-gap / parameters.segment_bases),
        share + (1.0 - share) * decay,
        (1.0 - share) + share * decay,
    )
}

/// **One sample's walk over a window's sites** (spec §3). `templates` is `K`; every site carries
/// that many templates' shares. Returns each site's genotype posterior and autozygous probability.
///
/// # Panics
///
/// If a site's likelihoods or shares do not match its allele count, or the autozygous share is
/// not strictly between 0 and 1.
#[must_use]
pub fn walk(sites: &[WalkSite<'_>], templates: usize, parameters: WalkParameters) -> WalkResult {
    walk_with_draw(sites, templates, parameters, None, None, None).0
}

/// **A sample's read link between two walk sites**: `L` between `one` and `other`, each with its
/// heterozygote's two alleles, the lower first; `None` when no fragment of the sample counts at both.
pub(crate) type ReadLinkOf<'a> =
    dyn Fn(usize, (usize, usize), usize, (usize, usize)) -> Option<f64> + 'a;

/// **What the reads give the haplotype draw** (`neighbour_evidence_reads_in_draw.md` §2.2): one
/// sample's read link between two walk sites, and the stream the orientations are drawn from.
pub(crate) struct DrawReads<'a> {
    /// The sample's read links.
    pub link: &'a ReadLinkOf<'a>,
    /// Whether the sample's reads at a walk site are kept out of the draw: one of its hidden
    /// records (§5) or of its records held out for the read check (§2.6). Such a site never joins a
    /// chain.
    pub kept_out: &'a dyn Fn(usize) -> bool,
    /// `HET_MASS_MIN`: the smallest heterozygous mass of an orientable record.
    pub het_mass_min: f64,
    /// The orientations' stream, read only where a read-linked group exists.
    pub uniforms: Uniforms,
}

/// **One haplotype pair at one site**: each copy's probability of carrying each allele, `A`
/// numbers a copy.
pub(crate) type HaplotypePair = [Vec<f64>; 2];

/// **The walk, and one state path drawn from it with the two haplotypes along it**
/// (`neighbour_evidence_phasing.md` §2.1). The walk is [`walk`]'s, unchanged. The path is drawn
/// backwards from the last site, each state given the one after it and the forward probabilities,
/// one uniform a site from `uniforms`, the states in the order autozygous `0..K` then pairs row by
/// row. Along the path, copy one carries allele `a` and copy two `b` with probability
/// `∝ q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})` in a pair state and `∝ q̃_k(a) · L({a, a})`, `a = b`, in
/// an autozygous one; each copy's haplotype is that joint's margin.
///
/// # Panics
///
/// As [`walk`]; and if every state at a site has zero weight — every genotype's likelihood zero
/// there, which a likelihood row normalised to its largest value never is.
#[must_use]
pub(crate) fn walk_and_draw(
    sites: &[WalkSite<'_>],
    templates: usize,
    parameters: WalkParameters,
    uniforms: &mut Uniforms,
    reads: Option<&mut DrawReads<'_>>,
) -> (WalkResult, Vec<HaplotypePair>) {
    let (result, drawn, _) =
        walk_with_draw(sites, templates, parameters, Some(uniforms), None, reads);
    (
        result,
        // PANIC-FREE: `walk_with_draw` draws whenever it is given a stream, as it is here.
        drawn.expect("a stream was given, so a path was drawn"),
    )
}

/// **[`walk_and_draw`], and `extra` more paths drawn from the same forward probabilities**
/// (`neighbour_evidence_phase_output.md` §2.3): the round's own path from `uniforms`, exactly as
/// [`walk_and_draw`] draws it, and the extra paths one after another from `extra_uniforms`, so
/// drawing them changes neither the round's path nor anything the round passes on. The extra paths
/// serve the phase's link test only.
///
/// # Panics
///
/// As [`walk_and_draw`].
#[must_use]
pub(crate) fn walk_and_draw_paths(
    sites: &[WalkSite<'_>],
    templates: usize,
    parameters: WalkParameters,
    uniforms: &mut Uniforms,
    extra_uniforms: &mut Uniforms,
    extra: usize,
    reads: Option<&mut DrawReads<'_>>,
) -> (WalkResult, Vec<HaplotypePair>, Vec<Vec<HaplotypePair>>) {
    let (result, drawn, extra_paths) = walk_with_draw(
        sites,
        templates,
        parameters,
        Some(uniforms),
        Some((extra_uniforms, extra)),
        reads,
    );
    (
        result,
        // PANIC-FREE: as in `walk_and_draw`.
        drawn.expect("a stream was given, so a path was drawn"),
        extra_paths,
    )
}

/// A state of the walk: autozygous on one template, or a pair of templates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum State {
    Autozygous(usize),
    Pair(usize, usize),
}

fn walk_with_draw(
    sites: &[WalkSite<'_>],
    templates: usize,
    parameters: WalkParameters,
    uniforms: Option<&mut Uniforms>,
    extra: Option<(&mut Uniforms, usize)>,
    mut reads: Option<&mut DrawReads<'_>>,
) -> (
    WalkResult,
    Option<Vec<HaplotypePair>>,
    Vec<Vec<HaplotypePair>>,
) {
    assert!(
        sites
            .windows(2)
            .all(|pair| pair[0].position <= pair[1].position),
        "a walk's sites must not descend in position"
    );
    let share = parameters.autozygous_share;
    assert!(
        share > 0.0 && share < 1.0 && templates > 0,
        "a walk needs an autozygous share strictly between 0 and 1 and at least one template: \
         got {share} and {templates}"
    );
    let k = templates;
    let k_f = k as f64;
    let mut genotype_pairs = GenotypePairs::default();
    let terms: Vec<SiteTerms> = sites
        .iter()
        .map(|site| {
            site_terms(
                site,
                k,
                parameters.mismatch,
                genotype_pairs.of(site.allele_count),
            )
        })
        .collect();

    // Forward: each site's state distribution given the sites up to it, normalised; the scales
    // multiply to the window's likelihood.
    let mut forward_autozygous: Vec<Vec<f64>> = Vec::with_capacity(sites.len());
    let mut forward_pairs: Vec<Vec<f64>> = Vec::with_capacity(sites.len());
    let mut log_likelihood = 0.0;
    for (index, site_terms) in terms.iter().enumerate() {
        let (mut autozygous, mut pairs) = if index == 0 {
            (
                vec![share / k_f; k],
                vec![(1.0 - share) / (k_f * k_f); k * k],
            )
        } else {
            let gap = (sites[index].position - sites[index - 1].position) as f64;
            let (keep, stay_autozygous, stay_other) = persistence(gap, &parameters);
            let previous_autozygous = &forward_autozygous[index - 1];
            let previous_pairs = &forward_pairs[index - 1];
            let mass_autozygous: f64 = previous_autozygous.iter().sum();
            let mass_pairs = lane_sum(previous_pairs);
            let (rows, columns) = margins(previous_pairs, k);
            let autozygous = previous_autozygous
                .iter()
                .map(|value| {
                    let within = keep * value + (1.0 - keep) * mass_autozygous / k_f;
                    stay_autozygous * within + (1.0 - stay_other) * mass_pairs / k_f
                })
                .collect();
            let pairs = pair_transition(
                previous_pairs,
                &rows,
                &columns,
                PairTransition::new(keep, mass_pairs, k_f, stay_other),
                (1.0 - stay_autozygous) * mass_autozygous / (k_f * k_f),
            );
            (autozygous, pairs)
        };
        for (value, emission) in autozygous.iter_mut().zip(&site_terms.autozygous) {
            *value *= emission;
        }
        for (value, emission) in pairs.iter_mut().zip(&site_terms.pairs) {
            *value *= emission;
        }
        let scale: f64 = autozygous.iter().sum::<f64>() + lane_sum(&pairs);
        let inverse = 1.0 / scale;
        for value in autozygous.iter_mut().chain(pairs.iter_mut()) {
            *value *= inverse;
        }
        log_likelihood += float::ln(scale);
        forward_autozygous.push(autozygous);
        forward_pairs.push(pairs);
    }

    // The round's own path from `uniforms`, then each extra path from `extra_uniforms`, one
    // uniform a site from the last site back. All are drawn in one sweep, so each site's totals
    // are summed once for every path.
    let mut streams: Vec<Vec<f64>> = Vec::new();
    let draws_own = uniforms.is_some();
    if let Some(uniforms) = uniforms {
        streams.push((0..sites.len()).map(|_| uniforms.next_uniform()).collect());
    }
    if let Some((extra_uniforms, count)) = extra {
        for _ in 0..count {
            streams.push(
                (0..sites.len())
                    .map(|_| extra_uniforms.next_uniform())
                    .collect(),
            );
        }
    }
    let mut paths: Vec<Vec<HaplotypePair>> = draw_state_paths(
        sites,
        &forward_autozygous,
        &forward_pairs,
        k,
        &parameters,
        &streams,
    )
    .iter()
    .map(|states| {
        haplotypes_along(
            sites,
            &terms,
            states,
            &mut genotype_pairs,
            reads.as_deref_mut(),
        )
    })
    .collect();
    let drawn = draws_own.then(|| paths.remove(0));
    let extra_paths = paths;

    // Backward, right to left, and each site's posterior as soon as its backward vector exists.
    let mut posteriors = vec![Vec::new(); sites.len()];
    let mut neighbour_predictions = vec![Vec::new(); sites.len()];
    let mut autozygous_probability = vec![0.0; sites.len()];
    let mut last_template_copy = vec![0.0; sites.len()];
    let mut backward_autozygous = vec![1.0; k];
    let mut backward_pairs = vec![1.0; k * k];
    let mut pairs_of_alleles = genotype_pairs;
    for index in (0..sites.len()).rev() {
        if index + 1 < sites.len() {
            let next = &terms[index + 1];
            let gap = (sites[index + 1].position - sites[index].position) as f64;
            let (keep, stay_autozygous, stay_other) = persistence(gap, &parameters);
            let weighted_autozygous: Vec<f64> = backward_autozygous
                .iter()
                .zip(&next.autozygous)
                .map(|(value, emission)| value * emission)
                .collect();
            let weighted_pairs: Vec<f64> = backward_pairs
                .iter()
                .zip(&next.pairs)
                .map(|(value, emission)| value * emission)
                .collect();
            let total_autozygous: f64 = weighted_autozygous.iter().sum();
            let total_pairs = lane_sum(&weighted_pairs);
            let (rows, columns) = margins(&weighted_pairs, k);
            backward_autozygous = weighted_autozygous
                .iter()
                .map(|value| {
                    let within = keep * value + (1.0 - keep) * total_autozygous / k_f;
                    stay_autozygous * within + (1.0 - stay_autozygous) * total_pairs / (k_f * k_f)
                })
                .collect();
            backward_pairs = pair_transition(
                &weighted_pairs,
                &rows,
                &columns,
                PairTransition::new(keep, total_pairs, k_f, stay_other),
                (1.0 - stay_other) * total_autozygous / k_f,
            );
            let total: f64 = backward_autozygous.iter().sum::<f64>() + lane_sum(&backward_pairs);
            let inverse = 1.0 / total;
            for value in backward_autozygous
                .iter_mut()
                .chain(backward_pairs.iter_mut())
            {
                *value *= inverse;
            }
        }
        let site = &sites[index];
        let site_terms = &terms[index];
        let alleles = site.allele_count;

        // The state posterior with this site's own emission divided back out, so the sample's
        // likelihood enters once, below (§3.6).
        //
        // **A state whose emission is exactly zero contributes zero**, not 0/0. It happens where the
        // reads rule out every genotype the state allows — a deep heterozygote's homozygous
        // likelihoods underflow to 0 at 30× — and then each of the emission's non-negative terms,
        // q̃ · q̃ · L(g), is zero, so the state adds nothing to any genotype whatever its weight.
        let without_emission = |both: f64, emission: f64| {
            if emission > 0.0 { both / emission } else { 0.0 }
        };
        let state_autozygous: Vec<f64> = (0..k)
            .map(|state| {
                without_emission(
                    forward_autozygous[index][state] * backward_autozygous[state],
                    site_terms.autozygous[state],
                )
            })
            .collect();
        let state_pairs: Vec<f64> = (0..k * k)
            .map(|state| {
                without_emission(
                    forward_pairs[index][state] * backward_pairs[state],
                    site_terms.pairs[state],
                )
            })
            .collect();
        let both_autozygous: f64 = (0..k)
            .map(|state| forward_autozygous[index][state] * backward_autozygous[state])
            .sum();
        let both_pairs = lane_dot(&forward_pairs[index], &backward_pairs);
        autozygous_probability[index] = both_autozygous / (both_autozygous + both_pairs);
        let last = k - 1;
        let with_last: f64 = forward_autozygous[index][last] * backward_autozygous[last]
            + (0..k)
                .map(|other| {
                    let row =
                        forward_pairs[index][last * k + other] * backward_pairs[last * k + other];
                    let column = if other == last {
                        0.0
                    } else {
                        forward_pairs[index][other * k + last] * backward_pairs[other * k + last]
                    };
                    row + column
                })
                .sum::<f64>();
        last_template_copy[index] = with_last / (both_autozygous + both_pairs);

        // The pair states' weight on each ordered allele pair, Q̃ᵀ · W · Q̃ (§3.6): first each
        // first-copy template's weight spread over the second copy's alleles (K² · A), then over
        // the first copy's (K · A²) — not K² · A² by summing every state over every pair.
        let blurred = &site_terms.blurred;
        let mut second_copy = vec![0.0; k * alleles];
        for first in 0..k {
            for second in 0..k {
                let weight = state_pairs[first * k + second];
                for b in 0..alleles {
                    second_copy[first * alleles + b] += weight * blurred[second * alleles + b];
                }
            }
        }
        let mut pair_weight = vec![0.0; alleles * alleles];
        for first in 0..k {
            for a in 0..alleles {
                let first_share = blurred[first * alleles + a];
                for b in 0..alleles {
                    pair_weight[a * alleles + b] += first_share * second_copy[first * alleles + b];
                }
            }
        }
        let mut autozygous_weight = vec![0.0; alleles];
        for state in 0..k {
            for a in 0..alleles {
                autozygous_weight[a] += state_autozygous[state] * blurred[state * alleles + a];
            }
        }
        let mut prediction: Vec<f64> = pairs_of_alleles
            .of(alleles)
            .iter()
            .map(|&(a, b)| {
                let mut weight = pair_weight[a * alleles + b];
                if a == b {
                    weight += autozygous_weight[a];
                } else {
                    weight += pair_weight[b * alleles + a];
                }
                weight
            })
            .collect();
        let mut posterior: Vec<f64> = prediction
            .iter()
            .zip(site.likelihoods)
            .map(|(weight, likelihood)| weight * likelihood)
            .collect();
        let total: f64 = posterior.iter().sum();
        for value in &mut posterior {
            *value /= total;
        }
        let total: f64 = prediction.iter().sum();
        for value in &mut prediction {
            *value /= total;
        }
        posteriors[index] = posterior;
        neighbour_predictions[index] = prediction;
    }

    (
        WalkResult {
            posteriors,
            neighbour_predictions,
            autozygous: autozygous_probability,
            last_template_copy,
            log_likelihood,
        },
        drawn,
        extra_paths,
    )
}

/// **Draws every path's states backwards from the forward probabilities** (phasing spec §2.1),
/// each state given the one after it. `streams` holds each path's uniforms, one a site from the
/// last site back. Each site's totals ([`SiteDraw`]) are summed once and serve every path.
fn draw_state_paths(
    sites: &[WalkSite<'_>],
    forward_autozygous: &[Vec<f64>],
    forward_pairs: &[Vec<f64>],
    k: usize,
    parameters: &WalkParameters,
    streams: &[Vec<f64>],
) -> Vec<Vec<State>> {
    let mut paths: Vec<Vec<State>> = vec![Vec::with_capacity(sites.len()); streams.len()];
    for (drawn, index) in (0..sites.len()).rev().enumerate() {
        let site = SiteDraw::new(&forward_autozygous[index], &forward_pairs[index], k);
        let persistence = (index + 1 < sites.len()).then(|| {
            persistence(
                (sites[index + 1].position - sites[index].position) as f64,
                parameters,
            )
        });
        for (path, stream) in paths.iter_mut().zip(streams) {
            let next = persistence.zip(path.last().copied());
            path.push(site.draw(next, stream[drawn]));
        }
    }
    for path in &mut paths {
        path.reverse();
    }
    paths
}

/// **One site's forward probabilities, arranged so a state can be drawn given the next in about
/// `K` steps** rather than a scan of all `K + K²` weights.
///
/// Given the next state, a state's weight is its forward probability times the transition, and
/// the transition is a share of every state of a kind, plus — towards a pair — a share of the next
/// pair's row and of its column, plus a share of the next state itself (spec §3.5). A uniform
/// first picks one of those parts by its mass, then a state inside it: by bisection over the
/// pairs' running totals for a share of every pair, by a scan of `K` for a row, a column or the
/// autozygous states. That draws from exactly the distribution the full scan does; the running
/// totals and the margins are summed once a site, for every path drawn through it.
struct SiteDraw<'a> {
    autozygous: &'a [f64],
    pairs: &'a [f64],
    k: usize,
    /// The pairs' running totals, row by row.
    running: Vec<f64>,
    rows: Vec<f64>,
    columns: Vec<f64>,
    autozygous_mass: f64,
}

impl<'a> SiteDraw<'a> {
    fn new(autozygous: &'a [f64], pairs: &'a [f64], k: usize) -> Self {
        let mut total = 0.0;
        let running = pairs
            .iter()
            .map(|value| {
                total += value;
                total
            })
            .collect();
        let (rows, columns) = margins(pairs, k);
        Self {
            autozygous,
            pairs,
            k,
            running,
            rows,
            columns,
            autozygous_mass: autozygous.iter().sum(),
        }
    }

    fn pair_mass(&self) -> f64 {
        self.running.last().copied().unwrap_or(0.0)
    }

    /// The state a uniform `u` draws here, given the next site's state and the persistence
    /// `(keep, stay autozygous, stay other)` between the two; `None` at the last site.
    fn draw(&self, next: Option<((f64, f64, f64), State)>, u: f64) -> State {
        let k = self.k;
        let k_f = k as f64;
        let any_template =
            |rest: f64| State::Autozygous(pick_at(self.autozygous.iter().copied(), rest));
        let any_pair = |rest: f64| {
            let position = pick_in_running(&self.running, rest);
            State::Pair(position / k, position % k)
        };
        match next {
            None => match choose_part(&[self.autozygous_mass, self.pair_mass()], u) {
                (0, rest) => any_template(rest),
                (_, rest) => any_pair(rest),
            },
            Some(((keep, stay_autozygous, stay_other), State::Pair(first, second))) => {
                let to_autozygous = (1.0 - stay_autozygous) / (k_f * k_f);
                let none_kept = stay_other * (1.0 - keep) * (1.0 - keep) / (k_f * k_f);
                let one_kept = stay_other * keep * (1.0 - keep) / k_f;
                let both_kept = stay_other * keep * keep;
                let masses = [
                    self.autozygous_mass * to_autozygous,
                    self.pair_mass() * none_kept,
                    self.rows[first] * one_kept,
                    self.columns[second] * one_kept,
                    self.pairs[first * k + second] * both_kept,
                ];
                match choose_part(&masses, u) {
                    (0, rest) => any_template(rest / to_autozygous),
                    (1, rest) => any_pair(rest / none_kept),
                    (2, rest) => {
                        let row = &self.pairs[first * k..][..k];
                        State::Pair(first, pick_at(row.iter().copied(), rest / one_kept))
                    }
                    (3, rest) => {
                        let column = self.pairs.iter().skip(second).step_by(k).copied();
                        State::Pair(pick_at(column, rest / one_kept), second)
                    }
                    _ => State::Pair(first, second),
                }
            }
            Some(((keep, stay_autozygous, stay_other), State::Autozygous(template))) => {
                let jumped = stay_autozygous * (1.0 - keep) / k_f;
                let to_pairs = (1.0 - stay_other) / k_f;
                let masses = [
                    self.autozygous_mass * jumped,
                    self.autozygous[template] * stay_autozygous * keep,
                    self.pair_mass() * to_pairs,
                ];
                match choose_part(&masses, u) {
                    (0, rest) => any_template(rest / jumped),
                    (1, _) => State::Autozygous(template),
                    (_, rest) => any_pair(rest / to_pairs),
                }
            }
        }
    }
}

/// **Which of a site's weight parts a uniform `u` falls in**, by their `masses` in order, and how
/// far into that part's mass it falls; a part of zero mass is never chosen.
///
/// # Panics
///
/// If every mass is zero.
fn choose_part(masses: &[f64], u: f64) -> (usize, f64) {
    let total: f64 = masses.iter().sum();
    assert!(total > 0.0, "a draw needs a positive weight");
    let target = u * total;
    let part = pick_at(masses.iter().copied(), target);
    let before: f64 = masses[..part].iter().sum();
    (part, (target - before).max(0.0))
}

/// **The two haplotypes along one drawn path** (phasing spec §2.1), each site's joint oriented by
/// the sample's reads where `reads` links it (reads-in-draw spec §2.2).
fn haplotypes_along(
    sites: &[WalkSite<'_>],
    terms: &[SiteTerms],
    states: &[State],
    genotype_pairs: &mut GenotypePairs,
    reads: Option<&mut DrawReads<'_>>,
) -> Vec<HaplotypePair> {
    let mut joints: Vec<Option<Vec<f64>>> = sites
        .iter()
        .zip(terms)
        .zip(states)
        .map(|((site, site_terms), &state)| {
            let alleles = site.allele_count;
            let mut by_alleles = vec![0.0; alleles * alleles];
            for (&(first, second), &likelihood) in
                genotype_pairs.of(alleles).iter().zip(site.likelihoods)
            {
                by_alleles[first * alleles + second] = likelihood;
                by_alleles[second * alleles + first] = likelihood;
            }
            joint_along(state, &site_terms.blurred, &by_alleles, alleles)
        })
        .collect();
    if let Some(reads) = reads {
        orient_by_reads(sites, &mut joints, reads);
    }
    sites
        .iter()
        .zip(&joints)
        .map(|(site, joint)| copies_of(joint.as_deref(), site.allele_count))
        .collect()
}

/// **The heterozygote a chain reads at a site**, its two alleles the lower first, or `None` when
/// the site is not orientable: its heterozygous mass is below `het_mass_min` (reads-in-draw spec
/// §2.2). The pair is the one with the most mass over both orders, ties to the first visited.
fn orientable_pair(joint: &[f64], alleles: usize, het_mass_min: f64) -> Option<(usize, usize)> {
    let mut heterozygous = 0.0;
    let mut best = None;
    let mut most = f64::NEG_INFINITY;
    for a in 0..alleles {
        for b in a + 1..alleles {
            let mass = joint[a * alleles + b] + joint[b * alleles + a];
            heterozygous += mass;
            if mass > most {
                most = mass;
                best = Some((a, b));
            }
        }
    }
    best.filter(|_| heterozygous >= het_mass_min)
}

/// **Draws the orientation of every read-linked group along one path and moves each grouped
/// site's heterozygous mass onto it** (reads-in-draw spec §2.2). A group is a run of consecutive
/// orientable sites, each linked to the next by the sample's fragments; a site outside every group
/// of two or more keeps its joint.
fn orient_by_reads(
    sites: &[WalkSite<'_>],
    joints: &mut [Option<Vec<f64>>],
    reads: &mut DrawReads<'_>,
) {
    let orientable: Vec<(usize, (usize, usize))> = joints
        .iter()
        .enumerate()
        .filter(|&(site, _)| !(reads.kept_out)(site))
        .filter_map(|(site, joint)| {
            let joint = joint.as_deref()?;
            orientable_pair(joint, sites[site].allele_count, reads.het_mass_min)
                .map(|pair| (site, pair))
        })
        .collect();
    let mut start = 0;
    let mut links = Vec::new();
    for index in 0..=orientable.len() {
        let link = (index > start && index < orientable.len())
            .then(|| {
                let (one, one_pair) = orientable[index - 1];
                let (other, other_pair) = orientable[index];
                (reads.link)(one, one_pair, other, other_pair)
            })
            .flatten();
        match link {
            Some(log_odds) => links.push(log_odds),
            None if index > start => {
                if index - start > 1 {
                    orient_group(
                        joints,
                        sites,
                        &orientable[start..index],
                        &links,
                        &mut reads.uniforms,
                    );
                }
                start = index;
                links.clear();
            }
            None => {}
        }
    }
}

/// **One group's chain** (reads-in-draw spec §2.2): own weights from each joint's two orders, the
/// link weight `e^{±L/2}` between neighbours, a forward pass normalised at each site, and a draw
/// backwards from the last; then each site's heterozygous mass moves onto the drawn order.
/// `links[i]` joins `group[i]` and `group[i + 1]`.
fn orient_group(
    joints: &mut [Option<Vec<f64>>],
    sites: &[WalkSite<'_>],
    group: &[(usize, (usize, usize))],
    links: &[f64],
    uniforms: &mut Uniforms,
) {
    debug_assert_eq!(links.len() + 1, group.len(), "one link between neighbours");
    // Own weights: [straight, swapped], straight putting the lower allele on copy one.
    let own: Vec<[f64; 2]> = group
        .iter()
        .map(|&(site, _)| {
            let alleles = sites[site].allele_count;
            // PANIC-FREE: only sites with a joint are orientable.
            let joint = joints[site]
                .as_deref()
                .expect("an orientable site has a joint");
            let mut weights = [0.0; 2];
            for a in 0..alleles {
                for b in a + 1..alleles {
                    weights[0] += joint[a * alleles + b];
                    weights[1] += joint[b * alleles + a];
                }
            }
            weights
        })
        .collect();
    let same = |log_odds: f64, one: usize, other: usize| {
        float::exp(if one == other { log_odds } else { -log_odds } / 2.0)
    };
    let mut forward: Vec<[f64; 2]> = Vec::with_capacity(group.len());
    for (index, weights) in own.iter().enumerate() {
        let mut next = *weights;
        if let Some(previous) = forward.last() {
            for (orientation, value) in next.iter_mut().enumerate() {
                *value *= (0..2)
                    .map(|before| previous[before] * same(links[index - 1], before, orientation))
                    .sum::<f64>();
            }
        }
        let total = next[0] + next[1];
        forward.push([next[0] / total, next[1] / total]);
    }
    let mut drawn = vec![0; group.len()];
    for index in (0..group.len()).rev() {
        let weights = match drawn.get(index + 1) {
            Some(&after) => [0, 1].map(|orientation| {
                forward[index][orientation] * same(links[index], orientation, after)
            }),
            None => forward[index],
        };
        drawn[index] = pick(&weights, uniforms.next_uniform());
    }
    for (&(site, _), &orientation) in group.iter().zip(&drawn) {
        let alleles = sites[site].allele_count;
        // PANIC-FREE: as above.
        let joint = joints[site]
            .as_mut()
            .expect("an orientable site has a joint");
        for a in 0..alleles {
            for b in a + 1..alleles {
                let (kept, emptied) = if orientation == 0 {
                    (a * alleles + b, b * alleles + a)
                } else {
                    (b * alleles + a, a * alleles + b)
                };
                joint[kept] += joint[emptied];
                joint[emptied] = 0.0;
            }
        }
    }
}

/// **The two copies' haplotypes at one site, given the state the path is in there**: the joint of
/// copy one's allele `a` and copy two's `b`, `q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})` in a pair state
/// and `q̃_k(a) · L({a, a})` on the diagonal in an autozygous one, normalised; each copy is its
/// margin. `blurred` holds the templates' blurred shares, `by_alleles` the likelihood as a
/// symmetric `A × A` matrix. Where every pair is ruled out each copy is flat.
#[cfg(test)]
fn copies_along(
    state: State,
    blurred: &[f64],
    by_alleles: &[f64],
    alleles: usize,
) -> HaplotypePair {
    copies_of(
        joint_along(state, blurred, by_alleles, alleles).as_deref(),
        alleles,
    )
}

/// **Each copy's haplotype, the margins of `joint`**; flat where there is no joint.
fn copies_of(joint: Option<&[f64]>, alleles: usize) -> HaplotypePair {
    let Some(joint) = joint else {
        // Every allele pair ruled out — the state's emission was zero (§3.6's guard).
        return [
            vec![1.0 / alleles as f64; alleles],
            vec![1.0 / alleles as f64; alleles],
        ];
    };
    let mut first_copy = vec![0.0; alleles];
    let mut second_copy = vec![0.0; alleles];
    for a in 0..alleles {
        for b in 0..alleles {
            let value = joint[a * alleles + b];
            first_copy[a] += value;
            second_copy[b] += value;
        }
    }
    [first_copy, second_copy]
}

/// **The joint of [`copies_along`] at one site**, normalised to sum to one; `None` where every
/// allele pair is ruled out.
fn joint_along(
    state: State,
    blurred: &[f64],
    by_alleles: &[f64],
    alleles: usize,
) -> Option<Vec<f64>> {
    let mut joint = vec![0.0; alleles * alleles];
    match state {
        State::Pair(first, second) => {
            for a in 0..alleles {
                for b in 0..alleles {
                    joint[a * alleles + b] = blurred[first * alleles + a]
                        * blurred[second * alleles + b]
                        * by_alleles[a * alleles + b];
                }
            }
        }
        State::Autozygous(template) => {
            for a in 0..alleles {
                joint[a * alleles + a] =
                    blurred[template * alleles + a] * by_alleles[a * alleles + a];
            }
        }
    }
    let total: f64 = joint.iter().sum();
    if total <= 0.0 {
        return None;
    }
    for value in &mut joint {
        *value /= total;
    }
    Some(joint)
}

/// **A pair state's transition, split by how many copies keep their template**: both
/// (`keep²`), one (`keep · (1 − keep)` spread over the `K` templates the other may jump to), or
/// neither (`(1 − keep)²` spread over all `K²` pairs, times the mass it comes from), all times the
/// chance of staying out of the autozygous state.
#[derive(Clone, Copy)]
struct PairTransition {
    stay_other: f64,
    both_kept: f64,
    one_kept: f64,
    none_kept: f64,
}

impl PairTransition {
    fn new(keep: f64, pair_mass: f64, k_f: f64, stay_other: f64) -> Self {
        Self {
            stay_other,
            both_kept: keep * keep,
            one_kept: keep * (1.0 - keep) / k_f,
            none_kept: (1.0 - keep) * (1.0 - keep) * pair_mass / (k_f * k_f),
        }
    }
}

/// **One step of the pair states** (spec §3.5), forward or backward: each pair's mass from
/// `values` (`k × k`, row by row) with its row and column sums, moved by `transition`, plus
/// `from_autozygous`, the share arriving at every pair from the autozygous states.
fn pair_transition(
    values: &[f64],
    rows: &[f64],
    columns: &[f64],
    transition: PairTransition,
    from_autozygous: f64,
) -> Vec<f64> {
    let PairTransition {
        stay_other,
        both_kept,
        one_kept,
        none_kept,
    } = transition;
    let column_terms: Vec<f64> = columns.iter().map(|column| one_kept * column).collect();
    let mut out = Vec::with_capacity(values.len());
    for (row, &row_sum) in values.chunks_exact(columns.len()).zip(rows) {
        let row_term = one_kept * row_sum + none_kept;
        out.extend(row.iter().zip(&column_terms).map(|(&value, column_term)| {
            stay_other * (both_kept * value + row_term + column_term) + from_autozygous
        }));
    }
    out
}

/// How many running totals [`lane_sum`] and [`lane_dot`] keep.
const LANES: usize = 8;

/// **A sum in [`LANES`] lanes**: each lane adds every eighth value, then the lanes are added in
/// order and the values past the last whole eight after them. The order is fixed by the length
/// alone, so the bits are the same on every machine and pool width, and the lanes vectorise where
/// one running total — each addition waiting for the one before — cannot (`float_portability.md`).
fn lane_sum(values: &[f64]) -> f64 {
    let (chunks, rest) = values.as_chunks::<LANES>();
    let mut lanes = [0.0; LANES];
    for chunk in chunks {
        for (lane, value) in lanes.iter_mut().zip(chunk) {
            *lane += value;
        }
    }
    lanes.iter().sum::<f64>() + rest.iter().sum::<f64>()
}

/// **[`lane_sum`] of the products of `one` and `other`**, pair by pair.
fn lane_dot(one: &[f64], other: &[f64]) -> f64 {
    let ((one_chunks, one_rest), (other_chunks, other_rest)) =
        (one.as_chunks::<LANES>(), other.as_chunks::<LANES>());
    let mut lanes = [0.0; LANES];
    for (a, b) in one_chunks.iter().zip(other_chunks) {
        for ((lane, a), b) in lanes.iter_mut().zip(a).zip(b) {
            *lane += a * b;
        }
    }
    let rest: f64 = one_rest.iter().zip(other_rest).map(|(a, b)| a * b).sum();
    lanes.iter().sum::<f64>() + rest
}

/// Row and column sums of a `k × k` matrix stored row by row, each summed in index order.
fn margins(matrix: &[f64], k: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rows = Vec::with_capacity(k);
    let mut columns = vec![0.0; k];
    for row in matrix.chunks_exact(k) {
        let mut row_total = 0.0;
        for (column, &value) in columns.iter_mut().zip(row) {
            row_total += value;
            *column += value;
        }
        rows.push(row_total);
    }
    (rows, columns)
}

#[cfg(test)]
mod tests;
