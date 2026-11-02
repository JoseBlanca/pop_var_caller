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

use super::random::{Uniforms, pick};
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
    walk_with_draw(sites, templates, parameters, None).0
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
) -> (WalkResult, Vec<HaplotypePair>) {
    let (result, drawn) = walk_with_draw(sites, templates, parameters, Some(uniforms));
    (
        result,
        // PANIC-FREE: `walk_with_draw` draws whenever it is given a stream, as it is here.
        drawn.expect("a stream was given, so a path was drawn"),
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
) -> (WalkResult, Option<Vec<HaplotypePair>>) {
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
            let mass_pairs: f64 = previous_pairs.iter().sum();
            let (rows, columns) = margins(previous_pairs, k);
            let autozygous = previous_autozygous
                .iter()
                .map(|value| {
                    let within = keep * value + (1.0 - keep) * mass_autozygous / k_f;
                    stay_autozygous * within + (1.0 - stay_other) * mass_pairs / k_f
                })
                .collect();
            let mut pairs = vec![0.0; k * k];
            for first in 0..k {
                for second in 0..k {
                    let within = keep * keep * previous_pairs[first * k + second]
                        + keep * (1.0 - keep) * (rows[first] + columns[second]) / k_f
                        + (1.0 - keep) * (1.0 - keep) * mass_pairs / (k_f * k_f);
                    pairs[first * k + second] = stay_other * within
                        + (1.0 - stay_autozygous) * mass_autozygous / (k_f * k_f);
                }
            }
            (autozygous, pairs)
        };
        for (value, emission) in autozygous.iter_mut().zip(&site_terms.autozygous) {
            *value *= emission;
        }
        for (value, emission) in pairs.iter_mut().zip(&site_terms.pairs) {
            *value *= emission;
        }
        let scale: f64 = autozygous.iter().sum::<f64>() + pairs.iter().sum::<f64>();
        for value in autozygous.iter_mut().chain(pairs.iter_mut()) {
            *value /= scale;
        }
        log_likelihood += float::ln(scale);
        forward_autozygous.push(autozygous);
        forward_pairs.push(pairs);
    }

    let drawn = uniforms.map(|uniforms| {
        draw_haplotypes(
            sites,
            &terms,
            &forward_autozygous,
            &forward_pairs,
            k,
            &parameters,
            &mut genotype_pairs,
            uniforms,
        )
    });

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
            let total_pairs: f64 = weighted_pairs.iter().sum();
            let (rows, columns) = margins(&weighted_pairs, k);
            backward_autozygous = weighted_autozygous
                .iter()
                .map(|value| {
                    let within = keep * value + (1.0 - keep) * total_autozygous / k_f;
                    stay_autozygous * within + (1.0 - stay_autozygous) * total_pairs / (k_f * k_f)
                })
                .collect();
            let mut pairs = vec![0.0; k * k];
            for first in 0..k {
                for second in 0..k {
                    let within = keep * keep * weighted_pairs[first * k + second]
                        + keep * (1.0 - keep) * (rows[first] + columns[second]) / k_f
                        + (1.0 - keep) * (1.0 - keep) * total_pairs / (k_f * k_f);
                    pairs[first * k + second] =
                        stay_other * within + (1.0 - stay_other) * total_autozygous / k_f;
                }
            }
            backward_pairs = pairs;
            let total: f64 =
                backward_autozygous.iter().sum::<f64>() + backward_pairs.iter().sum::<f64>();
            for value in backward_autozygous
                .iter_mut()
                .chain(backward_pairs.iter_mut())
            {
                *value /= total;
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
        let both_pairs: f64 = (0..k * k)
            .map(|state| forward_pairs[index][state] * backward_pairs[state])
            .sum();
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
    )
}

/// The backward draw and the re-estimation of [`walk_and_draw`].
#[allow(clippy::too_many_arguments)] // the forward walk's parts the draw reads, and the stream
fn draw_haplotypes(
    sites: &[WalkSite<'_>],
    terms: &[SiteTerms],
    forward_autozygous: &[Vec<f64>],
    forward_pairs: &[Vec<f64>],
    k: usize,
    parameters: &WalkParameters,
    genotype_pairs: &mut GenotypePairs,
    uniforms: &mut Uniforms,
) -> Vec<HaplotypePair> {
    let k_f = k as f64;
    let mut weights = vec![0.0; k + k * k];
    let mut states = Vec::with_capacity(sites.len());
    for index in (0..sites.len()).rev() {
        match states.last() {
            None => {
                weights[..k].copy_from_slice(&forward_autozygous[index]);
                weights[k..].copy_from_slice(&forward_pairs[index]);
            }
            Some(&next) => {
                let gap = (sites[index + 1].position - sites[index].position) as f64;
                let (keep, stay_autozygous, stay_other) = persistence(gap, parameters);
                match next {
                    State::Pair(next_first, next_second) => {
                        for (weight, value) in
                            weights[..k].iter_mut().zip(&forward_autozygous[index])
                        {
                            *weight = value * (1.0 - stay_autozygous) / (k_f * k_f);
                        }
                        for first in 0..k {
                            for second in 0..k {
                                let mut transition = (1.0 - keep) * (1.0 - keep) / (k_f * k_f);
                                if first == next_first {
                                    transition += keep * (1.0 - keep) / k_f;
                                }
                                if second == next_second {
                                    transition += keep * (1.0 - keep) / k_f;
                                }
                                if first == next_first && second == next_second {
                                    transition += keep * keep;
                                }
                                weights[k + first * k + second] = forward_pairs[index]
                                    [first * k + second]
                                    * stay_other
                                    * transition;
                            }
                        }
                    }
                    State::Autozygous(next_template) => {
                        for (template, (weight, value)) in weights[..k]
                            .iter_mut()
                            .zip(&forward_autozygous[index])
                            .enumerate()
                        {
                            let mut transition = (1.0 - keep) / k_f;
                            if template == next_template {
                                transition += keep;
                            }
                            *weight = value * stay_autozygous * transition;
                        }
                        for (weight, value) in weights[k..].iter_mut().zip(&forward_pairs[index]) {
                            *weight = value * (1.0 - stay_other) / k_f;
                        }
                    }
                }
            }
        }
        let picked = pick(&weights, uniforms.next_uniform());
        states.push(if picked < k {
            State::Autozygous(picked)
        } else {
            State::Pair((picked - k) / k, (picked - k) % k)
        });
    }
    states.reverse();

    sites
        .iter()
        .zip(terms)
        .zip(&states)
        .map(|((site, site_terms), &state)| {
            let alleles = site.allele_count;
            let mut by_alleles = vec![0.0; alleles * alleles];
            for (&(first, second), &likelihood) in
                genotype_pairs.of(alleles).iter().zip(site.likelihoods)
            {
                by_alleles[first * alleles + second] = likelihood;
                by_alleles[second * alleles + first] = likelihood;
            }
            copies_along(state, &site_terms.blurred, &by_alleles, alleles)
        })
        .collect()
}

/// **The two copies' haplotypes at one site, given the state the path is in there**: the joint of
/// copy one's allele `a` and copy two's `b`, `q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})` in a pair state
/// and `q̃_k(a) · L({a, a})` on the diagonal in an autozygous one, normalised; each copy is its
/// margin. `blurred` holds the templates' blurred shares, `by_alleles` the likelihood as a
/// symmetric `A × A` matrix. Where every pair is ruled out each copy is flat.
fn copies_along(
    state: State,
    blurred: &[f64],
    by_alleles: &[f64],
    alleles: usize,
) -> HaplotypePair {
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
        // Every allele pair ruled out — the state's emission was zero (§3.6's guard).
        return [
            vec![1.0 / alleles as f64; alleles],
            vec![1.0 / alleles as f64; alleles],
        ];
    }
    let mut first_copy = vec![0.0; alleles];
    let mut second_copy = vec![0.0; alleles];
    for a in 0..alleles {
        for b in 0..alleles {
            let value = joint[a * alleles + b] / total;
            first_copy[a] += value;
            second_copy[b] += value;
        }
    }
    [first_copy, second_copy]
}

/// Row and column sums of a `k × k` matrix stored row by row.
fn margins(matrix: &[f64], k: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rows = vec![0.0; k];
    let mut columns = vec![0.0; k];
    for first in 0..k {
        for second in 0..k {
            let value = matrix[first * k + second];
            rows[first] += value;
            columns[second] += value;
        }
    }
    (rows, columns)
}

#[cfg(test)]
mod tests;
