use super::*;
use crate::neighbour_evidence::records::background_shares;

/// One window from the reference file, with what the prototype computed for it.
struct ReferenceCase {
    /// The sampled templates; the walk adds the background when `background` says so.
    templates: usize,
    background: bool,
    parameters: WalkParameters,
    sites: Vec<(u64, usize, Vec<f64>, Vec<f64>)>,
    log_likelihood: f64,
    posteriors: Vec<Vec<f64>>,
}

/// The windows `make_reference_walks.py` ran through the prototype's `hmm.haplotypes`
/// (`doc/devel/ng/research/neighbour_evidence/`), parsed.
fn reference_cases() -> Vec<ReferenceCase> {
    let text = include_str!("../testdata/reference_walks.txt");
    let mut cases: Vec<ReferenceCase> = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let numbers = |from: usize| -> Vec<f64> {
            fields[from..]
                .iter()
                .map(|field| field.parse().expect("a number"))
                .collect()
        };
        match fields[0] {
            "case" => cases.push(ReferenceCase {
                templates: fields[1].parse().expect("a template count"),
                background: fields[7] == "1",
                parameters: WalkParameters {
                    segment_bases: fields[3].parse().expect("a segment length"),
                    run_bases: fields[4].parse().expect("a run length"),
                    autozygous_share: fields[5].parse().expect("a share"),
                    mismatch: fields[6].parse().expect("a mismatch rate"),
                },
                sites: Vec::new(),
                log_likelihood: 0.0,
                posteriors: Vec::new(),
            }),
            "site" => {
                let case = cases.last_mut().expect("a site follows its case");
                let alleles: usize = fields[2].parse().expect("an allele count");
                let values = numbers(3);
                let genotypes = alleles * (alleles + 1) / 2;
                let shares_end = genotypes + case.templates * alleles;
                let mut shares = values[genotypes..shares_end].to_vec();
                if case.background {
                    let mut background = Vec::new();
                    background_shares(&values[shares_end..], &mut background);
                    shares.extend(background);
                }
                case.sites.push((
                    fields[1].parse().expect("a position"),
                    alleles,
                    values[..genotypes].to_vec(),
                    shares,
                ));
            }
            "expect" => {
                cases.last_mut().expect("a case").log_likelihood =
                    fields[1].parse().expect("a log-likelihood");
            }
            "posterior" => cases
                .last_mut()
                .expect("a case")
                .posteriors
                .push(numbers(1)),
            other => panic!("unexpected line kind {other}"),
        }
    }
    cases
}

/// **The walk reproduces the prototype** on four windows — biallelic; two to four alleles with
/// long gaps; two templates, five-allele sites and a high mismatch rate; two templates and the
/// background template built from the records' `AF` (§3.3) — to 1e-9 in every posterior and in the
/// log-likelihood.
#[test]
fn the_walk_reproduces_the_python_prototype() {
    let cases = reference_cases();
    assert_eq!(cases.len(), 4, "the reference file holds four windows");
    for (number, case) in cases.iter().enumerate() {
        let sites: Vec<WalkSite<'_>> = case
            .sites
            .iter()
            .map(|(position, alleles, likelihoods, shares)| WalkSite {
                position: *position,
                allele_count: *alleles,
                likelihoods,
                template_shares: shares,
            })
            .collect();
        let templates = case.templates + usize::from(case.background);
        let result = walk(&sites, templates, case.parameters);
        assert!(
            (result.log_likelihood - case.log_likelihood).abs() < 1e-9,
            "case {number}: log-likelihood {} against {}",
            result.log_likelihood,
            case.log_likelihood
        );
        for (site, (got, expected)) in result.posteriors.iter().zip(&case.posteriors).enumerate() {
            for (g, e) in got.iter().zip(expected) {
                assert!(
                    (g - e).abs() < 1e-9,
                    "case {number}, site {site}: {got:?} against {expected:?}"
                );
            }
        }
    }
}

/// **The pairs come out in VCF order**: `0/0, 0/1, 1/1, 0/2, 1/2, 2/2`.
#[test]
fn the_genotype_pairs_are_in_vcf_order() {
    let mut pairs = GenotypePairs::default();
    assert_eq!(
        pairs.of(3),
        &[(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]
    );
}

/// **A sample with no reads at a site takes its templates' allele**: flat likelihood in the middle
/// of a window where every template carries the alternative and the sample's reads say it is
/// homozygous alternative on both sides.
#[test]
fn a_site_without_reads_takes_the_templates_allele() {
    let confident_alternative = [1e-6, 1e-3, 1.0];
    let flat = [1.0, 1.0, 1.0];
    let alternative_templates = [0.0, 1.0, 0.0, 1.0];
    let sites: Vec<WalkSite<'_>> = [
        (1_000_u64, &confident_alternative),
        (2_000, &flat),
        (3_000, &confident_alternative),
    ]
    .iter()
    .map(|(position, likelihoods)| WalkSite {
        position: *position,
        allele_count: 2,
        likelihoods: &likelihoods[..],
        template_shares: &alternative_templates,
    })
    .collect();
    let parameters = WalkParameters {
        segment_bases: 300_000.0,
        run_bases: 1_000_000.0,
        autozygous_share: 0.8,
        mismatch: 0.01,
    };
    let result = walk(&sites, 2, parameters);
    let middle = &result.posteriors[1];
    assert!(middle[2] > 0.95, "the site without reads: {middle:?}");
}

/// **A long gap returns the state to its long-run share, it does not swap it** (spec §3.5): a
/// sample confidently autozygous at one site — homozygous where its two templates differ — and
/// with no reads at a site a gigabase on, is autozygous there with about the long-run share, 0.8,
/// not its complement.
#[test]
fn a_long_gap_returns_the_autozygous_state_to_its_share() {
    let homozygous = [1.0, 1e-9, 1e-9];
    let flat = [1.0, 1.0, 1.0];
    // Template 0 carries the reference, template 1 the alternative.
    let shares = [1.0, 0.0, 0.0, 1.0];
    let sites = [
        WalkSite {
            position: 1_000,
            allele_count: 2,
            likelihoods: &homozygous,
            template_shares: &shares,
        },
        WalkSite {
            position: 1_000_001_000,
            allele_count: 2,
            likelihoods: &flat,
            template_shares: &shares,
        },
    ];
    let parameters = WalkParameters {
        segment_bases: 300_000.0,
        run_bases: 1_000_000.0,
        autozygous_share: 0.8,
        mismatch: 0.01,
    };
    let result = walk(&sites, 2, parameters);
    assert!(
        (result.autozygous[1] - 0.8).abs() < 1e-6,
        "after a gigabase: {}",
        result.autozygous[1]
    );
}

/// **A deep heterozygote whose homozygous likelihoods are exactly zero gets a posterior, not NaN**:
/// at 30× those likelihoods underflow the window's `f32` store, every autozygous state's emission
/// is then 0, and dividing it back out was 0/0 (found on the founder simulator's 30× cell,
/// 2026-11-02).
#[test]
fn a_zero_emission_state_contributes_nothing() {
    let heterozygous = [0.0, 1.0, 0.0];
    let shares = [1.0, 0.0, 0.0, 1.0];
    let sites = [
        WalkSite {
            position: 1_000,
            allele_count: 2,
            likelihoods: &heterozygous,
            template_shares: &shares,
        },
        WalkSite {
            position: 2_000,
            allele_count: 2,
            likelihoods: &heterozygous,
            template_shares: &shares,
        },
    ];
    let parameters = WalkParameters {
        segment_bases: 300_000.0,
        run_bases: 1_000_000.0,
        autozygous_share: 0.95,
        mismatch: 0.01,
    };
    let result = walk(&sites, 2, parameters);
    for posterior in &result.posteriors {
        assert!(
            posterior.iter().all(|value| value.is_finite()),
            "{posterior:?}"
        );
        assert!((posterior[1] - 1.0).abs() < 1e-12, "{posterior:?}");
    }
}

/// **A path follows the only template that fits** (phasing spec §8.3): with no mismatch, a sample
/// homozygous for the alternative at every site, one template carrying the alternative and the
/// other carrying the reference, every drawn state copies the first, so both copies carry the
/// alternative everywhere — at any seed.
#[test]
fn a_drawn_path_follows_the_only_template_that_fits() {
    use crate::neighbour_evidence::random::{StreamPurpose, Uniforms};
    let likelihoods = [0.0, 0.0, 1.0];
    let shares = [0.0, 1.0, 1.0, 0.0];
    let sites: Vec<WalkSite<'_>> = (0..40)
        .map(|index| WalkSite {
            position: 1 + index * 1_000,
            allele_count: 2,
            likelihoods: &likelihoods,
            template_shares: &shares,
        })
        .collect();
    let parameters = WalkParameters {
        segment_bases: 30_000.0,
        run_bases: 100_000.0,
        autozygous_share: 0.5,
        mismatch: 0.0,
    };
    for seed in 0..20 {
        let mut uniforms = Uniforms::new(StreamPurpose::Path, seed, 0, 1);
        let (_, haplotypes) = walk_and_draw(&sites, 2, parameters, &mut uniforms);
        for pair in &haplotypes {
            for copy in pair {
                assert!((copy[1] - 1.0).abs() < 1e-12, "seed {seed}: {copy:?}");
            }
        }
    }
}

/// **The re-estimation is the owner's subtraction as a special case** (phasing spec §8.3): copy one
/// on a homozygous-reference template, copy two on a template even between the alleles, and reads
/// that allow only the heterozygote give copy one the reference and copy two the alternative, each
/// with probability `1 − ε/2`.
#[test]
fn the_re_estimation_subtracts_a_known_haplotype() {
    let mismatch = 0.01;
    let blurred = [1.0 - mismatch / 2.0, mismatch / 2.0, 0.5, 0.5];
    let only_heterozygous = [0.0, 1.0, 1.0, 0.0];
    let [first, second] = copies_along(State::Pair(0, 1), &blurred, &only_heterozygous, 2);
    assert!((first[0] - (1.0 - mismatch / 2.0)).abs() < 1e-12);
    assert!((second[1] - (1.0 - mismatch / 2.0)).abs() < 1e-12);
    // Autozygous on the reference template, with reads that allow only the alternative homozygote:
    // the blur alone carries it.
    let only_alternative = [0.0, 0.0, 0.0, 1.0];
    let [first, second] = copies_along(State::Autozygous(0), &blurred, &only_alternative, 2);
    assert!((first[1] - 1.0).abs() < 1e-12 && (second[1] - 1.0).abs() < 1e-12);
}

/// **At three alleles a path still follows the only template that fits** (code review 2026-11-02:
/// nothing tested phasing beyond two alleles): a sample homozygous for allele 2, one template
/// carrying allele 2 and one carrying allele 0, no mismatch — both copies carry allele 2 at every
/// site, the genotype table's six genotypes read in VCF order.
#[test]
fn a_drawn_path_follows_the_only_template_that_fits_at_three_alleles() {
    use crate::neighbour_evidence::random::{StreamPurpose, Uniforms};
    // 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
    let likelihoods = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let shares = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0];
    let sites: Vec<WalkSite<'_>> = (0..30)
        .map(|index| WalkSite {
            position: 1 + index * 1_000,
            allele_count: 3,
            likelihoods: &likelihoods,
            template_shares: &shares,
        })
        .collect();
    let parameters = WalkParameters {
        segment_bases: 30_000.0,
        run_bases: 100_000.0,
        autozygous_share: 0.5,
        mismatch: 0.0,
    };
    for seed in 0..10 {
        let mut uniforms = Uniforms::new(StreamPurpose::Path, seed, 0, 1);
        let (_, haplotypes) = walk_and_draw(&sites, 2, parameters, &mut uniforms);
        for pair in &haplotypes {
            for copy in pair {
                assert_eq!(copy.len(), 3);
                assert!((copy[2] - 1.0).abs() < 1e-12, "seed {seed}: {copy:?}");
            }
        }
    }
}

/// **The re-estimation at three alleles**: copy one on a template carrying allele 1, copy two on a
/// template even over the three, and reads allowing only the genotype 1/2 — copy two carries
/// allele 2 with probability `(1/3)·(1 − ε·2/3) / ((1/3)·(1 − ε·2/3) + (1/3)·(ε/3))`, given copy
/// one's allele 1, and allele 1 when copy one drew allele 2 instead; computed here.
#[test]
fn the_re_estimation_at_three_alleles() {
    let mismatch: f64 = 0.03;
    let carrying_one = [mismatch / 3.0, 1.0 - 2.0 * mismatch / 3.0, mismatch / 3.0];
    let even = [1.0 / 3.0; 3];
    let blurred: Vec<f64> = carrying_one.iter().chain(&even).copied().collect();
    // Symmetric 3×3 likelihood: only the pair {1, 2}.
    let mut by_alleles = [0.0; 9];
    by_alleles[3 + 2] = 1.0;
    by_alleles[2 * 3 + 1] = 1.0;
    let [first, second] = copies_along(State::Pair(0, 1), &blurred, &by_alleles, 3);
    let one_then_two = carrying_one[1] * even[2];
    let two_then_one = carrying_one[2] * even[1];
    let total = one_then_two + two_then_one;
    assert!((first[1] - one_then_two / total).abs() < 1e-12);
    assert!((second[2] - one_then_two / total).abs() < 1e-12);
    assert!((second[1] - two_then_one / total).abs() < 1e-12);
    assert!(first[0].abs() < 1e-12 && second[0].abs() < 1e-12);
}
