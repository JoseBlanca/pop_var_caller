use super::*;
use crate::neighbour_evidence::random::StreamPurpose;
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
        let (_, haplotypes) = walk_and_draw(&sites, 2, parameters, &mut uniforms, None);
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
        let (_, haplotypes) = walk_and_draw(&sites, 2, parameters, &mut uniforms, None);
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

/// **Extra paths change nothing the round passes on** (phase output spec §2.3): the walk and the
/// round's own path are [`walk_and_draw`]'s exactly, bit for bit, and the extra paths are a
/// function of their own stream — the same stream gives the same paths, another gives others.
#[test]
fn extra_paths_leave_the_round_unchanged() {
    use crate::neighbour_evidence::random::{StreamPurpose, Uniforms};
    // A sample even between the alleles, two templates of opposite alleles: every path is a
    // genuine draw, so two streams disagree somewhere.
    let likelihoods = [0.2, 1.0, 0.2];
    let shares = [0.9, 0.1, 0.1, 0.9];
    let sites: Vec<WalkSite<'_>> = (0..60)
        .map(|index| WalkSite {
            position: 1 + index * 2_000,
            allele_count: 2,
            likelihoods: &likelihoods,
            template_shares: &shares,
        })
        .collect();
    let parameters = WalkParameters {
        segment_bases: 30_000.0,
        run_bases: 100_000.0,
        autozygous_share: 0.3,
        mismatch: 0.01,
    };
    let mut alone = Uniforms::new(StreamPurpose::Path, 5, 2, 3);
    let (walked, own) = walk_and_draw(&sites, 2, parameters, &mut alone, None);
    let mut with = Uniforms::new(StreamPurpose::Path, 5, 2, 3);
    let mut extra = Uniforms::new(StreamPurpose::ExtraPaths, 5, 2, 3);
    let (walked_again, own_again, paths) =
        walk_and_draw_paths(&sites, 2, parameters, &mut with, &mut extra, 4, None);
    assert_eq!(walked, walked_again);
    assert_eq!(own, own_again);
    assert_eq!(paths.len(), 4);
    let mut extra_again = Uniforms::new(StreamPurpose::ExtraPaths, 5, 2, 3);
    let mut with_again = Uniforms::new(StreamPurpose::Path, 5, 2, 3);
    let (_, _, paths_again) = walk_and_draw_paths(
        &sites,
        2,
        parameters,
        &mut with_again,
        &mut extra_again,
        4,
        None,
    );
    assert_eq!(paths, paths_again);
    assert!(
        paths.iter().any(|path| *path != own),
        "an extra path is a draw of its own"
    );
}

/// Walk sites carrying only their allele counts and positions, for the orientation chain, which
/// reads nothing else of them.
fn bare_sites(allele_counts: &[usize]) -> Vec<WalkSite<'static>> {
    allele_counts
        .iter()
        .enumerate()
        .map(|(index, &allele_count)| WalkSite {
            position: 1 + index as u64 * 100,
            allele_count,
            likelihoods: &[],
            template_shares: &[],
        })
        .collect()
}

/// A biallelic joint over `(0,0), (0,1), (1,0), (1,1)`.
fn biallelic(joint: [f64; 4]) -> Option<Vec<f64>> {
    Some(joint.to_vec())
}

/// Whether a biallelic site's joint was moved onto *straight*, copy one carrying allele 0.
fn straight(joint: &Option<Vec<f64>>) -> bool {
    let joint = joint.as_deref().expect("a joint");
    assert!(
        joint[1] == 0.0 || joint[2] == 0.0,
        "a grouped site carries one order: {joint:?}"
    );
    joint[2] == 0.0
}

/// One read fragment's `L` at a base error of 1 in 100, the error floor (read phase spec §2.5).
fn one_fragment() -> f64 {
    crate::read_phase::link::weight(40, 40)
}

/// **Runs the chain over `joints` with every consecutive pair linked by `log_odds`**, the hidden
/// sites `hidden`, from the orientation stream keyed by `seed`; returns the joints written.
fn oriented(
    sites: &[WalkSite<'_>],
    joints: &[Option<Vec<f64>>],
    log_odds: f64,
    hidden: &[usize],
    seed: u64,
) -> Vec<Option<Vec<f64>>> {
    let link = |_: usize, _: (usize, usize), _: usize, _: (usize, usize)| Some(log_odds);
    let is_hidden = |site: usize| hidden.contains(&site);
    let mut reads = DrawReads {
        link: &link,
        kept_out: &is_hidden,
        het_mass_min: 0.5,
        uniforms: Uniforms::new(StreamPurpose::Orientation, seed, 0, 1),
    };
    let mut joints = joints.to_vec();
    orient_by_reads(sites, &mut joints, &mut reads);
    joints
}

/// Asserts `hits` of `draws` lies within five standard deviations of the binomial with `p`.
fn assert_binomial(hits: usize, draws: usize, p: f64, what: &str) {
    let n = draws as f64;
    let spread = 5.0 * (n * p * (1.0 - p)).sqrt();
    assert!(
        (hits as f64 - n * p).abs() <= spread,
        "{what}: {hits} of {draws}, expected {:.1} ± {spread:.1}",
        n * p
    );
}

/// **Reads alone orient two sites whose templates are alike** (reads-in-draw spec §6): with equal
/// own weights and one fragment's `L` ≈ 3.9 saying cis, the two take the same orientation with
/// probability `e^L / (1 + e^L)` ≈ 0.98, and the first takes either with probability one half.
#[test]
fn reads_alone_orient_two_sites_alike_in_their_templates() {
    let sites = bare_sites(&[2, 2]);
    let alike = biallelic([0.0, 0.5, 0.5, 0.0]);
    let log_odds = one_fragment();
    let draws = 4_000;
    let (mut same, mut first_straight) = (0, 0);
    for seed in 0..draws {
        let joints = oriented(&sites, &[alike.clone(), alike.clone()], log_odds, &[], seed);
        same += usize::from(straight(&joints[0]) == straight(&joints[1]));
        first_straight += usize::from(straight(&joints[0]));
    }
    let p = float::exp(log_odds) / (1.0 + float::exp(log_odds));
    assert!((p - 0.98).abs() < 0.005, "one fragment's odds: {p}");
    assert_binomial(same, draws as usize, p, "same orientation");
    assert_binomial(first_straight, draws as usize, 0.5, "first straight");
}

/// **Two fragments do not overturn confident, differing templates; three do** (reads-in-draw spec
/// §2.3): the first site pinned straight, the second's templates favouring *swapped* by
/// `(0.995/0.005)²` ≈ 39,600, and reads saying cis. The second turns straight with probability
/// `e^L / (e^L + 39,600)`: about 0.06 at two fragments, 0.75 at three.
#[test]
fn three_fragments_overturn_confident_templates_and_two_do_not() {
    let sites = bare_sites(&[2, 2]);
    let pinned = biallelic([0.0, 1.0, 1e-30, 0.0]);
    let favour = float::powi(0.995 / 0.005, 2);
    let swapped = biallelic([0.0, 1.0 / (1.0 + favour), favour / (1.0 + favour), 0.0]);
    let draws = 4_000;
    for (fragments, low, high) in [(2.0, 0.03, 0.09), (3.0, 0.70, 0.80)] {
        let log_odds = fragments * one_fragment();
        let p = float::exp(log_odds) / (float::exp(log_odds) + favour);
        assert!(low < p && p < high, "{fragments} fragments: {p}");
        let turned = (0..draws)
            .filter(|&seed| {
                let joints = oriented(
                    &sites,
                    &[pinned.clone(), swapped.clone()],
                    log_odds,
                    &[],
                    seed,
                );
                assert!(straight(&joints[0]), "the pinned site stays straight");
                straight(&joints[1])
            })
            .count();
        assert_binomial(turned, draws as usize, p, "turned by reads");
    }
}

/// **A hidden site and one below the heterozygous mass are never in a group, and do not split it**
/// (reads-in-draw spec §2.2, §5): both keep their joints exactly, and the sites either side are
/// linked across them.
#[test]
fn hidden_and_homozygous_sites_stay_out_of_a_group() {
    let sites = bare_sites(&[2, 2, 2, 2]);
    let alike = biallelic([0.0, 0.5, 0.5, 0.0]);
    let below = biallelic([0.51, 0.245, 0.245, 0.0]);
    let joints = [alike.clone(), alike.clone(), below.clone(), alike.clone()];
    let mut same = 0;
    let draws: u64 = 200;
    for seed in 0..draws {
        let written = oriented(&sites, &joints, 30.0, &[1], seed);
        assert_eq!(written[1], alike, "the hidden site keeps its joint");
        assert_eq!(written[2], below, "the site below one half keeps its joint");
        same += usize::from(straight(&written[0]) == straight(&written[3]));
    }
    assert_eq!(
        same as u64, draws,
        "the sites either side are linked across"
    );
}

/// **A site with no link keeps its joint** — a group of one draws nothing — and **no links at all
/// leave every joint as it was, bit for bit, without reading the stream**.
#[test]
fn an_unlinked_site_keeps_its_joint() {
    let sites = bare_sites(&[2, 2]);
    let joints = vec![
        biallelic([0.1, 0.3, 0.5, 0.1]),
        biallelic([0.0, 0.7, 0.3, 0.0]),
    ];
    let link = |_: usize, _: (usize, usize), _: usize, _: (usize, usize)| None;
    let is_hidden = |_: usize| false;
    let mut reads = DrawReads {
        link: &link,
        kept_out: &is_hidden,
        het_mass_min: 0.5,
        uniforms: Uniforms::new(StreamPurpose::Orientation, 3, 0, 1),
    };
    let mut written = joints.clone();
    orient_by_reads(&sites, &mut written, &mut reads);
    assert_eq!(written, joints);
    let mut fresh = Uniforms::new(StreamPurpose::Orientation, 3, 0, 1);
    assert_eq!(reads.uniforms.next_uniform(), fresh.next_uniform());
}

/// **The written joint keeps every homozygous mass and moves the whole heterozygous mass onto the
/// drawn order** (reads-in-draw spec §2.2), at a three-allele site too.
#[test]
fn the_written_joint_keeps_homozygous_mass() {
    let sites = bare_sites(&[3, 2]);
    // Over (a, b) row by row: the 1/2 heterozygote carries most of the mass, both orders.
    let three = Some(vec![
        0.05, 0.02, 0.03, //
        0.04, 0.06, 0.30, //
        0.01, 0.40, 0.09,
    ]);
    let two = biallelic([0.1, 0.4, 0.4, 0.1]);
    for seed in 0..20 {
        let written = oriented(&sites, &[three.clone(), two.clone()], 5.0, &[], seed);
        let (before, after) = (three.as_deref().unwrap(), written[0].as_deref().unwrap());
        for a in 0..3 {
            assert_eq!(after[a * 3 + a], before[a * 3 + a], "homozygote {a}");
            for b in a + 1..3 {
                let pair = before[a * 3 + b] + before[b * 3 + a];
                let moved = after[a * 3 + b] + after[b * 3 + a];
                assert!((moved - pair).abs() < 1e-15, "pair {a}/{b}");
                assert!(after[a * 3 + b] == 0.0 || after[b * 3 + a] == 0.0);
            }
        }
    }
}

/// **At a three-allele site the link reads the heterozygote with the most mass**, its lower allele
/// first: a `1/2` heterozygote is read as alleles 1 and 2, never as the reference and another
/// (reads-in-draw spec §5).
#[test]
fn a_one_two_heterozygote_is_read_as_alleles_one_and_two() {
    let three = [
        0.0, 0.05, 0.05, //
        0.05, 0.0, 0.40, //
        0.05, 0.40, 0.0,
    ];
    assert_eq!(orientable_pair(&three, 3, 0.5), Some((1, 2)));
    assert_eq!(orientable_pair(&three, 3, 1.000_001), None);
    let sites = bare_sites(&[3, 2]);
    let seen = std::cell::RefCell::new(Vec::new());
    let link = |one: usize, one_pair: (usize, usize), other: usize, other_pair: (usize, usize)| {
        seen.borrow_mut().push((one, one_pair, other, other_pair));
        Some(1.0)
    };
    let is_hidden = |_: usize| false;
    let mut reads = DrawReads {
        link: &link,
        kept_out: &is_hidden,
        het_mass_min: 0.5,
        uniforms: Uniforms::new(StreamPurpose::Orientation, 1, 0, 1),
    };
    let mut joints = vec![Some(three.to_vec()), biallelic([0.0, 0.5, 0.5, 0.0])];
    orient_by_reads(&sites, &mut joints, &mut reads);
    assert_eq!(*seen.borrow(), vec![(0, (1, 2), 1, (0, 1))]);
}

/// **Reads that link nothing leave the draw bit for bit as without reads** (reads-in-draw spec
/// §1.3): the walk, the round's path and the extra paths.
#[test]
fn reads_that_link_nothing_change_no_bit() {
    let likelihoods = [0.2, 1.0, 0.2];
    let shares = [0.9, 0.1, 0.1, 0.9, 0.5, 0.5];
    let sites: Vec<WalkSite<'_>> = (0..60)
        .map(|index| WalkSite {
            position: 1 + index * 2_000,
            allele_count: 2,
            likelihoods: &likelihoods,
            template_shares: &shares,
        })
        .collect();
    let parameters = WalkParameters {
        segment_bases: 30_000.0,
        run_bases: 100_000.0,
        autozygous_share: 0.3,
        mismatch: 0.01,
    };
    let draw = |reads: Option<&mut DrawReads<'_>>| {
        let mut path = Uniforms::new(StreamPurpose::Path, 5, 2, 3);
        let mut extra = Uniforms::new(StreamPurpose::ExtraPaths, 5, 2, 3);
        walk_and_draw_paths(&sites, 3, parameters, &mut path, &mut extra, 4, reads)
    };
    let link = |_: usize, _: (usize, usize), _: usize, _: (usize, usize)| None;
    let is_hidden = |_: usize| false;
    let mut reads = DrawReads {
        link: &link,
        kept_out: &is_hidden,
        het_mass_min: 0.5,
        uniforms: Uniforms::new(StreamPurpose::Orientation, 5, 2, 3),
    };
    let without = draw(None);
    let with = draw(Some(&mut reads));
    assert_eq!(without.0, with.0);
    assert_eq!(without.1, with.1);
    assert_eq!(without.2, with.2);
}

/// **Reads linking every site leave the walk and the genotypes alone** (reads-in-draw spec §2.1):
/// only the haplotypes change, and they do.
#[test]
fn reads_change_the_haplotypes_and_not_the_walk() {
    let likelihoods = [0.05, 1.0, 0.05];
    let shares = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    let sites: Vec<WalkSite<'_>> = (0..40)
        .map(|index| WalkSite {
            position: 1 + index * 100,
            allele_count: 2,
            likelihoods: &likelihoods,
            template_shares: &shares,
        })
        .collect();
    let parameters = WalkParameters {
        segment_bases: 30_000.0,
        run_bases: 100_000.0,
        autozygous_share: 0.3,
        mismatch: 0.01,
    };
    let link = |_: usize, _: (usize, usize), _: usize, _: (usize, usize)| Some(20.0);
    let is_hidden = |_: usize| false;
    let mut reads = DrawReads {
        link: &link,
        kept_out: &is_hidden,
        het_mass_min: 0.5,
        uniforms: Uniforms::new(StreamPurpose::Orientation, 9, 0, 2),
    };
    let mut path = Uniforms::new(StreamPurpose::Path, 9, 0, 2);
    let (walked, without) = walk_and_draw(&sites, 3, parameters, &mut path, None);
    let mut path = Uniforms::new(StreamPurpose::Path, 9, 0, 2);
    let (walked_again, with) = walk_and_draw(&sites, 3, parameters, &mut path, Some(&mut reads));
    assert_eq!(walked, walked_again);
    // Templates even between the alleles: without reads both copies are the same margin.
    assert!(without.iter().all(|[one, two]| one == two));
    assert_ne!(with, without, "the reads changed the haplotypes");
    // With every pair linked at cis, every heterozygous site carries one orientation.
    let first_copy_reference: Vec<bool> = with.iter().map(|[one, _]| one[0] > 0.5).collect();
    assert!(
        first_copy_reference
            .iter()
            .all(|&straight| straight == first_copy_reference[0]),
        "{first_copy_reference:?}"
    );
}
