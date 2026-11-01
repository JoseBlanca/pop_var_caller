use super::*;

/// One window from the reference file, with what the prototype computed for it.
struct ReferenceCase {
    templates: usize,
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
                case.sites.push((
                    fields[1].parse().expect("a position"),
                    alleles,
                    values[..genotypes].to_vec(),
                    values[genotypes..].to_vec(),
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

/// **The walk reproduces the prototype** on three windows — biallelic; two to four alleles with
/// long gaps; two templates, five-allele sites and a high mismatch rate — to 1e-9 in every
/// posterior and in the log-likelihood.
#[test]
fn the_walk_reproduces_the_python_prototype() {
    let cases = reference_cases();
    assert_eq!(cases.len(), 3, "the reference file holds three windows");
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
        let result = walk(&sites, case.templates, case.parameters);
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
