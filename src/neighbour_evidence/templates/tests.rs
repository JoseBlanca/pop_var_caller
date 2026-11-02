use super::*;
use crate::neighbour_evidence::records::SampleColumns;

const HOMOZYGOUS_REFERENCE: ([u32; 3], [f64; 3]) = ([0, 30, 60], [1.0, 0.0, 0.0]);
const HOMOZYGOUS_ALTERNATIVE: ([u32; 3], [f64; 3]) = ([60, 30, 0], [0.0, 0.0, 1.0]);

/// A window of `records` records where each sample is homozygous for the allele its pattern
/// names — `true` alternative — or absent (`None`).
/// Which records a sample carries the alternative at, or `None` for a sample written `.`.
type Pattern = Option<fn(usize) -> bool>;

fn window(patterns: &[Pattern], records: usize) -> Vec<SiteRecord> {
    (0..records)
        .map(|record| {
            let columns: Vec<Option<SampleColumns<'_>>> = patterns
                .iter()
                .map(|pattern| {
                    pattern.map(|carries| {
                        let (phred, posteriors) = if carries(record) {
                            &HOMOZYGOUS_ALTERNATIVE
                        } else {
                            &HOMOZYGOUS_REFERENCE
                        };
                        SampleColumns {
                            phred_likelihoods: phred,
                            genotype_posteriors: posteriors,
                        }
                    })
                })
                .collect();
            SiteRecord::new(1 + record as u64 * 100, &[0.5], &columns)
        })
        .collect()
}

fn even(record: usize) -> bool {
    record.is_multiple_of(2)
}

fn third(record: usize) -> bool {
    record.is_multiple_of(3)
}

/// **A sample copies the samples that carry its alleles**: samples 0–2 carry one pattern, 3–4
/// another; sample 0's two templates are 1 and 2, sample 3's first is 4.
#[test]
fn a_sample_copies_the_samples_that_carry_its_alleles() {
    let records = window(
        &[Some(even), Some(even), Some(even), Some(third), Some(third)],
        30,
    );
    assert_eq!(choose_templates(&records, 0, 2), Some(vec![1, 2]));
    assert_eq!(choose_templates(&records, 3, 1), Some(vec![4]));
}

/// **Ties break by sample order, and a sample is never its own template.**
#[test]
fn ties_break_by_sample_order_and_never_choose_the_sample() {
    let records = window(&[Some(even as fn(usize) -> bool) as Pattern; 5], 10);
    assert_eq!(choose_templates(&records, 2, 3), Some(vec![0, 1, 3]));
    assert_eq!(choose_templates(&records, 2, 16), Some(vec![0, 1, 3, 4]));
}

/// **An absent sample shares no record, so it ranks last**: after one that agrees and even after
/// one that disagrees everywhere, whose rate is at least measured.
#[test]
fn an_absent_sample_ranks_last() {
    fn never(_: usize) -> bool {
        false
    }
    fn always(_: usize) -> bool {
        true
    }
    let records = window(&[Some(never), Some(always), None, Some(never)], 10);
    assert_eq!(choose_templates(&records, 0, 3), Some(vec![3, 1, 2]));
}

const FLAT: ([u32; 3], [f64; 3]) = ([0, 0, 0], [0.5, 0.0, 0.5]);

/// A window of 400 records where `column(sample, record)` gives each sample's `PL` and `GP`.
fn window_of(
    samples: usize,
    column: impl Fn(usize, usize) -> &'static ([u32; 3], [f64; 3]),
) -> Vec<SiteRecord> {
    (0..400)
        .map(|record| {
            let columns: Vec<Option<SampleColumns<'_>>> = (0..samples)
                .map(|sample| {
                    let (phred, posteriors) = column(sample, record);
                    Some(SampleColumns {
                        phred_likelihoods: phred,
                        genotype_posteriors: posteriors,
                    })
                })
                .collect();
            SiteRecord::new(1 + record as u64 * 100, &[0.5], &columns)
        })
        .collect()
}

fn carrier(carries: bool) -> &'static ([u32; 3], [f64; 3]) {
    if carries {
        &HOMOZYGOUS_ALTERNATIVE
    } else {
        &HOMOZYGOUS_REFERENCE
    }
}

/// **A shallow relative is judged where it has reads**: sample 1 agrees with sample 0 wherever it
/// has reads, which is one record in five; sample 2 has reads everywhere and disagrees at one record
/// in ten. Counting the records without reads as disagreement would put sample 2 first.
#[test]
fn a_shallow_relative_is_not_taken_for_a_distant_sample() {
    let records = window_of(3, |sample, record| match sample {
        0 => carrier(record % 2 == 0),
        1 if record % 5 == 0 => carrier(record % 2 == 0),
        1 => &FLAT,
        _ => carrier((record % 2 == 0) != (record % 10 == 1)),
    });
    assert_eq!(choose_templates(&records, 0, 2), Some(vec![1, 2]));
}

/// **Between two that disagree about as little, the one with more data**: sample 1 agrees
/// everywhere it has reads, one record in four; sample 2 has reads everywhere and disagrees at one
/// record in a hundred. Sample 3, with reads everywhere and disagreeing at one in ten, is further than
/// both.
#[test]
fn a_near_tie_goes_to_the_sample_with_more_data() {
    let records = window_of(4, |sample, record| match sample {
        0 => carrier(record % 2 == 0),
        1 if record % 4 == 0 => carrier(record % 2 == 0),
        1 => &FLAT,
        2 => carrier((record % 2 == 0) != (record % 100 == 1)),
        _ => carrier((record % 2 == 0) != (record % 10 == 1)),
    });
    assert_eq!(choose_templates(&records, 0, 3), Some(vec![2, 1, 3]));
}

/// **A sample with no informative record in the window has no templates**, so it is not
/// re-called there (§3.3).
#[test]
fn a_sample_without_informative_records_has_no_templates() {
    let flat_phred = [0_u32; 3];
    let prior = [0.25, 0.5, 0.25];
    let records: Vec<SiteRecord> = (0..5)
        .map(|record| {
            let columns = [
                Some(SampleColumns {
                    phred_likelihoods: &flat_phred,
                    genotype_posteriors: &prior,
                }),
                Some(SampleColumns {
                    phred_likelihoods: &HOMOZYGOUS_REFERENCE.0,
                    genotype_posteriors: &HOMOZYGOUS_REFERENCE.1,
                }),
            ];
            SiteRecord::new(1 + record * 100, &[0.5], &columns)
        })
        .collect();
    assert_eq!(choose_templates(&records, 0, 16), None);
    assert_eq!(choose_templates(&records, 1, 16), Some(vec![0]));
}

/// **More than [`AGREEMENT_RECORDS`] informative records are compared on that many, spread over
/// the window**: a sample agreeing with sample 1 over the window's first half and with sample 2
/// over its second still ties them, where taking the first 500 records would pick sample 1.
#[test]
fn the_records_compared_are_spread_over_the_window() {
    fn first_half(record: usize) -> bool {
        record < 1_000
    }
    fn second_half(record: usize) -> bool {
        record >= 1_000
    }
    fn always(_: usize) -> bool {
        true
    }
    let records = window(&[Some(always), Some(first_half), Some(second_half)], 2_000);
    let templates = choose_templates(&records, 0, 2).expect("informative records");
    assert_eq!(templates, vec![1, 2], "a tie, broken by sample order");
    let records = window(&[Some(always), Some(second_half), Some(first_half)], 2_000);
    assert_eq!(choose_templates(&records, 0, 2), Some(vec![1, 2]));
}

/// **A phasing round's haplotype choice takes each sample's better copy, one a sample** (phasing
/// spec §2.3): sample 0 is homozygous for the alternative; sample 1 carries it on copy 1 only,
/// sample 2 on neither, sample 3 on both. Samples 1 and 3 tie at no disagreement and come first in
/// sample order, each with the copy that carries the alternative (copy 0 for sample 3, on the tie
/// between its two); sample 2 offers its copy 0 last.
#[test]
fn a_round_takes_each_samples_better_haplotype_once() {
    let records = window_of(4, |_, _| &HOMOZYGOUS_ALTERNATIVE);
    let compared = compared_records(&records, 0);
    assert_eq!(compared.len(), 400.min(AGREEMENT_RECORDS));
    let homozygous_alternative = [0.0, 0.0, 1.0];
    let reference = [1.0, 0.0];
    let alternative = [0.0, 1.0];
    let haplotype = |other: usize, copy: usize, _: usize| -> &[f64] {
        match (other, copy) {
            (1, 1) | (3, _) => &alternative,
            _ => &reference,
        }
    };
    let chosen = choose_haplotypes(
        &records,
        0,
        &compared,
        &|_| &homozygous_alternative,
        &haplotype,
        3,
    );
    assert_eq!(chosen, vec![(1, 1), (3, 0), (2, 0)]);
    let chosen = choose_haplotypes(
        &records,
        0,
        &compared,
        &|_| &homozygous_alternative,
        &haplotype,
        1,
    );
    assert_eq!(chosen, vec![(1, 1)]);
}
