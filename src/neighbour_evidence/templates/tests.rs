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

/// [`choose_templates`]' samples, without the stretches each won.
fn chosen_samples(
    records: &[SiteRecord],
    sample: usize,
    max_templates: usize,
) -> Option<Vec<usize>> {
    choose_templates(records, sample, max_templates)
        .map(|chosen| chosen.into_iter().map(|template| template.source).collect())
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
    assert_eq!(chosen_samples(&records, 0, 2), Some(vec![1, 2]));
    assert_eq!(chosen_samples(&records, 3, 1), Some(vec![4]));
}

/// **Ties break by sample order, and a sample is never its own template.**
#[test]
fn ties_break_by_sample_order_and_never_choose_the_sample() {
    let records = window(&[Some(even as fn(usize) -> bool) as Pattern; 5], 10);
    assert_eq!(chosen_samples(&records, 2, 3), Some(vec![0, 1, 3]));
    assert_eq!(chosen_samples(&records, 2, 16), Some(vec![0, 1, 3, 4]));
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
    assert_eq!(chosen_samples(&records, 0, 3), Some(vec![3, 1, 2]));
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

/// **A shallow relative is judged where it has reads** (spec §3.3): sample 1 agrees with sample 0
/// wherever it has reads, one record in five; samples 2 and 3 have reads everywhere and each
/// disagrees at one record in four — sample 2 in the even stretches, sample 3 in the odd ones (a
/// stretch here is two records). Sample 1 wins every stretch it has a read in, 80 of 200: there it
/// beats whichever deep sample disagrees, as it could not if its records without reads counted as
/// disagreement. The deep samples win every stretch sample 1 has no read in — 120 for both — and
/// one in two of the others, so each ranks above it.
#[test]
fn a_shallow_relative_wins_the_stretches_it_has_reads_in() {
    let records = window_of(4, |sample, record| match sample {
        0 => carrier(record % 2 == 0),
        1 if record % 5 == 0 => carrier(record % 2 == 0),
        1 => &FLAT,
        2 => carrier((record % 2 == 0) != (record % 4 == 1)),
        _ => carrier((record % 2 == 0) != (record % 4 == 3)),
    });
    let won: Vec<(usize, usize)> = choose_templates(&records, 0, 3)
        .expect("informative records")
        .into_iter()
        .map(|template| (template.source, template.stretches_won))
        .collect();
    assert_eq!(won, vec![(2, 160), (3, 160), (1, 80)]);
}

/// One read of the reference, one of the alternative: shares 0.79 against 0.21.
const ONE_READ_REFERENCE: ([u32; 3], [f64; 3]) = ([0, 3, 13], [0.62, 0.34, 0.04]);
const ONE_READ_ALTERNATIVE: ([u32; 3], [f64; 3]) = ([13, 3, 0], [0.04, 0.34, 0.62]);
/// Reads that leave the genotype at the cohort's frequency, 0.5: shares even.
const UNSURE: ([u32; 3], [f64; 3]) = ([3, 0, 3], [0.25, 0.5, 0.25]);

/// **An unsure candidate is not taken for a close one** (spec §3.3, 2026-11-05): sample 0 has one
/// read a record, as a 1× F2 plant; sample 1 is certain and carries sample 0's allele at three
/// records in four, as a parent does; sample 2's shares sit at the cohort's frequency everywhere,
/// as a 1× sibling's do where it has no reads. Half the summed share difference ranked sample 2
/// first (0.29 against 0.36); against what the frequencies predict, sample 2 is no closer than an
/// unrelated sample and sample 1 is.
#[test]
fn an_unsure_candidate_is_not_taken_for_a_close_one() {
    let records = window_of(3, |sample, record| {
        let alternative = record % 2 == 0;
        match sample {
            0 if alternative => &ONE_READ_ALTERNATIVE,
            0 => &ONE_READ_REFERENCE,
            1 => carrier(alternative != (record % 4 == 1)),
            _ => &UNSURE,
        }
    });
    assert_eq!(chosen_samples(&records, 0, 2), Some(vec![1, 2]));
}

/// **Between two that agree as well where both are judged, the one with more data** (spec §3.3):
/// sample 1 agrees everywhere it has reads, one record in four; sample 2 has reads everywhere and
/// disagrees at one record in a hundred; sample 3 at one in ten. In a stretch where sample 1 has a
/// read, the deep samples agree there too and their score rests on twice the expected
/// disagreement, so they win it; sample 1 ranks last, after sample 3, which disagrees more but
/// carries data at every record the walk can copy.
#[test]
fn a_near_tie_goes_to_the_sample_with_more_data() {
    let records = window_of(4, |sample, record| match sample {
        0 => carrier(record % 2 == 0),
        1 if record % 4 == 0 => carrier(record % 2 == 0),
        1 => &FLAT,
        2 => carrier((record % 2 == 0) != (record % 100 == 1)),
        _ => carrier((record % 2 == 0) != (record % 10 == 1)),
    });
    assert_eq!(chosen_samples(&records, 0, 3), Some(vec![2, 3, 1]));
}

/// **The samples a copy copies win the stretches they match** (spec §3.3, owner 2026-11-05): sample
/// 0 copies sample 1 over the window's first half and sample 2 over its second, as an F2 plant
/// copies its parents either side of a crossover; each of those disagrees at every record of its
/// other half, a rate of ½ over the window. Samples 3–5 each disagree at one record in three
/// everywhere, a rate of ⅓, so the whole window's score ranks all three above samples 1 and 2. In
/// each two-record stretch the matching parent and the one of samples 3–5 that agrees there win;
/// each parent wins its 100 stretches, each of samples 3–5 about 67, and the parents are the two
/// templates.
#[test]
fn the_samples_a_copy_copies_win_the_stretches_they_match() {
    let records = window_of(6, |sample, record| {
        let target = record % 2 == 0;
        match sample {
            0 => carrier(target),
            1 => carrier(target != (record >= 200)),
            2 => carrier(target != (record < 200)),
            _ => carrier(target != (record % 3 == sample - 3)),
        }
    });
    let chosen = choose_templates(&records, 0, 5).expect("informative records");
    let won: Vec<(usize, usize)> = chosen
        .iter()
        .map(|template| (template.source, template.stretches_won))
        .collect();
    assert_eq!(won[..2], [(1, 100), (2, 100)]);
    for &(sample, stretches) in &won[2..] {
        assert!(
            (66..=67).contains(&stretches),
            "sample {sample} won {stretches}"
        );
    }
}

/// A [`Disagreements`] over `expected`, each candidate's observed disagreements one a record.
fn disagreements(expected: &[f64], candidates: &[&[Option<f64>]]) -> Disagreements {
    Disagreements {
        expected: expected.to_vec(),
        observed: candidates.iter().flat_map(|c| c.iter().copied()).collect(),
    }
}

/// **A stretch ends once it holds [`STRETCH_INFORMATION`]** of expected disagreement, counting
/// only records some candidate is informative at; the last holds what is left.
#[test]
fn a_stretch_ends_once_it_holds_enough_to_judge_by() {
    let expected = [0.5, 0.5, 0.3, 0.3, 0.3, 0.2, 0.4];
    let everywhere = [Some(0.0); 7];
    assert_eq!(
        disagreements(&expected, &[&everywhere]).stretches(),
        vec![0..2, 2..6, 6..7]
    );
    let mut gap = everywhere;
    gap[1] = None;
    assert_eq!(
        disagreements(&expected, &[&gap]).stretches(),
        vec![0..4, 4..7]
    );
}

/// **A stretch's winners are the [`STRETCH_WINNERS`] best and any tied with the last of them**; a
/// candidate not judged there does not win it. Each landing weight is the stretches won plus
/// [`LANDING_WEIGHT_BASE`].
#[test]
fn ties_with_the_last_winner_win_too() {
    let judged = disagreements(
        &[1.0],
        &[
            &[Some(0.0)],
            &[Some(0.5)],
            &[Some(0.5)],
            &[Some(1.0)],
            &[None],
        ],
    )
    .judged();
    let won: Vec<usize> = judged.iter().map(|one| one.stretches_won).collect();
    assert_eq!(won, vec![1, 1, 1, 0, 0]);
    let template = ChosenTemplate {
        source: 0,
        stretches_won: 3,
    };
    assert!((template.landing_weight() - (3.0 + LANDING_WEIGHT_BASE)).abs() < 1e-15);
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
    assert_eq!(chosen_samples(&records, 0, 16), None);
    assert_eq!(chosen_samples(&records, 1, 16), Some(vec![0]));
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
    let templates = chosen_samples(&records, 0, 2).expect("informative records");
    assert_eq!(templates, vec![1, 2], "a tie, broken by sample order");
    let records = window(&[Some(always), Some(second_half), Some(first_half)], 2_000);
    assert_eq!(chosen_samples(&records, 0, 2), Some(vec![1, 2]));
}

/// [`choose_haplotypes`]' haplotypes, without the stretches each won.
fn haplotype_sources(chosen: &[ChosenTemplate<(usize, usize)>]) -> Vec<(usize, usize)> {
    chosen.iter().map(|template| template.source).collect()
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
    let homozygous_alternative = vec![vec![0.0, 0.0, 1.0]; compared.len()];
    let reference = [1.0, 0.0];
    let alternative = [0.0, 1.0];
    let haplotype = |other: usize, copy: usize, _: usize, out: &mut Vec<f64>| {
        out.clear();
        out.extend_from_slice(match (other, copy) {
            (1, 1) | (3, _) => &alternative,
            _ => &reference,
        });
    };
    let chosen = choose_haplotypes(
        &records,
        0,
        &compared,
        &homozygous_alternative,
        &haplotype,
        3,
    );
    assert_eq!(haplotype_sources(&chosen), vec![(1, 1), (3, 0), (2, 0)]);
    let chosen = choose_haplotypes(
        &records,
        0,
        &compared,
        &homozygous_alternative,
        &haplotype,
        1,
    );
    assert_eq!(haplotype_sources(&chosen), vec![(1, 1)]);
}

/// **The haplotype choice at three alleles** (code review 2026-11-02): sample 0's previous genotype
/// is 2/2, so a haplotype lacks nothing only where it carries allele 2. Sample 1's copy 0 carries
/// allele 2 and its copy 1 allele 0; sample 2 carries allele 1 on both. Sample 1 offers copy 0 and
/// comes first.
#[test]
fn the_haplotype_choice_reads_every_allele() {
    use crate::neighbour_evidence::records::SampleColumns;
    let informative = [60_u32, 60, 60, 30, 30, 0];
    let posteriors = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let records: Vec<SiteRecord> = (0..50)
        .map(|record| {
            let columns: Vec<Option<SampleColumns<'_>>> = (0..3)
                .map(|_| {
                    Some(SampleColumns {
                        phred_likelihoods: &informative,
                        genotype_posteriors: &posteriors,
                    })
                })
                .collect();
            SiteRecord::new(1 + record * 100, &[0.2, 0.3], &columns)
        })
        .collect();
    let compared = compared_records(&records, 0);
    let previous = vec![posteriors.to_vec(); compared.len()];
    let allele = |carried: usize| -> [f64; 3] {
        let mut shares = [0.0; 3];
        shares[carried] = 1.0;
        shares
    };
    let (two, zero, one) = (allele(2), allele(0), allele(1));
    let haplotype = |other: usize, copy: usize, _: usize, out: &mut Vec<f64>| {
        out.clear();
        out.extend_from_slice(match (other, copy) {
            (1, 0) => &two,
            (1, 1) => &zero,
            _ => &one,
        });
    };
    let chosen = choose_haplotypes(&records, 0, &compared, &previous, &haplotype, 2);
    assert_eq!(haplotype_sources(&chosen), vec![(1, 0), (2, 0)]);
}
