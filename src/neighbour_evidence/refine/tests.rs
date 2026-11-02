use super::*;
use crate::neighbour_evidence::templates::MAX_TEMPLATES;
use crate::neighbour_evidence::test_cohort::{CohortShape, test_cohort};

const SHAPE: CohortShape = CohortShape {
    samples: 30,
    founders: 4,
    records: 3_000,
    spacing: 300,
    segment_bases: 300_000.0,
    inbred: true,
    depth: 1.0,
    read_error: 0.01,
    seed: 5,
    silent_from: None,
};

fn argmax(values: &[f64]) -> usize {
    (0..values.len()).fold(0, |best, index| {
        if values[index] > values[best] {
            index
        } else {
            best
        }
    })
}

/// **§4.2's rule**: evidence on both sides within reach; one side where the contig ends within
/// reach on the other; a record at the same position counting on neither side.
#[test]
fn a_flat_sample_needs_evidence_on_both_sides_or_a_contig_end() {
    assert!(has_evidence_nearby(&[100, 300], 200, 150, None));
    assert!(!has_evidence_nearby(&[100, 300], 200, 50, None));
    assert!(!has_evidence_nearby(&[100], 200, 150, None));
    assert!(has_evidence_nearby(&[100], 200, 150, Some(300)));
    assert!(!has_evidence_nearby(&[200], 200, 150, Some(1_000)));
    assert!(!has_evidence_nearby(&[300], 50, 100, None));
    assert!(has_evidence_nearby(&[300], 50, 300, None));
    assert!(!has_evidence_nearby(&[], 50, 300, Some(100)));
}

/// **On a drawn inbred cohort at one read a record the pass removes errors and calls the sites
/// without reads**: 30 samples from 4 founders over 900 kb.
///
/// - At records where a sample has a read, the pass is wrong at least a tenth less often than
///   calling (on seed 5: 334 against 445 of 57,042). The cut is small because this cohort's only
///   calling errors are single wrong reads, and a template's allele, blurred by the 1% mismatch
///   rate, can outweigh a read only when many templates agree.
/// - At records where it has none — which calling writes `./.`, its prior's guess wrong about a
///   quarter of the time — the pass calls at least half, those whose `GQ` reaches 10, and is wrong
///   at most once in a hundred calls (on seed 5: 43 wrong of 18,883 called of 32,958, against the
///   guess's 8,871; without the `GQ` floor it called 32,941 and was wrong at 1,970).
#[test]
fn the_pass_fixes_errors_and_fills_sites_without_reads_on_an_inbred_cohort() {
    let cohort = test_cohort(SHAPE);
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
    let refined = refine_window(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    )
    .records;
    let (mut calling_errors, mut pass_errors, mut with_reads) = (0, 0, 0);
    let (mut without_reads, mut filled, mut filled_wrong, mut guessed_wrong) = (0, 0, 0, 0);
    for (index, record) in cohort.records.iter().enumerate() {
        for sample in 0..SHAPE.samples {
            let posterior = &refined[index].posteriors[sample * 3..sample * 3 + 3];
            let truth = cohort.truth[index][sample];
            match record.row(sample) {
                SampleRow::Informative => {
                    with_reads += 1;
                    calling_errors += usize::from(cohort.called[index][sample] != truth);
                    pass_errors += usize::from(argmax(posterior) != truth);
                }
                SampleRow::Flat => {
                    without_reads += 1;
                    guessed_wrong += usize::from(cohort.called[index][sample] != truth);
                    if refined[index].outcomes[sample] == Outcome::Called {
                        filled += 1;
                        filled_wrong += usize::from(argmax(posterior) != truth);
                    }
                }
                SampleRow::Absent => unreachable!("the drawn cohort writes every sample"),
            }
        }
    }
    assert!(
        pass_errors * 10 < calling_errors * 9,
        "with reads: calling {calling_errors} wrong, the pass {pass_errors}, of {with_reads}"
    );
    assert!(
        filled * 2 >= without_reads && filled_wrong * 100 <= filled,
        "without reads: {filled} of {without_reads} called, {filled_wrong} wrong; \
         the prior's guess {guessed_wrong} wrong"
    );
}

/// **A sample without reads that is probably copying the background is withheld** (§4.2): with the
/// limit, every no-read call copies the background with probability at most the limit, and each
/// call the limit removes is one that copied it more; the walk itself does not change. 10 samples
/// from 8 founders, so some haplotypes have no other sample to copy.
#[test]
fn a_no_read_call_copying_the_background_is_withheld() {
    let shape = CohortShape {
        samples: 10,
        founders: 8,
        ..SHAPE
    };
    let cohort = test_cohort(shape);
    let refine = |no_read| {
        let mut fits: Vec<SampleFit> = (0..shape.samples).map(|_| SampleFit::new(0.95)).collect();
        refine_window(
            &cohort.records,
            0..cohort.records.len(),
            None,
            &mut fits,
            MAX_TEMPLATES,
            no_read,
        )
        .records
    };
    let limited = refine(NoReadRule::default());
    let unlimited = refine(NoReadRule {
        max_background_copy: None,
        ..NoReadRule::default()
    });
    let mut removed = 0;
    for ((record, with_limit), without) in cohort.records.iter().zip(&limited).zip(&unlimited) {
        assert_eq!(with_limit.posteriors, without.posteriors);
        for sample in 0..shape.samples {
            let background = with_limit.background_copy[sample];
            match (with_limit.outcomes[sample], without.outcomes[sample]) {
                (Outcome::Called, Outcome::Called) => {}
                (Outcome::Withheld, Outcome::Called) => {
                    assert!(background > NO_READ_MAX_BACKGROUND_COPY);
                    removed += 1;
                }
                (same, other) => assert_eq!(same, other),
            }
            if with_limit.outcomes[sample] == Outcome::Called
                && record.row(sample) == SampleRow::Flat
            {
                assert!(background <= NO_READ_MAX_BACKGROUND_COPY);
            }
        }
    }
    assert!(removed > 0, "the limit withheld nothing");
}

/// **A sample that stops having reads is withheld once no evidence brackets the site** (§4.2):
/// sample 0 has no reads from record 1,000 on, so every one of its records from there is `./.`
/// with the pass's posterior as its `GP`, while the samples around it are still called. Before
/// record 1,000 every record where it has reads is called, and most where it has none.
#[test]
fn a_sample_without_evidence_on_both_sides_is_withheld() {
    let cohort = test_cohort(CohortShape {
        silent_from: Some((0, 1_000)),
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
    let refined = refine_window(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    )
    .records;
    for record in &refined[1_000..] {
        assert_eq!(record.outcomes[0], Outcome::Withheld);
        let total: f64 = record.posteriors[..3].iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "a withheld sample still has a posterior"
        );
    }
    let called_before = refined[..999]
        .iter()
        .filter(|record| record.outcomes[0] == Outcome::Called)
        .count();
    assert!(
        called_before > 800,
        "{called_before} of 999 called before the reads stop"
    );
    assert!(
        cohort.records[..999]
            .iter()
            .zip(&refined[..999])
            .filter(|(record, _)| record.row(0) == SampleRow::Informative)
            .all(|(_, refined)| refined.outcomes[0] == Outcome::Called)
    );
    assert!(
        cohort.records[1_000..]
            .iter()
            .zip(&refined[1_000..])
            .filter(|(record, _)| record.row(1) == SampleRow::Informative)
            .all(|(_, refined)| refined.outcomes[1] == Outcome::Called)
    );
}

/// **Only the core is returned, and a margin's record is not written**: the solve of a window's
/// middle third returns that third, record for record what the same samples' walks give there.
#[test]
fn only_the_core_is_returned() {
    let cohort = test_cohort(CohortShape {
        records: 600,
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
    let whole = refine_window(
        &cohort.records,
        0..600,
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    );
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
    let middle = refine_window(
        &cohort.records,
        200..400,
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    );
    assert_eq!(middle.records.len(), 200);
    assert_eq!(middle.records, whole.records[200..400]);
}

/// **The result does not depend on the thread count** (§2.3): one thread against four, bit for
/// bit.
#[test]
fn the_result_does_not_depend_on_the_thread_count() {
    let cohort = test_cohort(CohortShape {
        records: 800,
        ..SHAPE
    });
    let solve = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool");
        pool.install(|| {
            let mut fits: Vec<SampleFit> =
                (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
            refine_window(
                &cohort.records,
                0..800,
                None,
                &mut fits,
                MAX_TEMPLATES,
                NoReadRule::default(),
            )
        })
    };
    assert_eq!(solve(1), solve(4));
}

/// **Each sample's fill check counts its fills and tries the rest where its evidence settles it**
/// (`FillCheck`): on a drawn inbred cohort at 3×, `filled` is the sample's no-read calls written,
/// the check runs at many sites, and the neighbours' prediction there disagrees with the settled
/// genotype at most once in twenty — the true fill error rate of this cohort is far lower.
#[test]
fn the_fill_check_counts_the_fills_and_tries_the_settled_sites() {
    let shape = CohortShape {
        depth: 3.0,
        ..SHAPE
    };
    let cohort = test_cohort(shape);
    let mut fits: Vec<SampleFit> = (0..shape.samples).map(|_| SampleFit::new(0.95)).collect();
    let solved = refine_window(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    );
    for (sample, check) in solved.checks.iter().enumerate() {
        let filled = cohort
            .records
            .iter()
            .zip(&solved.records)
            .filter(|(record, refined)| {
                record.row(sample) == SampleRow::Flat && refined.outcomes[sample] == Outcome::Called
            })
            .count();
        assert_eq!(check.filled as usize, filled, "sample {sample}");
    }
    let checked: u32 = solved.checks.iter().map(|check| check.checked).sum();
    let disagreed: u32 = solved.checks.iter().map(|check| check.disagreed).sum();
    assert!(
        checked > 10_000 && disagreed * 20 <= checked,
        "{disagreed} of {checked} checked sites disagreed"
    );
}

/// **A fill check's lower bound separates chance from a poor sample**: one wrong in 30 is a share
/// of 3.3% but plausibly under 2%; ten in 100 is not.
#[test]
fn a_few_wrong_sites_do_not_make_a_poor_sample() {
    let check = |checked, disagreed| FillCheck {
        filled: 0,
        checked,
        disagreed,
    };
    assert!(check(30, 1).wrong_share_at_least() < 0.02);
    assert!(check(100, 10).wrong_share_at_least() > 0.02);
    assert!(check(0, 0).wrong_share_at_least() < 0.0);
}

/// Wrong genotypes at records where samples have reads, and at records where they have none and
/// the pass called them, from a solve's posteriors.
fn errors(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    refined: &[RefinedRecord],
    samples: usize,
) -> (usize, usize, usize) {
    let (mut with_reads, mut filled, mut filled_wrong) = (0, 0, 0);
    for (index, record) in cohort.records.iter().enumerate() {
        for sample in 0..samples {
            let posterior = &refined[index].posteriors[sample * 3..sample * 3 + 3];
            let wrong = argmax(posterior) != cohort.truth[index][sample];
            match record.row(sample) {
                SampleRow::Informative => with_reads += usize::from(wrong),
                SampleRow::Flat if refined[index].outcomes[sample] == Outcome::Called => {
                    filled += 1;
                    filled_wrong += usize::from(wrong);
                }
                _ => {}
            }
        }
    }
    (with_reads, filled, filled_wrong)
}

/// **Phasing removes errors an outbred cohort's genotype templates cannot** (phasing spec §1): 30
/// outbred samples from 4 founders at 3×, solved in one round and in [`PHASING_ROUNDS`]. A
/// heterozygous template tells a copy only "either allele"; its haplotypes say which. The rounds
/// leave at most half the wrong genotypes where samples have reads (on seed 5: 867 against 2,795).
#[test]
fn phasing_removes_errors_on_an_outbred_cohort() {
    let shape = CohortShape {
        inbred: false,
        depth: 3.0,
        records: 1_500,
        ..SHAPE
    };
    let cohort = test_cohort(shape);
    let solve = |phasing| {
        let mut fits: Vec<SampleFit> = (0..shape.samples).map(|_| SampleFit::new(0.05)).collect();
        refine_window_in_rounds(
            &cohort.records,
            0..cohort.records.len(),
            None,
            &mut fits,
            MAX_TEMPLATES,
            NoReadRule::default(),
            phasing,
        )
        .records
    };
    let (one_round, _, _) = errors(&cohort, &solve(Phasing::none()), shape.samples);
    let (phased, _, _) = errors(&cohort, &solve(Phasing::rounds_for(0)), shape.samples);
    assert!(
        phased * 2 <= one_round,
        "wrong where samples have reads: one round {one_round}, phased {phased}"
    );
}

/// **Phasing gives the same result at any thread count** (phasing spec §5): the draws are keyed by
/// window, sample and round, and the rounds are synchronous.
#[test]
fn phasing_does_not_depend_on_the_thread_count() {
    let cohort = test_cohort(CohortShape {
        records: 600,
        inbred: false,
        ..SHAPE
    });
    let solve = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool");
        pool.install(|| {
            let mut fits: Vec<SampleFit> =
                (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
            refine_window_in_rounds(
                &cohort.records,
                100..500,
                None,
                &mut fits,
                MAX_TEMPLATES,
                NoReadRule::default(),
                Phasing::rounds_for(17),
            )
        })
    };
    assert_eq!(solve(1), solve(4));
}

/// **With phasing the fill check runs at hidden records only** (phasing spec §4.1): no sample is
/// checked at more than [`HIDDEN_RECORDS_PER_SAMPLE`] records, and most samples at some.
#[test]
fn phasing_checks_fills_at_hidden_records_only() {
    let cohort = test_cohort(CohortShape {
        depth: 6.0,
        records: 1_000,
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.95)).collect();
    let solved = refine_window_in_rounds(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
        Phasing::rounds_for(3),
    );
    assert!(
        solved
            .checks
            .iter()
            .all(|check| check.checked as usize <= HIDDEN_RECORDS_PER_SAMPLE)
    );
    let checked = solved
        .checks
        .iter()
        .filter(|check| check.checked > 0)
        .count();
    assert!(checked * 2 > SHAPE.samples, "{checked} samples checked");
}

/// **The phased window matches the Python reference to 1e-9** (phasing spec §8.1): the window
/// `phaser_reference.py` drew and solved — 8 samples, 300 records, some flat and some absent — with
/// the same random stream, state order, hidden records and haplotype choice, every sample's fit
/// frozen at the reference's grid point. Every written posterior agrees, so a disagreement is a bug,
/// not the seed.
#[test]
fn the_phased_window_matches_the_python_reference() {
    use crate::neighbour_evidence::records::SampleColumns;
    let text = include_str!("../testdata/phased_window.txt");
    let mut header: Vec<&str> = Vec::new();
    let mut inputs: Vec<(u64, f64, Vec<Option<(Vec<u32>, Vec<f64>)>>)> = Vec::new();
    let mut expected: Vec<(usize, usize, Vec<f64>)> = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields[0] {
            "window" => header = fields,
            "record" => inputs.push((
                fields[1].parse().expect("a position"),
                fields[2].parse().expect("a frequency"),
                fields[3..]
                    .iter()
                    .map(|column| {
                        (*column != ".").then(|| {
                            let (pl, gp) = column.split_once('|').expect("PL|GP");
                            (
                                pl.split(',').map(|x| x.parse().expect("a PL")).collect(),
                                gp.split(',').map(|x| x.parse().expect("a GP")).collect(),
                            )
                        })
                    })
                    .collect(),
            )),
            "posterior" => expected.push((
                fields[1].parse().expect("a record"),
                fields[2].parse().expect("a sample"),
                fields[3..]
                    .iter()
                    .map(|x| x.parse().expect("a probability"))
                    .collect(),
            )),
            _ => {}
        }
    }
    let samples: usize = header[1].parse().expect("samples");
    let share: f64 = header[3].parse().expect("a share");
    let chosen: usize = header[4].parse().expect("a grid point");
    let rounds: usize = header[5].parse().expect("rounds");
    let window_key: u64 = header[6].parse().expect("a key");
    let records: Vec<SiteRecord> = inputs
        .iter()
        .map(|(position, frequency, columns)| {
            let columns: Vec<Option<SampleColumns<'_>>> = columns
                .iter()
                .map(|column| {
                    column.as_ref().map(|(pl, gp)| SampleColumns {
                        phred_likelihoods: pl,
                        genotype_posteriors: gp,
                    })
                })
                .collect();
            SiteRecord::new(*position, &[*frequency], &columns)
        })
        .collect();
    let mut fits: Vec<SampleFit> = (0..samples)
        .map(|_| SampleFit::frozen_at(share, chosen))
        .collect();
    let solved = refine_window_in_rounds(
        &records,
        0..records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
        Phasing { rounds, window_key },
    );
    assert!(
        expected.len() > 1_000,
        "the fixture holds {} posteriors",
        expected.len()
    );
    let mut worst: f64 = 0.0;
    for (record, sample, posterior) in &expected {
        let refined = &solved.records[*record];
        assert_ne!(
            refined.outcomes[*sample],
            Outcome::Unchanged,
            "record {record}, sample {sample}"
        );
        let ours = &refined.posteriors[sample * 3..sample * 3 + 3];
        for (a, b) in ours.iter().zip(posterior) {
            worst = worst.max((a - b).abs());
        }
    }
    assert!(
        worst < 1e-9,
        "largest difference from the reference: {worst:e}"
    );
}
