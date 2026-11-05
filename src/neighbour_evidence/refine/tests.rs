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
/// leave at most half the wrong genotypes where samples have reads (on seed 5: 1,015 against 2,795, with up to 50 records a sample hidden for the fill check).
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
/// checked at more than [`HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW`] records, and most samples at some.
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
            .all(|check| check.checked as usize <= HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW)
    );
    let checked = solved
        .checks
        .iter()
        .filter(|check| check.checked > 0)
        .count();
    assert!(checked * 2 > SHAPE.samples, "{checked} samples checked");
}

/// **A record no walk uses, inserted mid-window, changes no other record's phased result**: the
/// window solve keeps record indices and site positions apart (they coincide in windows where every
/// record is a site, as in every other test here). From the code review of 2026-11-02.
#[test]
fn a_non_site_record_changes_no_other_phased_record() {
    use crate::neighbour_evidence::records::SampleColumns;
    let cohort = test_cohort(CohortShape {
        records: 600,
        inbred: false,
        ..SHAPE
    });
    let pl = [0_u32, 10, 40];
    let gp = [0.9, 0.09, 0.01];
    let mut columns: Vec<Option<SampleColumns<'_>>> = vec![None; SHAPE.samples];
    columns[0] = Some(SampleColumns {
        phred_likelihoods: &pl,
        genotype_posteriors: &gp,
    });
    let at = 50;
    let extra = SiteRecord::new(cohort.records[at].position, &[0.1], &columns);
    assert!(!extra.is_site());
    let mut with_extra = cohort.records.clone();
    with_extra.insert(at, extra);
    let solve = |records: &[SiteRecord]| {
        let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
        refine_window_in_rounds(
            records,
            0..records.len(),
            None,
            &mut fits,
            MAX_TEMPLATES,
            NoReadRule::default(),
            Phasing::rounds_for(17),
        )
    };
    let plain = solve(&cohort.records);
    let extended = solve(&with_extra);
    let mut kept = extended.records.clone();
    kept.remove(at);
    assert_eq!(plain.records, kept);
    assert_eq!(plain.checks, extended.checks);
}

/// **A phased window's posteriors, to the last bit** (code review 2026-11-02, float portability):
/// the 1e-9 reference test cannot see a one-bit difference between platforms' arithmetic, which
/// would change the written `GP`. This hash of every posterior's bits was recorded on macOS arm64;
/// the same test on Linux x86-64 must give it too. A deliberate change to the pass's arithmetic
/// re-records it.
#[test]
fn a_phased_windows_posteriors_are_the_same_bits_everywhere() {
    let cohort = test_cohort(CohortShape {
        records: 400,
        inbred: false,
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
    let solved = refine_window_in_rounds(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
        Phasing::rounds_for(window_key(3, 1_000)),
    );
    // FNV-1a over every posterior's bits, record by record.
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for record in &solved.records {
        for value in &record.posteriors {
            for byte in value.to_bits().to_le_bytes() {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    assert_eq!(
        hash, RECORDED_POSTERIOR_BITS,
        "the posteriors' bits are {hash:#018x}"
    );
}

/// The hash [`a_phased_windows_posteriors_are_the_same_bits_everywhere`] expects.
const RECORDED_POSTERIOR_BITS: u64 = 0x5a7f_cfa5_860e_1c2f;

/// **The fit runs in round 1 only** (phasing spec §2.4; code review 2026-11-02, a fit in every
/// round went unnoticed): two rounds and five leave every sample's fit in the same state, since
/// rounds 2 on read the parameters round 1 chose.
#[test]
fn only_round_one_fits() {
    let cohort = test_cohort(CohortShape {
        records: 500,
        inbred: false,
        ..SHAPE
    });
    let fits_after = |rounds: usize| {
        let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
        let _ = refine_window_in_rounds(
            &cohort.records,
            0..cohort.records.len(),
            None,
            &mut fits,
            MAX_TEMPLATES,
            NoReadRule::default(),
            Phasing::new(rounds, 9),
        );
        format!("{fits:?}")
    };
    assert_eq!(fits_after(2), fits_after(5));
}

/// **The fill check compares the neighbours' prediction, not the posterior, with the settled
/// genotype** (code review 2026-11-02: feeding it the posterior went unnoticed): every sample is
/// homozygous for the reference everywhere, except sample 0, whose deep reads settle the
/// alternative homozygote at one record. Its neighbours there predict the reference, so its check
/// counts exactly one disagreement; the posterior, which holds its reads, would count none.
#[test]
fn the_fill_check_compares_the_neighbours_prediction() {
    use crate::neighbour_evidence::records::SampleColumns;
    let reference = ([0_u32, 30, 60], [1.0, 0.0, 0.0]);
    let alternative = ([90_u32, 45, 0], [0.0, 0.0, 1.0]);
    let samples = 6;
    let odd = 100;
    let records: Vec<SiteRecord> = (0..200)
        .map(|record| {
            let columns: Vec<Option<SampleColumns<'_>>> = (0..samples)
                .map(|sample| {
                    let (pl, gp) = if sample == 0 && record == odd {
                        &alternative
                    } else {
                        &reference
                    };
                    Some(SampleColumns {
                        phred_likelihoods: pl,
                        genotype_posteriors: gp,
                    })
                })
                .collect();
            SiteRecord::new(1 + record as u64 * 200, &[0.1], &columns)
        })
        .collect();
    let mut fits: Vec<SampleFit> = (0..samples).map(|_| SampleFit::new(0.95)).collect();
    let solved = refine_window(
        &records,
        0..records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
    );
    assert!(solved.checks[0].checked > 100, "{:?}", solved.checks[0]);
    assert_eq!(solved.checks[0].disagreed, 1, "{:?}", solved.checks[0]);
}

/// **At most [`HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW`] records a sample are hidden, and at most
/// [`MAX_HIDDEN_SHARE`] of its informative core records** (code review 2026-11-02: neither limit was
/// reached by a test): a window where every record is informative hides 50 of 3,000 and 20 of 400.
#[test]
fn hidden_records_respect_both_limits() {
    for (records, expected) in [(3_000, HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW), (400, 20)] {
        let cohort = test_cohort(CohortShape {
            records,
            depth: 30.0,
            ..SHAPE
        });
        let layout = WindowLayout::new(&cohort.records, &(0..records));
        let hidden = draw_hidden_records(&cohort.records, &layout, 0, 5);
        assert_eq!(hidden.len(), expected, "{records} records");
        assert!(hidden.windows(2).all(|pair| pair[0] < pair[1]));
    }
}

/// **A phased window of two or three samples runs**, every template slot that exists filled, and
/// writes finite posteriors (code review 2026-11-02: no phased test below 30 samples).
#[test]
fn a_phased_window_of_a_few_samples_runs() {
    for samples in [2, 3] {
        let cohort = test_cohort(CohortShape {
            samples,
            founders: 2,
            records: 400,
            inbred: false,
            ..SHAPE
        });
        let mut fits: Vec<SampleFit> = (0..samples).map(|_| SampleFit::new(0.05)).collect();
        let solved = refine_window_in_rounds(
            &cohort.records,
            0..cohort.records.len(),
            None,
            &mut fits,
            MAX_TEMPLATES,
            NoReadRule::default(),
            Phasing::rounds_for(1),
        );
        assert!(
            solved
                .records
                .iter()
                .all(|record| record.posteriors.iter().all(|value| value.is_finite()))
        );
        assert!(
            solved
                .records
                .iter()
                .any(|record| record.outcomes.contains(&Outcome::Called))
        );
    }
}

/// **A window cannot be solved in no rounds.**
#[test]
#[should_panic(expected = "at least one round")]
fn phasing_refuses_zero_rounds() {
    let _ = Phasing::new(0, 0);
}

/// **A window's key tells windows apart by contig and by position**, so two windows never share
/// draws (code review 2026-11-02: replacing the stage's key with a constant went unnoticed).
#[test]
fn window_keys_differ_by_contig_and_position() {
    let keys = [
        window_key(0, 1),
        window_key(0, 2_000_001),
        window_key(1, 1),
        window_key(1, 2_000_001),
        window_key(24, 248_956_422),
    ];
    for (index, key) in keys.iter().enumerate() {
        assert!(!keys[index + 1..].contains(key), "{keys:?}");
    }
}

/// **The written phase against the truth** (phase output spec §7.4): on an outbred cohort of four
/// founders at 4×, the window phases about half the true heterozygotes (3,797 of 7,496 when
/// written), every block holds at least two records and is named by its first one's position, and
/// the links inside blocks are right (0 of 3,341 wrong when written; at most 1 in 100 allowed).
#[test]
fn the_written_phase_follows_the_true_copies() {
    let cohort = test_cohort(CohortShape {
        records: 800,
        inbred: false,
        depth: 4.0,
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
    let solved = refine_window_in_rounds(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
        Phasing::rounds_for(11),
    );
    let (mut heterozygous, mut phased, mut links, mut flips) = (0, 0, 0, 0);
    for sample in 0..SHAPE.samples {
        let mut previous: Option<(u64, bool)> = None;
        let mut block_sizes: std::collections::HashMap<u64, usize> = Default::default();
        for (index, record) in solved.records.iter().enumerate() {
            let [copy_0, copy_1] = cohort.truth_copies[index][sample];
            heterozygous += usize::from(copy_0 != copy_1);
            let Some(phase) = record.phases[sample] else {
                continue;
            };
            assert_ne!(
                phase.copy_0_allele, phase.copy_1_allele,
                "only heterozygotes are phased"
            );
            assert!(phase.phase_set <= cohort.records[index].position);
            *block_sizes.entry(phase.phase_set).or_default() += 1;
            if block_sizes[&phase.phase_set] == 1 {
                assert_eq!(
                    phase.phase_set, cohort.records[index].position,
                    "PS is the first POS"
                );
            }
            phased += 1;
            // Right where the record is truly heterozygous and copy 0 carries copy 0's allele.
            let right = copy_0 != copy_1 && phase.copy_0_allele == copy_0;
            if let Some((block, was_right)) = previous
                && block == phase.phase_set
            {
                links += 1;
                flips += usize::from(right != was_right);
            }
            previous = Some((phase.phase_set, right));
        }
        assert!(
            block_sizes.values().all(|&size| size >= 2),
            "no block of one"
        );
    }
    eprintln!(
        "{phased} of {heterozygous} true heterozygotes phased; {flips} flips in {links} links"
    );
    assert!(
        phased * 5 > heterozygous * 2,
        "{phased} of {heterozygous} phased"
    );
    assert!(flips * 100 <= links, "{flips} flips in {links} links");
}

/// **A phased window's written phase is the same bits everywhere** (phase output spec §7.5; code
/// review 2026-11-02: the posteriors' hash never reads the extra paths or the consensus): the
/// window of [`a_phased_windows_posteriors_are_the_same_bits_everywhere`] at 4×, FNV-1a over each
/// sample's `[copy 0's allele, copy 1's allele, phase set]` at each record as little-endian `u64`s,
/// `u64::MAX` three times where unphased. Recorded on macOS arm64 at 19102f1e and unchanged by the
/// review's fixes; to be checked on Linux x86-64.
#[test]
fn a_phased_windows_phase_is_the_same_bits_everywhere() {
    let cohort = test_cohort(CohortShape {
        records: 400,
        inbred: false,
        depth: 4.0,
        ..SHAPE
    });
    let mut fits: Vec<SampleFit> = (0..SHAPE.samples).map(|_| SampleFit::new(0.05)).collect();
    let solved = refine_window_in_rounds(
        &cohort.records,
        0..cohort.records.len(),
        None,
        &mut fits,
        MAX_TEMPLATES,
        NoReadRule::default(),
        Phasing::rounds_for(window_key(3, 1_000)),
    );
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    let mut phased = 0;
    for record in &solved.records {
        for phase in &record.phases {
            let words = match phase {
                Some(phase) => {
                    phased += 1;
                    [
                        phase.copy_0_allele as u64,
                        phase.copy_1_allele as u64,
                        phase.phase_set,
                    ]
                }
                None => [u64::MAX; 3],
            };
            for byte in words.iter().flat_map(|word| word.to_le_bytes()) {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    assert!(phased > 1_000, "{phased} phased columns");
    assert_eq!(
        (hash, phased),
        (RECORDED_PHASE_BITS, RECORDED_PHASED_COLUMNS),
        "the phase's bits are {hash:#018x} over {phased} phased columns"
    );
}

/// The hash and the count [`a_phased_windows_phase_is_the_same_bits_everywhere`] expects.
const RECORDED_PHASE_BITS: u64 = 0x2e55_296e_634e_c8a5;
const RECORDED_PHASED_COLUMNS: usize = 2_610;

/// **More rounds than the link test's masks hold are refused when the phasing is set up**, not in
/// the middle of a window: 9 rounds would draw 8 × 9 = 72 paths, more than 64; 8 rounds draw 63.
#[test]
#[should_panic(expected = "72 paths")]
fn more_rounds_than_the_link_test_holds_are_refused() {
    let _ = Phasing::new(8, 0);
    let _ = Phasing::new(9, 0);
}

/// **The draw's link is positive where the reads say cis** (reads-in-draw spec §2.2; review M2):
/// fragments that show the lower allele at both records, or the upper at both, give a positive `L`,
/// which the chain reads as the same orientation; fragments that pair a lower with an upper give a
/// negative one. At a `1/2` heterozygote the lower allele is 1, not the reference.
#[test]
fn the_draws_link_is_positive_for_cis() {
    let allele = |fragment: u32, allele: u16| FragmentAllele {
        sample: 0,
        fragment,
        allele,
        error_phred: 30,
    };
    // Two fragments, each the lower allele at both records or the upper at both: cis.
    let cis_one = [allele(0, 0), allele(1, 1)];
    let cis_other = [allele(0, 0), allele(1, 1)];
    let cis = draw_link(&cis_one, (0, 1), &cis_other, (0, 1)).expect("two fragments count");
    assert!(cis > 0.0, "{cis}");
    // The same fragments with the second record's alleles swapped: trans, the same size.
    let trans_other = [allele(0, 1), allele(1, 0)];
    let trans = draw_link(&cis_one, (0, 1), &trans_other, (0, 1)).expect("two fragments count");
    assert_eq!(trans, -cis);
    // A 1/2 heterozygote beside a 0/1: allele 1 there is the lower, so 1 with 0 is cis.
    let one_two = [allele(0, 1), allele(1, 2)];
    let zero_one = [allele(0, 0), allele(1, 1)];
    let multi = draw_link(&one_two, (1, 2), &zero_one, (0, 1)).expect("two fragments count");
    assert_eq!(multi, cis);
    // No fragment at both records: no link.
    assert_eq!(draw_link(&cis_one, (0, 1), &[allele(7, 0)], (0, 1)), None);
    assert_eq!(draw_link(&[], (0, 1), &cis_other, (0, 1)), None);
}

/// **About one in twenty of a sample's core records with fragments is held out for the read check,
/// never a hidden one, and the draw depends on the window and the sample alone** (reads-in-draw
/// spec §2.6). Sample 0 has a fragment at every record of a 3,000-record window, sample 1 at every
/// other one; records 0–99 are hidden.
#[test]
fn read_check_records_are_one_in_twenty_and_never_hidden() {
    let records = 3_000;
    let cohort = test_cohort(CohortShape {
        records,
        depth: 30.0,
        ..SHAPE
    });
    let layout = WindowLayout::new(&cohort.records, &(0..records));
    let fragments: Vec<RecordFragments> = (0..records)
        .map(|index| {
            let entry = |sample: u32| FragmentAllele {
                sample,
                fragment: index as u32,
                allele: 0,
                error_phred: 30,
            };
            let mut entries = vec![entry(0)];
            if index % 2 == 0 {
                entries.push(entry(1));
            }
            RecordFragments::from_entries(entries)
        })
        .collect();
    let hidden: Vec<usize> = (0..100).collect();
    let with_fragments = |sample: usize| {
        layout
            .core_sites
            .iter()
            .map(|&site| layout.record_of_site[site])
            .filter(|&index| sample == 0 || index % 2 == 0)
            .filter(|index| !hidden.contains(index))
            .count()
    };
    for sample in [0, 1] {
        let held_out = draw_read_check_records(&fragments, &layout, sample, &hidden, 5);
        assert!(held_out.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(held_out.iter().all(|index| !hidden.contains(index)));
        assert!(
            held_out
                .iter()
                .all(|&index| !fragments[index].of_sample(sample as u32).is_empty())
        );
        // Binomial(n, 1/20) within five standard deviations.
        let n = with_fragments(sample) as f64;
        let expected = n * READ_CHECK_SHARE;
        let spread = 5.0 * (n * READ_CHECK_SHARE * (1.0 - READ_CHECK_SHARE)).sqrt();
        assert!(
            (held_out.len() as f64 - expected).abs() <= spread,
            "sample {sample}: {} held out of {n}",
            held_out.len()
        );
        assert_eq!(
            held_out,
            draw_read_check_records(&fragments, &layout, sample, &hidden, 5)
        );
        assert_ne!(
            held_out,
            draw_read_check_records(&fragments, &layout, sample, &hidden, 6),
            "another window draws other records"
        );
    }
}
