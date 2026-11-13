use super::*;
use crate::neighbour_evidence::test_cohort::{CohortShape, test_cohort};
use crate::types::{GenomePosition, Position};

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a ploidy")
}

/// **`auto` runs on a diploid cohort of at least ten samples with median inbreeding at least
/// 0.5; `on` wherever it can; `off` never** — and every decision says why (§5).
#[test]
fn the_switch_and_the_cohort_decide() {
    let inbred = [0.9; 12];
    let decision = decide(NeighbourEvidenceSwitch::Auto, &inbred, diploid());
    assert!(decision.ran);
    assert_eq!(decision.why, "12 samples, median inbreeding 0.90");
    assert!(!decide(NeighbourEvidenceSwitch::Auto, &inbred[..9], diploid()).ran);
    assert!(!decide(NeighbourEvidenceSwitch::Auto, &[0.1; 12], diploid()).ran);
    assert!(decide(NeighbourEvidenceSwitch::On, &[0.1; 3], diploid()).ran);
    assert!(!decide(NeighbourEvidenceSwitch::On, &[0.9], diploid()).ran);
    assert!(!decide(NeighbourEvidenceSwitch::Off, &inbred, diploid()).ran);
    let tetraploid = Ploidy::try_new(4).expect("a ploidy");
    let decision = decide(NeighbourEvidenceSwitch::On, &inbred, tetraploid);
    assert!(!decision.ran);
    assert_eq!(decision.why, "ploidy 4; the pass models two copies");
    // The median of an even count is the mean of the middle two: 0.45 here.
    assert!(
        !decide(
            NeighbourEvidenceSwitch::Auto,
            &[[0.0; 5], [0.9; 5]].concat(),
            diploid()
        )
        .ran
    );
}

/// **The header line names the decision, and its parameters only when the pass ran** (§4.1).
#[test]
fn the_header_line_states_the_decision() {
    let decision = decide(NeighbourEvidenceSwitch::Auto, &[0.98; 57], diploid());
    assert_eq!(
        decision.as_header_value(),
        "<ran=yes,why=\"57 samples, median inbreeding 0.98\",templates=16,mismatch=0.01,\
         core=2000000,margin=1000000,window_bytes=32000000000,rounds=5,haplotype_templates=16,hidden_per_window=50,\
         hidden_share=0.05,extra_paths=8,lean_min=0.5,path_agreement=0.8,join_records=3,\
         alignment_costs=1/2,read_link_min=5.29,read_break_min=2.20,read_error_floor=0.01,\
         read_gq_min=30,read_reach=1000,read_cap=32>"
    );
    let decision = decide(NeighbourEvidenceSwitch::Auto, &[0.0; 2], diploid());
    assert_eq!(
        decision.as_header_value(),
        "<ran=no,why=\"2 samples, median inbreeding 0.00; fewer than 10 samples\">"
    );
}

fn place(contig: u32, position: u64) -> RecordPlace {
    RecordPlace {
        at: GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        },
        is_repeat_tract: false,
    }
}

/// The drawn cohort's lines on two contigs with a non-PASS line among them, run through a stage
/// with windows small enough that each contig holds several, on `threads` threads.
/// The lines put in, the lines handed back with their places, and the pass's counts.
type StageRun = (Vec<String>, Vec<(RecordPlace, Vec<u8>)>, RewriteTally);

fn run_stage(threads: usize) -> StageRun {
    let cohort = test_cohort(CohortShape {
        samples: 20,
        founders: 4,
        records: 2_000,
        spacing: 500,
        segment_bases: 300_000.0,
        inbred: true,
        depth: 1.0,
        read_error: 0.01,
        seed: 9,
        silent_from: None,
    });
    let mut decision = decide(NeighbourEvidenceSwitch::On, &[0.95; 20], diploid());
    decision.sizes = WindowSizes {
        core_bases: 300_000,
        margin_bases: 300_000,
        margin_floor_bases: 300_000,
        budget_bytes: crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES,
    };
    let mut inputs: Vec<(RecordPlace, String)> = Vec::new();
    for contig in 0..2 {
        for (index, line) in cohort.lines.iter().enumerate() {
            let line = line.replace("chr1", &format!("chr{}", contig + 1));
            let position = 1 + index as u64 * 500;
            let line = if index == 700 {
                line.replace("PASS", "lowDepth")
            } else {
                line
            };
            inputs.push((place(contig, position), line));
        }
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a pool");
    let (output, tally) = pool.install(|| {
        let mut stage = NeighbourEvidenceStage::new(&decision, vec![1_000_000, 1_000_000]);
        let mut output = Vec::new();
        for (place, line) in &inputs {
            stage
                .push(
                    *place,
                    line.clone().into_bytes(),
                    RecordFragments::default(),
                    &mut output,
                )
                .expect("within the window budget");
        }
        stage.finish(&mut output).expect("within the window budget");
        (output, stage.tally())
    });
    (
        inputs.into_iter().map(|(_, line)| line).collect(),
        output,
        tally,
    )
}

/// **Every line comes out once, in the order it went in, unless the pass dropped it**; a
/// non-PASS line comes out byte for byte; most samples without reads are called; and one thread
/// and four write the same bytes (§8.7).
#[test]
fn every_line_comes_out_once_in_order() {
    let (inputs, output, tally) = run_stage(4);
    assert_eq!(
        output.len() as u64 + tally.records_dropped,
        inputs.len() as u64
    );
    let positions: Vec<(u32, u64)> = output
        .iter()
        .map(|(place, _)| (place.at.contig.0, place.at.position.0))
        .collect();
    assert!(
        positions.windows(2).all(|pair| pair[0] < pair[1]),
        "in order, each once"
    );
    for contig in [0, 1] {
        let unchanged = &output
            .iter()
            .find(|(place, _)| place.at.contig.0 == contig && place.at.position.0 == 1 + 700 * 500)
            .expect("the non-PASS line is written")
            .1;
        assert_eq!(
            String::from_utf8_lossy(unchanged),
            inputs[contig as usize * 2_000 + 700]
        );
    }
    let no_calls_before: usize = inputs.iter().map(|line| line.matches("./.").count()).sum();
    let no_calls_after: usize = output
        .iter()
        .map(|(_, line)| String::from_utf8_lossy(line).matches("./.").count())
        .sum();
    assert!(
        no_calls_after * 2 < no_calls_before,
        "{no_calls_before} no-calls before the pass, {no_calls_after} after"
    );
    assert_eq!(tally.sites, 2 * 1_999);
    let (_, one_thread, _) = run_stage(1);
    assert_eq!(output, one_thread);
}

/// The outbred cohort the join tests run: 20 samples of 4 founders at 4×, 2,000 records.
fn join_cohort() -> crate::neighbour_evidence::test_cohort::TestCohort {
    test_cohort(CohortShape {
        samples: 20,
        founders: 4,
        records: 2_000,
        spacing: 500,
        segment_bases: 300_000.0,
        inbred: false,
        depth: 4.0,
        read_error: 0.01,
        seed: 13,
        silent_from: None,
    })
}

/// **The join cohort through a stage**, record `index` at `position(index)` on one contig long
/// enough for all, windows of `sizes`, on `threads` threads; the lines handed back.
fn join_stage(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    position: impl Fn(usize) -> u64 + Sync,
    sizes: WindowSizes,
    threads: usize,
) -> Vec<(RecordPlace, Vec<u8>)> {
    join_stage_and_cores(cohort, position, sizes, threads).0
}

/// [`join_stage`], and the first position of every core the stage solved.
fn join_stage_and_cores(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    position: impl Fn(usize) -> u64 + Sync,
    sizes: WindowSizes,
    threads: usize,
) -> (Vec<(RecordPlace, Vec<u8>)>, Vec<u64>) {
    let mut decision = decide(NeighbourEvidenceSwitch::On, &[0.05; 20], diploid());
    decision.sizes = sizes;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a pool");
    pool.install(|| {
        let mut stage = NeighbourEvidenceStage::new(&decision, vec![2_000_000]);
        let mut output = Vec::new();
        for (index, line) in cohort.lines.iter().enumerate() {
            let at = position(index);
            // The drawn line's POS is `1 + 500·index`; it is moved to `at`.
            let line = line.replacen(&format!("\t{}\t", 1 + index * 500), &format!("\t{at}\t"), 1);
            stage
                .push(
                    place(0, at),
                    line.into_bytes(),
                    RecordFragments::default(),
                    &mut output,
                )
                .expect("within the window budget");
        }
        stage.finish(&mut output).expect("within the window budget");
        (output, stage.core_starts.clone())
    })
}

/// What the join tests count over every sample's phased columns.
#[derive(Debug, Default)]
struct JoinScore {
    links: usize,
    flips: usize,
    /// Links between two records in different window cores, and of them the wrong ones.
    crossing: usize,
    crossing_flips: usize,
    /// The most window cores one phase set reaches into.
    most_cores_spanned: usize,
}

/// **Scores a join stage's phase against the truth**: every column ends in `PS`, every `PS` is its
/// set's first position, and the links are counted; `index_of(position)` is the drawn record at a
/// position, `core_of(position)` the window core it is written in.
fn score_joins(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    output: &[(RecordPlace, Vec<u8>)],
    index_of: impl Fn(u64) -> usize,
    core_of: impl Fn(u64) -> usize,
) -> JoinScore {
    let mut score = JoinScore::default();
    for sample in 0..20 {
        // (position, right, PS) of every phased column of the sample.
        let mut phased: Vec<(u64, bool, u64)> = Vec::new();
        for (place, line) in output {
            let line = String::from_utf8_lossy(line);
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(fields[8].ends_with(":PS"), "{}", fields[8]);
            let column: Vec<&str> = fields[9 + sample].split(':').collect();
            let Some(genotype) = column[0].split_once('|') else {
                assert_eq!(*column.last().expect("PS"), ".");
                continue;
            };
            let position = place.at.position.0;
            let [copy_0, _] = cohort.truth_copies[index_of(position)][sample];
            let phase_set: u64 = column.last().expect("PS").parse().expect("a PS");
            phased.push((position, genotype.0 == copy_0.to_string(), phase_set));
        }
        let mut spans: std::collections::BTreeMap<u64, (u64, usize, usize)> = Default::default();
        for &(position, _, phase_set) in &phased {
            let span =
                spans
                    .entry(phase_set)
                    .or_insert((position, core_of(position), core_of(position)));
            span.2 = core_of(position);
        }
        for (phase_set, (first, first_core, last_core)) in &spans {
            assert_eq!(phase_set, first, "a PS is its set's first position");
            score.most_cores_spanned = score.most_cores_spanned.max(last_core - first_core + 1);
        }
        for pair in phased.windows(2) {
            let ((before, was_right, before_set), (after, right, phase_set)) = (pair[0], pair[1]);
            if before_set == phase_set {
                score.links += 1;
                score.flips += usize::from(was_right != right);
                if core_of(before) != core_of(after) {
                    score.crossing += 1;
                    score.crossing_flips += usize::from(was_right != right);
                }
            }
        }
    }
    eprintln!("{score:?}");
    score
}

const JOIN_WINDOWS: WindowSizes = WindowSizes {
    core_bases: 200_000,
    margin_bases: 200_000,
    margin_floor_bases: 200_000,
    budget_bytes: crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES,
};

/// **Blocks join across windows, and the phase stays right across them** (phase output spec
/// §2.4, §7.3): the join cohort on 1 Mb, windows of 200 kb with 200-kb margins. Some phase sets
/// run across a window's boundary, and the links that cross one are right — a join written with
/// its block reversed turns every one of them (code review 2026-11-02: 29 of 29). When written,
/// 9 links wrong in 5,494, and 29 crossed a boundary, none wrong.
#[test]
fn blocks_join_across_windows() {
    let cohort = join_cohort();
    let output = join_stage(&cohort, |index| 1 + index as u64 * 500, JOIN_WINDOWS, 4);
    let score = score_joins(
        &cohort,
        &output,
        |position| ((position - 1) / 500) as usize,
        |position| ((position - 1) / 200_000) as usize,
    );
    assert!(score.crossing >= 10, "{score:?}");
    assert!(score.crossing_flips * 10 <= score.crossing, "{score:?}");
    assert!(score.flips * 100 <= score.links, "{score:?}");
}

/// **The phase is the same at one, four and eight threads** across windows that join (spec
/// §7.5; code review 2026-11-02: the older one-against-four test writes no phase).
#[test]
fn the_phase_is_the_same_at_any_thread_count() {
    let cohort = join_cohort();
    let at = |threads| {
        join_stage(
            &cohort,
            |index| 1 + index as u64 * 500,
            JOIN_WINDOWS,
            threads,
        )
    };
    let one = at(1);
    assert!(
        one.iter().any(|(_, line)| line.contains(&b'|')),
        "the run writes phase"
    );
    assert_eq!(one, at(4));
    assert_eq!(one, at(8));
}

/// **A phase set carries through short windows** (spec §2.4, §7.3): records every 500 bases with
/// a dense stretch of 400 every 100 bases, and a budget that shrinks the windows over it — cores
/// shorter than the margin after them, which then reaches back into the core before. The last
/// written block must carry its records from window to window, and the links that cross a core
/// boundary stay right: when written, 17 cores, 95 crossing links none wrong, and a phase set
/// reaching into 3 cores.
#[test]
fn a_phase_set_carries_through_short_windows() {
    let cohort = join_cohort();
    // Records 800 to 1,199 every 100 bases, the rest every 500.
    let position = |index: usize| -> u64 {
        let index = index as u64;
        match index {
            0..800 => 1 + index * 500,
            800..1_200 => 400_001 + (index - 800) * 100,
            _ => 440_001 + (index - 1_200) * 500,
        }
    };
    let positions: Vec<u64> = (0..2_000).map(position).collect();
    let sizes = WindowSizes {
        core_bases: 60_000,
        margin_bases: 60_000,
        margin_floor_bases: 20_000,
        // Small enough that the dense records shrink the windows, large enough that the floor
        // fits them.
        budget_bytes: 2_500_000,
    };
    let (output, core_starts) = join_stage_and_cores(&cohort, position, sizes, 4);
    let score = score_joins(
        &cohort,
        &output,
        |at| positions.binary_search(&at).expect("a drawn position"),
        |at| core_starts.partition_point(|&start| start <= at) - 1,
    );
    eprintln!("cores start at {core_starts:?}");
    assert!(core_starts.len() >= 6, "{core_starts:?}");
    assert!(score.most_cores_spanned >= 3, "{score:?}");
    assert!(score.crossing_flips * 10 <= score.crossing, "{score:?}");
    // Margins of 20 to 60 kb hold far less than a copying length; the phase inside them is
    // poorer than in full-sized windows (59 wrong in 5,693 when written).
    assert!(score.flips * 50 <= score.links, "{score:?}");
}

/// **The records the pass drops are those no sample carries the alternative at**: the drawn cohort
/// writes a record at every drawn site, and at most of the rare ones no founder carries the
/// alternative, so the pass drops them as calling would have. Moved here from the removed
/// `refine-genotypes` command's tests (read phase spec §5).
#[test]
fn the_records_kept_are_those_some_sample_carries() {
    let cohort = test_cohort(CohortShape {
        samples: 20,
        founders: 4,
        records: 1_500,
        spacing: 600,
        segment_bases: 300_000.0,
        inbred: true,
        depth: 1.0,
        read_error: 0.01,
        seed: 13,
        silent_from: None,
    });
    let carried = cohort
        .truth
        .iter()
        .filter(|genotypes| genotypes.iter().any(|&genotype| genotype > 0))
        .count();
    let decision = decide(NeighbourEvidenceSwitch::On, &[0.95; 20], diploid());
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![1_000_000]);
    let mut output = Vec::new();
    for (index, line) in cohort.lines.iter().enumerate() {
        let position = 1 + index as u64 * 600;
        stage
            .push(
                place(0, position),
                line.clone().into_bytes(),
                RecordFragments::default(),
                &mut output,
            )
            .expect("within the window budget");
    }
    stage.finish(&mut output).expect("within the window budget");
    let kept = output.len();
    assert!(
        kept * 100 >= carried * 95 && kept * 100 <= carried * 105,
        "{} records in, {kept} out; {carried} carry the alternative",
        cohort.lines.len()
    );
}

/// The join cohort through a stage in one window, `fragments(index)` beside record `index`.
fn join_stage_with_reads(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    fragments: impl Fn(usize) -> RecordFragments,
) -> Vec<(RecordPlace, Vec<u8>)> {
    join_stage_drawing(cohort, fragments, true)
}

/// [`join_stage_with_reads`], the fragments orienting the rounds' draws only if `reads_in_draw`.
fn join_stage_drawing(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
    fragments: impl Fn(usize) -> RecordFragments,
    reads_in_draw: bool,
) -> Vec<(RecordPlace, Vec<u8>)> {
    let decision = NeighbourEvidenceDecision {
        reads_in_draw,
        ..decide(NeighbourEvidenceSwitch::On, &[0.05; 20], diploid())
    };
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![2_000_000]);
    let mut output = Vec::new();
    for (index, line) in cohort.lines.iter().enumerate() {
        stage
            .push(
                place(0, 1 + index as u64 * 500),
                line.clone().into_bytes(),
                fragments(index),
                &mut output,
            )
            .expect("within the window budget");
    }
    stage.finish(&mut output).expect("within the window budget");
    output
}

/// One sample's written `GT` and `PS` at each output line, by the drawn record's index.
fn phase_of(output: &[(RecordPlace, Vec<u8>)], sample: usize) -> Vec<(usize, String, String)> {
    output
        .iter()
        .filter_map(|(place, line)| {
            let text = String::from_utf8_lossy(line).into_owned();
            let columns: Vec<&str> = text.split('\t').collect();
            let keys: Vec<&str> = columns[8].split(':').collect();
            let values: Vec<&str> = columns[9 + sample].split(':').collect();
            let value = |key: &str| {
                keys.iter()
                    .position(|name| *name == key)
                    .and_then(|at| values.get(at))
                    .map(|value| (*value).to_owned())
            };
            Some((
                ((place.at.position.get() - 1) / 500) as usize,
                value("GT")?,
                value("PS")?,
            ))
        })
        .collect()
}

/// **The reads a line carries reach the pass's link test** (read phase spec §2.7): two records a
/// sample's population phase links are given two fragments that say the opposite relation, and the
/// block breaks between them — the fragments having gone from the stage's held lines, by site, to
/// the window's solve and the sample's link test. Without them the two are written in one block.
#[test]
fn fragments_a_line_carries_reach_the_pass_and_break_a_link_they_contradict() {
    let cohort = join_cohort();
    let without = join_stage_with_reads(&cohort, |_| RecordFragments::default());
    // A sample's first two consecutive records written in one phase set.
    let (sample, first, second) = (0..20)
        .find_map(|sample| {
            let phased: Vec<(usize, String, String)> = phase_of(&without, sample)
                .into_iter()
                .filter(|(_, genotype, _)| genotype.contains('|'))
                .collect();
            phased
                .windows(2)
                .find(|pair| pair[0].2 == pair[1].2)
                .map(|pair| (sample, pair[0].clone(), pair[1].clone()))
        })
        .expect("the population phases some sample's records together");
    let alleles = |genotype: &str| -> (u16, u16) {
        let (copy_0, copy_1) = genotype.split_once('|').expect("phased");
        (
            copy_0.parse().expect("an allele"),
            copy_1.parse().expect("an allele"),
        )
    };
    let ((first_0, first_1), (second_0, second_1)) = (alleles(&first.1), alleles(&second.1));
    // Two fragments, each pairing copy 0's allele at one record with copy 1's at the other.
    let fragments = |index: usize| {
        let entries: &[(u32, u16)] = if index == first.0 {
            &[(1, first_0), (2, first_1)]
        } else if index == second.0 {
            &[(1, second_1), (2, second_0)]
        } else {
            &[]
        };
        RecordFragments::from_entries(
            entries
                .iter()
                .map(
                    |&(fragment, allele)| crate::read_phase::fragments::FragmentAllele {
                        sample: sample as u32,
                        fragment,
                        allele,
                        error_phred: 30,
                    },
                )
                .collect(),
        )
    };
    let with = join_stage_with_reads(&cohort, fragments);
    let after = phase_of(&with, sample);
    let at = |index: usize| {
        after
            .iter()
            .find(|(record, _, _)| *record == index)
            .map(|(_, _, phase_set)| phase_set.clone())
    };
    let (first_after, second_after) = (at(first.0), at(second.0));
    assert!(
        first_after != second_after || first_after.as_deref() == Some("."),
        "sample {sample}: records {} and {} must no longer be phased together, and are in {:?}",
        first.0,
        second.0,
        first_after
    );
}

/// Fragments over every heterozygous column of the join cohort: fragments `2k` and `2k + 1` cover
/// records `k` and `k + 1`; about one in five shows the upper allele, with error Phreds from 10 to
/// 49, so links of both signs and of every strength reach the link test.
fn mixed_fragments(
    cohort: &crate::neighbour_evidence::test_cohort::TestCohort,
) -> impl Fn(usize) -> RecordFragments + Sync {
    let heterozygotes: Vec<Vec<Option<(u16, u16)>>> = cohort
        .lines
        .iter()
        .map(|line| {
            line.split('\t')
                .skip(9)
                .map(|column| {
                    let (x, y) =
                        crate::vcf::line_fields::parse_genotype(column.split(':').next()?)?;
                    (x != y).then_some((x as u16, y as u16))
                })
                .collect()
        })
        .collect();
    move |index: usize| {
        let mut entries = Vec::new();
        for (sample, heterozygote) in heterozygotes[index].iter().enumerate() {
            let Some((lower, upper)) = *heterozygote else {
                continue;
            };
            for k in [Some(index), index.checked_sub(1)].into_iter().flatten() {
                for fragment in [2 * k, 2 * k + 1] {
                    let allele = if (sample + index + fragment).is_multiple_of(5) {
                        upper
                    } else {
                        lower
                    };
                    entries.push(crate::read_phase::fragments::FragmentAllele {
                        sample: sample as u32,
                        fragment: fragment as u32,
                        allele,
                        error_phred: (10 + (fragment * 7 + sample) % 40) as u8,
                    });
                }
            }
        }
        RecordFragments::from_entries(entries)
    }
}

/// **With reads, the pass writes the same bytes at any thread count** (review Mi5): every
/// heterozygous column carries fragments of both signs, the stage runs in pools of one and eight
/// threads, and the output is identical — and differs from the output without reads, and from the
/// output with reads kept out of the rounds' draws, so both read paths are what ran.
#[test]
fn the_pass_with_reads_writes_the_same_bytes_at_any_thread_count() {
    let cohort = join_cohort();
    let fragments = mixed_fragments(&cohort);
    let at = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool")
            .install(|| join_stage_with_reads(&cohort, &fragments))
    };
    let one = at(1);
    assert_ne!(
        one,
        join_stage_with_reads(&cohort, |_| RecordFragments::default()),
        "the fragments reach the link test"
    );
    assert_ne!(
        one,
        join_stage_drawing(&cohort, &fragments, false),
        "the fragments reach the rounds' draws (reads-in-draw spec §6)"
    );
    assert_eq!(one, at(8));
}

/// **A column reads newly phase has `GQ` 30 or more** (read phase spec §2.5; review r1): the link
/// test's `GQ` comes from the refined posterior, which the writer's `GQ` can only raise. With
/// fragments over every heterozygous column, every column written phased with the reads and
/// unphased without them is at `GQ` 30 or more. The cohort has phased-by-population heterozygotes
/// below 30 next to ones the reads link, so a link test that ignored `GQ` would phase some. The
/// fragments stay out of the rounds' draws: there they change the population's phase, which may
/// phase a heterozygote below 30 (reads-in-draw spec §2.1), and this test is of the link test.
#[test]
fn a_column_reads_newly_phase_has_gq_30_or_more() {
    let cohort = join_cohort();
    let without = join_stage_with_reads(&cohort, |_| RecordFragments::default());
    let with = join_stage_drawing(&cohort, mixed_fragments(&cohort), false);
    // The cohort's lines carry `GP` and no `GQ` key, so the written `GQ` is read from `GP` as the
    // pass computes it. The written `GP` is renormalised, which can only raise it, so this is at
    // least the `GQ` the link test saw.
    let genotype_quality = |output: &[(RecordPlace, Vec<u8>)], sample: usize, index: usize| {
        output.iter().find_map(|(place, line)| {
            (((place.at.position.get() - 1) / 500) as usize == index).then(|| {
                let text = String::from_utf8_lossy(line).into_owned();
                let columns: Vec<&str> = text.split('\t').collect();
                let at = columns[8].split(':').position(|key| key == "GP")?;
                let posterior: Vec<f64> = columns[9 + sample]
                    .split(':')
                    .nth(at)?
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<_, _>>()
                    .ok()?;
                // Renormalised, as the written values are rounded.
                let total: f64 = posterior.iter().sum();
                let posterior: Vec<f64> = posterior.iter().map(|value| value / total).collect();
                let (_, quality) = crate::calling::quality::score_best_genotype(&posterior);
                Some(crate::vcf::encode::written_genotype_quality(quality))
            })?
        })
    };
    let (mut newly_phased, mut low_quality_heterozygotes) = (0, 0);
    for sample in 0..20 {
        let before = phase_of(&without, sample);
        for (index, genotype, _) in phase_of(&with, sample) {
            let quality = genotype_quality(&with, sample, index);
            let heterozygote = genotype
                .split(['/', '|'])
                .collect::<Vec<_>>()
                .windows(2)
                .any(|pair| pair[0] != pair[1]);
            if heterozygote && quality.is_some_and(|quality| quality < 30) {
                low_quality_heterozygotes += 1;
            }
            let was_phased = before
                .iter()
                .any(|(at, earlier, _)| *at == index && earlier.contains('|'));
            if genotype.contains('|') && !was_phased {
                newly_phased += 1;
                assert!(
                    quality.is_some_and(|quality| quality >= 30),
                    "sample {sample}, record {index}: phased by reads at GQ {quality:?}"
                );
            }
        }
    }
    assert!(
        newly_phased > 0,
        "the reads phase something the population did not"
    );
    assert!(
        low_quality_heterozygotes > 0,
        "the cohort has heterozygotes below GQ 30"
    );
}

/// **A cohort whose smallest window does not fit the budget stops the stage with an error**
/// saying what the window needs (owner, 2026-11-13): no line of that window is handed back, and the
/// error names its records, its samples and its bytes, over the budget.
#[test]
fn a_window_over_budget_at_the_floor_stops_the_stage() {
    let cohort = join_cohort();
    let mut decision = decide(NeighbourEvidenceSwitch::On, &[0.05; 20], diploid());
    decision.sizes = WindowSizes {
        core_bases: 200_000,
        margin_bases: 200_000,
        margin_floor_bases: 200_000,
        budget_bytes: 100_000,
    };
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![2_000_000]);
    let mut output = Vec::new();
    let mut refused = None;
    for (index, line) in cohort.lines.iter().enumerate() {
        if let Err(over) = stage.push(
            place(0, 1 + index as u64 * 500),
            line.clone().into_bytes(),
            RecordFragments::default(),
            &mut output,
        ) {
            refused = Some(over);
            break;
        }
    }
    let over = refused.expect("the first window is over budget");
    assert!(output.is_empty());
    assert_eq!((over.contig, over.position, over.samples), (0, 1, 20));
    assert!(
        over.records > 0 && over.bytes > over.budget && over.budget == 100_000,
        "{over:?}"
    );
}
