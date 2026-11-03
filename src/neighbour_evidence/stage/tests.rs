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
         core=2000000,margin=1000000,rounds=5,haplotype_templates=16,hidden_per_window=50,\
         hidden_share=0.05,extra_paths=8,lean_min=0.5,path_agreement=1,join_records=3,\
         alignment_costs=1/2>"
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
        budget_sample_records: 120_000_000,
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
            stage.push(
                *place,
                line.clone().into_bytes(),
                RecordFragments::default(),
                &mut output,
            );
        }
        stage.finish(&mut output);
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
            stage.push(
                place(0, at),
                line.into_bytes(),
                RecordFragments::default(),
                &mut output,
            );
        }
        stage.finish(&mut output);
        output
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
    budget_sample_records: 120_000_000,
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
        budget_sample_records: 9_000,
    };
    // The cores the planner writes, as the stage plans them.
    let mut core_starts = Vec::new();
    let mut start = 0;
    while start < positions.len() {
        let plan =
            crate::neighbour_evidence::windows::plan_window(&positions, start, true, 20, &sizes)
                .expect("a whole contig plans");
        core_starts.push(positions[plan.core_start]);
        start = plan.core_end;
    }
    let output = join_stage(&cohort, position, sizes, 4);
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
        stage.push(
            place(0, position),
            line.clone().into_bytes(),
            RecordFragments::default(),
            &mut output,
        );
    }
    stage.finish(&mut output);
    let kept = output.len();
    assert!(
        kept * 100 >= carried * 95 && kept * 100 <= carried * 105,
        "{} records in, {kept} out; {carried} carry the alternative",
        cohort.lines.len()
    );
}
