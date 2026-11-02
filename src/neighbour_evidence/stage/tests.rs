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
         hidden_share=0.05,extra_paths=8,lean_min=0.5,path_agreement=1>"
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
            stage.push(*place, line.clone().into_bytes(), &mut output);
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

/// **Blocks join across windows, and the phase stays right across them** (phase output spec
/// §2.4, §7.3): an outbred cohort at 4× on a 1-Mb contig, windows of 200 kb with 200-kb margins.
/// Every `PS` is the position of its block's first record; some blocks run across a window's
/// boundary; and the links inside blocks — across boundaries too — are right: when written, 9
/// wrong in 5,494, and 29 links crossed a boundary, none wrong.
#[test]
fn blocks_join_across_windows() {
    let cohort = test_cohort(CohortShape {
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
    });
    let mut decision = decide(NeighbourEvidenceSwitch::On, &[0.05; 20], diploid());
    decision.sizes = WindowSizes {
        core_bases: 200_000,
        margin_bases: 200_000,
        margin_floor_bases: 200_000,
        budget_sample_records: 120_000_000,
    };
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![1_000_000]);
    let mut output = Vec::new();
    for (index, line) in cohort.lines.iter().enumerate() {
        stage.push(
            place(0, 1 + index as u64 * 500),
            line.clone().into_bytes(),
            &mut output,
        );
    }
    stage.finish(&mut output);

    let (mut links, mut flips, mut crossing, mut crossing_flips) = (0, 0, 0, 0);
    for sample in 0..20 {
        // (position, right, PS) of every phased column of the sample.
        let mut phased: Vec<(u64, bool, u64)> = Vec::new();
        for (place, line) in &output {
            let line = String::from_utf8_lossy(line);
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(fields[8].ends_with(":PS"), "{}", fields[8]);
            let column: Vec<&str> = fields[9 + sample].split(':').collect();
            let Some(genotype) = column[0].split_once('|') else {
                assert_eq!(*column.last().expect("PS"), ".");
                continue;
            };
            let position = place.at.position.0;
            let [copy_0, _] = cohort.truth_copies[((position - 1) / 500) as usize][sample];
            let block: u64 = column.last().expect("PS").parse().expect("a PS");
            phased.push((position, genotype.0 == copy_0.to_string(), block));
        }
        let mut first_of: std::collections::HashMap<u64, u64> = Default::default();
        for &(position, _, block) in &phased {
            first_of.entry(block).or_insert(position);
        }
        assert!(
            first_of.iter().all(|(block, first)| block == first),
            "every PS is its block's first position"
        );
        for pair in phased.windows(2) {
            let ((before, was_right, before_block), (after, right, block)) = (pair[0], pair[1]);
            if before_block == block {
                links += 1;
                flips += usize::from(was_right != right);
                // Window cores start at 1 and every 200 kb after.
                if (before - 1) / 200_000 != (after - 1) / 200_000 {
                    crossing += 1;
                    crossing_flips += usize::from(was_right != right);
                }
            }
        }
    }
    eprintln!(
        "{flips} flips in {links} links; {crossing} links cross a window boundary, \
         {crossing_flips} of them flips"
    );
    assert!(crossing > 0, "some block runs across a window's boundary");
    assert!(flips * 100 <= links, "{flips} flips in {links} links");
}
