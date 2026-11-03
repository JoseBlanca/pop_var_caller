use super::*;
use crate::types::{GenomePosition, Position};

fn place(contig: u32, position: u64) -> RecordPlace {
    RecordPlace {
        at: GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        },
        is_repeat_tract: false,
    }
}

/// A line at `position` with `FILTER` `filter` and one `GT:GQ` column per sample.
fn line(position: u64, filter: &str, columns: &[&str]) -> Vec<u8> {
    format!(
        "chr1\t{position}\t.\tA\tG,T\t50\t{filter}\t.\tGT:GQ\t{}",
        columns.join("\t")
    )
    .into_bytes()
}

/// Fragments as `(sample, fragment, allele)`, all at Phred 30.
fn fragments(entries: &[(u32, u32, u16)]) -> RecordFragments {
    RecordFragments::from_entries(
        entries
            .iter()
            .map(|&(sample, fragment, allele)| FragmentAllele {
                sample,
                fragment,
                allele,
                error_phred: 30,
            })
            .collect(),
    )
}

/// Every line through a two-sample stage; the lines that come out, as text.
fn phased(lines: Vec<(u64, Vec<u8>, RecordFragments)>) -> (Vec<String>, ReadPhaseTally) {
    let mut stage = ReadPhaseStage::new(2);
    let mut out = Vec::new();
    for (position, text, fragments) in lines {
        stage.push(place(0, position), text, fragments, &mut out);
    }
    stage.finish(&mut out);
    (
        out.into_iter()
            .map(|(_, line)| String::from_utf8(line).expect("text"))
            .collect(),
        stage.tally(),
    )
}

/// The `GT:GQ:PS` column of one sample on one written line.
fn column(written: &str, sample: usize) -> &str {
    written.split('\t').nth(9 + sample).expect("a column")
}

/// **Two heterozygotes two good fragments link are phased together** (§2.6): the first `a|b`, its
/// lower allele on copy 0; the second oriented by the reads — here trans, so its upper allele is
/// on copy 0 — and both carry `PS` the first's `POS`. The other sample, without fragments, gets
/// `PS` `.` and its `GT` unchanged; every line's `FORMAT` ends with `PS`.
#[test]
fn two_heterozygotes_linked_by_two_fragments_share_a_phase_set() {
    let (written, tally) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "0/1:40"]),
            fragments(&[(0, 1, 0), (0, 2, 1)]),
        ),
        (
            400,
            line(400, "PASS", &["0/2:40", "0/0:40"]),
            fragments(&[(0, 1, 2), (0, 2, 0)]),
        ),
    ]);
    assert_eq!(column(&written[0], 0), "0|1:40:100");
    assert_eq!(column(&written[1], 0), "2|0:40:100");
    assert_eq!(column(&written[0], 1), "0/1:40:.");
    assert_eq!(column(&written[1], 1), "0/0:40:.");
    assert!(written.iter().all(|line| line.contains("\tGT:GQ:PS\t")));
    assert_eq!(tally.columns_phased, 2);
    assert_eq!(tally.blocks, 1);
    assert_eq!(tally.links, 1);
}

/// **A record with `GQ` 29, or not `PASS`, is never read-linked** (§2.5): the fragments are the
/// same as above, and nothing is phased.
#[test]
fn a_record_below_gq_30_or_not_pass_is_never_linked() {
    for (filter, quality) in [("PASS", "29"), ("lowDepth", "40")] {
        let (written, tally) = phased(vec![
            (
                100,
                line(100, "PASS", &["0/1:40", "./.:."]),
                fragments(&[(0, 1, 0), (0, 2, 1)]),
            ),
            (
                400,
                line(400, filter, &[&format!("0/1:{quality}"), "./.:."]),
                fragments(&[(0, 1, 1), (0, 2, 0)]),
            ),
        ]);
        assert_eq!(column(&written[0], 0), "0/1:40:.", "{filter} {quality}");
        assert_eq!(column(&written[1], 0), format!("0/1:{quality}:."));
        assert_eq!(tally.columns_phased, 0);
    }
}

/// **One fragment is not a link**, and a record that links to nothing is a block of one, written
/// unphased.
#[test]
fn one_fragment_is_not_a_link_and_a_block_of_one_is_unphased() {
    let (written, tally) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0)]),
        ),
        (
            400,
            line(400, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0)]),
        ),
    ]);
    assert_eq!(column(&written[0], 0), "0/1:40:.");
    assert_eq!(column(&written[1], 0), "0/1:40:.");
    assert_eq!((tally.columns_phased, tally.blocks), (0, 0));
}

/// **An eligible record without a fragment shared with another is skipped, not a gap**: records 1
/// and 3 are linked by two fragments across record 2, which has none, and are phased together;
/// record 2 is not.
#[test]
fn an_eligible_record_without_shared_fragments_is_skipped_not_a_gap() {
    let (written, _) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0)]),
        ),
        (200, line(200, "PASS", &["0/1:40", "./.:."]), fragments(&[])),
        (
            300,
            line(300, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0)]),
        ),
    ]);
    assert_eq!(column(&written[0], 0), "0|1:40:100");
    assert_eq!(column(&written[1], 0), "0/1:40:.");
    assert_eq!(column(&written[2], 0), "0|1:40:100");
}

/// **A fragment over a non-consecutive pair ends a block it contradicts** (§2.6): records 1–2 and
/// 2–3 are each linked cis by two fragments, which makes 1–3 cis; two other fragments over 1 and 3
/// alone say trans, at −7.8 nats against it. The block ends before 3, which links to nothing else
/// and is written unphased.
#[test]
fn a_fragment_over_a_skipped_record_ends_a_block_it_contradicts() {
    let (written, tally) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0), (0, 5, 0), (0, 6, 0)]),
        ),
        (
            300,
            line(300, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0)]),
        ),
        (
            500,
            line(500, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 3, 0), (0, 4, 0), (0, 5, 1), (0, 6, 1)]),
        ),
    ]);
    assert_eq!(column(&written[0], 0), "0|1:40:100");
    assert_eq!(column(&written[1], 0), "0|1:40:100");
    assert_eq!(column(&written[2], 0), "0/1:40:.");
    assert_eq!(tally.links_refused_by_an_earlier_record, 1);
}

/// **A line is held until the stream is more than `2R` past it, lines leave in the order they
/// came, and a contig change releases everything**; a line with nothing to phase still gains `PS`.
#[test]
fn lines_are_held_two_reaches_and_leave_in_order() {
    let mut stage = ReadPhaseStage::new(2);
    let mut out = Vec::new();
    stage.push(
        place(0, 100),
        line(100, "PASS", &["0/1:40", "./.:."]),
        RecordFragments::default(),
        &mut out,
    );
    stage.push(
        place(0, 100 + 2 * READ_REACH),
        line(100 + 2 * READ_REACH, "PASS", &["0/0:40", "./.:."]),
        RecordFragments::default(),
        &mut out,
    );
    assert!(out.is_empty(), "exactly 2R past is not more than 2R");
    stage.push(
        place(0, 101 + 2 * READ_REACH),
        line(101 + 2 * READ_REACH, "PASS", &["0/0:40", "./.:."]),
        RecordFragments::default(),
        &mut out,
    );
    assert_eq!(out.len(), 1);
    stage.push(
        place(1, 5),
        line(5, "PASS", &["0/0:40", "./.:."]),
        RecordFragments::default(),
        &mut out,
    );
    assert_eq!(out.len(), 3, "a new contig releases the old one's lines");
    stage.finish(&mut out);
    let positions: Vec<(u32, u64)> = out
        .iter()
        .map(|(place, _)| (place.at.contig.0, place.at.position.get()))
        .collect();
    assert_eq!(positions, vec![(0, 100), (0, 2_100), (0, 2_101), (1, 5)],);
    assert!(
        out.iter()
            .all(|(_, line)| String::from_utf8_lossy(line).ends_with(":.\t./.:.:."))
    );
}

/// **A column that dropped trailing fields is padded, so `PS` lands in its own place.**
#[test]
fn a_short_column_is_padded_before_ps() {
    let written = with_phase(
        b"chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:7\t./.",
        &[None, None],
    );
    assert_eq!(
        String::from_utf8(written).expect("text"),
        "chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:PS\t0/1:40:7:.\t./.:.:.:."
    );
}
