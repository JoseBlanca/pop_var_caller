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

/// **An eligible record whose fragments reach no other eligible record is skipped, not a gap**:
/// records 1 and 3 are linked by two fragments across record 2, whose one fragment is seen
/// nowhere else, and are phased together; record 2 is not. Kills "phaseable when it has any
/// fragment" in place of "when a fragment of it is at two eligible records" (review M4).
#[test]
fn an_eligible_record_whose_fragments_reach_no_other_is_skipped_not_a_gap() {
    let (written, _) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0)]),
        ),
        (
            200,
            line(200, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 9, 0)]),
        ),
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

/// **A fragment is remembered until every line it reaches has settled** (review M3): lines at 100
/// and 700 share two fragments; lines at 1,500 and 1,750 push the stream past both, so both settle
/// inside `push`, not in `finish`. The line at 700 settles when 1,750 arrives, 1,650 bases after
/// the fragments were first seen: forgotten after `R` instead of `2R`, they would no longer make
/// it phaseable, and the line at 100 would be a block of one.
#[test]
fn a_fragment_is_remembered_until_every_line_it_reaches_has_settled() {
    let (written, _) = phased(vec![
        (
            100,
            line(100, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0)]),
        ),
        (
            700,
            line(700, "PASS", &["0/1:40", "./.:."]),
            fragments(&[(0, 1, 0), (0, 2, 0)]),
        ),
        (
            1_500,
            line(1_500, "PASS", &["0/0:40", "./.:."]),
            fragments(&[]),
        ),
        (
            1_750,
            line(1_750, "PASS", &["0/0:40", "./.:."]),
            fragments(&[]),
        ),
    ]);
    assert_eq!(column(&written[0], 0), "0|1:40:100");
    assert_eq!(column(&written[1], 0), "0|1:40:100");
}

/// **A block that grows across `2R` is written whole, and its first line leaves phased before the
/// last arrives** (review M3): lines at 100, 900, 1,700 and 2,500, each linked to the one before
/// by its own two fragments. The line at 100 is released when 2,500 arrives, by which time 900 has
/// joined it; every line shares `PS` 100.
#[test]
fn a_block_that_grows_across_two_reaches_is_written_whole() {
    let mut stage = ReadPhaseStage::new(2);
    let mut out = Vec::new();
    let linked = [
        (100, [1, 2], None),
        (900, [3, 4], Some([1, 2])),
        (1_700, [5, 6], Some([3, 4])),
    ];
    for (position, next, previous) in linked {
        let mut entries: Vec<(u32, u32, u16)> = next.iter().map(|&f| (0, f, 0)).collect();
        entries.extend(previous.iter().flatten().map(|&f| (0, f, 0)));
        stage.push(
            place(0, position),
            line(position, "PASS", &["0/1:40", "./.:."]),
            fragments(&entries),
            &mut out,
        );
    }
    assert!(out.is_empty());
    stage.push(
        place(0, 2_500),
        line(2_500, "PASS", &["0/1:40", "./.:."]),
        fragments(&[(0, 5, 0), (0, 6, 0)]),
        &mut out,
    );
    assert_eq!(
        out.len(),
        1,
        "the line at 100 leaves when the stream is past 2,100"
    );
    stage.finish(&mut out);
    let written: Vec<String> = out
        .into_iter()
        .map(|(_, line)| String::from_utf8(line).expect("text"))
        .collect();
    for line in &written {
        assert_eq!(column(line, 0), "0|1:40:100", "{line}");
    }
    assert_eq!(stage.tally().blocks, 1);
    assert_eq!(stage.tally().links, 3);
}

/// **A contig change forgets every fragment count** (review Mi22): fragment numbers restart on
/// every contig, so fragment 1 on the second contig is a different read from fragment 1 on the
/// first. On contig 0, a heterozygote carries fragment 1. On contig 1, heterozygotes at 100 and
/// 300 share fragments 5 and 6, and the one at 200 between them has fragment 1, seen nowhere else
/// on contig 1, so it is skipped and 100 and 300 are phased together. Had contig 0's count of
/// fragment 1 carried over, the record at 200 would count as phaseable, link to neither
/// neighbour, and split them.
#[test]
fn a_contig_change_forgets_fragment_counts() {
    let mut stage = ReadPhaseStage::new(2);
    let mut out = Vec::new();
    for (contig, position, entries) in [
        (0, 100, vec![(0, 1, 0)]),
        (1, 100, vec![(0, 5, 0), (0, 6, 0)]),
        (1, 200, vec![(0, 1, 0)]),
        (1, 300, vec![(0, 5, 0), (0, 6, 0)]),
    ] {
        stage.push(
            place(contig, position),
            line(position, "PASS", &["0/1:40", "./.:."]),
            fragments(&entries),
            &mut out,
        );
    }
    stage.finish(&mut out);
    let written: Vec<String> = out
        .into_iter()
        .map(|(_, line)| String::from_utf8(line).expect("text"))
        .collect();
    assert_eq!(column(&written[1], 0), "0|1:40:100");
    assert_eq!(column(&written[2], 0), "0/1:40:.");
    assert_eq!(column(&written[3], 0), "0|1:40:100");
}

/// **Eligibility reads `GT` and `GQ` wherever `FORMAT` puts them**: a haploid `GT`, a no-call, a
/// missing `GQ`, `GQ` 29 and a homozygote are not eligible; alleles of two digits are.
#[test]
fn eligibility_reads_haploid_no_call_and_missing_gq_as_not_eligible() {
    let columns = [
        "1:40", "./.:40", "0/1:.", "0/1:29", "1/1:40", "12/10:40", "0/1",
    ];
    let eligible = eligible_columns(&line(5, "PASS", &columns), columns.len());
    assert_eq!(
        eligible,
        vec![None, None, None, None, None, Heterozygote::of(10, 12), None]
    );
    let reordered = b"chr1\t5\t.\tA\tG\t50\tPASS\t.\tGQ:GT\t40:0/1";
    assert_eq!(eligible_columns(reordered, 1), vec![Heterozygote::of(0, 1)]);
    assert_eq!(
        eligible_columns(&line(5, "PASS", &["0/1:40"]), 2),
        vec![None, None],
        "a line whose column count is not the run's"
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
        String::from_utf8(written.expect("a line")).expect("text"),
        "chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:PS\t0/1:40:7:.\t./.:.:.:."
    );
}

/// **`GT` is rewritten where `FORMAT` puts it, and an existing `PS` key is reused** (review M8):
/// with `GT` second, the phased `GT` replaces it and not the first value; with `PS` already in
/// `FORMAT`, no second `PS` is appended and the value goes in its place.
#[test]
fn gt_is_rewritten_at_its_place_and_an_existing_ps_is_reused() {
    let phase = Some((
        Heterozygote::of(0, 1).expect("a heterozygote"),
        BlockPlacement {
            block: 0,
            upper_on_copy_0: true,
            phase_set: 77,
        },
    ));
    let written = with_phase(
        b"chr1\t5\t.\tA\tG\t50\tPASS\t.\tGQ:GT:PS:DP\t40:0/1:.:7\t40:0/1:.:7",
        &[phase, None],
    );
    assert_eq!(
        String::from_utf8(written.expect("a line")).expect("text"),
        "chr1\t5\t.\tA\tG\t50\tPASS\t.\tGQ:GT:PS:DP\t40:1|0:77:7\t40:0/1:.:7"
    );
}

/// **A line the stage cannot read is written as it came and counted** (review Mi22), so the run's
/// report can warn; a line without sample columns is written as it came and not counted.
#[test]
fn a_line_the_stage_cannot_read_is_written_as_it_came_and_counted() {
    let mut stage = ReadPhaseStage::new(1);
    let mut out = Vec::new();
    stage.push(
        place(0, 5),
        b"chr1\t5\t\xff".to_vec(),
        RecordFragments::default(),
        &mut out,
    );
    stage.finish(&mut out);
    assert_eq!(out[0].1, b"chr1\t5\t\xff".to_vec());
    assert_eq!(stage.tally().lines_not_read, 1);
    let sites_only: &[u8] = b"chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT";
    assert_eq!(with_phase(sites_only, &[]).as_deref(), Some(sites_only));
}
