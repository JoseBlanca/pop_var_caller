//! An aligned-reads reader with no file behind it: a scripted list, handed back in order.

use std::io;

use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;

use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::types::GenomeRegion;

/// A fixed list of records, yielded in the order it was given.
///
/// **What it yields is undecoded** — a [`NoodlesRawAlignedRead`], which is a SAM flag and a
/// mapping quality readable without unpacking anything else. The name no longer says so, so
/// this does. It matters most on *this* arm: the records are handed in already built, so it
/// would be easy to assume they arrive converted. They do not — this arm goes through the
/// same lines above it as a record freshly read from a file, which is the whole reason it can
/// serve as the oracle.
///
/// **Permanent, not a test fixture.** The forget rule — the one comparison that decides
/// which kept reads may be dropped — is the only part of the cursor design that can lose
/// reads *silently*: a rule that drops a read it should have kept produces a wrong genotype,
/// not a crash. The first attempt at it lost 3,830 of 236,081 loci while all 1,471 unit
/// tests passed. So the rule is built and driven from a scripted list, where "what should
/// this region return?" is answerable by scanning the same list by hand, before any indexed
/// file can hide a defect in it (spec §6, §11).
///
/// # What `begin_region` does here, and why it is not a search
///
/// Nothing is indexed, so there is nothing to look up: the reader rewinds to the start of
/// the list. Every region therefore sees every record, in order, and the layer above decides
/// which ones overlap — which is precisely the linear scan the BAM arm's index query has to
/// agree with. Slow by construction, and that is the point: a scan cannot skip a record the
/// index would have missed.
///
/// The records are held as written. Callers are expected to supply them in position order,
/// because that is what a coordinate-sorted file yields and what the layer above assumes;
/// this reader does not sort them, so a scripted list that is out of order is a way to drive
/// the order guard rather than a mistake this type will correct.
///
/// # A reader can be scripted to break, and that is what makes the fatal paths testable
///
/// [`with_failure_at_read`](Self::with_failure_at_read) marks a read at which this reader hands
/// back an `Err` instead of a record — a truncated block, a read cut off part-way, whatever a
/// real file does when it breaks. [`with_failing_seek`](Self::with_failing_seek) breaks the
/// *reposition* instead, which is the other way a reader can fail and a different fatal route
/// above.
///
/// **What this buys is narrower than "the chain carries a fault", because that was already
/// covered.** Two tests carry real faults up the whole chain on real inputs:
/// `open_bam.rs`'s `t10_a_truncated_file_fails_once_and_then_refuses_later_regions` truncates an
/// indexed BAM mid-walk, and `cursor.rs`'s
/// `a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short` drives a
/// reference fetch off the end of a contig. Both fail if a layer swallows the fault.
///
/// What neither of them can see is **which** [`ReadFilterError`](crate::ng::read::filtering)
/// a fault is charged to: both match `Err(_)`, so swapping `Source` and `Reference` at their
/// call sites leaves the whole suite green. A scripted fault is a fault whose kind the script
/// chose, so the test can assert the charge — measured: swapping either one fails exactly one
/// test in the tree, the scripted one.
#[derive(Debug)]
pub(crate) struct InMemoryAlignedReadsReader {
    /// The header the records' `reference_sequence_id`s are resolved against.
    header: sam::Header,
    /// The script, in the order it will be handed back.
    records: Vec<RecordBuf>,
    /// How far through `records` the reader is — the index of the next one to hand back.
    /// Reset by [`begin_region`](Self::begin_region).
    next_index: usize,
    /// The index into `records` of the read that fails instead of being handed back, or `None`
    /// for a reader that never fails. Counted in the same space as `next_index` — a read of the
    /// script, **not** a reference coordinate.
    failing_read_index: Option<usize>,
    /// Which reposition fails — counted from the first, or `None` for a reader whose seeks all
    /// succeed. See [`with_failing_seek_at`](Self::with_failing_seek_at).
    failing_seek_index: Option<usize>,
    /// How many repositions have been asked for, so `failing_seek_index` can name one.
    seeks_asked_for: usize,
}

impl InMemoryAlignedReadsReader {
    pub(crate) fn new(header: sam::Header, records: Vec<RecordBuf>) -> Self {
        Self {
            header,
            records,
            next_index: 0,
            failing_read_index: None,
            failing_seek_index: None,
            seeks_asked_for: 0,
        }
    }

    /// Fail at read `read_index` of the script, instead of handing back a record.
    ///
    /// **A read at or past the end of the script is meaningful, not a mistake**: it is the
    /// truncated file, which breaks exactly where it should have said it was finished. So the
    /// fault is decided before the script is consulted, and it fires at the first read that
    /// reaches it *or past it* — otherwise a fault scripted beyond the last record would be
    /// accepted in silence and exercise nothing, which is the one thing a fault-injection knob
    /// must not do.
    ///
    /// **The failure is not consumed**: reading again fails again, and
    /// [`begin_region`](Self::begin_region) does not clear it. A file that cannot be read stays
    /// unreadable, which is the condition the layers above are built to stop on. Nothing above
    /// ever asks twice — the filter fuses on the first one — so this is about the reader telling
    /// the truth rather than about a path anyone walks.
    pub(crate) fn with_failure_at_read(mut self, read_index: usize) -> Self {
        self.failing_read_index = Some(read_index);
        self
    }

    /// Fail the **reposition** rather than a read, at reposition `seek_index` counting from the
    /// first.
    ///
    /// The other way a reader can break, and a different fatal route: on a BAM,
    /// [`begin_region`](Self::begin_region) runs an index query, so a corrupt index fails the
    /// *move* and no read is ever attempted. A caller that swallowed that would answer the
    /// region from wherever the reader happened to be left — a plausible, silently short answer,
    /// which is the condition `CursorError::AfterFailure` exists to make loud one layer later.
    ///
    /// **Positional rather than all-or-nothing, and that is what makes it useful.** A reader
    /// whose *first* seek fails has served nothing, so a caller cannot tell "left exactly as it
    /// was" from "was never anywhere" — which is precisely the state a failed reposition must
    /// preserve. Failing a *later* seek is what reaches it.
    ///
    /// A reader whose seek fails has not moved, so `next_index` is left where it was.
    pub(crate) fn with_failing_seek_at(mut self, seek_index: usize) -> Self {
        self.failing_seek_index = Some(seek_index);
        self
    }

    pub(crate) fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Rewind to the start of the script.
    ///
    /// The region is accepted and ignored: this reader finds nothing, and the overlap test
    /// belongs to the layer above (`aligned_reads_reader/mod.rs`, the contract). Taking it anyway
    /// keeps the arm's shape identical to the ones that *will* use it, so the enum's
    /// delegation is uniform and a later arm cannot quietly need a different signature.
    pub(crate) fn begin_region(&mut self, _region: GenomeRegion) -> io::Result<()> {
        let this_seek = self.seeks_asked_for;
        self.seeks_asked_for += 1;
        if self.failing_seek_index == Some(this_seek) {
            // The reader has not moved: a failed seek leaves the position it had, which is what
            // makes swallowing the error produce a wrong answer rather than an empty one.
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("the script is set to fail reposition {this_seek}"),
            ));
        }
        self.next_index = 0;
        Ok(())
    }

    /// Hand back the next record, cloning it into the caller's reused buffer.
    ///
    /// **The clone is deliberate, and it is a real clone.** The script is replayed by every
    /// region, so handing the stored record out by move would empty the list after one pass
    /// and make the second region's answer depend on the first's — the exact class of bug
    /// this reader exists to catch elsewhere.
    ///
    /// It also allocates, every time. `RecordBuf` derives `Clone`, and a derived `Clone` gets
    /// the default `clone_from` — `*self = source.clone()` — which drops the destination's
    /// buffers rather than filling them. So this arm's cost shape is the **opposite** of a
    /// file arm's, where noodles' `read_record_buf` refills in place. That is acceptable on a
    /// reader with no file behind it and it is not something to generalise from; if the churn
    /// ever matters to a harness driving long scripts, the fix is a field-wise copy through
    /// `RecordBuf`'s `*_mut()` accessors, and it should wait for a measurement that asks.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawAlignedRead) -> io::Result<bool> {
        // Cleared before anything else, and never set: records come out raw, and a stale
        // group left in the reused buffer would attribute this record to the previous one's
        // read group without a word.
        buf.read_group = None;

        // Before the script is consulted, so a fault scripted at or past its end is the
        // truncated file rather than a clean stop — see `with_failure_at_read`.
        //
        // **Clamped to the end of the script, and neither `==` nor a bare `>=` is enough.**
        // `next_index` stops advancing once the script runs out, so it never reaches a fault
        // scripted beyond the last record: under either of those the fault would be unreachable
        // and accepted in silence, which is the one thing a fault-injection knob must not do.
        // Clamping makes every past-the-end fault fire at the end, where the file breaks.
        // `next_index` does not move: a file that cannot be read stays unreadable.
        if self
            .failing_read_index
            .is_some_and(|failing| self.next_index >= failing.min(self.records.len()))
        {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                // Nothing above this reads the `ErrorKind` — `io::Error::other` re-kinds the
                // outer error to `Other` before any caller sees it — so the message is what an
                // operator gets, and it names the read.
                format!("the script is set to fail at read {}", self.next_index),
            ));
        }

        let Some(record) = self.records.get(self.next_index) else {
            return Ok(false);
        };
        self.next_index += 1;
        buf.record.clone_from(record);
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::{bam_header, matching_contigs, read_named};
    use crate::ng::types::{ContigId, Position};

    fn reader(records: Vec<RecordBuf>) -> InMemoryAlignedReadsReader {
        InMemoryAlignedReadsReader::new(bam_header(&matching_contigs()), records)
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// Drain the reader into the names it yielded — how most tests below state their
    /// expectation, because a name identifies a record across two independently-produced
    /// streams, which is what the region-query oracle needs.
    ///
    /// **Names alone are not enough**, and [`the_whole_record_survives_the_clone`] is why: a
    /// clone carrying the name and garbage everywhere else would satisfy every assertion made
    /// through this helper. Use [`drain_records`] where the record's *content* is the point.
    fn drain(reader: &mut InMemoryAlignedReadsReader) -> Vec<String> {
        drain_records(reader)
            .iter()
            .map(|record| String::from_utf8_lossy(record.name().expect("named")).into_owned())
            .collect()
    }

    /// Drain the reader into the whole records it yielded.
    fn drain_records(reader: &mut InMemoryAlignedReadsReader) -> Vec<RecordBuf> {
        let mut buf = NoodlesRawAlignedRead::default();
        let mut records = Vec::new();
        while reader
            .read_next(&mut buf)
            .expect("this script has no scripted fault")
        {
            records.push(buf.record.clone());
        }
        records
    }

    /// **A reader is unarmed unless a fault is scripted onto it**, which nothing states
    /// separately because it does not need to: [`drain`] `expect`s every read, so a `new` that
    /// armed its reader by accident fails this test and eight others in this module before it
    /// reaches anything that scripted a fault deliberately. Verified — setting
    /// `failing_read_index: Some(0)` in `new` fails 60 tests.
    #[test]
    fn the_script_comes_back_in_the_order_it_was_given() {
        let mut reader = reader(vec![
            read_named("a", 0, 10),
            read_named("b", 0, 20),
            read_named("c", 0, 30),
        ]);
        reader.begin_region(region(1, 100)).expect("positions");

        assert_eq!(drain(&mut reader), ["a", "b", "c"]);
    }

    /// **The property that makes this a usable oracle**: every region sees the whole script,
    /// so what a region *should* return is answerable by scanning the same list by hand. A
    /// reader that consumed its records would make the second region's answer depend on the
    /// first's.
    #[test]
    fn every_region_replays_the_whole_script_from_the_start() {
        let mut reader = reader(vec![read_named("a", 0, 10), read_named("b", 0, 20)]);

        for region in [region(1, 100), region(50, 60), region(1, 5)] {
            reader.begin_region(region).expect("positions");
            assert_eq!(
                drain(&mut reader),
                ["a", "b"],
                "after begin_region({region:?})"
            );
        }
    }

    /// A reader drained without being repositioned stays drained, rather than silently
    /// starting over — which would make a caller that read past the end see the script twice.
    #[test]
    fn reading_past_the_end_keeps_returning_nothing() {
        let mut reader = reader(vec![read_named("a", 0, 10)]);
        reader.begin_region(region(1, 100)).expect("positions");

        assert_eq!(drain(&mut reader), ["a"]);
        let mut buf = NoodlesRawAlignedRead::default();
        assert!(!reader.read_next(&mut buf).expect("ends"));
        assert!(!reader.read_next(&mut buf).expect("still ends"));
    }

    /// An empty script is a legal one, and answers immediately rather than by index panic.
    #[test]
    fn an_empty_script_yields_nothing() {
        let mut reader = reader(Vec::new());
        reader.begin_region(region(1, 100)).expect("positions");

        assert!(drain(&mut reader).is_empty());
    }

    /// **The stale-read-group guard, from the one direction that can fail.** The buffer is
    /// reused across reads and across readers, so a reader that did not clear the field
    /// would hand this record out wearing the previous one's read group — and
    /// `RawAlignedRead::decode` would accept it, because the field is populated.
    #[test]
    fn a_record_comes_out_raw_with_no_read_group_stamped() {
        let mut reader = reader(vec![read_named("a", 0, 10)]);
        reader.begin_region(region(1, 100)).expect("positions");

        // A group left over from an earlier pass through some other source.
        // Both fields spelled out: the spread would fill exactly one and would
        // silently absorb a third if the buffer ever grows one — on the very
        // literal whose subject is the field being spread past.
        let mut buf = NoodlesRawAlignedRead {
            record: RecordBuf::default(),
            read_group: Some(crate::ng::types::ReadGroupId(7)),
        };
        assert!(reader.read_next(&mut buf).expect("reads"));

        assert!(
            buf.read_group.is_none(),
            "the reader must clear the reused buffer's read group; resolving it belongs to \
             the layer above",
        );
    }

    /// The script is handed back as given, **not sorted**. A caller supplying records out of
    /// position order is driving the order guard above this layer, and this reader silently
    /// reordering them would disarm exactly that test.
    #[test]
    fn an_out_of_order_script_is_not_quietly_sorted() {
        let mut reader = reader(vec![read_named("late", 0, 90), read_named("early", 0, 10)]);
        reader.begin_region(region(1, 100)).expect("positions");

        assert_eq!(drain(&mut reader), ["late", "early"]);
    }

    /// **Spec §7's first named case: a region abandoned half-way.** A caller that stops
    /// pulling and moves elsewhere leaves this reader mid-script, and the next region must
    /// still see all of it. A reader that only rewound when it had been drained would hand
    /// the next region a truncated script — fewer records, no error — and this arm is the
    /// oracle the forget rule will be judged against, so the result is not a red test but a
    /// *wrong oracle*.
    ///
    /// Found by mutation: `if self.next_index >= self.records.len() { self.next_index = 0 }`
    /// passed the whole suite.
    #[test]
    fn a_reader_abandoned_part_way_through_still_rewinds() {
        let mut reader = reader(vec![
            read_named("a", 0, 10),
            read_named("b", 0, 20),
            read_named("c", 0, 30),
        ]);
        reader.begin_region(region(1, 100)).expect("positions");

        // One record, then the caller moves on.
        let mut buf = NoodlesRawAlignedRead::default();
        assert!(reader.read_next(&mut buf).expect("reads"));

        reader.begin_region(region(50, 60)).expect("repositions");
        assert_eq!(
            drain(&mut reader),
            ["a", "b", "c"],
            "a reader abandoned after one record did not rewind",
        );
    }

    /// **The clone has to carry the whole record, not just enough to be recognised.** Every
    /// other test here identifies records by name, so a clone that dropped `alignment_start`
    /// — the one field the forget rule compares — would satisfy all of them. Both mutations
    /// were run and survived before this test existed.
    #[test]
    fn the_whole_record_survives_the_clone() {
        let script = vec![
            read_named("a", 0, 10),
            read_named("b", 1, 20),
            read_named("c", 0, 30),
        ];
        let mut reader = reader(script.clone());
        reader.begin_region(region(1, 100)).expect("positions");

        assert_eq!(drain_records(&mut reader), script);
    }

    /// The buffer is reused across every read of a pass, so the guard has to hold on every
    /// read — not only the first. A reader that cleared once per pass would leave the
    /// *second* record wearing whatever the caller last had, which is the shape of the bug
    /// `NoodlesRawAlignedRead`'s own doc records having been bitten by.
    #[test]
    fn every_record_of_a_pass_comes_out_with_no_read_group() {
        let mut reader = reader(vec![
            read_named("a", 0, 10),
            read_named("b", 0, 20),
            read_named("c", 0, 30),
        ]);
        reader.begin_region(region(1, 100)).expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        let mut read = 0;
        loop {
            // Re-stamped between every read, as a caller resolving read groups would.
            buf.read_group = Some(crate::ng::types::ReadGroupId(7));
            if !reader.read_next(&mut buf).expect("reads") {
                break;
            }
            read += 1;
            assert!(
                buf.read_group.is_none(),
                "record {read} of the pass kept a read group the reader must have cleared",
            );
        }
        assert_eq!(read, 3, "the pass did not reach every record");
    }

    /// **The scripted fault fires at its own read, and not before.**
    ///
    /// The "and not before" half is what makes this able to fail: a reader that failed on
    /// every read — the shape the deleted `ErroringSource` double had — would satisfy an
    /// assertion that only looked for an error, and would then be useless for driving a fault
    /// that arrives *mid-walk*, which is the interesting one.
    #[test]
    fn a_scripted_fault_fires_at_its_own_read_and_not_before() {
        let mut reader = reader(vec![
            read_named("a", 0, 10),
            read_named("b", 0, 20),
            read_named("c", 0, 30),
        ])
        .with_failure_at_read(1);
        reader.begin_region(region(1, 100)).expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        assert!(
            reader.read_next(&mut buf).expect("the first read is clean"),
            "the reader must hand back the reads before the scripted fault",
        );
        assert_eq!(
            String::from_utf8_lossy(buf.record.name().expect("named")),
            "a",
        );

        let error = reader
            .read_next(&mut buf)
            .expect_err("the second read is the scripted fault");
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
        assert!(
            // "at read 1" rather than "read 1", which also matches "read 10".
            error.to_string().contains("at read 1"),
            "the message must say which read failed, got {error}",
        );
    }

    /// **A fault does not heal.** Reading again fails again, and repositioning does not clear
    /// it — a file that cannot be read stays unreadable, which is the condition every layer
    /// above is built to stop on. A reader that consumed its scripted fault would let a caller
    /// walk on past a broken file and answer the next region short.
    #[test]
    fn a_scripted_fault_is_not_consumed_by_reading_or_by_repositioning() {
        let mut reader =
            reader(vec![read_named("a", 0, 10), read_named("b", 0, 20)]).with_failure_at_read(0);
        reader.begin_region(region(1, 100)).expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        assert!(reader.read_next(&mut buf).is_err(), "it fails at once");
        assert!(reader.read_next(&mut buf).is_err(), "and again");

        reader.begin_region(region(50, 60)).expect("repositions");
        assert!(
            reader.read_next(&mut buf).is_err(),
            "repositioning cleared the scripted fault",
        );
    }

    /// **A fault at a non-zero read survives a rewind in the right place**, which the test
    /// above cannot tell: at read 0 "rewound to the start" and "rewound into the fault" are the
    /// same observation. Here they are not — the first read must come back clean a second time,
    /// and the fault must fire on the second read again.
    ///
    /// Kills a rewind that shifts the fault relative to the script rather than leaving it
    /// where it was.
    #[test]
    fn a_scripted_fault_survives_a_rewind_at_the_same_read() {
        let mut reader =
            reader(vec![read_named("a", 0, 10), read_named("b", 0, 20)]).with_failure_at_read(1);
        reader.begin_region(region(1, 100)).expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        assert!(reader.read_next(&mut buf).expect("the first read is clean"));
        assert!(reader.read_next(&mut buf).is_err(), "the scripted fault");

        reader.begin_region(region(50, 60)).expect("repositions");
        assert!(
            reader
                .read_next(&mut buf)
                .expect("the first read is clean again"),
            "the rewind did not return to the start of the script",
        );
        assert_eq!(
            String::from_utf8_lossy(buf.record.name().expect("named")),
            "a",
        );
        assert!(
            reader.read_next(&mut buf).is_err(),
            "the fault moved relative to the rewound script",
        );
    }

    /// **A fault scripted at or past the end of the script is the truncated file**, and it must
    /// fail rather than report the clean end of input it looks like from inside.
    ///
    /// Two shapes, and both have caught an implementation. Checking the script *before* the
    /// fault turns the at-the-end case into `Ok(false)`; comparing the fault for **equality**
    /// with the read index makes anything further past the end unreachable, because the index
    /// stops advancing at the end of the script — so the fault would be accepted in silence and
    /// exercise nothing. A truncated file reported as a finished one is a silently short answer.
    #[test]
    fn a_fault_scripted_at_or_past_the_end_of_the_script_fails_rather_than_ending_cleanly() {
        let mut buf = NoodlesRawAlignedRead::default();

        // Exactly at the end: the file breaks where it should have said it was finished.
        let mut at_the_end = reader(vec![read_named("a", 0, 10)]).with_failure_at_read(1);
        at_the_end.begin_region(region(1, 100)).expect("positions");
        assert!(
            at_the_end.read_next(&mut buf).expect("the one record"),
            "reads"
        );
        assert!(
            at_the_end.read_next(&mut buf).is_err(),
            "the reader reported a clean end of input where the script says it breaks",
        );

        // Well past the end — the case an equality comparison makes unreachable.
        let mut past_the_end = reader(vec![read_named("a", 0, 10)]).with_failure_at_read(5);
        past_the_end
            .begin_region(region(1, 100))
            .expect("positions");
        assert!(
            past_the_end.read_next(&mut buf).expect("the one record"),
            "reads"
        );
        assert!(
            past_the_end.read_next(&mut buf).is_err(),
            "a fault scripted well past the end of the script never fired",
        );

        // The degenerate end: a file that broke before its first record, not one that had none.
        let mut empty = reader(Vec::new()).with_failure_at_read(0);
        empty.begin_region(region(1, 100)).expect("positions");
        assert!(
            empty.read_next(&mut buf).is_err(),
            "an empty script with a fault at read 0 reported a clean end of input",
        );
    }

    /// **A scripted seek failure breaks the reposition, and no read is attempted.** The other
    /// way a reader can break: on a BAM `begin_region` runs an index query, so a corrupt index
    /// fails the *move* rather than a read.
    #[test]
    fn a_scripted_seek_failure_breaks_the_reposition() {
        let mut reader = reader(vec![read_named("a", 0, 10)]).with_failing_seek_at(0);

        let error = reader
            .begin_region(region(1, 100))
            .expect_err("the scripted seek failure");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(
            reader
                .read_next(&mut NoodlesRawAlignedRead::default())
                .is_ok(),
            "a failed seek breaks the move, not the reads — the reader simply has not moved",
        );
    }

    /// `header()` hands back **this reader's** header, not merely one of the right shape: the
    /// `@SQ` list is what a decode resolves a record's `reference_sequence_id` against, so a
    /// header with the right number of contigs under different names would put every read on
    /// the wrong chromosome. Asserting the count alone let that mutation through.
    #[test]
    fn the_header_is_the_one_the_reader_was_built_with() {
        let reader = reader(vec![read_named("a", 0, 10)]);

        let names: Vec<String> = reader
            .header()
            .reference_sequences()
            .keys()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect();
        let expected: Vec<String> = matching_contigs()
            .iter()
            .map(|(name, _, _)| (*name).to_owned())
            .collect();
        assert_eq!(names, expected);
    }
}
