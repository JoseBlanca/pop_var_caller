//! A record reader with no file behind it: a scripted list, handed back in order.

use std::io;

use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;

use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::types::GenomeRegion;

/// A fixed list of records, yielded in the order it was given.
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
#[derive(Debug)]
pub(crate) struct InMemoryRecordReader {
    /// The header the records' `reference_sequence_id`s are resolved against.
    header: sam::Header,
    /// The script, in the order it will be handed back.
    records: Vec<RecordBuf>,
    /// How far through `records` the reader is — the index of the next one to hand back.
    /// Reset by [`begin_region`](Self::begin_region).
    next_index: usize,
}

impl InMemoryRecordReader {
    pub(crate) fn new(header: sam::Header, records: Vec<RecordBuf>) -> Self {
        Self {
            header,
            records,
            next_index: 0,
        }
    }

    pub(crate) fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Rewind to the start of the script.
    ///
    /// The region is accepted and ignored: this reader finds nothing, and the overlap test
    /// belongs to the layer above (`record_reader/mod.rs`, the contract). Taking it anyway
    /// keeps the arm's shape identical to the ones that *will* use it, so the enum's
    /// delegation is uniform and a later arm cannot quietly need a different signature.
    pub(crate) fn begin_region(&mut self, _region: GenomeRegion) -> io::Result<()> {
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

    fn reader(records: Vec<RecordBuf>) -> InMemoryRecordReader {
        InMemoryRecordReader::new(bam_header(&matching_contigs()), records)
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
    fn drain(reader: &mut InMemoryRecordReader) -> Vec<String> {
        drain_records(reader)
            .iter()
            .map(|record| String::from_utf8_lossy(record.name().expect("named")).into_owned())
            .collect()
    }

    /// Drain the reader into the whole records it yielded.
    fn drain_records(reader: &mut InMemoryRecordReader) -> Vec<RecordBuf> {
        let mut buf = NoodlesRawAlignedRead::default();
        let mut records = Vec::new();
        while reader
            .read_next(&mut buf)
            .expect("an in-memory read cannot fail")
        {
            records.push(buf.record.clone());
        }
        records
    }

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
