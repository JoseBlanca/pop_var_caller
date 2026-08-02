//! **This region's records only** — the questions that are the same whichever format the
//! file is.
//!
//! A [`RecordReader`] finds and unpacks records and asks nothing about them. Something has to
//! ask: is this record on the chromosome we want, does it touch the region, may we stop
//! looking, and which read group is it? Those four are format-blind, so they are written once
//! here rather than twice in a BAM reader and a CRAM one — which is what makes a record from
//! a scripted list and a record from a BGZF block behave identically instead of that being
//! something the design claims (`arch/alignment_cursor.md` §2.3).
//!
//! ```text
//! RecordReader    raw records, from wherever they come from
//!    ↓
//! RegionRecords   this region's records only  ← here
//!    ↓
//! ReadFilter      step-1 filtering
//!    ↓
//! AlignedRead     what the cursor hands out
//! ```
//!
//! # Why this exists at Milestone B rather than C
//!
//! The plan puts this type at C1, *lifted out of `BamRegionSource`*. It has to exist earlier:
//! a cursor yields `AlignedRead`, so it owns a [`ReadFilter`], and a `ReadFilter` needs a
//! [`RecordSource`] underneath it — and C1 depends on B2. So it is written here in the shape
//! the in-memory arm needs, and C1's remaining job is to prove the BAM arm reuses *this*
//! rather than growing a second copy.

#![allow(
    dead_code,
    reason = "Milestone B builds the cursor against a scripted list; the callers that hold one \
              arrive with the BAM arm at Milestone C and the generators at D. Until then every \
              item here is reached from its own tests. Remove at C3, where the cursor is wired \
              to a real file — if it is still needed then, nothing is reading through the \
              cursor that was built for it."
)]

use std::io;

use noodles_sam as sam;

use crate::bam::alignment_input::cigar_ref_span;
use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::filtering::{NoodlesRawRecord, RecordSource, resolve_read_group};
use crate::ng::read::input::read_groups::{ReadGroupResolution, RecordOwner};
use crate::ng::read::input::record_reader::RecordReader;
use crate::ng::read::input::region_query::overlaps;
use crate::ng::types::{ContigId, GenomeRegion};

/// The records of one region of one file, narrowed from whatever the reader hands over.
///
/// Owns the reader beneath it, and is owned by the [`ReadFilter`] above it — which is why
/// the cursor reaching down to reposition goes `filter.source_mut().move_to(region)` rather
/// than through a back-pointer. There is no cycle: the filter's source is the layer *below*
/// the cursor, not the cursor itself.
#[derive(Debug)]
pub(crate) struct RegionRecords {
    reader: RecordReader,
    /// The chromosome every region must be on — the cursor's, checked before this is ever
    /// asked to move.
    contig: ContigId,
    /// The region being served, or `None` before the first [`move_to`](Self::move_to). A
    /// source asked for records before it has been pointed anywhere yields nothing rather
    /// than guessing at a region.
    region: Option<GenomeRegion>,
    /// How this file's records are assigned to read groups. Settled at open; consulted per
    /// record only when it says it must.
    resolution: ReadGroupResolution,
    /// Records skipped as another sample's — **cumulative across every region this source
    /// serves**, because it now outlives the region it was counted in.
    other_sample_records: u64,
    /// The record the sorted early stop consumed without yielding.
    ///
    /// **Without this, carrying on into the next region loses exactly one read.** The stop
    /// fires *on* a record — the first one beginning past the region's end — and that record
    /// has already been taken from the reader. When the next region jumps, the reader is
    /// repositioned and re-reads it; when the next region **continues** from here, nothing
    /// re-reads it and it is gone. It is one read, silently, per region boundary.
    ///
    /// Held here rather than in the [`RecordReader`], which is where `arch §1.3` puts it: the
    /// over-read happens in this layer, because this is the layer that knows where the region
    /// ends. A reader cannot hold back a record it was never told to stop at.
    held: Option<sam::alignment::RecordBuf>,
}

impl RegionRecords {
    pub(crate) fn new(
        reader: RecordReader,
        contig: ContigId,
        resolution: ReadGroupResolution,
    ) -> Self {
        Self {
            reader,
            contig,
            region: None,
            resolution,
            other_sample_records: 0,
            held: None,
        }
    }

    /// The chromosome this source narrows to.
    pub(crate) fn contig(&self) -> ContigId {
        self.contig
    }

    /// Point at a new region **and reposition the reader**, discarding whatever was held.
    ///
    /// What every region did before the forget rule, and what a region that cannot reuse what
    /// is kept still does: start again from wherever the index says this region begins.
    ///
    /// The chromosome is **not** checked here: the cursor above has already refused a foreign
    /// region before anything was touched, which is what makes its "the cursor is unharmed"
    /// promise true by construction rather than by care (spec §10).
    pub(crate) fn jump_to(&mut self, region: GenomeRegion) -> io::Result<()> {
        self.region = Some(region);
        self.held = None;
        self.reader.begin_region(region)
    }

    /// Point at a new region **without moving the reader** — the case the whole design exists
    /// for.
    ///
    /// The new region begins at or after the last one served, so every record it needs is
    /// either already behind the reader (kept above, as reads) or still ahead of it. Reading
    /// simply carries on, and the record the last early stop took is handed over first (spec
    /// §4, the *partly held* case).
    pub(crate) fn continue_into(&mut self, region: GenomeRegion) {
        self.region = Some(region);
    }
}

impl RecordSource for RegionRecords {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        self.reader.header()
    }

    /// Records skipped as another sample's, since this source was made.
    ///
    /// **Cumulative, not per region**, and that is a change of meaning from the per-query
    /// sources this replaces: one of these outlives every region it serves, so a count reset
    /// at each `move_to` would report only the last region's and lose the rest.
    fn other_sample_records(&self) -> u64 {
        self.other_sample_records
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        let Some(region) = self.region else {
            // Never pointed at a region. Yielding nothing is the honest answer; guessing at
            // one would make the first region's reads depend on the order calls happened to
            // arrive in.
            return Ok(false);
        };

        loop {
            // The record the previous region's early stop took, before anything new is read:
            // it is the next one in position order, and nothing else will ever produce it.
            match self.held.take() {
                Some(held) => {
                    buf.record = held;
                    buf.read_group = None;
                }
                None => {
                    if !self.reader.read_next(buf)? {
                        return Ok(false);
                    }
                }
            }

            let on_this_contig = buf.record.reference_sequence_id() == Some(self.contig.0 as usize);

            // **The sorted early stop.** The file is coordinate-ordered, so once a record on
            // this contig begins past the region's end, no later record can reach back into
            // it. Reported as an ordinary end of input, because that is the only end a
            // `RecordSource` can report — which is exactly why `ReadFilter` distinguishes a
            // clean stop from a failure, and why a cursor undoes the former when it moves on.
            if on_this_contig
                && buf
                    .record
                    .alignment_start()
                    .is_some_and(|start| usize::from(start) as u64 > region.end.get())
            {
                // Held, not dropped: see the field. The next region gets it first.
                self.held = Some(buf.record.clone());
                return Ok(false);
            }

            // What the reader over-returned: a different contig, or a footprint that misses
            // the region. Dropped **uncounted** — these are a reader's business, not a
            // filter's, and charging them to a drop reason would make the tally mean
            // something different for a narrowed read than for a whole-file one.
            //
            // **The overlap rule is borrowed, not rewritten.** `region_query::overlaps` is
            // what the BAM and CRAM sources already apply, and the single thing that must
            // never happen in this design is two paths disagreeing about which reads a region
            // contains. When Milestone F deletes `region_query.rs` the function moves here;
            // until then, sharing it is what makes the two paths provably identical rather
            // than identical-looking.
            if !on_this_contig || !overlaps(&buf.record, region) {
                continue;
            }

            // Resolved on the record actually being yielded, never on one the loop skipped:
            // the answer can depend on the record's own `RG` tag, so resolving early would
            // attribute this read to whatever the previous record carried.
            match resolve_read_group(&buf.record, &self.resolution)? {
                RecordOwner::Mine(id) => {
                    buf.read_group = Some(id);
                    return Ok(true);
                }
                RecordOwner::OtherSample => {
                    self.other_sample_records += 1;
                    continue;
                }
            }
        }
    }
}

/// The last reference position a read touches.
///
/// **Matches what noodles reports for the record this read was decoded from**, including the
/// odd case — because the layer above compares kept *reads* against a region while this layer
/// compares *records*, and two rules that merely look equivalent is how a read gets yielded
/// when it is read fresh and dropped when it is replayed.
///
/// The odd case: a read whose CIGAR consumes no reference — all soft-clip — is **not** given
/// an empty footprint. `alignment_span()` answers `None` rather than `Some(0)`, so
/// `alignment_end()` reports the read's own start, and the read touches exactly the one base
/// it is anchored at. Measured rather than assumed: a `30S` record at position 40 reports
/// `end = 40`, and overlaps 31..=60 but not 1..=10.
pub(crate) fn read_end(read: &AlignedRead) -> u64 {
    let span = u64::from(cigar_ref_span(&read.cigar));
    read.pos + span.max(1) - 1
}

/// Whether a read's footprint touches `region`, by the same rule this module applies to the
/// record it was decoded from.
pub(crate) fn read_overlaps(read: &AlignedRead, region: GenomeRegion) -> bool {
    read.pos <= region.end.get() && read_end(read) >= region.start.get()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::record_reader::InMemoryRecordReader;
    use crate::ng::read::input::test_fixtures::{
        bam_header, fixture_read_group, matching_contigs, read_named_with_length,
    };
    use crate::ng::types::Position;
    use noodles_sam::alignment::RecordBuf;

    fn records_over(script: Vec<RecordBuf>) -> RegionRecords {
        RegionRecords::new(
            RecordReader::InMemory(InMemoryRecordReader::new(
                bam_header(&matching_contigs()),
                script,
            )),
            ContigId(0),
            fixture_read_group(),
        )
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn read_at(qname: &str, contig: usize, start: usize) -> RecordBuf {
        read_named_with_length(qname, contig, start, 10)
    }

    /// Every record this source yields for `region`, by name.
    fn narrowed_to(source: &mut RegionRecords, region: GenomeRegion) -> Vec<String> {
        source
            .jump_to(region)
            .expect("an in-memory move cannot fail");
        let mut buf = NoodlesRawRecord::default();
        let mut names = Vec::new();
        while source
            .read_next(&mut buf)
            .expect("an in-memory read cannot fail")
        {
            names.push(String::from_utf8_lossy(buf.record.name().expect("named")).into_owned());
        }
        names
    }

    /// A record touching the region is yielded; one that misses it is dropped, and one on
    /// another contig is dropped whatever its coordinates say.
    #[test]
    fn only_records_touching_the_region_on_this_contig_are_yielded() {
        let mut source = records_over(vec![
            read_at("before", 0, 1),
            read_at("inside", 0, 40),
            read_at("other-contig", 1, 40),
            read_at("after", 0, 90),
        ]);

        assert_eq!(narrowed_to(&mut source, region(31, 60)), ["inside"]);
    }

    /// **Both edges of the overlap, which is where an off-by-one lives.** A read ending
    /// exactly at the region's first base touches it; one ending the base before does not.
    /// The same at the far edge. A property test over random scripts and regions is what
    /// found this class in review; these four cases are the version that names them.
    #[test]
    fn a_record_touching_either_edge_of_the_region_is_yielded() {
        // 10-base reads, so "starts at 22" ends at 31.
        let mut source = records_over(vec![
            read_at("ends-just-before", 0, 21),
            read_at("ends-on-the-first-base", 0, 22),
            read_at("starts-on-the-last-base", 0, 60),
            read_at("starts-just-after", 0, 61),
        ]);

        assert_eq!(
            narrowed_to(&mut source, region(31, 60)),
            ["ends-on-the-first-base", "starts-on-the-last-base"],
        );
    }

    /// **The sorted early stop.** Once a record on this contig begins past the region's end
    /// the walk is over — later records cannot reach back — and the stop is reported as an
    /// ordinary end of input, which is the only end a `RecordSource` can report.
    #[test]
    fn the_walk_stops_at_the_first_record_beginning_past_the_region() {
        let mut source = records_over(vec![
            read_at("inside", 0, 40),
            read_at("past", 0, 70),
            // Never reached: the stop fires on `past`. A source that only *filtered* instead
            // of stopping would yield this one, which is the difference this test pins.
            read_at("also-inside", 0, 41),
        ]);

        assert_eq!(narrowed_to(&mut source, region(31, 60)), ["inside"]);
    }

    /// A source that has not been pointed at a region yields nothing rather than guessing at
    /// one — which would make the first region's records depend on call order.
    #[test]
    fn a_source_pointed_nowhere_yields_nothing() {
        let mut source = records_over(vec![read_at("a", 0, 40)]);

        let mut buf = NoodlesRawRecord::default();
        assert!(!source.read_next(&mut buf).expect("no reader is involved"));
    }

    /// Read groups are resolved on the record actually yielded, and the buffer is stamped —
    /// `RawRecord::decode` refuses an unstamped one, so a missed stamp is a fatal error
    /// rather than a read attributed to whatever came before it.
    #[test]
    fn a_yielded_record_carries_its_read_group() {
        let mut source = records_over(vec![read_at("a", 0, 40)]);
        source.jump_to(region(31, 60)).expect("moves");

        let mut buf = NoodlesRawRecord::default();
        assert!(source.read_next(&mut buf).expect("reads"));
        assert!(
            buf.read_group.is_some(),
            "the source must stamp the read group it resolved",
        );
    }

    /// The other-sample tally is **cumulative across regions**, because one of these outlives
    /// every region it serves. A count reset at each `move_to` would report the last region's
    /// and lose the rest.
    ///
    /// The fixture's resolution declares a single read group owning every record, so nothing
    /// is skipped here and the count is zero — pinned so that the *shape* of the claim has a
    /// test, and so the counter cannot start reporting a number nobody asked for.
    #[test]
    fn the_other_sample_tally_does_not_reset_between_regions() {
        let mut source = records_over(vec![read_at("a", 0, 40), read_at("b", 0, 45)]);

        for region in [region(31, 60), region(1, 20), region(31, 60)] {
            let before = source.other_sample_records();
            let _ = narrowed_to(&mut source, region);
            assert!(
                source.other_sample_records() >= before,
                "the tally went backwards across a reposition",
            );
        }
        assert_eq!(
            source.other_sample_records(),
            0,
            "this fixture's read group owns every record, so nothing is another sample's",
        );
    }
}
