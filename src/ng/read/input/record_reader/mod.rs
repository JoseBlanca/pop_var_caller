//! **Where a cursor's records come from** — one arm per format, and one contract they all
//! honour.
//!
//! A [`RecordReader`] finds records and unpacks them. That is *all* it does: it runs no
//! contig test, no overlap test, no read-group resolution and no early stop. Those are the
//! same questions whatever the file format is, so they are asked once, above this, by the
//! layer that turns records into reads (`arch/alignment_cursor.md` §2.3). A record from a
//! scripted list and a record freshly decoded from a file therefore go through *the same
//! lines* above this point, which is what makes them behave identically — rather than that
//! being something this module claims.
//!
//! (A read *replayed* from what the cursor kept is a different matter, and spec §5's sentence
//! is about **reads**, not records: a kept read sits **above** the filter, so replaying it
//! skips decode and filtering entirely — arch §2.3. It goes through *fewer* lines, not the
//! same ones. Nothing at this layer keeps anything.)
//!
//! # The contract, stated once — which is a place to look, not a mechanism
//!
//! Nothing below enforces these on a new arm: they are prose, and prose does not fail a
//! build. What will hold two arms together is the oracle the plan builds next — a run of
//! regions through one cursor, compared against a plain scan of the same file (spec §11.3) —
//! and, when the second arm lands, a shared test harness driving both. Until then this list
//! is where to check an arm against, and its being only a list is worth knowing.
//!
//! - **[`begin_region`](RecordReader::begin_region)** positions the reader for a region on
//!   the chromosome it is reading. Cheap when the reader is already there, a seek when it is
//!   not.
//! - **[`read_next`](RecordReader::read_next)** fills the caller's buffer with the next
//!   record **in position order** and answers whether it filled one. It is the shape
//!   [`RecordSource::read_next`](crate::ng::read::filtering::RecordSource::read_next)
//!   already has, because the layer above this is a `RecordSource`.
//!
//!   **One buffer per pass, but not one allocation per pass, and the difference is per
//!   arm.** The caller keeps a single buffer and hands the same `&mut` to every call —
//!   that part is the contract. Whether filling it *reallocates* is the arm's business: a
//!   file arm refills in place through noodles' `read_record_buf`, which is true reuse,
//!   while the in-memory arm clones, because `RecordBuf` derives `Clone` and a derived
//!   `Clone` gets the default `clone_from` (`*self = source.clone()`) rather than one that
//!   keeps the destination's buffers. Fine on an arm with no file behind it; **not a
//!   property to generalise from** — it is the opposite of what the file arms do.
//! - **Nothing is unpacked ahead of demand.** A block is decompressed only when a
//!   `read_next` needs a record inside it, so a caller that pulls one record and moves on
//!   has unpacked at most one block.
//! - **Nothing is kept across regions.** What a cursor keeps is *its* reads — decoded and
//!   filtered, held one level up — so a read is transformed once rather than once per region
//!   that returns it (spec §5). A reader holds only its position, and (from Milestone C) the
//!   single record the sorted early stop consumes without yielding.
//! - **Records come out raw.** In particular `read_group` is cleared, never stamped: the
//!   buffer is reused, and a reader that left the previous record's group in place would
//!   attribute this record to it silently. Clearing turns that into a refusal at
//!   [`RawRecord::decode`](crate::ng::read::filtering::RawRecord::decode).
//!
//! # What is here so far
//!
//! Only [`InMemoryRecordReader`]. The BAM arm lands in Milestone C and the CRAM arm in
//! Milestone E — and the order is deliberate: the forget rule that decides which kept reads
//! may be dropped is the one part of this design that can lose reads *silently*, so it is
//! built and tested against a reader with no file behind it before any real input can hide a
//! defect in it (spec §6, and the plan's principles).

// Covers this module and `in_memory` below it, so the reason is stated once.
#![allow(
    dead_code,
    reason = "Milestone A is types only: nothing constructs a RecordReader until the cursor \
              does at B1, so every item here is reachable from its own tests and from nowhere \
              else. Remove this attribute at B1 — if it is still needed then, the cursor is \
              not using the reader it was built for. `expect` would be the self-removing \
              choice and does not work here: the tests *do* use these items, so the lint \
              fires for the lib target and not the test target, and `--all-targets` then \
              reports the expectation unfulfilled — an error under this repo's -D warnings."
)]

pub(crate) mod in_memory;

use std::io;

use noodles_sam as sam;

use crate::ng::read::filtering::NoodlesRawRecord;
use crate::ng::types::GenomeRegion;

pub(crate) use in_memory::InMemoryRecordReader;

/// Finds records and unpacks them, one variant per place they can come from.
///
/// **An enum rather than a trait**, and the reason is not taste: the set of formats is
/// closed, and a trait would add a type parameter through four layers — it would reach
/// `PileupGenerator`, which already carries three, and then the walker, the generator and
/// the dispatcher. The dynamic-dispatch alternative pays a virtual call per record, about a
/// million per run (spec §5).
#[derive(Debug)]
pub(crate) enum RecordReader {
    /// A fixed list of records with no file behind it. **Permanent, not test-only**: it is
    /// what lets the forget rule be driven from a scripted list and compared against a plain
    /// scan of the same list, which is the oracle the first attempt at this feature did not
    /// have.
    InMemory(InMemoryRecordReader),
}

impl RecordReader {
    /// The header of the file these records came from — the `@SQ` list a decode resolves
    /// contig references against.
    pub(crate) fn header(&self) -> &sam::Header {
        match self {
            Self::InMemory(reader) => reader.header(),
        }
    }

    /// Position for `region`, which must be on the chromosome this reader is reading.
    ///
    /// The reader does **not** remember the region: it is not what decides which records
    /// overlap it, and a reader that filtered on its own would be a second copy of the test
    /// the layer above already makes.
    pub(crate) fn begin_region(&mut self, region: GenomeRegion) -> io::Result<()> {
        match self {
            Self::InMemory(reader) => reader.begin_region(region),
        }
    }

    /// Fill `buf` with the next record in position order. `Ok(false)` at the end, after
    /// which `buf` holds a stale record the caller must not read.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        match self {
            Self::InMemory(reader) => reader.read_next(buf),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::{bam_header, matching_contigs, read_named};
    use crate::ng::types::{ContigId, Position};

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// Every method the enum offers reaches its arm.
    ///
    /// **`begin_region` is repositioned after a drain, not before one**, which is the whole
    /// difference between this test and a version that cannot fail: a freshly built reader
    /// is already at the start, so an enum whose `begin_region` never reached the arm would
    /// look identical. Verified by mutation — replacing the enum's body with `Ok(())` leaves
    /// the earlier version green and fails this one.
    #[test]
    fn the_enum_forwards_every_contract_method_to_its_arm() {
        let header = bam_header(&matching_contigs());
        let mut reader = RecordReader::InMemory(InMemoryRecordReader::new(
            header.clone(),
            vec![read_named("a", 0, 10), read_named("b", 0, 20)],
        ));

        assert_eq!(
            reader.header().reference_sequences().len(),
            header.reference_sequences().len(),
        );

        let mut buf = NoodlesRawRecord::default();
        reader.begin_region(region(0, 1, 100)).expect("positions");
        let mut read_to_the_end = |reader: &mut RecordReader| {
            let mut records = 0;
            while reader
                .read_next(&mut buf)
                .expect("an in-memory read cannot fail")
            {
                records += 1;
            }
            records
        };
        assert_eq!(read_to_the_end(&mut reader), 2);

        // Drained. Only a `begin_region` that reaches the arm can produce records again.
        reader.begin_region(region(0, 1, 100)).expect("repositions");
        assert_eq!(
            read_to_the_end(&mut reader),
            2,
            "the enum's begin_region did not reach its arm: the reader stayed drained",
        );
    }
}
