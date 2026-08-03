//! **Where a cursor's records come from** — one arm per format, and one contract they all
//! honour.
//!
//! An [`AlignedReadsReader`] finds reads and unpacks them. That is *all* it does: it runs no
//! contig test, no overlap test, no read-group resolution and no early stop. Those are the
//! same questions whatever the file format is, so they are asked once, above this, by the
//! layer that turns records into reads (`arch/alignment_cursor.md` §2.3). A record from a
//! scripted list and a record freshly decoded from a file therefore go through *the same
//! lines* above this point, which is what makes them behave identically — rather than that
//! being something this module claims.
//!
//! **What every arm yields is undecoded**, and the type name no longer says so, so the module
//! does: a [`NoodlesRawAlignedRead`] is a SAM flag and a mapping quality readable without
//! unpacking anything else. Nothing here builds an [`AlignedRead`](crate::ng::read::AlignedRead)
//! — the conversion happens above, and only for the reads that clear the flag/MAPQ filter,
//! which is the whole point of the split (`spec/read_filtering_stages.md` §2). Each arm's own
//! doc repeats it, because that is where a caller meets one.
//!
//! **"Reads" in the type name, "records" in the contract below, and the difference is on
//! purpose.** A caller of this module wants *reads*, so that is what the type promises. A
//! record is the *encoding* of one — the line in the file, the entry in a CRAM container — and
//! the contract below is about finding and unpacking those, which is the only level at which
//! the arms differ from each other. Where this module says "record" it means the thing on
//! disk; where it says "read" it means what a caller gets back.
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
//! - **[`begin_region`](AlignedReadsReader::begin_region)** positions the reader for a region on
//!   the chromosome it is reading. Cheap when the reader is already there, a seek when it is
//!   not.
//!
//!   **It *positions*; it must never *bound*.** After it, reading on must yield every record
//!   from that point to the end of the chromosome — not just the ones inside the region it
//!   was handed. The region says *where to start*, and nothing else. This is not a style
//!   preference: the cursor above serves a forward region by carrying on reading without
//!   repositioning at all (`RegionRawAlignedReads::continue_into`), so a reader that had quietly
//!   stopped at the previous region's end would lose every record past it, silently, for
//!   every region after the first.
//!
//!   **The trap this names is specific and easy to fall into at Milestone C.** noodles'
//!   `query()` returns an iterator bounded by the interval it was given, which satisfies
//!   every other word of this contract and breaks this one. A BAM arm must either query a
//!   range that runs to the end of the chromosome, or re-query when it runs out — and either
//!   way it owes this contract a test.
//! - **[`read_next`](AlignedReadsReader::read_next)** fills the caller's buffer with the next
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
//! - **Records come out raw, and `read_group` is cleared rather than left stale.** The buffer
//!   is reused, so a reader that left the previous record's group in place would attribute
//!   this record to it silently; clearing turns that into a refusal at
//!   [`RawAlignedRead::decode`](crate::ng::read::aligned_read::RawAlignedRead::decode).
//!
//!   **The CRAM arm stamps it instead of clearing it, and that is the one exception in this
//!   contract** — a CRAM stores the read group as a container-level number, so it is decided
//!   at decode and travels with the record ([`CramAlignedReadsReader`]'s own doc says why).
//!   Stated here because this list is where a new arm's author checks their arm, and a
//!   universal with a live counterexample is how an arm ends up clearing a field it should
//!   have carried.
//!
//! # What is here
//!
//! All three arms: [`InMemoryAlignedReadsReader`], [`BamAlignedReadsReader`] and
//! [`CramAlignedReadsReader`], plus [`container`], which holds the CRAM decode's unit of work
//! and is not an arm. The order they were built in was deliberate — the forget rule that
//! decides which kept reads may be dropped is the one part of this design that can lose reads
//! *silently*, so it was built and tested against a reader with no file behind it before any
//! real input could hide a defect in it (spec §6, and the plan's principles).

// Covers this module and `in_memory` below it, so the reason is stated once.
// Narrowed after Milestone C: the BAM arm is live, and what is left unused is the in-memory
// arm's constructor, which the cursor's own tests reach and nothing else does. The blanket
// version of this said "remove at B1" and had outlived that by two milestones — a trigger
// nobody checks is worse than no trigger.
#![allow(
    dead_code,
    reason = "the in-memory arm is reached from tests and from the differential harness; it \
              is permanent, not scaffolding (see the enum's own doc)"
)]

mod bam;
mod container;
mod cram;
mod in_memory;

use std::io;

use noodles_sam as sam;

use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::types::GenomeRegion;

pub(crate) use bam::BamAlignedReadsReader;
pub(crate) use cram::CramAlignedReadsReader;
pub(crate) use in_memory::InMemoryAlignedReadsReader;

/// Finds reads and unpacks them, one variant per place they can come from.
///
/// **What it yields is undecoded** — a [`NoodlesRawAlignedRead`], not an
/// [`AlignedRead`](crate::ng::read::AlignedRead). The name says *reads* because that is what a
/// caller wants; it does not say *raw*, so this does (`spec/read_filtering_stages.md` §6).
///
/// **An enum rather than a trait**, and the reason is not taste: the set of formats is
/// closed, and a trait would add a type parameter through four layers — it would reach
/// `PileupGenerator`, which already carries three, and then the walker, the generator and
/// the dispatcher. The dynamic-dispatch alternative pays a virtual call per record, about a
/// million per run (spec §5).
#[derive(Debug)]
pub(crate) enum AlignedReadsReader {
    /// A fixed list of records with no file behind it. **Permanent, not test-only**: it is
    /// what lets the forget rule be driven from a scripted list and compared against a plain
    /// scan of the same list, which is the oracle the first attempt at this feature did not
    /// have.
    InMemory(InMemoryAlignedReadsReader),
    /// A BAM, read through its index.
    Bam(BamAlignedReadsReader),
    /// A CRAM, read through its `.crai` — container by container.
    Cram(CramAlignedReadsReader),
}

impl AlignedReadsReader {
    /// The header of the file these records came from — the `@SQ` list a decode resolves
    /// contig references against.
    pub(crate) fn header(&self) -> &sam::Header {
        match self {
            Self::InMemory(reader) => reader.header(),
            Self::Bam(reader) => reader.header(),
            Self::Cram(reader) => reader.header(),
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
            Self::Bam(reader) => reader.begin_region(region),
            Self::Cram(reader) => reader.begin_region(region),
        }
    }

    /// Fill `buf` with the next record in position order. `Ok(false)` at the end, after
    /// which `buf` holds a stale record the caller must not read.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawAlignedRead) -> io::Result<bool> {
        match self {
            Self::InMemory(reader) => reader.read_next(buf),
            Self::Bam(reader) => reader.read_next(buf),
            Self::Cram(reader) => reader.read_next(buf),
        }
    }

    /// Records this reader dropped as another sample's, before they were built.
    ///
    /// **Only the CRAM arm has any**, and that is the one asymmetry in this contract. A CRAM
    /// decides who owns a record while decoding a container, because the read group is a number
    /// the container carries rather than a tag on the record — so a foreign record is dropped
    /// there, and the layer above never sees it to step over. Every other arm hands its records
    /// up untouched and the layer above does the skipping and the counting.
    ///
    /// The consequence for the number is real and is stated where it is counted: this one is
    /// container-granular and can run ahead of where a walk has reached, while the BAM path's is
    /// exact.
    pub(crate) fn other_sample_records(&self) -> u64 {
        match self {
            Self::InMemory(_) | Self::Bam(_) => 0,
            Self::Cram(reader) => reader.other_sample_records(),
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
        let mut reader = AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
            header.clone(),
            vec![read_named("a", 0, 10), read_named("b", 0, 20)],
        ));

        assert_eq!(
            reader.header().reference_sequences().len(),
            header.reference_sequences().len(),
        );

        let mut buf = NoodlesRawAlignedRead::default();
        reader.begin_region(region(0, 1, 100)).expect("positions");
        let mut read_to_the_end = |reader: &mut AlignedReadsReader| {
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
