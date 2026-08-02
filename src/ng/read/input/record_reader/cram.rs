//! Records from a CRAM: the `.crai` says which container to start in, and reading carries on
//! from there.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::ng::read::filtering::NoodlesRawRecord;
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::region_query::{DecodedContainer, decode_container_at};
use crate::ng::types::GenomeRegion;

/// A CRAM reader that stays where it is between regions.
///
/// # What is different about CRAM, in one paragraph
///
/// A BAM is read record by record; a CRAM is read **container by container**. A container
/// holds on the order of ten thousand records and must be decompressed and decoded whole
/// before any one of them can be looked at. So this reader's unit of work is a container: it
/// decodes one, serves its records in order, and moves to the next.
///
/// # The one thing this must get right
///
/// [`begin_region`](Self::begin_region) **positions; it never bounds.** After it, reading on
/// yields every record from that point to the end of the chromosome — not just the ones inside
/// the region it was handed. That is the contract in `record_reader/mod.rs`, and it is
/// load-bearing rather than tidy: the cursor above serves a forward region by *not
/// repositioning at all*, so a reader that had quietly stopped at the previous region's end
/// would lose every record past it, silently, for every region after the first.
///
/// **This is where the per-region CRAM source and this one differ most**, and it is not a
/// re-layering: `CramRegionSource` stops its `.crai` walk at the first container beginning past
/// the region's end, and filters each container's records against the region. Both belong to a
/// reader that answers one region and is then thrown away. Neither can be here.
///
/// # One record, one read group — the documented exception
///
/// Every other arm hands up a record with no read group attached, and the layer above resolves
/// it from the record's `RG` tag. This arm attaches it, because a CRAM does not have an `RG`
/// tag to resolve: it stores the read group as a **number**, an index into the header's `@RG`
/// list. Deciding it at decode is what lets every auxiliary tag be dropped — 6.2 MiB per open
/// file, measured, the largest single item in what an open file cost — and re-inflating a
/// number into a string so the layer above can parse it back would be perverse. So the answer
/// travels with the record, and `RegionRecords` uses it rather than asking again.
pub(crate) struct CramRecordReader {
    reader: cram::io::Reader<File>,
    /// Parsed once at open and shared, never re-read per region.
    header: Arc<sam::Header>,
    /// The reference bases decoding consults, for this cursor's chromosome. Cheap to clone —
    /// it is internally shared — and cloned per decode, as noodles requires.
    repository: fasta::Repository,
    /// **This contig's** `.crai` entries, in file order, grouped once at open.
    entries: Arc<[cram::crai::Record]>,
    /// A copy of the file's, settled at open. Owned, so this reader carries no lifetime.
    resolution: ReadGroupResolution,
    path: Arc<Path>,
    /// The next `.crai` entry to decode. Set by [`begin_region`](Self::begin_region) and
    /// advanced by reading; never reset by reading, so a walk carries on to the contig's end.
    next_entry: usize,
    /// The container being served, and how far into it we have got.
    ///
    /// **One container, and it is not a cache.** It is what is currently being read: a
    /// container is decoded, drained, and dropped when the next one is decoded. What survives
    /// between regions is the cursor's *reads*, one layer up (spec §5).
    container: Option<DecodedContainer>,
    served: usize,
    /// The offset last decoded, so a container that appears under several `.crai` entries — one
    /// per slice — is decoded once rather than served twice.
    last_decoded_offset: Option<u64>,
    /// Records dropped at decode as another sample's, summed over the containers this reader
    /// has decoded. Container-granular; see [`DecodedContainer::other_sample_records`].
    other_sample_records: u64,
    /// Latched when the entries run out or a read fails. Cleared by the next
    /// [`begin_region`](Self::begin_region) — and **only** by it, so a reader that has run off
    /// the end of the contig stays quiet until it is repositioned rather than restarting.
    done: bool,
}

/// Hand-written because noodles' reader is not `Debug`, and so the output says what identifies
/// this reader rather than dumping a parsed index.
impl std::fmt::Debug for CramRecordReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramRecordReader")
            .field("path", &self.path)
            .field("entries", &self.entries.len())
            .field("next_entry", &self.next_entry)
            .field("done", &self.done)
            .finish_non_exhaustive()
    }
}

impl CramRecordReader {
    pub(crate) fn new(
        reader: cram::io::Reader<File>,
        header: Arc<sam::Header>,
        repository: fasta::Repository,
        entries: Arc<[cram::crai::Record]>,
        resolution: ReadGroupResolution,
        path: Arc<Path>,
    ) -> Self {
        Self {
            reader,
            header,
            repository,
            entries,
            resolution,
            path,
            next_entry: 0,
            container: None,
            served: 0,
            last_decoded_offset: None,
            other_sample_records: 0,
            // Nothing has been asked for yet, so there is nothing to read. `begin_region` is
            // what starts this reader, exactly as it restarts it.
            done: true,
        }
    }

    pub(crate) fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Records dropped at decode as another sample's — see the field.
    pub(crate) fn other_sample_records(&self) -> u64 {
        self.other_sample_records
    }

    /// Point the reader at the first container that can hold a record reaching `region.start`,
    /// and leave it able to read to the end of the contig.
    ///
    /// **Positioning early is always safe; positioning late loses reads.** The layer above
    /// discards records that do not overlap, so starting a container or two before the region
    /// costs a decode and nothing else. Starting after it drops reads silently. Every choice
    /// below leans the safe way.
    pub(crate) fn begin_region(&mut self, region: GenomeRegion) -> io::Result<()> {
        self.next_entry = self.first_entry_reaching(region.start.get());
        self.container = None;
        self.served = 0;
        self.last_decoded_offset = None;
        self.done = false;
        Ok(())
    }

    /// The first `.crai` entry whose container can hold a record reaching `start`.
    ///
    /// **A binary search and then a walk back**, and both halves are needed. The entries are
    /// sorted by the position a container *starts* at, so a binary search finds the first one
    /// starting at or after `start` — but a container beginning *before* `start` can still hold
    /// records reaching into it, and those records are ours to serve. So the search is followed
    /// by stepping back over every earlier container whose span reaches `start`.
    ///
    /// Walking from entry 0 instead would be correct and is what the per-region source does. It
    /// cannot be done here: that source is built once per region, while this one is
    /// repositioned on every jump of a walk that may make millions of them, and a prefix scan of
    /// a many-container contig would be paid at each.
    fn first_entry_reaching(&self, start: u64) -> usize {
        let container_start = |entry: &cram::crai::Record| {
            entry
                .alignment_start()
                .map_or(0, |at| usize::from(at) as u64)
        };
        let mut at = self
            .entries
            .partition_point(|entry| container_start(entry) < start);
        // Back over containers that begin earlier but reach `start`. A container with no
        // recorded start cannot be reasoned about, so it is kept rather than skipped.
        while at > 0 {
            let previous = &self.entries[at - 1];
            let begins = container_start(previous);
            let span = previous.alignment_span() as u64;
            let reaches = span == 0 || begins + span > start;
            if !reaches {
                break;
            }
            at -= 1;
        }
        at
    }

    /// The next record, raw except for its read group — see the type's doc for why that one
    /// field travels with it.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        loop {
            if let Some(container) = &self.container
                && self.served < container.len()
            {
                let i = self.served;
                self.served += 1;
                // Rebuilt into the caller's buffer rather than moved in, so the buffer's
                // allocations are reused across the whole walk.
                container.fill_record(i, &mut buf.record);
                buf.read_group = Some(container.read_group(i));
                return Ok(true);
            }

            match self.decode_next_container() {
                Ok(true) => {}
                Ok(false) => {
                    self.done = true;
                    return Ok(false);
                }
                Err(error) => {
                    self.done = true;
                    return Err(error);
                }
            }
        }
    }

    /// Decode the next container of this contig. `Ok(false)` once the entries run out.
    ///
    /// **No region is consulted.** Which containers hold records worth having is the caller's
    /// question, asked above; this walks the contig's entries in order until they are spent.
    fn decode_next_container(&mut self) -> io::Result<bool> {
        loop {
            let Some(entry) = self.entries.get(self.next_entry) else {
                return Ok(false);
            };
            self.next_entry += 1;
            let offset = entry.offset();

            // A container may hold several slices, and each slice is its own `.crai` entry
            // sharing the container's offset. Decoding it once per entry would serve every one
            // of its records again — caught loudly by the order guard above rather than
            // silently inflating depth, but wrong either way.
            //
            // **⚠ Untested, and the fixture is why.** noodles' writer puts one slice in each
            // container, so every `.crai` it produces has distinct offsets and this branch
            // never runs — deleting it leaves the whole suite green, which was checked. Files
            // written by samtools do hold multi-slice containers, so the guard is needed for
            // real input and cannot be exercised by input this project can build. Testing it
            // needs a fixture builder that can write several slices per container, or a
            // committed CRAM from another writer.
            if self.last_decoded_offset == Some(offset) {
                continue;
            }
            self.last_decoded_offset = Some(offset);

            let Some(container) = decode_container_at(
                &mut self.reader,
                &self.header,
                &self.repository,
                &self.resolution,
                offset,
            )?
            else {
                // End of stream reached through the index — nothing further.
                return Ok(false);
            };
            self.other_sample_records += container.other_sample_records();
            let has_records = container.len() > 0;
            self.container = Some(container);
            self.served = 0;
            if has_records {
                return Ok(true);
            }
            // A container of nothing but another sample's reads. Keep walking rather than
            // reporting end of input.
        }
    }
}
