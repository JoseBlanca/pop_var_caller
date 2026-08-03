//! Records from a CRAM: the `.crai` says which container to start in, and reading carries on
//! from there.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::read::input::aligned_reads_reader::container::{
    DecodedContainer, decode_container_at,
};
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::types::GenomeRegion;

/// A CRAM reader that stays where it is between regions.
///
/// **What it yields is undecoded** — a [`NoodlesRawAlignedRead`], which is a SAM flag and a
/// mapping quality readable without unpacking anything else. The name no longer says so, so
/// this does. "Undecoded" here means *not converted into ng's own read*: a CRAM container
/// still has to be decompressed and decoded into `RecordBuf`s before any record can be looked
/// at, which is a different sense of the word and the reason to be explicit. Nothing here
/// builds an [`AlignedRead`](crate::ng::read::AlignedRead).
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
/// the region it was handed. That is the contract in `aligned_reads_reader/mod.rs`, and it is
/// load-bearing rather than tidy: the cursor above serves a forward region by *not
/// repositioning at all*, so a reader that had quietly stopped at the previous region's end
/// would lose every record past it, silently, for every region after the first.
///
/// **This is where the per-region CRAM source and this one differed most**, and it was not a
/// re-layering: that source stopped its `.crai` walk at the first container beginning past the
/// region's end, and filtered each container's records against the region. Both behaviours
/// belong to a reader that answers one region and is then thrown away, and neither can be
/// here — which is why Milestone F deleting that source retired those rules rather than
/// leaving them for this arm to inherit.
///
/// # Each record carries its own read group — the documented exception
///
/// Every other arm hands up a record with no read group attached, and the layer above resolves
/// it from the record's `RG` tag. This arm attaches it, because a CRAM does not have an `RG`
/// tag to resolve: it stores the read group as a **number**, an index into the header's `@RG`
/// list. Deciding it at decode is what lets every auxiliary tag be dropped — 6.2 MiB per open
/// file, measured, the largest single item in what an open file cost — and re-inflating a
/// number into a string so the layer above can parse it back would be perverse. So the answer
/// travels with the record, and `RegionRecords` uses it rather than asking again.
///
/// **A sample is not one read group, and this is the arm where getting that wrong would be
/// worst.** One individual sequenced as several libraries declares an `@RG` for each, and every
/// one of them is ours — the reads are all kept, each carrying the library it came from,
/// because that is the grouping a per-chemistry error model is fitted on. The *only* thing that
/// makes a record foreign is belonging to another **sample**. Since this is the one arm that
/// decides the question itself, a mistake here would drop a whole library before any other
/// layer could see it, and the result would look like a sample sequenced less deeply rather
/// than like a fault.
pub(crate) struct CramAlignedReadsReader {
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
impl std::fmt::Debug for CramAlignedReadsReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramAlignedReadsReader")
            .field("path", &self.path)
            .field("entries", &self.entries.len())
            .field("next_entry", &self.next_entry)
            .field("done", &self.done)
            .finish_non_exhaustive()
    }
}

impl CramAlignedReadsReader {
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
        self.next_entry = first_entry_reaching(&self.entries, region.start.get());
        // **What is held is what is being served, and it belongs to where the reader was.**
        // Dropping it here is what makes a reposition a reposition: a reader that kept it
        // would serve the rest of the old container before reaching the new position, which
        // is another part of the contig's reads under this region's name.
        self.container = None;
        self.served = 0;
        self.last_decoded_offset = None;
        self.done = false;
        Ok(())
    }

    /// The next record, raw except for its read group — see the type's doc for why that one
    /// field travels with it.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawAlignedRead) -> io::Result<bool> {
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

/// The first `.crai` entry whose container can hold a record reaching `start`.
///
/// **A binary search and then a walk back**, and both halves are needed. The entries are
/// sorted by the position a container *starts* at, so a binary search finds the first one
/// starting at or after `start` — but a container beginning *before* `start` can still hold
/// records reaching into it, and those records are ours to serve. So the search is followed
/// by stepping back over every earlier container whose span reaches `start`.
///
/// Walking from entry 0 instead would be correct and is what the per-region source did. It
/// cannot be done here: that source was built once per region, while this one is repositioned
/// on every jump of a walk that may make millions of them, and a prefix scan of a
/// many-container contig would be paid at each.
///
/// A free function over the entries rather than a method, so the rule can be stated against a
/// hand-built index — a `CramAlignedReadsReader` needs a real file behind it, and the property this
/// has to get right is about the index alone.
fn first_entry_reaching(entries: &[cram::crai::Record], start: u64) -> usize {
    let container_start = |entry: &cram::crai::Record| {
        entry
            .alignment_start()
            .map_or(0, |at| usize::from(at) as u64)
    };
    let mut at = entries.partition_point(|entry| container_start(entry) < start);
    // Back over containers that begin earlier but reach `start`. A container with no
    // recorded start cannot be reasoned about, so it is kept rather than skipped.
    while at > 0 {
        let previous = &entries[at - 1];
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A `.crai` entry for a container covering `[start, start + span - 1]`.
    fn entry(start: usize, span: usize) -> cram::crai::Record {
        cram::crai::Record::new(Some(0), noodles_core::Position::new(start), span, 0, 0, 0)
    }

    /// Ten containers of 100 bases each, back to back from base 1.
    fn ten_containers() -> Vec<cram::crai::Record> {
        (0..10).map(|i| entry(1 + i * 100, 100)).collect()
    }

    /// **The container-level skip — and this is the assertion that discriminates.** A region
    /// deep in a contig must be reached by *landing* on its container, not by walking every
    /// earlier one. The reads would be identical either way, so nothing but an explicit
    /// assertion on where the walk starts can see the difference — and the difference is a
    /// prefix scan of the contig's whole index, paid once per region, ~10⁶ times.
    ///
    /// Inherited from the per-region source's `the_crai_walk_skips_containers_that_end_before
    /// _the_region`, which asserted the same thing about a walk that started at entry 0.
    #[test]
    fn containers_ending_before_the_region_are_skipped_rather_than_walked() {
        let entries = ten_containers();

        // Base 750 is inside the eighth container (601..=700 is the seventh, 701..=800 the
        // eighth), so the seven before it must be stepped over, not visited.
        assert_eq!(first_entry_reaching(&entries, 750), 7);
        assert_eq!(first_entry_reaching(&entries, 1), 0);
        assert_eq!(
            first_entry_reaching(&entries, 10_000),
            entries.len(),
            "past every container, the walk starts past the end and reads nothing",
        );
    }

    /// **The walk back, which is the half that loses reads if it is missing.** A container
    /// beginning *before* the region can still hold a record reaching into it — so landing on
    /// the first container starting at or after the region drops exactly those reads,
    /// silently. Positioning early costs a decode; positioning late costs reads.
    #[test]
    fn a_container_beginning_earlier_but_reaching_the_region_is_stepped_back_to() {
        // The second container starts at 101 and reaches to 400 — a long read, or a container
        // whose records were written out of span order.
        let entries = vec![entry(1, 50), entry(101, 300), entry(401, 100)];

        assert_eq!(
            first_entry_reaching(&entries, 350),
            1,
            "the container starting at 101 reaches 400, so it holds records touching 350",
        );
        assert_eq!(
            first_entry_reaching(&entries, 401),
            2,
            "and one that ends at 400 does not reach 401",
        );
    }

    /// A container the index gives no span for cannot be reasoned about, so it is kept rather
    /// than skipped — the safe direction, since positioning early only costs a decode.
    #[test]
    fn a_container_with_no_recorded_span_is_kept_rather_than_skipped() {
        let entries = vec![entry(1, 0), entry(101, 100)];

        assert_eq!(first_entry_reaching(&entries, 150), 0);
    }

    /// An empty index has nothing to reach, whatever is asked of it.
    #[test]
    fn an_empty_index_reaches_nothing() {
        assert_eq!(first_entry_reaching(&[], 1), 0);
        assert_eq!(first_entry_reaching(&[], 10_000), 0);
    }
}
