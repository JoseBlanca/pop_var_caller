//! Records from a CRAM: the `.crai` says which container to start in, and reading carries on
//! from there.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles_cram as cram;
use noodles_sam as sam;

use crate::read::aligned_read::NoodlesRawAlignedRead;
use crate::read::input::aligned_reads_reader::container::{DecodedContainer, container_at};
use crate::read::input::aligned_reads_reader::prefetch::SharedDecode;
use crate::read::input::read_groups::ReadGroupResolution;
use crate::ref_seq::SliceReference;
use crate::types::GenomeRegion;

/// A CRAM reader that stays where it is between regions.
///
/// **What it yields is undecoded** — a [`NoodlesRawAlignedRead`], which is a SAM flag and a
/// mapping quality readable without unpacking anything else. The name no longer says so, so
/// this does. "Undecoded" here means *not converted into ng's own read*: a CRAM container
/// still has to be decompressed and decoded into `RecordBuf`s before any record can be looked
/// at, which is a different sense of the word and the reason to be explicit. Nothing here
/// builds an [`AlignedRead`](crate::read::AlignedRead).
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
/// here — which is why deleting that source retired those rules rather than
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
/// travels with the record, and `RegionRawAlignedReads` uses it rather than asking again.
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
    /// **The reference bases this decode reads through**, minted by the cursor factory and
    /// owned for this reader's life (`alignment_cursor.md` §10 point 1). One slice's declared
    /// span at a time — this reader holds no chromosome and no repository.
    ///
    /// Its own reader, not the cursor's: a reference reader is an open file position and a
    /// resident window, and this project gives one to every consumer of bases rather than
    /// sharing. Sharing the cursor's would mean the walk's `evict_before` could drop bases a
    /// decode is about to want — costing a re-read — and would make a same-thread re-entry
    /// possible, which `WindowedRefSeq`'s lock answers with a hang rather than a panic.
    ///
    /// **`Box<dyn …>` beside an enum whose own doc refuses dynamic dispatch, and the numbers
    /// are why.** That enum is on the per-record path — about a million virtual calls a run,
    /// which is what it declines. This reader is reached **once per slice**, and a slice holds
    /// on the order of ten thousand records, so the same run makes on the order of a hundred
    /// such calls. Making it generic instead would put a type parameter through
    /// `AlignedReadsReader` and all 72 sites that name it, for a field only this arm has.
    /// `+ Send` keeps the cursor `Send` (`a_sample_cursor_is_send_in_both_arms`).
    reference_reader: Box<dyn SliceReference + Send>,
    /// **This contig's** `.crai` entries, in file order, grouped once at open.
    entries: Arc<ContigCraiEntries>,
    /// A copy of the file's, settled at open. Owned, so this reader carries no lifetime.
    resolution: ReadGroupResolution,
    path: Arc<Path>,
    /// The next `.crai` entry to decode. Set by [`begin_region`](Self::begin_region) and
    /// advanced by reading; never reset by reading, so a walk carries on to the contig's end.
    next_entry: usize,
    /// **The decoded containers this open file's readers share, and the thread that fills them
    /// ahead of the walk** — see [`SharedDecode`].
    ///
    /// Both belong to the open file rather than to this reader, and that is the point of each:
    /// two readers walking one file's interleaved regions take a container from the cache
    /// instead of each inflating it, and the prefetcher puts the next one there before either
    /// asks for it.
    shared: SharedDecode,
    /// The container being served, and how far into it we have got.
    ///
    /// **One container, and this reader does not own it.** It is what is currently being read:
    /// a container is taken from the shared cache, drained, and let go when the next is taken.
    /// What survives between regions is the cursor's *reads*, one layer up (spec §5).
    container: Option<Arc<DecodedContainer>>,
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
        reference_reader: Box<dyn SliceReference + Send>,
        entries: Arc<ContigCraiEntries>,
        resolution: ReadGroupResolution,
        path: Arc<Path>,
        shared: SharedDecode,
    ) -> Self {
        Self {
            reader,
            header,
            reference_reader,
            entries,
            resolution,
            path,
            shared,
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

    #[cfg_attr(
        not(test),
        expect(dead_code, reason = "only tests ask a reader for its header")
    )]
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
        self.next_entry = self.entries.first_reaching(region.start.get());
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
                // allocations are reused across the whole walk. **Both halves** — the record
                // and the read group the container decided — because this arm is the one that
                // knows the group, and splitting the two left a caller able to take a record
                // and not know to ask for the rest.
                container.fill_raw_read(i, buf);
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

    /// Take the next container of this contig — from the shared cache, or by decoding it.
    /// `Ok(false)` once the entries run out.
    ///
    /// **No region is consulted.** Which containers hold records worth having is the caller's
    /// question, asked above; this walks the contig's entries in order until they are spent.
    fn decode_next_container(&mut self) -> io::Result<bool> {
        loop {
            let Some(entry) = self.entries.entries().get(self.next_entry) else {
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

            let Some(container) = container_at(
                &self.shared.containers,
                &mut self.reader,
                &self.header,
                &self.resolution,
                offset,
                &*self.reference_reader,
            )?
            else {
                // End of stream reached through the index — nothing further.
                return Ok(false);
            };
            // **Charged here, so a container the cache served counts exactly as a container
            // this reader decoded.** The tally is per reader and container-granular (see the
            // field): a reader that skipped the charge because the other one had already
            // decoded the container would silently report half the foreign records it met.
            self.other_sample_records += container.other_sample_records();
            // **Said once the container is in hand, not before.** The walk is about to spend
            // thousands of `next_locus` calls inside this container, which is the time the
            // prefetch thread has to decode the one after it. Asked earlier — before this
            // container was decoded — the hint would be overwritten by this reader's own claim
            // and buy nothing.
            self.ask_for_the_next_container();
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

impl CramAlignedReadsReader {
    /// Tell the prefetcher which container this reader expects to want after the one it has.
    ///
    /// **The next *distinct* offset**, because a container holding several slices appears under
    /// one `.crai` entry per slice and the reader steps over the repeats — hinting one of those
    /// would name the container already in hand and waste the thread's turn.
    ///
    /// A hint is only ever a hint: it is dropped if the walk repositions elsewhere, and nothing
    /// waits for it. What it must not do is send this reader's *current* container, which is why
    /// the search starts at `next_entry` and skips `last_decoded_offset`.
    fn ask_for_the_next_container(&self) {
        let Some(prefetcher) = &self.shared.prefetcher else {
            return;
        };
        let Some(next) = self.entries.entries()[self.next_entry.min(self.entries.len())..]
            .iter()
            .map(cram::crai::Record::offset)
            .find(|offset| Some(*offset) != self.last_decoded_offset)
        else {
            return;
        };
        // Already decoded, or already being decoded by somebody — the thread's turn is better
        // spent asleep than re-deriving what the cache will hand over anyway.
        if self.shared.containers.is_spoken_for(next) {
            return;
        }
        prefetcher.wants(next);
    }
}

/// **One contig's `.crai` entries, with how far the containers up to each one reach.**
///
/// The entries are in file order, which on a coordinate-sorted file is ascending by the
/// position a container *starts* at. Where a container *ends* has no such order: one long
/// or spliced read makes its container reach far past containers that start after it. So
/// positioning cannot be read off the starts alone, and [`reach_end`](Self::reach_end)
/// holds, for each entry, the furthest any container up to and including it reaches. That
/// running maximum never decreases, which is what lets one binary search answer
/// [`first_reaching`](Self::first_reaching).
///
/// **8 bytes an entry**, on top of the 56 a `crai::Record` costs: 0.9 MB for a whole-genome
/// tomato CRAM's 112,140 slices. Built once when the file is opened, from the `.crai` on
/// disk; nothing of it is written anywhere.
#[derive(Debug, Default)]
pub(crate) struct ContigCraiEntries {
    entries: Box<[cram::crai::Record]>,
    /// `reach_end[i]`: one past the last base any of `entries[..=i]` reaches.
    ///
    /// An entry with no recorded span counts as reaching only its own first base here. The
    /// span it does not record is handled by [`first_reaching`](Self::first_reaching) instead,
    /// because counting it as reaching everything would make every later jump on the contig
    /// decode from that entry onward.
    reach_end: Box<[u64]>,
}

impl ContigCraiEntries {
    /// Group one contig's entries, which must be in file order.
    pub(crate) fn new(entries: Vec<cram::crai::Record>) -> Self {
        let mut furthest = 0u64;
        let reach_end = entries
            .iter()
            .map(|entry| {
                let span = (entry.alignment_span() as u64).max(1);
                furthest = furthest.max(container_start(entry) + span);
                furthest
            })
            .collect();
        Self {
            entries: entries.into_boxed_slice(),
            reach_end,
        }
    }

    /// The entries, in file order.
    pub(crate) fn entries(&self) -> &[cram::crai::Record] {
        &self.entries
    }

    /// The first entry whose container can hold a record reaching `start`.
    ///
    /// **One binary search on the running reach, then a step back over entries with no
    /// recorded span.** Every entry before the answer reaches no further than `start - 1`,
    /// by the running maximum; an entry the index gives no span for cannot be reasoned
    /// about, so one immediately before the answer is kept rather than skipped.
    ///
    /// Until 2026-11-28 this was a binary search on the starts and a walk back that stopped
    /// at the first container not reaching `start`. A container before that one could still
    /// reach it, and its long read was then never served after a jump (code-quality review
    /// 2026-11-26, io B2).
    ///
    /// Walking from entry 0 instead would be correct and is what the per-region source did.
    /// It cannot be done here: this reader is repositioned on every jump of a walk that may
    /// make millions of them, and a prefix scan of a many-container contig would be paid at
    /// each.
    pub(crate) fn first_reaching(&self, start: u64) -> usize {
        let mut at = self.reach_end.partition_point(|&end| end <= start);
        while at > 0 && self.entries[at - 1].alignment_span() == 0 {
            at -= 1;
        }
        at
    }

    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }

    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// The position a `.crai` entry's container starts at, 0 where it records none.
fn container_start(entry: &cram::crai::Record) -> u64 {
    entry
        .alignment_start()
        .map_or(0, |at| usize::from(at) as u64)
}

/// [`ContigCraiEntries::first_reaching`] over hand-built entries, for the tests.
#[cfg(test)]
fn first_entry_reaching(entries: &[cram::crai::Record], start: u64) -> usize {
    ContigCraiEntries::new(entries.to_vec()).first_reaching(start)
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

    /// **A container that reaches the region from behind a nearer one that does not.** The
    /// first entry starts at 1,000 and reaches 50,999 (one long or spliced read); the second
    /// covers 1,100 to 1,399. Walking back from the first entry starting past 20,000 stops at
    /// the second, which does not reach 20,000, so the first was never decoded and its long
    /// read was lost.
    #[test]
    fn a_container_reaching_the_region_behind_one_that_does_not_is_stepped_back_to() {
        let entries = vec![entry(1_000, 50_000), entry(1_100, 300)];

        assert_eq!(first_entry_reaching(&entries, 20_000), 0);
        assert_eq!(
            first_entry_reaching(&entries, 51_000),
            entries.len(),
            "and past every container's reach, nothing is decoded",
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
