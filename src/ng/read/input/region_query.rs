//! Serving one region: the BAM and CRAM region-query record sources, and the
//! order guard over the filtered stream.
//!
//! Both sources are index-driven siblings of the whole-file
//! [`BamRecordSource`](crate::ng::read::filtering::BamRecordSource) /
//! [`CramRecordSource`](crate::ng::read::filtering::CramRecordSource): they
//! query the **already-parsed** index for candidate chunks, seek a pooled
//! reader, and drop records the index over-returned — chunk-edge slop, which is
//! a reader concern and so is dropped **uncounted**, never charged to a filter
//! drop reason. There is no new trait: BAM and CRAM are two containers for one
//! idea, not competing implementations to bake off
//! (`doc/devel/ng/arch/alignment_file.md` §4).
//!
//! Two things here fail *quietly* rather than loudly, which is why each is
//! built and committed on its own with an independent oracle:
//!
//! - **A missed chunk edge is a wrong genotype, not a crash.** The region
//!   query's oracle is the existing whole-file source: an indexed query must
//!   return exactly what a full linear scan filtered to the same region returns
//!   (spec §7, T5).
//! - **A guard that never fires looks exactly like a guard that works.** The
//!   order check is mutation-verified — removing it must let a planted
//!   regression through (T4a).
//!
//! The order guard's state lives in the per-region iterator, never on the
//! handle, so querying region B and then region A is a new forward scan rather
//! than a spurious regression (spec §3.2).

use std::fs::File;
use std::io;
use std::io::SeekFrom;
use std::iter::FusedIterator;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position as RecordPosition;
use noodles_cram as cram;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::{Flags, MappingQuality};

use crate::bam::index_preflight::AlignmentIndex;
use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::filtering::{
    NoodlesRawRecord, ReadFilterError, RecordSource, resolve_read_group,
};
use crate::ng::read::input::read_groups::{ReadGroupResolution, RecordOwner};
use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position};

use super::AlignmentFileError;

/// The reads of one region of a BAM, as a [`RecordSource`] the step-1 filter
/// can drive.
///
/// Index-driven: the already-parsed index gives candidate chunks, each is
/// seeked in turn, and records are read one at a time into the caller's reused
/// buffer. The index answers at *bin* granularity, so it over-returns at the
/// edges — those records are dropped **uncounted**, because they are not reads
/// the filter rejected, they are reads the index was never asked about.
///
/// Owns its reader rather than borrowing one, so the whole chain can be built
/// by value and the reader handed back to the pool when the stream ends
/// ([`Self::into_reader`]).
pub(crate) struct BamRegionSource {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    /// Shared with the `AlignmentFile`: parsed once, never per query, and
    /// **owned** rather than borrowed so this source does not tie the stream it
    /// serves to the file's lifetime (arch `locus_generation_pileup.md` §2.2).
    header: Arc<sam::Header>,
    chunks: std::vec::IntoIter<Chunk>,
    /// Where the chunk being scanned ends, or `None` to take the next chunk.
    current_chunk_end: Option<bgzf::VirtualPosition>,
    /// The region's contig as a `ref_id`. Equal to `ContigId` by the open
    /// gate, which is what makes this a cast rather than a lookup.
    target_reference_sequence_id: usize,
    region: GenomeRegion,
    /// A copy of the `AlignmentFile`'s, settled at open and never rebuilt per
    /// query. Held by value, not behind an `Arc`: cloning one is a tag and a
    /// `u32` for the common single-read-group file, and one atomic increment
    /// otherwise (`ReadGroupResolution`). Owned, so this source carries no
    /// lifetime.
    resolution: ReadGroupResolution,
    /// Records skipped as another sample's — see `RecordSource::other_sample_records`.
    other_sample_records: u64,
    /// Latched once the scan passes the region end or the chunks run out.
    done: bool,
}

/// What a region query resolves to before any reader is involved: which
/// contig, and which chunks of the file to scan.
///
/// Split out from the source so **every fallible part of setting up a query
/// happens before a reader is taken from the pool**. Taking the handle
/// transfers the obligation to return it; a `?` after that point would lose a
/// reader silently (`BorrowedReader::take`).
pub(crate) struct RegionPlan {
    chunks: Vec<Chunk>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
}

impl BamRegionSource {
    /// Resolve the region and query the index for its chunks.
    ///
    /// The index query happens here, once per region — an in-memory lookup on
    /// the index parsed at open, never a re-parse (spec §3.3). No reader is
    /// touched, so a failure costs nothing to recover from.
    pub(crate) fn plan(
        header: &sam::Header,
        index: &AlignmentIndex,
        region: GenomeRegion,
        path: &Path,
    ) -> Result<RegionPlan, AlignmentFileError> {
        let invalid_region = || AlignmentFileError::Region { region };

        let target_reference_sequence_id =
            usize::try_from(region.contig.get()).map_err(|_| invalid_region())?;
        if target_reference_sequence_id >= header.reference_sequences().len() || region.is_empty() {
            return Err(invalid_region());
        }

        let interval = interval_of(region).ok_or_else(invalid_region)?;
        let chunks = match index {
            AlignmentIndex::BamCsi(index) => index.query(target_reference_sequence_id, interval),
            AlignmentIndex::BamBai(index) => index.query(target_reference_sequence_id, interval),
            // `open` pairs a `.crai` only with a CRAM, and CRAM goes through
            // `CramRegionSource`. Reaching here means the handle was built
            // wrong, which is this module's bug and not the caller's — so it
            // must not masquerade as bad input.
            AlignmentIndex::Crai(_) => unreachable!(
                "a .crai index reached the BAM region source; \
                 open pairs .crai only with CRAM"
            ),
        }
        // The region was checked against the contig table above, so a failure
        // here is the index itself being unusable, not the query being silly.
        .map_err(|source| AlignmentFileError::IndexQuery {
            path: path.to_path_buf(),
            region,
            source,
        })?;

        Ok(RegionPlan {
            chunks,
            target_reference_sequence_id,
            region,
        })
    }

    /// Attach a reader to a plan. **Infallible**, which is what makes it safe
    /// to call after the pooled handle has been taken.
    pub(crate) fn new(
        reader: bam::io::Reader<bgzf::io::Reader<File>>,
        header: Arc<sam::Header>,
        plan: RegionPlan,
        resolution: ReadGroupResolution,
    ) -> Self {
        Self {
            reader,
            header,
            chunks: plan.chunks.into_iter(),
            current_chunk_end: None,
            target_reference_sequence_id: plan.target_reference_sequence_id,
            region: plan.region,
            resolution,
            other_sample_records: 0,
            done: false,
        }
    }

    /// Give the reader back, for the pool.
    pub(crate) fn into_reader(self) -> bam::io::Reader<bgzf::io::Reader<File>> {
        self.reader
    }
}

/// The region as a noodles query interval. `None` if the 1-based bounds cannot
/// be represented (a `0` start, or a value past `usize`).
fn interval_of(region: GenomeRegion) -> Option<noodles_core::region::Interval> {
    let start = noodles_core::Position::new(usize::try_from(region.start.get()).ok()?)?;
    let end = noodles_core::Position::new(usize::try_from(region.end.get()).ok()?)?;
    Some((start..=end).into())
}

impl RecordSource for BamRegionSource {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        &self.header
    }

    fn other_sample_records(&self) -> u64 {
        self.other_sample_records
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        buf.read_group = None;
        loop {
            // 1. Land inside a chunk, seeking to its start if we have just
            //    taken it.
            let chunk_end = match self.current_chunk_end {
                Some(end) => end,
                None => {
                    let Some(chunk) = self.chunks.next() else {
                        self.done = true;
                        return Ok(false);
                    };
                    self.reader
                        .get_mut()
                        .seek(chunk.start())
                        .inspect_err(|_| self.done = true)?;
                    self.current_chunk_end = Some(chunk.end());
                    chunk.end()
                }
            };

            // 2. The previous read may have carried us past this chunk's end.
            if self.reader.get_ref().virtual_position() >= chunk_end {
                self.current_chunk_end = None;
                continue;
            }

            // 3. One record, into the caller's reused buffer.
            let bytes_read = self
                .reader
                .read_record_buf(&self.header, &mut buf.record)
                .inspect_err(|_| self.done = true)?;
            if bytes_read == 0 {
                // End of file inside a chunk — try the next one.
                self.current_chunk_end = None;
                continue;
            }

            let on_target_contig =
                buf.record.reference_sequence_id() == Some(self.target_reference_sequence_id);

            // 4. The sorted early stop: once a read *on the target contig*
            //    starts past the region, nothing later in this scan can
            //    overlap, because the file is coordinate-sorted and the open
            //    gate proved it claims to be.
            if on_target_contig
                && let Some(start) = buf.record.alignment_start()
                && usize::from(start) as u64 > self.region.end.get()
            {
                self.done = true;
                return Ok(false);
            }

            // 5. Drop what the index over-returned — a different contig in a
            //    chunk that straddles one, or a footprint that misses the
            //    region. Uncounted: these are a reader's business, not a
            //    filter's, and charging them to a `DropReason` would make the
            //    tally mean something different for an indexed read than for a
            //    whole-file one.
            if !on_target_contig || !overlaps(&buf.record, self.region) {
                continue;
            }

            // Resolved here, on the record that is actually being yielded: the
            // read group can depend on the record's own `RG` tag, so resolving
            // before the read — or on one the loop above skipped — would
            // attribute this read to whatever the previous record carried. It
            // also means the records this query discards never pay the lookup.
            match resolve_read_group(&buf.record, &self.resolution)? {
                RecordOwner::Mine(id) => {
                    buf.read_group = Some(id);
                    return Ok(true);
                }
                // Another sample's read, in a file this sample shares. Counted
                // apart from every drop — it says nothing about how this read
                // group behaved — and skipped before the filter sees it.
                RecordOwner::OtherSample => {
                    self.other_sample_records += 1;
                    continue;
                }
            }
        }
    }
}

/// Which `.crai` entries a CRAM region query has to walk, resolved before any
/// reader is involved. The sibling of [`RegionPlan`], and split out for the
/// same reason: nothing fallible may happen after a reader leaves the pool.
pub(crate) struct CramRegionPlan {
    /// Just this contig's entries, in file order — see
    /// [`crate::ng::read::input::open_bam::group_crai_by_contig`] for why the
    /// grouping is done once at open rather than searched per query.
    entries: Arc<[cram::crai::Record]>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
}

/// One CRAM container, decoded, kept for the *next* query.
///
/// **Why this exists.** CRAM decodes at container granularity, and a container holds thousands
/// of records covering a wide span of the reference, so consecutive region queries land in the
/// same container over and over. Decoding it per query costs the block decompression *and*
/// noodles' per-slice reference MD5 verification, which hashes the slice's whole reference span
/// on every `Slice::records` call — measured at 35% of a per-locus STR walk over a CRAM
/// (`doc/devel/reports/reviews/perf_ng-ssr-cohort-stutter_2026-07-26.md`, H2). Holding the last
/// container decoded through a given reader turns that per-query cost into a per-container one.
///
/// It travels with the pooled [`ReaderHandle`](super::open_bam::ReaderKind), so it is **per
/// worker**: each concurrent caller holds its own reader and therefore its own container, and no
/// lock is added to the query path. That is the shape `arch/alignment_file.md` §7 anticipated.
///
/// **The key is the `.crai` offset**, which identifies a container within the file. Reusing on a
/// stale key would serve another region's reads as this region's — silently wrong depth, not a
/// crash — so the offset comparison is the whole safety argument, and
/// `a_stale_container_is_not_reused` pins it.
/// ## Why it is stored packed rather than as records
///
/// A container is ~10⁴ records, and a `RecordBuf` per record was ~680 bytes across seven separate
/// allocations — against ~270 bytes of read that anything actually consumes (name, sequence,
/// qualities, CIGAR). The other 60% was `Vec` headers, capacity slack and per-allocation
/// overhead, ~72,000 allocations per open file, resident for the whole run.
///
/// So the bytes go into two flat buffers, [`payload`](Self::payload) and
/// [`cigar_ops`](Self::cigar_ops), and each record becomes a fixed-size [`RecordIndex`] naming
/// its slices of them. A query materialises a `RecordBuf` **only for the records it actually
/// serves**, straight into the caller's reused buffer — which is also how the per-query
/// `record.clone()` this replaces disappeared, so the query path now allocates nothing per read.
pub(crate) struct DecodedContainer {
    /// The `.crai` offset the records were decoded from.
    offset: u64,
    /// One entry per record of the container, in file order — *not* filtered to any region,
    /// because the next query filters against its own.
    index: Vec<RecordIndex>,
    /// Every record's name, sequence and quality scores, back to back. Sliced by the offsets in
    /// [`RecordIndex`]; meaningless on its own.
    payload: Vec<u8>,
    /// Every record's CIGAR operations, back to back. Separate from [`payload`](Self::payload)
    /// because an `Op` is not a byte and packing it into one would mean encoding and decoding it.
    cigar_ops: Vec<Op>,
}

/// Where one record's bytes are, and every scalar field of it that anything reads.
///
/// Fixed size, so the whole index is one allocation. The scalars are the noodles types rather
/// than raw integers: they are all `Copy`, so nothing is gained by unpacking them, and keeping
/// them means [`DecodedContainer::fill`] hands them back without conversion.
#[derive(Clone, Copy)]
struct RecordIndex {
    /// Where the record sits on the reference — the only thing a query needs before deciding
    /// whether to serve it at all.
    ///
    /// A record's footprint cannot change after it is decoded, but the region it is asked about
    /// changes with every query — and a container is re-filtered on every query that reuses it
    /// (one per locus, so ~5,500 per sample per chromosome on the STR path against a container of
    /// ~10⁴ records). Recomputing it there meant a fresh CIGAR walk per record per query, which a
    /// sampling profile put at 9.4% of a cohort run.
    footprint: Footprint,
    /// Which read group this record belongs to, resolved **once per decode**.
    ///
    /// The read group is the *only* thing anything reads out of a record's auxiliary tags
    /// (`resolve_read_group`, and `AlignedRead` carries no tag field), so resolving it while the
    /// tags are in hand is what lets them be dropped — 6.2 MiB per open file, measured, the
    /// single largest item in what an open file cost. Keeping a whole tag map to reach one field
    /// of it was the waste.
    owner: RecordOwner,
    flags: Flags,
    mapping_quality: Option<MappingQuality>,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<RecordPosition>,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<RecordPosition>,
    template_length: i32,
    /// `None` for a record with no name at all, which is not the same as an empty one.
    name: Option<Span>,
    sequence: Span,
    quality_scores: Span,
    /// Indexes [`DecodedContainer::cigar_ops`]; the others index `payload`.
    cigar: Span,
}

/// A half-open range into one of the container's flat buffers, in elements.
///
/// `u32` rather than `usize`: a container is bounded by the CRAM writer's records-per-slice
/// (10,000 in every file this project has met), so neither buffer approaches 4 GiB, and halving
/// the index's width matters more than the headroom does.
#[derive(Clone, Copy)]
struct Span {
    start: u32,
    end: u32,
}

impl Span {
    /// **Refuses rather than wraps.** A container whose payload passed 4 GiB would silently
    /// truncate these offsets and hand back another record's bytes — wrong reads, not a crash.
    /// No CRAM this project has met comes within three orders of magnitude of it (10,000 records
    /// a container, ~270 bytes a record), but the spec's record count is a 32-bit field, so the
    /// bound is the index's and has to be checked rather than assumed.
    fn new(start: usize, end: usize) -> io::Result<Self> {
        match (u32::try_from(start), u32::try_from(end)) {
            (Ok(start), Ok(end)) => Ok(Self { start, end }),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "a CRAM container holds more than 4 GiB of read data, which this reader's \
                 per-container index cannot address",
            )),
        }
    }

    fn range(self) -> std::ops::Range<usize> {
        self.start as usize..self.end as usize
    }
}

impl DecodedContainer {
    /// Append one decoded record: its bytes move into the flat
    /// buffers and its scalars into the index, and the `RecordBuf` itself is dropped by the
    /// caller.
    ///
    /// The transient `RecordBuf` is still built — noodles' conversion has no field-selective
    /// form, so the allocations happen either way. What changes is that they do not *survive*:
    /// the churn is unchanged, the residency is not.
    fn push(&mut self, record: &RecordBuf, owner: RecordOwner) -> io::Result<()> {
        let name = record
            .name()
            .map(|name| self.append_bytes(AsRef::<[u8]>::as_ref(name)))
            .transpose()?;
        let sequence = self.append_bytes(record.sequence().as_ref())?;
        let quality_scores = self.append_bytes(record.quality_scores().as_ref())?;
        let cigar_start = self.cigar_ops.len();
        self.cigar_ops.extend_from_slice(record.cigar().as_ref());
        let cigar = Span::new(cigar_start, self.cigar_ops.len())?;

        self.index.push(RecordIndex {
            footprint: Footprint::of(record),
            owner,
            flags: record.flags(),
            mapping_quality: record.mapping_quality(),
            reference_sequence_id: record.reference_sequence_id(),
            alignment_start: record.alignment_start(),
            mate_reference_sequence_id: record.mate_reference_sequence_id(),
            mate_alignment_start: record.mate_alignment_start(),
            template_length: record.template_length(),
            name,
            sequence,
            quality_scores,
            cigar,
        });
        Ok(())
    }

    /// Give back every byte the buffers over-reserved while growing.
    fn shrink_to_fit(&mut self) {
        self.index.shrink_to_fit();
        self.payload.shrink_to_fit();
        self.cigar_ops.shrink_to_fit();
    }

    fn append_bytes(&mut self, bytes: &[u8]) -> io::Result<Span> {
        let start = self.payload.len();
        self.payload.extend_from_slice(bytes);
        Span::new(start, self.payload.len())
    }

    /// Rebuild record `i` **into `out`**, reusing its allocations.
    ///
    /// Every owned field is cleared and refilled rather than replaced, so a stream that serves a
    /// million reads through one buffer allocates for the longest read it meets and nothing
    /// after — the same buffer-reuse property the BAM path gets from noodles'
    /// `read_record_buf`, which the CRAM path used to give up by moving whole records in.
    fn fill(&self, i: usize, out: &mut RecordBuf) {
        let entry = &self.index[i];

        *out.flags_mut() = entry.flags;
        *out.reference_sequence_id_mut() = entry.reference_sequence_id;
        *out.alignment_start_mut() = entry.alignment_start;
        *out.mapping_quality_mut() = entry.mapping_quality;
        *out.mate_reference_sequence_id_mut() = entry.mate_reference_sequence_id;
        *out.mate_alignment_start_mut() = entry.mate_alignment_start;
        *out.template_length_mut() = entry.template_length;

        match entry.name {
            Some(span) => {
                let name = out.name_mut().get_or_insert_with(Default::default);
                name.clear();
                name.extend_from_slice(&self.payload[span.range()]);
            }
            None => *out.name_mut() = None,
        }

        let sequence = out.sequence_mut().as_mut();
        sequence.clear();
        sequence.extend_from_slice(&self.payload[entry.sequence.range()]);

        let quality_scores = out.quality_scores_mut().as_mut();
        quality_scores.clear();
        quality_scores.extend_from_slice(&self.payload[entry.quality_scores.range()]);

        let cigar = out.cigar_mut().as_mut();
        cigar.clear();
        cigar.extend_from_slice(&self.cigar_ops[entry.cigar.range()]);

        // The container carries no auxiliary tags — they are dropped at decode, once the read
        // group has been resolved from them. Cleared rather than assumed empty because the buffer
        // is the caller's and outlives any one record: if it ever reached here carrying another
        // record's tags, those tags would silently become this record's.
        out.data_mut().clear();
    }
}

/// One record's reference footprint — what [`overlaps`] needs and all it needs.
///
/// `start`/`end` are 1-based inclusive, matching [`GenomeRegion`]. An unmapped or
/// position-less record gets `reference_sequence_id: None` and never overlaps anything, which is
/// exactly what the `RecordBuf`-shaped test it replaces did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Footprint {
    reference_sequence_id: Option<usize>,
    start: u64,
    end: u64,
}

impl Footprint {
    /// The footprint of `record`, one CIGAR walk.
    fn of(record: &sam::alignment::RecordBuf) -> Self {
        match (record.alignment_start(), record.alignment_end()) {
            (Some(first), Some(last)) => Self {
                reference_sequence_id: record.reference_sequence_id(),
                start: usize::from(first) as u64,
                end: usize::from(last) as u64,
            },
            _ => Self {
                reference_sequence_id: None,
                start: 0,
                end: 0,
            },
        }
    }

    /// Whether this footprint is on `contig` and touches `region`. The same test [`overlaps`]
    /// applies, against pre-resolved coordinates.
    fn overlaps(&self, contig: usize, region: GenomeRegion) -> bool {
        self.reference_sequence_id == Some(contig)
            && self.start <= region.end.get()
            && self.end >= region.start.get()
    }
}

/// The reads of one region of a CRAM.
///
/// CRAM decodes at **container** granularity — a container holds slices whose
/// records are decoded together against the reference — so this walks the
/// `.crai` for containers on the target contig, decodes each (or reuses the one
/// [`DecodedContainer`] this reader already holds), and buffers the overlapping
/// records. Non-overlapping ones are dropped **uncounted**, as in the BAM
/// source, and for the same reason.
///
/// **The `.crai` walk starts from a binary search, not from entry 0.**
/// Production rescans the index from the beginning on every call
/// (`segment_reader.rs`'s `fetch_mapped_reads`), which for a late contig in a
/// many-contig `.crai` is a repeated O(n) prefix scan — paid once per locus,
/// ~10⁶ times on the STR path. The `.crai` is sorted by
/// `(reference_sequence_id, alignment_start)`, so the contig's first entry is
/// one binary search away, and the forward walk from there is bounded by the
/// container-level early stop below (spec §3.3).
pub(crate) struct CramRegionSource {
    reader: cram::io::Reader<File>,
    /// Shared with the `AlignmentFile`, for the reason its BAM sibling's is
    /// ([`BamRegionSource`]): owned, so the stream this serves outlives the
    /// borrow that built it.
    header: Arc<sam::Header>,
    /// The reference bases slice decoding consults. Cheap to clone — it is
    /// internally shared — and cloned per decode, as noodles requires.
    repository: fasta::Repository,
    /// This contig's `.crai` entries, in file order.
    entries: Arc<[cram::crai::Record]>,
    /// Cursor into [`Self::entries`].
    next_index_record: usize,
    /// Indices into [`Self::container`] of the records overlapping this query, drained one per
    /// `read_next` before the next container is touched. Indices rather than records: the bytes
    /// stay in the container until a record is actually served.
    pending: std::vec::IntoIter<u32>,
    /// The offset last decoded **in this query**, so a container that appears
    /// once per slice is drained once. Distinct from `container`'s key: this one
    /// exists to stop the same records being yielded twice *within* a query,
    /// while `container` exists to stop them being *decoded* twice *across*
    /// queries. Resetting this per query is what makes a new query re-serve a
    /// container it already holds.
    last_decoded_offset: Option<u64>,
    /// The container this reader last decoded, reused when the walk asks for it
    /// again ([`DecodedContainer`]). Moved in from the pooled handle and moved
    /// back out by [`Self::into_parts`].
    container: Option<DecodedContainer>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
    /// A copy of the `AlignmentFile`'s, settled at open and never rebuilt per
    /// query. Held by value, not behind an `Arc`: cloning one is a tag and a
    /// `u32` for the common single-read-group file, and one atomic increment
    /// otherwise (`ReadGroupResolution`). Owned, so this source carries no
    /// lifetime.
    resolution: ReadGroupResolution,
    /// Records skipped as another sample's — see `RecordSource::other_sample_records`.
    other_sample_records: u64,
    done: bool,
}

impl CramRegionSource {
    /// Pick out the target contig's `.crai` entries. **O(1)** — the grouping
    /// was done once, at open (`open_bam::group_crai_by_contig`, which also
    /// explains why this is not a binary search).
    pub(crate) fn plan(
        header: &sam::Header,
        crai_by_contig: &[Arc<[cram::crai::Record]>],
        region: GenomeRegion,
    ) -> Result<CramRegionPlan, AlignmentFileError> {
        let invalid_region = || AlignmentFileError::Region { region };

        let target_reference_sequence_id =
            usize::try_from(region.contig.get()).map_err(|_| invalid_region())?;
        if target_reference_sequence_id >= header.reference_sequences().len() || region.is_empty() {
            return Err(invalid_region());
        }

        let entries = crai_by_contig
            .get(target_reference_sequence_id)
            .cloned()
            .ok_or_else(invalid_region)?;

        Ok(CramRegionPlan {
            entries,
            target_reference_sequence_id,
            region,
        })
    }

    /// Attach a reader to a plan. **Infallible**, so it is safe after the
    /// pooled handle has been taken.
    pub(crate) fn new(
        reader: cram::io::Reader<File>,
        header: Arc<sam::Header>,
        repository: fasta::Repository,
        plan: CramRegionPlan,
        resolution: ReadGroupResolution,
        container: Option<DecodedContainer>,
    ) -> Self {
        Self {
            reader,
            header,
            repository,
            entries: plan.entries,
            next_index_record: 0,
            pending: Vec::new().into_iter(),
            last_decoded_offset: None,
            container,
            target_reference_sequence_id: plan.target_reference_sequence_id,
            region: plan.region,
            resolution,
            other_sample_records: 0,
            done: false,
        }
    }

    /// Give the reader and the decoded container back, for the pool.
    pub(crate) fn into_parts(self) -> (cram::io::Reader<File>, Option<DecodedContainer>) {
        (self.reader, self.container)
    }

    /// Walk the `.crai` to the next container that could overlap, decode it,
    /// and buffer the records that do. `Ok(false)` once the walk is done.
    fn refill(&mut self) -> io::Result<bool> {
        loop {
            let Some(record) = self.entries.get(self.next_index_record).cloned() else {
                return Ok(false);
            };
            self.next_index_record += 1;

            debug_assert_eq!(
                record.reference_sequence_id(),
                Some(self.target_reference_sequence_id),
                "the grouping at open put only this contig's entries here"
            );

            // Container-level narrowing. Stopping once a container *starts*
            // past the region is what keeps a tiny locus near the start of a
            // large contig from walking that contig's whole index tail.
            if let Some(start) = record.alignment_start() {
                let container_start = usize::from(start) as u64;
                if container_start > self.region.end.get() {
                    return Ok(false);
                }
                let span = record.alignment_span() as u64;
                if span > 0 && container_start + span - 1 < self.region.start.get() {
                    continue;
                }
            }

            // `.crai` offsets are *container* positions and a container may
            // hold several slices, so a multi-slice container appears as
            // several entries sharing one offset. Decoding it once per entry
            // would surface every record again — caught loudly by the order
            // guard rather than silently inflating depth, but wrong either way.
            if self.last_decoded_offset == Some(record.offset()) {
                continue;
            }
            self.last_decoded_offset = Some(record.offset());

            // Decode only if this reader is not already holding this container. The two steps are
            // sequenced rather than written as one `match` because decoding needs `&mut self`
            // while reading the held records needs `&self`.
            let offset = record.offset();
            if self
                .container
                .as_ref()
                .is_none_or(|held| held.offset != offset)
            {
                let Some(container) = self.decode_container_at(offset)? else {
                    // End of stream reached through the index — nothing further.
                    return Ok(false);
                };
                self.container = Some(container);
            }
            let held = self
                .container
                .as_ref()
                .expect("just decoded, or already held at this offset");

            // Indices, not records: the bytes stay in the container and only the records this
            // query actually serves are ever rebuilt, one at a time, into the caller's buffer.
            let overlapping: Vec<u32> = held
                .index
                .iter()
                .enumerate()
                .filter(|(_, entry)| {
                    entry
                        .footprint
                        .overlaps(self.target_reference_sequence_id, self.region)
                })
                .map(|(i, _)| i as u32)
                .collect();

            if !overlapping.is_empty() {
                self.pending = overlapping.into_iter();
                return Ok(true);
            }
            // The container was in range but held nothing overlapping — keep
            // walking rather than reporting end of input.
        }
    }

    /// Seek to a container and decode every slice into owned records, each paired with the read
    /// group it belongs to. `Ok(None)` at end of stream (`read_container` reads 0 — the EOF
    /// marker).
    fn decode_container_at(&mut self, offset: u64) -> io::Result<Option<DecodedContainer>> {
        self.reader.seek(SeekFrom::Start(offset))?;

        let mut container = cram::io::reader::Container::default();
        if self.reader.read_container(&mut container)? == 0 {
            return Ok(None);
        }

        // **The auxiliary tags do not survive this function.** A `RecordBuf` owns every
        // `@RG`/`MD`/`NM`/… field the record carried, and this container is held for the life of
        // the file — so those fields would be resident for the whole run, once per open file.
        // Measured on the tomato cohort they are 6.2 MiB of the 12.7 MiB an open file costs: the
        // single largest item (`examples/dhat_ng_open_files.rs`).
        //
        // Nothing in ng reads a tag except `resolve_read_group`, and only the `RG` one.
        // `AlignedRead`, the only thing a record ever becomes, carries no tag field, so no other
        // consumer can be reading one. Resolving the read group **here**, while the tags are
        // still in hand, is therefore what lets every tag be dropped — on both arms, not just on
        // the single-read-group one.
        //
        // **One behaviour change, and it is only visible on a malformed file.** A record whose
        // `RG` is unreadable used to fail when that record was *served*; it now fails when its
        // container is decoded, which can be earlier and can involve records outside the queried
        // region. The condition was already fatal either way (`ReadGroupResolution::PerRecord`) —
        // this reaches it sooner. `a_cram_with_several_read_groups_resolves_each_record_by_its_tag`
        // (`input/mod.rs`) covers the well-formed path through here.
        let compression_header = container.compression_header()?;
        let mut decoded = DecodedContainer {
            offset,
            index: Vec::new(),
            payload: Vec::new(),
            cigar_ops: Vec::new(),
        };
        for slice in container.slices() {
            let slice = slice?;
            // The decoded block data and the borrowed records live only within this block;
            // copying each record's bytes into the container's own buffers here keeps the result
            // independent of those borrows.
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            for record in slice.records(
                self.repository.clone(),
                &self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )? {
                let record = RecordBuf::try_from_alignment_record(&self.header, &record)?;
                let owner = resolve_read_group(&record, &self.resolution)?;
                decoded.push(&record, owner)?;
            }
        }

        // The three buffers grew by doubling, so each can be holding up to twice what it needs —
        // and this container is then resident for the life of the file. Measured, the slack was
        // 1.6 MiB of 5.4 MiB per open file, so paying one copy per container to drop it is
        // strongly worth it: a container is decoded once and read from thousands of times.
        decoded.shrink_to_fit();

        Ok(Some(decoded))
    }
}

impl RecordSource for CramRegionSource {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        &self.header
    }

    fn other_sample_records(&self) -> u64 {
        self.other_sample_records
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        buf.read_group = None;
        loop {
            if let Some(index) = self.pending.next() {
                let held = self
                    .container
                    .as_ref()
                    .expect("pending indices only exist while their container is held");
                // The read group was resolved at decode, while the record still carried its `RG`
                // tag — which is what lets the tags be dropped (see `decode_container_at`).
                // Copied out before the record is rebuilt so the borrow ends here.
                let owner = held.index[index as usize].owner;
                match owner {
                    RecordOwner::Mine(id) => {
                        // Rebuilt into the caller's buffer rather than moved in, so the buffer's
                        // allocations are reused across the whole stream (`DecodedContainer::fill`).
                        held.fill(index as usize, &mut buf.record);
                        buf.read_group = Some(id);
                        return Ok(true);
                    }
                    RecordOwner::OtherSample => {
                        self.other_sample_records += 1;
                        continue;
                    }
                }
            }

            match self.refill() {
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
}

/// Whether a record's reference footprint touches the region. Shared by both
/// sources so the two containers cannot disagree about what "overlap" means —
/// which is precisely what the BAM/CRAM parity oracle (T8) would otherwise be
/// left to discover.
fn overlaps(record: &sam::alignment::RecordBuf, region: GenomeRegion) -> bool {
    match (record.alignment_start(), record.alignment_end()) {
        (Some(first), Some(last)) => {
            usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get()
        }
        _ => false,
    }
}

/// Proves, while streaming, that the reads really do arrive in genome order.
///
/// Everything downstream is built on that, so it is **checked, not trusted**:
/// the last emitted [`GenomePosition`] is kept, and a read whose key is
/// **strictly less** than it is a hard [`AlignmentFileError::OutOfOrderRead`].
/// There is no tolerance, no warning-and-continue, and no silent re-sort — an
/// unsorted file is a fatal input error the user must see.
///
/// Four things fix the semantics (spec §3.2):
///
/// - **The key is `(contig, position)`**, so a position regression within a
///   contig and a contig-order regression are one violation, not two checks.
///   "Contig order" means the reference's, because the open gate proved
///   `ref_id == ContigId`.
/// - **Equal keys are legal.** Several reads may start at the same base; only a
///   decrease is rejected.
/// - **It sits on the *filtered* stream**, so the guarantee is about exactly
///   the reads the caller sees. Dropped reads cannot break the monotonicity of
///   what survives.
/// - **The state lives here, in the per-region iterator, never on the handle.**
///   A caller is entitled to query region B and then region A; that is a new
///   forward scan, not a regression. Carrying the last key on the handle would
///   turn legitimate random access into a spurious error.
///
/// This and the `@HD SO` check at open are complementary rather than redundant:
/// the header check is cheap and fails before any work, while this one catches
/// the file that *claims* to be sorted and is not — which the header check
/// structurally cannot see.
pub(crate) struct OrderVerified<I> {
    inner: I,
    /// The last key emitted, or `None` before the first read.
    last: Option<GenomePosition>,
    path: Arc<Path>,
    /// Set on clean end of input or after an error is yielded, so the iterator
    /// is fused: the first `Err` surfaces once and then `None`.
    done: bool,
}

impl<I> OrderVerified<I> {
    pub(crate) fn new(inner: I, path: Arc<Path>) -> Self {
        Self {
            inner,
            last: None,
            path,
            done: false,
        }
    }

    /// Unwrap, so the stream that owns this can reclaim what it lent.
    pub(crate) fn into_inner(self) -> I {
        self.inner
    }
}

/// The region reader for whichever container the file is.
///
/// A [`RecordSource`] that delegates to the per-format source. BAM and CRAM are
/// two containers for one idea, so this is an enum rather than a trait — and it
/// is what lets everything above the source be written once, generic over
/// neither format nor container (`arch/alignment_file.md` §4).
pub(crate) enum RegionSource {
    Bam(BamRegionSource),
    Cram(CramRegionSource),
}

impl RegionSource {
    /// Give the reader back, in the shape the pool stores it — with the decoded container a CRAM
    /// query is holding, so the next query through this reader can reuse it.
    pub(crate) fn into_parts(self) -> (super::open_bam::ReaderKind, Option<DecodedContainer>) {
        match self {
            Self::Bam(source) => (super::open_bam::ReaderKind::Bam(source.into_reader()), None),
            Self::Cram(source) => {
                let (reader, container) = source.into_parts();
                (super::open_bam::ReaderKind::Cram(reader), container)
            }
        }
    }
}

impl RecordSource for RegionSource {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        match self {
            Self::Bam(source) => source.header(),
            Self::Cram(source) => source.header(),
        }
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        match self {
            Self::Bam(source) => source.read_next(buf),
            Self::Cram(source) => source.read_next(buf),
        }
    }

    /// Forwarded, not defaulted. The trait's default is `0`, which for a
    /// delegating wrapper is a silent wrong answer rather than an absent one:
    /// the reads really were skipped, and the tally would simply lose them.
    fn other_sample_records(&self) -> u64 {
        match self {
            Self::Bam(source) => source.other_sample_records(),
            Self::Cram(source) => source.other_sample_records(),
        }
    }
}

/// Where a read sits in the genome. Sound as a cross-file comparison key only
/// because the open gate proved this file's `ref_id`s are the reference's
/// `ContigId`s.
fn key_of(read: &AlignedRead) -> GenomePosition {
    GenomePosition {
        // PANIC-FREE: `ref_id` comes from a 32-bit field in both BAM and CRAM,
        // so it fits by construction. Checked rather than `as`-cast because a
        // wrapped value would collapse two contigs onto one key and silently
        // disarm this guard — the one thing it must never do. Matches
        // `filtering.rs`'s handling of the same value.
        contig: ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32")),
        position: Position(read.pos),
    }
}

impl<I> Iterator for OrderVerified<I>
where
    I: Iterator<Item = Result<AlignedRead, ReadFilterError>>,
{
    type Item = Result<AlignedRead, AlignmentFileError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match self.inner.next() {
            None => {
                self.done = true;
                None
            }
            Some(Err(error)) => {
                self.done = true;
                Some(Err(AlignmentFileError::Filter(error)))
            }
            Some(Ok(read)) => {
                let current = key_of(&read);
                if let Some(previous) = self.last
                    && current < previous
                {
                    self.done = true;
                    return Some(Err(AlignmentFileError::OutOfOrderRead {
                        path: self.path.to_path_buf(),
                        previous,
                        current,
                    }));
                }
                self.last = Some(current);
                Some(Ok(read))
            }
        }
    }
}

impl<I> FusedIterator for OrderVerified<I> where
    I: Iterator<Item = Result<AlignedRead, ReadFilterError>>
{
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use tempfile::TempDir;

    use super::*;
    use crate::bam::index_preflight::load_alignment_index;
    use crate::ng::read::filtering::BamRecordSource;
    use crate::ng::read::input::open_bam::group_crai_by_contig;
    use crate::ng::read::input::test_fixtures::{bam_header, indexed_bam, read_named};
    use crate::ng::types::ReadGroupId;
    use crate::ng::types::{ContigId, Position};

    /// Two contigs long enough that the reads below span **several BGZF
    /// blocks**, so the index returns more than one chunk and the seek-per-chunk
    /// path is genuinely exercised. A fixture small enough to fit one block
    /// would make the oracle agree for the wrong reason.
    const CONTIGS: [(&str, usize); 2] = [("chr1", 100_000), ("chr2", 50_000)];

    fn contig_specs() -> Vec<(&'static str, usize, Option<&'static str>)> {
        CONTIGS.iter().map(|(n, l)| (*n, *l, None)).collect()
    }

    /// A coordinate-sorted spread of reads over both contigs: one every 20 bp,
    /// plus deliberate pile-ups so several reads share a start.
    fn spread_of_reads() -> Vec<RecordBuf> {
        let mut records = Vec::new();
        for (reference_sequence_id, (_, length)) in CONTIGS.iter().enumerate() {
            let mut start = 1;
            while start + 10 < *length {
                records.push(read_named(
                    &format!("r{reference_sequence_id}_{start}"),
                    reference_sequence_id,
                    start,
                ));
                // Every 50th read is doubled at the same start, so equal
                // positions are part of the fixture rather than an edge case
                // the oracle never sees.
                if start % 1000 == 1 {
                    records.push(read_named(
                        &format!("r{reference_sequence_id}_{start}_b"),
                        reference_sequence_id,
                        start,
                    ));
                }
                start += 20;
            }
        }

        // The shapes `overlaps_region`'s `_ => false` arm exists for, which a
        // spread of ordinary reads would never produce: a placed record with no
        // alignment start, and an entirely unplaced one. Both must be dropped,
        // and neither is a *filter* drop.
        records.push(unmapped_but_placed("unmapped_placed", 0));
        records.push(unplaced("unplaced"));
        records
    }

    /// A record naming a contig but carrying no alignment start — legal in SAM,
    /// and it has no footprint to overlap anything.
    fn unmapped_but_placed(qname: &str, reference_sequence_id: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(reference_sequence_id)
            .set_flags(noodles_sam::alignment::record::Flags::UNMAPPED)
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    /// A record on no contig at all — the unplaced reads that sit at a sorted
    /// file's tail.
    fn unplaced(qname: &str) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_flags(noodles_sam::alignment::record::Flags::UNMAPPED)
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    /// The header comes back behind an `Arc` because that is what a region
    /// source now takes — the same shape `AlignmentFile` hands out. It still
    /// derefs to `&sam::Header` wherever a plan or the oracle wants one.
    fn fixture() -> (TempDir, std::path::PathBuf, Arc<sam::Header>) {
        let header = bam_header(&contig_specs());
        let (dir, path) = indexed_bam(&header, &spread_of_reads());
        (dir, path, Arc::new(header))
    }

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **The oracle: a whole-file linear scan, filtered to the region.**
    ///
    /// Independent of the code under test in the dimension that matters here:
    /// it reaches every record through the whole-file `BamRecordSource` and
    /// never seeks, so nothing about chunks, seeking or the early stop is
    /// shared. The overlap rule is retyped rather than called, which catches a
    /// typo but not a shared misunderstanding of what "overlap" means — that
    /// part is pinned by the explicit boundary cases below instead.
    fn reads_by_linear_scan(
        path: &Path,
        header: &sam::Header,
        region: GenomeRegion,
    ) -> Vec<String> {
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        reader.read_header().expect("read header");
        let mut source = BamRecordSource::new(
            reader,
            header.clone(),
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );

        let target = region.contig.get() as usize;
        let mut buf = NoodlesRawRecord::default();
        let mut names = Vec::new();

        while source.read_next(&mut buf).expect("read") {
            let record = &buf.record;
            if record.reference_sequence_id() != Some(target) {
                continue;
            }
            let (Some(first), Some(last)) = (record.alignment_start(), record.alignment_end())
            else {
                continue;
            };
            let overlaps = usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get();
            if overlaps {
                names.push(read_name(record));
            }
        }
        names
    }

    /// The same question through the **indexed** path.
    fn reads_by_index(path: &Path, header: &Arc<sam::Header>, region: GenomeRegion) -> Vec<String> {
        let index = load_alignment_index(path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let plan = BamRegionSource::plan(header, &index, region, path).expect("plan the query");
        let mut source = BamRegionSource::new(
            reader,
            Arc::clone(header),
            plan,
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );

        let mut buf = NoodlesRawRecord::default();
        let mut names = Vec::new();
        while source.read_next(&mut buf).expect("read") {
            names.push(read_name(&buf.record));
        }
        names
    }

    fn read_name(record: &RecordBuf) -> String {
        String::from_utf8_lossy(record.name().expect("named").as_ref()).into_owned()
    }

    /// **T5 — the indexed query equals a full linear scan of the same region.**
    ///
    /// This is the step whose failure is silent: a chunk edge handled wrongly
    /// drops reads that were really there, and the caller gets a plausible
    /// answer with a wrong genotype rather than a crash. Nothing the indexed
    /// path could assert about *itself* would catch that, so it is checked
    /// against an independent implementation over a range of regions —
    /// interior, boundary, empty, whole-contig, and the second contig.
    #[test]
    fn t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns() {
        let (_dir, path, header) = fixture();

        let regions = [
            region(0, 1, 100),          // the very start
            region(0, 1, 1),            // a single base
            region(0, 40_000, 40_100),  // deep in the interior
            region(0, 21, 21),          // exactly one read's start
            region(0, 30, 35),          // covered only by an overlapping read
            region(0, 1, 100_000),      // the whole contig
            region(0, 99_990, 100_000), // the far end
            region(1, 1, 500),          // the second contig
            region(1, 25_000, 25_050),  // its interior
            region(1, 49_000, 50_000),  // its far end
            region(0, 50_001, 50_001),  // a position with no read starting
        ];

        for region in regions {
            let expected = reads_by_linear_scan(&path, &header, region);
            let actual = reads_by_index(&path, &header, region);

            assert_eq!(
                actual,
                expected,
                "indexed query disagreed with the linear scan for contig {} [{}, {}]",
                region.contig.get(),
                region.start.get(),
                region.end.get()
            );
        }
    }

    /// The oracle is only worth something if it can fail. A region the fixture
    /// genuinely covers must return reads — otherwise T5 could be comparing two
    /// empty vectors and calling it agreement.
    #[test]
    fn the_oracle_is_not_vacuous() {
        let (_dir, path, header) = fixture();

        let interior = reads_by_linear_scan(&path, &header, region(0, 40_000, 40_100));
        assert!(
            interior.len() >= 5,
            "the fixture should cover this region several times over, got {}",
            interior.len()
        );

        let empty = reads_by_linear_scan(&path, &header, region(0, 99_999, 100_000));
        assert!(
            empty.len() < interior.len(),
            "and the regions must not all return the same thing"
        );
    }

    /// Records with no footprint — placed-but-unmapped, and unplaced — are
    /// dropped by the region query, and dropped **uncounted**: they are not
    /// reads the filter rejected, they are reads the index was never asked
    /// about. T5 alone cannot show this, because the oracle drops them by the
    /// same rule and the two would agree however the arm behaved.
    #[test]
    fn records_without_a_footprint_never_surface_from_a_region_query() {
        let (_dir, path, header) = fixture();

        for region in [region(0, 1, 100_000), region(1, 1, 50_000)] {
            let names = reads_by_index(&path, &header, region);
            assert!(
                !names
                    .iter()
                    .any(|n| n.starts_with("unmapped_placed") || n.starts_with("unplaced")),
                "a footprint-less record surfaced from {names:?}"
            );
        }
    }

    /// Reads are returned in coordinate order, which the whole chain above
    /// assumes and the merge one layer up depends on.
    #[test]
    fn the_indexed_query_returns_reads_in_coordinate_order() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            Arc::clone(&header),
            BamRegionSource::plan(&header, &index, region(0, 1, 100_000), &path).expect("plan"),
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );

        let mut buf = NoodlesRawRecord::default();
        let mut previous = 0u64;
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            let start = usize::from(buf.record.alignment_start().expect("mapped")) as u64;
            assert!(start >= previous, "{start} came after {previous}");
            previous = start;
            seen += 1;
        }
        assert!(seen > 100, "the scan should have covered the contig");
    }

    /// **The early stop.** Once a read on the target contig starts past the
    /// region, the sort guarantees nothing later overlaps — so the scan must
    /// end rather than read on to the end of the file.
    ///
    /// The assertion that matters is the *unconsumed chunks*, not `done`:
    /// `done` latches both on the early stop and on running out of chunks, so
    /// a test that only checked it would pass with the early stop deleted. This
    /// failure changes no answers, only work — a whole-file scan per region —
    /// which is exactly why it needs its own discriminating assertion.
    #[test]
    fn the_scan_stops_once_it_passes_the_region() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            Arc::clone(&header),
            BamRegionSource::plan(&header, &index, region(0, 1, 100), &path).expect("plan"),
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );

        // Drain, then confirm the source latched `done` rather than running on.
        let mut buf = NoodlesRawRecord::default();
        while source.read_next(&mut buf).expect("read") {}
        assert!(
            source.done,
            "the early stop must latch, so a drained query does not resume"
        );
        assert!(
            source.chunks.len() > 0,
            "the scan consumed every chunk for a 100 bp region of a 100 kb \
             contig — it read on instead of stopping early"
        );
    }

    /// The source owns its reader so the pool can have it back; a stream that
    /// ended must be able to hand over a reader that still works, or every
    /// query after the first would open a fresh one.
    #[test]
    fn the_reader_survives_the_query_and_comes_back_usable() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            Arc::clone(&header),
            BamRegionSource::plan(&header, &index, region(0, 1, 100), &path).expect("plan"),
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );
        let mut buf = NoodlesRawRecord::default();
        while source.read_next(&mut buf).expect("read") {}

        // Hand it back, then drive a second query with the very same reader.
        let reader = source.into_reader();
        let mut second = BamRegionSource::new(
            reader,
            Arc::clone(&header),
            BamRegionSource::plan(&header, &index, region(1, 1, 200), &path).expect("plan"),
            ReadGroupResolution::Sole(ReadGroupId(0)),
        );
        let mut seen = 0;
        while second.read_next(&mut buf).expect("read") {
            seen += 1;
        }
        assert!(
            seen > 0,
            "the returned reader must still be able to seek and read"
        );
    }

    #[test]
    fn a_region_naming_a_contig_the_file_does_not_have_is_an_error() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        assert!(matches!(
            BamRegionSource::plan(&header, &index, region(9, 1, 100), &path),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    // -----------------------------------------------------------------
    // The order guard (C3) — T4a..T4d
    // -----------------------------------------------------------------

    /// A `AlignedRead` at a genome position. Only `ref_id` and `pos` matter to
    /// the guard, so the rest is minimal — the guard is a pure adapter over the
    /// filtered stream and never looks at sequence or CIGAR.
    fn read_at(qname: &str, ref_id: usize, pos: u64) -> AlignedRead {
        AlignedRead {
            qname: qname.as_bytes().to_vec(),
            flag: 0,
            ref_id,
            pos,
            mapq: 60,
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    fn planted_path() -> Arc<Path> {
        Arc::from(Path::new("/data/planted.bam"))
    }

    /// An iterator that yields `None` and then **keeps going** — what a
    /// non-fused inner looks like. `std::iter::from_fn` cannot model this,
    /// because it fuses itself.
    struct Resuming {
        items: Vec<Option<Result<AlignedRead, ReadFilterError>>>,
        next_index: usize,
    }

    impl Iterator for Resuming {
        type Item = Result<AlignedRead, ReadFilterError>;

        fn next(&mut self) -> Option<Self::Item> {
            let item = self.items.get_mut(self.next_index)?.take();
            self.next_index += 1;
            item
        }
    }

    /// Drive the guard over a planted stream and collect what a caller would
    /// see: the read names that surfaced, and the error if one did.
    fn through_the_guard(
        reads: Vec<AlignedRead>,
    ) -> (Vec<String>, Option<AlignmentFileError>, usize) {
        let stream = reads.into_iter().map(Ok);
        let mut guard = OrderVerified::new(stream, planted_path());

        let mut names = Vec::new();
        let mut error = None;
        for item in &mut guard {
            match item {
                Ok(read) => names.push(String::from_utf8_lossy(&read.qname).into_owned()),
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }
        // How many items the iterator yields *after* the error — must be zero,
        // because the guard fuses.
        let after = guard.count();
        (names, error, after)
    }

    /// **T4a — a planted position regression is fatal.**
    ///
    /// The header can say `SO:coordinate` and the data can still be unsorted;
    /// that is the case the open gate structurally cannot see. The error names
    /// both keys, so the message can say *where* the file breaks rather than
    /// only that it does.
    ///
    /// Mutation-verified: deleting the comparison in `OrderVerified::next`
    /// makes this test fail, which is the whole point — a guard that never
    /// fires looks exactly like a guard that works.
    #[test]
    fn t4a_a_position_regression_within_a_contig_is_an_error() {
        let (names, error, after) = through_the_guard(vec![
            read_at("a", 0, 100),
            read_at("b", 0, 200),
            read_at("c", 0, 150), // backwards
            read_at("d", 0, 300),
        ]);

        assert_eq!(names, vec!["a", "b"], "the reads before the break flow");
        match error.expect("the regression must be fatal") {
            AlignmentFileError::OutOfOrderRead {
                path,
                previous,
                current,
            } => {
                assert_eq!(previous.position.get(), 200);
                assert_eq!(current.position.get(), 150);
                // Asserted because `path` is a constructor argument: a wiring
                // mistake at C4 would produce a correct-looking error naming
                // the wrong file.
                assert_eq!(path, PathBuf::from("/data/planted.bam"));
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
        assert_eq!(after, 0, "and the iterator fuses — no read after the error");
    }

    /// **T4b — a contig-order regression is the same violation.** The key is
    /// `(contig, position)`, so a read on an earlier contig after a later one
    /// is caught by the same comparison, even though its *position* went up.
    #[test]
    fn t4b_a_contig_order_regression_is_the_same_error() {
        let (names, error, _) = through_the_guard(vec![
            read_at("a", 1, 100),
            read_at("b", 0, 900), // earlier contig, later position
        ]);

        assert_eq!(names, vec!["a"]);
        match error.expect("a contig regression must be fatal") {
            AlignmentFileError::OutOfOrderRead {
                previous, current, ..
            } => {
                assert_eq!(previous.contig.get(), 1);
                assert_eq!(current.contig.get(), 0);
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
    }

    /// **T4c — equal keys are legal.** Several reads may start at the same
    /// base; only a strict decrease is a violation. A guard using `<=` would
    /// reject every pile-up in every real file.
    #[test]
    fn t4c_reads_sharing_a_start_position_are_not_a_regression() {
        let (names, error, _) = through_the_guard(vec![
            read_at("a", 0, 100),
            read_at("b", 0, 100),
            read_at("c", 0, 100),
            read_at("d", 0, 101),
        ]);

        assert_eq!(names, vec!["a", "b", "c", "d"]);
        assert!(error.is_none(), "equal positions are ordinary, not a fault");
    }

    /// **T4d — the check does not span queries.** A caller may query region B
    /// and then region A; the second is a new forward scan, not a regression.
    /// This is why the state lives in the iterator and never on the handle —
    /// carrying it there would turn legitimate random access into a spurious
    /// error.
    #[test]
    fn t4d_querying_a_later_region_then_an_earlier_one_is_not_a_regression() {
        let (later, error, _) = through_the_guard(vec![read_at("b1", 0, 9000)]);
        assert_eq!(later, vec!["b1"]);
        assert!(error.is_none());

        // A *separate* guard, as a second query gets — starting at position 10,
        // far behind where the first one ended.
        let (earlier, error, _) = through_the_guard(vec![read_at("a1", 0, 10)]);
        assert_eq!(earlier, vec!["a1"]);
        assert!(
            error.is_none(),
            "a new region query starts a new scan; it cannot regress against \
             a previous one"
        );

        // And the fresh guard is *armed*, not merely permissive — a guard that
        // had been silently disabled would also pass the assertion above.
        let (_, error, _) = through_the_guard(vec![read_at("a1", 0, 10), read_at("a2", 0, 5)]);
        assert!(
            matches!(error, Some(AlignmentFileError::OutOfOrderRead { .. })),
            "the second query's own guard still catches its own regression"
        );
    }

    /// **The other half of the fuse.** `FusedIterator` is claimed without
    /// requiring `I: FusedIterator`, so the latch has to hold against an inner
    /// that returns `None` and then resumes — otherwise a read could reappear
    /// after the stream was reported ended. Only the *error* path's fuse was
    /// covered before; deleting `done = true` from the clean-exhaustion arm
    /// left every test green.
    #[test]
    fn the_guard_stays_ended_even_if_what_it_wraps_resumes() {
        let inner = Resuming {
            items: vec![
                Some(Ok(read_at("a", 0, 1))),
                None,
                Some(Ok(read_at("b", 0, 2))),
            ],
            next_index: 0,
        };
        let mut guard = OrderVerified::new(inner, planted_path());

        assert!(matches!(guard.next(), Some(Ok(_))));
        assert!(guard.next().is_none(), "the inner reported end of input");
        assert!(
            guard.next().is_none(),
            "and the guard stays ended — the latch holds, not just the inner"
        );
    }

    /// An error from the filter below is wrapped and fuses the stream, rather
    /// than being mistaken for a clean end of input.
    #[test]
    fn a_filter_error_is_wrapped_and_fuses_the_guard() {
        let stream = vec![
            Ok(read_at("a", 0, 1)),
            Err(ReadFilterError::Source(io::Error::other("planted"))),
            Ok(read_at("b", 0, 2)),
        ]
        .into_iter();
        let mut guard = OrderVerified::new(stream, planted_path());

        assert!(matches!(guard.next(), Some(Ok(_))));
        assert!(matches!(
            guard.next(),
            Some(Err(AlignmentFileError::Filter(_)))
        ));
        assert!(
            guard.next().is_none(),
            "fused: the read after the error is not reachable"
        );
    }

    // -----------------------------------------------------------------
    // The CRAM `.crai` walk (C5)
    // -----------------------------------------------------------------

    /// A `.crai` entry for `contig` covering `[start, start + span - 1]`.
    fn crai_index(records: Vec<cram::crai::Record>) -> cram::crai::Index {
        records
    }

    /// A `.crai` entry for `contig` covering `[start, start + span - 1]`, at
    /// `offset` — which the tests use to tell entries apart.
    fn crai_entry(
        contig: Option<usize>,
        start: usize,
        span: usize,
        offset: u64,
    ) -> cram::crai::Record {
        cram::crai::Record::new(
            contig,
            noodles_core::Position::new(start),
            span,
            offset,
            0,
            0,
        )
    }

    fn three_contig_header() -> sam::Header {
        bam_header(&[
            ("chr1", 100, None),
            ("chr2", 200, None),
            ("chr3", 300, None),
        ])
    }

    /// **The fix for production's O(n) prefix rescan** — and, more importantly,
    /// for an ordering assumption that does not hold. A query on a late contig
    /// must reach that contig's entries without walking (or mis-searching) the
    /// earlier ones.
    ///
    /// Hand-built rather than read from a file, because the file fixture cannot
    /// hold more than one contig (`test_fixtures::indexed_cram`).
    #[test]
    fn the_crai_is_grouped_so_each_contig_sees_only_its_own_entries() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(0), 51, 50, 200),
            crai_entry(Some(1), 1, 100, 300),
            crai_entry(Some(1), 101, 100, 400),
            crai_entry(Some(2), 1, 300, 500),
        ]);
        let grouped = group_crai_by_contig(&index, 3);
        let header = three_contig_header();

        for (contig, expected_offsets) in [
            (0u32, vec![100u64, 200]),
            (1, vec![300, 400]),
            (2, vec![500]),
        ] {
            let plan =
                CramRegionSource::plan(&header, &grouped, region(contig, 1, 10)).expect("plan");
            let offsets: Vec<u64> = plan.entries.iter().map(|e| e.offset()).collect();
            assert_eq!(offsets, expected_offsets, "contig {contig}");
        }
    }

    /// **The ordering assumption a binary search would have made, violated by
    /// noodles itself.** Within a multi-reference slice, `fs::index` sorts by
    /// `Option<usize>` — and `None < Some(0)` — so the *unplaced* entry is
    /// emitted first. A `partition_point` over that is unspecified, and the
    /// walk would then meet a foreign entry and report end-of-input, losing
    /// every read of the region with no error at all.
    ///
    /// Grouping assumes nothing about order, so an interleaved index is fine.
    #[test]
    fn an_unplaced_entry_before_the_placed_ones_does_not_hide_a_contig() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(None, 1, 0, 150),
            crai_entry(Some(1), 1, 100, 300),
        ]);
        let grouped = group_crai_by_contig(&index, 3);

        let plan = CramRegionSource::plan(&three_contig_header(), &grouped, region(1, 1, 10))
            .expect("plan");
        let offsets: Vec<u64> = plan.entries.iter().map(|e| e.offset()).collect();
        assert_eq!(
            offsets,
            vec![300],
            "contig 1's entry must be found even though an unplaced entry \
             precedes it — the case that silently returned nothing before"
        );

        // And the unplaced entry belongs to no contig at all.
        for (contig, entries) in grouped.iter().enumerate() {
            assert!(
                entries.iter().all(|e| e.reference_sequence_id().is_some()),
                "an unplaced entry leaked into contig {contig}"
            );
        }
    }

    /// A contig with no entries yields an empty walk rather than someone
    /// else's entries.
    #[test]
    fn a_contig_absent_from_the_crai_has_an_empty_walk() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(2), 1, 300, 500),
        ]);
        let grouped = group_crai_by_contig(&index, 3);

        let plan = CramRegionSource::plan(&three_contig_header(), &grouped, region(1, 1, 10))
            .expect("plan");
        assert!(plan.entries.is_empty());
    }

    /// **The CRAM container-level early stop, with a discriminating
    /// assertion.** Draining a region near the start of a long contig must
    /// leave most of the `.crai` unwalked; a test that only checked the reads
    /// were right would pass with the stop deleted, because deleting it changes
    /// the *work*, not the answer — a whole-file decode per region, ~10⁶ times.
    ///
    /// Driven through `CramRegionSource` directly, since the cursor is not
    /// visible once the source is wrapped in the filter and the order guard.
    #[test]
    fn the_crai_walk_stops_once_it_passes_the_region() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        const CONTIG_LENGTH: usize = 400_000;
        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, 30_000);

        let index = load_alignment_index(&cram_path).expect("load crai");
        let repository =
            crate::bam::alignment_input::build_fasta_repository(&fasta).expect("repository");
        let header_for_plan = bam_header(&[("chr1", CONTIG_LENGTH, None)]);
        let grouped = group_crai_by_contig(
            match &index {
                AlignmentIndex::Crai(crai) => crai,
                _ => panic!("a .cram must carry a .crai"),
            },
            1,
        );
        assert!(
            grouped[0].len() > 1,
            "the fixture must span several containers, or there is nothing to \
             stop short of — got {}",
            grouped[0].len()
        );

        let mut reader = cram::io::Reader::new(File::open(&cram_path).expect("open"));
        let parsed_header = reader.read_header().expect("read header");

        let plan =
            CramRegionSource::plan(&header_for_plan, &grouped, region(0, 1, 200)).expect("plan");
        let total_entries = plan.entries.len();
        let mut source = CramRegionSource::new(
            reader,
            Arc::new(parsed_header),
            repository,
            plan,
            ReadGroupResolution::Sole(ReadGroupId(0)),
            None,
        );

        let mut buf = NoodlesRawRecord::default();
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            seen += 1;
        }

        assert!(seen > 0, "the region is covered");
        assert!(
            source.next_index_record < total_entries,
            "the walk consumed all {total_entries} entries for a 200 bp region \
             of a 400 kb contig — it decoded to the end instead of stopping"
        );
    }

    /// **The decoded-container cache must never serve a stale container's reads**
    /// ([`DecodedContainer`]). A reader that has just decoded the container holding one region is
    /// then asked for a region in a *different* container, and then for the first region again;
    /// every answer must equal what a reader with a cold cache returns for that region alone.
    ///
    /// This is the load-bearing test for the cache, because the failure is silent: reusing on a
    /// stale key does not crash, it reports another part of the contig's reads as this region's
    /// depth. A cache that ignored the `.crai` offset passes every "the reads are right" test on
    /// the *first* query and fails here on the second.
    #[test]
    fn a_stale_container_is_not_reused() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::open_bam::AlignmentFile;
        use crate::ng::reference_info::{ReferenceSource, read_reference_info};

        const CONTIG_LENGTH: usize = 400_000;
        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, 30_000);
        // `Fasta`, not `Fai`: a CRAM needs the reference *bases* to decode against, so `open`
        // rejects a reference that cannot supply a path to them.
        let reference = crate::ng::read::input::reference::OpenReference::from(
            read_reference_info(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .expect("the fixture reference reads"),
        );

        // Two regions far enough apart to be in different containers (the fixture writes ~30k
        // records per container over 400 kb).
        let early = region(0, 1_000, 1_200);
        let late = region(
            0,
            (CONTIG_LENGTH as u64 * 3) / 4,
            (CONTIG_LENGTH as u64 * 3) / 4 + 200,
        );

        let open = || {
            AlignmentFile::open(
                &cram_path,
                &reference,
                ReadFilterConfig::default(),
                false,
                crate::ng::read::input::test_fixtures::fixture_read_group(),
            )
            .expect("the fixture CRAM opens")
        };

        // What each region looks like to a reader that has never decoded anything.
        let cold_early = positions_of(&open(), early, reference.info());
        let cold_late = positions_of(&open(), late, reference.info());
        assert!(!cold_early.is_empty(), "the early region is covered");
        assert_ne!(
            cold_early, cold_late,
            "the two regions must return different reads, or the test cannot \
             distinguish a stale cache from a correct one"
        );

        // One file, so all three queries share the pool — and therefore the cache.
        let file = open();
        let warm_early = positions_of(&file, early, reference.info());
        let warm_late = positions_of(&file, late, reference.info());
        let warm_early_again = positions_of(&file, early, reference.info());

        assert_eq!(warm_early, cold_early, "first query, cache cold");
        assert_eq!(
            warm_late, cold_late,
            "the second query must decode its own container, not reuse the first's"
        );
        assert_eq!(
            warm_early_again, cold_early,
            "returning to the first region must not serve the late container's reads"
        );
        assert_eq!(
            file.readers_opened(),
            1,
            "all three queries went through one pooled reader, so they shared one cache"
        );
    }

    /// Every read position a region query yields, through the whole real chain (pool → source →
    /// filter → order guard) — the shape a caller sees, so the cache is exercised where it lives.
    fn positions_of(
        file: &Arc<crate::ng::read::input::open_bam::AlignmentFile>,
        span: GenomeRegion,
        reference: &crate::ng::reference_info::ReferenceInfo,
    ) -> Vec<u64> {
        let bases = crate::ng::ref_seq::WindowedRefSeq::new(
            reference
                .fasta_path
                .clone()
                .expect("the fixture reference is a FASTA"),
            reference.contig_list(),
        );
        file.reads_in_region(span, bases)
            .expect("the query plans")
            .map(|read| read.expect("no fatal read error").pos)
            .collect()
    }

    /// The span-based skip: a container that ends *before* the region must be
    /// stepped over without being decoded. Observable as the cursor advancing
    /// past entries that produced no reads.
    #[test]
    fn the_crai_walk_skips_containers_that_end_before_the_region() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        const CONTIG_LENGTH: usize = 400_000;
        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, 30_000);

        let index = load_alignment_index(&cram_path).expect("load crai");
        let repository =
            crate::bam::alignment_input::build_fasta_repository(&fasta).expect("repository");
        let header_for_plan = bam_header(&[("chr1", CONTIG_LENGTH, None)]);
        let grouped = group_crai_by_contig(
            match &index {
                AlignmentIndex::Crai(crai) => crai,
                _ => panic!("a .cram must carry a .crai"),
            },
            1,
        );

        let mut reader = cram::io::Reader::new(File::open(&cram_path).expect("open"));
        let parsed_header = reader.read_header().expect("read header");

        // A region in the *last* stretch of the contig: every earlier container
        // ends before it and must be skipped rather than decoded.
        let late = (CONTIG_LENGTH as u64 * 3) / 4;
        let plan = CramRegionSource::plan(&header_for_plan, &grouped, region(0, late, late + 200))
            .expect("plan");
        let mut source = CramRegionSource::new(
            reader,
            Arc::new(parsed_header),
            repository,
            plan,
            ReadGroupResolution::Sole(ReadGroupId(0)),
            None,
        );

        let mut buf = NoodlesRawRecord::default();
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            seen += 1;
        }

        assert!(seen > 0, "the late region is covered");
        assert!(
            source.last_decoded_offset.is_some(),
            "at least one container was decoded"
        );
    }

    #[test]
    fn a_cram_region_naming_a_contig_the_header_lacks_is_an_error() {
        let grouped = group_crai_by_contig(&crai_index(vec![crai_entry(Some(0), 1, 50, 0)]), 3);
        assert!(matches!(
            CramRegionSource::plan(&three_contig_header(), &grouped, region(9, 1, 10)),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    /// The CRAM planner rejects an inverted or empty region, as the BAM one
    /// does — the two must not disagree about what a valid region is.
    #[test]
    fn a_cram_planner_rejects_an_inverted_region() {
        let grouped = group_crai_by_contig(&crai_index(vec![crai_entry(Some(0), 1, 50, 0)]), 3);
        assert!(matches!(
            CramRegionSource::plan(&three_contig_header(), &grouped, region(0, 500, 100)),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    #[test]
    fn an_inverted_region_is_an_error() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        assert!(matches!(
            BamRegionSource::plan(&header, &index, region(0, 500, 100), &path),
            Err(AlignmentFileError::Region { .. })
        ));
    }
}
