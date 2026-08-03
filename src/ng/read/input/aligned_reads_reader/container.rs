//! **One CRAM container, decoded** — and the one body that turns a container into records.
//!
//! CRAM decodes at container granularity: a container holds thousands of records covering a
//! wide span of the reference, and it must be decompressed and decoded whole before any single
//! record can be looked at. So this is the unit the CRAM arm reads in, and
//! [`decode_container_at`] is how one becomes records.
//!
//! It lived in `region_query.rs` until Milestone F deleted that module. Nothing about it was
//! the per-region query's: a container is a container whoever is walking them.
//!
//! ## Why it is stored packed rather than as records
//!
//! A container is ~10⁴ records, and a `RecordBuf` per record was ~680 bytes across seven
//! separate allocations — against ~270 bytes of read that anything actually consumes (name,
//! sequence, qualities, CIGAR). The other 60% was `Vec` headers, capacity slack and
//! per-allocation overhead, ~72,000 allocations per open file, resident for as long as the
//! container is being served.
//!
//! So the bytes go into two flat buffers, [`payload`](DecodedContainer::payload) and
//! [`cigar_ops`](DecodedContainer::cigar_ops), and each record becomes a fixed-size
//! [`PackedReadEntry`] naming its slices of them. A record is materialised as a `RecordBuf`
//! **only when it is actually served**, straight into the caller's reused buffer — which is
//! also how the per-record `clone()` this replaced disappeared, so serving a read allocates
//! nothing.

use std::fs::File;
use std::io;
use std::io::SeekFrom;

use noodles_core::Position as RecordPosition;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::{Flags, MappingQuality};

use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::read::input::read_groups::{ReadGroupResolution, RecordOwner};
use crate::ng::types::ReadGroupId;

/// One CRAM container, decoded and held in two flat buffers.
///
/// **It is what is currently being served, not a cache.** The reader above decodes one, hands
/// its records out in order, and drops it when the next one is decoded. What survives between
/// regions is the cursor's *reads*, one layer up (spec §5).
pub(crate) struct DecodedContainer {
    /// Records of this container that belong to **another sample**, counted at decode and not
    /// otherwise kept.
    ///
    /// **They are counted here rather than stepped over later, and that is a real difference
    /// from the BAM path.** A foreign record is dropped the moment its read group is known, so
    /// nothing of it is built and nothing of it is stored — which is the saving. The price is
    /// that this number is *container-granular*: a container holds around ten thousand records
    /// while a region is a few hundred bases, so it can run ahead of where a walk has actually
    /// reached. The BAM path's is exact, because it steps over records one at a time. Both are
    /// honest counts of foreign records met; they are not the same quantity.
    other_sample_records: u64,
    /// One entry per record of the container, in file order — *not* filtered to any region.
    /// Which records a region wants is decided above this, by the layer that knows the region.
    index: Vec<PackedReadEntry>,
    /// Every record's name, sequence and quality scores, back to back. Sliced by the offsets in
    /// [`PackedReadEntry`]; meaningless on its own.
    payload: Vec<u8>,
    /// Every record's CIGAR operations, back to back. Separate from [`payload`](Self::payload)
    /// because an `Op` is not a byte and packing it into one would mean encoding and decoding it.
    cigar_ops: Vec<Op>,
}

/// One read in its packed form: where its bytes are in the container's flat buffers, and every
/// scalar field of it that anything reads.
///
/// Fixed size, so the whole index is one allocation. The scalars are the noodles types rather
/// than raw integers: they are all `Copy`, so nothing is gained by unpacking them, and keeping
/// them means [`DecodedContainer::fill_raw_read`] hands them back without conversion.
#[derive(Clone, Copy)]
struct PackedReadEntry {
    /// Which read group this record belongs to, resolved **once per decode** and before the
    /// record was built.
    ///
    /// A [`ReadGroupId`] rather than a [`RecordOwner`], because a container keeps only this
    /// sample's records: a foreign one is counted and dropped where its group is decided, so
    /// "belongs to another sample" is not a state anything here can be in.
    ///
    /// The read group is also the *only* thing anything reads out of a record's auxiliary tags
    /// (`AlignedRead` carries no tag field), so deciding it here is what lets them be dropped —
    /// 6.2 MiB per open file, measured, the single largest item in what an open file cost.
    owner: ReadGroupId,
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
    /// How many records this container held that belong to another sample — counted at decode
    /// and not otherwise kept. See the field.
    pub(crate) fn other_sample_records(&self) -> u64 {
        self.other_sample_records
    }

    /// How many of this sample's records the container holds.
    pub(crate) fn len(&self) -> usize {
        self.index.len()
    }

    // `read_group(i)` was here and went with B2. It existed for exactly one caller — the CRAM
    // arm stamping the buffer on the line after `fill_raw_read` — and `fill_raw_read` sets both
    // halves itself now, so nothing outside this file needs to ask which group an entry has.

    /// Append one decoded record: its bytes move into the flat buffers and its scalars into the
    /// index, and the `RecordBuf` itself is dropped by the caller.
    ///
    /// The transient `RecordBuf` is still built — noodles' conversion has no field-selective
    /// form, so the allocations happen either way. What changes is that they do not *survive*:
    /// the churn is unchanged, the residency is not.
    fn push(&mut self, record: &RecordBuf, owner: ReadGroupId) -> io::Result<()> {
        let name = record
            .name()
            .map(|name| self.append_bytes(AsRef::<[u8]>::as_ref(name)))
            .transpose()?;
        let sequence = self.append_bytes(record.sequence().as_ref())?;
        let quality_scores = self.append_bytes(record.quality_scores().as_ref())?;
        let cigar_start = self.cigar_ops.len();
        self.cigar_ops.extend_from_slice(record.cigar().as_ref());
        let cigar = Span::new(cigar_start, self.cigar_ops.len())?;

        self.index.push(PackedReadEntry {
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

    /// Rebuild read `i` **into `out`**, reusing its allocations — **both halves of it**.
    ///
    /// A raw aligned read is a record *and* the read group it belongs to, and this sets both.
    /// It used to take a bare `RecordBuf` and fill only the record, leaving its one caller to
    /// stamp the group on the next line — so a function named for a raw aligned read filled
    /// half of one, and a second caller could have taken the record and not known to ask for
    /// the rest.
    ///
    /// **The read group comes from here because on CRAM it is decided at decode.** A CRAM
    /// stores it as a container-level number rather than a per-record `RG` tag, so it is
    /// resolved once while the container is decoded and travels with the entry — this arm is
    /// the documented exception to the readers' "records come out raw, read group cleared"
    /// contract (`aligned_reads_reader/mod.rs`). Setting it here rather than at the call site
    /// puts the exception in one place instead of two.
    ///
    /// Every owned field is cleared and refilled rather than replaced, so a walk that serves a
    /// million reads through one buffer allocates for the longest read it meets and nothing
    /// after — the same buffer-reuse property the BAM arm gets from noodles'
    /// `read_record_buf`.
    pub(crate) fn fill_raw_read(&self, i: usize, raw_read: &mut NoodlesRawAlignedRead) {
        let entry = &self.index[i];

        // Destructured, not reached field by field: this function's doc promises it fills
        // **both halves**, and a third field added to `NoodlesRawAlignedRead` would otherwise
        // compile silently here and leave that promise vouching for something false.
        let NoodlesRawAlignedRead { record, read_group } = raw_read;
        *read_group = Some(entry.owner);
        let out = record;

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

/// Seek to a container and decode it into a [`DecodedContainer`].
///
/// `Ok(None)` at end of stream (`read_container` reads 0 — the EOF marker).
pub(crate) fn decode_container_at(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    resolution: &ReadGroupResolution,
    offset: u64,
) -> io::Result<Option<DecodedContainer>> {
    reader.seek(SeekFrom::Start(offset))?;

    let mut container = cram::io::reader::Container::default();
    if reader.read_container(&mut container)? == 0 {
        return Ok(None);
    }

    // **The auxiliary tags do not survive this function.** A `RecordBuf` owns every
    // `@RG`/`MD`/`NM`/… field the record carried, and this container is held for as long as it
    // is being served — so those fields would be resident alongside it, once per open file.
    // Measured on the tomato cohort they are 6.2 MiB of the 12.7 MiB an open file costs: the
    // single largest item (`examples/dhat_ng_open_files.rs`).
    //
    // Nothing in ng reads a tag except `resolve_read_group`, and only the `RG` one.
    // `AlignedRead`, the only thing a record ever becomes, carries no tag field, so no other
    // consumer can be reading one. Resolving the read group **here**, while the tags are
    // still in hand, is therefore what lets every tag be dropped.
    //
    // **One behaviour change, and it is only visible on a malformed file.** A record whose
    // `RG` is unreadable fails when its container is decoded rather than when that record is
    // served, which can be earlier and can involve records outside the region asked for. The
    // condition was already fatal either way (`ReadGroupResolution::PerRecord`) — this reaches
    // it sooner. `a_cursor_keeps_every_read_group_of_its_sample_not_just_one`
    // (`input/mod.rs`) covers the well-formed path through here.
    let compression_header = container.compression_header()?;
    let mut decoded = DecodedContainer {
        other_sample_records: 0,
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
            repository.clone(),
            header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
        )? {
            // **Who owns this record is decided before anything is built.** A CRAM stores
            // the read group as a number, so asking costs an index lookup into the header
            // and no allocation — where building the `RecordBuf` below copies the name,
            // the bases, the qualities, the CIGAR and every auxiliary tag. A record
            // belonging to another sample is therefore counted and dropped without any of
            // that ever happening. On a single-read-group file the question is not asked
            // at all: there is one group a record could be in.
            let RecordOwner::Mine(owner) = owner_of_cram_record(&record, resolution)? else {
                decoded.other_sample_records += 1;
                continue;
            };
            let record = RecordBuf::try_from_alignment_record(header, &record)?;
            decoded.push(&record, owner)?;
        }
    }

    // The three buffers grew by doubling, so each can be holding up to twice what it needs —
    // and this container is then resident for as long as it is served. Measured, the slack was
    // 1.6 MiB of 5.4 MiB per open file, so paying one copy per container to drop it is
    // strongly worth it: a container is decoded once and read from thousands of times.
    decoded.shrink_to_fit();

    Ok(Some(decoded))
}

/// Which read group a **borrowed** CRAM record belongs to, without building anything.
///
/// A CRAM does not store a read-group string per record: it stores a *number*, an index into
/// the header's `@RG` list, and noodles turns that number into the declared name on demand
/// without copying it. So this question can be answered while the record is still borrowed
/// from the decompressed block — before the copy that turns it into a `RecordBuf`.
///
/// **A file declaring one read group is not asked at all.** Every record in it belongs to that
/// group whatever it says, which is also what lets a file re-headered without rewriting its
/// records be read ([`ReadGroupResolution::Sole`]).
fn owner_of_cram_record(
    record: &cram::Record<'_>,
    resolution: &ReadGroupResolution,
) -> io::Result<RecordOwner> {
    use noodles_sam::alignment::Record as _;
    use noodles_sam::alignment::record::data::field::{Tag, Value};

    if resolution.every_record_is_mine() {
        return resolution
            .owner_of(None)
            .map_err(|unresolved| io::Error::new(io::ErrorKind::InvalidData, unresolved));
    }

    let data = record.data();
    let name = match data.get(&Tag::READ_GROUP).transpose()? {
        Some(Value::String(name)) => Some(name),
        // A non-string `RG` is a malformed record, not an absent tag; reading it as absent
        // would report a broken file as a missing tag.
        Some(_) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "a CRAM record has an RG field that is not a string",
            ));
        }
        None => None,
    };
    resolution
        .owner_of(name.as_ref().map(AsRef::<[u8]>::as_ref))
        .map_err(|unresolved| io::Error::new(io::ErrorKind::InvalidData, unresolved))
}
