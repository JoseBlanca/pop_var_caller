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
    ///
    /// **The `?`s here are gated on a 4 GiB payload and are therefore not reachable from a
    /// test** — building the input costs 4 GiB of memory. The refusal they propagate is pinned
    /// where the truncation would happen instead, on [`Span::new`], which a test can call
    /// directly. Recorded so a coverage audit does not re-derive it.
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
    /// `read_record_buf`. **The buffer really does arrive with a history**: one is refilled for
    /// a whole pass, so every read after the first carries the previous one, and the clears are
    /// what stands between them.
    ///
    /// # Panics
    ///
    /// If `i` is not an index of this container — a caller bug, since `len` is what bounds it.
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
    // **One buffer for every record in this container, not one per record.** The conversion
    // below is `RecordBuf::default()` followed by `try_clone_from_alignment_record`, and a
    // fresh destination builds each field from capacity zero — the name, the CIGAR, and a
    // quality-score push loop that doubles 0 → 4 → … → 256 for a 150-base read. Reused, the
    // clone `clear()`s and refills buffers that are already the right size, so the growth is
    // paid once per container rather than once per record.
    let mut record_buf = RecordBuf::default();
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
            record_buf.try_clone_from_alignment_record(header, &record)?;
            decoded.push(&record_buf, owner)?;
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

#[cfg(test)]
mod tests {
    //! **The packing round trip, driven without a CRAM.**
    //!
    //! [`decode_container_at`] is the only thing that builds a [`DecodedContainer`] in
    //! production, and it needs a real file. But the part worth testing is not the decode — it is
    //! [`DecodedContainer::push`] flattening a record into the two buffers and
    //! [`DecodedContainer::fill_raw_read`] rebuilding it, which is where a record can lose a
    //! field or inherit the previous one's. Both are reachable from here, so these drive the
    //! round trip directly and the CRAM walks in `open_bam.rs` cover the decode.
    //!
    //! **Four things had never been checked** — found reviewing the step that made
    //! `fill_raw_read` fill both halves of a raw aligned read, and deferred to here:
    //!
    //! - a record with no name, which is not the same as one with an empty name;
    //! - a record whose sequence, qualities and CIGAR are all empty;
    //! - the clear-and-refill promise the function's own doc makes;
    //! - `out.data_mut().clear()`, whose doc names its own silent failure mode and whose
    //!   deletion left the whole suite green.
    //!
    //! **The buffer already has a history, and an earlier draft of this note said otherwise.**
    //! It claimed all four were latent because every caller passes a fresh buffer, and that the
    //! step moving the filtering loop would be what first hands this function a used one. That is
    //! false: `ReadFilter` refills **one** `NoodlesRawAlignedRead` for a whole pass, so every read
    //! after the first already arrives carrying the previous one. Measured by instrumenting this
    //! function and running the CRAM cursor walk — read 0 arrives empty and every read after it
    //! arrives with 30 bases in place.
    //!
    //! So what was missing is not the condition but the check. **No test anywhere compares a
    //! served read's bases against an independent expectation**, which is why deleting
    //! `sequence.clear()` grows sequences past 300,000 bases while the rest of the suite passes —
    //! the CRAM-versus-BAM oracle included. A regression in these clears would corrupt production
    //! reads today, silently, and only this module would catch it.

    use super::*;
    use crate::ng::types::ReadGroupId;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    fn empty_container() -> DecodedContainer {
        DecodedContainer {
            other_sample_records: 0,
            index: Vec::new(),
            payload: Vec::new(),
            cigar_ops: Vec::new(),
        }
    }

    /// A container holding `records`, each attributed to the read group beside it — the same
    /// path `decode_container_at` takes, minus the file.
    fn container_of(records: &[(RecordBuf, ReadGroupId)]) -> DecodedContainer {
        let mut container = empty_container();
        for (record, owner) in records {
            container
                .push(record, *owner)
                .expect("the fixture is small enough to address");
        }
        container
    }

    /// A fully-populated record, so a field dropped in the round trip has somewhere to show up.
    /// Every scalar is given a *distinct* value: a packing that transposed two of them would
    /// survive a fixture that left them at their defaults.
    fn full_record(name: &str, bases: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::from(0x0043))
            .set_reference_sequence_id(1usize)
            .set_alignment_start(
                RecordPosition::try_from(37usize).expect("a valid 1-based position"),
            )
            .set_mapping_quality(MappingQuality::new(57).expect("a mapping quality in range"))
            .set_mate_reference_sequence_id(0usize)
            .set_mate_alignment_start(
                RecordPosition::try_from(91usize).expect("a valid 1-based position"),
            )
            .set_template_length(-249)
            .set_cigar([Op::new(Kind::Match, bases.len())].into_iter().collect())
            .set_sequence(Sequence::from(bases.to_vec()))
            .set_quality_scores(QualityScores::from(vec![31u8; bases.len()]))
            .build()
    }

    /// Read `i` of the container, rebuilt into a fresh buffer of its own.
    fn raw_read_from(container: &DecodedContainer, i: usize) -> NoodlesRawAlignedRead {
        let mut raw_read = NoodlesRawAlignedRead::default();
        container.fill_raw_read(i, &mut raw_read);
        raw_read
    }

    /// **Every field survives the flattening**, which is what makes the packed form a
    /// representation rather than a lossy summary.
    ///
    /// The whole record is compared, not a field list: a scalar added to `PackedReadEntry` and
    /// forgotten in either direction fails here without anyone remembering to extend an
    /// assertion.
    #[test]
    fn a_record_and_its_read_group_survive_the_round_trip_through_the_packed_form() {
        let record = full_record("r", b"ACGTACGTAC");
        let container = container_of(&[(record.clone(), ReadGroupId(3))]);

        let raw_read = raw_read_from(&container, 0);
        assert_eq!(raw_read.record, record);
        assert_eq!(
            raw_read.read_group,
            Some(ReadGroupId(3)),
            "the read group is the other half of a raw aligned read, and this function sets both",
        );
    }

    /// **A record with no name comes back with no name**, into a buffer that held one.
    ///
    /// `None` and an empty name are different states, and the packed form keeps them apart
    /// (`PackedReadEntry::name` is an `Option<Span>`). The buffer is the caller's and outlives
    /// any one record, so a rebuild that only handled the `Some` arm would leave this record
    /// wearing the previous one's name — a read attributed to another read's identifier, which
    /// is exactly the class of fault the reused buffer invites.
    #[test]
    fn an_unnamed_record_does_not_inherit_the_name_left_in_the_buffer() {
        let mut unnamed = full_record("ignored", b"ACGTACGTAC");
        *unnamed.name_mut() = None;
        let container = container_of(&[
            (full_record("named", b"ACGTACGTAC"), ReadGroupId(0)),
            (unnamed, ReadGroupId(0)),
        ]);

        // One buffer, both records — which is the only way this can fail.
        let mut raw_read = NoodlesRawAlignedRead::default();
        container.fill_raw_read(0, &mut raw_read);
        assert_eq!(
            raw_read.record.name().map(|name| name.to_vec()),
            Some(b"named".to_vec()),
        );

        container.fill_raw_read(1, &mut raw_read);
        assert_eq!(
            raw_read.record.name(),
            None,
            "the unnamed record kept the name the buffer was carrying",
        );
    }

    /// An **empty** name is not an absent one, and the round trip must not collapse the two.
    #[test]
    fn an_empty_name_stays_empty_rather_than_becoming_absent() {
        let mut empty_name = full_record("ignored", b"ACGTACGTAC");
        *empty_name.name_mut() = Some(Vec::new().into());
        let container = container_of(&[(empty_name, ReadGroupId(0))]);

        assert_eq!(
            raw_read_from(&container, 0)
                .record
                .name()
                .map(|name| name.to_vec()),
            Some(Vec::new()),
            "an empty name was reported as no name at all",
        );
    }

    /// **A record whose sequence, qualities and CIGAR are all empty round-trips**, into a buffer
    /// that held a full record.
    ///
    /// A zero-length byte range is a legal one and copying it is a no-op — which means the
    /// *clears* are the only thing standing between this record and the previous one's contents.
    /// A record with no sequence is what a SAM line carrying `*` for SEQ decodes to, so this is a
    /// real shape rather than a constructed one.
    #[test]
    fn an_empty_record_inherits_no_bases_qualities_or_cigar_from_the_one_before_it() {
        let empty = RecordBuf::builder()
            .set_name("empty")
            .set_flags(Flags::from(0x0004))
            .set_reference_sequence_id(1usize)
            .set_alignment_start(
                RecordPosition::try_from(5usize).expect("a valid 1-based position"),
            )
            .build();
        let container = container_of(&[
            (full_record("full", b"ACGTACGTAC"), ReadGroupId(0)),
            (empty, ReadGroupId(0)),
        ]);

        let mut raw_read = NoodlesRawAlignedRead::default();
        container.fill_raw_read(0, &mut raw_read);
        assert_eq!(raw_read.record.sequence().as_ref().len(), 10);

        container.fill_raw_read(1, &mut raw_read);
        assert!(
            raw_read.record.sequence().as_ref().is_empty(),
            "the empty record kept the previous one's bases",
        );
        assert!(
            raw_read.record.quality_scores().as_ref().is_empty(),
            "the empty record kept the previous one's base qualities",
        );
        assert!(
            raw_read.record.cigar().as_ref().is_empty(),
            "the empty record kept the previous one's CIGAR",
        );
    }

    /// **The clear-and-refill promise the function's own doc makes**, stated as a test: a record
    /// served after a longer one through *one* buffer must keep no tail of the longer one.
    ///
    /// This is the property that makes the reuse safe, and nothing exercised it. It is not a
    /// future condition either — a walk hits it on its second record, because one buffer serves a
    /// whole pass (see this module's doc).
    #[test]
    fn a_shorter_record_keeps_no_tail_of_the_longer_one_before_it() {
        let container = container_of(&[
            (full_record("long", b"ACGTACGTACGTACGTACGT"), ReadGroupId(0)),
            (full_record("short", b"TTTT"), ReadGroupId(0)),
        ]);

        let mut raw_read = NoodlesRawAlignedRead::default();
        container.fill_raw_read(0, &mut raw_read);
        container.fill_raw_read(1, &mut raw_read);

        assert_eq!(raw_read.record.sequence().as_ref(), b"TTTT");
        assert_eq!(raw_read.record.quality_scores().as_ref(), &[31u8; 4]);
        assert_eq!(
            raw_read.record.name().map(|name| name.to_vec()),
            Some(b"short".to_vec()),
        );
        assert_eq!(
            raw_read.record.cigar().as_ref(),
            [Op::new(Kind::Match, 4)],
            "the short record kept the long one's CIGAR operations",
        );
    }

    /// **Serving a second read reuses the first one's allocations**, which is the property the
    /// packed form exists for and the one every other test here is blind to.
    ///
    /// The tests above pin that reuse is *safe* — no content survives. None of them pins that
    /// reuse *happens*: replacing each field wholesale instead of clearing and refilling is
    /// behaviourally identical, leaves the whole suite green, and deletes the reason a
    /// `RecordBuf` per record was abandoned in the first place (~680 bytes across seven
    /// allocations, ~72,000 allocations per open file).
    ///
    /// Capacities are read through the `_mut()` accessors because that is where the owned
    /// buffers are; they must not change between the two serves.
    #[test]
    fn serving_a_second_read_reuses_the_first_reads_allocations() {
        let container = container_of(&[
            (
                full_record("a-long-read-name", b"ACGTACGTACGTACGTACGT"),
                ReadGroupId(0),
            ),
            (full_record("short", b"TTTT"), ReadGroupId(0)),
        ]);

        let mut raw_read = NoodlesRawAlignedRead::default();
        container.fill_raw_read(0, &mut raw_read);
        let grown = (
            raw_read.record.sequence_mut().as_mut().capacity(),
            raw_read.record.quality_scores_mut().as_mut().capacity(),
            raw_read.record.cigar_mut().as_mut().capacity(),
            raw_read
                .record
                .name_mut()
                .as_mut()
                .expect("the first record is named")
                .capacity(),
        );

        container.fill_raw_read(1, &mut raw_read);

        assert_eq!(raw_read.record.sequence_mut().as_mut().capacity(), grown.0);
        assert_eq!(
            raw_read.record.quality_scores_mut().as_mut().capacity(),
            grown.1
        );
        assert_eq!(raw_read.record.cigar_mut().as_mut().capacity(), grown.2);
        assert_eq!(
            raw_read
                .record
                .name_mut()
                .as_mut()
                .expect("the second record is named")
                .capacity(),
            grown.3,
            "the name was replaced rather than refilled, so the reuse is gone",
        );
    }

    /// The container drops every tag at decode, once the read group has been resolved out of
    /// them, so a rebuilt record has none to restore. The clear is therefore about the *buffer*,
    /// not the entry: if one ever arrived here carrying another record's tags, those tags would
    /// silently become this record's — and the read group is precisely what is read out of them
    /// elsewhere, so the fault would attribute a read to another library.
    ///
    /// Its own doc names that failure mode; deleting the line left the whole suite green.
    #[test]
    fn a_buffer_carrying_auxiliary_tags_comes_back_with_none() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let container = container_of(&[(full_record("r", b"ACGTACGTAC"), ReadGroupId(1))]);

        // Tags in the buffer: the one kind of history the CRAM arm never produces itself, since
        // it drops every tag at decode — so only the clear stands between them and this record.
        let mut raw_read = NoodlesRawAlignedRead::default();
        raw_read.record.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(b"someone-else".to_vec().into()),
        );
        assert!(
            raw_read.record.data().get(&Tag::READ_GROUP).is_some(),
            "the fixture must actually plant a tag, or this test pins nothing",
        );

        container.fill_raw_read(0, &mut raw_read);
        assert!(
            raw_read.record.data().get(&Tag::READ_GROUP).is_none(),
            "the rebuilt record kept an auxiliary tag the buffer was carrying",
        );
        assert_eq!(
            raw_read.read_group,
            Some(ReadGroupId(1)),
            "and the read group comes from the entry, never from a leftover tag",
        );
    }

    /// **Every *scalar* is read from entry `i` too, not just the byte spans.**
    ///
    /// The fixtures above pin the spans per-entry, because their records differ in name and in
    /// bases — but their scalars are identical, so a rebuild that took one of *those* from the
    /// wrong entry went unseen. Measured: reading `flags`, `reference_sequence_id`,
    /// `mapping_quality`, `mate_reference_sequence_id`, `mate_alignment_start` or
    /// `template_length` from `index[0]` left the **whole 2,856-test suite green**. Only
    /// `alignment_start` died, and only in the real-CRAM walks.
    ///
    /// That is this module's own failure mode one field further along: a read reporting another
    /// read's position, mate, mapping quality or template length, silently and without a panic.
    /// One `let entry = &self.index[i]` binding feeds all seven today, so no single-line slip
    /// produces it — the exposure is to the refactor C2 performs, which is why this module was
    /// asked for before it.
    #[test]
    fn every_scalar_is_read_from_the_entry_asked_for() {
        let first = full_record("first", b"ACGTACGTAC");
        let mut second = full_record("second", b"ACGTACGTAC");
        *second.flags_mut() = Flags::from(0x0099);
        *second.reference_sequence_id_mut() = Some(4);
        *second.alignment_start_mut() = RecordPosition::try_from(1_201usize).ok();
        *second.mapping_quality_mut() = MappingQuality::new(12);
        *second.mate_reference_sequence_id_mut() = Some(5);
        *second.mate_alignment_start_mut() = RecordPosition::try_from(1_777usize).ok();
        *second.template_length_mut() = 604;

        let container = container_of(&[
            (first.clone(), ReadGroupId(0)),
            (second.clone(), ReadGroupId(0)),
        ]);

        assert_eq!(raw_read_from(&container, 0).record, first);
        assert_eq!(raw_read_from(&container, 1).record, second);
    }

    /// **`Span::new` refuses an offset it cannot address rather than truncating it.**
    ///
    /// Truncation would not fail — it would build a span pointing into *another record's* bytes
    /// and hand those back as this record's, which is the silent wrong read the check exists to
    /// prevent and which its own doc names. Measured: replacing the checked conversion with a
    /// bare `as u32` left the whole suite green.
    ///
    /// The 4 GiB payload that would reach this through `push` is not constructible in a test, so
    /// the refusal is pinned where the truncation would happen rather than where it propagates.
    #[test]
    fn a_span_past_the_index_width_is_refused_rather_than_truncated() {
        let past_u32 = u32::MAX as usize + 1;
        assert!(Span::new(0, past_u32).is_err());
        assert!(Span::new(past_u32, past_u32 + 10).is_err());
        // …and an addressable span is still built, so the refusal is not simply always-on.
        assert_eq!(Span::new(3, 7).expect("addressable").range(), 3..7);
    }

    /// **`shrink_to_fit` gives back the slack the buffers grew by**, which nothing could see.
    ///
    /// The three buffers grow by doubling, so each can hold up to twice what it needs — and the
    /// container is then resident for as long as it is served. `decode_container_at` quantifies
    /// what the call buys (1.6 MiB of 5.4 MiB per open file, measured) and the module doc makes
    /// residency the whole reason the packed form exists, yet deleting all three calls left the
    /// suite green. To a reader who has not read that comment the call looks like dead code.
    ///
    /// Capacities are asserted as strict *decreases* rather than equal to the length, because
    /// `Vec::shrink_to_fit` is documented as possibly leaving more — an equality would be flaky
    /// under a different allocator.
    #[test]
    fn shrinking_gives_back_the_slack_the_buffers_grew_by() {
        let mut container = container_of(&[(full_record("r", b"ACGTACGTAC"), ReadGroupId(0))]);
        container.payload.reserve(4096);
        container.cigar_ops.reserve(4096);
        container.index.reserve(4096);
        let grown = container.payload.capacity();

        container.shrink_to_fit();

        assert!(container.payload.capacity() < grown);
        assert!(container.cigar_ops.capacity() < 4096);
        assert!(container.index.capacity() < 4096);
    }

    /// Each entry carries **its own** read group, so a container serving several does not stamp
    /// them all with the first one's. Asserting one group would pass on a `fill_raw_read` that
    /// ignored `i` entirely.
    #[test]
    fn each_record_is_stamped_with_its_own_read_group() {
        let container = container_of(&[
            (full_record("a", b"ACGTACGTAC"), ReadGroupId(2)),
            (full_record("b", b"ACGTACGTAC"), ReadGroupId(7)),
        ]);

        assert_eq!(
            raw_read_from(&container, 0).read_group,
            Some(ReadGroupId(2))
        );
        assert_eq!(
            raw_read_from(&container, 1).read_group,
            Some(ReadGroupId(7))
        );
    }

    /// `len` counts this sample's records, and `other_sample_records` is a separate quantity
    /// that packing never touches.
    #[test]
    fn a_container_counts_the_records_it_packed_and_charges_none_elsewhere() {
        let container = container_of(&[
            (full_record("a", b"ACGTACGTAC"), ReadGroupId(0)),
            (full_record("b", b"ACGTACGTAC"), ReadGroupId(0)),
        ]);

        assert_eq!(container.len(), 2);
        assert_eq!(
            container.other_sample_records(),
            0,
            "packing a record must not charge it to another sample",
        );
    }

    /// The flat buffers are shared, so entry *i* must read **its own** slices rather than the
    /// whole payload — the failure a single-record fixture cannot see.
    #[test]
    fn records_read_their_own_slices_of_the_shared_buffers() {
        let container = container_of(&[
            (full_record("first", b"AAAAAAAAAA"), ReadGroupId(0)),
            (full_record("second", b"CCCCCCCCCC"), ReadGroupId(0)),
            (full_record("third", b"GGGGGGGGGG"), ReadGroupId(0)),
        ]);

        assert_eq!(
            raw_read_from(&container, 0).record.sequence().as_ref(),
            b"AAAAAAAAAA"
        );
        assert_eq!(
            raw_read_from(&container, 1).record.sequence().as_ref(),
            b"CCCCCCCCCC"
        );
        assert_eq!(
            raw_read_from(&container, 2).record.sequence().as_ref(),
            b"GGGGGGGGGG"
        );
    }
}
