mod header;
pub mod records;

mod tag_streams;

use self::records::TagPolicy;
use std::{borrow::Cow, io, sync::Arc};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::Record as _};

use self::{
    header::read_header,
    records::{ExternalDataReaders, Records},
};
use super::read_block_as;
use crate::{
    Record, calculate_normalized_sequence_digest,
    container::{
        CompressionHeader, ReferenceSequenceContext,
        block::{self, ContentType},
        slice::Header,
    },
    io::BitReader,
    record::Feature,
};

/// A container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
pub struct Slice<'c> {
    header: Header,
    src: &'c [u8],
}

/// Which reference bases a slice needs, in the three shapes a slice can take.
///
/// **Three states rather than an `Option`, because the two "no single window" cases want
/// opposite things of a caller**, and folding them together is a mistake that decodes silently:
/// an unmapped slice needs no bases at all, while one spanning several reference sequences needs
/// *one window per sequence* and will otherwise reach for whole sequences from a repository.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceExtent {
    /// Every record sits on one reference sequence, and these are the first and last positions
    /// they touch — 1-based and inclusive. Fetch this and call
    /// [`Slice::records_over_window`].
    OneSequence {
        /// Index into the header's reference sequences.
        reference_sequence_id: usize,
        /// First position any record touches.
        start: Position,
        /// Last position any record touches.
        end: Position,
    },
    /// The records sit on **several** reference sequences, and the slice header names none of
    /// them: it carries no id, no start and no span. Which sequences and which stretches is
    /// answered by [`Slice::record_extents`], which decodes the records without needing any
    /// bases; the windows then go to [`Slice::records_over_windows`].
    SeveralSequences,
    /// No record is placed on a reference. Nothing is reconstructed against one, so no bases are
    /// needed and no repository is consulted.
    Unmapped,
}

/// One reference sequence's window, for a slice whose records span several
/// ([`Slice::records_over_windows`]).
#[derive(Clone, Copy, Debug)]
pub struct SequenceWindow<'w> {
    /// Which reference sequence these bases are, as an index into the header's list.
    pub reference_sequence_id: usize,
    /// The bases, starting at [`start`](Self::start).
    pub bases: &'w [u8],
    /// The position `bases[0]` is, 1-based.
    pub start: Position,
}

/// The stretch of one reference sequence a slice's records actually touch — what
/// [`Slice::record_extents`] reports, and what a caller fetches to build a [`SequenceWindow`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SequenceExtent {
    /// Which reference sequence, as an index into the header's list.
    pub reference_sequence_id: usize,
    /// First position any of this slice's records touches on it.
    pub start: Position,
    /// Last position any of them touches.
    pub end: Position,
}

/// Where a decode gets the reference bases it rebuilds each read from.
#[derive(Clone, Copy)]
enum ReferenceBases<'w> {
    /// Whole sequences, from the repository passed alongside — upstream's behaviour, and what
    /// [`Slice::records`] and [`Slice::records_discarding_tags`] use.
    Repository,
    /// One window, for a slice all of whose records sit on one reference sequence.
    OneWindow(&'w [u8], Position),
    /// One window per reference sequence, for a slice whose records span several.
    WindowPerSequence(&'w [SequenceWindow<'w>]),
    /// None at all: the records are decoded but no reference is attached to them, so their
    /// coordinates are readable and their bases are not ([`Slice::record_extents`]).
    Unresolved,
}

/// Whether a decode pass links each record to its mate.
///
/// Mate resolution walks the whole slice and names every unnamed record, which the extents pass
/// has no use for — it reads coordinates and drops the records.
#[derive(Clone, Copy, Eq, PartialEq)]
enum MateResolution {
    Resolve,
    Skip,
}

impl<'c> Slice<'c> {
    /// **Which reference bases this slice needs, before any of them are fetched.**
    ///
    /// Read from the slice header, which is parsed before any block is decoded — so a caller can
    /// fetch exactly what is wanted rather than holding a whole sequence. See
    /// [`ReferenceExtent`] for what each answer obliges the caller to do.
    ///
    /// *(Was `reference_span`, returning an `Option` that folded the unmapped and
    /// several-sequences cases together. They are not the same case, and ng shipped a decode
    /// that panicked on the second because of the fold.)*
    pub fn reference_extent(&self) -> ReferenceExtent {
        match self.header.reference_sequence_context() {
            ReferenceSequenceContext::Some(context) => ReferenceExtent::OneSequence {
                reference_sequence_id: context.reference_sequence_id(),
                start: context.alignment_start(),
                end: context.alignment_end(),
            },
            ReferenceSequenceContext::Many => ReferenceExtent::SeveralSequences,
            ReferenceSequenceContext::None => ReferenceExtent::Unmapped,
        }
    }

    /// **Which sequences this slice's records touch, and how much of each — reading no reference
    /// bases at all.**
    ///
    /// For a slice whose header says only "several reference sequences"
    /// ([`ReferenceExtent::SeveralSequences`]) this is the only way to learn what to fetch: the
    /// header carries no id, no start and no span, and the answer is in the records.
    ///
    /// **It costs a decode of the records, and that is the whole cost of supporting such a
    /// slice.** It is not a second decode of the *blocks*: `decode_blocks` has already run, and
    /// its output is what both this and [`records_over_windows`](Self::records_over_windows)
    /// read. What makes the pass possible is that a record's reference sequence, its start and
    /// its CIGAR come out of the file, and only its **bases** are rebuilt against a reference —
    /// so nothing here needs one.
    ///
    /// One entry per reference sequence any *mapped* record sits on, in ascending id order, each
    /// running from the first to the last position touched on it. An unmapped slice yields an
    /// empty `Vec`.
    ///
    /// **The records are dropped rather than returned**, deliberately: nothing was resolved for
    /// them, so asking one for its bases would panic, and a function that hands back such a
    /// record invites exactly that.
    pub fn record_extents<'h: 'c, 'ch: 'c>(
        &self,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<SequenceExtent>> {
        let records = self.read_records(
            fasta::Repository::default(),
            ReferenceBases::Unresolved,
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            if tag_streams::tags_can_be_left_unread(compression_header) {
                TagPolicy::SkipStreams
            } else {
                TagPolicy::DiscardValues
            },
            MateResolution::Skip,
        )?;

        let mut extents: Vec<SequenceExtent> = Vec::new();
        for record in &records {
            if record.bam_flags.is_unmapped() || record.cram_flags.sequence_is_missing() {
                continue;
            }
            let (Some(reference_sequence_id), Some(start), Some(end)) = (
                record.reference_sequence_id,
                record.alignment_start,
                record.alignment_end(),
            ) else {
                continue;
            };
            match extents
                .iter_mut()
                .find(|extent| extent.reference_sequence_id == reference_sequence_id)
            {
                Some(extent) => {
                    extent.start = extent.start.min(start);
                    extent.end = extent.end.max(end);
                }
                None => extents.push(SequenceExtent {
                    reference_sequence_id,
                    start,
                    end,
                }),
            }
        }
        extents.sort_by_key(|extent| extent.reference_sequence_id);
        Ok(extents)
    }

    /// The slice's records, decoded against **one window per reference sequence** — for a slice
    /// whose records span several ([`ReferenceExtent::SeveralSequences`]). Tags are dropped, as
    /// [`records_discarding_tags`](Self::records_discarding_tags).
    ///
    /// `windows` is what [`record_extents`](Self::record_extents) asked for, fetched: one entry
    /// per reference sequence, each carrying its bases and the position the first of them is.
    /// Order does not matter and extra entries are harmless.
    ///
    /// **A record whose sequence is missing from `windows`, or whose span its window does not
    /// cover, is refused with both spans named** — never decoded against whatever bases happen
    /// to be to hand. That refusal is the only guard here: a slice spanning several sequences
    /// carries no reference MD5 to check a window against, where the single-sequence path has
    /// one.
    pub fn records_over_windows<'h: 'c, 'ch: 'c>(
        &self,
        windows: &'c [SequenceWindow<'c>],
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        self.read_records(
            fasta::Repository::default(),
            ReferenceBases::WindowPerSequence(windows),
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            if tag_streams::tags_can_be_left_unread(compression_header) {
                TagPolicy::SkipStreams
            } else {
                TagPolicy::DiscardValues
            },
            MateResolution::Resolve,
        )
    }

    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    #[allow(clippy::type_complexity)]
    pub fn decode_blocks(
        &self,
    ) -> io::Result<(Cow<'c, [u8]>, Vec<(block::ContentId, Cow<'c, [u8]>)>)> {
        let mut src = self.src;

        let block = read_block_as(&mut src, ContentType::CoreData)?;
        let core_data_src = block.decode()?;

        let external_data_block_count = self.header.block_count() - 1;
        let external_data_srcs = (0..external_data_block_count)
            .map(|_| {
                let block = read_block_as(&mut src, ContentType::ExternalData)?;
                block.decode().map(|src| (block.content_id, src))
            })
            .collect::<io::Result<_>>()?;

        Ok((core_data_src, external_data_srcs))
    }

    /// Reads and returns a list of raw records in this slice.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram::{self as cram, io::reader::Container};
    /// use noodles_fasta as fasta;
    ///
    /// let data = [];
    /// let mut reader = cram::io::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     let compression_header = container.compression_header()?;
    ///
    ///     for result in container.slices() {
    ///         let slice = result?;
    ///
    ///         let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
    ///
    ///         let records = slice.records(
    ///             fasta::Repository::default(),
    ///             &header,
    ///             &compression_header,
    ///             &core_data_src,
    ///             &external_data_srcs,
    ///         )?;
    ///
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        self.read_records(
            reference_sequence_repository,
            ReferenceBases::Repository,
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            TagPolicy::Keep,
            MateResolution::Resolve,
        )
    }

    /// The slice's records, **with every auxiliary tag decoded and then thrown away**.
    ///
    /// A CRAM's tags are interleaved with everything else in the same encoded streams, so they
    /// cannot be left unread — reading past them is what keeps the streams in step. What can be
    /// skipped is turning each one into a value and storing it, which is a heap allocation per
    /// record with a non-empty tag set and, on a whole-genome CRAM, most of what decoding a
    /// record costs.
    ///
    /// **A record from here reports no tags at all.** `data()` is empty, and a caller that asks
    /// it for a tag gets `None` rather than an error — so this is only for a caller that reads
    /// none of them. The one exception is the read group, which a CRAM stores as a number
    /// rather than a tag and which `data().get(&Tag::READ_GROUP)` answers from that number; it
    /// is unaffected.
    pub fn records_discarding_tags<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        self.read_records(
            reference_sequence_repository,
            ReferenceBases::Repository,
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            if tag_streams::tags_can_be_left_unread(compression_header) {
                TagPolicy::SkipStreams
            } else {
                TagPolicy::DiscardValues
            },
            MateResolution::Resolve,
        )
    }

    /// The slice's records, decoded against a **window** of the reference rather than a whole
    /// contig — and with the auxiliary tags dropped, as
    /// [`records_discarding_tags`](Self::records_discarding_tags).
    ///
    /// `window` is the bases from `window_start` onward, and it must cover every position
    /// [`reference_span`](Self::reference_span) named or the decode is refused. Nothing else
    /// about it is assumed: it may start anywhere and may be longer than the span.
    ///
    /// **This is what lets a caller reading a genome in order hold megabases instead of a
    /// chromosome.** A CRAM stores a mapped read as differences from the reference, so the
    /// bases under the read have to be in hand — but only those, and a coordinate-ordered walk
    /// knows which they are one slice ahead. On a human chromosome 1 the difference between
    /// this and a whole contig is 249 MB of resident memory.
    ///
    /// Falls back to `repository` for a slice that names no single reference sequence — an
    /// unmapped slice needs none, and a multi-reference slice resolves its bases per record.
    pub fn records_over_window<'h: 'c, 'ch: 'c>(
        &self,
        window: &'c [u8],
        window_start: Position,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        self.read_records(
            reference_sequence_repository,
            ReferenceBases::OneWindow(window, window_start),
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            if tag_streams::tags_can_be_left_unread(compression_header) {
                TagPolicy::SkipStreams
            } else {
                TagPolicy::DiscardValues
            },
            MateResolution::Resolve,
        )
    }

    fn read_records<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        bases: ReferenceBases<'c>,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
        tag_policy: TagPolicy,
        mate_resolution: MateResolution,
    ) -> io::Result<Vec<Record<'c>>> {
        let core_data_reader = BitReader::new(core_data_src);

        let mut external_data_readers = ExternalDataReaders::new();

        for (block_content_id, src) in external_data_srcs {
            external_data_readers.insert(*block_content_id, src);
        }

        let reference_sequence_context = self.header.reference_sequence_context();
        let initial_id = self.header.record_counter();

        let mut reader = Records::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            initial_id,
            tag_policy,
        );

        #[cfg(feature = "perf-counters")]
        let phase_started = std::time::Instant::now();

        // **Nothing is resolved for the extents pass**, on either path: it reads coordinates,
        // and a slice-wide fetch here would consult the repository for a single-sequence slice
        // and defeat the whole point of the pass.
        let slice_reference_sequence = if matches!(bases, ReferenceBases::Unresolved) {
            None
        } else {
            get_slice_reference_sequence(
                &reference_sequence_repository,
                match bases {
                    ReferenceBases::OneWindow(window, start) => Some((window, start)),
                    _ => None,
                },
                header,
                compression_header,
                &self.header,
                external_data_srcs,
            )?
        };

        #[cfg(feature = "perf-counters")]
        {
            crate::perf::record_slice_phase(
                crate::perf::SlicePhase::Reference,
                phase_started.elapsed().as_nanos() as u64,
            );
            crate::perf::record_slice_records(self.header.record_count());
        }

        let substitution_matrix = compression_header.preservation_map().substitution_matrix();

        #[cfg(feature = "perf-counters")]
        let phase_started = std::time::Instant::now();

        let mut records = vec![Record::default(); self.header.record_count()];

        #[cfg(feature = "perf-counters")]
        {
            crate::perf::record_slice_phase(
                crate::perf::SlicePhase::Allocate,
                phase_started.elapsed().as_nanos() as u64,
            );
        }

        #[cfg(feature = "perf-counters")]
        let phase_started = std::time::Instant::now();

        for record in &mut records {
            reader.read_record(record)?;

            record.header = Some(header);

            if !record.bam_flags.is_unmapped() && !record.cram_flags.sequence_is_missing() {
                record.reference_sequence = if reference_sequence_context.is_many() {
                    match bases {
                        // **The multi-sequence path, and the only one that needs a window per
                        // record.** The slice header names no sequence, so each record answers
                        // for itself, and its window is looked up by the sequence it names.
                        ReferenceBases::WindowPerSequence(windows) => {
                            Some(window_for_record(windows, header, record)?)
                        }
                        // No reference at all: the extents pass, which reads coordinates.
                        ReferenceBases::Unresolved => None,
                        // Whole sequences from the repository — upstream's behaviour. An *empty*
                        // repository panics here rather than erroring, which is upstream's
                        // `expect`; a caller that cannot serve whole sequences uses
                        // `records_over_windows` instead.
                        ReferenceBases::Repository | ReferenceBases::OneWindow(..) => {
                            get_record_reference_sequence(
                                &reference_sequence_repository,
                                header,
                                record,
                            )?
                        }
                    }
                } else {
                    slice_reference_sequence.clone()
                };

                record.substitution_matrix = substitution_matrix.clone();
            }
        }

        #[cfg(feature = "perf-counters")]
        {
            crate::perf::record_slice_phase(
                crate::perf::SlicePhase::ReadRecords,
                phase_started.elapsed().as_nanos() as u64,
            );
        }

        #[cfg(feature = "perf-counters")]
        let phase_started = std::time::Instant::now();

        if mate_resolution == MateResolution::Resolve {
            resolve_mates(&mut records)?;
        }

        #[cfg(feature = "perf-counters")]
        {
            crate::perf::record_slice_phase(
                crate::perf::SlicePhase::ResolveMates,
                phase_started.elapsed().as_nanos() as u64,
            );
        }

        Ok(records)
    }
}

pub fn read_slice<'c>(src: &mut &'c [u8]) -> io::Result<Slice<'c>> {
    let header = read_header(src)?;
    Ok(Slice { header, src })
}

fn resolve_mates(records: &mut [Record]) -> io::Result<()> {
    let mut mate_indices: Vec<_> = records
        .iter()
        .enumerate()
        .map(|(i, record)| record.mate_distance.map(|len| i + len + 1))
        .collect();

    for i in 0..records.len() {
        let record = &mut records[i];

        if record.name.is_none() {
            let name = record.id.to_string().into_bytes();
            record.name = Some(Cow::from(name));
        }

        if mate_indices[i].is_none() {
            continue;
        }

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let mid = j + 1;
            let (left, right) = records.split_at_mut(mid);

            let record = &mut left[j];
            let mate = &mut right[mate_index - mid];
            set_mate(record, mate);

            if mate.name.is_none() {
                mate.name = record.name.clone();
            }

            j = mate_index;
        }

        let (left, right) = records.split_at_mut(j);
        let record = &mut right[0];
        let mate = &mut left[i];
        set_mate(record, mate);

        // "The TLEN field is positive for the leftmost segment of the template, negative for the
        // rightmost, and the sign for any middle segment is undefined. If segments cover the same
        // coordinates then the choice of which is leftmost and rightmost is arbitrary..."
        let template_length = calculate_template_length(record, mate);
        records[i].template_length = template_length;

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let record = &mut records[mate_index];
            record.template_length = -template_length;
            mate_indices[j] = None;
            j = mate_index;
        }
    }

    Ok(())
}

fn set_mate(record: &mut Record, mate: &mut Record) {
    set_mate_chunk(
        &mut record.bam_flags,
        &mut record.mate_reference_sequence_id,
        &mut record.mate_alignment_start,
        mate.bam_flags,
        mate.reference_sequence_id,
        mate.alignment_start,
    );
}

fn calculate_template_length(record: &Record, mate: &Record) -> i32 {
    calculate_template_length_chunk(
        record.alignment_start,
        record.read_length,
        &record.features,
        mate.alignment_start,
        mate.read_length,
        &mate.features,
    )
}

fn set_mate_chunk(
    record_bam_flags: &mut sam::alignment::record::Flags,
    record_mate_reference_sequence_id: &mut Option<usize>,
    record_mate_alignment_start: &mut Option<Position>,
    mate_bam_flags: sam::alignment::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
) {
    if mate_bam_flags.is_reverse_complemented() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_flags.is_unmapped() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
    }

    *record_mate_reference_sequence_id = mate_reference_sequence_id;
    *record_mate_alignment_start = mate_alignment_start;
}

// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.4.9 "TLEN"
fn calculate_template_length_chunk(
    record_alignment_start: Option<Position>,
    record_read_length: usize,
    record_features: &[Feature],
    mate_alignment_start: Option<Position>,
    mate_read_length: usize,
    mate_features: &[Feature],
) -> i32 {
    use crate::record::calculate_alignment_span;

    fn alignment_end(
        alignment_start: Option<Position>,
        read_length: usize,
        features: &[Feature],
    ) -> Option<Position> {
        alignment_start.and_then(|start| {
            let span = calculate_alignment_span(read_length, features);
            let end = usize::from(start) + span - 1;
            Position::new(end)
        })
    }

    let Some(start) = record_alignment_start
        .min(mate_alignment_start)
        .map(usize::from)
    else {
        return 0;
    };

    let record_alignment_end =
        alignment_end(record_alignment_start, record_read_length, record_features);
    let mate_alignment_end = alignment_end(mate_alignment_start, mate_read_length, mate_features);

    let end = record_alignment_end
        .max(mate_alignment_end)
        .map(usize::from)
        .expect("invalid end position");

    // "...the absolute value of TLEN equals the distance between the mapped end of the template
    // and the mapped start of the template, inclusively..."
    let len = if start > end {
        start - end + 1
    } else {
        end - start + 1
    };

    i32::try_from(len).expect("invalid template length")
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ReferenceSequence<'c> {
    Embedded {
        reference_start: Position,
        sequence: &'c [u8],
    },
    External {
        sequence: Arc<fasta::record::Sequence>,
    },
    /// Bases the caller supplied, starting at `reference_start` — see
    /// [`Slice::records_over_window`]. Indexed exactly as `Embedded` is; kept apart from it
    /// because one comes out of the file and the other out of the caller, and a reader of this
    /// enum should be able to tell which.
    Window {
        reference_start: Position,
        sequence: &'c [u8],
    },
}

/// The refusal a window that does not cover the slice earns. Spelled out rather than a bare
/// "out of range", because the two spans are what a caller needs to see to fix the fetch.
fn window_too_small(
    slice_start: Position,
    slice_end: Position,
    window_start: Position,
    window_len: usize,
) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "the reference window given covers {}..{} but this slice spans {}..{}",
            usize::from(window_start),
            usize::from(window_start) + window_len.saturating_sub(1),
            usize::from(slice_start),
            usize::from(slice_end),
        ),
    )
}

fn get_slice_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    window: Option<(&'c [u8], Position)>,
    header: &sam::Header,
    compression_header: &CompressionHeader,
    slice_header: &Header,
    external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
) -> io::Result<Option<ReferenceSequence<'c>>> {
    let reference_sequence_context = slice_header.reference_sequence_context();

    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let external_reference_sequence_is_required = compression_header
        .preservation_map()
        .external_reference_sequence_is_required();

    let embedded_reference_bases_block_content_id =
        slice_header.embedded_reference_bases_block_content_id();

    if external_reference_sequence_is_required {
        // A caller-supplied window answers the same question the repository would, from
        // memory the caller already holds. It has to cover the span the slice declares —
        // checked here rather than trusted, because getting it wrong decodes every read
        // against the wrong bases and says nothing.
        if let Some((bases, window_start)) = window {
            let slice_start = context.alignment_start();
            let slice_end = context.alignment_end();
            let offset = usize::from(slice_start)
                .checked_sub(usize::from(window_start))
                .ok_or_else(|| {
                    window_too_small(slice_start, slice_end, window_start, bases.len())
                })?;
            let span = usize::from(slice_end) - usize::from(slice_start);
            if offset + span >= bases.len() {
                return Err(window_too_small(
                    slice_start,
                    slice_end,
                    window_start,
                    bases.len(),
                ));
            }

            if let Some(expected_md5) = slice_header.reference_md5() {
                validate_sequence(&bases[offset..=offset + span], expected_md5)?;
            }

            return Ok(Some(ReferenceSequence::Window {
                reference_start: window_start,
                sequence: bases,
            }));
        }

        let reference_sequence_name = header
            .reference_sequences()
            .get_index(context.reference_sequence_id())
            .map(|(name, _)| name)
            .expect("invalid slice reference sequence ID");

        let sequence = reference_sequence_repository
            .get(reference_sequence_name)
            .transpose()?
            .expect("invalid slice reference sequence name");

        // § 8.5 "Slice header block" (2024-09-04): "MD5sums should not be validated if the stored
        // checksum is all-zero."
        if let Some(expected_md5) = slice_header.reference_md5() {
            let interval = context.alignment_start()..=context.alignment_end();
            let subsequence = &sequence[interval];
            validate_sequence(subsequence, expected_md5)?;
        }

        Ok(Some(ReferenceSequence::External { sequence }))
    } else if let Some(block_content_id) = embedded_reference_bases_block_content_id {
        let sequence = external_data_srcs
            .iter()
            .find(|(id, _)| *id == block_content_id)
            .map(|(_, src)| src)
            .expect("invalid block content ID");

        Ok(Some(ReferenceSequence::Embedded {
            reference_start: context.alignment_start(),
            sequence,
        }))
    } else {
        Ok(None)
    }
}

/// The window one record of a several-sequence slice is rebuilt against.
///
/// **Every failure here is an error, never a fallback.** The bases a read is reconstructed from
/// decide what that read *says*, so serving it the wrong window — or the right window a base
/// short — produces a read that is wrong and complains about nothing. A several-sequence slice
/// carries no reference MD5 either, so this check is the whole guard.
fn window_for_record<'c>(
    windows: &'c [SequenceWindow<'c>],
    header: &sam::Header,
    record: &Record<'c>,
) -> io::Result<ReferenceSequence<'c>> {
    let name_of = |id: usize| {
        header
            .reference_sequences()
            .get_index(id)
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .unwrap_or_else(|| format!("reference sequence {id}"))
    };

    let reference_sequence_id = record.reference_sequence_id.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "a mapped record in a slice spanning several reference sequences names none",
        )
    })?;

    let window = windows
        .iter()
        .find(|window| window.reference_sequence_id == reference_sequence_id)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "no reference window was given for {}, which this slice holds records on",
                    name_of(reference_sequence_id)
                ),
            )
        })?;

    let (Some(start), Some(end)) = (record.alignment_start, record.alignment_end()) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "a mapped record in a slice spanning several reference sequences has no position",
        ));
    };

    let covers = usize::from(start) >= usize::from(window.start)
        && usize::from(end) <= usize::from(window.start) + window.bases.len() - 1;
    if !covers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "the reference window given for {} covers {}..{} but a record of this slice \
                 spans {}..{}",
                name_of(reference_sequence_id),
                usize::from(window.start),
                usize::from(window.start) + window.bases.len().saturating_sub(1),
                usize::from(start),
                usize::from(end),
            ),
        ));
    }

    Ok(ReferenceSequence::Window {
        reference_start: window.start,
        sequence: window.bases,
    })
}

fn get_record_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    record: &Record<'c>,
) -> io::Result<Option<ReferenceSequence<'c>>> {
    if record.bam_flags.is_unmapped() {
        return Ok(None);
    }

    let reference_sequence_name = record
        .reference_sequence(header)
        .transpose()?
        .map(|(name, _)| name)
        .expect("invalid reference sequence ID");

    let sequence = reference_sequence_repository
        .get(reference_sequence_name)
        .transpose()?
        .expect("invalid reference sequence name");

    Ok(Some(ReferenceSequence::External { sequence }))
}

fn validate_sequence(sequence: &[u8], expected_checksum: &[u8; 16]) -> io::Result<()> {
    let actual_checksum = calculate_normalized_sequence_digest(sequence);

    if &actual_checksum == expected_checksum {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sequence checksum mismatch: expected {expected_checksum:?}, got {actual_checksum:?}"
            ),
        ))
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;
    use crate::record::Flags;

    #[test]
    fn test_resolve_mates() -> io::Result<()> {
        let mut records = vec![
            Record {
                id: 1,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 2,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(8),
                mate_distance: Some(1),
                ..Default::default()
            },
            Record {
                id: 3,
                ..Default::default()
            },
            Record {
                id: 4,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(13),
                ..Default::default()
            },
        ];

        resolve_mates(&mut records)?;

        let name_1 = b"1".as_bstr();

        assert_eq!(records[0].name(), Some(name_1));
        assert_eq!(
            records[0].mate_reference_sequence_id,
            records[1].reference_sequence_id
        );
        assert_eq!(records[0].mate_alignment_start, records[1].alignment_start);
        assert_eq!(records[0].template_length, 12);

        assert_eq!(records[1].name(), Some(name_1));
        assert_eq!(
            records[1].mate_reference_sequence_id,
            records[3].reference_sequence_id
        );
        assert_eq!(records[1].mate_alignment_start, records[3].alignment_start);
        assert_eq!(records[1].template_length, -12);

        let name_3 = b"3".as_bstr();
        assert_eq!(records[2].name(), Some(name_3));

        assert_eq!(records[3].name(), Some(name_1));
        assert_eq!(
            records[3].mate_reference_sequence_id,
            records[0].reference_sequence_id
        );
        assert_eq!(records[3].mate_alignment_start, records[0].alignment_start);
        assert_eq!(records[3].template_length, -12);

        Ok(())
    }

    #[test]
    fn test_calculate_template_length() {
        use sam::alignment::record::Flags;

        // --> -->
        let record = Record {
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

        // --> <--
        // This is the example given in _Sequence Alignment/Map Format Specification_ (2021-06-03)
        // § 1.4.9 "TLEN" (footnote 14).
        let record = Record {
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

        // <-- -->
        let record = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

        // <-- <--
        let record = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

        // No alignment start position.
        let record = Record::default();
        assert_eq!(calculate_template_length(&record, &record), 0);
    }
}
