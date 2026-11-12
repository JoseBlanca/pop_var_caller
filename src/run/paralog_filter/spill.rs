//! **The spill — the finished VCF lines, parked on disk between the calling pass and the
//! verdict.**
//!
//! A record's *verdict* is what the filter decides about it: kept, dropped, or tagged as
//! better explained by a hidden duplication than by a variant. That decision needs the run's
//! *cut* — the likelihood-ratio threshold the operator's target false-discovery rate resolves
//! to — and the cut is a property of the whole run, computed from the distribution of every
//! record's score. So no record's verdict is known until the last one has been called. Pass
//! one therefore writes each finished record's **line** — the bytes the encoder produced — to a
//! file beside the output, together with what the scorer reads for each sample. Pass
//! two reads that file and scores it; pass three reads it again and writes the VCF, rewriting
//! only `FILTER` and `INFO`, and only on the records the verdict touches.
//!
//! **The line and not the record** (spec §3.4). A record is ten fields behind a constructor
//! that asserts every parity between them ([`VcfRecord::new`](crate::vcf::VcfRecord)); a
//! codec over it is ten field codecs plus a reconstruction that must satisfy every one of those
//! assertions, and a field the record gains later is a field the codec drops in silence. The
//! line is the record's own bytes, so byte-identity for every record the verdict does not touch
//! holds by construction, and nothing about the record's shape has to be mirrored here.
//!
//! **The layout is spec §3.4's, field for field.** An entry is two zstd frames, each behind its
//! compressed length:
//!
//! ```text
//! entry :=
//!   body_bytes       varint    -- the compressed length of the frame below
//!   body             one zstd frame, holding:
//!     contig           varint    -- the three fields the writer's ordering check reads
//!     position         varint
//!     is_repeat_tract  u8        -- 0 or 1
//!     line_length      varint
//!     line             bytes     -- the record's line, no newline
//!     sample_count     varint    -- the run's sample count, dense
//!     per sample:
//!       gc_fraction    4 bytes   -- an f32's bits, little-endian; NaN = absent
//!       mean_depth     4 bytes   -- an f32's bits, little-endian
//!       ref_reads      varint    -- AD[0]
//!       alt_reads      varint    -- every alternative's reads, summed
//!   fragments_bytes  varint    -- the compressed length of the frame below; 0 = no fragments,
//!                                 and no frame follows
//!   fragments        one zstd frame, holding:
//!     count            varint
//!     per entry, sorted by sample and then fragment:
//!       sample_step    varint    -- this entry's sample less the previous entry's (0 at the first)
//!       fragment_step  varint    -- this entry's fragment less the previous entry's of the same
//!                                   sample, or the fragment itself at a sample's first entry
//!       allele         varint    -- numbered as the record numbers its alleles
//!       error_phred    u8        -- at most LARGEST_ERROR_PHRED
//! ```
//!
//! **Why compressed.** Uncompressed, an entry is the record's VCF line as text plus its
//! fragments, and on 2,176 tomato samples that came to about 80 kB a called locus — a projected
//! 25 to 30 TB of spill for a 783 Mb genome, against 6.8 TB for the same cohort before the
//! fragments were carried. The
//! line is a genotype column per sample, most of them alike, so it shrinks several-fold under
//! zstd. The fragments are written as steps because they are sorted: a sample's fragment numbers
//! are read indices, three or four varint bytes each, while the step from one to the next is
//! mostly one byte — and zstd alone does not find that.
//!
//! **Measured on the first 6,000 entries of a real spill** (2,176 tomato samples, the start of
//! chromosome 1, 1.28 GB uncompressed), compressing each frame on its own on one thread:
//!
//! | frame | uncompressed | level 1 | level 3 |
//! |---|---|---|---|
//! | body | 105.0 kB an entry | 26.4 kB (4.0×) | 25.6 kB (4.1×) |
//! | fragments, numbers as they are | 107.9 kB | 86.9 kB (1.2×) | 38.8 kB (2.8×) |
//! | fragments, numbers as steps | 74.3 kB | 15.7 kB (6.9×) | 15.4 kB (7.0×) |
//!
//! So an entry goes from about 213 kB to about 42 kB, five times smaller. Level 1 compressed at
//! about 600 MB/s and inflated at about 1,100 MB/s a thread; level 3 gains 3% and compresses at
//! about 350 MB/s, which is why the level is 1 ([`ZSTD_COMPRESSION_LEVEL`]).
//!
//! **Why each entry on its own, and two frames to an entry.** A frame per entry is compressed by
//! whichever thread encodes the entry ([`EncodedSpillEntry::of`]) — the threads that called the
//! loci, where the run is wired for it — so the thread that writes in genome order still only
//! copies bytes. A stream compressor over the whole file would put every byte through that one
//! thread instead. The fragments are their own frame so that pass two, which never reads them,
//! steps over them by their compressed length without inflating them. Each frame records its
//! uncompressed size and a checksum: the size bounds what the reader allocates before it inflates
//! anything, and the checksum turns damage inside a frame into a refusal rather than a plausible
//! wrong record.
//!
//! **The fragments are the record's reads that reach another record** (read phase spec §2.4),
//! which the line never carries. Only pass three reads them, handing them to the writer beside
//! the line; pass two, which reads every entry to score it, steps over the frame by its length.
//! **They travel beside an entry, never inside it**: [`SpillEntry`] has no fragments field, the
//! writer takes them as their own argument, and only [`SpillReaderWithFragments`] gives them back,
//! as the second half of each item. A pass three switched to pass two's reader does not compile,
//! rather than writing every line without its reads.
//! At the top of the range they are the larger part of an entry — 3,000 samples × 32 read alleles
//! × about 6 bytes, against a line of about 150 kB (spec §4.1) — so skipping them unread is what
//! keeps pass two's reading at what it was.
//!
//! **Every record's samples carry both signals** (spec §3.2): coverage, and the split between
//! reads carrying the reference allele and reads carrying any other — repeat tracts included,
//! since the STR path spells a segment whole and writes SNPs inside it as part of an STR record's
//! alleles (owner, 2026-10-10). The tract flag is there for the writer's ordering rule alone: a tract may share
//! a position with the generic locus owning its anchor base.
//!
//! Every field inside a frame is either fixed-width, length-prefixed, or self-delimiting — a varint
//! ends at the first byte whose top bit is clear — and each frame is behind its length, so an entry
//! ends where the next begins and there is no outer frame. **No manifest, no layout check, no
//! version**: the file is created, read
//! twice and deleted inside one run, and its only reader is the process that wrote it. That is
//! what separates it from the psp, which is a user-facing artifact that has to interpret
//! itself.
//!
//! **`NaN` means *this sample has no usable window here*, and it has to survive the file**
//! (spec §6 trap 4). Both floats are written as bit patterns and read back as bit patterns,
//! never compared with `==`, so an absent pair arrives absent rather than as a zero that would
//! make a sample with no evidence look like one with average coverage.
//!
//! **The encoder destructures its input exhaustively, and matches both row shapes**, which is
//! the codebase's way of making a struct that gains a field fail to compile rather than lose it
//! quietly
//! ([`cohort_merge`](crate::run::cohort_merge)'s `render`,
//! `var_calling::types`).

use std::io::{self, BufRead, Read, Write};
use std::iter::FusedIterator;

use thiserror::Error;

use super::WindowCoverage;
use crate::psp::varint::{MAX_VARINT_BYTES, VarintError, decode_u64_leb128, encode_u64_leb128};
use crate::read_phase::capture::LARGEST_ERROR_PHRED;
use crate::read_phase::fragments::{FragmentAllele, RecordFragments};
use crate::types::{ContigId, Position};

/// The field names the decoder reports, one per field of the layout above.
///
/// Kept together so that a field renamed on [`SpillEntry`] or [`SpilledSample`] is renamed
/// once, and so no two call sites can spell the same field differently — a failure message
/// naming a field that no longer exists sends the reader somewhere there is nothing to find.
mod field {
    pub const BODY_BYTES: &str = "body_bytes";
    pub const BODY: &str = "body";
    pub const FRAGMENTS: &str = "fragments";
    pub const CONTIG: &str = "contig";
    pub const POSITION: &str = "position";
    pub const IS_REPEAT_TRACT: &str = "is_repeat_tract";
    pub const LINE_LENGTH: &str = "line_length";
    pub const LINE: &str = "line";
    pub const SAMPLE_COUNT: &str = "sample_count";
    pub const GC_FRACTION: &str = "gc_fraction";
    pub const MEAN_DEPTH: &str = "mean_depth";
    pub const REF_READS: &str = "ref_reads";
    pub const ALT_READS: &str = "alt_reads";
    pub const FRAGMENTS_BYTES: &str = "fragments_bytes";
    pub const FRAGMENT_COUNT: &str = "fragment_count";
    pub const FRAGMENT_SAMPLE: &str = "fragment_sample";
    pub const FRAGMENT: &str = "fragment";
    pub const FRAGMENT_ALLELE: &str = "fragment_allele";
    pub const FRAGMENT_ERROR_PHRED: &str = "fragment_error_phred";
}

/// The most per-sample entries reserved up front from a decoded count, so that a corrupt count
/// cannot make the decoder allocate a run's worth of memory before the bytes that would fail it
/// have even been read: a count of 4,026,531,840 would otherwise reserve 4,026,531,840 × 16
/// bytes, about 64 GB.
///
/// **8,192 samples, so 128 KiB reserved at the cap** (a [`SpilledSample`] is 16 bytes). The
/// largest cohort spec §4 contemplates is three thousand samples, which reserves its exact
/// count in one allocation and never grows. A run between 8,192 and [`MAX_SAMPLES`] decodes
/// correctly and pays `Vec` growth — this is a reservation policy, not a limit.
const MAX_SAMPLES_RESERVED_UP_FRONT: usize = 8192;

/// The largest sample count the decoder will believe.
///
/// **The reservation cap above bounds what is reserved and not what is built.** Without a limit
/// on the count itself, a corrupt one keeps the loop decoding until the file runs out, so the
/// memory held is the rest of the file rather than one entry — the same failure
/// [`MAX_LINE_BYTES`] closes on the other variable-length field, and against spec §5's one
/// entry in hand at a time.
///
/// **A million samples, so 16 MB of `SpilledSample` at the ceiling.** Spec §4's largest cohort
/// is three thousand, so this leaves room for three hundred times it — a count past here is
/// corruption, not a cohort.
const MAX_SAMPLES: usize = 1_000_000;

/// The longest record line the decoder will read.
///
/// A line is the eight fixed columns plus one genotype column per sample. At spec §4's three
/// thousand samples those columns are a few hundred kilobytes, and 16 MiB leaves room for a
/// cohort five times larger with a hundred bytes a sample — so a longer length is corruption,
/// not a record.
///
/// **Without a ceiling the field is unbounded in the worst way**: the length is read, then
/// every byte the file still holds is appended before the short-read check can fire, so on a
/// spill the size of a genome's VCF one corrupt byte turns a decode error into an
/// out-of-memory kill of a run that has already finished calling.
const MAX_LINE_BYTES: usize = 16 * 1024 * 1024;

/// The most bytes one sample's row inflates to: two floats and two `u32` varints of at most five
/// bytes each.
const MAX_SAMPLE_BYTES: usize = 4 + 4 + 5 + 5;

/// The largest body the decoder will inflate: the head's four fields at their widest, the longest
/// line, and the most samples at their widest row. **About 35 MB**, so a frame that claims more is
/// damage — and a frame records its own size, so this is checked before anything is allocated
/// for it.
const MAX_BODY_BYTES: usize =
    3 * MAX_VARINT_BYTES + 1 + MAX_LINE_BYTES + MAX_SAMPLES * MAX_SAMPLE_BYTES;

/// The longest compressed frame the decoder will read for something that inflates to at most
/// `inflated` bytes.
///
/// **Above zstd's own bound on what compressing that many bytes can produce**
/// (`ZSTD_COMPRESSBOUND`, a little over `inflated + inflated / 256`), so every frame the writer can
/// make is under it, and a length past it is damage that would otherwise have the reader hold
/// the rest of the file — the failure [`MAX_LINE_BYTES`] closes for the line.
const fn max_frame_bytes(inflated: usize) -> usize {
    inflated + inflated / 128 + 4096
}

/// The level each frame is compressed at.
///
/// **Level 1, zstd's fastest of its ordinary levels**, because the frames are compressed while
/// the run is calling and every entry passes through it. On a real spill it gives within 3% of
/// level 3's size at about twice its speed (the table at the top of this module).
const ZSTD_COMPRESSION_LEVEL: i32 = 1;

/// One record on the spill: the line the writer would have written, and what the scorer needs
/// beside it.
///
/// `contig`, `position` and `is_repeat_tract` are the fields the writer's ordering check reads,
/// carried so that pass three can run that same check without rebuilding the record
/// ([`place_of`](crate::vcf::writer), spec §6 trap 6).
#[derive(Clone, Debug)]
pub struct SpillEntry {
    /// Which contig the record is written on.
    pub contig: ContigId,
    /// The position the record is *written* at, after the padding rule has moved it.
    pub position: Position,
    /// Whether the record is a repeat tract. Part of the ordering rule: a tract may share a
    /// position with the generic locus that owns its anchor base, and nothing else may.
    pub is_repeat_tract: bool,
    /// The record's VCF line as the encoder produced it, with no trailing newline.
    pub line: Vec<u8>,
    /// What each sample showed, one entry per sample in the run's order — dense, so a sample
    /// that covered nothing still has its place.
    pub samples: Vec<SpilledSample>,
}

/// What one sample showed at a record: its window, and the split between reads carrying the
/// reference allele and reads carrying any other.
#[derive(Clone, Copy, Debug)]
pub struct SpilledSample {
    /// The sample's window at this locus, or the `NaN` pair where it has none.
    pub window: WindowCoverage,
    /// Reads supporting the reference allele — `AD[0]`.
    pub ref_reads: u32,
    /// **Every alternative's reads, summed.** Pooling is what lets a multiallelic site use this
    /// path: the collapsed-duplication story asks whether the non-reference share sits at some
    /// whole number of copies over the total, and two copies of three carrying two *different*
    /// non-reference bases give the same two-thirds as two copies carrying one (spec §3.2).
    pub alt_reads: u32,
}

/// What can go wrong reading or writing the spill.
///
/// **Every variant names the field or the operation it failed on**, because the file has no
/// self-description to fall back on: a decoder that says only "truncated" leaves the reader
/// nowhere to look. The path is not here — the writer and the reader are generic over their
/// stream — and is added by the run, which wraps these into a
/// [`RunError`](crate::run::RunError) naming the file (spec §5). The record's ordinal is
/// the caller's too; the sample's index is only this module's, so it is here.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SpillError {
    /// The sink refused the bytes.
    #[error("the paralog spill could not be written")]
    Write {
        /// What the sink said.
        source: io::Error,
    },

    /// The sink refused to flush.
    #[error("the paralog spill could not be flushed to disk")]
    Flush {
        /// What the sink said.
        source: io::Error,
    },

    /// zstd refused to compress one of an entry's frames. Not reachable from any input the
    /// encoder is given — zstd refuses only a context it could not allocate — but carried rather
    /// than unwrapped, so that a refusal reaches the run as an error naming the spill.
    #[error("the paralog spill's {field} could not be compressed: {reason}")]
    Compress {
        /// The frame being compressed.
        field: &'static str,
        /// What zstd said.
        reason: &'static str,
    },

    /// One of an entry's frames is not a frame this writer makes: zstd refused it, it does not
    /// record its own size, or it inflated to a size other than the one it records. Damage
    /// inside a frame lands here, which is what the frame's checksum is for.
    #[error("the paralog spill's {field} is not a readable zstd frame: {reason}")]
    NotAFrame {
        /// The frame that could not be inflated.
        field: &'static str,
        /// What was wrong with it.
        reason: &'static str,
    },

    /// The source failed while a field was being decoded.
    #[error("the paralog spill could not be read while decoding its {field}")]
    Read {
        /// The field being decoded when the read failed.
        field: &'static str,
        /// What the source said.
        source: io::Error,
    },

    /// The stream ended inside an entry. A spill that ends at an entry boundary is a
    /// complete spill; one that ends mid-entry is a partial final write or a framing bug.
    #[error("the paralog spill ends inside a record, while reading its {field}")]
    Truncated {
        /// The field being decoded when the bytes ran out.
        field: &'static str,
    },

    /// A variable-length integer did not encode a `u64` — it ran past the 10-byte cap, or its
    /// last byte carried bits the `u64` cannot hold. Corruption, since nothing this codec
    /// writes can produce one.
    #[error("the paralog spill's {field} is not a valid variable-length integer")]
    OverlongVarint {
        /// The field whose encoding was over-long.
        field: &'static str,
    },

    /// A decoded value does not fit the field it belongs to.
    #[error("the paralog spill's {field} holds {value}, which does not fit the field")]
    OutOfRange {
        /// The field that overflowed.
        field: &'static str,
        /// What was decoded.
        value: u64,
    },

    /// A length-prefixed section held bytes past the last entry its own count declared: its
    /// length and its count disagree, so one of them is damaged.
    #[error(
        "the paralog spill's {field} section holds {unread} byte(s) past its last entry; its \
         length and its count disagree"
    )]
    SectionLongerThanItsEntries {
        /// The section whose length and count disagree.
        field: &'static str,
        /// The bytes left over.
        unread: u64,
    },

    /// The file ended at an entry boundary holding a different number of entries than were
    /// written — fewer, if it lost its tail; more, if something else wrote to it.
    ///
    /// **Nothing in the file itself says how many it should hold** — a spill cut on a boundary
    /// is byte for byte the prefix of a longer one — so the count comes from the writer, through
    /// [`SpillReader::new`]. Without it a lost tail reads back as a complete, shorter spill,
    /// and pass three writes a VCF short by its last records with no panic and no message.
    #[error("the paralog spill holds {read} record(s) where {expected} were written")]
    TheWrongNumberOfRecords {
        /// How many the writer counted.
        expected: u64,
        /// How many the reader found.
        read: u64,
    },

    /// A field written as `0` or `1` came back as something else.
    #[error("the paralog spill's {field} holds the byte {byte}, and only 0 or 1 are written")]
    NotABoolean {
        /// The field that is neither true nor false.
        field: &'static str,
        /// The byte found there.
        byte: u8,
    },

    /// One sample of a record could not be decoded. **Carries which sample**: at a cohort of
    /// three thousand, a failure that names only the field points at three thousand places at
    /// once.
    #[error("the paralog spill could not be read at sample {index} of a record")]
    InSample {
        /// The sample's index in the run's sample order.
        index: usize,
        /// What went wrong there.
        #[source]
        source: Box<SpillError>,
    },
}

/// **Writes entries to the spill, one at a time.**
///
/// Each entry is encoded into a scratch buffer that is reused across entries and then written
/// out, so the writer holds one entry's bytes at a time however many records the run has and
/// however many samples it carries.
///
/// **Wrap the sink in a [`BufWriter`](std::io::BufWriter).** [`Self::append`] issues one
/// `write_all` per entry, so an unbuffered `File` costs a write syscall per called record —
/// the same requirement [`SpillReader`]'s [`BufRead`] bound makes unskippable.
#[must_use = "a SpillWriter that is not finished may leave its last entries unflushed"]
pub struct SpillWriter<W: Write> {
    sink: W,
    /// The reused per-entry encode buffer. It grows to the largest entry and stays there —
    /// one entry's worth, never the stream's.
    scratch: Vec<u8>,
    entries_written: u64,
}

impl<W: Write> SpillWriter<W> {
    /// Start writing to `sink`.
    pub fn new(sink: W) -> Self {
        Self {
            sink,
            scratch: Vec::new(),
            entries_written: 0,
        }
    }

    /// Append one entry, with the record's fragments beside it.
    ///
    /// # Errors
    ///
    /// If the entry cannot be compressed, or the sink refuses the bytes. An entry the sink
    /// refused is not counted.
    pub fn append(
        &mut self,
        entry: &SpillEntry,
        fragments: &RecordFragments,
    ) -> Result<(), SpillError> {
        self.scratch.clear();
        encode_entry(entry, fragments, &mut self.scratch)?;
        self.sink
            .write_all(&self.scratch)
            .map_err(|source| SpillError::Write { source })?;
        self.entries_written += 1;
        Ok(())
    }

    /// Append one entry that was encoded elsewhere ([`EncodedSpillEntry::of`]) — the same bytes
    /// [`append`](Self::append) would have written.
    ///
    /// # Errors
    ///
    /// If the sink refuses the bytes. An entry the sink refused is not counted.
    pub fn append_encoded(&mut self, entry: &EncodedSpillEntry) -> Result<(), SpillError> {
        self.sink
            .write_all(&entry.bytes)
            .map_err(|source| SpillError::Write { source })?;
        self.entries_written += 1;
        Ok(())
    }

    #[cfg(test)]
    /// How many entries have been appended.
    #[must_use]
    pub fn entries_written(&self) -> u64 {
        self.entries_written
    }

    /// Flush the sink and hand it back.
    ///
    /// **The flush is not automatic.** The writer has no `Drop`, so a writer dropped rather
    /// than finished leaves whatever the sink itself does with unflushed bytes — a
    /// `BufWriter` flushes and discards any error, an unbuffered sink has nothing to flush.
    /// **And a spill that lost its tail on an entry boundary reads back as a complete, shorter
    /// spill**: nothing in the file says how many entries it should hold, so the reader cannot
    /// tell. Call this on every path that ends a run, and compare [`Self::entries_written`]
    /// against what the reader yields when the file's lifecycle can carry the number across.
    ///
    /// # Errors
    ///
    /// If the flush fails.
    pub fn finish(mut self) -> Result<W, SpillError> {
        self.sink
            .flush()
            .map_err(|source| SpillError::Flush { source })?;
        Ok(self.sink)
    }
}

/// **Reads entries back from the spill, one at a time**, in the order they were written —
/// which is genome order, and is what makes pass three's output a function of pass one's.
///
/// Takes a [`BufRead`] rather than a [`Read`] because the variable-length integers are read a
/// byte at a time: buffered, that is a bounds check per byte; unbuffered it would be a read
/// syscall per byte.
///
/// **A failure is terminal.** A stream that failed once is not positioned at an entry boundary
/// any more, so anything decoded after it is bytes from the middle of an entry read as if they
/// were the start of one — which decodes, plausibly, into a record that was never written.
pub struct SpillReader<R: BufRead> {
    source: R,
    /// Set by the first failure; after it, the reader yields nothing.
    failed: bool,
    /// How many entries the writer counted. **An argument rather than a setting**: a reader
    /// that could be built without it would make a short file undetectable, and forgetting a
    /// call is silent where forgetting an argument does not compile.
    entries_expected: u64,
    entries_read: u64,
    /// Whether each entry's fragments are read, or stepped over unread (pass two).
    fragments: FragmentsSection,
    /// What inflating a frame needs, reused across entries.
    inflating: Inflating,
}

/// **What a reader holds to inflate frames**: one decompression context and two buffers, reused
/// across entries, so a reader allocates for its largest entry once rather than for every one.
struct Inflating {
    context: zstd::zstd_safe::DCtx<'static>,
    /// One frame's compressed bytes, as read off the stream.
    compressed: Vec<u8>,
    /// The same frame, inflated.
    inflated: Vec<u8>,
}

impl Inflating {
    fn new() -> Self {
        Self {
            context: zstd::zstd_safe::DCtx::create(),
            compressed: Vec::new(),
            inflated: Vec::new(),
        }
    }
}

/// What a reader does with each entry's fragments section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FragmentsSection {
    /// Step over it by its length: pass two's reader.
    Skip,
    /// Decode it: pass three's.
    Read,
}

impl<R: BufRead> SpillReader<R> {
    /// Start reading from `source`, which must be positioned at an entry boundary, expecting
    /// `entries_expected` of them. Each entry's fragments are stepped over unread; only
    /// [`SpillReaderWithFragments`] gives them back.
    ///
    /// **The count is required because a spill that lost its tail on an entry boundary is byte
    /// for byte the prefix of a complete one**, so nothing in the bytes can tell the two apart.
    /// It has to come from the side that wrote them, and a reader that could be built without
    /// it would make the shortfall silent again.
    pub fn new(source: R, entries_expected: u64) -> Self {
        Self::over(source, entries_expected, FragmentsSection::Skip)
    }

    fn over(source: R, entries_expected: u64, fragments: FragmentsSection) -> Self {
        Self {
            source,
            failed: false,
            entries_expected,
            entries_read: 0,
            fragments,
            inflating: Inflating::new(),
        }
    }

    /// The next entry, `None` at a clean end of stream, and `None` for ever after a failure.
    ///
    /// # Errors
    ///
    /// If the stream fails, ends inside an entry, holds a value no entry can carry, or — where
    /// `expecting` was told how many were written — ends after fewer than that. **A
    /// stream that ends *inside* an entry is an error, not an end**, and so is one that ends
    /// on a boundary too early.
    pub fn next_entry(&mut self) -> Option<Result<SpillEntry, SpillError>> {
        self.next_decoded()
            .map(|decoded| decoded.map(|(entry, _)| entry))
    }

    /// The next entry and its fragments — empty unless this reader reads them.
    fn next_decoded(&mut self) -> Option<Result<(SpillEntry, RecordFragments), SpillError>> {
        if self.failed {
            return None;
        }
        let at_an_entry_boundary_and_out_of_bytes = match self.source.fill_buf() {
            Ok(buffered) => buffered.is_empty(),
            Err(source) => {
                self.failed = true;
                return Some(Err(SpillError::Read {
                    field: field::BODY_BYTES,
                    source,
                }));
            }
        };
        if at_an_entry_boundary_and_out_of_bytes {
            if self.entries_read != self.entries_expected {
                self.failed = true;
                return Some(Err(SpillError::TheWrongNumberOfRecords {
                    expected: self.entries_expected,
                    read: self.entries_read,
                }));
            }
            return None;
        }
        let decoded = decode_entry(&mut self.source, self.fragments, &mut self.inflating);
        match &decoded {
            Ok(_) => self.entries_read += 1,
            Err(_) => self.failed = true,
        }
        Some(decoded)
    }
}

impl<R: BufRead> Iterator for SpillReader<R> {
    type Item = Result<SpillEntry, SpillError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_entry()
    }
}

impl<R: BufRead> FusedIterator for SpillReader<R> {}

/// **Pass three's reader**: each entry with the record's fragments beside it, by type — the one
/// reader that decodes them (read phase spec §4.1). Otherwise [`SpillReader`]: the same count, the
/// same refusals, the same terminal failure.
pub struct SpillReaderWithFragments<R: BufRead>(SpillReader<R>);

impl<R: BufRead> SpillReaderWithFragments<R> {
    /// As [`SpillReader::new`], decoding each entry's fragments.
    pub fn new(source: R, entries_expected: u64) -> Self {
        Self(SpillReader::over(
            source,
            entries_expected,
            FragmentsSection::Read,
        ))
    }
}

impl<R: BufRead> Iterator for SpillReaderWithFragments<R> {
    type Item = Result<(SpillEntry, RecordFragments), SpillError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next_decoded()
    }
}

impl<R: BufRead> FusedIterator for SpillReaderWithFragments<R> {}

/// **One entry, already turned into the bytes the spill stores** — so that the encoding can be
/// done on the thread that called the locus, and the thread that writes only copies bytes.
///
/// At 2,169 samples an entry is about 50 kB, most of it the record's VCF line; building that
/// line and these bytes was done on the run's one writing thread while every calling thread
/// waited for the next round.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodedSpillEntry {
    bytes: Vec<u8>,
}

impl EncodedSpillEntry {
    /// Encode `entry`, with the record's fragments beside it, as [`SpillWriter::append`] would
    /// write it.
    ///
    /// **The compression happens here**, on the thread that calls this, which is what makes it
    /// worth calling ahead of the writing thread.
    ///
    /// # Errors
    ///
    /// If zstd refuses to compress a frame ([`SpillError::Compress`]).
    pub fn of(entry: &SpillEntry, fragments: &RecordFragments) -> Result<Self, SpillError> {
        let mut bytes = Vec::new();
        encode_entry(entry, fragments, &mut bytes)?;
        Ok(Self { bytes })
    }
}

/// **What a thread holds to compress frames**: one compression context set to
/// [`ZSTD_COMPRESSION_LEVEL`], and two buffers, reused for every entry the thread encodes.
struct Deflating {
    context: zstd::zstd_safe::CCtx<'static>,
    /// One frame's uncompressed bytes, as the field encoders write them.
    raw: Vec<u8>,
    /// The same frame, compressed.
    compressed: Vec<u8>,
}

impl Deflating {
    fn new() -> Self {
        use zstd::zstd_safe::CParameter;
        let mut context = zstd::zstd_safe::CCtx::create();
        // PANIC-FREE: each of the three is a parameter zstd defines with a value in its range —
        // a level between 1 and 22, and two flags — so none can be refused.
        for parameter in [
            CParameter::CompressionLevel(ZSTD_COMPRESSION_LEVEL),
            CParameter::ContentSizeFlag(true),
            CParameter::ChecksumFlag(true),
        ] {
            context
                .set_parameter(parameter)
                .expect("a zstd parameter in its range is accepted");
        }
        Self {
            context,
            raw: Vec::new(),
            compressed: Vec::new(),
        }
    }

    /// Compress what is in `raw` into one frame, and append the frame to `out` behind its length.
    fn append_frame(&mut self, field: &'static str, out: &mut Vec<u8>) -> Result<(), SpillError> {
        self.compressed.clear();
        self.compressed
            .reserve(zstd::zstd_safe::compress_bound(self.raw.len()));
        self.context
            .compress2(&mut self.compressed, &self.raw)
            .map_err(|code| SpillError::Compress {
                field,
                reason: zstd::zstd_safe::get_error_name(code),
            })?;
        encode_u64_leb128(self.compressed.len() as u64, out);
        out.extend_from_slice(&self.compressed);
        Ok(())
    }
}

thread_local! {
    /// **One compressor a thread**: an entry is encoded on whichever thread called its locus, and
    /// a zstd context is about a megabyte of tables at this level, too much to build per entry.
    static DEFLATING: std::cell::RefCell<Deflating> = std::cell::RefCell::new(Deflating::new());
}

/// Append one entry's bytes to `out`, in the layout at the top of this module: the body's frame,
/// then the fragments' frame, or a length of 0 where the record has none.
///
/// # Errors
///
/// If zstd refuses to compress either frame.
fn encode_entry(
    entry: &SpillEntry,
    fragments: &RecordFragments,
    out: &mut Vec<u8>,
) -> Result<(), SpillError> {
    DEFLATING.with_borrow_mut(|deflating| {
        deflating.raw.clear();
        encode_body(entry, &mut deflating.raw);
        deflating.append_frame(field::BODY, out)?;
        if fragments.is_empty() {
            encode_u64_leb128(0, out);
            return Ok(());
        }
        deflating.raw.clear();
        encode_fragments(fragments, &mut deflating.raw);
        deflating.append_frame(field::FRAGMENTS, out)
    })
}

/// Append the body's fields to `out`, uncompressed: the record's head, its line and its samples.
///
/// **The entry and each sample are destructured exhaustively**, so a field added to either
/// struct stops the build here instead of vanishing from the file.
fn encode_body(entry: &SpillEntry, out: &mut Vec<u8>) {
    let SpillEntry {
        contig,
        position,
        is_repeat_tract,
        line,
        samples,
    } = entry;

    encode_u64_leb128(u64::from(contig.get()), out);
    encode_u64_leb128(position.get(), out);
    out.push(u8::from(*is_repeat_tract));
    encode_u64_leb128(line.len() as u64, out);
    out.extend_from_slice(line);
    encode_u64_leb128(samples.len() as u64, out);

    for sample in samples {
        let SpilledSample {
            window,
            ref_reads,
            alt_reads,
        } = sample;

        encode_window(window, out);
        encode_u64_leb128(u64::from(*ref_reads), out);
        encode_u64_leb128(u64::from(*alt_reads), out);
    }
}

/// Append the fragments' fields to `out`, uncompressed: their count, then each one.
///
/// The sample and the fragment are written as steps from the entry before (the layout at the top
/// of this module).
fn encode_fragments(fragments: &RecordFragments, out: &mut Vec<u8>) {
    encode_u64_leb128(fragments.len() as u64, out);
    let mut previous = Steps::default();
    for entry in fragments.entries() {
        let FragmentAllele {
            sample,
            fragment,
            allele,
            error_phred,
        } = *entry;
        let (sample_step, fragment_step) = previous.step_to(sample, fragment);
        encode_u64_leb128(u64::from(sample_step), out);
        encode_u64_leb128(u64::from(fragment_step), out);
        encode_u64_leb128(u64::from(allele), out);
        out.push(error_phred);
    }
}

/// **The sample and fragment the next fragment entry is written as a step from**: the previous
/// entry's, or sample 0 and no fragment before the first.
#[derive(Default)]
struct Steps {
    sample: u32,
    /// The previous entry's fragment, or `None` at the start of a sample, where the fragment is
    /// written whole.
    fragment: Option<u32>,
}

impl Steps {
    /// The steps that reach `sample` and `fragment` from here, and move here to them.
    ///
    /// # Panics
    ///
    /// If `sample` and `fragment` come before where this is. **They cannot**: a
    /// [`RecordFragments`] sorts its entries by sample and then fragment when it is built, and
    /// the encoder walks them in that order.
    fn step_to(&mut self, sample: u32, fragment: u32) -> (u32, u32) {
        let sample_step = sample
            .checked_sub(self.sample)
            .expect("a record's fragments are sorted by sample");
        if sample_step > 0 {
            self.fragment = None;
        }
        let fragment_step = match self.fragment {
            None => fragment,
            Some(previous) => fragment
                .checked_sub(previous)
                .expect("a sample's fragments are sorted by fragment"),
        };
        self.sample = sample;
        self.fragment = Some(fragment);
        (sample_step, fragment_step)
    }

    /// The sample and fragment `sample_step` and `fragment_step` reach from here, and move here to
    /// them.
    ///
    /// # Errors
    ///
    /// Where either passes the end of a `u32`, which the encoder cannot write.
    fn take(&mut self, sample_step: u64, fragment_step: u64) -> Result<(u32, u32), SpillError> {
        let sample = u64::from(self.sample).saturating_add(sample_step);
        let sample = u32::try_from(sample).map_err(|_| SpillError::OutOfRange {
            field: field::FRAGMENT_SAMPLE,
            value: sample,
        })?;
        if sample != self.sample {
            self.fragment = None;
        }
        let fragment = self
            .fragment
            .map_or(0, u64::from)
            .saturating_add(fragment_step);
        let fragment = u32::try_from(fragment).map_err(|_| SpillError::OutOfRange {
            field: field::FRAGMENT,
            value: fragment,
        })?;
        self.sample = sample;
        self.fragment = Some(fragment);
        Ok((sample, fragment))
    }
}

/// The window pair, by bit pattern, so an absent one survives the round trip.
fn encode_window(window: &WindowCoverage, out: &mut Vec<u8>) {
    let WindowCoverage {
        gc_fraction,
        mean_depth,
    } = window;

    out.extend_from_slice(&gc_fraction.to_bits().to_le_bytes());
    out.extend_from_slice(&mean_depth.to_bits().to_le_bytes());
}

/// Read one entry from `source`, which is positioned at its first byte, and its fragments —
/// empty when `fragments` says to step over them.
fn decode_entry<R: BufRead>(
    source: &mut R,
    fragments: FragmentsSection,
    inflating: &mut Inflating,
) -> Result<(SpillEntry, RecordFragments), SpillError> {
    read_frame(
        source,
        field::BODY_BYTES,
        field::BODY,
        max_frame_bytes(MAX_BODY_BYTES),
        &mut inflating.compressed,
    )?;
    inflating.inflate(field::BODY, MAX_BODY_BYTES)?;
    let mut body = inflating.inflated.as_slice();
    let entry = decode_body(&mut body)?;
    if !body.is_empty() {
        return Err(SpillError::SectionLongerThanItsEntries {
            field: field::BODY,
            unread: body.len() as u64,
        });
    }
    let fragments = read_fragments(source, fragments, entry.samples.len(), inflating)?;
    Ok((entry, fragments))
}

impl Inflating {
    /// Inflate the frame in `compressed` into `inflated`, refusing one that claims more than `max`
    /// bytes before anything is allocated for it.
    fn inflate(&mut self, field: &'static str, max: usize) -> Result<(), SpillError> {
        let size = match zstd::zstd_safe::get_frame_content_size(&self.compressed) {
            Ok(Some(size)) => size,
            Ok(None) => {
                return Err(SpillError::NotAFrame {
                    field,
                    reason: "the frame does not record its size",
                });
            }
            Err(_) => {
                return Err(SpillError::NotAFrame {
                    field,
                    reason: "the bytes do not begin a zstd frame",
                });
            }
        };
        if size > max as u64 {
            return Err(SpillError::OutOfRange { field, value: size });
        }
        self.inflated.clear();
        // The size is under `max`, so this fits a `usize`.
        self.inflated.reserve(size as usize);
        let written = self
            .context
            .decompress(&mut self.inflated, &self.compressed)
            .map_err(|code| SpillError::NotAFrame {
                field,
                reason: zstd::zstd_safe::get_error_name(code),
            })?;
        if written as u64 != size {
            return Err(SpillError::NotAFrame {
                field,
                reason: "the frame inflated to a size other than the one it records",
            });
        }
        Ok(())
    }
}

/// Read a frame's length, as `length_field`, and then that many bytes into `out`.
///
/// The length is refused before a byte of the frame is read if it is past `max`, so the bytes
/// held are bounded by that and not by what the file still holds.
fn read_frame<R: BufRead>(
    source: &mut R,
    length_field: &'static str,
    frame_field: &'static str,
    max: usize,
    out: &mut Vec<u8>,
) -> Result<(), SpillError> {
    let length = read_usize(source, length_field)?;
    if length > max {
        return Err(SpillError::OutOfRange {
            field: length_field,
            value: length as u64,
        });
    }
    out.clear();
    let bytes_read = source
        .by_ref()
        .take(length as u64)
        .read_to_end(out)
        .map_err(|source| SpillError::Read {
            field: frame_field,
            source,
        })?;
    if bytes_read != length {
        return Err(SpillError::Truncated { field: frame_field });
    }
    Ok(())
}

/// Read the body's fields from `source`, which holds the inflated body.
fn decode_body<R: BufRead>(source: &mut R) -> Result<SpillEntry, SpillError> {
    let contig = read_u32(source, field::CONTIG)?;
    let position = read_varint(source, field::POSITION)?;
    let is_repeat_tract = read_bool(source, field::IS_REPEAT_TRACT)?;

    let line = read_line(source)?;

    let sample_count = read_usize(source, field::SAMPLE_COUNT)?;
    if sample_count > MAX_SAMPLES {
        return Err(SpillError::OutOfRange {
            field: field::SAMPLE_COUNT,
            value: sample_count as u64,
        });
    }
    let reserve = sample_count.min(MAX_SAMPLES_RESERVED_UP_FRONT);
    let mut samples = Vec::with_capacity(reserve);
    for index in 0..sample_count {
        let sample = decode_sample(source).map_err(|source| SpillError::InSample {
            index,
            source: Box::new(source),
        })?;
        samples.push(sample);
    }

    Ok(SpillEntry {
        contig: ContigId(contig),
        position: Position(position),
        is_repeat_tract,
        line,
        samples,
    })
}

/// The largest fragments frame the decoder will inflate: 3,000 samples × 32 read alleles × at
/// most 16 bytes is 1.5 MiB, so 64 MiB is far above any real entry and still bounds what one
/// corrupt frame can make the reader hold, as [`MAX_LINE_BYTES`] does for the line.
const MAX_FRAGMENTS_BYTES: usize = 64 * 1024 * 1024;

/// The fewest bytes one fragment entry encodes to: four fields of at least one byte each.
const MIN_FRAGMENT_ENTRY_BYTES: usize = 4;

/// **Read the fragments' frame, or step over it**: its length, then — when `fragments` says
/// read — the frame, inflated and decoded, and otherwise that many bytes consumed unread, which
/// is pass two's case. A length of 0 is a record with no fragments.
///
/// A decoded entry is refused when its sample is not one of the entry's `sample_count`, its
/// allele does not fit a `u16`, or its Phred is above [`LARGEST_ERROR_PHRED`]: none of those can
/// be written, so each is damage, and read as data it would change or drop phase silently.
fn read_fragments<R: BufRead>(
    source: &mut R,
    fragments: FragmentsSection,
    sample_count: usize,
    inflating: &mut Inflating,
) -> Result<RecordFragments, SpillError> {
    if fragments == FragmentsSection::Skip {
        let length = read_usize(source, field::FRAGMENTS_BYTES)?;
        if length > max_frame_bytes(MAX_FRAGMENTS_BYTES) {
            return Err(SpillError::OutOfRange {
                field: field::FRAGMENTS_BYTES,
                value: length as u64,
            });
        }
        let skipped = io::copy(&mut source.by_ref().take(length as u64), &mut io::sink()).map_err(
            |source| SpillError::Read {
                field: field::FRAGMENTS,
                source,
            },
        )?;
        if skipped != length as u64 {
            return Err(SpillError::Truncated {
                field: field::FRAGMENTS,
            });
        }
        return Ok(RecordFragments::default());
    }
    read_frame(
        source,
        field::FRAGMENTS_BYTES,
        field::FRAGMENTS,
        max_frame_bytes(MAX_FRAGMENTS_BYTES),
        &mut inflating.compressed,
    )?;
    if inflating.compressed.is_empty() {
        return Ok(RecordFragments::default());
    }
    inflating.inflate(field::FRAGMENTS, MAX_FRAGMENTS_BYTES)?;
    decode_fragments(&inflating.inflated, sample_count)
}

/// Decode the inflated fragments' frame, `section`, of an entry with `sample_count` samples.
fn decode_fragments(section: &[u8], sample_count: usize) -> Result<RecordFragments, SpillError> {
    let mut bytes = section;
    let count = read_usize(&mut bytes, field::FRAGMENT_COUNT)?;
    // Every entry takes at least four bytes, so a count the section cannot hold is refused
    // before anything is reserved for it.
    if count > section.len() / MIN_FRAGMENT_ENTRY_BYTES {
        return Err(SpillError::OutOfRange {
            field: field::FRAGMENT_COUNT,
            value: count as u64,
        });
    }
    let mut entries = Vec::with_capacity(count);
    let mut previous = Steps::default();
    for _ in 0..count {
        let sample_step = read_varint(&mut bytes, field::FRAGMENT_SAMPLE)?;
        let fragment_step = read_varint(&mut bytes, field::FRAGMENT)?;
        let (sample, fragment) = previous.take(sample_step, fragment_step)?;
        if sample as usize >= sample_count {
            return Err(SpillError::OutOfRange {
                field: field::FRAGMENT_SAMPLE,
                value: u64::from(sample),
            });
        }
        let allele = read_varint(&mut bytes, field::FRAGMENT_ALLELE)?;
        let allele = u16::try_from(allele).map_err(|_| SpillError::OutOfRange {
            field: field::FRAGMENT_ALLELE,
            value: allele,
        })?;
        let error_phred =
            read_byte(&mut bytes, field::FRAGMENT_ERROR_PHRED)?.ok_or(SpillError::Truncated {
                field: field::FRAGMENT_ERROR_PHRED,
            })?;
        if error_phred > LARGEST_ERROR_PHRED {
            return Err(SpillError::OutOfRange {
                field: field::FRAGMENT_ERROR_PHRED,
                value: u64::from(error_phred),
            });
        }
        entries.push(FragmentAllele {
            sample,
            fragment,
            allele,
            error_phred,
        });
    }
    if !bytes.is_empty() {
        return Err(SpillError::SectionLongerThanItsEntries {
            field: field::FRAGMENTS,
            unread: bytes.len() as u64,
        });
    }
    Ok(RecordFragments::from_entries(entries))
}

/// Read the record's line: a length, then that many bytes.
///
/// The length is refused before a byte of it is read if it is past [`MAX_LINE_BYTES`], so the
/// bytes held are bounded by that constant and not by what the file still holds.
fn read_line<R: BufRead>(source: &mut R) -> Result<Vec<u8>, SpillError> {
    let line_length = read_usize(source, field::LINE_LENGTH)?;
    if line_length > MAX_LINE_BYTES {
        return Err(SpillError::OutOfRange {
            field: field::LINE_LENGTH,
            value: line_length as u64,
        });
    }
    let mut line = Vec::new();
    let bytes_read = source
        .by_ref()
        .take(line_length as u64)
        .read_to_end(&mut line)
        .map_err(|source| SpillError::Read {
            field: field::LINE,
            source,
        })?;
    if bytes_read != line_length {
        return Err(SpillError::Truncated { field: field::LINE });
    }
    Ok(line)
}

/// Read one sample's four numbers.
fn decode_sample<R: BufRead>(source: &mut R) -> Result<SpilledSample, SpillError> {
    let window = decode_window(source)?;
    let ref_reads = read_u32(source, field::REF_READS)?;
    let alt_reads = read_u32(source, field::ALT_READS)?;
    Ok(SpilledSample {
        window,
        ref_reads,
        alt_reads,
    })
}

/// Read the window pair back from its bit patterns.
fn decode_window<R: BufRead>(source: &mut R) -> Result<WindowCoverage, SpillError> {
    let gc_fraction = f32::from_bits(read_f32_bits(source, field::GC_FRACTION)?);
    let mean_depth = f32::from_bits(read_f32_bits(source, field::MEAN_DEPTH)?);
    Ok(WindowCoverage {
        gc_fraction,
        mean_depth,
    })
}

/// Read one variable-length integer, naming `field` in whatever goes wrong.
fn read_varint<R: BufRead>(source: &mut R, field: &'static str) -> Result<u64, SpillError> {
    let mut bytes = [0u8; MAX_VARINT_BYTES];
    let mut length = 0;
    loop {
        let byte = read_byte(source, field)?.ok_or(SpillError::Truncated { field })?;
        bytes[length] = byte;
        length += 1;
        if byte < 0x80 || length == MAX_VARINT_BYTES {
            break;
        }
    }
    // A ten-byte encoding's last byte contributes `data << 63`, so only its lowest bit fits a
    // `u64` and the primitive drops the rest without complaint. Refuse it here rather than
    // absorb it: everything else this decoder cannot represent is refused, and the primitive
    // is the psp's, shared with the `.psp` reader.
    if length == MAX_VARINT_BYTES && bytes[MAX_VARINT_BYTES - 1] > 0x01 {
        return Err(SpillError::OverlongVarint { field });
    }
    decode_u64_leb128(&bytes[..length])
        .map(|(value, _consumed)| value)
        .map_err(|error| match error {
            VarintError::Overflow => SpillError::OverlongVarint { field },
            VarintError::Truncated => SpillError::Truncated { field },
        })
}

/// Read one variable-length integer that has to fit a `u32`.
fn read_u32<R: BufRead>(source: &mut R, field: &'static str) -> Result<u32, SpillError> {
    let value = read_varint(source, field)?;
    u32::try_from(value).map_err(|_| SpillError::OutOfRange { field, value })
}

/// Read one variable-length integer that has to fit a `usize` — a count or a length.
fn read_usize<R: BufRead>(source: &mut R, field: &'static str) -> Result<usize, SpillError> {
    let value = read_varint(source, field)?;
    usize::try_from(value).map_err(|_| SpillError::OutOfRange { field, value })
}

/// Read one byte written as `0` or `1`.
fn read_bool<R: BufRead>(source: &mut R, field: &'static str) -> Result<bool, SpillError> {
    match read_byte(source, field)?.ok_or(SpillError::Truncated { field })? {
        0 => Ok(false),
        1 => Ok(true),
        byte => Err(SpillError::NotABoolean { field, byte }),
    }
}

/// Read an `f32`'s four little-endian bytes as the bit pattern they are, so that a `NaN`
/// arrives as the same `NaN` rather than as whatever an arithmetic round trip would leave.
fn read_f32_bits<R: BufRead>(source: &mut R, field: &'static str) -> Result<u32, SpillError> {
    let mut bytes = [0u8; 4];
    match source.read_exact(&mut bytes) {
        Ok(()) => Ok(u32::from_le_bytes(bytes)),
        Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => {
            Err(SpillError::Truncated { field })
        }
        Err(source) => Err(SpillError::Read { field, source }),
    }
}

/// The next byte, or `None` at the end of the stream.
fn read_byte<R: BufRead>(source: &mut R, field: &'static str) -> Result<Option<u8>, SpillError> {
    let buffered = source
        .fill_buf()
        .map_err(|source| SpillError::Read { field, source })?;
    let Some(&byte) = buffered.first() else {
        return Ok(None);
    };
    source.consume(1);
    Ok(Some(byte))
}

#[cfg(test)]
mod tests;
