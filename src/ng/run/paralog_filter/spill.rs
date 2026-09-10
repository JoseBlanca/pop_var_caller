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
//! that asserts every parity between them ([`VcfRecord::new`](crate::ng::vcf::VcfRecord)); a
//! codec over it is ten field codecs plus a reconstruction that must satisfy every one of those
//! assertions, and a field the record gains later is a field the codec drops in silence. The
//! line is the record's own bytes, so byte-identity for every record the verdict does not touch
//! holds by construction, and nothing about the record's shape has to be mirrored here.
//!
//! **The layout is spec §3.4's, field for field:**
//!
//! ```text
//! entry :=
//!   contig           varint    -- the three fields the writer's ordering check reads
//!   position         varint
//!   is_repeat_tract  u8        -- and it also selects the row shape below
//!   line_length      varint
//!   line             bytes     -- the record's line, no newline
//!   sample_count     varint    -- the run's sample count, dense
//!   per sample, when is_repeat_tract is 0:
//!     gc_fraction    4 bytes   -- an f32's bits, little-endian; NaN = absent
//!     mean_depth     4 bytes   -- an f32's bits, little-endian
//!     ref_reads      varint    -- AD[0]
//!     alt_reads      varint    -- every alternative's reads, summed
//!   per sample, when it is 1:
//!     gc_fraction    4 bytes
//!     mean_depth     4 bytes
//! ```
//!
//! **Why the row has two shapes** (spec §3.2, §3.7). Everything but a repeat tract is scored on
//! coverage *and* allele balance; a tract is scored on coverage alone, because slippage moves
//! reads between its length alleles and the split the scorer's model expects is not the split it
//! would see. Writing that abstention as `0, 0` would put it in the same two fields a real
//! measurement lives in — and production's scorer drops a sample at zero total reads, which would
//! then silently empty every tract's score and leave every tract unfiltered while the file looked
//! correct (spec §6 trap 1). A tract's row carries no read counts at all, so there is no zero to
//! misread. It is also eight bytes a sample instead of ten or more.
//!
//! **One flag does both jobs.** `is_repeat_tract` is in the entry for the writer's ordering rule
//! — a tract may share a position with the generic locus owning its anchor base — and the scoring
//! rule is the same question, so the row shape is read off it rather than off a second flag. One
//! byte less a record, and one fewer pair of fields that can disagree.
//!
//! Every field is either fixed-width, length-prefixed, or self-delimiting — a varint ends at
//! the first byte whose top bit is clear — so an entry ends where the next begins and there is
//! no outer frame. **No manifest, no layout check, no version**: the file is created, read
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
//! ([`cohort_merge`](crate::ng::run::cohort_merge)'s `render`,
//! [`var_calling::types`](crate::var_calling::types)).

use std::io::{self, BufRead, Read, Write};
use std::iter::FusedIterator;

use thiserror::Error;

use super::WindowCoverage;
use crate::ng::types::{ContigId, Position};
use crate::psp::errors::VarintError;
use crate::psp::varint::{MAX_VARINT_BYTES, decode_u64_leb128, encode_u64_leb128};

/// The field names the decoder reports, one per field of the layout above.
///
/// Kept together so that a field renamed on [`SpillEntry`] or [`SpilledSample`] is renamed
/// once, and so no two call sites can spell the same field differently — a failure message
/// naming a field that no longer exists sends the reader somewhere there is nothing to find.
mod field {
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

/// One record on the spill: the line the writer would have written, and what the scorer needs
/// beside it.
///
/// `contig`, `position` and `is_repeat_tract` are the fields the writer's ordering check reads,
/// carried so that pass three can run that same check without rebuilding the record
/// ([`place_of`](crate::ng::vcf::writer), spec §6 trap 6). Which shape [`SpilledSamples`] takes
/// is settled once by whoever fills the entry, from `is_repeat_tract` and nothing else — not from
/// the record's span, and not from its alleles (spec §3.2). The decision is made by
/// [`entry_for`](super::entry_for), in the sink that fills the entry.
///
/// **The flag and the row shape must agree**, since the file stores only the flag and the reader
/// builds the rows from it. An entry whose two halves disagree would be written as one thing and
/// read back as the other, so the writer refuses it.
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
    /// What each sample showed, in the shape the record's span allows.
    pub samples: SpilledSamples,
}

/// **What the samples of one record carry, and it depends on how many bases the record covers.**
///
/// Spec §3.2 draws its line at the locus's span rather than at its alleles: a record occupying
/// **one** genomic position has a reference/non-reference split at that base, whatever its
/// alleles are, and a record occupying **more** does not have one place to take that split at.
///
/// **The wide variant carries no read counts, and that is the whole point of the type.** The
/// filter declining to read a tract's alleles is an *abstention*; a sample that observed nothing
/// is a *measurement*. Written as `0, 0` in the fields a measurement lives in, the two are
/// indistinguishable — and production's scorer drops a sample at zero total reads
/// ([`calibrate.rs:173`](../../../../src/var_calling/paralog_filter/calibrate.rs)), which would
/// then silently empty every wide locus's score and leave every repeat tract unfiltered while
/// the file looked correct (spec §6 trap 1). With two shapes there is no zero to misread.
#[derive(Clone, Debug)]
pub enum SpilledSamples {
    /// Anything that is not a repeat tract — a SNP of any allele count, an insertion, a
    /// deletion. Both signals reach the scorer, and a sample at zero total reads is absent from
    /// the score, as in production.
    GenericLocus(Vec<GenericLocusSample>),
    /// A repeat tract. Coverage alone, and **no sample is skipped**, because there is no read
    /// count that could be zero.
    RepeatTract(Vec<RepeatTractSample>),
}

impl SpilledSamples {
    /// How many samples the record carries, whichever shape it took.
    ///
    /// **The run's sample count, always** — both variants are dense, so a sample that covered
    /// nothing still has its place.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::GenericLocus(samples) => samples.len(),
            Self::RepeatTract(samples) => samples.len(),
        }
    }

    /// Whether the record carries no samples at all, which is a run with no samples.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The window pair for one sample, whichever shape the record took.
    ///
    /// This is the one thing every record has for every sample, so the coverage half of the
    /// score reads it without asking which variant it is looking at.
    #[must_use]
    pub fn window(&self, index: usize) -> Option<WindowCoverage> {
        match self {
            Self::GenericLocus(samples) => samples.get(index).map(|sample| sample.window),
            Self::RepeatTract(samples) => samples.get(index).map(|sample| sample.window),
        }
    }
}

/// What one sample showed at a record covering a single base: its window, and the split between
/// reads carrying the reference base and reads carrying anything else.
#[derive(Clone, Copy, Debug)]
pub struct GenericLocusSample {
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

/// What one sample showed at a record covering more than one base: its window, and nothing else.
///
/// The sample has reads and its call rests on them; what it has no single base to give is a
/// reference/non-reference split. See [`SpilledSamples`] for why that is a missing field rather
/// than a zero.
#[derive(Clone, Copy, Debug)]
pub struct RepeatTractSample {
    /// The sample's window at this locus, or the `NaN` pair where it has none.
    pub window: WindowCoverage,
}

/// What can go wrong reading or writing the spill.
///
/// **Every variant names the field or the operation it failed on**, because the file has no
/// self-description to fall back on: a decoder that says only "truncated" leaves the reader
/// nowhere to look. The path is not here — the writer and the reader are generic over their
/// stream — and is added by the run, which wraps these into a
/// [`RunError`](crate::ng::run::RunError) naming the file (spec §5). The record's ordinal is
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

    /// A field written as `0` or `1` came back as something else.
    #[error("the paralog spill's {field} holds the byte {byte}, and only 0 or 1 are written")]
    NotABoolean {
        /// The field that is neither true nor false.
        field: &'static str,
        /// The byte found there.
        byte: u8,
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

    /// A record's sample shape disagrees with its `is_repeat_tract` flag.
    ///
    /// **The flag is what the file stores**, and the reader builds the row shape from it (spec
    /// §3.4), so an entry whose two halves disagree would be written as one thing and read back
    /// as the other — a tract's samples silently gaining read counts they never measured, or a
    /// generic locus's being dropped. Refused at the writer, where the caller that built them
    /// can still be blamed.
    #[error(
        "the record at contig {contig} position {position} has is_repeat_tract={is_repeat_tract} \
         but the other sample shape, and the flag is what the file stores"
    )]
    SampleShapeDisagreesWithTheTractFlag {
        /// The record's contig.
        contig: u32,
        /// The record's written position.
        position: u64,
        /// What the entry's flag said.
        is_repeat_tract: bool,
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

    /// Append one entry.
    ///
    /// # Errors
    ///
    /// If the entry's sample shape disagrees with its tract flag — the flag is what the file
    /// stores, so a mismatch would be written as the flag and read back as the other shape; or
    /// if the sink refuses the bytes. An entry the sink refused is not counted.
    pub fn append(&mut self, entry: &SpillEntry) -> Result<(), SpillError> {
        // **Matched, not `matches!`**, so a third row shape cannot slip through as "not a
        // tract". The spec already defers a tract-aware allele term that would want one, and a
        // boolean test would have compiled unchanged and quietly written the wrong shape.
        let carries_tract_rows = match entry.samples {
            SpilledSamples::RepeatTract(_) => true,
            SpilledSamples::GenericLocus(_) => false,
        };
        if entry.is_repeat_tract != carries_tract_rows {
            return Err(SpillError::SampleShapeDisagreesWithTheTractFlag {
                contig: entry.contig.get(),
                position: entry.position.get(),
                is_repeat_tract: entry.is_repeat_tract,
            });
        }
        self.scratch.clear();
        encode_entry(entry, &mut self.scratch);
        self.sink
            .write_all(&self.scratch)
            .map_err(|source| SpillError::Write { source })?;
        self.entries_written += 1;
        Ok(())
    }

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
}

impl<R: BufRead> SpillReader<R> {
    /// Start reading from `source`, which must be positioned at an entry boundary, expecting
    /// `entries_expected` of them.
    ///
    /// **The count is required because a spill that lost its tail on an entry boundary is byte
    /// for byte the prefix of a complete one**, so nothing in the bytes can tell the two apart.
    /// It has to come from the side that wrote them, and a reader that could be built without
    /// it would make the shortfall silent again.
    pub fn new(source: R, entries_expected: u64) -> Self {
        Self {
            source,
            failed: false,
            entries_expected,
            entries_read: 0,
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
        if self.failed {
            return None;
        }
        let at_an_entry_boundary_and_out_of_bytes = match self.source.fill_buf() {
            Ok(buffered) => buffered.is_empty(),
            Err(source) => {
                self.failed = true;
                return Some(Err(SpillError::Read {
                    field: field::CONTIG,
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
        let decoded = decode_entry(&mut self.source);
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

/// Append one entry's bytes to `out`, in the layout at the top of this module.
///
/// **The entry and each sample are destructured exhaustively**, so a field added to either
/// struct stops the build here instead of vanishing from the file.
fn encode_entry(entry: &SpillEntry, out: &mut Vec<u8>) {
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

    match samples {
        SpilledSamples::GenericLocus(samples) => {
            for sample in samples {
                let GenericLocusSample {
                    window,
                    ref_reads,
                    alt_reads,
                } = sample;

                encode_window(window, out);
                encode_u64_leb128(u64::from(*ref_reads), out);
                encode_u64_leb128(u64::from(*alt_reads), out);
            }
        }
        SpilledSamples::RepeatTract(samples) => {
            for sample in samples {
                let RepeatTractSample { window } = sample;

                encode_window(window, out);
            }
        }
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

/// Read one entry from `source`, which is positioned at its first byte.
fn decode_entry<R: BufRead>(source: &mut R) -> Result<SpillEntry, SpillError> {
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
    // **The tract flag chooses the row shape** (spec §3.2, §3.4). Reading it rather than a
    // second flag is what makes a mismatch unrepresentable in the file: there is only one thing
    // to be wrong.
    //
    // **⚑ A third row shape needs a wider discriminant in the file, and the compiler will not
    // say so.** One bit distinguishes two shapes and no more, so a third `SpilledSamples`
    // variant — spec §8's deferred tract-aware allele term is the candidate — would be written
    // by an encoder that must then also widen this, while `if is_repeat_tract` goes on
    // compiling and silently decodes it as a generic locus. The writer's own check is a `match`
    // for that reason; this one cannot be, so it is a comment instead.
    let samples = if is_repeat_tract {
        let mut samples = Vec::with_capacity(reserve);
        for index in 0..sample_count {
            let sample =
                decode_repeat_tract_sample(source).map_err(|source| SpillError::InSample {
                    index,
                    source: Box::new(source),
                })?;
            samples.push(sample);
        }
        SpilledSamples::RepeatTract(samples)
    } else {
        let mut samples = Vec::with_capacity(reserve);
        for index in 0..sample_count {
            let sample =
                decode_generic_locus_sample(source).map_err(|source| SpillError::InSample {
                    index,
                    source: Box::new(source),
                })?;
            samples.push(sample);
        }
        SpilledSamples::GenericLocus(samples)
    };

    Ok(SpillEntry {
        contig: ContigId(contig),
        position: Position(position),
        is_repeat_tract,
        line,
        samples,
    })
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

/// Read one sample's four numbers, at a record covering a single base.
fn decode_generic_locus_sample<R: BufRead>(
    source: &mut R,
) -> Result<GenericLocusSample, SpillError> {
    let window = decode_window(source)?;
    let ref_reads = read_u32(source, field::REF_READS)?;
    let alt_reads = read_u32(source, field::ALT_READS)?;
    Ok(GenericLocusSample {
        window,
        ref_reads,
        alt_reads,
    })
}

/// Read one sample's window, at a record covering more than one base — there is nothing else.
fn decode_repeat_tract_sample<R: BufRead>(source: &mut R) -> Result<RepeatTractSample, SpillError> {
    Ok(RepeatTractSample {
        window: decode_window(source)?,
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
