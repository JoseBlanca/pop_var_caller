//! **The spill — the finished VCF lines, parked on disk between the calling pass and the
//! verdict.**
//!
//! A record's *verdict* is what the filter decides about it: kept, dropped, or tagged as
//! better explained by a hidden duplication than by a variant. That decision needs the run's
//! *cut* — the likelihood-ratio threshold the operator's target false-discovery rate resolves
//! to — and the cut is a property of the whole run, computed from the distribution of every
//! record's score. So no record's verdict is known until the last one has been called. Pass
//! one therefore writes each finished record's **line** — the bytes the encoder produced — to a
//! file beside the output, together with the four numbers per sample the scorer reads. Pass
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
//!   is_repeat_tract  u8
//!   is_biallelic_snp u8        -- settled once by whoever fills the entry
//!   line_length      varint
//!   line             bytes     -- the record's line, no newline
//!   sample_count     varint    -- the run's sample count, dense
//!   per sample:
//!     gc_fraction    4 bytes   -- an f32's bits, little-endian; NaN = absent
//!     mean_depth     4 bytes   -- an f32's bits, little-endian
//!     ref_reads      varint    -- AD[0]
//!     alt_reads      varint    -- AD[1] on a biallelic SNP, else 0
//! ```
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
//! **The encoder destructures its input exhaustively**, which is the codebase's way of making
//! a struct that gains a field fail to compile rather than lose it quietly
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
    pub const IS_BIALLELIC_SNP: &str = "is_biallelic_snp";
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
/// count in one allocation and never grows. A run above 8,192 decodes correctly and pays `Vec`
/// growth — this is a reservation policy, not a limit.
const MAX_SAMPLES_RESERVED_UP_FRONT: usize = 8192;

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
/// ([`place_of`](crate::ng::vcf::writer), spec §6 trap 6). `is_biallelic_snp` is settled once by
/// whoever fills the entry, from the record's own alleles *before* the padding base is added —
/// testing the written `REF`/`ALT` would call a one-base deletion a two-base SNP (spec §6 trap
/// 2). Nothing in this module makes that decision; the sink that fills the entry does, and it
/// arrives with the run wiring.
///
/// **A repeat tract is never a biallelic SNP.** Spec §3.2 scores biallelic SNP records on
/// coverage and allele balance and every other record — repeat tracts among them — on coverage
/// alone, so the two flags have three meaningful combinations and not four. The writer and the
/// reader both refuse the fourth.
#[derive(Clone, Debug)]
pub struct SpillEntry {
    /// Which contig the record is written on.
    pub contig: ContigId,
    /// The position the record is *written* at, after the padding rule has moved it.
    pub position: Position,
    /// Whether the record is a repeat tract. Part of the ordering rule: a tract may share a
    /// position with the generic locus that owns its anchor base, and nothing else may.
    pub is_repeat_tract: bool,
    /// Whether the record is a biallelic SNP, which decides whether the scorer sees each
    /// sample's allele counts or only its coverage (spec §3.2).
    pub is_biallelic_snp: bool,
    /// The record's VCF line as the encoder produced it, with no trailing newline.
    pub line: Vec<u8>,
    /// One entry per sample of the run, in the run's sample order — dense, so a sample that
    /// covered nothing still has its place.
    pub per_sample: Vec<SpilledSample>,
}

/// What one sample showed at one record, as the scorer reads it: its window, and the reads
/// behind the two alleles.
#[derive(Clone, Copy, Debug)]
pub struct SpilledSample {
    /// The sample's window at this locus, or the `NaN` pair where it has none.
    pub window: WindowCoverage,
    /// Reads supporting the reference allele — `AD[0]`.
    pub ref_reads: u32,
    /// Reads supporting the alternative — `AD[1]` where the *entry* is a biallelic SNP, `0` on
    /// every other record, because there is no single alternative for the score to read.
    pub alt_reads: u32,
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

    /// A record is marked as both a repeat tract and a biallelic SNP, which spec §3.2 excludes.
    #[error(
        "the record at contig {contig} position {position} is marked as both a repeat tract \
         and a biallelic SNP, and a tract is scored on coverage alone"
    )]
    TractMarkedAsABiallelicSnp {
        /// The record's contig.
        contig: u32,
        /// The record's written position.
        position: u64,
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
    /// If the entry is marked as both a repeat tract and a biallelic SNP, which spec §3.2
    /// excludes; or if the sink refuses the bytes. An entry the sink refused is not counted.
    pub fn append(&mut self, entry: &SpillEntry) -> Result<(), SpillError> {
        if entry.is_repeat_tract && entry.is_biallelic_snp {
            return Err(SpillError::TractMarkedAsABiallelicSnp {
                contig: entry.contig.get(),
                position: entry.position.get(),
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
}

impl<R: BufRead> SpillReader<R> {
    /// Start reading from `source`, which must be positioned at an entry boundary.
    pub fn new(source: R) -> Self {
        Self {
            source,
            failed: false,
        }
    }

    /// The next entry, `None` at a clean end of stream, and `None` for ever after a failure.
    ///
    /// # Errors
    ///
    /// If the stream fails, ends inside an entry, or holds a value no entry can carry. **A
    /// stream that ends *inside* an entry is an error, not an end.**
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
            return None;
        }
        let decoded = decode_entry(&mut self.source);
        if decoded.is_err() {
            self.failed = true;
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
        is_biallelic_snp,
        line,
        per_sample,
    } = entry;

    encode_u64_leb128(u64::from(contig.get()), out);
    encode_u64_leb128(position.get(), out);
    out.push(u8::from(*is_repeat_tract));
    out.push(u8::from(*is_biallelic_snp));
    encode_u64_leb128(line.len() as u64, out);
    out.extend_from_slice(line);
    encode_u64_leb128(per_sample.len() as u64, out);

    for sample in per_sample {
        let SpilledSample {
            window,
            ref_reads,
            alt_reads,
        } = sample;
        let WindowCoverage {
            gc_fraction,
            mean_depth,
        } = window;

        out.extend_from_slice(&gc_fraction.to_bits().to_le_bytes());
        out.extend_from_slice(&mean_depth.to_bits().to_le_bytes());
        encode_u64_leb128(u64::from(*ref_reads), out);
        encode_u64_leb128(u64::from(*alt_reads), out);
    }
}

/// Read one entry from `source`, which is positioned at its first byte.
fn decode_entry<R: BufRead>(source: &mut R) -> Result<SpillEntry, SpillError> {
    let contig = read_u32(source, field::CONTIG)?;
    let position = read_varint(source, field::POSITION)?;
    let is_repeat_tract = read_bool(source, field::IS_REPEAT_TRACT)?;
    let is_biallelic_snp = read_bool(source, field::IS_BIALLELIC_SNP)?;
    if is_repeat_tract && is_biallelic_snp {
        return Err(SpillError::TractMarkedAsABiallelicSnp { contig, position });
    }

    let line = read_line(source)?;

    let sample_count = read_usize(source, field::SAMPLE_COUNT)?;
    let mut per_sample = Vec::with_capacity(sample_count.min(MAX_SAMPLES_RESERVED_UP_FRONT));
    for index in 0..sample_count {
        let sample = decode_sample(source).map_err(|source| SpillError::InSample {
            index,
            source: Box::new(source),
        })?;
        per_sample.push(sample);
    }

    Ok(SpillEntry {
        contig: ContigId(contig),
        position: Position(position),
        is_repeat_tract,
        is_biallelic_snp,
        line,
        per_sample,
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

/// Read one sample's four numbers.
fn decode_sample<R: BufRead>(source: &mut R) -> Result<SpilledSample, SpillError> {
    let gc_fraction = f32::from_bits(read_f32_bits(source, field::GC_FRACTION)?);
    let mean_depth = f32::from_bits(read_f32_bits(source, field::MEAN_DEPTH)?);
    let ref_reads = read_u32(source, field::REF_READS)?;
    let alt_reads = read_u32(source, field::ALT_READS)?;
    Ok(SpilledSample {
        window: WindowCoverage {
            gc_fraction,
            mean_depth,
        },
        ref_reads,
        alt_reads,
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
