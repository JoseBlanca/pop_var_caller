//! The **psp block**: a run of consecutive records over one span of reference, and the rule
//! that decides where one ends.
//!
//! A block is the unit a reader can start at. Everything that makes that possible is here:
//! the cut is a **grid on the reference coordinate**, so every sample cuts at the same
//! coordinates; a block never crosses a contig, so a reader starting at one knows which
//! chromosome it is on from the block itself; and every block opens with the three facts a
//! reader needs before it can read a record.
//!
//! ```text
//! block payload = contig-id | first-position | record-count | record 0 | record 1 | ...
//!                 └──────────── the block head ───────────┘   └── each with its own head ──┘
//! ```
//!
//! **The block head is to a block what [`RecordHead`] is to a record** — the fixed fields in
//! front, taken before anything else is decided. It sits *inside* the unit that gets
//! compressed rather than in front of it: spec §3.2 asks that a block be self-contained, and a
//! block whose opening facts lived outside its own bytes would not be.
//!
//! **Two halves, and the file is read in this order.** [`BlockBuilder`] turns a stream of
//! records into block payloads, cutting them on the grid; [`BlockCompressor`] turns one payload
//! into a whole on-disk block — a four-byte length and then one zstd frame — with the
//! compressor's look-back window capped at what the file declares, which is what unties the
//! block's size from the reader's memory. **Nothing here writes a file** (Milestone F) and
//! nothing reads one back through a rolling buffer (Milestone D3).
//!
//! Design authority: `doc/devel/ng/spec/psp_file_format.md` §3.2 (a block is self-contained),
//! §4.1 (the cut rule and what it buys), §4.2 (the window is declared and the reader honours
//! it), §8 (the traps), and `doc/devel/ng/arch/psp_file_format.md` §1.
//!
//! [`RecordHead`]: crate::ng::psp::RecordHead

use std::num::NonZeroU64;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::psp::header::{MAX_LOOK_BACK_WINDOW_LOG, MIN_LOOK_BACK_WINDOW_LOG, Manifest};
use crate::ng::psp::record::{RecordEncodeError, RecordEncoder, RecordLayout, RecordLayoutError};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};
use crate::psp::errors::VarintError;
use crate::psp::varint::{decode_u64_leb128, encode_u64_leb128};

// ---------------------------------------------------------------------
// The block head
// ---------------------------------------------------------------------

/// The declared name of a block head's first field, for the messages a refusal carries.
///
/// **These three are not manifest fields.** The manifest declares how a *record's* fields are
/// encoded (spec §4.5); a block's opening three are the container's own framing, fixed by the
/// format version the way the file header's magic and length prefix are. They are named here
/// so a refusal says which of them broke, and so the name a message carries cannot drift from
/// the name this file uses.
///
/// *`record.rs` keeps its field names in a table with a `const fn` indexer, so the name a
/// message carries and the name the manifest declares cannot drift apart. **That purpose does
/// not carry over here**, because there is no manifest entry for these to drift from — so they
/// are three plain constants and the table would buy nothing.*
const BLOCK_CONTIG_ID: &str = "contig-id";
/// See [`BLOCK_CONTIG_ID`].
const BLOCK_FIRST_POSITION: &str = "first-position";
/// See [`BLOCK_CONTIG_ID`].
const BLOCK_RECORD_COUNT: &str = "record-count";

/// What a reader learns about a block before it reads a record: which contig it is on, where
/// its first record starts, and how many records it holds.
///
/// **These three are what make a block self-contained** (spec §3.2), which is the goal a
/// reader starting at an arbitrary block rests on. The contig is here because a block never
/// crosses one. The first position is the base the block's first record's position offset is
/// measured from, and it restarts here rather than continuing from the previous block. The
/// record count is what lets a reader say *the block ended where it should have* rather than
/// *the bytes ran out*.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use]
pub struct BlockHead {
    /// Which contig every record in the block sits on.
    pub contig: ContigId,
    /// Where the block's first record starts — and the base its position offset is measured
    /// from, so that record's offset is zero.
    pub first_position: Position,
    /// How many records the block holds.
    ///
    /// **Never zero, and the type says so.** A block exists because a record went into it, so
    /// a head claiming otherwise is one no writer should be able to build — and while this was
    /// a `u64` guarded by a doc comment, [`encode`](Self::encode) would happily write a head
    /// [`decode`](Self::decode) then refused.
    pub record_count: NonZeroU64,
}

impl BlockHead {
    /// Append the block's opening three fields to `out`.
    ///
    /// It cannot fail: all three are unbounded variable-length integers, and a count of zero
    /// has no representation in the type.
    pub fn encode(&self, out: &mut Vec<u8>) {
        // Destructured with no `..`: **a field added to the head is a compile error here**, at
        // the one place that decides what a file actually carries. The struct literals in
        // `decode` catch a new field on the way in; nothing but this catches it on the way out.
        let Self {
            contig,
            first_position,
            record_count,
        } = self;
        encode_u64_leb128(u64::from(contig.get()), out);
        encode_u64_leb128(first_position.get(), out);
        encode_u64_leb128(record_count.get(), out);
    }

    /// Read a block's opening three fields, and say how many bytes they took.
    ///
    /// **A buffer that stops inside them is [`BlockHeadDecodeError::Truncated`], not damage**,
    /// and the two are different instructions to a streaming reader: the first means
    /// *decompress more of this block and try again*, the second means *the file is damaged*.
    /// At the very start of a block a short buffer is the ordinary state of affairs rather than
    /// an exceptional one, which is why the split matters here as much as it does in a record.
    ///
    /// **A caller that already holds the whole block wants [`BlockRecords::split`] instead**,
    /// which converts the first class into the second: with every byte in hand, *fetch more*
    /// is a retry that never ends.
    pub fn decode(bytes: &[u8]) -> Result<DecodedBlockHead, BlockHeadDecodeError> {
        let mut reader = BlockHeadReader::new(bytes);
        let contig = ContigId(reader.read_u32(BLOCK_CONTIG_ID)?);
        let first_position = Position(reader.read_varint(BLOCK_FIRST_POSITION)?);
        let record_count = reader.read_non_zero(BLOCK_RECORD_COUNT)?;
        Ok(DecodedBlockHead {
            head: Self {
                contig,
                first_position,
                record_count,
            },
            head_bytes: reader.bytes_read(),
        })
    }
}

/// A block's opening three fields, read back, and how many bytes they took.
///
/// **A named pair rather than a tuple**, for the reason `record.rs` gives at
/// [`DecodedRecordBody`]: `let (head, _) = …` is shorter than using the second value, and the
/// second value is where the block's first record begins.
///
/// [`DecodedRecordBody`]: crate::ng::psp::DecodedRecordBody
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use]
pub struct DecodedBlockHead {
    /// What the block's opening fields said.
    pub head: BlockHead,
    /// How many bytes the head occupied — where the block's first record begins.
    pub head_bytes: usize,
}

/// A cursor over a block head's bytes that reports running out instead of indexing past the
/// end.
///
/// **The same shape as `record.rs`'s `FieldReader`, and deliberately a second one.** That type
/// is hard-wired to `RecordDecodeError`; making it generic over the fault would put a trait
/// call or a callback on the path that decodes about twenty million records a second, to serve
/// a parser that runs once per block — roughly a hundred and sixty times per sample against
/// several million. What the duplication has to carry over is the *invariant*, which is
/// `bytes_read <= bytes.len()`: it is what makes the slicing below total, and it is why every
/// method that advances checks the bytes are there first.
struct BlockHeadReader<'a> {
    bytes: &'a [u8],
    bytes_read: usize,
}

impl<'a> BlockHeadReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            bytes_read: 0,
        }
    }

    fn bytes_read(&self) -> usize {
        self.bytes_read
    }

    fn truncated(&self, field: &'static str) -> BlockHeadDecodeError {
        BlockHeadDecodeError::Truncated {
            field,
            bytes_in: self.bytes_read,
        }
    }

    fn malformed(&self, field: &'static str, reason: String) -> BlockHeadDecodeError {
        BlockHeadDecodeError::Malformed {
            field,
            bytes_in: self.bytes_read,
            reason,
        }
    }

    /// One variable-length integer, read through production's codec.
    ///
    /// **[`VarintError`] is matched exhaustively and not through a wildcard**, because this is
    /// the one place the two reader instructions are told apart: a variant added to it later
    /// has to be a compile error here rather than land in the damage class by default. The
    /// overflow's wording is this format's rather than the codec's, for the same reason
    /// `record.rs` rewords it — "the 10-byte cap" is the codec's vocabulary, not a psp's.
    fn read_varint(&mut self, field: &'static str) -> Result<u64, BlockHeadDecodeError> {
        match decode_u64_leb128(&self.bytes[self.bytes_read..]) {
            Ok((value, used)) => {
                self.bytes_read += used;
                Ok(value)
            }
            Err(VarintError::Truncated) => Err(self.truncated(field)),
            Err(VarintError::Overflow) => Err(self.malformed(
                field,
                "a variable-length integer longer than any 64-bit value needs".to_string(),
            )),
        }
    }

    fn read_u32(&mut self, field: &'static str) -> Result<u32, BlockHeadDecodeError> {
        let value = self.read_varint(field)?;
        u32::try_from(value).map_err(|_| {
            self.malformed(
                field,
                format!("{value}, which is past the {} this field holds", u32::MAX),
            )
        })
    }

    fn read_non_zero(&mut self, field: &'static str) -> Result<NonZeroU64, BlockHeadDecodeError> {
        let value = self.read_varint(field)?;
        NonZeroU64::new(value).ok_or_else(|| {
            self.malformed(
                field,
                "a block holding no records; a block exists because a record went into it"
                    .to_string(),
            )
        })
    }
}

/// Why a block's opening fields could not be read.
///
/// **Two classes, and a streaming reader branches on them** — the same split
/// [`RecordDecodeError`] makes for a record, and for the same reason: one says *fetch more
/// bytes and retry*, the other says *this file is damaged*, and a fault put in the wrong class
/// makes a reader either reject a good block or retry for ever on a bad one.
///
/// **There is no third class, and the rule that licenses that is a versioning one.** A record
/// can say *upgrade the reader* because the manifest describes its fields, so a reader meets an
/// unfamiliar one knowing it is unfamiliar. A block head is container framing that no manifest
/// describes, so a reader has no way to learn a fourth opening field exists — which means
/// **the block head can only change in a major version**, and a file whose major version this
/// reader does not know is refused before any block is read (spec §3.1). Widening the head in a
/// minor version would need a third class here first.
///
/// [`RecordDecodeError`]: crate::ng::psp::RecordDecodeError
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BlockHeadDecodeError {
    /// The bytes ran out while this field was being read, and what it declared was possible.
    #[error("the block's {field} runs past the bytes it was given, {bytes_in} bytes in")]
    Truncated {
        field: &'static str,
        /// How far into the block head the reader had got.
        bytes_in: usize,
    },
    /// The bytes were there and cannot mean what they say.
    #[error("the block's {field}, {bytes_in} bytes in, is unreadable: {reason}")]
    Malformed {
        field: &'static str,
        bytes_in: usize,
        reason: String,
    },
}

/// A whole block payload split into its head and the records behind it.
///
/// **For a caller holding a block entire** — a writer's own tests, and the compressed-block
/// round trip Milestone D2 adds. The streaming reader of Milestone D3 holds a rolling buffer
/// rather than a whole block and goes to [`BlockHead::decode`] directly.
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub struct BlockRecords<'a> {
    /// What the block's own head said.
    pub head: BlockHead,
    /// The records, still encoded, from the first record's first byte.
    ///
    /// **Not checked against [`head`](Self::head)'s record count.** Splitting reads the head
    /// and bounds nothing behind it; whether the block holds the number of records it declares
    /// is only known once they have been walked, so it is the walk that owes that check.
    pub records: &'a [u8],
}

impl<'a> BlockRecords<'a> {
    /// Split a whole block payload into its head and the records behind it.
    ///
    /// **A short head is damage here, never `Truncated`.** The caller holds the block entire,
    /// so no quantity of further bytes changes the answer, and a reader handed the retry class
    /// about a block already wholly in memory retries for ever. That is the C2 review's Blocker
    /// one level up: `record.rs` converts the same class at the same kind of boundary, once the
    /// container's length is known (`RecordDecodeError::inside_a_bounded_body`).
    pub fn split(payload: &'a [u8]) -> Result<Self, BlockHeadDecodeError> {
        let decoded = BlockHead::decode(payload).map_err(|fault| match fault {
            BlockHeadDecodeError::Truncated { field, bytes_in } => {
                BlockHeadDecodeError::Malformed {
                    field,
                    bytes_in,
                    reason: format!(
                        "it runs past the {} bytes this whole block holds",
                        payload.len()
                    ),
                }
            }
            damage => damage,
        })?;
        Ok(Self {
            head: decoded.head,
            records: &payload[decoded.head_bytes..],
        })
    }
}

// ---------------------------------------------------------------------
// The cut
// ---------------------------------------------------------------------

/// Which block a coordinate belongs to: its contig, and which cell of the coordinate grid it
/// falls in.
///
/// **The grid is the point of the cut rule.** A block ends when a position crosses into the
/// next multiple of the genomic block size, which is not the same thing as a block ending once
/// it has covered that many bases: a grid makes every sample cut at the *same* coordinates, so
/// a cohort reader stepping across a region touches one aligned block per sample rather than
/// one in some samples and two in others (spec §4.1). A running count would not align, and
/// losing that would have been an accident rather than a decision.
///
/// *Named for the grid and not for the block, because the cell **contains** the block: a type
/// called `BlockCell` reads as one more level down — blocks contain records, records contain
/// fields — which is the containment backwards.*
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct GridCell {
    contig: ContigId,
    /// Which multiple of the genomic block size the coordinate falls in — an index along the
    /// grid, not a coordinate and not a byte count.
    cell_index: u64,
}

/// The block being built: what its head will say, once it is known how many records went in.
#[derive(Debug, Clone, Copy)]
struct OpenBlock {
    grid_cell: GridCell,
    first_position: Position,
    record_count: NonZeroU64,
}

impl OpenBlock {
    fn head(self) -> BlockHead {
        BlockHead {
            contig: self.grid_cell.contig,
            first_position: self.first_position,
            record_count: self.record_count,
        }
    }
}

/// Turns a stream of records into block payloads, cutting them on the coordinate grid.
///
/// **It hands back bytes, not a file.** One payload is exactly what a compressor is given
/// (Milestone D2) and what the block index will point at (Milestone F1); nothing here
/// compresses or writes.
///
/// # What closes a block
///
/// - **The grid**, always: a record whose start falls in a different multiple of
///   `genomic_block_size_bp` from the open block's opens a new one.
/// - **A contig change**, always: a block never crosses a contig (spec §3.2).
/// - **The byte ceiling**, when the writer declared one: a block that has already reached it is
///   closed before the next record goes in. It is a secondary rule because a span-cut block has
///   a variable size in bytes — at three hundred reads a position a fully covered 100 kb span
///   is a large thing to hold while writing (spec §4.1; spec §12 question 2 is still open on
///   what the value should be, which is why the default is no ceiling at all).
///
/// # What it refuses
///
/// **Records must arrive in coordinate order, and the cut is where that would otherwise be
/// laundered.** [`RecordEncoder`] refuses a record starting before the base its offset is
/// measured from — but that base *resets* at every cut, so a record going backwards far enough
/// to land in an earlier grid cell would open a new block and be written without complaint, at
/// a coordinate the index would then seek to wrongly. This type checks the order before it
/// decides anything, which is why the check lives here and not in the encoder.
///
/// **A refused record leaves the builder exactly as it was.** Nothing is half-written and no
/// block half-closed, so a caller that reports a refusal and carries on produces the file it
/// would have produced had the refused record never been offered — measured over 40,000 random
/// record streams with 91,247 refusals interleaved, byte-identical to replaying only the
/// records that were accepted (the D1 review's fuzzing).
#[derive(Debug)]
#[must_use]
pub struct BlockBuilder {
    // ---- what lives for the whole file ----
    genomic_block_size_bp: Bp,
    block_byte_ceiling: Option<u32>,
    encoder: RecordEncoder,
    /// Where the next record is laid down while a cut is still undecided. Swapped with
    /// `records` when the cut commits, so both buffers stay warm and neither is reallocated
    /// per block.
    next_records: Vec<u8>,
    /// The block just closed: its head, then its records. Handed out by reference and
    /// overwritten by the next close.
    closed_block_payload: Vec<u8>,
    /// Where the last accepted record sits, so the order check has something to compare
    /// against across a cut — where the encoder's own base has been reset and cannot see it.
    last_accepted_region: Option<GenomeRegion>,

    // ---- what lives for one block, and is replaced whole at every cut ----
    /// The records of the block being built, each already behind its own head.
    records: Vec<u8>,
    open_block: Option<OpenBlock>,
}

impl BlockBuilder {
    /// A builder cutting on `genomic_block_size_bp`, optionally closing a block early once it
    /// holds `block_byte_ceiling` bytes of records.
    ///
    /// **Both zeroes are refused rather than acted on**, and the argument is one `header.rs`
    /// already makes for the grid: `Manifest`'s fields are public, so a caller can build one
    /// that never met the header's validation. A zero grid has no cells to divide by; a zero
    /// ceiling is reached by every block before its second record, so it gives every position
    /// its own block — one index entry and one compressed frame each.
    pub fn new(
        genomic_block_size_bp: Bp,
        block_byte_ceiling: Option<u32>,
    ) -> Result<Self, BlockCutRuleError> {
        if genomic_block_size_bp.get() == 0 {
            return Err(BlockCutRuleError::ZeroGenomicBlockSize);
        }
        if block_byte_ceiling == Some(0) {
            return Err(BlockCutRuleError::ZeroBlockByteCeiling);
        }
        Ok(Self {
            genomic_block_size_bp,
            block_byte_ceiling,
            // Replaced before the first record is written: a builder that has been handed none
            // has no block, and so no base to measure one from.
            encoder: RecordEncoder::for_block(Position(0)),
            next_records: Vec::new(),
            closed_block_payload: Vec::new(),
            last_accepted_region: None,
            records: Vec::new(),
            open_block: None,
        })
    }

    /// A builder driven by a file's declared cut rule, which is the only rule a writer
    /// extending that file may use — the manifest is fixed when the file is created and an
    /// append does not rewrite it (spec §6.4).
    ///
    /// **The declared field layout is checked here too, and it has to be**: an append writes
    /// records with whatever layout this build encodes, into a file that already declares one.
    /// A manifest naming different fields, or the same fields differently encoded, is a file
    /// this writer cannot extend — and left unchecked it would be extended anyway, with the
    /// added records unreadable under the header the file keeps.
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, BlockCutRuleError> {
        // Destructured with no `..`, and every field named: a field added to the manifest is
        // then a compile error here rather than a setting a writer silently fails to honour.
        let Manifest {
            genomic_block_size_bp,
            block_byte_ceiling,
            // Milestone D2's: it configures the compressor, which this builder does not own.
            look_back_window_log: _,
            // Checked below rather than ignored — `RecordLayout::from_manifest` is the one
            // reader of this field, and it refuses a layout this build cannot write.
            fields: _,
        } = manifest;
        RecordLayout::from_manifest(manifest)
            .map_err(|source| BlockCutRuleError::UnsupportedRecordLayout { source })?;
        Self::new(*genomic_block_size_bp, *block_byte_ceiling)
    }

    /// Lay `record` down, and hand back the block it closed if it closed one.
    ///
    /// The returned bytes are one whole block payload — its head, then its records — and they
    /// stay valid until the next call. `None` means the record joined the block being built and
    /// nothing is ready yet; [`finish`](Self::finish) closes the last one.
    pub fn push(
        &mut self,
        record: &SampleLocusObservations,
    ) -> Result<Option<&[u8]>, BlockWriteError> {
        self.check_order(record.region)?;
        let grid_cell = self.grid_cell_of(record.region);

        let Some(open) = self.open_block else {
            // The file's first record. The encoder is positioned before anything is written,
            // so a record the codec refuses leaves no block open behind it.
            self.encoder.start_block(record.region.start);
            if let Err(refused) = self.encoder.encode_record(record, &mut self.records) {
                self.records.clear();
                return Err(refused.into());
            }
            self.open_block_at(grid_cell, record.region);
            return Ok(None);
        };

        if grid_cell == open.grid_cell && !self.has_reached_byte_ceiling() {
            // Straight into the live block's buffer, so the rollback is here: the encoder
            // promises `out` is untouched when it refuses, and truncating makes this block's
            // bytes not depend on that promise holding.
            let before = self.records.len();
            if let Err(refused) = self.encoder.encode_record(record, &mut self.records) {
                self.records.truncate(before);
                return Err(refused.into());
            }
            self.count_one_more_record(record.region);
            return Ok(None);
        }

        // A cut. The record goes down first, in the buffer that will become the new block's,
        // so that a record the codec refuses leaves the open block open and loses nothing.
        self.next_records.clear();
        let resume_at = self.encoder.measured_from().position();
        self.encoder.start_block(record.region.start);
        if let Err(refused) = self.encoder.encode_record(record, &mut self.next_records) {
            self.encoder.start_block(resume_at);
            self.next_records.clear();
            return Err(refused.into());
        }

        self.closed_block_payload.clear();
        open.head().encode(&mut self.closed_block_payload);
        self.closed_block_payload.extend_from_slice(&self.records);

        std::mem::swap(&mut self.records, &mut self.next_records);
        self.open_block_at(grid_cell, record.region);
        Ok(Some(&self.closed_block_payload))
    }

    /// Close the block being built, if there is one. `None` when every record pushed has
    /// already been handed back inside a block, which is also what a builder that was never
    /// pushed to returns.
    ///
    /// **It consumes the builder**, which is the one thing that cannot be got wrong afterwards:
    /// a builder that could be pushed to after closing would put the closed block's records
    /// into the next one, and a builder that could be closed twice would put the last block in
    /// the file twice. Both were reachable while it took `&mut self`, and neither was caught by
    /// a test. The type is what closes them, not a test.
    ///
    /// It returns owned bytes where [`push`](Self::push) lends them — once per file, against a
    /// borrow per block.
    pub fn finish(mut self) -> Option<Vec<u8>> {
        let open = self.open_block.take()?;
        self.closed_block_payload.clear();
        open.head().encode(&mut self.closed_block_payload);
        self.closed_block_payload.extend_from_slice(&self.records);
        Some(self.closed_block_payload)
    }

    /// Which block a record's start belongs to. **Its start and not its end**: a record widened
    /// by a deletion may reach past its own block, and it must still belong to the cell its
    /// start falls in — a span is sample-dependent, so a cut taken from the end would put the
    /// same reference span in different blocks in different samples, which is exactly what the
    /// grid exists to prevent. A reader learns each record's span from its own head.
    fn grid_cell_of(&self, region: GenomeRegion) -> GridCell {
        GridCell {
            contig: region.contig,
            cell_index: region.start.get() / self.genomic_block_size_bp.get(),
        }
    }

    /// Whether the open block has already reached the declared byte ceiling.
    ///
    /// **It measures the records laid down, not the block head in front of them**, and it is
    /// checked before the next record rather than after the last — so a block may pass the
    /// ceiling by one record. That is what the rule costs: the alternative decides a record's
    /// fate from a length nothing knows until the record has been encoded.
    fn has_reached_byte_ceiling(&self) -> bool {
        self.block_byte_ceiling
            .is_some_and(|ceiling| self.records.len() >= ceiling as usize)
    }

    /// Open a block on `region`'s grid cell, with the record that opened it already counted.
    fn open_block_at(&mut self, grid_cell: GridCell, region: GenomeRegion) {
        self.open_block = Some(OpenBlock {
            grid_cell,
            first_position: region.start,
            record_count: NonZeroU64::MIN,
        });
        self.last_accepted_region = Some(region);
    }

    /// Count one more record into the block that is open.
    ///
    /// **The block must be open, and saying so is the point**: the version that absorbed a
    /// missing one with `if let` would have left a head declaring one record fewer than the
    /// block holds, which a reader meets as the bytes running out.
    fn count_one_more_record(&mut self, region: GenomeRegion) {
        let open = self
            .open_block
            .as_mut()
            .expect("a record is only counted into a block that is already open");
        // `NonZeroU64` has no `+=`; a block reaching 2^64 records is not a case worth an error
        // path, and saturating there costs nothing that could ever be reached.
        open.record_count = open.record_count.saturating_add(1);
        self.last_accepted_region = Some(region);
    }

    /// Refuse a record that does not follow the one before it along the reference.
    ///
    /// Contigs are visited in ascending order and never revisited; within one contig a record
    /// must not start before the record before it. **Two records starting on the same base are
    /// allowed** — a repeat tract and a generic locus can begin together — because that is not
    /// out of order.
    fn check_order(&self, offered: GenomeRegion) -> Result<(), BlockWriteError> {
        let Some(previous) = self.last_accepted_region else {
            return Ok(());
        };
        if offered.contig != previous.contig {
            return if offered.contig > previous.contig {
                Ok(())
            } else {
                Err(BlockWriteError::ContigOutOfOrder {
                    previous: previous.contig,
                    offered: offered.contig,
                })
            };
        }
        if offered.start < previous.start {
            return Err(BlockWriteError::OutOfOrder { previous, offered });
        }
        Ok(())
    }
}

/// Why a builder could not be made: the cut rule it was handed cannot cut.
///
/// **Its own type rather than a variant of [`BlockWriteError`]**, because the two operations
/// fail in disjoint ways: a builder that exists can never raise these, and building one can
/// never raise the others. Milestone F3 maps each to a `PspWriteError` that also names the
/// file, and a variant that cannot occur at a call site is how that mapping acquires an
/// `unreachable!`.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BlockCutRuleError {
    /// The grid has no cells to cut on.
    #[error(
        "the genomic block size is zero; the cut is a grid on the coordinate and a zero grid \
         has no cells"
    )]
    ZeroGenomicBlockSize,
    /// The byte ceiling is zero, which every block reaches before its second record.
    #[error(
        "the block byte ceiling is zero; a ceiling no block can stay under gives every record \
         a block of its own"
    )]
    ZeroBlockByteCeiling,

    /// The file declares a record layout this build does not write, so its blocks cannot be
    /// extended with records this build produces.
    #[error("the file declares a record layout this writer cannot honour: {source}")]
    UnsupportedRecordLayout {
        #[source]
        source: RecordLayoutError,
    },
}

/// Why a record could not be laid down in a block.
///
/// **Every variant is a record the writer was handed that the format cannot hold**, not an
/// internal fault. The writer of Milestone F3 turns each into a [`PspWriteError`] that also
/// names the file, which is the one thing a cohort gathering sixty samples at once needs and
/// this type does not know — including a `ContigOutOfOrder` of its own, since
/// `PspWriteError::OutOfOrder` renders a sentence about positions that is wrong for a contig
/// revisited.
///
/// [`PspWriteError`]: crate::ng::psp::PspWriteError
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BlockWriteError {
    /// A record the record codec cannot lay down: an empty region, a region at the coordinate
    /// ceiling, or a body longer than a head can describe.
    #[error(transparent)]
    RecordRefused(#[from] RecordEncodeError),

    /// A record that starts before the record before it, on the same contig.
    ///
    /// **Refused here rather than left to the record encoder**, whose own check cannot see it:
    /// the base a record's offset is measured from resets at every cut, so a record going back
    /// far enough to land in an earlier grid cell would open a new block and be written without
    /// complaint.
    #[error("a record at {offered} starts before the previous record at {previous}")]
    OutOfOrder {
        previous: GenomeRegion,
        offered: GenomeRegion,
    },

    /// A record on a contig the writer has already finished with, or on one that comes before
    /// the contig it is writing. Blocks are indexed in genomic order, and a contig visited
    /// twice gives two runs of blocks a seek cannot choose between.
    #[error(
        "a record on contig {offered_contig} after a record on contig {previous_contig}; \
         contigs are written in ascending order and never revisited",
        offered_contig = offered.get(),
        previous_contig = previous.get()
    )]
    ContigOutOfOrder {
        previous: ContigId,
        offered: ContigId,
    },
}

// ---------------------------------------------------------------------
// Compressing a block
// ---------------------------------------------------------------------

/// The level every psp block is compressed at.
///
/// **Production's number** (`psp::block::ZSTD_COMPRESSION_LEVEL`), kept because every byte
/// figure the specs quote was measured at it and because nothing here has re-measured the
/// curve.
///
/// **It is not declared in the manifest, and does not need to be**: a zstd frame carries
/// everything a decoder needs, so no reader is ever driven by the level — which is the
/// difference between it and the genomic block size or the look-back window (spec §4.5). But
/// *nobody is driven by it* is not the same as *nobody needs to know it*: goal 5's
/// byte-identity check, an append that must match the bytes already in a file, and every byte
/// figure the specs quote all rest on the level a file was written at. So a compressor keeps
/// its own ([`BlockCompressor::compression_level`]), and **recording it in the header's
/// writer-provenance parameters is Milestone F3's** — that is where a writer's exposed knobs
/// go, and nothing writes a header yet.
pub const ZSTD_COMPRESSION_LEVEL: i32 = 9;

/// How many bytes stand in front of each compressed block, saying how long it is.
///
/// **A byte count and nothing else.** The measuring prototype's framing carries the
/// *uncompressed* length beside it; this does not, and spec §8 says why — "a block header that
/// carries its uncompressed length is a temptation to allocate it", and the whole design is
/// that a reader never sizes a buffer from a block. What the count is for is the opposite: it
/// bounds how many compressed bytes belong to *this* block, so a reader feeding a decoder can
/// stop at the frame's end without trusting the frame to tell it.
pub const COMPRESSED_BLOCK_LENGTH_BYTES: usize = 4;

/// Compresses one block payload at a time, with the look-back window capped at what the file
/// declares.
///
/// **The cap is the whole design in one parameter.** Without it zstd sizes its window from the
/// data it is given, so a reader would have to hold a whole block's worth of history to resolve
/// a back-reference — which ties the block size to the reader's memory, and untying those two
/// is what this format exists for (spec §1). With it, a block can be as large as compression
/// and a small index want, while what a reader must hold is the declared window and nothing
/// more.
///
/// **One compressor for the life of a writer, not one per block.** zstd allocates its working
/// space and tables when the context is made; production found the per-call setup of a fresh
/// encoder to be the dominant allocator traffic in its own writer, and kept one alive for the
/// same reason.
///
/// **What it hands back is a whole on-disk block**: a four-byte length and then one zstd frame.
pub struct BlockCompressor {
    zstd: zstd::bulk::Compressor<'static>,
    look_back_window_log: u8,
    compression_level: i32,
    /// The zstd frame alone, and it has to be its own buffer.
    ///
    /// **`compress_to_buffer` writes from the buffer's *start*, not from its end** — its
    /// `WriteBuf` for `Vec<u8>` hands zstd `as_mut_ptr()` and then `set_len`. So a length
    /// prefix cannot be reserved in front and filled in afterwards: zstd overwrites it, and
    /// what a reader meets is a four-byte length where the frame's magic should be. Measured,
    /// on the first run of this file's round-trip test.
    zstd_frame_scratch: Vec<u8>,
    /// The length and the frame together — what a caller is handed.
    on_disk_block: Vec<u8>,
}

impl std::fmt::Debug for BlockCompressor {
    /// `zstd::bulk::Compressor` is not `Debug`, so this is written by hand — and it
    /// destructures `Self` with no `..`, so a setting added to the compressor is a compile
    /// error here rather than one that silently stops being reportable.
    ///
    /// **Buffer lengths, never buffer contents**: what they hold is a whole block, and a
    /// `{:?}` of a compressor mid-file would otherwise run to megabytes.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            zstd: _,
            look_back_window_log,
            compression_level,
            zstd_frame_scratch,
            on_disk_block,
        } = self;
        f.debug_struct("BlockCompressor")
            .field("look_back_window_log", look_back_window_log)
            .field("compression_level", compression_level)
            .field("zstd_frame_scratch_bytes", &zstd_frame_scratch.len())
            .field("on_disk_block_bytes", &on_disk_block.len())
            .finish_non_exhaustive()
    }
}

impl BlockCompressor {
    /// A compressor for the file `manifest` describes, at the level every psp is written at.
    ///
    /// **This is the constructor a writer uses, and the reason is that the window has to come
    /// from the same place the reader will read it from.** A compressor handed a window by hand
    /// can be given one the file does not declare, and a file whose frames need more window
    /// than its own manifest promises is one every reader in the cohort refuses, block by
    /// block, with a zstd error code — precisely the unactionable failure spec §4.2 says the
    /// declaration exists to prevent.
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, BlockCompressError> {
        // Destructured with no `..`, and every field named: a field added to the manifest that
        // a compressor would have to honour is a compile error here. Before this existed,
        // adding one flagged three sites in the crate and none of them was the compressor.
        let Manifest {
            look_back_window_log,
            // `BlockBuilder`'s, both of them: they decide where a block ends, not how it is
            // compressed.
            genomic_block_size_bp: _,
            block_byte_ceiling: _,
            // The record codec's, honoured by `RecordEncoder`.
            fields: _,
        } = manifest;
        Self::with_level(*look_back_window_log, ZSTD_COMPRESSION_LEVEL)
    }

    /// A compressor for a file declaring `look_back_window_log`, at the level every psp is
    /// written at.
    ///
    /// **[`from_manifest`](Self::from_manifest) is what a writer should use.** This is for a
    /// caller that has a window and no manifest — the tests below, and Milestone D3's own
    /// fixtures.
    pub fn new(look_back_window_log: u8) -> Result<Self, BlockCompressError> {
        Self::with_level(look_back_window_log, ZSTD_COMPRESSION_LEVEL)
    }

    /// The same, at a level of the caller's choosing. **For measuring what the level costs**,
    /// which nothing in this project has done on ng's own records; a writer uses
    /// [`from_manifest`](Self::from_manifest).
    pub fn with_level(look_back_window_log: u8, level: i32) -> Result<Self, BlockCompressError> {
        if !(MIN_LOOK_BACK_WINDOW_LOG..=MAX_LOOK_BACK_WINDOW_LOG).contains(&look_back_window_log) {
            return Err(BlockCompressError::WindowLogOutOfRange {
                look_back_window_log,
            });
        }
        // **Checked rather than left to zstd**, which clamps silently: measured before this
        // existed, a level of 100 produced exactly the level-22 bytes and a level of −131,073
        // produced a block 169 times the shipped level's. Either is a file written at a level
        // nobody asked for and nothing records.
        let allowed = zstd::compression_level_range();
        if !allowed.contains(&level) {
            return Err(BlockCompressError::LevelOutOfRange {
                level,
                lowest: *allowed.start(),
                highest: *allowed.end(),
            });
        }
        let mut zstd = zstd::bulk::Compressor::new(level)
            .map_err(|source| BlockCompressError::zstd("creating the compressor", source))?;
        // The cap that makes the block size and the reader's memory independent.
        zstd.set_parameter(zstd::zstd_safe::CParameter::WindowLog(u32::from(
            look_back_window_log,
        )))
        .map_err(|source| BlockCompressError::zstd("capping the look-back window", source))?;
        // **Off deliberately.** With it on, zstd writes the payload's length into the frame
        // header, and a reader that meets it is one `reserve` away from allocating a whole
        // block — the temptation spec §8 names. It also costs bytes we do not need.
        zstd.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
            .map_err(|source| BlockCompressError::zstd("clearing the content size", source))?;
        // A frame checksum, as production's writer sets: it is what turns a damaged block into
        // a refusal rather than into records that decode to plausible values.
        zstd.include_checksum(true)
            .map_err(|source| BlockCompressError::zstd("enabling the frame checksum", source))?;
        Ok(Self {
            zstd,
            look_back_window_log,
            compression_level: level,
            zstd_frame_scratch: Vec::new(),
            on_disk_block: Vec::new(),
        })
    }

    /// The look-back window this compressor caps at, as the exponent the file declares.
    pub fn look_back_window_log(&self) -> u8 {
        self.look_back_window_log
    }

    /// The level this compressor writes at. **Nothing in a psp records it** — see
    /// [`ZSTD_COMPRESSION_LEVEL`] for why, and for whose job recording it is.
    pub fn compression_level(&self) -> i32 {
        self.compression_level
    }

    /// Compress one block payload into a whole on-disk block: its length, then its frame.
    ///
    /// The bytes stay valid until the next call, the way [`BlockBuilder::push`] lends the
    /// payload that comes in here.
    ///
    /// **The same payload always gives the same bytes**, whether this compressor has seen a
    /// thousand blocks before it or none — which is what lets goal 5's byte-identity check across
    /// worker counts mean anything.
    ///
    /// **A frame can be longer than the payload it came from**, which is why the reservation is
    /// zstd's own bound and not the payload's length: a block holding one small record is
    /// incompressible, and zstd's framing and checksum then cost more than they save.
    pub fn compress(&mut self, payload: &[u8]) -> Result<&[u8], BlockCompressError> {
        self.zstd_frame_scratch.clear();
        // `compress_to_buffer` writes into the spare capacity from the buffer's start, so the
        // worst case has to be there before it is called; the reservation is kept between
        // blocks, so only a larger block than any before it allocates.
        self.zstd_frame_scratch
            .reserve(zstd::zstd_safe::compress_bound(payload.len()));
        self.zstd
            .compress_to_buffer(payload, &mut self.zstd_frame_scratch)
            .map_err(|source| BlockCompressError::zstd("compressing a block", source))?;

        let frame_bytes = self.zstd_frame_scratch.len();
        let declared =
            u32::try_from(frame_bytes).map_err(|_| BlockCompressError::FrameTooLong {
                frame_bytes,
                payload_bytes: payload.len(),
            })?;
        self.on_disk_block.clear();
        self.on_disk_block
            .reserve(COMPRESSED_BLOCK_LENGTH_BYTES + frame_bytes);
        self.on_disk_block
            .extend_from_slice(&declared.to_le_bytes());
        self.on_disk_block
            .extend_from_slice(&self.zstd_frame_scratch);
        Ok(&self.on_disk_block)
    }
}

/// What is at the front of a run of bytes a reader has pulled from a file.
///
/// **Three states and not a length**, because the distinction is the one the whole module
/// types: a reader that has some of a block does something different from a reader that has all
/// of it, and a bare `Option<usize>` says neither. Before this existed the function handed back
/// a length taken from four untrusted bytes with nothing bounding it — measured under fuzzing,
/// it returned a length longer than the slice it was given 656,785 times, the largest
/// 4,294,967,299 from eight bytes in — and the module's own walk sliced on it directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use]
pub enum CompressedBlockAt<'a> {
    /// The whole block is here: its zstd frame, and how many bytes to advance to reach the
    /// next block.
    Whole {
        zstd_frame: &'a [u8],
        block_bytes: usize,
    },
    /// The length is here and the block is not all of it. **`block_bytes` is what the four
    /// bytes *declared*, not a fact** — a reader pulls that many and finds out.
    NotAllHere { block_bytes: usize },
    /// Fewer than [`COMPRESSED_BLOCK_LENGTH_BYTES`] bytes: not even a length yet.
    NoLengthYet,
}

/// What is at the front of `bytes`, at a block boundary.
///
/// **Nothing is allocated and no frame is touched** — the four-byte length is read and that is
/// all. It is what a reader advances by to reach the next block, and what tells a decoder where
/// this block's bytes stop.
pub fn compressed_block_at(bytes: &[u8]) -> CompressedBlockAt<'_> {
    let Some(declared) = bytes.get(..COMPRESSED_BLOCK_LENGTH_BYTES) else {
        return CompressedBlockAt::NoLengthYet;
    };
    let declared: [u8; COMPRESSED_BLOCK_LENGTH_BYTES] = declared
        .try_into()
        .expect("a slice of exactly the length prefix's width is an array of it");
    let frame_bytes = u32::from_le_bytes(declared) as usize;
    // `usize` is at least 32 bits wide everywhere this builds, and the declaration is a `u32`,
    // so the sum is exact on a 64-bit target and checked rather than assumed on a 32-bit one.
    let Some(block_bytes) = COMPRESSED_BLOCK_LENGTH_BYTES.checked_add(frame_bytes) else {
        return CompressedBlockAt::NotAllHere {
            block_bytes: usize::MAX,
        };
    };
    match bytes.get(COMPRESSED_BLOCK_LENGTH_BYTES..block_bytes) {
        Some(zstd_frame) => CompressedBlockAt::Whole {
            zstd_frame,
            block_bytes,
        },
        None => CompressedBlockAt::NotAllHere { block_bytes },
    }
}

/// Why a block could not be compressed.
///
/// **Every variant is the writer's own configuration or its own output**, not a fault in data
/// it was handed — unlike every other error in this module, which is about a record or a file a
/// caller supplied.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum BlockCompressError {
    /// The look-back window the file declares is outside what zstd takes. `header.rs` refuses
    /// such a manifest; this is the same rule where the parameter is set, because `Manifest`'s
    /// fields are public and a caller can build one that never met that check.
    #[error(
        "a look-back window of 2^{look_back_window_log} bytes; zstd takes between \
         2^{MIN_LOOK_BACK_WINDOW_LOG} and 2^{MAX_LOOK_BACK_WINDOW_LOG}"
    )]
    WindowLogOutOfRange { look_back_window_log: u8 },

    /// A compression level outside what zstd takes. **Refused rather than clamped**, which is
    /// what zstd does silently: a level of 100 gives exactly the level-22 bytes, and a level of
    /// −131,073 gives a block 169 times the size, with nothing anywhere saying so.
    #[error("a compression level of {level}; zstd takes between {lowest} and {highest}")]
    LevelOutOfRange {
        level: i32,
        lowest: i32,
        highest: i32,
    },

    /// A compressed frame longer than the four-byte length in front of it can describe.
    ///
    /// **Not reachable from a block this format produces** — it needs more than four gibibytes
    /// of compressed output from one block — and refused rather than truncated for the reason
    /// the record head's own byte count is: a length that means something else is a block
    /// nothing could skip.
    #[error(
        "a compressed block of {frame_bytes} bytes from a {payload_bytes}-byte payload, longer \
         than the {} a block's length can describe",
        u32::MAX
    )]
    FrameTooLong {
        frame_bytes: usize,
        payload_bytes: usize,
    },

    /// zstd refused. **It names what was being done**, because zstd's own message is a code.
    #[error("zstd failed while {while_doing}")]
    Zstd {
        while_doing: &'static str,
        #[source]
        source: std::io::Error,
    },
}

impl BlockCompressError {
    fn zstd(while_doing: &'static str, source: std::io::Error) -> Self {
        Self::Zstd {
            while_doing,
            source,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::psp::header::{
        DEFAULT_BLOCK_BYTE_CEILING, DEFAULT_GENOMIC_BLOCK_SIZE_BP, DEFAULT_LOOK_BACK_WINDOW_LOG,
    };
    use crate::ng::psp::record::{OffsetBase, RecordLayout, decode_record, record_fields};
    use crate::ng::types::{ReadGroupId, SummedLogError};
    use proptest::prelude::*;

    /// The grid these tests cut on. **A fixture of its own and not the shipped default**: every
    /// coordinate below is chosen against 100,000, and spec §4.1 records that the default may
    /// yet move to 1,000 kb — at which point a change of mind about the default should fail the
    /// one test that is about the default, not five that are about the cut.
    const A_GRID: Bp = Bp(100_000);

    /// One ordinary record: a covered position whose reads all agreed with the reference.
    ///
    /// **The body is here to be round-tripped, not inspected.** The builder cuts on
    /// coordinates and never reads inside a record, so nothing in this file can tell one
    /// record's payload from another's — and no fixture could make it able to: over a
    /// four-letter alphabet, two one-base records must sometimes carry the same base. Whether
    /// a decode can return the wrong record's payload is `record.rs`'s property test, which
    /// generates arbitrary bytes. *An earlier version of this comment claimed the opposite and
    /// nothing held it; the claim was not merely unheld but unachievable.*
    fn a_record(contig: u32, start: u64, span: u64) -> SampleLocusObservations {
        debug_assert!(span > 0, "a record covers at least one reference base");
        let bases: Box<[u8]> = (0..span)
            .map(|offset| b"ACGT"[((start + offset) % 4) as usize])
            .collect::<Vec<_>>()
            .into_boxed_slice();
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(start + span - 1),
            },
            reference_bases: bases.clone(),
            observations: vec![SequenceObservation {
                bases,
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 2,
                q_sum: SummedLogError::from_steps(-4_096),
                mapq_sum: 180,
                mapq_sum_sq: 10_800,
                placed_left: 1,
                chain_ids: Vec::new(),
            }],
            reads_without_observation: 1,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// A region covering no base, which the record codec refuses.
    fn a_record_over_no_base(contig: u32, start: u64) -> SampleLocusObservations {
        debug_assert!(
            start > 0,
            "an empty region is built by stepping `end` back one"
        );
        let mut empty = a_record(contig, start, 1);
        empty.region.end = Position(start - 1);
        empty
    }

    /// Push every record and collect the block payloads, so a test can say what the cut did
    /// without repeating the drive loop.
    fn cut(
        mut builder: BlockBuilder,
        records: &[SampleLocusObservations],
    ) -> Result<Vec<Vec<u8>>, BlockWriteError> {
        let mut blocks = Vec::new();
        for record in records {
            if let Some(closed) = builder.push(record)? {
                blocks.push(closed.to_vec());
            }
        }
        if let Some(last) = builder.finish() {
            blocks.push(last);
        }
        Ok(blocks)
    }

    /// A block payload walked back: what its head said, and every record in it.
    #[derive(Debug, Clone, PartialEq)]
    struct WalkedBlock {
        head: BlockHead,
        records: Vec<SampleLocusObservations>,
    }

    /// Walk one block payload back, each position rebuilt from the block's own first position
    /// and the offsets since.
    ///
    /// **This is what the cut is checked against, and it takes nothing from the builder**: it
    /// starts from the block's declared first position, which is all a reader beginning at that
    /// block would have.
    fn walk(payload: &[u8]) -> WalkedBlock {
        let found = BlockRecords::split(payload).expect("the block head reads");
        let layout = RecordLayout::as_this_build_writes_it();
        let mut measured_from = OffsetBase::at_block_start(found.head.first_position);
        let mut at = 0usize;
        let mut records = Vec::new();
        while at < found.records.len() {
            let decoded = decode_record(
                &found.records[at..],
                found.head.contig,
                measured_from,
                &layout,
            )
            .unwrap_or_else(|refused| panic!("record {} reads: {refused}", records.len()));
            at += decoded.record_bytes;
            measured_from = OffsetBase::after(&decoded.head);
            records.push(decoded.record);
        }
        assert_eq!(
            at,
            found.records.len(),
            "the walk consumed the block exactly"
        );
        assert_eq!(
            records.len() as u64,
            found.head.record_count.get(),
            "the block held the number of records its head declared"
        );
        WalkedBlock {
            head: found.head,
            records,
        }
    }

    /// Every block's records, in the order the blocks were cut.
    fn walk_all(blocks: &[Vec<u8>]) -> Vec<SampleLocusObservations> {
        blocks
            .iter()
            .flat_map(|block| walk(block).records)
            .collect()
    }

    /// How many bytes of records a block payload holds, past its own head.
    fn record_bytes_in(payload: &[u8]) -> usize {
        BlockRecords::split(payload)
            .expect("the head reads")
            .records
            .len()
    }

    /// What one `a_record(0, n, 1)` costs on the wire, measured rather than guessed — so a
    /// ceiling derived from it cannot stop firing because a record changed size.
    fn one_records_bytes() -> u32 {
        let alone = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[a_record(0, 100, 1)],
        )
        .expect("in order");
        u32::try_from(record_bytes_in(&alone[0])).expect("one record is a small number of bytes")
    }

    fn a_head(contig: u32, first_position: u64, record_count: u64) -> BlockHead {
        BlockHead {
            contig: ContigId(contig),
            first_position: Position(first_position),
            record_count: NonZeroU64::new(record_count).expect("a block holds a record"),
        }
    }

    // -----------------------------------------------------------------
    // The block head
    // -----------------------------------------------------------------

    #[test]
    fn a_block_head_round_trips_and_says_how_many_bytes_it_took() {
        let head = a_head(7, 90_600_000, 24_881);
        let mut bytes = Vec::new();
        head.encode(&mut bytes);
        let read = BlockHead::decode(&bytes).expect("it reads back");
        assert_eq!(read.head, head);
        assert_eq!(read.head_bytes, bytes.len());
    }

    /// **The bytes a file carries.** A block head is container framing, not a manifest field,
    /// so nothing in a file says how it is laid out — which makes this array the format itself.
    /// A reordering or a change of encoding has to fail here and be a version bump.
    #[test]
    fn a_block_head_is_these_exact_bytes() {
        let mut bytes = Vec::new();
        a_head(1, 300, 2).encode(&mut bytes);
        assert_eq!(
            bytes,
            vec![
                0x01, // contig-id 1
                0xac, 0x02, // first-position 300
                0x02, // record-count 2
            ]
        );
    }

    /// A block head cut short is `Truncated` at every cut and never damage: a streaming reader
    /// fetches more bytes for the first and refuses the file for the second, and at the start
    /// of a block a short buffer is the ordinary state of affairs.
    ///
    /// **The field and the offset are asserted too.** They are the whole reason those two
    /// members exist, and hard-coding both to constants left the suite green.
    #[test]
    fn a_block_head_cut_short_is_truncated_at_every_cut_and_never_malformed() {
        let mut whole = Vec::new();
        a_head(300, 90_600_000, 24_881).encode(&mut whole);
        // contig 300 is 2 bytes, first-position 90,600,000 is 4, record-count 24,881 is 3.
        assert_eq!(
            whole.len(),
            9,
            "the cuts below are placed against this layout; it is {whole:02x?}"
        );
        let field_and_offset = |cut_at: usize| match cut_at {
            0..=1 => (BLOCK_CONTIG_ID, 0usize),
            2..=5 => (BLOCK_FIRST_POSITION, 2),
            _ => (BLOCK_RECORD_COUNT, 6),
        };

        for cut_at in 0..whole.len() {
            match BlockHead::decode(&whole[..cut_at]) {
                Err(BlockHeadDecodeError::Truncated { field, bytes_in }) => {
                    assert_eq!(
                        (field, bytes_in),
                        field_and_offset(cut_at),
                        "a head cut at {cut_at} bytes"
                    );
                }
                other => panic!("a head of {cut_at} bytes gave {other:?}"),
            }
        }
        assert!(
            BlockHead::decode(&whole).is_ok(),
            "and the whole head reads"
        );
    }

    /// **The other half of the retry rule, and the half that was missing.** A head cut short
    /// says *fetch more bytes*; a head whose variable-length integer no 64-bit value could hold
    /// says *this file is damaged*, however many bytes arrive. A reader handed the second as
    /// the first retries for ever on a block that will never read.
    #[test]
    fn a_varint_no_u64_could_hold_is_malformed_and_never_truncated() {
        for (field, before) in [
            (BLOCK_CONTIG_ID, Vec::new()),
            (BLOCK_FIRST_POSITION, vec![0x01u8]),
            (BLOCK_RECORD_COUNT, vec![0x01u8, 0x01u8]),
        ] {
            for trailing in [0usize, 4_096] {
                let mut bytes = before.clone();
                bytes.extend_from_slice(&[0x80u8; 12]);
                bytes.extend(std::iter::repeat_n(0x00u8, trailing));

                match BlockHead::decode(&bytes) {
                    Err(BlockHeadDecodeError::Malformed {
                        field: broke,
                        bytes_in,
                        reason,
                    }) => {
                        assert_eq!(broke, field);
                        assert_eq!(bytes_in, before.len());
                        assert!(
                            reason.contains("longer than any 64-bit value"),
                            "got {reason}"
                        );
                    }
                    other => panic!("a {field} of twelve continuation bytes gave {other:?}"),
                }
            }
        }
    }

    /// A block claiming no records is damage. Nothing this builder writes can say it — the
    /// count's type has no zero — so a file that does is not one this writer produced.
    #[test]
    fn a_block_claiming_no_records_is_refused() {
        let mut bytes = Vec::new();
        a_head(1, 300, 1).encode(&mut bytes);
        let last = bytes.len() - 1;
        bytes[last] = 0;

        match BlockHead::decode(&bytes) {
            Err(BlockHeadDecodeError::Malformed { field, reason, .. }) => {
                assert_eq!(field, BLOCK_RECORD_COUNT);
                assert!(reason.contains("no records"), "got {reason}");
            }
            other => panic!("expected Malformed, got {other:?}"),
        }
    }

    /// A contig id past what a contig id is, is damage rather than a narrowed number.
    #[test]
    fn a_contig_id_too_large_for_its_field_is_refused_rather_than_narrowed() {
        let mut bytes = Vec::new();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut bytes);
        encode_u64_leb128(300, &mut bytes);
        encode_u64_leb128(2, &mut bytes);

        match BlockHead::decode(&bytes) {
            Err(BlockHeadDecodeError::Malformed { field, reason, .. }) => {
                assert_eq!(field, BLOCK_CONTIG_ID);
                assert!(reason.contains(&u32::MAX.to_string()), "got {reason}");
            }
            other => panic!("expected Malformed, got {other:?}"),
        }
    }

    /// **A whole block whose head stops early is damage, never `Truncated`.** The caller holds
    /// the block entire, so *fetch more bytes* is a retry that never ends — the C2 review's
    /// Blocker one level up, and the reason `split` exists as more than a two-line forward.
    #[test]
    fn a_short_whole_block_is_malformed_and_never_truncated() {
        let mut whole = Vec::new();
        a_head(300, 90_600_000, 24_881).encode(&mut whole);

        for cut_at in 0..whole.len() {
            match BlockRecords::split(&whole[..cut_at]) {
                Err(BlockHeadDecodeError::Malformed { reason, .. }) => {
                    assert!(
                        reason.contains(&format!("{cut_at} bytes this whole block holds")),
                        "at {cut_at} bytes, got {reason}"
                    );
                }
                other => panic!("a whole block of {cut_at} bytes gave {other:?}"),
            }
        }
        assert!(
            BlockRecords::split(&whole).is_ok(),
            "and the whole head splits"
        );
    }

    /// Damage stays damage: `split` passes a `Malformed` through rather than restating it, so
    /// the field and the reason a reader acts on are the head's own.
    #[test]
    fn split_passes_a_damaged_head_through_unchanged() {
        let mut bytes = Vec::new();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut bytes);
        encode_u64_leb128(300, &mut bytes);
        encode_u64_leb128(2, &mut bytes);

        assert_eq!(
            BlockRecords::split(&bytes).expect_err("a contig id no contig id could be"),
            BlockHead::decode(&bytes).expect_err("the same fault")
        );
    }

    // -----------------------------------------------------------------
    // The cut
    // -----------------------------------------------------------------

    /// **The cut is a grid on the coordinate, not a running total** — the property every other
    /// one here rests on. Two records 99,998 bases apart share a block when both fall in the
    /// same multiple of 100,000; two records one base apart are in different blocks when the
    /// second crosses the multiple.
    #[test]
    fn the_cut_is_a_grid_on_the_coordinate_and_not_a_running_total() {
        let blocks = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[a_record(0, 1, 1), a_record(0, 99_999, 1)],
        )
        .expect("both records are in order");
        assert_eq!(
            blocks.len(),
            1,
            "99,998 bases apart, both under 100,000: one block"
        );

        let blocks = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[a_record(0, 99_999, 1), a_record(0, 100_000, 1)],
        )
        .expect("both records are in order");
        assert_eq!(blocks.len(), 2, "one base apart, across the multiple: two");
        assert_eq!(walk(&blocks[0]).head.first_position, Position(99_999));
        assert_eq!(walk(&blocks[1]).head.first_position, Position(100_000));
    }

    /// **A record whose span reaches past a grid multiple belongs to the cell its *start* falls
    /// in.** A span is sample-dependent — a deletion widens a locus in one sample and not in
    /// another — so a cut taken from the record's end would put the same reference span in
    /// different blocks in different samples, which is the one thing the grid exists to prevent
    /// (spec §4.1). Every file would still read back self-consistently, so only a cohort read
    /// across samples would ever notice.
    #[test]
    fn a_record_that_straddles_a_grid_multiple_is_cut_by_its_start_and_not_its_end() {
        let straddling = a_record(0, 99_998, 5);
        assert_eq!(
            straddling.region.end,
            Position(100_002),
            "the fixture must reach past the multiple, or it tests nothing"
        );

        let blocks = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[straddling, a_record(0, 99_999, 1)],
        )
        .expect("in order");

        assert_eq!(
            blocks.len(),
            1,
            "both records start under 100,000, so both belong to the same block"
        );
        assert_eq!(walk(&blocks[0]).head.first_position, Position(99_998));
    }

    /// **Every sample cuts at the same coordinates**, which is what the grid buys and what a
    /// running total would lose. Two samples covering the same region with entirely different
    /// records — including, on the sparse one, a record widened across a multiple — give blocks
    /// whose first positions fall in the same grid cells, and the same number of blocks.
    #[test]
    fn two_samples_with_different_records_cut_at_the_same_grid_cells() {
        let dense: Vec<_> = (0..170)
            .map(|index| a_record(0, 90_000 + index * 1_000, 1))
            .collect();
        let sparse: Vec<_> = [
            (90_500u64, 1u64),
            // Widened across the 100,000 multiple, which is where a cut taken from the end
            // would put this sample's boundary somewhere the dense sample has none.
            (99_998, 5),
            (150_000, 1),
            (250_000, 1),
            (250_001, 1),
        ]
        .into_iter()
        .map(|(start, span)| a_record(0, start, span))
        .collect();

        let cells_of = |blocks: &[Vec<u8>]| -> Vec<u64> {
            blocks
                .iter()
                .map(|block| walk(block).head.first_position.get() / A_GRID.get())
                .collect()
        };

        let a_builder = || BlockBuilder::new(A_GRID, None).expect("a grid");
        let dense_blocks = cut(a_builder(), &dense).expect("in order");
        let sparse_blocks = cut(a_builder(), &sparse).expect("in order");

        assert_eq!(cells_of(&dense_blocks), vec![0, 1, 2]);
        assert_eq!(cells_of(&sparse_blocks), vec![0, 1, 2]);
        assert_eq!(
            dense_blocks.len(),
            sparse_blocks.len(),
            "the two samples cut the same number of blocks over the same span"
        );
    }

    /// A block never crosses a contig, even when both records fall in the same grid cell — the
    /// property that lets a reader starting at a block know its chromosome from the block.
    #[test]
    fn a_block_never_crosses_a_contig_even_inside_one_grid_cell() {
        let blocks = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[a_record(0, 500, 1), a_record(1, 500, 1)],
        )
        .expect("ascending contigs are in order");

        assert_eq!(blocks.len(), 2);
        assert_eq!(walk(&blocks[0]).head.contig, ContigId(0));
        assert_eq!(walk(&blocks[1]).head.contig, ContigId(1));
    }

    /// The byte ceiling closes a block early, and **the same run without one stays whole** — so
    /// this cannot pass because the grid happened to cut in the same places.
    #[test]
    fn the_byte_ceiling_closes_a_block_early_and_nothing_else_does() {
        let records: Vec<_> = (0..20).map(|index| a_record(0, 100 + index, 1)).collect();

        let whole =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");
        assert_eq!(whole.len(), 1, "the grid alone leaves these twenty in one");

        let cut_up = cut(
            BlockBuilder::new(A_GRID, Some(60)).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert!(
            cut_up.len() > 1,
            "a 60-byte ceiling over {} bytes of records must close early",
            record_bytes_in(&whole[0])
        );
        assert_eq!(
            walk_all(&cut_up),
            walk_all(&whole),
            "and the records are the same records, however they were cut"
        );
    }

    /// **The ceiling closes a block once its records have *reached* it, not once they have
    /// passed it** — the one byte at which the two rules differ, and the boundary that decides
    /// what a writer's memory budget actually buys.
    #[test]
    fn the_ceiling_closes_a_block_when_the_records_reach_it_and_not_one_record_later() {
        let records: Vec<_> = (0..8).map(|index| a_record(0, 100 + index, 1)).collect();
        let one_record = one_records_bytes();
        let a_run = |ceiling: u32| {
            cut(
                BlockBuilder::new(A_GRID, Some(ceiling)).expect("a grid"),
                &records,
            )
            .expect("in order")
        };

        assert_eq!(
            a_run(one_record).len(),
            records.len(),
            "a ceiling of exactly one record's {one_record} bytes is reached by that one record"
        );
        assert_eq!(
            a_run(one_record + 1).len(),
            records.len() / 2,
            "one byte more and it takes two records to reach it"
        );
    }

    /// **The ceiling measures the records laid down and not the block head in front of them.**
    /// A ceiling set to what two records cost closes after the second; were the head counted
    /// against it, it would close after the first.
    #[test]
    fn the_ceiling_measures_the_records_and_not_the_head_in_front_of_them() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();
        let two_records = one_records_bytes() * 2;

        let blocks = cut(
            BlockBuilder::new(A_GRID, Some(two_records)).expect("a grid"),
            &records,
        )
        .expect("in order");

        assert_eq!(blocks.len(), 3, "six records, two to a block");
        for block in &blocks {
            assert_eq!(walk(block).records.len(), 2, "two records a block, not one");
        }
    }

    /// **`None` means no ceiling at all, not a large one.** Spec §12 question 2 is open on what
    /// a ceiling should be, so the value that has *not* been chosen is the one to pin: a grid
    /// cell holding two megabytes of records is one block, and stays one however deep the
    /// sample. A "sensible" fallback added later is exactly what this refuses.
    #[test]
    fn a_ceiling_of_none_never_closes_a_block_however_large_it_grows() {
        let records: Vec<_> = (1..90_000u64).map(|at| a_record(0, at, 1)).collect();
        let blocks =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");

        let bytes: usize = blocks.iter().map(|block| record_bytes_in(block)).sum();
        assert!(
            bytes > 1 << 20,
            "the fixture must pass any plausible fallback to say anything; it is {bytes} bytes"
        );
        assert_eq!(blocks.len(), 1, "one grid cell, no ceiling: one block");
    }

    /// **The oracle for the whole cut: every record comes back, once, in order, wherever the
    /// blocks fell.** Records over three contigs, three grid cells each, one record per cell
    /// widened across the cell's own upper boundary, and a byte ceiling that fires inside them.
    #[test]
    fn every_record_comes_back_once_and_in_order_however_the_blocks_fell() {
        let mut records = Vec::new();
        for contig in 0..3u32 {
            for cell in 0..3u64 {
                for step in 0..7u64 {
                    records.push(a_record(
                        contig,
                        cell * A_GRID.get() + 1 + step * 3_000,
                        1 + step % 4,
                    ));
                }
                // The last record of each cell reaches past the cell's own end, so a cut taken
                // from a record's end would put it in the next cell.
                records.push(a_record(contig, (cell + 1) * A_GRID.get() - 2, 5));
            }
        }

        let by_grid =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");
        assert_eq!(by_grid.len(), 9, "three contigs of three grid cells each");
        assert_eq!(walk_all(&by_grid), records);

        // A ceiling measured from what the blocks actually hold, so it must fire whatever a
        // record's size turns out to be. Guessing one is how a test that proves nothing about
        // the ceiling passes: an earlier draft guessed 200 bytes and never fired.
        let smallest = by_grid
            .iter()
            .map(|block| record_bytes_in(block))
            .min()
            .expect("nine blocks");
        let ceiling = u32::try_from(smallest / 3).expect("a third of a block's records");
        assert!(ceiling >= 1, "a ceiling of zero is refused");

        let blocks = cut(
            BlockBuilder::new(A_GRID, Some(ceiling)).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert!(
            blocks.len() > by_grid.len(),
            "a {ceiling}-byte ceiling over blocks of at least {smallest} bytes must cut further \
             than the grid's {} blocks; it gave {}",
            by_grid.len(),
            blocks.len()
        );
        assert_eq!(walk_all(&blocks), records);
    }

    /// Every block's head describes its own records: the contig they are all on, the position
    /// its first record starts at, and how many it holds.
    #[test]
    fn a_blocks_head_describes_the_records_behind_it() {
        let records: Vec<_> = (0..6u32)
            .map(|index| a_record(index / 3, 100 + u64::from(index), 1))
            .collect();
        let blocks = cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records)
            .expect("ascending contigs");

        let mut seen = 0usize;
        for block in &blocks {
            let walked = walk(block);
            assert_eq!(walked.head.record_count.get(), walked.records.len() as u64);
            assert_eq!(walked.head.first_position, walked.records[0].region.start);
            for record in &walked.records {
                assert_eq!(record.region.contig, walked.head.contig);
            }
            seen += walked.records.len();
        }
        assert_eq!(seen, records.len());
    }

    // -----------------------------------------------------------------
    // What it refuses
    // -----------------------------------------------------------------

    /// **The defect the cut would otherwise launder.** A record going backwards far enough to
    /// land in an earlier grid cell resets the encoder's own base, so the record codec cannot
    /// see it at all: without this check the pair below is written without complaint, and the
    /// index built over it seeks wrongly.
    #[test]
    fn a_record_that_goes_backwards_across_a_grid_cell_is_refused() {
        let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
        builder.push(&a_record(0, 250_000, 1)).expect("the first");

        match builder.push(&a_record(0, 150_000, 1)) {
            Err(BlockWriteError::OutOfOrder { previous, offered }) => {
                assert_eq!(previous.start, Position(250_000));
                assert_eq!(offered.start, Position(150_000));
            }
            other => panic!("expected OutOfOrder, got {other:?}"),
        }
    }

    /// And one going backwards inside a block is refused too — here rather than by the encoder,
    /// so the message names both regions the way the file's own error will.
    #[test]
    fn a_record_that_goes_backwards_inside_a_block_is_refused() {
        let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
        builder.push(&a_record(0, 500, 1)).expect("the first");

        let refused = builder
            .push(&a_record(0, 400, 1))
            .expect_err("a record before the one before it");
        assert_eq!(
            refused.to_string(),
            "a record at contig 0:400-400 starts before the previous record at contig 0:500-500"
        );
    }

    /// Two records starting on the same base are not out of order — a repeat tract and a
    /// generic locus can begin together — so the check must not refuse them.
    #[test]
    fn two_records_starting_on_the_same_base_are_accepted() {
        let blocks = cut(
            BlockBuilder::new(A_GRID, None).expect("a grid"),
            &[a_record(0, 500, 3), a_record(0, 500, 1)],
        )
        .expect("the same start twice is not out of order");
        assert_eq!(walk_all(&blocks).len(), 2);
    }

    /// A contig already left is refused, and so is one that goes backwards. Blocks are indexed
    /// in genomic order, and a contig visited twice gives two runs of blocks a seek cannot
    /// choose between.
    ///
    /// **The message is rendered, not only the fields**: it is the only refusal a writer sees
    /// for a contig fault, and its whole content is which contig came after which.
    #[test]
    fn a_contig_already_left_or_one_that_goes_backwards_is_refused() {
        for second in [0u32, 1] {
            let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
            builder.push(&a_record(0, 500, 1)).expect("contig 0");
            builder.push(&a_record(2, 500, 1)).expect("contig 2");

            match builder.push(&a_record(second, 500, 1)) {
                Err(refused @ BlockWriteError::ContigOutOfOrder { previous, offered }) => {
                    assert_eq!(previous, ContigId(2));
                    assert_eq!(offered, ContigId(second));
                    assert_eq!(
                        refused.to_string(),
                        format!(
                            "a record on contig {second} after a record on contig 2; contigs are \
                             written in ascending order and never revisited"
                        )
                    );
                }
                other => panic!("contig {second} after contig 2 gave {other:?}"),
            }
        }
    }

    /// **A refused record leaves the builder exactly as it was.** A run with three kinds of
    /// refusal offered before every record — one of them where a cut is about to happen — gives
    /// byte-identical blocks to the same run without them.
    #[test]
    fn a_refused_record_changes_nothing_a_later_record_can_see() {
        // Eight records on contig 1 that cross a grid boundary, with a ceiling that also
        // fires, so the refusals land next to both kinds of cut.
        let good: Vec<_> = (0..8)
            .map(|index| a_record(1, 99_000 + index * 500, 1))
            .collect();

        let without =
            cut(BlockBuilder::new(A_GRID, Some(40)).expect("a grid"), &good).expect("in order");
        assert!(
            without.len() >= 3,
            "the fixture must cut more than once for this to say anything; it cut into {} blocks",
            without.len()
        );

        let mut interrupted = BlockBuilder::new(A_GRID, Some(40)).expect("a grid");
        let mut with = Vec::new();
        assert!(
            interrupted.push(&good[0]).expect("the first").is_none(),
            "the first record opens a block rather than closing one"
        );
        for record in &good[1..] {
            assert!(matches!(
                interrupted.push(&a_record(1, 1, 1)),
                Err(BlockWriteError::OutOfOrder { .. })
            ));
            assert!(matches!(
                interrupted.push(&a_record(0, 300_000, 1)),
                Err(BlockWriteError::ContigOutOfOrder { .. })
            ));
            // Start 300,000 is a cell of its own, so this refusal reaches the codec through the
            // cut path — the one the rollback of the encoder's base defends.
            assert!(matches!(
                interrupted.push(&a_record_over_no_base(1, 300_000)),
                Err(BlockWriteError::RecordRefused(
                    RecordEncodeError::EmptyRegion { .. }
                ))
            ));
            if let Some(closed) = interrupted.push(record).expect("in order") {
                with.push(closed.to_vec());
            }
        }
        if let Some(last) = interrupted.finish() {
            with.push(last);
        }

        assert_eq!(with, without, "the refusals left no trace in the bytes");
    }

    /// **A record the codec refuses *inside* an open block leaves that block untouched** — the
    /// one refusal path with no second buffer in front of it. The block's own count is what is
    /// at risk: a record counted or half-written before it is refused leaves a head declaring
    /// one more record than the block holds, which a reader meets as bytes running out.
    #[test]
    fn a_codec_refusal_inside_the_open_block_leaves_it_untouched() {
        let good: Vec<_> = (0..4).map(|index| a_record(0, 500 + index, 1)).collect();
        let without =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &good).expect("in order");
        assert_eq!(
            without.len(),
            1,
            "the fixture is one block, so no cut can absorb the refusal"
        );

        let mut interrupted = BlockBuilder::new(A_GRID, None).expect("a grid");
        assert!(interrupted.push(&good[0]).expect("the first").is_none());
        for record in &good[1..] {
            // Same contig, same grid cell, and after the last accepted start: the cut path is
            // not taken, so this is the codec refusing with the block open.
            assert!(matches!(
                interrupted.push(&a_record_over_no_base(0, 600)),
                Err(BlockWriteError::RecordRefused(
                    RecordEncodeError::EmptyRegion { .. }
                ))
            ));
            assert!(interrupted.push(record).expect("in order").is_none());
        }

        assert_eq!(
            vec![interrupted.finish().expect("a block")],
            without,
            "the refusals left no trace in the bytes"
        );
    }

    /// A grid with no cells is refused rather than divided by, and a ceiling no block can stay
    /// under is refused rather than acted on. `header.rs` refuses both in a file's manifest;
    /// this is the same rule where the arithmetic happens, because `Manifest`'s fields are
    /// public and a caller can build one that never met that check.
    #[test]
    fn a_cut_rule_that_cannot_cut_is_refused() {
        let refused = BlockBuilder::new(Bp(0), None).expect_err("a zero grid has no cells");
        assert_eq!(refused, BlockCutRuleError::ZeroGenomicBlockSize);
        assert!(refused.to_string().contains("zero grid has no cells"));

        let refused =
            BlockBuilder::new(A_GRID, Some(0)).expect_err("a zero ceiling closes every block");
        assert_eq!(refused, BlockCutRuleError::ZeroBlockByteCeiling);
        assert!(refused.to_string().contains("a block of its own"));

        let mut manifest = a_manifest();
        manifest.block_byte_ceiling = Some(0);
        assert_eq!(
            BlockBuilder::from_manifest(&manifest).expect_err("the same rule through a manifest"),
            BlockCutRuleError::ZeroBlockByteCeiling
        );
    }

    /// A refusal on the very first record leaves no block behind it, so the next record still
    /// opens the file's first block.
    #[test]
    fn a_refusal_on_the_first_record_leaves_no_block_open() {
        let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
        assert!(builder.push(&a_record_over_no_base(0, 500)).is_err());
        let blocks = cut(builder, &[a_record(0, 500, 1)]).expect("in order");
        assert_eq!(blocks.len(), 1);
        assert_eq!(walk(&blocks[0]).head.first_position, Position(500));
        assert_eq!(walk(&blocks[0]).records.len(), 1);
    }

    /// A builder that was never pushed to has no block to close.
    ///
    /// **Closing twice, and pushing after closing, are not tested because they cannot be
    /// written**: closing consumes the builder. Both were live hazards while it did not — the
    /// second is what leaked one block's records into the next.
    #[test]
    fn finishing_a_builder_that_was_never_pushed_to_gives_nothing() {
        let never_used = BlockBuilder::new(A_GRID, None).expect("a grid");
        assert!(never_used.finish().is_none());
    }

    // -----------------------------------------------------------------
    // The manifest, and the values that ship
    // -----------------------------------------------------------------

    fn a_manifest() -> Manifest {
        Manifest::as_this_build_writes_it()
    }

    /// The manifest is where a file's cut rule lives, so a builder made from one cuts the way
    /// the numbers in it say — **both of them**, the grid and the ceiling.
    #[test]
    fn a_builder_from_a_manifest_cuts_on_the_manifests_grid_and_ceiling() {
        let mut manifest = a_manifest();
        manifest.genomic_block_size_bp = Bp(1_000);
        let blocks = cut(
            BlockBuilder::from_manifest(&manifest).expect("a grid"),
            &[a_record(0, 999, 2), a_record(0, 1_000, 1)],
        )
        .expect("in order");
        assert_eq!(blocks.len(), 2, "a 1,000 bp grid cuts at 1,000");

        manifest.genomic_block_size_bp = A_GRID;
        manifest.block_byte_ceiling = Some(1);
        let records: Vec<_> = (0..4).map(|index| a_record(0, 100 + index, 1)).collect();
        let blocks = cut(
            BlockBuilder::from_manifest(&manifest).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert_eq!(
            blocks.len(),
            records.len(),
            "a one-byte ceiling from the manifest gives each record a block of its own"
        );
    }

    /// A file declaring a record layout this build does not write cannot be extended with
    /// records this build produces, so the builder that would extend it is refused rather than
    /// producing records unreadable under the header the file keeps.
    #[test]
    fn a_manifest_declaring_a_layout_this_writer_cannot_honour_is_refused() {
        let mut manifest = a_manifest();
        manifest.fields.truncate(manifest.fields.len() - 1);

        let refused = BlockBuilder::from_manifest(&manifest)
            .expect_err("a manifest missing a field this build writes");
        assert!(
            matches!(refused, BlockCutRuleError::UnsupportedRecordLayout { .. }),
            "got {refused:?}"
        );
        assert!(
            refused.to_string().contains("cannot honour"),
            "got {refused}"
        );
    }

    /// The values a writer takes when nothing says otherwise, pinned where a change to one is a
    /// decision taken here rather than five cut tests that stop passing. Spec §4.1 records that
    /// the grid is "not an optimum" and that 1,000 kb is live, and spec §12 question 2 leaves
    /// the ceiling open — so both are expected to move, and neither should move by accident.
    #[test]
    fn the_values_a_writer_ships_with_are_these() {
        assert_eq!(DEFAULT_GENOMIC_BLOCK_SIZE_BP, Bp(100_000));
        assert_eq!(DEFAULT_BLOCK_BYTE_CEILING, None);
        let manifest = a_manifest();
        assert_eq!(
            manifest.genomic_block_size_bp,
            DEFAULT_GENOMIC_BLOCK_SIZE_BP
        );
        assert_eq!(manifest.block_byte_ceiling, DEFAULT_BLOCK_BYTE_CEILING);
        assert_eq!(manifest.look_back_window_log, DEFAULT_LOOK_BACK_WINDOW_LOG);
        assert_eq!(manifest.fields, record_fields());
    }

    // -----------------------------------------------------------------
    // Compressing a block
    // -----------------------------------------------------------------

    /// Why a decoder would not inflate a frame. **Two states and not one**, because a test that
    /// cannot tell them apart proves nothing: asking for a window below zstd's own minimum is
    /// refused at the parameter, before the frame is looked at, so an assertion of "this
    /// decoder refuses it" is satisfied by any bytes at all. Measured — that is exactly what
    /// half of the first window test did before this type existed.
    #[derive(Debug, PartialEq, Eq)]
    enum NotInflated {
        /// zstd would not accept the window cap itself.
        TheCapWasRefused,
        /// The cap was accepted and the frame was refused under it.
        TheFrameWasRefused,
    }

    /// Decompress a whole zstd frame with the decoder's window capped at `window_log_max`,
    /// which is what a reader configured from the file's declaration would do.
    ///
    /// **It grows its output rather than sizing it from the frame**, which is the property the
    /// design turns on: the frame carries no content size, so there is nothing to size from
    /// even if this wanted to. It is not the reader — that streams into a rolling buffer at
    /// Milestone D3 — it is just enough to ask zstd whether the window it was promised is the
    /// window the frame needs.
    fn decompress_capped_at(frame: &[u8], window_log_max: u32) -> Result<Vec<u8>, NotInflated> {
        let mut decoder = zstd::zstd_safe::DCtx::create();
        decoder
            .set_parameter(zstd::zstd_safe::DParameter::WindowLogMax(window_log_max))
            .map_err(|_| NotInflated::TheCapWasRefused)?;

        let mut out: Vec<u8> = Vec::new();
        let mut input = zstd::zstd_safe::InBuffer::around(frame);
        loop {
            let at = out.len();
            out.reserve(64 * 1024);
            let mut output = zstd::zstd_safe::OutBuffer::around_pos(&mut out, at);
            let hint = decoder
                .decompress_stream(&mut output, &mut input)
                .map_err(|_| NotInflated::TheFrameWasRefused)?;
            if hint == 0 {
                break;
            }
            if input.pos == frame.len() && output.pos() == at {
                return Err(NotInflated::TheFrameWasRefused);
            }
        }
        Ok(out)
    }

    /// The zstd frame of a whole on-disk block, and a check that the block is exactly what its
    /// own length prefix declares.
    fn frame_of(block: &[u8]) -> &[u8] {
        match compressed_block_at(block) {
            CompressedBlockAt::Whole {
                zstd_frame,
                block_bytes,
            } => {
                assert_eq!(block_bytes, block.len(), "the block is what it declares");
                zstd_frame
            }
            other => panic!("a whole block was expected; got {other:?}"),
        }
    }

    /// Enough records to fill more than one look-back window, so a frame's declared window is
    /// the one it was capped at rather than one zstd narrowed to fit a small payload.
    ///
    /// **The slack is one record, not sixty-four.** An oversized fixture is what hid the
    /// `compress_bound` finding: every compression test built forty records or more, so nothing
    /// ever reached a block zstd cannot shrink.
    fn records_past_a_window(window_log: u8) -> Vec<SampleLocusObservations> {
        let want = 1usize << window_log;
        let one = one_records_bytes() as usize;
        (1..=(want / one + 1) as u64)
            .map(|at| a_record(0, at, 1))
            .collect()
    }

    /// One block payload holding `records.len()` records, on a grid that never cuts.
    fn one_payload(records: &[SampleLocusObservations]) -> Vec<u8> {
        let payloads = cut(
            BlockBuilder::new(Bp(u64::MAX), None).expect("a grid"),
            records,
        )
        .expect("in order");
        assert_eq!(payloads.len(), 1, "the grid must not cut this fixture");
        payloads.into_iter().next().expect("one block")
    }

    #[test]
    fn a_block_round_trips_through_the_compressor_record_for_record() {
        let records: Vec<_> = (0..40).map(|index| a_record(0, 100 + index, 1)).collect();
        let payloads = cut(
            BlockBuilder::new(A_GRID, Some(200)).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert!(
            payloads.len() > 1,
            "more than one block, or this proves less"
        );

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let mut back = Vec::new();
        for payload in &payloads {
            let block = compressor.compress(payload).expect("a block compresses");
            let inflated =
                decompress_capped_at(frame_of(block), u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG))
                    .expect("it inflates under the window it declares");
            assert_eq!(&inflated, payload, "the payload came back byte for byte");
            back.extend(walk(&inflated).records);
        }
        assert_eq!(
            back, records,
            "and every record came back through the block"
        );
    }

    /// **The declared window is the window the frame needs, at the settings a writer ships
    /// with.** A decoder allowed exactly what the file declares inflates the block; one allowed
    /// a single exponent less refuses *the frame*, not the cap.
    ///
    /// Two things make this test able to fail, and it could not before either was here. The
    /// payload has to exceed the window, or zstd narrows the frame's own declaration to fit and
    /// the assertion says nothing. And the window has to be above zstd's minimum, or "one
    /// exponent less" is refused at the parameter and the frame is never looked at — which is
    /// what the first version of this test did, at a window of 2^10.
    #[test]
    fn the_shipped_settings_cap_a_payload_larger_than_the_window() {
        let payload = one_payload(&records_past_a_window(DEFAULT_LOOK_BACK_WINDOW_LOG));
        assert!(
            payload.len() > (1usize << DEFAULT_LOOK_BACK_WINDOW_LOG),
            "the payload must exceed the window; it is {} bytes against {}",
            payload.len(),
            1usize << DEFAULT_LOOK_BACK_WINDOW_LOG
        );

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let frame = frame_of(compressor.compress(&payload).expect("it compresses")).to_vec();

        assert!(
            decompress_capped_at(&frame, u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG)).is_ok(),
            "a decoder allowed exactly what the file declares must inflate it"
        );
        assert_eq!(
            decompress_capped_at(&frame, u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG) - 1),
            Err(NotInflated::TheFrameWasRefused),
            "and one exponent less must be refused by the frame, not by the cap"
        );
    }

    /// **⚠ A frame's own declared window is never *wider* than the file's, and is often
    /// narrower.** zstd narrows it to fit a payload smaller than the cap, so an ordinary block
    /// under a 32 kB declaration can write a frame needing 1 kB. Spec §4.2 calls a mismatch
    /// between our declaration and the frame's "a corruption worth detecting rather than
    /// tolerating" — **and the check that implements that has to be `≤`, not `=`**, or Milestone
    /// D3 rejects almost every block of every file this writer produces.
    #[test]
    fn a_frames_own_window_is_never_wider_than_the_file_declares() {
        let a_few = vec![a_record(0, 1, 1), a_record(0, 2, 1), a_record(0, 3, 1)];
        let small = one_payload(&a_few);
        let large = one_payload(&records_past_a_window(DEFAULT_LOOK_BACK_WINDOW_LOG));

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let narrowed = frame_of(compressor.compress(&small).expect("it compresses")).to_vec();
        let at_the_cap = frame_of(compressor.compress(&large).expect("it compresses")).to_vec();

        // The declared window is whatever the smallest cap that still inflates it is.
        let window_needed_by = |frame: &[u8]| -> u32 {
            (MIN_LOOK_BACK_WINDOW_LOG..=DEFAULT_LOOK_BACK_WINDOW_LOG)
                .map(u32::from)
                .find(|cap| decompress_capped_at(frame, *cap).is_ok())
                .expect("some cap at or under the file's inflates it")
        };

        assert!(
            window_needed_by(&narrowed) < u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG),
            "a small payload's frame declares less than the file does"
        );
        assert_eq!(
            window_needed_by(&at_the_cap),
            u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG),
            "and a payload past the cap declares exactly the file's"
        );
    }

    /// **The same payload gives the same bytes, whatever the compressor has seen before.** A
    /// writer keeps one compressor for the life of a file; if any state carried from one block
    /// to the next, a file's bytes would depend on how its blocks were scheduled, and goal 5's
    /// byte-identity check across worker counts would mean nothing.
    #[test]
    fn compressing_a_block_does_not_depend_on_the_blocks_before_it() {
        let records: Vec<_> = (0..60).map(|index| a_record(0, 100 + index, 1)).collect();
        let payloads = cut(
            BlockBuilder::new(A_GRID, Some(120)).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert!(
            payloads.len() >= 3,
            "several blocks, or history cannot show"
        );

        let mut long_lived =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let in_sequence: Vec<Vec<u8>> = payloads
            .iter()
            .map(|payload| {
                long_lived
                    .compress(payload)
                    .expect("it compresses")
                    .to_vec()
            })
            .collect();

        for (index, payload) in payloads.iter().enumerate() {
            let mut fresh =
                BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
            assert_eq!(
                fresh.compress(payload).expect("it compresses"),
                in_sequence[index].as_slice(),
                "block {index} compressed alone and in sequence"
            );
        }
    }

    /// A file's declared window reaches the compressor, so two files declaring different
    /// windows over the same records are not the same file.
    #[test]
    fn the_declared_window_reaches_the_bytes() {
        let payload = one_payload(&records_past_a_window(MIN_LOOK_BACK_WINDOW_LOG));

        let at = |window_log: u8| {
            BlockCompressor::new(window_log)
                .expect("a valid window")
                .compress(&payload)
                .expect("it compresses")
                .to_vec()
        };
        assert_ne!(
            at(MIN_LOOK_BACK_WINDOW_LOG),
            at(MIN_LOOK_BACK_WINDOW_LOG + 4),
            "a window four exponents wider must reach the bytes"
        );
    }

    /// **The level reaches the bytes too, and nothing else in the suite could see it.** Every
    /// other test that mentions the level compares two compressors at the *same* level, so a
    /// constructor that ignored its argument — or one that forced every caller to 9 — passed
    /// all of them.
    #[test]
    fn the_chosen_level_reaches_the_bytes() {
        let payload = one_payload(&records_past_a_window(DEFAULT_LOOK_BACK_WINDOW_LOG));
        let at = |level: i32| {
            let mut compressor = BlockCompressor::with_level(DEFAULT_LOOK_BACK_WINDOW_LOG, level)
                .expect("a valid level");
            assert_eq!(compressor.compression_level(), level);
            compressor.compress(&payload).expect("it compresses").len()
        };

        let cheap = at(1);
        let shipped = at(ZSTD_COMPRESSION_LEVEL);
        let dear = at(19);
        assert!(
            dear < shipped && shipped < cheap,
            "a higher level must give a smaller block: level 1 {cheap}, level \
             {ZSTD_COMPRESSION_LEVEL} {shipped}, level 19 {dear}"
        );
        assert_eq!(
            at(ZSTD_COMPRESSION_LEVEL),
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG)
                .expect("a window")
                .compress(&payload)
                .expect("it compresses")
                .len(),
            "and the constructor a writer uses writes at the shipped level"
        );
    }

    /// **A damaged frame is refused or harmless, never decoded into different records.** The
    /// frame checksum is what makes that true.
    ///
    /// **The claim is exactly that, and not "every damaged byte is refused".** Flipping a whole
    /// byte is refused every time; flipping a single *bit* is not — measured over one frame,
    /// nine flips of 552 were accepted and inflated to the payload unchanged, because they land
    /// in bits the format does not use. None inflated to anything different, which is the
    /// property that matters and the one asserted here.
    #[test]
    fn damage_to_a_frame_is_refused_or_harmless_and_never_silently_different() {
        let payload = one_payload(&(0..40).map(|i| a_record(0, 100 + i, 1)).collect::<Vec<_>>());
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let frame = frame_of(compressor.compress(&payload).expect("it compresses")).to_vec();

        let mut refused = 0usize;
        let mut harmless = 0usize;
        for at in 0..frame.len() {
            for bit in 0..8u32 {
                let mut damaged = frame.clone();
                damaged[at] ^= 1 << bit;
                match decompress_capped_at(&damaged, u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG)) {
                    Err(NotInflated::TheFrameWasRefused) => refused += 1,
                    Err(NotInflated::TheCapWasRefused) => {
                        panic!("the cap is this test's, not the frame's")
                    }
                    Ok(inflated) => {
                        assert_eq!(
                            inflated, payload,
                            "bit {bit} of byte {at} was accepted and gave different records"
                        );
                        harmless += 1;
                    }
                }
            }
        }
        assert_eq!(
            refused + harmless,
            frame.len() * 8,
            "every bit of the frame was tried"
        );
        assert!(
            refused > harmless * 10,
            "damage should overwhelmingly be refused; {refused} refused against {harmless} \
             harmless over {} bytes",
            frame.len()
        );
    }

    /// **A frame carries no uncompressed length, and that is deliberate.** Spec §8 lists a
    /// block header that says how large it inflates to among the traps — "a temptation to
    /// allocate it" — and the whole design is that a reader never sizes a buffer from a block.
    #[test]
    fn a_frame_says_nothing_about_how_large_it_inflates_to() {
        let payload = one_payload(&(0..40).map(|i| a_record(0, 100 + i, 1)).collect::<Vec<_>>());
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let block = compressor.compress(&payload).expect("it compresses");

        match zstd::zstd_safe::get_frame_content_size(frame_of(block)) {
            Ok(None) => {}
            other => panic!(
                "the frame must not declare the {} bytes it inflates to; it said {other:?}",
                payload.len()
            ),
        }
    }

    /// **A block zstd cannot shrink still compresses.** A block holding one small record is
    /// incompressible — zstd's framing and checksum cost more than they save — so the buffer a
    /// frame is written into has to be zstd's own worst-case bound and not the payload's
    /// length. Nothing reached this before: every compression fixture built forty records or
    /// more.
    #[test]
    fn a_payload_zstd_cannot_shrink_still_compresses() {
        let one_record = one_payload(&[a_record(0, 1, 1)]);
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let block = compressor.compress(&one_record).expect("it compresses");
        assert!(
            frame_of(block).len() > one_record.len(),
            "the fixture must be incompressible for this to say anything: {} payload bytes \
             against a {}-byte frame",
            one_record.len(),
            frame_of(block).len()
        );
        assert_eq!(
            decompress_capped_at(frame_of(block), u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG))
                .expect("it inflates"),
            one_record
        );

        // And an empty payload, which a block never holds but a caller can hand over.
        let block = compressor.compress(&[]).expect("nothing compresses too");
        assert!(
            decompress_capped_at(frame_of(block), u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG))
                .expect("it inflates")
                .is_empty()
        );
    }

    /// **The framing a file carries, pinned without pinning zstd's own bytes.** The length is
    /// four little-endian bytes and the frame's magic follows immediately; a round-trip test
    /// alone cannot see either, because a writer and a reader that both changed to big-endian
    /// would agree with each other and with nothing else.
    #[test]
    fn an_on_disk_block_is_a_little_endian_length_then_a_zstd_frame() {
        let payload = one_payload(&(0..40).map(|i| a_record(0, 100 + i, 1)).collect::<Vec<_>>());
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let block = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let frame_bytes = block.len() - COMPRESSED_BLOCK_LENGTH_BYTES;
        assert_eq!(
            &block[..COMPRESSED_BLOCK_LENGTH_BYTES],
            &u32::try_from(frame_bytes)
                .expect("a small frame")
                .to_le_bytes(),
            "the first four bytes are the frame's length, little-endian"
        );
        assert!(
            frame_bytes < 256,
            "the fixture's frame must be under 256 bytes, or a big-endian length would agree \
             with a little-endian one; it is {frame_bytes}"
        );
        assert_eq!(
            &block[COMPRESSED_BLOCK_LENGTH_BYTES..COMPRESSED_BLOCK_LENGTH_BYTES + 4],
            &[0x28, 0xb5, 0x2f, 0xfd],
            "and a zstd frame begins immediately after it"
        );
    }

    /// The length in front of a block says where the next one starts, and reading it touches
    /// nothing else — no decompression, no allocation.
    #[test]
    fn a_blocks_length_prefix_says_where_the_next_block_starts() {
        let records: Vec<_> = (0..40).map(|index| a_record(0, 100 + index, 1)).collect();
        let payloads = cut(
            BlockBuilder::new(A_GRID, Some(120)).expect("a grid"),
            &records,
        )
        .expect("in order");
        assert!(payloads.len() >= 3);

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let mut file = Vec::new();
        let mut lengths = Vec::new();
        for payload in &payloads {
            let block = compressor.compress(payload).expect("it compresses");
            lengths.push(block.len());
            file.extend_from_slice(block);
        }

        let mut at = 0usize;
        let mut walked = Vec::new();
        while at < file.len() {
            let CompressedBlockAt::Whole {
                zstd_frame,
                block_bytes,
            } = compressed_block_at(&file[at..])
            else {
                panic!("every block of a whole file is whole");
            };
            walked.push(block_bytes);
            assert_eq!(
                decompress_capped_at(zstd_frame, u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG))
                    .expect("it inflates"),
                payloads[walked.len() - 1]
            );
            at += block_bytes;
        }
        assert_eq!(walked, lengths, "every block was found by its own length");
        assert_eq!(at, file.len(), "and the walk consumed the file exactly");
    }

    /// **A declared length is not a fact, and the three states say which.** A reader pulling
    /// bytes from a file meets all three: too few for a length at all, a length whose block has
    /// not arrived, and a whole block. Before these were separate the function handed back the
    /// declaration alone, so a block truncated to half its length still reported its full one
    /// and the caller sliced on it.
    #[test]
    fn a_length_is_a_declaration_until_the_block_behind_it_is_there() {
        let payload = one_payload(&(0..40).map(|i| a_record(0, 100 + i, 1)).collect::<Vec<_>>());
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let whole = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        for cut_at in 0..COMPRESSED_BLOCK_LENGTH_BYTES {
            assert_eq!(
                compressed_block_at(&whole[..cut_at]),
                CompressedBlockAt::NoLengthYet,
                "at {cut_at} bytes there is not even a length"
            );
        }
        for cut_at in COMPRESSED_BLOCK_LENGTH_BYTES..whole.len() {
            assert_eq!(
                compressed_block_at(&whole[..cut_at]),
                CompressedBlockAt::NotAllHere {
                    block_bytes: whole.len()
                },
                "at {cut_at} bytes the length has arrived and the block has not"
            );
        }
        assert!(matches!(
            compressed_block_at(&whole),
            CompressedBlockAt::Whole { .. }
        ));

        // And a length no buffer could satisfy is `NotAllHere`, never a slice.
        let mut hostile = u32::MAX.to_le_bytes().to_vec();
        hostile.extend_from_slice(&[0u8; 4]);
        assert_eq!(
            compressed_block_at(&hostile),
            CompressedBlockAt::NotAllHere {
                block_bytes: COMPRESSED_BLOCK_LENGTH_BYTES + u32::MAX as usize
            }
        );
    }

    /// A look-back window zstd does not take is refused where the parameter is set, not left to
    /// zstd to report as a code. `header.rs` refuses the same values in a manifest.
    #[test]
    fn a_look_back_window_outside_what_zstd_takes_is_refused() {
        for window_log in [MIN_LOOK_BACK_WINDOW_LOG - 1, MAX_LOOK_BACK_WINDOW_LOG + 1] {
            match BlockCompressor::new(window_log) {
                Err(BlockCompressError::WindowLogOutOfRange {
                    look_back_window_log,
                }) => assert_eq!(look_back_window_log, window_log),
                other => panic!("2^{window_log} gave {other:?}"),
            }
        }
        assert!(BlockCompressor::new(MIN_LOOK_BACK_WINDOW_LOG).is_ok());
        assert!(BlockCompressor::new(MAX_LOOK_BACK_WINDOW_LOG).is_ok());
    }

    /// **A level zstd does not take is refused rather than clamped**, which is what zstd does
    /// silently — and silently is the problem: a level of 100 gives exactly the level-22 bytes,
    /// so a caller measuring the level's cost would read the same number twice and conclude the
    /// level had stopped mattering.
    #[test]
    fn a_compression_level_outside_what_zstd_takes_is_refused() {
        let allowed = zstd::compression_level_range();
        for level in [*allowed.start() - 1, *allowed.end() + 1, i32::MIN, i32::MAX] {
            match BlockCompressor::with_level(DEFAULT_LOOK_BACK_WINDOW_LOG, level) {
                Err(BlockCompressError::LevelOutOfRange {
                    level: refused,
                    lowest,
                    highest,
                }) => {
                    assert_eq!(refused, level);
                    assert_eq!((lowest, highest), (*allowed.start(), *allowed.end()));
                }
                other => panic!("level {level} gave {other:?}"),
            }
        }
        assert!(
            BlockCompressor::with_level(DEFAULT_LOOK_BACK_WINDOW_LOG, *allowed.start()).is_ok()
        );
        assert!(BlockCompressor::with_level(DEFAULT_LOOK_BACK_WINDOW_LOG, *allowed.end()).is_ok());
    }

    /// **A compressor made from a manifest caps at the window that manifest declares.** That is
    /// the only way a writer can be sure the frames it writes fit the promise its own header
    /// makes: a window taken from anywhere else gives a file every reader in the cohort refuses,
    /// block by block, with a zstd error code.
    #[test]
    fn a_compressor_from_a_manifest_caps_at_the_window_the_manifest_declares() {
        let mut manifest = a_manifest();
        manifest.look_back_window_log = MIN_LOOK_BACK_WINDOW_LOG + 2;
        let compressor = BlockCompressor::from_manifest(&manifest).expect("a valid window");
        assert_eq!(
            compressor.look_back_window_log(),
            manifest.look_back_window_log
        );
        assert_eq!(compressor.compression_level(), ZSTD_COMPRESSION_LEVEL);

        let payload = one_payload(&records_past_a_window(manifest.look_back_window_log));
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a valid window");
        let frame = frame_of(compressor.compress(&payload).expect("it compresses")).to_vec();
        assert!(
            decompress_capped_at(&frame, u32::from(manifest.look_back_window_log)).is_ok(),
            "a reader configured from the same manifest inflates it"
        );

        // And a manifest whose window is outside what zstd takes is refused here, not left to
        // produce a file with a window nobody can honour.
        manifest.look_back_window_log = MAX_LOOK_BACK_WINDOW_LOG + 1;
        assert!(matches!(
            BlockCompressor::from_manifest(&manifest),
            Err(BlockCompressError::WindowLogOutOfRange { .. })
        ));
    }

    /// The compressor a writer gets is the one every measurement was taken with.
    #[test]
    fn a_writers_compressor_uses_the_shipped_level() {
        assert_eq!(ZSTD_COMPRESSION_LEVEL, 9);
        let compressor = BlockCompressor::from_manifest(&a_manifest()).expect("a window");
        assert_eq!(
            compressor.look_back_window_log(),
            DEFAULT_LOOK_BACK_WINDOW_LOG
        );
        assert_eq!(compressor.compression_level(), ZSTD_COMPRESSION_LEVEL);
    }

    /// Every refusal names what it was doing and the number whoever sees it must act on.
    #[test]
    fn a_compressors_refusals_say_what_broke_and_what_the_bounds_are() {
        let window = BlockCompressor::new(MAX_LOOK_BACK_WINDOW_LOG + 1)
            .expect_err("a window zstd does not take");
        assert_eq!(
            window.to_string(),
            format!(
                "a look-back window of 2^{} bytes; zstd takes between 2^{MIN_LOOK_BACK_WINDOW_LOG} \
                 and 2^{MAX_LOOK_BACK_WINDOW_LOG}",
                MAX_LOOK_BACK_WINDOW_LOG + 1
            )
        );

        let allowed = zstd::compression_level_range();
        let level = BlockCompressor::with_level(DEFAULT_LOOK_BACK_WINDOW_LOG, *allowed.end() + 1)
            .expect_err("a level zstd does not take");
        assert_eq!(
            level.to_string(),
            format!(
                "a compression level of {}; zstd takes between {} and {}",
                *allowed.end() + 1,
                allowed.start(),
                allowed.end()
            )
        );

        let zstd_said = BlockCompressError::zstd(
            "compressing a block",
            std::io::Error::other("Destination buffer is too small"),
        );
        assert_eq!(
            zstd_said.to_string(),
            "zstd failed while compressing a block"
        );
        assert_eq!(
            std::error::Error::source(&zstd_said)
                .expect("zstd's own account is the cause")
                .to_string(),
            "Destination buffer is too small"
        );
    }

    // -----------------------------------------------------------------
    // The two round-trip laws, over the whole value range
    // -----------------------------------------------------------------

    proptest! {
        /// Any block head reads back as itself and says exactly how many bytes it took.
        #[test]
        fn a_block_head_round_trips_for_any_contig_position_and_count(
            contig in proptest::num::u32::ANY,
            first_position in proptest::num::u64::ANY,
            record_count in 1u64..=u64::MAX,
        ) {
            let head = a_head(contig, first_position, record_count);
            let mut bytes = Vec::new();
            head.encode(&mut bytes);
            let read = BlockHead::decode(&bytes).expect("it reads back");
            prop_assert_eq!(read.head, head);
            prop_assert_eq!(read.head_bytes, bytes.len());
        }

        /// **Records in, the same records out — whatever the grid and the ceiling are** — and
        /// every block holds one contig's worth of one grid cell. The spans are drawn wide
        /// enough to straddle a multiple against the smallest grids, which is the class the
        /// hand-written fixtures reach at only one coordinate each.
        #[test]
        fn every_record_comes_back_unchanged_for_any_grid_and_ceiling(
            grid in 1u64..=1_000u64,
            ceiling in proptest::option::of(1u32..=200u32),
            steps in proptest::collection::vec((0u64..=300u64, 1u64..=40u64), 1..40),
        ) {
            let mut start = 1u64;
            let mut records = Vec::new();
            for (step, span) in steps {
                start += step;
                records.push(a_record(0, start, span));
            }

            let blocks = cut(
                BlockBuilder::new(Bp(grid), ceiling).expect("a grid"),
                &records,
            )
            .expect("the starts never go backwards");

            for block in &blocks {
                let walked = walk(block);
                prop_assert_eq!(walked.head.record_count.get(), walked.records.len() as u64);
                for record in &walked.records {
                    prop_assert_eq!(
                        record.region.start.get() / grid,
                        walked.head.first_position.get() / grid,
                        "every record in a block falls in the block's own grid cell"
                    );
                }
            }
            prop_assert_eq!(walk_all(&blocks), records);
        }
    }
}
