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
//! **Three halves, and the file is read in this order.** [`BlockBuilder`] turns a stream of
//! records into block payloads, cutting them on the grid. [`BlockCompressor`] turns one payload
//! into a whole on-disk block — a four-byte length and then one zstd frame — with the
//! compressor's look-back window capped at what the file declares, which is what unties the
//! block's size from the reader's memory. [`BlockStream`] reads a run of them back a record at a
//! time, holding two 16 kB buffers and nothing that grows with the block.
//!
//! **Nothing here writes or opens a file.** Milestone F owns the header, the index and the
//! footer; this module is handed bytes and hands bytes back, which is what lets a caller start a
//! reader at any block it has an offset for.
//!
//! Design authority: `doc/devel/ng/spec/psp_file_format.md` §3.2 (a block is self-contained),
//! §4.1 (the cut rule and what it buys), §4.2 (the window is declared and the reader honours
//! it), §4.4 (the reader's two buffers), §5.1 (the reader's contract), §8 (the traps), and
//! `doc/devel/ng/arch/psp_file_format.md` §1.
//!
//! [`RecordHead`]: crate::ng::psp::RecordHead

use std::num::NonZeroU64;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::psp::chain_ids::{LiveSet, LiveSetReader};
use crate::ng::psp::header::{MAX_LOOK_BACK_WINDOW_LOG, MIN_LOOK_BACK_WINDOW_LOG, Manifest};
use crate::ng::psp::record::{
    OffsetBase, RecordDecodeError, RecordEncodeError, RecordEncoder, RecordHead, RecordLayout,
    RecordLayoutError, decode_the_body_of, read_record_head,
};
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
/// **A grid cell holding no records produces no block.** The rule above decides where a block
/// *ends*; it never asks for one per 100 kb of reference. A sample covering two cells ninety
/// apart writes two blocks, so a thin sample pays no index entry and no compressed frame for
/// reference it did not cover.
///
/// **And blocks that are merely *small* are not merged, which is a ruling rather than an
/// omission.** Spec §4.1 offers a second secondary rule — accumulate across empty stretches so a
/// patchy sample gets one large block instead of several thin ones, each compressing from cold —
/// and spec §12 question 3 says it ships, leaving only its threshold open. **The owner ruled
/// against it on 2026-08-27: merging would complicate the alignment between samples.** That is
/// the property the grid exists for: every sample's boundaries fall on the same coordinates, so
/// a cohort reader stepping across a region knows which block of each sample holds a position.
/// Merge, and one sample's block may begin ninety cells earlier than its neighbour's, so the
/// block holding a given position differs from sample to sample and a reader wanting one
/// position decodes from far behind it.
///
/// **⚠ What merging would have saved was never measured, here or in the spec.** An earlier
/// version of this comment priced it at "about 7 %" and cited spec §4.1's 17.557 bytes a record
/// against 16.444 — but those two are the 100 kb and 1,000 kb rows of the *same*, non-merging
/// cut rule, so they price a change of block size and say nothing about merging. No writer that
/// accumulates across empty stretches has been built or measured. The ruling stands on the
/// alignment argument above, which does not need a number; if someone later wants the price,
/// they have to build the merging writer to get it.
///
/// **Spec §4.1 and §12 question 3 record the ruling** as of 2026-08-27; before that they said the
/// rule shipped and only its threshold was open, which is what this comment used to flag.
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
            // The file's first record. The encoder resets nothing until the record is known to
            // be writable, so a record the codec refuses leaves no block open behind it.
            if let Err(refused) = self
                .encoder
                .encode_record_starting_a_block(record, &mut self.records)
            {
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
        //
        // **And nothing here puts the encoder back.** `encode_record_starting_a_block` makes
        // every refusal before it resets anything, so a refused record leaves the open block's
        // running differences exactly as they were. ⚠ There used to be a `start_block(resume_at)`
        // on this path, which restored the coordinate base and would have thrown away the live
        // set that the open block still needed — a reset is not undoable once a second difference
        // rides on it.
        self.next_records.clear();
        if let Err(refused) = self
            .encoder
            .encode_record_starting_a_block(record, &mut self.next_records)
        {
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

// ---------------------------------------------------------------------
// Reading a block back, a record at a time
// ---------------------------------------------------------------------

/// How many compressed bytes a reader pulls from the file at a time.
///
/// **16 kB, and it is an optimum rather than a floor** (spec §4.4). **⚠ The numbers below are
/// the measuring prototype's reader, not this one** — nothing has yet timed `BlockStream`, and
/// Milestone H is where that happens. Re-measured there on a quiet machine over a walk of
/// 7.69 M records: 4 kB takes 0.149 s and holds 233 kB an open sample,
/// 16 kB takes 0.143 s and holds 257 kB, 64 kB takes 0.161 s and holds 353 kB, and 256 kB takes
/// 0.200 s. Going up costs both time and memory; going down costs a little time and saves a
/// little memory. Why the curve turns is not established.
///
/// **It is the reader's choice and not the file's**, which is the point: nothing a reader holds
/// is a function of the block size, and untying those two is what this format exists for.
pub const READ_CHUNK_BYTES: usize = 16 * 1024;

/// How many decompressed bytes a reader keeps in front of the parser. See [`READ_CHUNK_BYTES`]
/// for the measurement; the two were swept together.
///
/// **A record larger than this is not an error**, and a fixed maximum record size is not a safe
/// assumption to bake in (spec §8): the buffer grows for one, and shrinks back to this at the
/// next block.
pub const ROLLING_BYTES: usize = 16 * 1024;

/// How far a reader's rolling buffer may grow for one record before it refuses the block.
///
/// **A reader's budget, not a maximum record size the format fixes** — and the name says buffer
/// rather than record for that reason. Spec §8 refuses the second outright ("many alleles, many
/// chain ids"), because a ceiling baked into the *format* makes a legitimate file unreadable
/// everywhere at once. But spec §1.1 puts an open sample at 500 kB, and on a **corrupt** block the
/// two do not both hold: no record parses, so nothing ever fits, and the buffer doubles until the
/// block's whole decompressed size is in it. That size is not bounded by the block's size on disk
/// — measured, **4,132 bytes on disk drove a reader to hold 67,125,248**.
///
/// **Why half a mebibyte.** Two numbers, both measured rather than estimated:
///
/// - **It is 10× the largest record this caller's own depth cap can produce.** Built to the top of
///   the committed range — three hundred reads a position, one observation each, with their
///   moments — a record encodes to **18,292 bytes** at a 50-base span and **48,693** at 150 bases.
///   (An earlier version of this comment said "roughly 30 kB" from estimation; the measured figure
///   is the first of those. Milestone E's chain ids add about 8 bytes a read on top, since
///   `encode_record_body` drops them today.)
/// - **It is what keeps the worst case inside the budget spec §1.1 states.** A caller holds one
///   reader per sample at up to several thousand; three thousand samples every one of which met a
///   damaged block at the same moment is **1,572,864,000 bytes**, against the 1.5 GB §1.1 gives
///   three thousand open samples — the same number to within 5 %. At 1 MiB — the value first
///   proposed — it would have been 3.07 GB, twice the whole budget, which is a weak bound for
///   something whose purpose is bounding hostile input. The arithmetic is a test, so it moves
///   when the constant does: `the_buffer_ceiling_is_priced_at_the_top_of_the_committed_cohort_range`.
///
/// **And it is raisable at run time**, which is the other half of why a ceiling is tolerable here:
/// [`BlockStream::with_a_buffer_ceiling`] takes one, and
/// [`BlockReadError::RecordLargerThanTheReaderAllows`] names it. That is the pattern spec §4.2
/// uses for a look-back window wider than a reader budgeted for — *the fix is a knob rather than a
/// rebuild*. **⚠ It was a `const` read straight out of the reader when this was first written, so
/// "raisable" meant recompiling**; the knob is what makes the sentence true.
///
/// Spec §4.4 carries the same two numbers, §7 the refusal's row, and §8 the trap they answer.
pub const ROLLING_BUFFER_CEILING_BYTES: usize = 512 * 1024;

const _: () = assert!(
    ROLLING_BUFFER_CEILING_BYTES > ROLLING_BYTES,
    "a ceiling at or under the buffer would refuse records the buffer already holds"
);

const _: () = assert!(
    ROLLING_BUFFER_CEILING_BYTES >= 64 * 1024,
    "a ceiling under 64 kB refuses records this caller's own depth cap produces: 48,693 bytes \
     at a 150-base span, and Milestone E's chain ids are still to come"
);

/// One record as a walk met it.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct StreamedRecord {
    /// Which block it came from — its contig, its first position, and how many records it holds.
    pub block: BlockHead,
    /// What the record's own head said, whether or not the body was built.
    pub head: RecordHead,
    /// The record, when the caller wanted it. **`None` is the skip**, and it is what the head
    /// exists for: the body was never decoded, only advanced past.
    pub record: Option<SampleLocusObservations>,
}

/// Reads records back out of a run of compressed blocks, holding nothing that grows with them.
///
/// **Two conditions make that true and only the first is the compressor's doing** (spec §5.1).
/// One: do not inflate a whole block — pull decompressed bytes out incrementally, bounded by
/// the declared window. Two: **do not accumulate what you inflated** — parse a record out of the
/// rolling buffer, hand it to the caller, keep nothing. Satisfying only the first moves the
/// memory rather than removing it: in production's cohort run the assembled per-sample columns
/// are the largest single mass of the heap, larger than the decompression buffers they came
/// from. This type hands each record over and retains only its two buffers, the decoder's own
/// state, and where it is.
///
/// **It reads forward from wherever `source` is.** Milestone F's index turns a coordinate into a
/// byte offset and seeks before handing the source over; this type knows nothing about files,
/// which is what lets a caller start one at an arbitrary block.
pub struct BlockStream<R> {
    source: R,
    decoder: zstd::zstd_safe::DCtx<'static>,
    look_back_window_log: u8,
    layout: RecordLayout,

    // ---- what lives for the whole stream: storage, budgets, and totals ----
    /// Compressed bytes pulled from the source, and how much of them has been fed to the
    /// decoder.
    compressed: Vec<u8>,
    compressed_at: usize,
    compressed_filled: usize,
    /// Decompressed bytes the parser works out of. **Storage, not state** — how much of it has
    /// been consumed is in the cursor, because that is per-block.
    rolling: Vec<u8>,

    /// How many times a record's parse ran out of bytes and was started again — **counted over
    /// the whole file and never reset at a block**, because a total is what the oracle compares.
    ///
    /// ⚠ The increment is at the top of the truncated arm, *before* the refill is asked for, so a
    /// stream that then refuses with `RecordRunsPastItsBlock` has counted one restart that never
    /// happened. That is why the name is not `retries_after_a_refill`, which it was: at most one
    /// per refused stream, and nothing pins the distinction.
    ///
    /// **Kept unconditionally, at eight bytes against a 500 kB budget.** It is the one thing
    /// about a restartable parse that cannot be seen from outside — a walk that comes out right
    /// is silent about whether the retry ever ran — and it was `#[cfg(test)]` until a review
    /// pointed out that this hides it from `tests/`, from benches, and from Milestone H, which
    /// is the one place this reader is to be timed.
    parses_restarted: u64,

    /// How many times a *block head*'s parse has been restarted after asking for more bytes.
    ///
    /// Counted separately because it answers a different question, and the answer is
    /// surprising: a block head arrives whole or not at all. zstd emits an internal block in
    /// one piece, so a refill lands *before* a head and never inside one — measured, over
    /// 108,746 head restarts the largest partial head ever seen was **zero bytes**.
    block_heads_restarted: u64,

    /// How far `rolling` may grow for one record. Defaults to
    /// [`ROLLING_BUFFER_CEILING_BYTES`]; [`BlockStream::with_a_buffer_ceiling`] sets it.
    buffer_ceiling: usize,

    /// Which reads are live, and the changes read out of each record's head to move between
    /// records.
    ///
    /// **Its own per-block state is inside it**, emptied by [`LiveSetReader::start_block`], which
    /// `begin_next_block` calls. It sits here rather than in [`BlockCursor`] because it owns
    /// scratch buffers that must survive a boundary; what must not survive one is the set.
    live_reads: LiveSetReader,

    /// Set once this reader has refused, and never cleared. **A stream that has refused is
    /// finished**: without a state that says so, what stopped a refused reader was only that its
    /// compressed buffer had been emptied — which lasts exactly as long as the whole file fits in
    /// one read. Measured on a 219,758-byte file of 150 blocks: after the first refusal the
    /// reader resynchronised at the next read-chunk boundary, took four arbitrary bytes for a
    /// block length, and carried on; on a smaller file it refused once and then handed back
    /// 5,681 further records and ended cleanly. Milestone F hands this reader a seeked `File`,
    /// where that is the ordinary case rather than the exotic one.
    refused: bool,

    // ---- what lives for one block, and is replaced whole at every boundary ----
    //
    // **Nothing per-block goes beside `cursor`.** A field here is initialised once in
    // `new` and never reset; the same field inside [`BlockCursor`] cannot be, because
    // `opening` and `between_blocks` are its only constructors and both literals are
    // exhaustive. The hand-written `Debug` below does make a field here a compile error —
    // but the fix it asks for is `field: _`, which is not the fix.
    cursor: BlockCursor,
}

/// Everything a reader **restarts when a block does** (spec §3.2).
///
/// **A struct rather than loose fields, for the reason D1 gave `RecordEncoder`'s own per-block
/// state**: spec §3.2 requires every running difference inside a block to restart, Milestone E
/// adds the chain-id live set, and a field added beside these and initialised once is one that
/// silently never resets — a file that then reads back wrong from every block's first record,
/// and plausibly wrong, because coverage is smooth. It is rebuilt by one assignment from
/// [`opening`](Self::opening), whose literal is exhaustive, so a field added here is a compile
/// error at the reset.
///
/// **The two buffers are not in here**, and that is the distinction: they are storage a reader
/// reuses across blocks, not state a record is parsed against.
///
/// **⚠ Not `Copy`, deliberately.** Milestone E's chain-id live set is a collection, and a `Copy`
/// struct cannot hold one — `error[E0204]` — so a coder meeting that would put the field on
/// [`BlockStream`] instead, where it compiles, passes every test, and is never reset. Measured: a
/// live set added there builds clean and all 197 `ng::psp` tests pass. Dropping `Copy` is what
/// keeps the field's natural home the one that forces the reset.
#[derive(Debug, Clone)]
struct BlockCursor {
    /// How much of the rolling buffer the parser has consumed. **Per-block, and it belongs
    /// here**: dropping its reset fails eight tests, so it is parse state rather than the
    /// buffer's own bookkeeping.
    rolling_at: usize,
    /// `None` until the block's head has been read out of the rolling buffer.
    block: Option<BlockHead>,
    /// How many records of this block are still to come.
    records_left: u64,
    /// The coordinate the next record's position offset is measured from.
    measured_from: OffsetBase,
    /// Compressed bytes of this block that have not reached the decoder yet.
    compressed_left: usize,
    /// True once this block's frame has been decompressed to its end.
    inflated: bool,
}

impl BlockCursor {
    /// A block whose compressed bytes have been counted and whose head has not been read.
    fn opening(compressed_left: usize) -> Self {
        Self {
            rolling_at: 0,
            block: None,
            records_left: 0,
            measured_from: OffsetBase::at_block_start(Position(0)),
            compressed_left,
            inflated: false,
        }
    }

    /// Before any block has been opened, and after one has been abandoned: nothing due, nothing
    /// left to inflate.
    ///
    /// **Written out rather than `..Self::opening(0)`**, so a field added to this type is a
    /// compile error at *both* constructors. With the update syntax it was an error at
    /// `opening` alone and this one silently inherited whatever `opening` chose — which for a
    /// field that must differ between "not started" and "started" is exactly the wrong default.
    fn between_blocks() -> Self {
        Self {
            rolling_at: 0,
            block: None,
            records_left: 0,
            measured_from: OffsetBase::at_block_start(Position(0)),
            compressed_left: 0,
            inflated: true,
        }
    }
}

impl<R> std::fmt::Debug for BlockStream<R> {
    /// `zstd::zstd_safe::DCtx` is not `Debug`, so this is written by hand — and it destructures
    /// `Self` with no `..`, so a field added to the reader is a compile error here rather than
    /// one that silently stops being reportable. **Buffer lengths, never buffer contents.**
    ///
    /// **⚠ Being named here is not the same as being reset.** This is the *reporting* guard, and
    /// the compile error it raises is answered by one line — `new_field: _`. A per-block field
    /// added beside `cursor` therefore compiles clean after that one line and is never reset at a
    /// block boundary; the guard that catches *that* is `BlockCursor`'s two exhaustive
    /// constructors, and the field has to be inside it to meet them.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            source: _,
            decoder: _,
            look_back_window_log,
            layout,
            compressed: _,
            compressed_at,
            compressed_filled,
            rolling,
            cursor,
            live_reads,
            parses_restarted,
            block_heads_restarted,
            buffer_ceiling,
            refused,
        } = self;
        f.debug_struct("BlockStream")
            .field("look_back_window_log", look_back_window_log)
            .field("unknown_declared_fields", &layout.unknown_field_count())
            .field("compressed_buffered", &(compressed_filled - compressed_at))
            .field("rolling_buffered", &(rolling.len() - cursor.rolling_at))
            .field("rolling_capacity", &rolling.capacity())
            .field("cursor", cursor)
            .field("live_reads", &live_reads.live().len())
            .field("parses_restarted", parses_restarted)
            .field("block_heads_restarted", block_heads_restarted)
            .field("buffer_ceiling", buffer_ceiling)
            .field("refused", refused)
            .finish_non_exhaustive()
    }
}

impl<R: std::io::Read> BlockStream<R> {
    /// A reader for a file whose manifest is `manifest`, reading forward from wherever `source`
    /// is positioned — which must be the start of a block.
    ///
    /// **The decoder's window comes from the manifest and nothing else.** zstd refuses a frame
    /// whose window exceeds what its decoder was configured for, so a reader that assumed a
    /// value would reject a perfectly good file with an error naming a zstd code (spec §4.2).
    pub fn new(source: R, manifest: &Manifest) -> Result<Self, BlockReadError> {
        // Destructured with no `..`: a manifest field a reader must honour is a compile error
        // here rather than a setting it silently ignores.
        let Manifest {
            look_back_window_log,
            // `BlockBuilder`'s: they decide where a *writer* ends a block, and a reader is told
            // where each one ends by the block's own framing.
            genomic_block_size_bp: _,
            block_byte_ceiling: _,
            // Read below, into the layout every record is decoded against.
            fields: _,
        } = manifest;
        let look_back_window_log = *look_back_window_log;
        if !(MIN_LOOK_BACK_WINDOW_LOG..=MAX_LOOK_BACK_WINDOW_LOG).contains(&look_back_window_log) {
            return Err(BlockReadError::WindowLogOutOfRange {
                look_back_window_log,
            });
        }
        let layout = RecordLayout::from_manifest(manifest)
            .map_err(|source| BlockReadError::UnsupportedRecordLayout { source })?;

        let mut decoder = zstd::zstd_safe::DCtx::create();
        decoder
            .set_parameter(zstd::zstd_safe::DParameter::WindowLogMax(u32::from(
                look_back_window_log,
            )))
            .map_err(|code| BlockReadError::zstd("configuring the decoder's window", code))?;

        Ok(Self {
            source,
            decoder,
            look_back_window_log,
            layout,
            compressed: vec![0u8; READ_CHUNK_BYTES],
            compressed_at: 0,
            compressed_filled: 0,
            rolling: Vec::with_capacity(ROLLING_BYTES),
            cursor: BlockCursor::between_blocks(),
            live_reads: LiveSetReader::new(),
            parses_restarted: 0,
            block_heads_restarted: 0,
            buffer_ceiling: ROLLING_BUFFER_CEILING_BYTES,
            refused: false,
        })
    }

    /// How many times a record's parse has been restarted after asking for more bytes — see
    /// the field for why a walk that comes out right does not answer this.
    /// The reads live at the record last handed back.
    ///
    /// **Exact after a skipped record too**, which is the whole point of putting the changes in
    /// the head: a caller walking with a predicate that declines most records still knows which
    /// reads are live at the ones it takes.
    #[must_use]
    pub fn live_reads(&self) -> &LiveSet {
        self.live_reads.live()
    }

    pub fn parses_restarted(&self) -> u64 {
        self.parses_restarted
    }

    /// How many times a block head's parse has been restarted after asking for more bytes.
    pub fn block_heads_restarted(&self) -> u64 {
        self.block_heads_restarted
    }

    /// The same reader, with a different ceiling on how far its rolling buffer may grow for one
    /// record.
    ///
    /// **This is the knob [`ROLLING_BUFFER_CEILING_BYTES`] is the default of.** A genuine record
    /// larger than the default makes a file unreadable until someone raises it, which is why the
    /// refusal names the number — the same shape spec §4.2 gives a look-back window wider than a
    /// reader budgeted for, where the fix is a setting rather than a rebuilt file.
    ///
    /// A ceiling at or under [`ROLLING_BYTES`] is refused: it would turn away records the buffer
    /// already holds without growing at all.
    pub fn with_a_buffer_ceiling(mut self, ceiling: usize) -> Result<Self, BlockReadError> {
        if ceiling <= ROLLING_BYTES {
            return Err(BlockReadError::BufferCeilingUnderTheBuffer {
                ceiling,
                buffer_bytes: ROLLING_BYTES,
            });
        }
        self.buffer_ceiling = ceiling;
        Ok(self)
    }

    /// How far this reader's rolling buffer may grow for one record.
    pub fn buffer_ceiling(&self) -> usize {
        self.buffer_ceiling
    }

    /// The look-back window this reader's decoder is configured for, as the file declared it.
    pub fn look_back_window_log(&self) -> u8 {
        self.look_back_window_log
    }

    /// How many bytes this reader's two buffers are holding right now.
    ///
    /// **This is the number the whole format is shaped around.** A caller holds one of these per
    /// sample for the length of a run, so what one costs is multiplied by the cohort size, and
    /// goal 1 puts the budget at 500 kB an open sample. What it does *not* include is zstd's own
    /// context — about 190 kB, and no buffer choice reaches it (spec §5.3) — so this is the part
    /// a reader controls.
    ///
    /// **It must not grow with the block, the depth, or the length of the genome.** It grows for
    /// exactly one thing: a record larger than [`ROLLING_BYTES`], and it goes back at the next
    /// block.
    pub fn buffered_bytes(&self) -> usize {
        self.compressed.capacity() + self.rolling.capacity()
    }

    /// The next record, built — the walk a caller takes when it wants everything.
    pub fn next_record(&mut self) -> Option<Result<StreamedRecord, BlockReadError>> {
        self.next_record_where(|_| true)
    }

    /// The next record, built only if `want` says so.
    ///
    /// **This is the whole of the skip, and it is the reader's decision rather than a separate
    /// call** (spec §6.2). A record the predicate declines costs its head and a pointer advance;
    /// its body is never decoded. **⚠ What that is worth was measured on the prototype's reader
    /// and not on this one**: over 7.69 M records of a tomato accession, a walk keeping one
    /// record in a hundred took 0.141 s against 0.29 s for one building every record. Milestone
    /// H times this reader.
    ///
    /// **The bytes still have to arrive either way.** Skipping saves building the record, not
    /// decompressing it: a block comes out of zstd sequentially and there is nothing to seek
    /// past.
    pub fn next_record_where<F>(
        &mut self,
        mut want: F,
    ) -> Option<Result<StreamedRecord, BlockReadError>>
    where
        F: FnMut(&RecordHead) -> bool,
    {
        if self.refused {
            return None;
        }
        loop {
            if self.cursor.records_left == 0 {
                match self.begin_next_block() {
                    Ok(true) => {}
                    Ok(false) => return None,
                    Err(refused) => return Some(Err(self.fail(refused))),
                }
            }
            let block = self
                .cursor
                .block
                .expect("a block is open once its head has been read");

            match read_record_head(
                &self.rolling[self.cursor.rolling_at..],
                block.contig,
                self.cursor.measured_from,
                &mut self.live_reads,
            ) {
                Ok(found) => {
                    let head = found.head;
                    let record_bytes = found.record_bytes;
                    // **The head is read once, and the body is built from what it located.**
                    // Reading it applied this record's chain-id changes, so parsing it a second
                    // time to reach the body would apply them twice.
                    let record = if want(&head) {
                        match decode_the_body_of(&found, self.live_reads.live(), &self.layout) {
                            Ok(decoded) => Some(decoded.record),
                            Err(refused) => {
                                return Some(Err(self.fail(BlockReadError::from_record(refused))));
                            }
                        }
                    } else {
                        None
                    };
                    self.cursor.rolling_at += record_bytes;
                    self.cursor.measured_from = OffsetBase::after(&head);
                    self.cursor.records_left -= 1;
                    if self.cursor.records_left == 0
                        && let Err(refused) = self.check_the_block_ended_here()
                    {
                        return Some(Err(self.fail(refused)));
                    }
                    return Some(Ok(StreamedRecord {
                        block,
                        head,
                        record,
                    }));
                }
                Err(RecordDecodeError::Truncated { .. }) => {
                    self.parses_restarted += 1;
                    // **Nothing but `pump` may go here.** The parse starts again from the
                    // record's first byte, against state this arm has not touched — which is
                    // what makes it restartable rather than half-advanced. Spec §8: "a parse
                    // that half-advances that state before failing corrupts every record after
                    // it, plausibly."
                    let refilled = self.pump();
                    match refilled {
                        Ok(true) => {}
                        Ok(false) => {
                            return Some(Err(self.fail(BlockReadError::RecordRunsPastItsBlock {
                                records_left: self.cursor.records_left,
                                bytes_left: self.rolling.len() - self.cursor.rolling_at,
                            })));
                        }
                        Err(refused) => return Some(Err(self.fail(refused))),
                    }
                }
                Err(damage) => return Some(Err(self.fail(BlockReadError::from_record(damage)))),
            }
        }
    }

    /// Read the next block's framing and its head, or report that the source is finished.
    ///
    /// **Everything a record is parsed against resets here** — the coordinate the position
    /// offsets are measured from, and how many records are still to come (spec §3.2). A block
    /// that carried state in from the one before it would read back wrong from its first record,
    /// and plausibly wrong, because coverage is smooth.
    fn begin_next_block(&mut self) -> Result<bool, BlockReadError> {
        let mut declared = [0u8; COMPRESSED_BLOCK_LENGTH_BYTES];
        if !self.read_exactly(&mut declared)? {
            return Ok(false);
        }
        // **One assignment, and it is the reset.** Everything a record is parsed against comes
        // from here (spec §3.2); a block that carried state in from the one before it would read
        // back wrong from its first record, plausibly, because coverage is smooth.
        self.cursor = BlockCursor::opening(u32::from_le_bytes(declared) as usize);
        // **The second running difference, and it lives inside its own type.** The set of reads
        // live restarts here too, which is what lets this reader begin at an arbitrary block:
        // with nothing live, the block's first record restates its whole set as arrivals.
        self.live_reads.start_block();
        // **Defensive, and measured to be so.** zstd takes consecutive frames on one context
        // without being told, so removing this line passes every test in the module — the only
        // states it would rescue are ones where the previous frame ended mid-way, and those
        // already end the stream (`fail`). It is here because a block boundary is where a
        // reader's state resets (spec §3.2), and leaving the decoder's out of that reset would
        // make the exception the thing a later coder has to remember.
        self.decoder
            .reset(zstd::zstd_safe::ResetDirective::SessionOnly)
            .map_err(|code| {
                BlockReadError::zstd("resetting the decoder for the next block", code)
            })?;

        self.rolling.clear();
        // Back to what a reader budgets for, so one enormous record does not leave every later
        // block holding its buffer: nothing a reader holds may be a function of the data.
        self.rolling.shrink_to(ROLLING_BYTES);

        loop {
            match BlockHead::decode(&self.rolling[self.cursor.rolling_at..]) {
                Ok(decoded) => {
                    self.cursor.rolling_at += decoded.head_bytes;
                    self.cursor.measured_from =
                        OffsetBase::at_block_start(decoded.head.first_position);
                    self.cursor.records_left = decoded.head.record_count.get();
                    self.cursor.block = Some(decoded.head);
                    return Ok(true);
                }
                Err(BlockHeadDecodeError::Truncated { field, bytes_in }) => {
                    self.block_heads_restarted += 1;
                    if !self.pump()? {
                        return Err(BlockReadError::BlockHeadRunsPastItsBlock { field, bytes_in });
                    }
                }
                Err(damage) => return Err(BlockReadError::DamagedBlockHead { source: damage }),
            }
        }
    }

    /// The block declared its last record, so nothing of it may be left — and the three ways
    /// something can be are three different faults.
    ///
    /// **They were one variant reporting `bytes_left: 0` for two of them**, which told whoever
    /// read the message that a block held nothing past its records and was refused anyway.
    fn check_the_block_ended_here(&mut self) -> Result<(), BlockReadError> {
        // The frame may not have been decompressed to its end; pump until it has.
        while self.cursor.rolling_at >= self.rolling.len() && self.pump()? {}
        let bytes_left = self.rolling.len() - self.cursor.rolling_at;
        if bytes_left > 0 {
            return Err(BlockReadError::BlockHoldsMoreThanItDeclared { bytes_left });
        }
        if self.cursor.compressed_left > 0 {
            // The length in front of the block covers more than the frame behind it — on a run
            // of blocks that reaches into the *next* one, and a reader that believed it would
            // hand back records from a block nobody asked for.
            return Err(BlockReadError::BlockDeclaresMoreBytesThanItsFrame {
                compressed_bytes_left: self.cursor.compressed_left,
            });
        }
        if !self.cursor.inflated {
            return Err(BlockReadError::BlockFrameDidNotEndWithItsRecords);
        }
        Ok(())
    }

    /// Decompress more of the block into the rolling buffer. `false` when the block's frame is
    /// finished and nothing more will arrive.
    fn pump(&mut self) -> Result<bool, BlockReadError> {
        if self.cursor.inflated {
            return Ok(false);
        }
        // Drop what the parser has already consumed before asking for more. **This is what keeps
        // the rolling buffer rolling**: without it it grows to the whole block.
        if self.cursor.rolling_at > 0 {
            self.rolling.drain(..self.cursor.rolling_at);
            self.cursor.rolling_at = 0;
        }
        // **The ceiling counts one record, and this is where that is true.** The drain above
        // leaves `rolling` holding exactly what the record in front of the parser has claimed so
        // far, so `len()` and "one record" are the same number — which is what lets the refusal
        // below name a record. `pump`'s other callers keep it true the same way:
        // `begin_next_block` clears the buffer before its own retry loop, and
        // `check_the_block_ended_here` pumps only once the parser has consumed everything.
        debug_assert_eq!(
            self.cursor.rolling_at, 0,
            "the ceiling counts one record, so nothing already consumed may still be held"
        );
        // ⚠ **`>=` rather than `>` is a one-byte convention, and no test distinguishes them.**
        // Measured: swapping it fails none of the 201 `ng::psp` tests. What it decides is whether
        // the buffer may ever *hold* exactly `buffer_ceiling` bytes or must stop one below, and
        // building an oracle for that would mean searching for the exact rolling length at the
        // deciding pump — a fixture pinned to zstd's emission sizes, to prove one byte on a
        // 512 kB budget. Recorded as an unpinned convention rather than given a test that would
        // break whenever the compressor's internals moved.
        if self.rolling.len() >= self.buffer_ceiling {
            // The buffer has grown as far as this reader allows and the record in front of it
            // still does not fit. **Refused rather than grown into**, so that a block nothing
            // can be parsed out of costs a bounded amount of memory rather than its own whole
            // decompressed size.
            let block = self.cursor.block.map_or(
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(0),
                    end: Position(0),
                },
                |head| GenomeRegion {
                    contig: head.contig,
                    start: head.first_position,
                    end: head.first_position,
                },
            );
            let records_read = self
                .cursor
                .block
                .map_or(0, |head| head.record_count.get() - self.cursor.records_left);
            return Err(BlockReadError::RecordLargerThanTheReaderAllows {
                block,
                records_read,
                allowed_bytes: self.buffer_ceiling,
            });
        }
        if self.rolling.len() >= self.rolling.capacity() {
            // One record needs more than the buffer holds. **Grow rather than fail**: spec §8
            // refuses a fixed maximum record size outright — "many alleles, many chain ids" —
            // and `begin_next_block` shrinks back at the next block.
            //
            // **The growth is bounded by [`ROLLING_BUFFER_CEILING_BYTES`], checked above.**
            // Without that ceiling a corrupt block grew this buffer to the block's whole
            // decompressed size, which is not bounded by its size on disk — measured, 4,132
            // bytes on disk against 67,125,248 held. Nothing is sized from a *declared* length
            // either way; the ceiling is what keeps a damaged file's cost bounded rather than
            // merely finite.
            self.rolling
                .reserve(self.rolling.capacity().max(ROLLING_BYTES));
        }

        if self.cursor.compressed_left == 0 {
            // The block's declared bytes are spent. If zstd has not finished the frame, the
            // length in front of the block is shorter than the frame behind it — and the answer
            // is *this block has no more to give*, which the caller turns into damage. **Checked
            // before the buffer, not inside it**: the next block's bytes are usually already
            // buffered, and the version that looked at the buffer first fed the decoder an empty
            // slice instead and left zstd to produce an error code nobody can act on.
            self.cursor.inflated = true;
            return Ok(false);
        }
        if self.compressed_at == self.compressed_filled && !self.fill_from_source()? {
            return Err(BlockReadError::FileEndsInsideABlock {
                compressed_bytes_left: self.cursor.compressed_left,
            });
        }

        let feeding =
            (self.compressed_filled - self.compressed_at).min(self.cursor.compressed_left);
        let mut input = zstd::zstd_safe::InBuffer {
            src: &self.compressed[self.compressed_at..self.compressed_at + feeding],
            pos: 0,
        };
        let at = self.rolling.len();
        let mut output = zstd::zstd_safe::OutBuffer::around_pos(&mut self.rolling, at);
        let hint = self
            .decoder
            .decompress_stream(&mut output, &mut input)
            .map_err(|code| BlockReadError::zstd("decompressing a block", code))?;
        let taken = input.pos;
        self.compressed_at += taken;
        self.cursor.compressed_left -= taken;
        if hint == 0 {
            self.cursor.inflated = true;
        }
        Ok(true)
    }

    /// Pull the next chunk of compressed bytes. `false` at the end of the source.
    fn fill_from_source(&mut self) -> Result<bool, BlockReadError> {
        self.compressed_at = 0;
        self.compressed_filled = 0;
        loop {
            match self.source.read(&mut self.compressed) {
                Ok(0) => return Ok(false),
                Ok(read) => {
                    self.compressed_filled = read;
                    return Ok(true);
                }
                Err(fault) if fault.kind() == std::io::ErrorKind::Interrupted => {}
                Err(source) => {
                    return Err(BlockReadError::Io {
                        while_doing: "reading a block's compressed bytes",
                        source,
                    });
                }
            }
        }
    }

    /// Read exactly `out.len()` bytes through the compressed buffer. `Ok(false)` when the source
    /// is finished and **nothing at all** was read, which is the ordinary end of a run of blocks;
    /// a source that ends part-way through is a truncated file.
    fn read_exactly(&mut self, out: &mut [u8]) -> Result<bool, BlockReadError> {
        let mut filled = 0usize;
        while filled < out.len() {
            if self.compressed_at == self.compressed_filled && !self.fill_from_source()? {
                if filled == 0 {
                    return Ok(false);
                }
                return Err(BlockReadError::FileEndsInsideABlockLength { bytes_read: filled });
            }
            let take = (self.compressed_filled - self.compressed_at).min(out.len() - filled);
            out[filled..filled + take]
                .copy_from_slice(&self.compressed[self.compressed_at..self.compressed_at + take]);
            self.compressed_at += take;
            filled += take;
        }
        Ok(true)
    }

    /// Mark the stream finished, so a caller that keeps asking after a refusal gets `None`
    /// rather than the same refusal for ever, or a record built from state a refusal left
    /// half-advanced.
    fn fail(&mut self, refused: BlockReadError) -> BlockReadError {
        self.refused = true;
        self.cursor = BlockCursor::between_blocks();
        // **The rolling buffer goes and the read chunk stays.** The rolling one may have grown
        // past what this reader budgeted for, and holding that until the caller drops the
        // reader is memory nobody asked for. The read chunk is the budget, so freeing it saves
        // nothing — and freeing it would make a refused reader stop for the wrong reason:
        // `read` into an empty buffer returns zero, which looks exactly like the end of the
        // source. What stops a refused reader is `refused`, and it has to be the only thing, or
        // the next coder to change either of these lines resurrects a reader that carries on
        // past damage.
        self.rolling = Vec::new();
        self.compressed_at = 0;
        self.compressed_filled = 0;
        refused
    }
}

/// Why a run of blocks could not be read back.
///
/// **Every variant is an input problem, not a bug**, and each carries a different instruction —
/// which is why [`RecordDecodeError`] and [`BlockHeadDecodeError`] are not folded into one class
/// here. A record that stopped early means *decompress more* to this reader and *the file is
/// damaged* to its caller, because by the time it reaches a caller the bytes have run out for
/// good.
///
/// [`RecordDecodeError`]: crate::ng::psp::RecordDecodeError
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum BlockReadError {
    /// The file's declared look-back window is outside what zstd takes, so no decoder can be
    /// configured for it.
    #[error(
        "the file declares a look-back window of 2^{look_back_window_log} bytes; zstd takes \
         between 2^{MIN_LOOK_BACK_WINDOW_LOG} and 2^{MAX_LOOK_BACK_WINDOW_LOG}"
    )]
    WindowLogOutOfRange { look_back_window_log: u8 },

    /// The file declares a record layout this reader does not read. **Upgrade the reader**; the
    /// file is not damaged.
    #[error("the file declares a record layout this reader cannot read: {source}")]
    UnsupportedRecordLayout {
        #[source]
        source: RecordLayoutError,
    },

    /// A record's bytes cannot mean what they say. The file is damaged.
    #[error("a record in this block is unreadable: {source}")]
    DamagedRecord {
        #[source]
        source: RecordDecodeError,
    },

    /// A record names something a later writer added. **Upgrade the reader**; the file is fine.
    #[error("a record in this block names something this reader does not know: {source}")]
    UnsupportedRecord {
        #[source]
        source: RecordDecodeError,
    },

    /// A record's bytes ran out at the end of its block, so no further bytes can arrive. **This
    /// is where the retry class becomes damage** — the same conversion `record.rs` makes for a
    /// field inside a bounded body, one level up.
    #[error(
        "a record runs past the end of its block, with {records_left} records still declared \
         and {bytes_left} bytes left"
    )]
    RecordRunsPastItsBlock {
        records_left: u64,
        bytes_left: usize,
    },

    /// The block's opening fields ran out before the block did.
    #[error("a block's {field} runs past the end of the block, {bytes_in} bytes in")]
    BlockHeadRunsPastItsBlock {
        field: &'static str,
        bytes_in: usize,
    },

    /// A block's opening fields cannot mean what they say.
    #[error("a block's head is unreadable: {source}")]
    DamagedBlockHead {
        #[source]
        source: BlockHeadDecodeError,
    },

    /// One record needs more of the rolling buffer than this reader allows.
    ///
    /// **A knob, not a rebuild**, which is why it is its own class: a genuine record this large
    /// is read by raising [`ROLLING_BUFFER_CEILING_BYTES`], the same instruction spec §4.2
    /// attaches to a look-back window wider than a reader budgeted for. On a damaged block —
    /// which is where it is actually expected — it is what stops the buffer growing to the
    /// block's whole decompressed size.
    #[error(
        "a record in the block starting at {block} needs more than the {allowed_bytes} bytes \
         this reader allows one record to hold; raise the ceiling to read it"
    )]
    RecordLargerThanTheReaderAllows {
        /// Which block it was in — the only thing that locates the record, since its own head
        /// never finished parsing.
        block: GenomeRegion,
        /// How many records of that block were handed over before this one.
        records_read: u64,
        allowed_bytes: usize,
    },

    /// A ceiling at or under the rolling buffer, which would refuse records the buffer holds
    /// without growing at all.
    #[error(
        "a buffer ceiling of {ceiling} bytes, which is not above the {buffer_bytes}-byte buffer \
         it is a ceiling on"
    )]
    BufferCeilingUnderTheBuffer { ceiling: usize, buffer_bytes: usize },

    /// The block held bytes past the last record its head declared.
    #[error("a block holds {bytes_left} bytes past the last record its head declared")]
    BlockHoldsMoreThanItDeclared { bytes_left: usize },

    /// The length in front of the block claims more compressed bytes than its frame used. On a
    /// run of blocks that length reaches into the next block, and a reader that believed it
    /// would hand back records from a block nobody asked for.
    #[error(
        "a block's length claims {compressed_bytes_left} compressed bytes more than its frame \
         used"
    )]
    BlockDeclaresMoreBytesThanItsFrame { compressed_bytes_left: usize },

    /// The block's records ran out while its frame had not finished: the two disagree about
    /// where the block ends.
    #[error("a block's records end before its compressed frame does")]
    BlockFrameDidNotEndWithItsRecords,

    /// The file ended inside a block's compressed bytes. **Refuse rather than read short**: a
    /// run that was killed part-way must not look like a sample covering less of the genome.
    #[error("the file ends with {compressed_bytes_left} compressed bytes of a block still due")]
    FileEndsInsideABlock { compressed_bytes_left: usize },

    /// The file ended inside the four bytes that say how long a block is.
    #[error("the file ends {bytes_read} bytes into a block's length")]
    FileEndsInsideABlockLength { bytes_read: usize },

    /// zstd refused: the frame is damaged, or it needs a larger window than this reader was
    /// configured for.
    ///
    /// **It names what was being done and what zstd said**, because zstd's own answer is a
    /// numeric code — and the two cases a caller most needs to tell apart, a corrupt frame and a
    /// window too small, differ *only* in that code. Resolved, they read "Restored data doesn't
    /// match checksum" and "Frame requires too much memory for decoding".
    #[error("zstd failed while {while_doing}: {said}")]
    Zstd {
        while_doing: &'static str,
        /// zstd's own account, resolved from its code.
        said: &'static str,
        code: usize,
    },

    /// Reading the file's bytes failed.
    #[error("the file could not be read while {while_doing}")]
    Io {
        while_doing: &'static str,
        #[source]
        source: std::io::Error,
    },
}

impl BlockReadError {
    fn zstd(while_doing: &'static str, code: usize) -> Self {
        Self::Zstd {
            while_doing,
            said: zstd::zstd_safe::get_error_name(code),
            code,
        }
    }

    /// Put a record's own fault in the class its instruction belongs to.
    ///
    /// **`Truncated` never arrives here**, and that is the point: the reader answers it by
    /// decompressing more, and turns it into
    /// [`RecordRunsPastItsBlock`](Self::RecordRunsPastItsBlock) only once the block has no more
    /// to give. If one ever did reach here it would be damage — no further bytes are coming —
    /// which is what the arm below says.
    fn from_record(fault: RecordDecodeError) -> Self {
        match fault {
            RecordDecodeError::Unsupported { .. } => Self::UnsupportedRecord { source: fault },
            RecordDecodeError::Malformed { .. } | RecordDecodeError::Truncated { .. } => {
                Self::DamagedRecord { source: fault }
            }
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
    use crate::pileup_record::ChainId;
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
        // **One reader across the block, because the chain-id changes carry state.** A fresh one
        // per record would meet the second record's departures against an empty set.
        let mut live_reads = LiveSetReader::new();
        live_reads.start_block();
        let mut at = 0usize;
        let mut records = Vec::new();
        while at < found.records.len() {
            let decoded = decode_record(
                &found.records[at..],
                found.head.contig,
                measured_from,
                &mut live_reads,
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

    /// **A grid cell with no records produces no block.** The cut is a rule for deciding where a
    /// block *ends*, not an instruction to emit one per 100 kb of reference: a block exists
    /// because a record went into it, so a sample covering two cells of a chromosome writes two
    /// blocks and not one per cell between them.
    ///
    /// The type carries the same statement — a block head's record count is a `NonZeroU64`, so a
    /// block holding none has no representation — and this is the behaviour that matters at the
    /// other end: a thin sample does not pay an index entry and a compressed frame for every
    /// stretch of reference it did not cover.
    ///
    /// **And its hardest-working line is `blocks.len() == 3`, which pins the owner's ruling
    /// against merging** (2026-08-27, recorded on [`OpenBlock`]): a writer that accumulated
    /// across the ninety empty cells here would emit **two** blocks, not three. That is the newer
    /// and more surprising of the two decisions, so it is worth saying which assertion carries
    /// it. ⚠ What the fixture pins is a merge rule with a *byte* threshold, as spec §4.1
    /// describes; one keyed on a cell count below 90 would survive this file.
    #[test]
    fn a_grid_cell_with_no_records_produces_no_block() {
        // Two records ninety cells apart on one contig, and one on the next contig.
        let records = vec![
            a_record(0, 1, 1),
            a_record(0, 90 * A_GRID.get() + 7, 1),
            a_record(1, 40 * A_GRID.get() + 3, 1),
        ];
        let blocks =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");

        assert_eq!(
            blocks.len(),
            3,
            "one block for each record's cell, and none for the ninety cells between them"
        );
        for block in &blocks {
            assert_eq!(walk(block).records.len(), 1);
        }
        assert_eq!(walk_all(&blocks), records);
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
        // The reader's ceiling on one record, which is a budget rather than a format rule — so
        // moving it is a decision taken here rather than a memory test that stops binding.
        assert_eq!(ROLLING_BUFFER_CEILING_BYTES, 512 * 1024);
        let manifest = a_manifest();
        assert_eq!(
            manifest.genomic_block_size_bp,
            DEFAULT_GENOMIC_BLOCK_SIZE_BP
        );
        assert_eq!(manifest.block_byte_ceiling, DEFAULT_BLOCK_BYTE_CEILING);
        assert_eq!(manifest.look_back_window_log, DEFAULT_LOOK_BACK_WINDOW_LOG);
        assert_eq!(manifest.fields, record_fields());
    }

    /// **The argument for the ceiling's value, as something that can fail.**
    ///
    /// The test above pins the number; this pins the reason, which is what the number is *for*.
    /// [`ROLLING_BUFFER_CEILING_BYTES`] is a per-open-sample figure, so what it costs is set by
    /// the cohort, and the cohort this caller is committed to runs to several thousand samples
    /// (`design_principles.md` §0). Spec §1.1 budgets an open sample at 500 kB and says that is
    /// what makes a run of that size fit.
    ///
    /// So the worst case — every sample meeting a damaged block at the same moment — has to come
    /// out at the budget rather than at a multiple of it. At the 1 MiB first proposed it came out
    /// at 3.07 GB, twice the whole budget, which is a weak bound for something whose only job is
    /// bounding hostile input.
    #[test]
    fn the_buffer_ceiling_is_priced_at_the_top_of_the_committed_cohort_range() {
        const SAMPLES_AT_THE_TOP_OF_THE_RANGE: usize = 3_000;
        const WHAT_SPEC_1_1_BUDGETS_AN_OPEN_SAMPLE: usize = 500_000;

        let worst_case = ROLLING_BUFFER_CEILING_BYTES * SAMPLES_AT_THE_TOP_OF_THE_RANGE;
        assert_eq!(worst_case, 1_572_864_000, "3,000 readers at the ceiling");
        assert!(
            worst_case
                <= WHAT_SPEC_1_1_BUDGETS_AN_OPEN_SAMPLE * SAMPLES_AT_THE_TOP_OF_THE_RANGE * 21 / 20,
            "the ceiling has to stay within a twentieth of spec §1.1's whole-cohort budget of \
             {} bytes, and {worst_case} does not",
            WHAT_SPEC_1_1_BUDGETS_AN_OPEN_SAMPLE * SAMPLES_AT_THE_TOP_OF_THE_RANGE
        );
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
        // **Varied records, not a repeated one.** A payload of identical records compresses to
        // the same size at every level, because there is nothing for a longer search to find:
        // measured, the window-spanning fixture used here before gave 91 bytes at both level 1
        // and level 9, so the middle of the ordering below could not be seen.
        let records: Vec<_> = (1..4_000u64)
            .map(|at| a_record(0, at, 1 + at % 23))
            .collect();
        let payload = one_payload(&records);
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
    // Reading a block back, a record at a time
    // -----------------------------------------------------------------

    /// A whole file of compressed blocks, and where each one begins — which is what Milestone
    /// F's index will hold and what lets a test start a reader at an arbitrary block.
    struct BlocksOnDisk {
        bytes: Vec<u8>,
        block_offsets: Vec<usize>,
        payloads: Vec<Vec<u8>>,
    }

    /// Cut `records` into blocks and compress each one, as a writer will.
    fn blocks_on_disk(
        records: &[SampleLocusObservations],
        grid: Bp,
        ceiling: Option<u32>,
    ) -> BlocksOnDisk {
        let manifest = a_manifest();
        let payloads =
            cut(BlockBuilder::new(grid, ceiling).expect("a grid"), records).expect("in order");
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let mut bytes = Vec::new();
        let mut block_offsets = Vec::new();
        for payload in &payloads {
            block_offsets.push(bytes.len());
            bytes.extend_from_slice(compressor.compress(payload).expect("it compresses"));
        }
        BlocksOnDisk {
            bytes,
            block_offsets,
            payloads,
        }
    }

    /// Read every record a reader over `bytes` yields, and refuse to loop for ever.
    fn stream_every_record(bytes: &[u8]) -> Result<Vec<StreamedRecord>, BlockReadError> {
        stream_records_where(bytes, |_| true)
    }

    fn stream_records_where<F>(
        bytes: &[u8],
        mut want: F,
    ) -> Result<Vec<StreamedRecord>, BlockReadError>
    where
        F: FnMut(&RecordHead) -> bool,
    {
        let manifest = a_manifest();
        let mut stream = BlockStream::new(bytes, &manifest).expect("a valid manifest");
        let mut met = Vec::new();
        while let Some(next) = stream.next_record_where(&mut want) {
            met.push(next?);
            assert!(
                met.len() < 1_000_000,
                "a reader that never ends is a reader that is looping"
            );
        }
        // **And it stays finished.** A caller that keeps asking after the end gets nothing more,
        // rather than the last record again or a refusal that repeats.
        assert!(stream.next_record().is_none());
        Ok(met)
    }

    /// The records a walk built, in order — for comparing against what was written.
    fn built(met: &[StreamedRecord]) -> Vec<SampleLocusObservations> {
        met.iter().filter_map(|one| one.record.clone()).collect()
    }

    #[test]
    fn a_stream_reads_back_every_record_that_was_written() {
        let records: Vec<_> = (0..200).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(300));
        assert!(
            on_disk.payloads.len() >= 4,
            "several blocks, or the reader never crosses one"
        );

        let met = stream_every_record(&on_disk.bytes).expect("it reads");
        assert_eq!(built(&met), records);
        assert_eq!(met.len(), records.len());
    }

    /// **The blocks a reader crosses are the blocks the writer cut**, and every record says
    /// which one it came from.
    #[test]
    fn every_record_says_which_block_it_came_from() {
        let records: Vec<_> = (0..200).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(300));
        let met = stream_every_record(&on_disk.bytes).expect("it reads");

        let mut blocks_met: Vec<BlockHead> = Vec::new();
        for one in &met {
            if blocks_met.last() != Some(&one.block) {
                blocks_met.push(one.block);
            }
            assert_eq!(one.head.region.contig, one.block.contig);
        }
        assert_eq!(blocks_met.len(), on_disk.payloads.len());
        for (index, block) in blocks_met.iter().enumerate() {
            let walked = walk(&on_disk.payloads[index]);
            assert_eq!(*block, walked.head, "block {index}");
            assert_eq!(
                met.iter().filter(|one| one.block == *block).count() as u64,
                block.record_count.get()
            );
        }
    }

    /// A stream over records on several contigs and several grid cells reads them all back, in
    /// the order they were written.
    #[test]
    fn a_stream_crosses_contigs_and_grid_cells() {
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
            }
        }
        let on_disk = blocks_on_disk(&records, A_GRID, None);
        assert_eq!(on_disk.payloads.len(), 9);
        assert_eq!(
            built(&stream_every_record(&on_disk.bytes).expect("it reads")),
            records
        );
    }

    /// **A record the caller declines costs its head and a pointer advance.** The heads a
    /// skipping walk meets are the heads a full walk meets, and the records it builds are the
    /// ones it asked for — byte for byte the same as a full decode of those.
    #[test]
    fn a_skipping_walk_meets_every_head_and_builds_only_what_it_wants() {
        let records: Vec<_> = (0..200).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(300));
        let whole = stream_every_record(&on_disk.bytes).expect("it reads");

        for (name, keep) in [
            ("every fourth", 4usize),
            ("every twentieth", 20),
            ("none at all", usize::MAX),
        ] {
            let mut index = 0usize;
            let met = stream_records_where(&on_disk.bytes, |_| {
                let wanted = keep != usize::MAX && index.is_multiple_of(keep);
                index += 1;
                wanted
            })
            .unwrap_or_else(|refused| panic!("{name}: {refused}"));

            assert_eq!(
                met.iter().map(|one| one.head).collect::<Vec<_>>(),
                whole.iter().map(|one| one.head).collect::<Vec<_>>(),
                "{name}: every head is met whatever is built"
            );
            let wanted: Vec<_> = records
                .iter()
                .enumerate()
                .filter(|(at, _)| keep != usize::MAX && at.is_multiple_of(keep))
                .map(|(_, record)| record.clone())
                .collect();
            assert_eq!(built(&met), wanted, "{name}: exactly what was asked for");
        }
    }

    /// **A reader can start at any block, and gets what a full read gives from there.** The
    /// offsets are what Milestone F's index will hold; nothing else about the file is needed.
    #[test]
    fn a_reader_starting_at_a_block_gets_the_tail_of_a_full_read() {
        let records: Vec<_> = (0..200).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(300));
        let whole = stream_every_record(&on_disk.bytes).expect("it reads");
        assert!(on_disk.block_offsets.len() >= 4);

        let mut skipped = 0usize;
        for (index, offset) in on_disk.block_offsets.iter().enumerate() {
            let from_here = stream_every_record(&on_disk.bytes[*offset..])
                .unwrap_or_else(|refused| panic!("from block {index}: {refused}"));
            assert_eq!(
                from_here,
                whole[skipped..],
                "block {index} at byte {offset}"
            );
            skipped += usize::try_from(walk(&on_disk.payloads[index]).head.record_count.get())
                .expect("a small count");
        }
        assert_eq!(skipped, records.len());
    }

    /// A record larger than the rolling buffer is read, not refused — spec §8 says a fixed
    /// maximum record size is not a safe assumption — and the buffer goes back to what a reader
    /// budgets for at the next block.
    #[test]
    fn a_record_larger_than_the_rolling_buffer_is_read_and_the_buffer_shrinks_back() {
        let mut enormous = a_record(0, 500, 1);
        enormous.observations = (0..4_000u32)
            .map(|read| SequenceObservation {
                bases: vec![b"ACGT"[(read % 4) as usize]; 24].into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read % 7),
                num_obs: 1,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(i64::from(read) + 1)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 1,
                chain_ids: Vec::new(),
            })
            .collect();
        let records = vec![enormous, a_record(0, 200_000, 1)];

        let payloads =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");
        assert!(
            payloads[0].len() > ROLLING_BYTES,
            "the record must exceed the rolling buffer, or this proves nothing: {} bytes \
             against {ROLLING_BYTES}",
            payloads[0].len()
        );

        let on_disk = blocks_on_disk(&records, A_GRID, None);
        let manifest = a_manifest();
        let mut stream = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let first = stream.next_record().expect("a record").expect("it reads");
        assert_eq!(first.record.as_ref(), Some(&records[0]));

        let second = stream.next_record().expect("a record").expect("it reads");
        assert_eq!(second.record.as_ref(), Some(&records[1]));
        assert_eq!(
            stream.buffered_bytes(),
            READ_CHUNK_BYTES + ROLLING_BYTES,
            "the buffer must shrink back at the next block"
        );
        assert!(stream.next_record().is_none());
    }

    /// **The rolling buffer rolls, and that is half of what the format is for.** Spec §5.1 says
    /// two conditions make a reader cheap and only the first is the compressor's doing: do not
    /// inflate the whole block, *and do not accumulate what you inflated*. Satisfying only the
    /// first moves the memory rather than removing it.
    ///
    /// So this reads a file whose blocks are many times the rolling buffer and asserts the
    /// reader never holds more than it budgeted for. **Nothing else in the module sees it**:
    /// removing the line that drops consumed bytes before pumping leaves every other test green,
    /// while the buffer grows to the whole block.
    #[test]
    fn a_reader_holds_its_two_buffers_and_never_the_block() {
        // Blocks far larger than the rolling buffer, over several grid cells.
        let records: Vec<_> = (1..12_000u64).map(|at| a_record(0, at, 1)).collect();
        let on_disk = blocks_on_disk(&records, Bp(4_000), None);
        let biggest = on_disk
            .payloads
            .iter()
            .map(|payload| payload.len())
            .max()
            .expect("some blocks");
        assert!(
            biggest > ROLLING_BYTES * 4,
            "a block must be several times the rolling buffer, or this proves nothing: \
             {biggest} bytes against {ROLLING_BYTES}"
        );

        let manifest = a_manifest();
        let mut stream = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let budget = READ_CHUNK_BYTES + ROLLING_BYTES;
        let mut most_held = stream.buffered_bytes();
        let mut back = Vec::new();
        while let Some(next) = stream.next_record() {
            // **Compared, not counted.** This walk crosses the retry arm thousands of times —
            // it is the module's densest exercise of it — and a retry that advanced the
            // coordinate before failing would give the right *number* of records at the wrong
            // positions, which a count cannot see.
            back.push(
                next.expect("it reads")
                    .record
                    .expect("every record is built"),
            );
            most_held = most_held.max(stream.buffered_bytes());
        }
        assert_eq!(back, records);
        assert!(
            most_held <= budget,
            "a reader over {biggest}-byte blocks held {most_held} bytes against a budget of \
             {budget}"
        );
    }

    /// **A reader is configured from the file's manifest and refuses what it cannot honour.**
    /// A window zstd does not take, and a record layout this build does not write, are both
    /// refused at the reader rather than block by block.
    #[test]
    fn a_reader_refuses_a_manifest_it_cannot_honour() {
        let mut manifest = a_manifest();
        manifest.look_back_window_log = MAX_LOOK_BACK_WINDOW_LOG + 1;
        assert!(matches!(
            BlockStream::new([].as_slice(), &manifest),
            Err(BlockReadError::WindowLogOutOfRange { .. })
        ));

        let mut manifest = a_manifest();
        manifest.fields.truncate(manifest.fields.len() - 1);
        assert!(matches!(
            BlockStream::new([].as_slice(), &manifest),
            Err(BlockReadError::UnsupportedRecordLayout { .. })
        ));

        let good = a_manifest();
        let stream = BlockStream::new([].as_slice(), &good).expect("a manifest this build wrote");
        assert_eq!(stream.look_back_window_log(), good.look_back_window_log);
    }

    /// **A window narrower than the file's is refused, and the message says both numbers.**
    /// A reader whose budget is smaller than what the file needs is the case spec §4.2 exists
    /// for: zstd's own answer is a code, and this one names the exponent.
    #[test]
    fn a_reader_whose_window_is_narrower_than_the_file_refuses_the_block() {
        let records = records_past_a_window(DEFAULT_LOOK_BACK_WINDOW_LOG);
        let payload = one_payload(&records);
        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let on_disk = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let mut narrower = a_manifest();
        narrower.look_back_window_log = DEFAULT_LOOK_BACK_WINDOW_LOG - 1;
        let mut stream = BlockStream::new(on_disk.as_slice(), &narrower).expect("a valid window");
        match stream.next_record() {
            Some(Err(BlockReadError::Zstd { while_doing, .. })) => {
                assert_eq!(while_doing, "decompressing a block");
            }
            other => panic!("expected a zstd refusal, got {other:?}"),
        }
        assert!(
            stream.next_record().is_none(),
            "and a refused stream stays refused rather than repeating itself"
        );
    }

    /// **An empty source is an empty read, not an error.** A run of blocks that holds none is
    /// what a writer that was handed no records produces.
    #[test]
    fn a_source_with_no_blocks_reads_nothing() {
        assert!(
            stream_every_record(&[])
                .expect("no blocks is not a fault")
                .is_empty()
        );
    }

    /// **A file cut short anywhere is refused, at every cut.** Spec §3.3's rule is that a run
    /// killed part-way must not look like a sample covering less of the genome, and the reader
    /// is where that starts: it refuses rather than stopping early and reporting success.
    #[test]
    fn a_file_cut_short_anywhere_inside_a_block_is_refused() {
        let records: Vec<_> = (0..40).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(120));
        assert!(on_disk.block_offsets.len() >= 3);
        let whole = stream_every_record(&on_disk.bytes).expect("it reads");

        for cut_at in 1..on_disk.bytes.len() {
            let short = &on_disk.bytes[..cut_at];
            let met = stream_every_record(short);
            match met {
                Err(_) => {}
                Ok(records_read) => {
                    // A cut exactly at a block boundary is a shorter *complete* file, which is
                    // the one case that is not damage — and it must be a prefix of the whole.
                    assert!(
                        on_disk.block_offsets.contains(&cut_at),
                        "a file cut at {cut_at} bytes, which is not a block boundary, must be \
                         refused"
                    );
                    assert_eq!(records_read, whole[..records_read.len()]);
                }
            }
        }
    }

    /// A block whose head declares more records than it holds is refused, and one that holds
    /// more than it declares is refused too. **Both are the count doing its job**: it is what
    /// lets a reader say *the block ended where it should have* rather than *the bytes ran out*.
    #[test]
    fn a_block_whose_record_count_disagrees_with_its_records_is_refused() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();
        let payload = one_payload(&records);
        let found = BlockRecords::split(&payload).expect("it splits");
        assert_eq!(found.head.record_count.get(), 6);

        for miscount in [5u64, 7] {
            let mut damaged = Vec::new();
            BlockHead {
                record_count: NonZeroU64::new(miscount).expect("not zero"),
                ..found.head
            }
            .encode(&mut damaged);
            damaged.extend_from_slice(found.records);

            let mut compressor =
                BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
            let on_disk = compressor
                .compress(&damaged)
                .expect("it compresses")
                .to_vec();
            let refused = stream_every_record(&on_disk)
                .expect_err("a block that miscounts its records is damaged");
            match (miscount, &refused) {
                (5, BlockReadError::BlockHoldsMoreThanItDeclared { .. }) => {}
                (7, BlockReadError::RecordRunsPastItsBlock { .. }) => {}
                _ => panic!("a count of {miscount} gave {refused}"),
            }
        }
    }

    /// A block whose head cannot mean what it says is refused as damage, not as a short read.
    #[test]
    fn a_block_whose_head_is_damaged_is_refused() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();
        let payload = one_payload(&records);
        let found = BlockRecords::split(&payload).expect("it splits");

        // A contig id no contig id could be, in front of the records that were there.
        let mut damaged = Vec::new();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut damaged);
        encode_u64_leb128(found.head.first_position.get(), &mut damaged);
        encode_u64_leb128(found.head.record_count.get(), &mut damaged);
        damaged.extend_from_slice(found.records);

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let on_disk = compressor
            .compress(&damaged)
            .expect("it compresses")
            .to_vec();
        assert!(matches!(
            stream_every_record(&on_disk),
            Err(BlockReadError::DamagedBlockHead { .. })
        ));
    }

    /// **A record's bytes reaching the end of its block is damage, not a retry.** Inside the
    /// reader a record that stops early means *decompress more*; once the block has no more to
    /// give, no quantity of further bytes changes the answer, and a reader that kept asking
    /// would never finish.
    #[test]
    fn a_record_that_runs_past_its_block_is_damage_and_not_a_retry() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();
        let payload = one_payload(&records);
        let found = BlockRecords::split(&payload).expect("it splits");

        // The same records with the last one's tail cut off, and the head still declaring six.
        let mut damaged = Vec::new();
        found.head.encode(&mut damaged);
        damaged.extend_from_slice(&found.records[..found.records.len() - 3]);

        let mut compressor =
            BlockCompressor::new(DEFAULT_LOOK_BACK_WINDOW_LOG).expect("a valid window");
        let on_disk = compressor
            .compress(&damaged)
            .expect("it compresses")
            .to_vec();

        let manifest = a_manifest();
        let mut stream = BlockStream::new(on_disk.as_slice(), &manifest).expect("a manifest");
        let mut met = 0usize;
        let refused = loop {
            match stream.next_record() {
                Some(Ok(_)) => {
                    met += 1;
                    assert!(met <= 6, "a reader that keeps going is one that is looping");
                }
                Some(Err(refused)) => break refused,
                None => panic!("a truncated block must refuse, not end"),
            }
        };
        assert!(
            matches!(refused, BlockReadError::RecordRunsPastItsBlock { .. }),
            "got {refused}"
        );
        assert_eq!(met, 5, "the five whole records were handed over first");
    }

    /// **A refused stream stays refused.** A caller that keeps asking after a refusal gets
    /// nothing more — not the same refusal for ever, and not a record built from state the
    /// refusal left half-advanced. Every way a stream can refuse is checked, because the state a
    /// refusal leaves behind differs by which one it was.
    #[test]
    fn a_stream_that_refuses_yields_nothing_afterwards() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();
        let payload = one_payload(&records);
        let found = BlockRecords::split(&payload).expect("it splits");
        let manifest = a_manifest();
        let compress = |bytes: &[u8]| {
            BlockCompressor::from_manifest(&manifest)
                .expect("a window")
                .compress(bytes)
                .expect("it compresses")
                .to_vec()
        };

        // A record cut short at the end of its block; a block declaring one record too few; a
        // block head no reader can believe; and a file that stops inside a block's bytes.
        let mut short_record = Vec::new();
        found.head.encode(&mut short_record);
        short_record.extend_from_slice(&found.records[..found.records.len() - 3]);

        let mut too_few = Vec::new();
        BlockHead {
            record_count: NonZeroU64::new(5).expect("not zero"),
            ..found.head
        }
        .encode(&mut too_few);
        too_few.extend_from_slice(found.records);

        let mut bad_head = Vec::new();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut bad_head);
        encode_u64_leb128(found.head.first_position.get(), &mut bad_head);
        encode_u64_leb128(found.head.record_count.get(), &mut bad_head);
        bad_head.extend_from_slice(found.records);

        let whole = compress(&payload);
        let files: [(&str, Vec<u8>); 4] = [
            ("a record cut short", compress(&short_record)),
            ("a block declaring too few records", compress(&too_few)),
            ("a block head that cannot be believed", compress(&bad_head)),
            (
                "a file that stops inside a block",
                whole[..whole.len() - 5].to_vec(),
            ),
        ];

        for (what, bytes) in files {
            let mut stream = BlockStream::new(bytes.as_slice(), &manifest).expect("a manifest");
            let mut refused = None;
            for _ in 0..records.len() + 4 {
                match stream.next_record() {
                    Some(Ok(_)) => {}
                    Some(Err(fault)) => {
                        refused = Some(fault);
                        break;
                    }
                    None => break,
                }
            }
            let refused = refused.unwrap_or_else(|| panic!("{what} must be refused"));
            for again in 0..3 {
                assert!(
                    stream.next_record().is_none(),
                    "{what}: after {refused}, ask {again} gave something back"
                );
            }
        }
    }

    /// **A block whose declared length is shorter than its own frame is refused, not fed past
    /// its end.** The decoder is handed at most the bytes the block declared: without that cap
    /// zstd consumes into whatever follows, and the count of bytes still due goes below zero.
    #[test]
    fn a_block_declaring_fewer_bytes_than_its_frame_is_refused() {
        let records: Vec<_> = (0..40).map(|index| a_record(0, 100 + index, 1)).collect();
        let payload = one_payload(&records);
        let manifest = a_manifest();
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let honest = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let frame_bytes = honest.len() - COMPRESSED_BLOCK_LENGTH_BYTES;
        let mut lying = honest.clone();
        lying[..COMPRESSED_BLOCK_LENGTH_BYTES].copy_from_slice(
            &u32::try_from(frame_bytes - 8)
                .expect("a small frame")
                .to_le_bytes(),
        );
        // Something after it, so a decoder fed past the block's end has bytes to run into.
        lying.extend_from_slice(&honest);

        let mut stream = BlockStream::new(lying.as_slice(), &manifest).expect("a manifest");
        let mut refused = None;
        for _ in 0..records.len() + 4 {
            match stream.next_record() {
                Some(Ok(_)) => {}
                Some(Err(fault)) => {
                    refused = Some(fault);
                    break;
                }
                None => break,
            }
        }
        let refused = refused.expect("a block that lies about its length must be refused");
        // **The class is this module's, not zstd's.** A truncated frame may decompress to
        // nothing at all, so which of the "this block has no more to give" refusals arrives
        // depends on how far the frame got — but a zstd error code is the one answer nobody can
        // act on, and it is what the version that fed the decoder an empty slice produced.
        assert!(
            matches!(
                refused,
                BlockReadError::BlockHeadRunsPastItsBlock { .. }
                    | BlockReadError::RecordRunsPastItsBlock { .. }
                    | BlockReadError::BlockHoldsMoreThanItDeclared { .. }
                    | BlockReadError::DamagedRecord { .. }
                    | BlockReadError::DamagedBlockHead { .. }
            ),
            "got {refused}"
        );
    }

    /// A source that hands back at most `most_bytes_a_read` bytes a read, and optionally fails with
    /// `Interrupted` every so often — which is what a pipe, a socket or a signalled file does,
    /// and what no fixture built from a `&[u8]` ever does.
    struct DribblingSource {
        bytes: Vec<u8>,
        at: usize,
        most_bytes_a_read: usize,
        interrupt_every: usize,
        reads: usize,
    }

    impl DribblingSource {
        fn new(bytes: Vec<u8>, most_bytes_a_read: usize) -> Self {
            Self::interrupted(bytes, most_bytes_a_read, 0)
        }

        fn interrupted(bytes: Vec<u8>, most_bytes_a_read: usize, interrupt_every: usize) -> Self {
            Self {
                bytes,
                at: 0,
                most_bytes_a_read,
                interrupt_every,
                reads: 0,
            }
        }
    }

    impl std::io::Read for DribblingSource {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            self.reads += 1;
            if self.interrupt_every > 0 && self.reads.is_multiple_of(self.interrupt_every) {
                return Err(std::io::Error::from(std::io::ErrorKind::Interrupted));
            }
            let take = (self.bytes.len() - self.at)
                .min(self.most_bytes_a_read)
                .min(out.len());
            out[..take].copy_from_slice(&self.bytes[self.at..self.at + take]);
            self.at += take;
            Ok(take)
        }
    }

    /// One record whose bases zstd cannot shrink, so a file of them is large on disk.
    ///
    /// **`a_record`'s bodies are a four-letter cycle and compress to almost nothing** — twelve
    /// thousand of them are 4,543 bytes on disk at the `Bp(200)` grid `a_file_of_several_read_chunks`
    /// uses, and 165 bytes at the 100 kb default: a quarter of one read chunk at best. A
    /// reader's behaviour past its first read cannot be tested with a fixture that never fills
    /// one. (An earlier version of this comment said 1,823 bytes, which is neither.)
    fn an_incompressible_record(start: u64) -> SampleLocusObservations {
        // A multiplicative hash over the position, taken a byte at a time: no run of it repeats
        // within a file, so zstd's match finder has nothing to work with.
        let mut noise = start.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let bases: Vec<u8> = (0..48)
            .map(|_| {
                noise = noise
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[((noise >> 33) & 3) as usize]
            })
            .collect();
        let mut bulky = a_record(0, start, 1);
        bulky.observations = (0..4u32)
            .map(|read| SequenceObservation {
                bases: bases[read as usize * 12..(read as usize + 1) * 12]
                    .to_vec()
                    .into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read),
                num_obs: 1 + read,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(noise as i64 & 0xffff)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: read % 2,
                chain_ids: Vec::new(),
            })
            .collect();
        bulky
    }

    /// A file of many blocks, larger than four read chunks — which every earlier fixture was
    /// not, and which is why several properties below could not fail before.
    fn a_file_of_several_read_chunks() -> (Vec<SampleLocusObservations>, BlocksOnDisk) {
        let records: Vec<_> = (1..3_000u64).map(an_incompressible_record).collect();
        let on_disk = blocks_on_disk(&records, Bp(200), None);
        assert!(
            on_disk.bytes.len() > READ_CHUNK_BYTES * 4,
            "the file must be several read chunks, or a refused reader stops for the wrong \
             reason; it is {} bytes",
            on_disk.bytes.len()
        );
        assert!(on_disk.block_offsets.len() > 10);
        (records, on_disk)
    }

    /// **A refused stream is finished, and on a file bigger than one read chunk that is a
    /// property rather than an accident.**
    ///
    /// Before this, nothing marked a reader as having refused: what stopped one was that its
    /// compressed buffer had been emptied, which lasts exactly as long as the whole file fits in
    /// a single 16 kB read. Past that the reader read on from wherever the source stood, took
    /// four arbitrary bytes for a block length, and carried on — measured by three review
    /// agents, one of which saw a file refuse after 5,980 records and then hand back **5,681
    /// more** before ending cleanly. Every earlier fixture here was a few hundred bytes, so the
    /// test named for this property could not fail.
    #[test]
    fn a_stream_that_refuses_stays_refused_on_a_file_larger_than_one_read() {
        let (_, on_disk) = a_file_of_several_read_chunks();

        // Damage inside the first block's frame, so the refusal comes early and most of the
        // file is still ahead of the reader.
        let first_frame_at = COMPRESSED_BLOCK_LENGTH_BYTES + 8;
        let mut damaged = on_disk.bytes.clone();
        damaged[first_frame_at] ^= 0xff;

        let manifest = a_manifest();
        // Both a source that hands over everything at once and one that dribbles: the second is
        // what a pipe or a socket does, and it is where a resynchronising reader lands exactly
        // on a block boundary and carries on as though nothing had happened.
        /// A source to try, made fresh so each arm starts from the file's first byte.
        type ASourceToTry = (&'static str, Box<dyn FnOnce() -> Box<dyn std::io::Read>>);
        let sources: [ASourceToTry; 2] = [
            ("a whole-file source", {
                let bytes = damaged.clone();
                Box::new(move || Box::new(std::io::Cursor::new(bytes)))
            }),
            ("a source that dribbles", {
                let bytes = damaged.clone();
                Box::new(move || Box::new(DribblingSource::new(bytes, 64)))
            }),
        ];
        for (what, source) in sources {
            let mut stream = BlockStream::new(source(), &manifest).expect("a manifest");
            let mut refused = None;
            for _ in 0..20_000 {
                match stream.next_record() {
                    Some(Ok(_)) => {}
                    Some(Err(fault)) => {
                        refused = Some(fault);
                        break;
                    }
                    None => break,
                }
            }
            let refused = refused.unwrap_or_else(|| panic!("{what}: damage must be refused"));
            for again in 0..5 {
                assert!(
                    stream.next_record().is_none(),
                    "{what}: after {refused}, ask {again} handed something back"
                );
            }
            assert_eq!(
                stream.buffered_bytes(),
                READ_CHUNK_BYTES,
                "{what}: a refused reader releases what it grew and keeps only its budget"
            );
        }
    }

    /// **A block whose length prefix claims more than its frame used is refused.** On a run of
    /// blocks that length reaches into the *next* block, so a reader that believed it would hand
    /// back records from a block nobody asked for — and the guard for it was held by nothing:
    /// dropping it left all 66 tests green while a two-block file read back as one.
    #[test]
    fn a_block_declaring_more_bytes_than_its_frame_is_refused() {
        let records: Vec<_> = (0..12).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(60));
        assert!(on_disk.block_offsets.len() >= 2, "two blocks at least");

        let whole = stream_every_record(&on_disk.bytes).expect("the honest file reads");
        assert_eq!(built(&whole), records);

        // Make the first block's length cover the second block too.
        let CompressedBlockAt::Whole { block_bytes, .. } = compressed_block_at(&on_disk.bytes)
        else {
            panic!("the first block is whole");
        };
        assert!(
            block_bytes < on_disk.bytes.len(),
            "the first block must not be the whole file"
        );
        let mut lying = on_disk.bytes.clone();
        lying[..COMPRESSED_BLOCK_LENGTH_BYTES].copy_from_slice(
            &u32::try_from(on_disk.bytes.len() - COMPRESSED_BLOCK_LENGTH_BYTES)
                .expect("a small file")
                .to_le_bytes(),
        );

        let refused = stream_every_record(&lying)
            .expect_err("a block whose length reaches into the next one is damaged");
        assert!(
            matches!(
                refused,
                BlockReadError::BlockDeclaresMoreBytesThanItsFrame { .. }
            ),
            "got {refused}"
        );
    }

    /// **`buffered_bytes` is the two buffers and nothing else.** The memory property rests on
    /// it, and a stub returning zero satisfied every check the module made.
    #[test]
    fn buffered_bytes_is_the_two_buffers() {
        let manifest = a_manifest();
        let fresh = BlockStream::new([].as_slice(), &manifest).expect("a manifest");
        assert_eq!(fresh.buffered_bytes(), READ_CHUNK_BYTES + ROLLING_BYTES);

        let (_, on_disk) = a_file_of_several_read_chunks();
        let mut stream = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let _ = stream.next_record().expect("a record").expect("it reads");
        assert_eq!(
            stream.buffered_bytes(),
            READ_CHUNK_BYTES + ROLLING_BYTES,
            "a reader mid-file holds exactly what it budgeted for"
        );
    }

    /// **A source that hands back one byte at a time reads the same records**, and one
    /// interrupted mid-read is retried rather than reported as a damaged store. Every other
    /// fixture is a `&[u8]`, which returns everything asked for and never fails — so neither
    /// path had a witness.
    #[test]
    fn a_dribbling_or_interrupted_source_reads_the_same_records() {
        let records: Vec<_> = (0..60).map(|index| a_record(0, 100 + index, 1)).collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(120));
        assert!(on_disk.block_offsets.len() >= 3);
        let manifest = a_manifest();

        for (what, most_bytes_a_read, interrupt_every) in [
            ("one byte at a time", 1usize, 0usize),
            ("seven bytes at a time", 7, 0),
            ("interrupted every third read", 64, 3),
            ("one byte at a time, interrupted every other read", 1, 2),
        ] {
            let source = DribblingSource::interrupted(
                on_disk.bytes.clone(),
                most_bytes_a_read,
                interrupt_every,
            );
            let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
            let mut back = Vec::new();
            while let Some(next) = stream.next_record() {
                let met = next.unwrap_or_else(|refused| panic!("{what}: {refused}"));
                back.push(met.record.expect("every record is built"));
            }
            assert_eq!(back, records, "{what}");
        }
    }

    /// **A record retried after a refill is parsed from its first byte, against state the retry
    /// did not touch.** This is spec §8's named trap — *"a parse that half-advances that state
    /// before failing corrupts every record after it, plausibly"* — and it is the property
    /// Milestone D4 exists to prove. Putting the defect in advances every later coordinate by
    /// one, which only a walk that compares the records it built can see.
    ///
    /// The source hands back one byte at a time, so **every record in the file is retried**,
    /// most of them many times.
    #[test]
    fn a_record_retried_after_a_refill_is_the_record_that_was_written() {
        let records: Vec<_> = (1..900u64).map(|at| a_record(0, at, 1)).collect();
        let on_disk = blocks_on_disk(&records, Bp(300), None);
        assert!(on_disk.block_offsets.len() >= 3, "several blocks");

        let manifest = a_manifest();
        let source = DribblingSource::new(on_disk.bytes.clone(), 1);
        let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
        let mut back = Vec::new();
        while let Some(next) = stream.next_record() {
            back.push(next.expect("it reads").record.expect("built"));
        }
        assert_eq!(
            back, records,
            "a retried record is the record that was written"
        );
    }

    /// **A block nothing can be parsed out of costs a bounded amount of memory, not its own
    /// decompressed size.**
    ///
    /// A psp block's decompressed size is not bounded by its size on disk. Before this ceiling
    /// existed a review agent built one that is **4,132 bytes on disk and drove a reader to hold
    /// 67,125,248 bytes** — 2,048 times the 32 kB the two buffers hold between them, and 131
    /// times what spec §1.1 gives a whole open sample — because no record parsed, so nothing
    /// ever "fitted", and the buffer doubled until the frame ran out. At a thousand open samples
    /// that is a run that dies on memory rather than one that reports a bad file.
    ///
    /// The fixture is the same shape: a block whose records are one long run of a single byte,
    /// which compresses to almost nothing and inflates far past the ceiling, with a record head
    /// in front of it declaring a body longer than anything that follows.
    ///
    /// **Its size is fixed rather than a multiple of the ceiling.** Written as `ceiling × 8` it
    /// would grow with the very knob whose purpose is to bound memory — at a raised ceiling of
    /// 1 GiB the fixture alone would be 8 GB.
    #[test]
    fn a_block_that_never_parses_costs_a_bounded_amount_of_memory() {
        const INFLATED_BYTES: usize = 8 * 1024 * 1024;
        const {
            assert!(
                INFLATED_BYTES > ROLLING_BUFFER_CEILING_BYTES * 4,
                "the fixture must inflate well past the ceiling to say anything — and it is a \
                 fixed 8 MiB, so raising the ceiling past 2 MiB means rethinking this test \
                 rather than rescaling it"
            )
        };
        let mut payload = Vec::with_capacity(INFLATED_BYTES);
        BlockHead {
            contig: ContigId(4),
            first_position: Position(90_000),
            record_count: NonZeroU64::MIN,
        }
        .encode(&mut payload);
        // position-offset 0, reference-span 1, non-reference-reads 0, and a body length larger
        // than anything that follows.
        payload.extend_from_slice(&[0x00, 0x01, 0x00]);
        encode_u64_leb128(u64::from(u32::MAX), &mut payload);
        // No chain-id departures and no arrivals, so the head is whole and what the reader
        // cannot find the end of is the body.
        payload.extend_from_slice(&[0x00, 0x00]);
        payload.resize(INFLATED_BYTES, b'A');

        let manifest = a_manifest();
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let on_disk = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();
        assert!(
            on_disk.len() * 64 < payload.len(),
            "the fixture must inflate far beyond its size on disk: {} bytes on disk against {} \
             inflated",
            on_disk.len(),
            payload.len()
        );

        let mut stream = BlockStream::new(on_disk.as_slice(), &manifest).expect("a manifest");
        let refused = match stream.next_record() {
            Some(Ok(_)) => panic!("nothing in this block is a record"),
            Some(Err(refused)) => refused,
            None => panic!("a block that never parses must be refused, not ended"),
        };

        // **What the refusal reports, not what the reader holds afterwards** — by then `fail`
        // has released the buffer, so a reading taken after the fact is the budget whatever
        // happened before it, and an assertion on it could not fail.
        let BlockReadError::RecordLargerThanTheReaderAllows {
            block,
            records_read,
            allowed_bytes,
        } = refused
        else {
            panic!("got {refused}");
        };
        assert_eq!(allowed_bytes, ROLLING_BUFFER_CEILING_BYTES);
        assert_eq!(
            block.contig,
            ContigId(4),
            "the refusal locates the block the record was in"
        );
        assert_eq!(block.start, Position(90_000));
        assert_eq!(records_read, 0, "it was the block's first record");
        assert_eq!(
            stream.buffered_bytes(),
            READ_CHUNK_BYTES,
            "and it releases what it grew once the block is refused"
        );
    }

    /// **⚠ The risk a ceiling creates: it must not refuse a file a writer legitimately produced.**
    ///
    /// At the shipped 100 kb grid and three hundred reads a position a block is about 1.76 MB
    /// decompressed (spec §4.1's 17.557 bytes a record), which is **larger than the ceiling** —
    /// and that has to be fine, because the buffer grows for one *record*, not for a block. No
    /// reader test decoded a block past the ceiling before this one, and the gap was live:
    /// checking the buffer's *capacity* rather than its length passes all the other tests while
    /// refusing a well-formed 2 MB block.
    #[test]
    fn a_well_formed_block_larger_than_the_ceiling_reads_back_every_record() {
        let records: Vec<_> = (1..14_000u64).map(an_incompressible_record).collect();
        let payload = one_payload(&records);
        assert!(
            payload.len() > ROLLING_BUFFER_CEILING_BYTES * 2,
            "the block must be well past the ceiling: {} bytes against {}",
            payload.len(),
            ROLLING_BUFFER_CEILING_BYTES
        );

        let manifest = a_manifest();
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let on_disk = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let mut stream = BlockStream::new(on_disk.as_slice(), &manifest).expect("a manifest");
        let mut back = Vec::new();
        let mut most_held = stream.buffered_bytes();
        while let Some(next) = stream.next_record() {
            back.push(
                next.expect("a well-formed block reads")
                    .record
                    .expect("built"),
            );
            most_held = most_held.max(stream.buffered_bytes());
        }
        assert_eq!(back, records);
        assert!(
            most_held <= READ_CHUNK_BYTES + ROLLING_BYTES,
            "a block {} bytes decompressed was read holding {most_held}",
            payload.len()
        );
    }

    /// One record whose encoded body is a given size, near enough: `observations` reads of
    /// sixteen bases each, which the encoder writes out one after another.
    fn a_wide_record(start: u64, observations: u32) -> SampleLocusObservations {
        let mut wide = a_record(0, start, 1);
        wide.observations = (0..observations)
            .map(|read| SequenceObservation {
                bases: vec![b"ACGT"[(read % 4) as usize]; 16].into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read % 7),
                num_obs: 1,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(i64::from(read) + 1)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 1,
                chain_ids: Vec::new(),
            })
            .collect();
        wide
    }

    /// **The ceiling counts what the buffer *holds*, not how large it has grown.** Once one
    /// record has pushed the rolling buffer's capacity up to the ceiling, every later record of
    /// the same block still has to be read — the buffer is then a large empty vector, and what it
    /// holds is one record's worth.
    ///
    /// This is the one mutation that survived the whole D4/D5 suite: `self.rolling.capacity() >=
    /// self.buffer_ceiling` in place of `self.rolling.len() >= …`. Under it a well-formed file
    /// whose records each need more than half the ceiling is refused with *"a record needs more
    /// than the … bytes this reader allows one to hold"* — a valid file reported as a record too
    /// large for the reader, which is the refusal that sends an operator to raise a limit that was
    /// never binding.
    #[test]
    fn a_block_whose_records_each_grow_the_buffer_to_the_ceiling_still_reads() {
        // Sized so one record needs more than half the ceiling: the buffer doubles from
        // ROLLING_BYTES, so holding one takes a capacity of exactly the ceiling.
        let records: Vec<_> = (1..=3u64)
            .map(|at| a_wide_record(at * 10, 12_000))
            .collect();
        let payload = one_payload(&records);
        let a_record_holds = payload.len() / records.len();
        assert!(
            (ROLLING_BUFFER_CEILING_BYTES / 2..ROLLING_BUFFER_CEILING_BYTES)
                .contains(&a_record_holds),
            "each record must be over half the ceiling and under it, or the buffer's capacity \
             never reaches the ceiling and the mutation this test exists for cannot fire: \
             {a_record_holds} bytes against a ceiling of {ROLLING_BUFFER_CEILING_BYTES}"
        );

        let manifest = a_manifest();
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let on_disk = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let mut stream = BlockStream::new(on_disk.as_slice(), &manifest).expect("a manifest");
        let mut back = Vec::new();
        while let Some(next) = stream.next_record() {
            back.push(
                next.unwrap_or_else(|refused| {
                    panic!("a well-formed record was refused: {refused}")
                })
                .record
                .expect("built"),
            );
        }
        assert_eq!(back, records);
    }

    /// A record genuinely larger than the ceiling is read by raising it, which is what makes a
    /// ceiling tolerable at all — and a ceiling under the buffer it caps is refused rather than
    /// accepted as a way of turning every record away.
    #[test]
    fn a_record_past_the_ceiling_is_read_by_raising_it() {
        let mut enormous = a_record(0, 500, 1);
        enormous.observations = (0..40_000u32)
            .map(|read| SequenceObservation {
                bases: vec![b"ACGT"[(read % 4) as usize]; 16].into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read % 7),
                num_obs: 1,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(i64::from(read) + 1)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 1,
                chain_ids: Vec::new(),
            })
            .collect();
        let records = vec![enormous];
        let payload = one_payload(&records);
        assert!(
            payload.len() > ROLLING_BUFFER_CEILING_BYTES,
            "one record must exceed the ceiling: {} bytes against {}",
            payload.len(),
            ROLLING_BUFFER_CEILING_BYTES
        );

        let manifest = a_manifest();
        let mut compressor = BlockCompressor::from_manifest(&manifest).expect("a window");
        let on_disk = compressor
            .compress(&payload)
            .expect("it compresses")
            .to_vec();

        let mut refused = BlockStream::new(on_disk.as_slice(), &manifest).expect("a manifest");
        assert!(matches!(
            refused.next_record(),
            Some(Err(BlockReadError::RecordLargerThanTheReaderAllows { .. }))
        ));

        let mut raised = BlockStream::new(on_disk.as_slice(), &manifest)
            .expect("a manifest")
            .with_a_buffer_ceiling(payload.len() * 2)
            .expect("a ceiling above the buffer");
        assert_eq!(raised.buffer_ceiling(), payload.len() * 2);
        let met = raised.next_record().expect("a record").expect("it reads");
        assert_eq!(met.record.as_ref(), Some(&records[0]));
        assert!(raised.next_record().is_none());

        let refused = BlockStream::new([].as_slice(), &manifest)
            .expect("a manifest")
            .with_a_buffer_ceiling(ROLLING_BYTES)
            .expect_err("a ceiling at the buffer refuses records the buffer already holds");
        assert!(
            matches!(refused, BlockReadError::BufferCeilingUnderTheBuffer { .. }),
            "got {refused}"
        );
    }

    /// The refusal names the number to raise, because raising it is the fix when the record is
    /// genuine — the same instruction spec §4.2 attaches to a look-back window wider than a
    /// reader budgeted for.
    #[test]
    fn a_record_past_the_readers_ceiling_names_the_number_to_raise() {
        let refused = BlockReadError::RecordLargerThanTheReaderAllows {
            block: GenomeRegion {
                contig: ContigId(7),
                start: Position(90_600_000),
                end: Position(90_600_000),
            },
            records_read: 41,
            allowed_bytes: 524_288,
        };
        assert_eq!(
            refused.to_string(),
            "a record in the block starting at contig 7:90600000-90600000 needs more than the \
             524288 bytes this reader allows one record to hold; raise the ceiling to read it"
        );
    }

    /// Every refusal says what it was and carries the number whoever sees it must act on.
    #[test]
    fn a_readers_refusals_name_what_broke() {
        assert_eq!(
            BlockReadError::RecordRunsPastItsBlock {
                records_left: 2,
                bytes_left: 7,
            }
            .to_string(),
            "a record runs past the end of its block, with 2 records still declared and 7 bytes \
             left"
        );
        assert_eq!(
            BlockReadError::FileEndsInsideABlock {
                compressed_bytes_left: 91,
            }
            .to_string(),
            "the file ends with 91 compressed bytes of a block still due"
        );
        assert_eq!(
            BlockReadError::WindowLogOutOfRange {
                look_back_window_log: 40,
            }
            .to_string(),
            format!(
                "the file declares a look-back window of 2^40 bytes; zstd takes between \
                 2^{MIN_LOOK_BACK_WINDOW_LOG} and 2^{MAX_LOOK_BACK_WINDOW_LOG}"
            )
        );

        let damaged = BlockReadError::from_record(RecordDecodeError::Malformed {
            field: "observation-count",
            bytes_in: 6,
            reason: "a count larger than the body".to_string(),
        });
        assert!(matches!(damaged, BlockReadError::DamagedRecord { .. }));
        assert!(std::error::Error::source(&damaged).is_some());

        let upgrade = BlockReadError::from_record(RecordDecodeError::Unsupported {
            field: "locus-kind",
            bytes_in: 6,
            tag: 9,
        });
        assert!(
            matches!(upgrade, BlockReadError::UnsupportedRecord { .. }),
            "a tag a later writer added is *upgrade the reader*, not damage"
        );
    }

    // -----------------------------------------------------------------
    // The chain ids ride in the head, so a skipped body cannot strand them
    // -----------------------------------------------------------------

    /// Coverage that looks like read pairs, as records a block builder will take: one record a
    /// position, each naming the reads live there, with the mates of a pair separated by an
    /// unsequenced hole so identifiers go live twice.
    fn records_naming_paired_reads(positions: u64) -> Vec<SampleLocusObservations> {
        const MATE: u64 = 12;
        const HOLE: u64 = 9;
        const PAIRS_A_POSITION: u64 = 2;
        let span = MATE * 2 + HOLE;
        (0..positions)
            .map(|at| {
                let live: Vec<ChainId> = (at.saturating_sub(span - 1)..=at)
                    .filter(|start| {
                        let into_the_pair = at - start;
                        into_the_pair < MATE || (MATE + HOLE..span).contains(&into_the_pair)
                    })
                    .flat_map(|start| {
                        (0..PAIRS_A_POSITION).map(move |which| start * PAIRS_A_POSITION + which)
                    })
                    .collect();
                let mut record = a_record(0, at + 1, 1);
                record.observations[0].chain_ids = live;
                record
            })
            .collect()
    }

    /// **A record the codec refuses at a block cut leaves the open block's live set alone.**
    ///
    /// ⚠ **This is a hazard the chain ids created and the shape of the writer removed.** The cut
    /// path used to reset the encoder to the new block, try the record, and — on a refusal — put
    /// the coordinate base back. A coordinate base can be put back; **a live set cannot**, and
    /// the open block still needed it. So `encode_record_starting_a_block` makes every refusal
    /// before it resets anything, and there is nothing to put back.
    ///
    /// The evidence is the bytes: a run written with refusals interleaved must equal the run
    /// written without them, over records that actually name reads — with empty chain-id lists
    /// every record's changes are two zero bytes and a lost set would not show.
    #[test]
    fn a_refusal_at_a_cut_leaves_the_open_blocks_live_reads_alone() {
        let records = records_naming_paired_reads(120);
        let clean =
            cut(BlockBuilder::new(Bp(20), None).expect("a grid"), &records).expect("in order");
        assert!(
            clean.len() >= 3,
            "the fixture has to cut for this to say anything"
        );

        let mut interrupted = BlockBuilder::new(Bp(20), None).expect("a grid");
        let mut with = Vec::new();
        for record in &records {
            // A region covering no base: refused by the codec, and at a cell of its own so the
            // refusal comes through the cut path.
            assert!(matches!(
                interrupted.push(&a_record_over_no_base(0, 900_000)),
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

        assert_eq!(
            with, clean,
            "a refused record at a cut must leave the open block's reads exactly as they were"
        );
    }

    /// **A record that straddles the reader's buffer is retried, and its reads come back once.**
    ///
    /// ⚠ **This is the regression test for a defect that shipped in this milestone twice, one
    /// level apart.** Answering a short read means *fetch more bytes and re-parse this record
    /// from its first byte*, so a parse that moved the live set before it could still fail meets
    /// its own changes a second time: an arrival into a set it is already in, refused as damage
    /// on a perfectly good file; or a departure position resolved against an already-shrunk set,
    /// which removes a different read and says nothing at all.
    ///
    /// **The suite could not see either.** Every fixture that reaches the retry path had empty
    /// chain-id lists, so each record's changes were the two bytes `0, 0` and applying them twice
    /// is applying them never. This one names reads, and its blocks are larger than the 16 kB
    /// rolling buffer, so records really do straddle it — asserted, not assumed.
    #[test]
    fn records_that_straddle_the_buffer_name_their_reads_exactly_once() {
        let records: Vec<_> = (1..1_400u64)
            .map(|at| {
                let mut record = an_incompressible_record(at);
                record.observations[0].chain_ids =
                    (at.saturating_sub(19)..=at).map(|id| id * 3).collect();
                record
            })
            .collect();
        let on_disk = blocks_on_disk(&records, Bp(400), None);
        let biggest = on_disk
            .payloads
            .iter()
            .map(Vec::len)
            .max()
            .expect("some blocks");
        assert!(
            biggest > ROLLING_BYTES * 2,
            "the blocks must be larger than the buffer, or nothing straddles it: {biggest} \
             against {ROLLING_BYTES}"
        );

        let manifest = a_manifest();
        // What a reader handed the whole file at once has live at each record.
        let mut whole = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let mut live_at_each = Vec::new();
        while let Some(next) = whole.next_record() {
            let _ = next.expect("a well-formed file reads");
            live_at_each.push(whole.live_reads().ids().to_vec());
        }
        assert_eq!(live_at_each.len(), records.len());

        // And the same file dribbled in, so records are retried on the way.
        let mut restarts_at_a_byte_a_read = 0u64;
        for most_bytes_a_read in [1usize, 3, 1_024] {
            let source = DribblingSource::new(on_disk.bytes.clone(), most_bytes_a_read);
            let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
            let mut at = 0usize;
            while let Some(next) = stream.next_record() {
                let _ = next.unwrap_or_else(|refused| {
                    panic!(
                        "at {most_bytes_a_read} bytes a read, record {at} of a well-formed file \
                         was refused: {refused}"
                    )
                });
                assert_eq!(
                    stream.live_reads().ids(),
                    live_at_each[at].as_slice(),
                    "at {most_bytes_a_read} bytes a read, record {at}'s reads"
                );
                at += 1;
            }
            assert_eq!(at, records.len());
            assert!(
                stream.parses_restarted() > 0,
                "at {most_bytes_a_read} bytes a read nothing was retried, so this schedule says \
                 nothing about a record that straddles the buffer"
            );
            if most_bytes_a_read == 1 {
                restarts_at_a_byte_a_read = stream.parses_restarted();
            }
        }
        assert!(
            restarts_at_a_byte_a_read > records.len() as u64,
            "a byte a read has to retry more often than there are records for this to be about \
             retries at all: {restarts_at_a_byte_a_read} over {} records",
            records.len()
        );
    }

    /// **A walk that skips almost every body still knows which reads are live.**
    ///
    /// This is the whole of why the chain ids' changes are in the record *head*. They carry
    /// state: a reader knows which reads are live only because it applied every arrival and
    /// departure since the block began. **Kept in the body, a reader that declined a record
    /// would never see that record's changes, so its set would go stale and every later record
    /// it did want would be wrong** — silently, because a stale set is still a plausible set
    /// (spec `psp_record_encoding.md` §6).
    ///
    /// So the oracle is: take one record in seven, and at each one the live set must be what a
    /// walk that built every record has there. The fixture is paired-end coverage, so most
    /// identifiers go live twice and a stale set drifts rather than merely lagging.
    #[test]
    fn a_walk_that_skips_almost_every_body_knows_which_reads_are_live() {
        let records = records_naming_paired_reads(300);
        let on_disk = blocks_on_disk(&records, Bp(90), None);
        assert!(
            on_disk.payloads.len() >= 3,
            "several blocks, so the set restarts inside the walk as well"
        );
        let manifest = a_manifest();

        // What a reader that builds everything has live at each record.
        let mut whole = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let mut live_at_each = Vec::new();
        while let Some(next) = whole.next_record() {
            let _ = next.expect("it reads");
            live_at_each.push(whole.live_reads().ids().to_vec());
        }
        assert_eq!(live_at_each.len(), records.len());
        assert!(
            live_at_each.iter().any(|live| live.len() > 20),
            "the fixture has to put reads live for this to say anything: the most at once was {}",
            live_at_each.iter().map(Vec::len).max().unwrap_or(0)
        );

        // And a reader that declines six records in seven.
        let mut skipping =
            BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let mut met = 0usize;
        let mut built = 0usize;
        while let Some(next) = skipping.next_record_where(|_| met.is_multiple_of(7)) {
            let streamed = next.expect("it reads");
            assert_eq!(
                skipping.live_reads().ids(),
                live_at_each[met].as_slice(),
                "record {met}: a skipping walk must have the reads a full walk has"
            );
            if streamed.record.is_some() {
                built += 1;
            }
            met += 1;
        }
        assert_eq!(met, records.len());
        assert!(
            built * 6 < met,
            "the walk has to skip most of the file to say anything: it built {built} of {met}"
        );
    }

    /// **And a reader that starts at an arbitrary block, skipping, agrees from there.**
    ///
    /// Milestone F's index hands a reader a block offset and nothing else. The two rules meet
    /// here: the set restarts at the block, and the skipping reader still applies every head.
    #[test]
    fn a_skipping_walk_from_any_block_has_the_reads_a_full_walk_has_there() {
        let records = records_naming_paired_reads(200);
        let on_disk = blocks_on_disk(&records, Bp(60), None);
        let manifest = a_manifest();

        let mut whole = BlockStream::new(on_disk.bytes.as_slice(), &manifest).expect("a manifest");
        let mut live_at_each = Vec::new();
        while let Some(next) = whole.next_record() {
            let _ = next.expect("it reads");
            live_at_each.push(whole.live_reads().ids().to_vec());
        }

        let mut skipped_records = 0usize;
        for (block, offset) in on_disk.block_offsets.iter().enumerate() {
            let mut from_here =
                BlockStream::new(&on_disk.bytes[*offset..], &manifest).expect("a manifest");
            let mut at = skipped_records;
            while let Some(next) = from_here.next_record_where(|_| false) {
                let _ = next.expect("it reads");
                assert_eq!(
                    from_here.live_reads().ids(),
                    live_at_each[at].as_slice(),
                    "block {block}, record {at}: starting here must give what reading through gives"
                );
                at += 1;
            }
            assert_eq!(
                at,
                records.len(),
                "block {block} read to the end of the file"
            );
            skipped_records += BlockRecords::split(&on_disk.payloads[block])
                .expect("the block head reads")
                .head
                .record_count
                .get() as usize;
        }
    }

    // -----------------------------------------------------------------
    // The restartable parse
    // -----------------------------------------------------------------

    /// Read a whole file through a source that hands over exactly `most_bytes_a_read` bytes a read, and say
    /// how many times a record's parse had to be retried on the way.
    fn stream_through_a_source_yielding_at_most(
        bytes: &[u8],
        most_bytes_a_read: usize,
    ) -> (Vec<SampleLocusObservations>, u64) {
        let manifest = a_manifest();
        let source = DribblingSource::new(bytes.to_vec(), most_bytes_a_read);
        let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
        let mut back = Vec::new();
        while let Some(next) = stream.next_record() {
            back.push(
                next.unwrap_or_else(|refused| {
                    panic!("at {most_bytes_a_read} bytes a read: {refused}")
                })
                .record
                .expect("every record is built"),
            );
        }
        (back, stream.parses_restarted())
    }

    /// **The oracle: a decode forced to refill at every possible boundary.**
    ///
    /// A reader meets the file as a stream of compressed bytes, so a refill can fall between any
    /// two of them — and a source that hands over `n` bytes a read puts the refills at every
    /// multiple of `n`. Running `n` from one byte to the whole file therefore places a refill at
    /// **every byte offset the reader can meet one at**, and every run must give the same
    /// records as a single-shot read.
    ///
    /// **⚠ And on blocks this small it retries nothing**, which is a fact about zstd rather than
    /// about the reader: it decodes in internal blocks and emits one whole, so a payload that
    /// fits a single emission is delivered in one piece however slowly its input arrived. The
    /// count is asserted below rather than left implicit, because the sibling test exists
    /// entirely because of it — and because a docstring here once claimed the opposite.
    ///
    /// What this sweep does hold is every *compressed-side* alignment: it was among the killers
    /// of three mutations, which is why it stays.
    #[test]
    fn a_decode_refilled_at_every_source_byte_boundary_gives_the_same_records() {
        // Small enough to run a schedule per byte, varied enough that records differ in width.
        let records: Vec<_> = (0..40)
            .map(|index| a_record(0, 100 + index * 3, 1 + index % 4))
            .collect();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(90));
        assert!(
            on_disk.block_offsets.len() >= 3,
            "several blocks, so refills land inside block heads as well as inside records"
        );

        let (whole, _) = stream_through_a_source_yielding_at_most(&on_disk.bytes, usize::MAX);
        assert_eq!(whole, records, "a single-shot read is the thing to match");
        let mut ever_retried = 0u64;
        for most_bytes_a_read in 1..=on_disk.bytes.len() {
            let (back, retries) =
                stream_through_a_source_yielding_at_most(&on_disk.bytes, most_bytes_a_read);
            assert_eq!(back, records, "reading {most_bytes_a_read} bytes a read");
            ever_retried += retries;
        }
        assert_eq!(
            ever_retried, 0,
            "blocks that fit one zstd emission are never straddled, whatever the input \
             schedule — the fact the next test exists for"
        );
    }

    /// **And the same sweep over blocks the rolling buffer cannot hold, which is where a record
    /// actually straddles one.**
    ///
    /// The sweep above moves the point at which *compressed* bytes arrive, and on small blocks
    /// that turns out not to move where a *record* is cut in half: zstd decodes in internal
    /// blocks and emits one whole, so a payload that fits a single emission is delivered in one
    /// piece however slowly its input arrived. Measured on the fixture above — **837 schedules,
    /// none of which retried a record even once**. A test that stopped there would have proved
    /// the walk works and nothing about restarting it.
    ///
    /// What straddles a record is the buffer running out, so the blocks here are larger than it.
    ///
    /// **⚠ The first version of this fixture was a single block.** Its 1,999 records all fell in
    /// cell 0 of the grid it named, so the sweep retried tens of thousands of times and never
    /// crossed a block boundary while doing it — and the test's name, and D4's report, were both
    /// plural about one block. The grid below cuts, and `payloads.len() >= 3` is what says so.
    ///
    /// Measured on the fixture as it now stands — three blocks, the largest 73,849 bytes against
    /// a 16 kB buffer:
    ///
    /// | bytes a read | records retried |
    /// |---:|---:|
    /// | 1 | 37,209 |
    /// | 17 | 2,200 |
    /// | 1,024 | 48 |
    /// | the whole file at once | 14 |
    ///
    /// **Even the single-shot read retries fourteen times**, because a payload larger than the
    /// buffer straddles it however fast the input arrives. That is the whole difference between
    /// this sweep and the one above.
    #[test]
    fn a_decode_of_blocks_larger_than_the_buffer_is_retried_and_still_exact() {
        let records: Vec<_> = (1..2_000u64).map(an_incompressible_record).collect();
        // A grid that cuts, so the sweep crosses block boundaries as well as buffer ones — the
        // first version put all 1,999 records in one cell and never crossed a block at all.
        let on_disk = blocks_on_disk(&records, Bp(700), None);
        let biggest = on_disk
            .payloads
            .iter()
            .map(|payload| payload.len())
            .max()
            .expect("some blocks");
        assert!(
            on_disk.payloads.len() >= 3,
            "several blocks, so a refill and a cut are different boundaries here"
        );
        assert!(
            biggest > ROLLING_BYTES * 2,
            "the blocks must be larger than the buffer, or no record straddles one: {biggest} \
             bytes against {ROLLING_BYTES}"
        );

        let mut most_retries = 0u64;
        for most_bytes_a_read in [1usize, 2, 3, 17, 251, 1024, READ_CHUNK_BYTES, usize::MAX] {
            let (back, retries) =
                stream_through_a_source_yielding_at_most(&on_disk.bytes, most_bytes_a_read);
            assert_eq!(back, records, "reading {most_bytes_a_read} bytes a read");
            most_retries = most_retries.max(retries);
        }
        // **The schedules have to actually retry, or the loop above proves nothing.**
        assert!(
            most_retries > records.len() as u64 / 10,
            "the sweep retried at most {most_retries} times over {} records, too few to have \
             crossed many of them",
            records.len()
        );
    }

    /// **A record larger than the rolling buffer is read at every refill schedule too**, and the
    /// buffer that grew for it goes back to what a reader budgets for at the next block.
    ///
    /// This is the case spec §8 refuses to bake a ceiling for — "many alleles, many chain ids" —
    /// met at the boundaries where it is hardest: a record that does not fit is retried until it
    /// does, so every retry is one the buffer has to survive.
    #[test]
    fn a_record_larger_than_the_buffer_survives_five_refill_schedules() {
        let mut enormous = a_record(0, 500, 1);
        enormous.observations = (0..2_600u32)
            .map(|read| SequenceObservation {
                bases: vec![b"ACGT"[(read % 4) as usize]; 12].into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read % 7),
                num_obs: 1,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(i64::from(read) + 1)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 1,
                chain_ids: Vec::new(),
            })
            .collect();
        let records = vec![enormous, a_record(0, 200_000, 1), a_record(0, 200_001, 2)];

        let payloads =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");
        assert!(
            payloads[0].len() > ROLLING_BYTES,
            "the record must exceed the rolling buffer: {} bytes against {ROLLING_BYTES}",
            payloads[0].len()
        );

        let on_disk = blocks_on_disk(&records, A_GRID, None);
        for most_bytes_a_read in [1usize, 3, 64, 1024, usize::MAX] {
            let (back, _) =
                stream_through_a_source_yielding_at_most(&on_disk.bytes, most_bytes_a_read);
            assert_eq!(back, records, "reading {most_bytes_a_read} bytes a read");
        }

        // And the buffer comes back down: the last block's records are ordinary.
        let manifest = a_manifest();
        let source = DribblingSource::new(on_disk.bytes.clone(), 7);
        let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
        while let Some(next) = stream.next_record() {
            let _ = next.expect("it reads");
        }
        assert_eq!(
            stream.buffered_bytes(),
            READ_CHUNK_BYTES + ROLLING_BYTES,
            "the buffer that grew for one record is back to the budget"
        );
    }

    /// **A block head is read after however many refills it takes, and comes out whole.**
    ///
    /// ⚠ *This test was called `a_refill_inside_a_block_head_is_retried` and claimed a one-byte
    /// schedule "refills inside every one of" the head's three integers. It does not, and it
    /// cannot: a reader clears its buffer at a block start and zstd's first emission delivers
    /// the whole internal block, so a refill always lands **before** a head and never inside
    /// one. Measured over 108,746 head restarts, the largest partial head ever seen was zero
    /// bytes.* What is real, and what this holds, is that the head's parse is **restarted** —
    /// many times — and that the fields it finally reads are the ones the writer put there.
    #[test]
    fn a_block_head_is_read_after_however_many_restarts_it_takes() {
        // A contig, first position and record count that each need several bytes, so a
        // one-byte-a-read schedule stops inside each of them.
        let records: Vec<_> = (0..300)
            .map(|index| a_record(300, 90_600_000 + index, 1))
            .collect();
        let on_disk = blocks_on_disk(&records, A_GRID, None);
        let head = walk(&on_disk.payloads[0]).head;
        assert_eq!(head.contig, ContigId(300));
        assert!(head.first_position.get() > 1 << 21, "a four-byte position");
        assert!(head.record_count.get() > 0x7f, "a multi-byte count");

        let manifest = a_manifest();
        for most_bytes_a_read in [1usize, 2, 3] {
            let source = DribblingSource::new(on_disk.bytes.clone(), most_bytes_a_read);
            let mut stream = BlockStream::new(source, &manifest).expect("a manifest");
            let mut back = Vec::new();
            while let Some(next) = stream.next_record() {
                back.push(next.expect("it reads").record.expect("built"));
            }
            assert_eq!(back, records, "reading {most_bytes_a_read} bytes a read");
            // **The witness.** Without it this is a walk that came out right, which is what the
            // test claimed to be about and was not.
            assert!(
                stream.block_heads_restarted() > 0,
                "reading {most_bytes_a_read} bytes a read never restarted a block head, so the schedule \
                 proves nothing about restarting one"
            );
        }
    }

    // -----------------------------------------------------------------
    // Every running difference resets at a block boundary
    // -----------------------------------------------------------------

    /// Records over three contigs and several grid cells, with the cells **far apart** — so a
    /// difference that was not reset lands a record about **86 million bases** from where it
    /// belongs (cells 0, 7, 41 and 900 of a 100 kb grid) rather than one base out, and a
    /// comparison sees it.
    fn records_whose_blocks_are_far_apart() -> Vec<SampleLocusObservations> {
        let mut records = Vec::new();
        for contig in 0..3u32 {
            for cell in [0u64, 7, 41, 900] {
                for step in 0..9u64 {
                    records.push(a_record(
                        contig,
                        cell * A_GRID.get() + 1 + step * 4_000,
                        1 + step % 5,
                    ));
                }
            }
        }
        records
    }

    /// **A block read on its own gives exactly what it gives in the middle of a file.**
    ///
    /// Spec §3.2: *"every running difference inside it restarts — the position difference, the
    /// coverage difference, the chain-id difference"*. A reader handed one block's bytes and
    /// nothing else has no history at all, so any difference still being measured from a previous
    /// block shows up as records at the wrong coordinates — and **plausibly** wrong, because
    /// coverage is smooth and a position difference that is slightly off still parses.
    ///
    /// **What it adds over D3's `a_reader_starting_at_a_block_gets_the_tail_of_a_full_read`**,
    /// which reads every *tail* of the same file: a byte range that both begins **and ends** at a
    /// block boundary — the shape Milestone F's index will hand a reader — and a fixture whose
    /// blocks are about 86 million bases apart, so a carried difference is not a small coordinate
    /// error. ⚠ It is not a strictly stronger test: a review could build no mutation that kills
    /// this one and spares D3's.
    ///
    /// Today there is one running difference, the position offset. **Milestone E adds the
    /// chain-id live set, which is the second**, and this test is what it will meet.
    #[test]
    fn a_block_read_alone_gives_what_it_gives_in_the_middle_of_a_file() {
        let records = records_whose_blocks_are_far_apart();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(150));
        assert!(
            on_disk.block_offsets.len() >= 12,
            "several blocks per contig, or a lone block is barely different from the file"
        );

        let whole = stream_every_record(&on_disk.bytes).expect("the file reads");
        assert_eq!(built(&whole), records);

        let mut taken = 0usize;
        for (index, offset) in on_disk.block_offsets.iter().enumerate() {
            let ends_at = on_disk
                .block_offsets
                .get(index + 1)
                .copied()
                .unwrap_or(on_disk.bytes.len());
            let alone = stream_every_record(&on_disk.bytes[*offset..ends_at])
                .unwrap_or_else(|refused| panic!("block {index} alone: {refused}"));

            let holds = usize::try_from(walk(&on_disk.payloads[index]).head.record_count.get())
                .expect("a small count");
            assert_eq!(
                alone,
                whole[taken..taken + holds],
                "block {index} read alone against the same block read in the file"
            );
            taken += holds;
        }
        assert_eq!(
            taken,
            records.len(),
            "every record was in exactly one block"
        );
    }

    /// **Restart equals sequential, from every block, at every refill schedule.**
    ///
    /// The plan's own oracle: reading from an arbitrary block gives what a full read gives from
    /// that point. D4 showed that where a reader refills is not where a writer cut, so the two
    /// boundaries are crossed together here — a stale difference that survived a block boundary
    /// and one that survived a buffer boundary are different defects.
    #[test]
    fn restart_equals_sequential_from_every_block_at_five_refill_schedules() {
        let records = records_whose_blocks_are_far_apart();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(150));
        let whole = stream_every_record(&on_disk.bytes).expect("the file reads");

        for most_bytes_a_read in [1usize, 5, 97, READ_CHUNK_BYTES, usize::MAX] {
            let mut taken = 0usize;
            for (index, offset) in on_disk.block_offsets.iter().enumerate() {
                let (from_here, _) = stream_through_a_source_yielding_at_most(
                    &on_disk.bytes[*offset..],
                    most_bytes_a_read,
                );
                assert_eq!(
                    from_here,
                    built(&whole[taken..]),
                    "from block {index} at {most_bytes_a_read} bytes a read"
                );
                taken += usize::try_from(walk(&on_disk.payloads[index]).head.record_count.get())
                    .expect("a small count");
            }
        }
    }

    /// **The writer's half: a block's bytes do not depend on the blocks written before it.**
    ///
    /// Self-containment has to hold on the way out as well as on the way in — a reader that
    /// restarts correctly is no use if the block it restarts at was written against a coordinate
    /// only the previous block knew. So every block is rebuilt from its own records alone and
    /// compared byte for byte with the one the whole run produced.
    #[test]
    fn a_block_is_written_the_same_alone_as_it_is_after_other_blocks() {
        let records = records_whose_blocks_are_far_apart();
        let on_disk = blocks_on_disk(&records, A_GRID, Some(150));

        let mut taken = 0usize;
        for (index, payload) in on_disk.payloads.iter().enumerate() {
            let walked = walk(payload);
            let holds = usize::try_from(walked.head.record_count.get()).expect("a small count");
            let its_own = &records[taken..taken + holds];

            // A builder that has seen nothing else, handed just this block's records. The grid
            // must not cut them further, so they are given a cell of their own.
            let alone = cut(
                BlockBuilder::new(Bp(u64::MAX), None).expect("a grid"),
                its_own,
            )
            .expect("in order");
            assert_eq!(
                alone.len(),
                1,
                "block {index}'s records are one block alone"
            );
            assert_eq!(
                alone[0], *payload,
                "block {index} written alone against written after {index} others"
            );
            taken += holds;
        }
        assert_eq!(taken, records.len());
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
