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
//! **What this file does not do yet.** Nothing here compresses (Milestone D2) and nothing
//! writes a file (Milestone F). [`BlockBuilder`] turns a stream of records into block
//! payloads — the bytes a compressor is handed — and stops there.
//!
//! Design authority: `doc/devel/ng/spec/psp_file_format.md` §3.2 (a block is self-contained),
//! §4.1 (the cut rule and what it buys), and `doc/devel/ng/arch/psp_file_format.md` §1.
//!
//! [`RecordHead`]: crate::ng::psp::RecordHead

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::psp::header::Manifest;
use crate::ng::psp::record::{RecordEncodeError, RecordEncoder};
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
pub struct BlockHead {
    /// Which contig every record in the block sits on.
    pub contig: ContigId,
    /// Where the block's first record starts — and the base its position offset is measured
    /// from, so that record's offset is zero.
    pub first_position: Position,
    /// How many records the block holds. **Never zero**: a block exists because a record went
    /// into it, and a reader refuses one claiming otherwise.
    pub record_count: u64,
}

impl BlockHead {
    /// Append the block's opening three fields to `out`.
    ///
    /// It cannot fail: all three are unbounded variable-length integers, and a block with no
    /// records is unreachable in [`BlockBuilder`] rather than refused here.
    pub fn encode(&self, out: &mut Vec<u8>) {
        encode_u64_leb128(u64::from(self.contig.get()), out);
        encode_u64_leb128(self.first_position.get(), out);
        encode_u64_leb128(self.record_count, out);
    }

    /// Read a block's opening three fields, and say how many bytes they took.
    ///
    /// **A buffer that stops inside them is [`BlockHeadError::Truncated`], not damage**, and
    /// the two are different instructions to a streaming reader: the first means *decompress
    /// more of this block and try again*, the second means *the file is damaged*. At the very
    /// start of a block a short buffer is the ordinary state of affairs rather than an
    /// exceptional one, which is why the split matters here as much as it does in a record.
    pub fn decode(bytes: &[u8]) -> Result<(Self, usize), BlockHeadError> {
        let mut at = 0usize;

        let contig_at = at;
        let contig = read_field(bytes, &mut at, BLOCK_CONTIG_ID)?;
        let contig = u32::try_from(contig).map_err(|_| BlockHeadError::Malformed {
            field: BLOCK_CONTIG_ID,
            bytes_in: contig_at,
            reason: format!("{contig} names no contig; a contig id is a 32-bit index"),
        })?;

        let first_position = read_field(bytes, &mut at, BLOCK_FIRST_POSITION)?;

        let count_at = at;
        let record_count = read_field(bytes, &mut at, BLOCK_RECORD_COUNT)?;
        if record_count == 0 {
            return Err(BlockHeadError::Malformed {
                field: BLOCK_RECORD_COUNT,
                bytes_in: count_at,
                reason: "a block holding no records; a block exists because a record went into it"
                    .to_string(),
            });
        }

        Ok((
            Self {
                contig: ContigId(contig),
                first_position: Position(first_position),
                record_count,
            },
            at,
        ))
    }
}

/// One variable-length integer of a block head, advancing `at` past it.
///
/// **The truncated case is separated from every other varint fault here**, because that is the
/// only place the two classes are told apart and a reader's whole retry loop hangs on it.
fn read_field(bytes: &[u8], at: &mut usize, field: &'static str) -> Result<u64, BlockHeadError> {
    match decode_u64_leb128(&bytes[*at..]) {
        Ok((value, used)) => {
            *at += used;
            Ok(value)
        }
        Err(VarintError::Truncated) => Err(BlockHeadError::Truncated {
            field,
            bytes_in: *at,
        }),
        Err(fault) => Err(BlockHeadError::Malformed {
            field,
            bytes_in: *at,
            reason: fault.to_string(),
        }),
    }
}

/// Why a block's opening fields could not be read.
///
/// **Two classes, and a streaming reader branches on them** — the same split
/// [`RecordDecodeError`] makes for a record, and for the same reason: one says *fetch more
/// bytes and retry*, the other says *this file is damaged*, and a fault put in the wrong class
/// makes a reader either reject a good block or retry for ever on a bad one. There is no third
/// class here because a block head names nothing a later writer could add to.
///
/// [`RecordDecodeError`]: crate::ng::psp::RecordDecodeError
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BlockHeadError {
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
/// **For a caller holding a block entire** — a writer's own tests, the parity oracle, a tool
/// dumping one block. The streaming reader of Milestone D3 never has a whole block in memory
/// and does not go through this.
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub struct BlockRecords<'a> {
    /// What the block's own head said.
    pub head: BlockHead,
    /// The records, still encoded, from the first record's first byte.
    pub records: &'a [u8],
}

impl<'a> BlockRecords<'a> {
    /// Split a whole block payload into its head and the records behind it.
    pub fn split(payload: &'a [u8]) -> Result<Self, BlockHeadError> {
        let (head, used) = BlockHead::decode(payload)?;
        Ok(Self {
            head,
            records: &payload[used..],
        })
    }
}

// ---------------------------------------------------------------------
// The cut
// ---------------------------------------------------------------------

/// Which block a coordinate belongs to: its contig, and which cell of the coordinate grid its
/// start falls in.
///
/// **The grid is the point of the cut rule.** A block ends when a position crosses into the
/// next multiple of the genomic block size, which is not the same thing as a block ending once
/// it has covered that many bases: a grid makes every sample cut at the *same* coordinates, so
/// a cohort reader stepping across a region touches one aligned block per sample rather than
/// one in some samples and two in others (spec §4.1). A running count would not align, and
/// losing that would have been an accident rather than a decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BlockCell {
    contig: ContigId,
    cell: u64,
}

/// The block being built: what its head will say, once it is known how many records went in.
#[derive(Debug, Clone, Copy)]
struct OpenBlock {
    cell: BlockCell,
    first_position: Position,
    record_count: u64,
}

impl OpenBlock {
    fn head(self) -> BlockHead {
        BlockHead {
            contig: self.cell.contig,
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
/// would have produced had the refused record never been offered.
#[derive(Debug)]
pub struct BlockBuilder {
    genomic_block_size_bp: Bp,
    block_byte_ceiling: Option<u32>,
    encoder: RecordEncoder,
    /// The records of the block being built, each already behind its own head.
    records: Vec<u8>,
    /// Where the next record is laid down while a cut is still undecided. Swapped with
    /// `records` when the cut commits, so both buffers stay warm and neither is reallocated
    /// per block.
    next_records: Vec<u8>,
    /// The block just closed: its head, then its records. Handed out by reference and
    /// overwritten by the next close.
    closed: Vec<u8>,
    open: Option<OpenBlock>,
    /// Where the last accepted record sits, so the order check has something to compare
    /// against across a cut — where the encoder's own base has been reset and cannot see it.
    last_accepted: Option<GenomeRegion>,
}

impl BlockBuilder {
    /// A builder cutting on `genomic_block_size_bp`, optionally closing a block early once it
    /// holds `block_byte_ceiling` bytes of records.
    ///
    /// **A zero grid is refused rather than divided by**: the cut is which multiple of this a
    /// coordinate falls in, and a grid with no cells has no answer. A header carrying one is
    /// refused when it is validated (`header.rs`); this is the same rule where the arithmetic
    /// happens, because `Manifest`'s fields are public and a caller can build one that never
    /// met that check.
    pub fn new(
        genomic_block_size_bp: Bp,
        block_byte_ceiling: Option<u32>,
    ) -> Result<Self, BlockWriteError> {
        if genomic_block_size_bp.get() == 0 {
            return Err(BlockWriteError::ZeroGenomicBlockSize);
        }
        Ok(Self {
            genomic_block_size_bp,
            block_byte_ceiling,
            // Replaced before the first record is written: a builder that has been handed none
            // has no block, and so no base to measure one from.
            encoder: RecordEncoder::for_block(Position(0)),
            records: Vec::new(),
            next_records: Vec::new(),
            closed: Vec::new(),
            open: None,
            last_accepted: None,
        })
    }

    /// A builder driven by a file's declared cut rule, which is the only rule a writer
    /// extending that file may use — the manifest is fixed when the file is created and an
    /// append does not rewrite it (spec §6.4).
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, BlockWriteError> {
        Self::new(manifest.genomic_block_size_bp, manifest.block_byte_ceiling)
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
        let cell = self.cell_of(record.region);

        let Some(open) = self.open else {
            // The file's first record. The encoder is positioned before anything is written,
            // so a record the codec refuses leaves no block open behind it.
            self.encoder.start_block(record.region.start);
            self.encoder.encode_record(record, &mut self.records)?;
            self.open_at(cell, record.region.start);
            self.accept(record.region);
            return Ok(None);
        };

        if cell == open.cell && !self.byte_ceiling_reached() {
            self.encoder.encode_record(record, &mut self.records)?;
            self.accept(record.region);
            return Ok(None);
        }

        // A cut. The record goes down first, in the buffer that will become the new block's,
        // so that a record the codec refuses leaves the open block open and loses nothing.
        self.next_records.clear();
        let resume_at = self.encoder.measured_from().position();
        self.encoder.start_block(record.region.start);
        if let Err(refused) = self.encoder.encode_record(record, &mut self.next_records) {
            self.encoder.start_block(resume_at);
            return Err(refused.into());
        }

        self.closed.clear();
        open.head().encode(&mut self.closed);
        self.closed.extend_from_slice(&self.records);

        std::mem::swap(&mut self.records, &mut self.next_records);
        self.open_at(cell, record.region.start);
        self.accept(record.region);
        Ok(Some(&self.closed))
    }

    /// Close the block being built, if there is one. `None` when every record pushed has
    /// already been handed back inside a block, which is also what a builder that was never
    /// pushed to returns.
    ///
    /// **It consumes the builder**, which is the one thing that cannot be got wrong afterwards:
    /// a builder that could be pushed to after closing would put the closed block's records
    /// into the next one, and a builder that could be closed twice would put the last block in
    /// the file twice. Measured on the version where it took `&mut self`: dropping the line
    /// that emptied the record buffer left all twenty-one tests green, because nothing pushed
    /// after closing. The type is what closes that, not a test.
    ///
    /// It returns owned bytes where [`push`](Self::push) lends them — once per file, against a
    /// borrow per block.
    pub fn finish(mut self) -> Option<Vec<u8>> {
        let open = self.open.take()?;
        self.closed.clear();
        open.head().encode(&mut self.closed);
        self.closed.extend_from_slice(&self.records);
        Some(self.closed)
    }

    /// Which block a record's start belongs to. **Its start and not its end**: a record widened
    /// by a deletion may reach past its own block, which costs nothing here — a reader learns
    /// each record's span from its head.
    fn cell_of(&self, region: GenomeRegion) -> BlockCell {
        BlockCell {
            contig: region.contig,
            cell: region.start.get() / self.genomic_block_size_bp.get(),
        }
    }

    /// Whether the open block has already reached the declared byte ceiling.
    ///
    /// **It measures the records laid down, not the block head in front of them**, and it is
    /// checked before the next record rather than after the last — so a block may pass the
    /// ceiling by one record. That is what the rule costs: the alternative decides a record's
    /// fate from a length nothing knows until the record has been encoded.
    fn byte_ceiling_reached(&self) -> bool {
        self.block_byte_ceiling
            .is_some_and(|ceiling| self.records.len() >= ceiling as usize)
    }

    fn open_at(&mut self, cell: BlockCell, first_position: Position) {
        self.open = Some(OpenBlock {
            cell,
            first_position,
            record_count: 0,
        });
    }

    fn accept(&mut self, region: GenomeRegion) {
        if let Some(open) = &mut self.open {
            open.record_count += 1;
        }
        self.last_accepted = Some(region);
    }

    /// Refuse a record that does not follow the one before it along the reference.
    ///
    /// Contigs are visited in ascending order and never revisited; within one contig a record
    /// must not start before the record before it. **Two records starting on the same base are
    /// allowed** — a repeat tract and a generic locus can begin together — because that is not
    /// out of order.
    fn check_order(&self, offered: GenomeRegion) -> Result<(), BlockWriteError> {
        let Some(previous) = self.last_accepted else {
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

/// Why a record could not be laid down in a block.
///
/// **Every variant is a record or a setting the writer was handed that the format cannot
/// hold**, not an internal fault. The writer of Milestone F3 turns each into a
/// [`PspWriteError`] that also names the file, which is the one thing a cohort gathering sixty
/// samples at once needs and this type does not know.
///
/// [`PspWriteError`]: crate::ng::psp::PspWriteError
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BlockWriteError {
    /// A record the record codec cannot lay down: an empty region, a region at the coordinate
    /// ceiling, or a body longer than a head can describe.
    #[error(transparent)]
    Record(#[from] RecordEncodeError),

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
        "a record on contig {} after a record on contig {}; contigs are written in ascending \
         order and never revisited",
        offered.get(), previous.get()
    )]
    ContigOutOfOrder {
        previous: ContigId,
        offered: ContigId,
    },

    /// The cut rule the builder was handed has no cells to cut on.
    #[error(
        "the genomic block size is zero; the cut is a grid on the coordinate and a zero grid \
         has no cells"
    )]
    ZeroGenomicBlockSize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::psp::header::{DEFAULT_GENOMIC_BLOCK_SIZE_BP, DEFAULT_LOOK_BACK_WINDOW_LOG};
    use crate::ng::psp::record::{OffsetBase, RecordLayout, decode_record, record_fields};
    use crate::ng::types::{ReadGroupId, SummedLogError};

    /// The 100 kb grid the format ships with.
    const A_GRID: Bp = DEFAULT_GENOMIC_BLOCK_SIZE_BP;

    /// One ordinary record: a covered position whose reads all agreed with the reference.
    ///
    /// **The reference bases differ between records of the same length**, so a walk that took
    /// one record's payload from another is visible rather than hidden by identical fixtures —
    /// which a review measured happening on the record codec's own run.
    fn a_record(contig: u32, start: u64, span: u64) -> SampleLocusObservations {
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

    /// Walk one block payload back: its head, then every record in it, each position rebuilt
    /// from the block's own first position and the offsets since.
    ///
    /// **This is what the cut is checked against, and it takes nothing from the builder**: it
    /// starts from the block's declared first position, which is all a reader beginning at that
    /// block would have.
    fn walk(payload: &[u8]) -> (BlockHead, Vec<SampleLocusObservations>) {
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
            found.head.record_count,
            "the block held the number of records its head declared"
        );
        (found.head, records)
    }

    /// Every block's records, in the order the blocks were cut.
    fn walk_all(blocks: &[Vec<u8>]) -> Vec<SampleLocusObservations> {
        blocks.iter().flat_map(|block| walk(block).1).collect()
    }

    /// How many bytes of records a block payload holds, past its own head.
    fn record_bytes_in(payload: &[u8]) -> usize {
        BlockRecords::split(payload)
            .expect("the head reads")
            .records
            .len()
    }

    // -----------------------------------------------------------------
    // The block head
    // -----------------------------------------------------------------

    #[test]
    fn a_block_head_round_trips_and_says_how_many_bytes_it_took() {
        let head = BlockHead {
            contig: ContigId(7),
            first_position: Position(90_600_000),
            record_count: 24_881,
        };
        let mut bytes = Vec::new();
        head.encode(&mut bytes);
        let (read, used) = BlockHead::decode(&bytes).expect("it reads back");
        assert_eq!(read, head);
        assert_eq!(used, bytes.len());
    }

    /// **The bytes a file carries.** A block head is container framing, not a manifest field,
    /// so nothing in a file says how it is laid out — which makes this array the format itself.
    /// A reordering or a change of encoding has to fail here and be a version bump.
    #[test]
    fn a_block_head_is_these_exact_bytes() {
        let mut bytes = Vec::new();
        BlockHead {
            contig: ContigId(1),
            first_position: Position(300),
            record_count: 2,
        }
        .encode(&mut bytes);
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
    #[test]
    fn a_block_head_cut_short_is_truncated_at_every_cut_and_never_malformed() {
        let mut whole = Vec::new();
        BlockHead {
            contig: ContigId(300),
            first_position: Position(90_600_000),
            record_count: 24_881,
        }
        .encode(&mut whole);
        assert!(
            whole.len() > 6,
            "every field of the fixture must be multi-byte, or the cuts miss the interesting \
             boundaries; it is {} bytes",
            whole.len()
        );

        for cut_at in 0..whole.len() {
            match BlockHead::decode(&whole[..cut_at]) {
                Err(BlockHeadError::Truncated { .. }) => {}
                other => panic!("a head of {cut_at} bytes gave {other:?}"),
            }
        }
        assert!(
            BlockHead::decode(&whole).is_ok(),
            "and the whole head reads"
        );
    }

    /// A block claiming no records is damage. Nothing this builder writes can say it — a block
    /// exists because a record went into it — so a file that does is not one this writer
    /// produced.
    #[test]
    fn a_block_claiming_no_records_is_refused() {
        let mut bytes = Vec::new();
        BlockHead {
            contig: ContigId(1),
            first_position: Position(300),
            record_count: 1,
        }
        .encode(&mut bytes);
        let last = bytes.len() - 1;
        bytes[last] = 0;

        match BlockHead::decode(&bytes) {
            Err(BlockHeadError::Malformed { field, reason, .. }) => {
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
            Err(BlockHeadError::Malformed { field, .. }) => assert_eq!(field, BLOCK_CONTIG_ID),
            other => panic!("expected Malformed, got {other:?}"),
        }
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
        assert_eq!(walk(&blocks[0]).0.first_position, Position(99_999));
        assert_eq!(walk(&blocks[1]).0.first_position, Position(100_000));
    }

    /// **Every sample cuts at the same coordinates**, which is what the grid buys and what a
    /// running total would lose. Two samples covering the same region with entirely different
    /// records give blocks whose first positions fall in the same grid cells.
    #[test]
    fn two_samples_with_different_records_cut_at_the_same_grid_cells() {
        let dense: Vec<_> = (0..170)
            .map(|index| a_record(0, 90_000 + index * 1_000, 1))
            .collect();
        let sparse: Vec<_> = [90_500u64, 150_000, 250_000, 250_001]
            .into_iter()
            .map(|start| a_record(0, start, 1))
            .collect();

        let cells_of = |blocks: &[Vec<u8>]| -> Vec<u64> {
            blocks
                .iter()
                .map(|block| walk(block).0.first_position.get() / A_GRID.get())
                .collect()
        };

        let a_builder = || BlockBuilder::new(A_GRID, None).expect("a grid");
        let dense_cells = cells_of(&cut(a_builder(), &dense).expect("in order"));
        let sparse_cells = cells_of(&cut(a_builder(), &sparse).expect("in order"));

        assert_eq!(dense_cells, vec![0, 1, 2]);
        assert_eq!(sparse_cells, vec![0, 1, 2]);
        assert_eq!(
            dense_cells, sparse_cells,
            "which cells a sample's blocks start in is a function of the coordinate alone"
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
        assert_eq!(walk(&blocks[0]).0.contig, ContigId(0));
        assert_eq!(walk(&blocks[1]).0.contig, ContigId(1));
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

    /// **A block may pass the ceiling by one record**, because the check runs before the next
    /// record rather than after the last: the alternative decides a record's fate from a length
    /// nothing knows until the record has been encoded. Pinned so a change is a decision rather
    /// than a surprise.
    #[test]
    fn a_block_reaches_the_ceiling_before_it_is_closed_so_it_may_pass_it_by_one_record() {
        let records: Vec<_> = (0..6).map(|index| a_record(0, 100 + index, 1)).collect();

        let blocks = cut(
            BlockBuilder::new(A_GRID, Some(1)).expect("a grid"),
            &records,
        )
        .expect("in order");

        assert_eq!(
            blocks.len(),
            records.len(),
            "a one-byte ceiling gives each record a block of its own, never a block of none"
        );
        for block in &blocks {
            assert_eq!(walk(block).1.len(), 1);
            assert!(
                record_bytes_in(block) > 1,
                "and each of those blocks is past the ceiling it was cut on"
            );
        }
    }

    /// **The oracle for the whole cut: every record comes back, once, in order, wherever the
    /// blocks fell.** Records over three contigs, three grid cells each, and a byte ceiling
    /// that fires inside them.
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
            }
        }

        let by_grid =
            cut(BlockBuilder::new(A_GRID, None).expect("a grid"), &records).expect("in order");
        assert_eq!(by_grid.len(), 9, "three contigs of three grid cells each");
        assert_eq!(walk_all(&by_grid), records);

        // A ceiling measured from what the blocks actually hold, so it must fire whatever a
        // record's size turns out to be. Guessing one is how a test that proves nothing about
        // the ceiling passes: an earlier draft of this guessed 200 bytes and never fired.
        let smallest = by_grid
            .iter()
            .map(|block| record_bytes_in(block))
            .min()
            .expect("nine blocks");
        let ceiling = u32::try_from(smallest / 3).expect("a third of a block's records");
        assert!(ceiling >= 1, "a ceiling of zero is refused by the header");

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
            let (head, in_block) = walk(block);
            assert_eq!(head.record_count, in_block.len() as u64);
            assert_eq!(head.first_position, in_block[0].region.start);
            for record in &in_block {
                assert_eq!(record.region.contig, head.contig);
            }
            seen += in_block.len();
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
        assert!(
            matches!(refused, BlockWriteError::OutOfOrder { .. }),
            "got {refused:?}"
        );
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
    #[test]
    fn a_contig_already_left_or_one_that_goes_backwards_is_refused() {
        for second in [0u32, 1] {
            let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
            builder.push(&a_record(0, 500, 1)).expect("contig 0");
            builder.push(&a_record(2, 500, 1)).expect("contig 2");

            match builder.push(&a_record(second, 500, 1)) {
                Err(BlockWriteError::ContigOutOfOrder { previous, offered }) => {
                    assert_eq!(previous, ContigId(2));
                    assert_eq!(offered, ContigId(second));
                }
                other => panic!("contig {second} after contig 2 gave {other:?}"),
            }
        }
    }

    /// **A refused record leaves the builder exactly as it was.** A run with three kinds of
    /// refusal offered in the middle of it — one of them where a cut is about to happen — gives
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
            // Offered before every record after the first, so a refusal that leaves a trace
            // shows wherever the cuts fall rather than only at one chosen point. The third
            // refusal is the one the codec makes, which happens *after* the cut decision — it
            // is the rollback that keeps it invisible.
            assert!(matches!(
                interrupted.push(&a_record(1, 1, 1)),
                Err(BlockWriteError::OutOfOrder { .. })
            ));
            assert!(matches!(
                interrupted.push(&a_record(0, 300_000, 1)),
                Err(BlockWriteError::ContigOutOfOrder { .. })
            ));
            assert!(matches!(
                interrupted.push(&a_record_over_no_base(1, 300_000)),
                Err(BlockWriteError::Record(
                    RecordEncodeError::EmptyRegion { .. }
                ))
            ));
            if let Some(closed) = interrupted.push(record).expect("in order") {
                with.push(closed.to_vec());
            }
        }
        if let Some(last) = interrupted.finish() {
            with.push(last.to_vec());
        }

        assert_eq!(with, without, "the refusals left no trace in the bytes");
    }

    /// A grid with no cells is refused rather than divided by.
    #[test]
    fn a_zero_genomic_block_size_is_refused() {
        let refused = BlockBuilder::new(Bp(0), None).expect_err("a zero grid has no cells");
        assert!(
            matches!(refused, BlockWriteError::ZeroGenomicBlockSize),
            "got {refused:?}"
        );
        assert!(refused.to_string().contains("zero grid has no cells"));
    }

    /// A refusal on the very first record leaves no block behind it, so the next record still
    /// opens the file's first block.
    #[test]
    fn a_refusal_on_the_first_record_leaves_no_block_open() {
        let mut builder = BlockBuilder::new(A_GRID, None).expect("a grid");
        assert!(builder.push(&a_record_over_no_base(0, 500)).is_err());
        assert!(
            builder.finish().is_none(),
            "a builder whose only record was refused has no block to close"
        );

        let mut again = BlockBuilder::new(A_GRID, None).expect("a grid");
        assert!(again.push(&a_record_over_no_base(0, 500)).is_err());
        let blocks = cut(again, &[a_record(0, 500, 1)]).expect("in order");
        assert_eq!(walk(&blocks[0]).0.first_position, Position(500));
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

    /// The manifest is where a file's cut rule lives, so a builder made from one cuts the way
    /// the numbers in it say.
    #[test]
    fn a_builder_from_a_manifest_cuts_on_the_manifests_grid() {
        let manifest = Manifest {
            genomic_block_size_bp: Bp(1_000),
            block_byte_ceiling: None,
            look_back_window_log: DEFAULT_LOOK_BACK_WINDOW_LOG,
            fields: record_fields(),
        };
        let blocks = cut(
            BlockBuilder::from_manifest(&manifest).expect("a grid"),
            &[a_record(0, 999, 1), a_record(0, 1_000, 1)],
        )
        .expect("in order");
        assert_eq!(blocks.len(), 2, "a 1,000 bp grid cuts at 1,000");
    }
}
