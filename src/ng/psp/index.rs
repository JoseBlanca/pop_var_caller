//! The block index: one entry per psp block, in genomic order, so a reader can start at
//! any block without reading what comes before it.
//!
//! **A flat vector, decoded whole at open** — production's shape (`src/psp/index.rs`), and
//! the large psp block is what makes it affordable to keep. Measured on a tomato
//! accession: **154 entries at 1 MB blocks against 1,674 in the equivalent `.psp`**, which
//! scales to roughly 14,000 entries and a few hundred kilobytes for a whole genome
//! (spec §3.3). At production's 5 kb blocks the same genome needs about 156,000 entries and
//! a 3.8 MB index *per open sample*, which is what the memory budget cannot afford — so a
//! large block removes the problem rather than solving it, and the coarse-index-and-chain
//! scheme `run_streaming.md` §7.2 asks for must not be built.

use crate::ng::types::{ContigId, GenomePosition, Position};
use crate::psp::varint::{decode_u64_leb128, encode_u64_leb128};

/// One psp block, as the index names it.
///
/// **And nothing more.** An earlier draft carried the largest non-reference support in the
/// block so a cohort scan could skip whole blocks without decompressing them; at a 100 kb
/// block essentially every block contains something that varies, so the field never fired
/// and it is gone (spec psp_record_encoding.md §2.4, spec §6.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BlockIndexEntry {
    /// The contig and the position of the block's first record.
    ///
    /// **One field, not a contig beside a position**, because [`GenomePosition`] is that
    /// pair and its ordering is contig-then-position — which is exactly the key a lookup
    /// turning a coordinate into a block searches on. A block never crosses a contig
    /// (spec §3.2), so one entry names one contig.
    pub first_position: GenomePosition,
    /// Byte offset of the block's first byte from the start of the file.
    pub offset: u64,
}

/// Bytes one entry can occupy at its smallest: two one-byte varints and the fixed offset.
///
/// **It exists to bound an allocation, not to describe a file.** The block count comes from the
/// footer, which is as damaged as the rest of a damaged file, so the reservation is capped
/// against what the buffer could actually hold.
///
/// Measured, by removing the cap: a footer claiming `u64::MAX` blocks reserves that many
/// [`BlockIndexEntry`], **24 bytes each**, and the decode dies inside `Vec::with_capacity` with
/// `capacity overflow`. **That is a panic on a damaged file**, which this module does not do —
/// a corrupt psp is data a run was handed, not a bug (see `PspReadError`). So the cap is what
/// keeps a hostile tail an error rather than a crash, and the memory is the lesser half of it.
const SMALLEST_ENTRY_BYTES: u64 = 1 + 1 + 8;

/// Bytes one entry can occupy at its largest: two five-byte varints and the fixed offset.
/// Used to size the writer's buffer once, so the worst case still allocates once.
const LARGEST_ENTRY_BYTES: usize = 5 + 5 + 8;

/// The index's own bytes: for each block in order, its contig, its first position, and the
/// byte offset it starts at.
///
/// **Contig and position are variable-length; the offset is a fixed eight bytes.** That is
/// production's split ([`src/psp/index.rs`](../../../../src/psp/index.rs)) and its reasoning
/// holds here: contig numbers and positions are small in most files and compress to one or two
/// bytes each, while an offset is large and arbitrary, so a varint would usually cost the same
/// eight bytes and sometimes nine.
///
/// **The index is not compressed.** It is read whole at open, before any block is touched, and
/// spec §3.3 puts a whole genome at roughly 14,000 entries — a few hundred kilobytes, paid once
/// per open sample.
pub fn encode_index(entries: &[BlockIndexEntry]) -> Vec<u8> {
    let mut out = Vec::with_capacity(entries.len() * LARGEST_ENTRY_BYTES);
    for entry in entries {
        encode_u64_leb128(u64::from(entry.first_position.contig.0), &mut out);
        encode_u64_leb128(entry.first_position.position.get(), &mut out);
        out.extend_from_slice(&entry.offset.to_le_bytes());
    }
    out
}

/// Read the index back, and refuse one a seek could not be driven by.
///
/// `expected_blocks` is the footer's block count: exactly that many entries are decoded and the
/// buffer must be exhausted afterwards, so a footer and an index that disagree is a refusal
/// rather than a short read.
///
/// **The order is checked here rather than left to a caller**, and that is the difference from
/// production, which says coordinate monotonicity is "the reader's responsibility (one layer
/// up)". Every use of this index is a search for the block holding a coordinate; on entries
/// that are not ordered, that search does not fail — **it silently returns the wrong block**,
/// and the records that come back are a real sample's records from the wrong place. The check
/// costs one comparison per entry, paid once per open file.
///
/// **Two blocks may share a first position, and the rule is written for it.** A block closed
/// early by the byte ceiling is followed by one that starts at the next record, and two records
/// may begin on the same base — a repeat tract and a generic locus do. So positions must be
/// **non-decreasing** while offsets must **strictly increase**: it is the offsets that say the
/// blocks are distinct, and two entries at one offset would be one block indexed twice.
pub fn decode_index(
    bytes: &[u8],
    expected_blocks: u64,
) -> Result<Vec<BlockIndexEntry>, IndexDecodeError> {
    let bounded = (bytes.len() as u64 / SMALLEST_ENTRY_BYTES).saturating_add(1);
    let mut entries = Vec::with_capacity(expected_blocks.min(bounded) as usize);

    let mut at = 0usize;
    for entry in 0..expected_blocks {
        let entry = entry as usize;
        let contig = take_varint(bytes, &mut at, entry, "contig")?;
        let contig = u32::try_from(contig).map_err(|_| IndexDecodeError::FieldTooLarge {
            entry,
            field: "contig",
            found: contig,
        })?;
        let position = take_varint(bytes, &mut at, entry, "first-position")?;
        let offset = take_offset(bytes, &mut at, entry)?;
        entries.push(BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(contig),
                position: Position(position),
            },
            offset,
        });
    }

    if at != bytes.len() {
        return Err(IndexDecodeError::TrailingBytes {
            trailing: bytes.len() - at,
            entries: entries.len(),
        });
    }

    for (entry, pair) in entries.windows(2).enumerate() {
        let (previous, offered) = (pair[0], pair[1]);
        if offered.first_position < previous.first_position {
            return Err(IndexDecodeError::OutOfOrder {
                entry: entry + 1,
                previous: previous.first_position,
                offered: offered.first_position,
            });
        }
        if offered.offset <= previous.offset {
            return Err(IndexDecodeError::OffsetNotAscending {
                entry: entry + 1,
                previous: previous.offset,
                offered: offered.offset,
            });
        }
    }

    Ok(entries)
}

/// The checksum the footer carries over the index's bytes.
///
/// **Production's function, called rather than copied**
/// ([`src/psp/index.rs`](../../../../src/psp/index.rs)): it is XXH3-64 truncated to its low 32
/// bits, which is the truncation zstd uses for its own frame checksum, and the reason to share
/// it is that there is then one XXH3 in the codebase rather than two that could differ. The
/// index is the one region of a psp no zstd frame checksum covers, which is why it has this.
pub fn checksum_index(bytes: &[u8]) -> u32 {
    crate::psp::index::checksum_index(bytes)
}

fn take_varint(
    bytes: &[u8],
    at: &mut usize,
    entry: usize,
    field: &'static str,
) -> Result<u64, IndexDecodeError> {
    let rest = bytes
        .get(*at..)
        .ok_or(IndexDecodeError::Truncated { entry, field })?;
    let (value, used) =
        decode_u64_leb128(rest).map_err(|_| IndexDecodeError::Truncated { entry, field })?;
    *at += used;
    Ok(value)
}

fn take_offset(bytes: &[u8], at: &mut usize, entry: usize) -> Result<u64, IndexDecodeError> {
    let field = "offset";
    let end = at
        .checked_add(8)
        .ok_or(IndexDecodeError::Truncated { entry, field })?;
    let slice = bytes
        .get(*at..end)
        .ok_or(IndexDecodeError::Truncated { entry, field })?;
    let offset = u64::from_le_bytes(
        slice
            .try_into()
            .expect("an eight-byte slice converts to an eight-byte array"),
    );
    *at = end;
    Ok(offset)
}

/// Why an index could not be read.
///
/// **Its own type, carrying no path**, like this module's other codec errors: the file's name
/// belongs to whoever opened the file, and F4's `open` is what dresses one of these as a
/// [`PspReadError`](super::PspReadError). Every variant means *this file is damaged, rebuild
/// it* — there is no variant here that a caller could act on by changing a setting.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum IndexDecodeError {
    /// The index ended in the middle of an entry the footer said was there.
    #[error("the index ends inside entry {entry}, before its {field}")]
    Truncated { entry: usize, field: &'static str },

    /// A contig number too large for a [`ContigId`]. **Refused rather than truncated**: a
    /// wrapped contig number indexes a real contig, and the seek would land on it.
    #[error("entry {entry} declares {field} {found}, which is too large for a contig number")]
    FieldTooLarge {
        entry: usize,
        field: &'static str,
        found: u64,
    },

    /// The footer's block count and the index's bytes disagree — bytes are left over once that
    /// many entries have been read.
    #[error("the index holds {trailing} bytes after the {entries} entries the footer declares")]
    TrailingBytes { trailing: usize, entries: usize },

    /// An entry starts before the one before it. **A search over these does not fail on
    /// disorder, it returns the wrong block**, which is why this is refused at open.
    #[error(
        "index entry {entry} starts at {offered:?}, before entry {} at {previous:?}",
        entry - 1
    )]
    OutOfOrder {
        entry: usize,
        previous: GenomePosition,
        offered: GenomePosition,
    },

    /// Two entries at the same byte offset, or offsets that go backwards. Positions may repeat
    /// — a block closed by the byte ceiling is followed by one starting at the same base — but
    /// the offsets are what say two entries are two blocks.
    #[error(
        "index entry {entry} is at byte {offered}, which does not follow entry {} at byte {previous}",
        entry - 1
    )]
    OffsetNotAscending {
        entry: usize,
        previous: u64,
        offered: u64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    /// Turning a coordinate into a block is a search over the index, and it has to order
    /// by contig before position — otherwise a position on a later contig would be found
    /// on an earlier one whenever the numbers happen to fall that way.
    #[test]
    fn entries_sort_by_contig_before_position() {
        let entry = |contig: u32, position: u64, offset: u64| BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(contig),
                position: Position(position),
            },
            offset,
        };
        let late_on_the_first_contig = entry(0, 90_000_000, 4_096);
        let early_on_the_second = entry(1, 1, 8_192);
        assert!(late_on_the_first_contig < early_on_the_second);

        let mut blocks = [early_on_the_second, late_on_the_first_contig];
        blocks.sort();
        assert_eq!(blocks[0], late_on_the_first_contig);
    }

    // -----------------------------------------------------------------
    // F1 — the index's own bytes
    // -----------------------------------------------------------------

    fn at(contig: u32, position: u64, offset: u64) -> BlockIndexEntry {
        BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(contig),
                position: Position(position),
            },
            offset,
        }
    }

    /// A tomato accession's whole index, as spec §3.3 measures it: 154 blocks. Round-tripped
    /// entry for entry, because a search over these is only as good as the coordinates in them.
    #[test]
    fn an_index_round_trips_entry_for_entry() {
        let entries: Vec<BlockIndexEntry> = (0..154)
            .map(|block| {
                at(
                    block / 20,
                    1 + u64::from(block) * 100_000,
                    4_096 + u64::from(block) * 65_536,
                )
            })
            .collect();
        let bytes = encode_index(&entries);
        let read_back = decode_index(&bytes, entries.len() as u64).expect("its own bytes decode");
        assert_eq!(read_back, entries);
    }

    /// An index with no blocks in it is legal and is not an error: a psp holding no records has
    /// one. **Its bytes are empty**, and a decode asking for none must take none rather than
    /// reading the first entry of whatever follows.
    #[test]
    fn an_empty_index_round_trips_and_takes_no_bytes() {
        let bytes = encode_index(&[]);
        assert!(bytes.is_empty());
        assert_eq!(decode_index(&bytes, 0), Ok(Vec::new()));
    }

    /// The extremes of every field, in one index: contig zero and `u32::MAX`, position 1 and
    /// `u64::MAX`, offset 0 and `u64::MAX`. **The widest values are where a varint's width and a
    /// fixed field's endianness disagree if either is wrong**, and they are cheap to state.
    #[test]
    fn the_widest_value_of_every_field_round_trips() {
        let entries = vec![
            at(0, 1, 0),
            at(1, u64::MAX / 2, 1),
            at(u32::MAX, u64::MAX, u64::MAX),
        ];
        let bytes = encode_index(&entries);
        assert_eq!(
            decode_index(&bytes, 3).expect("the extremes decode"),
            entries
        );
    }

    /// **Two blocks may share a first position**, and the index must accept it.
    ///
    /// Measured on `BlockBuilder`, not assumed: a byte ceiling closes a block after its first
    /// record, and two records may begin on the same base — a repeat tract and a generic locus
    /// do — so the block that opens next starts where the one just closed did. A rule demanding
    /// strictly increasing positions would refuse a file this writer produces.
    #[test]
    fn two_blocks_starting_on_the_same_base_are_accepted() {
        let entries = vec![at(0, 500, 4_096), at(0, 500, 8_192)];
        let bytes = encode_index(&entries);
        assert_eq!(decode_index(&bytes, 2).expect("a legal file"), entries);
    }

    /// A position that goes backwards is refused, and the message names the entry and both
    /// coordinates. **A search over a disordered index does not fail — it returns the wrong
    /// block**, and the records that come back are real records from the wrong place.
    #[test]
    fn an_index_whose_positions_go_backwards_is_refused() {
        for (previous, offered) in [
            (at(0, 900, 4_096), at(0, 800, 8_192)),
            (at(1, 10, 4_096), at(0, 10, 8_192)),
        ] {
            let bytes = encode_index(&[previous, offered]);
            let refused = decode_index(&bytes, 2).expect_err("backwards must be refused");
            assert_eq!(
                refused,
                IndexDecodeError::OutOfOrder {
                    entry: 1,
                    previous: previous.first_position,
                    offered: offered.first_position,
                }
            );
        }
    }

    /// Offsets must strictly increase: two entries at one offset are one block indexed twice,
    /// and an offset that goes backwards points a seek behind where it has already read.
    ///
    /// **This is the half positions cannot carry**, because positions are allowed to repeat.
    #[test]
    fn an_index_whose_offsets_repeat_or_go_backwards_is_refused() {
        for (previous, offered) in [
            (at(0, 500, 4_096), at(0, 600, 4_096)),
            (at(0, 500, 8_192), at(0, 600, 4_096)),
        ] {
            let bytes = encode_index(&[previous, offered]);
            let refused = decode_index(&bytes, 2).expect_err("offsets must ascend");
            assert_eq!(
                refused,
                IndexDecodeError::OffsetNotAscending {
                    entry: 1,
                    previous: previous.offset,
                    offered: offered.offset,
                }
            );
        }
    }

    /// The footer's block count and the index's bytes must agree in both directions.
    ///
    /// Too few entries for the bytes leaves a tail: **that is the case a lenient decoder would
    /// read as a shorter file**, and every block past the count would be unreachable with
    /// nothing said. Too many runs off the end and is a truncation.
    #[test]
    fn a_block_count_that_disagrees_with_the_bytes_is_refused_both_ways() {
        let entries = vec![
            at(0, 1, 4_096),
            at(0, 100_001, 8_192),
            at(0, 200_001, 16_384),
        ];
        let entries_left_over = vec![entries[2]];
        let bytes = encode_index(&entries);

        let refused = decode_index(&bytes, 2).expect_err("a tail must be refused");
        match refused {
            IndexDecodeError::TrailingBytes { trailing, entries } => {
                assert_eq!(entries, 2);
                // The third entry's own bytes, computed rather than recalled: a contig varint,
                // a three-byte position varint for 200,001, and the eight-byte offset.
                assert_eq!(trailing, encode_index(&entries_left_over).len());
                assert_eq!(trailing, 12);
            }
            other => panic!("expected trailing bytes, got {other}"),
        }

        let refused = decode_index(&bytes, 4).expect_err("a short read must be refused");
        assert_eq!(
            refused,
            IndexDecodeError::Truncated {
                entry: 3,
                field: "contig",
            }
        );
    }

    /// Every byte prefix of a real index is refused, and **none of them panics**. A file cut
    /// short anywhere inside its index is damage, not a shorter index.
    #[test]
    fn every_truncation_of_an_index_is_refused_without_panicking() {
        let entries = vec![
            at(0, 1, 4_096),
            at(3, 100_001, 8_192),
            at(3, 200_001, 16_384),
        ];
        let whole = encode_index(&entries);
        // Ten bytes, then twelve, then twelve: a contig varint, a position varint that is one
        // byte at position 1 and three at 100,001 and 200,001, and the fixed eight-byte offset.
        assert_eq!(whole.len(), 34);

        let mut refused_as_truncated = 0;
        for cut in 0..whole.len() {
            match decode_index(&whole[..cut], entries.len() as u64) {
                Err(IndexDecodeError::Truncated { .. }) => refused_as_truncated += 1,
                Err(other) => panic!("a cut at {cut} gave {other}"),
                Ok(read) => panic!("a cut at {cut} decoded {} entries", read.len()),
            }
        }
        assert_eq!(
            refused_as_truncated,
            whole.len(),
            "every cut short of the whole is a truncation"
        );
    }

    /// A contig number too large for a [`ContigId`] is refused rather than wrapped. **A wrapped
    /// number names a real contig**, so the seek would succeed and land on the wrong chromosome.
    #[test]
    fn a_contig_number_too_large_for_the_type_is_refused_rather_than_wrapped() {
        let mut bytes = Vec::new();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut bytes);
        encode_u64_leb128(1, &mut bytes);
        bytes.extend_from_slice(&4_096u64.to_le_bytes());
        assert_eq!(
            decode_index(&bytes, 1).expect_err("it does not fit a contig number"),
            IndexDecodeError::FieldTooLarge {
                entry: 0,
                field: "contig",
                found: u64::from(u32::MAX) + 1,
            }
        );
    }

    /// **A footer claiming more blocks than any file could hold must not size a buffer from that
    /// number.** The count comes from the tail of a file that may be damaged or hostile, so the
    /// reservation is bounded by the bytes actually present.
    ///
    /// Measured rather than argued: with the cap removed, this same input panics inside
    /// `Vec::with_capacity` with `capacity overflow` — `u64::MAX` entries of 24 bytes is more
    /// than a `usize` of bytes can express. **A panic is the wrong answer to a damaged file**,
    /// so this test fails loudly on that mutation rather than merely reporting more memory.
    #[test]
    fn a_block_count_from_a_damaged_footer_does_not_size_the_allocation() {
        let bytes = encode_index(&[at(0, 1, 4_096)]);
        let refused = decode_index(&bytes, u64::MAX).expect_err("the file cannot hold them");
        assert!(matches!(
            refused,
            IndexDecodeError::Truncated { entry: 1, .. }
        ));
    }

    /// The checksum is over the index's bytes, and a single changed byte moves it.
    ///
    /// **It is production's function**, so this test is about the wiring rather than the hash:
    /// the one region of a psp no zstd frame checksum covers is this one.
    #[test]
    fn the_checksum_moves_when_any_byte_of_the_index_does() {
        let entries = vec![at(0, 1, 4_096), at(0, 100_001, 8_192)];
        let bytes = encode_index(&entries);
        let pristine = checksum_index(&bytes);
        assert_eq!(pristine, checksum_index(&encode_index(&entries)));

        let mut changed = 0;
        for byte in 0..bytes.len() {
            let mut damaged = bytes.clone();
            damaged[byte] ^= 0x01;
            if checksum_index(&damaged) != pristine {
                changed += 1;
            }
        }
        assert_eq!(
            changed,
            bytes.len(),
            "flipping a bit in any of the {} bytes must move the checksum",
            bytes.len()
        );
    }

    /// The bytes are the format, so their layout is stated once against a literal rather than
    /// only against this module's own decoder. **A change here is a format change**; from
    /// Milestone F it costs a version.
    #[test]
    fn one_entry_encodes_to_these_exact_bytes() {
        let bytes = encode_index(&[at(3, 300, 4_096)]);
        assert_eq!(
            bytes,
            vec![
                3, // contig 3, one varint byte
                0xAC, 0x02, // position 300, two varint bytes
                0x00, 0x10, 0, 0, 0, 0, 0, 0, // offset 4,096, eight little-endian bytes
            ]
        );
    }
}
