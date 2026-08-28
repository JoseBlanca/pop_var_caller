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
use crate::psp::errors::VarintError;
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
    ///
    /// **`block_offset` and not `offset`**, because unqualified *offset* already means a
    /// base-pair distance inside this module: a record's first wire field is its
    /// `position-offset`, the number of bases since the record before it. Both are `u64` and
    /// they would sit on one struct. Production spells this same value `block_offset`, and the
    /// two neighbours in the footer that mean bytes qualify it too (`index_offset`,
    /// `trailer_offset`).
    pub block_offset: u64,
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

/// Bytes one entry can occupy at its largest, used to size the writer's buffer once.
///
/// **A contig number is a `u32` and a position is a `u64`**, so their variable-length forms run
/// to five bytes and **ten**, not five and five: 23 bytes an entry, not 18.
///
/// ⚠ The 18 was here first, carried over from production, where all four of the index's fields
/// genuinely are `u32` ([`src/psp/index.rs`](../../../../src/psp/index.rs)). ng's position is
/// not, and `the_widest_value_of_every_field_round_trips` already builds the very entry the old
/// comment said could not exist. Nothing was corrupt — a `Vec` grows — but a constant named for
/// a bound has to hold it.
const LARGEST_ENTRY_BYTES: usize = 5 + 10 + 8;

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
        // Destructured with no `..`: **a field added to the entry is a compile error here**
        // rather than a field silently left out of every index this build writes. Reading the
        // fields one at a time compiled clean through exactly that mutation, and the round trip
        // stayed green because the fixtures zero-filled the new field.
        let BlockIndexEntry {
            first_position,
            block_offset,
        } = *entry;
        encode_u64_leb128(u64::from(first_position.contig.0), &mut out);
        encode_u64_leb128(first_position.position.get(), &mut out);
        out.extend_from_slice(&block_offset.to_le_bytes());
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
/// **What it does not check, and where those live.** A contig number is checked against the
/// width of a [`ContigId`] and **not** against the header's contig list, so an index may name a
/// contig the file does not have; and a byte offset is not checked against the file's length,
/// because nothing here knows it. Both belong to `open`, which has the header and the file size
/// — this function's guarantee is that the entries are well-formed and ordered, not that they
/// point anywhere real.
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

    let mut byte_cursor = 0usize;
    for entry_number in 0..expected_blocks {
        let entry_number = entry_number as usize;
        let contig = take_varint(
            bytes,
            &mut byte_cursor,
            entry_number,
            IndexEntryField::Contig,
        )?;
        let contig = u32::try_from(contig).map_err(|_| IndexDecodeError::ContigNumberTooLarge {
            entry_number,
            found: contig,
        })?;
        let position = take_varint(
            bytes,
            &mut byte_cursor,
            entry_number,
            IndexEntryField::FirstPosition,
        )?;
        let block_offset = take_block_offset(bytes, &mut byte_cursor, entry_number)?;
        entries.push(BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(contig),
                position: Position(position),
            },
            block_offset,
        });
    }

    if byte_cursor != bytes.len() {
        return Err(IndexDecodeError::TrailingBytes {
            trailing_bytes: bytes.len() - byte_cursor,
            entries_read: entries.len(),
        });
    }

    // **Every consecutive pair, and the numbers are carried rather than computed.** An earlier
    // version worked the earlier entry's number out as `entry - 1` inside the message template,
    // which panics with `attempt to subtract with overflow` when rendering an error whose entry
    // is 0 — a panic in `Display`, on the path that reports a damaged file, from the very type
    // whose contract is that a damaged file is never a panic. `decode_index` cannot produce
    // entry 0 here, but the variants are public and constructible, and this file's own tests
    // build them by hand.
    for (earlier, pair) in entries.windows(2).enumerate() {
        let (previous, offered) = (pair[0], pair[1]);
        let (previous_entry, entry_number) = (earlier, earlier + 1);
        if offered.first_position < previous.first_position {
            return Err(IndexDecodeError::OutOfOrder {
                entry_number,
                previous_entry,
                previous: previous.first_position,
                offered: offered.first_position,
            });
        }
        if offered.block_offset <= previous.block_offset {
            return Err(IndexDecodeError::OffsetNotAscending {
                entry_number,
                previous_entry,
                previous: previous.block_offset,
                offered: offered.block_offset,
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
/// it is that there is then one XXH3 in the codebase rather than two that could differ.
///
/// **Why the index carries one when the header and the footer do not.** All four regions
/// outside the blocks are uncompressed and so are outside zstd's own per-frame checksum, but
/// the other three carry their own framing that damage shows up in: the header is length-
/// prefixed with a sentinel that must follow it, and the footer ends in a magic. The index has
/// neither — it is a bare run of entries whose only cross-check is the block count in the
/// footer, and a bit flipped inside an offset produces a perfectly well-formed index pointing
/// somewhere else.
///
/// *(This comment previously said the index is the one region no frame checksum covers. It is
/// not; that is true of production's layout, not of ng's.)*
pub fn checksum_index(bytes: &[u8]) -> u32 {
    crate::psp::index::checksum_index(bytes)
}

fn take_varint(
    bytes: &[u8],
    at: &mut usize,
    entry: usize,
    field: IndexEntryField,
) -> Result<u64, IndexDecodeError> {
    let rest = bytes.get(*at..).ok_or(IndexDecodeError::Truncated {
        entry_number: entry,
        field,
    })?;
    // **The two damages a varint has are two errors, not one.** `Truncated` means the bytes
    // ran out; `Overflow` means there are more continuation bytes than any `u64` can carry, and
    // the bytes have *not* run out. Collapsing them printed "the index ends inside entry 0"
    // about a buffer with nine bytes still unread, which sends whoever holds the file looking
    // for a cut that is not there.
    let (value, used) = decode_u64_leb128(rest).map_err(|damage| match damage {
        VarintError::Overflow => IndexDecodeError::Overlong {
            entry_number: entry,
            field,
        },
        _ => IndexDecodeError::Truncated {
            entry_number: entry,
            field,
        },
    })?;
    *at += used;
    Ok(value)
}

fn take_block_offset(bytes: &[u8], at: &mut usize, entry: usize) -> Result<u64, IndexDecodeError> {
    let field = IndexEntryField::BlockOffset;
    // **The width is stated once, as a type.** It used to be written three times — a
    // `checked_add(8)`, a slice, and an `expect` that the slice was eight bytes — and changing
    // only the first to 4 compiled and turned eight tests into *panics* inside a module whose
    // contract is that a damaged psp is an error and never a panic.
    let chunk: &[u8; 8] =
        bytes
            .get(*at..)
            .and_then(<[u8]>::first_chunk)
            .ok_or(IndexDecodeError::Truncated {
                entry_number: entry,
                field,
            })?;
    *at += chunk.len();
    Ok(u64::from_le_bytes(*chunk))
}

/// Which field of an entry a fault was found in.
///
/// **A closed set rather than a `&'static str`**, so a message cannot name a field the format
/// does not have, and so the three are spelled in one place. They are the three an entry holds,
/// in wire order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexEntryField {
    Contig,
    FirstPosition,
    BlockOffset,
}

impl std::fmt::Display for IndexEntryField {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.write_str(match self {
            IndexEntryField::Contig => "contig",
            IndexEntryField::FirstPosition => "first position",
            IndexEntryField::BlockOffset => "block offset",
        })
    }
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
    #[error("the index ends inside entry {entry_number}, before its {field}")]
    Truncated {
        entry_number: usize,
        field: IndexEntryField,
    },

    /// A variable-length integer with no end in sight: more continuation bytes than any `u64`
    /// can carry.
    ///
    /// **Not a truncation.** The bytes are there and cannot mean a number, which is a different
    /// fault with a different instruction — telling whoever holds the file that it ends where it
    /// does not sends them looking for a cut that is not there.
    #[error("entry {entry_number}'s {field} runs past any number a u64 can hold")]
    Overlong {
        entry_number: usize,
        field: IndexEntryField,
    },

    /// A contig number too large for a [`ContigId`]. **Refused rather than truncated**: a
    /// wrapped contig number indexes a real contig, and the seek would land on it.
    ///
    /// **Named for the one field that can raise it**, because it is the only one: the position
    /// is a `u64` and so is its type, and the offset is read as eight fixed bytes.
    #[error("entry {entry_number} declares contig number {found}, which is not a contig number")]
    ContigNumberTooLarge { entry_number: usize, found: u64 },

    /// The footer's block count and the index's bytes disagree — bytes are left over once that
    /// many entries have been read.
    #[error(
        "the index holds {trailing_bytes} bytes after the {entries_read} entries the footer \
         declares"
    )]
    TrailingBytes {
        trailing_bytes: usize,
        entries_read: usize,
    },

    /// An entry starts before the one before it. **A search over these does not fail on
    /// disorder, it returns the wrong block**, which is why this is refused at open.
    ///
    /// Both entry numbers are carried rather than one computed from the other: `entry - 1` in
    /// the message template panics when rendering an error whose entry is 0.
    #[error(
        "index entry {entry_number} starts at contig {}, position {}, before entry \
         {previous_entry} at contig {}, position {}",
        offered.contig.0, offered.position.get(),
        previous.contig.0, previous.position.get()
    )]
    OutOfOrder {
        entry_number: usize,
        previous_entry: usize,
        previous: GenomePosition,
        offered: GenomePosition,
    },

    /// Two entries at the same byte offset, or offsets that go backwards. Positions may repeat
    /// — a block closed by the byte ceiling is followed by one starting at the same base — but
    /// the offsets are what say two entries are two blocks.
    #[error(
        "index entry {entry_number} is at byte {offered}, which does not follow entry \
         {previous_entry} at byte {previous}"
    )]
    OffsetNotAscending {
        entry_number: usize,
        previous_entry: usize,
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
        let late_on_the_first_contig = index_entry(0, 90_000_000, 4_096);
        let early_on_the_second = index_entry(1, 1, 8_192);
        assert!(late_on_the_first_contig < early_on_the_second);

        let mut blocks = [early_on_the_second, late_on_the_first_contig];
        blocks.sort();
        assert_eq!(blocks[0], late_on_the_first_contig);
    }

    // -----------------------------------------------------------------
    // F1 — the index's own bytes
    // -----------------------------------------------------------------

    fn index_entry(contig: u32, position: u64, block_offset: u64) -> BlockIndexEntry {
        BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(contig),
                position: Position(position),
            },
            block_offset,
        }
    }

    /// A tomato accession's whole index, as spec §3.3 measures it: 154 blocks. Round-tripped
    /// entry for entry, because a search over these is only as good as the coordinates in them.
    #[test]
    fn an_index_round_trips_entry_for_entry() {
        let entries: Vec<BlockIndexEntry> = (0..154)
            .map(|block| {
                index_entry(
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
            index_entry(0, 1, 0),
            index_entry(1, u64::MAX / 2, 1),
            index_entry(u32::MAX, u64::MAX, u64::MAX),
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
        let entries = vec![index_entry(0, 500, 4_096), index_entry(0, 500, 8_192)];
        let bytes = encode_index(&entries);
        assert_eq!(decode_index(&bytes, 2).expect("a legal file"), entries);
    }

    /// A position that goes backwards is refused **wherever in the index it happens**, and the
    /// message names both entries and both coordinates.
    ///
    /// ⚠ **The fixture is eight entries and the break is walked across all seven pairs,
    /// deliberately.** The first version of this test used two entries, so the scan over
    /// consecutive pairs ran exactly one iteration in every test that could fail — and
    /// `.take(1)` on that scan left all thirteen tests green while the decoder accepted an index
    /// whose *third* entry went backwards. On a 154-entry index that is 152 of 153 pairs
    /// unchecked.
    #[test]
    fn a_position_that_goes_backwards_anywhere_in_the_index_is_refused() {
        // **Starting at 1,001 leaves room below every entry.** Starting at 1 made the broken
        // entry *equal* to its predecessor at the first break position, and equal positions are
        // legal — so the case the loop exists to cover was not a step back at all.
        let ascending = |n: u64| index_entry(0, 1_001 + n * 100, 4_096 + n * 4_096);
        for break_at in 1..8usize {
            let mut entries: Vec<BlockIndexEntry> = (0..8).map(ascending).collect();
            // Only the position moves: the offsets stay ascending, so this test can fail on the
            // position rule and on nothing else.
            entries[break_at].first_position.position = Position(1);
            let bytes = encode_index(&entries);
            assert_eq!(
                decode_index(&bytes, entries.len() as u64)
                    .expect_err("a step back must be refused"),
                IndexDecodeError::OutOfOrder {
                    entry_number: break_at,
                    previous_entry: break_at - 1,
                    previous: entries[break_at - 1].first_position,
                    offered: entries[break_at].first_position,
                }
            );
        }
    }

    /// A contig that goes backwards is the same refusal, and it is a separate fixture because
    /// the coordinate that moves is the other half of the key.
    #[test]
    fn an_index_whose_contigs_go_backwards_is_refused() {
        let entries = vec![
            index_entry(0, 10, 4_096),
            index_entry(1, 10, 8_192),
            index_entry(0, 10, 12_288),
        ];
        let bytes = encode_index(&entries);
        assert_eq!(
            decode_index(&bytes, 3).expect_err("a contig going backwards"),
            IndexDecodeError::OutOfOrder {
                entry_number: 2,
                previous_entry: 1,
                previous: entries[1].first_position,
                offered: entries[2].first_position,
            }
        );
    }

    /// **A position resets when a new contig begins, and that is what every real psp looks
    /// like.** Chromosome 2 starts at base 1 again, far below where chromosome 1 ended.
    ///
    /// ⚠ The 154-entry round-trip fixture does not cover this: its positions rise across contig
    /// boundaries as well as within them, so a decoder that compared raw positions and ignored
    /// contigs would pass it. Measured — adding such a rule left all thirteen tests green while
    /// refusing this index.
    #[test]
    fn a_position_that_resets_at_a_new_contig_is_accepted() {
        let entries = vec![
            index_entry(0, 1, 4_096),
            index_entry(0, 90_000_000, 8_192),
            index_entry(1, 1, 12_288),
            index_entry(1, 500, 16_384),
            index_entry(2, 1, 20_480),
        ];
        let bytes = encode_index(&entries);
        assert_eq!(
            decode_index(&bytes, entries.len() as u64).expect("a real multi-contig index"),
            entries
        );
    }

    /// Offsets must strictly increase: two entries at one offset are one block indexed twice,
    /// and an offset that goes backwards points a seek behind where it has already read.
    ///
    /// **This is the half positions cannot carry**, because positions are allowed to repeat.
    /// Walked across the pairs for the reason given on the position test.
    #[test]
    fn an_index_whose_offsets_repeat_or_go_backwards_is_refused() {
        for break_at in 1..8usize {
            for backwards in [false, true] {
                let mut entries: Vec<BlockIndexEntry> = (0..8)
                    .map(|n| index_entry(0, 1 + n * 100, 4_096 + n * 4_096))
                    .collect();
                let previous = entries[break_at - 1].block_offset;
                entries[break_at].block_offset = if backwards { previous - 1 } else { previous };
                let bytes = encode_index(&entries);
                assert_eq!(
                    decode_index(&bytes, entries.len() as u64).expect_err("offsets must ascend"),
                    IndexDecodeError::OffsetNotAscending {
                        entry_number: break_at,
                        previous_entry: break_at - 1,
                        previous,
                        offered: entries[break_at].block_offset,
                    }
                );
            }
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
            index_entry(0, 1, 4_096),
            index_entry(0, 100_001, 8_192),
            index_entry(0, 200_001, 16_384),
        ];
        let entries_left_over = vec![entries[2]];
        let bytes = encode_index(&entries);

        let refused = decode_index(&bytes, 2).expect_err("a tail must be refused");
        match refused {
            IndexDecodeError::TrailingBytes {
                trailing_bytes,
                entries_read,
            } => {
                assert_eq!(entries_read, 2);
                // The third entry's own bytes, computed rather than recalled: a contig varint,
                // a three-byte position varint for 200,001, and the eight-byte offset.
                assert_eq!(trailing_bytes, encode_index(&entries_left_over).len());
                assert_eq!(trailing_bytes, 12);
            }
            other => panic!("expected trailing bytes, got {other}"),
        }

        let refused = decode_index(&bytes, 4).expect_err("a short read must be refused");
        assert_eq!(
            refused,
            IndexDecodeError::Truncated {
                entry_number: 3,
                field: IndexEntryField::Contig,
            }
        );
    }

    /// Every byte prefix of a real index is refused, and **none of them panics**. A file cut
    /// short anywhere inside its index is damage, not a shorter index.
    #[test]
    fn every_truncation_of_an_index_is_refused_without_panicking() {
        let entries = vec![
            index_entry(0, 1, 4_096),
            index_entry(3, 100_001, 8_192),
            index_entry(3, 200_001, 16_384),
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
            IndexDecodeError::ContigNumberTooLarge {
                entry_number: 0,
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
        let bytes = encode_index(&[index_entry(0, 1, 4_096)]);
        let refused = decode_index(&bytes, u64::MAX).expect_err("the file cannot hold them");
        assert!(matches!(
            refused,
            IndexDecodeError::Truncated {
                entry_number: 1,
                ..
            }
        ));
    }

    /// The checksum is over the index's bytes: a golden value, then a bit flipped in each byte.
    ///
    /// ⚠ **The golden value is the half that pins which hash it is.** Without it this test
    /// asserted only avalanche, which any hash satisfies — swapping the body for FNV-1a left all
    /// thirteen tests green while every psp this build wrote carried a checksum no other build
    /// would agree with. The number below is XXH3-64 truncated to its low 32 bits, which is what
    /// production computes and what zstd uses for its own frames.
    #[test]
    fn the_checksum_moves_when_any_byte_of_the_index_does() {
        let entries = vec![index_entry(0, 1, 4_096), index_entry(0, 100_001, 8_192)];
        let bytes = encode_index(&entries);
        let pristine = checksum_index(&bytes);
        // XXH3-64 over these 22 bytes, truncated to its low 32 bits. Swapping the body for
        // FNV-1a gives 1,094,245,170 here and left every other assertion in this test green.
        assert_eq!(pristine, 683_841_834);
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

    /// A variable-length integer with no end is **not** a truncation, and does not say the file
    /// ends where it does not.
    ///
    /// Ten continuation bytes is more than any `u64` can carry. Measured before the split: this
    /// 19-byte buffer, with nine bytes still unread, was reported as *"the index ends inside
    /// entry 0, before its contig"* — which sends whoever holds the file looking for a cut that
    /// is not there. No test reached this varint error at all.
    #[test]
    fn a_variant_length_integer_with_no_end_is_refused_as_overlong_not_truncated() {
        let mut bytes = vec![0x80u8; 10];
        bytes.push(0x01);
        bytes.extend_from_slice(&4_096u64.to_le_bytes());
        assert_eq!(bytes.len(), 19);
        assert_eq!(
            decode_index(&bytes, 1).expect_err("no u64 is that long"),
            IndexDecodeError::Overlong {
                entry_number: 0,
                field: IndexEntryField::Contig,
            }
        );
    }

    /// **Every refusal renders**, including at entry 0, and none of them panics while doing it.
    ///
    /// ⚠ Two of these used to work the earlier entry's number out as `entry - 1` inside the
    /// message template. Rendering one whose entry is 0 panicked with *attempt to subtract with
    /// overflow* — a panic inside `Display`, on the path that reports a damaged file, from the
    /// one type whose whole contract is that a damaged file is never a panic. `decode_index`
    /// cannot produce entry 0 there, but these variants are public and constructible, and this
    /// very test builds them by hand.
    #[test]
    fn every_refusal_renders_at_entry_zero_without_panicking() {
        let somewhere = GenomePosition {
            contig: ContigId(0),
            position: Position(1),
        };
        let rendered: Vec<String> = vec![
            IndexDecodeError::Truncated {
                entry_number: 0,
                field: IndexEntryField::Contig,
            },
            IndexDecodeError::Overlong {
                entry_number: 0,
                field: IndexEntryField::FirstPosition,
            },
            IndexDecodeError::ContigNumberTooLarge {
                entry_number: 0,
                found: u64::MAX,
            },
            IndexDecodeError::TrailingBytes {
                trailing_bytes: 0,
                entries_read: 0,
            },
            IndexDecodeError::OutOfOrder {
                entry_number: 0,
                previous_entry: 0,
                previous: somewhere,
                offered: somewhere,
            },
            IndexDecodeError::OffsetNotAscending {
                entry_number: 0,
                previous_entry: 0,
                previous: 0,
                offered: 0,
            },
        ]
        .into_iter()
        .map(|refusal| refusal.to_string())
        .collect();

        assert_eq!(rendered.len(), 6);
        for message in &rendered {
            assert!(!message.is_empty());
            // A struct dump in a user-facing message is what `{:?}` on a coordinate gives.
            assert!(
                !message.contains("GenomePosition"),
                "a message must read as a sentence, not a struct dump: {message}"
            );
        }
    }

    /// **The messages are the contract**, so they are pinned rather than left to whatever
    /// `thiserror` last rendered.
    ///
    /// ⚠ Nothing pinned them before: cutting `Truncated`'s message down to *"the index is
    /// damaged"*, or mislabelling the position as the contig, left all thirteen tests green.
    /// Only `contig` was ever asserted; the other two field names appeared nowhere.
    #[test]
    fn each_refusal_says_what_is_wrong_and_where() {
        assert_eq!(
            IndexDecodeError::Truncated {
                entry_number: 4,
                field: IndexEntryField::FirstPosition,
            }
            .to_string(),
            "the index ends inside entry 4, before its first position"
        );
        assert_eq!(
            IndexDecodeError::Overlong {
                entry_number: 4,
                field: IndexEntryField::BlockOffset,
            }
            .to_string(),
            "entry 4's block offset runs past any number a u64 can hold"
        );
        assert_eq!(
            IndexDecodeError::ContigNumberTooLarge {
                entry_number: 2,
                found: 4_294_967_296,
            }
            .to_string(),
            "entry 2 declares contig number 4294967296, which is not a contig number"
        );
        assert_eq!(
            IndexDecodeError::TrailingBytes {
                trailing_bytes: 12,
                entries_read: 2,
            }
            .to_string(),
            "the index holds 12 bytes after the 2 entries the footer declares"
        );
        assert_eq!(
            IndexDecodeError::OutOfOrder {
                entry_number: 3,
                previous_entry: 2,
                previous: GenomePosition {
                    contig: ContigId(1),
                    position: Position(900),
                },
                offered: GenomePosition {
                    contig: ContigId(1),
                    position: Position(800),
                },
            }
            .to_string(),
            "index entry 3 starts at contig 1, position 800, before entry 2 at contig 1, \
             position 900"
        );
        assert_eq!(
            IndexDecodeError::OffsetNotAscending {
                entry_number: 3,
                previous_entry: 2,
                previous: 8_192,
                offered: 4_096,
            }
            .to_string(),
            "index entry 3 is at byte 4096, which does not follow entry 2 at byte 8192"
        );
    }

    /// Each field of an entry is named in the message by the same words in every refusal.
    #[test]
    fn every_field_of_an_entry_has_one_spelling() {
        let spellings: Vec<String> = [
            IndexEntryField::Contig,
            IndexEntryField::FirstPosition,
            IndexEntryField::BlockOffset,
        ]
        .iter()
        .map(|field| field.to_string())
        .collect();
        assert_eq!(spellings, ["contig", "first position", "block offset"]);
    }

    /// **Both width constants, against the bytes the encoder actually produces.**
    ///
    /// ⚠ Without this, neither can fail. `LARGEST_ENTRY_BYTES` was 18 for the whole of F1 — five
    /// bytes for the contig, five for the position — and nothing noticed, because a reservation
    /// that is too small only makes a `Vec` grow. But a position is a `u64`, so its
    /// variable-length form runs to ten bytes and the widest entry is 23. The number was carried
    /// over from production, where all four of the index's fields genuinely are `u32`.
    #[test]
    fn the_width_constants_are_the_widths_the_encoder_produces() {
        let widest = encode_index(&[index_entry(u32::MAX, u64::MAX, u64::MAX)]);
        assert_eq!(
            widest.len(),
            23,
            "five varint bytes, ten, and the fixed eight"
        );
        assert!(
            LARGEST_ENTRY_BYTES >= widest.len(),
            "the reservation is sized by a bound the widest entry breaks: \
             {LARGEST_ENTRY_BYTES} against {}",
            widest.len()
        );

        let smallest = encode_index(&[index_entry(0, 0, 0)]);
        assert_eq!(
            smallest.len(),
            10,
            "one varint byte, one, and the fixed eight"
        );
        assert!(
            SMALLEST_ENTRY_BYTES <= smallest.len() as u64,
            "the allocation bound assumes entries no smaller than any entry can be"
        );
    }

    /// The bytes are the format, so their layout is stated once against a literal rather than
    /// only against this module's own decoder. **A change here is a format change**; from
    /// Milestone F it costs a version.
    #[test]
    fn one_entry_encodes_to_these_exact_bytes() {
        let bytes = encode_index(&[index_entry(3, 300, 4_096)]);
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
