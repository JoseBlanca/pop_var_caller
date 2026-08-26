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

use crate::ng::types::GenomePosition;

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
}
