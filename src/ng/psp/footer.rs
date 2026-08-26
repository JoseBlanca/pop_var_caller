//! The footer: the fixed tail whose presence is what says the file is complete.
//!
//! It carries where the index and the trailer are, how many blocks the file holds, a
//! checksum over the index, and a magic **placed last**, so a four-byte read at
//! end-of-file rejects a truncated or foreign file before anything else is touched.
//!
//! **A file with no valid footer is refused rather than read short**, and that is the whole
//! reason it exists: it is the only signal that distinguishes a completed file from a run
//! that was killed, and a caller reading one short would silently get a sample that stops
//! in the middle of a chromosome (spec §3.3, goal 3).
//!
//! Production's 32-byte tail is the model (`src/psp/trailer.rs` — which production calls a
//! *trailer*, a word this module uses for something else entirely). This one is wider,
//! because it has to locate a section production has no equivalent of: the trailer, the
//! writer's closing payload.

/// Bytes the footer occupies: five `u64` offsets and counts, the index checksum, and the
/// magic.
pub const FOOTER_BYTES: usize = 5 * 8 + 4 + FOOTER_MAGIC.len();

/// The magic at the very end of a finished file — `NGPE`, ng psp end.
///
/// **Last, so a four-byte read at end-of-file rejects a truncated or foreign file before
/// anything else is touched**, and different from the head magic so that a truncation which
/// happened to copy the head's bytes would not pass the tail check.
pub const FOOTER_MAGIC: [u8; 4] = *b"NGPE";

/// The fixed tail of a finished file.
///
/// Field order is the wire order, and the magic follows them; a decoder checks the magic
/// before it believes any of the offsets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Footer {
    /// Where the block index starts. **Also where the blocks end**, which is what makes an
    /// append cheap: a writer reopening a file truncates here and carries on (spec §3).
    pub index_offset: u64,
    /// How many bytes the block index occupies.
    pub index_bytes: u64,
    /// Where the trailer starts. **The index sits before it deliberately**, so replacing a
    /// trailer means truncating here and writing forward, leaving the blocks and the index
    /// untouched (spec §3, §6.5).
    pub trailer_offset: u64,
    /// How many bytes the trailer occupies. Zero is legal — the payload may be empty.
    pub trailer_bytes: u64,
    /// How many psp blocks the file holds, and therefore how many entries the index has.
    pub n_blocks: u64,
    /// Checksum over the index's encoded bytes.
    pub index_checksum: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The footer is read by seeking that many bytes back from the end of the file, so its
    /// width is part of the format. **This adds up the fields' own wire widths** rather than
    /// asking for `size_of::<Footer>()`: a sixth field of four bytes lands in the struct's
    /// tail padding, so `size_of` does not move while the bytes on disk do — and every reader
    /// then seeks four bytes into the middle of the footer it was looking for.
    #[test]
    fn the_footer_constant_is_the_width_of_the_fields_it_stands_for() {
        let footer = a_footer();
        let on_the_wire = footer.index_offset.to_le_bytes().len()
            + footer.index_bytes.to_le_bytes().len()
            + footer.trailer_offset.to_le_bytes().len()
            + footer.trailer_bytes.to_le_bytes().len()
            + footer.n_blocks.to_le_bytes().len()
            + footer.index_checksum.to_le_bytes().len()
            + FOOTER_MAGIC.len();
        assert_eq!(FOOTER_BYTES, on_the_wire);
        assert_eq!(FOOTER_BYTES, 48);
    }

    /// Production's tail is 32 bytes and this one is wider, because it locates a section
    /// production has no equivalent of. Stated as prose against a literal rather than as
    /// arithmetic on production's constant: `FOOTER_BYTES` is not a port of that number, so a
    /// change to production's tail should not fail an ng test.
    #[test]
    fn the_footer_is_two_offsets_wider_than_productions_tail() {
        assert_eq!(FOOTER_BYTES, 32 + 2 * 8);
    }

    /// A footer holds only offsets and counts, so it can be read into a fixed buffer and
    /// copied out without an allocation.
    #[test]
    fn a_footer_is_plain_data_and_holds_no_allocation() {
        let footer = a_footer();
        let copied = footer;
        assert_eq!(copied, footer);
        assert!(std::mem::size_of::<Footer>() <= FOOTER_BYTES);
    }

    fn a_footer() -> Footer {
        Footer {
            index_offset: 4_096,
            index_bytes: 2_464,
            trailer_offset: 6_560,
            trailer_bytes: 0,
            n_blocks: 154,
            index_checksum: 0xDEAD_BEEF,
        }
    }
}
