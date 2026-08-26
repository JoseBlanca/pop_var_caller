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
pub const FOOTER_BYTES: usize = 5 * 8 + 4 + 4;

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
    /// width is part of the format rather than a consequence of the struct's layout. This
    /// pins the two together: a field added without widening the constant would make every
    /// reader seek to the wrong place.
    #[test]
    fn the_footer_is_forty_eight_bytes_and_wider_than_productions() {
        assert_eq!(FOOTER_BYTES, 48);
        assert_eq!(
            FOOTER_BYTES,
            crate::psp::trailer::TRAILER_BYTES + 16,
            "two more offsets than production's tail, because this format has a trailer \
             section to locate and production has none"
        );
    }

    /// A footer holds only offsets and counts, so it can be read into a fixed buffer and
    /// copied out without an allocation.
    #[test]
    fn a_footer_is_plain_data_and_holds_no_allocation() {
        let footer = Footer {
            index_offset: 4_096,
            index_bytes: 2_464,
            trailer_offset: 6_560,
            trailer_bytes: 0,
            n_blocks: 154,
            index_checksum: 0xDEAD_BEEF,
        };
        let copied = footer;
        assert_eq!(copied, footer);
        assert!(std::mem::size_of::<Footer>() <= FOOTER_BYTES);
    }
}
