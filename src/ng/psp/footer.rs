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

/// The footer's bytes: the five counts and offsets little-endian in field order, then the
/// checksum, then the magic.
///
/// **Fixed width and fixed order**, which is what lets a reader seek [`FOOTER_BYTES`] back from
/// the end of a file and read the whole thing in one go without knowing anything else about it.
///
/// **⚠ This will encode a footer [`decode_footer`] refuses**, and deliberately: the reader's
/// rule that the index ends exactly where the trailer begins is not checked here, because this
/// module's own tests must be able to *write* a footer that breaks it in order to prove the
/// reader refuses one. So the obligation is real and it is not this function's — **it belongs to
/// F3's `finish`**, which is what computes the two offsets from what it has actually written,
/// and which must not be able to produce a file its own reader rejects.
pub fn encode_footer(footer: &Footer) -> [u8; FOOTER_BYTES] {
    // Destructured with no `..`: **a field added to the footer is a compile error here** rather
    // than a field silently left out of every file this build writes, and the footer is the one
    // structure whose width is also a seek distance.
    let Footer {
        index_offset,
        index_bytes,
        trailer_offset,
        trailer_bytes,
        n_blocks,
        index_checksum,
    } = *footer;

    let mut bytes = [0u8; FOOTER_BYTES];
    let mut at = 0;
    for value in [
        index_offset,
        index_bytes,
        trailer_offset,
        trailer_bytes,
        n_blocks,
    ] {
        bytes[at..at + 8].copy_from_slice(&value.to_le_bytes());
        at += 8;
    }
    bytes[at..at + 4].copy_from_slice(&index_checksum.to_le_bytes());
    at += 4;
    bytes[at..].copy_from_slice(&FOOTER_MAGIC);
    bytes
}

/// Read the footer back, refusing one that cannot describe a file.
///
/// **The magic is checked before any offset is believed**, which is the whole reason it sits
/// last: a file that is not an ng psp, or one a killed writer never finished, is rejected on
/// four bytes rather than on arithmetic over numbers that mean nothing.
///
/// **What it checks is what a footer can check about itself**: the magic, that the sections it
/// locates do not overflow the address space, and that the index and the trailer abut — the
/// index ends exactly where the trailer begins, which is the layout spec §3 fixes and what
/// makes replacing a trailer a truncate-and-write-forward (spec §6.5).
///
/// **What it does not check, and where those live.** Nothing here knows the file's length, so
/// nothing here can say the trailer ends where the footer begins, nor that the index lies inside
/// the file at all; that is `open`'s, which has the length. Nor does it check the block count
/// against the index's byte length: `decode_index` already requires the bytes to hold exactly
/// that many entries and no more, which is the exact form of the same check.
pub fn decode_footer(bytes: &[u8; FOOTER_BYTES]) -> Result<Footer, FooterDecodeError> {
    let mut magic = [0u8; 4];
    magic.copy_from_slice(&bytes[FOOTER_BYTES - FOOTER_MAGIC.len()..]);
    if magic != FOOTER_MAGIC {
        return Err(FooterDecodeError::NotAFooter {
            found: magic,
            expected: FOOTER_MAGIC,
        });
    }

    // **PANIC-FREE: both `expect`s below are unreachable by type, and that is the difference
    // from `index.rs`.** This function's argument is a fixed-size `[u8; FOOTER_BYTES]`, so every
    // window taken from it has a length the compiler knows; the index's decoder worked on a
    // slice whose length was a runtime fact, which is why the same shape there was replaced by
    // `first_chunk` after a mutation turned eight of its tests into panics.
    let mut counts = [0u64; 5];
    for (number, slot) in counts.iter_mut().enumerate() {
        let at = number * 8;
        *slot = u64::from_le_bytes(
            bytes[at..at + 8]
                .try_into()
                .expect("an eight-byte window of a fixed array is eight bytes"),
        );
    }
    let [
        index_offset,
        index_bytes,
        trailer_offset,
        trailer_bytes,
        n_blocks,
    ] = counts;
    let index_checksum = u32::from_le_bytes(
        bytes[5 * 8..5 * 8 + 4]
            .try_into()
            .expect("a four-byte window of a fixed array is four bytes"),
    );

    let index_ends = index_offset.checked_add(index_bytes).ok_or(
        FooterDecodeError::SectionEndIsPastAnyFile {
            section: FileSection::BlockIndex,
            offset: index_offset,
            bytes: index_bytes,
        },
    )?;
    trailer_offset.checked_add(trailer_bytes).ok_or(
        FooterDecodeError::SectionEndIsPastAnyFile {
            section: FileSection::Trailer,
            offset: trailer_offset,
            bytes: trailer_bytes,
        },
    )?;
    if index_ends != trailer_offset {
        return Err(FooterDecodeError::IndexDoesNotEndWhereTheTrailerBegins {
            index_ends,
            trailer_starts: trailer_offset,
        });
    }

    Ok(Footer {
        index_offset,
        index_bytes,
        trailer_offset,
        trailer_bytes,
        n_blocks,
        index_checksum,
    })
}

/// Which of the two sections a footer locates.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileSection {
    BlockIndex,
    Trailer,
}

impl std::fmt::Display for FileSection {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.write_str(match self {
            FileSection::BlockIndex => "block index",
            FileSection::Trailer => "trailer",
        })
    }
}

/// Why a footer could not be read.
///
/// **Its own type, carrying no path**, like this module's other codec errors: the file's name
/// belongs to whoever opened the file, and F4's `open` is what dresses one of these as a
/// [`PspReadError`](super::PspReadError).
///
/// **The first variant is not the same instruction as the others.** No magic means *this file
/// was never finished — rebuild it*, which is goal 3 and the everyday case of a killed pileup;
/// the rest mean *this file is damaged*. `open` maps them apart.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FooterDecodeError {
    /// The last four bytes are not this format's tail magic: the writer never finished, or the
    /// file is not an ng psp at all.
    #[error(
        "the file does not end with {}; it was never finished, or it is not an ng psp (its \
         last bytes are {found:02x?})",
        String::from_utf8_lossy(expected.as_slice()).escape_debug()
    )]
    NotAFooter { found: [u8; 4], expected: [u8; 4] },

    /// A section's offset plus its length is a number no file can reach.
    #[error(
        "the footer puts the {section} at byte {offset} with {bytes} bytes, which is past any \
         file; this psp is damaged"
    )]
    SectionEndIsPastAnyFile {
        section: FileSection,
        offset: u64,
        bytes: u64,
    },

    /// The index does not end where the trailer begins. **They abut by construction** (spec §3),
    /// and a gap or an overlap means one of the two offsets is not what the writer wrote.
    #[error(
        "the block index ends at byte {index_ends} but the trailer starts at byte \
         {trailer_starts}; in a whole psp the trailer begins where the index ends, so this one \
         is damaged"
    )]
    IndexDoesNotEndWhereTheTrailerBegins {
        index_ends: u64,
        trailer_starts: u64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------
    // F2 — the footer's own bytes
    // -----------------------------------------------------------------

    /// A footer round-trips field for field, including a trailer of no bytes.
    #[test]
    fn a_footer_round_trips_field_for_field() {
        for footer in [
            // A trailer that holds something — what every finished psp has.
            a_footer(),
            // And one that is empty, which is legal too.
            Footer {
                trailer_bytes: 0,
                ..a_footer()
            },
            // A file holding no blocks at all: the index is empty and the trailer begins where
            // it began.
            Footer {
                index_offset: 4_096,
                index_bytes: 0,
                trailer_offset: 4_096,
                trailer_bytes: 0,
                n_blocks: 0,
                index_checksum: 0,
            },
        ] {
            let bytes = encode_footer(&footer);
            assert_eq!(decode_footer(&bytes).expect("its own bytes decode"), footer);
        }
    }

    /// The widest value of every field, so a byte that is dropped or transposed shows up.
    #[test]
    fn the_widest_value_of_every_field_round_trips() {
        let footer = Footer {
            index_offset: u64::MAX - 2,
            index_bytes: 1,
            trailer_offset: u64::MAX - 1,
            trailer_bytes: 1,
            n_blocks: u64::MAX,
            index_checksum: u32::MAX,
        };
        let bytes = encode_footer(&footer);
        assert_eq!(decode_footer(&bytes).expect("the extremes decode"), footer);
    }

    /// **The bytes are the format**, so the layout is stated once against a literal rather than
    /// only against this module's own decoder. A change here is a format change.
    #[test]
    fn a_footer_encodes_to_these_exact_bytes() {
        let bytes = encode_footer(&Footer {
            index_offset: 4_096,
            index_bytes: 2_464,
            trailer_offset: 6_560,
            trailer_bytes: 0,
            n_blocks: 154,
            index_checksum: 0xDEAD_BEEF,
        });
        assert_eq!(
            bytes,
            [
                0x00, 0x10, 0, 0, 0, 0, 0, 0, // index at byte 4,096
                0xA0, 0x09, 0, 0, 0, 0, 0, 0, // 2,464 bytes of index
                0xA0, 0x19, 0, 0, 0, 0, 0,
                0, // trailer at byte 6,560 — where the index ends
                0, 0, 0, 0, 0, 0, 0, 0, // an empty trailer
                0x9A, 0, 0, 0, 0, 0, 0, 0, // 154 blocks
                0xEF, 0xBE, 0xAD, 0xDE, // the index's checksum
                b'N', b'G', b'P', b'E', // and the magic, last
            ]
        );
    }

    /// **The magic sits in the last four bytes**, which is what lets a reader reject a foreign
    /// or unfinished file on a four-byte read at end-of-file before touching anything else.
    #[test]
    fn the_magic_is_the_last_four_bytes() {
        let bytes = encode_footer(&a_footer());
        assert_eq!(&bytes[FOOTER_BYTES - 4..], &FOOTER_MAGIC);
        assert_eq!(&FOOTER_MAGIC, b"NGPE");
    }

    /// A file that does not end with the magic is refused, and **the message says the two things
    /// that are worth telling apart**: a run that was killed before it finished, and a file that
    /// was never an ng psp.
    ///
    /// This is the refusal goal 3 rests on. A psp read short is a sample that stops in the middle
    /// of a chromosome and says nothing about it.
    #[test]
    fn a_file_that_does_not_end_with_the_magic_is_refused() {
        for wrong in [
            *b"\0\0\0\0", // a writer killed before it wrote anything here
            *b"PSPE",     // production's tail — the everyday wrong file
            *b"NGP\n",    // this format's *head* magic, which a truncation could leave
            *b"EPGN",     // the right bytes the wrong way round
        ] {
            let mut bytes = encode_footer(&a_footer());
            bytes[FOOTER_BYTES - 4..].copy_from_slice(&wrong);
            assert_eq!(
                decode_footer(&bytes).expect_err("that is not a footer"),
                FooterDecodeError::NotAFooter {
                    found: wrong,
                    expected: FOOTER_MAGIC,
                }
            );
        }
    }

    /// **The magic is checked before any offset is believed.** A footer whose numbers are
    /// nonsense *and* whose magic is wrong must be refused for the magic, because that is the
    /// answer a caller can act on — rebuild the file — where the arithmetic would only say it is
    /// damaged.
    #[test]
    fn the_magic_is_checked_before_the_offsets_are_believed() {
        let mut bytes = encode_footer(&Footer {
            index_offset: u64::MAX,
            index_bytes: u64::MAX,
            ..a_footer()
        });
        bytes[FOOTER_BYTES - 4..].copy_from_slice(b"PSPE");
        assert!(matches!(
            decode_footer(&bytes),
            Err(FooterDecodeError::NotAFooter { .. })
        ));
    }

    /// A section whose offset plus length is past any address is refused rather than wrapped.
    /// **Wrapping would give a small number that looks like a legal offset.**
    #[test]
    fn a_section_that_runs_past_any_address_is_refused() {
        let index = encode_footer(&Footer {
            index_offset: u64::MAX,
            index_bytes: 1,
            ..a_footer()
        });
        assert_eq!(
            decode_footer(&index).expect_err("that overflows"),
            FooterDecodeError::SectionEndIsPastAnyFile {
                section: FileSection::BlockIndex,
                offset: u64::MAX,
                bytes: 1,
            }
        );

        let trailer = encode_footer(&Footer {
            index_offset: 0,
            index_bytes: u64::MAX,
            trailer_offset: u64::MAX,
            trailer_bytes: 2,
            ..a_footer()
        });
        assert_eq!(
            decode_footer(&trailer).expect_err("that overflows"),
            FooterDecodeError::SectionEndIsPastAnyFile {
                section: FileSection::Trailer,
                offset: u64::MAX,
                bytes: 2,
            }
        );
    }

    /// The index ends exactly where the trailer begins. A gap or an overlap means one of the two
    /// offsets is not what the writer wrote — and **replacing a trailer truncates at that
    /// offset**, so a gap would leave stale bytes inside the file and an overlap would eat the
    /// end of the index.
    #[test]
    fn an_index_and_trailer_that_do_not_abut_are_refused() {
        for (index_bytes, trailer_offset, ends) in [
            (2_464u64, 6_561u64, 6_560u64), // a one-byte gap
            (2_464, 6_559, 6_560),          // a one-byte overlap
            (0, 6_560, 4_096),              // an empty index that does not reach the trailer
        ] {
            let bytes = encode_footer(&Footer {
                index_offset: 4_096,
                index_bytes,
                trailer_offset,
                ..a_footer()
            });
            assert_eq!(
                decode_footer(&bytes).expect_err("they must abut"),
                FooterDecodeError::IndexDoesNotEndWhereTheTrailerBegins {
                    index_ends: ends,
                    trailer_starts: trailer_offset,
                }
            );
        }
    }

    /// **When both sections overflow, the index is the one reported** — the checks run in wire
    /// order, so the first fault a reader meets is the first fault it names.
    ///
    /// Nothing pinned this: swapping the two checks changed which section a damaged footer was
    /// reported against, and every test still passed, because no fixture made both fail at once.
    #[test]
    fn the_first_section_to_overflow_is_the_one_named() {
        let bytes = encode_footer(&Footer {
            index_offset: u64::MAX,
            index_bytes: 3,
            trailer_offset: u64::MAX,
            trailer_bytes: 4,
            ..a_footer()
        });
        assert_eq!(
            decode_footer(&bytes).expect_err("both overflow"),
            FooterDecodeError::SectionEndIsPastAnyFile {
                section: FileSection::BlockIndex,
                offset: u64::MAX,
                bytes: 3,
            }
        );
    }

    /// **No 48 bytes make this panic**, which is the module's standing rule and a class it has
    /// shipped once already: a corrupt psp is data a run was handed, not a bug.
    ///
    /// The argument is exactly 48 bytes, so the whole input domain is reachable. Half the draws
    /// carry the real magic, because the arithmetic past the magic check is only reached by
    /// those and a uniform draw would almost never get there.
    #[test]
    fn no_forty_eight_bytes_make_the_decoder_panic() {
        // A cheap deterministic generator: no dependency, and the same draws every run, which is
        // what makes a failure reproducible from the seed alone.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let (mut accepted, mut refused) = (0u32, 0u32);
        for draw in 0..20_000 {
            let mut bytes = [0u8; FOOTER_BYTES];
            for chunk in bytes.chunks_mut(8) {
                let word = next().to_le_bytes();
                chunk.copy_from_slice(&word[..chunk.len()]);
            }
            if draw % 2 == 0 {
                bytes[FOOTER_BYTES - 4..].copy_from_slice(&FOOTER_MAGIC);
            }
            match decode_footer(&bytes) {
                Ok(_) => accepted += 1,
                Err(_) => refused += 1,
            }
        }
        assert_eq!(accepted + refused, 20_000);
        assert!(
            refused > 0 && accepted < 20_000,
            "the draws must reach both outcomes: {accepted} accepted, {refused} refused"
        );
    }

    /// Every refusal renders as a sentence naming what a caller must act on.
    #[test]
    fn each_refusal_says_what_is_wrong() {
        assert_eq!(
            FooterDecodeError::NotAFooter {
                found: *b"PSPE",
                expected: FOOTER_MAGIC,
            }
            .to_string(),
            "the file does not end with NGPE; it was never finished, or it is not an ng psp \
             (its last bytes are [50, 53, 50, 45])"
        );
        assert_eq!(
            FooterDecodeError::SectionEndIsPastAnyFile {
                section: FileSection::Trailer,
                offset: 18_446_744_073_709_551_615,
                bytes: 2,
            }
            .to_string(),
            "the footer puts the trailer at byte 18446744073709551615 with 2 bytes, which is \
             past any file; this psp is damaged"
        );
        assert_eq!(
            FooterDecodeError::IndexDoesNotEndWhereTheTrailerBegins {
                index_ends: 6_560,
                trailer_starts: 6_561,
            }
            .to_string(),
            "the block index ends at byte 6560 but the trailer starts at byte 6561; in a whole \
             psp the trailer begins where the index ends, so this one is damaged"
        );
    }

    /// Both sections are named by the same words wherever they appear.
    #[test]
    fn every_section_has_one_spelling() {
        assert_eq!(FileSection::BlockIndex.to_string(), "block index");
        assert_eq!(FileSection::Trailer.to_string(), "trailer");
    }

    /// **Every single-byte change to a well-formed footer is either refused or decodes to a
    /// different footer — none is silently ignored.** A byte nothing reads would be a byte a
    /// writer could leave uninitialised, which is how two files that should be identical stop
    /// being so (spec §7's worker-count invariance).
    #[test]
    fn no_byte_of_the_footer_is_ignored() {
        let pristine = encode_footer(&a_footer());
        // ⚠ **Compared against what the pristine bytes decode to, not against the footer they
        // were encoded from.** Written the other way, a decoder that masked a byte off still
        // passed whenever the fixture's value had no zeros there — measured: masking the top
        // byte of `n_blocks` was caught, because 154 has zeros there, and masking the top byte
        // of `index_checksum` was not, because 0xDEADBEEF does not.
        let as_decoded = decode_footer(&pristine).expect("the fixture is a valid footer");
        let mut refused = 0;
        for byte in 0..FOOTER_BYTES {
            let mut damaged = pristine;
            damaged[byte] ^= 0xFF;
            match decode_footer(&damaged) {
                Err(_) => refused += 1,
                Ok(read_back) => assert_ne!(
                    read_back, as_decoded,
                    "byte {byte} changed and the footer decoded to the same values"
                ),
            }
        }
        assert!(refused > 0, "some damage must be refused outright");
    }

    /// The footer is read by seeking that many bytes back from the end of the file, so its
    /// width is part of the format. **This adds up the fields' own wire widths** rather than
    /// asking for `size_of::<Footer>()`: a sixth field of four bytes lands in the struct's
    /// tail padding, so `size_of` does not move while the bytes on disk do — and every reader
    /// then seeks four bytes into the middle of the footer it was looking for.
    #[test]
    fn the_footer_constant_is_the_width_of_the_fields_it_stands_for() {
        // **Destructured, not read field by field.** Reading fields off a value is unaffected by
        // the struct gaining one, so this test passed with a seventh field that reached no file:
        // the encoder's own no-`..` destructure forces a new field to be *mentioned*, and
        // rustc's suggested repair — `trailer_checksum: _` — compiles clean and writes nothing.
        // This is the one place the constant is tied to the field *set* rather than to a value.
        let Footer {
            index_offset,
            index_bytes,
            trailer_offset,
            trailer_bytes,
            n_blocks,
            index_checksum,
        } = a_footer();
        let on_the_wire = index_offset.to_le_bytes().len()
            + index_bytes.to_le_bytes().len()
            + trailer_offset.to_le_bytes().len()
            + trailer_bytes.to_le_bytes().len()
            + n_blocks.to_le_bytes().len()
            + index_checksum.to_le_bytes().len()
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

    /// A tomato accession's footer, as spec §3.3 measures it: 154 blocks, and a trailer that
    /// **holds something**.
    ///
    /// ⚠ It carried `trailer_bytes: 0` at first, and so did every other fixture that reached an
    /// `Ok` — the round-trip shapes, the widest-value shape, and the three non-abut cases that
    /// spread this one. So **no well-formed footer with a real trailer was ever decoded**, which
    /// is the shape every finished psp has: the trailer is the writer's closing payload (spec
    /// §3.4). Two wrong abut rules passed the whole suite because of it — one folding
    /// `trailer_bytes` into the comparison, one checking only when the trailer is empty.
    fn a_footer() -> Footer {
        Footer {
            index_offset: 4_096,
            index_bytes: 2_464,
            trailer_offset: 6_560,
            trailer_bytes: 1_288,
            n_blocks: 154,
            index_checksum: 0xDEAD_BEEF,
        }
    }
}
