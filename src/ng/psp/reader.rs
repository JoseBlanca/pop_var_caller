//! Opening a psp: the footer, the index and the header — and no block at all.
//!
//! **What opening costs is what a cohort pays per sample before it reads anything** (spec §6.2),
//! multiplied by the cohort size, so it is deliberately three reads at fixed places rather than
//! a walk: the fixed tail, the index it points at, and the plain-text header at the front.
//!
//! **A file with no valid footer is refused rather than read short** (spec §3.3, goal 3). That
//! is the only thing distinguishing a run that was killed from a sample that genuinely covers
//! less of the genome, and reading one short would hand a caller a chromosome that stops in the
//! middle with nothing said.
//!
//! **The record walks a caller asks for start here and are built in [`super::walk`]**, which is
//! where everything that can inflate a frame lives. This file's part of a walk is the index
//! lookup, the seek, and the bound that says where the blocks end; that it names no block-
//! decoding code at all is what `the_opener_cannot_reach_any_block_decoding_code` reads its
//! imports to check.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use super::footer::{FOOTER_BYTES, Footer, decode_footer};
use super::header::{HEAD_MAGIC, Header, MAX_LOOK_BACK_WINDOW_LOG};
use super::index::BlockIndexEntry;
use super::walk::{self, RecordIter, SelectiveRecordIter};
// **The head through the module root, which is where its public path already is.** Naming the
// type a predicate is shown is not importing a decoder — but a path into the record module is
// what `the_opener_cannot_reach_any_block_decoding_code` forbids this file, comments included,
// so the type comes through `super` and this note spells no such path.
use super::{PspReadError, RecordHead, index, read_header_from};
use crate::ng::types::GenomePosition;

/// How large a compressor look-back window this reader will hold, unless told otherwise.
///
/// **256 kB, derived from the spec's own figures rather than chosen.** The budget is 500 kB an
/// open sample (spec §1.1, §7). Of that, zstd's decoder floor is about 190 kB (§5.3) and the
/// reader's two buffers are 16 kB each (§4.4). **The 190 kB already contains a 32 kB window** —
/// §5.3 says so in its first line — so a wider window costs only the difference, and the room
/// for one is `500 − 190 − 16 − 16 + 32 = 310 kB`. 2^18 is the largest power of two under that.
///
/// ⚠ **An earlier version of this arithmetic charged the window twice** and cited §5.3 for the
/// buffer figure, which is §4.4's. It reached 278 kB, and 2^18 is the largest power of two under
/// that too — a wrong derivation for a right number, which is the kind that survives review by
/// looking finished.
///
/// It is **eight times the window this build's writer produces** (2^15 = 32 kB), so an ordinary
/// psp opens with room to spare, and a file written with a deliberately wide window is refused
/// with a number the operator can act on rather than a zstd error code. Spec §7's table names the
/// two ways out: *raise the budget, or rewrite the file*.
///
/// **⚠ Two things this is not.** It is arithmetic on the spec's figures, not a measurement of
/// this reader; and the budget it is carved out of does not account for the block index, which
/// `open` also holds — at spec §3.3's ~14,000 entries for a whole genome and 24 bytes an entry
/// that is about 336 kB, more than half the budget on its own. **Milestone H4 measures what an
/// open sample actually costs**, and is where both of those are settled.
pub const DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES: u64 = 1 << 18;

/// A finished psp, open: its header, its block index, and where everything in it is.
///
/// **No block has been touched** — that is the property spec §6.2 fixes and the reason a cohort
/// can hold thousands of these open. Records come from [`records`](Self::records) and its two
/// siblings, each of which turns a starting point into a [`RecordIter`] over the file bounded at
/// the end of its blocks; the decoding itself is [`super::BlockStream`]'s.
#[derive(Debug)]
pub struct PspReader {
    path: PathBuf,
    file: File,
    header: Header,
    footer: Footer,
    block_index: Vec<BlockIndexEntry>,
    /// The ceiling this reader was opened under. **Kept because an open reader should be able to
    /// say what it holds**: at several thousand open samples the budget is the number that
    /// multiplies, and a caller tuning it needs to read back what it actually got.
    window_budget_bytes: u64,
    /// How far a walk's rolling buffer may grow for **one record**, which is a different
    /// ceiling from the window above and is spent only while a walk is running.
    record_buffer_ceiling_bytes: usize,
}

impl PspReader {
    /// Open a finished psp: footer, then index, then header. **No block is touched.**
    ///
    /// The order is the one spec §6.2 fixes, and it is forced by the layout: the footer is the
    /// only part at a known place, and it is what says where the index is.
    pub fn open(path: &Path) -> Result<Self, PspReadError> {
        Self::open_with_a_look_back_window_budget(path, DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES)
    }

    /// The same, with a different ceiling on the look-back window this reader will hold.
    ///
    /// **This is the knob [`DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES`] is the default of**, and it
    /// exists because spec §4.2 makes the fix for a too-wide window *a setting rather than a
    /// rebuilt file*.
    pub fn open_with_a_look_back_window_budget(
        path: &Path,
        window_budget_bytes: u64,
    ) -> Result<Self, PspReadError> {
        let mut file = File::open(path).map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "opening the file",
            source,
        })?;
        let file_bytes = file
            .metadata()
            .map_err(|source| PspReadError::Io {
                path: path.to_path_buf(),
                while_doing: "measuring the file",
                source,
            })?
            .len();

        let footer = Self::read_footer(path, &mut file, file_bytes)?;
        let blocks = Self::read_index(path, &mut file, &footer, file_bytes)?;
        // **From the handle already open**, rather than opening the file a second time: a
        // cohort opens thousands of these, and a second `open(2)` per sample buys nothing.
        let for_the_header = file.try_clone().map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "reading the header",
            source,
        })?;
        let (header, header_bytes) = read_header_from(for_the_header, path)?;

        // **Every block must start inside the blocks**, which begin where the header ends and
        // end where the index begins. An offset is where a reader would seek: one pointing into
        // the header would inflate the header, and one at or past the index would inflate the
        // index.
        //
        // ⚠ The lower bound was missing. An entry at byte 0 opened, and the refusal was deferred
        // to a corrupt block at read time — the wrong instruction, arriving after a cohort had
        // committed to the sample.
        // **And the index itself must start after the header**, which is a rule about the file
        // and not about any entry — so the loop below cannot carry it: **on an empty index there
        // are no entries to check.**
        //
        // ⚠ That is not hypothetical. A footer saying `index_offset = 4, index_bytes = 0,
        // n_blocks = 0` passed every check this reader made, and `PspWriter::append` — which
        // truncates at exactly this offset — then cut a 3,742-byte psp down to 109 bytes and
        // reported success. The four-byte bound `read_and_check_the_footer` applies is what 48
        // bytes can say on their own; this is what the header's length adds, and it subsumes it.
        if footer.index_offset < header_bytes as u64 {
            return Err(PspReadError::damaged(
                path,
                format!(
                    "the footer puts the block index at byte {}, which is inside the                      {header_bytes}-byte header",
                    footer.index_offset
                ),
            ));
        }
        for entry in &blocks {
            let inside = (header_bytes as u64..footer.index_offset).contains(&entry.block_offset);
            if !inside {
                return Err(PspReadError::damaged(
                    path,
                    format!(
                        "the index puts a block at byte {}; the blocks run from byte \
                         {header_bytes} to byte {}",
                        entry.block_offset, footer.index_offset
                    ),
                ));
            }
        }

        let needed_bytes = 1u64 << header.manifest.look_back_window_log;
        if needed_bytes > window_budget_bytes {
            return Err(PspReadError::WindowTooLarge {
                path: path.to_path_buf(),
                needed_bytes,
                allowed_bytes: window_budget_bytes,
            });
        }

        Ok(Self {
            path: path.to_path_buf(),
            file,
            header,
            footer,
            block_index: blocks,
            window_budget_bytes,
            record_buffer_ceiling_bytes: walk::DEFAULT_RECORD_BUFFER_CEILING_BYTES,
        })
    }

    /// The same reader, with a different ceiling on how much of a walk's rolling buffer **one
    /// record** may have.
    ///
    /// **This is the knob [`PspReadError::RecordLargerThanTheReaderAllows`] names**, and spec
    /// §7 is why it exists: a genuine record can be larger than any fixed budget — §8 refuses
    /// to fix a maximum record size in the format — so that refusal's instruction is *raise the
    /// ceiling*, and an instruction with nothing to turn is not an instruction. It is the
    /// record-shaped sibling of [`open_with_a_look_back_window_budget`](Self::open_with_a_look_back_window_budget).
    ///
    /// A ceiling at or under the rolling buffer itself is refused **here, where the setting is
    /// made**, rather than at the record that would have tripped over it.
    /// **The rule is the walk's, not this file's**: the buffer the ceiling is a ceiling on
    /// belongs to the walk, and `reader.rs` may not name the module that owns it.
    pub fn with_a_record_buffer_ceiling(mut self, ceiling: usize) -> Result<Self, PspReadError> {
        walk::check_a_record_buffer_ceiling(&self.path, ceiling)?;
        self.record_buffer_ceiling_bytes = ceiling;
        Ok(self)
    }

    /// How much of a walk's rolling buffer one record may have, under this reader.
    pub fn record_buffer_ceiling_bytes(&self) -> usize {
        self.record_buffer_ceiling_bytes
    }

    /// The file's own account of how it was written, the manifest included.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// One entry per block, in genomic order — the cheap survey, with no block decompressed.
    ///
    /// **`block_index` and not `blocks`**: every other document about this format calls this
    /// section the *block index* (spec §3.3), and a caller reading `psp.blocks()` beside
    /// `psp.records()` reasonably expects the blocks themselves.
    pub fn block_index(&self) -> &[BlockIndexEntry] {
        &self.block_index
    }

    /// Where each section of the file is.
    pub fn footer(&self) -> &Footer {
        &self.footer
    }

    /// The look-back window ceiling this reader was opened under.
    pub fn look_back_window_budget_bytes(&self) -> u64 {
        self.window_budget_bytes
    }

    /// The writer's closing payload: one seek and one read, and it may be empty.
    ///
    /// **Opaque here.** The container stores bytes and hands them back; what is in them is the
    /// writer's business (spec §3.4).
    pub fn trailer(&mut self) -> Result<Vec<u8>, PspReadError> {
        let mut payload = vec![0u8; self.footer.trailer_bytes as usize];
        if payload.is_empty() {
            return Ok(payload);
        }
        self.seek_to(self.footer.trailer_offset, "seeking to the trailer")?;
        self.file
            .read_exact(&mut payload)
            .map_err(|source| PspReadError::Io {
                path: self.path.clone(),
                while_doing: "reading the trailer",
                source,
            })?;
        Ok(payload)
    }

    /// Every record in the file, from its first block.
    ///
    /// **This is where a psp starts costing more than an open file**: until now nothing has
    /// been decompressed, and the walk this returns holds two buffers and a decoder's state
    /// (spec §5.1). Nothing it holds is a function of the block size, the depth, or the length
    /// of the genome.
    pub fn records(&mut self) -> Result<RecordIter<'_>, PspReadError> {
        // **A sample with no records has no blocks**, and the walk over it is empty rather than
        // refused: `finish` writes such a file and `open` accepts it. Starting at the index puts
        // the bound below at zero bytes, which ends the walk before it reads anything.
        let at = self
            .block_index
            .first()
            .map_or(self.footer.index_offset, |first| first.block_offset);
        self.walk_from(0, at)
    }

    /// Every record in the file, building only the bodies `want` asks for.
    ///
    /// **This is the shape the cohort's first pass uses** (spec §6.2), and it is the whole-file
    /// case of [`RecordIter::building_only_where`] — a walk from a coordinate takes a predicate the same
    /// way. A record the predicate declines still arrives, in order, with its head; what it does
    /// not carry is a body.
    ///
    /// ```no_run
    /// # use pop_var_caller::ng::psp::PspReader;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// # let mut psp = PspReader::open(std::path::Path::new("a.psp"))?;
    /// for found in psp.records_where(|head| head.non_reference_reads > 0)? {
    ///     let found = found?;
    ///     if let Some(record) = found.record {
    ///         // only the records the predicate wanted were built
    ///         let _ = record;
    ///     }
    /// }
    /// # Ok(()) }
    /// ```
    pub fn records_where<F>(&mut self, want: F) -> Result<SelectiveRecordIter<'_, F>, PspReadError>
    where
        F: FnMut(&RecordHead) -> bool,
    {
        Ok(self.records()?.building_only_where(want))
    }

    /// Every record from one block onwards, named by its ordinal in
    /// [`block_index`](Self::block_index).
    ///
    /// **The block-level entry point [`records_from`](Self::records_from) is built on**
    /// (spec §6.2), and the one a caller that has already searched the index itself wants.
    pub fn records_from_block(&mut self, block: usize) -> Result<RecordIter<'_>, PspReadError> {
        let Some(entry) = self.block_index.get(block).copied() else {
            return Err(PspReadError::NoSuchBlock {
                path: self.path.clone(),
                ordinal_asked_for: block as u64,
                blocks_in_the_file: self.block_index.len() as u64,
            });
        };
        self.walk_from(block as u64, entry.block_offset)
    }

    /// Every record from the block holding `at` onwards.
    ///
    /// **`records_from` exists because callers think in coordinates** (spec §6.2). The index
    /// turns one into a block with a single binary search over the entries `open` already
    /// holds; no block is searched, because the entries carry only where each block *starts*.
    ///
    /// **⚠ The walk starts at a block's first record, not at `at`.** A reader cannot start
    /// mid-block (spec §1.2), so the records that come back begin at or before the coordinate
    /// asked for — usually well before it, since a block spans 100 kb by default. A caller that
    /// wants only records from `at` onwards drops the ones in front itself, which costs a head
    /// each and no body.
    ///
    /// **⚠ Which block, and it is one earlier than it looks.** The index says where each block
    /// *starts* and nothing about where its records end, and a block's last record may begin on
    /// the very base the next block begins on — a byte ceiling closes a block after a record and
    /// the record after it may start on the same base (`index.rs`). So the block whose first
    /// position equals `at` is **not** necessarily the first one holding a record at `at`, and
    /// the walk enters the block *before* the first block that starts at or after `at`. The cost
    /// is one extra block, and only when `at` falls exactly on a block's first position; the
    /// alternative loses records silently.
    ///
    /// **⚠ And it selects on where records *start*, not on what they span.** A record that
    /// begins in the block before the one chosen and reaches past `at` — a deletion is the case
    /// that does this — is not in the walk. **This is not an overlap query and the format
    /// cannot make it one**: an index entry carries a block's first position and nothing else,
    /// and spec §3.3 removed the only field that could have said how far a block's records
    /// reach. Production's index keeps a `last_pos` for exactly this and ng's dropped it
    /// deliberately. A caller that needs the records covering a coordinate starts a block
    /// earlier and looks at their spans.
    pub fn records_from(&mut self, at: GenomePosition) -> Result<RecordIter<'_>, PspReadError> {
        if self.block_index.is_empty() {
            return self.records();
        }
        // The entries are non-decreasing in `first_position` (`index.rs` refuses an index that
        // is not, at open), so this is the first block that starts at or after `at`.
        let at_or_after = self
            .block_index
            .partition_point(|entry| entry.first_position < at);
        // **And the walk enters the one before it**, because that block's records run up to and
        // including the next block's first position — so it is the earliest block that can hold
        // a record starting at `at`, and the index cannot say whether it does.
        //
        // ⚠ This was `<=` and then one back, which reads as *the last block starting at or
        // before `at`* and is wrong twice over: on a run of blocks sharing one first position it
        // enters the run at its **end**, and even without a run it skips a previous block whose
        // last record starts exactly at `at`. Both lose records with no error — the shape spec
        // §3.3 goal 3 exists against — and the fixture that found it is
        // `a_walk_from_a_position_two_blocks_share_starts_at_the_first_of_them`.
        //
        // **A coordinate in front of every block starts at the first**, rather than being
        // refused: it is what a caller asking for a whole contig's records from position 1
        // writes, and the first block of that contig may well start further in.
        self.records_from_block(at_or_after.saturating_sub(1))
    }

    /// The walk itself: hand the file, bounded at the end of its blocks, to `walk.rs`.
    ///
    /// **The seek and the bound are all this reader contributes.** Everything that inflates a
    /// frame lives in `walk.rs` and is deliberately not reachable from here — see this module's
    /// own `the_opener_cannot_reach_any_block_decoding_code`.
    fn walk_from(
        &mut self,
        first_block: u64,
        block_offset: u64,
    ) -> Result<RecordIter<'_>, PspReadError> {
        // Destructured with no `..`: a field added to the reader has to be considered here.
        let Self {
            path,
            file,
            header,
            footer,
            block_index: _,
            window_budget_bytes: _,
            record_buffer_ceiling_bytes,
        } = self;
        walk::walk_from(
            path,
            file,
            &header.manifest,
            walk::WalkStart {
                blocks_end: footer.index_offset,
                block_offset,
                first_block,
                record_buffer_ceiling_bytes: *record_buffer_ceiling_bytes,
            },
        )
    }

    /// Read and check the fixed tail — [`read_and_check_the_footer`], with this reader's own
    /// handle.
    fn read_footer(path: &Path, file: &mut File, file_bytes: u64) -> Result<Footer, PspReadError> {
        read_and_check_the_footer(path, file, file_bytes)
    }
}

/// Read a psp's fixed tail and apply the checks 48 bytes cannot make about themselves.
///
/// **Shared between `open` and the trailer replacement**, which are the two operations that
/// start from the footer. It was copied into the second and the copy carried the first's
/// reasoning with it, which was already untrue of the copy — the G3 review's finding.
///
/// **A file too short to hold a footer is incomplete, not damaged** — it is what a writer killed
/// before it finished leaves, which is the everyday case.
///
/// The two rules below are the ones that need the file's **length**, which the footer does not
/// carry: that the sections end exactly where the footer begins, and that the index does not
/// start inside the header's magic. Everything the 48 bytes can check about themselves —
/// the index ending exactly where the trailer begins among it — is [`decode_footer`]'s.
pub(super) fn read_and_check_the_footer(
    path: &Path,
    file: &mut File,
    file_bytes: u64,
) -> Result<Footer, PspReadError> {
    let footer_bytes = FOOTER_BYTES as u64;
    let Some(footer_at) = file_bytes.checked_sub(footer_bytes) else {
        return Err(PspReadError::Incomplete {
            path: path.to_path_buf(),
        });
    };
    file.seek(SeekFrom::Start(footer_at))
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "seeking to the footer",
            source,
        })?;
    let mut tail = [0u8; FOOTER_BYTES];
    file.read_exact(&mut tail)
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "reading the footer",
            source,
        })?;

    let footer = decode_footer(&tail).map_err(|refused| match refused {
        // **No tail magic has two readings and they need telling apart**: a psp whose writer
        // never finished, and a file that was never a psp. The head magic is what separates
        // them, and it is read only on this path — the everyday open pays nothing for it.
        super::footer::FooterDecodeError::NotAFooter { .. } => match read_the_head_magic_of(file) {
            Some(found) if found != HEAD_MAGIC => PspReadError::NotAnNgPsp {
                path: path.to_path_buf(),
                found,
                expected: HEAD_MAGIC,
            },
            _ => PspReadError::Incomplete {
                path: path.to_path_buf(),
            },
        },
        // **`reason` names the section and the cause says what was wrong with it.** An
        // earlier version formatted the decoder's own sentence into `reason` as well, so a
        // caller printing the chain saw it twice.
        damaged => PspReadError::damaged_by(path, "the footer does not decode", damaged.into()),
    })?;

    // **What the footer could not check about itself, because it does not know the length.**
    // The trailer must end exactly where the footer begins: a file whose sections stop short
    // has bytes nothing accounts for, and one that runs past has sections overlapping the
    // footer.
    let sections_end = footer
        .trailer_offset
        .checked_add(footer.trailer_bytes)
        .ok_or_else(|| {
            PspReadError::damaged(path, "the trailer's end is past any address".to_string())
        })?;
    if sections_end != footer_at {
        return Err(PspReadError::damaged(
            path,
            format!(
                "the footer says the file's sections end at byte {sections_end}, but the \
                     footer itself begins at byte {footer_at}"
            ),
        ));
    }
    if footer.index_offset < HEAD_MAGIC.len() as u64 {
        return Err(PspReadError::damaged(
            path,
            format!(
                "the footer puts the block index at byte {}, which is inside the header",
                footer.index_offset
            ),
        ));
    }
    Ok(footer)
}

/// The file's first four bytes, for [`read_and_check_the_footer`]'s one failure path.
fn read_the_head_magic_of(file: &mut File) -> Option<[u8; 4]> {
    file.seek(SeekFrom::Start(0)).ok()?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic).ok()?;
    Some(magic)
}

impl PspReader {
    /// Read the index, check it against the footer's checksum, and decode it.
    fn read_index(
        path: &Path,
        file: &mut File,
        footer: &Footer,
        file_bytes: u64,
    ) -> Result<Vec<BlockIndexEntry>, PspReadError> {
        // **What actually bounds this allocation, and it is not the check below.**
        //
        // `decode_footer` has already proved `index_offset + index_bytes == trailer_offset`, and
        // `read_footer` has already proved `trailer_offset + trailer_bytes` is exactly where the
        // footer begins, which is `file_bytes - 48`. Together those give
        // `index_offset + index_bytes <= file_bytes` before this function is entered — so the
        // length can never exceed the file, and the buffer below can never be sized past it.
        //
        // ⚠ **The check that follows therefore cannot fire, and it is kept anyway.** It was
        // written believing it was the bound; it is not, and a reader who assumed it was would be
        // trusting the wrong line. It stays because the bound it restates lives in two other
        // functions with nothing in this signature saying so, and a rule this cheap should not
        // depend on both of them continuing to hold. Replacing its body with `unreachable!()`
        // leaves the suite green, which is how it was found.
        let index_end = footer
            .index_offset
            .checked_add(footer.index_bytes)
            .ok_or_else(|| {
                PspReadError::damaged(
                    path,
                    "the block index's end is past any address".to_string(),
                )
            })?;
        if index_end > file_bytes {
            return Err(PspReadError::damaged(
                path,
                format!(
                    "the footer puts a {}-byte block index at byte {} of a {file_bytes}-byte file",
                    footer.index_bytes, footer.index_offset
                ),
            ));
        }
        file.seek(SeekFrom::Start(footer.index_offset))
            .map_err(|source| PspReadError::Io {
                path: path.to_path_buf(),
                while_doing: "seeking to the block index",
                source,
            })?;
        let mut bytes = vec![0u8; footer.index_bytes as usize];
        file.read_exact(&mut bytes)
            .map_err(|source| PspReadError::Io {
                path: path.to_path_buf(),
                while_doing: "reading the block index",
                source,
            })?;

        // **The checksum before the decode.** The index is the one region of a psp no zstd frame
        // checksum covers and which carries no framing of its own, so this is what says its bytes
        // are the bytes that were written — and a damaged index that still happens to decode
        // would seek a reader to the wrong block without failing.
        let found = index::checksum_index(&bytes);
        if found != footer.index_checksum {
            return Err(PspReadError::damaged(
                path,
                format!(
                    "the block index checksums to {found:#010x}; the footer says {:#010x}",
                    footer.index_checksum
                ),
            ));
        }

        let blocks = index::decode_index(&bytes, footer.n_blocks).map_err(|refused| {
            PspReadError::damaged_by(path, "the block index does not decode", refused.into())
        })?;

        Ok(blocks)
    }

    fn seek_to(&mut self, at: u64, while_doing: &'static str) -> Result<(), PspReadError> {
        self.file
            .seek(SeekFrom::Start(at))
            .map(|_| ())
            .map_err(|source| PspReadError::Io {
                path: self.path.clone(),
                while_doing,
                source,
            })
    }
}

const _: () = assert!(
    DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES < (1u64 << MAX_LOOK_BACK_WINDOW_LOG),
    "a budget at or above the format's widest window could never refuse anything"
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::psp::footer::{FOOTER_MAGIC, encode_footer};
    use crate::ng::psp::header::DEFAULT_LOOK_BACK_WINDOW_LOG;
    use crate::ng::psp::writer::PspWriter;
    use crate::ng::psp::writer::tests_support::{
        a_file, a_finished_psp, a_header, a_record, a_sample, bytes_of, rewrite,
    };

    // -----------------------------------------------------------------
    // What opening gives
    // -----------------------------------------------------------------

    /// Opening gives the header, the block list and the trailer — everything the file says about
    /// itself.
    #[test]
    fn opening_gives_the_header_the_blocks_and_the_trailer() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");

        assert_eq!(psp.header().sample, "SRR7279481");
        assert_eq!(psp.header().contigs.len(), 2);
        assert_eq!(psp.block_index().len() as u64, psp.footer().n_blocks);
        assert!(
            psp.block_index().len() >= 8,
            "two contigs of four grid cells"
        );
        assert_eq!(
            psp.trailer().expect("the trailer reads"),
            b"a per-sample summary"
        );

        // The blocks come back in genomic order, which is what a seek searches on.
        let mut sorted = psp.block_index().to_vec();
        sorted.sort();
        assert_eq!(sorted, psp.block_index());
    }

    /// **Opening decompresses no block**, shown rather than asserted: every byte of the blocks
    /// region is overwritten with rubbish and the file still opens and still says the same things
    /// about itself. Verified in review — making `open` decompress the first block fails this
    /// test and only this one.
    ///
    /// ⚠ **It was called `opening_touches_no_block`, and it does not show that.** An `open` that
    /// read every byte of the blocks region into a buffer without decompressing passes it, and
    /// *touching* is the cost spec §6.2 is about. The stronger statement is
    /// `the_opener_cannot_reach_any_block_decoding_code`, which reads it off the module's
    /// imports; showing that nothing is *read* needs `open` drivable over an arbitrary source,
    /// which it is not.
    #[test]
    fn opening_decompresses_no_block() {
        let (_dir, path) = a_finished_psp();
        let intact = PspReader::open(&path).expect("it opens");
        let (blocks_start, blocks_end) = (
            intact.block_index()[0].block_offset,
            intact.footer().index_offset,
        );
        let expected_blocks = intact.block_index().to_vec();
        drop(intact);

        let mut bytes = bytes_of(&path);
        assert!(blocks_end > blocks_start, "the file holds blocks");
        for byte in &mut bytes[blocks_start as usize..blocks_end as usize] {
            *byte = 0xA5;
        }
        rewrite(&path, &bytes);

        let mut wrecked = PspReader::open(&path).expect("open must not read a block");
        assert_eq!(wrecked.block_index(), expected_blocks);
        assert_eq!(wrecked.header().sample, "SRR7279481");
        assert_eq!(
            wrecked.trailer().expect("still reads"),
            b"a per-sample summary"
        );
    }

    // -----------------------------------------------------------------
    // What is refused
    // -----------------------------------------------------------------

    /// **A file whose writer never finished is refused, not read short.** This is goal 3, and
    /// the file here is not empty — it holds a header and blocks, which is exactly why reading
    /// it short would be so easy and so wrong.
    #[test]
    fn a_file_with_no_footer_is_refused_as_incomplete() {
        let (_dir, path) = a_file();
        {
            let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            // Dropped without `finish`.
        }
        let refused = PspReader::open(&path).expect_err("an unfinished psp must be refused");
        assert!(
            matches!(refused, PspReadError::Incomplete { .. }),
            "got {refused}"
        );
        assert_eq!(
            refused.to_string(),
            format!(
                "{} has no valid footer — the writer did not finish",
                path.display()
            )
        );
    }

    /// A file too short to hold even a footer is incomplete, and does not panic on the
    /// subtraction that finds the footer.
    #[test]
    fn a_file_shorter_than_a_footer_is_refused_without_panicking() {
        for length in [0usize, 1, FOOTER_BYTES - 1] {
            let (_dir, path) = a_file();
            rewrite(&path, &vec![0u8; length]);
            let refused = PspReader::open(&path).expect_err("too short to be a psp");
            assert!(
                matches!(refused, PspReadError::Incomplete { .. }),
                "at {length} bytes, got {refused}"
            );
        }
    }

    /// **A file that was never an ng psp is told apart from one whose writer was killed**, and
    /// the head magic is what separates them. Both lack the tail magic; only one can be rebuilt
    /// by re-running the pileup.
    #[test]
    fn a_foreign_file_is_refused_as_foreign_and_not_as_unfinished() {
        let (_dir, path) = a_file();
        // Production's own psp head magic: the everyday wrong file, since both use `.psp`.
        let mut bytes = b"PSP\n".to_vec();
        bytes.extend_from_slice(&[0u8; 200]);
        rewrite(&path, &bytes);
        let refused = PspReader::open(&path).expect_err("that is not an ng psp");
        match refused {
            PspReadError::NotAnNgPsp {
                found, expected, ..
            } => {
                assert_eq!(found, *b"PSP\n");
                assert_eq!(expected, HEAD_MAGIC);
            }
            other => panic!("expected a foreign file, got {other}"),
        }
    }

    /// A footer whose sections do not end where the footer begins is damage, not incompleteness:
    /// the bytes on disk disagree with each other.
    #[test]
    fn a_footer_that_does_not_account_for_the_file_is_refused_as_damaged() {
        let (_dir, path) = a_finished_psp();
        let mut bytes = bytes_of(&path);
        let tail_at = bytes.len() - FOOTER_BYTES;
        let mut footer =
            crate::ng::psp::footer::decode_footer(&bytes[tail_at..].try_into().unwrap())
                .expect("a footer");
        // Claim a trailer one byte shorter than it is: the sections then stop short of the tail.
        footer.trailer_bytes -= 1;
        bytes[tail_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &bytes);

        let refused = PspReader::open(&path).expect_err("the sections must account for the file");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("sections end at byte"),
            "got {refused}"
        );
    }

    /// **A damaged index is caught by the checksum, before it is decoded.** A flipped bit inside
    /// an offset gives a perfectly well-formed index pointing somewhere else, so nothing but the
    /// checksum can see it.
    #[test]
    fn an_index_that_does_not_match_its_checksum_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut bytes = bytes_of(&path);
        let footer = crate::ng::psp::footer::decode_footer(
            &bytes[bytes.len() - FOOTER_BYTES..].try_into().unwrap(),
        )
        .expect("a footer");
        // One bit, inside the index, in a byte that is part of an offset.
        bytes[footer.index_offset as usize + 4] ^= 0x01;
        rewrite(&path, &bytes);

        let refused = PspReader::open(&path).expect_err("a damaged index must be refused");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("checksums to"),
            "got {refused}"
        );
    }

    /// An index naming a block outside the blocks is refused: the offset is where a reader would
    /// seek, and seeking into the index would decompress the index.
    #[test]
    fn an_index_pointing_outside_the_blocks_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut bytes = bytes_of(&path);
        let tail_at = bytes.len() - FOOTER_BYTES;
        let footer = crate::ng::psp::footer::decode_footer(&bytes[tail_at..].try_into().unwrap())
            .expect("a footer");

        // Rebuild the index with its last block moved into the index's own bytes, and restamp
        // the checksum so the checksum test is not what fires.
        let mut entries = index::decode_index(
            &bytes
                [footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize],
            footer.n_blocks,
        )
        .expect("the index reads");
        // **Exactly at `index_offset`, not past it.** The old fixture used `+ 1`, so relaxing
        // the bound from `>=` to `>` — a block starting precisely where the index does — left
        // every test green.
        let last = entries.len() - 1;
        entries[last].block_offset = footer.index_offset;
        let rebuilt = index::encode_index(&entries);
        assert_eq!(rebuilt.len() as u64, footer.index_bytes, "same width");
        bytes[footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize]
            .copy_from_slice(&rebuilt);
        let mut footer = footer;
        footer.index_checksum = index::checksum_index(&rebuilt);
        bytes[tail_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &bytes);

        let refused = PspReader::open(&path).expect_err("a block outside the blocks");
        assert!(
            refused.to_string().contains("the blocks run from byte"),
            "got {refused}"
        );
    }

    /// **Both sides of "the sections end where the footer begins."**
    ///
    /// ⚠ Only one side was tested: the fixture claimed a trailer one byte *shorter* than it was.
    /// Weakening the rule from `!=` to `<` kept all fifteen tests green — and **62 of the 384
    /// single-bit flips in a real psp's footer stopped being refused and started opening**. One
    /// character removed the property the whole format rests on, and the suite said nothing.
    #[test]
    fn sections_that_stop_short_or_run_past_the_footer_are_both_refused() {
        for (nudge, what) in [(-1i64, "stop short of"), (1, "run past")] {
            let (_dir, path) = a_finished_psp();
            let mut bytes = bytes_of(&path);
            let tail_at = bytes.len() - FOOTER_BYTES;
            let mut footer =
                crate::ng::psp::footer::decode_footer(&bytes[tail_at..].try_into().unwrap())
                    .expect("a footer");
            footer.trailer_bytes = footer.trailer_bytes.wrapping_add(nudge as u64);
            bytes[tail_at..].copy_from_slice(&encode_footer(&footer));
            rewrite(&path, &bytes);

            let refused =
                PspReader::open(&path).expect_err("the sections must account for the file");
            assert!(
                refused.to_string().contains("sections end at byte"),
                "sections that {what} the footer must be refused: {refused}"
            );
        }
    }

    /// A block offset inside the header is refused **at open**, not deferred to a corrupt block
    /// at read time.
    ///
    /// ⚠ The range check bounded offsets above and not below, so an entry at byte 0 opened. The
    /// wrong instruction, and it arrives after a cohort has committed to the sample.
    #[test]
    fn an_index_pointing_into_the_header_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut bytes = bytes_of(&path);
        let tail_at = bytes.len() - FOOTER_BYTES;
        let mut footer =
            crate::ng::psp::footer::decode_footer(&bytes[tail_at..].try_into().unwrap())
                .expect("a footer");
        let mut entries = index::decode_index(
            &bytes
                [footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize],
            footer.n_blocks,
        )
        .expect("the index reads");
        entries[0].block_offset = 0;
        let rebuilt = index::encode_index(&entries);
        assert_eq!(rebuilt.len() as u64, footer.index_bytes, "same width");
        bytes[footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize]
            .copy_from_slice(&rebuilt);
        footer.index_checksum = index::checksum_index(&rebuilt);
        bytes[tail_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &bytes);

        let refused = PspReader::open(&path).expect_err("byte 0 is the header, not a block");
        assert!(
            refused.to_string().contains("the blocks run from byte"),
            "got {refused}"
        );
    }

    /// A footer declaring an index larger than the file is refused **before a buffer for it
    /// exists**, which is the discipline `read_header` states for the header's own length.
    #[test]
    fn an_index_longer_than_the_file_does_not_size_a_buffer() {
        let (_dir, path) = a_finished_psp();
        let mut bytes = bytes_of(&path);
        let tail_at = bytes.len() - FOOTER_BYTES;
        let mut footer =
            crate::ng::psp::footer::decode_footer(&bytes[tail_at..].try_into().unwrap())
                .expect("a footer");
        // Keep the sections abutting so it is this rule that fires, not that one.
        footer.index_bytes = u64::MAX / 2;
        footer.trailer_offset = footer.index_offset.wrapping_add(footer.index_bytes);
        footer.trailer_bytes = 0;
        bytes[tail_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &bytes);

        let refused = PspReader::open(&path).expect_err("that index cannot be in that file");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
    }

    /// The trailer reads the same twice: the seek is not an accident of where the cursor
    /// happened to be.
    ///
    /// ⚠ Every other test calls it once, when the cursor is already in the right place — so a
    /// missing seek would return the footer's bytes on the second call and nothing would notice.
    #[test]
    fn the_trailer_reads_the_same_twice() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("it opens");
        let first = psp.trailer().expect("the trailer reads");
        let second = psp.trailer().expect("and reads again");
        assert_eq!(first, b"a per-sample summary");
        assert_eq!(first, second);
    }

    /// **`open` reaches no block-decoding code at all** — a stronger statement than the
    /// rubbish-bytes test, which shows only that opening does not *decode* a block.
    ///
    /// Read off the module's own imports: `reader.rs` names neither the block module nor the
    /// record module, so anything that inflates a frame has to be reached through `walk`, which
    /// is one named seam a reviewer can hold in view — rather than from anywhere in this file.
    ///
    /// **It is why the record walk is a file of its own.** Milestone G put `records` and its
    /// siblings on this type, and the walk they hand back needs a `BlockStream` — so building it
    /// here would have added the very import this test forbids.
    ///
    /// ⚠ **What survives G1 is the seam, not "no block code is reachable from here"**, which is
    /// what an earlier version of this comment claimed: `walk::walk_from` builds a `BlockStream`,
    /// and this file calls it. The weaker statement is the true one and is the one worth having.
    ///
    /// ⚠ **The braced import forms are on the list too.** `use super::{block::BlockStream, …}`
    /// contains neither `super::block` nor `psp::block`, compiles, and passed this test verbatim
    /// until `block::` and `record::` were added — a form `rustfmt` produces on its own.
    #[test]
    fn the_opener_cannot_reach_any_block_decoding_code() {
        let source = include_str!("reader.rs");
        let before_tests = source
            .split("#[cfg(test)]")
            .next()
            .expect("the file has a non-test half");
        for forbidden in [
            "psp::block",
            "psp::record",
            "super::block",
            "super::record",
            // The braced forms, which none of the four above contain.
            "block::",
            "record::",
        ] {
            assert!(
                !before_tests.contains(forbidden),
                "opening must not be able to reach {forbidden}: a block is decompressed only \
                 by a walk, never by an open"
            );
        }
    }

    /// **A window wider than this reader budgeted for is its own refusal, and the message names
    /// both numbers** — because the fix is a setting rather than a rebuilt file (spec §4.2).
    #[test]
    fn a_window_wider_than_the_budget_is_refused_by_name() {
        let (_dir, path) = a_file();
        let mut header = a_header(1_000);
        header.manifest.look_back_window_log = 22; // 4 MiB
        let writer = PspWriter::create(&path, header).expect("a wide window is writable");
        let _ = writer.finish(&[]).expect("it finishes");

        let refused = PspReader::open(&path).expect_err("wider than the default budget");
        match refused {
            PspReadError::WindowTooLarge {
                needed_bytes,
                allowed_bytes,
                ..
            } => {
                assert_eq!(needed_bytes, 1 << 22);
                assert_eq!(allowed_bytes, DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES);
            }
            other => panic!("expected a window refusal, got {other}"),
        }

        // And raising the budget opens it, which is what makes the refusal a setting.
        let raised = PspReader::open_with_a_look_back_window_budget(&path, 1 << 22)
            .expect("a raised budget opens it");
        assert_eq!(raised.header().manifest.look_back_window_log, 22);
    }

    /// The window this build writes opens under the default budget with room to spare.
    #[test]
    fn the_window_this_build_writes_is_well_inside_the_default_budget() {
        let written = 1u64 << DEFAULT_LOOK_BACK_WINDOW_LOG;
        assert_eq!(written, 32 * 1024);
        assert_eq!(DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES, 256 * 1024);
        assert_eq!(
            DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES / written,
            8,
            "the budget is eight times the window this build produces"
        );
    }

    /// A psp with no records opens, and says so: no blocks, and a trailer that is still there.
    #[test]
    fn a_psp_with_no_records_opens() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let _ = writer.finish(b"nothing was found").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("an empty psp is a psp");
        assert!(psp.block_index().is_empty());
        assert_eq!(psp.footer().n_blocks, 0);
        assert_eq!(psp.trailer().expect("it reads"), b"nothing was found");
    }

    /// An empty trailer reads back as no bytes rather than as an error.
    #[test]
    fn an_empty_trailer_reads_back_empty() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        writer.push(&a_record(0, 1, 1)).expect("one record");
        let _ = writer.finish(&[]).expect("it finishes");

        let mut psp = PspReader::open(&path).expect("it opens");
        assert!(psp.trailer().expect("it reads").is_empty());
    }

    /// A file that does not exist names itself and what was being done to it.
    #[test]
    fn a_missing_file_names_itself() {
        let (dir, path) = a_file();
        drop(dir);
        let refused = PspReader::open(&path).expect_err("it is not there");
        assert!(matches!(refused, PspReadError::Io { .. }), "got {refused}");
        assert!(
            refused.to_string().contains("opening the file"),
            "got {refused}"
        );
    }

    /// **Every truncation of a finished psp is refused, and none panics.** A file cut anywhere
    /// is either incomplete or damaged — never quietly readable.
    ///
    /// **Every byte, not every sixteenth** (Milestone H2). A killed writer stops at whatever byte
    /// the kernel had taken, and one cut in sixteen leaves fifteen of every sixteen stopping
    /// points unvisited — including most of the two- and four-byte fields a section boundary is
    /// made of. Measured: **3,742 cuts on this fixture, 0.19 s**, against 234 cuts before. The
    /// count is asserted below rather than left to this sentence.
    #[test]
    fn every_truncation_of_a_finished_psp_is_refused_without_panicking() {
        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let mut incomplete = 0;
        let mut damaged = 0;
        for cut in 0..whole.len() {
            rewrite(&path, &whole[..cut]);
            match PspReader::open(&path) {
                Err(PspReadError::Incomplete { .. }) => incomplete += 1,
                Err(PspReadError::Damaged { .. }) => damaged += 1,
                Err(PspReadError::NotAnNgPsp { .. }) => incomplete += 1,
                Err(other) => panic!("a cut at {cut} gave {other}"),
                Ok(_) => panic!("a cut at {cut} opened"),
            }
        }
        assert!(
            incomplete > 0 && damaged == 0,
            "a truncated psp has lost its footer, so every cut is incomplete: \
             {incomplete} incomplete, {damaged} damaged"
        );
        assert_eq!(
            incomplete,
            whole.len(),
            "every cut has to be accounted for, or the sweep is proving less than it counts"
        );
    }

    /// The tail magic alone is not enough: a file that ends with it but holds nothing else is
    /// refused rather than opened.
    #[test]
    fn a_file_that_is_only_a_footer_magic_is_refused() {
        let (_dir, path) = a_file();
        let mut bytes = vec![0u8; FOOTER_BYTES];
        bytes[FOOTER_BYTES - 4..].copy_from_slice(&FOOTER_MAGIC);
        rewrite(&path, &bytes);
        let refused = PspReader::open(&path).expect_err("that describes no file");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
    }
}
