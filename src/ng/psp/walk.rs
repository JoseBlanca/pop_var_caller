//! Walking a psp's records: seek to a block, and stream forward to where the blocks end.
//!
//! **A file away from [`super::reader`], and the split is the guarantee.** Opening a psp must
//! touch no block (spec §6.2) — that is the cost a cohort pays per open sample, multiplied by
//! the cohort size — and the way that is held is structural: `reader.rs` names nothing that can
//! inflate a frame, and one of its own tests reads its imports to say so. The walk needs
//! [`BlockStream`], so it lives here instead of weakening that.
//!
//! **What this file adds to `BlockStream` is two things and no decoding.** One: the source is
//! bounded at the end of the blocks, so a stream that would otherwise read the index, the
//! trailer and the footer as further blocks stops cleanly instead. Two: a block-read failure is
//! put in the class whose instruction fits it — the file is damaged, the reader is too old, or
//! a limit needs raising (spec §7) — and named with the block ordinal a caller's own
//! [`super::PspReader::block_index`] uses.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use super::PspReadError;
use super::block::{BlockReadError, BlockStream, ROLLING_BYTES, StreamedRecord};
use super::chain_ids::LiveSet;
use super::header::Manifest;
use super::record::RecordHead;

/// The default ceiling on how much of a walk's rolling buffer **one record** may have.
///
/// **Named here rather than taken from `block.rs` by the reader**, because `reader.rs` may not
/// name the module that owns the buffer — see this file's own header.
pub(super) const DEFAULT_RECORD_BUFFER_CEILING_BYTES: usize =
    super::block::ROLLING_BUFFER_CEILING_BYTES;

/// Whether a walk could honour this ceiling on one record's share of its rolling buffer.
///
/// **The rule lives beside the buffer it is about.** A ceiling at or under the buffer itself
/// would turn away records the buffer already holds without growing at all, so it is refused —
/// and refused where the setting is made rather than at the record that would have tripped it.
pub(super) fn check_a_record_buffer_ceiling(
    path: &Path,
    ceiling: usize,
) -> Result<(), PspReadError> {
    if ceiling <= ROLLING_BYTES {
        return Err(PspReadError::RecordBufferCeilingTooSmall {
            path: path.to_path_buf(),
            ceiling,
            buffer_bytes: ROLLING_BYTES,
        });
    }
    Ok(())
}

/// Where a walk begins, where it ends, and which block it begins at.
///
/// **Named fields rather than three `u64` parameters.** `blocks_end` and `block_offset` are both
/// byte offsets into the same file, and swapping them bounds the walk at zero bytes — a walk
/// that yields nothing, which no error and no failed seek would report. A psp read as covering
/// nothing is the failure spec §3.3 goal 3 is written against, and an argument order is no way
/// to arrive at it.
pub(super) struct WalkStart {
    /// The file's index offset: the first byte that is no longer a block.
    pub blocks_end: u64,
    /// Where the block to start at begins.
    pub block_offset: u64,
    /// The ordinal of the block at `block_offset`, carried only so that a failure can name the
    /// block it happened in.
    pub first_block: u64,
    /// How much of the rolling buffer one record may have.
    pub record_buffer_ceiling_bytes: usize,
}

/// Seek to a block and hand back a walk that ends where the blocks end.
///
/// **The bound is what stops the reader walking into the index.** A [`BlockStream`] reads
/// blocks until its source runs out, and a psp does not end with its blocks. Handed the whole
/// file, it would take the index's first four bytes for a block length and try to inflate what
/// followed. Bounding the source turns the blocks' end into an end of file, which is the one
/// thing a stream already stops cleanly on.
pub(super) fn walk_from<'a>(
    path: &'a Path,
    file: &'a mut File,
    manifest: &Manifest,
    start: WalkStart,
) -> Result<RecordIter<'a>, PspReadError> {
    let WalkStart {
        blocks_end,
        block_offset,
        first_block,
        record_buffer_ceiling_bytes,
    } = start;
    // `open` proved every block offset lies in the blocks, and the empty-file case is given the
    // index's own offset, so this never saturates; it is written this way because an underflow
    // here would be a sixteen-exabyte read rather than an error.
    let blocks_bytes = blocks_end.saturating_sub(block_offset);
    // ⚠ **This arm has no test, and the H2 review could not find a fixture that reaches it.**
    // `Seek` on a regular file fails only for an offset `SeekFrom::Start(u64)` cannot express,
    // and every offset arriving here has already been bounded inside the file by `open`. The
    // write-only descriptor that gives the *read* below a genuine `EBADF` does not help: seeking
    // a write-only descriptor succeeds, which is exactly why that trick works for the read. It is
    // written out rather than left to `?` on a bare `io::Error` so that the class is right the
    // day something does reach it — but do not read it as covered.
    file.seek(SeekFrom::Start(block_offset))
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "seeking to a block",
            source,
        })?;
    let stream = BlockStream::new(file.take(blocks_bytes), manifest).map_err(|refused| {
        // **The two variants `BlockStream::new` can return, named rather than caught by `_`.**
        // Everything else `BlockReadError` holds comes out of a walk, not out of building one.
        match refused {
            // **Upgrade the reader, not rebuild the file**: the manifest declares a field
            // encoding this build does not know.
            BlockReadError::UnsupportedRecordLayout { .. } => {
                PspReadError::UnsupportedRecordEncoding {
                    path: path.to_path_buf(),
                    source: refused,
                }
            }
            // A look-back window outside the format's range, which `open` has already refused by
            // parsing the header (`header.rs` bounds it against `MIN..=MAX` at parse) — so this
            // arm is unreachable through a `PspReader`. It is written out rather than left to
            // `unreachable!` because the day the header stops checking is not the day to learn
            // that a reader panics on a bad manifest.
            BlockReadError::WindowLogOutOfRange {
                look_back_window_log,
            } => PspReadError::damaged(
                path,
                format!(
                    "the file declares a look-back window of 2^{look_back_window_log}, which \
                     is outside the format's range"
                ),
            ),
            other => PspReadError::damaged(
                path,
                format!("the file's manifest cannot drive a reader: {other}"),
            ),
        }
    })?;
    let stream = stream
        .with_a_buffer_ceiling(record_buffer_ceiling_bytes)
        .map_err(|_| PspReadError::RecordBufferCeilingTooSmall {
            path: path.to_path_buf(),
            ceiling: record_buffer_ceiling_bytes,
            buffer_bytes: ROLLING_BYTES,
        })?;
    Ok(RecordIter {
        path,
        stream,
        first_block,
    })
}

/// A walk over a psp's records, from some block to the end of the blocks.
///
/// **Lazy, and it retains no record it has handed over** (arch §4.1). What it holds is the
/// compressed read buffer, the rolling decompressed buffer, the decompressor's state and the
/// record being built — and none of those is a function of the block size, which is goal 1 and
/// the reason the format exists.
///
/// **It borrows the reader that made it**, so a psp has one walk at a time and the file's cursor
/// belongs to that walk while it lives. The reader's own seeks — the trailer's — cannot then be
/// interleaved with a walk's, and no second file handle is opened per sample.
///
/// **An iterator that fails yields `Err` once and then `None`**, and never a half-built record
/// (spec §6.7).
#[derive(Debug)]
#[must_use = "a walk that is not iterated reads nothing"]
pub struct RecordIter<'a> {
    /// **Borrowed rather than cloned**: it is only read when something fails, and a cohort
    /// makes one of these per open sample.
    path: &'a Path,
    stream: BlockStream<std::io::Take<&'a mut File>>,
    /// The ordinal of the block this walk started at, so a failure can name the block it is in
    /// the way a caller's own index names it.
    first_block: u64,
}

impl<'a> RecordIter<'a> {
    /// Which reads are live at the record last handed back.
    ///
    /// **The chain ids the residual observation is derived from** (spec psp_chain_id_encoding
    /// §5): the ids a record does not list are the live set minus the ones it does, so a caller
    /// reconstructing them needs this beside the record.
    pub fn live_reads(&self) -> &LiveSet {
        self.stream.live_reads()
    }

    /// How many blocks this walk has opened, the one it is inside included.
    ///
    /// **Not how many it has finished**, so a walk that has handed back one record of the first
    /// block already answers 1. The same name as [`BlockStream::blocks_begun`], which is what it
    /// forwards: one number, one name.
    pub fn blocks_begun(&self) -> u64 {
        self.stream.blocks_begun()
    }

    /// The ordinal of the block the record last handed back came from, into
    /// [`super::PspReader::block_index`].
    ///
    /// **Nothing else names it.** Two blocks may share a first position (`index.rs`), so the
    /// coordinate a record carries does not identify the block it was in; and a walk started by
    /// [`super::PspReader::records_from`] never learns which block the search chose. This is
    /// that ordinal, and it is why the walk carries the one it began at.
    ///
    /// Before the first record it answers the ordinal the walk began at.
    pub fn current_block(&self) -> u64 {
        self.first_block + self.stream.blocks_begun().saturating_sub(1)
    }

    /// The same walk, building only the records `want` asks for.
    ///
    /// **Here rather than only on [`super::PspReader`]** so that every entry point gets the skip:
    /// `records_from(at)?.building_only_where(…)` is what a cohort reading one region of every sample
    /// writes, and spec §6.2's `records_where` is the whole-file case of it.
    pub fn building_only_where<F>(self, want: F) -> SelectiveRecordIter<'a, F>
    where
        F: FnMut(&RecordHead) -> bool,
    {
        SelectiveRecordIter { walk: self, want }
    }

    /// One record out of the stream, with its failure put in the class whose instruction fits.
    ///
    /// **The one place a record is stepped**, so [`Iterator`] and [`SelectiveRecordIter`] cannot drift:
    /// the only difference between a full walk and a selective one is the predicate handed here.
    fn step(
        &mut self,
        want: impl FnMut(&RecordHead) -> bool,
    ) -> Option<Result<StreamedRecord, PspReadError>> {
        match self.stream.next_record_where(want) {
            Some(Ok(record)) => Some(Ok(record)),
            // **Classified after the stream has been asked**, so the block count names the block
            // the failure happened in rather than the one before it.
            Some(Err(failed)) => Some(Err(refuse(
                self.path,
                self.first_block,
                self.stream.blocks_begun(),
                failed,
            ))),
            None => None,
        }
    }
}

/// Put a walk's refusal in the class whose instruction fits it, and name the block.
///
/// **Three instructions, not one** (spec §7): the file is damaged, the reader is too old, or a
/// limit needs raising. Folding them together is what [`PspReadError::CorruptBlock`]'s own doc
/// warns against, because *rebuild the file* is wrong advice for two of the three.
///
/// **Exhaustive, with no `_` arm, and that is the point.** A variant added to [`BlockReadError`]
/// tomorrow would otherwise become *the file is corrupt* silently, with no compile error at the
/// place that has to choose. It is a free function rather than a method so that G2's selective
/// walk can share it whatever shape it takes.
///
/// `first_block` is the ordinal the walk began at and `blocks_begun` how many it has opened
/// since, the one it is inside included.
fn refuse(
    path: &Path,
    first_block: u64,
    blocks_begun: u64,
    refused: BlockReadError,
) -> PspReadError {
    let path = path.to_path_buf();
    // The block the walk was inside: one per block begun, less the one not yet finished.
    let inside = first_block + blocks_begun.saturating_sub(1);
    // The four bytes that introduce a block are read before the block is counted, so a fault in
    // them is about a block that never began — one past the last that did, and therefore one
    // past the last entry the index has.
    let never_began = first_block + blocks_begun;
    match refused {
        BlockReadError::RecordLargerThanTheReaderAllows { allowed_bytes, .. } => {
            PspReadError::RecordLargerThanTheReaderAllows {
                path,
                block: inside,
                allowed_bytes,
                source: refused,
            }
        }
        BlockReadError::UnsupportedRecord { .. }
        | BlockReadError::UnsupportedRecordLayout { .. } => {
            PspReadError::UnsupportedRecordEncoding {
                path,
                source: refused,
            }
        }
        // A ceiling under the buffer is the caller's mistake, not the file's, and it has the
        // same instruction as `NoSuchBlock`. Unreachable from a `PspReader`, which refuses such
        // a ceiling where it is set — but *rebuild the file* would be the wrong answer here.
        BlockReadError::BufferCeilingUnderTheBuffer {
            ceiling,
            buffer_bytes,
        } => PspReadError::RecordBufferCeilingTooSmall {
            path,
            ceiling,
            buffer_bytes,
        },
        BlockReadError::Io {
            while_doing,
            source,
        } => PspReadError::Io {
            path,
            while_doing,
            source,
        },
        // The file ended inside the length in front of a block, so the block it names is the one
        // that never began.
        BlockReadError::FileEndsInsideABlockLength { .. } => PspReadError::CorruptBlock {
            path,
            block: never_began,
            source: refused,
        },
        // Everything else is the file disagreeing with itself *inside* a block: a frame that
        // will not inflate, a record running past its block, a block holding more than it
        // declared.
        BlockReadError::WindowLogOutOfRange { .. }
        | BlockReadError::DamagedRecord { .. }
        | BlockReadError::RecordRunsPastItsBlock { .. }
        | BlockReadError::BlockHeadRunsPastItsBlock { .. }
        | BlockReadError::DamagedBlockHead { .. }
        | BlockReadError::BlockHoldsMoreThanItDeclared { .. }
        | BlockReadError::BlockDeclaresMoreBytesThanItsFrame { .. }
        | BlockReadError::BlockFrameDidNotEndWithItsRecords
        | BlockReadError::FileEndsInsideABlock { .. }
        | BlockReadError::Zstd { .. } => PspReadError::CorruptBlock {
            path,
            block: inside,
            source: refused,
        },
    }
}

impl Iterator for RecordIter<'_> {
    type Item = Result<StreamedRecord, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.step(|_| true)
    }
}

/// A walk that builds only the records a predicate asks for.
///
/// **The skip is the reader's decision rather than a separate call** (spec §6.2). Every record
/// opens with the head of spec §4.3, so a walk hands the caller the head and the caller says
/// whether it wants the body: a record the predicate declines costs its head and a pointer
/// advance, and comes back with [`StreamedRecord::record`] as `None`.
///
/// **⚠ It is not a filter, and reading it as one is the mistake to avoid.** Every record of every
/// block still arrives, in order; what the predicate decides is whether the *body* was built. A
/// caller that wants only the records it kept drops the rest, and the skip has already saved it
/// the work of building them.
///
/// **The bytes still have to arrive either way.** Skipping saves building the record, not
/// decompressing it: a block comes out of zstd sequentially and there is nothing to seek past.
///
/// **⚠ A declined body is never checked, so this is a weaker reader of damage than a full
/// walk.** The two agreements between a record's head and its body — the declared body length
/// against the bytes the body used, and the head's non-reference read count against the body's —
/// are made in [`super::record::decode_the_body_of`], which a declined record never reaches.
/// **A walk that came back without an `Err` says the file's framing held, not that the records
/// it skipped are sound**, and a caller re-walking to build them may meet damage the first pass
/// did not report. Measured on a three-record block of 102 payload bytes, every byte flipped in
/// turn: a full walk refuses 93 of them, and **a walk declining every body accepts 72 of those
/// 93** — about three in four. `a_declining_walk_accepts_damage_a_full_walk_refuses` is that
/// measurement.
///
/// **What the skip is worth, on this reader.** Keeping one record in a hundred: **3.038× on
/// tomato at 10.3 reads a record and 2.869× on HG002 at 280.0** — depth costs about 5 % of it,
/// not the collapse that was feared when the chain ids' live-set changes joined the head
/// (`reports/implementations/ng_psp_h5_2026-08-30.md`). Those stores were converted from a
/// production `.psp`, whose chain-id column is a fraction of ng's, so they read high.
///
/// **Re-taken 2026-09-04 on a store ng wrote itself** — 8,105,483 loci over tomato SRR7279481's
/// whole genome at 9.7 reads a record, with ng's own chain ids and the head this build writes:
/// **2.930×**, stepping over 99.0 % of the body bytes. It is 3.6 % below the converted store's
/// reading at a comparable depth, which is the direction a fuller head predicts — though the two
/// corpora are not identical, so the gap is not cleanly the chain ids' doing.
#[must_use = "a walk that is not iterated reads nothing"]
pub struct SelectiveRecordIter<'a, F> {
    walk: RecordIter<'a>,
    want: F,
}

/// **Written out rather than derived, because a derived one needs `F: Debug` and a closure is
/// not** — so the derive was inert for every predicate a caller can actually pass, and the type
/// looked printable in rustdoc while being unprintable at every call site. The predicate is the
/// one field a reader cannot be shown, so it is named and skipped.
impl<F> std::fmt::Debug for SelectiveRecordIter<'_, F> {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Destructured with no `..`: a field added here has to be considered.
        let Self { walk, want: _ } = self;
        out.debug_struct("SelectiveRecordIter")
            .field("walk", walk)
            .field("want", &"<predicate>")
            .finish()
    }
}

impl<'a, F> SelectiveRecordIter<'a, F> {
    /// The walk underneath, for everything a full walk can be asked.
    ///
    /// **An accessor rather than a forwarded method each.** The three below are the questions a
    /// selective walk is asked often enough to be worth spelling; anything added to
    /// [`RecordIter`] later — Milestone H4 and H5 are both about numbers a walk could be asked
    /// for — is reachable through this without also having to be added here.
    pub fn walk(&self) -> &RecordIter<'a> {
        &self.walk
    }

    /// Which reads are live at the record last handed back — see [`RecordIter::live_reads`].
    ///
    /// **Exact after a declined record too**, which is the whole point of putting the chain-id
    /// changes in the head rather than in the body: a caller walking with a predicate that
    /// declines most records still knows which reads are live at the ones it takes.
    pub fn live_reads(&self) -> &LiveSet {
        self.walk.live_reads()
    }

    /// How many blocks this walk has opened — see [`RecordIter::blocks_begun`].
    pub fn blocks_begun(&self) -> u64 {
        self.walk.blocks_begun()
    }

    /// The ordinal of the block the record last handed back came from — see
    /// [`RecordIter::current_block`].
    pub fn current_block(&self) -> u64 {
        self.walk.current_block()
    }
}

impl<F: FnMut(&RecordHead) -> bool> Iterator for SelectiveRecordIter<'_, F> {
    type Item = Result<StreamedRecord, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        // `&mut F` is itself `FnMut`, so the predicate is borrowed for the step rather than
        // moved: a walk keeps its own across every record.
        self.walk.step(&mut self.want)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::psp::block::{
        COMPRESSED_BLOCK_LENGTH_BYTES, ROLLING_BUFFER_CEILING_BYTES, ROLLING_BYTES,
    };
    use crate::ng::psp::footer::{FOOTER_BYTES, encode_footer};
    use crate::ng::psp::writer::PspWriter;
    use crate::ng::psp::writer::tests_support::{
        a_file, a_finished_psp, a_header, a_record, a_sample, bytes_of, footer_of, rewrite,
        wreck_the_block,
    };
    use crate::ng::psp::{DamageFound, PspReadError, PspReader};
    use crate::ng::types::{ContigId, GenomePosition, Position};
    // -----------------------------------------------------------------

    // -----------------------------------------------------------------
    // A read that fails, rather than a file that ends
    // -----------------------------------------------------------------

    /// The kernel's *bad file descriptor* refusal, which is what a read on a write-only
    /// descriptor gives. Named so the number is not one of two bare `9`s in this milestone — the
    /// other is `SIGKILL`, in `writer.rs`.
    const EBADF: i32 = 9;

    /// **A `read(2)` that fails part-way through a block reaches the caller as
    /// [`PspReadError::Io`], naming the file and what was being done to it.**
    ///
    /// This is the one branch of [`refuse`] — the function that puts a block-level failure in
    /// the class whose instruction fits it — that no test reached before this step, and the
    /// reason is that the obvious way to break a walk does not break it in this way:
    /// **truncating a file under an open reader gives an end of file, not an error**, and a
    /// stream that runs out of source is a *truncated file*, which is a different class with a
    /// different instruction.
    ///
    /// The failure is produced here without `unsafe`, without closing a descriptor out from
    /// under a `File`, and without a platform trick: **the walk is handed a descriptor opened
    /// write-only on a real psp.** Seeking succeeds, and every `read(2)` on it fails with
    /// `EBADF`, the kernel's *bad file descriptor* refusal — deterministic, and the same on
    /// Linux and macOS.
    ///
    /// **What it holds is that the class survives the trip.** A read failure must not arrive as
    /// `CorruptBlock` — *rebuild the file* is the wrong instruction for a disk that went away,
    /// and the file may be perfectly sound.
    #[test]
    fn a_read_that_fails_part_way_through_a_block_is_an_io_error_and_not_damage() {
        let (_dir, path) = a_finished_psp();
        let footer = footer_of(&bytes_of(&path));
        // The real reader, only to learn where the blocks are and to prove this file is sound.
        let mut sound_reader = PspReader::open(&path).expect("the fixture opens");
        assert!(
            sound_reader.records().expect("a walk").next().is_some(),
            "the fixture has to hold a record, or a failing read has nothing to fail during"
        );
        let first_block_starts_at = sound_reader
            .block_index()
            .first()
            .expect("a block")
            .block_offset;

        // **Write-only, and deliberately not truncating**: the bytes on disk stay a valid psp,
        // so the only thing wrong is that this descriptor cannot be read from.
        let mut write_only_handle = std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .expect("a psp can be opened for writing");
        // **The manifest comes from the file being walked**, not from a fixture builder that
        // happens to agree with it: a walk configured from one file and pointed at another is a
        // wrong test that still fails plausibly.
        let manifest = crate::ng::psp::read_header(&path)
            .expect("the fixture's header reads")
            .manifest;
        let start = WalkStart {
            blocks_end: footer.index_offset,
            block_offset: first_block_starts_at,
            first_block: 0,
            record_buffer_ceiling_bytes: DEFAULT_RECORD_BUFFER_CEILING_BYTES,
        };
        let mut walk = walk_from(&path, &mut write_only_handle, &manifest, start)
            .expect("building a walk reads nothing, so it succeeds even here");
        let refused = walk
            .next()
            .expect("the walk reports the failure rather than ending")
            .expect_err("a descriptor that cannot be read from cannot yield a record");
        match refused {
            PspReadError::Io {
                path: named_path,
                while_doing,
                source,
            } => {
                assert_eq!(
                    named_path, path,
                    "the refusal names the file it was reading"
                );
                assert_eq!(
                    while_doing, "reading a block's compressed bytes",
                    "the refusal says what was being done to the file, in the words `block.rs` \
                     sets when a source read fails — pinned, because an empty phrase would pass \
                     a check that it is merely non-empty"
                );
                assert_eq!(
                    source.raw_os_error(),
                    Some(EBADF),
                    "the cause is the kernel's own refusal, not a manufactured error: {source}"
                );
            }
            other => panic!(
                "a read failure must not arrive as damage — the file is sound and the \
                 instruction is not `rebuild it`; got {other}"
            ),
        }
    }

    /// The contrast the test above rests on: **the same file cut short gives an end of file,
    /// which is damage, not an `Io` failure.** Without this pair, `Io` and `CorruptBlock` could
    /// be swapped and both tests would still pass.
    #[test]
    fn a_file_that_ends_inside_a_block_is_damage_and_not_an_io_error() {
        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let footer = footer_of(&whole);
        let first_block_starts_at = PspReader::open(&path)
            .expect("the fixture opens")
            .block_index()
            .first()
            .expect("a block")
            .block_offset;
        let manifest = crate::ng::psp::read_header(&path)
            .expect("the fixture's header reads")
            .manifest;
        // **One byte off the end of the last block**, so the cut is certainly *inside* a block.
        // A cut that happens to land on a block boundary ends the walk cleanly and correctly —
        // the stream read whole blocks and then ran out — which is not the case this contrast
        // is about, and is how the first version of this test failed.
        let cut = usize::try_from(
            footer
                .index_offset
                .checked_sub(1)
                .expect("the blocks do not end at byte zero"),
        )
        .expect("the fixture is small enough to index");
        assert!(
            cut > usize::try_from(first_block_starts_at).expect("the first block starts early"),
            "the cut has to land after the first block begins"
        );
        rewrite(&path, &whole[..cut]);
        let mut file = std::fs::File::open(&path).expect("the cut file opens");
        let start = WalkStart {
            blocks_end: footer.index_offset,
            block_offset: first_block_starts_at,
            first_block: 0,
            record_buffer_ceiling_bytes: DEFAULT_RECORD_BUFFER_CEILING_BYTES,
        };
        let refused = walk_from(&path, &mut file, &manifest, start)
            .expect("a walk")
            .find_map(Result::err)
            .expect("a walk over a file that ends inside a block must not end cleanly");
        assert!(
            matches!(refused, PspReadError::CorruptBlock { .. }),
            "a source that runs out is a truncated file, not a failed read; got {refused}"
        );
    }

    /// Every record the sample was written from comes back, in the order it was pushed, equal
    /// field for field.
    ///
    /// **The expectation is the fixture, not the reader.** `a_sample` is the list of records
    /// `a_finished_psp` pushes, so this compares the walk against what was written rather than
    /// against what the walk itself decided.
    #[test]
    fn every_record_written_comes_back_in_order() {
        let (_dir, path) = a_finished_psp();
        let written = a_sample();
        assert_eq!(
            written.len(),
            40,
            "two contigs of four grid cells of five records"
        );

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let read: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .map(|found| found.expect("a finished psp walks"))
            .collect();

        assert_eq!(read.len(), written.len());
        for (found, expected) in read.iter().zip(&written) {
            assert_eq!(
                found.record.as_ref(),
                Some(expected),
                "the record at {:?}",
                expected.region
            );
        }
    }

    /// **A walk ends where the blocks end, and reads none of what follows them.**
    ///
    /// A psp does not end with its blocks: the index, the trailer and the footer come after,
    /// and a reader handed the whole file takes the index's first four bytes for a block
    /// length. So the walk must end cleanly — `None`, not `Err` — having opened exactly as many
    /// blocks as the index names.
    #[test]
    fn a_walk_ends_at_the_last_block_and_never_reads_the_index() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let blocks = psp.block_index().len() as u64;
        assert_eq!(blocks, 8, "two contigs of four grid cells");

        let mut walk = psp.records().expect("the walk starts");
        let mut records = 0;
        for found in walk.by_ref() {
            let _ = found.expect("nothing past the blocks is read as one");
            records += 1;
        }
        assert_eq!(records, 40);
        assert_eq!(
            walk.blocks_begun(),
            blocks,
            "the walk opened every block the index names, and nothing else"
        );
    }

    /// **`records_from` on a block's own first position enters the block *before* it**, and
    /// holds every record of the block asked about.
    ///
    /// The block before it is entered because the index says where blocks start and nothing
    /// about where their records end: that block's last record may begin on the very base this
    /// one begins on (`index.rs`), and no entry can say whether it does. Checked for every block
    /// in the file, not for one.
    #[test]
    fn records_from_a_blocks_first_position_enters_the_block_before_it() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();
        assert_eq!(entries.len(), 8);

        for (ordinal, entry) in entries.iter().enumerate() {
            let expected = ordinal.saturating_sub(1);
            let mut walk = psp.records_from(entry.first_position).expect("it starts");
            let first = walk
                .next()
                .expect("a block holds at least one record")
                .expect("a finished psp walks");
            assert_eq!(
                walk.current_block(),
                expected as u64,
                "the walk from block {ordinal}'s first position began in block {}",
                walk.current_block()
            );
            let region = first.record.expect("the body was built").region;
            assert_eq!(
                (region.contig, region.start),
                (
                    entries[expected].first_position.contig,
                    entries[expected].first_position.position
                ),
                "and at that block's first record"
            );
        }
    }

    /// **Two blocks may share a first position, and a walk from it must enter the first of
    /// them.** A byte ceiling closes a block after a record, and the record after it may begin
    /// on the same base (`index.rs` §"Two blocks may share a first position"); asking for that
    /// base must not skip the earlier block.
    ///
    /// ⚠ This is the fixture that found the G1 review's Blocker. The search was *the last block
    /// starting at or before `at`*, which enters a run of equal positions at its **end**: the
    /// walk came back with two records where the file holds three, with no error.
    #[test]
    fn a_walk_from_a_position_two_blocks_share_starts_at_the_first_of_them() {
        let (_dir, path) = a_file();
        let mut header = a_header(1_000_000);
        // One block a record, which is what puts two entries on one base.
        header.manifest.block_byte_ceiling = Some(1);
        let mut writer = PspWriter::create(&path, header).expect("a header");
        for record in [
            a_record(0, 500, 1),
            a_record(0, 500, 10),
            a_record(0, 700, 1),
        ] {
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();
        assert_eq!(entries.len(), 3, "one block a record");
        assert_eq!(
            entries[0].first_position, entries[1].first_position,
            "the fixture must put two blocks on the same base"
        );

        let asked = GenomePosition {
            contig: ContigId(0),
            position: Position(500),
        };
        let starts: Vec<_> = psp
            .records_from(asked)
            .expect("the walk starts")
            .map(|found| {
                found
                    .expect("a finished psp walks")
                    .record
                    .expect("a body")
                    .region
                    .start
            })
            .collect();
        assert_eq!(
            starts,
            vec![Position(500), Position(500), Position(700)],
            "both records at 500 are in a walk from 500"
        );
    }

    /// **A coordinate past every block gives the last block's records, not an empty walk.** It
    /// follows from the contract — the walk starts at a block's first record, not at the
    /// coordinate — and a caller scanning forward from beyond a sample's coverage gets records
    /// in front of it rather than nothing.
    #[test]
    fn a_coordinate_past_every_block_gives_the_last_blocks_records() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();
        let last = entries.last().expect("the fixture has blocks");

        let beyond = GenomePosition {
            contig: last.first_position.contig,
            position: Position(last.first_position.position.get() + 1_000_000),
        };
        let mut walk = psp.records_from(beyond).expect("the walk starts");
        let first = walk
            .next()
            .expect("the last block's records")
            .expect("a finished psp walks");
        assert_eq!(
            walk.current_block(),
            entries.len() as u64 - 1,
            "the walk began in the last block"
        );
        assert_eq!(
            first.record.expect("a body").region.start,
            last.first_position.position
        );
        assert_eq!(
            walk.count() + 1,
            5,
            "five records in the fixture's last block"
        );
    }

    /// **`live_reads` answers about the record just handed back**, and the answer moves as reads
    /// arrive and leave. It is what the residual observation's chain ids are derived from (spec
    /// psp_chain_id_encoding §5).
    ///
    /// ⚠ **This needs a fixture of its own, and the first version of it did not have one.**
    /// Written against `a_finished_psp`, whose every observation carries an empty `chain_ids`,
    /// the live set is empty at all 40 records — so an accessor answering about the wrong
    /// record, or about the block's opening set, passes exactly as the right one does. That is
    /// the second of Milestone F's four shapes: a fixture where several wrong implementations
    /// agree.
    #[test]
    fn live_reads_names_the_reads_live_at_the_record_just_yielded() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        for (at, ids) in [(1u64, vec![5u64, 9]), (101, vec![9, 12]), (201, vec![12])] {
            let mut record = a_record(0, at, 1);
            record.observations[0].chain_ids = ids;
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(psp.block_index().len(), 1, "one block, so one live set");
        let mut walk = psp.records().expect("the walk starts");
        let expected: [&[u64]; 3] = [&[5, 9], &[9, 12], &[12]];
        for want in expected {
            let _ = walk
                .next()
                .expect("a record")
                .expect("a finished psp walks");
            assert_eq!(walk.live_reads().ids(), want);
        }
        assert!(walk.next().is_none());
    }

    /// **A coordinate inside a block starts at that block's first record, not at the
    /// coordinate.** A reader cannot start mid-block (spec §1.2), so what comes back begins
    /// before what was asked for — and that is the contract, not a defect.
    #[test]
    fn a_coordinate_inside_a_block_starts_at_that_blocks_first_record() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();

        for (ordinal, entry) in entries.iter().enumerate() {
            // 50 bases in: the fixture puts its records 100 bases apart, so this lands inside
            // the block and past its first record without reaching the second.
            let asked = GenomePosition {
                contig: entry.first_position.contig,
                position: Position(entry.first_position.position.get() + 50),
            };
            if let Some(next) = entries.get(ordinal + 1) {
                assert!(
                    next.first_position > asked,
                    "the fixture must keep {asked:?} inside block {ordinal}"
                );
            }
            let first = psp
                .records_from(asked)
                .expect("the walk starts")
                .next()
                .expect("a block holds at least one record")
                .expect("a finished psp walks")
                .record
                .expect("the body was built")
                .region;
            assert_eq!(
                (first.contig, first.start),
                (entry.first_position.contig, entry.first_position.position),
                "asking for {asked:?} began at {first:?}"
            );
            assert!(
                first.start.get() < asked.position.get(),
                "the walk began before the coordinate asked for"
            );
        }
    }

    /// **A coordinate in front of every block starts at the first block**, rather than being
    /// refused: it is what asking for a whole contig from position 1 looks like.
    #[test]
    fn a_coordinate_before_every_block_starts_at_the_first() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let first_block = psp.block_index()[0].first_position;
        assert!(
            first_block.position.get() > 0,
            "the fixture's first record is not at position 0"
        );

        let ahead = GenomePosition {
            contig: first_block.contig,
            position: Position(0),
        };
        let walked = psp.records_from(ahead).expect("the walk starts").count();
        assert_eq!(walked, 40, "everything, as `records` would give");
    }

    /// **Restart equals sequential, through the file's own surface.** Starting at block *n*
    /// gives exactly the tail a full walk gives from that block — which is the property every
    /// running difference resetting at a block boundary exists for (spec §3.2), tested here
    /// against a real file rather than a buffer.
    ///
    /// The block each record belongs to is taken from the walk's own count of blocks opened,
    /// not from comparing coordinates: two blocks may legitimately share a first position
    /// (`index.rs`), so a coordinate does not identify one.
    #[test]
    fn a_walk_from_any_block_is_the_tail_of_a_walk_from_the_first() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let blocks = psp.block_index().len();

        let mut whole = Vec::new();
        {
            let mut walk = psp.records().expect("the walk starts");
            while let Some(found) = walk.next() {
                let found = found.expect("a finished psp walks");
                whole.push((walk.blocks_begun() - 1, found.record.expect("a body")));
            }
        }
        assert_eq!(whole.len(), 40);

        for ordinal in 0..blocks {
            let expected: Vec<_> = whole
                .iter()
                .filter(|(block, _)| *block >= ordinal as u64)
                .map(|(_, record)| record.clone())
                .collect();
            assert!(
                !expected.is_empty(),
                "block {ordinal} holds records in the full walk"
            );
            let found: Vec<_> = psp
                .records_from_block(ordinal)
                .expect("the walk starts")
                .map(|found| found.expect("a finished psp walks").record.expect("a body"))
                .collect();
            assert_eq!(found, expected, "restarting at block {ordinal}");
        }
    }

    /// **A record that begins in one block and reaches into the next is not in a walk that
    /// starts at the next.** `records_from` selects on where records *start*, not on what they
    /// span, and the index cannot make it otherwise: an entry carries a block's first position
    /// and nothing else, spec §3.3 having removed the only field that could say how far a
    /// block's records reach.
    ///
    /// The fixture is a deletion-shaped record — 300 bases from position 900, so it covers
    /// 1,100 — followed by records in the next grid cell. Both halves are asserted: that the
    /// file really does hold a record covering 1,100, and that asking for 1,100 does not
    /// return it.
    #[test]
    fn a_record_spanning_into_the_next_block_is_not_in_a_walk_that_starts_there() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        for record in [
            a_record(0, 900, 300),
            a_record(0, 1_050, 1),
            a_record(0, 1_150, 1),
        ] {
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(
            psp.block_index().len(),
            2,
            "the record at 1,050 crosses into the next grid cell and cuts a block"
        );

        // Half one: the file holds a record that starts at 900 and covers 1,100.
        let spanning = psp
            .records()
            .expect("the walk starts")
            .next()
            .expect("the first record")
            .expect("a finished psp walks")
            .record
            .expect("a body")
            .region;
        assert_eq!(spanning.start, Position(900));
        assert!(
            spanning.end >= Position(1_100),
            "the fixture's first record covers 1,100: {spanning:?}"
        );

        // Half two: asking for 1,100 does not return it.
        let asked = GenomePosition {
            contig: ContigId(0),
            position: Position(1_100),
        };
        let starts: Vec<_> = psp
            .records_from(asked)
            .expect("the walk starts")
            .map(|found| {
                found
                    .expect("a finished psp walks")
                    .record
                    .expect("a body")
                    .region
                    .start
            })
            .collect();
        assert_eq!(
            starts,
            vec![Position(1_050), Position(1_150)],
            "the walk begins at the second block's first record"
        );
    }

    /// A block ordinal the file does not have is refused, and the refusal says how many there
    /// are — it is the caller's mistake, and a cohort holding thousands of these open must lose
    /// one sample rather than the run.
    #[test]
    fn a_block_ordinal_the_file_does_not_have_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let blocks = psp.block_index().len();

        let refused = psp
            .records_from_block(blocks)
            .expect_err("one past the last block is not a block");
        match refused {
            PspReadError::NoSuchBlock {
                ordinal_asked_for,
                blocks_in_the_file,
                ..
            } => {
                assert_eq!(ordinal_asked_for, blocks as u64);
                assert_eq!(blocks_in_the_file, blocks as u64);
            }
            other => panic!("got {other}"),
        }
        // And the last block that does exist still walks.
        assert!(psp.records_from_block(blocks - 1).is_ok());
    }

    /// A sample with no records walks to nothing, from either entry point — it is a file
    /// `finish` writes and `open` accepts, so neither walk may refuse it.
    #[test]
    fn a_psp_with_no_records_walks_to_nothing() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let _ = writer.finish(b"a summary of nothing").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert!(psp.block_index().is_empty());
        assert_eq!(psp.records().expect("the walk starts").count(), 0);
        let anywhere = GenomePosition {
            contig: ContigId(0),
            position: Position(1_000),
        };
        assert_eq!(
            psp.records_from(anywhere).expect("the walk starts").count(),
            0
        );
        // **And by ordinal it is a refusal, not an empty walk**: block 0 of a file with no
        // blocks is a block the file does not have, which is the caller's mistake.
        assert!(matches!(
            psp.records_from_block(0),
            Err(PspReadError::NoSuchBlock {
                blocks_in_the_file: 0,
                ..
            })
        ));
    }

    /// **A block that will not inflate is named, and the walk ends there.** The records of the
    /// blocks before it are handed over first — a corrupt block costs the file from that point,
    /// not from the beginning — and the error names the block by the ordinal
    /// `block_index` uses.
    #[test]
    fn a_block_that_will_not_inflate_names_itself_and_ends_the_walk() {
        let (_dir, path) = a_finished_psp();
        let wrecked_block = 2_usize;
        wreck_the_block(&path, wrecked_block);

        let mut psp = PspReader::open(&path).expect("the index and footer are untouched");
        let mut walk = psp.records().expect("the walk starts");
        let mut good = 0;
        let refused = loop {
            match walk.next() {
                Some(Ok(_)) => good += 1,
                Some(Err(refused)) => break refused,
                None => panic!("a wrecked block was walked past"),
            }
        };
        assert_eq!(good, 10, "five records in each of the two blocks before it");
        match &refused {
            PspReadError::CorruptBlock { block, .. } => {
                assert_eq!(*block, wrecked_block as u64);
            }
            other => panic!("got {other}"),
        }
        // **The cause is kept**, so a caller walking the chain reaches what zstd said rather
        // than a sentence about it.
        let cause = std::error::Error::source(&refused).expect("the block reader's own error");
        assert!(
            cause.downcast_ref::<BlockReadError>().is_some(),
            "got {cause}"
        );
        assert!(walk.next().is_none(), "a refused walk is finished");
    }

    /// **A refusal that has a codec's error underneath it keeps it, and one that does not says
    /// so.** Four different causes reached `Damaged` as one string until this walk was built;
    /// `source()` is what tells a caller which rule broke.
    #[test]
    fn a_damaged_file_carries_the_codec_error_underneath_when_there_is_one() {
        // A rule `open` checks itself, having the file's length: no error underneath.
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        footer.trailer_bytes += 1;
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);
        let refused = PspReader::open(&path).expect_err("the sections no longer end at the footer");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
        assert!(
            std::error::Error::source(&refused).is_none(),
            "no codec refused: `open` measured the file itself"
        );

        // A rule the block index's own decoder checks: its error is kept as the cause.
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        let index_at = footer.index_offset as usize;
        let index_ends = index_at + footer.index_bytes as usize;
        // A varint whose continuation bit never clears: the first entry's contig runs off the
        // end of the index.
        whole[index_at..index_ends].fill(0x80);
        footer.index_checksum = crate::ng::psp::index::checksum_index(&whole[index_at..index_ends]);
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);

        let refused = PspReader::open(&path).expect_err("the index does not decode");
        // **The cause is matched, not downcast.** Telling the footer from the block index is
        // what a caller deciding whether a re-run would help has to do, and a typed cause is
        // what lets it.
        match &refused {
            PspReadError::Damaged {
                reason,
                source: Some(DamageFound::BlockIndex(_)),
                ..
            } => {
                assert_eq!(reason, "the block index does not decode");
            }
            other => panic!("got {other}"),
        }
        // And `reason` does not repeat what the cause says, so a chain reads once.
        let cause = std::error::Error::source(&refused).expect("the index decoder's own error");
        assert_ne!(
            cause.to_string(),
            refused.to_string(),
            "a chain that prints the same sentence twice reads as a broken printer"
        );
    }

    /// **A record larger than this reader lets one record hold is refused by name**, not as
    /// damage. Spec §7 keeps this class apart from a corrupt block because the two arrive at the
    /// same line and want opposite instructions: raise the ceiling, or rebuild the file. A
    /// genuine record can be this large — §8 refuses to fix a maximum record size in the format
    /// — and this test writes one.
    #[test]
    fn a_record_larger_than_the_reader_allows_is_refused_by_name() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        // Two observations of a 400 kb span: the reference bases and the observed bases each
        // reach the record's body, so the body is past the 512 kB one record may hold.
        writer
            .push(&a_record(0, 1, 400_000))
            .expect("a long record is legal");
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut walk = psp.records().expect("the walk starts");
        let refused = walk
            .next()
            .expect("the walk reaches the record")
            .expect_err("and cannot hold it");
        match &refused {
            PspReadError::RecordLargerThanTheReaderAllows {
                block,
                allowed_bytes,
                ..
            } => {
                assert_eq!(*block, 0);
                assert_eq!(*allowed_bytes, ROLLING_BUFFER_CEILING_BYTES);
            }
            other => panic!("got {other}"),
        }
        assert!(walk.next().is_none(), "a refused walk is finished");
    }

    /// **A file whose manifest names a field this build does not read there is *upgrade the
    /// reader*, not *rebuild the file*.** The bytes are not damaged; they are another build's.
    /// Folding this into the corrupt-block class would send whoever met it to re-run a pileup
    /// that would produce the very same file.
    ///
    /// **The file is made by patching a finished psp's header text**, because this writer
    /// refuses such a manifest at `create` — which is the right behaviour and leaves no other
    /// way to build the file a *different* build would have written. The replacement is the
    /// same length as what it replaces, so the header's declared body length, every block
    /// offset and the footer are all still true.
    #[test]
    fn a_manifest_this_build_cannot_read_asks_for_a_newer_reader() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let declared = b"position-offset";
        let renamed = b"position-offsex";
        assert_eq!(declared.len(), renamed.len(), "the header keeps its length");
        let at = whole
            .windows(declared.len())
            .position(|window| window == declared)
            .expect("the manifest names the record's first field");
        assert!(
            whole[at + 1..]
                .windows(declared.len())
                .all(|window| window != declared),
            "the name appears once, so patching it patches the manifest and nothing else"
        );
        whole[at..at + renamed.len()].copy_from_slice(renamed);
        rewrite(&path, &whole);

        // The header itself is well formed, so opening works — this is not a damaged file.
        let mut psp = PspReader::open(&path).expect("the header parses");
        let refused = psp
            .records()
            .expect_err("the records cannot be read against this manifest");
        assert!(
            matches!(refused, PspReadError::UnsupportedRecordEncoding { .. }),
            "got {refused}"
        );
        let cause = std::error::Error::source(&refused).expect("the layout's own error");
        assert!(
            cause.downcast_ref::<BlockReadError>().is_some(),
            "got {cause}"
        );
    }

    /// **A fault in the four bytes that introduce a block names the block that never began**,
    /// which is one past the last entry the index has — and so is *not* an ordinal into
    /// `block_index()`.
    ///
    /// Two junk bytes between the last block and the index put the walk's bound two bytes into
    /// what the stream reads as the next block's length. The footer's own offsets are moved with
    /// them, so the file still opens and the fault is genuinely in the blocks.
    #[test]
    fn a_fault_in_the_bytes_introducing_a_block_names_one_past_the_last() {
        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        let blocks = footer.n_blocks;
        let index_at = footer.index_offset as usize;

        let mut patched = whole[..index_at].to_vec();
        patched.extend_from_slice(&[0x00, 0x00]);
        patched.extend_from_slice(&whole[index_at..footer_at]);
        footer.index_offset += 2;
        footer.trailer_offset += 2;
        patched.extend_from_slice(&encode_footer(&footer));
        rewrite(&path, &patched);

        let mut psp = PspReader::open(&path).expect("the sections still add up");
        let mut walk = psp.records().expect("the walk starts");
        let mut good = 0;
        let refused = loop {
            match walk.next() {
                Some(Ok(_)) => good += 1,
                Some(Err(refused)) => break refused,
                None => panic!("the two junk bytes were walked past"),
            }
        };
        assert_eq!(good, 40, "every real record is handed over first");
        match &refused {
            PspReadError::CorruptBlock { block, .. } => {
                assert_eq!(*block, blocks, "one past the last block the index names");
                assert!(
                    psp.block_index().get(*block as usize).is_none(),
                    "and therefore not an ordinal into `block_index`"
                );
            }
            other => panic!("got {other}"),
        }
    }

    /// **The record-buffer ceiling is a knob, and the refusal names it.** A record this reader
    /// will not hold at the default ceiling is read once the ceiling is raised — which is what
    /// makes *raise the ceiling* an instruction rather than a sentence (spec §7).
    #[test]
    fn a_record_over_the_ceiling_is_read_once_the_ceiling_is_raised() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        writer
            .push(&a_record(0, 1, 400_000))
            .expect("a long record is legal");
        let _ = writer.finish(b"").expect("it finishes");

        let mut raised = PspReader::open(&path)
            .expect("a finished psp opens")
            .with_a_record_buffer_ceiling(4 * 1024 * 1024)
            .expect("a ceiling above the buffer");
        assert_eq!(raised.record_buffer_ceiling_bytes(), 4 * 1024 * 1024);
        let record = raised
            .records()
            .expect("the walk starts")
            .next()
            .expect("the record")
            .expect("and it is read")
            .record
            .expect("a body");
        assert_eq!(record.region.start, Position(1));
        assert_eq!(record.reference_bases.len(), 400_000);
    }

    /// **A ceiling at or under the buffer it is a ceiling on is refused where it is set**, not
    /// at the record that would have tripped over it — it would turn away records the buffer
    /// already holds without growing at all.
    #[test]
    fn a_record_buffer_ceiling_under_the_buffer_is_refused_at_the_setting() {
        let (_dir, path) = a_finished_psp();
        let refused = PspReader::open(&path)
            .expect("a finished psp opens")
            .with_a_record_buffer_ceiling(ROLLING_BYTES)
            .expect_err("a ceiling at the buffer is no ceiling");
        match refused {
            PspReadError::RecordBufferCeilingTooSmall {
                ceiling,
                buffer_bytes,
                ..
            } => {
                assert_eq!(ceiling, ROLLING_BYTES);
                assert_eq!(buffer_bytes, ROLLING_BYTES);
            }
            other => panic!("got {other}"),
        }
    }

    /// **A walk that started at block 4 and fails in block 5 says 5.** Every other failing-walk
    /// test starts at block 0, where the walk's own count and the file's ordinal are the same
    /// number — so nothing there separates the two, and deleting the `first_block` addend left
    /// the whole suite green.
    #[test]
    fn a_walk_from_a_later_block_names_the_failing_block_by_its_own_ordinal() {
        let (_dir, path) = a_finished_psp();
        let wrecked_block = 5_usize;
        wreck_the_block(&path, wrecked_block);

        let mut psp = PspReader::open(&path).expect("the index and footer are untouched");
        let mut walk = psp.records_from_block(4).expect("the walk starts");
        let mut good = 0;
        let refused = loop {
            match walk.next() {
                Some(Ok(_)) => good += 1,
                Some(Err(refused)) => break refused,
                None => panic!("a wrecked block was walked past"),
            }
        };
        assert_eq!(good, 5, "the five records of block 4");
        match &refused {
            PspReadError::CorruptBlock { block, .. } => {
                assert_eq!(*block, wrecked_block as u64, "the block's absolute ordinal");
            }
            other => panic!("got {other}"),
        }
    }

    /// A finished psp holding exactly one record of the given kind.
    fn a_one_block_psp(
        kind: crate::ng::locus_generation::LocusKind,
    ) -> (tempfile::TempDir, std::path::PathBuf) {
        let (dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        let mut record = a_record(0, 1, 1);
        record.kind = kind;
        writer.push(&record).expect("in order");
        let _ = writer.finish(b"").expect("it finishes");
        (dir, path)
    }

    /// Where that file's one block sits, and what it holds once inflated.
    fn the_only_block(path: &Path) -> (usize, Vec<u8>) {
        let whole = bytes_of(path);
        let psp = PspReader::open(path).expect("a finished psp opens");
        assert_eq!(psp.block_index().len(), 1, "one block");
        let at = psp.block_index()[0].block_offset as usize;
        let ends = psp.footer().index_offset as usize;
        let declared = u32::from_le_bytes(
            whole[at..at + COMPRESSED_BLOCK_LENGTH_BYTES]
                .try_into()
                .expect("four bytes"),
        ) as usize;
        let frame_at = at + COMPRESSED_BLOCK_LENGTH_BYTES;
        assert_eq!(frame_at + declared, ends, "the block runs up to the index");
        let payload = zstd::decode_all(&whole[frame_at..frame_at + declared])
            .expect("the block's own frame inflates");
        (at, payload)
    }

    /// Put `payload` back in the file as its one block, recompressed, with the footer's offsets
    /// moved to wherever the new frame ends — so the file still opens and only the bytes under
    /// test have changed.
    fn rewrite_the_only_block(path: &Path, block_at: usize, payload: &[u8]) {
        use crate::ng::psp::block::BlockCompressor;

        let header = PspReader::open(path)
            .expect("a finished psp opens")
            .header()
            .clone();
        let mut compressor =
            BlockCompressor::from_manifest(&header.manifest).expect("the file's own compressor");
        // `compress` hands back the whole on-disk block: its four-byte length, then its frame.
        let reframed = compressor
            .compress(payload)
            .expect("it compresses")
            .to_vec();

        let whole = bytes_of(path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        let mut patched = Vec::new();
        patched.extend_from_slice(&whole[..block_at]);
        patched.extend_from_slice(&reframed);
        let index_at = footer.index_offset as usize;
        let new_index_at = patched.len() as u64;
        let shift = new_index_at as i64 - footer.index_offset as i64;
        patched.extend_from_slice(&whole[index_at..footer_at]);
        footer.index_offset = new_index_at;
        footer.trailer_offset = (footer.trailer_offset as i64 + shift) as u64;
        patched.extend_from_slice(&encode_footer(&footer));
        rewrite(path, &patched);
    }

    /// A finished psp of one record whose locus-kind tag no build knows — what a **newer
    /// writer's** file looks like to this reader, framed correctly and unreadable in its head.
    ///
    /// **The byte is located rather than hard-coded**: the same record is written twice, once
    /// `Generic` and once `SsrBundle`, and the two decompressed payloads differ in exactly the
    /// kind tag. The block is then recompressed and the footer's offsets moved with it, so the
    /// file still opens and only the record's kind is beyond this build.
    fn a_psp_whose_record_this_build_cannot_read() -> (tempfile::TempDir, std::path::PathBuf) {
        use crate::ng::locus_generation::LocusKind;

        let (dir_a, path_a) = a_one_block_psp(LocusKind::Generic);
        let (_dir_b, path_b) = a_one_block_psp(LocusKind::SsrBundle);
        let (block_at, mut generic) = the_only_block(&path_a);
        let (_, bundle) = the_only_block(&path_b);
        assert_eq!(generic.len(), bundle.len(), "each kind tag is one byte");
        let differ: Vec<usize> = (0..generic.len())
            .filter(|&index| generic[index] != bundle[index])
            .collect();
        assert_eq!(differ.len(), 1, "exactly one byte differs");
        let kind_at = differ[0];
        assert_eq!((generic[kind_at], bundle[kind_at]), (0, 2), "the kind tags");

        // A tag no kind uses.
        generic[kind_at] = 7;
        rewrite_the_only_block(&path_a, block_at, &generic);
        (dir_a, path_a)
    }

    /// A finished psp of one repeat-tract record whose **body** this build cannot read, its head
    /// entirely sound.
    ///
    /// **The fault is the stored motif emptied to no bases**, which `Motif::new` refuses where
    /// the motif is minted. It is one byte's value and not a length, so the body is exactly as
    /// long as its head declares and everything before the motif reads normally — which is what
    /// makes this a body-only fault rather than a damaged record.
    ///
    /// **The byte is located from the end**: a tract's motif and its two flanks are the body's
    /// last three fields, in that order, and this record's flanks are empty — so the last five
    /// bytes of the block are the motif's length, its two bases, and a zero length for each
    /// flank. The assertion below is what fails if that stops being true.
    fn a_psp_whose_record_body_this_build_cannot_read() -> (tempfile::TempDir, std::path::PathBuf) {
        use crate::ng::locus_generation::{LocusKind, SsrDetail};
        use crate::ng::types::Motif;

        let (dir, path) = a_one_block_psp(LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
            left_flank: Box::from(&b""[..]),
            right_flank: Box::from(&b""[..]),
        }));
        let (block_at, mut payload) = the_only_block(&path);
        let at_the_motif_length = payload.len() - 5;
        assert_eq!(
            &payload[at_the_motif_length..],
            &[2, b'A', b'T', 0, 0],
            "the block ends with the tract's motif and its two empty flanks"
        );

        // A motif of no bases: read as a motif, refused as one.
        payload[at_the_motif_length] = 0;
        rewrite_the_only_block(&path, block_at, &payload);
        (dir, path)
    }

    /// **A record naming a locus kind this build does not know is *upgrade the reader*, not
    /// *rebuild the file*.** The file is another build's, not a damaged one — and re-running a
    /// pileup would produce the very same bytes.
    ///
    /// It is the record-level arm of the *upgrade the reader* class.
    /// `a_manifest_this_build_cannot_read_asks_for_a_newer_reader` reaches the same class
    /// through the manifest, which is a different arm a hundred lines earlier.
    #[test]
    fn a_record_naming_a_locus_kind_this_build_does_not_know_asks_for_a_newer_reader() {
        let (_dir, path) = a_psp_whose_record_this_build_cannot_read();
        let mut psp = PspReader::open(&path).expect("the header and index are untouched");
        let mut walk = psp.records().expect("the walk starts");
        let refused = walk
            .next()
            .expect("the walk reaches the record")
            .expect_err("and cannot read its kind");
        assert!(
            matches!(refused, PspReadError::UnsupportedRecordEncoding { .. }),
            "got {refused}"
        );
    }

    /// **A declined record's body is never decoded, not merely dropped.** That claim is what the
    /// selective walk is *for*, and only a body a full walk **cannot** read can hold it: this
    /// file's one record is a repeat tract whose stored motif has been emptied, so a walk that
    /// builds it is refused where the motif is minted and a walk that declines it walks past
    /// without ever looking.
    ///
    /// ⚠ **The fault has to be in the body, and the locus kind stopped qualifying.** This test
    /// used to corrupt the kind tag, which was the body's last field until it moved into the
    /// head — where an unknown tag now refuses a *declining* walk too, because a reader meets it
    /// while deciding whether it wants the record at all.
    ///
    /// ⚠ **Every other test of the skip checks only that the body came back `None`**, which a
    /// decode-then-discard implementation satisfies exactly as well.
    #[test]
    fn a_declined_records_body_is_never_decoded() {
        let (_dir, path) = a_psp_whose_record_body_this_build_cannot_read();

        // The half that gives the other half its meaning: building this body does fail.
        let refused = PspReader::open(&path)
            .expect("the header and index are untouched")
            .records()
            .expect("the walk starts")
            .next()
            .expect("the walk reaches the record")
            .expect_err("and the emptied motif cannot be minted");
        // The displayed message stops at the block; the field it could not read is in the cause
        // chain under it, which is what says the fault was the motif and not framing.
        assert!(
            format!("{refused:?}").contains("repeat-motif"),
            "the refusal names the field it could not read: {refused:?}"
        );

        let mut psp = PspReader::open(&path).expect("the header and index are untouched");
        let mut walk = psp.records_where(|_| false).expect("the walk starts");
        let found = walk
            .next()
            .expect("the record still arrives")
            .expect("its body was never decoded, so nothing in it can be refused");
        assert!(found.record.is_none());
        assert_eq!(found.head.region.start, Position(1));
        assert!(walk.next().is_none(), "and the walk ends cleanly");
    }

    /// **A corrupt psp is an input, not a bug** (spec §6.7): bytes flipped inside the blocks are
    /// refused, never panicked on, and never followed by another record.
    #[test]
    fn a_psp_with_damaged_blocks_walks_without_panicking() {
        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let (blocks_start, blocks_end) = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            (
                psp.block_index()[0].block_offset as usize,
                psp.footer().index_offset as usize,
            )
        };
        assert!(blocks_end - blocks_start > 64, "there are blocks to damage");

        // A seeded xorshift, so a failure is reproducible from the source alone.
        let mut state: u64 = 0x5eed_1234_9abc_def0;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let (mut refused, mut opened) = (0u32, 0u32);
        for _ in 0..400 {
            let mut bytes = whole.clone();
            let flips = 1 + (next() % 4) as usize;
            for _ in 0..flips {
                let at = blocks_start + (next() as usize % (blocks_end - blocks_start));
                bytes[at] = (next() % 256) as u8;
            }
            rewrite(&path, &bytes);
            let Ok(mut psp) = PspReader::open(&path) else {
                continue;
            };
            opened += 1;
            let mut walk = psp.records().expect("the walk starts");
            let mut broke = false;
            while let Some(found) = walk.next() {
                match found {
                    Ok(_) => assert!(!broke, "a record arrived after a refusal"),
                    Err(_) => {
                        broke = true;
                        assert!(walk.next().is_none(), "a refused walk is finished");
                        break;
                    }
                }
            }
            if broke {
                refused += 1;
            }
        }
        assert!(opened > 0);
        assert!(refused > 0, "the damage must reach the walk");
    }

    /// **A block index a reader cannot trust.** Its checksum is a CRC, so anyone who edits a psp
    /// can restore it; `open` then accepts offsets that point at the middle of a block, and the
    /// walk reads four arbitrary bytes as a block length. Every entry point, no panic, no hang.
    #[test]
    fn a_hostile_block_index_walks_without_panicking() {
        use crate::ng::psp::index::{checksum_index, encode_index};

        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let (entries, blocks_start, blocks_end) = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            (
                psp.block_index().to_vec(),
                psp.block_index()[0].block_offset,
                psp.footer().index_offset,
            )
        };

        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let mut opened = 0u32;
        for _ in 0..300 {
            // Offsets scattered inside the blocks, kept strictly increasing so the index
            // decoder accepts them: this is what a hand-edited psp looks like.
            let mut moved = entries.clone();
            let span = blocks_end - blocks_start;
            let mut offsets: Vec<u64> = (0..moved.len())
                .map(|_| blocks_start + next() % span)
                .collect();
            offsets.sort_unstable();
            offsets.dedup();
            while offsets.len() < moved.len() {
                offsets.push(offsets.last().copied().unwrap_or(blocks_start) + 1);
            }
            for (entry, at) in moved.iter_mut().zip(&offsets) {
                entry.block_offset = *at;
            }

            let index_bytes = encode_index(&moved);
            let footer_at = whole.len() - FOOTER_BYTES;
            let mut footer = footer_of(&whole);
            let old_index_at = footer.index_offset as usize;
            let old_index_ends = old_index_at + footer.index_bytes as usize;

            let mut patched = Vec::new();
            patched.extend_from_slice(&whole[..old_index_at]);
            patched.extend_from_slice(&index_bytes);
            patched.extend_from_slice(&whole[old_index_ends..footer_at]);
            let shift = index_bytes.len() as i64 - footer.index_bytes as i64;
            footer.index_bytes = index_bytes.len() as u64;
            footer.trailer_offset = (footer.trailer_offset as i64 + shift) as u64;
            footer.index_checksum = checksum_index(&index_bytes);
            patched.extend_from_slice(&encode_footer(&footer));
            rewrite(&path, &patched);

            let Ok(mut psp) = PspReader::open(&path) else {
                continue;
            };
            opened += 1;
            for ordinal in 0..psp.block_index().len() {
                let mut walk = psp
                    .records_from_block(ordinal)
                    .expect("an ordinal the index has");
                while let Some(found) = walk.next() {
                    if found.is_err() {
                        assert!(walk.next().is_none(), "a refused walk is finished");
                        break;
                    }
                }
            }
        }
        assert!(opened > 0, "some hostile index must reach a walk");
    }

    /// **A coordinate on a contig the file does not carry starts at the last block**, which is
    /// the far end of the same rule: nothing starts at or after it, so the walk enters the block
    /// before — the last one there is. Pinned so that changing it is a deliberate act.
    #[test]
    fn a_coordinate_on_a_contig_the_file_does_not_have_starts_at_the_last_block() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let last = *psp.block_index().last().expect("blocks");
        assert!(
            psp.header().contigs.len() < 9,
            "contig 9 must be one the file does not carry"
        );

        let beyond = GenomePosition {
            contig: ContigId(9),
            position: Position(1),
        };
        let starts: Vec<_> = psp
            .records_from(beyond)
            .expect("the walk starts")
            .map(|found| found.expect("a finished psp walks").head.region.start)
            .collect();
        assert_eq!(starts.len(), 5, "the last block's records, and no others");
        assert_eq!(starts[0], last.first_position.position);
    }

    // -----------------------------------------------------------------
    // The head-driven skip (G2)
    // -----------------------------------------------------------------

    /// **A declined record still arrives; what it does not carry is a body.** `records_where`
    /// is not a filter, and a caller reading it as one would take the walk's length for the
    /// number of records it kept.
    ///
    /// The predicate keeps every other record, so the two halves are the same size and neither
    /// count could stand in for the other: 40 records, 20 built.
    #[test]
    fn a_declined_record_arrives_with_its_head_and_no_body() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");

        let mut seen = 0usize;
        let mut built = 0usize;
        let mut heads = Vec::new();
        let walk = psp
            .records_where(|_| {
                seen += 1;
                seen % 2 == 1
            })
            .expect("the walk starts");
        for found in walk {
            let found = found.expect("a finished psp walks");
            heads.push(found.head.region);
            if found.record.is_some() {
                built += 1;
            }
        }
        assert_eq!(heads.len(), 40, "every record's head arrives");
        assert_eq!(built, 20, "and half of them were built");
    }

    /// **The records a predicate keeps are the records a full walk gives**, field for field.
    ///
    /// This is Milestone C3's property through the file: every difference a body carries
    /// restarts at the record, so a skipped body strands nothing — and the failure it guards
    /// against is silent, because the records after a skipped one decode plausibly and wrong.
    #[test]
    fn the_records_a_predicate_keeps_match_a_full_walk() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let whole: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .map(|found| found.expect("a finished psp walks").record.expect("a body"))
            .collect();
        assert_eq!(whole.len(), 40);

        // Every third record, so the kept ones are spread across blocks rather than bunched.
        let mut at = 0usize;
        let kept: Vec<_> = psp
            .records_where(|_| {
                at += 1;
                at.is_multiple_of(3)
            })
            .expect("the walk starts")
            .filter_map(|found| found.expect("a finished psp walks").record)
            .collect();
        let expected: Vec<_> = whole
            .iter()
            .enumerate()
            .filter(|(index, _)| (index + 1) % 3 == 0)
            .map(|(_, record)| record.clone())
            .collect();
        assert_eq!(kept.len(), 13, "40 records, every third");
        assert_eq!(kept, expected);
    }

    /// **The live set is exact after a declined record too**, which is why the chain-id changes
    /// ride in the record's head and not in its body (spec psp_record_encoding §6).
    ///
    /// The fixture's middle record is declined, and the ids that arrive and depart at it still
    /// move — so a caller keeping one record in a hundred can still derive the residual
    /// observation at the ones it keeps.
    #[test]
    fn the_live_set_is_exact_after_a_declined_record() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        for (at, ids) in [(1u64, vec![5u64, 9]), (101, vec![9, 12]), (201, vec![12])] {
            let mut record = a_record(0, at, 1);
            record.observations[0].chain_ids = ids;
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut at = 0usize;
        // Decline the middle record, which is where 5 departs and 12 arrives.
        let mut walk = psp
            .records_where(|_| {
                at += 1;
                at != 2
            })
            .expect("the walk starts");
        let expected: [(&[u64], bool); 3] = [(&[5, 9], true), (&[9, 12], false), (&[12], true)];
        for (live, built) in expected {
            let found = walk
                .next()
                .expect("a record")
                .expect("a finished psp walks");
            assert_eq!(found.record.is_some(), built);
            assert_eq!(walk.live_reads().ids(), live);
        }
        assert!(walk.next().is_none());
    }

    /// **A predicate that declines everything still walks every block and every record.** The
    /// bytes have to arrive either way: skipping saves building the record, not decompressing
    /// it, because a block comes out of zstd sequentially and there is nothing to seek past.
    #[test]
    fn a_predicate_that_declines_everything_still_reads_every_block() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let blocks = psp.block_index().len() as u64;

        let mut walk = psp.records_where(|_| false).expect("the walk starts");
        let mut seen = 0;
        for found in walk.by_ref() {
            let found = found.expect("a finished psp walks");
            assert!(found.record.is_none(), "nothing was wanted");
            seen += 1;
        }
        assert_eq!(seen, 40);
        assert_eq!(walk.blocks_begun(), blocks);
    }

    /// **The predicate spec §6.2 writes is the one that works** — the cohort's first pass,
    /// keeping the records where a read showed something other than the reference.
    ///
    /// ⚠ **The fixture has to separate the run or this proves nothing.** `a_record`'s
    /// observation is a copy of its own reference bases, so on `a_finished_psp` every head reads
    /// zero, the predicate is constant-false, and a reader whose `non_reference_reads` did not
    /// come from the file at all would pass. Every third record here carries a read that
    /// disagrees.
    #[test]
    fn the_cohorts_first_pass_predicate_keeps_the_records_where_something_varied() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let mut varies = Vec::new();
        for step in 0..12u64 {
            let mut record = a_record(0, 1 + step * 100, 1);
            let differs = step.is_multiple_of(3);
            if differs {
                record.observations[0].bases = Box::new(*b"N");
            }
            varies.push(differs);
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");
        assert!(
            varies.iter().any(|it| *it) && varies.iter().any(|it| !it),
            "the head has to separate the run, or this proves nothing"
        );

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let built: Vec<bool> = psp
            .records_where(|head| head.non_reference_reads > 0)
            .expect("the walk starts")
            .map(|found| {
                let found = found.expect("a finished psp walks");
                assert_eq!(
                    found.head.non_reference_reads > 0,
                    found.record.is_some(),
                    "the head the predicate saw and the body it got must agree"
                );
                found.record.is_some()
            })
            .collect();
        assert_eq!(built, varies, "exactly the records where something varied");
    }

    /// **A walk from a coordinate takes a predicate too**, which is what a cohort reading one
    /// region of every sample writes. `records_where` is the whole-file case of the same thing.
    #[test]
    fn a_walk_from_a_coordinate_takes_a_predicate() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();
        // Inside the last block, so the walk covers that block and no more.
        let asked = GenomePosition {
            contig: entries[7].first_position.contig,
            position: Position(entries[7].first_position.position.get() + 50),
        };
        let mut at = 0usize;
        let found: Vec<_> = psp
            .records_from(asked)
            .expect("the walk starts")
            .building_only_where(|_| {
                at += 1;
                at == 1
            })
            .map(|found| {
                let found = found.expect("a finished psp walks");
                (found.head.region.start, found.record.is_some())
            })
            .collect();
        assert_eq!(found.len(), 5, "the last block's five records");
        assert!(found[0].1, "the first was wanted");
        assert!(
            found[1..].iter().all(|(_, built)| !built),
            "and none of the rest was built"
        );
        assert_eq!(found[0].0, entries[7].first_position.position);
    }

    /// **A selective walk refuses the way a full one does**, in the same class and naming the
    /// same block — the classification is one function, not one per walk shape.
    #[test]
    fn a_selective_walk_refuses_a_wrecked_block_the_same_way() {
        let (_dir, path) = a_finished_psp();
        let wrecked_block = 2_usize;
        wreck_the_block(&path, wrecked_block);

        let mut psp = PspReader::open(&path).expect("the index and footer are untouched");
        let mut walk = psp.records_where(|_| false).expect("the walk starts");
        let mut good = 0;
        let refused = loop {
            match walk.next() {
                Some(Ok(_)) => good += 1,
                Some(Err(refused)) => break refused,
                None => panic!("a wrecked block was walked past"),
            }
        };
        assert_eq!(good, 10, "the heads of the two blocks before it");
        match &refused {
            PspReadError::CorruptBlock { block, .. } => assert_eq!(*block, wrecked_block as u64),
            other => panic!("got {other}"),
        }
        assert!(walk.next().is_none(), "a refused walk is finished");
    }

    /// **A walk that declines every body accepts damage a full walk refuses**, and this measures
    /// how much.
    ///
    /// A record's two self-consistency checks — the body length its head declared against the
    /// bytes the body used, and the head's non-reference read count against the body's — are made
    /// while the body is decoded, so a declined record never reaches them. That is inherent and
    /// correct: **you cannot check a body you did not decode.** What it means for a caller is
    /// that a clean selective walk says the file's *framing* held, and no more.
    ///
    /// The fixture is one block of three records; every byte of its decompressed payload is
    /// flipped in turn, the block recompressed, and the file walked twice — once building every
    /// body, once building none.
    #[test]
    fn a_declining_walk_accepts_damage_a_full_walk_refuses() {
        use crate::ng::psp::block::BlockCompressor;

        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000_000)).expect("a header");
        for at in [1u64, 101, 201] {
            writer.push(&a_record(0, at, 4)).expect("in order");
        }
        let _ = writer.finish(b"").expect("it finishes");

        let whole = bytes_of(&path);
        let (block_at, blocks_end, header) = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            assert_eq!(psp.block_index().len(), 1, "one block");
            (
                psp.block_index()[0].block_offset as usize,
                psp.footer().index_offset as usize,
                psp.header().clone(),
            )
        };
        let frame_at = block_at + COMPRESSED_BLOCK_LENGTH_BYTES;
        let payload =
            zstd::decode_all(&whole[frame_at..blocks_end]).expect("the block's frame inflates");
        assert!(payload.len() > 32, "there is a payload to damage");

        let mut compressor =
            BlockCompressor::from_manifest(&header.manifest).expect("the file's own compressor");
        let footer = footer_of(&whole);
        let footer_at = whole.len() - FOOTER_BYTES;

        let (mut both_refuse, mut only_the_full_walk) = (0usize, 0usize);
        for byte in 0..payload.len() {
            let mut damaged = payload.clone();
            damaged[byte] ^= 0xff;
            let reframed = compressor
                .compress(&damaged)
                .expect("it compresses")
                .to_vec();

            let mut footer = footer;
            let mut patched = whole[..block_at].to_vec();
            patched.extend_from_slice(&reframed);
            let moved = patched.len() as i64 - footer.index_offset as i64;
            patched.extend_from_slice(&whole[footer.index_offset as usize..footer_at]);
            footer.index_offset = (footer.index_offset as i64 + moved) as u64;
            footer.trailer_offset = (footer.trailer_offset as i64 + moved) as u64;
            patched.extend_from_slice(&encode_footer(&footer));
            rewrite(&path, &patched);

            let Ok(mut psp) = PspReader::open(&path) else {
                continue;
            };
            let full_refused = psp
                .records()
                .expect("the walk starts")
                .any(|found| found.is_err());
            let declining_refused = psp
                .records_where(|_| false)
                .expect("the walk starts")
                .any(|found| found.is_err());
            match (full_refused, declining_refused) {
                (true, true) => both_refuse += 1,
                (true, false) => only_the_full_walk += 1,
                (false, true) => panic!("a declining walk refused where a full walk did not"),
                (false, false) => {}
            }
        }

        assert!(
            only_the_full_walk > both_refuse,
            "most of the damage a full walk catches is in the bodies a predicate skips: \
             {only_the_full_walk} against {both_refuse} over {} payload bytes",
            payload.len()
        );
    }

    /// **A selective walk names its block by the file's own ordinal**, like a full one — the
    /// addend that makes it absolute is invisible on a walk from block 0, which is where every
    /// other test of it starts.
    #[test]
    fn a_selective_walk_from_a_later_block_names_that_blocks_own_ordinal() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut walk = psp
            .records_from_block(6)
            .expect("the walk starts")
            .building_only_where(|_| false);
        assert_eq!(walk.current_block(), 6, "before the first record");
        let _ = walk
            .next()
            .expect("a record")
            .expect("a finished psp walks");
        assert_eq!(walk.current_block(), 6, "the block the record came from");
        assert_eq!(walk.blocks_begun(), 1, "one block opened so far");
    }

    /// **`records_where` refuses at the walk, not at the record**, and the predicate is never
    /// shown a head from a file whose record encoding this build cannot read.
    #[test]
    fn records_where_refuses_a_manifest_this_build_cannot_read() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let declared = b"position-offset";
        let at = whole
            .windows(declared.len())
            .position(|window| window == declared)
            .expect("the manifest names the record's first field");
        whole[at..at + declared.len()].copy_from_slice(b"position-offsex");
        rewrite(&path, &whole);

        let mut psp = PspReader::open(&path).expect("the header parses");
        let mut called = 0usize;
        let refused = psp
            .records_where(|_| {
                called += 1;
                true
            })
            .expect_err("the records cannot be read against this manifest");
        assert!(
            matches!(refused, PspReadError::UnsupportedRecordEncoding { .. }),
            "got {refused}"
        );
        assert_eq!(called, 0, "the predicate was never shown a head");
    }

    /// **A sample with no records walks selectively to nothing**, rather than being refused —
    /// the same contract `records` carries, through the predicate-taking entry point.
    #[test]
    fn records_where_on_a_psp_with_no_records_is_empty() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert!(psp.block_index().is_empty(), "no blocks");
        let mut walk = psp.records_where(|_| true).expect("the walk starts");
        assert!(walk.next().is_none());
        assert_eq!(walk.blocks_begun(), 0);
        assert_eq!(walk.current_block(), 0);
    }

    /// **A corrupt psp is an input on the selective walk too** (spec §6.7). The skip branch
    /// advances the cursor on a length the file supplied without ever decoding a body, so it is
    /// the branch a damaged file reaches furthest into — and no sweep reached it before.
    #[test]
    fn a_psp_with_damaged_blocks_walks_selectively_without_panicking() {
        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let (blocks_start, blocks_end) = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            (
                psp.block_index()[0].block_offset as usize,
                psp.footer().index_offset as usize,
            )
        };

        // A seeded xorshift, so a failure is reproducible from the source alone.
        let mut state: u64 = 0x5eed_1234_9abc_def0;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let (mut refused, mut opened) = (0u32, 0u32);
        for _ in 0..400 {
            let mut bytes = whole.clone();
            let flips = 1 + (next() % 4) as usize;
            for _ in 0..flips {
                let at = blocks_start + (next() as usize % (blocks_end - blocks_start));
                bytes[at] = (next() % 256) as u8;
            }
            rewrite(&path, &bytes);
            let Ok(mut psp) = PspReader::open(&path) else {
                continue;
            };
            opened += 1;
            let mut walk = psp.records_where(|_| false).expect("the walk starts");
            let mut broke = false;
            while let Some(found) = walk.next() {
                match found {
                    Ok(found) => {
                        assert!(!broke, "a record arrived after a refusal");
                        assert!(found.record.is_none(), "nothing was wanted");
                    }
                    Err(_) => {
                        broke = true;
                        assert!(walk.next().is_none(), "a refused walk is finished");
                        break;
                    }
                }
            }
            if broke {
                refused += 1;
            }
        }
        assert!(opened > 0);
        assert!(refused > 0, "the damage must reach the selective walk");
    }
}
