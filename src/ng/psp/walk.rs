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
use super::block::{BlockReadError, BlockStream, StreamedRecord};
use super::chain_ids::LiveSet;
use super::header::Manifest;

/// Seek to a block and hand back a walk that ends where the blocks end.
///
/// `blocks_end` is the file's index offset: the first byte that is no longer a block.
/// `first_block` is the ordinal of the block at `block_offset`, carried only so that a failure
/// can name the block it happened in.
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
    blocks_end: u64,
    block_offset: u64,
    first_block: u64,
) -> Result<RecordIter<'a>, PspReadError> {
    // `open` proved every block offset lies in the blocks, and the empty-file case is given the
    // index's own offset, so this never saturates; it is written this way because an underflow
    // here would be a sixteen-exabyte read rather than an error.
    let blocks_bytes = blocks_end.saturating_sub(block_offset);
    file.seek(SeekFrom::Start(block_offset))
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "seeking to a block",
            source,
        })?;
    let stream = BlockStream::new(file.take(blocks_bytes), manifest).map_err(|refused| {
        match refused {
            // **Upgrade the reader, not rebuild the file**: the manifest declares a field
            // encoding this build does not know.
            BlockReadError::UnsupportedRecordLayout { .. } => {
                PspReadError::UnsupportedRecordEncoding {
                    path: path.to_path_buf(),
                    source: refused,
                }
            }
            // The only other way `new` fails is a look-back window outside the format's range,
            // which `open` has already refused by parsing the header — so this arm is
            // unreachable through a `PspReader`. It is written out rather than left to
            // `unreachable!` because the day the header stops checking is not the day to learn
            // that a reader panics on a bad manifest.
            other => PspReadError::damaged_by(
                path,
                format!("the file's manifest cannot drive a reader: {other}"),
                other,
            ),
        }
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

impl RecordIter<'_> {
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
    /// block already answers 1. Added to the ordinal the walk started at, it names the block a
    /// record came from — which nothing else can, since two blocks may share a first position
    /// (`index.rs`).
    pub fn blocks_read(&self) -> u64 {
        self.stream.blocks_begun()
    }

    /// Put a walk's refusal in the class whose instruction fits it, and name the block.
    ///
    /// **Three instructions, not one** (spec §7): the file is damaged, the reader is too old,
    /// or a limit needs raising. Folding them together is what [`PspReadError::CorruptBlock`]'s
    /// own doc warns against, because *rebuild the file* is wrong advice for two of the three.
    fn refuse(&self, refused: BlockReadError) -> PspReadError {
        let path = self.path.to_path_buf();
        let begun = self.stream.blocks_begun();
        let block = match refused {
            // The four bytes that introduce a block are read before the block is counted, so a
            // fault in them is about a block that never began — one past the last that did.
            BlockReadError::FileEndsInsideABlockLength { .. } => self.first_block + begun,
            _ => self.first_block + begun.saturating_sub(1),
        };
        match refused {
            BlockReadError::RecordLargerThanTheReaderAllows { allowed_bytes, .. } => {
                PspReadError::RecordLargerThanTheReaderAllows {
                    path,
                    block,
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
            BlockReadError::Io {
                while_doing,
                source,
            } => PspReadError::Io {
                path,
                while_doing,
                source,
            },
            // Everything else is the file disagreeing with itself: a frame that will not
            // inflate, a record running past its block, a block holding more than it declared.
            damage => PspReadError::CorruptBlock {
                path,
                block,
                source: damage,
            },
        }
    }
}

impl Iterator for RecordIter<'_> {
    type Item = Result<StreamedRecord, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.stream.next_record() {
            Some(Ok(record)) => Some(Ok(record)),
            // **Classified after the stream has been asked**, so the block count names the block
            // the failure happened in rather than the one before it.
            Some(Err(refused)) => Some(Err(self.refuse(refused))),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::psp::block::{COMPRESSED_BLOCK_LENGTH_BYTES, ROLLING_BUFFER_CEILING_BYTES};
    use crate::ng::psp::footer::{FOOTER_BYTES, decode_footer, encode_footer};
    use crate::ng::psp::writer::PspWriter;
    use crate::ng::psp::writer::tests_support::{
        a_file, a_finished_psp, a_header, a_record, a_sample, bytes_of, rewrite,
    };
    use crate::ng::psp::{PspReadError, PspReader};
    use crate::ng::types::{ContigId, GenomePosition, Position};
    // -----------------------------------------------------------------

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
            walk.blocks_read(),
            blocks,
            "the walk opened every block the index names, and nothing else"
        );
    }

    /// **`records_from` on a block's own first position starts that block.** Checked for every
    /// block in the file, not for one.
    #[test]
    fn records_from_a_blocks_first_position_starts_that_block() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let entries = psp.block_index().to_vec();
        assert!(entries.len() >= 8);

        for entry in &entries {
            let first = psp
                .records_from(entry.first_position)
                .expect("the walk starts")
                .next()
                .expect("a block holds at least one record")
                .expect("a finished psp walks");
            let region = first.record.expect("the body was built").region;
            assert_eq!(
                (region.contig, region.start),
                (entry.first_position.contig, entry.first_position.position),
                "the walk from {:?} began at {region:?}",
                entry.first_position
            );
        }
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
                whole.push((walk.blocks_read() - 1, found.record.expect("a body")));
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
                asked_for,
                blocks: n,
                ..
            } => {
                assert_eq!(asked_for, blocks as u64);
                assert_eq!(n, blocks as u64);
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
    }

    /// **A block that will not inflate is named, and the walk ends there.** The records of the
    /// blocks before it are handed over first — a corrupt block costs the file from that point,
    /// not from the beginning — and the error names the block by the ordinal
    /// `block_index` uses.
    #[test]
    fn a_block_that_will_not_inflate_names_itself_and_ends_the_walk() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let wrecked_block = 2usize;
        let (at, ends) = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            let entries = psp.block_index();
            assert!(entries.len() > wrecked_block + 1);
            (
                entries[wrecked_block].block_offset as usize,
                entries[wrecked_block + 1].block_offset as usize,
            )
        };
        // Everything after the four-byte length: zstd is handed a frame with no magic, which it
        // refuses outright rather than inflating to something plausible.
        whole[at + COMPRESSED_BLOCK_LENGTH_BYTES..ends].fill(0xff);
        rewrite(&path, &whole);

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
        let mut footer = decode_footer(
            &whole[footer_at..]
                .try_into()
                .expect("the file ends with a footer"),
        )
        .expect("a finished file's footer reads");
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
        let mut footer = decode_footer(
            &whole[footer_at..]
                .try_into()
                .expect("the file ends with a footer"),
        )
        .expect("a finished file's footer reads");
        let index_at = footer.index_offset as usize;
        let index_ends = index_at + footer.index_bytes as usize;
        // A varint whose continuation bit never clears: the first entry's contig runs off the
        // end of the index.
        whole[index_at..index_ends].fill(0x80);
        footer.index_checksum = crate::ng::psp::index::checksum_index(&whole[index_at..index_ends]);
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);

        let refused = PspReader::open(&path).expect_err("the index does not decode");
        assert!(
            matches!(refused, PspReadError::Damaged { .. }),
            "got {refused}"
        );
        let cause = std::error::Error::source(&refused).expect("the index decoder's own error");
        assert!(
            cause
                .downcast_ref::<crate::ng::psp::index::IndexDecodeError>()
                .is_some(),
            "got {cause}"
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
}
