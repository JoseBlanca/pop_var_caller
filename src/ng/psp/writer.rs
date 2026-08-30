//! Writing a psp: create a file, push records into it, and finish it.
//!
//! **Finishing is what makes the file readable at all.** Before [`PspWriter::finish`] there is
//! no footer, and every reader refuses a file without one — which is exactly what should happen
//! to a run that was killed, and is goal 3 of the format (spec §6.3). A `PspWriter` dropped
//! without finishing therefore leaves a file no reader will touch, and that is the intended
//! outcome rather than a leak.
//!
//! ```text
//! create   → header
//! push     → blocks, one compressed frame each, as the cut rule closes them
//! finish   → the last block, the index, the trailer, the footer — then flushed and synced
//! ```

use std::fs::File;
use std::io::{BufWriter, Seek, Write};
use std::path::{Path, PathBuf};

use super::block::{BlockBuilder, BlockCompressor, BlockHead, BlockWriteError};
use super::footer::{FOOTER_BYTES, Footer, encode_footer};
use super::header::Header;
use super::index::BlockIndexEntry;
use super::{PspWriteError, footer, index};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion};

/// What one finished file cost, handed back by [`PspWriter::finish`].
///
/// **Counts of what was written, not of what was offered**: a record the writer refused is not
/// here, because it is not in the file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use]
pub struct WriteStats {
    /// How many records reached the file.
    pub records: u64,
    /// How many blocks they were cut into — and therefore how many entries the index holds.
    pub blocks: u64,
    /// The finished file's length in bytes, header and footer included.
    pub bytes: u64,
}

/// Writes one psp, from its header to its footer.
///
/// **The manifest is fixed at [`create`](Self::create) and cannot change afterwards** (spec
/// §6.3): a field's encoding, the genomic block size and the look-back window are all decided
/// before the first record and recorded in the header, because a writer that could change one
/// half-way through would produce a file no reader could interpret without re-reading the header
/// per block.
#[derive(Debug)]
#[must_use]
pub struct PspWriter {
    /// Which file is being written. **Every error this type raises names it**, because a cohort
    /// gathering sixty samples at once raises them from sixty writers at once.
    path: PathBuf,
    out: BufWriter<File>,
    /// Where the next byte will land — the file's length so far. Advanced only by a write that
    /// returned, so it always describes bytes this writer has handed to the buffer.
    written: u64,
    /// **An `Option` because [`BlockBuilder::finish`] consumes the builder**, which is its own
    /// guard: a builder that could be closed twice would put the last block in the file twice.
    /// It is `Some` for the whole of this writer's life, because the only thing that takes it is
    /// [`PspWriter::finish`], which consumes the writer in the same breath.
    builder: Option<BlockBuilder>,
    compressor: BlockCompressor,
    /// One entry per block closed so far, in the order they were written.
    index: Vec<BlockIndexEntry>,
    records: u64,
    /// Why this writer can no longer produce a whole file, once something has gone wrong that
    /// it cannot undo.
    ///
    /// **The block builder hands a closed block over and reopens in the same call**, so once
    /// [`push`](Self::push) holds a payload the builder will never offer it again. Everything
    /// after that point — reading the head back, compressing, writing — can still fail, and a
    /// failure there means those records exist nowhere. Without this flag `finish` then wrote a
    /// **valid** file: footer, matching index checksum, entries in ascending order, and a
    /// thousand bases missing out of the middle of a contig with nothing saying so. That is
    /// worse than the unreadable stump a killed run leaves, because every reader accepts it.
    spent: Option<&'static str>,
}

impl PspWriter {
    /// Create a psp and write its header.
    ///
    /// **⚠ An existing file at `path` is truncated**, the way `File::create` truncates — so
    /// creating over a finished psp destroys it, and there is no footer left to say a sample was
    /// ever there. That is the conventional meaning of *create* and it is what the pipeline
    /// wants when it re-runs a sample, but it is not stated anywhere in the spec, and the spec
    /// *does* warn about the milder destruction `append` causes (§6.4: "write to a new path and
    /// rename if that matters"). A caller that must not destroy an existing psp checks for one
    /// first.
    ///
    /// **Nothing touches the filesystem until the header and the manifest have been accepted.**
    /// The header is encoded, and the cut rule and compressor built from its manifest, before
    /// the file is created — so a header this writer cannot honour leaves no file behind at all,
    /// rather than an empty one that every reader then refuses for a different reason.
    pub fn create(path: &Path, mut header: Header) -> Result<Self, PspWriteError> {
        // **The compression level is recorded because it is a setting, and goal 4 is that
        // settings are recorded.** A reader needs nothing from it — zstd decodes any level — but
        // `append` must match bytes already in a file, and without this it would have no way to
        // learn what those bytes were written at. `block.rs`'s own doc assigns this to F3 by
        // name.
        header.writer.parameters.insert(
            "zstd-compression-level".to_string(),
            crate::ng::psp::header::ParameterValue::Integer(i64::from(
                super::block::ZSTD_COMPRESSION_LEVEL,
            )),
        );
        let header_bytes = header.encode()?;
        let builder = BlockBuilder::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedManifest {
                path: path.to_path_buf(),
                source: source.into(),
            }
        })?;
        let compressor = BlockCompressor::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedManifest {
                path: path.to_path_buf(),
                source: source.into(),
            }
        })?;

        let file = File::create(path).map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "creating the file",
            source,
        })?;
        let mut writer = Self {
            path: path.to_path_buf(),
            out: BufWriter::new(file),
            written: 0,
            builder: Some(builder),
            compressor,
            index: Vec::new(),
            records: 0,
            spent: None,
        };
        writer.put(&header_bytes, "writing the header")?;
        Ok(writer)
    }

    /// Reopen a finished psp and add records to it.
    ///
    /// **The footer says where the blocks end, so appending is truncating at the index offset and
    /// carrying on** (spec §6.4): the old index, trailer and footer are discarded, the blocks and
    /// the header are not touched, and `finish` writes a whole new index made of the entries that
    /// were already there and the ones this writer adds.
    ///
    /// **The header is not rewritten, so the manifest is the file's.** A field's encoding, the
    /// grid the blocks cut on and the look-back window are all already declared, and appended
    /// records must use them — so a manifest this build cannot honour is
    /// [`UnsupportedManifest`](PspWriteError::UnsupportedManifest) here rather than a file whose
    /// added records no reader can interpret under the header it keeps.
    ///
    /// **Coordinate order runs across the seam.** The last record already in the file is found by
    /// walking the last block's heads — the bodies are never built — and the first appended record
    /// is checked against it.
    ///
    /// **⚠ A file being appended to has no footer while this writer holds it open**, exactly like
    /// a new one, so an append interrupted half-way leaves a file every reader refuses. That is
    /// right — but note it has *lost the trailer the file had*, which is not recoverable. Spec
    /// §6.4: **write to a new path and rename if that matters to the caller.**
    pub fn append(path: &Path) -> Result<Self, PspWriteError> {
        let reopen = |source: super::PspReadError| PspWriteError::Reopen {
            path: path.to_path_buf(),
            source,
        };
        // **Opened as a reader first, and every check `open` makes is made.** An append writes a
        // fresh index and footer onto whatever it finds, so a file the reader would refuse must
        // not be extended — the same lesson `replace_trailer` learned from its own review, and
        // here the index is needed anyway.
        let mut psp = super::PspReader::open(path).map_err(reopen)?;
        let header = psp.header().clone();
        let blocks_end = psp.footer().index_offset;
        let index = psp.block_index().to_vec();
        // **The manifest is checked before the file is walked, and the order is the finding.**
        // `BlockBuilder::from_manifest` is what refuses a layout this build cannot write, and the
        // walk below needs the same layout to read — so walking first made a manifest this writer
        // cannot honour arrive as a *reader's* refusal wrapped in `Reopen`, which is the wrong
        // class for what spec §6.4 calls out by name.
        let mut builder = BlockBuilder::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedManifest {
                path: path.to_path_buf(),
                source: source.into(),
            }
        })?;
        let last_record = Self::find_the_last_record_in(&mut psp).map_err(reopen)?;
        drop(psp);
        if let Some(region) = last_record {
            builder = builder.continuing_after(region);
        }
        // **At the level the file records, not this build's** (spec goal 4). `create` writes it
        // into the header's parameters for exactly this: an append that used another level would
        // put blocks in one file compressed two ways, with nothing saying so.
        let compressor =
            Self::build_the_compressor_the_header_records(&header).map_err(|source| {
                PspWriteError::UnsupportedManifest {
                    path: path.to_path_buf(),
                    source,
                }
            })?;

        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .map_err(|source| PspWriteError::Io {
                path: path.to_path_buf(),
                while_doing: "reopening the file to append to it",
                source,
            })?;
        // **The truncation is what removes the footer**, and from here until `finish` writes a new
        // one every reader refuses the file — which is the same state a new writer leaves and is
        // goal 3.
        file.set_len(blocks_end)
            .map_err(|source| PspWriteError::Io {
                path: path.to_path_buf(),
                while_doing: "truncating the file at its block index",
                source,
            })?;
        file.seek(std::io::SeekFrom::Start(blocks_end))
            .map_err(|source| PspWriteError::Io {
                path: path.to_path_buf(),
                while_doing: "seeking to the end of the blocks",
                source,
            })?;

        Ok(Self {
            path: path.to_path_buf(),
            out: BufWriter::new(file),
            // **Where the next byte lands, which is where the blocks ended.** Every offset the
            // new index and footer carry is measured from this, so a zero here would put every
            // appended block at an address inside the file that was already there.
            written: blocks_end,
            builder: Some(builder),
            compressor,
            index,
            // **Not the records already in the file**, which nothing counts: `WriteStats` says
            // what this writer wrote, and the blocks it inherited are in the index it inherited.
            records: 0,
            spent: None,
        })
    }

    /// The region of the last record already in a file, or `None` if it holds none.
    ///
    /// **Heads only.** The walk declines every body, so this costs the last block's framing and
    /// no record building — and the last block is the only one read.
    fn find_the_last_record_in(
        psp: &mut super::PspReader,
    ) -> Result<Option<GenomeRegion>, super::PspReadError> {
        let Some(last_block) = psp.block_index().len().checked_sub(1) else {
            return Ok(None);
        };
        let mut last = None;
        for found in psp
            .records_from_block(last_block)?
            .building_only_where(|_| false)
        {
            last = Some(found?.head.region);
        }
        Ok(last)
    }

    /// Build a compressor at the window the manifest declares and the level the header records.
    ///
    /// **A header that records no level gives this build's**, which is what a file written before
    /// the parameter existed looks like. **Every other shape is refused**, including one recorded
    /// as a string or a float — ⚠ until the G4 review those fell into the same arm as *absent*,
    /// so a file recording `"1"` was appended to at level 9 with nothing said, which is exactly
    /// the file §2.4 says must not exist.
    fn build_the_compressor_the_header_records(
        header: &Header,
    ) -> Result<BlockCompressor, super::ManifestRefusal> {
        use crate::ng::psp::header::ParameterValue;

        let level = match header.writer.parameters.get("zstd-compression-level") {
            Some(ParameterValue::Integer(recorded)) => {
                i32::try_from(*recorded).map_err(|_| super::ManifestRefusal::LevelPastAnyLevel {
                    recorded: *recorded,
                })?
            }
            // **Absent is the one shape that falls back**, and only because it is what a file
            // written before this parameter existed looks like. Any other shape is a level
            // recorded and then ignored, which is the file §2.4 says must not exist.
            Some(other) => {
                return Err(super::ManifestRefusal::UnreadableLevel {
                    recorded: format!("{other:?}"),
                });
            }
            None => super::block::ZSTD_COMPRESSION_LEVEL,
        };
        Ok(BlockCompressor::with_level(
            header.manifest.look_back_window_log,
            level,
        )?)
    }

    /// Lay one record down.
    ///
    /// **Coordinate order is enforced** (spec §6.3): a record starting before the one before it,
    /// or on a contig already finished with, is refused rather than written — a file that breaks
    /// the order seeks wrongly instead of failing, because the index and every seek rest on it.
    ///
    /// **A record refused for its coordinates leaves the file exactly as it was**, and the
    /// writer stays usable: the cut rule guarantees that for its own state, and nothing here
    /// writes a byte until the builder has handed back a whole block.
    ///
    /// ⚠ **Every other failure is unrecoverable, and marks the writer so.** Once the builder has
    /// handed a closed block over it has already reopened, so nothing can offer that block
    /// again; a failure while reading its head back, compressing it or writing it means those
    /// records exist nowhere. The writer records why and [`finish`](Self::finish) then refuses.
    ///
    /// An earlier version said an I/O failure here "is terminal and is meant to be", reasoning
    /// that the file would keep no footer. **That holds only while the failure persists.** A
    /// transient one — one full device that empties, one refused write — left `finish` free to
    /// write a footer over a file with a thousand bases missing from the middle of a contig,
    /// which every reader accepts.
    pub fn push(&mut self, record: &SampleLocusObservations) -> Result<(), PspWriteError> {
        let closed = self
            .builder
            .as_mut()
            .expect("the builder is taken only by `finish`, which consumes the writer")
            .push(record)
            .map_err(|refused| match refused {
                BlockWriteError::OutOfOrder { previous, offered } => PspWriteError::OutOfOrder {
                    path: self.path.clone(),
                    previous,
                    offered,
                },
                other => PspWriteError::RecordRefused {
                    path: self.path.clone(),
                    source: other,
                },
            })?;
        self.records += 1;
        let Some(payload) = closed else {
            return Ok(());
        };
        // **From here on a failure is unrecoverable and is recorded as such.** The builder has
        // already closed this block and reopened; nothing can hand it back.
        let head = Self::decode_the_head_of(payload, &self.path).inspect_err(|_| {
            self.spent = Some("a block it had just built could not be read back");
        })?;
        // **Destructured so the compressed block is written where it lies.** `compress` hands
        // back a borrow of the compressor's own buffer, and copying it out was the price of
        // borrowing `self` twice — one allocation the size of a compressed block, per block, for
        // nothing. The fields it needs are disjoint from the compressor's, and naming them says
        // so to the borrow checker.
        let Self {
            path,
            out,
            written,
            compressor,
            index,
            ..
        } = self;
        let block = match compressor.compress(payload) {
            Ok(block) => block,
            Err(source) => {
                self.spent = Some("a block could not be compressed");
                return Err(PspWriteError::BlockRefused {
                    path: self.path.clone(),
                    source,
                });
            }
        };
        let entry = BlockIndexEntry {
            first_position: GenomePosition {
                contig: head.contig,
                position: head.first_position,
            },
            block_offset: *written,
        };
        let put = out
            .write_all(block)
            .map(|()| *written += block.len() as u64)
            .map_err(|source| PspWriteError::Io {
                path: path.clone(),
                while_doing: "writing a block",
                source,
            });
        match put {
            Ok(()) => {
                index.push(entry);
                Ok(())
            }
            Err(refused) => {
                self.spent = Some("a block could not be written");
                Err(refused)
            }
        }
    }

    /// Write the last block, the index, the trailer and the footer, then make the file durable.
    ///
    /// **It consumes the writer**, which is what says in the type system that a file with no
    /// footer is not a file: there is no way to hold a writer that has finished, and no way to
    /// produce a readable psp without calling this (arch §4.2, §5).
    ///
    /// **⚠ Durability is three steps in this order, and it is easy to get wrong** (spec §6.3):
    /// flush the format, *then* surface the buffered writer's errors, *then* sync. A `BufWriter`
    /// dropped without the middle step can swallow a failed flush, and a truncated footer on a
    /// billions-of-records file looks exactly like an interrupted run.
    ///
    /// **What it writes, it reads back before believing it.** The index and the footer are
    /// decoded by the very functions a reader will use, and a failure is refused here rather
    /// than left on disk — see [`Self::check_the_index_and_footer_read_back`].
    pub fn finish(mut self, trailer: &[u8]) -> Result<WriteStats, PspWriteError> {
        // **A writer that lost a block cannot produce a whole file, and must not produce a file
        // that looks whole.** See the `spent` field for what that looked like before this
        // existed.
        if let Some(why) = self.spent {
            return Err(PspWriteError::WouldNotBeReadable {
                path: self.path.clone(),
                reason: format!("records were lost earlier in the walk: {why}"),
                // The loss happened earlier and its own error went to the caller then.
                source: None,
            });
        }
        let last_block = self
            .builder
            .take()
            .expect("the builder is taken only here, and `finish` consumes the writer")
            .finish();
        if let Some(payload) = last_block {
            let head = Self::decode_the_head_of(&payload, &self.path)?;
            let block = self
                .compressor
                .compress(&payload)
                .map_err(|source| PspWriteError::BlockRefused {
                    path: self.path.clone(),
                    source,
                })?
                .to_vec();
            self.put_block(&head, &block)?;
        }

        let index_offset = self.written;
        let index_bytes = index::encode_index(&self.index);
        let footer = Footer {
            index_offset,
            index_bytes: index_bytes.len() as u64,
            trailer_offset: index_offset + index_bytes.len() as u64,
            trailer_bytes: trailer.len() as u64,
            n_blocks: self.index.len() as u64,
            index_checksum: index::checksum_index(&index_bytes),
        };
        self.check_the_index_and_footer_read_back(&index_bytes, &footer)?;

        self.put(&index_bytes, "writing the block index")?;
        self.put(trailer, "writing the trailer")?;
        self.put(&encode_footer(&footer), "writing the footer")?;

        // **Every byte of the file is accounted for by a section the footer names.** The blocks
        // end where the index starts, the index and the trailer have the lengths the footer
        // carries, and the footer is the last 48 bytes — so the counter and the footer must
        // agree, and if they ever do not, every offset the footer holds is measured from a
        // different file than the one on disk.
        //
        // **Carried forward from the F3 and F4 reviews, and re-measured here**: a defect that
        // adds a fifth section between the trailer and the footer fails **82 tests with this
        // assertion and 66 without it**, so it accounts for sixteen of them. ⚠ The F-era figure
        // was *1 failing test to 7*, on a suite a third this size; a number measured on another
        // tree is not a fact about this one, and the G4 review caught it still standing here.
        //
        // It is a `debug_assert` because it is arithmetic on numbers this function itself
        // produced — a release build gains nothing by re-deriving them. **What it buys is where
        // the failure lands**: at the writer, naming the disagreement, rather than at sixty-six
        // readers.
        debug_assert_eq!(
            self.written,
            footer.trailer_offset + footer.trailer_bytes + FOOTER_BYTES as u64,
            "the bytes written and the sections the footer names must be the same file"
        );

        let stats = WriteStats {
            records: self.records,
            blocks: self.index.len() as u64,
            bytes: self.written,
        };

        // **The durability steps, in the one order that surfaces every failure** (spec §6.3).
        //
        // ⚠ `into_inner` *is* the flush, and it is the flush that surfaces. An earlier version
        // called `flush()` first and then `into_inner()`, which reads like the spec's three
        // steps but is not: a successful `flush()` empties the buffer, so the flush inside
        // `into_inner` then has nothing to write and cannot fail. The error arm was
        // unreachable — deleting it changed no behaviour under mutation, which is how it was
        // found.
        let file = self.out.into_inner().map_err(|failed| PspWriteError::Io {
            path: self.path.clone(),
            while_doing: "flushing the finished file",
            source: failed.into_error(),
        })?;
        file.sync_all().map_err(|source| PspWriteError::Io {
            path: self.path.clone(),
            while_doing: "syncing the finished file",
            source,
        })?;

        Ok(stats)
    }

    /// **A writer must not be able to produce a file its own reader refuses**, so the two
    /// structures with rules of their own are read back with the very functions that will read
    /// them from disk.
    ///
    /// This is the obligation F1 and F2 both routed here rather than solving in their codecs:
    /// neither could assert its rule at the encoder, because each module's own tests must be
    /// able to write the bytes that prove the reader refuses them. `finish` is the one place
    /// where the bytes being written are meant to be readable, so it is the one place the
    /// assertion belongs.
    ///
    /// It costs one index decode and one 48-byte footer decode per file — against a walk that
    /// wrote every block in it.
    fn check_the_index_and_footer_read_back(
        &self,
        index_bytes: &[u8],
        footer: &Footer,
    ) -> Result<(), PspWriteError> {
        index::decode_index(index_bytes, footer.n_blocks).map_err(|source| {
            PspWriteError::WouldNotBeReadable {
                path: self.path.clone(),
                reason: "the block index it wrote does not decode".to_string(),
                source: Some(source.into()),
            }
        })?;
        // **The checksum too**, because a reader checks it before it decodes: a footer carrying
        // a checksum of something other than the index is a file `open` refuses, and computing
        // it over the trailer instead sailed through this guard before the check existed.
        let found = index::checksum_index(index_bytes);
        if found != footer.index_checksum {
            return Err(PspWriteError::WouldNotBeReadable {
                path: self.path.clone(),
                reason: format!(
                    "the footer carries {:#010x} as the block index's checksum; the index it \
                     wrote checksums to {found:#010x}",
                    footer.index_checksum
                ),
                // A rule the writer checks itself: no decoder can see the two together.
                source: None,
            });
        }
        footer::decode_footer(&encode_footer(footer)).map_err(|source| {
            PspWriteError::WouldNotBeReadable {
                path: self.path.clone(),
                reason: "the footer it wrote does not decode".to_string(),
                source: Some(source.into()),
            }
        })?;
        Ok(())
    }

    /// The head of a block this writer has just built, read back out of its own bytes.
    ///
    /// **The index entry comes from the bytes going to disk, not from a number kept beside
    /// them**, so the index describes the file rather than the writer's intentions.
    fn decode_the_head_of(payload: &[u8], path: &Path) -> Result<BlockHead, PspWriteError> {
        BlockHead::decode(payload)
            .map(|decoded| decoded.head)
            .map_err(|source| PspWriteError::WouldNotBeReadable {
                path: path.to_path_buf(),
                reason: "a block it just built does not decode".to_string(),
                source: Some(source.into()),
            })
    }

    /// Write one compressed block and give it its index entry.
    fn put_block(&mut self, head: &BlockHead, block: &[u8]) -> Result<(), PspWriteError> {
        let entry = BlockIndexEntry {
            first_position: GenomePosition {
                contig: head.contig,
                position: head.first_position,
            },
            block_offset: self.written,
        };
        self.put(block, "writing a block")?;
        self.index.push(entry);
        Ok(())
    }

    /// Hand bytes to the buffer and count them.
    fn put(&mut self, bytes: &[u8], while_doing: &'static str) -> Result<(), PspWriteError> {
        self.out
            .write_all(bytes)
            .map_err(|source| PspWriteError::Io {
                path: self.path.clone(),
                while_doing,
                source,
            })?;
        self.written += bytes.len() as u64;
        Ok(())
    }
}

/// Fixtures shared by this module's tests and by the reader's.
///
/// **Here rather than duplicated**, because a reader's tests need a psp that a writer wrote, and
/// two copies of a fixture drift: the one thing worse than a fixture that cannot fail is two
/// fixtures that disagree about what a sample looks like.
#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::psp::header::{
        ContigIdentity, FORMAT_VERSION, ParameterValue, ReferenceIdentity, WriterProvenance,
    };
    use crate::ng::types::{Bp, ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};
    use std::io::Read as _;

    /// A tomato-shaped header, with the block grid small enough that a handful of records cut
    /// several blocks.
    pub(crate) fn a_header(genomic_block_size_bp: u64) -> Header {
        let mut header = Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            contigs: vec![
                ContigIdentity {
                    name: "SL4.0ch01".to_string(),
                    length: 90_863_682,
                    md5: Some([0x1b; 16]),
                },
                ContigIdentity {
                    name: "SL4.0ch02".to_string(),
                    length: 53_473_368,
                    md5: Some([0x2c; 16]),
                },
            ],
            writer: WriterProvenance {
                tool: "ng".to_string(),
                version: "0.1.0".to_string(),
                subcommand: "pileup".to_string(),
                input_alignments: vec!["SRR7279481.cram".to_string()],
                input_reference: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                command_line: "ng pileup --sample SRR7279481".to_string(),
                parameters: std::collections::BTreeMap::from([(
                    "depth-cap".to_string(),
                    ParameterValue::Integer(300),
                )]),
                created: "2026-08-28T00:00:00Z"
                    .parse()
                    .expect("a valid RFC 3339 stamp"),
            },
            manifest: crate::ng::psp::header::Manifest::as_this_build_writes_it(),
        };
        header.manifest.genomic_block_size_bp = Bp(genomic_block_size_bp);
        header
    }

    pub(crate) fn a_record(contig: u32, start: u64, span: u64) -> SampleLocusObservations {
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

    /// Records across two contigs, cutting several blocks on a 1 kb grid.
    pub(crate) fn a_sample() -> Vec<SampleLocusObservations> {
        let mut records = Vec::new();
        for contig in 0..2u32 {
            for block in 0..4u64 {
                for step in 0..5u64 {
                    records.push(a_record(contig, 1 + block * 1_000 + step * 100, 1));
                }
            }
        }
        records
    }

    /// A finished psp holding [`a_sample`], cut on a 1 kb grid, and where it is.
    ///
    /// **One fixture for the writer's tests, the reader's and the walk's.** Three copies of
    /// "what a finished file looks like" drift, and the one thing worse than a fixture that
    /// cannot fail is two that disagree about what a sample is.
    pub(crate) fn a_finished_psp() -> (tempfile::TempDir, PathBuf) {
        let (dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        for record in a_sample() {
            writer.push(&record).expect("in order");
        }
        let _ = writer.finish(b"a per-sample summary").expect("it finishes");
        (dir, path)
    }

    /// The footer of a psp read whole into memory.
    ///
    /// **The three-line preamble every test that damages a file writes**: slice the fixed tail,
    /// widen it to an array, decode it. Here rather than in each of them, for the reason
    /// [`a_finished_psp`] is.
    pub(crate) fn footer_of(bytes: &[u8]) -> crate::ng::psp::footer::Footer {
        let tail: [u8; crate::ng::psp::footer::FOOTER_BYTES] = bytes
            [bytes.len() - crate::ng::psp::footer::FOOTER_BYTES..]
            .try_into()
            .expect("the file is at least a footer long");
        crate::ng::psp::footer::decode_footer(&tail).expect("a finished file's footer reads")
    }

    /// Fill one block's compressed frame with `0xff`, leaving the four-byte length in front of
    /// it, so the block cannot inflate.
    ///
    /// **`0xff` rather than a flipped bit**, and deliberately: zstd refuses a frame with no magic
    /// outright, where a flipped bit inside a frame can inflate to something plausible and make
    /// the test depend on which byte was chosen.
    pub(crate) fn wreck_the_block(path: &Path, ordinal: usize) {
        let mut whole = bytes_of(path);
        let (at, ends) = {
            let psp = crate::ng::psp::PspReader::open(path).expect("a finished psp opens");
            let entries = psp.block_index();
            let at = entries[ordinal].block_offset as usize;
            let ends = entries
                .get(ordinal + 1)
                .map_or(psp.footer().index_offset as usize, |next| {
                    next.block_offset as usize
                });
            (at, ends)
        };
        whole[at + crate::ng::psp::COMPRESSED_BLOCK_LENGTH_BYTES..ends].fill(0xff);
        rewrite(path, &whole);
    }

    /// Replace a file's contents wholesale — how a test lays down a psp it has damaged.
    pub(crate) fn rewrite(path: &Path, bytes: &[u8]) {
        use std::io::Write as _;
        let mut file = File::create(path).expect("the file is writable");
        file.write_all(bytes).expect("it writes");
    }

    pub(crate) fn a_file() -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("SRR7279481.psp");
        (dir, path)
    }

    pub(crate) fn bytes_of(path: &Path) -> Vec<u8> {
        let mut bytes = Vec::new();
        File::open(path)
            .expect("the file exists")
            .read_to_end(&mut bytes)
            .expect("it reads");
        bytes
    }
}

#[cfg(test)]
mod tests {
    use super::tests_support::{
        a_file, a_finished_psp, a_header, a_record, a_sample, bytes_of, footer_of, rewrite,
        wreck_the_block,
    };
    use super::*;
    use crate::ng::psp::NotReadable;
    use crate::ng::psp::PspReader;
    use crate::ng::psp::footer::{FOOTER_BYTES, FOOTER_MAGIC, decode_footer};
    use crate::ng::psp::index::decode_index;
    use crate::ng::types::{ContigId, Position};

    // -----------------------------------------------------------------
    // What a finished file is
    // -----------------------------------------------------------------

    /// **A finished file ends with the footer, and everything the footer says is where it says.**
    /// The index decodes, its checksum matches, the trailer is the payload handed to `finish`,
    /// and the block count is the number of entries.
    #[test]
    fn a_finished_file_is_laid_out_the_way_its_footer_says() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a writable header");
        let records = a_sample();
        for record in &records {
            writer.push(record).expect("the fixture is in order");
        }
        let stats = writer.finish(b"a per-sample summary").expect("it finishes");

        let bytes = bytes_of(&path);
        assert_eq!(stats.bytes, bytes.len() as u64, "the count is the file");
        assert_eq!(stats.records, records.len() as u64);
        assert!(stats.blocks >= 8, "two contigs of four grid cells each");

        let tail: [u8; FOOTER_BYTES] = bytes[bytes.len() - FOOTER_BYTES..]
            .try_into()
            .expect("the file is at least a footer long");
        let footer = decode_footer(&tail).expect("a finished file has a readable footer");
        assert_eq!(footer.n_blocks, stats.blocks);

        let index = &bytes
            [footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize];
        assert_eq!(
            crate::ng::psp::index::checksum_index(index),
            footer.index_checksum
        );
        let entries = decode_index(index, footer.n_blocks).expect("the index reads");
        assert_eq!(entries.len() as u64, stats.blocks);

        let trailer = &bytes[footer.trailer_offset as usize
            ..(footer.trailer_offset + footer.trailer_bytes) as usize];
        assert_eq!(trailer, b"a per-sample summary");

        // The footer starts exactly where the trailer ends: nothing is left over.
        assert_eq!(
            footer.trailer_offset + footer.trailer_bytes,
            (bytes.len() - FOOTER_BYTES) as u64
        );
    }

    /// **Every index entry points at a block whose head names the coordinate the entry
    /// claims.** The offsets are read back out of the file, each is checked to start a whole
    /// compressed block, and the block is inflated and its head decoded so the coordinate can be
    /// compared.
    ///
    /// ⚠ **This test used to stop at "it starts a whole block", while its name promised the
    /// rest.** Shifting every entry's position 100 bases past its block's true head survived all
    /// twelve tests — an index that seeks to the wrong block is exactly the silent-wrong-answer
    /// failure the whole ordering apparatus exists to prevent, and the one test named for it
    /// never decoded a head.
    #[test]
    fn every_index_entry_points_at_the_block_it_names() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a writable header");
        for record in a_sample() {
            writer.push(&record).expect("the fixture is in order");
        }
        let _ = writer.finish(&[]).expect("it finishes");

        let bytes = bytes_of(&path);
        let tail: [u8; FOOTER_BYTES] = bytes[bytes.len() - FOOTER_BYTES..].try_into().unwrap();
        let footer = decode_footer(&tail).expect("a footer");
        let entries = decode_index(
            &bytes
                [footer.index_offset as usize..(footer.index_offset + footer.index_bytes) as usize],
            footer.n_blocks,
        )
        .expect("the index reads");
        assert!(!entries.is_empty(), "the fixture must produce blocks");

        let mut decompressor = zstd::zstd_safe::DCtx::create();
        for entry in &entries {
            let at = entry.block_offset as usize;
            assert!(
                at < footer.index_offset as usize,
                "a block offset must land in the blocks, not in the index"
            );
            let frame = match crate::ng::psp::block::compressed_block_at(&bytes[at..]) {
                crate::ng::psp::block::CompressedBlockAt::Whole { zstd_frame, .. } => zstd_frame,
                other => panic!("entry at byte {at} does not start a whole block: {other:?}"),
            };
            // `decompress` writes into spare capacity, so the room has to be there first.
            let mut payload = Vec::with_capacity(1 << 20);
            decompressor
                .decompress(&mut payload, frame)
                .expect("a block this writer wrote inflates");
            let head = crate::ng::psp::block::BlockHead::decode(&payload)
                .expect("a block opens with its head")
                .head;
            assert_eq!(
                (head.contig, head.first_position),
                (entry.first_position.contig, entry.first_position.position),
                "the entry at byte {at} names a coordinate the block does not start at"
            );
        }
    }

    // -----------------------------------------------------------------
    // What is refused
    // -----------------------------------------------------------------

    /// A record that goes backwards is refused, and **the file is untouched by the refusal**:
    /// pushing the records that follow still produces a file with exactly the accepted ones.
    #[test]
    fn a_record_that_goes_backwards_is_refused_and_costs_the_file_nothing() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a writable header");
        writer.push(&a_record(0, 1_000, 1)).expect("the first");
        let refused = writer
            .push(&a_record(0, 900, 1))
            .expect_err("that goes backwards");
        assert!(
            matches!(refused, PspWriteError::OutOfOrder { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("SRR7279481.psp"),
            "the message must name the file: {refused}"
        );
        writer.push(&a_record(0, 1_100, 1)).expect("and on we go");
        let stats = writer.finish(&[]).expect("it finishes");
        assert_eq!(stats.records, 2, "the refused record is not in the file");
    }

    /// A contig already written and revisited is refused: two runs of blocks on one contig give
    /// a seek nothing to choose between.
    #[test]
    fn a_contig_that_comes_back_is_refused() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a writable header");
        writer.push(&a_record(1, 10, 1)).expect("contig 1");
        let refused = writer
            .push(&a_record(0, 10, 1))
            .expect_err("contig 0 comes before contig 1");
        // **And it carries the block builder's own account.** The variant has no `reason` of
        // its own: the cause *is* the detail, so a wiring that dropped it would leave a caller
        // with "a record could not be written" and nothing more.
        match refused {
            PspWriteError::RecordRefused { source, .. } => assert!(
                matches!(source, BlockWriteError::ContigOutOfOrder { .. }),
                "got {source}"
            ),
            other => panic!("got {other}"),
        }
    }

    /// **A header this writer cannot honour leaves no file behind at all.** An empty contig list
    /// is refused, and the path it was handed does not exist afterwards — rather than an empty
    /// file that every reader then refuses for a different reason.
    #[test]
    fn a_header_that_cannot_be_written_leaves_no_file() {
        let (_dir, path) = a_file();
        let mut header = a_header(1_000);
        header.contigs.clear();
        let refused = PspWriter::create(&path, header).expect_err("an empty contig list");
        assert!(
            matches!(refused, PspWriteError::InvalidHeaderField { .. }),
            "got {refused}"
        );
        assert!(!path.exists(), "no file must have been created");
    }

    // -----------------------------------------------------------------
    // What an unfinished file is
    // -----------------------------------------------------------------

    /// **A writer dropped without finishing leaves a file with no footer**, which is the whole
    /// of goal 3: a killed run must be refused, not read short.
    ///
    /// The file has a header and blocks in it — it is not empty, which is exactly why reading it
    /// short would be so easy and so wrong.
    #[test]
    fn a_writer_dropped_without_finishing_leaves_a_file_with_no_footer() {
        let (_dir, path) = a_file();
        {
            let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            // Dropped here, deliberately, without `finish`.
        }
        let bytes = bytes_of(&path);
        assert!(
            bytes.len() > FOOTER_BYTES,
            "the file holds a header and blocks: {} bytes",
            bytes.len()
        );
        assert_ne!(
            &bytes[bytes.len() - 4..],
            &FOOTER_MAGIC,
            "an unfinished file must not end with the footer magic"
        );
        // And its header still reads, which is what `read_header` exists for (spec §6.6).
        let header = crate::ng::psp::read_header(&path).expect("the header survives");
        assert_eq!(header.sample, "SRR7279481");
    }

    /// The environment variable that turns this test binary into the writer the test kills. Its
    /// **value is the path the child writes its psp to**.
    ///
    /// **The child is this binary, run again.** A killed writer has to be a real process, and the
    /// only program in this tree that writes a psp is this test suite.
    #[cfg(unix)]
    const WRITE_THE_PSP_HERE: &str = "NG_PSP_KILLED_WRITER_TARGET";

    /// The file the child touches **when it stops pushing records of its own accord**. The
    /// parent requires it to be absent, which is how *the child was still writing when it was
    /// killed* is checked rather than assumed.
    #[cfg(unix)]
    fn where_the_child_says_it_stopped_pushing(psp: &Path) -> PathBuf {
        psp.with_extension("stopped-pushing")
    }

    /// How long the child keeps pushing before it gives up. Far longer than the parent's own
    /// wait, so the kill always arrives first; bounded only so that an orphan left by a parent
    /// that died itself does not write for ever.
    #[cfg(unix)]
    const THE_CHILD_PUSHES_FOR: std::time::Duration = std::time::Duration::from_secs(60);

    /// The name of the test the child is told to run. **Kept beside the test it names**, because
    /// a rename that misses it makes the child run nothing at all.
    #[cfg(unix)]
    const THE_TEST_THE_CHILD_RUNS: &str = "ng::psp::writer::tests::a_writer_killed_before_finishing_leaves_a_file_every_reader_refuses";

    /// How long the parent waits for those bytes before failing rather than hanging, and how
    /// often it looks.
    #[cfg(unix)]
    const WAIT_FOR_THE_CHILD: std::time::Duration = std::time::Duration::from_secs(6);
    #[cfg(unix)]
    const LOOK_EVERY: std::time::Duration = std::time::Duration::from_millis(10);

    /// The signal [`std::process::Child::kill`] sends. Asserted, so an ordinary exit cannot pass
    /// for a kill.
    #[cfg(unix)]
    const SIGKILL: i32 = 9;

    /// **A writer killed with `SIGKILL` before `finish` leaves a file every reader refuses**, and
    /// its header still reads — which is Milestone H2's oracle and spec §6.3's goal 3.
    ///
    /// ⚠ **This is not the same test as
    /// [`a_writer_dropped_without_finishing_leaves_a_file_with_no_footer`].** Dropping a
    /// `PspWriter` in-process runs `BufWriter`'s own `Drop`, which **flushes** — so the file on
    /// disk is everything the writer had produced. `SIGKILL` runs no destructor, so **the blocks
    /// still sitting in the buffer are lost**, and the file is shorter than the same work dropped
    /// would leave.
    ///
    /// ⚠ **What the kill does *not* produce here is a cut inside a block, and the first version
    /// of this doc claimed it did.** `PspWriter` hands each finished block to the `BufWriter` in
    /// one `write_all`, and the buffer is the default 8 kB while this fixture's blocks average
    /// about 57 bytes — so a flush can only ever land where a `write_all` ended, which is a block
    /// boundary. Measured on the file the child leaves: **233 whole blocks, the chain ending
    /// exactly at the file's last byte.** The mid-block state is real and is worth covering, and
    /// what covers it is
    /// [`super::super::reader::tests::every_truncation_of_a_finished_psp_is_refused_without_panicking`],
    /// which manufactures every cut by truncation. **This test's own contribution is the
    /// mechanism end to end on a real process**, plus the `Incomplete` refusal.
    ///
    /// The block-edge property is asserted below rather than left to this paragraph, so that a
    /// change to the writer's buffering says the premise moved instead of quietly making the
    /// claim true again.
    ///
    /// **What is asserted, and why each matters:**
    ///
    /// - the child really was killed by signal 9, so no destructor ran — without this the test
    ///   could be passing on an ordinary exit and proving nothing about the buffer;
    /// - blocks reached disk, so a reader that ignored the missing footer would hand back real
    ///   records from this file — the *sample that covered less of the genome* spec §6.3's goal 3
    ///   forbids;
    /// - `PspReader::open` refuses it as [`PspReadError::Incomplete`], not as damage — the
    ///   instruction is *the run was interrupted, rebuild it*;
    /// - `read_header` still succeeds (spec §6.6), because a tool that reports what a
    ///   half-written file was going to be needs it.
    #[cfg(unix)]
    #[test]
    fn a_writer_killed_before_finishing_leaves_a_file_every_reader_refuses() {
        use std::os::unix::process::ExitStatusExt as _;

        // The child arm. Reached only in the re-executed copy of this binary.
        if let Ok(psp_path) = std::env::var(WRITE_THE_PSP_HERE) {
            let psp_path = PathBuf::from(psp_path);
            let mut writer = PspWriter::create(&psp_path, a_header(1_000)).expect("a header");
            // **It pushes until it is killed, rather than a fixed count**, and that is a
            // correction rather than a preference. The first version pushed 80,000 records and
            // then slept, and the H2 review measured what that did: the whole workload takes
            // about 5 ms while the parent kills at about 12 ms, so in **25 runs out of 25** the
            // child had *finished* and was asleep when the signal arrived. A writer that has
            // stopped writing is not what this test is about, and it made the test a race —
            // adding a `finish` call after the loop failed 40 of 40 standalone runs and passed
            // 3 of 10 under the parallel suite.
            //
            // Four bases apart on the 1 kb grid, so a block closes every 250 records and the
            // parent has a boundary to catch it at whenever it looks.
            // **Bounded by the coordinates as well as by the clock.** Positions are four bases
            // apart and the header's two contigs are about 90 and 53 megabases, so a per-contig
            // cap of twenty million records keeps every position inside its contig. Forty million
            // pushes is roughly 2.5 s of work against a parent that kills at about 12 ms — a
            // margin of some two hundred, where the version this replaces *lost* the race.
            const RECORDS_PER_CONTIG: u64 = 20_000_000;
            let give_up_at = std::time::Instant::now() + THE_CHILD_PUSHES_FOR;
            let mut at = 0u64;
            while std::time::Instant::now() < give_up_at && at < 2 * RECORDS_PER_CONTIG {
                let contig = u32::try_from(at / RECORDS_PER_CONTIG).expect("two contigs");
                let position = 1 + (at % RECORDS_PER_CONTIG) * 4;
                writer
                    .push(&a_record(contig, position, 1))
                    .expect("in order");
                at += 1;
            }
            // **Only reached if the kill never came**, and the marker is what says so. Its
            // presence tells the parent that the file it is about to judge is not the file of an
            // interrupted writer but of one that stopped on its own.
            std::fs::write(
                where_the_child_says_it_stopped_pushing(&psp_path),
                b"stopped",
            )
            .expect("the marker is written");
            // **It panics rather than returning**, because this branch is inside the test it is
            // spawned to re-run: returning would be a *pass*. A process that reaches here was
            // never killed, which means it was not the child the parent spawned — the usual
            // cause being the environment variable already set in the shell.
            panic!(
                "this process took the writer branch and pushed for {THE_CHILD_PUSHES_FOR:?} \
                 without being killed: {WRITE_THE_PSP_HERE} was set in an environment this test \
                 did not create"
            );
        }

        let (_dir, path) = a_file();
        let mut child = std::process::Command::new(
            std::env::current_exe().expect("the test binary knows its own path"),
        )
        .args(["--exact", THE_TEST_THE_CHILD_RUNS, "--nocapture"])
        .env(WRITE_THE_PSP_HERE, &path)
        .stdout(std::process::Stdio::null())
        // **Kept, not discarded.** Every way the child can fail — a stale `--exact` name after a
        // rename, a panic in `create`, a panic in `push` — reaches the parent only as an
        // ordinary exit, and the assertion that fires first is the one about SIGKILL. Without
        // the child's own words that refusal reads as *a destructor ran*, which sends the reader
        // into `BufWriter::drop` instead of at the four characters that went stale.
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("re-executing this test binary");

        // Wait until blocks have reached disk, so the kill lands on a file worth refusing rather
        // than on an empty one. Bounded, so a child that never writes fails the test instead of
        // hanging it.
        //
        // **The threshold is the header's own length, not a multiple of the footer.** An earlier
        // version used eight footers — 384 bytes — against a header of 3,136, so a child killed
        // before it had finished writing its *header* cleared the guard and then failed three
        // assertions later, inside `read_header`, pointing at the header reader.
        let enough_bytes = {
            let (_scratch_dir, sound) = a_finished_psp();
            let (_header, header_bytes) = crate::ng::psp::read_header_and_its_length(&sound)
                .expect("a finished psp's header reads");
            u64::try_from(header_bytes).expect("a header is a small number of bytes")
        };
        let mut bytes_on_disk = 0;
        let mut why_the_file_could_not_be_measured = None;
        let rounds = WAIT_FOR_THE_CHILD.as_millis() / LOOK_EVERY.as_millis();
        for _ in 0..rounds {
            match std::fs::metadata(&path) {
                Ok(found) => bytes_on_disk = found.len(),
                Err(source) => why_the_file_could_not_be_measured = Some(source),
            }
            if bytes_on_disk > enough_bytes {
                break;
            }
            std::thread::sleep(LOOK_EVERY);
        }
        child
            .kill()
            .expect("the child is still alive and can be killed");
        let status = child.wait().expect("the killed child is reaped");
        let mut what_the_child_said = String::new();
        if let Some(mut complaints) = child.stderr.take() {
            use std::io::Read as _;
            let _ = complaints.read_to_string(&mut what_the_child_said);
        }

        // **This assertion cannot tell its causes apart, so it names all of them.** It is the
        // first to fire whenever the child exits on its own, whatever the reason.
        assert_eq!(
            status.signal(),
            Some(SIGKILL),
            "the child exited on its own instead of dying by SIGKILL (exit code {:?}, \
             {bytes_on_disk} bytes written, file not measurable: \
             {why_the_file_could_not_be_measured:?}). Either it never ran this test — \
             `THE_TEST_THE_CHILD_RUNS` repeats this test's own name and goes stale on a rename, \
             and a libtest binary whose filter matches nothing exits 0 — or it ran and exited, \
             which would run the flush this test exists to prevent. The child said: \
             {what_the_child_said}",
            status.code()
        );
        assert!(
            bytes_on_disk > enough_bytes,
            "the child wrote {bytes_on_disk} bytes before the kill; with no blocks on disk this \
             file is refused for being empty rather than for being unfinished"
        );
        // **The child was still pushing when the signal arrived.** That is the difference between
        // interrupting a writer and killing one that had already stopped, and it is the property
        // the first version of this test claimed and did not have.
        assert!(
            !where_the_child_says_it_stopped_pushing(&path).exists(),
            "the child stopped pushing of its own accord before the kill, so the file it left is \
             not an interrupted writer's — it is one that finished its work and waited"
        );

        let bytes = bytes_of(&path);
        assert_ne!(
            bytes.get(bytes.len().saturating_sub(4)..),
            Some(&FOOTER_MAGIC[..]),
            "a killed writer cannot have written a footer"
        );

        // **Where the kill actually cut.** Walk the length-prefixed block chain forward from the
        // header's end: if it lands exactly on the file's last byte, every block on disk is
        // whole, and the kill lost buffered blocks rather than splitting one.
        let (_, header_bytes) = crate::ng::psp::read_header_and_its_length(&path)
            .expect("the header survives the kill");
        let mut at = header_bytes;
        let mut whole_blocks = 0;
        while at + crate::ng::psp::COMPRESSED_BLOCK_LENGTH_BYTES <= bytes.len() {
            let declared: [u8; crate::ng::psp::COMPRESSED_BLOCK_LENGTH_BYTES] = bytes
                [at..at + crate::ng::psp::COMPRESSED_BLOCK_LENGTH_BYTES]
                .try_into()
                .expect("a slice of exactly the prefix's width");
            let next = at
                + crate::ng::psp::COMPRESSED_BLOCK_LENGTH_BYTES
                + u32::from_le_bytes(declared) as usize;
            if next > bytes.len() {
                break;
            }
            at = next;
            whole_blocks += 1;
        }
        assert!(
            whole_blocks > 0,
            "the kill has to leave whole blocks on disk, or this file is a header and nothing else"
        );
        assert_eq!(
            at,
            bytes.len(),
            "the kill left a partial block on disk. That is not a defect — it is the state this \
             test's doc used to claim it produced and does not. If the writer's buffering has \
             changed so that a block can now exceed the `BufWriter`, update the ⚠ paragraph \
             above: {whole_blocks} whole blocks, chain ending at {at} of {} bytes",
            bytes.len()
        );
        let refused = PspReader::open(&path).expect_err("a file with no footer is not readable");
        assert!(
            matches!(refused, crate::ng::psp::PspReadError::Incomplete { .. }),
            "the instruction is `the run was interrupted, rebuild it` and not `the file is \
             damaged`; got {refused}"
        );
        // Spec §6.6: the half-written file still has a header, and this is what reads it.
        let header = crate::ng::psp::read_header(&path).expect("the header survives the kill");
        assert_eq!(header.sample, "SRR7279481");
    }

    /// A file with no records at all still finishes, and its index is empty rather than absent.
    #[test]
    fn a_sample_with_no_records_still_finishes() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let stats = writer.finish(b"nothing was found").expect("it finishes");
        assert_eq!(stats.records, 0);
        assert_eq!(stats.blocks, 0);

        let bytes = bytes_of(&path);
        let tail: [u8; FOOTER_BYTES] = bytes[bytes.len() - FOOTER_BYTES..].try_into().unwrap();
        let footer = decode_footer(&tail).expect("an empty psp is still a psp");
        assert_eq!(footer.n_blocks, 0);
        assert_eq!(footer.index_bytes, 0);
        assert_eq!(footer.trailer_bytes, b"nothing was found".len() as u64);
    }

    /// An empty trailer is legal, and the footer says so.
    #[test]
    fn an_empty_trailer_is_legal() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        writer.push(&a_record(0, 1, 1)).expect("one record");
        let _ = writer.finish(&[]).expect("it finishes");

        let bytes = bytes_of(&path);
        let tail: [u8; FOOTER_BYTES] = bytes[bytes.len() - FOOTER_BYTES..].try_into().unwrap();
        let footer = decode_footer(&tail).expect("a footer");
        assert_eq!(footer.trailer_bytes, 0);
        assert_eq!(
            footer.trailer_offset,
            footer.index_offset + footer.index_bytes
        );
    }

    /// **A failed flush is surfaced, not swallowed** — the durability contract spec §6.3 calls
    /// easy to get wrong.
    ///
    /// `/dev/full` accepts an open and fails every write with `ENOSPC`. The fixture is well under
    /// `BufWriter`'s 8 KiB buffer, so nothing reaches the device until `into_inner` flushes
    /// inside `finish` — which pins exactly the step that surfaces.
    ///
    /// ⚠ **F3 recorded this as untestable without a failing file descriptor and routed it to
    /// H2.** It was testable here; the review found the route. Linux only — macOS has no
    /// `/dev/full`, and the container this project builds in is Linux.
    #[cfg(target_os = "linux")]
    #[test]
    fn a_failed_flush_is_surfaced_rather_than_swallowed() {
        let full = Path::new("/dev/full");
        if !full.exists() {
            return;
        }
        let mut writer = PspWriter::create(full, a_header(1_000))
            .expect("/dev/full accepts an open and buffers the header");
        writer.push(&a_record(0, 1, 1)).expect("it buffers");
        let refused = writer
            .finish(b"a per-sample summary")
            .expect_err("a device with no room must not report a finished file");
        assert!(matches!(refused, PspWriteError::Io { .. }), "got {refused}");
        assert!(
            refused.to_string().contains("flushing the finished file"),
            "the message must name the step that failed: {refused}"
        );
    }

    // -----------------------------------------------------------------
    // The guard that reads back what it is about to write
    // -----------------------------------------------------------------

    /// **The guard fires on an index the reader would refuse**, driven directly because the
    /// writer cannot currently produce one: `BlockBuilder` enforces coordinate order, so every
    /// index `finish` builds is already sorted.
    ///
    /// That is exactly why the guard is worth having and worth testing this way — it is a net
    /// under a defect somewhere else, and a net nothing ever falls into is untested by
    /// construction. Removing the call changed no test until this one existed.
    #[test]
    fn the_readable_check_refuses_an_index_a_reader_would_not_take() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        writer.push(&a_record(0, 1, 1)).expect("one record");

        // Two entries whose offsets go backwards — the shape `decode_index` refuses.
        let entry = |position: u64, block_offset: u64| BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(0),
                position: Position(position),
            },
            block_offset,
        };
        writer.index = vec![entry(1, 8_192), entry(2, 4_096)];

        let index_bytes = index::encode_index(&writer.index);
        let footer = Footer {
            index_offset: 100,
            index_bytes: index_bytes.len() as u64,
            trailer_offset: 100 + index_bytes.len() as u64,
            trailer_bytes: 0,
            n_blocks: writer.index.len() as u64,
            index_checksum: index::checksum_index(&index_bytes),
        };
        let refused = writer
            .check_the_index_and_footer_read_back(&index_bytes, &footer)
            .expect_err("a reader would refuse that index");
        // **The index decoder's own account is kept.** Replacing the cause with `None` left the
        // whole suite green until this looked at it (the G3 review).
        match &refused {
            PspWriteError::WouldNotBeReadable {
                source: Some(NotReadable::BlockIndex(_)),
                reason,
                ..
            } => assert_eq!(reason, "the block index it wrote does not decode"),
            other => panic!("got {other}"),
        }

        // **And through `finish`, which is what wires the guard in.** Called directly the guard
        // could be right while nothing called it: deleting the call from `finish` left the
        // assertion above green.
        let through_finish = writer
            .finish(&[])
            .expect_err("finish must not write a file it cannot read");
        assert!(
            matches!(through_finish, PspWriteError::WouldNotBeReadable { .. }),
            "got {through_finish}"
        );
        assert!(
            refused.to_string().contains("SRR7279481.psp"),
            "the message must name the file: {refused}"
        );
        assert!(
            refused.to_string().contains("block index"),
            "the message must say which structure: {refused}"
        );
    }

    /// The same guard on the footer half: sections that do not abut.
    #[test]
    fn the_readable_check_refuses_a_footer_a_reader_would_not_take() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let index_bytes = index::encode_index(&[]);
        let footer = Footer {
            index_offset: 100,
            index_bytes: 0,
            // A gap: the trailer does not start where the index ends.
            trailer_offset: 200,
            trailer_bytes: 0,
            n_blocks: 0,
            index_checksum: index::checksum_index(&index_bytes),
        };
        let refused = writer
            .check_the_index_and_footer_read_back(&index_bytes, &footer)
            .expect_err("a reader would refuse that footer");
        assert!(
            refused.to_string().contains("footer"),
            "the message must say which structure: {refused}"
        );
    }

    /// **A writer that lost a block refuses to finish**, rather than writing a file every reader
    /// accepts with records missing out of the middle of it.
    ///
    /// The loss is forced the way the review forced it: the writer's own `spent` marker is set,
    /// standing in for a compressor refusal, a head that would not read back, or a transient
    /// write error. What matters is what `finish` then does.
    ///
    /// ⚠ Before this, `finish` returned `Ok` and produced a file with a valid footer, a
    /// checksum-matching index, entries in ascending order — and a thousand bases missing from
    /// the middle of a contig. That is worse than the unreadable stump a killed run leaves,
    /// because every reader takes it.
    #[test]
    fn a_writer_that_lost_a_block_refuses_to_finish() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        for record in a_sample() {
            writer.push(&record).expect("in order");
        }
        writer.spent = Some("a block could not be written");

        let refused = writer
            .finish(b"a per-sample summary")
            .expect_err("a walk that lost records must not report a finished file");
        // **No cause here, and that is right**: nothing decoded, the loss happened earlier and
        // its own error went to the caller then.
        assert!(
            matches!(
                refused,
                PspWriteError::WouldNotBeReadable { source: None, .. }
            ),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("records were lost"),
            "the message must say what happened: {refused}"
        );

        // And the file on disk is not a finished psp: no footer, so every reader refuses it.
        let bytes = bytes_of(&path);
        assert_ne!(&bytes[bytes.len() - 4..], &FOOTER_MAGIC);
    }

    /// **The guard runs before anything is written, and that ordering is the point.**
    ///
    /// Moved after the three writes it still returns `Err` — and leaves a file ending in the
    /// footer magic on disk, which by this module's own contract is a *finished* file. A
    /// `finish` that fails must leave nothing a reader will take.
    #[test]
    fn a_finish_that_refuses_leaves_no_finished_file_behind() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        writer.push(&a_record(0, 1, 1)).expect("one record");

        // An index the reader refuses: offsets that go backwards.
        let entry = |position: u64, block_offset: u64| BlockIndexEntry {
            first_position: GenomePosition {
                contig: ContigId(0),
                position: Position(position),
            },
            block_offset,
        };
        writer.index = vec![entry(1, 8_192), entry(2, 4_096)];

        let refused = writer
            .finish(&[])
            .expect_err("a reader would refuse that index");
        assert!(
            matches!(refused, PspWriteError::WouldNotBeReadable { .. }),
            "got {refused}"
        );
        let bytes = bytes_of(&path);
        assert!(
            bytes.len() < FOOTER_BYTES || bytes[bytes.len() - 4..] != FOOTER_MAGIC,
            "a refused finish must not leave a file ending in the footer magic"
        );
    }

    /// The footer's checksum must be over the index the writer actually wrote.
    ///
    /// Computing it over the trailer instead sailed through the guard until this was asserted —
    /// a file `open` refuses at the checksum, produced by a `finish` that returned `Ok`.
    #[test]
    fn the_readable_check_refuses_a_checksum_over_the_wrong_bytes() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let index_bytes = index::encode_index(&[]);
        let footer = Footer {
            index_offset: 100,
            index_bytes: 0,
            trailer_offset: 100,
            trailer_bytes: 0,
            n_blocks: 0,
            index_checksum: index::checksum_index(b"something else entirely"),
        };
        let refused = writer
            .check_the_index_and_footer_read_back(&index_bytes, &footer)
            .expect_err("that checksum is of other bytes");
        assert!(refused.to_string().contains("checksum"), "got {refused}");
    }

    /// The compression level the writer used reaches the file, because it is a setting and
    /// goal 4 is that settings are recorded — and because `append` must match bytes already
    /// written and has no other way to learn what level produced them.
    #[test]
    fn the_compression_level_reaches_the_header() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let _ = writer.finish(&[]).expect("it finishes");

        let header = crate::ng::psp::read_header(&path).expect("the header reads");
        assert_eq!(
            header.writer.parameters.get("zstd-compression-level"),
            Some(&crate::ng::psp::header::ParameterValue::Integer(i64::from(
                crate::ng::psp::block::ZSTD_COMPRESSION_LEVEL
            ))),
            "the level this build compresses at must be in the file"
        );
    }

    // -----------------------------------------------------------------
    // The same sample twice
    // -----------------------------------------------------------------

    /// **The same records give the same bytes**, which is what spec §7's worker-count invariance
    /// rests on — the header's timestamp is the one field allowed to differ, and this fixture
    /// holds it fixed so the rest can be compared byte for byte.
    #[test]
    fn the_same_records_written_twice_give_the_same_bytes() {
        let write = || {
            let (dir, path) = a_file();
            let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            let _ = writer.finish(b"the same payload").expect("it finishes");
            let bytes = bytes_of(&path);
            drop(dir);
            bytes
        };
        assert_eq!(write(), write());
    }

    /// The block grid decides how many blocks a sample cuts, and the index follows it. **A
    /// coarser grid gives fewer blocks and a shorter index**, which is the whole memory argument
    /// of the format (spec §3.3).
    #[test]
    fn a_coarser_block_grid_gives_fewer_index_entries() {
        let blocks_at = |grid: u64| {
            let (_dir, path) = a_file();
            let mut writer = PspWriter::create(&path, a_header(grid)).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            writer.finish(&[]).expect("it finishes").blocks
        };
        let fine = blocks_at(1_000);
        let coarse = blocks_at(1_000_000);
        assert!(
            coarse < fine,
            "a 1 Mb grid must cut fewer blocks than a 1 kb one: {coarse} against {fine}"
        );
        // Two contigs, and a block never crosses one, so the floor is two.
        assert_eq!(coarse, 2);
    }

    // -----------------------------------------------------------------
    // Appending to a finished file (G4)
    // -----------------------------------------------------------------

    /// **A file appended to reads back as one file**: every record that was there and every
    /// record added, in order, from one walk.
    ///
    /// The comparison is against the two record lists concatenated, field for field — an append
    /// that lost the old records, or wrote the new ones somewhere a walk does not reach, fails
    /// here rather than in a count.
    #[test]
    fn an_appended_file_reads_back_as_one_file() {
        let (_dir, path) = a_finished_psp();
        let already_there = a_sample();
        let added: Vec<_> = (0..6u64)
            .map(|step| a_record(1, 60_000 + step * 100, 1))
            .collect();

        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        for record in &added {
            writer
                .push(record)
                .expect("past the last record already there");
        }
        let stats = writer
            .finish(b"a summary of both halves")
            .expect("it finishes");
        // **The three fields are measured on two different populations**, and each is pinned:
        // records is this writer's, blocks and bytes are the file's.
        assert_eq!(
            stats.records, 6,
            "what this writer wrote, not what the file holds"
        );
        assert_eq!(
            stats.blocks, 9,
            "eight already there and one block appended"
        );
        assert_eq!(stats.bytes, bytes_of(&path).len() as u64, "the whole file");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let read: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .map(|found| found.expect("a finished psp walks").record.expect("a body"))
            .collect();
        let expected: Vec<_> = already_there.iter().chain(&added).cloned().collect();
        assert_eq!(read, expected);
        assert_eq!(
            psp.trailer().expect("it reads"),
            b"a summary of both halves"
        );
    }

    /// **The blocks and the header that were already there are byte-identical afterwards**, and
    /// the appended blocks start exactly where the old index did.
    #[test]
    fn appending_leaves_the_header_and_the_old_blocks_untouched() {
        let (_dir, path) = a_finished_psp();
        let before = bytes_of(&path);
        let footer = footer_of(&before);
        let blocks_end = footer.index_offset as usize;
        let old_blocks = before[..blocks_end].to_vec();
        let old_entries = {
            let psp = PspReader::open(&path).expect("it opens");
            psp.block_index().to_vec()
        };

        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        writer.push(&a_record(1, 60_000, 1)).expect("in order");
        let _ = writer.finish(b"").expect("it finishes");

        let after = bytes_of(&path);
        assert_eq!(
            after[..blocks_end],
            old_blocks[..],
            "the header and every old block"
        );
        let psp = PspReader::open(&path).expect("it opens");
        assert_eq!(
            &psp.block_index()[..old_entries.len()],
            &old_entries[..],
            "the old index entries come back unchanged"
        );
        assert_eq!(
            psp.block_index()[old_entries.len()].block_offset,
            blocks_end as u64,
            "and the first appended block starts where the old index did"
        );
    }

    /// **Coordinate order runs across the seam** (spec §6.4): a first appended record that
    /// precedes the last one already in the file is refused, on the position and on the contig.
    ///
    /// ⚠ **A builder that started blank would accept both** — it is the seed taken from the last
    /// block's heads that makes this fail, and nothing else in the file could.
    #[test]
    fn a_record_that_precedes_the_seam_is_refused() {
        for (contig, at, what) in [
            (1u32, 1u64, "an earlier position on the same contig"),
            (0, 9_000, "an earlier contig"),
        ] {
            let (_dir, path) = a_finished_psp();
            let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
            let refused = writer
                .push(&a_record(contig, at, 1))
                .expect_err("that is behind the seam");
            assert!(
                matches!(
                    refused,
                    PspWriteError::OutOfOrder { .. } | PspWriteError::RecordRefused { .. }
                ),
                "{what}: got {refused}"
            );
        }
    }

    /// **A record at exactly the last record's position is accepted**, because the rule is *must
    /// not precede* and two records may begin on one base.
    #[test]
    fn a_record_at_the_seams_own_position_is_accepted() {
        let (_dir, path) = a_finished_psp();
        let last = {
            let mut psp = PspReader::open(&path).expect("it opens");
            psp.records()
                .expect("the walk starts")
                .last()
                .expect("the fixture has records")
                .expect("a finished psp walks")
                .head
                .region
        };
        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        writer
            .push(&a_record(last.contig.0, last.start.get(), 1))
            .expect("the same base is not before it");
        let _ = writer.finish(b"").expect("it finishes");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(psp.records().expect("the walk starts").count(), 41);
    }

    /// **A file with no footer cannot be appended to**, because nothing says where its blocks
    /// end — it is what a killed run leaves, and the answer is to re-run it.
    #[test]
    fn a_file_with_no_footer_cannot_be_appended_to() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        writer.push(&a_record(0, 1, 1)).expect("one record");
        drop(writer);

        let refused = PspWriter::append(&path).expect_err("there is no footer");
        // ⚠ **`Reopen { .. }` alone cannot tell what this test names from three other things** —
        // it is satisfied by `NotAnNgPsp`, by `Damaged` and by `Io` too. That is G3's M2, in a
        // test written after that finding was fixed.
        assert!(
            matches!(
                refused,
                PspWriteError::Reopen {
                    source: crate::ng::psp::PspReadError::Incomplete { .. },
                    ..
                }
            ),
            "got {refused}"
        );
    }

    /// **An append is refused on a manifest this writer cannot honour** (spec §6.4), because the
    /// header is not rewritten: the added records would have to use encodings the file declares
    /// and this build does not write, and left unchecked they would be added anyway — unreadable
    /// under the header the file keeps.
    #[test]
    fn an_append_is_refused_on_a_manifest_this_writer_cannot_honour() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let declared = b"position-offset";
        let at = whole
            .windows(declared.len())
            .position(|window| window == declared)
            .expect("the manifest names the record's first field");
        whole[at..at + declared.len()].copy_from_slice(b"position-offsex");
        rewrite(&path, &whole);

        let refused = PspWriter::append(&path).expect_err("this build cannot write that layout");
        match refused {
            PspWriteError::UnsupportedManifest { source, .. } => {
                assert!(
                    matches!(source, crate::ng::psp::ManifestRefusal::CutRule(_)),
                    "got {source}"
                )
            }
            other => panic!("got {other}"),
        }
        assert_eq!(bytes_of(&path), whole, "and the file is exactly as it was");
    }

    /// **An append dropped before it finishes leaves a file every reader refuses** — and the
    /// trailer the file had is gone with the index it followed (spec §6.4).
    #[test]
    fn an_append_dropped_before_finishing_leaves_a_file_no_reader_accepts() {
        let (_dir, path) = a_finished_psp();
        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        writer.push(&a_record(1, 60_000, 1)).expect("in order");
        drop(writer);

        let refused = PspReader::open(&path).expect_err("there is no footer any more");
        assert!(
            matches!(refused, crate::ng::psp::PspReadError::Incomplete { .. }),
            "got {refused}"
        );
    }

    /// **An append writes at the level the file records**, not at this build's — so a file's
    /// blocks are all compressed the same way, which is goal 4.
    ///
    /// ⚠ **The fixture has to make the two levels give different bytes, and the first one did
    /// not.** It compared an appended block against a fresh block at the level the file records,
    /// which *is* this build's — so an append that ignored the record entirely passed. Here the
    /// file's recorded level is patched to 1 in its own header text, which is one character wide
    /// and moves no offset, and the payload is varied enough that 1 and 9 differ. **Both halves
    /// are asserted**: that the levels separate, and that the appended block is the recorded
    /// level's.
    #[test]
    fn an_append_compresses_at_the_level_the_file_records() {
        use crate::ng::psp::header::ParameterValue;
        use crate::ng::psp::{
            BlockCompressor, COMPRESSED_BLOCK_LENGTH_BYTES, ZSTD_COMPRESSION_LEVEL,
        };

        // A record whose bases do not repeat trivially, so the levels have something to disagree
        // about.
        let varied = || {
            let mut record = a_record(1, 60_000, 1);
            let bases: Vec<u8> = (0..8_000u64)
                .map(|i| b"ACGT"[((i * 7 + i / 3 + i * i / 5) % 4) as usize])
                .collect();
            record.reference_bases = bases.clone().into_boxed_slice();
            record.observations[0].bases = bases.into_boxed_slice();
            record
        };

        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        assert_eq!(
            ZSTD_COMPRESSION_LEVEL, 9,
            "or the fixture patches the wrong text"
        );
        let recorded = b"zstd-compression-level = 9";
        let patched_to = b"zstd-compression-level = 1";
        let at = whole
            .windows(recorded.len())
            .position(|window| window == recorded)
            .expect("the header records the level it was written at");
        // **The same width, so no offset in the file moves** — the whole reason for a one-digit
        // level here.
        whole[at..at + patched_to.len()].copy_from_slice(patched_to);
        rewrite(&path, &whole);
        {
            let psp = PspReader::open(&path).expect("the patched file still opens");
            assert_eq!(
                psp.header().writer.parameters.get("zstd-compression-level"),
                Some(&ParameterValue::Integer(1)),
                "the fixture must record a level this build does not use"
            );
        }

        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        writer.push(&varied()).expect("in order");
        let _ = writer.finish(b"").expect("it finishes");
        let appended = {
            let whole = bytes_of(&path);
            let psp = PspReader::open(&path).expect("it opens");
            let at = psp.block_index().last().expect("blocks").block_offset as usize;
            whole[at..psp.footer().index_offset as usize].to_vec()
        };

        // The same record, alone, in a fresh file — which this build writes at its own level.
        let (_dir2, fresh) = a_file();
        let mut writer = PspWriter::create(&fresh, a_header(1_000)).expect("a header");
        writer.push(&varied()).expect("one record");
        let _ = writer.finish(b"").expect("it finishes");
        let (payload, at_this_builds_level) = {
            let whole = bytes_of(&fresh);
            let psp = PspReader::open(&fresh).expect("it opens");
            let at = psp.block_index()[0].block_offset as usize;
            let block = whole[at..psp.footer().index_offset as usize].to_vec();
            let payload = zstd::decode_all(&block[COMPRESSED_BLOCK_LENGTH_BYTES..])
                .expect("the block inflates");
            (payload, block)
        };

        let mut at_one =
            BlockCompressor::with_level(a_header(1_000).manifest.look_back_window_log, 1)
                .expect("a level zstd takes");
        let at_one = at_one.compress(&payload).expect("it compresses").to_vec();
        assert_ne!(
            at_one, at_this_builds_level,
            "the fixture must separate the levels, or it proves nothing"
        );
        assert_eq!(
            appended, at_one,
            "the appended block is compressed at the level the file records"
        );
    }

    /// **A footer that puts the block index inside the header is refused, and the file is left
    /// whole.**
    ///
    /// ⚠ **This is the G4 review's Blocker, and it is G3's one operation over.** `append`
    /// truncates at `index_offset`, so this rule is all that stands between a footer of nonsense
    /// and a psp reduced to a stump — and **the per-entry rule cannot carry it: on an empty index
    /// there are no entries to check.** A footer saying the index is at byte 4 and holds nothing
    /// passed `PspReader::open`, and `append` then cut a 3,742-byte psp down to 109 bytes and
    /// returned `Ok`.
    #[test]
    fn a_footer_that_puts_the_block_index_inside_the_header_is_refused() {
        use crate::ng::psp::footer::Footer;

        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let mut crafted = whole[..whole.len() - FOOTER_BYTES].to_vec();
        let footer = Footer {
            index_offset: 4,
            index_bytes: 0,
            trailer_offset: 4,
            trailer_bytes: (crafted.len() - 4) as u64,
            n_blocks: 0,
            index_checksum: crate::ng::psp::checksum_index(&[]),
        };
        crafted.extend_from_slice(&crate::ng::psp::encode_footer(&footer));
        rewrite(&path, &crafted);

        let refused = PspWriter::append(&path).expect_err("the index is inside the header");
        match refused {
            PspWriteError::Reopen {
                source: crate::ng::psp::PspReadError::Damaged { reason, .. },
                ..
            } => assert!(reason.contains("inside the"), "got {reason}"),
            other => panic!("got {other}"),
        }
        assert_eq!(bytes_of(&path), crafted, "and every byte is still there");
    }

    /// **The seam is the last record in the file, not the last block's first.**
    ///
    /// ⚠ **The first version of the seam test could not tell those apart.** Both records it
    /// offered preceded every record in the last block, so it separated *seeded* from *not
    /// seeded* and nothing finer — and an implementation keeping the **first** record of the last
    /// block passed all 381 tests. A record falling between the two is what closes it.
    #[test]
    fn a_record_inside_the_last_blocks_span_is_refused() {
        let (_dir, path) = a_finished_psp();
        let (first_of_last, last_of_last) = {
            let mut psp = PspReader::open(&path).expect("it opens");
            let n = psp.block_index().len();
            let regions: Vec<_> = psp
                .records_from_block(n - 1)
                .expect("the walk starts")
                .building_only_where(|_| false)
                .map(|found| found.expect("a finished psp walks").head.region)
                .collect();
            (regions[0], *regions.last().expect("the block has records"))
        };
        assert!(
            first_of_last.start < last_of_last.start,
            "the fixture's last block must hold more than one position, or this proves nothing"
        );

        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        let refused = writer
            .push(&a_record(
                last_of_last.contig.0,
                first_of_last.start.get() + 1,
                1,
            ))
            .expect_err("that is behind the last record, though inside the last block");
        assert!(
            matches!(refused, PspWriteError::OutOfOrder { .. }),
            "got {refused}"
        );
    }

    /// **A psp that holds no records is appendable**, and the first appended record is checked
    /// against nothing: `finish` writes such a file, `open` accepts it, and there is no block to
    /// walk for a seam.
    #[test]
    fn a_psp_that_holds_no_records_is_appendable() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let empty = writer.finish(b"no coverage").expect("it finishes");
        assert_eq!(empty.blocks, 0, "the fixture must hold no blocks");

        let mut writer = PspWriter::append(&path).expect("a psp with no blocks is appendable");
        // The seam is nothing, so even the file's lowest coordinate is legal.
        writer
            .push(&a_record(0, 1, 1))
            .expect("nothing precedes it");
        let stats = writer.finish(b"one record now").expect("it finishes");
        assert_eq!((stats.records, stats.blocks), (1, 1));

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(psp.records().expect("the walk starts").count(), 1);
        assert_eq!(psp.trailer().expect("it reads"), b"one record now");
    }

    /// **A recorded compression level this build cannot use is refused, not written at another
    /// one** — and the refusal names the level the file records, not a substitute for it.
    ///
    /// ⚠ **It used to name a substitute.** Every level outside an `i32` became `i32::MAX` before
    /// the check saw it, so a file recording −5,000,000,000 was refused for recording
    /// 2,147,483,647 — a number that is not in the file. And nothing reached the path at all.
    #[test]
    fn an_append_is_refused_on_a_recorded_level_zstd_will_not_take() {
        use crate::ng::psp::header::ParameterValue;

        for recorded in [99i64, -5_000_000_000, 5_000_000_000] {
            let mut header = a_header(1_000);
            header.writer.parameters.insert(
                "zstd-compression-level".to_string(),
                ParameterValue::Integer(recorded),
            );
            let refused = PspWriter::build_the_compressor_the_header_records(&header)
                .expect_err("zstd does not take that level");
            assert!(
                matches!(
                    refused,
                    crate::ng::psp::ManifestRefusal::Compressor(_)
                        | crate::ng::psp::ManifestRefusal::LevelPastAnyLevel { .. }
                ),
                "got {refused}"
            );
            // **The refusal names the level the file records**, not a substitute for it — which
            // is the half that used to be wrong: every level outside an `i32` became `i32::MAX`
            // before the range check saw it.
            assert!(
                refused.to_string().contains(&recorded.to_string()),
                "the refusal must name the level the file records; got {refused}"
            );
        }
    }

    /// **A level recorded in a shape this writer cannot read is refused, not ignored.**
    ///
    /// A setting recorded and then silently replaced is what §2.4 argues is worse than one not
    /// recorded — and *absent* is the one case that legitimately falls back, because it is what a
    /// file written before the parameter existed looks like.
    #[test]
    fn a_level_recorded_in_a_shape_this_writer_cannot_read_is_refused() {
        use crate::ng::psp::header::ParameterValue;

        for shape in [
            ParameterValue::String("1".to_string()),
            ParameterValue::Float(1.0),
            ParameterValue::Boolean(true),
        ] {
            let mut header = a_header(1_000);
            header
                .writer
                .parameters
                .insert("zstd-compression-level".to_string(), shape.clone());
            let refused = PspWriter::build_the_compressor_the_header_records(&header)
                .expect_err("a level recorded and unreadable is not a level to guess at");
            assert!(
                matches!(
                    refused,
                    crate::ng::psp::ManifestRefusal::UnreadableLevel { .. }
                ),
                "for {shape:?}, got {refused}"
            );
        }

        // And absent is the compatibility case, which still falls back.
        let mut header = a_header(1_000);
        header.writer.parameters.remove("zstd-compression-level");
        PspWriter::build_the_compressor_the_header_records(&header)
            .expect("a file that records none");
    }

    /// **The seam is found from the last block alone.** `records_from_block` walks to the end of
    /// the file, so starting anywhere earlier gives the same answer at the cost of every block
    /// before it — a change nothing else would notice.
    #[test]
    fn the_seam_is_found_by_reading_only_the_last_block() {
        let (_dir, path) = a_finished_psp();
        let mut psp = PspReader::open(&path).expect("it opens");
        let blocks = psp.block_index().len();
        let seen = psp
            .records_from_block(blocks - 1)
            .expect("the walk starts")
            .building_only_where(|_| false)
            .count();
        assert_eq!(seen, 5, "the last block's records, not the file's");
        assert!(
            seen < 40,
            "or the walk is reading blocks the seam does not need"
        );
    }

    /// **An append stopped part way leaves a file no reader accepts**, at every stopping point
    /// from the truncation onwards — the shape `replace_trailer` was given at G3, on the more
    /// destructive of the two operations.
    #[test]
    fn an_append_stopped_part_way_leaves_a_file_no_reader_accepts() {
        let (_dir, path) = a_finished_psp();
        let blocks_end = footer_of(&bytes_of(&path)).index_offset as usize;
        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        for step in 0..6u64 {
            writer
                .push(&a_record(1, 60_000 + step * 100, 1))
                .expect("in order");
        }
        let _ = writer.finish(b"whole").expect("it finishes");
        let whole = bytes_of(&path);

        for stopped_at in blocks_end..whole.len() {
            rewrite(&path, &whole[..stopped_at]);
            assert!(
                PspReader::open(&path).is_err(),
                "a file stopped at byte {stopped_at} of {} opened",
                whole.len()
            );
        }
        rewrite(&path, &whole);
        assert!(PspReader::open(&path).is_ok(), "and the complete one opens");
    }

    /// **A second append inherits the first's index, not the original's** — and the file grows a
    /// block a time, which §6's trade-off note is about.
    #[test]
    fn a_second_append_inherits_the_first_ones_index() {
        let (_dir, path) = a_finished_psp();
        for round in 0..2u64 {
            let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
            writer
                .push(&a_record(1, 60_000 + round * 1_000, 1))
                .expect("in order");
            let stats = writer.finish(b"round").expect("it finishes");
            assert_eq!(stats.records, 1, "what this writer wrote");
            assert_eq!(stats.blocks, 9 + round, "and every block in the file");
        }
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(psp.records().expect("the walk starts").count(), 42);
    }

    /// **A psp whose last block cannot inflate is not extended**, and the file is left whole.
    ///
    /// This is the one check `append` makes that `open` does not: the seam walk reads the last
    /// block, so a damaged one is met before a byte is truncated. §2.1's "every check `open`
    /// makes comes free" is therefore only half of it.
    #[test]
    fn a_psp_whose_last_block_cannot_inflate_is_not_extended() {
        let (_dir, path) = a_finished_psp();
        let blocks = {
            let psp = PspReader::open(&path).expect("it opens");
            psp.block_index().len()
        };
        wreck_the_block(&path, blocks - 1);
        let wrecked = bytes_of(&path);

        let refused = PspWriter::append(&path).expect_err("the last block does not inflate");
        assert!(
            matches!(refused, PspWriteError::Reopen { .. }),
            "got {refused}"
        );
        assert_eq!(
            bytes_of(&path),
            wrecked,
            "and the file is exactly as it was"
        );
    }

    /// **A record appended into the grid cell the old last block covered opens a second block
    /// for that cell**, rather than reopening the first — which is legal, and is one of the two
    /// ways two blocks come to share a first position (`index.rs`).
    ///
    /// ⚠ **`continuing_after`'s doc claims both halves and nothing asserted either.** It is the
    /// shape G1's Blocker was about, so it is worth pinning from the writer's side too.
    #[test]
    fn a_record_appended_into_the_old_last_cell_opens_a_second_block_for_it() {
        let (_dir, path) = a_finished_psp();
        let (blocks_before, last_entry) = {
            let psp = PspReader::open(&path).expect("it opens");
            (
                psp.block_index().len(),
                *psp.block_index().last().expect("blocks"),
            )
        };
        // The fixture's grid is 1 kb and its last block starts at 3,001; 3,500 is the same cell.
        let same_cell = last_entry.first_position.position.get() + 499;
        assert_eq!(
            same_cell / 1_000,
            last_entry.first_position.position.get() / 1_000
        );

        let mut writer = PspWriter::append(&path).expect("a finished psp is appendable");
        writer
            .push(&a_record(last_entry.first_position.contig.0, same_cell, 1))
            .expect("past the seam");
        let _ = writer.finish(b"").expect("it finishes");

        let psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(
            psp.block_index().len(),
            blocks_before + 1,
            "the appended record opened a block of its own"
        );
        assert_eq!(
            psp.block_index()
                .last()
                .expect("blocks")
                .first_position
                .position
                .get(),
            same_cell
        );
    }
}
