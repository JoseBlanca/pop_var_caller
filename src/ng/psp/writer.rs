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
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use super::block::{BlockBuilder, BlockCompressor, BlockHead, BlockWriteError};
use super::footer::{Footer, encode_footer};
use super::header::Header;
use super::index::BlockIndexEntry;
use super::{PspWriteError, footer, index};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::GenomePosition;

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
                reason: source.to_string(),
            }
        })?;
        let compressor = BlockCompressor::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedManifest {
                path: path.to_path_buf(),
                reason: source.to_string(),
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
                    reason: other.to_string(),
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
        let block = match self.compressor.compress(payload) {
            Ok(block) => block.to_vec(),
            Err(source) => {
                self.spent = Some("a block could not be compressed");
                return Err(PspWriteError::BlockRefused {
                    path: self.path.clone(),
                    reason: source.to_string(),
                });
            }
        };
        self.put_block(&head, &block)
            .inspect_err(|_| self.spent = Some("a block could not be written"))
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
                    reason: source.to_string(),
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
                reason: format!("the block index it wrote reads back as: {source}"),
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
            });
        }
        footer::decode_footer(&encode_footer(footer)).map_err(|source| {
            PspWriteError::WouldNotBeReadable {
                path: self.path.clone(),
                reason: format!("the footer it wrote reads back as: {source}"),
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
                reason: format!("a block it just built reads back as: {source}"),
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
    use super::tests_support::{a_file, a_header, a_record, a_sample, bytes_of};
    use super::*;
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
        assert!(
            matches!(refused, PspWriteError::RecordRefused { .. }),
            "got {refused}"
        );
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
        assert!(
            matches!(refused, PspWriteError::WouldNotBeReadable { .. }),
            "got {refused}"
        );

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
        assert!(
            matches!(refused, PspWriteError::WouldNotBeReadable { .. }),
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
}
