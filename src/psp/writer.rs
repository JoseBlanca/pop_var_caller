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

use super::block::{BlockBuilder, BlockCompressError, BlockCompressor, BlockHead, BlockWriteError};
use super::footer::{FOOTER_BYTES, Footer, encode_footer};
use super::header::Header;
use super::index::BlockIndexEntry;
use super::{PspWriteError, footer, index};
use crate::locus_generation::SampleLocusObservations;
use crate::types::GenomePosition;

/// The header-parameter key under which a writer records the zstd level it compressed at.
///
/// **Named because more than one site must spell it identically**: `create` writes it, and
/// anything comparing a header it built against a header read back has to add it. A rename that
/// reached only some of them would look like a store defect rather than a typo.
pub const ZSTD_COMPRESSION_LEVEL_KEY: &str = "zstd-compression-level";

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
    /// How many of those bytes are the closing payload the caller supplied — the trailer.
    ///
    /// **Handed back because the caller has no other way to say what its own payload cost in
    /// the file**, and reopening the finished psp to read the footer would be a second answer
    /// to a question this call already knows. For a sample walk the payload is that sample's
    /// census (`psp_census_pair.md` §3.1), so this is how much of a psp the parameters fit
    /// reads and the rest is what it does not.
    pub trailer_bytes: u64,
}

/// **How many closed blocks may be waiting to be compressed before the walk has to wait.**
///
/// Two, not more: a block is the walk's own record buffer copied out, so every slot costs a
/// block payload of resident memory, and the walk closes one block per genomic block (100 kb by
/// default) while compressing one costs a few per cent of the walk. One in the queue and one
/// under the compressor is already enough for the compressor never to be the reason the walk
/// stops; a deeper queue would buy nothing and cost payloads.
const BLOCKS_AWAITING_COMPRESSION: usize = 2;

/// One closed block, on its way to the compressor.
struct BlockToCompress {
    head: BlockHead,
    /// The block payload — its head, then its records. A buffer recycled through
    /// [`BlockCompressionLine::spare_payloads`], not a fresh allocation per block.
    payload: Vec<u8>,
    /// Where the compressed block is to be written. Recycled the same way.
    frame: Vec<u8>,
}

/// One compressed block, on its way back to the writing thread.
struct CompressedBlock {
    head: BlockHead,
    /// The payload buffer, to be filled again.
    payload: Vec<u8>,
    /// `Ok` holds the whole on-disk block — its four-byte length, then its zstd frame. `Err`
    /// is why it could not be compressed, carried back rather than raised on the compressing
    /// thread, so that the walk still owns every failure it reports.
    frame: Result<Vec<u8>, BlockCompressError>,
}

/// **Compressing a psp's blocks on a thread of its own, while the walk carries on.**
///
/// A closed block is compressed and written and nothing waits for it: no record after it
/// depends on its bytes, and its index entry is not built until it is written. So the only
/// thing tying compression to the walk is *order* — the file's blocks and its index entries
/// must land in the order the walk closed them — and one compressing thread with FIFO queues
/// keeps that by construction. **The bytes are the same bytes**: [`BlockCompressor::compress`]
/// is a pure function of one payload, holding no state between blocks, so which thread ran it
/// cannot be read out of the file.
///
/// **What it costs in memory** is `BLOCKS_AWAITING_COMPRESSION` block payloads plus the frames
/// they compress to, all recycled rather than reallocated — on the default 100 kb grid, a few
/// megabytes at one sample. What it costs in threads is one, for as long as one psp is open.
struct BlockCompressionLine {
    to_compress: crossbeam_channel::Sender<BlockToCompress>,
    /// **Unbounded on purpose, and it is what makes this deadlock-free.** The walk can block
    /// handing a block over; the compressor must never block handing one back, or the two
    /// would wait on each other. It is bounded in fact by the sending side, which never has
    /// more than `BLOCKS_AWAITING_COMPRESSION` blocks in flight.
    compressed: crossbeam_channel::Receiver<CompressedBlock>,
    /// Joined when the writer is dropped, so no compressing thread outlives its psp.
    worker: Option<std::thread::JoinHandle<()>>,
    /// How many blocks the compressor has been handed and not yet given back.
    in_flight: usize,
    /// Payload and frame buffers that have come back, ready to be filled again.
    spare_payloads: Vec<Vec<u8>>,
    spare_frames: Vec<Vec<u8>>,
}

impl std::fmt::Debug for BlockCompressionLine {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("BlockCompressionLine")
            .field("in_flight", &self.in_flight)
            .finish_non_exhaustive()
    }
}

impl BlockCompressionLine {
    /// Start the compressing thread, which owns `compressor` from here on.
    fn start(compressor: BlockCompressor) -> Self {
        let (to_compress, jobs) = crossbeam_channel::bounded(BLOCKS_AWAITING_COMPRESSION);
        let (done, compressed) = crossbeam_channel::unbounded();
        let worker = std::thread::Builder::new()
            .name("psp-block-compression".to_string())
            .spawn(move || compress_blocks(compressor, &jobs, &done))
            .expect("a thread for compressing psp blocks");
        Self {
            to_compress,
            compressed,
            worker: Some(worker),
            in_flight: 0,
            spare_payloads: Vec::new(),
            spare_frames: Vec::new(),
        }
    }

    /// Hand one closed block over, having first written every block that has come back.
    ///
    /// The `&'static str` in the error is what [`PspWriter::spent`] records: which of the two
    /// unrecoverable things happened to a block whose records exist nowhere else.
    fn hand_over(
        &mut self,
        head: BlockHead,
        payload: &[u8],
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
    ) -> Result<(), (&'static str, PspWriteError)> {
        self.write_what_came_back(path, out, written, index)?;
        while self.in_flight >= BLOCKS_AWAITING_COMPRESSION {
            self.write_one(path, out, written, index)?;
        }
        let mut buffer = self.spare_payloads.pop().unwrap_or_default();
        buffer.clear();
        buffer.extend_from_slice(payload);
        let frame = self.spare_frames.pop().unwrap_or_default();
        self.to_compress
            .send(BlockToCompress {
                head,
                payload: buffer,
                frame,
            })
            .expect("the compressing thread outlives every block handed to it");
        self.in_flight += 1;
        Ok(())
    }

    /// Write every block already compressed, without waiting for any that is not.
    fn write_what_came_back(
        &mut self,
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
    ) -> Result<(), (&'static str, PspWriteError)> {
        while let Ok(block) = self.compressed.try_recv() {
            self.in_flight -= 1;
            self.write(block, path, out, written, index)?;
        }
        Ok(())
    }

    /// Wait for the next block to come back, and write it.
    fn write_one(
        &mut self,
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
    ) -> Result<(), (&'static str, PspWriteError)> {
        let block = self
            .compressed
            .recv()
            .expect("a block handed over comes back");
        self.in_flight -= 1;
        self.write(block, path, out, written, index)
    }

    /// Wait for every block still in flight, and write them all in order.
    fn write_all(
        &mut self,
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
    ) -> Result<(), (&'static str, PspWriteError)> {
        while self.in_flight > 0 {
            self.write_one(path, out, written, index)?;
        }
        Ok(())
    }

    fn write(
        &mut self,
        block: CompressedBlock,
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
    ) -> Result<(), (&'static str, PspWriteError)> {
        let CompressedBlock {
            head,
            payload,
            frame,
        } = block;
        self.spare_payloads.push(payload);
        let frame = match frame {
            Ok(frame) => frame,
            Err(source) => {
                return Err((
                    "a block could not be compressed",
                    PspWriteError::BlockRefused {
                        path: path.to_path_buf(),
                        source,
                    },
                ));
            }
        };
        let put = PspWriter::put_block(path, out, written, index, &head, &frame);
        self.spare_frames.push(frame);
        put.map_err(|error| ("a block could not be written", error))
    }
}

/// **Joined, not detached.** A psp that is finished, or one whose writer was dropped
/// half-written, must leave no thread behind — a cohort walked in one process would otherwise
/// accumulate one per sample.
impl Drop for BlockCompressionLine {
    fn drop(&mut self) {
        // Replacing the sender closes the one the worker is waiting on, which ends its loop.
        let (closed, _) = crossbeam_channel::bounded(0);
        self.to_compress = closed;
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

/// The compressing thread's whole life: take a block, compress it, hand both buffers back.
///
/// **A failure travels back down the channel rather than panicking here**, so that the walk
/// still reports every failure of its own file — a compressing thread that died would leave
/// the walk waiting on a block that never comes.
fn compress_blocks(
    mut compressor: BlockCompressor,
    jobs: &crossbeam_channel::Receiver<BlockToCompress>,
    done: &crossbeam_channel::Sender<CompressedBlock>,
) {
    while let Ok(job) = jobs.recv() {
        let BlockToCompress {
            head,
            payload,
            mut frame,
        } = job;
        let compressed = match compressor.compress(&payload) {
            Ok(bytes) => {
                frame.clear();
                frame.extend_from_slice(bytes);
                Ok(frame)
            }
            Err(source) => Err(source),
        };
        if done
            .send(CompressedBlock {
                head,
                payload,
                frame: compressed,
            })
            .is_err()
        {
            return;
        }
    }
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
    /// **Compression, on a thread of its own.** See [`BlockCompressionLine`].
    line: BlockCompressionLine,
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
    /// *does* warn about the milder destruction an append would cause (§6.4, not implemented:
    /// "write to a new path and rename if that matters"). A caller that must not destroy an
    /// existing psp checks for one first.
    ///
    /// **Nothing touches the filesystem until the header and the manifest have been accepted.**
    /// The header is encoded, and the cut rule and compressor built from its manifest, before
    /// the file is created — so a header this writer cannot honour leaves no file behind at all,
    /// rather than an empty one that every reader then refuses for a different reason.
    pub fn create(path: &Path, mut header: Header) -> Result<Self, PspWriteError> {
        // **The compression level is recorded because it is a setting, and goal 4 is that
        // settings are recorded.** A reader needs nothing from it — zstd decodes any level — but
        // byte-identity checks and the size figures the specs quote rest on it, and an append
        // (spec §6.4, not implemented) would need it to match bytes already in a file.
        // `block.rs`'s own doc assigns this to the writer by name.
        header.writer.parameters.insert(
            ZSTD_COMPRESSION_LEVEL_KEY.to_string(),
            crate::psp::header::ParameterValue::Integer(i64::from(
                super::block::ZSTD_COMPRESSION_LEVEL,
            )),
        );
        let header_bytes = header.encode()?;
        let builder = BlockBuilder::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedHeader {
                path: path.to_path_buf(),
                source: source.into(),
            }
        })?;
        let compressor = BlockCompressor::from_manifest(&header.manifest).map_err(|source| {
            PspWriteError::UnsupportedHeader {
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
            line: BlockCompressionLine::start(compressor),
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
        // **Destructured because the block goes to the compressing thread and the blocks that
        // have come back from it go to the file, and the two touch disjoint fields.** Naming
        // them says so to the borrow checker.
        //
        // **Handing a block over is not writing it**, so a failure raised here can belong to a
        // block closed several records ago. That costs nothing a reader can see — the file
        // stops at the same place either way, and `finish` refuses to seal it — and it is why
        // the line hands back the sentence [`spent`](Self::spent) records rather than the call
        // site guessing which of the two things went wrong.
        let Self {
            path,
            out,
            written,
            line,
            index,
            ..
        } = self;
        match line.hand_over(head, payload, path, out, written, index) {
            Ok(()) => Ok(()),
            Err((why, error)) => {
                self.spent = Some(why);
                Err(error)
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
            // **Destructured for the reason `push` destructures**: the compressing line and the
            // file are disjoint fields, and naming them says so.
            let Self {
                path,
                out,
                written,
                line,
                index,
                ..
            } = &mut self;
            line.hand_over(head, &payload, path, out, written, index)
                .map_err(|(_, error)| error)?;
        }
        // **Every block that is still in flight is written before the index is built**, because
        // the index describes them and the footer describes the index. This is the one place
        // the walk waits for the compressor, and by then there is nothing else left to do.
        {
            let Self {
                path,
                out,
                written,
                line,
                index,
                ..
            } = &mut self;
            line.write_all(path, out, written, index)
                .map_err(|(_, error)| error)?;
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
        // **Carried forward from two earlier reviews, and re-measured here**: a defect that
        // adds a fifth section between the trailer and the footer fails **82 tests with this
        // assertion and 66 without it**, so it accounts for sixteen of them. ⚠ The earlier figure
        // was *1 failing test to 7*, on a suite a third this size; a number measured on another
        // tree is not a fact about this one, and a review caught it still standing here.
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
            trailer_bytes: footer.trailer_bytes,
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
    /// This is the obligation the index and footer codecs both route here rather than solving:
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
    ///
    /// **The entry is pushed only if the write returned**, so the index never claims a block the
    /// file does not hold — and its offset is where the block went, taken before the write moves
    /// the counter.
    ///
    /// **The fields are named rather than reached through `self`**, and that is not a style
    /// choice: both callers hold `block` as a borrow of the compressor's own buffer, so a method
    /// taking `&mut self` could not be given it without copying the block out first — one
    /// allocation the size of a compressed block, per block, for nothing. Naming the disjoint
    /// fields says to the borrow checker that the compressor is not among them.
    fn put_block(
        path: &Path,
        out: &mut BufWriter<File>,
        written: &mut u64,
        index: &mut Vec<BlockIndexEntry>,
        head: &BlockHead,
        block: &[u8],
    ) -> Result<(), PspWriteError> {
        let entry = BlockIndexEntry {
            first_position: GenomePosition {
                contig: head.contig,
                position: head.first_position,
            },
            block_offset: *written,
        };
        out.write_all(block).map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "writing a block",
            source,
        })?;
        *written += block.len() as u64;
        index.push(entry);
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
    use crate::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::psp::header::{
        ContigIdentity, FORMAT_VERSION, ParameterValue, ReferenceIdentity, WriterProvenance,
    };
    use crate::types::{Bp, ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};
    use std::io::Read as _;

    /// A tomato-shaped header, with the block grid small enough that a handful of records cut
    /// several blocks.
    pub(crate) fn a_header(genomic_block_size_bp: u64) -> Header {
        let contigs = vec![
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
        ];
        let mut header = Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            segmentation_inputs: crate::psp::segmentation_section::segmentation_inputs_for_tests(
                &contigs,
            ),
            contigs,
            read_groups: crate::psp::header::read_groups_for_tests(),
            observation_reach_ceiling_bp: Bp(4_000),
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
            manifest: crate::psp::header::Manifest::as_this_build_writes_it(),
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
            chain_id_base: 0,
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(start + span - 1),
            },
            reference_bases: bases.clone(),
            observations: vec![SequenceObservation {
                flank_deletion: crate::locus_generation::FlankDeletion::NONE,
                held_flank: None,
                bases,
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 2,
                q_sum: SummedLogError::from_steps(-4_096),
                mapq_sum: 180,
                mapq_sum_sq: 10_800,
                placed_left: 1,
                num_low_quality: 0,
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
    pub(crate) fn footer_of(bytes: &[u8]) -> crate::psp::footer::Footer {
        let tail: [u8; crate::psp::footer::FOOTER_BYTES] = bytes
            [bytes.len() - crate::psp::footer::FOOTER_BYTES..]
            .try_into()
            .expect("the file is at least a footer long");
        crate::psp::footer::decode_footer(&tail).expect("a finished file's footer reads")
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
            let psp = crate::psp::PspReader::open(path).expect("a finished psp opens");
            let entries = psp.block_index();
            let at = entries[ordinal].block_offset as usize;
            let ends = entries
                .get(ordinal + 1)
                .map_or(psp.footer().index_offset as usize, |next| {
                    next.block_offset as usize
                });
            (at, ends)
        };
        whole[at + crate::psp::COMPRESSED_BLOCK_LENGTH_BYTES..ends].fill(0xff);
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
        a_file, a_finished_psp, a_header, a_record, a_sample, bytes_of, rewrite,
    };
    use super::*;
    use crate::psp::NotReadable;
    use crate::psp::PspReader;
    use crate::psp::footer::{FOOTER_BYTES, FOOTER_MAGIC, decode_footer};
    use crate::psp::index::decode_index;
    use crate::types::{ContigId, Position};

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
            crate::psp::index::checksum_index(index),
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
            let frame = match crate::psp::block::compressed_block_at(&bytes[at..]) {
                crate::psp::block::CompressedBlockAt::Whole { zstd_frame, .. } => zstd_frame,
                other => panic!("entry at byte {at} does not start a whole block: {other:?}"),
            };
            // `decompress` writes into spare capacity, so the room has to be there first.
            let mut payload = Vec::with_capacity(1 << 20);
            decompressor
                .decompress(&mut payload, frame)
                .expect("a block this writer wrote inflates");
            let head = crate::psp::block::BlockHead::decode(&payload)
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
        let header = crate::psp::read_header(&path).expect("the header survives");
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
    const THE_TEST_THE_CHILD_RUNS: &str =
        "psp::writer::tests::a_writer_killed_before_finishing_leaves_a_file_every_reader_refuses";

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
    /// its header still reads — which is spec §6.3's goal 3.
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
            // then slept, and a review measured what that did: the whole workload takes
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
            let (_header, header_bytes) = crate::psp::read_header_and_its_length(&sound)
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
        let (_, header_bytes) =
            crate::psp::read_header_and_its_length(&path).expect("the header survives the kill");
        let mut at = header_bytes;
        let mut whole_blocks = 0;
        while at + crate::psp::COMPRESSED_BLOCK_LENGTH_BYTES <= bytes.len() {
            let declared: [u8; crate::psp::COMPRESSED_BLOCK_LENGTH_BYTES] = bytes
                [at..at + crate::psp::COMPRESSED_BLOCK_LENGTH_BYTES]
                .try_into()
                .expect("a slice of exactly the prefix's width");
            let next = at
                + crate::psp::COMPRESSED_BLOCK_LENGTH_BYTES
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
            matches!(refused, crate::psp::PspReadError::Incomplete { .. }),
            "the instruction is `the run was interrupted, rebuild it` and not `the file is \
             damaged`; got {refused}"
        );
        // Spec §6.6: the half-written file still has a header, and this is what reads it.
        let header = crate::psp::read_header(&path).expect("the header survives the kill");
        assert_eq!(header.sample, "SRR7279481");
    }

    // -----------------------------------------------------------------
    // Worker-count invariance
    // -----------------------------------------------------------------

    /// The block grid the sharded sample is written on.
    ///
    /// ⚠ **Named because the mid-block precondition divides by it**, and the two used to be two
    /// literals that had to agree with nothing making them. Measured: writing the fixture on a
    /// 40 bp grid instead gives every record its own block — 1,000 instead of 40 — so every shard
    /// boundary lands on a grid line, the exact condition the precondition exists to refuse, **and
    /// the test still passed**. The filter is not conservative: a position off the 1 kb grid can
    /// sit exactly on a finer one.
    const GRID_BP_OF_THE_SHARDED_SAMPLE: u64 = 1_000;

    /// How many records each of the two contigs holds.
    const RECORDS_PER_CONTIG_IN_THE_SHARDED_SAMPLE: u64 = 500;

    /// How far apart consecutive records sit. **Not a divisor of the grid times any small worker
    /// count**, which is what keeps some shard boundary off the grid — see [`shard_boundaries`].
    const BASES_BETWEEN_RECORDS: u64 = 40;

    /// How many contigs the sharded sample spans, which is what the fixture header declares.
    const CONTIGS_IN_THE_SHARDED_SAMPLE: u64 = 2;

    /// Where one run's shards begin and end, as record ordinals into the whole sample.
    ///
    /// **Deliberately uneven.** An even split of two equal contigs puts the two-worker boundary
    /// exactly on the contig's first record — a grid line — and the four-worker boundaries on
    /// grid lines too, since 250 records at 40 bases is 10,000, a whole number of 1 kb cells. A
    /// sweep whose every boundary sits on a grid line cannot tell a cut that follows the
    /// coordinate from one that follows the shards, which is the thing it exists to tell apart.
    /// The `+ shard` skew is what moves them off.
    fn shard_boundaries(workers: u64) -> Vec<u64> {
        let records = CONTIGS_IN_THE_SHARDED_SAMPLE * RECORDS_PER_CONTIG_IN_THE_SHARDED_SAMPLE;
        let mut boundaries = vec![0];
        boundaries
            .extend((1..workers).map(|shard| (shard * records / workers + shard).min(records)));
        boundaries.push(records);
        boundaries.dedup();
        boundaries
    }

    /// The contig and position a record ordinal lands on.
    fn the_record_ordinal_lands_at(ordinal: u64) -> (u32, u64) {
        let per_contig = RECORDS_PER_CONTIG_IN_THE_SHARDED_SAMPLE;
        let contig = ordinal / per_contig;
        // **Checked, rather than left to `try_from`**, which only fails past `u32::MAX`: an
        // ordinal past the fixture would otherwise return contig 2 and surface much later, as a
        // record naming a contig the header does not declare.
        assert!(
            contig < CONTIGS_IN_THE_SHARDED_SAMPLE,
            "the fixture header declares {CONTIGS_IN_THE_SHARDED_SAMPLE} contigs, and ordinal \
             {ordinal} is past them"
        );
        (
            u32::try_from(contig).expect("a contig index fits a u32"),
            1 + (ordinal % per_contig) * BASES_BETWEEN_RECORDS,
        )
    }

    /// One sample's records, produced by `workers` threads each owning a contiguous span of the
    /// genome, then put back in coordinate order.
    ///
    /// **The threads are real and the spans are uneven**, so the interleaving in time genuinely
    /// varies with `workers`; what must not vary is the sequence that comes out.
    fn a_sample_gathered_by(workers: u64) -> Vec<SampleLocusObservations> {
        let boundaries = shard_boundaries(workers);
        let shards: Vec<(u64, u64)> = boundaries
            .windows(2)
            .map(|pair| (pair[0], pair[1]))
            .collect();
        let gathered: Vec<Vec<SampleLocusObservations>> = std::thread::scope(|threads| {
            let running: Vec<_> = shards
                .iter()
                .map(|&(from, to)| {
                    threads.spawn(move || {
                        (from..to)
                            .map(|at| {
                                let (contig, position) = the_record_ordinal_lands_at(at);
                                a_record(contig, position, 1)
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            running
                .into_iter()
                .map(|thread| thread.join().expect("a worker did not panic"))
                .collect()
        });
        // The shards are contiguous and in order, so concatenating them *is* the coordinate-order
        // merge a run would do.
        gathered.into_iter().flatten().collect()
    }

    /// Write a sample; give back the file's bytes **and how many blocks it cut**.
    ///
    /// **One writer, not three.** The block count used to be read by writing the same sample a
    /// second time and throwing the bytes away, which is the drift `a_finished_psp`'s own doc
    /// warns about — three spellings of *what a finished file looks like* in one module.
    fn the_file_written_from(records: &[SampleLocusObservations]) -> (Vec<u8>, u64) {
        let (_dir, path) = a_file();
        let mut writer =
            PspWriter::create(&path, a_header(GRID_BP_OF_THE_SHARDED_SAMPLE)).expect("a header");
        for record in records {
            writer.push(record).expect("in order");
        }
        let stats = writer.finish(b"a per-sample summary").expect("it finishes");
        (bytes_of(&path), stats.blocks)
    }

    /// **One sample gathered at 1, 2, 4, 8 and 16 workers gives byte-identical files** — spec §7's
    /// byte-identity claim.
    ///
    /// **What this can and cannot prove, stated plainly, because the difference decides what the
    /// test is worth.** ng's writer is serial: it is handed records in coordinate order and the
    /// producing is what a run parallelises. So byte-identity across worker counts rests on two
    /// things, and only the first is this file's doing:
    ///
    /// 1. **the writer is a function of the record sequence alone** — nothing about *when* records
    ///    arrive, or in what batches, reaches the bytes. **That half is already held** by
    ///    [`the_same_records_written_twice_give_the_same_bytes`], which is where a leaked
    ///    timestamp, a hash iteration order or an address would show;
    /// 2. **the sharding gives back the same sequence** whatever the worker count, which is the
    ///    run's job and is reproduced here so the claim is end to end rather than assumed.
    ///
    /// **And the block cut must not follow the shards.** The grid is a function of the coordinate
    /// (spec §4.1), so a shard boundary falling mid-block must not close one — which is only
    /// tested if some boundary actually falls mid-block. That is asserted below, not hoped for.
    #[test]
    fn one_sample_gathered_at_any_worker_count_gives_byte_identical_files() {
        let by_one = a_sample_gathered_by(1);
        let (expected, blocks) = the_file_written_from(&by_one);

        // **The fixture has to cut several blocks**, or "the cut does not follow the shards" is a
        // claim about a file with one block in it. It says nothing about *where* the boundaries
        // fall — that is the per-worker-count check below, and conflating the two is how this
        // message read before.
        assert!(
            blocks > 1,
            "the sample cuts {blocks} block; a one-block file makes \"the cut does not follow \
             the shards\" vacuous, whatever the boundaries do"
        );

        // **Counts that do not divide the records, as well as ones that do.** 1,000 records
        // divide evenly by 2, 4, 8 and 16, so the truncating branch of `records / workers` was
        // never taken with a remainder; 3, 5 and 7 take it.
        for workers in [2u64, 3, 4, 5, 7, 8, 16] {
            // **A shard boundary that falls mid-block**, without which a writer that closed a
            // block at every shard boundary would produce the same file anyway.
            let boundaries = shard_boundaries(workers);
            let mid_block_boundaries = boundaries[1..boundaries.len() - 1]
                .iter()
                .filter(|at| {
                    let (_contig, position) = the_record_ordinal_lands_at(**at);
                    (position - 1) % GRID_BP_OF_THE_SHARDED_SAMPLE != 0
                })
                .count();
            assert!(
                mid_block_boundaries > 0,
                "at {workers} workers every shard boundary lands on a grid line, so this round \
                 cannot tell a coordinate cut from a shard cut"
            );

            let gathered = a_sample_gathered_by(workers);
            assert_eq!(
                gathered.len(),
                by_one.len(),
                "at {workers} workers the sharding produced a different number of records"
            );
            assert!(
                gathered == by_one,
                "at {workers} workers the sharding produced a different record sequence"
            );
            assert!(
                the_file_written_from(&gathered).0 == expected,
                "at {workers} workers the file differs from the one-worker file"
            );
        }
    }

    /// **The header's timestamp is the only thing allowed to differ, and it differs only inside
    /// the header** — the second half of spec §7's sentence, which nothing held.
    ///
    /// The invariance claim is *byte-identical apart from any timestamp in the header*. The sweep
    /// above holds the first half by keeping `created` fixed; this holds the second, by writing
    /// the same records twice with two different stamps and requiring every differing byte to lie
    /// inside the header.
    ///
    /// **Why it is worth a test of its own.** The timestamp is written into the header's TOML, so
    /// a stamp of a different rendered *width* moves the header's length — and with it every
    /// offset in the footer, the index, and the blocks' positions in the file. A run whose files
    /// differed past the header would be one where nothing downstream could compare two samples
    /// byte for byte, which is what §7 offers.
    #[test]
    fn two_runs_differing_only_in_their_timestamp_differ_only_inside_the_header() {
        let write_stamped = |created: &str| {
            let (dir, path) = a_file();
            let mut header = a_header(GRID_BP_OF_THE_SHARDED_SAMPLE);
            header.writer.created = created.parse().expect("a valid RFC 3339 stamp");
            let mut writer = PspWriter::create(&path, header).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            let _ = writer.finish(b"a per-sample summary").expect("it finishes");
            let bytes = bytes_of(&path);
            drop(dir);
            bytes
        };
        // Two stamps of the same rendered width — **and that is a precondition, not a
        // description**: see the sibling test below, which measures what a wider one does.
        let early = write_stamped("2026-08-28T00:00:00Z");
        let later = write_stamped("2026-08-30T11:22:33Z");
        assert_eq!(
            early.len(),
            later.len(),
            "two stamps of the same width must not change the file's length"
        );

        let (_header, header_bytes) = {
            let (dir, path) = a_file();
            rewrite(&path, &early);
            let read = crate::psp::read_header_and_its_length(&path).expect("the header reads");
            drop(dir);
            read
        };
        let differing: Vec<usize> = early
            .iter()
            .zip(&later)
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(at, _)| at)
            .collect();
        assert!(
            !differing.is_empty(),
            "two different stamps have to change something, or this test proves nothing"
        );
        assert!(
            differing.iter().all(|at| *at < header_bytes),
            "the two files differ outside the {header_bytes}-byte header, at {:?} — a timestamp \
             is the one field spec §7 allows to differ, and only inside the header",
            differing
                .iter()
                .filter(|at| **at >= header_bytes)
                .take(8)
                .collect::<Vec<_>>()
        );
    }

    /// ⚠ **Spec §7's byte-identity holds only while every stamp renders to the same width, and
    /// nothing enforces that.** This test pins the limit rather than the promise.
    ///
    /// `created` is a `toml::value::Datetime`, which carries an optional sub-second fraction and
    /// an offset, so `2026-08-30T11:22:33.5Z` is two characters wider than
    /// `2026-08-28T00:00:00Z`. The stamp goes into the header's TOML, so those two characters
    /// move the header's length — **and with it every offset in the footer, the index, and where
    /// the blocks sit in the file**. Measured here: the two files differ in length, and the
    /// difference is exactly the two characters.
    ///
    /// **Why this matters beyond the test.** ng has no production header writer yet — every
    /// `WriterProvenance` built under `src/` is in a test. A future `pileup` rendering the
    /// clock with `to_rfc3339`, which prints sub-second digits only when they are non-zero, would
    /// produce stamps of *varying* width run to run, and §7's promise — the thing that lets
    /// anything downstream compare two samples byte for byte — would quietly stop holding.
    ///
    /// **The fix is not this test's to make.** Either the writer normalises `created` to whole
    /// seconds in UTC, or it refuses a stamp that does not render at a fixed width, or spec §7
    /// gains the words "of fixed rendered width". All three are design decisions and are raised
    /// rather than taken.
    #[test]
    fn a_wider_timestamp_moves_every_offset_in_the_file_and_that_is_the_limit_of_the_promise() {
        let write_stamped = |created: &str| {
            let (dir, path) = a_file();
            let mut header = a_header(GRID_BP_OF_THE_SHARDED_SAMPLE);
            header.writer.created = created.parse().expect("a valid RFC 3339 stamp");
            let mut writer = PspWriter::create(&path, header).expect("a header");
            for record in a_sample() {
                writer.push(&record).expect("in order");
            }
            let _ = writer.finish(b"a per-sample summary").expect("it finishes");
            let bytes = bytes_of(&path);
            drop(dir);
            bytes
        };
        let narrow = write_stamped("2026-08-28T00:00:00Z");
        let wider = write_stamped("2026-08-30T11:22:33.5Z");
        assert_eq!(
            wider.len(),
            narrow.len() + 2,
            "the wider stamp is two characters wider, and the file grows by exactly that"
        );

        let (_header, header_bytes) = {
            let (dir, path) = a_file();
            rewrite(&path, &narrow);
            let read = crate::psp::read_header_and_its_length(&path).expect("the header reads");
            drop(dir);
            read
        };
        // **And the difference does not stay inside the header**, which is the whole point: the
        // shorter file's own bytes past its header no longer match the wider file's.
        let past_the_header = narrow
            .iter()
            .zip(&wider)
            .skip(header_bytes)
            .any(|(a, b)| a != b);
        assert!(
            past_the_header,
            "a wider stamp has to move the bytes past the header — if it did not, spec §7 would \
             hold at any width and this test would be recording a limit that is not there"
        );
    }

    /// **The shard split never produces an empty or backwards shard**, even asked for more
    /// workers than there are records — the case `.min` and `.dedup` in [`shard_boundaries`] exist
    /// for and that no worker count in the sweep reaches.
    #[test]
    fn the_shard_split_stays_ordered_with_more_workers_than_records() {
        let records = CONTIGS_IN_THE_SHARDED_SAMPLE * RECORDS_PER_CONTIG_IN_THE_SHARDED_SAMPLE;
        for workers in [1u64, 2, 999, records, records + 1, 4 * records] {
            let boundaries = shard_boundaries(workers);
            assert!(
                boundaries.windows(2).all(|pair| pair[0] < pair[1]),
                "at {workers} workers a shard may be small but never empty or backwards: \
                 {boundaries:?}"
            );
            assert_eq!(
                boundaries.first(),
                Some(&0),
                "the split starts at the first record"
            );
            assert_eq!(
                boundaries.last(),
                Some(&records),
                "and ends at the last, or records go unwritten"
            );
        }
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
    /// ⚠ **This was once recorded as untestable without a failing file descriptor.** It was
    /// testable here; a review found the route. Linux only — macOS has no
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
        // whole suite green until a review looked at it.
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
    /// goal 4 is that settings are recorded.
    #[test]
    fn the_compression_level_reaches_the_header() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let _ = writer.finish(&[]).expect("it finishes");

        let header = crate::psp::read_header(&path).expect("the header reads");
        assert_eq!(
            header.writer.parameters.get("zstd-compression-level"),
            Some(&crate::psp::header::ParameterValue::Integer(i64::from(
                crate::psp::block::ZSTD_COMPRESSION_LEVEL
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
    // Abandoned writers, and a footer a reader must refuse
    // -----------------------------------------------------------------

    /// **A writer abandoned with blocks still being compressed must not hang, and must leave no
    /// thread behind.**
    ///
    /// Since compression moved to a thread of its own, dropping a writer has to close the
    /// channel that thread is waiting on and join it. Getting that wrong has two shapes and
    /// neither shows up as a failing assertion elsewhere: a join on a thread whose channel was
    /// never closed hangs the run at the drop, and a detached thread leaks — one per sample, in
    /// a cohort walked in a single process, each holding a block payload.
    ///
    /// So this abandons twenty writers in turn, each with several blocks closed and unwritten,
    /// in one process. A hang fails it by timing out; a leak is what the twenty are for. The
    /// existing drop tests cannot reach either: they push one record, which closes no block, so
    /// nothing is ever in flight.
    #[test]
    fn writers_abandoned_with_blocks_in_flight_neither_hang_nor_pile_up() {
        let (_dir, path) = a_file();
        for round in 0..20u64 {
            let mut writer =
                PspWriter::create(&path, a_header(1_000)).expect("a header for each round");
            // Ten genomic blocks at a thousand bases each, so blocks are closed and handed over
            // rather than only built.
            for step in 0..10u64 {
                writer
                    .push(&a_record(0, 1 + step * 1_000 + round, 1))
                    .expect("each record is past the last");
            }
            // No `finish`: this is what a killed run leaves, and it is the path that has to
            // close the channel rather than wait on it.
            drop(writer);
        }

        // The last file is what a killed run leaves — headed, unsealed, and refused by a reader.
        let refused = PspReader::open(&path).expect_err("nothing sealed the last one");
        assert!(
            matches!(refused, crate::psp::PspReadError::Incomplete { .. }),
            "an unsealed file is refused, not opened: {refused:?}",
        );
    }

    /// **A footer that puts the block index inside the header is refused at open.**
    ///
    /// **The per-entry rule cannot carry it: on an empty index there are no entries to check.** A
    /// footer saying the index is at byte 4 and holds nothing once passed `PspReader::open`, and
    /// the append writer this tree had then (deleted 2026-11-28, psp spec §6.4) cut a 3,742-byte
    /// psp down to 109 bytes and returned `Ok`.
    #[test]
    fn a_footer_that_puts_the_block_index_inside_the_header_is_refused() {
        use crate::psp::footer::Footer;

        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let mut crafted = whole[..whole.len() - FOOTER_BYTES].to_vec();
        let footer = Footer {
            index_offset: 4,
            index_bytes: 0,
            trailer_offset: 4,
            trailer_bytes: (crafted.len() - 4) as u64,
            n_blocks: 0,
            index_checksum: crate::psp::checksum_index(&[]),
        };
        crafted.extend_from_slice(&crate::psp::encode_footer(&footer));
        rewrite(&path, &crafted);

        let refused = PspReader::open(&path).expect_err("the index is inside the header");
        match refused {
            crate::psp::PspReadError::Damaged { reason, .. } => {
                assert!(reason.contains("inside the"), "got {reason}")
            }
            other => panic!("got {other}"),
        }
    }
}
