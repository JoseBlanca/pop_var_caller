//! **The inflated BGZF blocks one open BAM's readers share**, so a block two of them read is
//! inflated once rather than once per reader.
//!
//! # Why a BAM has several readers at all
//!
//! A walk drives more than one locus generator over a sample — the SNP/indel one, the
//! repeat-tract one and the STR bundle one — and each holds its own cursor, hence its own
//! [`BamAlignedReadsReader`](super::BamAlignedReadsReader), on purpose: their regions
//! interleave, and one cursor between them would tie their lifetimes together (`walker.rs`'s
//! `generic_path_generators`). Interleaving is exactly the case where they all pass through
//! the same stretch of the file, so before this cache **each reader inflated every block it
//! passed**, and a block of a stretch all three walk was inflated three times. Profiled on
//! `generate-psps` over 5 simulated samples at 30x, inflating was 10,147 samples of 84,427
//! (12 in 100) with all three readers and 6,862 with the bundle reader off.
//!
//! The CRAM arm solved the same problem first, one level up — it shares decoded *containers*
//! ([`SharedContainers`](super::SharedContainers)). This is the BAM counterpart, at the BAM's unit
//! of decompression, and it borrows that design's shape: one cache per open file, a fixed number
//! of slots given up least recently used first, keyed by where the unit starts in the file.
//!
//! # What is shared and what is not
//!
//! **Only the inflate is shared.** Each reader keeps its own place — the block it is in and how
//! far into it — and its own record decoding, so readers stay as independent as before; they
//! simply find the bytes already unpacked when another reader has just been there.
//!
//! # Why the block layer is ours rather than noodles'
//!
//! noodles' BGZF reader cannot be handed a block it did not inflate: it owns one buffer and its
//! `seek` re-inflates unconditionally (noodles-bgzf 0.47.0, `src/io/reader.rs:175-186`). So this
//! module reads the frames itself and hands noodles' BAM reader a plain byte stream
//! ([`SharedBgzfReader`]), which is all `bam::io::Reader::read_record_buf` asks of it.
//!
//! **Reading a frame mirrors noodles' `io/reader/frame.rs` check for check**, with the same
//! error kinds and messages, so a corrupt or truncated file fails as it did: a bad header is
//! `invalid BGZF header`, a bad checksum `block data checksum mismatch`, a `BSIZE` too small to
//! be a frame `invalid frame size`, a file ending inside a frame `UnexpectedEof`, and a file
//! ending *between* frames is the end of the stream. Two cases differ, both only for the better,
//! and neither can arise from a well-formed file and its index: a frame claiming more than
//! 64 KiB of data is refused as `InvalidData`, and a seek to an offset past the end of its block
//! as `InvalidInput`, where noodles would panic indexing past its buffer in both.

use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::{Arc, Mutex};

use noodles_bgzf::VirtualPosition;

/// How many inflated blocks one open file keeps — see [`SharedBgzfBlocks`] for why 32.
const CACHE_SLOTS: usize = 32;

/// The largest amount of data a BGZF block may hold (the format's `ISIZE` limit).
const MAX_BLOCK_DATA: usize = 1 << 16;

/// A BGZF frame's fixed header: the gzip header, `XLEN`, and the `BC` subfield holding `BSIZE`.
const FRAME_HEADER_SIZE: usize = 18;

/// A gzip trailer: the data's CRC32 and its length.
const FRAME_TRAILER_SIZE: usize = 8;

/// One BGZF block, inflated, and where it sits in the file.
pub(crate) struct InflatedBlock {
    /// Where this block's compressed frame starts — the first half of a virtual position.
    compressed_offset: u64,
    /// Where the next frame starts, which is where reading on goes.
    next_compressed_offset: u64,
    /// The block's data. Never empty except for the block that stands for the end of the
    /// stream: an empty block in the middle of a file is stepped over, as noodles steps over it.
    data: Box<[u8]>,
}

impl InflatedBlock {
    /// A block holding nothing, standing at `offset`: where a reader starts, and what the end of
    /// the stream reads as.
    fn empty_at(offset: u64) -> Self {
        Self {
            compressed_offset: offset,
            next_compressed_offset: offset,
            data: Box::new([]),
        }
    }
}

/// How often a reader found its block already inflated, and how often it had to inflate one.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BlockCacheCounts {
    /// Blocks inflated from the file — the work this cache exists to reduce.
    pub(crate) blocks_inflated: u64,
    /// Requests answered from a block already held.
    pub(crate) blocks_reused: u64,
}

/// Everything behind the lock: the file the blocks come from, the slots, and the scratch an
/// inflate needs.
struct BlockCacheState {
    /// The one descriptor every reader of this file inflates through. A reader holds no file of
    /// its own: what it needs from the file is blocks, and those come from here.
    file: File,
    /// The slots, least recently used first, each keyed by its block's compressed offset.
    held: Vec<Arc<InflatedBlock>>,
    /// One compressed frame, reused for every inflate.
    frame: Vec<u8>,
    /// libdeflate's decompressor state, reused rather than allocated per block.
    decompressor: libdeflater::Decompressor,
    counts: BlockCacheCounts,
}

/// The inflated blocks of one open BAM, shared by every reader on it.
///
/// **One per open file is one per walk.** `SampleReads` is not `Clone` and opens its own
/// [`AlignmentFile`](crate::read::input::AlignmentFile) per sample, so nothing shares this
/// between two samples or between two workers.
///
/// # Thirty-two slots, because the bundle reader trails the others by dozens of blocks
///
/// How far apart the readers are, counted in blocks, is what sets the slot count — a block can
/// only be shared if it is still held when the reader furthest behind reaches it. Measured with
/// `generate-psps` on one simulated sample at 30x (3 Mb, 4,114 distinct blocks, 12,334 block
/// requests from the three readers):
///
/// | slots | 1 | 4 | 8 | 16 | 32 | 64 | 128 |
/// |---|---:|---:|---:|---:|---:|---:|---:|
/// | blocks inflated, STR bundles on | 10,026 | 7,763 | 6,864 | 5,510 | **4,572** | 4,184 | 4,114 |
/// | blocks inflated, STR bundles off | 6,018 | 4,129 | — | 4,114 | — | 4,114 | — |
///
/// The generic and repeat-tract readers stay within four blocks of each other: with the
/// bundles off, four slots inflate 4,129 blocks of the 4,114 there are. The bundle reader is
/// the one that trails. It is only moved when the next bundle's region arrives, and a cursor
/// serves a region ahead of it by reading forward rather than seeking (`cursor.rs`, the forget
/// rule) — so it then reads every block from the last bundle to this one, blocks the other two
/// passed as they walked that stretch. Its lag is the gap between bundles, dozens of blocks at
/// 30x. Thirty-two slots inflate 4,572 blocks where the readers used to inflate 12,334 between
/// them: 94 in 100 of the inflates sharing can save, against 99 in 100 for 64 slots at twice the
/// memory. On HG002 at 30x (`call-from-alignments` over the tandem-repeat tier, one sample) the
/// same 32 slots inflate 16,433 of 46,391 requested blocks, where unbounded slots would inflate
/// 15,579.
///
/// **The distance is in bases, so the slots needed grow with depth.** A block holds about
/// 64 KiB of records — some 700 bases of genome at 30x, 70 at 300x, 7 kb at 3x. At 3x a few
/// slots catch everything. At 300x — arithmetic, not measured — the generic and repeat-tract
/// readers' gap would be some 40 blocks, more than 32 slots hold, so part of their sharing is
/// lost as well as most of the bundle reader's. Nothing is ever worse than before the cache:
/// a block not held is inflated, exactly as every reader used to.
///
/// **Memory: at most 32 blocks of at most 64 KiB, 2 MiB, per open BAM**, plus the block each
/// reader is in when the cache has already let it go, and one compressed frame of scratch. A
/// reader no longer holds a noodles BGZF reader — a 64 KiB block buffer and a frame buffer —
/// which pays back up to 128 KiB per reader. The
/// slots fill only as blocks are read and are released when the file is closed. A run holding
/// every sample open at once, as a cohort `call-from-alignments` does, holds up to 2 MiB per
/// sample: 2 GiB at a thousand samples.
///
/// # The lock is held through an inflate, and that costs nothing here
///
/// A sample's readers belong to one walk, and a walk is driven by one thread at a time, so they
/// never ask for blocks at the same moment and the lock is never contended. Holding it through
/// the inflate is what keeps this simple: no block is ever inflated twice because two readers
/// raced for it, with no claim protocol to get right (the CRAM cache needs one because a
/// prefetch thread decodes beside the walk). Should readers of one file ever run on two threads
/// at once, they would take turns inflating — slower, never wrong.
pub(crate) struct SharedBgzfBlocks {
    state: Mutex<BlockCacheState>,
}

impl SharedBgzfBlocks {
    /// Open `path` for inflating its blocks. Nothing is read yet.
    pub(crate) fn open(path: &Path) -> io::Result<Self> {
        Ok(Self {
            state: Mutex::new(BlockCacheState {
                file: File::open(path)?,
                held: Vec::with_capacity(CACHE_SLOTS + 1),
                frame: Vec::new(),
                decompressor: libdeflater::Decompressor::new(),
                counts: BlockCacheCounts::default(),
            }),
        })
    }

    /// How many blocks were inflated and how many requests were answered without inflating.
    #[cfg(test)]
    pub(crate) fn counts(&self) -> BlockCacheCounts {
        self.lock().counts
    }

    /// The compressed offsets of the blocks held, least recently used first.
    #[cfg(test)]
    fn held_offsets(&self) -> Vec<u64> {
        self.lock()
            .held
            .iter()
            .map(|block| block.compressed_offset)
            .collect()
    }

    /// The first non-empty block at or after `compressed_offset`: held if it is, inflated if not.
    ///
    /// At the end of the stream this is an empty block standing at `compressed_offset`, which is
    /// not kept — it would give up a slot for nothing.
    fn block_at(&self, compressed_offset: u64) -> io::Result<Arc<InflatedBlock>> {
        let mut state = self.lock();
        if let Some(at) = state
            .held
            .iter()
            .position(|block| block.compressed_offset == compressed_offset)
        {
            // Most recently used goes last, so the slot given up next is the one longest unasked.
            let block = state.held.remove(at);
            state.held.push(Arc::clone(&block));
            state.counts.blocks_reused += 1;
            return Ok(block);
        }

        let block = Arc::new(state.inflate_at(compressed_offset)?);
        if !block.data.is_empty() {
            state.counts.blocks_inflated += 1;
            // Asked for an empty block, this inflated the one after it, which may already be
            // held under its own offset. One slot per offset: the fresh copy replaces it.
            if let Some(at) = state
                .held
                .iter()
                .position(|held| held.compressed_offset == block.compressed_offset)
            {
                state.held.remove(at);
            } else if state.held.len() >= CACHE_SLOTS {
                state.held.remove(0);
            }
            state.held.push(Arc::clone(&block));
        }
        Ok(block)
    }

    /// **A poisoned lock is recovered from rather than propagated.** What it holds is inflated
    /// blocks, a descriptor that is repositioned before every read, and scratch — nothing a
    /// panicking thread could have left half-written for the next one to trust.
    fn lock(&self) -> std::sync::MutexGuard<'_, BlockCacheState> {
        self.state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }
}

impl BlockCacheState {
    /// Inflate the first non-empty block at or after `compressed_offset`, stepping over empty
    /// ones as noodles does.
    fn inflate_at(&mut self, compressed_offset: u64) -> io::Result<InflatedBlock> {
        self.file.seek(SeekFrom::Start(compressed_offset))?;
        let mut frame_start = compressed_offset;
        loop {
            if !read_frame_into(&mut self.file, &mut self.frame)? {
                return Ok(InflatedBlock::empty_at(frame_start));
            }
            let next_frame_start = frame_start + self.frame.len() as u64;
            let data = inflate_frame(&self.frame, &mut self.decompressor)?;
            if !data.is_empty() {
                return Ok(InflatedBlock {
                    compressed_offset: frame_start,
                    next_compressed_offset: next_frame_start,
                    data,
                });
            }
            frame_start = next_frame_start;
        }
    }
}

/// Read one whole BGZF frame into `frame`. `false` when the stream ends before a frame header
/// does — noodles' end of stream, and the same rule: a header cut short is the end, a frame cut
/// short is an error.
fn read_frame_into(file: &mut File, frame: &mut Vec<u8>) -> io::Result<bool> {
    frame.resize(FRAME_HEADER_SIZE, 0);
    match file.read_exact(frame) {
        Ok(()) => {}
        Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(error) => return Err(error),
    }
    // `BSIZE` is the header's last two bytes: the frame's total size, minus one.
    let frame_size = usize::from(u16::from_le_bytes([frame[16], frame[17]])) + 1;
    if frame_size < FRAME_HEADER_SIZE + FRAME_TRAILER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid frame size",
        ));
    }
    frame.resize(frame_size, 0);
    file.read_exact(&mut frame[FRAME_HEADER_SIZE..])?;
    Ok(true)
}

/// Check one frame's header, inflate its data, and check the data against the trailer's CRC32.
fn inflate_frame(
    frame: &[u8],
    decompressor: &mut libdeflater::Decompressor,
) -> io::Result<Box<[u8]>> {
    const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
    const DEFLATE: u8 = 0x08;
    const FEXTRA: u8 = 0x04;
    let header = &frame[..FRAME_HEADER_SIZE];
    let is_bgzf = header[0..2] == GZIP_MAGIC
        && header[2] == DEFLATE
        && header[3] == FEXTRA
        && header[10..12] == [0x06, 0x00] // XLEN
        && header[12..14] == *b"BC" // the BGZF subfield
        && header[14..16] == [0x02, 0x00]; // its length
    if !is_bgzf {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }

    let trailer_start = frame.len() - FRAME_TRAILER_SIZE;
    let compressed = &frame[FRAME_HEADER_SIZE..trailer_start];
    let trailer = &frame[trailer_start..];
    let expected_crc32 = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let data_length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
    if data_length > MAX_BLOCK_DATA {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BGZF block claims more than 64 KiB of data",
        ));
    }

    let mut data = vec![0u8; data_length].into_boxed_slice();
    decompressor
        .deflate_decompress(compressed, &mut data)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    let mut crc = libdeflater::Crc::new();
    crc.update(&data);
    if crc.sum() != expected_crc32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ));
    }
    Ok(data)
}

/// One reader's place in a BAM's uncompressed stream: the block it is in and how far into it,
/// with the blocks themselves drawn from the file's [`SharedBgzfBlocks`].
///
/// A drop-in for the noodles BGZF reader as far as a BAM record reader is concerned: it is
/// `Read` and `BufRead`, it seeks to a [`VirtualPosition`], and it reports one with the same
/// rule noodles uses — inside a block it is *(this block, offset)*; at a block's end it is
/// *(next block, 0)*. That rule matters, because the BAM reader compares it against the end of
/// an index chunk to know when the chunk is done.
pub(crate) struct SharedBgzfReader {
    blocks: Arc<SharedBgzfBlocks>,
    block: Arc<InflatedBlock>,
    /// How far into `block.data` this reader has read.
    position_in_block: usize,
}

impl SharedBgzfReader {
    /// A reader at the start of the file, as a freshly opened noodles reader is.
    pub(crate) fn new(blocks: Arc<SharedBgzfBlocks>) -> Self {
        Self {
            blocks,
            block: Arc::new(InflatedBlock::empty_at(0)),
            position_in_block: 0,
        }
    }

    /// Move to `position`: into the block it names, at the offset it names.
    pub(crate) fn seek(&mut self, position: VirtualPosition) -> io::Result<()> {
        let block = self.blocks.block_at(position.compressed())?;
        let offset = usize::from(position.uncompressed());
        if offset > block.data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "virtual position lies past the end of its BGZF block",
            ));
        }
        self.block = block;
        self.position_in_block = offset;
        Ok(())
    }

    /// Where this reader is, by noodles' rule (see the type).
    pub(crate) fn virtual_position(&self) -> VirtualPosition {
        let (compressed, uncompressed) = if self.position_in_block < self.block.data.len() {
            (self.block.compressed_offset, self.position_in_block)
        } else {
            (self.block.next_compressed_offset, 0)
        };
        // An offset inside a block is below 64 KiB, and a compressed offset read from a file is
        // below 2^48 on any file noodles' own reader would have opened.
        VirtualPosition::try_from((compressed, uncompressed as u16))
            .expect("a block's offsets fit in a virtual position")
    }
}

impl BufRead for SharedBgzfReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.position_in_block >= self.block.data.len() {
            self.block = self.blocks.block_at(self.block.next_compressed_offset)?;
            self.position_in_block = 0;
        }
        Ok(&self.block.data[self.position_in_block..])
    }

    fn consume(&mut self, amount: usize) {
        self.position_in_block = (self.position_in_block + amount).min(self.block.data.len());
    }
}

impl Read for SharedBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let amount = available.len().min(buf.len());
        buf[..amount].copy_from_slice(&available[..amount]);
        self.consume(amount);
        Ok(amount)
    }

    /// The common case — the bytes are all in the current block — is one copy.
    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        if let Some(available) = self
            .block
            .data
            .get(self.position_in_block..self.position_in_block + buf.len())
        {
            buf.copy_from_slice(available);
            self.position_in_block += buf.len();
            return Ok(());
        }
        let mut rest = buf;
        while !rest.is_empty() {
            match self.read(rest)? {
                0 => {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "failed to fill whole buffer",
                    ));
                }
                amount => rest = &mut rest[amount..],
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use noodles_bgzf as bgzf;

    /// The bytes of a BGZF stream of `blocks` blocks, each holding `block_length` bytes of a
    /// pattern that differs from block to block, then the end-of-file marker (an empty block).
    fn bgzf_bytes(blocks: usize, block_length: usize) -> Vec<u8> {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        for block in 0..blocks {
            let bytes: Vec<u8> = (0..block_length)
                .map(|i| ((block * 31 + i * 7) % 251) as u8)
                .collect();
            writer.write_all(&bytes).expect("write");
            writer.flush().expect("one block per flush");
        }
        writer.finish().expect("finish")
    }

    /// `bytes` written to a file. The directory is returned to keep the file alive.
    fn file_of(bytes: &[u8]) -> (tempfile::TempDir, std::path::PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("blocks.gz");
        std::fs::write(&path, bytes).expect("write");
        (dir, path)
    }

    /// Where each frame of a BGZF stream starts, read off the frames' own `BSIZE` fields.
    fn frame_offsets(bytes: &[u8]) -> Vec<u64> {
        let mut offsets = Vec::new();
        let mut at = 0;
        while at < bytes.len() {
            offsets.push(at as u64);
            at += usize::from(u16::from_le_bytes([bytes[at + 16], bytes[at + 17]])) + 1;
        }
        offsets
    }

    fn shared(path: &Path) -> Arc<SharedBgzfBlocks> {
        Arc::new(SharedBgzfBlocks::open(path).expect("open"))
    }

    fn noodles_reader(path: &Path) -> bgzf::io::Reader<File> {
        bgzf::io::Reader::new(File::open(path).expect("open"))
    }

    /// Read both to the end in pieces that straddle block boundaries, checking after every piece
    /// that they returned the same bytes and stand at the same virtual position.
    fn assert_reads_alike(ours: &mut SharedBgzfReader, theirs: &mut bgzf::io::Reader<File>) {
        let mut piece = [0u8; 333];
        loop {
            assert_eq!(ours.virtual_position(), theirs.virtual_position());
            let read_ours = ours.read(&mut piece).expect("ours reads");
            let ours_bytes = piece[..read_ours].to_vec();
            let read_theirs = theirs.read(&mut piece).expect("noodles reads");
            assert_eq!(ours_bytes, piece[..read_theirs]);
            if read_ours == 0 {
                assert_eq!(ours.virtual_position(), theirs.virtual_position());
                return;
            }
        }
    }

    /// **The positions are the contract that matters**, not only the bytes: the BAM reader
    /// ends an index chunk by comparing a virtual position against the chunk's end, so a reader
    /// that reported *(block, end)* where noodles reports *(next block, 0)* would read one record
    /// too many or too few at every chunk boundary.
    #[test]
    fn reading_through_matches_noodles_in_bytes_and_in_positions() {
        let (_dir, path) = file_of(&bgzf_bytes(5, 1_000));
        assert_reads_alike(
            &mut SharedBgzfReader::new(shared(&path)),
            &mut noodles_reader(&path),
        );
    }

    /// An empty block in the middle of a stream — here the end-of-file marker of the first of
    /// two concatenated streams — is stepped over, as noodles steps over it, both when reading
    /// on into it and when seeking straight to it.
    #[test]
    fn an_empty_block_mid_stream_is_stepped_over_as_noodles_steps_over_it() {
        let first = bgzf_bytes(2, 500);
        let mut bytes = first.clone();
        bytes.extend(bgzf_bytes(2, 700));
        let (_dir, path) = file_of(&bytes);
        assert_reads_alike(
            &mut SharedBgzfReader::new(shared(&path)),
            &mut noodles_reader(&path),
        );

        // The first stream's end-of-file marker is its last frame.
        let empty_block = *frame_offsets(&first).last().expect("frames");
        let position = VirtualPosition::try_from((empty_block, 0)).expect("fits");
        let mut ours = SharedBgzfReader::new(shared(&path));
        let mut theirs = noodles_reader(&path);
        ours.seek(position).expect("ours seeks");
        theirs.seek(position).expect("noodles seeks");
        assert_reads_alike(&mut ours, &mut theirs);
    }

    /// A seek into the middle of every block lands where noodles lands — including at a block's
    /// very end (offset 2,000 of a 2,000-byte block), which reads on from the next block.
    #[test]
    fn seeking_lands_where_noodles_lands() {
        let bytes = bgzf_bytes(4, 2_000);
        let (_dir, path) = file_of(&bytes);
        let blocks = shared(&path);
        let offsets = frame_offsets(&bytes);
        // The last frame is the end-of-file marker, which holds no data to land inside.
        for &block in &offsets[..offsets.len() - 1] {
            for within in [0u16, 1, 999, 1_999, 2_000] {
                let position = VirtualPosition::try_from((block, within)).expect("fits");
                let mut ours = SharedBgzfReader::new(Arc::clone(&blocks));
                let mut theirs = noodles_reader(&path);
                ours.seek(position).expect("ours seeks");
                theirs.seek(position).expect("noodles seeks");
                assert_reads_alike(&mut ours, &mut theirs);
            }
        }
    }

    /// A seek to the end-of-file marker reads nothing and stands at the end of the file, as
    /// noodles does. A seek to the file's length is checked on ours alone: noodles keeps serving
    /// the block it held before, the stale read this module exists partly to avoid.
    #[test]
    fn seeking_to_the_end_of_the_stream_reads_nothing() {
        let bytes = bgzf_bytes(2, 1_000);
        let (_dir, path) = file_of(&bytes);
        let end_marker = *frame_offsets(&bytes).last().expect("frames");
        let end_of_file = VirtualPosition::try_from((bytes.len() as u64, 0)).expect("fits");

        let mut ours = SharedBgzfReader::new(shared(&path));
        let mut theirs = noodles_reader(&path);
        let at_marker = VirtualPosition::try_from((end_marker, 0)).expect("fits");
        ours.seek(at_marker).expect("ours seeks");
        theirs.seek(at_marker).expect("noodles seeks");
        assert_eq!(ours.virtual_position(), end_of_file);
        assert_reads_alike(&mut ours, &mut theirs);

        ours.seek(end_of_file)
            .expect("seeking to the end is allowed");
        assert_eq!(ours.virtual_position(), end_of_file);
        assert_eq!(ours.read(&mut [0u8; 16]).expect("reads"), 0);

        let past_the_end = VirtualPosition::try_from((bytes.len() as u64, 1)).expect("fits");
        assert_eq!(
            ours.seek(past_the_end)
                .expect_err("no byte 1 at the end")
                .kind(),
            io::ErrorKind::InvalidInput,
        );
    }

    /// `read_exact` across a block boundary takes the slow path and must return the same bytes
    /// as noodles; asking for more than the stream holds fails as noodles fails.
    #[test]
    fn read_exact_across_a_block_boundary_matches_noodles() {
        let (_dir, path) = file_of(&bgzf_bytes(3, 1_000));
        let mut ours = SharedBgzfReader::new(shared(&path));
        let mut theirs = noodles_reader(&path);
        let near_the_end_of_the_first = VirtualPosition::try_from((0, 900)).expect("fits");
        ours.seek(near_the_end_of_the_first).expect("ours seeks");
        theirs
            .seek(near_the_end_of_the_first)
            .expect("noodles seeks");

        // 100 bytes of the first block and 1,000 of the second and 200 of the third.
        let mut ours_bytes = [0u8; 1_300];
        let mut theirs_bytes = [0u8; 1_300];
        ours.read_exact(&mut ours_bytes).expect("ours reads");
        theirs.read_exact(&mut theirs_bytes).expect("noodles reads");
        assert_eq!(ours_bytes, theirs_bytes);
        assert_eq!(ours.virtual_position(), theirs.virtual_position());

        // 800 bytes are left; ask for 801.
        let ours_error = ours.read_exact(&mut [0u8; 801]).expect_err("ours runs out");
        let theirs_error = theirs
            .read_exact(&mut [0u8; 801])
            .expect_err("noodles runs out");
        assert_eq!(ours_error.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(ours_error.kind(), theirs_error.kind());
    }

    /// A block reached by stepping over an empty one is held once, not once per route to it:
    /// the cache keeps one slot per offset.
    #[test]
    fn a_block_reached_past_an_empty_one_is_held_once() {
        let first = bgzf_bytes(1, 500);
        let mut bytes = first.clone();
        bytes.extend(bgzf_bytes(1, 700));
        let (_dir, path) = file_of(&bytes);
        let empty_block = *frame_offsets(&first).last().expect("frames");
        let after_it = first.len() as u64;
        let blocks = shared(&path);

        let mut reader = SharedBgzfReader::new(Arc::clone(&blocks));
        reader
            .seek(VirtualPosition::try_from((after_it, 0)).expect("fits"))
            .expect("seeks");
        reader
            .seek(VirtualPosition::try_from((empty_block, 0)).expect("fits"))
            .expect("seeks");

        assert_eq!(
            reader.virtual_position(),
            VirtualPosition::try_from((after_it, 0)).expect("fits")
        );
        assert_eq!(blocks.held_offsets(), vec![after_it]);
    }

    /// **The point of the cache**: a second reader over the same file finds every block the
    /// first one inflated, and inflates none itself.
    #[test]
    fn a_second_reader_takes_the_blocks_the_first_inflated() {
        let (_dir, path) = file_of(&bgzf_bytes(5, 1_000));
        let blocks = shared(&path);
        let read_all = |blocks: &Arc<SharedBgzfBlocks>| {
            let mut everything = Vec::new();
            SharedBgzfReader::new(Arc::clone(blocks))
                .read_to_end(&mut everything)
                .expect("reads");
            everything
        };

        let first = read_all(&blocks);
        assert_eq!(blocks.counts().blocks_inflated, 5);
        let second = read_all(&blocks);

        assert_eq!(first, second);
        assert_eq!(
            blocks.counts(),
            BlockCacheCounts {
                blocks_inflated: 5,
                blocks_reused: 5,
            },
            "the second reader inflated blocks the first had already inflated",
        );
    }

    /// Past [`CACHE_SLOTS`] blocks the one given up is the one longest unasked: after reading
    /// one more block than there are slots, the first block has to be inflated again and the
    /// last one does not.
    #[test]
    fn the_least_recently_used_block_is_the_one_given_up() {
        let bytes = bgzf_bytes(CACHE_SLOTS + 1, 100);
        let (_dir, path) = file_of(&bytes);
        let blocks = shared(&path);
        let offsets = frame_offsets(&bytes);
        let mut reader = SharedBgzfReader::new(Arc::clone(&blocks));
        reader.read_to_end(&mut Vec::new()).expect("reads");
        let inflated = blocks.counts().blocks_inflated;
        assert_eq!(inflated, (CACHE_SLOTS + 1) as u64);

        let last = offsets[CACHE_SLOTS];
        reader
            .seek(VirtualPosition::try_from((last, 0)).expect("fits"))
            .expect("seeks");
        assert_eq!(
            blocks.counts().blocks_inflated,
            inflated,
            "the most recent block was given up"
        );

        reader
            .seek(VirtualPosition::try_from((offsets[0], 0)).expect("fits"))
            .expect("seeks");
        assert_eq!(
            blocks.counts().blocks_inflated,
            inflated + 1,
            "the oldest block was still held, so something else was given up in its place",
        );
    }

    /// Read `bytes` to the end with both readers and hand back both errors.
    fn both_errors(bytes: &[u8]) -> (io::Error, io::Error) {
        let (_dir, path) = file_of(bytes);
        let ours = SharedBgzfReader::new(shared(&path))
            .read_to_end(&mut Vec::new())
            .expect_err("ours must refuse this file");
        let theirs = noodles_reader(&path)
            .read_to_end(&mut Vec::new())
            .expect_err("noodles refuses this file");
        (ours, theirs)
    }

    fn assert_fails_alike(bytes: &[u8]) {
        let (ours, theirs) = both_errors(bytes);
        assert_eq!(ours.kind(), theirs.kind());
        assert_eq!(ours.to_string(), theirs.to_string());
    }

    /// **A corrupt block fails as it did before the cache**: a checksum that does not match the
    /// data, a damaged header and an impossible frame size each give the error noodles gives.
    #[test]
    fn a_corrupt_block_fails_as_noodles_fails() {
        let clean = bgzf_bytes(3, 1_000);
        let second_block = frame_offsets(&clean)[1] as usize;
        let third_block = frame_offsets(&clean)[2] as usize;

        // The CRC32 is the first half of the second block's trailer.
        let mut bad_checksum = clean.clone();
        bad_checksum[third_block - 8] ^= 0xff;
        assert_fails_alike(&bad_checksum);
        assert_eq!(
            both_errors(&bad_checksum).0.to_string(),
            "block data checksum mismatch"
        );

        // The gzip magic number.
        let mut bad_header = clean.clone();
        bad_header[second_block] = 0;
        assert_fails_alike(&bad_header);

        // A `BSIZE` smaller than a header and a trailer.
        let mut bad_size = clean.clone();
        bad_size[second_block + 16] = 4;
        bad_size[second_block + 17] = 0;
        assert_fails_alike(&bad_size);
    }

    /// A file cut short inside a frame is an error, as it is for noodles; one cut short between
    /// frames is simply the end of the stream, as it is for noodles.
    #[test]
    fn a_truncated_file_fails_where_noodles_fails() {
        let clean = bgzf_bytes(3, 1_000);
        let third_block = frame_offsets(&clean)[2] as usize;

        let mut cut_inside_a_frame = clean.clone();
        cut_inside_a_frame.truncate(third_block + 30);
        let (ours, theirs) = both_errors(&cut_inside_a_frame);
        assert_eq!(ours.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(ours.kind(), theirs.kind());

        let mut cut_between_frames = clean;
        cut_between_frames.truncate(third_block);
        let (_dir, path) = file_of(&cut_between_frames);
        assert_reads_alike(
            &mut SharedBgzfReader::new(shared(&path)),
            &mut noodles_reader(&path),
        );
    }

    /// An error leaves nothing behind: the corrupt block is not cached, so asking again fails
    /// again rather than serving what a failed inflate left.
    #[test]
    fn a_block_that_failed_to_inflate_is_not_kept() {
        let mut bytes = bgzf_bytes(2, 1_000);
        let second_block = frame_offsets(&bytes)[1] as usize;
        bytes[second_block - 8] ^= 0xff;
        let (_dir, path) = file_of(&bytes);
        let blocks = shared(&path);
        let position = VirtualPosition::try_from((0, 0)).expect("fits");
        for _ in 0..2 {
            let mut reader = SharedBgzfReader::new(Arc::clone(&blocks));
            assert!(reader.seek(position).is_err());
        }
        assert_eq!(blocks.counts(), BlockCacheCounts::default());
    }
}
