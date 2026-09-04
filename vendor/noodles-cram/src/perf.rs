//! **Per-codec counters for a CRAM read**, added to this vendored copy so a run can say which
//! compression method its time went into. Not part of noodles: nothing upstream depends on it,
//! and it exists to answer "is a fork of the decompressors worth building?".
//!
//! Each counter is a relaxed atomic, so a threaded read stays correct while costing an
//! uncontended increment per block. The block count and the byte counts are exact; the
//! nanosecond figures carry one `Instant::now()` pair per block.

use std::sync::atomic::{AtomicU64, Ordering};

/// One row per CRAM compression method, in the numbering of the CRAM v3 specification.
pub const METHOD_NAMES: [&str; 9] = [
    "none",
    "gzip",
    "bzip2",
    "lzma",
    "rans_4x8",
    "rans_nx16",
    "arithmetic",
    "fqzcomp",
    "name_tokenizer",
];

const METHODS: usize = 9;

macro_rules! counters {
    ($name:ident) => {
        static $name: [AtomicU64; METHODS] = [
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
        ];
    };
}

counters!(BLOCKS);
counters!(COMPRESSED_BYTES);
counters!(UNCOMPRESSED_BYTES);
counters!(NANOS);

pub(crate) fn record(method: usize, compressed: usize, uncompressed: usize, nanos: u64) {
    BLOCKS[method].fetch_add(1, Ordering::Relaxed);
    COMPRESSED_BYTES[method].fetch_add(compressed as u64, Ordering::Relaxed);
    UNCOMPRESSED_BYTES[method].fetch_add(uncompressed as u64, Ordering::Relaxed);
    NANOS[method].fetch_add(nanos, Ordering::Relaxed);
}

/// What one compression method cost over the blocks decoded so far.
#[derive(Clone, Copy, Debug, Default)]
pub struct MethodCost {
    pub blocks: u64,
    pub compressed_bytes: u64,
    pub uncompressed_bytes: u64,
    pub nanos: u64,
}

/// The counters as they stand, one row per method, indexed as [`METHOD_NAMES`].
pub fn read_counters() -> [MethodCost; METHODS] {
    std::array::from_fn(|method| MethodCost {
        blocks: BLOCKS[method].load(Ordering::Relaxed),
        compressed_bytes: COMPRESSED_BYTES[method].load(Ordering::Relaxed),
        uncompressed_bytes: UNCOMPRESSED_BYTES[method].load(Ordering::Relaxed),
        nanos: NANOS[method].load(Ordering::Relaxed),
    })
}

/// Zero every counter, so a timed pass measures only itself.
pub fn reset_counters() {
    for method in 0..METHODS {
        BLOCKS[method].store(0, Ordering::Relaxed);
        COMPRESSED_BYTES[method].store(0, Ordering::Relaxed);
        UNCOMPRESSED_BYTES[method].store(0, Ordering::Relaxed);
        NANOS[method].store(0, Ordering::Relaxed);
    }
}

// ---------------------------------------------------------------------
// rANS 4x8: order, and the split between building the decode tables and decoding
// ---------------------------------------------------------------------

static RANS_ORDER0_BLOCKS: AtomicU64 = AtomicU64::new(0);
static RANS_ORDER1_BLOCKS: AtomicU64 = AtomicU64::new(0);
static RANS_TABLE_NANOS: AtomicU64 = AtomicU64::new(0);
static RANS_DECODE_NANOS: AtomicU64 = AtomicU64::new(0);
static RANS_ORDER0_BYTES: AtomicU64 = AtomicU64::new(0);
static RANS_ORDER1_BYTES: AtomicU64 = AtomicU64::new(0);

pub(crate) fn record_rans(order_one: bool, bytes: usize, table_nanos: u64, decode_nanos: u64) {
    if order_one {
        RANS_ORDER1_BLOCKS.fetch_add(1, Ordering::Relaxed);
        RANS_ORDER1_BYTES.fetch_add(bytes as u64, Ordering::Relaxed);
    } else {
        RANS_ORDER0_BLOCKS.fetch_add(1, Ordering::Relaxed);
        RANS_ORDER0_BYTES.fetch_add(bytes as u64, Ordering::Relaxed);
    }
    RANS_TABLE_NANOS.fetch_add(table_nanos, Ordering::Relaxed);
    RANS_DECODE_NANOS.fetch_add(decode_nanos, Ordering::Relaxed);
}

/// How the rANS 4x8 blocks split by order, and where their time went.
#[derive(Clone, Copy, Debug, Default)]
pub struct RansCost {
    pub order_0_blocks: u64,
    pub order_1_blocks: u64,
    pub order_0_bytes: u64,
    pub order_1_bytes: u64,
    pub table_nanos: u64,
    pub decode_nanos: u64,
}

pub fn read_rans_counters() -> RansCost {
    RansCost {
        order_0_blocks: RANS_ORDER0_BLOCKS.load(Ordering::Relaxed),
        order_1_blocks: RANS_ORDER1_BLOCKS.load(Ordering::Relaxed),
        order_0_bytes: RANS_ORDER0_BYTES.load(Ordering::Relaxed),
        order_1_bytes: RANS_ORDER1_BYTES.load(Ordering::Relaxed),
        table_nanos: RANS_TABLE_NANOS.load(Ordering::Relaxed),
        decode_nanos: RANS_DECODE_NANOS.load(Ordering::Relaxed),
    }
}

pub fn reset_rans_counters() {
    RANS_ORDER0_BLOCKS.store(0, Ordering::Relaxed);
    RANS_ORDER1_BLOCKS.store(0, Ordering::Relaxed);
    RANS_ORDER0_BYTES.store(0, Ordering::Relaxed);
    RANS_ORDER1_BYTES.store(0, Ordering::Relaxed);
    RANS_TABLE_NANOS.store(0, Ordering::Relaxed);
    RANS_DECODE_NANOS.store(0, Ordering::Relaxed);
}

// ---------------------------------------------------------------------
// Slice::records — where the record-decode layer's time goes
// ---------------------------------------------------------------------

static SLICE_ALLOC_NANOS: AtomicU64 = AtomicU64::new(0);
static SLICE_REFERENCE_NANOS: AtomicU64 = AtomicU64::new(0);
static SLICE_READ_NANOS: AtomicU64 = AtomicU64::new(0);
static SLICE_TAG_NANOS: AtomicU64 = AtomicU64::new(0);
static SLICE_MATE_NANOS: AtomicU64 = AtomicU64::new(0);
static SLICE_RECORDS: AtomicU64 = AtomicU64::new(0);

pub(crate) fn record_slice_phase(phase: SlicePhase, nanos: u64) {
    match phase {
        SlicePhase::Allocate => &SLICE_ALLOC_NANOS,
        SlicePhase::Reference => &SLICE_REFERENCE_NANOS,
        SlicePhase::ReadRecords => &SLICE_READ_NANOS,
        SlicePhase::ReadTags => &SLICE_TAG_NANOS,
        SlicePhase::ResolveMates => &SLICE_MATE_NANOS,
    }
    .fetch_add(nanos, Ordering::Relaxed);
}

pub(crate) fn record_slice_records(records: usize) {
    SLICE_RECORDS.fetch_add(records as u64, Ordering::Relaxed);
}

/// One stage of turning a decoded slice into records.
#[derive(Clone, Copy, Debug)]
pub enum SlicePhase {
    /// Building the `Vec<Record>` the slice's records are read into.
    Allocate,
    /// Fetching the slice's reference sequence, and verifying its digest.
    Reference,
    /// Reading every record out of the decoded blocks — tags excluded.
    ReadRecords,
    /// The auxiliary tags, within the above.
    ReadTags,
    /// Linking each record to its mate.
    ResolveMates,
}

/// What each stage of `Slice::records` cost over the slices decoded so far.
#[derive(Clone, Copy, Debug, Default)]
pub struct SliceCost {
    pub records: u64,
    pub allocate_nanos: u64,
    pub reference_nanos: u64,
    pub read_nanos: u64,
    pub tag_nanos: u64,
    pub resolve_mates_nanos: u64,
}

pub fn read_slice_counters() -> SliceCost {
    SliceCost {
        records: SLICE_RECORDS.load(Ordering::Relaxed),
        allocate_nanos: SLICE_ALLOC_NANOS.load(Ordering::Relaxed),
        reference_nanos: SLICE_REFERENCE_NANOS.load(Ordering::Relaxed),
        read_nanos: SLICE_READ_NANOS.load(Ordering::Relaxed),
        tag_nanos: SLICE_TAG_NANOS.load(Ordering::Relaxed),
        resolve_mates_nanos: SLICE_MATE_NANOS.load(Ordering::Relaxed),
    }
}

pub fn reset_slice_counters() {
    SLICE_RECORDS.store(0, Ordering::Relaxed);
    SLICE_ALLOC_NANOS.store(0, Ordering::Relaxed);
    SLICE_REFERENCE_NANOS.store(0, Ordering::Relaxed);
    SLICE_READ_NANOS.store(0, Ordering::Relaxed);
    SLICE_TAG_NANOS.store(0, Ordering::Relaxed);
    SLICE_MATE_NANOS.store(0, Ordering::Relaxed);
}
