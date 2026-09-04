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
