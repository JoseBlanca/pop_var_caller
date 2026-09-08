//! **Decoding one CRAM's containers a little ahead of the walk, on a thread of its own.**
//!
//! # What it is for
//!
//! Reading a CRAM record is not a step, it is a step that sometimes stops to inflate three and a
//! half megabytes. The walk pulls reads one at a time, and every ten thousandth of those calls
//! reaches the bottom of the stack and decodes a whole container before it can return the one
//! read that was asked for. Profiled on 10 Mb of tomato chromosome 1, that is **4.9 seconds of a
//! 41.4-second walk**, and the call tree shows exactly where it sits — under
//! `PileupGenerator::next_locus`, not beside it:
//!
//! ```text
//! GeneratorSet::next_locus            "give me the next position"
//! └ PreparedSampleReads::next           … which needs another read
//!   └ AlignmentCursor::next_read
//!     └ RegionRawAlignedReads::read_next
//!       └ Slice::decode_blocks_inner    ← a whole container, decoded here
//! ```
//!
//! Nothing about building a position needs the *next* container, and which container is next is
//! knowable: the `.crai` lists them in file order and the walk goes forward through it. So a
//! thread can decode it while the generator is still working through the one before, and put it
//! where the readers already look — the shared [`SharedContainers`] cache.
//!
//! # Why this is a hint and not a queue
//!
//! The prefetcher holds **one wanted offset**, not a backlog. A reader that takes a container
//! overwrites the hint with the next offset it expects to want, and if the prefetcher had not yet
//! picked up the previous hint, that hint is simply replaced. This is the right shape rather than
//! a lazy one: work the walk has already passed is worthless, and a queue would have to be drained
//! of it. One cell can hold nothing stale.
//!
//! # What it may not do
//!
//! **It may not decode a container twice.** The cache's claim mechanism is what enforces that:
//! the prefetcher takes a claim exactly as a reader does, and a reader arriving at a container
//! the prefetcher is midway through waits for it rather than starting a second decode. Without
//! that, a prefetcher that fell behind would double the decode it was meant to hide.
//!
//! **It may not evict what a reader is about to want.** That is why the cache has four slots
//! rather than two — two for the readers, two to run ahead into.
//!
//! **It may not outlive its file.** The thread is joined when the [`AlignmentFile`] is dropped,
//! so a cohort walked in one process leaves none behind.
//!
//! [`AlignmentFile`]: crate::ng::read::input::AlignmentFile

use std::fs::File;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};

use noodles_cram as cram;
use noodles_sam as sam;

use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::ref_seq::RawRefSeq;

use super::container::{SharedContainers, container_at};

/// **What an open CRAM lends every reader it mints**: the decoded containers they share, and
/// the thread that decodes the next one while they work.
///
/// One value rather than two arguments because they are one thing — a reader that had the cache
/// and not the prefetcher would take from a cache nothing fills ahead, and one with the
/// prefetcher and a cache of its own would have the thread filling a cache it never reads.
///
/// **The prefetcher is not inside [`SharedContainers`], and that is deliberate**: the thread
/// holds the cache, so a cache that held the thread would be a cycle and the thread would never
/// be joined.
#[derive(Clone)]
pub(crate) struct SharedDecode {
    /// The containers this file's readers take from.
    pub(crate) containers: Arc<SharedContainers>,
    /// The thread that fills them ahead of the walk, where there is one.
    pub(crate) prefetcher: Option<Arc<ContainerPrefetcher>>,
}

/// **Neither the cache nor the thread is printable, and a derived `Debug` would try.** What
/// identifies this is whether a prefetcher is running, which is the only thing a reader's own
/// `Debug` could usefully say about it.
impl std::fmt::Debug for SharedDecode {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SharedDecode")
            .field("prefetching", &self.prefetcher.is_some())
            .finish_non_exhaustive()
    }
}

/// The offset a reader would like decoded next, and the thread that decodes it.
///
/// See the module documentation for what it is for and what it may not do.
pub(crate) struct ContainerPrefetcher {
    wanted: Arc<Wanted>,
    /// Joined on drop, so no thread outlives the file it reads.
    worker: Option<std::thread::JoinHandle<()>>,
}

/// The one-cell hint, and the flag that ends the thread.
struct Wanted {
    offset: Mutex<Option<u64>>,
    changed: Condvar,
    stopping: AtomicBool,
}

impl std::fmt::Debug for ContainerPrefetcher {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ContainerPrefetcher")
            .finish_non_exhaustive()
    }
}

impl ContainerPrefetcher {
    /// Start the thread. It owns `reader` and `reference` — **its own descriptor and its own
    /// reference reader**, because both are positions rather than values and cannot be shared
    /// with a reader that is seeking elsewhere.
    pub(crate) fn start(
        shared: Arc<SharedContainers>,
        mut reader: cram::io::Reader<File>,
        header: Arc<sam::Header>,
        resolution: ReadGroupResolution,
        reference: Box<dyn RawRefSeq + Send>,
    ) -> Self {
        let wanted = Arc::new(Wanted {
            offset: Mutex::new(None),
            changed: Condvar::new(),
            stopping: AtomicBool::new(false),
        });
        let mine = Arc::clone(&wanted);
        let worker = std::thread::Builder::new()
            .name("cram-container-prefetch".to_string())
            .spawn(move || {
                while let Some(offset) = mine.next_wanted() {
                    // **The result is dropped on purpose.** A prefetch that fails is not this
                    // thread's to report: the reader that reaches the same container will decode
                    // it itself and raise the same failure where it can name the walk's position.
                    // Storing it in the cache is the whole of the success case too — the handle
                    // returned here is a second reference to what the cache now holds.
                    let _ = container_at(
                        &shared,
                        &mut reader,
                        &header,
                        &resolution,
                        offset,
                        &*reference,
                    );
                }
            })
            .expect("a thread for prefetching CRAM containers");
        Self {
            wanted,
            worker: Some(worker),
        }
    }

    /// Ask for `offset` to be decoded, replacing whatever was wanted before.
    ///
    /// Never blocks and never fails: a hint the thread does not get to is a hint the walk has
    /// overtaken, and losing it costs nothing.
    pub(crate) fn wants(&self, offset: u64) {
        let mut wanted = self
            .wanted
            .offset
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        *wanted = Some(offset);
        drop(wanted);
        self.wanted.changed.notify_one();
    }
}

impl Wanted {
    /// The next offset to decode, or `None` once the file is being dropped.
    fn next_wanted(&self) -> Option<u64> {
        let mut wanted = self
            .offset
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        loop {
            if self.stopping.load(Ordering::Acquire) {
                return None;
            }
            if let Some(offset) = wanted.take() {
                return Some(offset);
            }
            wanted = self
                .changed
                .wait(wanted)
                .unwrap_or_else(|poisoned| poisoned.into_inner());
        }
    }
}

/// **Joined, not detached.** A file that is closed must leave no thread behind — a cohort walked
/// in one process would otherwise accumulate one per sample, each holding a CRAM descriptor and a
/// reference reader.
impl Drop for ContainerPrefetcher {
    fn drop(&mut self) {
        self.wanted.stopping.store(true, Ordering::Release);
        self.wanted.changed.notify_all();
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}
