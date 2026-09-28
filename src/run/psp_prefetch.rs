//! **Reading the cohort's psps ahead of the merge, on a thread of its own**, so that the disk
//! is working while the pool is computing.
//!
//! Without this a calling run over stored files alternates. The merge covers a round — every
//! sample draws its next records, which on a cold file is a read that waits on the disk — and
//! only then builds and genotypes it, while the disk sits idle. Measured on 2,169 tomato
//! accessions, psps on one USB spinning disk, 32 regions of 500 bases in flight: about six
//! seconds reading at 90–115 MB/s, then five to seven with the disk doing nothing, over and
//! over. Half the wall clock was one resource waiting for the other.
//!
//! **What this does is fill the operating system's page cache, and nothing else.** A thread
//! reads each file's bytes a stretch of genome ahead of where the merge is working and throws
//! them away; when the merge's own reader reaches them they come from memory. The merge, the
//! readers and every byte of output are untouched — a run with this off and one with it on
//! read the same bytes in the same order through the same code, so the VCF cannot differ.
//!
//! **It is a hint, and it fails as one.** A file that will not open or read here is left
//! alone for the rest of the run; the merge's own reader meets the same fault and reports it
//! with the file's name, which is where a failure belongs.
//!
//! ## How far ahead, and what that costs
//!
//! The cost is page cache: memory the kernel gives back under pressure, not memory this process
//! holds. It is bounded by a byte budget ([`DEFAULT_PSP_PREFETCH_BUDGET_BYTES`],
//! `--psp-prefetch-bytes`), **shared between the files in proportion to their size**, so each
//! file's share covers about the same stretch of genome whatever its depth — on the tomato
//! cohort the files run from 4 GB to 40 GB over the same genome.
//!
//! **Bytes and not a stretch of genome**, which the first version used, turning the budget into
//! bases by the cohort's average bytes a base. On the tomato cohort, with an 8 GiB budget, that
//! version read about 44 GB in its first five minutes and the merge got only 17% further than
//! with no reading ahead at all: the look-ahead held far more than the budget, and the thread's
//! reading competed with the merge's own. The genome is not evenly dense, so an average turns a
//! byte budget into a stretch that can hold many times it. A byte count cannot overshoot.
//!
//! **Each file is refilled only once what is read ahead of the merge has fallen to half its
//! share**, and then read up to the whole of it. Topping every file up by one round's worth
//! would be a round's few tens of kilobytes per file — the seek-bound pattern this exists to
//! replace. Waiting for half means every refill is half a file's share, in one sequential run.
//!
//! What the range of cohorts gets from the default of 8 GiB:
//!
//! - **thousands of samples**: 2,169 files share it, about 3.8 MB each, so a refill is about
//!   1.9 MB of a file. On a spinning disk more is faster; the flag is how.
//! - **one sample**: the whole budget is that file's, read 4 GiB at a time.
//! - **any depth**: a deeper file is larger, so it gets a larger share and the stretch of genome
//!   stays about the same.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;

use crate::psp::PspReader;
use crate::types::GenomePosition;

/// The page cache a run may fill ahead of the merge when nobody says otherwise: 8 GiB.
///
/// **Sized for the seek, not for the memory.** A refill is half a file's share, and on a
/// spinning disk a read shorter than a few megabytes spends more of its time moving the head
/// than transferring. At 2,169 files 8 GiB gives refills of about 1.9 MB; a cohort of that size
/// on a machine with memory to spare does better with more (`--psp-prefetch-bytes`).
pub const DEFAULT_PSP_PREFETCH_BUDGET_BYTES: u64 = 8 << 30;

/// The coarsest and finest spacing, in file bytes, of the positions this thread keeps from each
/// file's block index.
///
/// **A coarse copy, because the whole index is the reader's and the reader is the merge's.**
/// Copying every entry would duplicate the index — at 2,169 files that is gigabytes. A file
/// keeps one entry per eighth of its share, within these bounds: fine enough that where the
/// merge has got to in the file is known to within an eighth of what is read ahead of it, and at
/// the 8 GiB default about 17 million entries, 400 MB, across the tomato cohort's 8.6 TB.
const MARK_SPACING_BYTES: std::ops::RangeInclusive<u64> = (256 << 10)..=(4 << 20);

/// How much one `read` asks for. Large enough that the kernel sees a sequential stream, small
/// enough that a stop request is noticed within one of them.
const READ_BYTES: usize = 1 << 20;

/// One file this thread reads ahead in.
struct FileAhead {
    path: PathBuf,
    /// Block first positions and byte offsets, one every eighth of `share` or so, in file order
    /// (which is genome order).
    marks: Vec<(GenomePosition, u64)>,
    /// Where the blocks end — nothing past it is a block.
    blocks_end: u64,
    /// This file's part of the budget: how many bytes past the merge it is read to.
    share: u64,
    /// The byte offset read to so far.
    read_to: u64,
    /// Set by the first fault; the file is left alone after it.
    given_up: bool,
}

impl FileAhead {
    fn of(path: PathBuf, psp: &PspReader, share: u64) -> Self {
        let spacing = (share / 8).clamp(*MARK_SPACING_BYTES.start(), *MARK_SPACING_BYTES.end());
        let mut marks: Vec<(GenomePosition, u64)> = Vec::new();
        for entry in psp.block_index() {
            let spaced = marks
                .last()
                .is_none_or(|&(_, at)| entry.block_offset >= at + spacing);
            if spaced {
                marks.push((entry.first_position, entry.block_offset));
            }
        }
        Self {
            path,
            marks,
            blocks_end: psp.footer().index_offset,
            share,
            read_to: 0,
            given_up: false,
        }
    }

    /// The offset of the last mark at or before `at` — where the block holding `at` begins or
    /// earlier — or the first block when `at` is before all of them.
    fn offset_holding(&self, at: GenomePosition) -> u64 {
        let after = self.marks.partition_point(|&(first, _)| first <= at);
        match after.checked_sub(1) {
            Some(mark) => self.marks[mark].1,
            None => self.marks.first().map_or(self.blocks_end, |&(_, at)| at),
        }
    }

    /// Where to read next, as `[from, to)`: nothing while more than half the share is read
    /// ahead of `merge_at`, and otherwise from where reading stopped up to a whole share past it.
    fn range_to_refill(&self, merge_at: GenomePosition) -> Option<(u64, u64)> {
        if self.given_up {
            return None;
        }
        let merge_from = self.offset_holding(merge_at);
        if self.read_to >= merge_from.saturating_add(self.share / 2) {
            return None;
        }
        let from = self.read_to.max(merge_from);
        let to = merge_from.saturating_add(self.share).min(self.blocks_end);
        (from < to).then_some((from, to))
    }

    /// Refill this file if it has fallen to half its share.
    fn refill(&mut self, merge_at: GenomePosition, buffer: &mut [u8], stop: &AtomicBool) {
        let Some((from, to)) = self.range_to_refill(merge_at) else {
            return;
        };
        if self.read_range(from, to, buffer, stop).is_err() {
            self.given_up = true;
        }
    }

    /// Read `[from, to)` and discard it. **Opened for each refill and closed after**, so the
    /// thread holds no descriptor between refills and the run's descriptor count is the one
    /// [`OpenPspCohort`](super::psp_caller::OpenPspCohort) checked.
    fn read_range(
        &mut self,
        from: u64,
        to: u64,
        buffer: &mut [u8],
        stop: &AtomicBool,
    ) -> std::io::Result<()> {
        let mut file = File::open(&self.path)?;
        file.seek(SeekFrom::Start(from))?;
        let mut at = from;
        while at < to && !stop.load(Ordering::Relaxed) {
            let want = usize::try_from(to - at).map_or(buffer.len(), |left| left.min(buffer.len()));
            let got = file.read(&mut buffer[..want])?;
            if got == 0 {
                break;
            }
            at += got as u64;
        }
        self.read_to = at;
        Ok(())
    }
}

/// What the merge and the reading thread share: where the merge is, and whether to stop.
struct Shared {
    merge_at: Mutex<Option<GenomePosition>>,
    moved: Condvar,
    stop: AtomicBool,
}

/// **The reading thread, running for as long as this value lives.** Dropping it stops the
/// thread and waits for it — within one [`READ_BYTES`] read.
pub struct PspPrefetch {
    shared: Arc<Shared>,
    thread: Option<JoinHandle<()>>,
}

impl PspPrefetch {
    /// Start reading `psps` ahead, filling at most about `budget_bytes` of page cache.
    ///
    /// **`None` when there is nothing to do**: a budget of zero, which is how the flag turns
    /// this off, or a cohort whose files hold no blocks.
    #[must_use]
    pub fn start(psps: &[PspReader], paths: &[PathBuf], budget_bytes: u64) -> Option<Self> {
        let block_bytes: Vec<u64> = psps
            .iter()
            .map(|psp| {
                let first = psp
                    .block_index()
                    .first()
                    .map_or(0, |entry| entry.block_offset);
                psp.footer().index_offset.saturating_sub(first)
            })
            .collect();
        let total: u64 = block_bytes.iter().sum();
        if budget_bytes == 0 || total == 0 {
            return None;
        }
        // **In proportion to the file's size**, and at least two reads' worth, so a share too
        // small to be worth a seek still reads something.
        let share_of = |bytes: u64| {
            u64::try_from(u128::from(budget_bytes) * u128::from(bytes) / u128::from(total))
                .unwrap_or(u64::MAX)
                .max(2 * READ_BYTES as u64)
        };
        let files: Vec<FileAhead> = psps
            .iter()
            .zip(paths)
            .zip(block_bytes)
            .map(|((psp, path), bytes)| FileAhead::of(path.clone(), psp, share_of(bytes)))
            .collect();
        let shared = Arc::new(Shared {
            merge_at: Mutex::new(None),
            moved: Condvar::new(),
            stop: AtomicBool::new(false),
        });
        let for_the_thread = Arc::clone(&shared);
        let thread = std::thread::Builder::new()
            .name("psp-prefetch".to_owned())
            .spawn(move || read_ahead(files, &for_the_thread))
            .ok()?;
        Some(Self {
            shared,
            thread: Some(thread),
        })
    }

    /// Tell the thread the merge has reached `at`. Cheap: one lock and one wake.
    pub fn merge_is_at(&self, at: GenomePosition) {
        let mut merge_at = self
            .shared
            .merge_at
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if merge_at.is_none_or(|was| was < at) {
            *merge_at = Some(at);
            self.shared.moved.notify_one();
        }
    }
}

impl Drop for PspPrefetch {
    fn drop(&mut self) {
        self.shared.stop.store(true, Ordering::Relaxed);
        {
            let _held = self
                .shared
                .merge_at
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            self.shared.moved.notify_one();
        }
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

/// The thread's body: wait for the merge to move, then refill every file that has fallen to
/// half its share, in sample order, one sequential run per file.
fn read_ahead(mut files: Vec<FileAhead>, shared: &Shared) {
    let mut buffer = vec![0_u8; READ_BYTES];
    let mut served: Option<GenomePosition> = None;
    loop {
        let merge_at = {
            let mut merge_at = shared
                .merge_at
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            while !shared.stop.load(Ordering::Relaxed) && *merge_at == served {
                merge_at = shared
                    .moved
                    .wait(merge_at)
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
            }
            if shared.stop.load(Ordering::Relaxed) {
                return;
            }
            *merge_at
        };
        served = merge_at;
        let Some(merge_at) = merge_at else { continue };
        for file in &mut files {
            if shared.stop.load(Ordering::Relaxed) {
                return;
            }
            file.refill(merge_at, &mut buffer, &shared.stop);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{ContigId, Position};

    fn at(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    /// A file with a mark every 100 bases, each 10 bytes on from the last, blocks ending at 50,
    /// and a share of 20 bytes.
    fn five_marks() -> FileAhead {
        FileAhead {
            path: PathBuf::new(),
            marks: (0..5)
                .map(|mark| (at(0, 1 + 100 * mark), 10 * mark))
                .collect(),
            blocks_end: 50,
            share: 20,
            read_to: 0,
            given_up: false,
        }
    }

    /// **Where the merge is in a file is the mark at or before it**: the block holding the
    /// merge's position begins there or later, so nothing the merge still needs is skipped.
    #[test]
    fn the_merge_is_placed_at_the_mark_at_or_before_it() {
        let file = five_marks();
        assert_eq!(file.offset_holding(at(0, 150)), 10);
        assert_eq!(file.offset_holding(at(0, 101)), 10);
        assert_eq!(file.offset_holding(at(0, 0)), 0);
        // A later contig is past every mark on this one.
        assert_eq!(file.offset_holding(at(1, 1)), 40);
    }

    /// **A refill reads up to a whole share past the merge, and only once half of it is left**,
    /// which is what keeps every refill one long read rather than a round's worth of short ones.
    #[test]
    fn a_file_is_refilled_only_once_it_has_fallen_to_half_its_share() {
        let mut file = five_marks();
        // Nothing read yet: a whole share from the merge's mark.
        assert_eq!(file.range_to_refill(at(0, 150)), Some((10, 30)));
        // Read to 30 with the merge at 10: 20 ahead, more than half of 20, so nothing.
        file.read_to = 30;
        assert_eq!(file.range_to_refill(at(0, 150)), None);
        // The merge at 20: 10 ahead, which is half, so still nothing.
        assert_eq!(file.range_to_refill(at(0, 250)), None);
        // The merge at 30: nothing ahead, so from where reading stopped to a share past it.
        assert_eq!(file.range_to_refill(at(0, 350)), Some((30, 50)));
        // Never past the blocks.
        file.read_to = 45;
        assert_eq!(file.range_to_refill(at(0, 450)), Some((45, 50)));
        file.read_to = 50;
        assert_eq!(file.range_to_refill(at(0, 450)), None);
    }

    /// **A file that failed once is left alone.**
    #[test]
    fn a_file_that_failed_is_not_read_again() {
        let mut file = five_marks();
        file.given_up = true;
        assert_eq!(file.range_to_refill(at(0, 150)), None);
    }
}
