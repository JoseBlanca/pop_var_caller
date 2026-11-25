//! **mimalloc's options as `pvc` runs it.** One setting: freed memory goes back to the system at
//! once instead of a second later.
//!
//! A crate of its own because it calls mimalloc's C functions, which needs `unsafe`, and the main
//! crate forbids it.

use libmimalloc_sys::{mi_option_get, mi_option_set_default, mi_option_t};

/// mimalloc's `mi_option_purge_delay`: how long, in milliseconds, freed memory is kept before it
/// is returned to the system. Its place in mimalloc's `mi_option_e` (v3 and v2 alike);
/// `libmimalloc-sys` does not export it as a constant.
const PURGE_DELAY: mi_option_t = 15;

/// **Return freed memory to the system at once** (mimalloc's purge delay 0 ms), unless the
/// environment says otherwise.
///
/// mimalloc 3 keeps memory it has been given back for 1,000 ms before returning it (`purge_delay`;
/// mimalloc 2 kept it 10 ms). Where a program frees a large structure and asks for another at
/// once — refinement does at the end of every window's rounds, freeing the haplotype store as it
/// builds the core's results — the process holds both for that second. Measured on refinement
/// alone (`doc/devel/reports/implementations/refine_memory_overhead_2026-11-25.md` §4, Linux,
/// 6 threads): resident memory 17 to 24% above the bytes the program held at the default, 6 to 7%
/// at 0 ms, the output the same. On a whole `pvc call-from-alignments` run — 63 tomato samples
/// over 8 Mb, refinement on, 6 threads, Linux — the peak resident memory fell from 6.8–7.0 GB to
/// 4.9 GB, for 3% more CPU time (most of it the kernel's, returning pages) and 2% more
/// wall time, the VCF byte for byte the same
/// (`doc/devel/reports/implementations/mimalloc_purge_delay_2026-11-25.md`).
///
/// **`MIMALLOC_PURGE_DELAY` set in the environment wins.** mimalloc reads its environment when
/// the process starts, before `main`, and marks an option it found there as set; this only
/// changes the default, which leaves such an option alone.
///
/// Call it first thing in `main`, before other threads start: mimalloc's option table is not
/// guarded against concurrent writes.
pub fn return_freed_memory_at_once() {
    // SAFETY: `mi_option_set_default` writes one entry of mimalloc's option table, in range, and
    // has no other precondition; called before any other thread exists, nothing reads it
    // concurrently.
    unsafe { mi_option_set_default(PURGE_DELAY, 0) }
}

/// The purge delay mimalloc uses now, in milliseconds.
pub fn purge_delay_ms() -> std::ffi::c_long {
    // SAFETY: reads one entry of mimalloc's option table, in range.
    unsafe { mi_option_get(PURGE_DELAY) }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The option index is the purge delay's: mimalloc 3's default is 1,000 ms before the call
    /// and 0 after it. Skipped if the environment sets the delay, which the call must not touch.
    #[test]
    fn sets_the_purge_delay_unless_the_environment_does() {
        let from_environment = std::env::var_os("MIMALLOC_PURGE_DELAY").is_some();
        let before = purge_delay_ms();
        return_freed_memory_at_once();
        if from_environment {
            assert_eq!(purge_delay_ms(), before);
        } else {
            assert_eq!(before, 1000, "not mimalloc 3's purge delay default");
            assert_eq!(purge_delay_ms(), 0);
        }
    }
}
