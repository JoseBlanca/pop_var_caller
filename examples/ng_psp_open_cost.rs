//! **What does one open psp cost in resident memory, and how does that cost grow with the
//! cohort?** — Milestone H4's measurement.
//!
//! ```text
//! cargo run --release --example ng_psp_open_cost -- <a store.ngpsp> --samples N [--rounds R]
//! ```
//!
//! # The number this is for
//!
//! The format is shaped around one requirement (spec `psp_file_format.md` §1.1): **an open file
//! costs no more than 500 kB of resident memory, and that figure does not grow with the block
//! size, the read depth or the length of the genome.** A cohort holds one open file per sample
//! for the whole run, so that cost is multiplied by the cohort size — which is why the plan asks
//! for the **slope**, not one point. A single reading at one sample count says nothing: it is
//! dominated by the process, the binary and the allocator's own arena.
//!
//! So this program opens the same store `N` times, advances every reader one record per round in
//! lockstep — the shape a cohort merge reads in — and reports the peak resident set. Run it at
//! several `N` and the fixed part cancels in the differences; the slope is the per-open-sample
//! cost.
//!
//! **⚠ Run it once per `N`, in a fresh process, and that is not a convenience.** The peak is read
//! from `VmHWM` in `/proc/self/status`, which is a **high-water mark for the life of the
//! process**. Two sample counts measured in one process would both report the larger, and the
//! slope would come back zero — a false pass on the one claim the format exists to make.
//!
//! # What it measures, and what it does not
//!
//! It opens **the same store N times**, which is what the measuring prototype did for the same
//! question: a reader's buffers are the same size whatever they carry, so this gives the
//! multiplier without needing the samples to differ. It is therefore a fact about the reader's
//! own cost and **not** about the variety of a real cohort.
//!
//! Linux only: `VmHWM` is a `/proc` field. The dev container is Linux and so is CI.

use std::path::PathBuf;

use pop_var_caller::ng::psp::{PspReader, RecordIter};

/// How many records each reader is advanced by default.
///
/// **Past the first block, which is what matters.** A reader's buffers reach their full size while
/// its first block inflates, so a round count in the hundreds measures steady state; more only
/// costs time. The figure is reported, so a reading is never quoted without it.
const DEFAULT_ROUNDS: usize = 1_000;

/// The peak resident set of this process, in kilobytes, from `VmHWM`.
///
/// **A high-water mark, not the current size** — which is what makes it the right number for a
/// peak and the wrong number to read twice in one process.
fn peak_resident_kib() -> Option<u64> {
    kib_from_status(&std::fs::read_to_string("/proc/self/status").ok()?)
}

/// The `VmHWM` line's kilobyte figure, split out from the file read so it can be tested.
fn kib_from_status(status: &str) -> Option<u64> {
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")?
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    })
}

fn main() {
    let mut args = std::env::args().skip(1);
    let store = PathBuf::from(
        args.next()
            .expect("usage: ng_psp_open_cost <a store.ngpsp> --samples N [--rounds R]"),
    );
    let mut samples = 1usize;
    let mut rounds = DEFAULT_ROUNDS;
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--samples" => {
                let value = args.next().expect("--samples needs a cohort size");
                samples = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--samples {value:?} is not a count: {why}"));
            }
            "--rounds" => {
                let value = args.next().expect("--rounds needs a record count");
                rounds = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--rounds {value:?} is not a count: {why}"));
            }
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(samples > 0, "--samples needs at least one open sample");

    let before_opening = peak_resident_kib().expect("this measurement is Linux-only");

    // **Opened first, all of them, before anything is walked.** That is the cost spec §6.2 says a
    // cohort pays per sample up front: the footer, the index and the header, and no block.
    let mut readers: Vec<PspReader> = (0..samples)
        .map(|which| {
            PspReader::open(&store).unwrap_or_else(|why| {
                panic!("opening {} for sample {which}: {why}", store.display())
            })
        })
        .collect();
    let after_opening = peak_resident_kib().expect("VmHWM");

    // **Disjoint mutable borrows, one walk per reader**, which is what lets a cohort hold them all
    // at once: `iter_mut` hands out one borrow per element, so N walks live together.
    let mut walks: Vec<RecordIter<'_>> = readers
        .iter_mut()
        .map(|reader| reader.records().expect("walking the store"))
        .collect();

    let mut records = 0u64;
    let mut checksum = 0u64;
    for _ in 0..rounds {
        let mut any = false;
        for walk in &mut walks {
            if let Some(found) = walk.next() {
                let found = found.expect("a record");
                records += 1;
                checksum = checksum.wrapping_add(found.head.region.start.0);
                any = true;
            }
        }
        if !any {
            break;
        }
    }
    let after_walking = peak_resident_kib().expect("VmHWM");

    println!("phase\tng-psp-open-cost");
    println!("store\t{}", store.display());
    println!("open-samples\t{samples}");
    println!("rounds\t{rounds}");
    println!("records-read\t{records}");
    println!("checksum\t{checksum}");
    println!("peak-rss-kib-before-opening\t{before_opening}");
    println!("peak-rss-kib-after-opening\t{after_opening}");
    println!("peak-rss-kib\t{after_walking}");
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **`VmHWM` and `VmRSS` differ by three characters and by everything that matters**: one is
    /// the peak for the life of the process and the other is the size right now. Reading the
    /// second would make every sample count report whatever the process happened to be holding at
    /// the moment it looked, and the slope would be noise.
    #[test]
    fn the_peak_is_taken_from_the_high_water_mark_and_not_the_current_size() {
        let status = "\
Name:\tng_psp_open_cost
VmPeak:\t 2410696 kB
VmSize:\t 2400000 kB
VmHWM:\t  487216 kB
VmRSS:\t   12345 kB
";
        assert_eq!(kib_from_status(status), Some(487_216));
    }

    #[test]
    fn a_status_without_the_field_is_none_rather_than_a_zero() {
        assert_eq!(kib_from_status("Name:\tsomething\nVmRSS:\t 12 kB\n"), None);
        assert_eq!(kib_from_status(""), None);
    }
}
