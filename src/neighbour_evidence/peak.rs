//! **How much memory refinement is expected to hold at its peak, and whether this run may hold
//! it** (`doc/devel/ng/spec/neighbour_evidence_pass.md` §6, memory).
//!
//! The window budget (`--refine-window-bytes`) counts what a window holds for every sample at every
//! record. The process holds more than that count, in two parts the budget leaves out on purpose
//! (`doc/devel/reports/implementations/refine_memory_overhead_2026-11-25.md`):
//!
//! - **the walks' working arrays, one set a thread**: about [`WALK_BYTES_A_RECORD_A_THREAD`] for
//!   every record the window walks. Counting them in the budget would make the windows, and so the
//!   output, depend on the thread count (§2.3);
//! - **the allocator's slack**: memory freed and not yet given back to the system, about
//!   [`ALLOCATOR_SLACK`] of what the program holds.
//!
//! So the expected peak is
//!
//! ```text
//! peak = (counted + 1.05 kB × walking threads × walked records) × (1 + 0.16)
//! ```
//!
//! where *walking threads* is the run's threads or its samples, whichever is fewer, since a thread
//! walks one sample at a time. It is estimated twice:
//!
//! - **before calling starts** ([`PeakEstimate::before_any_window`]), from the budget, the samples,
//!   the threads and the longest contig, so that a run that cannot fit stops before hours of
//!   calling. *counted* is the budget, or what a window of the longest contig could hold at
//!   [`RECORDS_A_KB_CEILING`] if that is less; *walked records* is what a window counted at the
//!   budget walks at [`COUNTED_BYTES_A_SAMPLE_WALKED_RECORD`], or what such a window could walk,
//!   if less. Both are upper bounds where the cohort is small and the window is set by its length;
//! - **at every window** ([`PeakEstimate::for_window`]), from what the budget counted for it and the
//!   records it walks, both known exactly once it is planned.
//!
//! **A run stops when the estimate is above [`STOP_SHARE`] of the memory it may use** (the smaller
//! of the machine's memory and its control group's limit, [`crate::memory_limit`]), unless
//! `--refine-ignore-memory-limit` was given. The constants were measured on Linux on a drawn cohort
//! of 300 samples at up to 6 threads and on 2,176 tomato samples at 48 threads; see each one.

use crate::memory_limit::MemoryLimit;
use crate::neighbour_evidence::windows::{NEAR_MARGIN_BASES, WindowSizes};

/// **The walks' working arrays, bytes a walked record a thread**: the templates' allele shares, the
/// sample's likelihoods, posteriors and drawn paths. *Measured* 2026-11-25 at 1.04 to 1.07 kB on the
/// memory probe's drawn cohort, at tomato25's mix of 2.8 alleles a record on average, 100- to
/// 400-kb cores, 1 to 6 threads (`refine_memory_overhead_2026-11-25.md` §3). It grows with the
/// alleles a record carries, about 136 bytes an allele (the shares of 16 templates and the
/// background, `f64`).
pub const WALK_BYTES_A_RECORD_A_THREAD: u64 = 1_050;

/// **What the budget counts, bytes a sample a walked record**, at a window filled to the budget
/// with thinned margins: 2,314 MB for 45,600 walked records of 300 samples on the drawn cohort
/// (169 bytes); at most 216 on tomato25's highest window (80 GB or less for about 170,000 walked
/// records of 2,176 samples). The lower figure gives more walked records for the same budget, so
/// the estimate before calling errs high.
pub const COUNTED_BYTES_A_SAMPLE_WALKED_RECORD: u64 = 170;

/// **What the budget counts, bytes a sample a record held**, for the estimate before calling of a
/// window whose length, not the budget, sets it: about 240 *measured* on 400 tomato samples (pass
/// spec §6), rounded up.
pub const HELD_BYTES_A_SAMPLE_RECORD: u64 = 250;

/// **What a record costs whatever the cohort**: its position, its line's own columns, its place.
/// *Not measured*: an allowance, so that the estimate before calling of a two-sample window is
/// not near zero. At 600,000 records it is 0.6 GB.
pub const HELD_BYTES_A_RECORD: u64 = 1_000;

/// **The densest records a kb the estimate before calling assumes**, where the window's length
/// sets its size: the drawn cohort's 150. tomato25's densest stretch holds about 90.
pub const RECORDS_A_KB_CEILING: u64 = 150;

/// **The allocator's slack, a share of what the program holds**: 16%, what tomato25's highest
/// window showed (103.0 GB measured against about 89 GB held, glibc, 48 threads), in the middle of
/// what mimalloc, pvc's allocator, held on the drawn cohort at its default settings (12 to 25%,
/// Linux, 1 to 6 threads). With mimalloc returning freed memory at once (`purge_delay` 0) it was
/// 6 to 7%.
pub const ALLOCATOR_SLACK: f64 = 0.16;

/// **The share of the memory limit above which the run stops**: 90%. The estimate is the middle of
/// what was measured; the peaks measured lay up to 8% above it (the slack's spread, 12 to 25%
/// around 16%) and the rest of what the program holds within 5% of the count either way, so a
/// run stopping at 90% keeps about 11% for those, and stops before the system would kill it.
pub const STOP_SHARE: f64 = 0.90;

/// **The expected peak of refinement**, and the parts it is made of.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeakEstimate {
    /// What the window budget counts for the window: its plan's count, or, before any window is
    /// planned, the most a window can count.
    pub counted_bytes: u64,
    /// The records the window walks, or before any window is planned, the most it may walk.
    pub walked_records: u64,
    /// The threads walking at once: the run's threads, or its samples if fewer.
    pub walking_threads: u64,
    /// `--refine-window-bytes`.
    pub budget_bytes: u64,
}

impl PeakEstimate {
    /// **Before any window is planned**: the most the run's windows may count and walk, from the
    /// window sizes, the samples, the threads and the longest contig.
    #[must_use]
    pub fn before_any_window(
        sizes: &WindowSizes,
        samples: usize,
        threads: usize,
        longest_contig: u64,
    ) -> Self {
        let samples = samples as u64;
        let span = (sizes.core_bases + 2 * sizes.margin_bases).min(longest_contig);
        // Every record within the core and the near margins is walked; beyond, at most the
        // thinning's records a kb, if it is denser than the records.
        let walked_by_length = match sizes.thinning {
            Some(thinning) => {
                let near = (sizes.core_bases + 2 * NEAR_MARGIN_BASES).min(span);
                let far_density = (u64::from(thinning.thousandths_per_kb.get()))
                    .min(RECORDS_A_KB_CEILING * 1_000);
                near * RECORDS_A_KB_CEILING / 1_000 + (span - near) * far_density / 1_000_000
            }
            None => span * RECORDS_A_KB_CEILING / 1_000,
        };
        let walked_by_budget =
            sizes.budget_bytes / (COUNTED_BYTES_A_SAMPLE_WALKED_RECORD * samples.max(1));
        let held_by_length = span * RECORDS_A_KB_CEILING / 1_000
            * (HELD_BYTES_A_RECORD + HELD_BYTES_A_SAMPLE_RECORD * samples);
        Self {
            counted_bytes: sizes.budget_bytes.min(held_by_length),
            walked_records: walked_by_budget.min(walked_by_length),
            walking_threads: (threads as u64).min(samples).max(1),
            budget_bytes: sizes.budget_bytes,
        }
    }

    /// **For a planned window**: what the budget counted for it, and the records its walks read.
    #[must_use]
    pub fn for_window(
        counted_bytes: u64,
        walked_records: u64,
        samples: usize,
        threads: usize,
        budget_bytes: u64,
    ) -> Self {
        Self {
            counted_bytes,
            walked_records,
            walking_threads: (threads as u64).min(samples as u64).max(1),
            budget_bytes,
        }
    }

    /// The walks' working arrays, every walking thread's.
    #[must_use]
    pub fn walk_bytes(&self) -> u64 {
        WALK_BYTES_A_RECORD_A_THREAD * self.walking_threads * self.walked_records
    }

    /// **The expected peak**: what the budget counts and the walks' arrays, with the allocator's
    /// slack on top.
    #[must_use]
    pub fn peak_bytes(&self) -> u64 {
        // A size to compare, to the byte: the rounding of a float is far below its precision.
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let peak =
            ((self.counted_bytes + self.walk_bytes()) as f64 * (1.0 + ALLOCATOR_SLACK)) as u64;
        peak
    }

    /// **The parts, for a person**: `70.0 GB for a window, 9.5 GB for 48 threads' working arrays
    /// (189,230 walked records), 16% for the allocator`.
    #[must_use]
    pub fn parts_in_words(&self) -> String {
        format!(
            "{} for a window as --refine-window-bytes counts it ({} allowed), {} for {} \
             thread{}' working arrays ({} walked records at {:.2} kB each a thread), and {:.0}% on \
             top for the memory allocator",
            gb(self.counted_bytes),
            gb(self.budget_bytes),
            gb(self.walk_bytes()),
            self.walking_threads,
            if self.walking_threads == 1 { "" } else { "s" },
            with_thousands(self.walked_records),
            WALK_BYTES_A_RECORD_A_THREAD as f64 / 1_000.0,
            ALLOCATOR_SLACK * 100.0
        )
    }
}

/// **What the memory check is given**: the limit, if it could be read, and whether
/// `--refine-ignore-memory-limit` asked the run to go on past it. The default checks nothing: a
/// decision made outside a command, such as a test's, carries no limit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MemoryGuard {
    /// The memory this process may use; `None` where it could not be read.
    pub limit: Option<MemoryLimit>,
    /// `--refine-ignore-memory-limit`.
    pub ignore_limit: bool,
}

/// **When the estimate was made**: before calling, or at a window about to be solved.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakMoment {
    /// Before any window is planned, from the budget, samples, threads and longest contig.
    BeforeCalling,
    /// At a planned window: its contig, counted from 0 in the header's order, and its core's first
    /// position.
    Window { contig: u32, position: u64 },
}

/// **Refinement is expected to need more of the memory limit than [`STOP_SHARE`]**: the estimate,
/// the limit and when it was made. Its message is what the run stops with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeakTooClose {
    pub estimate: PeakEstimate,
    pub limit: MemoryLimit,
    pub moment: PeakMoment,
    /// The largest `--refine-window-bytes` the estimate fits within [`STOP_SHARE`] of the limit at
    /// this run's threads; `None` where the threads' arrays alone do not fit.
    pub budget_that_fits: Option<u64>,
    /// The most threads that fit at this budget, where fewer than the run's but at least one.
    pub threads_that_fit: Option<u64>,
}

impl std::fmt::Display for PeakTooClose {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let peak = self.estimate.peak_bytes();
        let when = match self.moment {
            PeakMoment::BeforeCalling => String::new(),
            PeakMoment::Window { contig, position } => {
                format!(" in its window from contig {contig} position {position}")
            }
        };
        write!(
            f,
            "genotype refinement is expected to need about {}{when}, {:.0}% of the {} this run \
             may use ({}), and pvc stops above {:.0}% so that the system does not kill the run part \
             of the way through. The estimate is {}. ",
            gb(peak),
            share(peak, self.limit.bytes) * 100.0,
            gb(self.limit.bytes),
            self.limit.source,
            STOP_SHARE * 100.0,
            self.estimate.parts_in_words()
        )?;
        let mut remedies = Vec::new();
        if let Some(budget) = self.budget_that_fits {
            remedies.push(format!(
                "lower --refine-window-bytes to {budget} ({}) or less",
                gb(budget)
            ));
        }
        if let Some(threads) = self.threads_that_fit {
            remedies.push(format!("run on {threads} threads or fewer (--threads)"));
        }
        remedies.push(
            "give --refine-ignore-memory-limit to run anyway, if this run's peak is known to be \
             lower"
                .to_owned(),
        );
        let last = remedies.pop().unwrap_or_default();
        if remedies.is_empty() {
            write!(f, "To go on, {last}.")?;
        } else {
            write!(f, "To go on, {}, or {last}.", remedies.join(", "))?;
        }
        if matches!(self.moment, PeakMoment::Window { .. }) {
            write!(
                f,
                " With the hidden-duplication filter on, the run keeps its scratch directory \
                 and its calls, and `pvc resume-calling` finishes it with the new options."
            )?;
        }
        Ok(())
    }
}

/// **What the check found**, where the run goes on.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PeakCheck {
    /// The limit could not be read: nothing to check against.
    NoLimit,
    /// The estimate is within [`STOP_SHARE`] of the limit; `share` is its share of the limit.
    Fits { share: f64 },
    /// Above it, and `--refine-ignore-memory-limit` asked the run to go on.
    Ignored(PeakTooClose),
}

/// **Checks an estimate against the limit**: goes on where it fits, or where the limit could not
/// be read or the person asked to go on past it.
///
/// `budget_at` gives the estimate the run would make with another budget, at the same threads;
/// it finds the budget to suggest.
///
/// # Errors
///
/// If the estimate is above [`STOP_SHARE`] of the limit and the limit is not ignored.
pub fn check(
    estimate: PeakEstimate,
    guard: MemoryGuard,
    moment: PeakMoment,
    budget_at: impl Fn(u64) -> PeakEstimate,
) -> Result<PeakCheck, PeakTooClose> {
    let Some(limit) = guard.limit else {
        return Ok(PeakCheck::NoLimit);
    };
    let peak = estimate.peak_bytes();
    if !above_the_stop(peak, limit.bytes) {
        return Ok(PeakCheck::Fits {
            share: share(peak, limit.bytes),
        });
    }
    let too_close = PeakTooClose {
        estimate,
        limit,
        moment,
        budget_that_fits: largest_budget_that_fits(estimate.budget_bytes, limit.bytes, budget_at),
        threads_that_fit: threads_that_fit(estimate, limit.bytes),
    };
    if guard.ignore_limit {
        Ok(PeakCheck::Ignored(too_close))
    } else {
        Err(too_close)
    }
}

/// Whether `peak` is above [`STOP_SHARE`] of `limit`.
fn above_the_stop(peak: u64, limit: u64) -> bool {
    share(peak, limit) > STOP_SHARE
}

/// `part` as a share of `whole`.
#[allow(clippy::cast_precision_loss)]
fn share(part: u64, whole: u64) -> f64 {
    part as f64 / whole.max(1) as f64
}

/// **The largest budget whose estimate is within the stop**, rounded down to a whole GB (to a tenth
/// of one below 1 GB), found by halving: the estimate does not fall as the budget grows. `None`
/// where even a budget of nothing is over it.
fn largest_budget_that_fits(
    budget: u64,
    limit: u64,
    budget_at: impl Fn(u64) -> PeakEstimate,
) -> Option<u64> {
    let fits = |bytes: u64| !above_the_stop(budget_at(bytes).peak_bytes(), limit);
    if !fits(0) {
        return None;
    }
    let (mut low, mut high) = (0, budget);
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if fits(middle) {
            low = middle;
        } else {
            high = middle;
        }
    }
    let step = if low >= 1_000_000_000 {
        1_000_000_000
    } else {
        100_000_000
    };
    let rounded = low / step * step;
    (rounded > 0).then_some(rounded)
}

/// **The most threads that fit at this budget**, where that is fewer than the run walks with and
/// at least one.
fn threads_that_fit(estimate: PeakEstimate, limit: u64) -> Option<u64> {
    let fits = |threads: u64| {
        !above_the_stop(
            PeakEstimate {
                walking_threads: threads,
                ..estimate
            }
            .peak_bytes(),
            limit,
        )
    };
    (1..estimate.walking_threads)
        .rev()
        .find(|&threads| fits(threads))
}

/// A count with its thousands marked: `189,240`.
fn with_thousands(number: u64) -> String {
    let digits = number.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

/// Bytes as a person reads them here: `92.2 GB`.
#[allow(clippy::cast_precision_loss)]
fn gb(bytes: u64) -> String {
    format!("{:.1} GB", bytes as f64 / 1e9)
}

/// **The line a run logs for an estimate that fits**, or whose limit could not be read.
#[must_use]
pub fn fits_in_words(
    estimate: &PeakEstimate,
    found: &PeakCheck,
    limit: Option<MemoryLimit>,
) -> String {
    let peak = gb(estimate.peak_bytes());
    match (found, limit) {
        (PeakCheck::Fits { share }, Some(limit)) => format!(
            "about {peak} at its peak, {:.0}% of the {} this run may use ({}); {}",
            share * 100.0,
            gb(limit.bytes),
            limit.source,
            estimate.parts_in_words()
        ),
        (PeakCheck::Ignored(too_close), _) => format!(
            "about {peak} at its peak, above {:.0}% of the {} this run may use ({}); going on, \
             as --refine-ignore-memory-limit asks; {}",
            STOP_SHARE * 100.0,
            gb(too_close.limit.bytes),
            too_close.limit.source,
            estimate.parts_in_words()
        ),
        _ => format!(
            "about {peak} at its peak, not checked: the memory this run may use could not be \
             read; {}",
            estimate.parts_in_words()
        ),
    }
}

#[cfg(test)]
mod tests;
