//! **A window's memory three ways**, for the refinement memory probes (`refine_replay`,
//! `refine_mem_probe`): the bytes the program asked for and still held at its peak, the resident
//! memory the kernel saw at its peak, and — printed by the probe beside them — what the window
//! budget counted.
//!
//! The allocator is `pvc`'s, mimalloc, wrapped to count the bytes asked of it; built with
//! `--no-default-features`, the system allocator, wrapped the same way (glibc's malloc on Linux),
//! so the two can be compared. Under `dhat-heap` dhat takes the allocator slot and the requested
//! bytes are not counted here.
//!
//! Environment: `MEMORY_SAMPLE_MS`, how often the resident memory is read (20 ms by default);
//! `MEMORY_TRACE`, a file each reading goes to as `seconds resident anonymous requested`, in bytes.
#![allow(dead_code, reason = "each probe uses its own part")]

use std::time::Duration;

#[cfg(not(feature = "dhat-heap"))]
use counting_allocator::CountingAllocator;
use counting_allocator::{ResidentSampler, kernel_peak, reset_kernel_peak};

#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static ALLOC: CountingAllocator<mimalloc::MiMalloc> = CountingAllocator::new(mimalloc::MiMalloc);

#[cfg(not(any(feature = "alloc-mimalloc", feature = "dhat-heap")))]
#[global_allocator]
static ALLOC: CountingAllocator<std::alloc::System> = CountingAllocator::new(std::alloc::System);

/// The allocator this build runs on.
pub const ALLOCATOR: &str = if cfg!(feature = "dhat-heap") {
    "dhat"
} else if cfg!(feature = "alloc-mimalloc") {
    "mimalloc"
} else {
    "system"
};

/// The bytes held now and the peak since the last reset, as the counting allocator saw them.
fn requested() -> Option<(u64, u64)> {
    #[cfg(not(feature = "dhat-heap"))]
    return Some((ALLOC.live_bytes() as u64, ALLOC.peak_bytes() as u64));
    #[cfg(feature = "dhat-heap")]
    None
}

fn requested_now() -> u64 {
    requested().map_or(0, |(now, _)| now)
}

fn reset_requested() {
    #[cfg(not(feature = "dhat-heap"))]
    ALLOC.reset_peak();
}

/// **One window's memory**, in bytes, from the end of the window before to the end of this one.
pub struct WindowPeaks {
    /// The most bytes the program held at once, and what it holds at the window's end.
    pub requested_peak: Option<u64>,
    pub requested_end: Option<u64>,
    /// The highest resident memory sampled, and its anonymous part (Linux).
    pub resident_peak: u64,
    pub anonymous_peak: Option<u64>,
    /// The kernel's own resident peak (Linux), which no sampling interval can miss.
    pub kernel_peak: Option<u64>,
}

impl WindowPeaks {
    /// The peaks as `key=value` fields in `unit` bytes (1e9 for GB, 1e6 for MB).
    pub fn describe(&self, unit: f64) -> String {
        let show = |value: Option<u64>| {
            value.map_or_else(
                || "-".to_string(),
                |bytes| format!("{:.2}", bytes as f64 / unit),
            )
        };
        format!(
            "requested_peak={} requested_end={} resident_peak={} anonymous_peak={} kernel_peak={}",
            show(self.requested_peak),
            show(self.requested_end),
            show(Some(self.resident_peak)),
            show(self.anonymous_peak),
            show(self.kernel_peak),
        )
    }
}

/// **Follows the memory window by window**: [`Self::take`] returns the peaks since the last call
/// and starts them again.
pub struct WindowMemory {
    sampler: ResidentSampler,
}

impl WindowMemory {
    pub fn start() -> Self {
        let every = std::env::var("MEMORY_SAMPLE_MS")
            .ok()
            .and_then(|text| text.parse().ok())
            .unwrap_or(20);
        let trace = std::env::var("MEMORY_TRACE").ok().map(|path| {
            Box::new(std::io::BufWriter::new(
                std::fs::File::create(path).expect("MEMORY_TRACE is a file that can be made"),
            )) as Box<dyn std::io::Write + Send>
        });
        reset_kernel_peak();
        reset_requested();
        Self {
            sampler: ResidentSampler::start(Duration::from_millis(every), trace, requested_now),
        }
    }

    pub fn take(&self) -> WindowPeaks {
        let requested = requested();
        let resident = self.sampler.take_peak();
        let peaks = WindowPeaks {
            requested_peak: requested.map(|(_, peak)| peak),
            requested_end: requested.map(|(now, _)| now),
            resident_peak: resident.total,
            anonymous_peak: resident.anonymous,
            kernel_peak: kernel_peak(),
        };
        reset_kernel_peak();
        reset_requested();
        peaks
    }
}
