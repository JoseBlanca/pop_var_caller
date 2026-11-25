//! **What a process asks of its allocator, against what it holds resident** — for the memory
//! probes in `examples/`, which compare both with what refinement's window budget counts.
//!
//! - [`CountingAllocator`] wraps a global allocator and keeps two numbers: the bytes the program
//!   holds now (asked for and not yet freed) and the most it held since the last
//!   [`CountingAllocator::reset_peak`]. These are the program's live data, whatever the allocator
//!   underneath keeps cached or fragmented.
//! - [`resident`] reads what the kernel says the process holds in memory now; on Linux it splits
//!   off the anonymous part (heap and stacks, not mapped files).
//! - [`ResidentSampler`] reads [`resident`] on a thread of its own every few milliseconds and keeps
//!   the highest value since its last reset, so a probe can report a resident peak per window on
//!   macOS too, where the kernel's own peak cannot be reset.
//!
//! The difference between the resident peak and the requested peak is what the allocator holds
//! beyond the program's live data: its caches, its per-thread heaps, and freed memory it has not
//! returned to the system.
use std::alloc::{GlobalAlloc, Layout};
use std::io::Write as _;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

/// **A global allocator that counts the bytes asked of it**, around `inner`.
///
/// Counts the sizes asked for, not what `inner` rounds them up to: a 33-byte request counts 33.
/// Every allocation and free adds one atomic update to a shared counter, and one more when it
/// raises the peak; on many threads that costs some speed, never correctness.
pub struct CountingAllocator<A> {
    inner: A,
    live: AtomicUsize,
    peak: AtomicUsize,
}

impl<A> CountingAllocator<A> {
    /// Wraps `inner`, counting from zero.
    pub const fn new(inner: A) -> Self {
        Self {
            inner,
            live: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
        }
    }

    /// The bytes held now: asked for and not yet freed.
    pub fn live_bytes(&self) -> usize {
        self.live.load(Relaxed)
    }

    /// The most bytes held at once since the start or the last [`Self::reset_peak`].
    pub fn peak_bytes(&self) -> usize {
        self.peak.load(Relaxed)
    }

    /// Starts the peak again from what is held now.
    pub fn reset_peak(&self) {
        self.peak.store(self.live.load(Relaxed), Relaxed);
    }

    fn add(&self, bytes: usize) {
        let now = self.live.fetch_add(bytes, Relaxed) + bytes;
        if now > self.peak.load(Relaxed) {
            self.peak.fetch_max(now, Relaxed);
        }
    }

    fn remove(&self, bytes: usize) {
        self.live.fetch_sub(bytes, Relaxed);
    }
}

// SAFETY: every call is passed to `inner` unchanged; only the counters are added.
unsafe impl<A: GlobalAlloc> GlobalAlloc for CountingAllocator<A> {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        // SAFETY: the caller's contract is `inner`'s.
        let pointer = unsafe { self.inner.alloc(layout) };
        if !pointer.is_null() {
            self.add(layout.size());
        }
        pointer
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        // SAFETY: as `alloc`.
        let pointer = unsafe { self.inner.alloc_zeroed(layout) };
        if !pointer.is_null() {
            self.add(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        // SAFETY: as `alloc`.
        unsafe { self.inner.dealloc(pointer, layout) };
        self.remove(layout.size());
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // SAFETY: as `alloc`.
        let moved = unsafe { self.inner.realloc(pointer, layout, new_size) };
        if !moved.is_null() {
            if new_size >= layout.size() {
                self.add(new_size - layout.size());
            } else {
                self.remove(layout.size() - new_size);
            }
        }
        moved
    }
}

/// **What the kernel says the process holds in memory now**, in bytes.
#[derive(Clone, Copy, Debug, Default)]
pub struct Resident {
    /// Every resident page.
    pub total: u64,
    /// Of those, the anonymous ones — heap and stacks, not mapped files. Linux only.
    pub anonymous: Option<u64>,
}

/// The process's resident memory now, or `None` where it cannot be read.
pub fn resident() -> Option<Resident> {
    resident_now()
}

#[cfg(target_os = "linux")]
fn resident_now() -> Option<Resident> {
    // `statm`: size, resident, shared (file-backed resident), … in pages.
    let statm = std::fs::read_to_string("/proc/self/statm").ok()?;
    let mut fields = statm.split_whitespace().skip(1);
    let resident: u64 = fields.next()?.parse().ok()?;
    let shared: u64 = fields.next()?.parse().ok()?;
    // SAFETY: `sysconf` reads a constant.
    let page = unsafe { libc::sysconf(libc::_SC_PAGESIZE) } as u64;
    Some(Resident {
        total: resident * page,
        anonymous: Some(resident.saturating_sub(shared) * page),
    })
}

#[cfg(target_os = "macos")]
fn resident_now() -> Option<Resident> {
    let mut info: libc::proc_taskinfo = unsafe { std::mem::zeroed() };
    let size = std::mem::size_of::<libc::proc_taskinfo>() as libc::c_int;
    // SAFETY: `info` is a `proc_taskinfo` of the size passed.
    let written = unsafe {
        libc::proc_pidinfo(
            libc::getpid(),
            libc::PROC_PIDTASKINFO,
            0,
            (&raw mut info).cast(),
            size,
        )
    };
    (written == size).then_some(Resident {
        total: info.pti_resident_size,
        anonymous: None,
    })
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn resident_now() -> Option<Resident> {
    None
}

/// On Linux, the kernel's own resident peak (`VmHWM`) since the start or the last
/// [`reset_kernel_peak`], in bytes.
pub fn kernel_peak() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let kb: u64 = status
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))?
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()
        .ok()?;
    Some(kb * 1024)
}

/// On Linux, starts the kernel's resident peak again from what the process holds now.
pub fn reset_kernel_peak() {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
}

/// The highest resident memory seen since the last reset, as [`Resident`]'s fields.
#[derive(Default)]
struct Highest {
    total: AtomicU64,
    anonymous: AtomicU64,
}

/// **Reads the resident memory every `every`** on a thread of its own, keeping the highest value
/// since the last [`Self::take_peak`]; stops when dropped. Called with a `trace`, it also writes
/// each reading as a line — seconds since the start, resident and anonymous bytes, and whatever
/// `extra` returns (a probe's requested bytes, say) — so a window's memory can be followed in time.
pub struct ResidentSampler {
    highest: Arc<Highest>,
    stop: Arc<AtomicBool>,
    thread: Option<JoinHandle<()>>,
}

impl ResidentSampler {
    /// Starts sampling.
    pub fn start(
        every: Duration,
        mut trace: Option<Box<dyn std::io::Write + Send>>,
        extra: fn() -> u64,
    ) -> Self {
        let highest = Arc::new(Highest::default());
        let stop = Arc::new(AtomicBool::new(false));
        let thread = {
            let (highest, stop) = (Arc::clone(&highest), Arc::clone(&stop));
            std::thread::spawn(move || {
                let started = Instant::now();
                while !stop.load(Relaxed) {
                    if let Some(now) = resident() {
                        highest.total.fetch_max(now.total, Relaxed);
                        let anonymous = now.anonymous.unwrap_or(0);
                        highest.anonymous.fetch_max(anonymous, Relaxed);
                        if let Some(trace) = trace.as_mut() {
                            let _ = writeln!(
                                trace,
                                "{:.3}\t{}\t{}\t{}",
                                started.elapsed().as_secs_f64(),
                                now.total,
                                anonymous,
                                extra()
                            );
                        }
                    }
                    std::thread::sleep(every);
                }
                if let Some(trace) = trace.as_mut() {
                    let _ = trace.flush();
                }
            })
        };
        Self {
            highest,
            stop,
            thread: Some(thread),
        }
    }

    /// The highest resident memory since the last call, then starts again from now.
    pub fn take_peak(&self) -> Resident {
        let now = resident().unwrap_or_default();
        let total = self.highest.total.swap(now.total, Relaxed);
        let anonymous = self
            .highest
            .anonymous
            .swap(now.anonymous.unwrap_or(0), Relaxed);
        Resident {
            total: total.max(now.total),
            anonymous: now.anonymous.map(|now| anonymous.max(now)),
        }
    }
}

impl Drop for ResidentSampler {
    fn drop(&mut self) {
        self.stop.store(true, Relaxed);
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}
