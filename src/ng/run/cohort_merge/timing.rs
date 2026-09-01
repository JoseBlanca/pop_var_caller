//! Where a merge's wall time goes, counted from inside it — **compiled only under
//! `--features merge-timing`**.
//!
//! **Timing that ships in the hot path changes what it measures**, so every counter here is a
//! no-op type without the feature: [`Stopwatch::start`] reads no clock, [`Counter::add`]
//! touches no memory, and the calls in the merge compile away. The merge's own source reads
//! the same in both builds, which is the point — a driver with the timing spliced in by hand
//! is a different driver.
//!
//! **What it is for.** The parallel merge gives back much less than its thread count suggests
//! and there are three candidate reasons: the round barrier, the fixed cost the merge pays per
//! building region per sample, and allocation. These counters tell them apart. Nothing here
//! decides anything about the merge's answer; the report is printed by whichever probe asked
//! for it.
//!
//! **The counters are global and the report is not reentrant.** One merge at a time, and
//! [`reset`] before it. A run that merges on several threads at once would sum them, which no
//! probe does.

use std::fmt;

/// A running total of nanoseconds, or of events — an atomic counter under the feature and
/// nothing at all without it.
#[cfg(feature = "merge-timing")]
#[derive(Debug)]
pub struct Counter(std::sync::atomic::AtomicU64);

/// A running total of nanoseconds, or of events — an atomic counter under the feature and
/// nothing at all without it.
#[cfg(not(feature = "merge-timing"))]
#[derive(Debug)]
pub struct Counter;

#[cfg(feature = "merge-timing")]
impl Counter {
    /// A counter at zero.
    pub const fn new() -> Self {
        Self(std::sync::atomic::AtomicU64::new(0))
    }

    /// Add `amount` to the total.
    #[inline]
    pub fn add(&self, amount: u64) {
        self.0
            .fetch_add(amount, std::sync::atomic::Ordering::Relaxed);
    }

    /// Raise the total to `amount` if `amount` is larger — the round's slowest builder.
    #[inline]
    pub fn raise_to(&self, amount: u64) {
        self.0
            .fetch_max(amount, std::sync::atomic::Ordering::Relaxed);
    }

    /// The total so far.
    pub fn get(&self) -> u64 {
        self.0.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Read the total and put it back to zero.
    pub fn take(&self) -> u64 {
        self.0.swap(0, std::sync::atomic::Ordering::Relaxed)
    }
}

#[cfg(not(feature = "merge-timing"))]
impl Counter {
    /// A counter at zero.
    pub const fn new() -> Self {
        Self
    }

    /// Add `amount` to the total — nothing, in a build without the feature.
    #[inline(always)]
    pub fn add(&self, _amount: u64) {}

    /// Raise the total to `amount` — nothing, in a build without the feature.
    #[inline(always)]
    pub fn raise_to(&self, _amount: u64) {}

    /// The total so far — always zero without the feature.
    pub fn get(&self) -> u64 {
        0
    }

    /// Read the total and put it back to zero — always zero without the feature.
    pub fn take(&self) -> u64 {
        0
    }
}

impl Default for Counter {
    fn default() -> Self {
        Self::new()
    }
}

/// A clock started at some point in the merge, read when the thing being timed ends.
///
/// **It reads no clock without the feature**, so a build that is not measuring pays neither the
/// `clock_gettime` nor the register pressure of holding an instant across the timed span.
#[cfg(feature = "merge-timing")]
#[derive(Debug, Clone, Copy)]
pub struct Stopwatch(std::time::Instant);

/// A clock started at some point in the merge, read when the thing being timed ends.
#[cfg(not(feature = "merge-timing"))]
#[derive(Debug, Clone, Copy)]
pub struct Stopwatch;

#[cfg(feature = "merge-timing")]
impl Stopwatch {
    /// Start timing now.
    #[inline]
    pub fn start() -> Self {
        Self(std::time::Instant::now())
    }

    /// Nanoseconds since [`start`](Self::start).
    #[inline]
    pub fn elapsed_nanos(self) -> u64 {
        u64::try_from(self.0.elapsed().as_nanos()).unwrap_or(u64::MAX)
    }

    /// Add the time since [`start`](Self::start) to `counter`.
    #[inline]
    pub fn add_to(self, counter: &Counter) {
        counter.add(self.elapsed_nanos());
    }
}

#[cfg(not(feature = "merge-timing"))]
impl Stopwatch {
    /// Start timing now — nothing, in a build without the feature.
    #[inline(always)]
    pub fn start() -> Self {
        Self
    }

    /// Nanoseconds since `start` — always zero without the feature.
    #[inline(always)]
    pub fn elapsed_nanos(self) -> u64 {
        0
    }

    /// Add the time since `start` to `counter` — nothing, in a build without the feature.
    #[inline(always)]
    pub fn add_to(self, _counter: &Counter) {}
}

/// How many rounds of building regions the merge ran.
pub static ROUNDS: Counter = Counter::new();
/// How many building regions were handed to a builder, over the whole merge.
pub static REGIONS: Counter = Counter::new();
/// How many of those returned without opening the walk, because no locus begins in them.
pub static REGIONS_WITH_NO_LOCUS: Counter = Counter::new();
/// Nanoseconds the organiser spent deciding what to drop and moving it out — `evict_before`.
pub static EVICT_NANOS: Counter = Counter::new();
/// Nanoseconds the organiser spent drawing every sample's reader forward — `cover`.
pub static COVER_NANOS: Counter = Counter::new();
/// How many sweeps over the cohort those covers took, summed.
pub static COVER_SWEEPS: Counter = Counter::new();
/// Nanoseconds the samples' own drawing was busy inside a cover, summed across threads.
///
/// **Against [`COVER_NANOS`] this is what a cover's threads bought.** The cover is spread over
/// the cohort, so its wall time is this divided by the threads that were there — when the
/// spreading works.
pub static COVER_BUSY_NANOS: Counter = Counter::new();
/// Nanoseconds between the round's builders being launched and the last of them finishing.
pub static ROUND_WALL_NANOS: Counter = Counter::new();
/// Nanoseconds the builders themselves were busy, summed across threads.
pub static BUILDER_BUSY_NANOS: Counter = Counter::new();
/// The slowest single builder of each round, summed over rounds — the floor a barrier imposes.
pub static SLOWEST_BUILDER_NANOS: Counter = Counter::new();
/// The slowest builder of the round now running. Reset by the organiser at each round's start.
pub static SLOWEST_IN_THIS_ROUND_NANOS: Counter = Counter::new();
/// Nanoseconds spent building the per-sample window a builder reads — inside builder time.
pub static WINDOW_NANOS: Counter = Counter::new();
/// Nanoseconds spent setting a region's walk up, before any locus is closed — inside builder
/// time.
pub static WALK_SETUP_NANOS: Counter = Counter::new();
/// Nanoseconds the organiser spent releasing loci in region order — `submit` and `drain_ready`.
///
/// **The parallel driver's alone**, since 2026-09-01. The serial cached driver used to charge
/// its per-region `Vec::extend` here; it no longer has one, because each locus now goes
/// straight to a sink where it is built, so a serial breakdown reads zero for this and the
/// time it used to hold is inside `BUILDER_BUSY_NANOS`.
pub static ORGANISE_NANOS: Counter = Counter::new();
/// Nanoseconds spent on each built locus **after it was assembled and before the next one is
/// closed** — inside builder time, and subtracted from it to leave the assembling alone.
///
/// **In a calling run this is the genotyping**, which is the one term
/// `run_streaming.md` §11 question 7 has no measurement for and which decides the pool
/// milestone's shape. In the merge's own collecting driver it is the push onto the outcome's
/// vector, which is nearly nothing — the counter measures whatever the sink is, not what it
/// is for.
///
/// **Read once per built locus**, which costs two `clock_gettime` calls a locus under
/// `--features merge-timing` and nothing at all without it. Only
/// `merge_cohort_handing_each_locus_over` adds to it; the parallel driver goes through the
/// collecting `build_region`, whose sink is untimed, so a parallel breakdown reads zero here.
pub static AFTER_ASSEMBLY_NANOS: Counter = Counter::new();
/// Nanoseconds of the whole merge, from the driver's first line to its last.
pub static MERGE_WALL_NANOS: Counter = Counter::new();

/// Every counter, back to zero — call before the merge that is to be measured.
pub fn reset() {
    for counter in [
        &ROUNDS,
        &REGIONS,
        &REGIONS_WITH_NO_LOCUS,
        &EVICT_NANOS,
        &COVER_NANOS,
        &COVER_SWEEPS,
        &COVER_BUSY_NANOS,
        &ROUND_WALL_NANOS,
        &BUILDER_BUSY_NANOS,
        &SLOWEST_BUILDER_NANOS,
        &SLOWEST_IN_THIS_ROUND_NANOS,
        &WINDOW_NANOS,
        &WALK_SETUP_NANOS,
        &ORGANISE_NANOS,
        &AFTER_ASSEMBLY_NANOS,
        &MERGE_WALL_NANOS,
    ] {
        counter.take();
    }
}

/// What one merge's counters say, in milliseconds, with the derived shares worked out.
///
/// **Everything here is the merge that just ran**, so a caller that wants a median of several
/// runs collects several of these rather than summing counters across them.
#[derive(Debug, Clone, Copy)]
pub struct Report {
    /// How many threads rayon had while the merge ran — the divisor the ideal build time uses.
    pub threads: usize,
    /// Rounds of building regions.
    pub rounds: u64,
    /// Building regions handed to a builder.
    pub regions: u64,
    /// Of those, how many held no locus and so opened no walk.
    pub regions_with_no_locus: u64,
    /// Sweeps over the cohort the covers took, summed.
    pub cover_sweeps: u64,
    /// The whole merge.
    pub merge_wall_ms: f64,
    /// Deciding what to evict and moving it out, on the organiser's thread.
    pub evict_ms: f64,
    /// Drawing the readers forward, on the organiser's thread.
    pub cover_ms: f64,
    /// The samples' own drawing inside those covers, summed across threads.
    pub cover_busy_ms: f64,
    /// Launching a round's builders and waiting for the last of them.
    pub round_wall_ms: f64,
    /// The builders' own work, summed across threads.
    pub builder_busy_ms: f64,
    /// The slowest builder of each round, summed over rounds.
    pub slowest_builder_ms: f64,
    /// Building the per-sample windows, inside the builders' work.
    pub window_ms: f64,
    /// Setting the walks up, inside the builders' work.
    pub walk_setup_ms: f64,
    /// Releasing loci in region order, on the organiser's thread.
    pub organise_ms: f64,
    /// What was done with each locus once it was assembled — the genotyping, in a calling
    /// run. **Inside [`Self::builder_busy_ms`]**, so the assembling alone is the difference
    /// ([`Self::assembling_loci_ms`]).
    pub after_assembly_ms: f64,
}

/// Read every counter into a [`Report`]. `threads` is what rayon had while the merge ran.
pub fn report(threads: usize) -> Report {
    let ms = |counter: &Counter| counter.get() as f64 / 1e6;
    Report {
        threads,
        rounds: ROUNDS.get(),
        regions: REGIONS.get(),
        regions_with_no_locus: REGIONS_WITH_NO_LOCUS.get(),
        cover_sweeps: COVER_SWEEPS.get(),
        merge_wall_ms: ms(&MERGE_WALL_NANOS),
        evict_ms: ms(&EVICT_NANOS),
        cover_ms: ms(&COVER_NANOS),
        cover_busy_ms: ms(&COVER_BUSY_NANOS),
        round_wall_ms: ms(&ROUND_WALL_NANOS),
        builder_busy_ms: ms(&BUILDER_BUSY_NANOS),
        slowest_builder_ms: ms(&SLOWEST_BUILDER_NANOS),
        window_ms: ms(&WINDOW_NANOS),
        walk_setup_ms: ms(&WALK_SETUP_NANOS),
        organise_ms: ms(&ORGANISE_NANOS),
        after_assembly_ms: ms(&AFTER_ASSEMBLY_NANOS),
    }
}

impl Report {
    /// **Assembling the loci, with what was done to each of them afterwards taken out** —
    /// the builders' own work minus [`Self::after_assembly_ms`].
    ///
    /// This and `after_assembly_ms` are the two halves of a run's builder time, and telling
    /// them apart is what says whether a calling run is spending its time putting loci
    /// together or genotyping them — the split spec §11 question 7 asks for.
    ///
    /// **The subtraction cannot go negative**: every nanosecond counted after a locus is
    /// assembled is counted inside the builder's own stopwatch, which starts before the
    /// region's walk opens and stops after its last locus is handed over.
    pub fn assembling_loci_ms(&self) -> f64 {
        self.builder_busy_ms - self.after_assembly_ms
    }

    /// The build phase's time if every round's builders had been spread perfectly over the
    /// threads — the builders' summed work divided by the threads that were there.
    pub fn perfectly_spread_ms(&self) -> f64 {
        self.builder_busy_ms / self.threads as f64
    }

    /// What the round's slowest builder costs over a perfect spread: the wait a barrier
    /// imposes even when nothing is wasted launching the builders.
    ///
    /// **This is the round barrier's own price**, and it is zero when every builder in a round
    /// takes the same time.
    pub fn stragglers_ms(&self) -> f64 {
        (self.slowest_builder_ms - self.perfectly_spread_ms()).max(0.0)
    }

    /// What the build phase cost beyond its slowest builder and beyond a perfect spread —
    /// launching the builders and collecting them, which is rayon's.
    pub fn launching_builders_ms(&self) -> f64 {
        (self.round_wall_ms - self.perfectly_spread_ms().max(self.slowest_builder_ms)).max(0.0)
    }

    /// The merge's time that no counter above accounts for.
    pub fn unaccounted_ms(&self) -> f64 {
        self.merge_wall_ms - self.evict_ms - self.cover_ms - self.round_wall_ms - self.organise_ms
    }
}

impl fmt::Display for Report {
    /// One merge's breakdown, as lines a probe can print unchanged.
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let share = |part: f64| 100.0 * part / self.merge_wall_ms.max(f64::MIN_POSITIVE);
        writeln!(
            out,
            "# rounds: {}, building regions: {} ({} held no locus), cover sweeps: {}",
            self.rounds, self.regions, self.regions_with_no_locus, self.cover_sweeps,
        )?;
        writeln!(out, "# threads while merging: {}", self.threads)?;
        writeln!(out, "part, ms, share_of_merge_%")?;
        writeln!(out, "whole merge, {:.1}, 100.0", self.merge_wall_ms)?;
        writeln!(
            out,
            "  drawing the readers forward (cover), {:.1}, {:.1}",
            self.cover_ms,
            share(self.cover_ms)
        )?;
        writeln!(
            out,
            "    of which the samples' own drawing (summed over threads), {:.1}, -",
            self.cover_busy_ms
        )?;
        writeln!(
            out,
            "    the same drawing spread perfectly over {} threads, {:.1}, {:.1}",
            self.threads,
            self.cover_busy_ms / self.threads as f64,
            share(self.cover_busy_ms / self.threads as f64)
        )?;
        writeln!(
            out,
            "  evicting what the round has passed, {:.1}, {:.1}",
            self.evict_ms,
            share(self.evict_ms)
        )?;
        writeln!(
            out,
            "  the round's builders, start to last one finished, {:.1}, {:.1}",
            self.round_wall_ms,
            share(self.round_wall_ms)
        )?;
        writeln!(
            out,
            "    of which the builders were busy (summed over threads), {:.1}, -",
            self.builder_busy_ms
        )?;
        writeln!(
            out,
            "    of which building the per-sample windows, {:.1}, -",
            self.window_ms
        )?;
        writeln!(
            out,
            "    of which setting each region's walk up, {:.1}, -",
            self.walk_setup_ms
        )?;
        writeln!(
            out,
            "    of which assembling the loci, {:.1}, -",
            self.assembling_loci_ms()
        )?;
        writeln!(
            out,
            "    of which what each locus was then handed to (a run: genotyping it), {:.1}, -",
            self.after_assembly_ms
        )?;
        writeln!(
            out,
            "    the same work spread perfectly over {} threads, {:.1}, {:.1}",
            self.threads,
            self.perfectly_spread_ms(),
            share(self.perfectly_spread_ms())
        )?;
        writeln!(
            out,
            "    waiting for each round's slowest builder, {:.1}, {:.1}",
            self.stragglers_ms(),
            share(self.stragglers_ms())
        )?;
        writeln!(
            out,
            "    launching and collecting the builders, {:.1}, {:.1}",
            self.launching_builders_ms(),
            share(self.launching_builders_ms())
        )?;
        writeln!(
            out,
            "  releasing loci in region order, {:.1}, {:.1}",
            self.organise_ms,
            share(self.organise_ms)
        )?;
        write!(
            out,
            "  everything else, {:.1}, {:.1}",
            self.unaccounted_ms(),
            share(self.unaccounted_ms())
        )
    }
}
