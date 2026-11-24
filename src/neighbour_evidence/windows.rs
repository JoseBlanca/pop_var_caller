//! **Where each window's core and margins end** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §2.2, §6).
//!
//! The pass cannot write a record until it has seen the records after it, so it holds them in
//! windows: a *core*, written when the window is solved, and a *margin* either side read for
//! context and written by the neighbouring windows. Cores are runs of consecutive records in the
//! order they arrived, so every record is written exactly once and in order. A window never
//! crosses a contig, and spans the gaps between BED regions like any other distance.
//!
//! **Sizes come from length and from a memory budget in bytes.** A window holds, for every record
//! it spans, the record's line, its numbers and its reads, and each sample's state in the solve;
//! the caller says what each record costs ([`WindowBytes`]). It starts at a 2-Mb core with 1-Mb
//! margins; when that would hold more bytes than the budget allows, the core shrinks first, down to
//! the margins' length, then core and margins shrink together, down to the margin floor of
//! 500 kb, and then the core alone, down to the core floor (`--refine-core-floor`, the margin
//! floor unless set lower). A window still over budget at the floors is planned at them and says
//! by how much it is over, and the stage stops the run there: margins shorter than the floor would
//! cut off the evidence a record's walk reaches for (§4.2's `S` is at most 300 kb). A shorter core
//! cuts off nothing, since every record's margins stay as long; it only makes more windows, each
//! solving its margins again.
//!
//! **It counts what the caller holds outside the window too** ([`WindowBytes::outside`]): the
//! records before the window it keeps in case a later window plans longer margins, and those read
//! past it while larger windows were ruled out. They hold no solve state, but their lines, numbers
//! and reads are memory the run holds while the window is solved. Uncounted, the code review of
//! 2026-11-24 estimated them at about a sixth on top of the count of a window at the 500-kb floor
//! on 2,176 tomato samples; on the stage tests' drawn cohort, without reads, they add up to 4%.
//!
//! **It counts the rounds' haplotype store or the core's results, whichever is larger, not both**
//! ([`StoreAndResults`]): the solve lets the store go before it makes the results.
//!
//! The budget does not count the walks' own arrays, which every thread holds one of: counting
//! them would make the windows, and so the output, depend on the thread count (§2.3).
//!
//! **Thinned margins** (`--refine-margin-density`, §6, thinned margins): beyond [`NEAR_MARGIN_BASES`] of the core,
//! a window's walks read only the records the stage marked as walked in thinned margins, at most so
//! many a kb. The solve holds no state for the others, and the stage holds their numbers and reads
//! compressed, so the caller hands over what each record saves where it is skipped
//! ([`WindowBytes::skipped`]) and what it holds packed instead ([`WindowBytes::packed`]); a window
//! pays the first and not the second within the near margins, and the other way beyond them.

/// **How large a window may be** (§2.2, §6).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowSizes {
    /// The core's length, in bases, before the budget shrinks it.
    pub core_bases: u64,
    /// Each margin's length, in bases, before the budget shrinks it.
    pub margin_bases: u64,
    /// The shortest a margin may become (and, with it, the core).
    pub margin_floor_bases: u64,
    /// The shortest a core may become once its margins are at their floor
    /// (`--refine-core-floor`); at or above the margin floor, the core stops there too.
    pub core_floor_bases: u64,
    /// The most bytes a window should hold, as [`WindowBytes`] counts them
    /// (`--refine-window-bytes`).
    pub budget_bytes: u64,
    /// How thinly the margins are walked beyond [`NEAR_MARGIN_BASES`] of the core
    /// (`--refine-margin-density`); `None` walks every record.
    pub thinning: Option<MarginThinning>,
}

/// **How thinly a window's far margins are walked** (§6, thinned margins): the most records a kb a
/// walk reads there, on average, evenly spaced, in thousandths of a record. A stretch sparser than
/// that is walked in full.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MarginThinning {
    pub thousandths_per_kb: std::num::NonZeroU32,
}

impl MarginThinning {
    /// `--refine-margin-density`'s value, records a kb to three decimals; `None` at 0, which walks
    /// every record.
    ///
    /// # Errors
    ///
    /// If `text` is not a number from 0 to 1,000,000 with at most three decimals.
    pub fn parse(text: &str) -> Result<Option<Self>, String> {
        let refuse = || {
            format!(
                "`{text}` is not a number of records a kb, 0 to 1000000, at most three decimals"
            )
        };
        let (whole, fraction) = text.split_once('.').unwrap_or((text, ""));
        if fraction.len() > 3
            || whole.is_empty() && fraction.is_empty()
            || !whole
                .chars()
                .chain(fraction.chars())
                .all(|c| c.is_ascii_digit())
        {
            return Err(refuse());
        }
        let whole: u32 = if whole.is_empty() {
            0
        } else {
            whole.parse().map_err(|_| refuse())?
        };
        let fraction: u32 = format!("{fraction:0<3}").parse().map_err(|_| refuse())?;
        if whole > 1_000_000 {
            return Err(refuse());
        }
        Ok(std::num::NonZeroU32::new(whole * 1_000 + fraction)
            .map(|thousandths_per_kb| Self { thousandths_per_kb }))
    }
}

impl std::fmt::Display for MarginThinning {
    /// Records a kb, as few decimals as it needs.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let thousandths = self.thousandths_per_kb.get();
        let fraction = format!("{:03}", thousandths % 1_000);
        let fraction = fraction.trim_end_matches('0');
        if fraction.is_empty() {
            write!(f, "{}", thousandths / 1_000)
        } else {
            write!(f, "{}.{fraction}", thousandths / 1_000)
        }
    }
}

/// **The default `--refine-margin-density`**, records a kb: on 2,176 tomato samples' densest
/// stretch it held 41% less memory and took half the time with the fill check's error unchanged,
/// and the simulated cohorts, sparser than it, are left as they were; at 1 a kb and below the
/// simulated cohorts begin to lose genotypes (`doc/devel/reports/implementations/margin_thinning_2026-11-23.md`).
pub const DEFAULT_MARGIN_DENSITY: &str = "5";

/// **How far from the core every record is walked**, thinning or not (§6, thinned margins): the stretch next to
/// the core says most about which samples a record's haplotypes follow. 50 kb, the near distance of
/// the margin-thinning experiment (`doc/devel/reports/implementations/margin_thinning_2026-11-23.md`).
pub const NEAR_MARGIN_BASES: u64 = 50_000;

/// **The default shortest margin**, and with it the shortest core (`--refine-window-floor`): §4.2's
/// `S` is at most 300 kb, and a margin shorter than the walk's reach cuts off evidence a record's
/// genotype is read from.
pub const DEFAULT_MARGIN_FLOOR_BASES: u64 = 500_000;

/// **The longest a margin is**, and so the most `--refine-window-floor` may be.
pub const DEFAULT_MARGIN_BASES: u64 = 1_000_000;

/// **The default window budget**: 32 GB. At about 250 bytes a sample a record (400 tomato samples,
/// pass spec §6), 2,176 samples fit a window of about 58,000 records, 2.2 Mb at 26 records a kb.
pub const DEFAULT_WINDOW_BYTES: u64 = 32_000_000_000;

/// **What each record costs a window, in bytes, as running totals** over the caller's records:
/// entry `i` is the total over records `0..i`, so each slice is one longer than the records.
/// A window pays `held` for every record it spans, margins included, and `core` again for every
/// record of its core — what the solve keeps only for the records it writes. With thinned margins,
/// `skipped` is what `held` counts for a record where the walks skip the record and a window
/// holding it does not walk it — its solve state, and its numbers and reads — zero for a record
/// they read; a window does not pay it beyond its near margins. `packed` is what such a record
/// holds there instead, its numbers and reads compressed; a window pays it only beyond its near
/// margins, and none where it is `None`.
/// `outside` is what a record holds while no window reads it; a window pays it for every record
/// the caller holds outside the window, before or after it. `None` where the caller holds only
/// what the window spans.
#[derive(Clone, Copy, Debug)]
pub struct WindowBytes<'a> {
    pub held: &'a [u64],
    pub core: &'a [u64],
    pub skipped: Option<&'a [u64]>,
    pub packed: Option<&'a [u64]>,
    pub outside: Option<&'a [u64]>,
    /// The two parts of `held` and `core` the solve never holds at once; `None` counts both.
    pub store_and_results: Option<StoreAndResults<'a>>,
}

/// **Two parts of a window's count the solve never holds at once**, as running totals like
/// [`WindowBytes`]'s: the rounds' haplotype store, part of `held` at every record, which the solve
/// lets go when the rounds end, and the core's written results, part of `core`, which it makes
/// after. A window is counted at the larger of the store it walks and the results its core
/// writes, not their sum (pass spec §6).
#[derive(Clone, Copy, Debug)]
pub struct StoreAndResults<'a> {
    /// What `held` counts for the rounds' store at each record.
    pub store: &'a [u64],
    /// What a record skipped beyond the near margins does not hold of `store` — its store, as part
    /// of `WindowBytes::skipped`; zero for a record walked there. `None` without thinning.
    pub store_skipped: Option<&'a [u64]>,
    /// What `core` counts for the written results at each record.
    pub results: &'a [u64],
}

impl Default for WindowSizes {
    fn default() -> Self {
        Self {
            core_bases: 2_000_000,
            margin_bases: DEFAULT_MARGIN_BASES,
            margin_floor_bases: DEFAULT_MARGIN_FLOOR_BASES,
            core_floor_bases: DEFAULT_MARGIN_FLOOR_BASES,
            budget_bytes: DEFAULT_WINDOW_BYTES,
            thinning: None,
        }
    }
}

/// **One window**, as indices into the caller's records of one contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowPlan {
    /// The first record of the left margin.
    pub context_start: usize,
    /// The core is `core_start..core_end`.
    pub core_start: usize,
    pub core_end: usize,
    /// One past the last record of the right margin.
    pub context_end: usize,
    /// The records every walk reads, thinning or not: the core and its near margins
    /// ([`NEAR_MARGIN_BASES`]), `near_start..near_end`. The whole window without thinning.
    pub near_start: usize,
    pub near_end: usize,
    /// The core's length and each margin's, in bases, after the budget.
    pub core_bases: u64,
    pub margin_bases: u64,
    /// The bytes the window holds, as [`WindowBytes`] counts them, the records held outside it
    /// included.
    pub bytes: u64,
    /// Of `bytes`, what the records held outside the window cost ([`WindowBytes::outside`]).
    pub outside_bytes: u64,
    /// The bytes it holds past the budget; zero within it.
    pub over_budget: u64,
}

impl WindowPlan {
    /// How many records the window holds, margins included.
    #[must_use]
    pub fn records(&self) -> usize {
        self.context_end - self.context_start
    }
}

/// **The window whose core starts at record `core_start`**, or `None` when the caller has not yet
/// read far enough to decide it.
///
/// `positions` are the `POS` of the contig's records the caller holds, in arrival order and so
/// non-decreasing; they must reach back at least `sizes.margin_bases` before `positions[core_start]`
/// (or to the contig's start). `bytes` are those records' costs. `contig_complete` says the contig
/// has no records after the last one given. The next window's core starts at the returned
/// `core_end`.
///
/// # Panics
///
/// If `core_start` is not an index into `positions`, or `bytes` is not one entry longer than it.
#[must_use]
pub fn plan_window(
    positions: &[u64],
    bytes: WindowBytes<'_>,
    core_start: usize,
    contig_complete: bool,
    sizes: &WindowSizes,
) -> Option<WindowPlan> {
    assert!(
        bytes.held.len() == positions.len() + 1
            && bytes.core.len() == positions.len() + 1
            && bytes
                .skipped
                .is_none_or(|skipped| skipped.len() == positions.len() + 1)
            && bytes
                .packed
                .is_none_or(|packed| packed.len() == positions.len() + 1)
            && bytes.store_and_results.is_none_or(|apart| {
                apart.store.len() == positions.len() + 1
                    && apart.results.len() == positions.len() + 1
                    && apart
                        .store_skipped
                        .is_none_or(|skipped| skipped.len() == positions.len() + 1)
            })
            && bytes
                .outside
                .is_none_or(|outside| outside.len() == positions.len() + 1),
        "one running total a record, and one before them"
    );
    debug_assert!(
        positions.windows(2).all(|pair| pair[0] <= pair[1]),
        "a contig's positions must not descend"
    );
    let start = positions[core_start];
    let last = *positions.last().expect("core_start indexes positions");
    let measure = |core_bases: u64, margin_bases: u64| {
        let context_start = positions[..core_start]
            .partition_point(|&position| position < start.saturating_sub(margin_bases));
        let core_end = core_start
            + positions[core_start..].partition_point(|&position| position < start + core_bases);
        let context_end =
            positions.partition_point(|&position| position < start + core_bases + margin_bases);
        let (near_start, near_end) = match (sizes.thinning, bytes.skipped) {
            (Some(_), Some(_)) => {
                let near = NEAR_MARGIN_BASES.min(margin_bases);
                (
                    context_start
                        + positions[context_start..core_start]
                            .partition_point(|&position| position < start.saturating_sub(near)),
                    positions.partition_point(|&position| position < start + core_bases + near),
                )
            }
            _ => (context_start, context_end),
        };
        // A record's skipped cost is at most what `held` counts for it, so the difference only
        // grows as the window reads further.
        let beyond_near = |totals: &[u64]| {
            totals[near_start] - totals[context_start] + totals[context_end] - totals[near_end]
        };
        let skipped = bytes.skipped.map_or(0, beyond_near);
        let packed = bytes.packed.map_or(0, beyond_near);
        // What the caller holds before the window and after it. Reading further only adds to it.
        let outside_bytes = bytes.outside.map_or(0, |outside| {
            outside[context_start] - outside[0] + outside[positions.len()] - outside[context_end]
        });
        // The store the walks hold and the results the core writes are never held together, so
        // the smaller is not paid. Reading further only adds: a record adds at least its store to
        // `held` and at least its result to `core`.
        let never_together = bytes.store_and_results.map_or(0, |apart| {
            let store = apart.store[context_end]
                - apart.store[context_start]
                - apart.store_skipped.map_or(0, beyond_near);
            let results = apart.results[core_end] - apart.results[core_start];
            store.min(results)
        });
        let held = bytes.held[context_end] - bytes.held[context_start] - skipped
            + packed
            + bytes.core[core_end]
            - bytes.core[core_start]
            - never_together;
        let held = held + outside_bytes;
        let plan = WindowPlan {
            context_start,
            core_start,
            core_end,
            context_end,
            near_start,
            near_end,
            core_bases,
            margin_bases,
            bytes: held,
            outside_bytes,
            over_budget: held.saturating_sub(sizes.budget_bytes),
        };
        // Every record the window reaches has been read once a later one has, or none will come.
        let read_far_enough = contig_complete || last >= start + core_bases + margin_bases;
        (plan, read_far_enough)
    };
    for (core_bases, margin_bases) in size_ladder(sizes) {
        let (plan, read_far_enough) = measure(core_bases, margin_bases);
        if plan.over_budget > 0 {
            // Already over on what has been read, and reading more only adds.
            continue;
        }
        if !read_far_enough {
            return None;
        }
        return Some(plan);
    }
    let (plan, read_far_enough) = measure(core_floor(sizes), sizes.margin_floor_bases);
    read_far_enough.then_some(plan)
}

/// The shortest core: the core floor, but no longer than the margin floor or the core itself.
fn core_floor(sizes: &WindowSizes) -> u64 {
    sizes
        .core_floor_bases
        .min(sizes.margin_floor_bases)
        .min(sizes.core_bases)
}

/// The sizes tried, largest first, in steps of a tenth of the margin floor: the core shrinking to
/// the margins' length, then the two together down to the margin floor, then the core alone down
/// to the core floor.
fn size_ladder(sizes: &WindowSizes) -> Vec<(u64, u64)> {
    let step = (sizes.margin_floor_bases / 10).max(1);
    let mut ladder = Vec::new();
    let mut core = sizes.core_bases;
    while core > sizes.margin_bases {
        ladder.push((core, sizes.margin_bases));
        core = core.saturating_sub(step);
    }
    let mut both = sizes.margin_bases.min(sizes.core_bases);
    while both >= sizes.margin_floor_bases {
        ladder.push((both, both));
        if both < sizes.margin_floor_bases + step {
            break;
        }
        both -= step;
    }
    // The core alone, below the margin floor: the margins keep the walk's reach, and a short core
    // only means more windows.
    let floor = core_floor(sizes);
    let mut core = sizes.margin_floor_bases.saturating_sub(step);
    while core > floor {
        if core < sizes.core_bases {
            ladder.push((core, sizes.margin_floor_bases));
        }
        core = core.saturating_sub(step);
    }
    if floor < sizes.margin_floor_bases {
        ladder.push((floor, sizes.margin_floor_bases));
    }
    ladder
}

#[cfg(test)]
mod tests;
