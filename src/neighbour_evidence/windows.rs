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
//! the margins' length, and then core and margins shrink together, down to the margin floor of
//! 500 kb. A window still over budget at the floor is planned at the floor and says by how much it
//! is over, and the stage stops the run there: margins shorter than the floor would cut off the
//! evidence a record's walk reaches for (§4.2's `S` is at most 300 kb).
//!
//! The budget does not count the walks' own arrays, which every thread holds one of: counting
//! them would make the windows, and so the output, depend on the thread count (§2.3).

/// **How large a window may be** (§2.2, §6).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowSizes {
    /// The core's length, in bases, before the budget shrinks it.
    pub core_bases: u64,
    /// Each margin's length, in bases, before the budget shrinks it.
    pub margin_bases: u64,
    /// The shortest a margin may become (and, with it, the core).
    pub margin_floor_bases: u64,
    /// The most bytes a window should hold, as [`WindowBytes`] counts them
    /// (`--refine-window-bytes`).
    pub budget_bytes: u64,
}

/// **The default window budget**: 32 GB. At about 250 bytes a sample a record (400 tomato samples,
/// pass spec §6), 2,176 samples fit a window of about 58,000 records, 2.2 Mb at 26 records a kb.
pub const DEFAULT_WINDOW_BYTES: u64 = 32_000_000_000;

/// **What each record costs a window, in bytes, as running totals** over the caller's records:
/// entry `i` is the total over records `0..i`, so each slice is one longer than the records.
/// A window pays `held` for every record it spans, margins included, and `core` again for every
/// record of its core — what the solve keeps only for the records it writes.
#[derive(Clone, Copy, Debug)]
pub struct WindowBytes<'a> {
    pub held: &'a [u64],
    pub core: &'a [u64],
}

impl Default for WindowSizes {
    fn default() -> Self {
        Self {
            core_bases: 2_000_000,
            margin_bases: 1_000_000,
            margin_floor_bases: 500_000,
            budget_bytes: DEFAULT_WINDOW_BYTES,
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
    /// The core's length and each margin's, in bases, after the budget.
    pub core_bases: u64,
    pub margin_bases: u64,
    /// The bytes the window holds, as [`WindowBytes`] counts them.
    pub bytes: u64,
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
        bytes.held.len() == positions.len() + 1 && bytes.core.len() == positions.len() + 1,
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
        let held = bytes.held[context_end] - bytes.held[context_start] + bytes.core[core_end]
            - bytes.core[core_start];
        let plan = WindowPlan {
            context_start,
            core_start,
            core_end,
            context_end,
            core_bases,
            margin_bases,
            bytes: held,
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
    let (plan, read_far_enough) = measure(sizes.margin_floor_bases, sizes.margin_floor_bases);
    read_far_enough.then_some(plan)
}

/// The sizes tried, largest first, in steps of a tenth of the margin floor: the core shrinking to
/// the margins' length, then the two together down to the floor.
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
    ladder
}

#[cfg(test)]
mod tests;
