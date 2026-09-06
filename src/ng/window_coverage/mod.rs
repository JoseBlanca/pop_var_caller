//! Each sample's read depth and GC fraction over the 500-base window centred on a position.
//!
//! One computation with two consumers, both fed by the same accumulator so that they cannot
//! drift apart (spec [`window_coverage.md`](../../../doc/devel/ng/spec/window_coverage.md) §1):
//!
//! - **a per-locus pair** — at every locus the caller writes, each sample's
//!   [`WindowCoverage`] at the locus's first base, handed to the hidden-duplication filter
//!   beside the record;
//! - **a per-sample histogram** — every emitted pair binned by GC and by depth
//!   ([`CoverageByGcHistogram`]), complete when the calling pass ends, from which that filter
//!   fits what one copy's depth looks like in this sample.
//!
//! The pair is the measurement and the histogram is the yardstick; the filter divides one by
//! the other, so both have to be the same quantity measured the same way. Nothing here scores
//! or fits anything.
//!
//! **The sliding window and its histogram fold are production's**
//! ([`SlidingWindowCoverageAccumulator`](crate::sample_summary::coverage::SlidingWindowCoverageAccumulator),
//! transcribed with its tests under the freeze rule that `src/sample_summary/` is not edited).
//! The copy is deliberate and is held to its original by a differential test
//! (`production_parity.rs`) — a later reader who is tempted to delete it and call
//! production instead should know that ng's version is about to diverge: a floor under how
//! many covered positions a window may speak for, and a depth bin width fitted to the sample,
//! arrive in the next two steps of the plan.

mod accumulator;

#[cfg(test)]
mod production_parity;

pub use accumulator::WindowCoverageAccumulator;

/// Window width in bases — the analysis window and the GC covariate window are the same 500
/// bases, production's `--gc-window-bp` default, inherited rather than re-measured (spec §3.3).
///
/// **The look-ahead a cover has to keep is half of this.** A centre at `p` is complete only
/// once the stream has passed `p + WINDOW_BP / 2`, and the cover that draws each sample forward
/// has to reach that far past a building region or every region's last centres finalise after
/// their builder has run — silently, as an absent sample at exactly the loci nearest each
/// region boundary (spec §3.3, §6 trap 3). Both sites derive their number from here.
pub const WINDOW_BP: u32 = 500;

/// GC bins, uniform over the closed unit interval `[0, 1]` (spec §3.4).
///
/// The depth axis has no constant here on purpose: its bin count and its width are marked soft
/// in the spec and are settled by measurement, in plan steps D1 and D2.
pub const GC_BINS: u32 = 50;

/// How wide the window is and how the histogram's cells are cut.
///
/// Production's `CoverageBinScheme` transcribed, under the name spec §3.6 gives it. Every
/// field must be positive (`depth_bin_width` finite and `> 0`); [`assert_valid`] asserts that
/// in both debug and release rather than silently using a scheme the caller never chose.
///
/// **There is no `Default`**, and the settled values live in [`WINDOW_BP`] and [`GC_BINS`]
/// instead: every field here changes what the run measures, so the configuration is spelled at
/// the site that builds the accumulator — plan step C2, inside the merge's observation cache —
/// and a misconfigured run fails there rather than quietly measuring something else.
///
/// [`assert_valid`]: WindowCoverageConfig::assert_valid
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowCoverageConfig {
    /// Window width in bases. The window centred on `p` spans the covered positions in
    /// `[p − window_bp / 2, p + window_bp / 2]`. Settled at [`WINDOW_BP`].
    pub window_bp: u32,
    /// Number of GC bins, uniform over the closed unit interval `[0, 1]`. Settled at
    /// [`GC_BINS`].
    pub gc_bins: u32,
    /// Width of one depth bin, in mean-depth units.
    pub depth_bin_width: f64,
    /// Number of regular depth bins; one overflow bin follows them.
    pub depth_bins: u32,
}

impl WindowCoverageConfig {
    /// Cells per GC row: the regular depth bins plus one overflow bin.
    fn depth_columns(&self) -> usize {
        self.depth_bins as usize + 1
    }

    /// Total histogram cells — `gc_bins * (depth_bins + 1)`.
    fn cell_count(&self) -> usize {
        self.gc_bins as usize * self.depth_columns()
    }

    /// Panic (in both debug and release) if any field is non-positive, or `depth_bin_width`
    /// is not finite. The configuration is chosen by the run, not by the data, so a bad one
    /// is a programmer error and fails loudly rather than quietly substituting another — the
    /// `NaN` width is the case worth naming, because it would send every window into depth
    /// bin 0 and produce a plausible-looking histogram from a configuration nobody chose.
    /// The boundary that supplies the configuration is plan step C2, where the merge's
    /// observation cache builds one accumulator per sample.
    pub(crate) fn assert_valid(&self) {
        assert!(self.window_bp >= 1, "window_bp must be >= 1");
        assert!(self.gc_bins >= 1, "gc_bins must be >= 1");
        assert!(self.depth_bins >= 1, "depth_bins must be >= 1");
        assert!(
            self.depth_bin_width.is_finite() && self.depth_bin_width > 0.0,
            "depth_bin_width must be finite and > 0, got {}",
            self.depth_bin_width,
        );
    }

    /// Map a window's `(GC fraction in [0, 1], mean depth >= 0)` to a row-major cell index.
    /// GC saturates into the last GC bin at 1.0; a depth at or above
    /// `depth_bins * depth_bin_width` lands in the overflow column. `mean_depth` is always
    /// `>= 0` (a sum of non-negative counts over a positive divisor), so `0` maps to depth
    /// bin `0`; the `as usize` cast also saturates a hypothetical negative to `0`.
    fn cell_index(&self, gc_fraction: f64, mean_depth: f64) -> usize {
        let gc_bins = self.gc_bins as usize;
        let gc_bin = ((gc_fraction * gc_bins as f64) as usize).min(gc_bins - 1);
        let depth_bins = self.depth_bins as usize;
        let depth_bin = ((mean_depth / self.depth_bin_width) as usize).min(depth_bins);
        gc_bin * self.depth_columns() + depth_bin
    }
}

/// One sample's window at one position: the mean read depth and the GC fraction over the
/// covered positions of the window centred there.
///
/// **Which position it belongs to is not in here.** The accumulator hands the pair back
/// beside its centre ([`WindowCoverageAccumulator::pop_ready`]), and the cache keys it by
/// that centre — spec §3.5's shape, so that a locus can carry one pair per sample without
/// repeating a coordinate every consumer already knows.
///
/// **From the next step of the plan on, both fields are `NaN` where the sample has no usable
/// window there** — a value, not an error. Nothing at this step produces one: a window always
/// holds at least its own centre, so neither quotient can be `NaN`. [`PartialEq`] compares bit
/// patterns *ahead* of that need rather than after it, because a derived one would make an
/// absent window unequal to itself the moment the floor lands, and the tests that would then
/// start failing are the ones that compare two runs — a failure that would arrive a step away
/// from its cause. Production's `LocusWindowCoverage`
/// ([`types.rs:233`](../../../src/var_calling/types.rs)) compares the same way and for the
/// same reason, and withholds `Eq`/`Hash` as this does: bitwise equality separates `+0.0` from
/// `-0.0`, which is right for comparing two runs and wrong for a lookup key.
#[derive(Debug, Clone, Copy)]
pub struct WindowCoverage {
    /// Fraction of the window's covered positions whose reference base is `G` or `C`.
    pub gc_fraction: f32,
    /// Mean read depth over the window's covered positions.
    pub mean_depth: f32,
}

impl PartialEq for WindowCoverage {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            gc_fraction,
            mean_depth,
        } = self;
        gc_fraction.to_bits() == other.gc_fraction.to_bits()
            && mean_depth.to_bits() == other.mean_depth.to_bits()
    }
}

/// One sample's finished coverage-by-GC histogram: a row-major `[gc_bin][depth_bin]` count
/// matrix over every window the sample finalised, plus the bin scheme needed to read a cell.
///
/// These are sufficient statistics, not a fitted model — the curve and the single-copy scale
/// are the filter's own work. **The four scheme fields are how a consumer turns a cell index
/// back into a depth**, so they are not decoration: a `depth_bin_width` that does not match the
/// one the cells were cut with makes every fitted single-copy depth wrong by a constant factor,
/// with nothing to signal it.
///
/// **Production's shape with two fields dropped** (spec §7): `n_skipped_tiles`, which the
/// sliding model can only ever report as zero, and `callable_positions`, which is the same
/// counter as `windows_folded` here — production reads it in one place,
/// `var_calling::diversity`, as the denominator of a cohort diversity estimate, a number ng
/// obtains elsewhere. Neither is read by the coverage model fit.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageByGcHistogram {
    /// Window width in bases — the GC covariate window, equal to the analysis window.
    pub window_bp: u32,
    /// Number of GC bins, uniform over `[0, 1]`: a window of GC fraction `g` lands in bin
    /// `min(floor(g * gc_bins), gc_bins - 1)`.
    pub gc_bins: u32,
    /// Width of one depth bin, in mean-depth units.
    pub depth_bin_width: f64,
    /// Number of regular depth bins. One overflow bin — for a mean depth at or above
    /// `depth_bins * depth_bin_width` — follows them, so each GC row holds `depth_bins + 1`
    /// cells.
    pub depth_bins: u32,
    /// How many finalised windows were folded in. Every covered position finalises exactly
    /// one window, so this also counts the sample's covered positions.
    pub windows_folded: u64,
    /// Row-major `[gc_bin][depth_bin]` counts. Length is exactly `gc_bins * (depth_bins + 1)`.
    ///
    /// **A cell saturates at `u32::MAX` rather than wrapping.** One window is folded per
    /// covered position, so a cell can only reach that on a reference above about 4.3 Gbp —
    /// larger than tomato or human, not larger than the plant genomes this caller is meant to
    /// take. A saturated cell under-reports; a wrapped one would report a near-empty cell
    /// where the single-copy peak is, and the filter's fit anchors on exactly that mode.
    pub counts: Vec<u32>,
}
