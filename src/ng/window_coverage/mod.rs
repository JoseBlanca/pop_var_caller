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
//! (`production_parity.rs`) — a later reader who is tempted to delete it and call production
//! instead should know that ng's version has already diverged twice. A window built from fewer
//! than [`WindowCoverageConfig::min_window_positions`] covered positions reports nothing, where
//! production's would report a number; and the depth bin width is fitted to each sample rather
//! than configured, so the two sides' histograms are cut on different axes and only the windows
//! are still comparable.

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
pub const GC_BINS: u32 = 50;

/// Depth bins, before the overflow bin that follows them — **soft, and confirmed or moved by
/// measurement at plan step D2** (spec §3.4).
///
/// Production bins depth at a fixed 0.5× in 2,000 bins to 1,000×, which is 400 kB a sample
/// against a 500 kB per-open-sample budget, and which serves one end of the depth axis at a
/// time: at three reads a position a 0.5× bin is a sixth of the single-copy peak's own scatter,
/// and at 300 reads the range has to reach 1,200× for a four-copy carrier. ng bins fewer and
/// scales the width to the sample instead — 50 × 401 × 4 bytes, **80.2 kB a sample**.
pub const DEPTH_BINS: u32 = 400;

/// How many windows a sample buffers before it sets its depth bin width — **soft, and confirmed
/// or moved by measurement at plan step D2** (spec §3.4). Until then, 10,000.
///
/// The width is the sample's own: the accumulator holds its first windows back, takes the median
/// of their mean depths, and sets the width so the bins span ten times that median. A sample
/// that finalises fewer than this many sets the width from what it has when the pass ends.
pub const DEPTH_SCALE_WINDOWS: u32 = 10_000;

/// How many times the sample's median window depth the regular bins are to span — **soft, and
/// confirmed or moved by measurement at plan step D2** (spec §3.4).
///
/// Ten leaves room above the single-copy peak for a four-copy carrier and below it for a sample
/// whose median sits above its own mode. What share of windows still overflow is D2's
/// measurement, on both benchmarks, against the fit's own rejection guard at a fifth.
///
/// **A multiple, not a count**, so that D2 can land between two whole numbers: the width the
/// bins are cut at is this times the sample's median, over [`DEPTH_BINS`].
pub const DEPTH_RANGE_IN_MEDIANS: f64 = 10.0;

/// How many covered positions a window must hold before it is allowed to speak — **soft, and
/// set by measurement at plan step D1** (spec §3.3's open question). Until then, 50 of 500.
///
/// A window built from a handful of positions looks exactly as confident as one built from five
/// hundred, and in ng the holes are not scattered: a repeat-tract region that emits no loci, an
/// analysed-region edge, a stretch no read reached. Below this many, the window comes back
/// absent instead.
pub const MIN_WINDOW_POSITIONS: u32 = 50;

/// How wide the window is, how many cells the histogram has, and how its depth axis is scaled.
///
/// Production's `CoverageBinScheme`, under the name spec §3.6 gives it, with three of ng's own
/// additions in place of production's fixed depth bin width: a floor under how thin a window may
/// be, the size of the sample the depth width is fitted from, and how far past that sample's
/// median the regular bins are to reach. Every field must be positive; [`assert_valid`] asserts
/// that in both debug and release rather than silently using a scheme the caller never chose.
///
/// **The width of a depth bin is not here, because it is not configuration.** It is fitted from
/// each sample's own depth as the pass runs (spec §3.4), so it arrives on the finished histogram
/// and nowhere else.
///
/// **There is no `Default`**, and the values live in the constants above instead: every field
/// here changes what the run measures, so the configuration is spelled at the site that builds
/// the accumulator — plan step C2, inside the merge's observation cache — and a misconfigured
/// run fails there rather than quietly measuring something else.
///
/// **The constants are not all of the same kind.** [`WINDOW_BP`] and [`GC_BINS`] are
/// production's, inherited and not re-argued here. [`MIN_WINDOW_POSITIONS`], [`DEPTH_BINS`],
/// [`DEPTH_SCALE_WINDOWS`] and [`DEPTH_RANGE_IN_MEDIANS`] are ng's own and are provisional: a
/// reader who finds a run's behaviour turning on one of them is looking at a value nobody has
/// yet defended with data, and plan steps D1 and D2 are where that happens.
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
    /// Number of regular depth bins; one overflow bin follows them. Provisionally
    /// [`DEPTH_BINS`].
    pub depth_bins: u32,
    /// How many windows to buffer before fitting this sample's depth bin width from their
    /// median. Provisionally [`DEPTH_SCALE_WINDOWS`].
    pub depth_scale_windows: u32,
    /// How many times the sample's fitted median the regular bins span. Provisionally
    /// [`DEPTH_RANGE_IN_MEDIANS`].
    pub depth_range_in_medians: f64,
    /// Fewest covered positions a window may be built from and still report a number; below
    /// this it comes back absent. Provisionally [`MIN_WINDOW_POSITIONS`] — see spec §3.3.
    ///
    /// **This is ng's one departure from production's window**, which lets a single covered
    /// position emit a window over itself alone.
    pub min_window_positions: u32,
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

    /// Panic (in both debug and release) if any field is non-positive, or
    /// `depth_range_in_medians` is not finite. The configuration is chosen by the run, not by
    /// the data, so a bad one is a programmer error and fails loudly rather than quietly
    /// substituting another.
    ///
    /// **The `NaN` case is the one worth naming.** `depth_range_in_medians` is the only float
    /// here, and a `NaN` one makes the fitted depth bin width `NaN` — which is exactly the value
    /// the fit reads as "this sample cannot be scaled", so a misconfigured run would cost every
    /// sample its histogram and say nothing.
    ///
    /// The boundary that supplies the configuration is plan step C2, where the merge's
    /// observation cache builds one accumulator per sample.
    pub(crate) fn assert_valid(&self) {
        assert!(self.window_bp >= 1, "window_bp must be >= 1");
        assert!(self.gc_bins >= 1, "gc_bins must be >= 1");
        assert!(self.depth_bins >= 1, "depth_bins must be >= 1");
        assert!(
            self.depth_scale_windows >= 1,
            "depth_scale_windows must be >= 1",
        );
        assert!(
            self.depth_range_in_medians.is_finite() && self.depth_range_in_medians > 0.0,
            "depth_range_in_medians must be finite and > 0, got {}",
            self.depth_range_in_medians,
        );
        assert!(
            self.min_window_positions >= 1,
            "min_window_positions must be >= 1",
        );
        // The widest a window ever gets is its own centre plus half a window either side, so a
        // floor above that many bases can never be met and every window of every sample would
        // come back absent — silently, which is the one thing this configuration's docs say it
        // will not do.
        let widest_window_in_bases = 2 * u64::from(self.window_bp / 2) + 1;
        assert!(
            u64::from(self.min_window_positions) <= widest_window_in_bases,
            "min_window_positions {} is more positions than a {}-base window can hold ({}), so \
             every window would be absent",
            self.min_window_positions,
            self.window_bp,
            widest_window_in_bases,
        );
    }

    /// Map a window's `(GC fraction in [0, 1], mean depth >= 0)` to a row-major cell index,
    /// against a depth bin width fitted to this sample.
    ///
    /// GC saturates into the last GC bin at 1.0; a depth at or above `depth_bins * width` lands
    /// in the overflow column. `mean_depth` is always `>= 0` (a sum of non-negative counts over
    /// a positive divisor), so `0` maps to depth bin `0`; the `as usize` cast also saturates a
    /// hypothetical negative to `0`.
    fn cell_index(&self, gc_fraction: f64, mean_depth: f64, depth_bin_width: f64) -> usize {
        let gc_bins = self.gc_bins as usize;
        let gc_bin = ((gc_fraction * gc_bins as f64) as usize).min(gc_bins - 1);
        let depth_bins = self.depth_bins as usize;
        let depth_bin = ((mean_depth / depth_bin_width) as usize).min(depth_bins);
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
/// **Both fields are `NaN` where the sample has no usable window there** — a value, not an
/// error, built and tested through [`absent`](Self::absent) and [`is_absent`](Self::is_absent)
/// so that the two numbers can only go missing together: a half-absent pair would divide a real
/// depth by a missing yardstick. [`PartialEq`] compares bit patterns for that reason — a derived
/// one would make an absent window unequal to itself and break every comparison of two runs.
/// Production's `LocusWindowCoverage`
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

impl WindowCoverage {
    /// The value for a position where the sample has no usable window: no covered record
    /// there, or a window under the configured floor.
    ///
    /// **Both numbers go absent together**, which is why this is built here rather than written
    /// as a pair of `f32::NAN`s at each site that needs one.
    #[must_use]
    pub fn absent() -> Self {
        Self {
            gc_fraction: f32::NAN,
            mean_depth: f32::NAN,
        }
    }

    /// Whether either number is missing.
    ///
    /// **Ask this rather than comparing against [`absent`](Self::absent).** Equality here is
    /// bitwise, which is what comparing two runs needs and what asking "is this missing" does
    /// not: a `NaN` that arrived through a codec or an `f64` round trip carries a different
    /// payload and compares unequal, while this accepts it.
    ///
    /// **Fail-closed on either field.** The two numbers are meant to go missing together, and
    /// the accumulator only ever produces them that way; but both fields are `pub`, so a depth
    /// without the GC fraction that says what depth to expect there is representable, and it is
    /// not a usable measurement. Answering "absent" for it is the safe half of the answer — and
    /// this is a question, so it reports rather than panicking.
    #[must_use]
    pub fn is_absent(self) -> bool {
        self.mean_depth.is_nan() || self.gc_fraction.is_nan()
    }
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
    /// Width of one depth bin, in mean-depth units — **fitted to this sample, not configured**.
    /// It is the sample's median window depth times `depth_range_in_medians`, over `depth_bins`,
    /// so the regular bins span ten times the median. Always finite and greater than zero: a
    /// sample whose median came out zero or worse has no histogram at all.
    pub depth_bin_width: f64,
    /// Number of regular depth bins. One overflow bin — for a mean depth at or above
    /// `depth_bins * depth_bin_width` — follows them, so each GC row holds `depth_bins + 1`
    /// cells.
    pub depth_bins: u32,
    /// How many windows were folded in. Every covered position finalises exactly one window,
    /// but a window under the configured floor comes back absent and trains no yardstick, so
    /// this counts the sample's covered positions **minus** the ones whose window was too
    /// sparse to speak.
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
