//! The centred sliding window, transcribed from production's
//! `SlidingWindowCoverageAccumulator` (`src/sample_summary/coverage.rs`) with its tests.
//!
//! Production is frozen, so this is a copy rather than a call. What changed in the copying is
//! only the vocabulary: ng's [`ContigId`] / [`Position`] coordinates instead of a pair of
//! `u32`s, the emitted pair carried beside its centre instead of holding a coordinate of its
//! own, longer names for three fields, and the two histogram fields spec §7 drops. The window
//! arithmetic — the two-pointer buffer, the finalisation frontier, the contig reset — is
//! unchanged, and the tests below are production's, which is what makes that claim checkable.

use std::collections::VecDeque;

use super::{CoverageByGcHistogram, SampleHistogram, WindowCoverage, WindowCoverageConfig};
use crate::ng::types::{ContigId, GenomePosition, Position};

/// A covered position retained in the sliding-window buffer.
#[derive(Debug, Clone, Copy)]
struct CoveredPosition {
    position: u64,
    depth: u32,
    is_gc: bool,
}

/// What one finalised window reports, before it is narrowed to the `f32` pair a locus carries.
///
/// A named pair rather than a `(f64, f64)`: the two are the same primitive and travel together
/// through the fold, and transposing them is silent in every sense the module has — the emitted
/// pair is built separately, so the per-locus half would stay right and the differential green,
/// while every sample's depth bin width would be fitted from a median **GC fraction** and come
/// out hundreds of times too small.
#[derive(Debug, Clone, Copy)]
struct WindowMeans {
    gc_fraction: f64,
    mean_depth: f64,
}

/// Where a sample stands with the depth bin width its histogram is cut on.
///
/// **Three states, not two, and the third is why this is an enum.** A sample whose held-back
/// windows yield no positive median cannot be given a width — and must not be given a *later*
/// one either: refitting from a later stretch of the genome would hand the filter a yardstick
/// built from part of the sample, with the windows before it folded nowhere and the count
/// agreeing with the cells. [`Unfittable`](Self::Unfittable) latches that, which also keeps the
/// held-back list bounded.
#[derive(Debug, Clone)]
enum DepthBinWidth {
    /// Still collecting windows to fit the width from. Absent windows never enter the list: a
    /// window too thin to report a depth cannot help decide what a depth means either.
    AwaitingWindows(Vec<WindowMeans>),
    /// Fitted, and finite and greater than zero.
    Fitted(f64),
    /// The windows that arrived had no positive median. This sample has no histogram, and
    /// nothing later can change that.
    Unfittable,
}

/// The median of the held-back windows' mean depths, or `0.0` if there are none.
///
/// **The upper of the two middles on an even count**, rather than their average: either is a
/// defensible median and this one needs no arithmetic on the values, so the width a sample gets
/// is a depth the sample actually reported.
///
/// Reordering the slice is why it is taken by `&mut`; the caller is finished with the order.
/// `total_cmp` rather than `partial_cmp`: it is total, so there is nothing to unwrap — and no
/// `NaN` can be here anyway, since an absent window is never held back.
fn median_depth(windows: &mut [WindowMeans]) -> f64 {
    if windows.is_empty() {
        return 0.0;
    }
    let middle = windows.len() / 2;
    windows.select_nth_unstable_by(middle, |a, b| a.mean_depth.total_cmp(&b.mean_depth));
    windows[middle].mean_depth
}

/// Accumulates, for every covered position, the mean depth and GC fraction of the **centred**
/// window of width `window_bp` over the covered positions in it — and folds each finalised
/// window into the per-sample histogram, so one pass serves both the per-locus pair and the
/// filter's model fit.
///
/// Fed one sample's covered positions in coordinate order through [`observe`](Self::observe).
/// A position `p` is **finalised** — its centred window can gain no more positions — once the
/// stream reaches `p + window_bp / 2`, at which point its [`WindowCoverage`] is queued for
/// [`pop_ready`](Self::pop_ready) and folded into the histogram. **Unless the window holds
/// fewer than `min_window_positions` distinct covered positions**, in which case it is queued
/// absent and folded nowhere: a window too thin to report a depth is also too thin to train the
/// yardstick a depth is judged against. [`finish`](Self::finish)
/// finalises the tail (windows truncated at the last observed position or a contig end) and
/// returns it together with the histogram in one consuming call, so the tail cannot be
/// silently dropped.
///
/// Memory is `O(window_bp)` **for a stream with one record per position**, which is what the
/// walk produces; the order rule below is only non-decreasing, and a stream that repeated one
/// position without advancing would buffer without bound. On top of that sits the histogram,
/// and — until the depth bin width is fitted — the `depth_scale_windows` windows held back to
/// fit it from. Sixteen bytes each, but the list grows by doubling, so at the configured 10,000
/// windows its capacity reaches 16,384 and the transient peaks at **262 kB**, not 160.
///
/// **The depth bin width is the sample's own, and nothing is folded until it is known.** The
/// first `depth_scale_windows` windows that clear the floor are held back; their median mean
/// depth sets the width so the regular bins span `depth_range_in_medians` times it; then those
/// windows and every later one are folded. Production's fixed 0.5× bin serves one end of the
/// depth axis at a time — at three reads a position it is a sixth of the single-copy peak's own
/// scatter, and at 300 reads the range has to reach 1,200× for a four-copy carrier — which is
/// what scaling per sample avoids.
///
/// A window never spans contigs — a contig change finalises the previous contig's tail and
/// resets the sliding state.
///
/// # Invariants (the two-pointer window state, held between `observe` calls)
///
/// The covered positions of the current contig are appended to `positions` in coordinate
/// order; it holds the suffix still needed — from the left edge of the next centre's window to
/// the newest observed position. Two cursors index into it, and the running sums cover a
/// contiguous prefix:
///
/// - `centre_offset` — index in `positions` of the next centre to finalise.
/// - `summed_offset` — one past the last entry folded into
///   `sum_depth`/`sum_gc`/`summed_positions`; the sums equal `positions[0..summed_offset]`.
/// - Ordering: `centre_offset <= summed_offset <= positions.len()`.
/// - When the shrink-left step pops the front (a position now left of every remaining centre's
///   window), it decrements **both** cursors in lockstep so they keep pointing at the same
///   logical entries — safe because a front is popped only when it sits **strictly** left of
///   the current centre (so it was already summed and already before the centre), keeping both
///   cursors `>= 1` at the decrement. The strictness is what makes this hold when several
///   entries share one position.
#[derive(Debug, Clone)]
pub struct WindowCoverageAccumulator {
    config: WindowCoverageConfig,
    /// Half-window in bases (`window_bp / 2`); the window centred on `p` spans the covered
    /// positions in `[p − half_window_bp, p + half_window_bp]`.
    half_window_bp: u64,
    /// Row-major `[gc_bin][depth_bin]` counts, length `config.cell_count()`. Every cell is
    /// zero until the depth bin width is fitted, because nothing can be folded before then.
    counts: Vec<u32>,
    /// The depth bin width fitted from this sample's own depth, and the windows held back to
    /// fit it from while it is not yet known.
    ///
    /// **Nothing is folded until it is known**, because a fold before the width is set has no
    /// bins to fold into. The held-back windows are the exact pairs the fold will use, so that a
    /// window folded after the fit lands in the cell it would have landed in before.
    ///
    /// **It is per sample and is never reset at a contig boundary.** A window never spans
    /// contigs, but the axis a sample's whole histogram is cut on is one axis; refitting it at
    /// each contig would give one matrix rows cut on several.
    depth_bin_width: DepthBinWidth,
    /// Windows folded into the histogram so far — one per covered position, minus the ones
    /// the floor silenced.
    windows_folded: u64,
    /// Windows the floor silenced so far. Together with `windows_folded` this is every window
    /// finalised, which is what lets `finish` tell "this sample reached nothing" from "this
    /// sample reached ground its windows were all too thin to speak for".
    windows_under_the_floor: u64,
    /// Contig currently being folded; `None` before the first position.
    contig: Option<ContigId>,
    /// Covered positions retained, in coordinate order (see the type's `# Invariants`).
    positions: VecDeque<CoveredPosition>,
    /// Index in `positions` of the next centre to finalise (see `# Invariants`).
    centre_offset: usize,
    /// One past the last entry folded into the running sums (see `# Invariants`).
    summed_offset: usize,
    /// Running sums over `positions[0..summed_offset]` — the current centre's window once the
    /// shrink-left step has removed the entries left of it.
    sum_depth: u64,
    sum_gc: u64,
    /// How many observations those sums cover — the divisor of both means, and production's.
    summed_positions: u64,
    /// How many **different** coordinates those observations sit on — the number the floor is
    /// judged against.
    ///
    /// The two differ only when a stream observes one coordinate more than once, which
    /// [`observe`](Self::observe)'s non-decreasing order rule permits. Five records at one base
    /// are one base's worth of evidence, and a floor that counted them as five would pass
    /// exactly the window it exists to silence.
    distinct_positions_summed: u64,
    /// Finalised windows awaiting the caller, oldest first.
    ready: VecDeque<(GenomePosition, WindowCoverage)>,
    /// The previous [`observe`](Self::observe)'s position, for the coordinate-order assert.
    last_observed: Option<GenomePosition>,
}

impl WindowCoverageAccumulator {
    /// Construct an empty accumulator. Panics on an invalid configuration, via
    /// [`WindowCoverageConfig::assert_valid`].
    pub fn new(config: WindowCoverageConfig) -> Self {
        config.assert_valid();
        Self {
            half_window_bp: u64::from(config.window_bp / 2),
            counts: vec![0; config.cell_count()],
            depth_bin_width: DepthBinWidth::AwaitingWindows(Vec::new()),
            config,
            windows_folded: 0,
            windows_under_the_floor: 0,
            contig: None,
            positions: VecDeque::new(),
            centre_offset: 0,
            summed_offset: 0,
            sum_depth: 0,
            sum_gc: 0,
            summed_positions: 0,
            distinct_positions_summed: 0,
            ready: VecDeque::new(),
            last_observed: None,
        }
    }

    /// Fold one covered reference position. `reference_base` is that position's reference base
    /// (any case); `depth` is the depth there.
    ///
    /// Positions must arrive in non-decreasing genome order, **and that is a release assert**.
    /// Production carries the same guard as a `debug_assert!`; this repository ships with debug
    /// assertions off, so there it would never run, and what it stops is not a panic but a wrong
    /// number: a position behind the frontier finalises no centre and joins the buffer out of
    /// place, so the windows around it average over the wrong positions, the finalised list stops
    /// ascending, and the binary search a builder reads it with is then searching unsorted ground.
    /// One comparison per covered position, which is what the cache's own two release asserts
    /// cost. An `N`
    /// reference base is **not** a covered position: it becomes no centre and contributes to
    /// no window's sums, but it still advances the finalisation frontier. Windows finalised by
    /// this position are queued for [`pop_ready`](Self::pop_ready).
    pub fn observe(
        &mut self,
        contig: ContigId,
        position: Position,
        reference_base: u8,
        depth: u32,
    ) {
        let here = GenomePosition { contig, position };
        assert!(
            self.last_observed.is_none_or(|previous| here >= previous),
            "window coverage observed out of order: {here:?} after {:?}",
            self.last_observed,
        );
        self.last_observed = Some(here);

        // Contig change: no window spans contigs, so drain the previous contig's remaining
        // centres (truncated) and reset the sliding state.
        if self.contig != Some(contig) {
            self.finalise_all();
            self.reset_contig();
            self.contig = Some(contig);
        }

        // A non-`N` covered position joins the buffer as a centre and a sum contributor; an
        // `N` position only advances the frontier below.
        if !reference_base.eq_ignore_ascii_case(&b'N') {
            self.positions.push_back(CoveredPosition {
                position: position.get(),
                depth,
                is_gc: matches!(reference_base.to_ascii_uppercase(), b'G' | b'C'),
            });
        }

        // Finalise every centre whose right edge `centre + half_window_bp` the frontier has
        // now reached: no later covered position can fall in its window, so it is complete.
        while self.centre_offset < self.positions.len()
            && self.positions[self.centre_offset]
                .position
                .saturating_add(self.half_window_bp)
                <= position.get()
        {
            self.finalise_centre();
        }
    }

    /// Pop one finalised window and the centre it belongs to, oldest first, or `None` if none
    /// is ready.
    pub fn pop_ready(&mut self) -> Option<(GenomePosition, WindowCoverage)> {
        self.ready.pop_front()
    }

    /// Finalise the stream and hand back the histogram.
    ///
    /// Every remaining centre is finalised (its window truncated at the last observed position
    /// or the contig end), and the still-queued windows — those not yet drained through
    /// [`pop_ready`](Self::pop_ready), plus the tail just finalised — come back with it.
    /// Consuming `self` and returning the tail is what stops a caller forgetting to drain.
    ///
    /// **A sample the width could not be fitted from comes back with the reason instead** — it
    /// finalised no window at all, or none that cleared the floor, or the windows it did have
    /// reported no positive median depth. The three are three different reports, and this is the
    /// only place the difference exists; a caller that does not act on it can say
    /// [`fitted`](SampleHistogram::fitted) and get an `Option` back.
    ///
    /// A sample that finalised fewer windows than `depth_scale_windows` fits its width here,
    /// from what it has.
    #[must_use = "the tail windows and the histogram are the pass's whole output; dropping \
                  them discards every window the accumulator has not already handed over"]
    pub fn finish(mut self) -> (Vec<(GenomePosition, WindowCoverage)>, SampleHistogram) {
        self.finalise_all();
        self.fit_depth_bin_width_and_fold_held_back();
        let tail: Vec<(GenomePosition, WindowCoverage)> = self.ready.into_iter().collect();
        // Exhaustive destructure: a field added to the configuration must be either carried
        // into the histogram or explicitly ignored here, rather than silently omitted.
        let WindowCoverageConfig {
            window_bp,
            gc_bins,
            depth_bins,
            // None of the three is a bin scheme. The floor decides which windows are folded,
            // the scale sample decides when the width is fitted, and the range decides how wide
            // it comes out; a consumer reading a cell needs none of them, and the width they
            // between them produced is carried below.
            min_window_positions: _,
            depth_scale_windows: _,
            depth_range_in_medians: _,
        } = self.config;
        let histogram = match self.depth_bin_width {
            DepthBinWidth::Fitted(depth_bin_width) => {
                SampleHistogram::Fitted(CoverageByGcHistogram {
                    window_bp,
                    gc_bins,
                    depth_bin_width,
                    depth_bins,
                    windows_folded: self.windows_folded,
                    windows_under_the_floor: self.windows_under_the_floor,
                    counts: self.counts,
                })
            }
            // Nothing was ever held back. Either no window was finalised at all, or every one
            // of them was refused by the floor — and those are the two the counter separates.
            DepthBinWidth::AwaitingWindows(_) if self.windows_under_the_floor == 0 => {
                SampleHistogram::NoWindowFinalised
            }
            DepthBinWidth::AwaitingWindows(_) => SampleHistogram::EveryWindowUnderTheFloor,
            DepthBinWidth::Unfittable => SampleHistogram::MedianDepthNotPositive,
        };
        (tail, histogram)
    }

    /// Finalise the centre at `centre_offset`: complete its window (extend the running sums
    /// right to `centre + half_window_bp`, shrink them left past `centre − half_window_bp`),
    /// emit the pair, and fold it into the histogram — unless the window holds fewer than
    /// `min_window_positions` distinct covered positions, in which case it emits the absent
    /// pair and folds nothing. `summed_positions >= 1` always — the centre lies in its own
    /// window.
    fn finalise_centre(&mut self) {
        let centre = self.positions[self.centre_offset].position;
        // `saturating_add` clamps a centre within `half_window_bp` of `u64::MAX`, which no
        // reference coordinate reaches; it matches the `saturating_sub` on the other edge,
        // where a centre inside the first half-window of a contig is real and common.
        let right_edge = centre.saturating_add(self.half_window_bp);
        let left_edge = centre.saturating_sub(self.half_window_bp);

        // Extend right: pull in every buffered position at or before the window's right edge
        // (all are already observed — the frontier is here).
        while self.summed_offset < self.positions.len()
            && self.positions[self.summed_offset].position <= right_edge
        {
            let entry = self.positions[self.summed_offset];
            // The buffer is sorted, so an entry starts a new coordinate exactly when it differs
            // from the one before it.
            if self.summed_offset == 0
                || self.positions[self.summed_offset - 1].position != entry.position
            {
                self.distinct_positions_summed += 1;
            }
            self.sum_depth += u64::from(entry.depth);
            self.sum_gc += u64::from(entry.is_gc);
            self.summed_positions += 1;
            self.summed_offset += 1;
        }
        // Shrink left: drop every buffered position before the window's left edge — also
        // before every later centre's window, so never needed again. The `-=` cannot
        // underflow: a popped front was summed by the extend step above and sits left of the
        // centre, so `summed_positions`, `summed_offset` and `centre_offset` are all `>= 1`
        // here (see the type's `# Invariants`).
        while let Some(front) = self.positions.front() {
            if front.position < left_edge {
                // The coordinate leaves when its last entry does. Entries sharing a coordinate
                // are always summed together — they compare equal against the right edge — so
                // the one behind the front is summed too whenever it exists.
                if self.positions.len() == 1 || self.positions[1].position != front.position {
                    self.distinct_positions_summed -= 1;
                }
                self.sum_depth -= u64::from(front.depth);
                self.sum_gc -= u64::from(front.is_gc);
                self.summed_positions -= 1;
                self.positions.pop_front();
                self.summed_offset -= 1;
                self.centre_offset -= 1;
            } else {
                break;
            }
        }

        // PANIC-FREE: `observe` sets `contig` before buffering any position, and
        // `reset_contig` deliberately leaves it set, so a centre is only ever finalised while
        // a contig is active.
        let at = GenomePosition {
            contig: self.contig.expect("centre finalised within a contig"),
            position: Position(centre),
        };
        let window = if self.distinct_positions_summed < u64::from(self.config.min_window_positions)
        {
            // **Too thin to speak, so it says nothing — and it trains nothing either.** The
            // spec (§3.3) settles that such a window comes back absent; that it is also kept
            // out of the histogram is this module's decision, and the reason is arithmetic:
            // `NaN` divided by the bin width casts to 0, so folding an absent window would pile
            // every too-sparse window into the lowest-GC, lowest-depth cell — and the yardstick
            // the filter fits would carry a spike of windows that had no depth to report.
            self.windows_under_the_floor += 1;
            WindowCoverage::absent()
        } else {
            let divisor = self.summed_positions as f64;
            let gc_fraction = self.sum_gc as f64 / divisor;
            let mean_depth = self.sum_depth as f64 / divisor;
            self.fold_or_hold_back(WindowMeans {
                gc_fraction,
                mean_depth,
            });
            WindowCoverage {
                gc_fraction: gc_fraction as f32,
                // The stored `f32` holds a *mean* depth, bounded by the per-position depth and
                // far under `f32`'s ~1.6e7 exact-integer limit, so the narrowing loses nothing
                // in practice.
                mean_depth: mean_depth as f32,
            }
        };
        // One push and one cursor advance, on both paths: `finalise_all` loops on that cursor,
        // so an exit that forgot it would hang the run rather than fail a test.
        self.ready.push_back((at, window));
        self.centre_offset += 1;
    }

    /// Fold one window into the histogram, or hold it back if the depth bin width is not yet
    /// fitted — and fit it, and fold everything held back, once enough windows have arrived.
    fn fold_or_hold_back(&mut self, window: WindowMeans) {
        match &mut self.depth_bin_width {
            DepthBinWidth::Fitted(width) => {
                let width = *width;
                self.fold(window, width);
            }
            DepthBinWidth::AwaitingWindows(held_back) => {
                held_back.push(window);
                if held_back.len() >= self.config.depth_scale_windows as usize {
                    self.fit_depth_bin_width_and_fold_held_back();
                }
            }
            // A sample already found unfittable folds nothing and holds nothing: it has no
            // histogram, and a second scale sample from a later stretch of the genome would
            // only produce a yardstick built from part of it.
            DepthBinWidth::Unfittable => {}
        }
    }

    /// Fit this sample's depth bin width from the windows held back, then fold them.
    ///
    /// **A sample whose held-back windows have no positive median gets no width and no
    /// histogram, and stays that way.** A width of zero would put every window in the overflow
    /// column and a negative or non-finite one is not a scale at all — so the sample is latched
    /// unfittable and [`finish`](Self::finish) returns `None`, the same answer as for a sample
    /// that finalised no window and for the same reason: nothing was measured that a yardstick
    /// could be built from.
    ///
    /// **Having none held back is not the same as having a bad median**, so an empty list leaves
    /// the state where it is: that is a sample which has not measured anything *yet*, and at
    /// [`finish`](Self::finish) it is a sample that never did.
    fn fit_depth_bin_width_and_fold_held_back(&mut self) {
        let DepthBinWidth::AwaitingWindows(held_back) = &mut self.depth_bin_width else {
            return;
        };
        if held_back.is_empty() {
            return;
        }
        let mut held_back = std::mem::take(held_back);
        let median = median_depth(&mut held_back);
        let width = median * self.config.depth_range_in_medians / f64::from(self.config.depth_bins);
        if !(width.is_finite() && width > 0.0) {
            self.depth_bin_width = DepthBinWidth::Unfittable;
            return;
        }
        self.depth_bin_width = DepthBinWidth::Fitted(width);
        for window in held_back {
            self.fold(window, width);
        }
    }

    /// Add one window to its cell.
    fn fold(&mut self, window: WindowMeans, depth_bin_width: f64) {
        let WindowMeans {
            gc_fraction,
            mean_depth,
        } = window;
        // Saturating rather than wrapping: a cell takes at most one count per covered
        // position, so it can only reach `u32::MAX` on a reference above about 4.3 Gbp. A
        // saturated cell under-reports; a wrapped one reports a near-empty cell where the
        // single-copy peak is, and the filter's fit anchors on exactly that mode.
        let index = self
            .config
            .cell_index(gc_fraction, mean_depth, depth_bin_width);
        let cell = &mut self.counts[index];
        *cell = cell.saturating_add(1);
        self.windows_folded += 1;
    }

    /// Finalise all remaining centres (truncated windows).
    fn finalise_all(&mut self) {
        while self.centre_offset < self.positions.len() {
            self.finalise_centre();
        }
    }

    /// Reset the per-contig sliding state, after draining a contig. `contig` is deliberately
    /// *not* cleared — the caller sets it to the new contig immediately after.
    fn reset_contig(&mut self) {
        self.positions.clear();
        self.centre_offset = 0;
        self.summed_offset = 0;
        self.sum_depth = 0;
        self.sum_gc = 0;
        self.summed_positions = 0;
        self.distinct_positions_summed = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `window_bp = 4` configuration (`half = 2`): the window centred on `p` spans covered
    /// positions in `[p − 2, p + 2]`. Coarse GC and depth bins keep the per-value asserts about
    /// the emitted pair, not about the histogram.
    fn sliding_config() -> WindowCoverageConfig {
        WindowCoverageConfig {
            window_bp: 4,
            gc_bins: 2,
            depth_bins: 40,
            // The floor off, so these keep testing production's window: at 1 every window
            // clears it, because a window always holds at least its own centre.
            min_window_positions: 1,
            // Larger than any fixture here, so the width is fitted once at `finish` from every
            // window the stream produced rather than from its first few.
            depth_scale_windows: 1_000,
            depth_range_in_medians: 10.0,
        }
    }

    /// Feed a `(position, reference base, depth)` stream on one contig and collect every
    /// emitted window (drained after each `observe`, plus the tail `finish` returns), together
    /// with the finished histogram.
    fn run_sliding_allowing_no_histogram(
        config: WindowCoverageConfig,
        contig: u32,
        stream: &[(u64, u8, u32)],
    ) -> (
        Vec<(GenomePosition, WindowCoverage)>,
        Option<CoverageByGcHistogram>,
    ) {
        let (windows, histogram) = run_sliding_with_the_reason(config, contig, stream);
        (windows, histogram.fitted())
    }

    /// The whole answer, including which silence it was.
    fn run_sliding_with_the_reason(
        config: WindowCoverageConfig,
        contig: u32,
        stream: &[(u64, u8, u32)],
    ) -> (Vec<(GenomePosition, WindowCoverage)>, SampleHistogram) {
        let mut accumulator = WindowCoverageAccumulator::new(config);
        let mut out = Vec::new();
        for &(position, base, depth) in stream {
            accumulator.observe(ContigId(contig), Position(position), base, depth);
            while let Some(window) = accumulator.pop_ready() {
                out.push(window);
            }
        }
        let (tail, histogram) = accumulator.finish();
        out.extend(tail);
        (out, histogram)
    }

    /// The same, for the fixtures that do produce a histogram — most of them.
    fn run_sliding(
        config: WindowCoverageConfig,
        contig: u32,
        stream: &[(u64, u8, u32)],
    ) -> (Vec<(GenomePosition, WindowCoverage)>, CoverageByGcHistogram) {
        let (windows, histogram) = run_sliding_allowing_no_histogram(config, contig, stream);
        (
            windows,
            histogram.expect("this fixture folds windows, so a depth width was fitted"),
        )
    }

    /// **The three silences are three different answers**, and this is where that is asserted:
    /// a sample the pass never reached, one whose every window the floor refused, and one whose
    /// windows reported no positive depth. Collapsing them would lose the middle one, which is a
    /// reading on the floor rather than a fault and is the likely case at one low-coverage
    /// sample.
    #[test]
    fn a_sample_without_a_histogram_says_which_silence_it_was() {
        let (_, nothing_reached) = run_sliding_with_the_reason(sliding_config(), 0, &[]);
        assert_eq!(nothing_reached, SampleHistogram::NoWindowFinalised);

        // Two positions 100 apart at a window of 10: each window holds only its own centre, so
        // both are one short of a floor of two.
        let silenced = [(1u64, b'G', 7u32), (101, b'G', 7)];
        let (_, all_refused) = run_sliding_with_the_reason(floored_config(2), 0, &silenced);
        assert_eq!(all_refused, SampleHistogram::EveryWindowUnderTheFloor);

        let zero_depth: Vec<(u64, u8, u32)> = (1..=5u64).map(|p| (p, b'G', 0u32)).collect();
        let (_, no_median) = run_sliding_with_the_reason(sliding_config(), 0, &zero_depth);
        assert_eq!(no_median, SampleHistogram::MedianDepthNotPositive);
    }

    /// The histogram says how many windows the floor silenced, so that a yardstick fitted from a
    /// small share of a sample's positions can be told from one fitted from all of them.
    ///
    /// Positions 1..=6 at a floor of 3, in a ten-base window: the first two and the last two
    /// hold fewer than three covered positions and are refused; the middle two hold all six.
    #[test]
    fn the_histogram_counts_the_windows_the_floor_silenced() {
        let stream: Vec<(u64, u8, u32)> = [1u64, 2, 3, 100, 200, 300]
            .into_iter()
            .map(|p| (p, b'G', 9u32))
            .collect();
        let (windows, histogram) = run_sliding(floored_config(3), 0, &stream);
        assert_eq!(windows.len(), 6);
        assert_eq!(
            histogram.windows_folded + histogram.windows_under_the_floor,
            windows.len() as u64,
            "every finalised window is either folded or counted as silenced",
        );
        assert_eq!(
            histogram.windows_folded, 3,
            "only the run of three clears a floor of 3"
        );
        assert_eq!(histogram.windows_under_the_floor, 3);
    }

    // -- production's sliding-window tests, transcribed --------------------------------

    /// A ramp `depth = position` over positions 1..=10, all `G`, makes every centred window's
    /// mean hand-computable — on a symmetric window over a linear ramp the interior means equal
    /// the centre's own depth, and the edges are truncated.
    ///
    /// The mean-depth half discriminates; the GC half cannot fail, because every base in the
    /// stream is `G`. Mixed-case and mixed-GC windows are covered by
    /// `lowercase_reference_bases_are_read_the_same_as_uppercase` below.
    #[test]
    fn sliding_ramp_means_match_hand_computed() {
        let stream: Vec<(u64, u8, u32)> = (1..=10u64).map(|p| (p, b'G', p as u32)).collect();
        let (windows, histogram) = run_sliding(sliding_config(), 0, &stream);

        // One window per covered position, in position order.
        let positions: Vec<u64> = windows.iter().map(|(at, _)| at.position.get()).collect();
        assert_eq!(positions, (1..=10).collect::<Vec<_>>());

        // Hand-computed centred means over [p-2, p+2] ∩ [1,10]:
        let expected = [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0];
        for ((at, window), &want) in windows.iter().zip(expected.iter()) {
            assert!(
                (window.mean_depth - want).abs() < 1e-6,
                "position {} mean {} != {want}",
                at.position.get(),
                window.mean_depth,
            );
            assert!(
                (window.gc_fraction - 1.0).abs() < 1e-6,
                "all-G stream has GC 1.0"
            );
        }
        assert_eq!(histogram.windows_folded, 10);
    }

    /// A uniform stream: every window's mean equals the constant depth and GC, so every sample
    /// lands in a single histogram cell — **and the cell is named**, because "one non-empty cell
    /// holding twenty" is satisfied by any positive width and any in-range index formula.
    #[test]
    fn sliding_uniform_all_one_cell() {
        let stream: Vec<(u64, u8, u32)> = (1..=20u64).map(|p| (p, b'C', 12)).collect();
        let (windows, histogram) = run_sliding(sliding_config(), 3, &stream);
        assert_eq!(windows.len(), 20);
        for (at, window) in &windows {
            assert!((window.mean_depth - 12.0).abs() < 1e-6);
            assert!((window.gc_fraction - 1.0).abs() < 1e-6);
            assert_eq!(at.contig, ContigId(3));
        }
        // Median 12, range 10, 40 bins → width 3.0; depth 12 → column 4; GC 1.0 → the last of
        // two GC bins, 1; each row holds 41 cells, so the cell is 1 × 41 + 4.
        assert_eq!(histogram.depth_bin_width, 3.0);
        assert_eq!(histogram.counts[1 * 41 + 4], 20);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 20);
    }

    /// An `N` reference base is not a centre and contributes to no window's sums. Positions
    /// 1..=5 all at depth 10; position 3 is `N`. The windows are over the four covered
    /// positions {1,2,4,5}; position 3 emits nothing.
    ///
    /// This fixture cannot check the other half of production's claim — that an `N` still
    /// advances the finalisation frontier — because every covered depth is 10, so deferring a
    /// finalisation changes neither a value nor the order.
    /// `a_centre_closes_when_the_frontier_reaches_its_right_edge` is what pins the frontier.
    #[test]
    fn sliding_n_position_is_excluded_from_every_window() {
        let stream = [
            (1u64, b'G', 10u32),
            (2, b'G', 10),
            (3, b'N', 999), // huge depth, must not affect any mean
            (4, b'G', 10),
            (5, b'G', 10),
        ];
        let (windows, histogram) = run_sliding(sliding_config(), 0, &stream);
        let positions: Vec<u64> = windows.iter().map(|(at, _)| at.position.get()).collect();
        assert_eq!(
            positions,
            vec![1, 2, 4, 5],
            "position 3 (N) emits no window"
        );
        for (_, window) in &windows {
            assert!(
                (window.mean_depth - 10.0).abs() < 1e-6,
                "the N position's depth 999 must not enter any mean (got {})",
                window.mean_depth,
            );
        }
        assert_eq!(histogram.windows_folded, 4);
    }

    /// Windows never span contigs: a contig change finalises the previous contig's tail
    /// (truncated) and starts the next fresh. Two contigs, each positions 1..=3 at a distinct
    /// uniform depth.
    #[test]
    fn sliding_window_does_not_span_contigs() {
        let stream = [
            (1u64, b'G', 10u32),
            (2, b'G', 10),
            (3, b'G', 10),
            // contig 1: a different depth
            (1, b'G', 20),
            (2, b'G', 20),
            (3, b'G', 20),
        ];
        let mut accumulator = WindowCoverageAccumulator::new(sliding_config());
        let mut out = Vec::new();
        for (i, &(position, base, depth)) in stream.iter().enumerate() {
            let contig = ContigId(if i < 3 { 0 } else { 1 });
            accumulator.observe(contig, Position(position), base, depth);
            while let Some(window) = accumulator.pop_ready() {
                out.push(window);
            }
        }
        let (tail, _histogram) = accumulator.finish();
        out.extend(tail);
        // Six windows, three per contig; each carries only its own contig's depth, with no
        // averaging across the boundary.
        assert_eq!(out.len(), 6);
        for (at, window) in &out {
            let expected = if at.contig == ContigId(0) { 10.0 } else { 20.0 };
            assert!(
                (window.mean_depth - expected).abs() < 1e-6,
                "contig {:?} position {} mean {} != {expected}",
                at.contig,
                at.position.get(),
                window.mean_depth,
            );
        }
    }

    /// The accumulator is deterministic: the same stream yields identical windows and
    /// histogram.
    ///
    /// Kept for parity with production, and it cannot fail for any defect this module can have
    /// today — the fold holds no map iteration order, no clock and no randomness. It would
    /// start discriminating if one were introduced.
    #[test]
    fn sliding_is_deterministic() {
        let stream: Vec<(u64, u8, u32)> =
            (1..=15u64).map(|p| (p, b'G', (p % 4) as u32 + 1)).collect();
        let first = run_sliding(sliding_config(), 0, &stream);
        let second = run_sliding(sliding_config(), 0, &stream);
        assert_eq!(first.0, second.0);
        assert_eq!(first.1, second.1);
    }

    /// An odd `window_bp` truncates the half: `window_bp = 5 → half = 2`, so the centred window
    /// is `[p-2, p+2]`, the same as `window_bp = 4` here. Pins the documented
    /// `half = window_bp / 2` convention against silent drift.
    #[test]
    fn sliding_odd_window_uses_floor_half() {
        let config = WindowCoverageConfig {
            window_bp: 5,
            ..sliding_config()
        };
        let stream: Vec<(u64, u8, u32)> = (1..=6u64).map(|p| (p, b'G', p as u32)).collect();
        let (windows, _) = run_sliding(config, 0, &stream);
        // half = 2 → position 3's window is [1,5] = {1,2,3,4,5}, mean 3.0.
        let (_, third) = windows
            .iter()
            .find(|(at, _)| at.position == Position(3))
            .unwrap();
        assert!((third.mean_depth - 3.0).abs() < 1e-6, "half=2 window");
    }

    /// An empty stream yields no windows and **no histogram at all**: with nothing measured
    /// there is no median to scale a depth axis by, and an all-zero histogram carrying an
    /// invented width would be a yardstick presented as if it had been fitted.
    #[test]
    fn sliding_empty_stream_yields_no_window_and_no_histogram() {
        let (windows, histogram) = run_sliding_allowing_no_histogram(sliding_config(), 0, &[]);
        assert!(windows.is_empty());
        assert!(histogram.is_none());
    }

    /// `window_bp = 1` → `half = 0`: the window is `[p, p]`, the centre alone — the degenerate
    /// boundary where the two-pointer extents collapse. Each window's mean is exactly that
    /// position's own depth.
    #[test]
    fn sliding_window_bp_one_emits_self_only_window() {
        let config = WindowCoverageConfig {
            window_bp: 1,
            ..sliding_config()
        };
        let stream: Vec<(u64, u8, u32)> = (1..=5u64).map(|p| (p, b'G', p as u32 * 10)).collect();
        let (windows, histogram) = run_sliding(config, 0, &stream);
        assert_eq!(windows.len(), 5);
        for (at, window) in &windows {
            assert!(
                (window.mean_depth - (at.position.get() * 10) as f32).abs() < 1e-6,
                "position {} window must contain only itself, got mean {}",
                at.position.get(),
                window.mean_depth,
            );
        }
        assert_eq!(histogram.windows_folded, 5);
    }

    /// Consecutive covered positions farther apart than `window_bp`: every window is a
    /// singleton and the buffer fully drains between centres — the shrink-left / full-turnover
    /// path. A stale sum carried across the gap would contaminate a singleton mean with a
    /// far-away depth.
    #[test]
    fn sliding_large_gaps_yield_singleton_windows() {
        // window_bp 4 → half 2; positions much more than 4 apart never share a window.
        let stream = [(10u64, b'G', 5u32), (100, b'G', 50), (1000, b'G', 500)];
        let (windows, _) = run_sliding(sliding_config(), 0, &stream);
        assert_eq!(windows.len(), 3);
        assert!((windows[0].1.mean_depth - 5.0).abs() < 1e-6);
        assert!((windows[1].1.mean_depth - 50.0).abs() < 1e-6);
        assert!((windows[2].1.mean_depth - 500.0).abs() < 1e-6);
    }

    /// A single covered position on the contig: `summed_positions == 1`, the window truncated
    /// on both sides — the divisor-is-one path in isolation, which the type's docs claim always
    /// holds.
    #[test]
    fn sliding_single_covered_position_emits_one_self_window() {
        let (windows, histogram) = run_sliding(sliding_config(), 7, &[(42u64, b'C', 9u32)]);
        assert_eq!(windows.len(), 1);
        let (at, window) = &windows[0];
        assert_eq!(at.position, Position(42));
        assert_eq!(at.contig, ContigId(7));
        assert!((window.mean_depth - 9.0).abs() < 1e-6);
        assert!((window.gc_fraction - 1.0).abs() < 1e-6);
        assert_eq!(histogram.windows_folded, 1);
    }

    /// `finish` returns the still-un-finalised tail together with the un-drained ready queue,
    /// so a caller that never calls `pop_ready` still receives every window — the tail cannot
    /// be silently dropped.
    #[test]
    fn sliding_finish_returns_the_undrained_tail() {
        let mut accumulator = WindowCoverageAccumulator::new(sliding_config());
        for p in 1..=5u64 {
            accumulator.observe(ContigId(0), Position(p), b'G', 10);
            // Deliberately do NOT drain `pop_ready` during the stream.
        }
        let (windows, histogram) = accumulator.finish();
        assert_eq!(windows.len(), 5, "every window returned, none dropped");
        assert_eq!(
            histogram
                .fitted()
                .expect("five windows were folded")
                .windows_folded,
            5,
        );
    }

    // -- ng's own: the boundaries the transcribed eleven leave to the differential -----
    //
    // Each of these was written because a mutation to the code it covers left all eleven
    // tests above green. The differential against production used to catch four of them too;
    // its histogram half stopped applying once the depth bin width became a per-sample fit,
    // so these are what survive that.

    /// The bin scheme reaches the histogram unaltered — the configured part echoed, the depth
    /// width fitted.
    ///
    /// Those four fields are how a consumer turns a cell index back into a depth, and nothing
    /// else in this module reads them: a cross-wired or constant echo produces a histogram whose
    /// cells mean something other than what they say, with no panic.
    #[test]
    fn finish_echoes_the_configured_bins_and_the_fitted_depth_width() {
        // Every field a different value: a configuration whose fields coincide cannot tell a
        // cross-wired echo from a right one.
        let config = WindowCoverageConfig {
            window_bp: 6,
            gc_bins: 3,
            depth_bins: 5,
            min_window_positions: 1,
            depth_scale_windows: 1_000,
            depth_range_in_medians: 7.0,
        };
        let (_, histogram) = run_sliding(config, 0, &[(1u64, b'G', 3u32)]);
        assert_eq!(histogram.window_bp, config.window_bp);
        assert_eq!(histogram.gc_bins, config.gc_bins);
        assert_eq!(histogram.depth_bins, config.depth_bins);
        // One window at depth 3 → median 3; 3 × 7 / 5 = 4.2.
        assert!(
            (histogram.depth_bin_width - 4.2).abs() < 1e-12,
            "the width is the median times the range over the bins, got {}",
            histogram.depth_bin_width,
        );
        assert_eq!(
            histogram.counts.len(),
            config.gc_bins as usize * (config.depth_bins as usize + 1),
        );
    }

    /// A window of `window_bp` 1 is its own centre alone, so its mean depth is exactly the
    /// depth fed in — which is what makes a cell index hand-computable.
    fn one_position_per_window_config() -> WindowCoverageConfig {
        WindowCoverageConfig {
            window_bp: 1,
            gc_bins: 2,
            depth_bins: 4,
            min_window_positions: 1,
            depth_scale_windows: 1_000,
            depth_range_in_medians: 1.0,
        }
    }

    /// A mean depth exactly on the top regular bin's edge belongs to the overflow column.
    ///
    /// The overflow column is the fit's own rejection guard — how many windows sat above the
    /// range — so merging it into the last regular bin would hide the thing it measures.
    ///
    /// **A uniform stream at a range of one median puts every window exactly on that edge**:
    /// the median is the common depth, the regular bins span one median, so the top edge *is*
    /// that depth. Nothing here has to be arranged to hit the boundary — the fit lands on it.
    #[test]
    fn mean_depth_on_the_top_bin_edge_lands_in_the_overflow_column() {
        let stream: Vec<(u64, u8, u32)> = (1..=5u64).map(|p| (p, b'A', 4u32)).collect();
        let (_, histogram) = run_sliding(one_position_per_window_config(), 0, &stream);
        // Median 4, range 1, 4 bins → width 1.0, top regular edge 4 × 1.0 = 4.0.
        assert_eq!(histogram.depth_bin_width, 1.0);
        // GC 0 → gc bin 0; depth 4.0 → overflow column 4; cell = 0 * 5 + 4.
        assert_eq!(histogram.counts[4], 5);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 5);
    }

    /// A GC fraction of exactly 1.0 saturates into the last GC bin rather than off the end.
    #[test]
    fn gc_fraction_of_one_saturates_into_the_last_gc_bin() {
        // Depths 1 and 3 on `G` bases far enough apart to be their own windows: median 3
        // (the upper of two), range 1, 4 bins → width 0.75.
        let stream = [(1u64, b'G', 1u32), (100, b'G', 3u32)];
        let (_, histogram) = run_sliding(one_position_per_window_config(), 0, &stream);
        assert_eq!(histogram.depth_bin_width, 0.75);
        // GC 1.0 → floor(1.0 × 2) = 2, clamped to the last GC bin, 1 — the saturation this
        // test exists for. Depth 1 → bin 1 (1 / 0.75 = 1.33), depth 3 → the overflow column 4.
        // Cells 1 × 5 + 1 = 6 and 1 × 5 + 4 = 9.
        assert_eq!(histogram.counts[6], 1);
        assert_eq!(histogram.counts[9], 1);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 2);
    }

    /// A soft-masked reference reads exactly as an unmasked one.
    ///
    /// Lowercase is what a real FASTA supplies over repeats — 227,170 bases of tomato SL4.0 —
    /// so a case-sensitive `G`/`C` test would report those windows as GC-free, and a
    /// case-sensitive `N` test would count masked ground as covered at whatever depth the
    /// walker put there.
    #[test]
    fn lowercase_reference_bases_are_read_the_same_as_uppercase() {
        let stream = [(1u64, b'g', 4u32), (2, b'n', 999), (3, b'c', 4)];
        let (windows, histogram) = run_sliding(sliding_config(), 0, &stream);
        let positions: Vec<u64> = windows.iter().map(|(at, _)| at.position.get()).collect();
        assert_eq!(
            positions,
            vec![1, 3],
            "lowercase n is not a covered position"
        );
        for (_, window) in &windows {
            assert!(
                (window.gc_fraction - 1.0).abs() < 1e-6,
                "lowercase g and c are GC, got {}",
                window.gc_fraction,
            );
            assert!(
                (window.mean_depth - 4.0).abs() < 1e-6,
                "the lowercase n's depth 999 must not enter any mean (got {})",
                window.mean_depth,
            );
        }
        assert_eq!(histogram.windows_folded, 2);
    }

    /// The same position observed twice is two centres, not one.
    ///
    /// The order rule is non-decreasing rather than increasing, so this is a legal stream, and
    /// a consumer keying a cache by centre will see the centre twice.
    #[test]
    fn a_repeated_position_emits_one_window_per_observation() {
        let stream = [(5u64, b'G', 10u32), (5, b'G', 20)];
        let (windows, histogram) = run_sliding(sliding_config(), 0, &stream);
        assert_eq!(windows.len(), 2, "each observation is its own centre");
        for (at, window) in &windows {
            assert_eq!(at.position, Position(5));
            assert!(
                (window.mean_depth - 15.0).abs() < 1e-6,
                "both observations lie in each other's window, got {}",
                window.mean_depth,
            );
        }
        assert_eq!(histogram.windows_folded, 2);
    }

    /// A centre closes as soon as the frontier **reaches** its right edge, not one base later.
    ///
    /// A repeat arriving at exactly `centre + half` is the only input under which that
    /// boundary is observable: close on time and position 1's window holds one observation of
    /// position 3, close a base late and it holds both.
    #[test]
    fn a_centre_closes_when_the_frontier_reaches_its_right_edge() {
        // half = 2, so position 1's window is [0, 3]; the repeat at 3 arrives after it closed.
        let stream = [(1u64, b'G', 10u32), (3, b'G', 20), (3, b'G', 40)];
        let (windows, _) = run_sliding(sliding_config(), 0, &stream);
        assert_eq!(windows.len(), 3);
        assert!(
            (windows[0].1.mean_depth - 15.0).abs() < 1e-6,
            "position 1 closes at the first observation of position 3, got {}",
            windows[0].1.mean_depth,
        );
    }

    /// Every folded window lands in exactly one cell — the invariant `windows_folded` and the
    /// cell counts jointly claim, and which the dropped `n_skipped_tiles` field used to make
    /// visible in production's tiled model.
    #[test]
    fn every_folded_window_lands_in_exactly_one_histogram_cell() {
        let stream: Vec<(u64, u8, u32)> = (1..=25u64)
            .map(|p| {
                (
                    p,
                    [b'A', b'G', b'C', b'T', b'N'][(p % 5) as usize],
                    (p * 3) as u32,
                )
            })
            .collect();
        let (_, histogram) = run_sliding(sliding_config(), 0, &stream);
        assert_eq!(
            u64::from(histogram.counts.iter().sum::<u32>()),
            histogram.windows_folded,
        );
    }

    // -- the two panicking contracts ---------------------------------------------------

    // -- the floor: ng's one departure from production's window -----------------------

    /// A configuration whose floor is `at_least`, over a window wide enough that a run of
    /// consecutive positions all share one window.
    fn floored_config(at_least: u32) -> WindowCoverageConfig {
        WindowCoverageConfig {
            min_window_positions: at_least,
            window_bp: 10,
            ..sliding_config()
        }
    }

    /// A window holding exactly the floor's worth of covered positions still speaks; one
    /// position fewer does not.
    ///
    /// Five consecutive positions at depth 10, `half = 5`, so every one of the five lies in
    /// every other's window and each window holds all five. At a floor of 5 that is enough and
    /// the mean is 10; at a floor of 6 it is one short and every window comes back absent.
    /// The two runs differ only in the floor, so nothing else can explain the difference.
    #[test]
    fn a_window_at_the_floor_speaks_and_one_position_short_of_it_does_not() {
        let stream: Vec<(u64, u8, u32)> = (1..=5u64).map(|p| (p, b'G', 10)).collect();

        let (at_the_floor, _) = run_sliding(floored_config(5), 0, &stream);
        assert_eq!(at_the_floor.len(), 5);
        for (at, window) in &at_the_floor {
            assert!(
                !window.is_absent(),
                "position {} holds all five covered positions, which meets a floor of 5",
                at.position.get(),
            );
            assert!((window.mean_depth - 10.0).abs() < 1e-6);
        }

        let (one_short, _) = run_sliding_allowing_no_histogram(floored_config(6), 0, &stream);
        assert_eq!(one_short.len(), 5, "an absent window is still emitted");
        for (at, window) in &one_short {
            assert!(
                window.is_absent(),
                "position {} holds five positions against a floor of 6, so it must be \
                 absent, got mean {}",
                at.position.get(),
                window.mean_depth,
            );
        }
    }

    /// An absent window is absent in **both** numbers, and compares equal to the absent value
    /// by bit pattern rather than by IEEE comparison — which is what lets a consumer store it,
    /// hand it on, and compare two runs.
    #[test]
    fn an_absent_window_is_absent_in_both_fields_and_compares_by_bits() {
        let (windows, _) =
            run_sliding_allowing_no_histogram(floored_config(6), 0, &[(1u64, b'G', 10u32)]);
        let (_, only) = &windows[0];
        assert!(only.gc_fraction.is_nan(), "GC fraction must be absent too");
        assert!(only.mean_depth.is_nan());
        assert_eq!(
            *only,
            WindowCoverage::absent(),
            "an absent window equals the absent value; a derived PartialEq would not",
        );
    }

    /// **An absent window trains no yardstick.** It is emitted, but it is not folded, so the
    /// histogram carries no cell for it and `windows_folded` does not count it.
    ///
    /// Without this, `NaN / depth_bin_width` casts to 0 and every window too sparse to speak
    /// would pile into the histogram's first cell — the one nearest where a fit looks for the
    /// single-copy peak.
    #[test]
    fn a_window_under_the_floor_is_emitted_but_not_folded() {
        // Two positions 100 apart at a window of 10: each window holds only its own centre, so
        // both are one short of a floor of 2.
        let stream = [(1u64, b'G', 7u32), (101, b'G', 7)];
        let (windows, histogram) = run_sliding_allowing_no_histogram(floored_config(2), 0, &stream);
        assert_eq!(windows.len(), 2, "both centres are still emitted");
        assert!(windows.iter().all(|(_, w)| w.is_absent()));
        assert!(
            histogram.is_none(),
            "no cell may take an absent window, and a sample whose every window was silenced \
             has no depth to scale a histogram by either",
        );

        // The same stream with the floor off folds both, which is what makes the assertions
        // above about the floor rather than about the fixture.
        let (windows, histogram) = run_sliding(floored_config(1), 0, &stream);
        assert!(windows.iter().all(|(_, w)| !w.is_absent()));
        assert_eq!(histogram.windows_folded, 2);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 2);
    }

    /// The floor counts covered positions, not reference span: an `N` inside the window does
    /// not help a window reach it.
    #[test]
    fn an_n_position_does_not_count_toward_the_floor() {
        // Three positions in one window of 10, the middle one `N`: two covered, floor of 3.
        let stream = [(1u64, b'G', 8u32), (2, b'N', 8), (3, b'G', 8)];
        let (windows, _) = run_sliding_allowing_no_histogram(floored_config(3), 0, &stream);
        assert_eq!(windows.len(), 2, "the N emits no window of its own");
        assert!(
            windows.iter().all(|(_, w)| w.is_absent()),
            "two covered positions do not meet a floor of three",
        );

        // The control: the same two covered positions do clear a floor of two, so the
        // assertion above is about the `N` not counting and not about the fixture being thin.
        let (windows, _) = run_sliding(floored_config(2), 0, &stream);
        assert!(windows.iter().all(|(_, w)| !w.is_absent()));
    }

    /// **The floor is decided per window, not once per stream.** One isolated position and a
    /// run of four, at a floor of 3: the lone window holds one position and must be absent, the
    /// four hold four each and must speak — in one stream, from one accumulator, so a verdict
    /// that carried forward from the first window would show.
    #[test]
    fn a_stream_holding_both_a_silent_window_and_a_speaking_one_is_scored_per_window() {
        let stream: Vec<(u64, u8, u32)> = std::iter::once((1u64, b'G', 9u32))
            .chain((100..=103u64).map(|p| (p, b'G', 9u32)))
            .collect();
        let (windows, histogram) = run_sliding(floored_config(3), 0, &stream);
        assert_eq!(windows.len(), 5);
        assert!(windows[0].1.is_absent(), "the lone position is one of one");
        for (at, window) in &windows[1..] {
            assert!(
                !window.is_absent(),
                "position {} sits in a run of four, which clears a floor of 3",
                at.position.get(),
            );
        }
        assert_eq!(
            histogram.windows_folded, 4,
            "only the four that spoke train the yardstick",
        );
    }

    /// **The floor counts distinct coordinates, not records.** One base observed five times is
    /// one base: the window over it is as thin as any single-position window, and must not
    /// clear a floor of five because five records arrived. The order rule is non-decreasing, so
    /// such a stream is legal, and production's own accumulator averages all five.
    #[test]
    fn a_repeated_position_does_not_clear_the_floor_on_its_own() {
        let stream: Vec<(u64, u8, u32)> = (0..5).map(|_| (7u64, b'G', 10u32)).collect();
        let (windows, histogram) = run_sliding_allowing_no_histogram(floored_config(5), 0, &stream);
        assert_eq!(
            windows.len(),
            5,
            "one window per observation, as production emits"
        );
        assert!(windows.iter().all(|(_, w)| w.is_absent()));
        assert!(
            histogram.is_none(),
            "nothing was folded, so nothing was scaled"
        );
    }

    /// **A window never spans contigs, and neither does the floor's count.** Five positions on
    /// each of two contigs: every window holds its own contig's five, so a floor of 6 silences
    /// all ten and a floor of 5 lets all ten speak. Ten positions in total, so a count that
    /// leaked across the boundary would clear the floor of 6.
    #[test]
    fn the_floor_counts_the_positions_of_its_own_contig() {
        fn run(
            floor: u32,
        ) -> (
            Vec<(GenomePosition, WindowCoverage)>,
            Option<CoverageByGcHistogram>,
        ) {
            let mut accumulator = WindowCoverageAccumulator::new(floored_config(floor));
            let mut out = Vec::new();
            for contig in 0..2u32 {
                for position in 1..=5u64 {
                    accumulator.observe(ContigId(contig), Position(position), b'G', 10);
                    while let Some(window) = accumulator.pop_ready() {
                        out.push(window);
                    }
                }
            }
            let (tail, histogram) = accumulator.finish();
            out.extend(tail);
            (out, histogram.fitted())
        }

        let (windows, histogram) = run(6);
        assert_eq!(windows.len(), 10);
        assert!(
            windows.iter().all(|(_, w)| w.is_absent()),
            "five per contig is short of six",
        );
        assert!(histogram.is_none(), "nothing folded, so nothing scaled");

        let (windows, histogram) = run(5);
        assert!(
            windows.iter().all(|(_, w)| !w.is_absent()),
            "five per contig meets five",
        );
        assert_eq!(
            histogram.expect("ten windows were folded").windows_folded,
            10,
        );
    }

    /// **A pair is absent if either number is.** The two are meant to go missing together and
    /// the accumulator only ever produces them that way, but the fields are `pub`, so the
    /// predicate is what holds the line: a real depth over a missing GC fraction is not a usable
    /// measurement.
    #[test]
    fn is_absent_is_true_when_only_one_field_is_missing() {
        let no_depth = WindowCoverage {
            gc_fraction: 0.4,
            mean_depth: f32::NAN,
        };
        let no_gc = WindowCoverage {
            gc_fraction: f32::NAN,
            mean_depth: 12.0,
        };
        assert!(no_depth.is_absent());
        assert!(
            no_gc.is_absent(),
            "a real depth over a missing GC fraction is not a usable pair",
        );
    }

    /// **Ask `is_absent`, not `== absent()`.** Bitwise equality is right for comparing two runs
    /// and wrong for asking whether a value is missing: a `NaN` that arrived through a codec or
    /// an `f64` round trip carries a different payload, so it compares unequal to this module's
    /// own absent value while still being absent.
    #[test]
    fn is_absent_accepts_any_nan_payload_but_equality_does_not() {
        let foreign = WindowCoverage {
            gc_fraction: f32::from_bits(0x7fc0_1234),
            mean_depth: f32::from_bits(0x7fc0_1234),
        };
        assert!(foreign.is_absent());
        assert!(
            foreign != WindowCoverage::absent(),
            "bitwise equality is payload-sensitive, which is why `is_absent` exists",
        );
    }

    #[test]
    #[should_panic(expected = "window_bp must be >= 1")]
    fn new_panics_on_zero_window_bp() {
        let config = WindowCoverageConfig {
            window_bp: 0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    #[test]
    #[should_panic(expected = "gc_bins must be >= 1")]
    fn new_panics_on_zero_gc_bins() {
        let config = WindowCoverageConfig {
            gc_bins: 0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    /// Zero is the floor's only value below the valid range: one is met by every window,
    /// since a window always holds at least its own centre.
    #[test]
    #[should_panic(expected = "min_window_positions must be >= 1")]
    fn new_panics_on_a_zero_floor() {
        let config = WindowCoverageConfig {
            min_window_positions: 0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    /// A floor above what the window is wide enough to hold would silence every window of every
    /// sample and report nothing about it — the failure a run tuning this number could most
    /// easily walk into.
    #[test]
    #[should_panic(expected = "is more positions than a 10-base window can hold")]
    fn new_panics_on_a_floor_no_window_could_reach() {
        // window_bp 10 spans [p-5, p+5] — eleven bases, so eleven positions at most.
        let config = WindowCoverageConfig {
            min_window_positions: 12,
            ..floored_config(1)
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    #[test]
    #[should_panic(expected = "depth_bins must be >= 1")]
    fn new_panics_on_zero_depth_bins() {
        let config = WindowCoverageConfig {
            depth_bins: 0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    #[test]
    #[should_panic(expected = "depth_scale_windows must be >= 1")]
    fn new_panics_on_a_zero_depth_scale_sample() {
        let config = WindowCoverageConfig {
            depth_scale_windows: 0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    #[test]
    #[should_panic(expected = "depth_range_in_medians must be finite and > 0")]
    fn new_panics_on_a_zero_depth_range() {
        let config = WindowCoverageConfig {
            depth_range_in_medians: 0.0,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    /// A `NaN` range makes the fitted width `NaN`, which is exactly the value the fit reads as
    /// "this sample cannot be scaled" — so without this guard a misconfigured run would cost
    /// every sample its histogram and say nothing.
    #[test]
    #[should_panic(expected = "depth_range_in_medians must be finite and > 0")]
    fn new_panics_on_a_nan_depth_range() {
        let config = WindowCoverageConfig {
            depth_range_in_medians: f64::NAN,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
    }

    // -- the depth scale: a bin width fitted to the sample ------------------------------

    /// The width is the median of the sample's own window depths, times the range, over the
    /// bins — and the median is the upper of the two middles on an even count.
    ///
    /// Five windows at depths 1 to 5, far enough apart to be their own windows: the median is
    /// 3, and with a range of 7 over 5 bins the width is 4.2. Fed in an order that is not
    /// sorted, so a "take the last" or "take the first" would give 5 or 4 instead.
    #[test]
    fn the_depth_width_is_fitted_from_the_median_of_the_sample_own_window_depths() {
        let config = WindowCoverageConfig {
            depth_bins: 5,
            depth_range_in_medians: 7.0,
            ..one_position_per_window_config()
        };
        let stream = [
            (1u64, b'G', 4u32),
            (100, b'G', 1),
            (200, b'G', 5),
            (300, b'G', 3),
            (400, b'G', 2),
        ];
        let (_, histogram) = run_sliding(config, 0, &stream);
        assert!(
            (histogram.depth_bin_width - 4.2).abs() < 1e-12,
            "median 3 × range 7 / 5 bins = 4.2, got {}",
            histogram.depth_bin_width,
        );
        assert_eq!(histogram.windows_folded, 5);
    }

    /// **A window held back for the fit lands in the cell it would have landed in had the width
    /// already been known.** The same stream is run twice: once with the width fitted from the
    /// first window, so every later window is folded as it is finalised, and once with a scale
    /// sample larger than the stream, so every window is held back and folded at `finish`.
    ///
    /// The depth is uniform, so the median — and therefore the width — is the same either way;
    /// the bases are not, so the histogram has more than one cell in it and a mis-binning of the
    /// held-back windows would show. **What a uniform depth cannot see** is a held-back fold
    /// that uses some other width, because that defect is present in both arms and cancels;
    /// `a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately` is what covers it.
    #[test]
    fn a_stream_gives_one_histogram_whether_its_windows_are_folded_early_or_late() {
        let bases = [b'G', b'A', b'C', b'T', b'G', b'G', b'A', b'T'];
        let stream: Vec<(u64, u8, u32)> = (1..=8u64)
            .map(|p| (p, bases[(p - 1) as usize], 12u32))
            .collect();

        let folded_as_they_come = WindowCoverageConfig {
            depth_scale_windows: 1,
            ..sliding_config()
        };
        let all_held_back = WindowCoverageConfig {
            depth_scale_windows: 1_000,
            ..sliding_config()
        };
        let (early_windows, early) = run_sliding(folded_as_they_come, 0, &stream);
        let (late_windows, late) = run_sliding(all_held_back, 0, &stream);

        assert_eq!(
            early_windows, late_windows,
            "the emitted pairs never depended on the fit"
        );
        assert_eq!(
            early, late,
            "and neither does the histogram they are folded into"
        );
        assert!(
            early.counts.iter().filter(|&&c| c > 0).count() > 1,
            "the fixture must fill more than one cell, or this test cannot see a mis-binning",
        );
    }

    /// A sample that finalises fewer windows than the scale sample asks for fits its width at
    /// `finish`, from what it has — the common case for a run over a few regions, and for the
    /// last sample of any cohort.
    #[test]
    fn a_stream_shorter_than_the_scale_sample_fits_its_width_at_the_end() {
        let config = WindowCoverageConfig {
            depth_bins: 5,
            depth_range_in_medians: 7.0,
            depth_scale_windows: 1_000_000,
            ..one_position_per_window_config()
        };
        let stream = [(1u64, b'G', 3u32), (100, b'G', 3), (200, b'G', 3)];
        let (windows, histogram) = run_sliding(config, 0, &stream);
        assert_eq!(windows.len(), 3);
        assert!(
            (histogram.depth_bin_width - 4.2).abs() < 1e-12,
            "median 3 × range 7 / 5 bins = 4.2 from three windows, got {}",
            histogram.depth_bin_width,
        );
        assert_eq!(
            histogram.windows_folded, 3,
            "all three were folded at the end"
        );
    }

    /// The median is the upper of the two middles on an even count, the middle on an odd one,
    /// and `0.0` when there is nothing to take a median of — the value that makes the width
    /// non-positive and so leaves the sample without a histogram.
    #[test]
    fn median_depth_returns_the_upper_middle_on_an_even_count() {
        fn means(pairs: &[(f64, f64)]) -> Vec<WindowMeans> {
            pairs
                .iter()
                .map(|&(gc_fraction, mean_depth)| WindowMeans {
                    gc_fraction,
                    mean_depth,
                })
                .collect()
        }
        assert_eq!(median_depth(&mut []), 0.0);
        assert_eq!(median_depth(&mut means(&[(0.5, 7.0)])), 7.0);
        // Even, and unsorted, so taking the first or the last would give a different answer.
        assert_eq!(median_depth(&mut means(&[(0.1, 9.0), (0.2, 1.0)])), 9.0);
        assert_eq!(
            median_depth(&mut means(&[
                (0.1, 4.0),
                (0.2, 1.0),
                (0.3, 8.0),
                (0.4, 2.0)
            ])),
            4.0,
        );
        // Odd.
        assert_eq!(
            median_depth(&mut means(&[(0.1, 5.0), (0.2, 1.0), (0.3, 3.0)])),
            3.0
        );
        // All equal: the answer cannot depend on which of the ties the selection lands on.
        assert_eq!(
            median_depth(&mut means(&[(0.1, 6.0), (0.9, 6.0), (0.3, 6.0)])),
            6.0
        );
        // The GC fraction is never the key, even when it orders the pairs the other way.
        assert_eq!(median_depth(&mut means(&[(9.0, 1.0), (1.0, 9.0)])), 9.0);
    }

    /// **A sample whose held-back windows have no positive median is carried absent, and stays
    /// absent.** Three of the first four windows read depth zero, so the median is zero and no
    /// width can be fitted — but the fourth reads 80 and the two after it read 30. Refitting
    /// from the remainder would hand the filter a yardstick built from the last two windows of
    /// six, with `windows_folded` reporting two and nothing to say the other four existed.
    #[test]
    fn a_failed_fit_does_not_start_a_second_scale_sample() {
        let config = WindowCoverageConfig {
            depth_scale_windows: 4,
            ..one_position_per_window_config()
        };
        let stream = [
            (1u64, b'G', 0u32),
            (100, b'G', 0),
            (200, b'G', 0),
            (300, b'G', 80),
            (400, b'G', 30),
            (500, b'G', 30),
        ];
        let (windows, histogram) = run_sliding_allowing_no_histogram(config, 0, &stream);
        assert_eq!(windows.len(), 6);
        assert!(
            windows.iter().all(|(_, w)| !w.is_absent()),
            "every window clears a floor of one, so every one of them speaks",
        );
        assert!(
            histogram.is_none(),
            "the fit failed on a median of zero, so this sample has no yardstick — a histogram \
             here would have been fitted from the tail of the sample and would count only it",
        );
    }

    /// The width is fitted from exactly `depth_scale_windows` windows, not one more.
    ///
    /// Two windows at depths 10 and 2 — median 10, the upper of the two — then a third at 4. A
    /// fit one window late would see `[10, 2, 4]` and take 4.
    #[test]
    fn the_scale_sample_is_exactly_as_many_windows_as_it_says() {
        let config = WindowCoverageConfig {
            depth_scale_windows: 2,
            ..one_position_per_window_config()
        };
        let stream = [(1u64, b'G', 10u32), (100, b'G', 2), (200, b'G', 4)];
        let (_, histogram) = run_sliding(config, 0, &stream);
        assert_eq!(
            histogram.depth_bin_width, 2.5,
            "median 10 of the first two × range 1 / 4 bins; a fit one window late gives 1.0",
        );
        assert_eq!(histogram.windows_folded, 3);
    }

    /// **The depth axis is the sample's, not the contig's.** The width is fitted on the first
    /// contig from windows at depth 10; the second contig is four times deeper. A width refitted
    /// at the contig boundary would come out 10.0 instead of 2.5, and one histogram's rows would
    /// be cut on two different axes.
    #[test]
    fn the_fitted_depth_width_survives_a_contig_change() {
        let config = WindowCoverageConfig {
            depth_scale_windows: 2,
            ..one_position_per_window_config()
        };
        let mut accumulator = WindowCoverageAccumulator::new(config);
        for (contig, depth) in [(0u32, 10u32), (1, 40)] {
            for step in 0..2u64 {
                accumulator.observe(ContigId(contig), Position(1 + step * 100), b'G', depth);
                while accumulator.pop_ready().is_some() {}
            }
        }
        let (_, histogram) = accumulator.finish();
        let histogram = histogram.fitted().expect("four windows were folded");
        assert_eq!(
            histogram.depth_bin_width, 2.5,
            "median 10 from the first contig × range 1 / 4 bins",
        );
        assert_eq!(
            histogram.windows_folded, 4,
            "both contigs' windows are folded"
        );
    }

    /// **Held back and folded immediately are the same fold, at depths that differ.** With a
    /// scale sample of two, the first two windows are held back and folded at the fit and the
    /// two after are folded as they finalise. Their depths differ, so a held-back window folded
    /// against any width but the fitted one lands in a different column — which the uniform-depth
    /// fixture of `a_stream_gives_one_histogram_whether_its_windows_are_folded_early_or_late`
    /// cannot see, because there a defect present in both of its arms cancels.
    #[test]
    fn a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately() {
        let config = WindowCoverageConfig {
            depth_scale_windows: 2,
            ..one_position_per_window_config()
        };
        // Held back: two windows at depth 10 → median 10 → width 10 × 1 / 4 = 2.5.
        // Folded live: depth 3 → column 1 (3 / 2.5 = 1.2); depth 40 → the overflow column 4.
        let stream = [
            (1u64, b'G', 10u32),
            (100, b'G', 10),
            (200, b'G', 3),
            (300, b'G', 40),
        ];
        let (_, histogram) = run_sliding(config, 0, &stream);
        assert_eq!(histogram.depth_bin_width, 2.5);
        // GC 1.0 → the last of two GC bins, 1; each row holds depth_bins + 1 = 5 cells. The two
        // held back at depth 10 land at 10 / 2.5 = 4, the overflow column, and the live 40 with
        // them.
        assert_eq!(histogram.counts[1 * 5 + 4], 3);
        assert_eq!(histogram.counts[1 * 5 + 1], 1, "the live window at depth 3");
        assert_eq!(histogram.counts.iter().sum::<u32>(), 4);
        assert_eq!(histogram.windows_folded, 4);
    }

    /// **A sample whose windows all report zero depth gets no histogram.** The median is zero,
    /// so the width would be zero — every window in the overflow column and a yardstick with no
    /// scale on it. The same answer as for a sample that finalised no window at all, and for the
    /// same reason: nothing was measured that a yardstick could be built from.
    #[test]
    fn a_sample_whose_windows_all_report_zero_depth_gets_no_histogram() {
        let stream: Vec<(u64, u8, u32)> = (1..=5u64).map(|p| (p, b'G', 0u32)).collect();
        let (windows, histogram) = run_sliding_allowing_no_histogram(sliding_config(), 0, &stream);
        assert_eq!(
            windows.len(),
            5,
            "the windows are still emitted, and they are not absent"
        );
        assert!(windows.iter().all(|(_, w)| !w.is_absent()));
        assert!(histogram.is_none());
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "out of order")]
    fn observe_panics_in_debug_on_a_position_that_goes_backwards() {
        let mut accumulator = WindowCoverageAccumulator::new(sliding_config());
        accumulator.observe(ContigId(0), Position(20), b'A', 1);
        accumulator.observe(ContigId(0), Position(5), b'A', 1);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "out of order")]
    fn observe_panics_in_debug_on_a_contig_that_goes_backwards() {
        let mut accumulator = WindowCoverageAccumulator::new(sliding_config());
        accumulator.observe(ContigId(1), Position(5), b'A', 1);
        accumulator.observe(ContigId(0), Position(5), b'A', 1);
    }

    /// The absent pair equals itself, and equality is over bit patterns rather than IEEE
    /// comparison. The floor is what produces one, and three tests above assert on pairs it
    /// produced; this is the statement of the convention itself, on a pair built by hand.
    #[test]
    fn an_absent_window_equals_itself_and_signed_zeroes_differ() {
        let absent = WindowCoverage {
            gc_fraction: f32::NAN,
            mean_depth: f32::NAN,
        };
        let also_absent = WindowCoverage {
            gc_fraction: f32::NAN,
            mean_depth: f32::NAN,
        };
        assert_eq!(absent, also_absent, "two absent windows must compare equal");
        assert_ne!(
            WindowCoverage {
                gc_fraction: 0.0,
                mean_depth: 0.0,
            },
            WindowCoverage {
                gc_fraction: -0.0,
                mean_depth: 0.0,
            },
            "equality is over bit patterns, not IEEE comparison",
        );
    }
}
