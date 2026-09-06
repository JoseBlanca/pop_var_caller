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

use super::{CoverageByGcHistogram, WindowCoverage, WindowCoverageConfig};
use crate::ng::types::{ContigId, GenomePosition, Position};

/// A covered position retained in the sliding-window buffer.
#[derive(Debug, Clone, Copy)]
struct CoveredPosition {
    position: u64,
    depth: u32,
    is_gc: bool,
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
/// position without advancing would buffer without bound.
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
    /// Row-major `[gc_bin][depth_bin]` counts, length `config.cell_count()`.
    counts: Vec<u32>,
    /// Windows folded into the histogram so far — one per covered position, minus the ones
    /// the floor silenced.
    windows_folded: u64,
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
            config,
            windows_folded: 0,
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
    /// Positions must arrive in non-decreasing genome order, which is debug-asserted — the
    /// same guard production carries, kept because each sample's records reach its own
    /// accumulator in that sample's own coordinate order. Whether the check should hold in
    /// release too is a question for the step that supplies the caller (plan step C2). An `N`
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
        debug_assert!(
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
    /// **The return shape changes at plan step A3**, where a depth bin width fitted to the
    /// sample means a sample that finalised no window has no histogram at all (spec §3.6).
    #[must_use = "the tail windows and the histogram are the pass's whole output; dropping \
                  them discards every window the accumulator has not already handed over"]
    pub fn finish(mut self) -> (Vec<(GenomePosition, WindowCoverage)>, CoverageByGcHistogram) {
        self.finalise_all();
        let tail: Vec<(GenomePosition, WindowCoverage)> = self.ready.into_iter().collect();
        // Exhaustive destructure: a field added to the configuration must be either carried
        // into the histogram or explicitly ignored here, rather than silently omitted. A2 and
        // A3 both add one.
        let WindowCoverageConfig {
            window_bp,
            gc_bins,
            depth_bin_width,
            depth_bins,
            // Not a bin scheme: the floor decides which windows are folded, and a consumer
            // reading a cell does not need it.
            min_window_positions: _,
        } = self.config;
        let histogram = CoverageByGcHistogram {
            window_bp,
            gc_bins,
            depth_bin_width,
            depth_bins,
            windows_folded: self.windows_folded,
            counts: self.counts,
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
            WindowCoverage::absent()
        } else {
            let divisor = self.summed_positions as f64;
            let gc_fraction = self.sum_gc as f64 / divisor;
            let mean_depth = self.sum_depth as f64 / divisor;
            // Saturating rather than wrapping: a cell takes at most one count per covered
            // position, so it can only reach `u32::MAX` on a reference above about 4.3 Gbp. A
            // saturated cell under-reports; a wrapped one reports a near-empty cell where the
            // single-copy peak is, and the filter's fit anchors on exactly that mode.
            let cell = &mut self.counts[self.config.cell_index(gc_fraction, mean_depth)];
            *cell = cell.saturating_add(1);
            self.windows_folded += 1;
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
            depth_bin_width: 1.0,
            depth_bins: 40,
            // The floor off, so these keep testing production's window: at 1 every window
            // clears it, because a window always holds at least its own centre.
            min_window_positions: 1,
        }
    }

    /// Feed a `(position, reference base, depth)` stream on one contig and collect every
    /// emitted window (drained after each `observe`, plus the tail `finish` returns), together
    /// with the finished histogram.
    fn run_sliding(
        config: WindowCoverageConfig,
        contig: u32,
        stream: &[(u64, u8, u32)],
    ) -> (Vec<(GenomePosition, WindowCoverage)>, CoverageByGcHistogram) {
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
    /// lands in a single histogram cell.
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
        let nonzero: Vec<u32> = histogram
            .counts
            .iter()
            .copied()
            .filter(|&c| c > 0)
            .collect();
        assert_eq!(nonzero, vec![20], "all 20 windows in one cell");
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
        let (tail, _) = accumulator.finish();
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

    /// An empty stream yields no windows and an empty histogram.
    #[test]
    fn sliding_empty_stream_is_empty() {
        let (windows, histogram) = run_sliding(sliding_config(), 0, &[]);
        assert!(windows.is_empty());
        assert_eq!(histogram.windows_folded, 0);
        assert!(histogram.counts.iter().all(|&c| c == 0));
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
        assert_eq!(histogram.windows_folded, 5);
    }

    // -- ng's own: the boundaries the transcribed eleven leave to the differential -----
    //
    // Each of these was written because a mutation to the code it covers left all eleven
    // tests above green. The differential against production catches four of them too, but
    // its histogram half stops applying once step A3 fits the depth bin width per sample,
    // so these are what survive that.

    /// The configuration reaches the histogram unaltered.
    ///
    /// The four scheme fields are how a consumer turns a cell index back into a depth, and
    /// nothing else in this module reads them — a cross-wired or constant echo produces a
    /// histogram whose cells mean something other than what they say, with no panic.
    #[test]
    fn finish_echoes_the_configured_bin_scheme() {
        // Every field a different value, and `depth_bin_width` deliberately not 1.0: a
        // configuration whose fields coincide cannot tell a cross-wired echo from a right one.
        let config = WindowCoverageConfig {
            window_bp: 6,
            gc_bins: 3,
            depth_bin_width: 0.7,
            depth_bins: 5,
            min_window_positions: 1,
        };
        let (_, histogram) = run_sliding(config, 0, &[(1u64, b'G', 3u32)]);
        assert_eq!(histogram.window_bp, config.window_bp);
        assert_eq!(histogram.gc_bins, config.gc_bins);
        assert_eq!(histogram.depth_bin_width, config.depth_bin_width);
        assert_eq!(histogram.depth_bins, config.depth_bins);
        assert_eq!(
            histogram.counts.len(),
            config.gc_bins as usize * (config.depth_bins as usize + 1),
        );
    }

    /// A mean depth exactly on the top regular bin's edge belongs to the overflow column.
    ///
    /// The overflow column is the fit's own rejection guard — how many windows sat above the
    /// range — so merging it into the last regular bin would hide the thing it measures.
    #[test]
    fn mean_depth_on_the_top_bin_edge_lands_in_the_overflow_column() {
        // window_bp 1 → each window is its centre alone, so the mean is that depth.
        // depth_bins 4 at width 1.0 → the top regular edge is 4.0.
        let config = WindowCoverageConfig {
            window_bp: 1,
            gc_bins: 2,
            depth_bin_width: 1.0,
            depth_bins: 4,
            min_window_positions: 1,
        };
        let (_, histogram) = run_sliding(config, 0, &[(1u64, b'A', 4u32)]);
        // GC 0 → gc bin 0; depth 4.0 → overflow column 4; cell = 0 * 5 + 4.
        assert_eq!(histogram.counts[4], 1);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 1);
    }

    /// A GC fraction of exactly 1.0 saturates into the last GC bin rather than off the end.
    #[test]
    fn gc_fraction_of_one_saturates_into_the_last_gc_bin() {
        let config = WindowCoverageConfig {
            window_bp: 1,
            gc_bins: 2,
            depth_bin_width: 1.0,
            depth_bins: 4,
            min_window_positions: 1,
        };
        let (_, histogram) = run_sliding(config, 0, &[(1u64, b'G', 1u32)]);
        // GC 1.0 → floor(1.0 * 2) = 2, clamped to 1; depth 1.0 → bin 1; cell = 1 * 5 + 1.
        assert_eq!(histogram.counts[6], 1);
        assert_eq!(histogram.counts.iter().sum::<u32>(), 1);
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
            window_bp: 10,
            gc_bins: 2,
            depth_bin_width: 1.0,
            depth_bins: 40,
            min_window_positions: at_least,
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

        let (one_short, _) = run_sliding(floored_config(6), 0, &stream);
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
        let (windows, _) = run_sliding(floored_config(6), 0, &[(1u64, b'G', 10u32)]);
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
        let (windows, histogram) = run_sliding(floored_config(2), 0, &stream);
        assert_eq!(windows.len(), 2, "both centres are still emitted");
        assert!(windows.iter().all(|(_, w)| w.is_absent()));
        assert_eq!(histogram.windows_folded, 0);
        assert_eq!(
            histogram.counts.iter().sum::<u32>(),
            0,
            "no cell may take an absent window",
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
        let (windows, _) = run_sliding(floored_config(3), 0, &stream);
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
        let (windows, histogram) = run_sliding(floored_config(5), 0, &stream);
        assert_eq!(
            windows.len(),
            5,
            "one window per observation, as production emits"
        );
        assert!(windows.iter().all(|(_, w)| w.is_absent()));
        assert_eq!(histogram.windows_folded, 0);
    }

    /// **A window never spans contigs, and neither does the floor's count.** Five positions on
    /// each of two contigs: every window holds its own contig's five, so a floor of 6 silences
    /// all ten and a floor of 5 lets all ten speak. Ten positions in total, so a count that
    /// leaked across the boundary would clear the floor of 6.
    #[test]
    fn the_floor_counts_the_positions_of_its_own_contig() {
        fn run(floor: u32) -> (Vec<(GenomePosition, WindowCoverage)>, CoverageByGcHistogram) {
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
            (out, histogram)
        }

        let (windows, histogram) = run(6);
        assert_eq!(windows.len(), 10);
        assert!(
            windows.iter().all(|(_, w)| w.is_absent()),
            "five per contig is short of six",
        );
        assert_eq!(histogram.windows_folded, 0);

        let (windows, histogram) = run(5);
        assert!(
            windows.iter().all(|(_, w)| !w.is_absent()),
            "five per contig meets five",
        );
        assert_eq!(histogram.windows_folded, 10);
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

    /// A `NaN` width is the case worth its own test: without the guard, `mean / NaN` is `NaN`,
    /// `as usize` yields 0, and every window silently lands in depth bin 0.
    #[test]
    #[should_panic(expected = "depth_bin_width must be finite and > 0")]
    fn new_panics_on_non_finite_depth_bin_width() {
        let config = WindowCoverageConfig {
            depth_bin_width: f64::NAN,
            ..sliding_config()
        };
        let _ = WindowCoverageAccumulator::new(config);
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
    /// comparison. Nothing in this module produces the absent pair yet — plan step A2's floor
    /// is what starts emitting it — so this is what keeps the convention from being an
    /// undefended sentence in a doc comment until then.
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
