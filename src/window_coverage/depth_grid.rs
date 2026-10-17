//! The fixed depth grid a sample's windows are folded on while the pass runs, and the re-cut that
//! turns it into the histogram the coverage fit reads (spec `window_coverage.md` §3.4).
//!
//! **Why two axes.** The fit wants a linear depth axis scaled to the sample's median window depth
//! — its mode search, its median and its scatter all assume equal-width bins — and that median is
//! only known once every window has been seen. The grid has to accept a window before then, from
//! a sample of any depth, so its bins are equal steps of *relative* depth: the same 2.33% of the
//! depth at one read a position as at a thousand. When the pass ends, [`DepthGrid::median_depth`]
//! reads the median off the grid and [`DepthGrid::recut`] spreads each grid bin's count over the
//! linear bins it overlaps.

use super::{DEPTH_GRID_BINS, DEPTH_GRID_HIGHEST_READS, DEPTH_GRID_LOWEST_READS};
use crate::float;

/// One sample's windows, counted on the fixed relative-depth grid: a row-major
/// `[gc_bin][grid_bin]` matrix with [`DEPTH_GRID_BINS`] regular bins and one overflow bin a row.
///
/// Grid bin `k` holds mean depths in `[lowest · step^k, lowest · step^(k+1))`, where `step` is the
/// ratio between consecutive edges; a depth below the lowest edge — which no window reports, a
/// covered position having at least one read — is counted in bin 0, and one at or above the
/// highest edge in the overflow bin.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct DepthGrid {
    gc_bins: usize,
    /// Row-major counts, `gc_bins * (DEPTH_GRID_BINS + 1)` long. A cell saturates at `u32::MAX`
    /// for the reason the histogram's cells do ([`super::CoverageByGcHistogram::counts`]).
    counts: Vec<u32>,
}

/// The natural log of the ratio between consecutive grid edges.
fn log_step() -> f64 {
    float::ln(DEPTH_GRID_HIGHEST_READS / DEPTH_GRID_LOWEST_READS) / f64::from(DEPTH_GRID_BINS)
}

/// The lower edge of grid bin `k`, in reads a position. `k == DEPTH_GRID_BINS` is the highest
/// edge, where the overflow bin begins.
fn lower_edge(k: usize) -> f64 {
    DEPTH_GRID_LOWEST_READS * float::exp(k as f64 * log_step())
}

impl DepthGrid {
    pub(super) fn new(gc_bins: u32) -> Self {
        let gc_bins = gc_bins as usize;
        Self {
            gc_bins,
            counts: vec![0; gc_bins * Self::columns()],
        }
    }

    /// How many cells the grid holds, for the memory test that prices it.
    #[cfg(test)]
    pub(super) fn cell_count(&self) -> usize {
        self.counts.len()
    }

    /// Cells per GC row: the regular grid bins plus the overflow bin.
    fn columns() -> usize {
        DEPTH_GRID_BINS as usize + 1
    }

    /// Count one window, already placed in its GC bin.
    pub(super) fn fold(&mut self, gc_bin: usize, mean_depth: f64) {
        let grid_bins = DEPTH_GRID_BINS as usize;
        let grid_bin = if mean_depth >= DEPTH_GRID_HIGHEST_READS {
            grid_bins
        } else if mean_depth <= DEPTH_GRID_LOWEST_READS {
            0
        } else {
            // The `min` guards the float rounding just under the highest edge, which can land a
            // depth the branch above did not catch one past the last regular bin.
            ((float::ln(mean_depth / DEPTH_GRID_LOWEST_READS) / log_step()) as usize)
                .min(grid_bins - 1)
        };
        let cell = &mut self.counts[gc_bin * Self::columns() + grid_bin];
        *cell = cell.saturating_add(1);
    }

    /// The median window depth over every window counted, in reads a position, or `None` if none
    /// was. **The upper of the two middles** on an even count, as the prefix median it replaces
    /// took; inside its grid bin it is placed by interpolating the bin's count evenly over the
    /// bin's log-depth span, so it is not confined to a bin edge.
    ///
    /// **A median in the overflow bin comes back as the highest edge**, which is a lower bound on
    /// it; the linear axis cut from it then puts more than half the windows past its top, and the
    /// fit refuses the sample for range, which is the right report for a sample deeper than the
    /// grid.
    pub(super) fn median_depth(&self) -> Option<f64> {
        let columns = Self::columns();
        let marginal: Vec<u64> = (0..columns)
            .map(|bin| {
                (0..self.gc_bins)
                    .map(|row| u64::from(self.counts[row * columns + bin]))
                    .sum()
            })
            .collect();
        let total: u64 = marginal.iter().sum();
        if total == 0 {
            return None;
        }
        // The 0-based rank of the upper middle, as `select_nth_unstable` at `len / 2` gave.
        let rank = total / 2;
        let mut below = 0u64;
        for (bin, &count) in marginal.iter().enumerate() {
            if rank < below + count {
                if bin == DEPTH_GRID_BINS as usize {
                    return Some(DEPTH_GRID_HIGHEST_READS);
                }
                let fraction_into_bin = ((rank - below) as f64 + 0.5) / count as f64;
                return Some(lower_edge(bin) * float::exp(fraction_into_bin * log_step()));
            }
            below += count;
        }
        unreachable!("the rank is below the total, so some bin holds it")
    }

    /// The same windows on a linear axis of `depth_bins` bins `depth_bin_width` wide, plus an
    /// overflow bin — the row-major layout [`super::CoverageByGcHistogram::counts`] carries.
    ///
    /// **Each grid bin's count is spread evenly over its depth span** and shared out to the linear
    /// bins in proportion to their overlap with it; whatever lies at or past `depth_bins *
    /// depth_bin_width`, and the grid's own overflow, goes to the linear overflow bin. The shares
    /// are fractional, so each row is rounded back to whole counts by largest remainder, ties to
    /// the lower bin: **every row keeps its exact total**, so the histogram still counts every
    /// window once.
    pub(super) fn recut(&self, depth_bin_width: f64, depth_bins: u32) -> Vec<u32> {
        assert!(
            depth_bin_width.is_finite() && depth_bin_width > 0.0,
            "a linear axis needs a finite positive bin width, got {depth_bin_width}",
        );
        let grid_columns = Self::columns();
        let grid_bins = DEPTH_GRID_BINS as usize;
        let linear_bins = depth_bins as usize;
        let linear_columns = linear_bins + 1;
        let top_of_the_axis = linear_bins as f64 * depth_bin_width;
        let mut recut = vec![0u32; self.gc_bins * linear_columns];
        let mut shares = vec![0f64; linear_columns];

        for row in 0..self.gc_bins {
            let grid_row = &self.counts[row * grid_columns..(row + 1) * grid_columns];
            let row_total: u64 = grid_row.iter().map(|&count| u64::from(count)).sum();
            if row_total == 0 {
                continue;
            }
            shares.fill(0.0);
            for (grid_bin, &count) in grid_row[..grid_bins].iter().enumerate() {
                if count == 0 {
                    continue;
                }
                let count = f64::from(count);
                let low = lower_edge(grid_bin);
                let high = lower_edge(grid_bin + 1);
                let span = high - low;
                let mut linear_bin = (low / depth_bin_width) as usize;
                while linear_bin < linear_bins {
                    let bin_low = linear_bin as f64 * depth_bin_width;
                    if bin_low >= high {
                        break;
                    }
                    let overlap = high.min(bin_low + depth_bin_width) - low.max(bin_low);
                    shares[linear_bin] += count * overlap / span;
                    linear_bin += 1;
                }
                if high > top_of_the_axis {
                    shares[linear_bins] += count * (high - low.max(top_of_the_axis)) / span;
                }
            }
            shares[linear_bins] += f64::from(grid_row[grid_bins]);

            let out = &mut recut[row * linear_columns..(row + 1) * linear_columns];
            round_keeping_the_total(&shares, row_total, out);
        }
        recut
    }
}

/// Round fractional `shares` to whole counts in `out` so that they add up to exactly `total`:
/// every share is floored, and the counts still missing go one each to the shares with the largest
/// fractional parts, ties to the lower index — so the rounding is deterministic.
fn round_keeping_the_total(shares: &[f64], total: u64, out: &mut [u32]) {
    let mut floored_sum = 0u64;
    let mut remainders: Vec<(usize, f64)> = Vec::with_capacity(shares.len());
    for (bin, &share) in shares.iter().enumerate() {
        let floor = share.floor().max(0.0);
        out[bin] = floor.min(f64::from(u32::MAX)) as u32;
        floored_sum += u64::from(out[bin]);
        remainders.push((bin, share - floor));
    }
    // Float sums of shares can land a hair either side of the integer total; whatever is owed is
    // at most one per bin.
    let missing = total.saturating_sub(floored_sum) as usize;
    remainders.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    for &(bin, _) in remainders.iter().take(missing) {
        out[bin] = out[bin].saturating_add(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fold_all(depths: &[f64]) -> DepthGrid {
        let mut grid = DepthGrid::new(1);
        for &depth in depths {
            grid.fold(0, depth);
        }
        grid
    }

    #[test]
    fn the_grid_steps_by_a_fixed_share_of_the_depth() {
        let step = float::exp(log_step()) - 1.0;
        assert!((step - 0.0233).abs() < 0.0001, "step {step}");
        assert!((lower_edge(DEPTH_GRID_BINS as usize) - DEPTH_GRID_HIGHEST_READS).abs() < 1e-6);
    }

    /// The median off the grid is within one grid step of the exact one, at three reads a
    /// position and at three hundred.
    #[test]
    fn the_median_is_within_one_grid_step_at_any_depth() {
        for scale in [3.0, 30.0, 300.0] {
            let depths: Vec<f64> = (0..1001)
                .map(|i| scale * (0.5 + i as f64 / 1000.0))
                .collect();
            let exact = depths[500];
            let median = fold_all(&depths).median_depth().unwrap();
            assert!(
                (median / exact - 1.0).abs() < 0.0233,
                "at {scale}: median {median} against {exact}",
            );
        }
    }

    #[test]
    fn an_empty_grid_has_no_median() {
        assert_eq!(DepthGrid::new(3).median_depth(), None);
    }

    #[test]
    fn a_median_past_the_grid_is_its_highest_edge() {
        let median = fold_all(&[20_000.0, 30_000.0, 5.0]).median_depth().unwrap();
        assert_eq!(median, DEPTH_GRID_HIGHEST_READS);
    }

    /// Re-cutting keeps every row's total, and a window lands within one grid step of where a
    /// direct linear fold would have put it.
    #[test]
    fn recutting_keeps_every_window_and_moves_none_further_than_a_grid_step() {
        let depths: Vec<f64> = (0..5000).map(|i| 10.0 + (i % 400) as f64 * 0.1).collect();
        let grid = fold_all(&depths);
        let width = 0.25;
        let recut = grid.recut(width, 400);
        assert_eq!(recut.iter().map(|&c| u64::from(c)).sum::<u64>(), 5000);
        let recut_mean: f64 = recut[..400]
            .iter()
            .enumerate()
            .map(|(bin, &count)| (bin as f64 + 0.5) * width * f64::from(count))
            .sum::<f64>()
            / 5000.0;
        let exact_mean = depths.iter().sum::<f64>() / 5000.0;
        assert!(
            (recut_mean / exact_mean - 1.0).abs() < 0.0233,
            "{recut_mean} against {exact_mean}",
        );
    }

    #[test]
    fn what_lies_past_the_linear_axis_goes_to_its_overflow_bin() {
        let grid = fold_all(&[6.0, 6.0, 50.0, 20_000.0]);
        // Top of the axis at 4 × 2.5 = 10 reads: the two at 50 and 20,000 are over it.
        let recut = grid.recut(2.5, 4);
        assert_eq!(recut, vec![0, 0, 2, 0, 2]);
    }

    #[test]
    fn rounding_keeps_the_total_and_breaks_ties_to_the_lower_bin() {
        let mut out = vec![0u32; 4];
        round_keeping_the_total(&[0.5, 0.5, 0.5, 0.5], 2, &mut out);
        assert_eq!(out, vec![1, 1, 0, 0]);
    }
}
