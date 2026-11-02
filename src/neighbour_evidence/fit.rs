//! **Each sample's walk parameters, fitted on the run's first windows and then frozen**
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §3.5).
//!
//! Three numbers are the sample's own. The share of its genome in autozygous runs, `F`, starts at
//! its inbreeding coefficient from the parameters file and is re-estimated as the mean probability,
//! over the sites walked, that the sample is autozygous there. The mean length of a stretch one copy
//! copies from one template, `S`, and the mean autozygous run, `R`, are picked from small grids by
//! the likelihood of the sample's evidence. The grid costs twelve walks a window, so it runs only
//! until the windows fitted hold [`FIT_RECORDS`] of the sample's informative records; from then on
//! all three are fixed for the rest of the run. Until then, each window adds its walks'
//! log-likelihoods to a running total a grid point, and the point with the largest total is the
//! one walked with — so the second fitting window refines the first window's choice rather than
//! replacing it.

use super::model::{WalkParameters, WalkSite, walk};

/// The grid of `S`, the mean stretch one copy copies from one template, in bases (§3.5).
pub const SEGMENT_GRID: [f64; 2] = [30_000.0, 300_000.0];

/// The grid of `R`, the mean autozygous run, in bases (§3.5).
pub const RUN_GRID: [f64; 6] = [
    10_000.0,
    30_000.0,
    100_000.0,
    300_000.0,
    1_000_000.0,
    3_000_000.0,
];

/// How many of a sample's informative records the fitting windows must hold before its parameters
/// are frozen (§3.5).
pub const FIT_RECORDS: usize = 5_000;

/// The bounds `F` is clamped to (§3.5): a walk needs both states possible.
pub const AUTOZYGOUS_SHARE_BOUNDS: (f64, f64) = (0.01, 0.99);

const GRID_POINTS: usize = SEGMENT_GRID.len() * RUN_GRID.len();

/// **One sample's walk parameters**, fitting or frozen.
#[derive(Clone, Debug)]
pub struct SampleFit {
    /// `F`, the long-run autozygous share, clamped to [`AUTOZYGOUS_SHARE_BOUNDS`].
    autozygous_share: f64,
    /// Each grid point's log-likelihood, summed over the windows fitted so far; segment-major.
    summed_log_likelihood: [f64; GRID_POINTS],
    /// The autozygous probabilities of every site walked at the chosen grid point in the fitting
    /// windows, summed, and how many sites that is — `F`'s re-estimate is their ratio.
    autozygous_probability_sum: f64,
    autozygous_sites: usize,
    /// The sample's informative records in the windows fitted so far.
    informative_records_fitted: usize,
    /// The grid point walked with: the one with the largest summed log-likelihood.
    chosen: usize,
    frozen: bool,
}

impl SampleFit {
    /// A sample not yet fitted, starting from its inbreeding coefficient in the parameters file.
    #[must_use]
    pub fn new(inbreeding_coefficient: f64) -> Self {
        Self {
            autozygous_share: clamp_share(inbreeding_coefficient),
            summed_log_likelihood: [0.0; GRID_POINTS],
            autozygous_probability_sum: 0.0,
            autozygous_sites: 0,
            informative_records_fitted: 0,
            chosen: 0,
            frozen: false,
        }
    }

    /// A sample frozen at grid point `chosen` (segment-major) and share `inbreeding_coefficient`,
    /// for tests that compare a walk with the prototype at fixed parameters.
    #[cfg(test)]
    pub(crate) fn frozen_at(inbreeding_coefficient: f64, chosen: usize) -> Self {
        Self {
            autozygous_share: clamp_share(inbreeding_coefficient),
            summed_log_likelihood: [0.0; GRID_POINTS],
            autozygous_probability_sum: 0.0,
            autozygous_sites: 0,
            informative_records_fitted: 0,
            chosen,
            frozen: true,
        }
    }

    /// Whether the fitting windows have held [`FIT_RECORDS`] of the sample's informative records,
    /// so the parameters no longer move.
    #[must_use]
    pub fn is_frozen(&self) -> bool {
        self.frozen
    }

    /// The parameters the sample is walked with now.
    #[must_use]
    pub fn parameters(&self, mismatch: f64) -> WalkParameters {
        WalkParameters {
            segment_bases: SEGMENT_GRID[self.chosen / RUN_GRID.len()],
            run_bases: RUN_GRID[self.chosen % RUN_GRID.len()],
            autozygous_share: self.autozygous_share,
            mismatch,
        }
    }

    /// **The parameters to walk this window with**, fitting on it first if the sample is not yet
    /// frozen: the twelve grid walks at the current `F`, their log-likelihoods added to the running
    /// totals, the best total chosen, and `F` re-estimated from that grid point's walk.
    /// `informative_records` is how many of the window's sites are informative records of the
    /// sample.
    pub fn fit_window(
        &mut self,
        sites: &[WalkSite<'_>],
        templates: usize,
        informative_records: usize,
        mismatch: f64,
    ) -> WalkParameters {
        if self.frozen || sites.is_empty() {
            return self.parameters(mismatch);
        }
        let mut autozygous_sums = [0.0; GRID_POINTS];
        for (point, (summed, autozygous_sum)) in self
            .summed_log_likelihood
            .iter_mut()
            .zip(&mut autozygous_sums)
            .enumerate()
        {
            let result = walk(
                sites,
                templates,
                WalkParameters {
                    segment_bases: SEGMENT_GRID[point / RUN_GRID.len()],
                    run_bases: RUN_GRID[point % RUN_GRID.len()],
                    autozygous_share: self.autozygous_share,
                    mismatch,
                },
            );
            *summed += result.log_likelihood;
            *autozygous_sum = result.autozygous.iter().sum();
        }
        // The largest total; on a tie, the earlier grid point.
        self.chosen = (0..GRID_POINTS).fold(0, |best, point| {
            if self.summed_log_likelihood[point] > self.summed_log_likelihood[best] {
                point
            } else {
                best
            }
        });
        self.autozygous_probability_sum += autozygous_sums[self.chosen];
        self.autozygous_sites += sites.len();
        self.autozygous_share =
            clamp_share(self.autozygous_probability_sum / self.autozygous_sites as f64);
        self.informative_records_fitted += informative_records;
        self.frozen = self.informative_records_fitted >= FIT_RECORDS;
        self.parameters(mismatch)
    }
}

fn clamp_share(share: f64) -> f64 {
    let (lowest, highest) = AUTOZYGOUS_SHARE_BOUNDS;
    if share.is_nan() {
        lowest
    } else {
        share.clamp(lowest, highest)
    }
}

#[cfg(test)]
mod tests;
