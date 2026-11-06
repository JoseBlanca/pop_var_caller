//! **Each sample's walk parameters, fitted on the run's first windows and then frozen**
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §3.5).
//!
//! Three numbers are the sample's own. The share of its genome in autozygous runs, `F`, is its
//! inbreeding coefficient from the parameters file, fitted over the whole run, and the walk keeps it:
//! re-estimated from the walk's own autozygous probabilities, it ran away to 0.83 on average in an
//! outbred cohort whose samples were all at 2×, where one soft template explains both copies about
//! as well as two, and the pass called heterozygotes homozygous (§3.5). The mean length of a
//! stretch one copy copies from one template, `S`, and the mean autozygous run, `R`, are picked from
//! small grids by the likelihood of the sample's evidence. The grid costs twelve walks a window, so
//! it runs only until the windows fitted hold [`FIT_RECORDS`] of the sample's informative records;
//! from then on `S` and `R` are fixed for the rest of the run. Until then, each window adds its walks'
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
    /// `F`, the long-run autozygous share: the parameters file's, clamped to
    /// [`AUTOZYGOUS_SHARE_BOUNDS`].
    autozygous_share: f64,
    /// Each grid point's log-likelihood, summed over the windows fitted so far; segment-major.
    summed_log_likelihood: [f64; GRID_POINTS],
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
            informative_records_fitted: 0,
            chosen: 0,
            frozen: false,
        }
    }

    /// Whether the fitting windows have held [`FIT_RECORDS`] of the sample's informative records,
    /// so the parameters no longer move.
    #[must_use]
    pub fn is_frozen(&self) -> bool {
        self.frozen
    }

    /// The parameters the sample is walked with now, with the pass's `mismatch` and
    /// `background_landing_share` (`WalkParameters`).
    #[must_use]
    pub fn parameters(
        &self,
        mismatch: f64,
        background_landing_share: Option<f64>,
    ) -> WalkParameters {
        WalkParameters {
            segment_bases: SEGMENT_GRID[self.chosen / RUN_GRID.len()],
            run_bases: RUN_GRID[self.chosen % RUN_GRID.len()],
            autozygous_share: self.autozygous_share,
            mismatch,
            background_landing_share,
        }
    }

    /// **The parameters to walk this window with**, fitting on it first if the sample is not yet
    /// frozen: the twelve grid walks at the sample's `F`, their log-likelihoods added to the running
    /// totals, and the best total chosen.
    /// `informative_records` is how many of the window's sites are informative records of the
    /// sample. `template_weights` is the walk's, as [`walk`] takes it.
    pub fn fit_window(
        &mut self,
        sites: &[WalkSite<'_>],
        templates: usize,
        template_weights: Option<&[f64]>,
        informative_records: usize,
        mismatch: f64,
        background_landing_share: Option<f64>,
    ) -> WalkParameters {
        if self.frozen || sites.is_empty() {
            return self.parameters(mismatch, background_landing_share);
        }
        for (point, summed) in self.summed_log_likelihood.iter_mut().enumerate() {
            let result = walk(
                sites,
                templates,
                template_weights,
                WalkParameters {
                    segment_bases: SEGMENT_GRID[point / RUN_GRID.len()],
                    run_bases: RUN_GRID[point % RUN_GRID.len()],
                    autozygous_share: self.autozygous_share,
                    mismatch,
                    background_landing_share,
                },
            );
            *summed += result.log_likelihood;
        }
        // The largest total; on a tie, the earlier grid point.
        self.chosen = (0..GRID_POINTS).fold(0, |best, point| {
            if self.summed_log_likelihood[point] > self.summed_log_likelihood[best] {
                point
            } else {
                best
            }
        });
        self.informative_records_fitted += informative_records;
        self.frozen = self.informative_records_fitted >= FIT_RECORDS;
        self.parameters(mismatch, background_landing_share)
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
