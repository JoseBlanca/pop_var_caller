//! **The hidden-duplication filter's statistics, copied from production into ng.**
//!
//! At each locus the filter asks which of two stories better explains what every sample
//! showed: a real variant that some of them carry, or two gene copies the reference
//! collapsed into one, piling their reads onto the same position. Four quantities answer
//! that, and this module is where they will all live. **One has landed so far**: the
//! per-sample coverage model, which says what read depth *one copy* would produce in a
//! sample at a window of a given GC content, so that dividing an observed depth by it
//! gives a copy number — one for a single copy, about two for a collapsed pair. The
//! other three arrive with the plan's next two steps: how much better one story fits
//! than the other at a locus, how common hidden duplications are in this run, and where
//! to cut. Nothing here reads a file, writes a record or knows what a VCF is; the run's
//! wiring will live in `crate::ng::run::paralog_filter` (plan step B1) and has not
//! landed.
//!
//! **Copied, not re-derived.** The model is production's
//! (`doc/devel/specs/hidden_paralog_filter.md`, reformulated for one sample by
//! `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`). Two of the three
//! files beside this one are production's source, line for line and byte for byte:
//! [`coverage_model`], 1,158 lines past its module header, and [`model_params`], 236
//! lines. The third, `copy_fidelity`, is ng's own and is the test that asserts the other
//! two — so an edit to *either* tree fails the build instead of drifting quietly. What
//! ng changes about the filter is where its numbers come from and when they are known,
//! not what it computes: `doc/devel/ng/spec/hidden_paralog_filter.md`.
//!
//! **This file is ng's own**, and holds nothing but its declarations, its re-exports and
//! the tests ng adds to production's. The constants a copy of this module would carry in
//! *production's* `mod.rs` are in [`model_params`] instead, precisely so that a file-level
//! guard can reach them.
//!
//! Landed so far (plan `doc/devel/ng/impl_plan/hidden_paralog_filter.md`, Milestone A):
//!
//! - **A1:** the constants and grids ([`ParalogModelParams`] in [`model_params`]), and
//!   the per-sample fit of what one copy's depth looks like
//!   ([`SingleCopyCoverageModel`] in [`coverage_model`]).

pub mod coverage_model;
pub mod model_params;

/// **ng's, not a copy** — the textual check that the copies beside it are still
/// production's, written from outside the files it checks (spec §1.2).
#[cfg(test)]
mod copy_fidelity;

// Production's own surface (`src/paralog/mod.rs`), name for name, so that the call sites
// arriving with the later steps resolve unchanged.
pub use coverage_model::{
    CoverageFitConfig, CoverageModelError, ModeMedianRatioBounds, SingleCopyCoverageModel,
};
pub use model_params::{
    DEFAULT_ALLELE_FREQ_PRIOR_POINTS, DEFAULT_CARRIER_COPY_NUMBERS, DEFAULT_CARRIER_FREQ_HI,
    DEFAULT_CARRIER_FREQ_LO, DEFAULT_CARRIER_FREQ_POINTS, DEFAULT_HOMALT_MIN_DEPTH,
    DEFAULT_HOMALT_VAF_THRESHOLD, DEFAULT_MAX_RELATIVE_COPY_NUMBER, DEFAULT_PSEUDOCOUNT_VAF,
    GridSpec, ParalogModelParams, SfsPriorSpec,
};

/// **ng's own tests, beside production's transcribed ones** — the cases the copied suite
/// does not reach. They live here rather than in [`model_params`] because that file is
/// production's byte for byte and may not gain a line.
#[cfg(test)]
mod tests {
    use super::*;

    /// **An infinite endpoint is the input that tests the finiteness guard; `NaN` is not.**
    /// The transcribed `grid_spec_new_validates_bounds` asserts only
    /// `GridSpec::new(f64::NAN, 0.6, 40)`, and `NAN < 0.6` is `false`, so that call is
    /// refused by the `lo < hi` clause whether or not `is_finite` is there at all —
    /// deleting both `is_finite` calls leaves the whole transcribed suite green. An
    /// infinite `hi` is what separates them: it passes `lo < hi` and would otherwise build
    /// a carrier-frequency grid with no finite upper end, so H2's marginalisation would
    /// integrate over `[0.004, ∞)`.
    #[test]
    fn grid_spec_new_refuses_an_infinite_endpoint() {
        assert!(GridSpec::new(0.004, f64::INFINITY, 40).is_none());
        assert!(GridSpec::new(f64::NEG_INFINITY, 0.6, 40).is_none());
        assert!(GridSpec::new(0.004, f64::NAN, 40).is_none());
    }
}
