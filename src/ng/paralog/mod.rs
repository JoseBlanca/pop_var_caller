//! **The hidden-duplication filter's statistics, copied from production into ng.**
//!
//! At each locus the filter asks which of two stories better explains what every sample
//! showed: a real variant that some of them carry, or two gene copies the reference
//! collapsed into one, piling their reads onto the same position. Four quantities answer
//! that, and **all four have landed.** The per-sample coverage model says what read depth
//! *one copy* would produce in a sample at a window of a given GC content, so that dividing
//! an observed depth by it gives a copy number — one for a single copy, about two for a
//! collapsed pair. The per-locus score then says how much better the collapsed-pair story
//! fits every sample's copy numbers and allele counts than the real-variant story does. The
//! prior says how common hidden duplications are in this run, estimated from the run's own
//! scores rather than assumed. And the calibration turns each score into a probability and
//! resolves the operator's target false-discovery rate to a cut.
//!
//! Nothing here reads a file, writes a record or knows what a VCF is; the run's wiring will
//! live in `crate::ng::run::paralog_filter` (plan step B1) and has not landed.
//!
//! **Copied, not re-derived.** The model is production's
//! (`doc/devel/specs/hidden_paralog_filter.md`, reformulated for one sample by
//! `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`). Three of the five
//! files beside this one are production's source, line for line and byte for byte:
//! [`coverage_model`], 1,157 lines past its module header; [`locus_score`], 797; and
//! [`model_params`], 236. The other two are ng's own and are what assert the copies —
//! `copy_fidelity` textually, so an edit to *either* tree fails the build instead of
//! drifting quietly, and `production_parity` numerically, by scoring randomised loci
//! through both implementations and comparing every number by bit pattern. What ng
//! changes about the filter is where its numbers come from and when they are known, not
//! what it computes: `doc/devel/ng/spec/hidden_paralog_filter.md`.
//!
//! **This file is ng's own**, and holds nothing but its declarations, its re-exports and
//! the tests ng adds to production's. Production keeps the model's constants in its own
//! `mod.rs`, whose declarations ng cannot share; ng puts them in [`model_params`] instead,
//! precisely so that a file-level guard can reach them — and it does, comparing from the
//! first item's `///` on each side rather than from the end of a shared module header.
//!
//! Landed so far (plan `doc/devel/ng/impl_plan/hidden_paralog_filter.md`, Milestone A):
//!
//! - **A1:** the constants and grids ([`ParalogModelParams`] in [`model_params`]), and
//!   the per-sample fit of what one copy's depth looks like
//!   ([`SingleCopyCoverageModel`] in [`coverage_model`]).
//! - **A2:** the per-locus score ([`score_locus_for_paralogy`] in [`locus_score`]) and its
//!   per-pass tables, with the differential against production's that is the port's real
//!   proof.
//! - **A3:** how common hidden duplications are in this run and where to cut them
//!   ([`ParalogPrior`] and [`ParalogFdrCurve`] in [`prior`], [`calibration::ParalogCalibration`]).

pub mod calibrate;
pub mod calibration;
pub mod coverage_model;
pub mod locus_score;
pub mod model_params;
pub mod prior;

/// **ng's, not a copy** — the textual check that the copies beside it are still
/// production's, written from outside the files it checks (spec §1.2).
#[cfg(test)]
mod copy_fidelity;

/// **ng's, not a copy** — the differential that proves the copied scorer computes what
/// production's does, on randomised inputs, by bit pattern (spec §10).
#[cfg(test)]
mod production_parity;

// Production's own surface (`src/paralog/mod.rs`), name for name, so that the call sites
// arriving with the later steps resolve unchanged.
pub use calibrate::calibrate_from_the_ratio_histogram;
pub use calibration::{CalibrationConfig, DEFAULT_FALLBACK_PARALOG_PRIOR, ParalogCalibration};
pub use coverage_model::{
    CoverageFitConfig, CoverageModelError, ModeMedianRatioBounds, SingleCopyCoverageModel,
};
pub use locus_score::{
    LocusObservations, ParalogScore, ParalogScorePrecompute, SampleObservation,
    score_locus_for_paralogy,
};
pub use model_params::{
    DEFAULT_ALLELE_FREQ_PRIOR_POINTS, DEFAULT_CARRIER_COPY_NUMBERS, DEFAULT_CARRIER_FREQ_HI,
    DEFAULT_CARRIER_FREQ_LO, DEFAULT_CARRIER_FREQ_POINTS, DEFAULT_HOMALT_MIN_DEPTH,
    DEFAULT_HOMALT_VAF_THRESHOLD, DEFAULT_MAX_RELATIVE_COPY_NUMBER, DEFAULT_PSEUDOCOUNT_VAF,
    GridSpec, ParalogModelParams, SfsPriorSpec,
};
pub use prior::{
    DEFAULT_LR_HISTOGRAM_BINS, DEFAULT_LR_HISTOGRAM_HI, DEFAULT_LR_HISTOGRAM_LO, EmConfig,
    ParalogFdrCurve, ParalogLrHistogram, ParalogPrior,
};

/// **ng's own tests, beside production's transcribed ones** — the cases the copied suite
/// does not reach. They live here rather than in [`model_params`] because that file is
/// production's byte for byte and may not gain a line.
#[cfg(test)]
mod tests {
    use super::*;

    /// **A locus is flagged when its tail false-discovery rate is at or below the target,
    /// and never on a ratio that is not a number.**
    ///
    /// ng's own, because production's tests for this type live in the same file as the
    /// spill-streaming driver ng does not copy. The `NaN` half is spec §6 trap 4's input
    /// side: a locus with no usable evidence never entered the histogram, so it contributed
    /// to neither π nor the curve, and classifying on it here would judge it by a rule the
    /// calibration never applied.
    #[test]
    fn a_locus_is_flagged_only_when_its_q_value_meets_the_target() {
        let mut histogram = ParalogLrHistogram::with_defaults();
        for step in 0..2_000 {
            // A tenth strongly positive, the rest well negative — a run with duplications
            // in it, so the curve has somewhere to put a cut.
            histogram.push(if step % 10 == 0 { 18.0 } else { -8.0 });
        }
        let prior = ParalogPrior::estimate(&histogram, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&histogram, &prior);
        // **Sweep the operator's knob rather than fixing it.** With one target, an
        // implementation that read a literal `0.01` in place of `self.target_fdr` — ignoring
        // the knob entirely — answers every probe exactly as the real one does.
        for target_fdr in [0.0, 0.001, 0.01, 0.5, 1.0] {
            let calibration = calibration::ParalogCalibration {
                prior,
                curve: curve.clone(),
                lr_threshold: curve.lr_threshold_for_fdr(target_fdr),
                target_fdr,
            };
            for lr in [-30.0, -8.0, 0.0, 8.0, 18.0, 40.0] {
                assert_eq!(
                    calibration.flags(lr),
                    curve.q_of_lr(lr) <= target_fdr,
                    "at a target of {target_fdr}, a locus at a ratio of {lr} must be flagged \
                     exactly when its tail FDR meets it"
                );
            }
        }

        let target_fdr = 0.01;
        let calibration = calibration::ParalogCalibration {
            prior,
            curve: curve.clone(),
            lr_threshold: curve.lr_threshold_for_fdr(target_fdr),
            target_fdr,
        };
        assert!(
            calibration.flags(18.0),
            "the duplicated end of this run must be flagged at all, or this test asserts \
             nothing about the rule it is checking"
        );
        assert!(!calibration.flags(-30.0), "the variant end must not be");

        // **The boundary is `<=`, not `<`**, and a target set exactly to a q-value the curve
        // produces is the only input that separates the two.
        let exactly_a_q_value_the_curve_produces = curve.q_of_lr(-8.0);
        let at_the_boundary = calibration::ParalogCalibration {
            prior,
            curve: curve.clone(),
            lr_threshold: None,
            target_fdr: exactly_a_q_value_the_curve_produces,
        };
        assert!(
            at_the_boundary.flags(-8.0),
            "a locus whose tail FDR is exactly the target is flagged: the rule is `<=`"
        );

        for unscorable in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                !calibration.flags(unscorable),
                "a ratio that is not a number entered neither pi nor the curve, so it is \
                 never flagged"
            );
            assert_eq!(
                calibration.posterior(unscorable),
                None,
                "and it has no posterior either"
            );
        }
    }

    /// **The recorded cut and the flag agree to within one histogram bin, and no closer.**
    ///
    /// `calibration.rs` documents `flags` as *"equivalent to `lr >= lr_threshold` by the
    /// curve's monotonicity"*. **That is true only up to the bin.** `lr_threshold_for_fdr`
    /// returns the crossing bin's *centre* while `flags` decides on the bin, so every ratio
    /// in the lower half of that bin is flagged while sitting below the recorded cut — a
    /// window half a bin wide, `0.05` on the likelihood-ratio axis at the shipped resolution
    /// of 2,000 bins over `[-100, 100]`.
    ///
    /// It matters because `lr_threshold` is what the run writes into the VCF header for
    /// provenance (spec §3.3): an operator reading *"records at or above this ratio were
    /// dropped"* is wrong for every record inside that half-bin. Nothing in ng read the field
    /// at all before this test, and step C3 will wire it into the header.
    #[test]
    fn the_recorded_cut_and_the_flag_agree_to_within_one_bin() {
        let mut histogram = ParalogLrHistogram::with_defaults();
        for step in 0..2_000 {
            histogram.push(if step % 10 == 0 { 18.0 } else { -8.0 });
        }
        let prior = ParalogPrior::estimate(&histogram, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&histogram, &prior);
        let target_fdr = 0.01;
        let calibration = calibration::ParalogCalibration {
            prior,
            curve: curve.clone(),
            lr_threshold: curve.lr_threshold_for_fdr(target_fdr),
            target_fdr,
        };
        let cut = calibration
            .lr_threshold
            .expect("this fixture's target is reachable, or the test asserts nothing");
        let half_a_bin = 0.5 * (100.0 - -100.0) / 2_000.0;

        let mut disagreements = 0usize;
        for step in -20..=20 {
            let lr = cut + 0.01 * f64::from(step);
            if calibration.flags(lr) != (lr >= cut) {
                disagreements += 1;
                assert!(
                    (lr - cut).abs() <= half_a_bin,
                    "the flag and the recorded cut may disagree only inside the crossing \
                     bin; at a ratio of {lr} they disagree {} away from a cut of {cut}",
                    (lr - cut).abs(),
                );
            }
        }
        assert!(
            disagreements > 0,
            "they do disagree inside the bin — if this fixture stops showing it, the claim \
             above is no longer being tested"
        );
        assert!(calibration.flags(cut), "the recorded cut is itself flagged");
        assert!(
            !calibration.flags(cut - 2.0 * half_a_bin - 1.0),
            "a ratio a whole bin below the cut is not"
        );
    }

    /// **The posterior is the prior's log-odds shifted by the likelihood ratio**, and is
    /// withheld where that shift is undefined.
    ///
    /// `σ(LR + logit π)`. Asserted against the formula rather than against stored values, so
    /// that it pins the relation and not one arithmetic path through it; and at a degenerate
    /// π the field is omitted rather than emitted saturated at 0 or 1.
    #[test]
    fn the_posterior_is_the_logistic_of_the_ratio_shifted_by_the_prior_odds() {
        let mut histogram = ParalogLrHistogram::with_defaults();
        for step in 0..1_000 {
            histogram.push(if step % 8 == 0 { 12.0 } else { -6.0 });
        }
        let prior = ParalogPrior::estimate(&histogram, &EmConfig::default());
        let pi = prior.prior_probability;
        assert!(
            pi > 0.0 && pi < 1.0,
            "this fixture must fit a usable pi, or the comparison below is vacuous; got {pi}"
        );
        let calibration = calibration::ParalogCalibration {
            prior,
            curve: ParalogFdrCurve::from_histogram(
                &histogram,
                &ParalogPrior::estimate(&histogram, &EmConfig::default()),
            ),
            lr_threshold: None,
            target_fdr: 0.01,
        };

        for lr in [-20.0, -1.0, 0.0, 1.0, 12.0, 35.0] {
            let log_odds = lr + (pi / (1.0 - pi)).ln();
            let expected = 1.0 / (1.0 + (-log_odds).exp());
            assert_eq!(
                calibration.posterior(lr).map(f64::to_bits),
                Some(expected.to_bits()),
                "at a ratio of {lr} the posterior must be the logistic of the ratio plus \
                 the prior's log-odds"
            );
        }

        for degenerate in [0.0, 1.0] {
            let mut with_a_degenerate_prior = calibration.clone();
            with_a_degenerate_prior.prior.prior_probability = degenerate;
            assert_eq!(
                with_a_degenerate_prior.posterior(3.0),
                None,
                "at a prior of {degenerate} the log-odds are undefined, so the field is \
                 omitted rather than saturated"
            );
        }
    }

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
