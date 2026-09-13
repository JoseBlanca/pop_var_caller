//! **Fit the duplication rate, build the false-discovery curve, and resolve the cut** — the four
//! lines that put the three copied pieces together.
//!
//! [`ParalogPrior::estimate`], [`ParalogFdrCurve::from_histogram`] and
//! [`ParalogFdrCurve::lr_threshold_for_fdr`] are production's, transcribed into
//! [`super::prior`] and compared with production's bit for bit. What is ng's own is here: the
//! **fallback**. An estimate that did not settle is replaced by the documented rate rather than
//! used, because an unconverged iterate is not distinguishable from a real estimate by its value
//! alone and would otherwise calibrate the whole run. `converged` stays `false`, so the run report
//! can say which happened.
//!
//! **This is ng's file, not a copy, and it was never textually guarded** — production keeps its
//! counterpart (`calibrate_from_histogram`) in `src/var_calling/paralog_filter/calibrate.rs`, below
//! four items ng deliberately does not port, so it fell outside the span the copy guard compared.
//! It is checked against production's frozen answers all the same:
//! `the_fallback_and_the_cut_agree_with_productions_bit_for_bit` in `production_parity` is the
//! oracle, beside the nine other differentials.
//!
//! It lives in this module rather than beside its caller because it is statistics: it reads a
//! histogram and a configuration and returns a calibration, and knows nothing about the spill, the
//! passes, or ng's I/O.

use super::{
    CalibrationConfig, ParalogCalibration, ParalogFdrCurve, ParalogLrHistogram, ParalogPrior,
};

/// **The rate, the curve and the cut, from a histogram of likelihood ratios.**
///
/// `target_fdr` is the operator's target false-discovery rate among the records removed. The
/// caller has already checked it is a probability — pass two takes it as a
/// [`TargetFdr`](crate::ng::run::paralog_filter::TargetFdr), which cannot hold anything else.
///
/// Pure: the same histogram and configuration give the same calibration.
pub fn calibrate_from_the_ratio_histogram(
    histogram: &ParalogLrHistogram,
    target_fdr: f64,
    config: &CalibrationConfig,
) -> ParalogCalibration {
    let estimated = ParalogPrior::estimate(histogram, &config.em);
    let prior = if estimated.converged {
        estimated
    } else {
        ParalogPrior {
            prior_probability: config.fallback_prior,
            converged: false,
        }
    };
    let curve = ParalogFdrCurve::from_histogram(histogram, &prior);
    let lr_threshold = curve.lr_threshold_for_fdr(target_fdr);
    ParalogCalibration {
        prior,
        curve,
        lr_threshold,
        target_fdr,
    }
}
