//! **Pass two: score every parked record, and turn the scores into a cut.**
//!
//! Pass one parked each finished record on the spill because no record's verdict can be reached
//! until every record has been scored — the rate at which hidden duplications occur in *this
//! run* is fitted from the run's own scores, and the false-discovery cut is a quantile of the
//! resulting probabilities. This is the pass that does both: it reads the spill once, scores
//! each record through [`ParalogScoringContext`], keeps every ratio, fits the rate, and resolves
//! the operator's target false-discovery rate to a likelihood ratio to cut at.
//!
//! **The ratios are kept, one `f64` a record, and that is deliberate.** Pass three could
//! re-score instead, and would then be applying the cut to numbers the cut was not computed
//! from if anything about the scoring drifted between the two passes. Keeping them makes the
//! two agree by construction. Eight bytes a record is forty megabytes at five million records,
//! and it does not grow with the cohort.
//!
//! **The failure this pass is built against is a `NaN` becoming a `0.0`** (spec §6 trap 4). A
//! record no sample could speak for is *unscored*: it is kept, never flagged, and it must not
//! enter the fitted rate. Its ratio is `NaN`, and the histogram refuses to fold a `NaN`. Turn it
//! into a zero anywhere along the way and it becomes a record voting "not a duplication" — so a
//! run where nothing could be scored would fit a duplication rate out of nothing at all and
//! produce a plausible file. Two things stand against it here: [`ParalogScoringContext::score`]
//! returns `NaN` rather than the scorer's neutral `0.0`, and the ratio pushed to the histogram
//! and the ratio pushed to the vector are **the same binding**, so no rule can apply to one and
//! not the other.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.3, §6 trap 4.

use crate::ng::paralog::{
    CalibrationConfig, ParalogCalibration, ParalogFdrCurve, ParalogLrHistogram, ParalogPrior,
    SampleObservation,
};

use super::{CohortSizeMismatch, ParalogScoringContext, SpillFile, SpillFileError};

/// **What pass two settles for pass three**: how common duplications are in this run, where the
/// cut falls, and every record's ratio in spill order.
///
/// Spec §3.7's `ParalogVerdicts`, with the count of scored records added — see
/// [`Self::records_scored`].
#[derive(Debug, Clone)]
pub struct ParalogVerdicts {
    /// The fitted rate, the false-discovery curve, and the resolved cut. Its `flags` and
    /// `posterior` are what pass three asks about each record.
    pub calibration: ParalogCalibration,
    /// One ratio a record, **in spill order**, `NaN` where the record was not scored. Pass three
    /// walks the spill again in the same order, so the *i*th entry it reads is this *i*th ratio.
    pub ratios: Vec<f64>,
    /// How many records the calibration actually rests on — the number of ratios that were
    /// folded, which is the number that are finite.
    ///
    /// **Kept because the two counts can only be told apart from outside.** A run over a million
    /// records whose coverage models were all rejected produces a million `NaN`s, fits the
    /// fallback rate, and writes a VCF that looks exactly like a run where the filter worked;
    /// the difference between "a million records scored" and "none" is this number, and spec
    /// §3.5's run-report line is where an operator sees it.
    pub records_scored: u64,
}

impl ParalogVerdicts {
    /// **What to tell the operator when the fitted duplication rate is not a fitted rate.**
    ///
    /// `Some` when the estimate did not settle — the iteration ran out, or there was nothing to
    /// fit — in which case the run used the documented fallback instead. Spec §3.3 asks for a
    /// warning; this is its words, and the run report prints them. `None` when the estimate
    /// converged, which is the ordinary case.
    #[must_use]
    pub fn why_the_paralog_rate_is_not_fitted(&self) -> Option<String> {
        if self.calibration.prior.converged {
            return None;
        }
        Some(format!(
            "the hidden-duplication rate could not be fitted from this run's {} scored \
             record(s), so the filter used the fallback rate of {:.6}; the records it removes \
             are calibrated against that rate and not against this run",
            self.records_scored, self.calibration.prior.prior_probability,
        ))
    }
}

/// What can go wrong in pass two.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PassTwoError {
    /// The parked records could not be read back.
    #[error("the records parked for the hidden-duplication filter could not be read back")]
    Spill(#[source] SpillFileError),
    /// A parked record does not carry the run's cohort.
    #[error("a parked record does not describe the run's cohort")]
    CohortSize(#[source] CohortSizeMismatch),
}

/// **Score every record pass one parked, and resolve the target false-discovery rate to a cut.**
///
/// One sequential read of the spill. Nothing is held across records but the ratios: each entry
/// is decoded, scored, and dropped, and the observation buffer is refilled rather than
/// reallocated.
///
/// `target_fdr` is the operator's `--paralog-fdr`. `config` carries the fitting knobs and the
/// fallback rate; [`CalibrationConfig::default`] is production's.
///
/// # Errors
///
/// If the spill cannot be read — including a spill that ends before the number of records pass
/// one counted — or if a parked record carries a different number of samples than the run has.
pub fn score_the_parked_records_and_resolve_the_cut(
    spill: &SpillFile,
    context: &ParalogScoringContext,
    target_fdr: f64,
    config: &CalibrationConfig,
) -> Result<ParalogVerdicts, PassTwoError> {
    let entries = spill.read().map_err(PassTwoError::Spill)?;

    // Pass one counted them, so the vector is sized once rather than doubling its way up to
    // millions of records.
    let mut ratios: Vec<f64> =
        Vec::with_capacity(usize::try_from(spill.entries_written()).unwrap_or(0));
    let mut histogram = ParalogLrHistogram::with_defaults();
    let mut observations: Vec<Option<SampleObservation>> =
        Vec::with_capacity(context.sample_count());

    for entry in entries {
        let entry = entry.map_err(|source| PassTwoError::Spill(spill.naming_this_file(source)))?;
        let ratio = context
            .score(&entry, &mut observations)
            .map_err(PassTwoError::CohortSize)?;
        // **One binding, pushed to both.** The histogram drops a non-finite ratio and the vector
        // keeps it, which is the whole difference between "did not enter the fit" and "kept,
        // never flagged" — but they are the same number, so no later edit can fold one value and
        // record another (spec §6 trap 4).
        histogram.push(ratio);
        ratios.push(ratio);
        // `entry` is dropped here: its line and its per-sample rows are not carried forward.
    }

    let records_scored = histogram.total();
    Ok(ParalogVerdicts {
        calibration: calibrate_from_the_ratio_histogram(&histogram, target_fdr, config),
        ratios,
        records_scored,
    })
}

/// **Fit the duplication rate, build the false-discovery curve, and resolve the cut.**
///
/// ng's own three lines over three copied pieces — the estimate, the curve and the threshold are
/// all [`crate::ng::paralog`]'s, transcribed from production and checked against it bit for bit.
/// What is ng's is the *fallback*: an estimate that did not settle is replaced by the documented
/// rate rather than used, because an unconverged iterate is not distinguishable from a real
/// estimate by its value alone and would silently calibrate the whole run. `converged` stays
/// `false` so the run report can say which happened.
///
/// It lives here rather than beside the copied statistics because its only caller is pass two,
/// and because the file it would join is compared byte for byte against production's and may not
/// gain a line. Production's counterpart is `calibrate_from_histogram`
/// (`src/var_calling/paralog_filter/calibrate.rs`), and the test
/// `the_fallback_and_the_cut_agree_with_productions_bit_for_bit` is what says the two do the same
/// thing.
fn calibrate_from_the_ratio_histogram(
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

#[cfg(test)]
mod tests;
