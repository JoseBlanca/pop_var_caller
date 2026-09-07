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
//! produce a plausible file. Three things stand against it, and the third was missing until the
//! step's own review: [`ParalogScoringContext::score`] returns `NaN` rather than the scorer's
//! neutral `0.0`; the ratio pushed to the histogram and the ratio pushed to the vector are **the
//! same binding**, so no rule can apply to one and not the other; and the verdict screens on
//! finiteness before it consults the curve — which matters more than it looks, because the
//! curve's own answer for a `NaN` sits *inside* a target of one in two by about
//! 1.5 × 10⁻¹⁴. `an_unscored_record_is_never_flagged_however_loose_the_target` is what holds that.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.3, §6 trap 4.

use crate::ng::paralog::{
    CalibrationConfig, DEFAULT_LR_HISTOGRAM_BINS, DEFAULT_LR_HISTOGRAM_HI, DEFAULT_LR_HISTOGRAM_LO,
    ParalogCalibration, ParalogLrHistogram, SampleObservation, calibrate_from_the_ratio_histogram,
};

use super::{CohortSizeMismatch, ParalogScoringContext, SpillFile, SpillFileError};

/// **The operator's target false-discovery rate**, checked once at the boundary.
///
/// The share of the records the filter removes that were really variants: `0.01` means about one
/// wrong removal in a hundred.
///
/// **Zero is not one of these, and that is the point.** `--paralog-fdr 0` means *do not run the
/// filter*, and it has to mean that somewhere a caller cannot forget. Handed a target of zero the
/// scoring does **not** remove nothing — a strongly duplicated record's tail false-discovery value
/// underflows to exactly zero, and zero is not above zero, so the most extreme records go. Until
/// this type refused it, the only thing standing between that and an operator who asked for no
/// filtering was a `> 0.0` comparison copied into two call sites.
/// [`WhatTheOperatorAskedFor::from_the_flags`](crate::ng::run::paralog_filter::WhatTheOperatorAskedFor::from_the_flags)
/// is where zero becomes *no filter*, once.
///
/// **A newtype because both ways of getting it wrong are silent.** Handed a bare `f64`, the
/// calibration compares it against the curve and asks no questions — measured on a three-record
/// spill: a target of `5`, which is what someone typing five for "five percent" gets, removes
/// **every** scored record; a negative one, or a `NaN` arriving through some upstream arithmetic,
/// removes none and reports no cut, which is exactly what the calibration reports when a target
/// is simply unreachable. So a run report cannot tell a wrong knob from a clean cohort. The crate
/// already draws this line for [`InbreedingF`](crate::ng::types::InbreedingF); this mirrors it.
///
/// The admitted range is `[0, 1)`, the same range the `--paralog-fdr` flag already parses, so the
/// flag and this type cannot disagree about what a target is.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct TargetFdr(f64);

impl TargetFdr {
    /// The only constructor. A target that is not a fraction strictly between `0` and `1` is
    /// refused rather than coerced.
    ///
    /// **Zero is refused here** — it means *no filter*, which is a different thing from a target,
    /// and it is spelled by the absence of one of these. See the type's own note.
    ///
    /// **`-0.0` is refused too, and it is the one that needs saying.** It compares equal to zero,
    /// so a sign-blind check would read it as *off* — but a run asking for a negative target has
    /// asked for something, and answering "the filter did not run" is the wrong reply to a
    /// mistake.
    ///
    /// # Errors
    ///
    /// If the value is not finite, is zero or negative — negative zero included — or is `1` or
    /// more.
    pub fn try_new(target: f64) -> Result<Self, NotATargetFdr> {
        if target.is_finite() && target.is_sign_positive() && target > 0.0 && target < 1.0 {
            Ok(Self(target))
        } else {
            Err(NotATargetFdr { given: target })
        }
    }

    /// The target as the copied calibration takes it.
    #[must_use]
    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// **What was offered as a target false-discovery rate is not one.**
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[error(
    "the target false-discovery rate is {given}; it must be a fraction of one, above 0 and below \
     1 — 0.01 means about one wrongly removed record in a hundred, not one in a hundredth. Use \
     exactly 0 to run without the filter"
)]
pub struct NotATargetFdr {
    /// The value that was offered.
    pub given: f64,
}

/// **The range and resolution the run's ratios were binned over.**
///
/// Recorded because the cut depends on it and nothing else in the output would say what it was: a
/// ratio past either edge is folded into that end's bin, so a run whose ratios all sit outside the
/// range leaves the curve one occupied bin to work with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LrHistogramShape {
    /// The lowest ratio the histogram distinguishes; anything below folds into the first bin.
    pub lowest_ratio: f64,
    /// The highest; anything above folds into the last bin.
    pub highest_ratio: f64,
    /// How many bins the range is cut into.
    pub bins: usize,
}

impl LrHistogramShape {
    /// The shipped range and resolution — production's, inherited.
    #[must_use]
    pub fn shipped() -> Self {
        Self {
            lowest_ratio: DEFAULT_LR_HISTOGRAM_LO,
            highest_ratio: DEFAULT_LR_HISTOGRAM_HI,
            bins: DEFAULT_LR_HISTOGRAM_BINS,
        }
    }

    /// **An empty histogram of this shape — the only way pass two builds one.**
    ///
    /// The shape is recorded on the verdicts and the histogram is what the ratios are folded
    /// into, and the two must be the same thing: a run reporting a range it did not use would
    /// mis-count how many ratios saturated and mis-describe the cut. Building the histogram
    /// *from* the shape is what stops them drifting — there is no second place where a range is
    /// written down. A mutation that changed the range at the old call site survived the whole
    /// suite, which is why this is a method rather than a comment.
    ///
    /// # Panics
    ///
    /// If the shape is not a range a histogram can be cut from — impossible for
    /// [`Self::shipped`], which is the only shape pass two builds.
    #[must_use]
    pub fn histogram(self) -> ParalogLrHistogram {
        ParalogLrHistogram::new(self.lowest_ratio, self.highest_ratio, self.bins)
            .expect("the likelihood-ratio histogram's shape is a valid range")
    }
}

/// **What pass two settles for pass three**: how common duplications are in this run, where the
/// cut falls, and every record's ratio in spill order.
///
/// Spec §3.7's `ParalogVerdicts`, with four additions the run report needs — see
/// [`Self::records_in_the_fit`], [`Self::ratios_outside_the_histogram`], [`Self::config`] and
/// [`Self::lr_histogram`].
#[derive(Debug)]
pub struct ParalogVerdicts {
    /// The fitted rate, the false-discovery curve, and the resolved cut. Its `flags` and
    /// `posterior` are what pass three asks about each record.
    pub calibration: ParalogCalibration,
    /// One ratio a record, **in spill order**, `NaN` where the record was not scored.
    ///
    /// **A ratio's only tie to its record is its position here.** Pass three walks the spill again
    /// and must take these front to back, in step with its own read — `zip`, not an index. Nothing
    /// in the type enforces that, and it is the one invariant this step cannot test because its
    /// consumer does not exist yet: a pass three that read the spill filtered, chunked or sharded
    /// would give every record its neighbour's verdict, with nothing about the file looking wrong.
    pub ratios: Vec<f64>,
    /// How many records the fitted rate actually rests on — the ratios that were folded, which
    /// are the finite ones.
    ///
    /// **Kept because the two counts can only be told apart from outside.** A run over a million
    /// records whose coverage models were all rejected produces a million `NaN`s, falls back, and
    /// writes a VCF that looks exactly like a run where the filter worked; the difference between
    /// "a million records in the fit" and "none" is this number.
    pub records_in_the_fit: u64,
    /// How many of those ratios landed past an end of the histogram and were folded into its end
    /// bin.
    ///
    /// **The score grows with the cohort**, because it is a sum over the samples that were
    /// weighed, while the histogram's range is fixed — measured on one duplication-shaped record,
    /// 24.2 at one sample, 156.3 at six, 1,663 at 63, against a range of ±100. Saturating is
    /// harmless while it happens to one class, whose probability is saturated anyway. It stops
    /// being harmless if a cohort is large enough that real variants *and* duplications both land
    /// past the same edge: they then share one bin and the operator's target has nothing left to
    /// move. Nothing acts on this count yet; the plan's D2 and D3 are the runs that report it on
    /// real data.
    pub ratios_outside_the_histogram: u64,
    /// The rate-fitting knobs the run used — the iteration's start, tolerance and cap, and the
    /// fallback rate.
    ///
    /// **Recorded because they are inherited constants that nothing has re-measured.** Spec §3.1
    /// says the model's are prototype-tuned on tomato2 and taken as given; a number in that
    /// position which no run ever prints is a number nobody checks.
    pub config: CalibrationConfig,
    /// The range and resolution the ratios were binned over.
    pub lr_histogram: LrHistogramShape,
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
            self.records_in_the_fit, self.calibration.prior.prior_probability,
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
/// `config` carries the fitting knobs and the fallback rate; [`CalibrationConfig::default`] is
/// production's.
///
/// # Errors
///
/// If the spill cannot be read — including one that ends before, or runs past, the number of
/// records pass one counted — or if a parked record carries a different number of samples than
/// the run has.
pub fn score_the_parked_records_and_resolve_the_cut(
    spill: &SpillFile,
    context: &ParalogScoringContext,
    target_fdr: TargetFdr,
    config: &CalibrationConfig,
) -> Result<ParalogVerdicts, PassTwoError> {
    let entries = spill.read().map_err(PassTwoError::Spill)?;

    // **The shape is written down once and the histogram is built from it**, so the range the
    // run reports and the range the ratios were folded into cannot be two different ranges.
    let lr_histogram = LrHistogramShape::shipped();
    let mut histogram = lr_histogram.histogram();

    // Pass one counted them, so the vector is sized once rather than doubling its way up to
    // millions of records — 40 MB at the five million spec §3.3 contemplates. On a target too
    // narrow to hold that count the reservation is skipped and the vector grows by doubling
    // instead, which costs time and never correctness.
    let mut ratios: Vec<f64> =
        Vec::with_capacity(usize::try_from(spill.entries_written()).unwrap_or(0));
    let mut observations: Vec<Option<SampleObservation>> =
        Vec::with_capacity(context.sample_count());

    for entry in entries {
        let entry = entry.map_err(|source| PassTwoError::Spill(spill.naming_this_file(source)))?;
        let ratio = context
            .score(&entry, &mut observations)
            .map_err(PassTwoError::CohortSize)?;
        // **One binding, pushed to both.** The histogram drops a non-finite ratio and the vector
        // keeps it, which is the whole difference between "did not enter the fit" and "kept,
        // never flagged" — but they are the same number, so no rule can apply to one and not the
        // other (spec §6 trap 4).
        histogram.push(ratio);
        ratios.push(ratio);
        // `entry` is dropped here: its line and its per-sample rows are not carried forward.
    }

    // The reader refuses a file holding a different number of records than pass one counted, in
    // either direction, so arriving here with a different length would be a defect in the reader
    // rather than in the file — and pass three pairs the two by position.
    debug_assert_eq!(
        ratios.len() as u64,
        spill.entries_written(),
        "pass three pairs each spilled record with the ratio at its own index, so the two counts \
         cannot differ"
    );

    let records_in_the_fit = histogram.total();
    let ratios_outside_the_histogram = ratios
        .iter()
        .filter(|ratio| {
            ratio.is_finite()
                && (**ratio < lr_histogram.lowest_ratio || **ratio > lr_histogram.highest_ratio)
        })
        .count() as u64;

    Ok(ParalogVerdicts {
        calibration: calibrate_from_the_ratio_histogram(&histogram, target_fdr.get(), config),
        ratios,
        records_in_the_fit,
        ratios_outside_the_histogram,
        config: *config,
        lr_histogram,
    })
}

#[cfg(test)]
mod tests;
