//! **What a run does after its calling pass, when the hidden-duplication filter is on.**
//!
//! Fit each sample's coverage model, score every parked record, resolve the operator's target to
//! a cut, then open the VCF and write the records the cut keeps. One function, called from both
//! subcommands with the same arguments, because the two must not be able to drift into filtering
//! differently — spec §1.1 goal 2 asks that direct mode and psp mode run the same three passes on
//! the same numbers, and a shared call is what makes that checkable rather than a convention.
//!
//! **The VCF is opened here, after the verdicts exist.** Nothing before this point could open it:
//! no record's fate is known until every record has been scored, and a half-written VCF beside a
//! run that has decided nothing is a file someone can mistake for the output. It also means the
//! header can state what the filter was calibrated to, which a dropped record's absence otherwise
//! leaves no trace of.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.1, §3.5, §3.6.

use std::path::{Path, PathBuf};

use crate::ng::paralog::{CalibrationConfig, CoverageFitConfig, ParalogModelParams};
use crate::ng::types::{InbreedingF, Ploidy};
use crate::ng::vcf::{HiddenParalogProvenance, VcfHeaderMetadata, VcfWriteError, VcfWriter};
use crate::ng::window_coverage::SampleHistogram;

use super::{
    CoverageFitConfigRefused, ParalogScoringContext, ParalogVerdicts, PassThreeError, PassTwoError,
    SpillFile, TargetFdr, WhatTheFilterDid, score_the_parked_records_and_resolve_the_cut,
    write_the_records_the_filter_kept,
};

/// **The two knobs the operator sets, together** — because they are one decision.
///
/// A target of zero means the filter does not run at all and the tag flag has nothing to act on;
/// above zero, the tag flag says whether a removed record leaves the file or stays in it on the
/// filter's id. Passing them as two arguments among eight let a call site swap them past the type
/// checker, since one is a target and the other a plain `bool`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhatTheOperatorAskedFor {
    /// The share of the removed records that may really have been variants — `--paralog-fdr`.
    pub target_fdr: TargetFdr,
    /// Keep a removed record on the `hiddenParalog` filter instead of leaving it out —
    /// `--paralog-filter-tag`. **What it is for**: a run that drops can be audited against one
    /// that tags, since the two files differ by exactly the tagged lines.
    pub tag_instead_of_dropping: bool,
}

/// **What the filter did, and everything the run report needs to say so.**
#[derive(Debug)]
pub struct FilteredRun {
    /// The records written, dropped, tagged and left unscored.
    pub did: WhatTheFilterDid,
    /// The fitted rate, the cut, and every record's ratio.
    pub verdicts: ParalogVerdicts,
    /// The per-sample coverage models — kept for
    /// [`why_no_model`](ParalogScoringContext::why_no_model), which spec §3.5's report line names
    /// each rejected sample from.
    pub scoring: ParalogScoringContext,
}

/// What can stop a run once its records are parked.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ParalogFilterError {
    /// The coverage fit was configured with knobs that do not describe a fit, so no sample could
    /// be fitted. **An operator's mistake, not a verdict about the cohort** — which is why it is
    /// an error here rather than sixty-three identical per-sample refusals.
    #[error("the hidden-duplication filter's coverage fit could not be configured")]
    TheCoverageFitIsNotConfigured(#[source] CoverageFitConfigRefused),
    /// The parked records could not be scored.
    #[error("the parked records could not be scored for hidden duplications")]
    Scoring(#[source] PassTwoError),
    /// The VCF could not be created. **The records are still parked at this point**, and the
    /// spill's own guard removes them as the run unwinds.
    #[error("the calls could not be written to {}", path.display())]
    OutputNotCreated {
        /// Where they were going.
        path: PathBuf,
        /// What the writer said.
        #[source]
        source: VcfWriteError,
    },
    /// A record could not be given its verdict or written.
    #[error("the filtered calls could not be written to {}", path.display())]
    Writing {
        /// Where they were going.
        path: PathBuf,
        /// What pass three said.
        #[source]
        source: PassThreeError,
    },
    /// The finished file could not be put in place.
    #[error("the calls could not be finished at {}", path.display())]
    CallsNotFinished {
        /// Where they were going.
        path: PathBuf,
        /// What the writer said.
        #[source]
        source: VcfWriteError,
    },
}

/// **Fit, score, calibrate and write** — passes two and three, and the fit between them.
///
/// `histograms` and `inbreeding` are both in the run's sample order and are the same length;
/// `metadata` is the header the run would have written, which this adds the filter's provenance
/// to before opening the file.
///
/// # Panics
///
/// If the histograms and the coefficients disagree in length — they come from one run indexed by
/// one sample order, so a mismatch is a wiring error rather than an input.
///
/// # Errors
///
/// If the fit is misconfigured, if the spill cannot be read, or if the VCF cannot be written.
pub fn fit_score_and_write_the_calls(
    spill: &SpillFile,
    histograms: Vec<SampleHistogram>,
    inbreeding: &[InbreedingF],
    asked_for: WhatTheOperatorAskedFor,
    output: &Path,
    metadata: VcfHeaderMetadata,
    ploidy: Ploidy,
) -> Result<FilteredRun, ParalogFilterError> {
    let WhatTheOperatorAskedFor {
        target_fdr,
        tag_instead_of_dropping,
    } = asked_for;
    // **Between the passes, and only here.** The histograms are complete when the calling pass
    // ends and are consumed by the fit, so each sample's bins are freed as its model is built.
    let scoring = ParalogScoringContext::new(
        histograms,
        inbreeding,
        &ParalogModelParams::default(),
        &CoverageFitConfig::default(),
    )
    .map_err(ParalogFilterError::TheCoverageFitIsNotConfigured)?;

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        spill,
        &scoring,
        target_fdr,
        &CalibrationConfig::default(),
    )
    .map_err(ParalogFilterError::Scoring)?;

    // **The header states what the run was calibrated to**, because a dropped record leaves no
    // trace of it anywhere else in the file.
    let metadata = metadata.the_hidden_duplication_filter_ran(HiddenParalogProvenance {
        target_fdr: verdicts.calibration.target_fdr,
        paralog_rate: verdicts.calibration.prior.prior_probability,
        lr_cut: verdicts.calibration.lr_threshold,
        rate_was_fitted: verdicts.calibration.prior.converged,
        samples_with_a_coverage_model: scoring.how_many_samples_have_a_coverage_model(),
        samples_in_the_run: scoring.sample_count(),
    });

    let mut writer = VcfWriter::create(output, metadata, ploidy).map_err(|source| {
        ParalogFilterError::OutputNotCreated {
            path: output.to_path_buf(),
            source,
        }
    })?;

    let did =
        write_the_records_the_filter_kept(spill, &verdicts, tag_instead_of_dropping, &mut writer)
            .map_err(|source| ParalogFilterError::Writing {
            path: output.to_path_buf(),
            source,
        })?;

    writer
        .finish()
        .map_err(|source| ParalogFilterError::CallsNotFinished {
            path: output.to_path_buf(),
            source,
        })?;

    Ok(FilteredRun {
        did,
        verdicts,
        scoring,
    })
}

/// **The lines spec §3.5 asks a run report for**, when the filter ran.
///
/// Words rather than fields, because this is what an operator reads: how many records went, how
/// much of the cohort the coverage evidence rested on and why the rest dropped out, the fitted
/// rate and whether it was fitted at all.
#[must_use]
pub fn what_to_tell_the_operator(run: &FilteredRun) -> Vec<String> {
    let FilteredRun {
        did,
        verdicts,
        scoring,
    } = run;
    let mut lines = vec![format!(
        "hidden-duplication filter: {} record(s) dropped, {} tagged, {} written; \
         {} were scored, {} could not be",
        did.dropped, did.tagged, did.written, verdicts.records_in_the_fit, did.unscored,
    )];

    lines.push(format!(
        "  duplication rate {:.6} ({}), cut at {}, over {} of {} sample(s) with a coverage model",
        verdicts.calibration.prior.prior_probability,
        if verdicts.calibration.prior.converged {
            "fitted from this run"
        } else {
            "the documented fallback — this run's could not be fitted"
        },
        match verdicts.calibration.lr_threshold {
            Some(cut) => format!("{cut:.4}"),
            None => "no ratio the target reaches".to_string(),
        },
        scoring.how_many_samples_have_a_coverage_model(),
        scoring.sample_count(),
    ));

    // **Named, not counted.** "The model rested on 61 of 63" and "on 61 of 63, and the two that
    // dropped out were nearly uncovered" are different statements, and only the second tells an
    // operator whether to look at those samples' reads.
    for (sample, why) in scoring.why_no_model().iter().enumerate() {
        if let Some(why) = why {
            lines.push(format!("  sample {sample} has no coverage model: {why}"));
        }
    }

    if verdicts.ratios_outside_the_histogram > 0 {
        lines.push(format!(
            "  {} scored record(s) sat past the ends of the ratio range [{}, {}] and were \
             folded into its end bins",
            verdicts.ratios_outside_the_histogram,
            verdicts.lr_histogram.lowest_ratio,
            verdicts.lr_histogram.highest_ratio,
        ));
    }

    lines
}

#[cfg(test)]
mod tests;
