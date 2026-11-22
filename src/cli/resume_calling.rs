//! **`pvc resume-calling`: finishes a calling run from the scratch directory it kept** — the
//! hidden-duplication filter's scoring and writing, refinement among them, without calling again.
//!
//! A run with the filter on parks every called record on a spill in its scratch directory and,
//! once calling has ended, keeps beside it what the filter's passes need
//! ([`crate::run::paralog_filter::resume`]). When such a run stops later on — refinement over its
//! memory budget, a full disk, a power cut — and kept the directory (`--keep-tmp-on-error`, or a
//! killed process), this command picks it up there: it scores the spill and writes the VCF as the
//! run would have, through the same function. **On a cohort of thousands that is days of calling
//! not done twice.**
//!
//! **What may change on resuming** is what the writing does and nothing that would change what was
//! called: refinement's window budget and smallest window, the share of agreeing paths phase asks
//! for, the threads, where the VCF and the log go. The VCF's header states the run's own command
//! line and the refinement settings this run used.
//!
//! **What a resumed run cannot report** is calling's own report, which the stopped run logged; it
//! logs the filter's, refinement's and phase's.
use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use thiserror::Error;

use crate::cli::calling_run::{self, CallingRunError};
use crate::cli::run_output::{self, LogFileError};
use crate::neighbour_evidence::phase::PathAgreement;
use crate::run::paralog_filter::resume::{self, ResumeError};
use crate::run::paralog_filter::{self, ParalogFilterError};
use crate::vcf::writer::is_standard_output;

/// Finish a calling run from the scratch directory it kept: the hidden-duplication filter's
/// scoring and writing, refinement among them, without calling again.
///
/// The directory is the `pvc_tmp_*` a run with the hidden-dup filter on made beside its VCF, or in
/// --tmp-dir, and kept when it stopped after calling had ended (--keep-tmp-on-error, or a killed
/// run). It is removed when this run finishes, and kept again if it fails with
/// --keep-tmp-on-error.
#[derive(Debug, Args)]
pub struct ResumeCallingArgs {
    /// The scratch directory the stopped run kept, holding its spill and `resume.json`.
    #[arg(long)]
    pub scratch_dir: PathBuf,

    /// Where to write the VCF; `-` writes it to standard output, uncompressed. A name ending in
    /// `.vcf.gz` or `.vcf.bgz` is written compressed (bgzip); any other is plain text.
    #[arg(long)]
    pub output: PathBuf,

    /// Threads to use; 0 uses every core. The output is the same at any number.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Replace the VCF and the log if they are already there.
    #[arg(long)]
    pub force: bool,

    /// Where to write this run's log. Defaults to `pvc_calls.log` beside the VCF.
    #[arg(long)]
    pub log: Option<PathBuf>,

    /// Keep the scratch directory again if this run stops on an error. A run that finishes
    /// removes it.
    #[arg(long)]
    pub keep_tmp_on_error: bool,

    /// The most memory, in bytes, one window of genotype refinement may hold, in place of the
    /// stopped run's. See `pvc call-from-psps --help`.
    #[arg(long, help_heading = "Advanced")]
    pub refine_window_bytes: Option<u64>,

    /// The shortest each margin of a refinement window may become, in bases, in place of the
    /// stopped run's. See `pvc call-from-psps --help`.
    #[arg(
        long,
        value_name = "BASES",
        value_parser = clap::value_parser!(u64).range(1..=crate::neighbour_evidence::windows::DEFAULT_MARGIN_BASES),
        help_heading = "Advanced"
    )]
    pub refine_window_floor: Option<u64>,

    /// The share of drawn paths that must agree to phase two neighbouring sites together, in
    /// place of the stopped run's. See `pvc call-from-psps --help`.
    #[arg(long, value_parser = calling_run::parse_path_agreement, help_heading = "Advanced")]
    pub phase_min_path_agreement: Option<PathAgreement>,
}

/// Why a resumed run could not finish.
#[derive(Debug, Error)]
pub enum ResumeCallingCliError {
    /// What every calling command can refuse.
    #[error(transparent)]
    Calling(#[from] CallingRunError),
    /// The log could not be opened.
    #[error(transparent)]
    Log(#[from] LogFileError),
    /// The output is already there and --force was not given.
    #[error("{} is already there; give --force to replace it", path.display())]
    OutputAlreadyThere {
        /// The output.
        path: PathBuf,
    },
    /// The scratch directory holds no resume point the run can be finished from.
    #[error("the run cannot be resumed from `{}`", directory.display())]
    NotResumable {
        /// The directory given.
        directory: PathBuf,
        /// Why.
        #[source]
        source: ResumeError,
    },
    /// The filter's passes failed.
    #[error("the hidden-duplication filter could not finish this run")]
    ParalogFilter {
        /// What failed.
        #[source]
        source: ParalogFilterError,
    },
}

/// Runs `pvc resume-calling`.
///
/// # Errors
///
/// If the output or the log cannot be written, if the directory holds no resume point this build
/// can read, or if scoring or writing fails.
pub fn run_resume_calling(args: &ResumeCallingArgs) -> Result<(), ResumeCallingCliError> {
    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }
    calling_run::refuse_an_output_that_cannot_be_written(&args.output)?;
    if !is_standard_output(&args.output) && args.output.exists() && !args.force {
        return Err(ResumeCallingCliError::OutputAlreadyThere {
            path: args.output.clone(),
        });
    }
    let started = Instant::now();
    let log = args
        .log
        .clone()
        .unwrap_or_else(|| run_output::beside(&args.output, "pvc_calls.log"));
    let _run_log = run_output::open_log(&log, args.force)?;

    let mut scratch = calling_run::kept_scratch_directory(
        &args.scratch_dir,
        calling_run::ScratchOnFailure::from_keep_flag(args.keep_tmp_on_error),
    )?;
    let not_resumable = |source| ResumeCallingCliError::NotResumable {
        directory: args.scratch_dir.clone(),
        source,
    };
    let point = resume::load(scratch.path()).map_err(not_resumable)?;
    let mut metadata = point.metadata;
    if let Some(decision) = metadata.neighbour_evidence().cloned() {
        let decision = match args.refine_window_bytes {
            Some(bytes) => decision.with_window_bytes(bytes),
            None => decision,
        };
        let decision = match args.refine_window_floor {
            Some(bases) => decision.with_window_floor(bases),
            None => decision,
        };
        let decision = match args.phase_min_path_agreement {
            Some(agreement) => decision.with_path_agreement(agreement),
            None => decision,
        };
        metadata = metadata.the_neighbour_evidence_decision(decision);
    }
    tracing::info!(
        "resuming the run `{}`: {} samples, {} records on its spill, from {}",
        metadata.command_line(),
        run_output::with_thousands(metadata.sample_names().len() as u64),
        run_output::with_thousands(point.spill.entries_written()),
        scratch.path().display()
    );

    let filtered = paralog_filter::fit_score_and_write_the_calls_explaining(
        &point.spill,
        point.histograms,
        &point.inbreeding,
        point.asked_for,
        &args.output,
        metadata,
        point.ploidy,
        None,
    )
    .map_err(|source| ResumeCallingCliError::ParalogFilter { source })?;

    for line in paralog_filter::what_to_tell_the_operator(&filtered) {
        match line.trim_start().strip_prefix("warning: ") {
            Some(warning) => tracing::warn!("{warning}"),
            None => tracing::debug!("report: {line}"),
        }
    }
    calling_run::report_neighbour_evidence(&filtered.finished, &args.output);
    calling_run::report_read_phase(&filtered.finished);
    tracing::info!(
        "done: {} records written to {} in {}",
        run_output::with_thousands(filtered.finished.records_written),
        calling_run::where_the_calls_went(&args.output),
        run_output::elapsed_text(started.elapsed())
    );
    // The spill is closed before its directory goes.
    drop(point.spill);
    scratch.the_run_succeeded();
    Ok(())
}

#[cfg(test)]
mod tests;
