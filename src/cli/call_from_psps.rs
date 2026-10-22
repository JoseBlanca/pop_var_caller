//! `call-from-psps` — psp mode's calling stage at the command line.
//!
//! **The same cohort, the same parameters, the same VCF** (spec `run_streaming.md` §12.3). What
//! differs from `call-from-alignments` is where the evidence comes from: this command opens one
//! stored file per sample instead of one alignment file per sample, and every stage above the
//! observations — the merge, the candidate selection, the genotyping, the VCF — is the body
//! both modes drive. The numbers each run scores with come from one place too
//! ([`calling_run`]), because two copies of that decision would be two
//! places for the modes to drift apart while both kept running.
//!
//! # The ground is the files', not a flag's
//!
//! **There is no `--regions` here** (spec §5.3). A psp records the ground its walk covered, and
//! a cohort whose files disagree about it is refused rather than called over the ground they
//! happen to share — so the analysed regions are read from the headers, and asking for them
//! again on the command line would only let a person contradict the files. To call over less
//! ground than the psps hold, walk less ground.
//!
//! # What this run can say about each sample, and what it cannot
//!
//! **A psp carries no count of what its walk kept or dropped.** Those tallies belong to the read
//! cursor of a walk that ran in another process, so where `call-from-alignments` reports each
//! sample's kept and dropped reads with the reasons, this reports what it *drew*: how many
//! stored loci it read out of each file and how deep they were, both measured as the records
//! went past. It also compares one thing nothing else compares — the read filters each walk
//! applied, which every psp records and the format never checks — and says so where the
//! cohort's files disagree (owner's ruling, 2026-09-04).

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Args;
use thiserror::Error;

use crate::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
use crate::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::cli::calling_run::{self, CallingRunError};
use crate::cli::generate_psps::PSP_FILE_EXTENSION;
use crate::cli::psp_inputs::{PspArgumentRefusal, psps_named};
use crate::cli::run_ground::{self, GroundError};
use crate::cli::run_output::{self, LogFileError};
use crate::fasta::ContigList;
use crate::parameter_estimation::joint::loci::ReferenceDigest;
use crate::read::input::read_groups::ReadGroups;
use crate::read::input::reference::OpenReference;
use crate::reference_info::{
    ReferenceCheck, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::region_typing::GenomeRegions;
use crate::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN;
use crate::run::explain::{ExplainRegions, ParalogExplanations, write_explanations};
use crate::run::paralog_filter::{self, CalledRecordSink, SpillingSink};
use crate::run::report::BoundsTheRunCalledUnder;
use crate::run::{OpenPspCohort, PspVariantCaller, RunError, RunReport, StoredCohortInputs};
use crate::types::InbreedingF;
use crate::vcf::writer::{VcfWriteError, VcfWriter, is_standard_output};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line.
pub const SUBCOMMAND: &str = "call-from-psps";

/// Default target false-discovery rate for the hidden-duplication filter — spec §3.6's.
///
/// **About one record in a hundred of those the filter removes was really a variant.** It was
/// `0` while the filter's scoring and writing passes did not exist, because that default would
/// have parked every record and produced no VCF; both passes are built now, so the spec's value
/// stands and the filter is on unless an operator turns it off with `--hidden-dup-fdr 0`.
const DEFAULT_PARALOG_FDR: f64 = 0.01;

/// Call a cohort of stored psps and write a VCF.
///
/// **`--parameters` and `--defaults` are one group, and exactly one of them is required** — the
/// same shape `call-from-alignments` uses, so that a run naming neither is told both answers.
#[derive(Debug, Args, Clone)]
#[command(group(
    clap::ArgGroup::new("where the numbers come from")
        .required(true)
        .args(["parameters", "defaults"])
))]
pub struct CallFromPspsArgs {
    /// The reference FASTA the psps' samples were aligned to. A `.fai` index is built beside it
    /// if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The repeat catalogue the psps were written with. Defaults to
    /// `<reference>.repeats.parquet`.
    ///
    /// A catalogue other than the psps' is refused.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// A psp file, or a directory whose `.psp` files are all read. Repeat the option.
    ///
    /// The VCF's sample columns follow the order given; a directory's files are taken in name
    /// order.
    #[arg(long = "psp", required = true, num_args = 1..)]
    pub psps: Vec<PathBuf>,

    /// Where to write the VCF; `-` writes it to standard output, uncompressed. A name ending in
    /// `.vcf.gz` or `.vcf.bgz` is written compressed (bgzip); any other is plain text.
    ///
    /// A file appears whole or not at all: it is written as `<output>.tmp` and renamed when the
    /// last record is on disk. Standard output cannot be taken back: a run that fails part-way
    /// has sent part of a VCF, and its exit status says it failed.
    #[arg(long)]
    pub output: PathBuf,

    /// A parameters file, as `pvc estimate-parameters` writes it for these psps.
    ///
    /// It is refused if it was fitted against another reference, or for other samples or read
    /// groups.
    #[arg(long)]
    pub parameters: Option<PathBuf>,

    /// Call with the numbers built into pvc instead of a parameters file.
    ///
    /// These assume no base-quality miscalibration, no contamination and no inbreeding, and use
    /// built-in repeat stutter rates.
    #[arg(long)]
    pub defaults: bool,

    /// Copies of each chromosome in every sample. Defaults to 2, or to the parameters file's.
    ///
    /// With a parameters file it may only repeat the file's ploidy, which the numbers were
    /// fitted at. At most 16.
    #[arg(long)]
    pub ploidy: Option<u8>,

    /// The site quality (QUAL) below which a site is not written. 0 writes every site.
    ///
    /// The default, 1, drops the sites whose artifact penalties outweigh their evidence, which
    /// removes most false calls at little cost in true ones.
    #[arg(
        long,
        default_value_t = crate::run::records::DEFAULT_MIN_SITE_QUALITY,
        value_parser = crate::cli::calling_run::parse_min_site_quality
    )]
    pub min_site_quality: f32,

    /// Threads to use; 0 uses every core. The output is the same at any number.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Replace the VCF and the log if they are already there. Without it, a run refuses before
    /// reading anything.
    #[arg(long)]
    pub force: bool,

    /// Where to write the run's log. Defaults to `pvc_calls.log` beside the VCF.
    ///
    /// The log holds the inputs by name, every setting's value, and the run's full report.
    #[arg(long)]
    pub log: Option<PathBuf>,

    /// The directory for the run's scratch files. Defaults to the VCF's directory, or the
    /// working directory when the VCF goes to standard output.
    ///
    /// The run makes a directory of its own inside it, `pvc_tmp_<random>`, and removes it when
    /// it ends, whether it succeeded or failed; a run that is killed, or stops on an internal
    /// error, leaves it behind. While the hidden-dup filter is on, every call waits there until
    /// the whole genome is called: several times the compressed VCF's size.
    #[arg(long)]
    pub tmp_dir: Option<PathBuf>,

    /// A locus where any sample has a read group deeper than this is not called, for any sample.
    ///
    /// Such depth, far above a sample's usual coverage, is almost always a collapsed repeat:
    /// many copies of a sequence mapped onto one copy in the reference. A psp holding more
    /// reads of a read group at a position is also thinned to this many as it is read. The
    /// default, 1,000, is about three times the depth of a 300× sample.
    #[arg(
        long,
        default_value_t = crate::run::depth_ceiling::DEFAULT_MAX_READ_GROUP_DEPTH,
        value_parser = clap::value_parser!(u32).range(1..),
        help_heading = "Advanced"
    )]
    pub max_depth: u32,

    /// The most alleles a SNP/indel locus is called over, the reference among them.
    ///
    /// Where more segregate, the worst supported are left out. Repeat (STR) loci are not
    /// affected: they have their own limit of 32.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        help_heading = "Advanced"
    )]
    pub max_generic_alleles: u16,

    /// The widest locus assembled, in reference bases — and so the longest deletion called.
    ///
    /// The run's log says how many loci were too wide and how long they were. Long reads want a
    /// larger value than short ones.
    #[arg(long, default_value_t = DEFAULT_MAX_COHORT_LOCUS_SPAN, help_heading = "Advanced")]
    pub max_locus_span: u32,

    /// The false-discovery rate allowed among the calls the hidden-dup filter removes. 0 turns
    /// the filter off.
    ///
    /// The filter removes calls better explained by a hidden duplication — two copies of a
    /// sequence collapsed onto one in the reference — than by a variant. While it is on, the
    /// run needs scratch space in --tmp-dir of several times the compressed VCF's size.
    #[arg(long, default_value_t = DEFAULT_PARALOG_FDR, help_heading = "Advanced")]
    pub hidden_dup_fdr: f64,

    /// Keep the calls the hidden-dup filter removes, marked `hiddenDup` in the FILTER column,
    /// instead of dropping them.
    #[arg(long, help_heading = "Advanced")]
    pub hidden_dup_tag: bool,

    /// The width, in bases, of the stretch of genome over which the hidden-dup filter estimates
    /// how common hidden duplications are. 0 uses one genome-wide rate.
    ///
    /// Each call's rate comes from the calls within half the window on either side, pulled toward
    /// the genome-wide rate where few calls are near. A research option: its default may change.
    #[arg(
        long,
        value_name = "BASES",
        default_value_t = paralog_filter::DEFAULT_HIDDEN_DUP_WINDOW_BASES,
        help_heading = "Advanced"
    )]
    pub hidden_dup_window: u64,

    /// The most reads of one read group kept at one position. Defaults to `--max-depth`.
    ///
    /// A psp holding more is thinned to this as it is read: each allele keeps its share of the
    /// reads, chosen by a hash so neighbouring positions keep the same reads. At equal values a
    /// thinned position is always dropped by `--max-depth`.
    #[arg(
        long,
        value_parser = clap::value_parser!(u32).range(1..),
        help_heading = "Developer"
    )]
    pub max_reads_per_position: Option<u32>,

    /// How much reference one round of locus building covers, in bases. Chosen from the cohort
    /// size when not given; it trades memory for speed and changes no call.
    #[arg(long, help_heading = "Developer")]
    pub locus_builder_regions_len: Option<u32>,

    /// How much of the psps to read ahead into the operating system's file cache, in bytes. 0
    /// turns it off. Changes no call.
    ///
    /// On a spinning disk with thousands of psps, more is faster.
    #[arg(
        long,
        default_value_t = crate::run::psp_prefetch::DEFAULT_PSP_PREFETCH_BUDGET_BYTES,
        help_heading = "Developer"
    )]
    pub psp_prefetch_bytes: u64,

    /// Explain every decision at the loci overlapping a BED, written as a table: candidates,
    /// genotype likelihoods, GT, GQ, QUAL and its penalties, and the record's fate. Changes no
    /// call.
    #[arg(
        long,
        num_args = 2,
        value_names = ["REGIONS_BED", "OUTPUT_TSV"],
        help_heading = "Developer"
    )]
    pub explain_loci: Option<Vec<PathBuf>>,
}

impl CallFromPspsArgs {
    /// The read cap stored records are thinned to: `--max-reads-per-position` where it was given,
    /// and `--max-depth` otherwise (`doc/devel/ng/spec/command_line.md` §5.4).
    #[must_use]
    pub fn read_cap(&self) -> u32 {
        self.max_reads_per_position.unwrap_or(self.max_depth)
    }
}

/// Everything that can stop a run over stored files, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CallFromPspsCliError {
    /// The `--explain-loci` regions could not be read.
    #[error("reading the --explain-loci regions in {path}")]
    ExplainRegions {
        /// The BED.
        path: PathBuf,
        /// What the BED reader hit.
        #[source]
        source: crate::regions::BedError,
    },

    /// The `--explain-loci` explanation could not be written.
    #[error("writing the --explain-loci explanation to {path}")]
    ExplanationNotWritten {
        /// The TSV.
        path: PathBuf,
        /// What the write hit.
        #[source]
        source: std::io::Error,
    },

    /// The hidden-duplication filter's target is not a false-discovery rate.
    #[error("--hidden-dup-fdr is not a false-discovery rate")]
    ParalogTargetIsNotAFraction {
        /// What the type said was wrong with it.
        #[source]
        source: crate::run::paralog_filter::NotATargetFdr,
    },

    /// The hidden-duplication filter could not finish the run.
    #[error("the hidden-duplication filter could not finish this run")]
    ParalogFilter {
        /// What the filter said.
        #[source]
        source: crate::run::paralog_filter::ParalogFilterError,
    },

    /// The calls could not be written, or the records could not be parked for the filter.
    #[error("the calls could not be written to {}", path.display())]
    CallsNotWritten {
        /// Where they were going.
        path: PathBuf,
        /// What pass one said.
        #[source]
        source: crate::run::paralog_filter::PassOneError,
    },
    /// The reference could not be read.
    #[error("reading the reference {}", path.display())]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The reference's FASTA could not be verified against its index.
    #[error("verifying the reference {} against its index", path.display())]
    ReferenceVerification {
        /// The FASTA.
        path: PathBuf,
        /// What the verification said.
        #[source]
        source: ReferenceInfoError,
    },

    /// A `--psp` naming a directory that could not be listed.
    #[error("the psps in {} could not be listed", path.display())]
    PspDirectory {
        /// The directory.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },

    /// A `--psp` naming a directory with no psp in it.
    ///
    /// **Refused rather than skipped**, because a directory typed by mistake and a directory
    /// whose walk has not finished look the same from here, and calling a cohort short of a
    /// sample is not something to do quietly.
    #[error(
        "--psp {} holds no .{PSP_FILE_EXTENSION} file; name the files themselves, or point at \
         the --output-dir a generate-psps run wrote",
        path.display()
    )]
    NoPspsInDirectory {
        /// The directory.
        path: PathBuf,
    },

    /// The ground or the catalog this run would call over could not be worked out.
    ///
    /// **Transparent**, so the sentence is the shared one every mode renders.
    #[error(transparent)]
    Ground(#[from] GroundError),

    /// A number, a path or a parameters file this run cannot call under.
    #[error(transparent)]
    Calling(#[from] CallingRunError),

    /// The output could not be created, written or renamed into place.
    #[error("the output {} could not be written", path.display())]
    Output {
        /// The VCF.
        path: PathBuf,
        /// What the writer said.
        #[source]
        source: VcfWriteError,
    },

    /// The VCF is already there, and the run was not asked to replace it.
    #[error(
        "{} is already there — give --force to replace it, or another --output",
        path.display()
    )]
    OutputAlreadyThere {
        /// The VCF.
        path: PathBuf,
    },

    /// The run's log could not be opened.
    #[error(transparent)]
    Log(#[from] LogFileError),

    /// The run itself was refused, or stopped.
    #[error("the run stopped")]
    Run {
        /// What stopped it.
        #[source]
        source: RunError,
    },
}

/// **The psps this run calls, with every directory expanded** — the rule every psp-taking command
/// shares ([`psps_named`]), with its refusals dressed in this command's words.
fn psps_named_by(args: &CallFromPspsArgs) -> Result<Vec<PathBuf>, CallFromPspsCliError> {
    psps_named(&args.psps).map_err(|refusal| match refusal {
        PspArgumentRefusal::Unlistable { path, source } => {
            CallFromPspsCliError::PspDirectory { path, source }
        }
        PspArgumentRefusal::Empty { path } => CallFromPspsCliError::NoPspsInDirectory { path },
    })
}

/// Call `--psp`'s cohort and write `--output`. Prints a summary when it finishes.
///
/// # Errors
///
/// Every way a run can be refused or can stop — see [`CallFromPspsCliError`]. **In order, and
/// all of it before a block is decoded**: the ploidy, the allele cap and the two output paths;
/// then the psps named, because the round width is chosen from how many there are; then that
/// width and the locus-span bound; then the reference; then every check of §6.2, across the
/// cohort and against this run.
pub fn run_call_from_psps(args: &CallFromPspsArgs) -> Result<(), CallFromPspsCliError> {
    // **Everything a person typed is judged before a byte is read**, the same order the other
    // two commands use.
    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    // **A target that is not a fraction is refused before anything else looks at it.** `clap`
    // parses any `f64`, so `7`, `inf` and `nan` all arrive here, and the type is what says which
    // of those is a false-discovery rate.
    let asked_for = paralog_filter::WhatTheOperatorAskedFor::from_the_flags(
        args.hidden_dup_fdr,
        args.hidden_dup_tag,
    )
    .map_err(|source| CallFromPspsCliError::ParalogTargetIsNotAFraction { source })?
    .map(|asked| asked.with_window(args.hidden_dup_window));

    let asked_ploidy = calling_run::ploidy_asked_for(args.ploidy)?;
    let candidate_selection = calling_run::candidate_selection_for(args.max_generic_alleles)?;
    calling_run::refuse_an_output_that_cannot_be_written(&args.output)?;
    calling_run::refuse_a_tmp_dir_that_is_not_a_directory(args.tmp_dir.as_deref())?;
    if !is_standard_output(&args.output) && args.output.exists() && !args.force {
        return Err(CallFromPspsCliError::OutputAlreadyThere {
            path: args.output.clone(),
        });
    }
    let paths = psps_named_by(args)?;
    // **The round width comes from the cohort's size**, and the cohort is the psps — so this
    // waits for the directories to be expanded, where direct mode can count its `--alignment`
    // flags straight away.
    let merge_parameters = calling_run::merge_parameters_for(
        args.max_locus_span,
        args.locus_builder_regions_len,
        paths.len(),
    )?;
    // **The log opens once the command line has been judged**, so a typing mistake leaves no log
    // behind, and before the reference is read, so everything that takes time is in it.
    let started = Instant::now();
    let log = args
        .log
        .clone()
        .unwrap_or_else(|| run_output::beside(&args.output, "pvc_calls.log"));
    let run_log = run_output::open_log(&log, args.force)?;

    // **The FASTA is read to the end before a psp is opened**: every file's contig table is
    // compared against the reference's per-contig checksums, and a reference whose verification
    // has not been joined carries none.
    let cache = std::sync::Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        args.reference.clone(),
        ReferenceCheck::VerifyAgainstIndex,
    )
    .map_err(|source| CallFromPspsCliError::Reference {
        path: args.reference.clone(),
        source,
    })?;
    let with_checksums = match verify {
        Some(handle) => {
            handle
                .join()
                .map_err(|source| CallFromPspsCliError::ReferenceVerification {
                    path: args.reference.clone(),
                    source,
                })?
        }
        None => std::sync::Arc::clone(&info),
    };
    let contigs: ContigList = info.contig_list();
    // **Read before any work**, so a BED that will not read stops the run at once.
    let explain = match args.explain_loci.as_deref() {
        Some([bed, tsv]) => Some((
            ExplainRegions::from_bed_path(bed, &contigs).map_err(|source| {
                CallFromPspsCliError::ExplainRegions {
                    path: bed.clone(),
                    source,
                }
            })?,
            tsv.clone(),
        )),
        Some(_) => unreachable!("clap takes exactly two values for --explain-loci"),
        None => None,
    };
    let reference = OpenReference::new(info);

    // **The cohort is opened before the segmentation is built, because the cohort is what says
    // what ground there is** (spec §5.3): the analysed regions come from the headers, and the
    // files have already been forced to agree about them.
    let cohort =
        OpenPspCohort::open(&paths).map_err(|source| CallFromPspsCliError::Run { source })?;
    let analysed = cohort.analysed_regions().clone();
    // **What counts as a repeat is the psps' own**, recorded by the walk that wrote them: a run
    // over stored loci has to cut the genome the way that walk did, so it has no option that
    // could say otherwise (`doc/devel/ng/spec/command_line.md` §4.4).
    // **And so is whether there were repeats to sort at all**: psps written with
    // `--generic-path-only` read no catalog, and this run reads none either.
    let segmentation = match run_ground::sorting_of_the_psps(cohort.segmentation_inputs()) {
        Some((setting, walked_with)) => {
            let criteria = walked_with.repeat_tract_criteria.clone();
            let catalog_path =
                run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
            run_ground::segments_cut_with(
                &catalog_path,
                &args.reference,
                setting,
                &criteria,
                run_ground::CriteriaSource::ThePspHeaders,
                &analysed,
                &with_checksums,
            )?
        }
        None => {
            run_ground::refuse_a_catalog_the_psps_never_read(args.catalog.as_deref())?;
            run_ground::segments_without_a_catalog(&analysed)
        }
    };
    let satellite_bases_on_the_generic_path = segmentation.satellite_bases_on_the_generic_path();
    say_what_this_run_does(args, &paths, cohort.read_groups(), &analysed, &log);

    let numbers = calling_run::run_parameters(
        args.parameters.as_deref(),
        args.ploidy,
        cohort.read_groups(),
        &with_checksums,
        asked_ploidy,
        segmentation.inputs().repeat_tract_criteria(),
    )?;
    let ploidy = numbers.parameters.ploidy();
    // **Captured before the parameters move into the caller.** The filter reads one coefficient
    // a sample, in the run's sample order, and the caller consumes the parameters it lives on.
    let inbreeding_by_sample: Vec<InbreedingF> = numbers
        .parameters
        .inbreeding_coefficient_by_sample()
        .to_vec();

    let digest = ReferenceDigest::of(&with_checksums)
        .map_err(|source| CallingRunError::ReferenceNotDigested { source })?;
    let parameters_file = numbers.parameters_file(
        cohort.read_groups(),
        &digest,
        segmentation.inputs().repeat_tract_criteria(),
    );
    // **The numbers go into the log before a read is decoded**, so a run that stops part-way has
    // already recorded what it was scoring with (`doc/devel/ng/spec/command_line.md` §7.3).
    run_log.record_parameters(&parameters_file.to_toml());

    let caller = PspVariantCaller::open(
        cohort,
        StoredCohortInputs {
            reference: &reference,
            reference_with_checksums: &with_checksums,
        },
        segmentation,
        numbers.parameters,
        calling_run::calling_loop_config_for_this_run()?,
        candidate_selection,
        merge_parameters,
    )
    .map_err(|source| CallFromPspsCliError::Run { source })?
    .with_psp_prefetch_budget(args.psp_prefetch_bytes)
    .with_max_reads_per_position(args.read_cap())
    .with_depth_ceiling(crate::run::depth_ceiling::DepthCeiling::at(args.max_depth));
    let caller = match &explain {
        Some((regions, _)) => caller.with_explain_loci(regions.clone()),
        None => caller,
    };
    let caller = caller.with_min_site_quality(
        crate::types::Phred::try_new(args.min_site_quality)
            .expect("`parse_min_site_quality` admits only valid qualities"),
    );
    let sample_names: Vec<String> = caller.sample_names().map(str::to_owned).collect();

    let read_groups = caller.read_groups().clone();
    let metadata = calling_run::header_for(
        &log,
        args.parameters.as_deref(),
        &args.reference,
        &contigs,
        &with_checksums,
        caller.sample_names().map(str::to_owned).collect(),
    )?;
    // **The sink is where the filter's flag lands.** With the filter off, this hands each record
    // straight to the writer — the call the run made before the filter existed, which is what
    // makes `--hidden-dup-fdr 0` byte-identical by construction rather than by hope. With it on, the
    // records are parked and **the VCF is not opened at all**: no record's fate is known until
    // every record has been scored, so a file opened now would stand beside a run that has
    // decided nothing.
    // **The run's own scratch directory, made only where something will be parked in it**, and
    // removed when `scratch` is dropped — at the end of this function, on its success or on any
    // error that returns through it.
    let scratch = match asked_for {
        Some(_) => Some(calling_run::scratch_directory_for(
            &args.output,
            args.tmp_dir.as_deref(),
        )?),
        None => None,
    };
    let mut sink = if let Some(scratch) = &scratch {
        // The spill's lines are encoded against the same contigs the header names, taken from
        // the metadata itself rather than re-derived.
        CalledRecordSink::ParkedOnTheSpill(SpillingSink::in_directory(
            scratch.path(),
            metadata.contigs().to_vec(),
            ploidy,
        ))
    } else {
        CalledRecordSink::StraightToTheVcf(
            VcfWriter::create(&args.output, metadata.clone(), ploidy).map_err(|source| {
                CallFromPspsCliError::Output {
                    path: args.output.clone(),
                    source,
                }
            })?,
        )
    };

    // **The record's spill entry is encoded on the thread that called its locus** — the preparer
    // is the sink's, taken once here — and the sink only appends it, in genome order.
    let preparer = sink.preparer();
    let (mut calling, mut stored) = caller
        .call_cohort_preparing_each_record(
            &SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior),
            &|record, window_coverage| preparer.prepare(record, window_coverage),
            &mut |prepared, _| sink.accept_prepared(prepared),
        )
        .map_err(|source| CallFromPspsCliError::Run { source })?;

    // **The spill has to stay alive until pass three has read it**, so the sink is taken apart
    // rather than consumed: `finish` would drop the file the next two passes need.
    let parked = match sink {
        CalledRecordSink::StraightToTheVcf(writer) => {
            writer
                .finish()
                .map_err(|source| CallFromPspsCliError::CallsNotWritten {
                    path: args.output.clone(),
                    source: crate::run::paralog_filter::PassOneError::Vcf(source),
                })?;
            None
        }
        CalledRecordSink::ParkedOnTheSpill(mut spilling) => {
            spilling
                .finish_parking()
                .map_err(|source| CallFromPspsCliError::CallsNotWritten {
                    path: args.output.clone(),
                    source: crate::run::paralog_filter::PassOneError::Spill(source),
                })?;
            Some(spilling)
        }
    };

    // **Passes two and three, through the one function both subcommands call**, so direct mode
    // and psp mode cannot drift into filtering differently (spec §1.1 goal 2).
    // The filter is asked for its verdicts on the records the explained loci wrote.
    let mut paralog_explanations = explain
        .as_ref()
        .map(|_| ParalogExplanations::wanting(&calling.explanations));
    let filtered = match &parked {
        Some(spilling) => Some(
            paralog_filter::fit_score_and_write_the_calls_explaining(
                spilling.spill(),
                // **Taken, not cloned.** `ParalogScoringContext::new` consumes the histograms so
                // each sample's bins are freed as its model is fitted; cloning here would defeat
                // that and hold a second copy of every sample's — about 80 kB each. Nothing else
                // reads the field afterwards.
                std::mem::take(&mut stored.window_coverage_histograms),
                &inbreeding_by_sample,
                // PANIC-FREE: the sink is the spill only where `asked_for` is `Some`, and this
                // arm is reached only through that sink.
                asked_for.expect("the records are parked only when a filter was asked for"),
                &args.output,
                metadata,
                ploidy,
                paralog_explanations.as_mut(),
            )
            .map_err(|source| CallFromPspsCliError::ParalogFilter { source })?,
        ),
        None => None,
    };

    if let Some(verdicts) = paralog_explanations {
        verdicts.add_to(&mut calling.explanations);
    }
    if let Some((_, tsv)) = &explain {
        write_explanations(
            tsv,
            &mut std::mem::take(&mut calling.explanations),
            &contigs,
            &sample_names,
        )
        .map_err(|source| CallFromPspsCliError::ExplanationNotWritten {
            path: tsv.clone(),
            source,
        })?;
    }

    // **What the filter did goes at the end of the report**, and is empty on a run with the
    // filter off — so an off run prints what it printed before the filter existed.
    let filter_lines = filtered
        .as_ref()
        .map(paralog_filter::what_to_tell_the_operator)
        .unwrap_or_default();
    calling_run::log_the_report(
        &args.output,
        &RunReport::of_a_stored_cohort(
            &calling,
            &stored,
            &contigs,
            &read_groups,
            &parameters_file,
            &analysed,
            BoundsTheRunCalledUnder {
                max_cohort_locus_span: merge_parameters.max_cohort_locus_span.get(),
                max_candidate_alleles: candidate_selection.max_candidate_alleles.get(),
            },
        )
        .with_satellite_bases_on_the_generic_path(satellite_bases_on_the_generic_path),
        &filter_lines,
    );
    let records = filtered
        .as_ref()
        .map_or(calling.records_written, |run| run.did.written);
    tracing::info!(
        "done: {} records written to {} in {}",
        run_output::with_thousands(records),
        calling_run::where_the_calls_went(&args.output),
        run_output::elapsed_text(started.elapsed())
    );
    tracing::info!("log: {}", log.display());
    Ok(())
}

/// **The run's first words**: on the terminal, what it is doing and on how much; in the log, the
/// inputs by name (`doc/devel/ng/spec/command_line.md` §7.2, §7.3).
fn say_what_this_run_does(
    args: &CallFromPspsArgs,
    psps: &[PathBuf],
    read_groups: &ReadGroups,
    analysed: &GenomeRegions,
    log: &Path,
) {
    let samples: BTreeSet<&str> = read_groups
        .iter()
        .map(|(_, group)| &*group.sample)
        .collect();
    let bases: u64 = analysed.iter().map(|region| region.len()).sum();
    let threads = rayon::current_num_threads();
    tracing::info!(
        "pvc call-from-psps: calling {} sample{} from psp files over {}, {threads} thread{}",
        samples.len(),
        if samples.len() == 1 { "" } else { "s" },
        run_output::size_of_ground(bases),
        if threads == 1 { "" } else { "s" },
    );
    tracing::info!("log: {}", log.display());

    tracing::debug!("input: reference {}", args.reference.display());
    let catalog = run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
    tracing::debug!("input: repeat catalogue {}", catalog.display());
    for psp in psps {
        tracing::debug!("input: psp {}", psp.display());
    }
    for (_, group) in read_groups.iter() {
        tracing::debug!(
            "input: sample {}, read group {}, from {}",
            group.sample,
            group.id,
            group.file.display()
        );
    }
    tracing::debug!("input: ground the psps were written over, {bases} bases");
    tracing::debug!(
        "output: {}",
        calling_run::where_the_calls_went(&args.output)
    );
}
