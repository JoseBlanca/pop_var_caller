//! The `call-from-alignments` subcommand: **call a cohort of alignment files and write a VCF**,
//! in one process, with nothing written to disk in between — and a log holding, among the rest,
//! the model parameters the run scored with.
//!
//! This is ng's direct mode (`doc/devel/ng/spec/run_streaming.md` §2 and §5.1) at the command
//! line. Every sample's file is opened and held for the whole run, every sample's reads are
//! walked at one shared frontier, each cohort locus is genotyped where it is built, and each
//! record is written as it is finished. **Nothing accumulates but the file's bytes** — no
//! per-sample artefact, no intermediate table.
//!
//! # What it needs, and what it decides for itself
//!
//! In: the reference, the tandem-repeat catalog built beside it, one alignment file per
//! sample, the stretch of genome to call over, and the model's numbers — either a parameters
//! file or `--defaults`.
//!
//! Out: **two files**. The VCF, and the run's log (`pvc_calls.log` beside it unless `--log` says
//! otherwise), which holds every number the run scored with as a parameters file between two
//! marker lines (`doc/devel/ng/spec/command_line.md` §7.3). A run is then reproducible from its
//! own output, a run that guessed its numbers is auditable in the same form as one that measured
//! them, and somebody who wants to change one number copies the block out and starts from a file
//! rather than from a document.
//!
//! **The catalog is not optional and it is not a repeat caller's input.** It is what says where
//! the repeat tracts are, so the run can route them; every locus is called down the SNP/indel
//! path today and a tract is counted as ground this caller has not built yet, which the run's
//! own counts say out loud.
//!
//! # Where this run is short, and it says so rather than being wrong
//!
//! **Repeat tracts and STR bundles are called.** A tract goes down its own path
//! — `select_ssr`, the stutter emission, and the same frequency loop — and its record carries
//! `STR`, `RU`, `PERIOD` and each called allele's `REPCN`
//! (`doc/devel/ng/impl_plan/calling_loop_ssr.md` Milestone C). An STR bundle — repeats too close
//! together for any to have clean flanks — is one locus over its hull, called with stutter in
//! each of its runs (`str_bundle_loci.md` §7); its record carries per-run `RU` and `PERIOD` and
//! each called allele's run counts in `REPCNB`, and the summary says how many bundles were
//! called.
//!
//! **The decode is parallel; the rest is one thread (Milestone E1, 2026-09-01).** Decoding
//! reads is 88% of a calling run at 63 samples, and the run sweeps every sample's reader
//! concurrently on the rayon pool while assembling and genotyping stay on one thread — the
//! output is identical at every thread count. Nothing is printed while it works.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use clap::Args;
use thiserror::Error;

use crate::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
use crate::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::cli::calling_run;
use crate::cli::run_ground::{self, GroundError};
use crate::cli::run_output::{self, LogFileError};
use crate::fasta::ContigList;
use crate::locus_generation::pileup::PileupGeneratorConfig;
use crate::parameter_estimation::joint::loci::ReferenceDigest;
use crate::read::ReadFilterConfig;
use crate::read::input::read_groups::{ReadGroupError, ReadGroups, build_read_groups};
use crate::read::input::reference::OpenReference;
use crate::reference_info::{
    ReferenceCheck, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::region_typing::{DEFAULT_MAX_STR_BUNDLE_LEN, DEFAULT_MAX_STR_LEN, GenomeRegions};
use crate::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN;
use crate::run::explain::{ExplainRegions, ParalogExplanations, write_explanations};
use crate::run::paralog_filter::{self, CalledRecordSink, SpillingSink};
use crate::run::report::BoundsTheRunCalledUnder;
use crate::run::{AlignedFilesVariantCaller, AlignmentInputs, RunError, RunReport};
use crate::types::InbreedingF;
use crate::types::MAX_MOTIF_LEN;
use crate::vcf::writer::{VcfWriteError, VcfWriter, is_standard_output};

/// Default target false-discovery rate for the hidden-duplication filter — spec §3.6's.
///
/// **About one record in a hundred of those the filter removes was really a variant.** It was
/// `0` while the filter's scoring and writing passes did not exist, because that default would
/// have parked every record and produced no VCF; both passes are built now, so the spec's value
/// stands and the filter is on unless an operator turns it off with `--hidden-dup-fdr 0`.
const DEFAULT_PARALOG_FDR: f64 = 0.01;

/// `call-from-alignments` arguments.
///
/// **`--parameters` and `--defaults` are one group, and exactly one of them is required.** A
/// group rather than a pair of conflicting flags so that a run naming neither is told both
/// answers rather than only the first: the message reads *one of --parameters, --defaults*.
#[derive(Debug, Args, Clone)]
#[command(group(
    clap::ArgGroup::new("where the numbers come from")
        .required(true)
        .args(["parameters", "defaults"])
))]
pub struct CallFromAlignmentsArgs {
    /// The reference FASTA every alignment file was made against. A `.fai` index is built
    /// beside it if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The repeat catalogue. Defaults to `<reference>.repeats.parquet`, where
    /// `pvc repeat-catalog` writes it.
    ///
    /// A catalogue built on another reference is refused.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// A sample's BAM or CRAM file, indexed. Repeat the option, once per file.
    ///
    /// Samples are named by their files' `@RG SM` tag, and are the VCF's sample columns in the
    /// order given. Two files with the same `SM` are one sample.
    #[arg(long = "alignment", required = true, num_args = 1..)]
    pub alignments: Vec<PathBuf>,

    /// Where to write the VCF; `-` writes it to standard output, uncompressed. A name ending in
    /// `.vcf.gz` or `.vcf.bgz` is written compressed (bgzip); any other is plain text.
    ///
    /// A file appears whole or not at all: it is written as `<output>.tmp` and renamed when the
    /// last record is on disk. Standard output cannot be taken back: a run that fails part-way
    /// has sent part of a VCF, and its exit status says it failed.
    #[arg(long)]
    pub output: PathBuf,

    /// A BED file of the regions to call. Without it, the whole genome.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// A parameters file, as `pvc estimate-parameters` writes it for this cohort.
    ///
    /// It is refused if it was fitted against another reference, or for other samples or read
    /// groups.
    #[arg(long)]
    pub parameters: Option<PathBuf>,

    /// Call with the numbers built into pvc instead of a parameters file.
    ///
    /// These assume no base-quality miscalibration, no contamination and no inbreeding, and use
    /// built-in repeat stutter rates. Fitted numbers come from psp mode, through
    /// `pvc estimate-parameters`.
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
    /// many copies of a sequence mapped onto one copy in the reference. Reading also keeps no
    /// more than this many reads of a read group at one position. The default, 1,000, is about
    /// three times the depth of a 300× sample.
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

    /// Build a `.csi` (BAM) or `.crai` (CRAM) index beside any alignment file that has none.
    #[arg(long, help_heading = "Advanced")]
    pub build_index_if_missing: bool,

    /// Call satellites with the SNP/indel model instead of leaving them uncalled. A satellite is
    /// a repeat longer than --max-str-len, or a cluster of repeats longer than
    /// --max-str-bundle-len.
    ///
    /// What this adds on satellites is mostly SNPs and short indels: the satellite's own length
    /// variation mostly stays uncalled, because reads that disagree about its length make loci
    /// wider than --max-locus-span, which are not called. The run report says how many satellite
    /// bases were handed over and how many loci were too wide.
    #[arg(long, help_heading = "Advanced")]
    pub call_satellites: bool,

    /// The fewest repeat-unit copies a stretch needs to be called as a repeat: six
    /// comma-separated numbers, for repeat units of 1 to 6 bases.
    ///
    /// Below these, a stretch is ordinary sequence and is called as SNPs and indels.
    #[arg(
        long,
        value_parser = crate::cli::parsers::parse_min_copies,
        default_value = "8,6,6,6,5,4",
        help_heading = "What counts as a repeat"
    )]
    pub min_copies: MinCopies,

    /// The shortest repeat unit called as a repeat, in bases. 1 includes homopolymers.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub min_period: u8,

    /// The longest repeat unit called as a repeat, in bases. At most 6.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub max_period: u8,

    /// A repeat longer than this, in bases, is a satellite and is not called.
    ///
    /// With 150 bp reads, a read spans a repeat plus a flank on each side only up to about
    /// 90 bp, so longer repeats cannot be genotyped from reads that cross them.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_len: u64,

    /// An STR bundle — repeats within 15 bases of each other, called as one locus — longer
    /// than this, in bases, is a satellite and is not called. 0 calls no bundle.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_BUNDLE_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_bundle_len: u64,

    /// How much of a repeat must match a perfect run of its unit, from 0 to 1. Below this it is
    /// ordinary sequence.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        value_parser = crate::cli::parsers::parse_min_purity,
        help_heading = "What counts as a repeat"
    )]
    pub min_purity: f32,

    /// The most reads of one read group kept at one position. Defaults to `--max-depth`.
    ///
    /// Set apart from `--max-depth`, a position above this is thinned to it, by a hash of the
    /// read names so neighbouring positions keep the same reads. At equal values a thinned
    /// position is always dropped by `--max-depth`.
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

impl CallFromAlignmentsArgs {
    /// The read cap the walk applies: `--max-reads-per-position` where it was given, and
    /// `--max-depth` otherwise (`doc/devel/ng/spec/command_line.md` §5.4).
    #[must_use]
    pub fn read_cap(&self) -> u32 {
        self.max_reads_per_position.unwrap_or(self.max_depth)
    }
}

/// Everything that can stop a run, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CallFromAlignmentsCliError {
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
    ///
    /// **Named apart from the read itself** because the read succeeds and the verification runs
    /// on a second thread: a run reaches this after its contig table is already in hand.
    #[error("verifying the reference {} against its index", path.display())]
    ReferenceVerification {
        /// The FASTA.
        path: PathBuf,
        /// What the verification said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The ground this run would walk could not be worked out — the reference's contigs, the
    /// BED, the catalog, or what this run counts as a repeat.
    ///
    /// **Transparent**, so the sentence a person reads is the one
    /// [`run_ground`] writes: these refusals belong to
    /// every mode and are shared with `generate-psps`, and a wrapper of its own here would
    /// make the same mistake render differently depending on which command was typed.
    #[error(transparent)]
    Ground(#[from] GroundError),

    /// A number, a path or a parameters file this run cannot call under.
    ///
    /// **Transparent**, for [`Ground`](Self::Ground)'s reason: these refusals belong to every
    /// calling run and are shared with `call-from-psps`, so the same mistake must not read
    /// differently depending on which command was typed.
    #[error(transparent)]
    Calling(#[from] calling_run::CallingRunError),

    /// The alignment files' read groups could not be read.
    #[error("the read groups of this run's alignment files could not be read")]
    ReadGroups {
        /// What the reader said.
        #[source]
        source: ReadGroupError,
    },

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

    /// The run itself stopped.
    #[error("the run stopped")]
    Run {
        /// What stopped it.
        #[source]
        source: RunError,
    },
}

/// This command's flags, as the shared ground assembly asks for them.
fn ground_request(args: &CallFromAlignmentsArgs) -> run_ground::GroundRequest<'_> {
    run_ground::GroundRequest {
        reference: &args.reference,
        catalog: args.catalog.as_deref(),
        regions: args.regions.as_deref(),
        routing: run_ground::RepeatRouting {
            min_copies: args.min_copies,
            min_period: args.min_period,
            max_period: args.max_period,
            max_str_len: args.max_str_len,
            max_str_bundle_len: args.max_str_bundle_len,
            min_purity: args.min_purity,
        },
        setting: crate::run::RepeatGroundSetting::from_options(args.call_satellites, false),
    }
}

/// Call `--alignment`'s cohort over `--regions` and write `--output`. Logs the run report when it
/// finishes.
///
/// # Errors
///
/// Every way a run can be refused or can stop — see [`CallFromAlignmentsCliError`].
///
/// **No file that could be mistaken for a finished one is left behind**, which is a weaker
/// claim than *nothing is left behind* and the accurate one: the VCF's bytes go to
/// `<output>.tmp` and are renamed into place only once the last record is on disk, so a run that
/// stops leaves its output path untouched — and leaves the `.tmp` beside it. Production's writer
/// removes it on an abort and ng's has no such path yet, so clearing it is the operator's.
pub fn run_call_from_alignments(
    args: &CallFromAlignmentsArgs,
) -> Result<(), CallFromAlignmentsCliError> {
    // **Everything a person typed is judged before a byte is read.** Reading the reference and
    // opening a cohort of CRAMs is minutes; a number or a path that was never going to work
    // should cost none of them.
    if args.threads > 0 {
        // A failure here means a pool is already built, which in this binary means a second
        // call in one process — the calls are unaffected, so it is not worth an error. Same
        // reasoning, and the same line, as `estimate-contamination`.
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
    .map_err(|source| CallFromAlignmentsCliError::ParalogTargetIsNotAFraction { source })?;

    let asked_ploidy = calling_run::ploidy_asked_for(args.ploidy)?;
    let merge_parameters = calling_run::merge_parameters_for(
        args.max_locus_span,
        args.locus_builder_regions_len,
        args.alignments.len(),
    )?;
    let candidate_selection = calling_run::candidate_selection_for(args.max_generic_alleles)?;
    calling_run::refuse_an_output_that_cannot_be_written(&args.output)?;
    calling_run::refuse_a_tmp_dir_that_is_not_a_directory(args.tmp_dir.as_deref())?;
    if !is_standard_output(&args.output) && args.output.exists() && !args.force {
        return Err(CallFromAlignmentsCliError::OutputAlreadyThere {
            path: args.output.clone(),
        });
    }
    // **The log opens once the command line has been judged**, so a typing mistake leaves no log
    // behind, and before the reference is read, so everything that takes time is in it.
    let started = Instant::now();
    let log = args
        .log
        .clone()
        .unwrap_or_else(|| run_output::beside(&args.output, "pvc_calls.log"));
    let run_log = run_output::open_log(&log, args.force)?;

    // **The FASTA is read to the end before a file is opened.** Two of the run's construction
    // checks compare per-contig checksums — that each sample was aligned to this reference, and
    // that the catalog was built on it — and a reference whose verification has not been joined
    // carries none.
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        args.reference.clone(),
        ReferenceCheck::VerifyAgainstIndex,
    )
    .map_err(|source| CallFromAlignmentsCliError::Reference {
        path: args.reference.clone(),
        source,
    })?;
    let with_checksums =
        match verify {
            Some(handle) => handle.join().map_err(|source| {
                CallFromAlignmentsCliError::ReferenceVerification {
                    path: args.reference.clone(),
                    source,
                }
            })?,
            None => Arc::clone(&info),
        };
    let contigs: ContigList = info.contig_list();
    // **Read before any work**, so a BED that will not read stops the run at once.
    let explain = match args.explain_loci.as_deref() {
        Some([bed, tsv]) => Some((
            ExplainRegions::from_bed_path(bed, &contigs).map_err(|source| {
                CallFromAlignmentsCliError::ExplainRegions {
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

    let ground = ground_request(args);
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    // **Kept because the segmentation takes the original**, and the report names the ground the
    // run undertook to speak for — which is the one thing a VCF cannot say and the whole reason
    // a run reports at all.
    let analysed_for_the_report = analysed.clone();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;
    let satellite_bases_on_the_generic_path = segmentation.satellite_bases_on_the_generic_path();

    let read_groups = build_read_groups(&args.alignments)
        .map_err(|source| CallFromAlignmentsCliError::ReadGroups { source })?;
    say_what_this_run_does(args, &read_groups, &analysed, &log);
    let numbers = calling_run::run_parameters(
        args.parameters.as_deref(),
        args.ploidy,
        &read_groups,
        &with_checksums,
        asked_ploidy,
        segmentation.inputs().repeat_tract_criteria(),
    )?;
    // **The run's ploidy is the parameters', not the flag's.** A supplied file states the
    // ploidy its numbers were fitted at, and the records' `GT` has to be enumerated at the same
    // one the model scored at; `run_parameters` has already refused a flag that disagrees with
    // it, so on either path this is the number the operator meant.
    let ploidy = numbers.parameters.ploidy();
    // **Captured before the parameters move into the caller.** The filter reads one coefficient
    // a sample, in the run's sample order, and the caller consumes the parameters it lives on.
    let inbreeding_by_sample: Vec<InbreedingF> = numbers
        .parameters
        .inbreeding_coefficient_by_sample()
        .to_vec();

    // **The parameters file is assembled here, before a read is decoded, and written after the
    // last one.** Spec §7 makes writing unconditional, and `ParametersFile::of_run` holds its
    // wiring checks in release — that this run's read-group table, its parameters and its
    // inbreeding estimates were all minted from the same inputs. Those are startup questions,
    // and asking them at startup is what stops a panic from discarding a cohort's calling work,
    // which `of_run`'s own note leaves to the driver. **Nothing about the file changes while the
    // run calls**: it records what the run was configured with, not what it found.
    let digest = ReferenceDigest::of(&with_checksums)
        .map_err(|source| calling_run::CallingRunError::ReferenceNotDigested { source })?;
    // **The census the numbers were fitted under, not this run's.** Direct mode never has one of
    // its own (`run_streaming.md` §2) — it reads its evidence from the alignment files, builds no
    // psp and runs no fit — so on the defaults path the file names no terms, which is how a run
    // with no census spells itself; on the supplied path it is what the file recorded
    // ([`TheRunsNumbers::census`]). **What this run counted as a repeat** (`parameters_file.md`
    // §3.9) is taken from the segmentation rather than rebuilt from the flags, so the record cannot
    // say one thing while the catalog was asked another.
    let parameters_file = numbers.parameters_file(
        &read_groups,
        &digest,
        segmentation.inputs().repeat_tract_criteria(),
    );
    // **The numbers go into the log before a read is decoded**, so a run that stops part-way has
    // already recorded what it was scoring with (`doc/devel/ng/spec/command_line.md` §7.3).
    run_log.record_parameters(&parameters_file.to_toml());

    let caller = AlignedFilesVariantCaller::open(
        AlignmentInputs {
            read_groups: &read_groups,
            reference: &reference,
            read_filters: ReadFilterConfig::default(),
            build_index_if_missing: args.build_index_if_missing,
            locus_generator_settings: PileupGeneratorConfig::with_max_reads_per_position(
                args.read_cap(),
            ),
            reference_with_checksums: &with_checksums,
        },
        segmentation,
        numbers.parameters,
        calling_run::calling_loop_config_for_this_run()?,
        candidate_selection,
        merge_parameters,
    )
    .map_err(|source| CallFromAlignmentsCliError::Run { source })?
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
                CallFromAlignmentsCliError::Output {
                    path: args.output.clone(),
                    source,
                }
            })?,
        )
    };

    let mut written = caller
        .call_cohort_handing_each_record_over(
            &SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior),
            &mut |record, window_coverage| sink.accept(record, window_coverage),
        )
        .map_err(|source| CallFromAlignmentsCliError::Run { source })?;

    // **The spill has to stay alive until pass three has read it**, so the sink is taken apart
    // rather than consumed: `finish` would drop the file the next two passes need.
    let parked = match sink {
        CalledRecordSink::StraightToTheVcf(writer) => {
            writer
                .finish()
                .map_err(|source| CallFromAlignmentsCliError::CallsNotWritten {
                    path: args.output.clone(),
                    source: crate::run::paralog_filter::PassOneError::Vcf(source),
                })?;
            None
        }
        CalledRecordSink::ParkedOnTheSpill(mut spilling) => {
            spilling.finish_parking().map_err(|source| {
                CallFromAlignmentsCliError::CallsNotWritten {
                    path: args.output.clone(),
                    source: crate::run::paralog_filter::PassOneError::Spill(source),
                }
            })?;
            Some(spilling)
        }
    };

    // **Passes two and three, through the one function both subcommands call**, so direct mode
    // and psp mode cannot drift into filtering differently (spec §1.1 goal 2).
    // The filter is asked for its verdicts on the records the explained loci wrote.
    let mut paralog_explanations = explain
        .as_ref()
        .map(|_| ParalogExplanations::wanting(&written.calling.explanations));
    let filtered = match &parked {
        Some(spilling) => Some(
            paralog_filter::fit_score_and_write_the_calls_explaining(
                spilling.spill(),
                // **Taken, not cloned.** `ParalogScoringContext::new` consumes the histograms so
                // each sample's bins are freed as its model is fitted; cloning here would defeat
                // that and hold a second copy of every sample's — about 80 kB each. Nothing else
                // reads the field afterwards.
                std::mem::take(&mut written.window_coverage_histograms),
                &inbreeding_by_sample,
                // PANIC-FREE: the sink is the spill only where `asked_for` is `Some`, and this
                // arm is reached only through that sink.
                asked_for.expect("the records are parked only when a filter was asked for"),
                &args.output,
                metadata,
                ploidy,
                paralog_explanations.as_mut(),
            )
            .map_err(|source| CallFromAlignmentsCliError::ParalogFilter { source })?,
        ),
        None => None,
    };

    if let Some(verdicts) = paralog_explanations {
        verdicts.add_to(&mut written.calling.explanations);
    }
    if let Some((_, tsv)) = &explain {
        write_explanations(
            tsv,
            &mut std::mem::take(&mut written.calling.explanations),
            &contigs,
            &sample_names,
        )
        .map_err(|source| CallFromAlignmentsCliError::ExplanationNotWritten {
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
        &RunReport::of(
            &written,
            &contigs,
            &read_groups,
            &parameters_file,
            &analysed_for_the_report,
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
        .map_or(written.calling.records_written, |run| run.did.written);
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
    args: &CallFromAlignmentsArgs,
    read_groups: &ReadGroups,
    analysed: &GenomeRegions,
    log: &Path,
) {
    let samples: BTreeSet<&str> = read_groups
        .iter()
        .map(|(_, group)| &*group.sample)
        .collect();
    let bases: u64 = analysed.iter().map(|region| region.len()).sum();
    let ground = match &args.regions {
        Some(bed) => format!("{} ({})", bed.display(), run_output::size_of_ground(bases)),
        None => format!("the whole genome ({})", run_output::size_of_ground(bases)),
    };
    let threads = rayon::current_num_threads();
    tracing::info!(
        "pvc call-from-alignments: calling {} sample{} over {ground}, {threads} thread{}",
        samples.len(),
        if samples.len() == 1 { "" } else { "s" },
        if threads == 1 { "" } else { "s" },
    );
    tracing::info!("log: {}", log.display());

    tracing::debug!("input: reference {}", args.reference.display());
    let catalog = run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
    tracing::debug!("input: repeat catalogue {}", catalog.display());
    if let Some(bed) = &args.regions {
        tracing::debug!("input: regions {} ({bases} bases)", bed.display());
    }
    for (_, group) in read_groups.iter() {
        tracing::debug!(
            "input: alignments {} — sample {}, read group {}",
            group.file.display(),
            group.sample,
            group.id
        );
    }
    tracing::debug!(
        "output: {}",
        calling_run::where_the_calls_went(&args.output)
    );
}

#[cfg(test)]
mod tests;
