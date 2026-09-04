//! The `call-from-alignments` subcommand: **call a cohort of alignment files and write a VCF**,
//! in one process, with nothing written to disk in between — and, beside that VCF, the model
//! parameters the run scored with.
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
//! Out: **two files**. The VCF, and a `<output stem>.parameters.toml` beside it holding every
//! number the run scored with. The second is not optional and there is no flag to switch it off
//! (`doc/devel/ng/spec/parameters_file.md` §7): a run is then reproducible from its own output,
//! a run that guessed its numbers is auditable in the same form as one that measured them, and
//! somebody who wants to change one number starts from a file rather than from a document. It
//! is also why a run whose `--parameters` file is the one it would write is refused — that is
//! the second command a person types, and it would destroy the edit they just made.
//!
//! **The catalog is not optional and it is not a repeat caller's input.** It is what says where
//! the repeat tracts are, so the run can route them; every locus is called down the SNP/indel
//! path today and a tract is counted as ground this caller has not built yet, which the run's
//! own counts say out loud.
//!
//! # Where this run is short, and it says so rather than being wrong
//!
//! **Repeat tracts are called, and repeat *clusters* are not.** A tract goes down its own path
//! — `select_ssr`, the stutter emission, and the same frequency loop — and its record carries
//! `STR`, `RU`, `PERIOD` and each called allele's `REPCN`
//! (`doc/devel/ng/impl_plan/calling_loop_ssr.md` Milestone C). A repeat cluster with no clean
//! flanks has no generator and no model: its ground is charged to *not built yet* and the
//! summary says how much of the run's ground that was.
//!
//! **The decode is parallel; the rest is one thread (Milestone E1, 2026-09-01).** Decoding
//! reads is 88% of a calling run at 63 samples, and the run sweeps every sample's reader
//! concurrently on the rayon pool while assembling and genotyping stay on one thread — the
//! output is identical at every thread count. Nothing is printed while it works.

use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::ng::calling::parameters_file::{ParametersFile, beside_the_vcf};
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::read_groups::{ReadGroupError, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::ng::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN;
use crate::ng::run::report::BoundsTheRunCalledUnder;
use crate::ng::run::{AlignedFilesVariantCaller, AlignmentInputs, RunError, RunReport};
use crate::ng::types::MAX_MOTIF_LEN;
use crate::ng::vcf::writer::{VcfWriteError, VcfWriter};
use crate::pop_var_caller_exp::calling_run;
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

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
    /// Reference FASTA — the one every alignment file was made against. A `.fai` is built
    /// beside it if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog, which says where the repeat tracts are. Build it first with
    /// `pop_var_caller_exp repeat-catalog --reference <reference>`; it is not optional.
    ///
    /// Defaults to `<reference>.repeats.parquet`, which is where that command writes it. A
    /// catalog built on another reference is refused: its coordinates would put every tract in
    /// the wrong place, genome-wide, with nothing to notice.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One alignment file per sample (BAM or CRAM, indexed). Repeat the flag.
    ///
    /// A sample is named by its files' `@RG SM` tag, and those names, in the order the flags
    /// were given, are the VCF's sample columns. Two files carrying one `SM` are one sample.
    #[arg(long = "alignment", required = true, num_args = 1..)]
    pub alignments: Vec<PathBuf>,

    /// Where to write the VCF. `.vcf.gz` or `.vcf.bgz` writes bgzf; anything else is plain
    /// text.
    ///
    /// The named file appears whole or not at all: the bytes go to `<output>.tmp` and are
    /// renamed into place only once the last record is on disk. A run that stops leaves that
    /// `.tmp` behind.
    #[arg(long)]
    pub output: PathBuf,

    /// BED of the stretch of genome to call over. Without it, every base of every contig.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// The model's numbers, as a parameters file — what a fit over this cohort wrote, or a copy
    /// of one edited by hand. It is refused if it was fitted against another reference or names
    /// read groups this cohort does not have.
    #[arg(long)]
    pub parameters: Option<PathBuf>,

    /// Call with the defaults compiled into this binary instead of a parameters file.
    ///
    /// Not the same claim as a fit. A defaults run assumes no base-quality calibration, no
    /// contamination and no inbreeding; its genotypes are what the reads alone say under those
    /// assumptions.
    #[arg(long)]
    pub defaults: bool,

    /// How many copies of each chromosome every sample carries. Two when nothing says
    /// otherwise.
    ///
    /// A parameters file states its own, and where one is given this flag may only agree with
    /// it: a file records the ploidy its numbers were fitted at, and calling at another one
    /// would score genotypes the fit never saw. A run that types a different number is refused
    /// rather than quietly called at the file's.
    ///
    /// At most 16 copies, which is what the read likelihood scores.
    #[arg(long)]
    pub ploidy: Option<u8>,

    /// Build a `.bai`/`.crai` beside any alignment file that has none.
    #[arg(long, help_heading = "Advanced")]
    pub build_index_if_missing: bool,

    /// The widest a locus may be, in reference bases, before the caller declines to assemble it.
    ///
    /// A deletion joins the positions it covers into one locus, so this is what decides how
    /// long a deletion can be and still be called — and everything else chained into that
    /// locus goes with it. The run report says how many loci were refused and how long they
    /// were; raise this and call again if their lengths cluster just above it. Long reads want
    /// a larger number than short ones.
    #[arg(long, default_value_t = DEFAULT_MAX_COHORT_LOCUS_SPAN, help_heading = "Advanced")]
    pub max_cohort_locus_span: u32,

    /// The most alleles a locus may be called over, the reference counted among them.
    ///
    /// Where more alleles segregate than this, the worst-evidenced are cut — and a sample whose
    /// own reads earned one of the cut sequences is set aside at that locus rather than called
    /// over a set that cannot hold what it carries. The run report says where that left nobody
    /// callable at all.
    #[arg(long, default_value_t = DEFAULT_MAX_CANDIDATE_ALLELES.get(), help_heading = "Advanced")]
    pub max_candidate_alleles: u16,

    /// How much reference one round of locus building covers, in bases. Chosen from the
    /// cohort's size when it is not given.
    ///
    /// The run advances in rounds: it draws every sample's observations over the next stretch
    /// of reference, then builds and calls the loci in it. A wider stretch means fewer rounds
    /// and more of each sample's reading done at once — which is what the threads overlap —
    /// and more observations held at once, which is what it costs in memory. Both scale with
    /// the cohort, so a single number cannot be right at three samples and at three thousand.
    /// Left unset, the run picks a width that holds one round's observations to a fixed
    /// budget: about 8,000 bases at sixty-three samples, 500 at a thousand.
    #[arg(long, help_heading = "Advanced")]
    pub cohort_locus_builder_regions_len: Option<u32>,

    /// How many threads to use. Zero means every core.
    ///
    /// What they parallelise is the reading: each round draws the cohort's samples at once,
    /// and everything after — building the loci, calling them, writing the VCF — stays on one
    /// thread in genome order, so the output does not depend on this number. The other two
    /// subcommands take the same flag; before this one did, a run could only be narrowed
    /// through `RAYON_NUM_THREADS`.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// The fewest motif copies a tract needs before this run treats it as a repeat: six
    /// comma-separated numbers, one per period 1 to 6. Any other count is refused.
    ///
    /// The default is the copy count at which each period starts to stutter, measured over a
    /// tomato archive on 2026-08-10. Below its floor, a tract is ordinary sequence and the
    /// SNP/indel caller handles it.
    #[arg(
        long,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_copies,
        default_value = "8,6,6,6,5,4",
        help_heading = "What counts as a repeat"
    )]
    pub min_copies: MinCopies,

    /// The shortest repeat unit this run treats as a repeat. 1 puts homopolymers on the
    /// repeat path.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub min_period: u8,

    /// The longest repeat unit this run treats as a repeat. Six is the longest the catalog
    /// holds, and the longest a motif can be.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub max_period: u8,

    /// A tract longer than this many bases is a satellite: neither caller speaks for it, and
    /// the run report counts its ground as refused.
    ///
    /// A round number at the read-length limit rather than a measured one — with 150 bp reads
    /// a read spans a tract plus an anchor each side only up to about 90 bp, so past 100 the
    /// repeat path has nothing to offer.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_len: u64,

    /// How much of a tract must match a perfect tiling of its motif, from 0 to 1. Below this
    /// the tract is ordinary sequence.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_purity,
        help_heading = "What counts as a repeat"
    )]
    pub min_purity: f32,
}

/// Everything that can stop a run, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CallFromAlignmentsCliError {
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

    /// The parameters file could not be written beside the output.
    ///
    /// **The VCF is already whole when this can happen**, since the parameters go to disk after
    /// the last record is renamed into place — so a run that reaches this has its calls and not
    /// its provenance, which is what the message has to let an operator work out.
    #[error(
        "the calls are complete and written to {}, but the parameters that produced them could \
         not be saved to {} — keep the VCF and re-run to recover its parameters",
        calls.display(),
        path.display()
    )]
    ParametersNotWritten {
        /// Where they would have gone.
        path: PathBuf,
        /// The VCF that is finished and on disk. **Named because the run still fails**: a
        /// `set -e` pipeline would otherwise throw away a complete, correctly-headed file, and
        /// nothing else in the message says it is there.
        calls: PathBuf,
        /// What the write said.
        #[source]
        source: std::io::Error,
    },

    /// The ground this run would walk could not be worked out — the reference's contigs, the
    /// BED, the catalog, or what this run counts as a repeat.
    ///
    /// **Transparent**, so the sentence a person reads is the one
    /// [`run_ground`](crate::pop_var_caller_exp::run_ground) writes: these refusals belong to
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
            min_purity: args.min_purity,
        },
    }
}

/// Call `--alignment`'s cohort over `--regions` and write `--output`. Prints a summary when it
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

    let asked_ploidy = calling_run::ploidy_asked_for(args.ploidy)?;
    let merge_parameters = calling_run::merge_parameters_for(
        args.max_cohort_locus_span,
        args.cohort_locus_builder_regions_len,
        args.alignments.len(),
    )?;
    let candidate_selection = calling_run::candidate_selection_for(args.max_candidate_alleles)?;
    calling_run::refuse_an_output_that_cannot_be_written(&args.output)?;
    calling_run::refuse_an_output_whose_parameters_file_is_this_run_s_input(
        &args.output,
        args.parameters.as_deref(),
    )?;

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
    let reference = OpenReference::new(info);

    let ground = ground_request(args);
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    // **Kept because the segmentation takes the original**, and the report names the ground the
    // run undertook to speak for — which is the one thing a VCF cannot say and the whole reason
    // a run reports at all.
    let analysed_for_the_report = analysed.clone();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;

    let read_groups = build_read_groups(&args.alignments)
        .map_err(|source| CallFromAlignmentsCliError::ReadGroups { source })?;
    let numbers = calling_run::run_parameters(
        args.parameters.as_deref(),
        args.ploidy,
        &read_groups,
        &with_checksums,
        asked_ploidy,
        &segmentation.inputs().repeat_tract_criteria,
    )?;
    // **The run's ploidy is the parameters', not the flag's.** A supplied file states the
    // ploidy its numbers were fitted at, and the records' `GT` has to be enumerated at the same
    // one the model scored at; `run_parameters` has already refused a flag that disagrees with
    // it, so on either path this is the number the operator meant.
    let ploidy = numbers.parameters.ploidy();

    // **The parameters file is assembled here, before a read is decoded, and written after the
    // last one.** Spec §7 makes writing unconditional, and `ParametersFile::of_run` holds its
    // wiring checks in release — that this run's read-group table, its parameters and its
    // inbreeding estimates were all minted from the same inputs. Those are startup questions,
    // and asking them at startup is what stops a panic from discarding a cohort's calling work,
    // which `of_run`'s own note leaves to the driver. **Nothing about the file changes while the
    // run calls**: it records what the run was configured with, not what it found.
    let digest = ReferenceDigest::of(&with_checksums)
        .map_err(|source| calling_run::CallingRunError::ReferenceNotDigested { source })?;
    let parameters_file = ParametersFile::of_run(
        &numbers.parameters,
        &read_groups,
        &numbers.reads_behind_each_calibration,
        &numbers.inbreeding_by_sample,
        &digest,
        // **The census the numbers were fitted under, not this run's.** Direct mode never has
        // one of its own (`run_streaming.md` §2) — it reads its evidence from the alignment
        // files, builds no psp and runs no fit — so on the defaults path this names no terms,
        // which is how a run with no census spells itself. On the supplied path it is what the
        // file recorded, because a run writing its parameters out again writes back the terms it
        // read; see [`TheRunsNumbers::census`].
        numbers.census.clone(),
        // **What this run counted as a repeat** (`parameters_file.md` §3.9) — taken from the
        // segmentation rather than rebuilt from the flags, so the record cannot say one thing
        // while the catalog was asked another.
        &segmentation.inputs().repeat_tract_criteria,
    );

    let caller = AlignedFilesVariantCaller::open(
        AlignmentInputs {
            read_groups: &read_groups,
            reference: &reference,
            read_filters: ReadFilterConfig::default(),
            build_index_if_missing: args.build_index_if_missing,
            locus_generator_settings: PileupGeneratorConfig::default(),
            reference_with_checksums: &with_checksums,
        },
        segmentation,
        numbers.parameters,
        calling_run::calling_loop_config_for_this_run()?,
        candidate_selection,
        merge_parameters,
    )
    .map_err(|source| CallFromAlignmentsCliError::Run { source })?;

    let metadata = calling_run::header_for(
        &args.output,
        &args.reference,
        &contigs,
        &with_checksums,
        caller.sample_names().map(str::to_owned).collect(),
    )?;
    let mut writer = VcfWriter::create(&args.output, metadata, ploidy).map_err(|source| {
        CallFromAlignmentsCliError::Output {
            path: args.output.clone(),
            source,
        }
    })?;

    let written = caller
        .call_cohort_handing_each_record_over(
            &SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior),
            &mut |record| writer.write_record(record),
        )
        .map_err(|source| CallFromAlignmentsCliError::Run { source })?;

    writer
        .finish()
        .map_err(|source| CallFromAlignmentsCliError::Output {
            path: args.output.clone(),
            source,
        })?;

    // **After the VCF is on disk, not before.** Spec §7's three purposes — a run reproducible
    // from its own output, a defaults run auditable, an edit that starts from something — are
    // all about a run that finished, and a parameters file standing beside a VCF that does not
    // exist would answer none of them. The file goes to a temporary in the same directory and is
    // renamed over its destination, so what lands is whole or is not there at all.
    // **A file already there is replaced, and the run says so.** Spec §7 makes writing
    // unconditional, so this is not a refusal — but a person who edited the parameters in place
    // and then re-ran to compare has just lost the edit, and a run that did that in silence was
    // measured doing exactly that. The `--parameters` route is refused outright
    // (`refuse_an_output_whose_parameters_file_is_this_run_s_input`); this is the other route to
    // the same loss, where no flag names the file.
    let would_write = beside_the_vcf(&args.output);
    if would_write.exists() {
        eprintln!(
            "note: replacing the parameters file already at {}",
            would_write.display()
        );
    }
    let parameters_at = parameters_file
        .write_beside_the_vcf(&args.output)
        .map_err(|source| CallFromAlignmentsCliError::ParametersNotWritten {
            path: would_write,
            calls: args.output.clone(),
            source,
        })?;

    calling_run::print_report(
        &args.output,
        &parameters_at,
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
        ),
    );
    Ok(())
}

#[cfg(test)]
mod tests;
