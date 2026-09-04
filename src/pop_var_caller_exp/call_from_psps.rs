//! `call-from-psps` — psp mode's calling stage at the command line.
//!
//! **The same cohort, the same parameters, the same VCF** (spec `run_streaming.md` §12.3). What
//! differs from `call-from-alignments` is where the evidence comes from: this command opens one
//! stored file per sample instead of one alignment file per sample, and every stage above the
//! observations — the merge, the candidate selection, the genotyping, the VCF — is the body
//! both modes drive. The numbers each run scores with come from one place too
//! ([`calling_run`](super::calling_run)), because two copies of that decision would be two
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

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::ng::calling::parameters_file::{ParametersFile, beside_the_vcf};
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
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
use crate::ng::run::{OpenPspCohort, PspVariantCaller, RunError, RunReport, StoredCohortInputs};
use crate::ng::types::MAX_MOTIF_LEN;
use crate::ng::vcf::writer::{VcfWriteError, VcfWriter};
use crate::pop_var_caller_exp::calling_run::{self, CallingRunError};
use crate::pop_var_caller_exp::generate_psps::PSP_FILE_EXTENSION;
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line.
pub const SUBCOMMAND: &str = "call-from-psps";

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
    /// Reference FASTA — the one every psp's samples were aligned to. A `.fai` is built beside
    /// it if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog the psps were walked against. Defaults to
    /// `<reference>.repeats.parquet`.
    ///
    /// **A cohort walked under another catalog is refused**, naming the field that differs:
    /// every psp records the catalog its segmentation came from, and calling over a different
    /// one would score stored loci against a routing that never produced them.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One psp per sample, or a directory holding them. Repeat the flag.
    ///
    /// A directory contributes every `.psp` file directly inside it, in name order. **The order
    /// makes no difference to the calls**: samples are matched to their parameters by name, and
    /// the VCF's sample columns follow the order given here.
    #[arg(long = "psp", required = true, num_args = 1..)]
    pub psps: Vec<PathBuf>,

    /// Where to write the VCF. `.vcf.gz` or `.vcf.bgz` writes bgzf; anything else is plain
    /// text.
    ///
    /// The named file appears whole or not at all: the bytes go to `<output>.tmp` and are
    /// renamed into place only once the last record is on disk. A run that stops leaves that
    /// `.tmp` behind.
    #[arg(long)]
    pub output: PathBuf,

    /// The model's numbers, as a parameters file. It is refused if it was fitted against
    /// another reference or does not name this cohort's samples and read groups.
    #[arg(long)]
    pub parameters: Option<PathBuf>,

    /// Call with the defaults compiled into this binary instead of a parameters file.
    ///
    /// Not the same claim as a fit. A defaults run assumes no base-quality calibration, no
    /// contamination and no inbreeding; its genotypes are what the stored reads alone say under
    /// those assumptions.
    #[arg(long)]
    pub defaults: bool,

    /// How many copies of each chromosome every sample carries. Two when nothing says
    /// otherwise; where a parameters file is given this flag may only agree with it.
    #[arg(long)]
    pub ploidy: Option<u8>,

    /// The widest a locus may be, in reference bases, before the caller declines to assemble it.
    #[arg(long, default_value_t = DEFAULT_MAX_COHORT_LOCUS_SPAN, help_heading = "Advanced")]
    pub max_cohort_locus_span: u32,

    /// The most alleles a locus may be called over, the reference counted among them.
    #[arg(long, default_value_t = DEFAULT_MAX_CANDIDATE_ALLELES.get(), help_heading = "Advanced")]
    pub max_candidate_alleles: u16,

    /// How much reference one round of locus building covers, in bases. Chosen from the
    /// cohort's size when it is not given.
    #[arg(long, help_heading = "Advanced")]
    pub cohort_locus_builder_regions_len: Option<u32>,

    /// How many threads to use. Zero means every core.
    ///
    /// The output does not depend on this number: what the threads parallelise is the reading,
    /// and everything after — building the loci, calling them, writing the VCF — stays on one
    /// thread in genome order.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// The fewest motif copies a tract needs before this run treats it as a repeat: six
    /// comma-separated numbers, one per period 1 to 6.
    ///
    /// **These must be what the psps were walked with**, or the cohort is refused: a run's
    /// segmentation has to be the one that produced the stored loci.
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

    /// The longest repeat unit this run treats as a repeat.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub max_period: u8,

    /// A tract longer than this many bases is a satellite: neither caller speaks for it.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_len: u64,

    /// How much of a tract must match a perfect tiling of its motif, from 0 to 1.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_purity,
        help_heading = "What counts as a repeat"
    )]
    pub min_purity: f32,
}

/// Everything that can stop a run over stored files, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CallFromPspsCliError {
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

    /// The parameters file could not be written beside the output.
    ///
    /// **The VCF is already whole when this can happen**, since the parameters go to disk after
    /// the last record is renamed into place — so a run that reaches this has its calls and not
    /// its provenance, which is what the message has to let an operator work out. **The same
    /// sentence direct mode writes**, because the loss and the repair are the same: a `set -e`
    /// pipeline that read only *could not be written* would throw away a complete,
    /// correctly-headed file.
    #[error(
        "the calls are complete and written to {}, but the parameters that produced them could \
         not be saved to {} — keep the VCF and re-run to recover its parameters",
        calls.display(),
        path.display()
    )]
    ParametersNotWritten {
        /// Where it would have gone.
        path: PathBuf,
        /// The VCF that is on disk.
        calls: PathBuf,
        /// What the write said.
        #[source]
        source: std::io::Error,
    },

    /// The run itself was refused, or stopped.
    #[error("the run stopped")]
    Run {
        /// What stopped it.
        #[source]
        source: RunError,
    },
}

/// This command's flags, as the shared ground assembly asks for them.
///
/// **`regions` is `None` and that is not an omission**: psp mode takes its analysed ground from
/// the files' own headers, so the only thing this request is used for is the catalog and the
/// routing criteria.
fn ground_request(args: &CallFromPspsArgs) -> run_ground::GroundRequest<'_> {
    run_ground::GroundRequest {
        reference: &args.reference,
        catalog: args.catalog.as_deref(),
        regions: None,
        routing: run_ground::RepeatRouting {
            min_copies: args.min_copies,
            min_period: args.min_period,
            max_period: args.max_period,
            max_str_len: args.max_str_len,
            min_purity: args.min_purity,
        },
    }
}

/// **The psps this run calls, with every directory expanded** — one entry a sample, in the
/// order they were given.
///
/// A directory contributes the `.psp` files directly inside it, **sorted by name**, so that two
/// runs naming one directory open the same cohort in the same order however the filesystem
/// answers. Sub-directories are not descended into: a psp's home is a `generate-psps`
/// `--output-dir`, which is flat.
///
/// # Errors
///
/// [`CallFromPspsCliError::PspDirectory`] for a directory that will not list, and
/// [`CallFromPspsCliError::NoPspsInDirectory`] for one holding no psp.
fn psps_named_by(args: &CallFromPspsArgs) -> Result<Vec<PathBuf>, CallFromPspsCliError> {
    let mut paths = Vec::with_capacity(args.psps.len());
    for named in &args.psps {
        if !named.is_dir() {
            paths.push(named.clone());
            continue;
        }
        let mut inside = Vec::new();
        for entry in
            std::fs::read_dir(named).map_err(|source| CallFromPspsCliError::PspDirectory {
                path: named.clone(),
                source,
            })?
        {
            let entry = entry.map_err(|source| CallFromPspsCliError::PspDirectory {
                path: named.clone(),
                source,
            })?;
            let path = entry.path();
            if path.extension().is_some_and(|it| it == PSP_FILE_EXTENSION) && path.is_file() {
                inside.push(path);
            }
        }
        if inside.is_empty() {
            return Err(CallFromPspsCliError::NoPspsInDirectory {
                path: named.clone(),
            });
        }
        inside.sort();
        paths.extend(inside);
    }
    Ok(paths)
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

    let asked_ploidy = calling_run::ploidy_asked_for(args.ploidy)?;
    let candidate_selection = calling_run::candidate_selection_for(args.max_candidate_alleles)?;
    calling_run::refuse_an_output_that_cannot_be_written(&args.output)?;
    calling_run::refuse_an_output_whose_parameters_file_is_this_run_s_input(
        &args.output,
        args.parameters.as_deref(),
    )?;
    let paths = psps_named_by(args)?;
    // **The round width comes from the cohort's size**, and the cohort is the psps — so this
    // waits for the directories to be expanded, where direct mode can count its `--alignment`
    // flags straight away.
    let merge_parameters = calling_run::merge_parameters_for(
        args.max_cohort_locus_span,
        args.cohort_locus_builder_regions_len,
        paths.len(),
    )?;

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
    let reference = OpenReference::new(info);

    // **The cohort is opened before the segmentation is built, because the cohort is what says
    // what ground there is** (spec §5.3): the analysed regions come from the headers, and the
    // files have already been forced to agree about them.
    let cohort =
        OpenPspCohort::open(&paths).map_err(|source| CallFromPspsCliError::Run { source })?;
    let ground = ground_request(args);
    let analysed = cohort.analysed_regions().clone();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;

    let numbers = calling_run::run_parameters(
        args.parameters.as_deref(),
        args.ploidy,
        cohort.read_groups(),
        &with_checksums,
        asked_ploidy,
        &segmentation.inputs().repeat_tract_criteria,
    )?;
    let ploidy = numbers.parameters.ploidy();

    let digest = ReferenceDigest::of(&with_checksums)
        .map_err(|source| CallingRunError::ReferenceNotDigested { source })?;
    let parameters_file = ParametersFile::of_run(
        &numbers.parameters,
        cohort.read_groups(),
        &numbers.reads_behind_each_calibration,
        &numbers.inbreeding_by_sample,
        &digest,
        numbers.census.clone(),
        &segmentation.inputs().repeat_tract_criteria,
    );

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
    .map_err(|source| CallFromPspsCliError::Run { source })?;

    let read_groups = caller.read_groups().clone();
    let metadata = calling_run::header_for(
        &args.output,
        &args.reference,
        &contigs,
        &with_checksums,
        caller.sample_names().map(str::to_owned).collect(),
    )?;
    let mut writer = VcfWriter::create(&args.output, metadata, ploidy).map_err(|source| {
        CallFromPspsCliError::Output {
            path: args.output.clone(),
            source,
        }
    })?;

    let (calling, stored) = caller
        .call_cohort_handing_each_record_over(
            &SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior),
            &mut |record| writer.write_record(record),
        )
        .map_err(|source| CallFromPspsCliError::Run { source })?;

    writer
        .finish()
        .map_err(|source| CallFromPspsCliError::Output {
            path: args.output.clone(),
            source,
        })?;

    // **After the VCF is on disk, not before** — direct mode's rule and its reasons: a
    // parameters file standing beside a VCF that does not exist answers none of spec §7's three
    // purposes.
    let would_write = beside_the_vcf(&args.output);
    if would_write.exists() {
        eprintln!(
            "note: replacing the parameters file already at {}",
            would_write.display()
        );
    }
    let parameters_at = parameters_file
        .write_beside_the_vcf(&args.output)
        .map_err(|source| CallFromPspsCliError::ParametersNotWritten {
            path: would_write,
            calls: args.output.clone(),
            source,
        })?;

    calling_run::print_report(
        &args.output,
        &parameters_at,
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
        ),
    );
    Ok(())
}
