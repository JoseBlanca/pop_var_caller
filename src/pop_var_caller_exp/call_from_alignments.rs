//! The `call-from-alignments` subcommand: **call a cohort of alignment files and write a VCF**,
//! in one process, with nothing written to disk in between.
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
//! file or `--defaults`. Out: one VCF.
//!
//! **The catalog is not optional and it is not a repeat caller's input.** It is what says where
//! the repeat tracts are, so the run can route them; every locus is called down the SNP/indel
//! path today and a tract is counted as ground this caller has not built yet, which the run's
//! own counts say out loud.
//!
//! # Where this run is short, and it says so rather than being wrong
//!
//! **Repeat tracts are analysed and not called.** Candidate selection at a tract is specified
//! and unbuilt (`doc/devel/ng/impl_plan/candidate_alleles_ssr.md`), so both tract generator
//! slots are refused as unbuilt and their ground is charged to *not built yet*. A run over
//! tract-rich ground is therefore short, not wrong, and the summary says by how much.
//!
//! **It is single-threaded.** Decoding reads is 94–97% of a calling run and it runs on one
//! thread; the parallel form exists and is not reached from here (the run driver's plan,
//! Milestone E, deferred 2026-09-01).

use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::calling::allele_candidates::CandidateSelectionConfig;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::ng::calling::inference::CallingLoopConfig;
use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::ng::calling::likelihood::MAX_PLOIDY_COPIES;
use crate::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::ng::calling::parameters_file::{
    DeclaredInbreeding, ParametersFile, ParametersFileError,
};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::parameter_estimation::joint::loci::{ReferenceDigest, SelectionError};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::read_groups::{ReadGroupError, ReadGroups, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfo, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::GenomeRegions;
use crate::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria, sibling_catalog_path,
};
use crate::ng::run::{
    AlignedFilesVariantCaller, AlignmentInputs, MergeParameters, RunError, Segmentation,
    WrittenCohort,
};
use crate::ng::types::{DomainError, Ploidy};
use crate::ng::vcf::header::{HeaderContig, HeaderMetadataError, VcfHeaderMetadata};
use crate::ng::vcf::writer::{VcfWriteError, VcfWriter};
use crate::pop_var_caller::common::current_command_line;
use crate::regions::{BedError, ContigBounds};

/// The ploidy a run assumes when it is not told.
///
/// **A property of the run and not of any fit** (`doc/devel/ng/spec/parameters_file.md` §3.2),
/// which is why it is a flag here rather than a number read out of a parameters file.
const DEFAULT_PLOIDY: u8 = 2;

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

    /// There is no repeat catalog to route this run's ground with.
    #[error(
        "no repeat catalog at {}; build one with `pop_var_caller_exp repeat-catalog --reference {}`",
        path.display(),
        reference.display()
    )]
    MissingCatalog {
        /// Where one was looked for.
        path: PathBuf,
        /// The reference it would be built from.
        reference: PathBuf,
    },

    /// The repeat catalog could not be used.
    #[error("the repeat catalog {} could not be used", path.display())]
    Catalog {
        /// The catalog.
        path: PathBuf,
        /// What the catalog said.
        #[source]
        source: RepeatCatalogError,
    },

    /// The BED naming the ground to call could not be read.
    #[error("the regions file {} could not be read", path.display())]
    Bed {
        /// The BED.
        path: PathBuf,
        /// What the parser said.
        #[source]
        source: BedError,
    },

    /// A ploidy of zero, or one no genotype table can enumerate.
    #[error("--ploidy {asked}")]
    Ploidy {
        /// What was asked for.
        asked: u8,
        /// Why it is not a ploidy.
        #[source]
        source: DomainError,
    },

    /// A ploidy past what the read likelihood scores.
    ///
    /// **This was a panic until 2026-09-01**, and reached after the whole cohort had been
    /// opened: `Ploidy::try_new` turns down zero and nothing else, and the copy-share table the
    /// read likelihood builds asserts at seventeen. A polyploid crop is an ordinary thing to
    /// call — sugarcane runs to about twelve copies — so a person who types twenty is told the
    /// ceiling rather than handed a backtrace.
    #[error(
        "--ploidy {asked}: this caller scores genotypes up to {limit} copies of the genome, and \
         a run at {asked} would have nothing to score them against"
    )]
    PloidyPastWhatIsScored {
        /// What was asked for.
        asked: u8,
        /// The most copies the read likelihood builds a table for.
        limit: usize,
    },

    /// `--output` names a directory rather than a file.
    #[error(
        "--output {} is a directory; give the path of the VCF file to write, such as {}/calls.vcf.gz",
        path.display(),
        path.display()
    )]
    OutputIsADirectory {
        /// What was typed.
        path: PathBuf,
    },

    /// The directory `--output` would be written into does not exist.
    ///
    /// **Refused before the reference is read**, because it is a typing mistake that would
    /// otherwise be discovered after every alignment file had been opened.
    #[error(
        "--output {} cannot be written: there is no directory {}",
        path.display(),
        directory.display()
    )]
    OutputDirectoryIsMissing {
        /// What was typed.
        path: PathBuf,
        /// The directory it would have gone into.
        directory: PathBuf,
    },

    /// `--ploidy` was typed and disagrees with the parameters file's own.
    ///
    /// **A file's numbers were fitted at one ploidy and mean nothing at another**
    /// (`doc/devel/ng/spec/parameters_file.md` §3.2), so a run that wants the other number
    /// wants a different fit. Refused rather than called at the file's, which would be a run
    /// silently ignoring what its operator typed.
    #[error(
        "--ploidy {asked} disagrees with {}, whose numbers were fitted at ploidy {in_the_file}; \
         drop the flag to call at the file's ploidy, or fit the cohort at {asked}",
        path.display()
    )]
    PloidyIsNotTheParametersFiles {
        /// The parameters file.
        path: PathBuf,
        /// What was typed.
        asked: u8,
        /// What the file states.
        in_the_file: u8,
    },

    /// A contig longer than a `u32` can name.
    ///
    /// **The analysed ground is resolved against contig lengths a `u32` carries**
    /// (`ContigBounds`), so a reference whose contig is longer than that cannot be resolved
    /// against — and narrowing it silently would compute the run's ground against a wrong
    /// length. No assembly in use has one; the refusal is what makes that a fact rather than
    /// an assumption.
    #[error(
        "contig {name} of {} is {length} bases, longer than the {limit} a run can address",
        reference.display()
    )]
    ContigTooLong {
        /// The reference.
        reference: PathBuf,
        /// The contig that is too long.
        name: String,
        /// Its length.
        length: u64,
        /// The longest a run can address.
        limit: u64,
    },

    /// The alignment files' read groups could not be read.
    #[error("the read groups of this run's alignment files could not be read")]
    ReadGroups {
        /// What the reader said.
        #[source]
        source: ReadGroupError,
    },

    /// The parameters file could not be read, or is not this run's.
    #[error("the parameters file {} could not be used", path.display())]
    Parameters {
        /// The file.
        path: PathBuf,
        /// What reading or binding it said.
        #[source]
        source: ParametersFileError,
    },

    /// The parameters file could not be opened at all.
    #[error("the parameters file {} could not be opened", path.display())]
    ParametersUnreadable {
        /// The file.
        path: PathBuf,
        /// What opening it said.
        #[source]
        source: std::io::Error,
    },

    /// A parameters file cannot be bound to a reference the run cannot digest.
    ///
    /// **A `.fai` describes a genome's geometry and holds no bases**, so a run driven from one
    /// has nothing to compare a file's reference binding against. This run reads the FASTA, so
    /// it is unreachable through the command line and is carried rather than unwrapped.
    #[error("this run's reference carries no digest, so a parameters file cannot be bound to it")]
    ReferenceNotDigested {
        /// What the digest said.
        #[source]
        source: SelectionError,
    },

    /// The run's shipped calling-loop settings would not validate.
    ///
    /// **Not reachable from any flag** — the settings are compiled in — so this is a defect in
    /// this binary rather than in what was typed. It is an error and not a panic because a
    /// message naming the setting is worth more than a backtrace.
    #[error("this binary's calling-loop settings are not runnable: {0}")]
    CallingLoopSettings(String),

    /// The header could not honestly state what this run is.
    #[error("the output's header could not be written")]
    Header {
        /// What the header refused.
        #[source]
        source: HeaderMetadataError,
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
    let asked_ploidy = ploidy_asked_for(args)?;
    refuse_an_output_that_cannot_be_written(args)?;

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

    let analysed = analysed_regions(args, &contigs)?;
    let segmentation = segments_over(args, &analysed, &with_checksums)?;

    let read_groups = build_read_groups(&args.alignments)
        .map_err(|source| CallFromAlignmentsCliError::ReadGroups { source })?;
    let parameters = run_parameters(args, &read_groups, &with_checksums, asked_ploidy)?;
    // **The run's ploidy is the parameters', not the flag's.** A supplied file states the
    // ploidy its numbers were fitted at, and the records' `GT` has to be enumerated at the same
    // one the model scored at; `run_parameters` has already refused a flag that disagrees with
    // it, so on either path this is the number the operator meant.
    let ploidy = parameters.ploidy();

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
        parameters,
        CallingLoopConfig::DEFAULT.validate().map_err(|source| {
            CallFromAlignmentsCliError::CallingLoopSettings(source.to_string())
        })?,
        CandidateSelectionConfig::DEFAULT,
        MergeParameters::DEFAULT,
    )
    .map_err(|source| CallFromAlignmentsCliError::Run { source })?;

    let metadata = header_for(args, &contigs, &with_checksums, &caller)?;
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

    report(args, &written);
    Ok(())
}

/// **The ploidy this run was asked for, judged against what the caller can score.**
///
/// Two refusals, and the second is the one that used to be a panic. `Ploidy::try_new` turns
/// down zero and nothing else, so seventeen copies is a value the type admits and the read
/// likelihood's copy-share table does not — it asserts, and a person who typed `--ploidy 20`
/// got a Rust backtrace naming a source file after the whole cohort had been opened. A
/// polyploid crop is an ordinary thing to call, so the number is judged here, before anything
/// is read, and the ceiling is named.
fn ploidy_asked_for(args: &CallFromAlignmentsArgs) -> Result<Ploidy, CallFromAlignmentsCliError> {
    let asked = args.ploidy.unwrap_or(DEFAULT_PLOIDY);
    let ploidy = Ploidy::try_new(asked)
        .map_err(|source| CallFromAlignmentsCliError::Ploidy { asked, source })?;
    if usize::from(asked) > MAX_PLOIDY_COPIES {
        return Err(CallFromAlignmentsCliError::PloidyPastWhatIsScored {
            asked,
            limit: MAX_PLOIDY_COPIES,
        });
    }
    Ok(ploidy)
}

/// **Somewhere to write must exist before a cohort is opened.**
///
/// The writer's own refusals are honest and arrive far too late: a missing directory is
/// discovered after the reference is read and every alignment file is opened, and an `--output`
/// naming a *directory* is discovered only after the last locus has been called — leaving the
/// in-flight `<output>.tmp` beside it. Both are typing mistakes, and both are visible from the
/// path alone.
///
/// **What it cannot check is permission**, which on every platform is answerable only by
/// writing. That failure still comes from the writer, and it still names the path.
fn refuse_an_output_that_cannot_be_written(
    args: &CallFromAlignmentsArgs,
) -> Result<(), CallFromAlignmentsCliError> {
    if args.output.is_dir() {
        return Err(CallFromAlignmentsCliError::OutputIsADirectory {
            path: args.output.clone(),
        });
    }
    let directory = args.output.parent().unwrap_or_else(|| Path::new("."));
    // An empty parent is what a bare file name gives, and it means the working directory.
    if !directory.as_os_str().is_empty() && !directory.is_dir() {
        return Err(CallFromAlignmentsCliError::OutputDirectoryIsMissing {
            path: args.output.clone(),
            directory: directory.to_path_buf(),
        });
    }
    Ok(())
}

/// The ground this run calls over: the BED it was given, or every base of every contig.
fn analysed_regions(
    args: &CallFromAlignmentsArgs,
    contigs: &ContigList,
) -> Result<GenomeRegions, CallFromAlignmentsCliError> {
    // **The narrowing is refused, not taken.** `ContigBounds` carries a `u32`, and a contig
    // past that would have its ground resolved against a wrong length with nothing to notice —
    // the rule `typed_regions.rs`'s `ContigTooLong` records, and the reason the one precedent
    // in this tree that casts is a test.
    let mut bounds: Vec<ContigBounds<'_>> = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        let length =
            u32::try_from(entry.length).map_err(|_| CallFromAlignmentsCliError::ContigTooLong {
                reference: args.reference.clone(),
                name: entry.name.clone(),
                length: entry.length,
                limit: u64::from(u32::MAX),
            })?;
        bounds.push(ContigBounds {
            name: &entry.name,
            length,
        });
    }
    match &args.regions {
        Some(bed) => GenomeRegions::from_bed_path(bed, &bounds).map_err(|source| {
            CallFromAlignmentsCliError::Bed {
                path: bed.clone(),
                source,
            }
        }),
        None => Ok(GenomeRegions::whole_contigs(&bounds)),
    }
}

/// The run's segments: the analysed ground cut into the stretches each generator owns, drawn
/// from the catalog.
fn segments_over(
    args: &CallFromAlignmentsArgs,
    analysed: &GenomeRegions,
    with_checksums: &ReferenceInfo,
) -> Result<Segmentation, CallFromAlignmentsCliError> {
    let path = args
        .catalog
        .clone()
        .unwrap_or_else(|| sibling_catalog_path(&args.reference));
    if !path.exists() {
        return Err(CallFromAlignmentsCliError::MissingCatalog {
            path,
            reference: args.reference.clone(),
        });
    }
    let criteria = StrRepeatCriteria::default();
    let catalog = RepeatCatalog::open_checking_against_reference(&path, with_checksums).map_err(
        |source| CallFromAlignmentsCliError::Catalog {
            path: path.clone(),
            source,
        },
    )?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))
        .map_err(|source| CallFromAlignmentsCliError::Catalog {
            path: path.clone(),
            source,
        })?;
    Segmentation::build(
        segments,
        analysed.clone(),
        catalog.header().clone(),
        criteria,
        path,
    )
    .map_err(|source| CallFromAlignmentsCliError::Run { source })
}

/// The numbers this run scores with — a supplied file, or the defaults compiled in.
///
/// **A supplied file is bound to this run at its own door**: it is refused if it was fitted
/// against another reference or names read groups this cohort does not have, naming the
/// position where the two lists diverge (`doc/devel/ng/spec/parameters_file.md` §6). No census
/// is compared against, because direct mode has none — §2.1 settles that this keeps the file's
/// warrants rather than demoting them.
fn run_parameters(
    args: &CallFromAlignmentsArgs,
    read_groups: &ReadGroups,
    with_checksums: &ReferenceInfo,
    ploidy: Ploidy,
) -> Result<RunParameters, CallFromAlignmentsCliError> {
    let Some(path) = &args.parameters else {
        return Ok(RunParameters::of_defaults(
            read_groups,
            ploidy,
            &DeclaredInbreeding::nothing_said(),
        ));
    };
    let text = std::fs::read_to_string(path).map_err(|source| {
        CallFromAlignmentsCliError::ParametersUnreadable {
            path: path.clone(),
            source,
        }
    })?;
    let file = ParametersFile::from_toml(&text).map_err(|source| {
        CallFromAlignmentsCliError::Parameters {
            path: path.clone(),
            source,
        }
    })?;
    let digest = ReferenceDigest::of(with_checksums)
        .map_err(|source| CallFromAlignmentsCliError::ReferenceNotDigested { source })?;
    let bound = file
        .to_run_parameters_for(&digest, read_groups, None)
        .map_err(|source| CallFromAlignmentsCliError::Parameters {
            path: path.clone(),
            source,
        })?;
    let parameters = bound.from_file.parameters;
    // **A flag that was typed may only agree with the file.** Spec §3.2 puts the ploidy in the
    // file so that a supplied one "cannot be paired with a run at a different ploidy without
    // saying so", and calling at the file's number while an operator typed another is exactly
    // that. Only a ploidy that was *typed* is compared: `--ploidy` is an `Option` for this
    // reason, so a tetraploid file is not refused for the flag's default being two.
    if let Some(asked) = args.ploidy
        && asked != parameters.ploidy().get()
    {
        return Err(CallFromAlignmentsCliError::PloidyIsNotTheParametersFiles {
            path: path.clone(),
            asked,
            in_the_file: parameters.ploidy().get(),
        });
    }
    Ok(parameters)
}

/// What the file's header states about the run that wrote it.
///
/// **The `##parametersFile` line is left off**, because this step writes no parameters file
/// beside its output and a header naming one that is not there would send a reader looking for
/// it. The run driver's plan step F2 is what writes the file and fills the line in.
fn header_for(
    args: &CallFromAlignmentsArgs,
    contigs: &ContigList,
    with_checksums: &ReferenceInfo,
    caller: &AlignedFilesVariantCaller,
) -> Result<VcfHeaderMetadata, CallFromAlignmentsCliError> {
    let header_contigs = contigs
        .entries
        .iter()
        .enumerate()
        .map(|(index, entry)| HeaderContig {
            name: entry.name.clone(),
            length: entry.length,
            md5: with_checksums
                .contigs
                .get(index)
                .and_then(|contig| contig.md5),
        })
        .collect();
    VcfHeaderMetadata::try_new(
        header_contigs,
        caller.sample_names().map(str::to_owned).collect(),
        current_command_line(),
        args.reference.display().to_string(),
        String::new(),
    )
    .map_err(|source| CallFromAlignmentsCliError::Header { source })
}

/// **What the run wrote, and what it could not speak for** — printed because the two questions a
/// person has beside a new VCF are how many records are in it and how much of the ground it
/// covers.
///
/// **The last two rows are what stops an empty file being read as an empty genome.** Every locus
/// goes down the SNP/indel path today and a repeat tract is charged to *not built yet*, so a run
/// over tract-rich ground is short rather than wrong — and a summary of six zeros with no reason
/// beside them cannot be told from a run that looked everywhere and found nothing. **Measured on
/// a 60-base `AT` tract at 24 reads a sample: every count was zero and the exit status was
/// success**, which is what this row now answers.
///
/// The run report proper — every refusal, every parameter that was defaulted rather than fitted,
/// the per-read-group filter tallies — is the run driver's plan step F3; this is the minimum a
/// command that writes a file has to say about it.
fn report(args: &CallFromAlignmentsArgs, written: &WrittenCohort) {
    println!("output\t{}", args.output.display());
    println!("samples\t{}", written.walk.per_sample.len());
    println!("records_written\t{}", written.records_written);
    // **Indented, because these two are the parts of the total above them and the four rows
    // below are not.** A flat list of counts reads as a partition, and only this pair is one.
    println!("  loci_called\t{}", written.loci_called());
    println!(
        "  loci_called_establishing_no_variant\t{}",
        written.loci_called_but_not_written
    );
    println!(
        "loci_too_wide_to_assemble\t{}",
        written.loci_too_wide_to_assemble.len()
    );
    println!(
        "loci_with_nobody_to_call\t{}",
        written.loci_with_nobody_to_call.len()
    );

    // **Every sample walks the same ground, so one sample's region tally is the run's.** The
    // loci each sample's walk emitted differ and are not this line's business.
    let Some(ground) = written.walk.per_sample.first().map(|walk| &walk.regions) else {
        return;
    };
    println!(
        "regions_walked\t{} of {}",
        ground.regions_handled, ground.regions_in
    );
    println!(
        "bases_not_called_repeat_tracts_this_caller_has_not_built\t{}",
        ground.unhandled_not_implemented_bp
    );
    println!(
        "bases_never_called_satellite\t{}",
        ground.unhandled_out_of_scope_bp
    );
}

#[cfg(test)]
mod tests;
