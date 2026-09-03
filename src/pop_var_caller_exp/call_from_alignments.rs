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

use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::calling::allele_candidates::{
    CandidateSelectionConfig, DEFAULT_MAX_CANDIDATE_ALLELES, MaxCandidateAlleles,
};
use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::ng::calling::inference::CallingLoopConfig;
use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::ng::calling::likelihood::MAX_PLOIDY_COPIES;
use crate::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::ng::calling::parameters_file::{
    CensusIdentity, DeclaredInbreeding, ParametersFile, ParametersFileError,
    ReadsBehindEachCalibration, beside_the_vcf,
};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::parameter_estimation::Estimate;
use crate::ng::parameter_estimation::joint::loci::{ReferenceDigest, SelectionError};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::read_groups::{ReadGroupError, ReadGroups, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfo, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies, SsrSegmentCriteria,
};
use crate::ng::region_typing::{DEFAULT_MAX_STR_LEN, GenomeRegions, TypedRegionConfig};
use crate::ng::repeat_catalog::{
    CriteriaRefusal, ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
    sibling_catalog_path,
};
use crate::ng::run::cohort_merge::{
    CohortLocusBuilderRegionsLen, DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    DEFAULT_MAX_COHORT_LOCUS_SPAN, MaxCohortLocusSpan,
};
use crate::ng::run::report::BoundsTheRunCalledUnder;
use crate::ng::run::{
    AlignedFilesVariantCaller, AlignmentInputs, MergeParameters, RunError, RunReport, Segmentation,
};
use crate::ng::tandem_repeat::{PeriodRange, PeriodRangeError};
use crate::ng::types::{Bp, DomainError, InbreedingF, MAX_MOTIF_LEN, Ploidy};
use crate::ng::vcf::header::{HeaderContig, HeaderMetadataError, VcfHeaderMetadata};
use crate::ng::vcf::writer::{VcfWriteError, VcfWriter};
use crate::pop_var_caller::common::current_command_line;
use crate::regions::{BedError, ContigBounds};

/// The ploidy a run assumes when it is not told.
///
/// **A property of the run and not of any fit** (`doc/devel/ng/spec/parameters_file.md` §3.2),
/// which is why it is a flag here rather than a number read out of a parameters file.
const DEFAULT_PLOIDY: u8 = 2;

/// **The tract-accuracy program's measurement switch for the per-locus slippage re-fit**
/// (lever L3, `doc/devel/ng/research/tract_accuracy_program_report.md`): how many re-fit
/// rounds a repeat tract may run, as an environment variable.
///
/// Absent or `0` is the frozen shipped default — no rounds, the loop untouched. A positive
/// whole number runs that many rounds at the configuration's default pull-backs
/// (`SlippageRefitConfig`: 50 pseudo-counts on the direction split and fall-off, 20 slipped
/// reads on the level); anything else is refused before a read is decoded, because a
/// measurement switch that fell back silently would report the frozen loop as the re-fitted
/// arm.
///
/// **An environment variable by design, not an oversight**: there is no parameters-file key
/// for this yet, deliberately — the arm is enabled per run while the program measures it, and
/// a *keep* verdict owes this switch proper parameters-file plumbing before the experiment
/// spelling is retired.
const NG_SLIPPAGE_REFIT_ROUNDS: &str = "NG_SLIPPAGE_REFIT_ROUNDS";

/// Zero both of the re-fit's pull-backs — the spec's **free** setting (HipSTR's), where the
/// locus's own reads set the numbers outright. `1` is the only accepted value; anything else
/// set is refused loudly. The same measurement-switch caveat as
/// [`NG_SLIPPAGE_REFIT_ROUNDS`] applies: an arm's spelling, owed real plumbing on a keep.
const NG_SLIPPAGE_REFIT_FREE: &str = "NG_SLIPPAGE_REFIT_FREE";

/// The calling-loop configuration this run asks for: the shipped defaults, with the slippage
/// re-fit's round count read once from [`NG_SLIPPAGE_REFIT_ROUNDS`] and its pull-backs
/// zeroed where [`NG_SLIPPAGE_REFIT_FREE`] asks for the free setting.
fn calling_loop_config_for_this_run()
-> Result<crate::ng::calling::inference::RunnableCallingLoopConfig, CallFromAlignmentsCliError> {
    let mut config = CallingLoopConfig::DEFAULT;
    match std::env::var(NG_SLIPPAGE_REFIT_ROUNDS) {
        Ok(rounds) => {
            config.slippage_refit.max_rounds = rounds.trim().parse().map_err(|error| {
                CallFromAlignmentsCliError::CallingLoopSettings(format!(
                    "{NG_SLIPPAGE_REFIT_ROUNDS} must be a whole number of re-fit rounds \
                     (0 is the frozen default), not {rounds:?}: {error}"
                ))
            })?;
        }
        Err(std::env::VarError::NotPresent) => {}
        Err(std::env::VarError::NotUnicode(value)) => {
            return Err(CallFromAlignmentsCliError::CallingLoopSettings(format!(
                "{NG_SLIPPAGE_REFIT_ROUNDS} is set but is not text: {value:?}"
            )));
        }
    }
    match std::env::var(NG_SLIPPAGE_REFIT_FREE) {
        Ok(value) if value.trim() == "1" => {
            config
                .slippage_refit
                .direction_and_fall_off_pull_back_pseudocounts = 0.0;
            config.slippage_refit.level_pull_back_slipped_reads = 0.0;
        }
        Ok(value) => {
            return Err(CallFromAlignmentsCliError::CallingLoopSettings(format!(
                "{NG_SLIPPAGE_REFIT_FREE} accepts only 1 (zero both pull-backs), \
                 not {value:?}"
            )));
        }
        Err(std::env::VarError::NotPresent) => {}
        Err(std::env::VarError::NotUnicode(value)) => {
            return Err(CallFromAlignmentsCliError::CallingLoopSettings(format!(
                "{NG_SLIPPAGE_REFIT_FREE} is set but is not text: {value:?}"
            )));
        }
    }
    config
        .validate()
        .map_err(|source| CallFromAlignmentsCliError::CallingLoopSettings(source.to_string()))
}

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

    /// The run's period range runs the wrong way round.
    ///
    /// Both ends are bounded to 1..=6 by clap, so this is reachable only by asking for a
    /// narrowest period wider than the widest.
    #[error("--min-period and --max-period do not make a range")]
    PeriodRange {
        /// Which way the range is wrong.
        #[source]
        source: PeriodRangeError,
    },

    /// This run asked for repeats the catalog was never built to hold.
    ///
    /// **Not a policy refusal**: the rows below the file's own floors were never written, so
    /// the request cannot be served at all. Either move the named flag back up, or build a
    /// catalog at floors low enough to answer it.
    #[error(
        "{flag} asks for repeats the catalog {} does not hold; raise it, or rebuild the catalog \
         at lower floors",
        path.display()
    )]
    RoutingBelowCatalog {
        /// The flag whose value put the request outside the file.
        flag: &'static str,
        /// The catalog that was asked.
        path: PathBuf,
        /// Which axis, with both numbers.
        #[source]
        source: RepeatCatalogError,
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

    /// A bound that is not a bound.
    #[error(
        "--max-cohort-locus-span {asked}: a locus covers at least one reference base, so a \
         bound of zero would refuse every locus there is"
    )]
    MaxCohortLocusSpanIsZero {
        /// What was asked for.
        asked: u32,
    },

    /// A round that covers no ground.
    #[error(
        "--cohort-locus-builder-regions-len {asked}: a round covers at least one reference \
         base, so a width of zero would never advance"
    )]
    CohortLocusBuilderRegionsLenIsZero {
        /// What was asked for.
        asked: u32,
    },

    /// An allele cap that is a refusal under another name.
    ///
    /// **Below two the reference is the only survivor** and every alternative becomes a
    /// truncation, so a locus carrying two obvious variants would lose both — which
    /// `candidate_alleles.md` §4.1 rules out.
    #[error(
        "--max-candidate-alleles {asked}: a cap counts the reference among the alleles, so it \
         needs at least {smallest} — the reference and one alternative"
    )]
    MaxCandidateAllelesTooSmall {
        /// What was asked for.
        asked: u16,
        /// The smallest cap that is a cap rather than a refusal.
        smallest: u16,
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

    /// The parameters file this run would write is the one it was handed.
    ///
    /// **Spec §7 invites the collision**: it tells a user to copy the file their run wrote and
    /// change a line, so `--parameters calls.parameters.toml --output calls.vcf.gz` is the
    /// natural next command and would write over the file that was just edited. Refused before
    /// anything is read, because what it destroys is the edit and the run would look ordinary.
    #[error(
        "--parameters {} is the file this run would write beside --output {}; \
         name the output differently, or copy the parameters somewhere else first",
        path.display(),
        output.display()
    )]
    ParametersWouldBeOverwritten {
        /// The file that was handed in.
        path: PathBuf,
        /// The output whose sibling it is.
        output: PathBuf,
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
/// How many observations one round may hold across the whole cohort, before the round's width
/// is narrowed to keep to it.
///
/// **The thing that costs memory is the product, not the width.** A round holds roughly one
/// observation per covered base per sample, so `width × samples` is what has to be bounded —
/// and bounding it means the observations a round holds are about the same at three samples
/// and at three thousand. At the ~600 bytes an observation runs to on the tomato benchmark,
/// half a million of them is about 300 MB.
const ROUND_OBSERVATION_BUDGET: u32 = 500_000;

/// The widest round the budget allows, in reference bases, whatever the cohort.
///
/// **A ceiling because the gain saturates, not because the memory does.** Measured on four
/// accessions over 400 kb of SL4.0: 3.41 s at 500 bases, 3.29 at 8,000, 3.18 at 32,000 and
/// 3.22 at 64,000 — the last two are level, and 64,000 costs 407 MB of peak resident against
/// 340. There is nothing to buy above this.
const WIDEST_ROUND: u32 = 16_000;

/// The round width a cohort of `samples` files gets when the command line names none.
///
/// **Why this is not one number.** The compiled-in default the merge carries
/// ([`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`], 500 bases) was chosen on the merge reading
/// pre-built `.psp` files, where a round costs a scan over records already in memory. A
/// calling run draws its records out of one CRAM per sample instead, and its rounds are where
/// that reading is overlapped across threads — so a narrow round pays a fan-out, a barrier
/// and a thread wake-up per sample for a few microseconds of work each time, and the waste
/// grows with the cohort.
///
/// Measured on the tomato benchmark, 63 accessions over the whole 8 Mb of SL4.0, 18 threads,
/// with the VCF byte-identical at every width:
///
/// | round width | wall | peak resident |
/// |---|---|---|
/// | 500 | 193.2 s | 1,593 MB |
/// | 2,000 | 157.8 s | 1,431 MB |
/// | 8,000 | 115.3 s | 1,811 MB |
///
/// The rule gives 7,936 bases at that cohort size, 500 at a thousand samples — which is the
/// merge's own default, so a cohort large enough to be memory-bound gets today's behaviour —
/// and the ceiling at anything under 32 samples.
///
/// **It never changes an answer.** Where a round's edge falls decides only which observations
/// are resident when a locus is built, not which loci exist or what they are called; the
/// widths above were checked against one another byte for byte.
fn round_width_for(samples: usize) -> CohortLocusBuilderRegionsLen {
    let samples = u32::try_from(samples).unwrap_or(u32::MAX).max(1);
    let width = (ROUND_OBSERVATION_BUDGET / samples)
        .clamp(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, WIDEST_ROUND);
    // PANIC-FREE: the clamp's lower bound is the merge's own default, which is non-zero.
    CohortLocusBuilderRegionsLen(
        NonZeroU32::new(width).expect("the clamp's floor is a non-zero constant"),
    )
}

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

    let asked_ploidy = ploidy_asked_for(args)?;
    let merge_parameters =
        MergeParameters {
            max_cohort_locus_span: MaxCohortLocusSpan(
                NonZeroU32::new(args.max_cohort_locus_span).ok_or(
                    CallFromAlignmentsCliError::MaxCohortLocusSpanIsZero {
                        asked: args.max_cohort_locus_span,
                    },
                )?,
            ),
            cohort_locus_builder_regions_len: match args.cohort_locus_builder_regions_len {
                Some(asked) => CohortLocusBuilderRegionsLen(NonZeroU32::new(asked).ok_or(
                    CallFromAlignmentsCliError::CohortLocusBuilderRegionsLenIsZero { asked },
                )?),
                None => round_width_for(args.alignments.len()),
            },
            ..MergeParameters::DEFAULT
        };
    let candidate_selection = CandidateSelectionConfig {
        max_candidate_alleles: MaxCandidateAlleles::new(args.max_candidate_alleles).ok_or(
            CallFromAlignmentsCliError::MaxCandidateAllelesTooSmall {
                asked: args.max_candidate_alleles,
                smallest: MaxCandidateAlleles::SMALLEST,
            },
        )?,
        ..CandidateSelectionConfig::DEFAULT
    };
    refuse_an_output_that_cannot_be_written(args)?;
    refuse_an_output_whose_parameters_file_is_this_run_s_input(args)?;

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
    // **Kept because the segmentation takes the original**, and the report names the ground the
    // run undertook to speak for — which is the one thing a VCF cannot say and the whole reason
    // a run reports at all.
    let analysed_for_the_report = analysed.clone();
    let segmentation = segments_over(args, &analysed, &with_checksums)?;

    let read_groups = build_read_groups(&args.alignments)
        .map_err(|source| CallFromAlignmentsCliError::ReadGroups { source })?;
    let numbers = run_parameters(
        args,
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
        .map_err(|source| CallFromAlignmentsCliError::ReferenceNotDigested { source })?;
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
        calling_loop_config_for_this_run()?,
        candidate_selection,
        merge_parameters,
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

    report(
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

/// **A run may not write its parameters file over the one it was handed.**
///
/// Spec §7 tells a user to copy the file their run wrote and change a line, so a re-run whose
/// supplied file and whose output share a stem writes over its own input — `calls.vcf.gz` and
/// `calls.parameters.toml` are exactly that pair. What it would destroy is whatever the person
/// changed by hand, and the numbers that came back would look ordinary, so this is refused
/// rather than done. `write_beside_the_vcf`'s own note leaves the choice to the driver; this is
/// the driver making it.
///
/// **Only a file that is there is compared**, so a mistyped `--parameters` gets the message
/// about the file not existing rather than an instruction to copy a file that does not exist.
///
/// **Both sides are resolved through the file system where it can resolve them**, which is what
/// makes the comparison see three things a textual one does not: `./calls.parameters.toml` and
/// `calls.parameters.toml` are one file; a symlink pointing at the file the run would write is
/// that file, and following it would have destroyed the target through the link; and on a
/// case-insensitive volume — macOS's default, and this project builds there — `CALLS.vcf.gz`'s
/// sibling is the same file as `calls.parameters.toml`. Where a path cannot be resolved, its
/// directory is resolved and the name compared as typed, which is the honest answer for a
/// destination that does not exist yet.
fn refuse_an_output_whose_parameters_file_is_this_run_s_input(
    args: &CallFromAlignmentsArgs,
) -> Result<(), CallFromAlignmentsCliError> {
    let Some(supplied) = &args.parameters else {
        return Ok(());
    };
    if !supplied.exists() {
        return Ok(());
    }
    let would_write = beside_the_vcf(&args.output);
    if resolved(&would_write) == resolved(supplied) {
        return Err(CallFromAlignmentsCliError::ParametersWouldBeOverwritten {
            path: supplied.clone(),
            output: args.output.clone(),
        });
    }
    Ok(())
}

/// A path as the file system knows it: fully resolved where it exists, and otherwise its
/// directory resolved with the name as typed. A name with no directory part is the working
/// directory's.
fn resolved(path: &Path) -> PathBuf {
    if let Ok(whole) = path.canonicalize() {
        return whole;
    }
    let directory = path.parent().unwrap_or_else(|| Path::new("."));
    let directory = if directory.as_os_str().is_empty() {
        Path::new(".")
    } else {
        directory
    };
    let resolved = directory
        .canonicalize()
        .unwrap_or_else(|_| directory.to_path_buf());
    match path.file_name() {
        Some(name) => resolved.join(name),
        None => resolved,
    }
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

/// **What this run counts as a repeat**, built from the five flags that say so.
///
/// The catalog is a source of *candidates*: it is deliberately built below every calling
/// floor, so that a caller can put its own line anywhere inside that gap by filtering the
/// file rather than re-scanning the genome (`repeat_catalog.md` §4.1). Until this existed the
/// run asked the file with [`StrRepeatCriteria::default()`], which *is* the file's own storage
/// floors — so everything the file held became an STR locus of the run, and on the human
/// benchmark that routed about seven times more reference to the repeat path than ng's own
/// calling floors would (`run_ssr_observations.md` §2).
///
/// **The flank floor is not a flag**, because it is the one axis a reader cannot move
/// downwards: the rows below the file's 15 bp were never written, so the request could not be
/// served. It comes from the conversion, which is where that reasoning lives
/// ([`StrRepeatCriteria::from`]).
///
/// The scan half of the [`TypedRegionConfig`] built here is unread: the conversion takes only
/// the classification rules and the satellite cap, and this run detects nothing — it reads
/// tracts a `repeat-catalog` run already found.
fn routing_criteria(
    args: &CallFromAlignmentsArgs,
) -> Result<StrRepeatCriteria, CallFromAlignmentsCliError> {
    // Both ends are already bounded to 1..=6 by clap, so the only way left to fail is a
    // range the wrong way round.
    let periods = PeriodRange::new(args.min_period, args.max_period)
        .map_err(|source| CallFromAlignmentsCliError::PeriodRange { source })?;
    Ok(StrRepeatCriteria::from(&TypedRegionConfig {
        max_str_len: Bp(args.max_str_len),
        criteria: SsrSegmentCriteria {
            periods,
            min_copies: args.min_copies,
            min_purity: args.min_purity,
            // Not a flag: the score floor gates the *scanner*'s output, and a catalog reader
            // has no scanner. `SsrSegmentCriteria::default()`'s 0 rejects nothing.
            ..SsrSegmentCriteria::default()
        },
        ..TypedRegionConfig::default()
    }))
}

/// Render a catalog failure, naming the flag to move when the failure is that this run asked
/// for more than the file holds.
///
/// **The refusal is real and it is not policy** — the rows below the file's floors were never
/// written, so the request cannot be served (`run_ssr_observations.md` §2.3) — but on its own
/// it names two numbers and no way to change either. A person who typed `--min-copies 3,3,3,3,3,3`
/// should be told that flag, not left to infer which of five knobs produced *"period 1: catalog
/// holds tracts of 5 copies and up, reader asked for 3"*.
///
/// **The flank floor has no flag and so has no arm here**: this run pins it at the catalog's
/// own, so a catalog built at a wider flank than 15 bp is a file that has to be rebuilt, which
/// is what the general catalog error already says.
fn catalog_error_naming_the_flag(
    source: RepeatCatalogError,
    path: &Path,
) -> CallFromAlignmentsCliError {
    // Exhaustive on the refusal on purpose: a new bounded axis must not silently inherit the
    // no-flag answer and leave a person hunting five knobs for the one they moved.
    let flag = match &source {
        RepeatCatalogError::CriteriaTooPermissive(refusal) => match refusal {
            CriteriaRefusal::CopyFloor { .. } => Some("--min-copies"),
            // Whichever end reaches outside what was built; `serves` checks the low end
            // first, so a range outside at both ends names `--min-period`.
            CriteriaRefusal::PeriodRange {
                built_min,
                wanted_min,
                ..
            } if wanted_min < built_min => Some("--min-period"),
            CriteriaRefusal::PeriodRange { .. } => Some("--max-period"),
            CriteriaRefusal::MinFlank { .. } => None,
        },
        _ => None,
    };
    match flag {
        Some(flag) => CallFromAlignmentsCliError::RoutingBelowCatalog {
            flag,
            path: path.to_path_buf(),
            source,
        },
        None => CallFromAlignmentsCliError::Catalog {
            path: path.to_path_buf(),
            source,
        },
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
    let criteria = routing_criteria(args)?;
    let catalog = RepeatCatalog::open_checking_against_reference(&path, with_checksums).map_err(
        |source| CallFromAlignmentsCliError::Catalog {
            path: path.clone(),
            source,
        },
    )?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))
        .map_err(|source| catalog_error_naming_the_flag(source, &path))?;
    Segmentation::build(
        segments,
        analysed.clone(),
        catalog.header().clone(),
        criteria,
        path,
    )
    .map_err(|source| CallFromAlignmentsCliError::Run { source })
}

/// **The numbers this run scores with, and the two things a file recording them needs beside
/// each.**
///
/// `RunParameters` keeps what *calling* reads — one bare multiplier a read group, one bare
/// coefficient a sample — and a file has to say how much data stood behind each and under what
/// warrant. Neither is recoverable from the parameters afterwards, so both travel from wherever
/// the numbers came from (`ParametersFile::of_run`'s own note says the same of the inbreeding
/// warrants).
struct TheRunsNumbers {
    /// What every locus is scored against.
    parameters: RunParameters,
    /// How many reads stood behind each read group's base-quality multiplier, in the run's dense
    /// read-group order.
    reads_behind_each_calibration: ReadsBehindEachCalibration,
    /// Each sample's inbreeding coefficient with its warrant, in the run's sample order.
    inbreeding_by_sample: Vec<Estimate<InbreedingF>>,
    /// **The census these numbers were fitted under, as the file that carried them named it** —
    /// naming no terms where they came from the defaults, which is how a run with no census
    /// spells itself.
    ///
    /// **Carried rather than restated, and `of_run`'s own contract is why**: *"a run that read a
    /// file fitted under other terms and writes its parameters out again has to write back the
    /// terms it read, not its own."* Direct mode never has a census of its own, so writing
    /// *there was none* would erase what the fit recorded — and a later psp run over the same
    /// cohort and the same census would then find a disagreement and demote every warrant to
    /// `supplied`, where reading the original file would have kept them. Provenance would
    /// degrade by one hop through direct mode, silently, which is exactly the divergence spec
    /// §2.1 exists to prevent.
    census: CensusIdentity,
}

/// The numbers this run scores with — a supplied file, or the defaults compiled in.
///
/// **A supplied file is bound to this run at its own door**: it is refused if it was fitted
/// against another reference or names read groups this cohort does not have, naming the
/// position where the two lists diverge (`doc/devel/ng/spec/parameters_file.md` §6). No census
/// is compared against, because direct mode has none — §2.1 settles that this keeps the file's
/// warrants rather than demoting them.
///
/// **A run that read a file writes back what it read.** The counts behind its multipliers and
/// the warrants on its coefficients belong to whichever run fitted them, and a run that dropped
/// them would write a file claiming its supplied numbers rest on nothing (spec §7).
fn run_parameters(
    args: &CallFromAlignmentsArgs,
    read_groups: &ReadGroups,
    with_checksums: &ReferenceInfo,
    ploidy: Ploidy,
    routing: &StrRepeatCriteria,
) -> Result<TheRunsNumbers, CallFromAlignmentsCliError> {
    let Some(path) = &args.parameters else {
        // **The defaults run's warrants are a pure function of what it was told**, which is what
        // lets these two be built beside the parameters rather than out of them:
        // `RunParameters::of_defaults` and `DeclaredInbreeding::of_each_sample` take the same
        // two arguments, so the coefficients here cannot disagree with the ones being scored.
        let inbreeding = DeclaredInbreeding::nothing_said();
        return Ok(TheRunsNumbers {
            parameters: RunParameters::of_defaults(read_groups, ploidy, &inbreeding),
            reads_behind_each_calibration: ReadsBehindEachCalibration::nothing_was_fitted(
                read_groups.len(),
            ),
            inbreeding_by_sample: inbreeding.of_each_sample(read_groups),
            // Nothing was fitted, so no census produced these numbers — and that is a fact
            // about them rather than a gap.
            census: CensusIdentity::of_a_run_with_no_census(),
        });
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
    // **A file written by a run that routed differently is reported and never refused** — the
    // owner's ruling, `run_ssr_observations.md` §2.3 and `parameters_file.md` §3.9. Which
    // stretches count as repeats is the user's to choose independently of where the numbers
    // came from, so this changes nothing about the run: the numbers are as warranted as they
    // ever were, and only the ground they are applied to has moved. What the operator is owed is
    // knowing it happened, because a tract this run admits and the fit's selection did not is
    // scored from strata fitted over other loci.
    if let Some(axis) = file.routing_disagreement(routing) {
        eprintln!(
            "note: {} was written by a run that routed on a different {axis}; calling on",
            path.display()
        );
    }
    let digest = ReferenceDigest::of(with_checksums)
        .map_err(|source| CallFromAlignmentsCliError::ReferenceNotDigested { source })?;
    let bound = file
        .to_run_parameters_for(&digest, read_groups, None)
        .map_err(|source| CallFromAlignmentsCliError::Parameters {
            path: path.clone(),
            source,
        })?;
    let from_file = bound.from_file;
    // **A flag that was typed may only agree with the file.** Spec §3.2 puts the ploidy in the
    // file so that a supplied one "cannot be paired with a run at a different ploidy without
    // saying so", and calling at the file's number while an operator typed another is exactly
    // that. Only a ploidy that was *typed* is compared: `--ploidy` is an `Option` for this
    // reason, so a tetraploid file is not refused for the flag's default being two.
    if let Some(asked) = args.ploidy
        && asked != from_file.parameters.ploidy().get()
    {
        return Err(CallFromAlignmentsCliError::PloidyIsNotTheParametersFiles {
            path: path.clone(),
            asked,
            in_the_file: from_file.parameters.ploidy().get(),
        });
    }
    Ok(TheRunsNumbers {
        parameters: from_file.parameters,
        reads_behind_each_calibration: ReadsBehindEachCalibration::as_a_file_recorded_them(
            from_file.reads_behind_each_calibration,
        ),
        inbreeding_by_sample: from_file.inbreeding_by_sample,
        census: file.fitted_from.census.clone(),
    })
}

/// What the file's header states about the run that wrote it.
///
/// **`##parametersFile` names the file this run writes beside its VCF**, by name and not by
/// path: the two are siblings by construction (`beside_the_vcf`), so a path would say the same
/// thing at greater length and would be wrong the moment somebody moved the pair. A reader who
/// wants to know what the genotypes rest on opens it and reads the line it opens with.
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
        beside_the_vcf(&args.output)
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned(),
    )
    .map_err(|source| CallFromAlignmentsCliError::Header { source })
}

/// **What the run has to say about itself**, printed when it finishes.
///
/// The report proper is [`RunReport`], which is where every rule about what a run owes a reader
/// lives; this adds the two paths — a person needs to know which files to open — and prints the
/// lines.
///
/// **The lines are the report's and not this function's**, so that what a run says is something
/// a test can hold. It was the one part of this command a mutation could change with the whole
/// suite still green.
fn report(calls: &Path, parameters_at: &Path, report: &RunReport<'_>) {
    println!("calls: {}", calls.display());
    println!("parameters: {}", parameters_at.display());
    for line in report.lines() {
        println!("{line}");
    }
}

#[cfg(test)]
mod tests;
