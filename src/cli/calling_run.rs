//! **What every calling run assembles, whichever the observations came from** — the flags a
//! person types turned into the numbers a run scores with, and the refusals that judge them.
//!
//! # Why this is a module and not two copies
//!
//! There are two calling commands, `call-from-alignments` and `call-from-psps`, and spec
//! §12.3 asks them for **the same VCF over the same cohort**. Everything between the command
//! line and the calling loop is therefore not two decisions that happen to agree — it is one
//! decision, and a second copy of it is a place for the two modes to drift apart while both
//! keep running. Three of the pieces here would drift silently:
//!
//! - **the parameters** ([`run_parameters`]): which numbers a run scores with, and the
//!   refusals that bind a supplied file to this cohort and this reference;
//! - **the calling loop's settings** ([`calling_loop_config_for_this_run`]): the measurement
//!   switches the tract-accuracy program reads from the environment, so a mode that did not
//!   read one would report the frozen loop as the switched arm;
//! - **the round width** ([`round_width_for`]): a run parameter chosen from the cohort's size,
//!   which changes no answer but which the byte-comparison oracle would notice at once.
//!
//! This module is [`run_ground`](super::run_ground)'s sibling: that one owns the ground a run
//! speaks for, this one owns the numbers it speaks with. Both were lifted out of
//! `call-from-alignments` when the second command needed them, and both keep their refusals
//! `#[error(transparent)]` at the commands, so one mistake reads the same however it was
//! reached.

use std::num::{NonZeroU32, NonZeroUsize};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::calling::allele_candidates::{CandidateSelectionConfig, MaxCandidateAlleles};
use crate::calling::inference::{CallingLoopConfig, DiscoveryMode, RunnableCallingLoopConfig};
use crate::calling::likelihood::MAX_PLOIDY_COPIES;
use crate::calling::parameters_file::{
    CensusIdentity, DeclaredInbreeding, FittedDensity, FittedDuplicatedPositions, ParametersFile,
    ParametersFileError, ReadsBehindEachCalibration, RepeatTractsLeftOut, SnpIndelFitAsReturned,
    SnpIndelFitStart,
};
use crate::calling::run_parameters::RunParameters;
use crate::cli::provenance::current_command_line;
use crate::cli::run_output;
use crate::fasta::ContigList;
use crate::neighbour_evidence::phase::PathAgreement;
use crate::neighbour_evidence::stage::{
    NeighbourEvidenceDecision, NeighbourEvidenceSwitch, decide,
};
use crate::parameter_estimation::Estimate;
use crate::parameter_estimation::joint::loci::{ReferenceDigest, SelectionError};
use crate::read::input::read_groups::ReadGroups;
use crate::reference_info::ReferenceInfo;
use crate::repeat_catalog::StrRepeatCriteria;
use crate::run::cohort_merge::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen,
    DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, MaxCohortLocusSpan,
};
use crate::run::{MergeParameters, RunReport};
use crate::types::{DomainError, InbreedingF, Ploidy};
use crate::vcf::FinishedVcf;
use crate::vcf::header::{HeaderContig, HeaderMetadataError, VcfHeaderMetadata};

/// Everything a calling run can refuse before it reads an observation, whichever mode it is.
///
/// **Carried `#[error(transparent)]` by both commands' own error enums**, so a mistake in what
/// a person typed renders in one wording however it was reached — the rule
/// [`GroundError`](super::run_ground::GroundError) already follows for the ground.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CallingRunError {
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
        /// The most copies the read likelihood scores.
        limit: usize,
    },

    /// A bound that is not a bound.
    #[error(
        "--max-locus-span {asked}: a locus covers at least one reference base, so a \
         bound of zero would refuse every locus there is"
    )]
    MaxCohortLocusSpanIsZero {
        /// What was asked for.
        asked: u32,
    },

    /// A round that covers no ground.
    #[error(
        "--locus-builder-regions-len {asked}: a round covers at least one reference \
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
        "--max-generic-alleles {asked}: a cap counts the reference among the alleles, so it \
         needs at least {smallest} — the reference and one alternative"
    )]
    MaxCandidateAllelesTooSmall {
        /// What was asked for.
        asked: u16,
        /// The smallest cap that is a cap rather than a refusal.
        smallest: u16,
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
    /// otherwise be discovered after every input file had been opened.
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

    /// `--tmp-dir` names something that is not a directory.
    ///
    /// Judged with the rest of the command line, because the scratch directory is otherwise first
    /// made after the reference has been read and every input opened.
    #[error("--tmp-dir {} is not a directory", path.display())]
    TmpDirIsNotADirectory {
        /// What was typed.
        path: PathBuf,
    },

    /// The run's scratch directory could not be made.
    #[error(
        "the run's scratch directory could not be made in {} — give --tmp-dir a directory it \
         can write to",
        directory.display()
    )]
    ScratchDirectory {
        /// Where it was to be made.
        directory: PathBuf,
        /// What the file system said.
        #[source]
        source: std::io::Error,
    },

    /// The scratch directory a resumed run was pointed at is not there.
    #[error(
        "`{}` is not a directory: point --scratch-dir at the scratch directory a failed run kept",
        directory.display()
    )]
    ScratchDirectoryNotFound {
        /// What was given.
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
    /// has nothing to compare a file's reference binding against. Both calling commands read
    /// the FASTA, so it is unreachable through either command line and is carried rather than
    /// unwrapped.
    #[error("this run's reference carries no digest, so a parameters file cannot be bound to it")]
    ReferenceNotDigested {
        /// What the digest said.
        #[source]
        source: SelectionError,
    },

    /// The output's header could not honestly state what this run is.
    #[error("the output's header could not be written")]
    Header {
        /// What the header refused.
        #[source]
        source: HeaderMetadataError,
    },

    /// The run's shipped calling-loop settings would not validate.
    ///
    /// **Not reachable from any flag** — the settings are compiled in — so this is a defect in
    /// this binary rather than in what was typed. It is an error and not a panic because a
    /// message naming the setting is worth more than a backtrace.
    #[error("this binary's calling-loop settings are not runnable: {0}")]
    CallingLoopSettings(String),
}

/// The ploidy a run assumes when it is not told.
///
/// **A property of the run and not of any fit** (`doc/devel/ng/spec/parameters_file.md` §3.2),
/// which is why it is a flag here rather than a number read out of a parameters file.
pub const DEFAULT_PLOIDY: u8 = 2;

/// **A measurement switch for STR allele discovery** (lever L7). The discovery pre-pass runs
/// at every STR segment by default, admitting the sequences one sample showed too often for
/// slippage to explain — 2 reads **and** 15% of that sample's spanning reads.
///
/// Absent or `1` is that default, [`DiscoveryMode::BeforeTheLoop`]; `0` switches it off, for
/// measuring. Anything else is refused before a read is decoded, because a measurement switch
/// that fell back silently would report one arm as the other.
///
/// **Not a user option, by the owner's decision (2026-11-23)**: discovery is fixed behaviour,
/// documented in `doc/user/methods.md`; this variable exists only to measure without it.
const NG_TRACT_DISCOVERY: &str = "NG_TRACT_DISCOVERY";

/// Default target false-discovery rate for the hidden-duplication filter — spec §3.6's.
///
/// **About one record in a hundred of those the filter removes was really a variant.** It was
/// `0` while the filter's scoring and writing passes did not exist, because that default would
/// have parked every record and produced no VCF; both passes are built now, so the spec's value
/// stands and the filter is on unless an operator turns it off with `--hidden-dup-fdr 0`.
const DEFAULT_PARALOG_FDR: f64 = 0.01;

// **The flags `call-from-alignments` and `call-from-psps` share**, declared once and flattened
// into both. Two structs rather than one because `--help` lists a command's flags in the order
// they are declared: the shared flags sit on either side of `--max-depth`, whose text differs
// between the commands, so each half is flattened where it stood. A new shared flag goes into
// one of these, and the tests' argument values come from `for_tests` below.
//
// (Plain comments rather than doc comments: clap reads a flattened struct's doc comment as the
// command's description.)
#[derive(Debug, clap::Args, Clone)]
pub struct CallingFlags {
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

    /// Charge a site whose variant reads carry more low base qualities than its reference reads,
    /// at this many times the test's tail, subtracted from QUAL (`BQPEN`). 0 leaves the test off
    /// (`doc/devel/ng/spec/low_base_quality.md`).
    #[arg(
        long,
        default_value_t = crate::calling::quality::artifact_correction::DEFAULT_LOW_BASE_QUALITY_WEIGHT,
        value_parser = crate::cli::calling_run::parse_low_base_quality_weight
    )]
    pub low_base_quality_weight: f64,

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
    /// it ends, whether it succeeded or failed, unless --keep-tmp-on-error; a run that is
    /// killed, or stops on an internal error, leaves it behind. While the hidden-dup filter is
    /// on, every call waits there until the whole genome is called: several times the
    /// compressed VCF's size.
    #[arg(long)]
    pub tmp_dir: Option<PathBuf>,

    /// Keep the run's scratch directory, and the hidden-dup filter's spill in it, when the run
    /// stops on an error. Its path is logged. A run that finishes removes it as usual. A run
    /// stopped by a refinement window over --refine-window-bytes, after calling, keeps it with or
    /// without this flag, so `pvc resume-calling` can finish it.
    ///
    /// The spill holds every called record, so on a large cohort the kept directory can take
    /// terabytes: remove it once it is no longer needed.
    #[arg(long)]
    pub keep_tmp_on_error: bool,
}

// The advanced flags both calling commands share; see `CallingFlags`.
#[derive(Debug, clap::Args, Clone)]
pub struct AdvancedCallingFlags {
    /// The most alleles a SNP/indel locus is called over, the reference among them.
    ///
    /// Where more segregate, the worst supported are left out. Repeat (STR) loci are not
    /// affected: they have their own limit of 32.
    #[arg(
        long,
        default_value_t = crate::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        help_heading = "Advanced"
    )]
    pub max_generic_alleles: u16,

    /// The widest locus assembled, in reference bases — and so the longest deletion called.
    ///
    /// The run's log says how many loci were too wide and how long they were. Long reads want a
    /// larger value than short ones.
    #[arg(long, default_value_t = crate::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN, help_heading = "Advanced")]
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
    /// the genome-wide rate where few calls are near.
    #[arg(
        long,
        value_name = "BASES",
        default_value_t = crate::run::paralog_filter::DEFAULT_HIDDEN_DUP_WINDOW_BASES,
        help_heading = "Advanced"
    )]
    pub hidden_dup_window: u64,

    /// Whether to re-call each sample's genotypes with the evidence of the samples it shares
    /// stretches of chromosome with, and of its own neighbouring sites. Off by default; `auto`
    /// runs it on a diploid cohort of at least 10 samples whose median inbreeding is at least
    /// 0.5, and `on` on any diploid cohort of two or more. The VCF header says whether it ran
    /// and why. Where it does not run, on a diploid cohort, phase is
    /// written from reads alone (`GT` `a|b` with a `PS`); the header's `##readBackedPhase=` line
    /// states its thresholds.
    #[arg(
        long,
        value_enum,
        default_value_t = RefineGenotypesArg::Off,
        help_heading = "Advanced"
    )]
    pub refine_genotypes: RefineGenotypesArg,

    /// The share of the neighbour-evidence pass's drawn paths that must agree on how two
    /// neighbouring heterozygous sites are phased for the VCF to phase them together. 1, all of
    /// them, is the strictest; lower values give longer phase blocks with more wrong links.
    /// The default, 0.8, was measured to keep wrong links under 1 in 200 from 12 samples up.
    #[arg(
        long,
        default_value = "0.8",
        value_parser = parse_path_agreement,
        help_heading = "Advanced"
    )]
    pub phase_min_path_agreement: crate::neighbour_evidence::phase::PathAgreement,

    /// The most memory, in bytes, one window of genotype refinement may hold: its records' lines,
    /// numbers and reads, and each sample's state while it is solved. Windows shrink to fit it,
    /// down to --refine-window-floor; a cohort whose smallest window still does not fit stops
    /// the run with an error saying how much it needs. Each refinement thread holds, besides, its
    /// walk's arrays, about 1 kB a record it walks, and the memory allocator holds about a sixth
    /// more: on 2,176 samples at 48 threads the process peaked at about 1.3 times the budget. The
    /// run estimates that peak before calling and at each window, and warns, without stopping,
    /// where it is above 85% of the memory the run may use (the machine's, or the job's limit
    /// where lower, such as a systemd unit's MemoryMax), naming a budget that would fit.
    #[arg(
        long,
        default_value_t = crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES,
        help_heading = "Advanced"
    )]
    pub refine_window_bytes: u64,

    /// The shortest each margin of a genotype refinement window may become when windows shrink to
    /// fit --refine-window-bytes, in bases; the core shrinks with them. Shorter margins hold less
    /// memory but may cut off evidence: a sample copies stretches of up to about 300 kb.
    #[arg(
        long,
        value_name = "BASES",
        default_value_t = crate::neighbour_evidence::windows::DEFAULT_MARGIN_FLOOR_BASES,
        value_parser = clap::value_parser!(u64).range(1..=crate::neighbour_evidence::windows::DEFAULT_MARGIN_BASES),
        help_heading = "Advanced"
    )]
    pub refine_window_floor: u64,

    /// The shortest a genotype refinement window's core may become, in bases, when a window with
    /// margins at --refine-window-floor is still over --refine-window-bytes. The margins keep their
    /// length, so no record loses the evidence around it; a shorter core only makes more windows,
    /// each solving its margins again. Defaults to --refine-window-floor.
    #[arg(
        long,
        value_name = "BASES",
        value_parser = clap::value_parser!(u64).range(1..=crate::neighbour_evidence::windows::DEFAULT_MARGIN_BASES),
        help_heading = "Advanced"
    )]
    pub refine_core_floor: Option<u64>,

    /// Walk at most this many records a kb, evenly spaced, in each genotype refinement window's
    /// margins farther than 50 kb from its core; 0 walks every record. Stretches no denser than
    /// that are walked in full. The core and the 50 kb next to it are always walked in full, so
    /// every record is refined as before but sees fewer of the far records when choosing which
    /// samples it resembles; the records not walked are held compressed. A window then holds
    /// less memory and takes less time: on 2,176 tomato samples, 5 a kb held 41% less and took
    /// half the time, with as many genotypes wrong.
    #[arg(
        long,
        value_name = "RECORDS",
        default_value = crate::neighbour_evidence::windows::DEFAULT_MARGIN_DENSITY,
        value_parser = parse_margin_density,
        help_heading = "Advanced"
    )]
    pub refine_margin_density: MarginDensity,

    /// The longest a refinement window's core and each margin may be, in bases, before the budget
    /// shrinks them (at most 1 Mb, the default margin; the default core is 2 Mb). For measurements
    /// on short contigs, where a default window would span the contig with no margin.
    #[arg(long, value_name = "BASES", hide = true, value_parser = clap::value_parser!(u64).range(1..=crate::neighbour_evidence::windows::DEFAULT_MARGIN_BASES))]
    pub refine_window_max: Option<u64>,

    /// The longest a refinement window's core may be, in bases, before the budget shrinks it. For
    /// measuring how core length alone changes the output.
    #[arg(long, value_name = "BASES", hide = true, value_parser = clap::value_parser!(u64).range(1..=crate::neighbour_evidence::windows::DEFAULT_MARGIN_BASES))]
    pub refine_core_max: Option<u64>,
}

#[cfg(test)]
impl CallingFlags {
    /// The shared flags at their defaults, for a test to override the few it is about.
    pub(crate) fn for_tests() -> Self {
        Self {
            ploidy: None,
            min_site_quality: crate::run::records::DEFAULT_MIN_SITE_QUALITY,
            low_base_quality_weight:
                crate::calling::quality::artifact_correction::DEFAULT_LOW_BASE_QUALITY_WEIGHT,
            threads: 0,
            force: false,
            log: None,
            tmp_dir: None,
            keep_tmp_on_error: false,
        }
    }
}

#[cfg(test)]
impl AdvancedCallingFlags {
    /// The advanced flags as the tests hand them to a command: the hidden-duplication filter off,
    /// refinement on `auto`, every drawn path required to agree before phasing, and the window
    /// settings at their defaults.
    pub(crate) fn for_tests() -> Self {
        Self {
            max_generic_alleles: crate::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES
                .get(),
            max_locus_span: crate::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN,
            hidden_dup_fdr: 0.0,
            hidden_dup_tag: false,
            hidden_dup_window: 0,
            refine_genotypes: RefineGenotypesArg::Auto,
            phase_min_path_agreement: crate::neighbour_evidence::phase::PathAgreement::ALL,
            refine_window_bytes: crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES,
            refine_window_floor: crate::neighbour_evidence::windows::DEFAULT_MARGIN_FLOOR_BASES,
            refine_core_floor: None,
            refine_margin_density: MarginDensity(None),
            refine_window_max: None,
            refine_core_max: None,
        }
    }
}

/// **The run's own settings, from the flags both calling commands share** — built once here so the
/// two commands cannot configure a run two ways. Fragments are gathered, the default; the command
/// turns that off once its header says the output reads none
/// ([`with_fragments_gathered`](crate::run::callers::AlignedFilesVariantCaller::with_fragments_gathered)).
pub fn calling_settings(
    max_depth: u32,
    explain: Option<&crate::run::explain::ExplainRegions>,
    min_site_quality: f32,
    low_base_quality_weight: f64,
) -> crate::run::callers::CallingSettings {
    crate::run::callers::CallingSettings {
        depth_ceiling: crate::run::depth_ceiling::DepthCeiling::at(max_depth),
        explain: explain.cloned(),
        min_site_quality: crate::types::Phred::try_new(min_site_quality)
            .expect("`parse_min_site_quality` admits only valid qualities"),
        low_base_quality_weight,
        ..crate::run::callers::CallingSettings::default()
    }
}

/// The calling-loop configuration this run asks for: the shipped defaults, with the discovery
/// pre-pass switched off where [`NG_TRACT_DISCOVERY`] asks for it.
pub fn calling_loop_config_for_this_run() -> Result<RunnableCallingLoopConfig, CallingRunError> {
    let mut config = CallingLoopConfig::DEFAULT;
    match std::env::var(NG_TRACT_DISCOVERY) {
        // The pre-pass has been the shipped default since the owner adopted L7 (2026-09-03),
        // so `1` restates the default and `0` is the measurement arm that switches it off.
        Ok(value) if value.trim() == "1" => {
            config.discovery.mode = DiscoveryMode::BeforeTheLoop;
        }
        Ok(value) if value.trim() == "0" => {
            config.discovery.mode = DiscoveryMode::Off;
        }
        Ok(value) => {
            return Err(CallingRunError::CallingLoopSettings(format!(
                "{NG_TRACT_DISCOVERY} accepts 1 (the discovery pre-pass, the shipped \
                 default) or 0 (off, the measurement arm), not {value:?}"
            )));
        }
        Err(std::env::VarError::NotPresent) => {}
        Err(std::env::VarError::NotUnicode(value)) => {
            return Err(CallingRunError::CallingLoopSettings(format!(
                "{NG_TRACT_DISCOVERY} is set but is not text: {value:?}"
            )));
        }
    }
    config
        .validate()
        .map_err(|source| CallingRunError::CallingLoopSettings(source.to_string()))
}

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
pub fn round_width_for(samples: usize) -> CohortLocusBuilderRegionsLen {
    let samples = u32::try_from(samples).unwrap_or(u32::MAX).max(1);
    let width = (ROUND_OBSERVATION_BUDGET / samples)
        .clamp(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, WIDEST_ROUND);
    // PANIC-FREE: the clamp's lower bound is the merge's own default, which is non-zero.
    CohortLocusBuilderRegionsLen(
        NonZeroU32::new(width).expect("the clamp's floor is a non-zero constant"),
    )
}

/// **How a round is cut up so that it holds what one region holds today** — the width of each
/// building region and how many of them are worked at once.
///
/// # Why the two numbers are one decision
///
/// A round holds every sample's observations over `regions_in_flight × width` bases, where the
/// one-region-at-a-time driver holds them over `width`. So asking for more regions in flight
/// without narrowing them multiplies the memory that [`round_width_for`] was written to bound —
/// silently, and by the number of cores the machine happens to have. **The rule here is that
/// the round holds the ground one region holds today**: `ground` is that width — the flag's if
/// one was typed, [`round_width_for`]'s otherwise — and this divides it.
///
/// # How narrow, and how many, reading stored files
///
/// **A run over stored files streams its rounds**
/// ([`stream`](crate::run::cohort_merge::stream)), holding
/// [`ROUNDS_HELD`](crate::run::cohort_merge::stream::ROUNDS_HELD) of them at once, so each round
/// gets that share of `ground` and the run holds about `ground`, as before.
///
/// **Up to [`REGIONS_A_WORKER`] regions for every thread**, so that a thread whose region finishes
/// early takes another rather than waiting for the round's slowest. A round cannot start until the
/// one before has finished, so with one region a thread the round lasts as long as its slowest
/// region; with four it lasts about as long as the average thread's share.
///
/// **No narrower than [`narrowest_round_region`]**, which keeps a region's work above the cost of
/// handing it out: a region's work grows with its width times the samples, so a one-sample cohort
/// keeps the merge's own 500 bases and a cohort of eight or more goes down to 125.
///
/// Measured on the tomato25 cohort's stored files at about three reads a position, 48 threads, the
/// VCF byte-identical at every setting (the regions as narrow as 500 bases against 125, the
/// ground unchanged):
///
/// | samples | ground | 500 bases | 125 bases | peak resident |
/// |---|---|---|---|---|
/// | 1 (4 Mb) | 16,000 | 5 s | 7 s | 1.5 GB either way |
/// | 8 (4 Mb) | 16,000 | 22–24 s | 20 s | 3.9–4.2 GB |
/// | 63 (4 Mb) | 7,936 | 42–43 s | 22–32 s | 1.3 GB |
/// | 2,176 (300 kb) | 24,000 | 69 s | 58 s | 21.1 against 21.6 GB |
///
/// # Reading alignments: one region a thread, no narrower than 500 bases
///
/// **The narrower regions do not pay there.** At 48 tomato accessions over 200 kb of SL4.0
/// chromosome 1, 48 threads, 125-base regions called in 76 s against 65 s at 500, with 30% more
/// CPU, the VCF the same. A region there decodes its own reads and pays a start-up for it, which
/// a region over stored files does not. So `call-from-alignments` keeps the rule measured on the
/// tomato benchmark, 1 to 63 accessions: peak resident up 4% to 13% and the calling pass up to
/// 5.5× faster than one region at a time.
///
/// # Large cohorts
///
/// **Above about a thousand samples, with no width typed, there is little to divide**:
/// `round_width_for` gives 500 bases of ground there, which this cuts into four regions of 125.
/// A large cohort keeps the memory it has and gives up most of the parallelism, rather than keeping
/// the parallelism and multiplying the memory by the core count. A run that can afford more ground
/// asks for it with `--locus-builder-regions-len`, as the 2,176-sample run above did.
#[must_use]
pub fn round_shape_for(
    ground: CohortLocusBuilderRegionsLen,
    samples: usize,
    read_from: ObservationsReadFrom,
) -> (
    CohortLocusBuilderRegionsLen,
    CohortLocusBuilderRegionsInFlight,
) {
    let (ground, floor, regions_a_worker) = match read_from {
        // The stream holds `ROUNDS_HELD` rounds at once, so each gets that share of the ground.
        ObservationsReadFrom::StoredFiles => (
            ground.get() / crate::run::cohort_merge::stream::ROUNDS_HELD as u32,
            narrowest_round_region(samples),
            REGIONS_A_WORKER,
        ),
        ObservationsReadFrom::Alignments => {
            (ground.get(), DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, 1)
        }
    };
    let workers = rayon::current_num_threads().max(1);
    // PANIC-FREE: the clamp's floor is one, so the division is by at least one, and `workers` is
    // at least one.
    let in_flight = usize::try_from(ground / floor)
        .unwrap_or(usize::MAX)
        .clamp(1, workers.saturating_mul(regions_a_worker));
    let width = (ground / u32::try_from(in_flight).unwrap_or(u32::MAX)).max(floor);
    (
        CohortLocusBuilderRegionsLen(
            NonZeroU32::new(width).expect("the max above is a non-zero constant"),
        ),
        CohortLocusBuilderRegionsInFlight(
            NonZeroUsize::new(in_flight).expect("the clamp's floor is one"),
        ),
    )
}

/// **Where a calling run reads its cohort's observations from**, which decides how its rounds are
/// cut ([`round_shape_for`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationsReadFrom {
    /// One alignment file per sample, decoded as the merge goes (`call-from-alignments`).
    Alignments,
    /// One stored psp file per sample, written beforehand (`call-from-psps`).
    StoredFiles,
}

/// **How many regions a round over stored files gives each thread, at most** ([`round_shape_for`]). At 2,176 tomato
/// samples over 300 kb, one region a thread called in 69 s, two in 61 s, four in 58 s and eight in
/// 57 s, the ground the same each time.
pub const REGIONS_A_WORKER: usize = 4;

/// **The work a building region carries at least, in sample-bases** — its width times the
/// cohort's samples ([`narrowest_round_region`]).
const SAMPLE_BASES_A_REGION: u32 = 1_000;

/// The narrowest region a round of a cohort of `samples` is cut into: [`SAMPLE_BASES_A_REGION`]
/// over the samples, between 125 bases and the merge's own 500
/// ([`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`]).
///
/// **Why it depends on the cohort.** Handing a region out costs about the same whatever is in it,
/// and what is in it grows with the samples. At one sample, 125-base regions called 4 Mb in 7 s
/// against 5 s at 500; at eight they called it in 20 s against 22–24. So one or two samples keep
/// 500 bases, four get 250, and eight or more 125. Between one and eight the rule is interpolated,
/// not measured.
#[must_use]
pub fn narrowest_round_region(samples: usize) -> u32 {
    let samples = u32::try_from(samples.max(1)).unwrap_or(u32::MAX);
    (SAMPLE_BASES_A_REGION / samples).clamp(
        NARROWEST_ROUND_REGION,
        DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    )
}

/// The narrowest a building region is cut, whatever the cohort ([`narrowest_round_region`]): at
/// 2,176 samples, 64-base regions called 300 kb in 57 s against 58 s at 125, no faster.
const NARROWEST_ROUND_REGION: u32 = 125;

/// **The fewest checked sites a sample's fill disagreement rate is stated from**: below it, the
/// table gives the counts and the log counts the sample as not measured. At 30 sites a rate of one
/// in ten is known to within about a factor of two.
pub const MIN_FILL_CHECK_SITES: u32 = 30;

/// **Reports what reads did to the phase** (read phase spec §2.6, §2.8): phase from reads alone where
/// it ran, and where the neighbour-evidence pass ran, the read check of its population links — in
/// total in the log, and per sample at debug level.
pub fn report_read_phase(finished: &FinishedVcf) {
    if let Some(tally) = finished.read_phase_tally {
        tracing::info!(
            "phase from reads alone: {} sample columns written phased, in {} blocks of two records or \
             more; {} links made by reads, {} refused because a read over an earlier record of the \
             block contradicted them",
            run_output::with_thousands(tally.columns_phased),
            run_output::with_thousands(tally.blocks),
            run_output::with_thousands(tally.links),
            run_output::with_thousands(tally.links_refused_by_an_earlier_record),
        );
        if tally.lines_not_read > 0 {
            tracing::warn!(
                "phase from reads alone: {} lines could not be read and were written without PS; \
                 every line comes from this program's encoder, so this is a defect to report",
                run_output::with_thousands(tally.lines_not_read),
            );
        }
    }
    if finished.pass_read_counts.is_empty() {
        return;
    }
    let mut total = crate::neighbour_evidence::phase::PassReadCounts::default();
    for (_, counts) in &finished.pass_read_counts {
        total += *counts;
    }
    tracing::info!(
        "reads in the neighbour-evidence pass: the read check — of the population links reads \
         reach at odds of {:.0} to 1 with GQ {} or more at both ends, and with one end among the \
         records whose reads the haplotype draw held out, {} of {} contradicted by the \
         reads, a check of the population links that does not come from the model, for links \
         under {} bases and not acted on; {} population links broken by reads, {} links made by \
         reads, {} heterozygotes the consensus did not lean on attached by reads",
        crate::float::exp(crate::read_phase::link::LINK_MIN),
        crate::read_phase::link::GQ_MIN,
        total.population_links_contradicted,
        total.population_links_reached,
        crate::read_phase::fragments::READ_REACH,
        total.population_links_broken,
        total.links_made_by_reads,
        total.records_attached_by_reads,
    );
    for (sample, counts) in &finished.pass_read_counts {
        tracing::debug!(
            "  {sample}: read check {} of {} population links contradicted; {} broken, {} made, {} \
             attached by reads",
            counts.population_links_contradicted,
            counts.population_links_reached,
            counts.population_links_broken,
            counts.links_made_by_reads,
            counts.records_attached_by_reads,
        );
    }
}

/// **Reports what the neighbour-evidence pass did**, when it ran: its counts in the log, and each
/// sample's fill check (`FillCheck`) — summarised in the log, and in full in a table beside the VCF,
/// `<vcf name>.fill_check.tsv`, unless the VCF went to standard output. A table that cannot be
/// written is a warning, not a failed run: the VCF is already whole.
pub fn report_neighbour_evidence(finished: &FinishedVcf, vcf_path: &Path) {
    let Some(tally) = finished.neighbour_evidence else {
        return;
    };
    tracing::info!(
        "neighbour-evidence pass: {} records re-called; {} genotypes called ({} of them where the \
         sample had no reads), {} left ./. where the sample had no reads, {} changed; {} alleles \
         and {} records dropped as no longer called; {} repeat counts unknown",
        run_output::with_thousands(tally.records_recalled),
        run_output::with_thousands(tally.called),
        run_output::with_thousands(tally.called_without_reads),
        run_output::with_thousands(tally.withheld),
        run_output::with_thousands(tally.genotypes_changed),
        run_output::with_thousands(tally.alleles_dropped),
        run_output::with_thousands(tally.records_dropped),
        run_output::with_thousands(tally.repeat_counts_unknown),
    );

    let mut rates: Vec<(f64, &str)> = finished
        .fill_checks
        .iter()
        .filter(|(_, check)| check.checked >= MIN_FILL_CHECK_SITES)
        .map(|(sample, check)| {
            (
                f64::from(check.disagreed) / f64::from(check.checked),
                sample.as_str(),
            )
        })
        .collect();
    rates.sort_by(|left, right| left.0.total_cmp(&right.0));
    let table = fill_check_table_path(vcf_path);
    let table_note = table.as_ref().map_or(String::new(), |path| {
        format!("; each sample in {}", path.display())
    });
    let unmeasured = finished.fill_checks.len() - rates.len();
    if rates.is_empty() {
        tracing::info!(
            "fill check: no sample had {MIN_FILL_CHECK_SITES} sites where its reads settle its \
             genotype and its neighbours would have filled it, so no fill error rate is \
             estimated{table_note}"
        );
    } else {
        let percent = |rate: f64| format!("{:.1}%", 100.0 * rate);
        let worst: Vec<String> = rates
            .iter()
            .rev()
            .take(3)
            .map(|&(rate, sample)| format!("{sample} {}", percent(rate)))
            .collect();
        tracing::info!(
            "fill check — each sample's no-read calls tried at records hidden from the phasing \
             rounds where its reads settle its genotype, the neighbours' prediction against the \
             reads: {} of {} samples checked \
             at {MIN_FILL_CHECK_SITES} sites or more, wrong at a median of {} of sites (from {} to \
             {}; highest {}); {unmeasured} with fewer sites, not estimated. A best case: a site no \
             read reached can be harder{table_note}",
            rates.len(),
            finished.fill_checks.len(),
            percent(rates[rates.len() / 2].0),
            percent(rates[0].0),
            percent(rates[rates.len() - 1].0),
            worst.join(", "),
        );
    }
    warn_of_poor_fills(&finished.fill_checks, table_note.as_str());
    if let Some(path) = table
        && let Err(error) = write_fill_check_table(&path, &finished.fill_checks)
    {
        tracing::warn!("fill check: could not write {}: {error}", path.display());
    }
}

/// **The fill error rate above which a run warns** (owner, 2026-11-02): over the whole cohort, or
/// for a sample whose rate is plausibly above it ([`FillCheck::wrong_share_at_least`](crate::neighbour_evidence::refine::FillCheck::wrong_share_at_least)).
pub const POOR_FILL_SHARE: f64 = 0.02;

/// **Warns when the pass's fills look poor**, so the user can look at why: the cohort's fills, all
/// samples' checks pooled, wrong at more than [`POOR_FILL_SHARE`]; or samples whose own rate is
/// above it beyond the chance of their few checked sites. On the founder simulator the pooled rule
/// warned on both outbred cohorts (estimated 4.2% and 5.5%) and on no inbred one; the per-sample rule
/// named 8 and 6 samples of 50 in the outbred cohorts and one sample in all the inbred ones, where a
/// plain share above 2% named up to 10 samples of 50 in a cohort whose true fill error was 1.5%. A
/// run where no sample was checked — at 1× — says nothing here.
fn warn_of_poor_fills(
    checks: &[(String, crate::neighbour_evidence::refine::FillCheck)],
    table_note: &str,
) {
    let (checked, disagreed) = checks.iter().fold((0_u64, 0_u64), |(n, d), (_, check)| {
        (n + u64::from(check.checked), d + u64::from(check.disagreed))
    });
    let poor: Vec<String> = checks
        .iter()
        .filter(|(_, check)| {
            check.checked >= MIN_FILL_CHECK_SITES && check.wrong_share_at_least() > POOR_FILL_SHARE
        })
        .map(|(sample, check)| format!("{sample} ({} of {})", check.disagreed, check.checked))
        .collect();
    let pooled = if checked > 0 {
        disagreed as f64 / checked as f64
    } else {
        0.0
    };
    if checked >= u64::from(MIN_FILL_CHECK_SITES) && pooled > POOR_FILL_SHARE {
        tracing::warn!(
            "the neighbour-evidence pass's fills look poor: where the samples' reads settle their \
             genotype, the neighbours alone predicted it wrong at {:.1}% of {} sites, above {:.0}%. \
             The pass is built for cohorts of inbred lines; in an outbred cohort, or one whose \
             samples share few haplotypes, its fills are less reliable. Consider \
             --refine-genotypes off{table_note}",
            100.0 * pooled,
            run_output::with_thousands(checked),
            100.0 * POOR_FILL_SHARE,
        );
    }
    if !poor.is_empty() {
        let shown = poor.len().min(10);
        tracing::warn!(
            "{} sample(s) whose fills look poor, wrong at more than {:.0}% of the sites checked \
             beyond what chance explains — wrong of checked: {}{}{table_note}",
            poor.len(),
            100.0 * POOR_FILL_SHARE,
            poor[..shown].join(", "),
            if poor.len() > shown { ", …" } else { "" },
        );
    }
}

/// `<vcf name without .vcf, .vcf.gz or .vcf.bgz>.fill_check.tsv` beside the VCF; `None` when the VCF
/// went to standard output.
fn fill_check_table_path(vcf_path: &Path) -> Option<PathBuf> {
    if vcf_path.as_os_str() == "-" {
        return None;
    }
    let name = vcf_path.file_name()?.to_string_lossy().into_owned();
    let lower = name.to_ascii_lowercase();
    let stem = [".vcf.gz", ".vcf.bgz", ".vcf"]
        .iter()
        .find(|suffix| lower.ends_with(*suffix))
        .map_or(name.as_str(), |suffix| &name[..name.len() - suffix.len()]);
    Some(vcf_path.with_file_name(format!("{stem}.fill_check.tsv")))
}

fn write_fill_check_table(
    path: &Path,
    checks: &[(String, crate::neighbour_evidence::refine::FillCheck)],
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "sample\tfilled\tsites_checked\tfills_wrong_there\twrong_share"
    )?;
    for (sample, check) in checks {
        let share = if check.checked >= MIN_FILL_CHECK_SITES {
            format!(
                "{:.4}",
                f64::from(check.disagreed) / f64::from(check.checked)
            )
        } else {
            ".".to_string()
        };
        writeln!(
            out,
            "{sample}\t{}\t{}\t{}\t{share}",
            check.filled, check.checked, check.disagreed
        )?;
    }
    out.flush()
}

/// **`--refine-genotypes`'s three values** (`neighbour_evidence_pass.md` §5, ruled §9 Q3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum RefineGenotypesArg {
    /// Run where the cohort is at least 10 samples with median inbreeding at least 0.5.
    Auto,
    /// Run wherever it can: two or more samples, diploid.
    On,
    /// Do not run. The default, by the owner's decision (2026-11-23).
    Off,
}

/// **Warn when reads phase records with the hidden-duplication filter off** (read phase spec
/// §2.5): reads from a collapsed duplication err together, every one carrying the other copy's
/// alleles, which the error model does not cover; the filter is what keeps such records out of
/// `PASS`. Phase from reads runs at ploidy two whether or not refinement does.
pub fn warn_of_read_phase_without_the_duplication_filter(filter_on: bool, ploidy: Ploidy) {
    if !filter_on && ploidy.get() == 2 {
        tracing::warn!(
            "the hidden-duplication filter is off (--hidden-dup-fdr 0), and phase is written from \
             reads: in a collapsed duplication every read carries the other copy's alleles, which \
             can link two records confidently and wrongly"
        );
    }
}

/// **`--phase-min-path-agreement`'s value**: a share in `(0, 1]` (phase output spec §9).
///
/// # Errors
///
/// A value that is not a number in `(0, 1]`.
pub fn parse_path_agreement(text: &str) -> Result<PathAgreement, String> {
    let share: f64 = text
        .parse()
        .map_err(|_| format!("`{text}` is not a number"))?;
    PathAgreement::try_new(share).ok_or_else(|| format!("`{text}` is not a share in (0, 1]"))
}

/// **`--refine-margin-density`'s parser**: records a kb to three decimals; 0, no thinning.
///
/// # Errors
///
/// As [`MarginThinning::parse`](crate::neighbour_evidence::windows::MarginThinning::parse).
pub fn parse_margin_density(text: &str) -> Result<MarginDensity, String> {
    crate::neighbour_evidence::windows::MarginThinning::parse(text).map(MarginDensity)
}

/// **`--refine-margin-density`'s value**: the far margins' thinning, or none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MarginDensity(pub Option<crate::neighbour_evidence::windows::MarginThinning>);

impl std::fmt::Display for MarginDensity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Some(thinning) => write!(f, "{thinning}"),
            None => write!(f, "0"),
        }
    }
}

/// **The refinement window options a calling run was given** (`--refine-window-bytes`,
/// `--refine-window-floor`, `--refine-core-floor`, `--refine-margin-density` and the hidden
/// `--refine-window-max` and `--refine-core-max`), as [`neighbour_evidence_decision`] reads them.
#[derive(Clone, Copy, Debug)]
pub struct RefineWindowArgs {
    pub window_bytes: u64,
    pub window_floor: u64,
    /// The window floor where not given.
    pub core_floor: Option<u64>,
    /// 0 walks every record.
    pub margin_density: MarginDensity,
    pub window_max: Option<u64>,
    pub core_max: Option<u64>,
}

/// **Whether this run uses the neighbour-evidence pass, and why**, from the option, each sample's
/// inbreeding coefficient in the run's sample order, and the ploidy — taken once, before the
/// header is written, and stated in it either way. `path_agreement` is
/// `--phase-min-path-agreement`, `windows` the window options.
#[must_use]
pub fn neighbour_evidence_decision(
    asked: RefineGenotypesArg,
    inbreeding_by_sample: &[InbreedingF],
    ploidy: Ploidy,
    path_agreement: PathAgreement,
    windows: RefineWindowArgs,
) -> NeighbourEvidenceDecision {
    let switch = match asked {
        RefineGenotypesArg::Auto => NeighbourEvidenceSwitch::Auto,
        RefineGenotypesArg::On => NeighbourEvidenceSwitch::On,
        RefineGenotypesArg::Off => NeighbourEvidenceSwitch::Off,
    };
    let inbreeding: Vec<f64> = inbreeding_by_sample.iter().map(|f| f.get()).collect();
    let decision = decide(switch, &inbreeding, ploidy)
        .with_path_agreement(path_agreement)
        .with_window_bytes(windows.window_bytes)
        .with_window_floor(windows.window_floor)
        .with_margin_thinning(windows.margin_density.0);
    // Not given, the core floor is the window floor, which `with_window_floor` has set.
    let decision = match windows.core_floor {
        Some(bases) => decision.with_core_floor(bases),
        None => decision,
    };
    let decision = match windows.window_max {
        Some(bases) => decision.with_window_max(bases),
        None => decision,
    };
    let decision = match windows.core_max {
        Some(bases) => decision.with_core_max(bases),
        None => decision,
    };
    tracing::info!(
        "neighbour-evidence pass: {} ({})",
        if decision.ran { "runs" } else { "off" },
        decision.why
    );
    decision
}

/// **Whether refinement's expected peak fits the memory this run may use**, judged before calling
/// starts (`neighbour_evidence_pass.md` §6, memory): from the window budget, the samples, the
/// threads and the longest contig, against the smaller of the machine's memory and the control
/// group's limit ([`crate::memory_limit`]). **It never stops the run**: above
/// [`WARN_SHARE`](crate::neighbour_evidence::peak::WARN_SHARE) of the limit it warns, saying the
/// estimate is an upper bound and what would fit; otherwise it logs the estimate. It hands the
/// limit and `scratch_dir`, the run's scratch directory where it keeps one, to the decision, so
/// that every window is estimated again once it is planned and its warning can name the command
/// that resumes the run. Does nothing where refinement does not run.
#[must_use]
pub fn check_refinement_memory(
    metadata: VcfHeaderMetadata,
    scratch_dir: Option<&Path>,
) -> VcfHeaderMetadata {
    let Some(decision) = metadata
        .neighbour_evidence()
        .filter(|decision| decision.ran)
        .cloned()
    else {
        return metadata;
    };
    let longest_contig = metadata
        .contigs()
        .iter()
        .map(|contig| contig.length)
        .max()
        .unwrap_or(0);
    let (decision, _) = refinement_memory_checked(
        decision,
        crate::memory_limit::of_this_process(),
        scratch_dir,
        rayon::current_num_threads(),
        longest_contig,
    );
    metadata.the_neighbour_evidence_decision(decision)
}

/// [`check_refinement_memory`]'s judgement, given the limit and the threads: the decision carrying
/// the limit, and the warning it logged, if it warned.
pub(crate) fn refinement_memory_checked(
    decision: NeighbourEvidenceDecision,
    limit: Option<crate::memory_limit::MemoryLimit>,
    scratch_dir: Option<&Path>,
    threads: usize,
    longest_contig: u64,
) -> (NeighbourEvidenceDecision, Option<String>) {
    use crate::neighbour_evidence::peak::{self, MemoryWatch, PeakCheck, PeakEstimate, PeakMoment};
    let samples = decision.inbreeding_by_sample.len();
    let estimate_at = |budget: u64| {
        PeakEstimate::before_any_window(
            &crate::neighbour_evidence::windows::WindowSizes {
                budget_bytes: budget,
                ..decision.sizes
            },
            samples,
            threads,
            longest_contig,
        )
    };
    let estimate = estimate_at(decision.sizes.budget_bytes);
    let found = peak::check(estimate, limit, PeakMoment::BeforeCalling, estimate_at);
    let warning = match found {
        PeakCheck::NearLimit(near) => {
            let warning = near.warning(scratch_dir);
            tracing::warn!("{warning}");
            Some(warning)
        }
        PeakCheck::Fits { .. } | PeakCheck::NoLimit => {
            tracing::info!(
                "refinement memory, before calling: {}",
                peak::fits_in_words(&estimate, &found, limit)
            );
            None
        }
    };
    let decision = NeighbourEvidenceDecision {
        memory_watch: MemoryWatch {
            limit,
            scratch_dir: scratch_dir.map(Path::to_path_buf),
        },
        ..decision
    };
    (decision, warning)
}

/// **The ploidy this run was asked for, judged against what the caller can score.**
///
/// Two refusals, and the second is the one that used to be a panic. `Ploidy::try_new` turns
/// down zero and nothing else, so seventeen copies is a value the type admits and the read
/// likelihood's copy-share table does not — it asserts, and a person who typed `--ploidy 20`
/// got a Rust backtrace naming a source file after the whole cohort had been opened. A
/// polyploid crop is an ordinary thing to call, so the number is judged here, before anything
/// is read, and the ceiling is named.
pub fn ploidy_asked_for(flag: Option<u8>) -> Result<Ploidy, CallingRunError> {
    let asked = flag.unwrap_or(DEFAULT_PLOIDY);
    let ploidy =
        Ploidy::try_new(asked).map_err(|source| CallingRunError::Ploidy { asked, source })?;
    if usize::from(asked) > MAX_PLOIDY_COPIES {
        return Err(CallingRunError::PloidyPastWhatIsScored {
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
pub fn refuse_an_output_that_cannot_be_written(output: &Path) -> Result<(), CallingRunError> {
    if output.is_dir() {
        return Err(CallingRunError::OutputIsADirectory {
            path: output.to_path_buf(),
        });
    }
    let directory = output.parent().unwrap_or_else(|| Path::new("."));
    // An empty parent is what a bare file name gives, and it means the working directory.
    if !directory.as_os_str().is_empty() && !directory.is_dir() {
        return Err(CallingRunError::OutputDirectoryIsMissing {
            path: output.to_path_buf(),
            directory: directory.to_path_buf(),
        });
    }
    Ok(())
}

/// Refuse a `--tmp-dir` that is not a directory, before anything is read.
///
/// # Errors
///
/// [`CallingRunError::TmpDirIsNotADirectory`].
pub fn refuse_a_tmp_dir_that_is_not_a_directory(
    tmp_dir: Option<&Path>,
) -> Result<(), CallingRunError> {
    match tmp_dir {
        Some(directory) if !directory.is_dir() => Err(CallingRunError::TmpDirIsNotADirectory {
            path: directory.to_path_buf(),
        }),
        _ => Ok(()),
    }
}

/// **The run's own scratch directory**, `pvc_tmp_<random>` inside `--tmp-dir`, or inside the
/// VCF's directory when no `--tmp-dir` was given (the working directory for a VCF on standard
/// output) (`doc/devel/ng/spec/command_line.md` §7.5).
///
/// **A directory of its own**, so several runs can share one `--tmp-dir` without colliding, and
/// removing it never touches anything the user put there.
///
/// What happens to it when the run ends is [`RunScratch`]'s, and `on_failure` says which way.
///
/// # Errors
///
/// [`CallingRunError::ScratchDirectory`] when it cannot be made.
pub fn scratch_directory_for(
    output: &Path,
    tmp_dir: Option<&Path>,
    on_failure: ScratchOnFailure,
) -> Result<RunScratch, CallingRunError> {
    let directory = match tmp_dir {
        Some(directory) => directory.to_path_buf(),
        None => match output.parent() {
            Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
            _ => PathBuf::from("."),
        },
    };
    let made = tempfile::Builder::new()
        .prefix("pvc_tmp_")
        .tempdir_in(&directory)
        .map_err(|source| CallingRunError::ScratchDirectory {
            directory: directory.clone(),
            source,
        })?;
    let made = made.keep();
    tracing::debug!("scratch directory: {}", made.display());
    Ok(RunScratch {
        directory: Some(made),
        on_failure,
        the_run_succeeded: false,
    })
}

/// **The scratch directory an earlier run kept, taken over by a run that resumes it**
/// (`pvc resume-calling`): removed when this run succeeds, kept when it fails and `on_failure`
/// asks for that, as the run that made it would have.
///
/// # Errors
///
/// [`CallingRunError::ScratchDirectoryNotFound`] when `directory` is not a directory.
pub fn kept_scratch_directory(
    directory: &Path,
    on_failure: ScratchOnFailure,
) -> Result<RunScratch, CallingRunError> {
    if !directory.is_dir() {
        return Err(CallingRunError::ScratchDirectoryNotFound {
            directory: directory.to_path_buf(),
        });
    }
    tracing::debug!("scratch directory, resumed: {}", directory.display());
    Ok(RunScratch {
        directory: Some(directory.to_path_buf()),
        on_failure,
        the_run_succeeded: false,
    })
}

/// What becomes of a run's scratch directory when the run ends in an error —
/// `--keep-tmp-on-error`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ScratchOnFailure {
    /// Removed, as on success. The default: a failed run with the hidden-duplication filter on
    /// can leave a spill several terabytes long on a cohort of thousands.
    Removed,
    /// Left on disk, with everything in it, and its path logged — for looking at what the run
    /// had made before it stopped.
    Kept,
    /// **Left on disk because a refinement window did not fit its budget**, whatever
    /// `--keep-tmp-on-error` said: over `--refine-window-bytes` even at the smallest window. That
    /// stop comes after calling has ended, so the directory holds everything `pvc resume-calling`
    /// needs, and removing it would throw away hours of calling to a setting the user can change.
    /// Set by [`RunScratch::keep_if_refinement_stopped`].
    KeptForResuming,
}

impl ScratchOnFailure {
    /// The choice `--keep-tmp-on-error` makes.
    #[must_use]
    pub fn from_keep_flag(keep_tmp_on_error: bool) -> Self {
        if keep_tmp_on_error {
            Self::Kept
        } else {
            Self::Removed
        }
    }
}

/// **A run's scratch directory and the one decision about it: removed or kept when the run
/// ends.**
///
/// The directory is removed when this is dropped, unless the run never said it succeeded and
/// [`ScratchOnFailure::Kept`] was asked for. A run says so with [`Self::the_run_succeeded`] as
/// its last step, so **every way out that is not that step counts as a failure**: an error
/// returned through `?` from anywhere in the run drops this without it.
///
/// **It is the one owner of what is inside.** The hidden-duplication filter's spill lives here
/// and does not remove itself (`SpillFile::in_directory`), so keeping the directory keeps the
/// spill, and removing it removes the spill. A panic in a release build aborts, so no `Drop`
/// runs and the directory stays whatever was asked; in a debug build it unwinds through here and
/// counts as a failure.
#[derive(Debug)]
pub struct RunScratch {
    /// `None` only once `Drop` has taken it.
    directory: Option<PathBuf>,
    on_failure: ScratchOnFailure,
    the_run_succeeded: bool,
}

impl RunScratch {
    /// Where the directory is.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.directory
            .as_ref()
            // PANIC-FREE: `directory` is `Some` from construction until `Drop` takes it.
            .expect("the scratch directory is there until it is dropped")
    }

    /// Record that the run finished and wrote its output, so the directory is removed whatever
    /// `--keep-tmp-on-error` said.
    pub fn the_run_succeeded(&mut self) {
        self.the_run_succeeded = true;
    }

    /// **Keep the directory if `error` is a refinement window over its budget**
    /// ([`ScratchOnFailure::KeptForResuming`]); any other error leaves the choice as it was.
    pub fn keep_if_refinement_stopped(&mut self, error: &(dyn std::error::Error + 'static)) {
        if refinement_stopped_over_budget(error) {
            self.on_failure = ScratchOnFailure::KeptForResuming;
        }
    }
}

/// **Whether `error`, or any error under it, is refinement stopping on a window over its budget**
/// even at its smallest. Refinement never stops on its memory estimate; it only warns.
#[must_use]
pub fn refinement_stopped_over_budget(error: &(dyn std::error::Error + 'static)) -> bool {
    let mut next = Some(error);
    while let Some(this) = next {
        if matches!(
            this.downcast_ref::<crate::vcf::writer::VcfWriteError>(),
            Some(crate::vcf::writer::VcfWriteError::RefinementOverBudget(_))
        ) {
            return true;
        }
        next = this.source();
    }
    false
}

impl Drop for RunScratch {
    fn drop(&mut self) {
        let Some(directory) = self.directory.take() else {
            return;
        };
        if !self.the_run_succeeded && self.on_failure != ScratchOnFailure::Removed {
            // A run that stopped after calling left what resuming needs; one that stopped before
            // left a spill nothing can finish.
            let resumable = crate::run::paralog_filter::resume::is_resumable(&directory);
            let why = if self.on_failure == ScratchOnFailure::KeptForResuming {
                "because a refinement window did not fit --refine-window-bytes and calling had finished"
            } else {
                "as --keep-tmp-on-error asked"
            };
            tracing::warn!(
                "the run did not finish; its scratch directory is kept, {why}: {} — remove it \
                 once it is no longer needed{}",
                directory.display(),
                if resumable {
                    format!(
                        ", or finish the run from it with `pvc resume-calling --scratch-dir {}`",
                        directory.display()
                    )
                } else {
                    String::new()
                }
            );
            return;
        }
        let path = directory;
        if let Err(source) = std::fs::remove_dir_all(&path) {
            // A directory left on disk, possibly terabytes of spill, is worth a line; a `Drop`
            // must not panic, and a second message while a panic unwinds would bury the first.
            if !std::thread::panicking() {
                tracing::warn!(
                    "the run's scratch directory `{}` could not be removed and is still on \
                     disk: {source}",
                    path.display()
                );
            }
        }
    }
}

/// Where the calls went, as a person reads it: the file, or `standard output`.
#[must_use]
pub fn where_the_calls_went(output: &Path) -> String {
    if crate::vcf::writer::is_standard_output(output) {
        "standard output".to_owned()
    } else {
        output.display().to_string()
    }
}

/// **The numbers this run scores with, and the two things a file recording them needs beside
/// each.**
///
/// `RunParameters` keeps what *calling* reads — one bare multiplier a read group, one bare
/// coefficient a sample — and a file has to say how much data stood behind each and under what
/// warrant. Neither is recoverable from the parameters afterwards, so both travel from wherever
/// the numbers came from (`ParametersFile::of_run`'s own note says the same of the inbreeding
/// warrants).
pub struct TheRunsNumbers {
    /// What every locus is scored against.
    pub parameters: RunParameters,
    /// How many reads stood behind each read group's base-quality multiplier, in the run's dense
    /// read-group order.
    pub reads_behind_each_calibration: ReadsBehindEachCalibration,
    /// Each sample's inbreeding coefficient with its warrant, in the run's sample order.
    pub inbreeding_by_sample: Vec<Estimate<InbreedingF>>,
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
    pub census: CensusIdentity,
    /// **How each start of the SNP/indel fit behind these numbers ended, as the file that carried
    /// them recorded it** — `None` on the defaults path, or a file that recorded none. Carried
    /// rather than dropped for the census's reason: a run writing its parameters out again writes
    /// back what it read.
    pub snp_indel_fit_starts: Option<Vec<SnpIndelFitStart>>,
    /// **What the data determined about the fitted density and the duplicated class, as the file
    /// that carried these numbers recorded it** — format 4's two informational tables
    /// (`parameter_prepass_joint_fit.md` §2.1.3); `None` on the defaults path, or where the file
    /// has none. No run reads them; carried so that a calling run does not silently drop them
    /// from the file it writes back.
    pub fitted_density: Option<FittedDensity>,
    pub duplicated_positions: Option<FittedDuplicatedPositions>,
    /// **Which tracts the fit behind the file left out of its stutter fit** — format 7's
    /// informational table; carried, like the two above, so a calling run does not drop it.
    pub repeat_tracts_left_out: Option<RepeatTractsLeftOut>,
    /// **The SNP/indel fit's own numbers as the fit behind the file returned them** — format 10's
    /// informational table; carried, like the three above, so a calling run does not drop it.
    pub snp_indel_fit: Option<SnpIndelFitAsReturned>,
}

impl TheRunsNumbers {
    /// **The parameters file this run writes** — what it scored with, and everything the file it
    /// read said about where those numbers came from: the counts behind each multiplier, each
    /// coefficient's warrant, the census, how each start of the fit behind them ended, and what
    /// the data determined about the fit's density and duplicated class.
    ///
    /// **The one door both calling commands write through**, so that neither can drop something
    /// it read: `ParametersFile::of_run` writes no fit starts and no informational tables, and a
    /// writer that forgot [`ParametersFile::with_snp_indel_fit_starts`] or
    /// [`ParametersFile::with_what_the_data_determine`] would compile and lose them.
    #[must_use]
    pub fn parameters_file(
        &self,
        read_groups: &ReadGroups,
        reference: &ReferenceDigest,
        repeat_routing: Option<&StrRepeatCriteria>,
    ) -> ParametersFile {
        ParametersFile::of_run_under(
            &self.parameters,
            read_groups,
            &self.reads_behind_each_calibration,
            &self.inbreeding_by_sample,
            reference,
            self.census.clone(),
            repeat_routing,
        )
        .with_snp_indel_fit_starts(self.snp_indel_fit_starts.clone())
        .with_what_the_data_determine(
            self.fitted_density.clone(),
            self.duplicated_positions.clone(),
        )
        .with_repeat_tracts_left_out(self.repeat_tracts_left_out.clone())
        .with_snp_indel_fit(self.snp_indel_fit.clone())
    }
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
pub fn run_parameters(
    supplied: Option<&Path>,
    ploidy_flag: Option<u8>,
    read_groups: &ReadGroups,
    with_checksums: &ReferenceInfo,
    ploidy: Ploidy,
    routing: Option<&StrRepeatCriteria>,
) -> Result<TheRunsNumbers, CallingRunError> {
    let Some(path) = supplied else {
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
            snp_indel_fit_starts: None,
            fitted_density: None,
            duplicated_positions: None,
            repeat_tracts_left_out: None,
            snp_indel_fit: None,
        });
    };
    let text =
        std::fs::read_to_string(path).map_err(|source| CallingRunError::ParametersUnreadable {
            path: path.to_path_buf(),
            source,
        })?;
    let file = ParametersFile::from_toml(&text).map_err(|source| CallingRunError::Parameters {
        path: path.to_path_buf(),
        source,
    })?;
    // **A file written by a run that routed differently is reported and never refused** — the
    // owner's ruling, `run_ssr_observations.md` §2.3 and `parameters_file.md` §3.9. Which
    // stretches count as repeats is the user's to choose independently of where the numbers
    // came from, so this changes nothing about the run: the numbers are as warranted as they
    // ever were, and only the ground they are applied to has moved. What the operator is owed is
    // knowing it happened, because a tract this run admits and the fit's selection did not is
    // scored from strata fitted over other loci.
    // A run that sorted no repeats (`--generic-path-only`) calls no repeat tract, so no routing
    // the file records can disagree with anything it does.
    if let Some(routing) = routing
        && let Some(axis) = file.routing_disagreement(routing)
    {
        tracing::warn!(
            "{} was written by a run that routed on a different {axis}; calling on",
            path.display()
        );
    }
    let digest = ReferenceDigest::of(with_checksums)
        .map_err(|source| CallingRunError::ReferenceNotDigested { source })?;
    let bound = file
        .to_run_parameters_for(&digest, read_groups, None)
        .map_err(|source| CallingRunError::Parameters {
            path: path.to_path_buf(),
            source,
        })?;
    let from_file = bound.from_file;
    // **A flag that was typed may only agree with the file.** Spec §3.2 puts the ploidy in the
    // file so that a supplied one "cannot be paired with a run at a different ploidy without
    // saying so", and calling at the file's number while an operator typed another is exactly
    // that. Only a ploidy that was *typed* is compared: `--ploidy` is an `Option` for this
    // reason, so a tetraploid file is not refused for the flag's default being two.
    if let Some(asked) = ploidy_flag
        && asked != from_file.parameters.ploidy().get()
    {
        return Err(CallingRunError::PloidyIsNotTheParametersFiles {
            path: path.to_path_buf(),
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
        snp_indel_fit_starts: file.fitted_from.snp_indel_fit_starts.clone(),
        fitted_density: file.ordinary_site_prior.fitted_density.clone(),
        duplicated_positions: file.duplicated_positions.clone(),
        repeat_tracts_left_out: file.repeat_tracts_left_out.clone(),
        snp_indel_fit: file.snp_indel_fit.clone(),
    })
}

/// **The two bounds a run assembles a locus under**, judged where they are typed.
///
/// The round width is the one number here that is *chosen* rather than checked: left unset it
/// comes from the cohort's size ([`round_width_for`]), which changes no answer and does change
/// where a round's edge falls. Both calling commands take it from here, so the same cohort holds
/// the same ground a round whichever way its observations were read; how that ground is cut into
/// regions depends on which (`read_from`, [`round_shape_for`]). Neither changes the output, which
/// spec §12.3 compares between the two byte for byte.
///
/// # Errors
///
/// [`CallingRunError::MaxCohortLocusSpanIsZero`] and
/// [`CallingRunError::CohortLocusBuilderRegionsLenIsZero`] — a bound of zero refuses every
/// locus there is, and a round of zero bases never advances.
pub fn merge_parameters_for(
    max_cohort_locus_span: u32,
    cohort_locus_builder_regions_len: Option<u32>,
    samples: usize,
    read_from: ObservationsReadFrom,
) -> Result<MergeParameters, CallingRunError> {
    Ok(MergeParameters {
        max_cohort_locus_span: MaxCohortLocusSpan(NonZeroU32::new(max_cohort_locus_span).ok_or(
            CallingRunError::MaxCohortLocusSpanIsZero {
                asked: max_cohort_locus_span,
            },
        )?),
        // **One decision, two numbers.** The ground a round holds is the width the flag names
        // or the one the cohort's size gives; `round_shape_for` then divides that ground
        // between the regions worked at once, so asking for a round costs the memory one
        // region already cost rather than that memory times the core count.
        ..{
            let ground = match cohort_locus_builder_regions_len {
                Some(asked) => CohortLocusBuilderRegionsLen(
                    NonZeroU32::new(asked)
                        .ok_or(CallingRunError::CohortLocusBuilderRegionsLenIsZero { asked })?,
                ),
                None => round_width_for(samples),
            };
            let (width, in_flight) = round_shape_for(ground, samples, read_from);
            MergeParameters {
                cohort_locus_builder_regions_len: width,
                cohort_locus_builder_regions_in_flight: in_flight,
                ..MergeParameters::DEFAULT
            }
        }
    })
}

/// **How many alleles a locus may be called over**, judged where it is typed.
///
/// # Errors
///
/// [`CallingRunError::MaxCandidateAllelesTooSmall`] for a cap below the reference plus one
/// alternative, which is a refusal of every locus under another name.
pub fn candidate_selection_for(
    max_candidate_alleles: u16,
) -> Result<CandidateSelectionConfig, CallingRunError> {
    Ok(CandidateSelectionConfig {
        max_candidate_alleles: MaxCandidateAlleles::new(max_candidate_alleles).ok_or(
            CallingRunError::MaxCandidateAllelesTooSmall {
                asked: max_candidate_alleles,
                smallest: MaxCandidateAlleles::SMALLEST,
            },
        )?,
        ..CandidateSelectionConfig::DEFAULT
    })
}

#[cfg(test)]
mod tests;

/// **The parameters a run recorded in its log, copied out to a file beside its VCF** — what a
/// person does to reuse them, done for a test. `log` is the run's `--log`, or `None` for the
/// default beside the VCF.
///
/// # Panics
///
/// When the log is not there or holds no parameters block: a test calls this after a run that
/// must have written both.
#[cfg(test)]
pub(crate) fn parameters_the_run_logged(output: &Path, log: Option<&Path>) -> PathBuf {
    let log = log.map_or_else(
        || crate::cli::run_output::beside(output, "pvc_calls.log"),
        Path::to_path_buf,
    );
    let text = std::fs::read_to_string(&log).expect("the run's log reads");
    let parameters =
        crate::cli::run_output::parameters_in(&text).expect("the log holds the run's parameters");
    let mut name = output.as_os_str().to_os_string();
    name.push(".logged-parameters.toml");
    let at = PathBuf::from(name);
    std::fs::write(&at, parameters).expect("the parameters are copied out");
    at
}

/// **What the VCF's header states about the run that wrote it**, whichever mode wrote it: among
/// the rest, the run's log by file name, and where its numbers came from
/// (`doc/devel/ng/spec/command_line.md` §7.6).
///
/// **The log by name and not by path**, because the two usually travel as a directory and a path
/// would be wrong the moment somebody moved the pair.
///
/// # Errors
///
/// [`CallingRunError::Header`] where the metadata is not something a header can state — a
/// sample named twice, or a contig list the writer will not take.
pub fn header_for(
    log: &Path,
    parameters: Option<&Path>,
    reference: &Path,
    contigs: &ContigList,
    with_checksums: &ReferenceInfo,
    samples: Vec<String>,
) -> Result<VcfHeaderMetadata, CallingRunError> {
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
        samples,
        current_command_line(),
        reference.display().to_string(),
        log.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned(),
        where_the_numbers_came_from(parameters),
    )
    .map_err(|source| CallingRunError::Header { source })
}

/// **Where a run's numbers came from, as the header says it**: `defaults`, or the parameters
/// file's name and the MD5 of its bytes — the name alone would say nothing once the file was
/// edited. A file that cannot be read again here is named without a digest; the run read it once
/// already, and a header is not the place to fail.
#[must_use]
pub fn where_the_numbers_came_from(parameters: Option<&Path>) -> String {
    use md5::{Digest, Md5};

    let Some(path) = parameters else {
        return "defaults".to_owned();
    };
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    match std::fs::read(path) {
        Ok(bytes) => {
            let hex = crate::reference_info::format_md5_hex(Md5::digest(&bytes).into());
            format!("{name} (md5 {hex})")
        }
        Err(_) => name.into_owned(),
    }
}

/// **What a calling pass produced, in either mode** — the part of the run's tail that differs
/// between the two commands.
pub enum CalledRun {
    /// Direct mode's: the calling, and what each sample's walk saw.
    Walked(crate::run::WrittenCohort),
    /// psp mode's: the calling, and what this run drew out of each stored file.
    Stored(
        crate::run::callers::CohortCallingTallies,
        crate::run::psp_caller::StoredCohortTallies,
    ),
}

impl CalledRun {
    fn calling_mut(&mut self) -> &mut crate::run::callers::CohortCallingTallies {
        match self {
            Self::Walked(written) => &mut written.calling,
            Self::Stored(calling, _) => calling,
        }
    }

    fn window_coverage_histograms_mut(
        &mut self,
    ) -> &mut Vec<crate::window_coverage::SampleHistogram> {
        match self {
            Self::Walked(written) => &mut written.window_coverage_histograms,
            Self::Stored(_, stored) => &mut stored.window_coverage_histograms,
        }
    }

    fn report<'a>(
        &'a self,
        contigs: &'a ContigList,
        read_groups: &'a ReadGroups,
        parameters_file: &'a ParametersFile,
        analysed: &crate::region_typing::GenomeRegions,
        bounds: crate::run::report::BoundsTheRunCalledUnder,
    ) -> RunReport<'a> {
        match self {
            Self::Walked(written) => RunReport::of(
                written,
                contigs,
                read_groups,
                parameters_file,
                analysed,
                bounds,
            ),
            Self::Stored(calling, stored) => RunReport::of_a_stored_cohort(
                calling,
                stored,
                contigs,
                read_groups,
                parameters_file,
                analysed,
                bounds,
            ),
        }
    }
}

/// **What the run's tail reads beside the calls**: everything both commands had built before
/// calling and hand on unchanged.
pub struct RunTail<'a> {
    /// The VCF's path, or `-`.
    pub output: &'a Path,
    /// The `--explain-loci` regions and the table they are written to, if asked for.
    pub explain: Option<&'a (crate::run::explain::ExplainRegions, PathBuf)>,
    /// The reference's contigs.
    pub contigs: &'a ContigList,
    /// The run's samples, in its sample order.
    pub sample_names: &'a [String],
    /// Each sample's inbreeding coefficient, in the run's sample order.
    pub inbreeding_by_sample: &'a [InbreedingF],
    /// What the hidden-duplication filter was asked for; `None` when it is off.
    pub asked_for: Option<crate::run::paralog_filter::WhatTheOperatorAskedFor>,
    /// The VCF's header.
    pub metadata: VcfHeaderMetadata,
    /// The run's ploidy.
    pub ploidy: Ploidy,
    /// The run-wide read-group table.
    pub read_groups: &'a ReadGroups,
    /// The parameters file the run logged.
    pub parameters_file: &'a ParametersFile,
    /// The ground the run analysed, as the report states it.
    pub analysed: &'a crate::region_typing::GenomeRegions,
    /// The bounds the run called under, as the report states them.
    pub bounds: crate::run::report::BoundsTheRunCalledUnder,
    /// The satellite bases `--call-satellites` handed to the generic path.
    pub satellite_bases_on_the_generic_path: u64,
    /// **For `call-from-psps --regions` only**: the BED the run called and what the psps cover,
    /// as one phrase for the report. `None` for every other run.
    pub regions_asked_for: Option<String>,
    /// When the run started, for the closing line.
    pub started: std::time::Instant,
    /// The run's log, named in the closing line.
    pub log: &'a Path,
}

/// What can stop the run's tail. Each command turns it into its own error, with the same text.
#[derive(Debug)]
pub enum RunTailError {
    /// The calls could not be written, or the records could not be parked for the filter.
    CallsNotWritten {
        /// Where they were going.
        path: PathBuf,
        /// What pass one said.
        source: crate::run::paralog_filter::PassOneError,
    },
    /// The hidden-duplication filter could not finish the run.
    ParalogFilter {
        /// What the filter said.
        source: crate::run::paralog_filter::ParalogFilterError,
    },
    /// The `--explain-loci` explanation could not be written.
    ExplanationNotWritten {
        /// The TSV.
        path: PathBuf,
        /// What the write hit.
        source: std::io::Error,
    },
}

/// **The tail both calling commands share, from the moment calling returns**: finish the sink,
/// run the hidden-duplication filter's passes two and three on what it parked, write the
/// `--explain-loci` table, log the run's report, and say how many records were written.
///
/// Only the last step, [`RunScratch::the_run_succeeded`], counts as success: a run that fails
/// here leaves `scratch` to be removed, or kept, as its flags say.
///
/// # Errors
///
/// [`RunTailError`], which each command turns into its own error.
pub fn finish_the_run<E: From<RunTailError>>(
    sink: crate::run::paralog_filter::CalledRecordSink,
    mut called: CalledRun,
    scratch: &mut Option<RunScratch>,
    tail: RunTail<'_>,
) -> Result<(), E> {
    use crate::run::paralog_filter::{self, CalledRecordSink, PassOneError};
    let RunTail {
        output,
        explain,
        contigs,
        sample_names,
        inbreeding_by_sample,
        asked_for,
        metadata,
        ploidy,
        read_groups,
        parameters_file,
        analysed,
        bounds,
        satellite_bases_on_the_generic_path,
        regions_asked_for,
        started,
        log,
    } = tail;
    // **The spill has to stay alive until pass three has read it**, so the sink is taken apart
    // rather than consumed: `finish` would drop the file the next two passes need.
    let (parked, finished_straight) = match sink {
        CalledRecordSink::StraightToTheVcf(writer) => {
            let finished = writer
                .finish()
                .map_err(|source| RunTailError::CallsNotWritten {
                    path: output.to_path_buf(),
                    source: PassOneError::Vcf(source),
                })?;
            (None, Some(finished))
        }
        CalledRecordSink::ParkedOnTheSpill(mut spilling) => {
            spilling
                .finish_parking()
                .map_err(|source| RunTailError::CallsNotWritten {
                    path: output.to_path_buf(),
                    source: PassOneError::Spill(source),
                })?;
            (Some(spilling), None)
        }
    };

    // **Passes two and three, through the one function both subcommands call**, so direct mode
    // and psp mode cannot drift into filtering differently (spec §1.1 goal 2).
    // The filter is asked for its verdicts on the records the explained loci wrote.
    let mut paralog_explanations = explain.map(|_| {
        crate::run::explain::ParalogExplanations::wanting(&called.calling_mut().explanations)
    });
    let filtered = match &parked {
        Some(spilling) => Some(
            paralog_filter::fit_score_and_write_the_calls_explaining(
                spilling.spill(),
                // **Taken, not cloned.** `ParalogScoringContext::new` consumes the histograms so
                // each sample's bins are freed as its model is fitted; cloning here would defeat
                // that and hold a second copy of every sample's — about 80 kB each. Nothing else
                // reads the field afterwards.
                std::mem::take(called.window_coverage_histograms_mut()),
                inbreeding_by_sample,
                // PANIC-FREE: the sink is the spill only where `asked_for` is `Some`, and this
                // arm is reached only through that sink.
                asked_for.expect("the records are parked only when a filter was asked for"),
                output,
                metadata,
                ploidy,
                paralog_explanations.as_mut(),
            )
            .map_err(|source| {
                // A stop on memory comes after calling: keep what resuming needs.
                if let Some(scratch) = scratch.as_mut() {
                    scratch.keep_if_refinement_stopped(&source);
                }
                RunTailError::ParalogFilter { source }
            })?,
        ),
        None => None,
    };

    if let Some(verdicts) = paralog_explanations {
        verdicts.add_to(&mut called.calling_mut().explanations);
    }
    if let Some((_, tsv)) = explain {
        let mut explanations = std::mem::take(&mut called.calling_mut().explanations);
        crate::run::explain::write_explanations(tsv, &mut explanations, contigs, sample_names)
            .map_err(|source| RunTailError::ExplanationNotWritten {
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
    log_the_report(
        output,
        &called
            .report(contigs, read_groups, parameters_file, analysed, bounds)
            .with_satellite_bases_on_the_generic_path(satellite_bases_on_the_generic_path)
            .with_regions_asked_for(regions_asked_for),
        &filter_lines,
    );
    // **Counted from the file**: the neighbour-evidence pass can drop a record whose alternatives
    // it left uncalled, so neither calling's count nor the filter's is the file's.
    let finished = finished_straight
        .or_else(|| filtered.as_ref().map(|run| run.finished.clone()))
        .expect("the calls were written by one of the two paths");
    report_neighbour_evidence(&finished, output);
    report_read_phase(&finished);
    let records = finished.records_written;
    tracing::info!(
        "done: {} records written to {} in {}",
        run_output::with_thousands(records),
        where_the_calls_went(output),
        run_output::elapsed_text(started.elapsed())
    );
    tracing::info!("log: {}", log.display());
    if let Some(scratch) = scratch.as_mut() {
        scratch.the_run_succeeded();
    }
    Ok(())
}

/// **The run report, written to the run's log**: the terminal says only how the run ended, and
/// the counts behind that are the log's (`doc/devel/ng/spec/command_line.md` §7.2). Each line is
/// prefixed `report:` so it can be found among the log's events; `filter_lines` are the
/// hidden-dup filter's, empty on a run with the filter off.
pub fn log_the_report(calls: &Path, report: &RunReport<'_>, filter_lines: &[String]) {
    tracing::debug!("report: calls: {}", where_the_calls_went(calls));
    for line in report.lines().iter().chain(filter_lines) {
        // **A warning in the report is a warning on the terminal too**, where the rest of the
        // report is not: it is the line a person must not miss.
        match line.trim_start().strip_prefix("warning: ") {
            Some(warning) => tracing::warn!("{warning}"),
            None => tracing::debug!("report: {line}"),
        }
    }
}

/// `--low-base-quality-weight`'s parser: a finite number, zero or more.
///
/// # Errors
///
/// A value that is not a number, is negative, or is not finite.
pub fn parse_low_base_quality_weight(value: &str) -> Result<f64, String> {
    let weight: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if weight.is_finite() && weight >= 0.0 {
        Ok(weight)
    } else {
        Err(format!(
            "`{value}` is not a weight: it must be zero or more, and finite"
        ))
    }
}

/// `--min-site-quality`'s parser: a finite quality, zero or more.
///
/// # Errors
///
/// A value that is not a number, is negative, or is not finite — said in the words the command
/// line prints.
pub fn parse_min_site_quality(value: &str) -> Result<f32, String> {
    let quality: f32 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    crate::types::Phred::try_new(quality)
        .map(|_| quality)
        .map_err(|_| {
            format!("`{value}` is not a site quality: it must be zero or more, and finite")
        })
}
