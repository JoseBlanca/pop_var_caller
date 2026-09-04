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

use std::num::NonZeroU32;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::calling::allele_candidates::{CandidateSelectionConfig, MaxCandidateAlleles};
use crate::ng::calling::inference::{CallingLoopConfig, DiscoveryMode, RunnableCallingLoopConfig};
use crate::ng::calling::likelihood::MAX_PLOIDY_COPIES;
use crate::ng::calling::parameters_file::{
    CensusIdentity, DeclaredInbreeding, ParametersFile, ParametersFileError,
    ReadsBehindEachCalibration, beside_the_vcf,
};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::parameter_estimation::Estimate;
use crate::ng::parameter_estimation::joint::loci::{ReferenceDigest, SelectionError};
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::ng::run::cohort_merge::{
    CohortLocusBuilderRegionsLen, DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, MaxCohortLocusSpan,
};
use crate::ng::run::{MergeParameters, RunReport};
use crate::ng::types::{DomainError, InbreedingF, Ploidy};
use crate::ng::vcf::header::{HeaderContig, HeaderMetadataError, VcfHeaderMetadata};
use crate::pop_var_caller::common::current_command_line;

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

/// **The tract-accuracy program's measurement switch for allele discovery** (lever L7): run
/// the discovery pre-pass at every repeat tract, admitting the tract sequences one sample
/// showed too often for slippage to explain — 2 reads **and** 15% of that sample's spanning
/// reads, the shipped bar.
///
/// Absent is the shipped default — no discovery, selection untouched. `1` asks for
/// [`DiscoveryMode::BeforeTheLoop`], the pre-pass inside candidate selection that milestone
/// E1's finding put there (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md`
/// §6.5); anything else set is refused before a read is decoded, for the reason the other
/// switches give: a measurement switch that fell back silently would report the plain
/// selection as the discovery arm.
///
/// **An environment variable by design, not an oversight**: there is no parameters-file key
/// for this yet, deliberately — the arm is enabled per run while the program measures it, and
/// a *keep* verdict owes this switch proper parameters-file plumbing before the experiment
/// spelling is retired.
const NG_TRACT_DISCOVERY: &str = "NG_TRACT_DISCOVERY";

/// The calling-loop configuration this run asks for: the shipped defaults, with the slippage
/// re-fit's round count read once from [`NG_SLIPPAGE_REFIT_ROUNDS`], its pull-backs
/// zeroed where [`NG_SLIPPAGE_REFIT_FREE`] asks for the free setting, and the discovery
/// pre-pass switched on where [`NG_TRACT_DISCOVERY`] asks for it.
pub fn calling_loop_config_for_this_run() -> Result<RunnableCallingLoopConfig, CallingRunError> {
    let mut config = CallingLoopConfig::DEFAULT;
    match std::env::var(NG_SLIPPAGE_REFIT_ROUNDS) {
        Ok(rounds) => {
            config.slippage_refit.max_rounds = rounds.trim().parse().map_err(|error| {
                CallingRunError::CallingLoopSettings(format!(
                    "{NG_SLIPPAGE_REFIT_ROUNDS} must be a whole number of re-fit rounds \
                     (0 is the frozen default), not {rounds:?}: {error}"
                ))
            })?;
        }
        Err(std::env::VarError::NotPresent) => {}
        Err(std::env::VarError::NotUnicode(value)) => {
            return Err(CallingRunError::CallingLoopSettings(format!(
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
            return Err(CallingRunError::CallingLoopSettings(format!(
                "{NG_SLIPPAGE_REFIT_FREE} accepts only 1 (zero both pull-backs), \
                 not {value:?}"
            )));
        }
        Err(std::env::VarError::NotPresent) => {}
        Err(std::env::VarError::NotUnicode(value)) => {
            return Err(CallingRunError::CallingLoopSettings(format!(
                "{NG_SLIPPAGE_REFIT_FREE} is set but is not text: {value:?}"
            )));
        }
    }
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
pub fn refuse_an_output_whose_parameters_file_is_this_run_s_input(
    output: &Path,
    supplied: Option<&Path>,
) -> Result<(), CallingRunError> {
    let Some(supplied) = supplied else {
        return Ok(());
    };
    if !supplied.exists() {
        return Ok(());
    }
    let would_write = beside_the_vcf(output);
    if resolved(&would_write) == resolved(supplied) {
        return Err(CallingRunError::ParametersWouldBeOverwritten {
            path: supplied.to_path_buf(),
            output: output.to_path_buf(),
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
    routing: &StrRepeatCriteria,
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
    if let Some(axis) = file.routing_disagreement(routing) {
        eprintln!(
            "note: {} was written by a run that routed on a different {axis}; calling on",
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
    })
}

/// **The two bounds a run assembles a locus under**, judged where they are typed.
///
/// The round width is the one number here that is *chosen* rather than checked: left unset it
/// comes from the cohort's size ([`round_width_for`]), which changes no answer and does change
/// where a round's edge falls. Both calling commands take it from here so that the same cohort
/// gets the same width whichever way its observations were read — spec §12.3 compares their
/// output byte for byte.
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
) -> Result<MergeParameters, CallingRunError> {
    Ok(MergeParameters {
        max_cohort_locus_span: MaxCohortLocusSpan(NonZeroU32::new(max_cohort_locus_span).ok_or(
            CallingRunError::MaxCohortLocusSpanIsZero {
                asked: max_cohort_locus_span,
            },
        )?),
        cohort_locus_builder_regions_len: match cohort_locus_builder_regions_len {
            Some(asked) => CohortLocusBuilderRegionsLen(
                NonZeroU32::new(asked)
                    .ok_or(CallingRunError::CohortLocusBuilderRegionsLenIsZero { asked })?,
            ),
            None => round_width_for(samples),
        },
        ..MergeParameters::DEFAULT
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

/// **What the VCF's header states about the run that wrote it**, whichever mode wrote it.
///
/// **One copy, because one of the decisions in it survived a mutation until a command-level
/// test existed**: `##parametersFile` names the file this run writes beside its VCF, by name
/// and not by path — the two are siblings by construction ([`beside_the_vcf`]), so a path would
/// say the same thing at greater length and would be wrong the moment somebody moved the pair.
/// Two copies of that would be two places to get it wrong and one place to fix it.
///
/// # Errors
///
/// [`CallingRunError::Header`] where the metadata is not something a header can state — a
/// sample named twice, or a contig list the writer will not take.
pub fn header_for(
    output: &Path,
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
        beside_the_vcf(output)
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned(),
    )
    .map_err(|source| CallingRunError::Header { source })
}

/// **What the run has to say about itself**, printed when it finishes — the two paths a person
/// needs to open, then the report's own lines.
///
/// **The lines are the report's and not this function's**, so that what a run says is something
/// a test can hold. It was the one part of `call-from-alignments` a mutation could change with
/// the whole suite still green.
pub fn print_report(calls: &Path, parameters_at: &Path, report: &RunReport<'_>) {
    println!("calls: {}", calls.display());
    println!("parameters: {}", parameters_at.display());
    for line in report.lines() {
        println!("{line}");
    }
}
