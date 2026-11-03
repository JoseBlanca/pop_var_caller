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
    ParametersFileError, ReadsBehindEachCalibration, SnpIndelFitStart,
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
/// How far it can be divided is set by the merge's own floor
/// ([`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`], 500 bases), below which a region is too
/// narrow to be worth a fan-out and a barrier, and by the pool: there is nothing to gain from
/// more regions in flight than there are threads to work them.
///
/// # What it gives, and where it gives nothing
///
/// | samples | ground | regions in flight | width |
/// |---|---|---|---|
/// | 8 | 16,000 | 18 | 888 |
/// | 63 | 7,936 | 15 | 529 |
/// | 1,000 | 500 | 1 | 500 |
///
/// **Above about a thousand samples it hands back one region and today's behaviour**, because
/// `round_width_for` is already at its floor there and there is no ground left to divide. That
/// is the honest degradation and it is the one this rule chooses: a large cohort keeps the
/// memory it has and gives up the parallelism, rather than keeping the parallelism and
/// multiplying the memory by the core count. Nothing at that end of the range has been
/// measured; what has been is 1 through 63 accessions on the tomato benchmark, where peak
/// resident rises between 4% and 13% and the calling pass falls by up to 5.5×.
#[must_use]
pub fn round_shape_for(
    ground: CohortLocusBuilderRegionsLen,
) -> (
    CohortLocusBuilderRegionsLen,
    CohortLocusBuilderRegionsInFlight,
) {
    let ground = ground.get();
    let floor = DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN;
    let workers = rayon::current_num_threads().max(1);
    // PANIC-FREE: `ground >= floor` because `round_width_for` clamps to it, so the division is
    // at least one, and `workers` is at least one.
    let in_flight = usize::try_from(ground / floor)
        .unwrap_or(usize::MAX)
        .clamp(1, workers);
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

/// **The fewest checked sites a sample's fill disagreement rate is stated from**: below it, the
/// table gives the counts and the log counts the sample as not measured. At 30 sites a rate of one
/// in ten is known to within about a factor of two.
pub const MIN_FILL_CHECK_SITES: u32 = 30;

/// **Reports what reads did to the phase** (read phase spec §2.6, §2.8): phase from reads alone where
/// it ran, and where the neighbour-evidence pass ran, the read check of its population links — in
/// total in the log, and per sample at debug level.
pub fn report_read_phase(finished: &FinishedVcf) {
    if let Some(tally) = finished.read_phase {
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
         reach at odds of 199 to 1 with GQ 30 or more at both ends, {} of {} contradicted by the \
         reads, a check of the population links that does not come from the model, for links \
         under 1 kb and not acted on; {} population links broken by reads, {} links made by reads, \
         {} heterozygotes the consensus did not lean on attached by reads",
        total.population_links_contradicted,
        total.population_links_reached,
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
        run_output::with_thousands(tally.sites),
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
    /// Do not run.
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

/// **Whether this run uses the neighbour-evidence pass, and why**, from the option, each sample's
/// inbreeding coefficient in the run's sample order, and the ploidy — taken once, before the
/// header is written, and stated in it either way. `path_agreement` is
/// `--phase-min-path-agreement`.
#[must_use]
pub fn neighbour_evidence_decision(
    asked: RefineGenotypesArg,
    inbreeding_by_sample: &[InbreedingF],
    ploidy: Ploidy,
    path_agreement: PathAgreement,
) -> NeighbourEvidenceDecision {
    let switch = match asked {
        RefineGenotypesArg::Auto => NeighbourEvidenceSwitch::Auto,
        RefineGenotypesArg::On => NeighbourEvidenceSwitch::On,
        RefineGenotypesArg::Off => NeighbourEvidenceSwitch::Off,
    };
    let inbreeding: Vec<f64> = inbreeding_by_sample.iter().map(|f| f.get()).collect();
    let decision = decide(switch, &inbreeding, ploidy).with_path_agreement(path_agreement);
    tracing::info!(
        "neighbour-evidence pass: {} ({})",
        if decision.ran { "runs" } else { "off" },
        decision.why
    );
    decision
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
/// output) (`doc/devel/ng/spec/command_line.md` §7.5). Removed when the returned value is
/// dropped.
///
/// **A directory of its own**, so several runs can share one `--tmp-dir` without colliding, and
/// removing it never touches anything the user put there.
///
/// # Errors
///
/// [`CallingRunError::ScratchDirectory`] when it cannot be made.
pub fn scratch_directory_for(
    output: &Path,
    tmp_dir: Option<&Path>,
) -> Result<tempfile::TempDir, CallingRunError> {
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
    tracing::debug!("scratch directory: {}", made.path().display());
    Ok(made)
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
            let (width, in_flight) = round_shape_for(ground);
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
            let digest = Md5::digest(&bytes);
            let hex: String = digest.iter().map(|byte| format!("{byte:02x}")).collect();
            format!("{name} (md5 {hex})")
        }
        Err(_) => name.into_owned(),
    }
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
