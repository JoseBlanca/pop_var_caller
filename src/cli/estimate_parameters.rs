//! `estimate-parameters` — fitting a cohort's numbers from its censuses, and writing them down.
//!
//! **This is the file a calling run scores with.** Until it exists a run has two choices and
//! neither is a fit: the constants compiled into the binary, or a parameters file somebody hands
//! it. This command produces one from the cohort's own data — the per-library sequencing error
//! rates and base-quality calibrations, the population's allele-frequency density and the
//! genotype prior seeded from it, each library's contamination, and the repeat-tract slippage
//! ladder.
//!
//! # What it reads, and what it only checks
//!
//! It reads each psp's **census** — the small object in the file's tail — and none of its
//! records. From each header it takes the ground the walk covered and the repeat-tract criteria it
//! cut that ground with (`psp_census_pair.md` §5, §6); the census itself is read from the trailer
//! the same header's footer points at, so nothing has to be paired with anything.
//!
//! # Why the reference and the catalog
//!
//! **A census stores a repeat tract by its index within its stratum and nothing else** — no
//! coordinate, no stratum — so the order has to be rebuilt by choosing the selection again, which
//! is a function of the seed, the reference, the analysed ground and the catalog. A selection
//! rebuilt from another reference or another catalog would index one stratum's tracts as
//! another's, which is a wrong answer rather than a failure, so nothing about the rebuild is
//! trusted:
//!
//! - **the reference and the catalog are compared with the psp headers first**, and refused as
//!   the wrong file, because that fix is to name the right one and it rebuilds nothing;
//! - **then every census's recorded settings are compared with the ones this run's selection
//!   records under**, and every sample that differs is named with the first setting that does and
//!   the command that regenerates it (`psp_census_pair.md` plan step C5).
//!
//! # The inbreeding coefficient, and the three things it can be
//!
//! **The file this writes is the only way a coefficient reaches a calling run** — the calling
//! commands take no flag for it — so the ladder is resolved here, once, and written with the
//! warrant that says which rung it came from
//! ([`DeclaredInbreeding::of_each_sample_over`](crate::calling::parameters_file::DeclaredInbreeding::of_each_sample_over)):
//!
//! - **`supplied`** — `--inbreeding` was given. It overrides the fit, because a user who knows how
//!   their material was bred knows it whatever the cohort size (owner, 2026-08-27).
//! - **`fitted_here`** — nothing was stated, and this cohort's fit measured each sample's
//!   homozygote excess: how much less heterozygous it is than the fit's own allele frequencies
//!   predict. **Circular, and the run says so** rather than hiding it — it is the only fitted
//!   source since the runs-of-homozygosity estimator was removed.
//! - **`defaulted`** — nothing was stated and nothing could be fitted, so zero. **A single-sample
//!   run lands here even though the fit ran**, because one genome's totals cannot identify an
//!   excess and it comes out zero whatever the truth.

use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;

use crate::calling::parameters_file::{DeclaredInbreeding, ParametersFile};
use crate::cli::generate_psps::PSP_FILE_EXTENSION;
use crate::cli::psp_inputs::{PspArgumentRefusal, psps_named};
use crate::cli::run_ground::{self, GroundError};
use crate::parameter_estimation::joint::fit::JointFitConfig;
use crate::parameter_estimation::joint::loci::{ReferenceDigest, SelectionError, UnambiguousRuns};
use crate::parameter_estimation::joint::ssr_fit::{DEFAULT_STRATA_AT_ONCE, SsrFitConfig};
use crate::parameter_estimation::progress::StageProgress;
use crate::reference_info::{
    ReferenceCheck, ReferenceInfoError, read_reference_observing_or_creating_fai,
};
use crate::repeat_catalog::RepeatCatalogHeader;
use crate::run::{
    CensusCohortError, CensusPlan, CensusSelection, CensusesToRegenerate, CohortFitError,
    OpenPspCohort, RunError, THE_COMMAND_THAT_REBUILDS_A_CENSUS, each_census_in_the_cohorts_psps,
    every_read_group_pooled, fit_a_cohort, fitted_inbreeding_of, parameters_file_of,
    parameters_from_the_fit, read_groups_of, the_censuses_as_one_cohort,
    what_the_heads_say_about_every_census_in_a_cohort,
    what_the_run_says_about_every_census_in_a_cohort,
};
use crate::types::{ContigId, InbreedingF, Ploidy};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line.
pub const SUBCOMMAND: &str = "estimate-parameters";

/// Fit a cohort's parameters from its censuses and write them as a parameters file.
#[derive(Debug, Args)]
pub struct EstimateParametersArgs {
    /// Reference FASTA — the one every census's samples were aligned to.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog the psps were walked against. Defaults to
    /// `<reference>.repeats.parquet`.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One psp per sample, or a directory holding them. Repeat the flag.
    ///
    /// **Its records are not read.** What the fit reads is the census in each file's trailer, and
    /// what it reads from each header is the ground the walk covered and the criteria it cut that
    /// ground with (`psp_census_pair.md` §5, §6).
    #[arg(long = "psp", required = true, num_args = 1..)]
    pub psps: Vec<PathBuf>,

    /// Where to write the parameters file.
    #[arg(long)]
    pub output: PathBuf,

    /// Overwrite the output if it is already there.
    #[arg(long)]
    pub force: bool,

    /// How many copies of each chromosome every sample carries.
    #[arg(long, default_value_t = 2)]
    pub ploidy: u8,

    /// The inbreeding coefficient to record for every sample, **overriding the fitted one**.
    ///
    /// **Omit it and the fit's own number is written**: each sample's homozygote excess, how much
    /// less heterozygous it is than the allele frequencies this cohort's own fit predicts. The
    /// file marks that `fitted_here`, and marks what is stated here `supplied`, so a reader can
    /// tell them apart.
    ///
    /// **A stated coefficient wins, and that is a ruling** (owner, 2026-08-27): a user who knows
    /// how their material was bred knows it whatever the cohort size. **It is also why this is an
    /// option and not a default of zero** — `--inbreeding 0` is *these plants are not inbred*,
    /// which is a claim, and saying nothing is *use what you measured*, which is a different one.
    #[arg(long)]
    pub inbreeding: Option<f64>,

    /// How many repeat-tract strata the fit works on at once; 1 uses the least memory.
    ///
    /// A stratum is every repeat tract with one motif length and one reference copy number, and
    /// the fit estimates each stratum's slippage on its own. **This decides how much memory that
    /// step needs.** Fitting one stratum holds one table, with a row for each tract and each
    /// sample with reads there, 91 numbers wide (one for each pair of allele lengths a sample can
    /// carry): for a 5,000-tract stratum, about 0.34 GiB at 100 samples and 7.4 GiB at 2,169.
    /// N strata at once hold up to N of those tables.
    ///
    /// At 1, the default, the fit works on one stratum at a time with every thread on its tracts.
    /// Above 1, N threads each fit a stratum of their own, and the other cores are idle during
    /// this step when N is below the core count. The step's first progress line prints the
    /// largest stratum's table, so a run can choose N from that number. No value changes a fitted
    /// number.
    #[arg(long, value_name = "N", default_value_t = DEFAULT_STRATA_AT_ONCE)]
    pub str_param_estimates_at_once: NonZeroUsize,

    /// Do not estimate contamination: every read group is taken as uncontaminated.
    ///
    /// Contamination is the share of a library's reads that came from another individual.
    /// Estimating it holds, at each position where the cohort varies, an expected allele count for
    /// every sample and four read counts for every read group; a run that estimates it prints what
    /// that takes on its `contamination:` progress line. With this option the parameters file
    /// carries no contamination table, which a calling run reads as "no read group is
    /// contaminated" — the same model it uses when nothing was found, not a fraction measured at
    /// zero. Use it when the samples are known to be clean or when contamination would not change
    /// the calls you need.
    #[arg(long)]
    pub skip_contamination: bool,
}

/// Everything that can stop an `estimate-parameters` run.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum EstimateParametersCliError {
    /// The reference could not be read.
    #[error("reading the reference {}", path.display())]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The reference holds no unambiguous stretch to select positions from.
    #[error("choosing this run's census positions")]
    CensusGround {
        /// What the selection said.
        #[source]
        source: SelectionError,
    },

    /// The ground could not be built.
    #[error("working out what ground this run covers")]
    Ground(#[from] GroundError),

    /// The selection could not be rebuilt.
    #[error("rebuilding the census selection these censuses were written against")]
    CensusNotPlanned {
        /// What the plan said.
        #[source]
        source: Box<RunError>,
    },

    /// A `--psp` naming a directory could not be listed.
    #[error("listing the psps in {}", path.display())]
    PspDirectory {
        /// The directory.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },

    /// A `--psp` directory holds no psp.
    #[error(
        "--psp {} holds no .{PSP_FILE_EXTENSION} file; name the files themselves, or point at \
         the --output-dir a generate-psps run wrote",
        path.display()
    )]
    NoPspsInDirectory {
        /// The directory.
        path: PathBuf,
    },

    /// The cohort's psps could not be opened as one cohort.
    #[error("opening the cohort of psps")]
    PspCohort {
        /// What the opener said.
        #[source]
        source: Box<RunError>,
    },

    /// `--reference` is not the reference the cohort's psps were walked against.
    ///
    /// **Refused before the census selection is rebuilt, and in words that blame the reference.**
    /// Left to the comparison of recorded settings, every census would be named as stale and the
    /// person told to regenerate — against this reference, which costs a quarter of an hour a
    /// sample at whole-genome scale and leaves every census refused again by the next fit with the
    /// right one.
    #[error(
        "{} is not the reference these psps were walked against; they name {walked_against}, so \
         run this fit again with that one, which regenerates nothing",
        path.display()
    )]
    WalkedAgainstAnotherReference {
        /// The FASTA this run was given.
        path: PathBuf,
        /// What the psps call the reference they were walked against — a file name and no
        /// directory, which is all a psp header holds.
        walked_against: String,
        /// The first psp that disagrees with it, and how.
        #[source]
        source: Box<RunError>,
    },

    /// The catalog this run reads is not the one the cohort's psps were walked with.
    ///
    /// **Refused before a segment is cut, for the reason a wrong reference is**: a census judged
    /// against a selection rebuilt from another catalog is named as stale, and regenerating it
    /// against that catalog fixes nothing. `difference` says which part of the catalog's header
    /// is not theirs.
    #[error(
        "the repeat catalog {} is not the one these psps were walked with: {difference}; run \
         this fit again with --catalog naming the one they were, which regenerates nothing",
        path.display()
    )]
    WalkedWithAnotherCatalog {
        /// The catalog this run read — the one `--catalog` named, or the one beside the reference.
        path: PathBuf,
        /// How it differs from the psps' own, as a clause.
        difference: String,
    },

    /// Some of the cohort's psps carry no census this build can fit from, or carry one recorded
    /// under settings this run does not use.
    ///
    /// **The message is the whole of what the user acts on**, so it is the error's own
    /// (`psp_census_pair.md` §4.3): a line a sample, and the command that rebuilds them.
    ///
    /// **Several lines, and no `#[source]`, where an error in this tree is usually one lowercase
    /// clause over a cause.** That is deliberate: this message is the product rather than a frame
    /// around a fault, and there is no underlying error to carry — the psps were read
    /// successfully and what they hold is the news.
    #[error("{report}")]
    CohortCannotBeFitted {
        /// The samples to regenerate, those whose psps would not read, and the command line.
        report: Box<CensusesToRegenerate>,
    },

    /// The cohort's censuses could not be read.
    #[error("reading the censuses in the psps")]
    Cohort {
        /// What the cohort said.
        #[source]
        source: Box<CensusCohortError>,
    },

    /// The cohort could not be fitted.
    #[error("fitting the cohort")]
    Fit {
        /// What the fit said.
        #[source]
        source: Box<CohortFitError>,
    },

    /// The ploidy or the inbreeding coefficient is not a value this run can use.
    #[error("{what} is not a value a run can be given: {value}")]
    NotAValue {
        /// Which flag.
        what: &'static str,
        /// What was typed.
        value: String,
    },

    /// The output is already there and `--force` was not given.
    #[error("{} is already there; pass --force to replace it", path.display())]
    OutputAlreadyThere {
        /// The file that would be replaced.
        path: PathBuf,
    },

    /// The output could not be written.
    #[error("writing {}", path.display())]
    Output {
        /// The path.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },
}

/// Fit `--psp`'s cohort and write `--output`.
///
/// # Errors
///
/// See [`EstimateParametersCliError`]. The output path is judged before a census is opened, so a
/// run that cannot write its answer says so before spending the fit.
pub fn run_estimate_parameters(
    args: &EstimateParametersArgs,
) -> Result<(), EstimateParametersCliError> {
    let (file, samples) = fit_and_assemble(args)?;

    let toml = file.to_toml();
    std::fs::write(&args.output, &toml).map_err(|source| EstimateParametersCliError::Output {
        path: args.output.clone(),
        source,
    })?;
    println!(
        "fitted {} sample{} and wrote {} bytes to {}",
        samples,
        if samples == 1 { "" } else { "s" },
        toml.len(),
        args.output.display(),
    );
    Ok(())
}

/// The SNP/indel fit's settings this run asked for: the defaults, at the run's ploidy, and with
/// contamination estimated unless `--skip-contamination` said not to.
fn ordinary_position_config(args: &EstimateParametersArgs, ploidy: Ploidy) -> JointFitConfig {
    JointFitConfig {
        ploidy,
        estimate_contamination: !args.skip_contamination,
        ..JointFitConfig::default()
    }
}

/// The repeat-tract fit's settings this run asked for: the defaults, with the number of strata
/// fitted at once from `--str-param-estimates-at-once`.
fn repeat_tract_config(args: &EstimateParametersArgs) -> SsrFitConfig {
    SsrFitConfig {
        strata_at_once: args.str_param_estimates_at_once,
        ..SsrFitConfig::default()
    }
}

/// The run itself, with the writing left to the caller — which is what lets a test read the file
/// rather than parse it back off disk.
fn fit_and_assemble(
    args: &EstimateParametersArgs,
) -> Result<(ParametersFile, usize), EstimateParametersCliError> {
    let ploidy =
        Ploidy::try_new(args.ploidy).map_err(|_| EstimateParametersCliError::NotAValue {
            what: "--ploidy",
            value: args.ploidy.to_string(),
        })?;
    let declared = match args.inbreeding {
        None => DeclaredInbreeding::nothing_said(),
        Some(coefficient) => DeclaredInbreeding::one_value_for_every_sample(
            InbreedingF::try_new(coefficient).map_err(|_| {
                EstimateParametersCliError::NotAValue {
                    what: "--inbreeding",
                    value: coefficient.to_string(),
                }
            })?,
        ),
    };
    // **Before a census is opened**, so a run that cannot write its answer does not spend a fit
    // finding out.
    if !args.force && args.output.exists() {
        return Err(EstimateParametersCliError::OutputAlreadyThere {
            path: args.output.clone(),
        });
    }

    let paths = psps_named_by(args)?;
    // **The psps are opened as a cohort, by the opener every psp-taking command shares**: it
    // reads each header and refuses a set that was not walked as one — two files naming one
    // sample, or walked over different ground, or under a different catalog or different
    // repeat-tract criteria (`psp_census_pair.md` §6). No block is decoded by any of it.
    let mut cohort =
        OpenPspCohort::open(&paths).map_err(|source| EstimateParametersCliError::PspCohort {
            source: Box::new(source),
        })?;
    // **Every psp of the cohort is judged, and only then does the run stop** (spec §4.1, §4.3).
    // A cohort with three stale censuses is a message about three samples and not about the first
    // one: rebuilding one census is a quarter of an hour, so a refusal that named them one at a
    // time would cost that wait once a stale sample, in series, to learn a job that fits in one
    // message.
    //
    // **And before the reference is opened**, which is what makes it immediate: judging a psp is
    // one seek and ten bytes, where reading a human reference is minutes. The third cause a
    // census can be stale for — recorded under other settings than this run's — is the one this
    // cannot see, and it is caught below, once the reference has been read and the selection
    // rebuilt.
    if let Some(report) = CensusesToRegenerate::of(
        what_the_heads_say_about_every_census_in_a_cohort(&mut cohort),
        || the_command_that_regenerates(args),
    ) {
        return Err(EstimateParametersCliError::CohortCannotBeFitted {
            report: Box::new(report),
        });
    }
    // **The censuses come out of the psps' trailers**, read lazily: their headers and directories
    // now, a section when the fit asks for it (spec §5).
    //
    // **Each one apart, and not yet a cohort.** Assembling them compares the samples with each
    // other, and a refusal from that names two samples and nothing to do. Judged against this
    // run's own settings once the selection is rebuilt, below, every stale sample is named and so
    // is the fix (plan step C5).
    let setup = StageProgress::begin(format!(
        "reading the censuses of {} psp(s) and the reference",
        paths.len()
    ));
    let censuses = each_census_in_the_cohorts_psps(&cohort).map_err(|source| {
        EstimateParametersCliError::Cohort {
            source: Box::new(source),
        }
    })?;

    // **Read with an observer**, because the selection has to know where the reference is
    // sequence at all: a position inside a run of `N` has no base to compare a read against.
    let mut callable = UnambiguousRuns::default();
    let with_checksums = std::sync::Arc::new(
        read_reference_observing_or_creating_fai(
            args.reference.clone(),
            ReferenceCheck::VerifyAgainstIndex,
            &mut callable,
        )
        .map_err(|source| EstimateParametersCliError::Reference {
            path: args.reference.clone(),
            source,
        })?,
    );
    let unambiguous = callable
        .into_selectable()
        .map_err(|source| EstimateParametersCliError::CensusGround { source })?;
    // **The reference against the psps, before anything is built from it.** A selection rebuilt
    // from another reference differs from every census, and the comparison below would name each
    // one as stale with the fix *regenerate* — the wrong fix, since rebuilt against this reference
    // every census would be refused again by the next fit with the right one. The catalog is
    // opened after this for the same reason: a catalog built on the right reference would be
    // refused against a wrong one, blaming the one file that is correct.
    cohort
        .refuse_a_reference_it_was_not_walked_against(&with_checksums)
        .map_err(
            |source| EstimateParametersCliError::WalkedAgainstAnotherReference {
                path: args.reference.clone(),
                walked_against: cohort
                    .the_reference_the_psps_were_walked_against()
                    .to_string(),
                source: Box::new(source),
            },
        )?;

    // **The ground and what it is cut with come out of the cohort's psps**, not off this
    // command line (spec §6). The five flags that used to supply the criteria were values a
    // person had to retype from a walk that had already recorded them; they are gone, and there
    // is nothing here to type wrongly.
    //
    // **They are the first psp's, and what makes one psp enough is the cohort opener**: it
    // compares all three settings across every file and refuses a cohort that disagrees, naming
    // the two samples and the field (step B3). Take that away and this line is the hazard spec §6
    // describes — a selection built from the first psp's settings that the rest cannot match.
    //
    // **And a census written under other criteria is still named**, below: its recorded settings
    // include the criteria its positions were chosen under, and they are compared with this
    // run's. The fit's own comparison of the positions kept catches a difference only where it
    // moves one of them — on tomato about 1 position in 400 is kept (spec §2) — so it is the
    // backstop and not the check.
    let analysed = cohort.analysed_regions().clone();
    let criteria = cohort.segmentation_inputs().repeat_tract_criteria.clone();
    let catalog_path = run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
    let catalog = run_ground::open_catalog(&catalog_path, &args.reference, &with_checksums)?;
    // **The catalog against the psps, before a segment is cut** (spec §6). Cut with the psps'
    // criteria, a catalog that is not theirs can fail as one too coarse to serve them — true, and
    // a message about the file's contents where the fault is which file was named.
    refuse_a_catalog_the_psps_were_not_walked_with(
        &catalog_path,
        catalog.header(),
        &cohort.segmentation_inputs().catalog,
    )?;
    let segmentation = run_ground::segments_cut_from(
        &catalog,
        &catalog_path,
        &criteria,
        run_ground::CriteriaSource::ThePspHeaders,
        &analysed,
    )?;
    let plan = CensusPlan::of_run(
        CensusSelection::SHIPPED,
        &catalog,
        &analysed,
        &unambiguous,
        &with_checksums,
        &segmentation.inputs().repeat_tract_criteria,
    )
    .map_err(|source| EstimateParametersCliError::CensusNotPlanned {
        source: Box::new(source),
    })?;
    drop(catalog);

    // **Every census against the settings this run records under, and every stale one named**
    // (spec §4.1, §4.2's third row). The reference and the catalog were just proven to be the
    // psps' own, so a census that differs from the run differs from its own psp or was written by
    // a build that chooses positions differently — and regenerating is the fix for both.
    //
    // **After the reference is read, where the head's refusal comes before it**, because the run's
    // settings include a digest of the reference's bases and of the positions it keeps. So a
    // cohort with one psp carrying no census and another recorded under a different budget is
    // refused for the first alone, before the reference is opened — and that costs the person
    // nothing, because the command that refusal prints covers every psp given and
    // `regenerate-census` rebuilds both: its own skip rule is all three causes, the selection
    // included, since it rebuilds the selection anyway (spec §8).
    if let Some(report) = CensusesToRegenerate::of(
        what_the_run_says_about_every_census_in_a_cohort(
            &cohort,
            &censuses,
            &plan.recording_terms(),
        ),
        || the_command_that_regenerates(args),
    ) {
        return Err(EstimateParametersCliError::CohortCannotBeFitted {
            report: Box::new(report),
        });
    }
    // **What assembling can still refuse is damage**: every census now records this run's
    // settings, so two cannot disagree on them, and two psps declaring one read group were refused
    // when the cohort was opened, from their headers.
    let mut evidence = the_censuses_as_one_cohort(censuses).map_err(|source| {
        EstimateParametersCliError::Cohort {
            source: Box::new(source),
        }
    })?;

    let contigs = std::sync::Arc::clone(&plan.contigs);
    let contig_of = move |name: &str| {
        contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|index| ContigId(index as u32))
    };
    setup.always(|into| {
        format!(
            "censuses and reference read, {into}; {} sample(s) holding {} read group(s)",
            evidence.len(),
            evidence.read_groups().len()
        )
    });
    let pooled = every_read_group_pooled(&evidence);
    let fit = fit_a_cohort(
        &mut evidence,
        &plan.loci,
        &contig_of,
        &pooled,
        &ordinary_position_config(args, ploidy),
        &repeat_tract_config(args),
    )
    .map_err(|source| EstimateParametersCliError::Fit {
        source: Box::new(source),
    })?;

    let samples = evidence.len();
    let read_groups = read_groups_of(&evidence, cohort.paths());
    // **One list, resolved once, handed to both the run and the file.** The ladder is
    // `DeclaredInbreeding`'s: what the operator stated, else what this fit measured, else the
    // default — joined to the run's samples by name, and carrying which of the three it was.
    let inbreeding =
        declared.of_each_sample_over(&read_groups, &fitted_inbreeding_of(&fit.generic.hom_excess));
    let parameters = parameters_from_the_fit(&fit, &evidence, &pooled, &inbreeding, ploidy);
    let terms = evidence
        .terms()
        .expect("a cohort of one or more samples records terms")
        .clone();
    let reference = ReferenceDigest::of(&with_checksums).map_err(|source| {
        EstimateParametersCliError::CensusNotPlanned {
            source: Box::new(RunError::CensusNotPlanned {
                source: Box::new(source),
            }),
        }
    })?;
    let file = parameters_file_of(
        &parameters,
        &read_groups,
        &inbreeding,
        &reference,
        &terms,
        &segmentation.inputs().repeat_tract_criteria,
    );
    Ok((file, samples))
}

/// **The command that rebuilds this run's stale censuses, in the words it was given** — spec
/// §4.3's last line, written so it can be copied rather than reconstructed.
///
/// **The `--psp` arguments are the user's own, not the expanded list.** A run over a directory of
/// sixty psps typed one path, and a report that answered with sixty is a report nobody copies.
///
/// **`--reference` and `--catalog` come too, because rebuilding a census needs them**: the census
/// stores a repeat tract by its index within its stratum, so the selection has to be built again,
/// and that reads the reference and the catalog (spec §8). `--catalog` is left out when this run
/// was not given one, since the command finds the same file beside the reference the same way.
fn the_command_that_regenerates(args: &EstimateParametersArgs) -> String {
    use std::fmt::Write as _;

    let mut line = String::from(THE_COMMAND_THAT_REBUILDS_A_CENSUS);
    let _ = write!(line, " --reference {}", as_typed(&args.reference));
    if let Some(catalog) = &args.catalog {
        let _ = write!(line, " --catalog {}", as_typed(catalog));
    }
    for psp in &args.psps {
        let _ = write!(line, " --psp {}", as_typed(psp));
    }
    line
}

/// **Refuse a catalog the psps were not walked with**, saying how it differs from theirs (spec §6).
///
/// **The comparison is the catalog header's own**
/// ([`RepeatCatalogHeader::first_difference`]), so this command and `regenerate-census` name the
/// same difference in the same words; what is this command's is the sentence around it, which says
/// what to do. `call-from-psps` makes the same comparison through
/// [`SegmentationInputs::first_difference`](crate::run::SegmentationInputs::first_difference),
/// which names only *the repeat catalog* — enough for a run that will not proceed either way, and
/// not enough for a person deciding which of two files they hold.
///
/// **The digest clause cannot fire here**, and the assertion says why: this run's catalog was
/// checked against this run's reference as it opened, the psps' catalog was checked against the
/// walk's reference when the walk opened it, and the check above proved those two references are
/// one. What can still differ is the contig *table*, which also carries the FASTA's line
/// geometry — the same bases wrapped at another width.
fn refuse_a_catalog_the_psps_were_not_walked_with(
    path: &Path,
    given: &RepeatCatalogHeader,
    walked_with: &RepeatCatalogHeader,
) -> Result<(), EstimateParametersCliError> {
    debug_assert_eq!(
        given.reference_md5, walked_with.reference_md5,
        "each catalog was checked against the reference of the run that opened it, and those two \
         references were just proven to be one, so a difference here would mean one of those \
         three checks did not run",
    );
    match given.first_difference(walked_with) {
        None => Ok(()),
        Some(difference) => Err(EstimateParametersCliError::WalkedWithAnotherCatalog {
            path: path.to_path_buf(),
            difference,
        }),
    }
}

/// A path as it has to appear on a command line that will be pasted into a shell.
///
/// **Quoted when it holds a space**, because the whole value of that line is that it is copied,
/// and a directory called `tomato run 3` pasted bare is three arguments and a different meaning.
fn as_typed(path: &std::path::Path) -> String {
    let shown = path.display().to_string();
    match shown.contains(char::is_whitespace) {
        true => format!("'{shown}'"),
        false => shown,
    }
}

/// **The psps this run fits, with every directory expanded** — the rule every psp-taking command
/// shares ([`psps_named`]), with its refusals dressed in this command's words.
fn psps_named_by(
    args: &EstimateParametersArgs,
) -> Result<Vec<PathBuf>, EstimateParametersCliError> {
    psps_named(&args.psps).map_err(|refusal| match refusal {
        PspArgumentRefusal::Unlistable { path, source } => {
            EstimateParametersCliError::PspDirectory { path, source }
        }
        PspArgumentRefusal::Empty { path } => {
            EstimateParametersCliError::NoPspsInDirectory { path }
        }
    })
}
