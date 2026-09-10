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
//! is a function of the seed, the reference, the analysed ground and the catalog. The rebuild is
//! checked against a digest every census carries before a tract is read: a selection rebuilt from
//! another reference or another catalog would index one stratum's tracts as another's, which is a
//! wrong answer rather than a failure.
//!
//! # What it leaves declared
//!
//! **The inbreeding coefficient**, which is fitted from a sample's own windowed genome histogram
//! — the other pre-pass route, not this one. `--inbreeding` states one for the whole cohort and
//! the file records it as supplied.

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::ng::calling::parameters_file::ParametersFile;
use crate::ng::parameter_estimation::joint::fit::JointFitConfig;
use crate::ng::parameter_estimation::joint::loci::{
    ReferenceDigest, SelectionError, UnambiguousRuns,
};
use crate::ng::parameter_estimation::joint::ssr_fit::SsrFitConfig;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoError, read_reference_observing_or_creating_fai,
};
use crate::ng::repeat_catalog::RepeatCatalog;
use crate::ng::run::{
    CensusCohortError, CensusPlan, CensusSelection, CensusesToRegenerate, CohortFitError,
    OpenPspCohort, RunError, every_census_in_the_cohorts_psps, every_read_group_pooled,
    fit_a_cohort, parameters_file_of, parameters_from_the_fit, read_groups_of,
    what_the_heads_say_about_every_census_in_a_cohort,
};
use crate::ng::types::{ContigId, InbreedingF, Ploidy};
use crate::pop_var_caller_exp::generate_psps::PSP_FILE_EXTENSION;
use crate::pop_var_caller_exp::psp_inputs::{PspArgumentRefusal, psps_named};
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

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

    /// The inbreeding coefficient to record for every sample.
    ///
    /// **Declared, never fitted on this route.** It is fitted from a sample's own windowed genome
    /// histogram, which is the other pre-pass route; the file records whatever is stated here as
    /// supplied, so a reader can tell it from a measurement.
    #[arg(long, default_value_t = 0.0)]
    pub inbreeding: f64,
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

    /// Some of the cohort's psps carry no census this build can fit from.
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
    let inbreeding = InbreedingF::try_new(args.inbreeding).map_err(|_| {
        EstimateParametersCliError::NotAValue {
            what: "--inbreeding",
            value: args.inbreeding.to_string(),
        }
    })?;
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
    // census can be stale for — built against a different set of loci — is the one this cannot
    // see, and it is caught later, by the fit's own digest.
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
    let mut evidence = every_census_in_the_cohorts_psps(&cohort).map_err(|source| {
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
    // **What happens if they are wrong anyway is weaker than it sounds.** The fit compares a
    // digest of the *kept generic positions* against the censuses' own (`fit_a_cohort`), so a
    // difference is caught only where it moves one of them — and on tomato about 1 position in
    // 400 is kept (spec §2), so a criterion that retypes a single short tract can pass that
    // digest and be re-indexed in silence. That is the argument for taking the criteria from the
    // psps rather than a reason to trust the net below.
    let analysed = cohort.analysed_regions().clone();
    let criteria = cohort.segmentation_inputs().repeat_tract_criteria.clone();
    let catalog_path = run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
    let segmentation = run_ground::segments_cut_with(
        &catalog_path,
        &args.reference,
        &criteria,
        run_ground::CriteriaSource::ThePspHeaders,
        &analysed,
        &with_checksums,
    )?;
    let catalog = RepeatCatalog::open_checking_against_reference(&catalog_path, &with_checksums)
        .map_err(|source| {
            EstimateParametersCliError::Ground(GroundError::Catalog {
                path: catalog_path.clone(),
                source,
            })
        })?;
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

    let contigs = std::sync::Arc::clone(&plan.contigs);
    let contig_of = move |name: &str| {
        contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|index| ContigId(index as u32))
    };
    let pooled = every_read_group_pooled(&evidence);
    let fit = fit_a_cohort(
        &mut evidence,
        &plan.loci,
        &contig_of,
        &pooled,
        &JointFitConfig {
            ploidy,
            ..JointFitConfig::default()
        },
        &SsrFitConfig::default(),
    )
    .map_err(|source| EstimateParametersCliError::Fit {
        source: Box::new(source),
    })?;

    let samples = evidence.len();
    let read_groups = read_groups_of(&evidence, cohort.paths());
    let stated: Vec<Estimate<InbreedingF>> = (0..samples)
        .map(|_| Estimate {
            value: inbreeding,
            provenance: Provenance::Supplied,
            observations: 0,
        })
        .collect();
    let parameters = parameters_from_the_fit(
        &fit,
        &evidence,
        &pooled,
        (0..samples).map(|_| inbreeding).collect(),
        ploidy,
    );
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
        &stated,
        &reference,
        &terms,
        &segmentation.inputs().repeat_tract_criteria,
    );
    Ok((file, samples))
}

/// **The command that rebuilds a psp's census**, named here until it exists.
///
/// Plan step D1 turns today's `generate-census` — which writes a census *file* beside a psp,
/// which nothing reads any more — into `regenerate-census`, which replaces the psp's trailer.
/// Until then this report names a command a person cannot yet run. **That is the right name to
/// print rather than the old one**: `generate-census` would leave them with a file this fit does
/// not read and a psp still stale. The constant is here so that D1 has one place to point it at
/// the real subcommand's own name.
const THE_COMMAND_THAT_REBUILDS_A_CENSUS: &str = "regenerate-census";

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
