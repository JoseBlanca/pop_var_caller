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
//! It reads the **censuses**. It does not read the psps — but each census's psp has to be beside
//! it, because a census names the psp it was built from and evidence from other reads is
//! otherwise indistinguishable from this run's. What is taken from each psp is its header: one
//! short read, for the digest the census names it by and for the ground the walk covered.
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
use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::ng::repeat_catalog::RepeatCatalog;
use crate::ng::run::{
    CensusCohortError, CensusPlan, CensusSelection, CohortFitError, RunError,
    every_read_group_pooled, fit_a_cohort, open_census_cohort, parameters_file_of,
    parameters_from_the_fit, read_groups_of,
};
use crate::ng::types::{ContigId, InbreedingF, Ploidy};
use crate::pop_var_caller_exp::generate_census::CENSUS_FILE_EXTENSION;
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

    /// One census per sample, or a directory holding them. Repeat the flag.
    ///
    /// **Each census's psp must be beside it**, under the same stem. It is not read: its header
    /// is, for the digest the census names it by and for the ground its walk covered.
    #[arg(long = "census", required = true, num_args = 1..)]
    pub censuses: Vec<PathBuf>,

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

    /// The fewest motif copies a tract needs before this run treats it as a repeat: six
    /// comma-separated numbers, one per period 1 to 6. **Give the values the psps were walked
    /// under** — they decide which stretches the selection may keep.
    #[arg(
        long,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_copies,
        default_value = "8,6,6,6,5,4",
        help_heading = "What counts as a repeat"
    )]
    pub min_copies: MinCopies,

    /// The shortest repeat unit this run treats as a repeat.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=crate::ng::types::MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub min_period: u8,

    /// The longest repeat unit this run treats as a repeat.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=crate::ng::types::MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub max_period: u8,

    /// A tract longer than this many bases is a satellite.
    #[arg(long, default_value_t = DEFAULT_MAX_STR_LEN, help_heading = "What counts as a repeat")]
    pub max_str_len: u64,

    /// The least share of a tract's bases that must match its motif exactly.
    #[arg(long, default_value_t = DEFAULT_MIN_PURITY, help_heading = "What counts as a repeat")]
    pub min_purity: f32,
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

    /// A `--census` naming a directory could not be listed.
    #[error("listing the censuses in {}", path.display())]
    CensusDirectory {
        /// The directory.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },

    /// A `--census` directory holds no census.
    #[error("{} holds no .census file", path.display())]
    NoCensusesInDirectory {
        /// The directory.
        path: PathBuf,
    },

    /// The cohort could not be opened.
    #[error("opening the cohort of censuses")]
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

/// Fit `--census`'s cohort and write `--output`.
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

    let paths = censuses_named_by(args)?;
    let mut open =
        open_census_cohort(&paths).map_err(|source| EstimateParametersCliError::Cohort {
            source: Box::new(source),
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

    let ground = run_ground::GroundRequest {
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
    };
    // **The ground is the psps' own**, not this command's: the censuses were written over it, and
    // a selection rebuilt over anything else is refused by the fit.
    let analysed = open.analysed_regions.clone();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;
    let catalog =
        RepeatCatalog::open_checking_against_reference(&ground.catalog_path(), &with_checksums)
            .map_err(|source| {
                EstimateParametersCliError::Ground(GroundError::Catalog {
                    path: ground.catalog_path(),
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
    let pooled = every_read_group_pooled(&open.evidence);
    let fit = fit_a_cohort(
        &mut open.evidence,
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

    let samples = open.evidence.len();
    let read_groups = read_groups_of(&open.evidence, &open.samples);
    let stated: Vec<Estimate<InbreedingF>> = (0..samples)
        .map(|_| Estimate {
            value: inbreeding,
            provenance: Provenance::Supplied,
            observations: 0,
        })
        .collect();
    let parameters = parameters_from_the_fit(
        &fit,
        &open.evidence,
        &pooled,
        (0..samples).map(|_| inbreeding).collect(),
        ploidy,
    );
    let terms = open
        .evidence
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

/// **The censuses this run fits, with every directory expanded** — one entry a sample, in the
/// order they were given, and a directory's contents sorted by name so two runs naming one
/// directory read the same cohort in the same order.
fn censuses_named_by(
    args: &EstimateParametersArgs,
) -> Result<Vec<PathBuf>, EstimateParametersCliError> {
    let mut paths = Vec::with_capacity(args.censuses.len());
    for named in &args.censuses {
        if !named.is_dir() {
            paths.push(named.clone());
            continue;
        }
        let mut inside = Vec::new();
        for entry in std::fs::read_dir(named).map_err(|source| {
            EstimateParametersCliError::CensusDirectory {
                path: named.clone(),
                source,
            }
        })? {
            let entry = entry.map_err(|source| EstimateParametersCliError::CensusDirectory {
                path: named.clone(),
                source,
            })?;
            let path = entry.path();
            if path
                .extension()
                .is_some_and(|it| it == CENSUS_FILE_EXTENSION)
                && path.is_file()
            {
                inside.push(path);
            }
        }
        if inside.is_empty() {
            return Err(EstimateParametersCliError::NoCensusesInDirectory {
                path: named.clone(),
            });
        }
        inside.sort();
        paths.extend(inside);
    }
    Ok(paths)
}
