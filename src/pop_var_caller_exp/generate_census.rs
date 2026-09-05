//! `generate-census` — building a census from psps that already exist.
//!
//! **A census is the small file a parameters fit reads.** It holds what one sample showed at a
//! fixed set of positions and repeat tracts chosen for the whole run — a few thousand loci out
//! of a genome — so that the fit can ask the same question of every sample and compare their
//! answers. A psp holds everything; a census holds the part the fit needs, in the shape it needs
//! it (`doc/devel/ng/spec/parameter_prepass_census_sites.md`).
//!
//! **`generate-psps` already writes one beside each psp, and this does not replace it.** This is
//! the second producer, for the three situations that cannot re-walk the reads
//! (`parameter_prepass_joint_records.md` §6.1): psps written before censuses existed, a census
//! lost or built under settings since changed, and a census wanted larger than the one on disk.
//! The two are held to producing the same file byte for byte, which is what says a psp carries
//! everything a census needs.
//!
//! **There is no `--regions`, for the reason `call-from-psps` has none.** Which ground was
//! walked is recorded in every psp, the cohort is refused unless the files agree about it, and
//! that agreed ground is what the selection is made over. **This is not a convenience**: the
//! digest of the analysed regions travels in every census as one of its recording terms, so a
//! census built over ground the psps were not walked over could not be pooled with the cohort's
//! others — and would say so hours later, at the fit.
//!
//! **The selection is the run's, not the sample's**, and the numbers behind it are
//! [`CensusSelection::SHIPPED`]: about two million positions, five thousand tracts a stratum,
//! and a seed that is a compiled-in constant rather than a clock. Two invocations of this
//! command over one cohort therefore keep the same positions, which is what lets a cohort be
//! built up a sample at a time.
//!
//! **One psp is open at a time while censuses are built.** The whole cohort is opened first,
//! because opening it is what checks that the files agree and what says what ground there is —
//! and then closed, before the first census. Holding a thousand psps open to read them one by
//! one would spend the memory psp mode exists to save.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;

use crate::ng::parameter_estimation::joint::census_file::write_census;
use crate::ng::parameter_estimation::joint::loci::{SelectionError, UnambiguousRuns};
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoError, read_reference_observing_or_creating_fai,
};
use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::ng::repeat_catalog::RepeatCatalog;
use crate::ng::run::report::{describe, plural};
use crate::ng::run::{
    CensusFromPspError, CensusPlan, CensusSelection, CensusTally, OpenPspCohort, RunError,
    census_from_psp,
};
use crate::ng::types::MAX_MOTIF_LEN;
use crate::pop_var_caller_exp::generate_psps::PSP_FILE_EXTENSION;
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line.
pub const SUBCOMMAND: &str = "generate-census";

/// The extension a census file carries, matching what `generate-psps` writes.
pub const CENSUS_FILE_EXTENSION: &str = "census";

/// Build each stored psp's census, without re-reading a single alignment file.
#[derive(Debug, Args)]
pub struct GenerateCensusArgs {
    /// Reference FASTA — the one every psp's samples were aligned to. A `.fai` is built beside
    /// it if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog the psps were walked against. Defaults to
    /// `<reference>.repeats.parquet`.
    ///
    /// **A cohort walked under another catalog is refused**, naming the field that differs: the
    /// catalog decides which stretches are repeat tracts, and a selection made over a different
    /// one would keep tracts these psps hold no observations of.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One psp per sample, or a directory holding them. Repeat the flag.
    ///
    /// A directory contributes every `.psp` file directly inside it, in name order.
    #[arg(long = "psp", required = true, num_args = 1..)]
    pub psps: Vec<PathBuf>,

    /// The directory the census files are written into, one `<sample>.census` per sample.
    /// Created if it does not exist.
    ///
    /// **Name the psps' own directory to put each census beside its psp**, which is where
    /// `generate-psps` leaves them and where somebody looking for the pair will look.
    #[arg(long)]
    pub output_dir: PathBuf,

    /// Overwrite census files that are already in `--output-dir`.
    ///
    /// Without it the run refuses as soon as it finds one, **before a psp is read**, so a cohort
    /// is never left half-replaced.
    #[arg(long)]
    pub force: bool,

    /// The fewest motif copies a tract needs before this run treats it as a repeat: six
    /// comma-separated numbers, one per period 1 to 6.
    ///
    /// **Give the same values the psps were walked under.** They decide which stretches the
    /// segmentation calls repeat tracts, and therefore which loci the selection may keep.
    #[arg(
        long,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_copies,
        default_value = "8,6,6,6,5,4",
        help_heading = "What counts as a repeat"
    )]
    pub min_copies: MinCopies,

    /// The shortest repeat unit this run treats as a repeat. 1 puts homopolymers on the repeat
    /// path.
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

    /// A tract longer than this many bases is a satellite, and no generator speaks for it.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_len: u64,

    /// The least share of a tract's bases that must match its motif exactly.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        help_heading = "What counts as a repeat"
    )]
    pub min_purity: f32,
}

/// Everything that can stop a `generate-census` run.
///
/// **The refusals come before the first psp is read**, in the order a person would want them: the
/// reference, the output directory and what is already in it, then the cohort.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum GenerateCensusCliError {
    /// The reference could not be read.
    #[error("reading the reference {}", path.display())]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The ground the censuses would be selected over could not be built.
    #[error("working out what ground this run covers")]
    Ground(#[from] GroundError),

    /// The reference holds no unambiguous stretch to select positions from.
    #[error("choosing this run's census positions")]
    CensusGround {
        /// What the selection said.
        #[source]
        source: SelectionError,
    },

    /// The selection itself was refused.
    #[error("planning this run's census")]
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
    #[error("{} holds no .psp file", path.display())]
    NoPspsInDirectory {
        /// The directory.
        path: PathBuf,
    },

    /// The cohort could not be opened, or its files do not agree.
    #[error("opening the cohort of psps")]
    Cohort {
        /// What the cohort said.
        #[source]
        source: Box<RunError>,
    },

    /// A sample's name cannot be the file name of its own census.
    #[error(
        "the sample name {sample:?} cannot be a file name, so its census has nowhere to go; \
         the name comes from the psp's own header"
    )]
    SampleNameNotAFileName {
        /// The name the psp declares.
        sample: String,
    },

    /// A census is already at the path this run would write, and `--force` was not given.
    #[error(
        "{} is already there; pass --force to replace the censuses in this directory",
        path.display()
    )]
    CensusAlreadyThere {
        /// The file that would be replaced.
        path: PathBuf,
    },

    /// The output directory could not be created or written in.
    #[error("writing to {}", path.display())]
    OutputDir {
        /// The path that failed.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },

    /// One sample's census could not be built from its psp.
    #[error("building {sample}'s census from {}", psp.display())]
    Build {
        /// The individual, as its psp names it.
        sample: String,
        /// The psp being read.
        psp: PathBuf,
        /// What the producer said.
        #[source]
        source: Box<CensusFromPspError>,
    },

    /// A census was built and could not be encoded.
    #[error("encoding {sample}'s census")]
    CensusNotEncoded {
        /// The individual.
        sample: String,
        /// What the encoder said.
        #[source]
        source: Box<crate::ng::parameter_estimation::joint::census::CensusError>,
    },
}

/// What one sample's census cost and holds.
#[derive(Debug, Clone)]
pub struct SampleCensusOutcome {
    /// The individual, as its psp's header names it.
    pub sample: String,
    /// The psp it was built from.
    pub psp: PathBuf,
    /// How many records that psp holds.
    pub records: u64,
    /// The census written.
    pub census: PathBuf,
    /// How large it is.
    pub census_bytes: u64,
    /// What went into it.
    pub tally: CensusTally,
}

impl SampleCensusOutcome {
    /// The one line that says what this sample produced — **shared by the progress note printed
    /// as the sample finishes and by the report at the end**, so the two cannot come to say
    /// different things about one census.
    #[must_use]
    pub fn line(&self) -> String {
        let mut line = String::new();
        let _ = write!(
            line,
            "{}: {} stored loci read from {}, census {} bytes at {}",
            self.sample,
            self.records,
            self.psp.display(),
            self.census_bytes,
            self.census.display(),
        );
        if self.tally.contributes_nothing() {
            // **Named rather than omitted.** A census that is all denominator is a legitimate
            // outcome — the walk covered ground the selection kept nothing in, or no read
            // reached what it did keep — and a run that said nothing about it would leave
            // somebody hunting for a file written exactly as asked.
            let _ = write!(
                line,
                "; no kept locus has a read in this file, so it contributes nothing to a fit",
            );
            return line;
        }
        let _ = write!(
            line,
            "; reads at {} of {} kept positions and {} of {} kept tracts",
            self.tally.positions_with_reads,
            self.tally.positions_kept,
            self.tally.tracts_with_reads,
            self.tally.tracts_kept,
        );
        line
    }
}

/// What a whole run produced.
#[derive(Debug, Clone)]
pub struct CensusReport {
    /// The ground the psps agree they were walked over.
    pub ground: String,
    /// How many bases of it there are.
    pub analysed_bases: u64,
    /// One entry a sample, in the order the psps were given.
    pub samples: Vec<SampleCensusOutcome>,
}

impl CensusReport {
    /// The report, one line at a time.
    #[must_use]
    pub fn lines(&self) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.samples.len() + 2);
        lines.push(format!(
            "built {} census{} over {} — {} bases analysed",
            self.samples.len(),
            if self.samples.len() == 1 { "" } else { "es" },
            self.ground,
            self.analysed_bases,
        ));
        lines.extend(self.samples.iter().map(SampleCensusOutcome::line));
        let empty = self
            .samples
            .iter()
            .filter(|sample| sample.tally.contributes_nothing())
            .count();
        if empty > 0 {
            lines.push(format!(
                "{} of {} sample{} put nothing into the fit",
                empty,
                self.samples.len(),
                plural(self.samples.len() as u64),
            ));
        }
        lines
    }
}

/// Build a census for every psp named, and print what each one holds.
///
/// # Errors
///
/// See [`GenerateCensusCliError`]. **Every refusal it can make comes before the first psp is
/// read**: the reference, the ground, the selection, the psp paths, the cohort's own agreement,
/// each sample's name, and whether a census is already where this run would write one.
pub fn run_generate_census(args: &GenerateCensusArgs) -> Result<(), GenerateCensusCliError> {
    let report = build_every_census(args)?;
    for line in report.lines() {
        println!("{line}");
    }
    Ok(())
}

/// The run itself, with the printing left to the caller — which is what lets a test read the
/// report rather than parse it back out of a stream.
fn build_every_census(args: &GenerateCensusArgs) -> Result<CensusReport, GenerateCensusCliError> {
    // **The reference is read with an observer**, because the selection has to know where the
    // genome is sequence at all: a position inside a run of `N` has no reference base to compare
    // a read against, and keeping one would put a permanent hole in every sample's records.
    let mut callable = UnambiguousRuns::default();
    let with_checksums = std::sync::Arc::new(
        read_reference_observing_or_creating_fai(
            args.reference.clone(),
            ReferenceCheck::VerifyAgainstIndex,
            &mut callable,
        )
        .map_err(|source| GenerateCensusCliError::Reference {
            path: args.reference.clone(),
            source,
        })?,
    );
    let unambiguous = callable
        .into_selectable()
        .map_err(|source| GenerateCensusCliError::CensusGround { source })?;
    let contigs = with_checksums.contig_list();

    let paths = psps_named_by(args)?;

    // **The cohort is opened to be agreed with, then closed.** Opening it is what refuses files
    // that were walked over different ground or that share an `@RG ID`, and what says what the
    // analysed regions are; holding it open while censuses are built would keep every psp of the
    // cohort open to read them one at a time.
    let (analysed, samples) = {
        let cohort =
            OpenPspCohort::open(&paths).map_err(|source| GenerateCensusCliError::Cohort {
                source: Box::new(source),
            })?;
        let samples: Vec<String> = cohort.sample_names().map(str::to_string).collect();
        (cohort.analysed_regions().clone(), samples)
    };

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
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;

    let catalog =
        RepeatCatalog::open_checking_against_reference(&ground.catalog_path(), &with_checksums)
            .map_err(|source| {
                GenerateCensusCliError::Ground(GroundError::Catalog {
                    path: ground.catalog_path(),
                    source,
                })
            })?;
    // **Chosen once, before the first psp is read**, because the selection is the run's: a run
    // that chose per sample would choose the same set N times, and one whose samples chose
    // differently could not be fitted at all.
    let plan = CensusPlan::of_run(
        CensusSelection::SHIPPED,
        &catalog,
        &analysed,
        &unambiguous,
        &with_checksums,
        &segmentation.inputs().repeat_tract_criteria,
    )
    .map_err(|source| GenerateCensusCliError::CensusNotPlanned {
        source: Box::new(source),
    })?;
    drop(catalog);

    std::fs::create_dir_all(&args.output_dir).map_err(|source| {
        GenerateCensusCliError::OutputDir {
            path: args.output_dir.clone(),
            source,
        }
    })?;
    // **Every name is judged and every collision found before any work**, so a cohort is never
    // left half-replaced by a refusal that arrives at the fortieth sample.
    let mut census_paths = Vec::with_capacity(samples.len());
    for sample in &samples {
        refuse_a_sample_name_that_is_not_a_file_name(sample)?;
        let path = args
            .output_dir
            .join(format!("{sample}.{CENSUS_FILE_EXTENSION}"));
        if !args.force && path.exists() {
            return Err(GenerateCensusCliError::CensusAlreadyThere { path });
        }
        census_paths.push(path);
    }

    let mut built = Vec::with_capacity(paths.len());
    for ((psp, sample), census_path) in paths.iter().zip(&samples).zip(&census_paths) {
        let mut produced = census_from_psp(psp, &plan, &segmentation).map_err(|source| {
            GenerateCensusCliError::Build {
                sample: sample.clone(),
                psp: psp.clone(),
                source: Box::new(source),
            }
        })?;
        let tally =
            produced
                .tally()
                .map_err(|source| GenerateCensusCliError::CensusNotEncoded {
                    sample: sample.clone(),
                    source: Box::new(source),
                })?;

        // **Written to a name of its own and renamed once whole**, the treatment
        // `generate-psps` gives both its files and for the same reason: a run stopped part-way
        // must leave neither a stump at a name a fit would open, nor a destroyed copy of the
        // census it was replacing. The scratch name carries this process's id, because a cohort
        // spread over invocations can have two of them in one directory at once.
        let while_writing = census_path.with_extension(format!(
            "{CENSUS_FILE_EXTENSION}.{}.partial",
            std::process::id()
        ));
        let write = || -> Result<u64, std::io::Error> {
            let mut file = std::fs::File::create(&while_writing)?;
            write_census(&produced.evidence, Some(produced.identity), &mut file)
                .map_err(std::io::Error::other)?;
            file.sync_all()?;
            std::fs::metadata(&while_writing).map(|it| it.len())
        };
        let census_bytes = match write() {
            Ok(bytes) => bytes,
            Err(source) => {
                let _ = std::fs::remove_file(&while_writing);
                return Err(GenerateCensusCliError::OutputDir {
                    path: census_path.clone(),
                    source,
                });
            }
        };
        std::fs::rename(&while_writing, census_path).map_err(|source| {
            GenerateCensusCliError::OutputDir {
                path: census_path.clone(),
                source,
            }
        })?;

        let outcome = SampleCensusOutcome {
            sample: sample.clone(),
            psp: psp.clone(),
            records: produced.identity.records,
            census: census_path.clone(),
            census_bytes,
            tally,
        };
        // **Said as each sample finishes, not only at the end**, and to stderr — so a shell
        // capturing the report gets the report and a person watching gets the progress, in the
        // same words, because both come from `SampleCensusOutcome::line`.
        eprintln!("{}", outcome.line());
        built.push(outcome);
    }

    Ok(CensusReport {
        ground: describe(&analysed, &contigs),
        analysed_bases: analysed.iter().map(|region| region.len()).sum(),
        samples: built,
    })
}

/// **The psps this run reads, with every directory expanded** — one entry a sample, in the order
/// they were given.
///
/// A directory contributes the `.psp` files directly inside it, **sorted by name**, so that two
/// runs naming one directory read the same cohort in the same order however the filesystem
/// answers.
fn psps_named_by(args: &GenerateCensusArgs) -> Result<Vec<PathBuf>, GenerateCensusCliError> {
    let mut paths = Vec::with_capacity(args.psps.len());
    for named in &args.psps {
        if !named.is_dir() {
            paths.push(named.clone());
            continue;
        }
        let mut inside = Vec::new();
        for entry in
            std::fs::read_dir(named).map_err(|source| GenerateCensusCliError::PspDirectory {
                path: named.clone(),
                source,
            })?
        {
            let entry = entry.map_err(|source| GenerateCensusCliError::PspDirectory {
                path: named.clone(),
                source,
            })?;
            let path = entry.path();
            if path.extension().is_some_and(|it| it == PSP_FILE_EXTENSION) && path.is_file() {
                inside.push(path);
            }
        }
        if inside.is_empty() {
            return Err(GenerateCensusCliError::NoPspsInDirectory {
                path: named.clone(),
            });
        }
        inside.sort();
        paths.extend(inside);
    }
    Ok(paths)
}

/// Refuse a sample whose name cannot be the file name of its own census.
///
/// **One normal path component and nothing else**: `@RG SM` is free header text and travels into
/// the psp, so without this a sample could name a path outside `--output-dir` or one that cannot
/// be created at all.
fn refuse_a_sample_name_that_is_not_a_file_name(
    sample: &str,
) -> Result<(), GenerateCensusCliError> {
    let mut components = Path::new(sample).components();
    match (components.next(), components.next()) {
        (Some(std::path::Component::Normal(only)), None) if only == sample => Ok(()),
        _ => Err(GenerateCensusCliError::SampleNameNotAFileName {
            sample: sample.to_string(),
        }),
    }
}
