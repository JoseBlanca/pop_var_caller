//! `regenerate-census` — rebuilding a psp's census out of the records that psp already holds.
//!
//! **A census is what a parameters fit reads.** It holds what one sample showed at a fixed set of
//! positions and repeat tracts chosen for the whole run — a few thousand loci out of a genome — so
//! that the fit can ask the same question of every sample and compare their answers. It lives in
//! the psp's trailer, written there by the walk that wrote the psp
//! (`doc/devel/ng/spec/psp_census_pair.md` §3).
//!
//! **This is the repair, for the psps a fit refuses.** `estimate-parameters` stops and names every
//! sample whose census it cannot use — one written before the census moved into the psp, one
//! written by another build, one recorded under settings this run does not use (spec §4) — and
//! tells the person to run this. It
//! rebuilds each named psp's census from that psp's own records, without opening a single
//! alignment file, and replaces the psp's trailer with it. **Nothing else in the file is
//! rewritten**: the header, the blocks and the index are the bytes they were.
//!
//! **⚠ A repair interrupted mid-write costs a re-walk of that one sample.** The tail is truncated
//! before the new census is written ([`replace_trailer`]), so between those two moments the psp has
//! no footer and no reader accepts it — and a psp that cannot be read cannot have its census
//! rebuilt from its records. The command this replaced could only destroy a cache file beside the
//! psp. Every other failure leaves the file byte for byte what it was, and says so.
//!
//! **Psps in, and nothing else to type.** The ground the walk covered and the criteria it cut that
//! ground with are in every psp's header, and the cohort is refused unless the files agree about
//! them, so there is no flag here that says what a repeat is (spec §6, §8). What is left is the
//! reference and the catalog, which the selection cannot be rebuilt without — and both are checked
//! against the psps before anything is written, because a census recorded against another
//! reference is one every later fit refuses, after this command has spent the rebuild producing
//! it — 0.25 to 3.25 s a sample on the fixtures spec §2 measures, and a quarter of an hour at 50×
//! human.
//!
//! **The selection is the run's, not the sample's**, and its numbers are
//! [`CensusSelection::SHIPPED`]: about two million positions, five thousand tracts a stratum, and
//! a seed that is a compiled-in constant rather than a clock. So two invocations over one cohort
//! keep the same positions, which is what lets a cohort be repaired a sample at a time — on as
//! many machines as there are (spec §8).
//!
//! **One psp open at a time while censuses are rebuilt.** The whole cohort is opened first,
//! because opening it is what checks that the files agree and what says what ground there is — and
//! then closed, before the first census. Holding a thousand psps open to read them one by one
//! would spend the memory psp mode exists to save.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;

use crate::ng::parameter_estimation::joint::census::CensusError;
use crate::ng::parameter_estimation::joint::census_file::write_census;
use crate::ng::parameter_estimation::joint::loci::{SelectionError, UnambiguousRuns};
use crate::ng::psp::{TrailerReplacementFailure, replace_trailer};
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoError, read_reference_observing_or_creating_fai,
};
use crate::ng::repeat_catalog::RepeatCatalogHeader;
use crate::ng::run::report::{describe, plural};
use crate::ng::run::{
    CensusFromPspError, CensusPlan, CensusSelection, CensusTally, OpenPspCohort, RunError,
    Segmentation, THE_COMMAND_THAT_REBUILDS_A_CENSUS, census_from_psp,
};
use crate::pop_var_caller_exp::psp_inputs::{PspArgumentRefusal, psps_named};
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line.
///
/// **Defined from the constant the library names it by, never the reverse.** Two refusals in
/// `ng::run` tell a person to run this command — the report `estimate-parameters` refuses a stale
/// cohort with, and the fit's own backstop — and `ng` imports nothing from this module outside its
/// tests, so pointing that constant at this name would be the first place it did.
pub const SUBCOMMAND: &str = THE_COMMAND_THAT_REBUILDS_A_CENSUS;

/// Rebuild each psp's census from the records it holds, and replace its trailer with it.
#[derive(Debug, Args)]
pub struct RegenerateCensusArgs {
    /// Reference FASTA — the one every psp's samples were aligned to. A `.fai` is built beside
    /// it if there is none.
    ///
    /// **Checked against the psps before anything is written.** The selection of positions is a
    /// function of the reference's bases, and a census rebuilt against another reference is one
    /// every fit refuses.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog the psps were walked against. Defaults to
    /// `<reference>.repeats.parquet`.
    ///
    /// **Checked against the catalog the psps record**, for the same reason as the reference: the
    /// catalog decides which stretches are repeat tracts, so a selection made over another one
    /// would keep tracts these psps hold no observations of.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One psp per sample, or a directory holding them. Repeat the flag.
    ///
    /// A directory contributes every `.psp` file directly inside it, in name order.
    #[arg(long = "psp", required = true, num_args = 1..)]
    pub psps: Vec<PathBuf>,
}

/// Everything that can stop a `regenerate-census` run.
///
/// **Most of them come before any psp is rewritten**, in the order they are met: the reference,
/// the ground it offers to select from, the psp paths, the cohort's own agreement, the reference
/// and the catalog against the psps, and the selection.
///
/// **Three arrive inside the loop, and by then earlier psps have been rewritten**
/// ([`Build`](Self::Build), [`CensusNotEncoded`](Self::CensusNotEncoded),
/// [`TrailerNotReplaced`](Self::TrailerNotReplaced)): a cohort of sixty that fails at the fortieth
/// leaves thirty-nine rebuilt, and the run prints no report — what says which those were is the
/// per-sample progress each one printed as it finished. **The command this replaced had the
/// property this does not**, and could: it judged every output path before doing any work, because
/// what it wrote was a separate file. Here the psp is the output. **Plan step D2 is what makes the
/// re-run cheap** — it skips the psps that need nothing, so a second run does only what is still
/// owed; until then it rebuilds all sixty.
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum RegenerateCensusCliError {
    /// The reference could not be read.
    #[error("reading the reference {}", path.display())]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The ground the censuses are selected over could not be built.
    #[error("working out what ground these psps cover")]
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

    /// `--reference` is not the reference the psps were walked against.
    ///
    /// **Refused before a trailer is replaced.** A census rebuilt against another reference keeps
    /// other positions, records that reference's digest, and is refused by every fit that reads
    /// it — after this run has spent the rebuild.
    #[error(
        "{} is not the reference these psps were walked against; they name {walked_against}, so \
         run this again with that one",
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

    /// The catalog this run reads is not the one the psps were walked with.
    #[error(
        "the repeat catalog {} is not the one these psps were walked with: {difference}; run \
         this again with --catalog naming the one they were",
        path.display()
    )]
    WalkedWithAnotherCatalog {
        /// The catalog this run read — the one `--catalog` named, or the one beside the reference.
        path: PathBuf,
        /// How it differs from the psps' own, as a clause.
        difference: String,
    },

    /// One sample's census could not be rebuilt from its psp.
    #[error("rebuilding {sample}'s census from {}", psp.display())]
    Build {
        /// The individual, as its psp names it.
        sample: String,
        /// The psp being read.
        psp: PathBuf,
        /// What the producer said.
        #[source]
        source: Box<CensusFromPspError>,
    },

    /// A census was rebuilt and could not be encoded.
    #[error("encoding {sample}'s census")]
    CensusNotEncoded {
        /// The individual.
        sample: String,
        /// What the encoder said.
        #[source]
        source: Box<CensusError>,
    },

    /// A psp would not take its new census.
    ///
    /// **What state that file is left in is the cause's to say, and it says it**
    /// ([`crate::ng::psp::FileAfterAFailedReplacement`]): a psp left byte for byte what it was
    /// needs this command
    /// run again, and one already cut back to its trailer has no footer, so that sample has to be
    /// rewritten before anything can read it. **This variant carries no copy of that verdict** —
    /// two fields that could disagree about the state of one file would be a way to tell somebody
    /// to retry a file that is torn.
    #[error("replacing {sample}'s census in {}", psp.display())]
    TrailerNotReplaced {
        /// The individual, as its psp names it.
        sample: String,
        /// The psp that would not take the write.
        psp: PathBuf,
        /// What the writer said, and what it left behind.
        #[source]
        source: Box<TrailerReplacementFailure>,
    },
}

/// What one sample's rebuild cost and holds.
#[derive(Debug, Clone)]
pub struct SampleCensusOutcome {
    /// The individual, as its psp's header names it.
    pub sample: String,
    /// The psp it was rebuilt from and written back into.
    pub psp: PathBuf,
    /// How many records that psp holds.
    pub records: u64,
    /// How large the census written into its trailer is.
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
            "{}: {} stored loci read from {}, census {} bytes written into its trailer",
            self.sample,
            self.records,
            self.psp.display(),
            self.census_bytes,
        );
        if self.tally.contributes_nothing() {
            // **Named rather than omitted.** A census that is all denominator is a legitimate
            // outcome — the walk covered ground the selection kept nothing in, or no read
            // reached what it did keep — and a run that said nothing about it would leave
            // somebody hunting for a psp rewritten exactly as asked.
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
            "regenerated {} census{} over {} — {} bases analysed",
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

/// Rebuild the census in every psp named, and print what each one holds.
///
/// # Errors
///
/// See [`RegenerateCensusCliError`].
pub fn run_regenerate_census(args: &RegenerateCensusArgs) -> Result<(), RegenerateCensusCliError> {
    let report = regenerate_every_census(args)?;
    for line in report.lines() {
        println!("{line}");
    }
    Ok(())
}

/// The run itself, with the printing left to the caller — which is what lets a test read the
/// report rather than parse it back out of a stream.
fn regenerate_every_census(
    args: &RegenerateCensusArgs,
) -> Result<CensusReport, RegenerateCensusCliError> {
    // **The reference is read with an observer**, because the selection has to know where the
    // genome is sequence at all: a position inside a run of `N` has no reference base to compare
    // a read against, and keeping one would put a permanent hole in every sample's records.
    let mut callable = UnambiguousRuns::default();
    let with_checksums = read_reference_observing_or_creating_fai(
        args.reference.clone(),
        ReferenceCheck::VerifyAgainstIndex,
        &mut callable,
    )
    .map_err(|source| RegenerateCensusCliError::Reference {
        path: args.reference.clone(),
        source,
    })?;
    let unambiguous = callable
        .into_selectable()
        .map_err(|source| RegenerateCensusCliError::CensusGround { source })?;
    let contigs = with_checksums.contig_list();

    let paths = psps_named_by(args)?;

    // **The cohort is opened to be agreed with, then closed.** Opening it refuses files walked
    // over different ground, under a different catalog or different criteria, or sharing an
    // `@RG ID`, and it is what says what the analysed regions and the criteria are. Holding it
    // open while censuses are rebuilt would keep every psp of the cohort open to read them one at
    // a time — and each one is about to be rewritten.
    let (analysed, criteria, catalog_of_the_psps, samples) = {
        let cohort =
            OpenPspCohort::open(&paths).map_err(|source| RegenerateCensusCliError::Cohort {
                source: Box::new(source),
            })?;
        // **The reference against the psps, before anything is built from it.** A census rebuilt
        // against another reference keeps other positions and records that reference's digest, so
        // every fit that read it would refuse the cohort — having spent this command's rebuild
        // first. `estimate-parameters` makes the same comparison before it fits (spec §6).
        cohort
            .refuse_a_reference_it_was_not_walked_against(&with_checksums)
            .map_err(
                |source| RegenerateCensusCliError::WalkedAgainstAnotherReference {
                    path: args.reference.clone(),
                    walked_against: cohort
                        .the_reference_the_psps_were_walked_against()
                        .to_string(),
                    source: Box::new(source),
                },
            )?;
        let settings = cohort.segmentation_inputs();
        (
            cohort.analysed_regions().clone(),
            settings.repeat_tract_criteria.clone(),
            settings.catalog.clone(),
            cohort
                .sample_names()
                .map(str::to_string)
                .collect::<Vec<_>>(),
        )
    };

    let catalog_path = run_ground::catalog_path_for(args.catalog.as_deref(), &args.reference);
    let catalog = run_ground::open_catalog(&catalog_path, &args.reference, &with_checksums)?;
    // **The catalog against the psps, before a segment is cut**, for the reason the reference is
    // checked: the census records what the catalog was built at, so one rebuilt against another
    // catalog is refused by every fit that reads it.
    refuse_a_catalog_the_psps_were_not_walked_with(
        &catalog_path,
        catalog.header(),
        &catalog_of_the_psps,
    )?;
    // **The criteria are the psps' own** (spec §6): they decide which stretches are repeat tracts
    // and therefore which loci the selection may keep, and this command has no flag that could
    // say anything else.
    let segmentation = run_ground::segments_cut_from(
        &catalog,
        &catalog_path,
        &criteria,
        run_ground::CriteriaSource::ThePspHeaders,
        &analysed,
    )?;
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
    .map_err(|source| RegenerateCensusCliError::CensusNotPlanned {
        source: Box::new(source),
    })?;
    drop(catalog);

    let mut rebuilt = Vec::with_capacity(paths.len());
    for (psp, sample) in paths.iter().zip(&samples) {
        let outcome = regenerate_one_census(psp, sample, &plan, &segmentation)?;
        // **Said as each sample finishes, not only at the end**, and to stderr — so a shell
        // capturing the report gets the report and a person watching gets the progress, in the
        // same words, because both come from `SampleCensusOutcome::line`.
        eprintln!("{}", outcome.line());
        rebuilt.push(outcome);
    }

    Ok(CensusReport {
        ground: describe(&analysed, &contigs),
        analysed_bases: analysed.iter().map(|region| region.len()).sum(),
        samples: rebuilt,
    })
}

/// **One psp: its records read, its census rebuilt, its trailer replaced** — and nothing else in
/// the file touched.
///
/// **The census is written with no pileup identity.** A census in a file of its own names the psp
/// it was built from, so that the pair can be checked; a census that *is* that psp's trailer has
/// nothing left to pair wrongly with (spec §3.1), and writing one would make this command's output
/// differ from the walk's for a field neither needs.
fn regenerate_one_census(
    psp: &Path,
    sample: &str,
    plan: &CensusPlan,
    segmentation: &Segmentation,
) -> Result<SampleCensusOutcome, RegenerateCensusCliError> {
    let mut produced = census_from_psp(psp, plan, segmentation).map_err(|source| {
        RegenerateCensusCliError::Build {
            sample: sample.to_string(),
            psp: psp.to_path_buf(),
            source: Box::new(source),
        }
    })?;
    let tally = produced
        .tally()
        .map_err(|source| RegenerateCensusCliError::CensusNotEncoded {
            sample: sample.to_string(),
            source: Box::new(source),
        })?;
    let records = produced.identity.records;

    let mut bytes = Vec::new();
    write_census(&produced.evidence, None, &mut bytes).map_err(|source| {
        RegenerateCensusCliError::CensusNotEncoded {
            sample: sample.to_string(),
            source: Box::new(source),
        }
    })?;
    let census_bytes = bytes.len() as u64;
    // **The tail is replaced in place, and the psp's own writer is what decides how.** It keeps
    // everything before the trailer, writes the new footer last, and says which side of that
    // write a failure fell on.
    replace_trailer(psp, &bytes).map_err(|source| {
        RegenerateCensusCliError::TrailerNotReplaced {
            sample: sample.to_string(),
            psp: psp.to_path_buf(),
            source: Box::new(source),
        }
    })?;

    Ok(SampleCensusOutcome {
        sample: sample.to_string(),
        psp: psp.to_path_buf(),
        records,
        census_bytes,
        tally,
    })
}

/// **The psps this run rebuilds, with every directory expanded** — the rule every psp-taking
/// command shares ([`psps_named`]), with its refusals dressed in this command's words.
fn psps_named_by(args: &RegenerateCensusArgs) -> Result<Vec<PathBuf>, RegenerateCensusCliError> {
    psps_named(&args.psps).map_err(|refusal| match refusal {
        PspArgumentRefusal::Unlistable { path, source } => {
            RegenerateCensusCliError::PspDirectory { path, source }
        }
        PspArgumentRefusal::Empty { path } => RegenerateCensusCliError::NoPspsInDirectory { path },
    })
}

/// **Refuse a catalog the psps were not walked with**, saying how it differs from theirs (spec §6).
///
/// **The comparison is the catalog header's own** ([`RepeatCatalogHeader::first_difference`]), so
/// this command and `estimate-parameters` name the same difference in the same words; what is this
/// command's is the sentence around it, which says what to do.
fn refuse_a_catalog_the_psps_were_not_walked_with(
    path: &Path,
    given: &RepeatCatalogHeader,
    walked_with: &RepeatCatalogHeader,
) -> Result<(), RegenerateCensusCliError> {
    match given.first_difference(walked_with) {
        None => Ok(()),
        Some(difference) => Err(RegenerateCensusCliError::WalkedWithAnotherCatalog {
            path: path.to_path_buf(),
            difference,
        }),
    }
}
