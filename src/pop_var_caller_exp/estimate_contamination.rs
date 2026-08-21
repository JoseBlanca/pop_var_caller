//! The `estimate-contamination` subcommand: how much of each sample's DNA came from a
//! different individual, measured over ordinary positions and written as JSON.
//!
//! **This is a side tool and not a step of the caller.** ng estimates contamination as one of
//! the parameters it fits before calling, and uses it internally; this subcommand exists so the
//! same estimator can be handed to somebody comparing methods, who has alignments and wants a
//! number per sample. Nothing here is on the calling path.
//!
//! # What it does, in order
//!
//! 1. Reads the reference and the tandem-repeat catalog built beside it.
//! 2. Decides which stretch of genome to look at — a BED if one is given, otherwise blocks
//!    drawn across the contigs to a total the run asks for.
//! 3. Picks a fixed set of **ordinary positions** in that stretch, the same set in every
//!    sample: repeat tracts are excluded, which is what the catalog is for here.
//! 4. Walks each alignment once, writing down what that sample showed at those positions.
//! 5. Fits the whole cohort at once — error rates, the population's allele frequencies, each
//!    sample's contamination fraction. **No genotype is called anywhere in this**: every
//!    position sums over both unknown genotypes, this sample's and the stray reads', each
//!    weighted by what the allele frequencies make likely and by how well it explains the
//!    reads.
//!
//! # Two things a reader of the output has to know
//!
//! **It needs a panel, and about a dozen samples is the floor.** Each position's allele
//! frequency is fitted from the samples in the run and no outside panel, so that a diverged
//! sample is not scored against the average of everyone else. A sample that supplies more than
//! half of its own fitted frequency is reported as *not identified* rather than given a number,
//! and how much it supplies is roughly `(axes + 1) / samples`. At the default four axes that
//! clears the bar from about twelve samples upwards; below it, expect refusals.
//!
//! **The fraction reads low.** On a drawn panel at 30 reads a position, a sample contaminated
//! at 3% comes back at about 2.6% while every clean sample comes back at 0.0000
//! (`doc/devel/ng/reports/census_depth_resolution_2026-08-16.md`). So the separation between a
//! contaminated sample and a clean one is what this is good at; the value itself is a floor on
//! the truth rather than an estimate of it.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use clap::Args;
use rayon::prelude::*;
use serde::Serialize;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::locus_generation::pileup::{
    PileupGenerator, PileupGeneratorConfig, PileupGeneratorConfigError,
};
use crate::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusGenerationError, SampleLocusObservationsIterator,
    UnhandledReason,
};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::joint::census::{
    CensusWriter, CohortCensusEvidence, DepthCap, ReadCap, SampleCensusEvidence, TermsDisagreement,
};
use crate::ng::parameter_estimation::joint::contamination::{
    ContaminationConfig, ContaminationEstimate, ContaminationSource, DEFAULT_COMPONENTS,
};
use crate::ng::parameter_estimation::joint::fit::{JointFitConfig, JointFitError, fit_jointly};
use crate::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLoci, ReferenceDigest, RegionSetDigest, SelectableRegions,
    SelectionError, SelectionTerms, UnambiguousRuns, select_kept_loci,
};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::SampleReads;
use crate::ng::read::input::read_groups::{ReadGroupError, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::left_align::LeftAlignPreparer;
use crate::ng::ref_seq::WindowedRefSeq;
use crate::ng::reference_info::{
    ReferenceInfo, ReferenceInfoError, ReferenceSource, read_reference_info_observing,
};
use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig};
use crate::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria, sibling_catalog_path,
};
use crate::ng::types::ReadGroupId;
use crate::ng::types::{ContigId, GenomeRegion, Position};
use crate::regions::{BedError, ContigBounds};

/// How much genome the run looks at when no BED is given.
///
/// **Ten megabases, and the binding quantity is segregating positions rather than bases.** The
/// estimate pools information over positions where the cohort varies, and its error against
/// that count flattens out around ten thousand of them; 10 Mb of a plant or vertebrate genome
/// carries several times that. More genome costs a proportional amount of walking for an
/// improvement the estimate stops showing.
const DEFAULT_DRAWN_BP: u64 = 10_000_000;

/// How long each drawn block is.
///
/// **Blocks and not scattered bases**, because reads are fetched by interval and a hundred
/// thousand one-base intervals would be a hundred thousand seeks. A hundred kilobases is also
/// long enough that the reads spanning its two ends are a rounding error against its middle.
const DEFAULT_BLOCK_BP: u64 = 100_000;

/// Contigs shorter than this are not drawn from. Scaffolds carry mismapped reads out of
/// proportion to their length, and a drawn block landing on one buys a position count that the
/// fit then has to judge as noise.
const SMALLEST_CONTIG_DRAWN_FROM: u64 = 1_000_000;

/// How many ordinary positions are kept, when the run does not say.
///
/// Two million is what both of this project's own runs used — the tomato cohort and the human
/// trio — so it is the count the reported behaviour of this estimator was measured at.
const DEFAULT_KEPT_POSITIONS: u64 = 2_000_000;

/// The seed for the block draw and for the choice of kept positions. **Part of the run's
/// identity and reported in the output**: two runs at the same seed over the same reference
/// look at the same bases.
const DEFAULT_SEED: u64 = 42;

/// The depth above which a position's reads stop being counted one by one. Copied from the
/// census walk this subcommand drives; it bounds the per-position list of disagreeing reads.
const DEPTH_CAP: DepthCap = DepthCap::new(124);

/// `estimate-contamination` arguments.
#[derive(Debug, Args, Clone)]
pub struct EstimateContaminationArgs {
    /// Reference FASTA — the one the alignments were made against. Its `.fai` must be beside
    /// it.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog. Defaults to `<reference>.repeats.parquet`, which is where
    /// `repeat-catalog` writes it.
    ///
    /// **Repeat tracts are not measured here; the catalog is what says where they are** so the
    /// kept positions can avoid them. A read crossing a repeat tract can be one repeat unit
    /// short of the DNA it came from, which imitates contamination closely enough that those
    /// positions are left out.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One alignment file per sample (BAM or CRAM, indexed). Repeat the flag.
    ///
    /// **A file holding more than one sample is refused.** The estimate is per sample and the
    /// walk is per sample, so two samples in one file would be pooled into one set of reads.
    #[arg(long = "alignment", required = true, num_args = 1..)]
    pub alignments: Vec<PathBuf>,

    /// Where to write the JSON.
    #[arg(long)]
    pub output: PathBuf,

    /// BED of the stretch of genome to look at. Without it, blocks are drawn to
    /// `--drawn-bp`.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// How much genome to draw when `--regions` is not given.
    #[arg(long, default_value_t = DEFAULT_DRAWN_BP)]
    pub drawn_bp: u64,

    /// How long each drawn block is.
    #[arg(long, default_value_t = DEFAULT_BLOCK_BP, help_heading = "Advanced")]
    pub block_bp: u64,

    /// How many ordinary positions to keep and record in every sample.
    #[arg(long, default_value_t = DEFAULT_KEPT_POSITIONS, help_heading = "Advanced")]
    pub kept_positions: u64,

    /// How many axes of variation each sample's own allele frequency is a straight line in.
    ///
    /// **Lower it on a small panel.** A sample must supply less than half of its own fitted
    /// frequency to be given a number, and the share it supplies is about
    /// `(axes + 1) / samples` — so four axes need about twelve samples, and two need about
    /// seven. Fewer axes model less population structure, and structure that is not modelled
    /// *hides* contamination rather than inventing it.
    #[arg(long, default_value_t = DEFAULT_COMPONENTS)]
    pub components: usize,

    /// The seed for the block draw and the choice of kept positions.
    #[arg(long, default_value_t = DEFAULT_SEED, help_heading = "Advanced")]
    pub seed: u64,

    /// How many threads to use. Zero means every core.
    ///
    /// **It bounds both halves of the run**: how many alignments are walked at once, and how
    /// wide the fit runs. The walks are what it also bounds the *memory* of — each holds its
    /// own view of the reference while it runs.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

/// Everything that can stop a run, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum EstimateContaminationCliError {
    #[error("the reference could not be read")]
    Reference(#[source] ReferenceInfoError),

    #[error("no .fai beside {}; build one with `samtools faidx`", .path.display())]
    MissingFai {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("no repeat catalog at {}; build one with `pop_var_caller_exp repeat-catalog`", .path.display())]
    MissingCatalog { path: PathBuf },

    #[error("the repeat catalog could not be used")]
    Catalog(#[source] RepeatCatalogError),

    #[error("the BED at {} could not be read", .path.display())]
    Bed {
        path: PathBuf,
        #[source]
        source: BedError,
    },

    /// **A stretch of genome that holds nothing to look at is a usage error, not an empty
    /// result.** It happens when a BED names only repeat tracts, or only stretches of `N`.
    #[error(
        "no ordinary positions were kept: the analysed stretch holds {analysed} bases, of \
         which none are outside a repeat tract and unambiguous"
    )]
    NothingToKeep { analysed: u64 },

    #[error("the positions to record could not be chosen")]
    Selection(#[source] SelectionError),

    #[error("the reference has more contigs than a run can address")]
    ContigCount,

    #[error("{} declares {samples} samples; this run is one sample per file", .path.display())]
    NotOneSample { path: PathBuf, samples: usize },

    #[error("the read groups of {} could not be read", .path.display())]
    ReadGroups {
        path: PathBuf,
        #[source]
        source: ReadGroupError,
    },

    #[error("{} could not be opened against this reference", .path.display())]
    OpenAlignment {
        path: PathBuf,
        #[source]
        source: crate::ng::read::input::IngestError,
    },

    #[error("the walk over {} stopped", .path.display())]
    Walk {
        path: PathBuf,
        #[source]
        source: LocusGenerationError,
    },

    #[error("the locus generators could not be built")]
    Generators(#[source] PileupGeneratorConfigError),

    /// The samples did not record the same thing, so they cannot be one cohort. **Naming the
    /// value they differ on is the whole point**: they all fail the same way, silently, and
    /// only the name says what to fix.
    #[error(
        "{} and {} recorded different things and cannot be pooled — they differ on {}",
        .0.first, .0.second, .0.field
    )]
    Cohort(TermsDisagreement),

    #[error("the cohort could not be fitted")]
    Fit(#[source] JointFitError),

    #[error("{} could not be written", .path.display())]
    Output {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------
// What comes out
// ---------------------------------------------------------------------

/// What one read group is called, so a fraction reported against it can be read.
#[derive(Debug, Clone, Serialize)]
pub struct LibraryName {
    /// The index into this sample's read-group table. **Minted per alignment file**, so it is
    /// unique within a sample and not across the run — the sample name and this together are
    /// what identify a library here.
    #[serde(skip)]
    pub read_group: ReadGroupId,
    /// The `@RG ID`, as the file declares it.
    pub id: String,
    /// The `@RG LB`, or what was filled in for it when the file declared none.
    pub library: String,
}

/// **One library's answer** — the fraction of *its* reads that came from another plant.
///
/// **A row per library and not per sample, because a library is what gets contaminated**
/// (2026-08-20). A second plant's DNA enters at library preparation or on the sequencing
/// machine, so two libraries made from one plant can differ, one carrying stray reads and the
/// other none. A sample sequenced from one library has one row and nothing is lost.
#[derive(Debug, Serialize)]
pub struct LibraryContamination {
    pub sample: String,
    /// The `@RG ID` and `@RG LB` this row is about.
    pub read_group: String,
    pub library: String,
    /// The share of this library's reads that came from another individual, or `null` when the
    /// panel cannot answer for it. **Reads low**: see this module's header.
    pub contamination: Option<f64>,
    /// Why there is no number, when there is none.
    pub not_identified: Option<String>,
    /// **Whether the number was fitted from this library's own reads or from every read the
    /// plant produced.** A sample sequenced from one library says `the whole sample's reads`,
    /// because there the two are the same reads and no claim about libraries is being made.
    pub fitted_from: &'static str,
    /// How much of its own fitted allele frequency this library's **sample** supplied. Above 0.5
    /// there is no estimate, because it would be a reading of the sample's own noise — and the
    /// refusal does not carry the value, so this is `null` for a refused sample rather than the
    /// number that refused it.
    pub own_frequency_share: Option<f64>,
    /// **How many of the panel's varying positions this library put a read on, and how many
    /// reads it put there.** This is what tells a library *measured clean* from a library
    /// *nobody could measure*: both come back near zero, and only these say which is which.
    pub markers_with_reads: Option<u64>,
    pub reads_on_markers: Option<u64>,
    /// How many of the kept positions the **sample** had a read at, over all of its libraries.
    pub positions_with_reads: u64,
}

/// What the panel as a whole looked like — **the number a threshold has to clear**, since the
/// floor is a property of this panel's depth and size rather than a constant anyone can quote.
#[derive(Debug, Serialize)]
pub struct PanelSummary {
    pub samples: usize,
    /// How many libraries the panel holds, over all its samples — the number of rows in
    /// `samples`, and equal to `samples` unless somebody was sequenced twice.
    pub libraries: usize,
    pub estimated: usize,
    pub not_identified: usize,
    pub markers: u64,
    pub median: Option<f64>,
    pub highest: Option<f64>,
}

/// What the run was, in enough detail to repeat it.
#[derive(Debug, Serialize)]
pub struct RunIdentity {
    pub version: String,
    pub reference: String,
    pub catalog: String,
    pub regions: String,
    pub analysed_bases: u64,
    pub kept_positions: usize,
    pub components: usize,
    pub seed: u64,
}

/// The whole file.
#[derive(Debug, Serialize)]
pub struct ContaminationReport {
    pub run: RunIdentity,
    pub panel: PanelSummary,
    /// **One row per library**, highest fraction first.
    pub samples: Vec<LibraryContamination>,
}

// ---------------------------------------------------------------------
// The driver
// ---------------------------------------------------------------------

/// Run the estimate and write the JSON.
pub fn run_estimate_contamination(
    args: &EstimateContaminationArgs,
) -> Result<(), EstimateContaminationCliError> {
    if args.threads > 0 {
        // A failure here means a pool is already built, which in this binary means a second
        // call in one process — the estimate is unaffected, so it is not worth an error.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    // Progress goes to stderr, so that redirecting stdout still leaves a person watching a long
    // run something to watch.
    let started = Instant::now();
    let prepared = prepare(args)?;
    eprintln!(
        "analysed {} bases ({}); kept {} ordinary positions, {:.0} s",
        prepared.analysed_bases,
        prepared.regions_described,
        prepared.kept.generic().len(),
        started.elapsed().as_secs_f64()
    );

    // **One sample a thread, and it is this subcommand's own arrangement rather than the
    // caller's.** How ng parallelises a real run is being settled elsewhere; what makes it
    // easy here is that a sample's walk shares nothing with any other sample's — each builds
    // its own view of the reference and its own generators inside its own thread, and hands
    // back only the evidence it wrote down. Those views are file-backed and single-threaded,
    // so they may not cross a thread boundary; they never do, because nothing outside the
    // walk ever holds one.
    //
    // **The cost is one open view of the reference a thread.** At a hundred samples on a
    // machine with many cores that is the thing to watch, and `--threads` is what bounds it.
    let done = AtomicUsize::new(0);
    let total = args.alignments.len();
    let walked: Vec<(SampleCensusEvidence, Vec<LibraryName>)> = args
        .alignments
        .par_iter()
        .map(|alignment| {
            let at = Instant::now();
            let (sample, names) = walk_one_sample(alignment, &prepared)?;
            // The order lines appear in is the order walks *finish*, which is not the order
            // the samples were given — so each line names its sample, and the count is how
            // many are done rather than which one this is.
            eprintln!(
                "walked {} ({} of {total} done), {:.0} s",
                sample.sample,
                done.fetch_add(1, Ordering::Relaxed) + 1,
                at.elapsed().as_secs_f64()
            );
            Ok((sample, names))
        })
        // **Collected in the order the alignments were given**, whatever order they finished
        // in: a cohort assembled in a racing order would fit the same numbers but report them
        // against a different arrangement of samples from one run to the next.
        .collect::<Result<Vec<_>, EstimateContaminationCliError>>()?;
    let library_names: BTreeMap<String, Vec<LibraryName>> = walked
        .iter()
        .map(|(sample, names)| (sample.sample.clone(), names.clone()))
        .collect();
    let cohort: Vec<SampleCensusEvidence> = walked.into_iter().map(|(sample, _)| sample).collect();
    eprintln!(
        "walked {total} samples in {:.0} s",
        started.elapsed().as_secs_f64()
    );
    let mut cohort =
        CohortCensusEvidence::new(cohort).map_err(EstimateContaminationCliError::Cohort)?;

    let config = JointFitConfig {
        edges: Arc::new(DepthBinEdges::for_census()),
        contamination: ContaminationConfig {
            components: args.components,
            ..ContaminationConfig::default()
        },
        ..JointFitConfig::default()
    };
    let at = Instant::now();
    let fit = fit_jointly(&mut cohort, &config).map_err(EstimateContaminationCliError::Fit)?;
    eprintln!(
        "fitted {} samples in {} passes, {:.0} s",
        args.alignments.len(),
        fit.passes,
        at.elapsed().as_secs_f64()
    );

    let report = assemble(args, &prepared, &library_names, &fit);
    report_to_the_terminal(&report);
    let json = serde_json::to_string_pretty(&report).expect("a report of numbers and names");
    std::fs::write(&args.output, json + "\n").map_err(|source| {
        EstimateContaminationCliError::Output {
            path: args.output.clone(),
            source,
        }
    })
}

/// What a person watching the run should be told before they open the file.
///
/// **The panel's own spread and not each sample's value**, because there is no threshold to
/// quote: the floor a contaminated sample has to stand above is a property of this panel's
/// depth, size and diversity. What is worth saying out loud is where the panel sits and how
/// much evidence it rests on.
fn report_to_the_terminal(report: &ContaminationReport) {
    let panel = &report.panel;
    eprintln!(
        "\n{} of {} libraries from {} samples estimated over {} varying positions, {} not \
         identified",
        panel.estimated, panel.libraries, panel.samples, panel.markers, panel.not_identified
    );
    match (panel.median, panel.highest) {
        (Some(median), Some(highest)) => eprintln!(
            "the panel's own spread: median {median:.4}, highest {highest:.4} ({})",
            report.samples.first().map_or_else(
                || "no library".to_string(),
                |row| format!("{}, library {}", row.sample, row.library)
            )
        ),
        _ => eprintln!("no sample could be estimated; the panel is too small or too uniform"),
    }
    // **The count of varying positions is what the fraction's precision rests on**, and below
    // about a thousand the published error of this estimator at 1% contamination is of the same
    // size as the thing being estimated. Said here because a reader of the file cannot be
    // expected to know where that line falls.
    if panel.estimated > 0 && panel.markers < THIN_MARKERS {
        eprintln!(
            "\nWARNING: {} varying positions is thin. Below about {THIN_MARKERS} the fraction is \
             dominated by noise — look at more genome (--regions or --drawn-bp), or expect these \
             numbers to move if you run it again over a different stretch.",
            panel.markers
        );
    }
    if panel.samples < SMALL_PANEL {
        eprintln!(
            "\nWARNING: {} samples is a small panel. Each position's allele frequency is fitted \
             from the samples in this run, so a sample supplying more than half of its own \
             frequency is refused rather than given a number — {} were. Lower --components, or \
             add samples.",
            panel.samples, panel.not_identified
        );
    }
}

/// Below this many varying positions, the fraction is noise. From the estimator's own error
/// against marker count: at a thousand markers its error at 1% contamination is 0.69 of the
/// value, and at ten thousand it is 0.11.
const THIN_MARKERS: u64 = 1_000;

/// Below this many samples, expect refusals at the default four axes: a sample supplies about
/// `(axes + 1) / samples` of its own fitted frequency and is refused above a half.
const SMALL_PANEL: usize = 12;

/// Everything the walks share: which positions to record, which regions to walk, and what to
/// record them under.
struct Prepared {
    fasta: PathBuf,
    info: Arc<ReferenceInfo>,
    contigs: Arc<ContigList>,
    index: Arc<noodles_fasta::fai::Index>,
    kept: CensusLoci,
    terms: SelectionTerms,
    /// Only the ordinary stretches. **Repeat tracts are not walked at all**, so no read is
    /// realigned against one — the catalog's whole job here is to say which stretches to skip.
    generic_regions: Vec<TypedRegion>,
    analysed_bases: u64,
    regions_described: String,
    catalog_path: PathBuf,
}

fn prepare(args: &EstimateContaminationArgs) -> Result<Prepared, EstimateContaminationCliError> {
    // ---- the reference, and where it is sequence at all ---------------------------------
    let mut callable = UnambiguousRuns::default();
    let info = Arc::new(
        read_reference_info_observing(
            ReferenceSource::Fasta {
                fasta: args.reference.clone(),
                fai: None,
            },
            &mut callable,
        )
        .map_err(EstimateContaminationCliError::Reference)?,
    );
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&args.reference).map_err(|source| {
        EstimateContaminationCliError::MissingFai {
            path: args.reference.clone(),
            source,
        }
    })?;
    let unambiguous = callable
        .into_selectable()
        .map_err(EstimateContaminationCliError::Selection)?;

    // ---- what stretch of genome to look at ------------------------------------------------
    let bounds: Vec<ContigBounds> = contigs
        .entries
        .iter()
        .map(|entry| {
            u32::try_from(entry.length)
                .map(|length| ContigBounds {
                    name: &entry.name,
                    length,
                })
                .map_err(|_| EstimateContaminationCliError::ContigCount)
        })
        .collect::<Result<_, _>>()?;
    let (analysed, regions_described) = match &args.regions {
        Some(bed) => {
            let spans = GenomeRegions::from_bed_path(bed, &bounds)
                .map_err(|source| EstimateContaminationCliError::Bed {
                    path: bed.clone(),
                    source,
                })?
                .iter()
                .collect();
            (
                SelectableRegions::new(spans).map_err(EstimateContaminationCliError::Selection)?,
                bed.display().to_string(),
            )
        }
        None => {
            let drawn = draw_blocks(&contigs, args.drawn_bp, args.block_bp, args.seed);
            let described = format!(
                "{} blocks of {} bp drawn at seed {}",
                drawn.len(),
                args.block_bp,
                args.seed
            );
            (
                SelectableRegions::new(drawn).map_err(EstimateContaminationCliError::Selection)?,
                described,
            )
        }
    };

    // ---- the catalog, and the positions to record -----------------------------------------
    let catalog_path = args
        .catalog
        .clone()
        .unwrap_or_else(|| sibling_catalog_path(&args.reference));
    if !catalog_path.exists() {
        return Err(EstimateContaminationCliError::MissingCatalog { path: catalog_path });
    }
    let catalog = RepeatCatalog::open_checking_against_reference(&catalog_path, &info)
        .map_err(EstimateContaminationCliError::Catalog)?;
    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());

    let terms = SelectionTerms {
        seed: args.seed,
        reference: ReferenceDigest::of(&info).map_err(EstimateContaminationCliError::Selection)?,
        analysed_regions: RegionSetDigest::of(&analysed),
        catalog_built_under: CatalogBuildSettings::of(&catalog),
        ssr_criteria: criteria.clone(),
        generic_target: args.kept_positions,
        // **No repeat tracts at all.** Contamination is fitted on ordinary positions; a cap of
        // zero keeps none, so nothing about a tract is recorded, walked or fitted.
        ssr_cap: 0,
    };
    let kept =
        select_kept_loci(&terms, &catalog, &analysed, &unambiguous).map_err(|e| match e {
            SelectionError::Catalog(inner) => EstimateContaminationCliError::Catalog(inner),
            other => EstimateContaminationCliError::Selection(other),
        })?;
    if kept.generic().is_empty() {
        return Err(EstimateContaminationCliError::NothingToKeep {
            analysed: analysed.total_length().get(),
        });
    }

    // The walk and the selection run over the same stretch: the analysed regions with the
    // ambiguous runs cut out, and then only the parts that are not repeat tracts.
    let masked = analysed.intersect(&unambiguous);
    let generic_regions: Vec<TypedRegion> = catalog
        .genome_segments(&criteria, ReadScope::Regions(masked.spans()))
        .map_err(EstimateContaminationCliError::Catalog)?
        .map(|item| item.map_err(EstimateContaminationCliError::Catalog))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|region| region.kind == RegionKind::Generic)
        .collect();

    Ok(Prepared {
        fasta: args.reference.clone(),
        info,
        contigs,
        index,
        kept,
        terms,
        generic_regions,
        analysed_bases: analysed.total_length().get(),
        regions_described,
        catalog_path,
    })
}

/// One sample: one walk over the ordinary stretches, filling that sample's census.
fn walk_one_sample(
    alignment: &Path,
    prepared: &Prepared,
) -> Result<(SampleCensusEvidence, Vec<LibraryName>), EstimateContaminationCliError> {
    let read_groups =
        build_read_groups(std::slice::from_ref(&alignment.to_path_buf())).map_err(|source| {
            EstimateContaminationCliError::ReadGroups {
                path: alignment.to_path_buf(),
                source,
            }
        })?;
    let sample = match read_groups.read_groups_per_sample() {
        [only] => only.clone(),
        other => {
            return Err(EstimateContaminationCliError::NotOneSample {
                path: alignment.to_path_buf(),
                samples: other.len(),
            });
        }
    };
    let reference = OpenReference::new(Arc::clone(&prepared.info));
    let sample_reads = SampleReads::open(
        &sample,
        &read_groups,
        &reference,
        ReadFilterConfig::default(),
        true,
    )
    .map_err(|source| EstimateContaminationCliError::OpenAlignment {
        path: alignment.to_path_buf(),
        source,
    })?;

    // Owned captures: the generator outlives these borrows and keeps the maker so it can mint a
    // fresh view of the reference whenever it evicts one.
    let accessor = {
        let fasta = prepared.fasta.clone();
        let contigs = Arc::clone(&prepared.contigs);
        let index = Arc::clone(&prepared.index);
        move || {
            WindowedRefSeq::with_shared_index(
                fasta.clone(),
                Arc::clone(&contigs),
                Arc::clone(&index),
            )
        }
    };
    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "file-backed and single-threaded, as in the census walk this is taken from"
    )]
    let shared = Arc::new(accessor());
    let generic_generator = PileupGenerator::new(
        Arc::clone(&shared),
        accessor.clone(),
        LeftAlignPreparer::with_default_normalizer(accessor()),
        PileupGeneratorConfig::default(),
    )
    .map_err(EstimateContaminationCliError::Generators)?;
    // **No repeat-tract generator.** No tract is kept, so nothing would read what one produced,
    // and building one would realign every read that crosses a tract for nothing.
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::OutOfScope),
        GeneratorSlot::Generator(Box::new(generic_generator)),
        GeneratorSlot::Unfilled(UnhandledReason::OutOfScope),
    );

    let contig_of = |name: &str| {
        prepared
            .contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|i| ContigId(i as u32))
    };
    let mut writer = CensusWriter::new(
        sample.sample.to_string(),
        &prepared.kept,
        sample.read_groups.clone(),
        &contig_of,
        prepared.terms.clone(),
        DepthBinEdges::for_census(),
        ReadCap(crate::ng::locus_generation::ssr::DEFAULT_SSR_MAX_READS_PER_LOCUS),
        DEPTH_CAP,
    );

    // **Every stretch handed to the walk is marked walked, whatever it holds.** Without this a
    // position no read reached is indistinguishable from a stretch the run never opened, since
    // the generator emits no locus where there is no read.
    for region in &prepared.generic_regions {
        writer.mark_walked(region.region);
    }

    let regions: Vec<Result<TypedRegion, RepeatCatalogError>> =
        prepared.generic_regions.iter().cloned().map(Ok).collect();
    let mut stream =
        SampleLocusObservationsIterator::new(regions.into_iter(), sample_reads, generators);
    for locus in &mut stream {
        let locus = locus.map_err(|source| EstimateContaminationCliError::Walk {
            path: alignment.to_path_buf(),
            source,
        })?;
        writer.add_locus(&locus);
    }
    // **What each read group is called, carried out beside the evidence.** The fraction is now
    // reported per read group, and an index into a table the reader does not have is not an
    // answer — so the library and platform names travel with it.
    let names = sample
        .read_groups
        .iter()
        .map(|id| {
            let group = read_groups.get(*id);
            LibraryName {
                read_group: *id,
                id: group.id.to_string(),
                library: group.library.value.to_string(),
            }
        })
        .collect();
    Ok((writer.finish(), names))
}

/// Turn the fit into the file.
fn assemble(
    args: &EstimateContaminationArgs,
    prepared: &Prepared,
    library_names: &BTreeMap<String, Vec<LibraryName>>,
    fit: &crate::ng::parameter_estimation::joint::fit::JointFit,
) -> ContaminationReport {
    let mut samples: Vec<LibraryContamination> = Vec::new();
    let mut estimated: Vec<f64> = Vec::new();
    let mut markers_seen = 0_u64;

    for (name, per_read_group) in &fit.contamination {
        let positions_with_reads = fit
            .rates
            .get(name)
            .map_or(0, |rates| rates.value.positions_with_reads);
        for (group, estimate) in per_read_group {
            // The walk recorded what each read group is called; an id with no name behind it
            // means the two lists disagree, and saying so beats printing a bare number.
            let named = library_names
                .get(name)
                .and_then(|names| names.iter().find(|entry| entry.read_group == *group));
            let (id, library) = named.map_or_else(
                || (format!("read group {group}"), "unnamed".to_string()),
                |entry| (entry.id.clone(), entry.library.clone()),
            );
            let row = match estimate {
                ContaminationEstimate::Estimated {
                    alpha,
                    source,
                    panel_markers,
                    markers_with_reads,
                    reads_on_markers,
                    leverage,
                } => {
                    estimated.push(*alpha);
                    markers_seen = *panel_markers;
                    LibraryContamination {
                        sample: name.clone(),
                        read_group: id,
                        library,
                        contamination: Some(*alpha),
                        not_identified: None,
                        fitted_from: match source {
                            ContaminationSource::ThisReadGroupsReads => "this library's reads",
                            ContaminationSource::TheWholeSamplesReads => "the whole sample's reads",
                        },
                        own_frequency_share: Some(*leverage),
                        markers_with_reads: Some(*markers_with_reads),
                        reads_on_markers: Some(*reads_on_markers),
                        positions_with_reads,
                    }
                }
                ContaminationEstimate::NotIdentified { reason } => LibraryContamination {
                    sample: name.clone(),
                    read_group: id,
                    library,
                    contamination: None,
                    not_identified: Some(reason.to_string()),
                    fitted_from: "nothing — no estimate was made",
                    own_frequency_share: None,
                    markers_with_reads: None,
                    reads_on_markers: None,
                    positions_with_reads,
                },
            };
            samples.push(row);
        }
    }

    estimated.sort_by(f64::total_cmp);
    let libraries = samples.len();
    let panel = PanelSummary {
        samples: fit.contamination.len(),
        libraries,
        estimated: estimated.len(),
        not_identified: libraries - estimated.len(),
        markers: markers_seen,
        median: estimated.get(estimated.len() / 2).copied(),
        highest: estimated.last().copied(),
    };

    // Highest first: what a reader of this file is looking for is which sample stands out.
    samples.sort_by(|a, b| match (a.contamination, b.contamination) {
        (Some(x), Some(y)) => y.total_cmp(&x),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => (&a.sample, &a.read_group).cmp(&(&b.sample, &b.read_group)),
    });

    ContaminationReport {
        run: RunIdentity {
            version: env!("CARGO_PKG_VERSION").to_string(),
            reference: args.reference.display().to_string(),
            catalog: prepared.catalog_path.display().to_string(),
            regions: prepared.regions_described.clone(),
            analysed_bases: prepared.analysed_bases,
            kept_positions: prepared.kept.generic().len(),
            components: args.components,
            seed: args.seed,
        },
        panel,
        samples,
    }
}

// ---------------------------------------------------------------------
// Drawing the stretch to look at
// ---------------------------------------------------------------------

/// Blocks spread across the contigs, to a total of `total_bp`.
///
/// **Each contig gets blocks in proportion to its length**, and within a contig they are spaced
/// evenly from an offset the seed chooses. Evenly rather than at random because two blocks
/// drawn independently can overlap or abut, and a run that asks for 10 Mb should look at 10 Mb;
/// the seed still moves every block, so two seeds see different bases.
///
/// The result is sorted, disjoint, and inside the contigs.
fn draw_blocks(contigs: &ContigList, total_bp: u64, block_bp: u64, seed: u64) -> Vec<GenomeRegion> {
    let block_bp = block_bp.max(1);
    let eligible: Vec<(usize, u64)> = contigs
        .entries
        .iter()
        .enumerate()
        .filter(|(_, entry)| entry.length >= SMALLEST_CONTIG_DRAWN_FROM.max(block_bp))
        .map(|(index, entry)| (index, entry.length))
        .collect();
    let drawable: u64 = eligible.iter().map(|(_, length)| *length).sum();
    if drawable == 0 {
        return Vec::new();
    }
    let wanted = total_bp.min(drawable);

    let mut out = Vec::new();
    for (index, length) in eligible {
        // This contig's share of the total, rounded to whole blocks and at least one wherever
        // the share reaches half a block.
        let share = (wanted as u128 * length as u128 / drawable as u128) as u64;
        let blocks = ((share + block_bp / 2) / block_bp).min(length / block_bp);
        if blocks == 0 {
            continue;
        }
        let stride = length / blocks;
        let offset = hash(seed, index as u64) % stride.max(1);
        for block in 0..blocks {
            let start = block * stride + offset;
            let end = (start + block_bp - 1).min(length - 1);
            if start > end {
                continue;
            }
            out.push(GenomeRegion {
                contig: ContigId(index as u32),
                start: Position(start),
                end: Position(end),
            });
        }
    }
    out.sort_by_key(|region| (region.contig.get(), region.start.get()));
    out
}

/// A seeded value for one contig. Splitmix64, so two adjacent contig indices do not give two
/// adjacent offsets.
fn hash(seed: u64, index: u64) -> u64 {
    let mut z = seed
        .wrapping_add(index.wrapping_mul(0x9E37_79B9_7F4A_7C15))
        .wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::ContigEntry;

    fn contigs(lengths: &[u64]) -> ContigList {
        ContigList {
            entries: lengths
                .iter()
                .enumerate()
                .map(|(index, length)| ContigEntry {
                    name: format!("chr{index}"),
                    length: *length,
                    md5: None,
                })
                .collect(),
        }
    }

    #[test]
    fn the_draw_reaches_what_was_asked_for_and_stays_inside_the_contigs() {
        let contigs = contigs(&[50_000_000, 30_000_000, 20_000_000]);
        let blocks = draw_blocks(&contigs, 10_000_000, 100_000, 42);
        let total: u64 = blocks.iter().map(|region| region.len()).sum();
        // Whole blocks per contig, so the total lands within one block of what was asked.
        assert!(
            (9_900_000..=10_100_000).contains(&total),
            "drew {total} bases of a wanted 10,000,000"
        );
        for region in &blocks {
            let length = contigs.entries[region.contig.get() as usize].length;
            assert!(region.end.get() < length, "{region:?} runs off its contig");
            assert!(region.start.get() <= region.end.get());
        }
    }

    #[test]
    fn the_drawn_blocks_are_disjoint_and_sorted() {
        let blocks = draw_blocks(&contigs(&[50_000_000, 30_000_000]), 10_000_000, 100_000, 7);
        for pair in blocks.windows(2) {
            let (first, second) = (pair[0], pair[1]);
            assert!(
                (first.contig.get(), first.start.get()) < (second.contig.get(), second.start.get()),
                "{first:?} and {second:?} are out of order"
            );
            if first.contig == second.contig {
                assert!(
                    first.end.get() < second.start.get(),
                    "{first:?} and {second:?} overlap, so a position would be kept twice"
                );
            }
        }
    }

    /// **Two seeds must look at different bases**, or the seed is not doing anything and two
    /// runs asked to disagree would agree.
    #[test]
    fn a_different_seed_draws_different_blocks() {
        let contigs = contigs(&[50_000_000, 30_000_000]);
        let first = draw_blocks(&contigs, 10_000_000, 100_000, 1);
        let second = draw_blocks(&contigs, 10_000_000, 100_000, 2);
        assert_eq!(first.len(), second.len(), "the same count either way");
        assert_ne!(first, second, "the seed moved no block");
    }

    /// A contig shorter than the floor contributes nothing: a scaffold's mismapped reads are
    /// not what this estimate should be made of.
    #[test]
    fn short_contigs_are_left_out() {
        let contigs = contigs(&[50_000_000, 900_000]);
        let blocks = draw_blocks(&contigs, 10_000_000, 100_000, 42);
        assert!(
            blocks.iter().all(|region| region.contig.get() == 0),
            "a contig below the floor was drawn from"
        );
    }
}
