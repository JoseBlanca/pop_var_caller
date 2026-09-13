//! **What each route to a census costs** — plan step B3 of
//! `doc/devel/ng/impl_plan/parameter_prepass_runs.md`.
//!
//! There are two ways to end up with every sample's census, and they are not two ways of doing
//! the same job:
//!
//! - **during the walk** — one pass over the alignment files builds the census as it goes and
//!   seals it into the psp's trailer. This is what every run does: `generate-psps` always builds
//!   a census and there is no flag to skip it (`psp_census_pair.md` §3.3).
//! - **afterwards** — the psp is walked with an empty trailer, and a second pass reads the stored
//!   records back and builds the census from them. **No command produces a psp in that state**;
//!   it is the state a psp written by an older build arrives in, and the second pass is the one
//!   `regenerate-census` makes to repair it (§4.2, §8). This harness produces it on purpose, to
//!   price that repair.
//!
//! **The two produce the same census, byte for byte**, and neither carries a pileup identity —
//! the psp-header digest and record count by which a census kept in a file of its own names the
//! file it came from. A census that *is* its psp's trailer has nothing left to pair wrongly with
//! (`psp_census_pair.md` §3), so no shipped route writes one and this harness writes one on
//! neither side. What separates the routes is therefore only what they cost. The module
//! `run::census_from_psp::the_two_producers_agree` holds them to that agreement on fixtures;
//! this harness is where it is checked on real reads, which is why **both routes leave a
//! `<sample>.census` file for the wrapping script to compare** even though neither writes such a
//! file in a real run — `regenerate-census` puts its census in the psp's trailer too. On the
//! during-the-walk route that file is the trailer copied out, **after the clock stops**, so the
//! route it would slow is not charged for it.
//!
//! This harness runs **one** route per process, so that a wrapper measuring peak resident memory
//! measures one route rather than the larger of two.
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example ng_census_route_cost -- \
//!     <during-the-walk|after-the-walk> <reference.fa> <catalog.parquet> <regions.bed> <cram-or-dir>
//! ```
//!
//! `NG_SAMPLES=n` walks the first `n` alignment files of a directory (default 1); `NG_REGIONS=n`
//! analyses the first `n` BED intervals (default 2); `NG_WORK=dir` says where the files land
//! (default `tmp/ng_census_route_cost/<route>`, inside the project).
//!
//! # The ground is the one a real run walks
//!
//! The segmentation comes from `run_ground::segments_over`, the same call `generate-psps` and
//! `call-from-psps` make. **That is not a detail**: asking the catalog directly with
//! `StrRepeatCriteria::default()` cuts the genome at the storage floors and produces several
//! times as many segments over the same ground, and a cost measured over that segmentation is a
//! cost no run pays.
//!
//! # What it prints
//!
//! One `route=... samples=... psp_bytes=... census_bytes=... seconds=...` line, and per sample
//! the two file sizes. Wall time here is the harness's own clock over the work; the peak
//! resident memory is the wrapper's, since a process cannot see its own high-water mark
//! portably.
//!
//! **`census_bytes` is the census exactly as a real run stores it** — what the walk seals into
//! the psp's trailer, and what `regenerate-census` writes back there. Both routes also leave a
//! `<sample>.census` in the work directory for the wrapping script's `cmp`; neither route writes
//! such a file in a real run.

use std::error::Error;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::cli::run_ground::{self, GroundRequest, RepeatRouting};
use pop_var_caller::locus_generation::pileup::PileupGeneratorConfig;
use pop_var_caller::parameter_estimation::joint::census_file::write_census;
use pop_var_caller::parameter_estimation::joint::loci::UnambiguousRuns;
use pop_var_caller::psp::{ParameterValue, PspReader, WriterProvenance};
use pop_var_caller::read::ReadFilterConfig;
use pop_var_caller::read::input::reference::OpenReference;
use pop_var_caller::reference_info::{ReferenceCheck, read_reference_observing_or_creating_fai};
use pop_var_caller::region_typing::DEFAULT_MAX_STR_LEN;
use pop_var_caller::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use pop_var_caller::repeat_catalog::RepeatCatalog;
use pop_var_caller::run::{
    CensusPlan, CensusSelection, SampleObservationGatherer, SampleWalkInputs, census_from_psp,
};

const SAMPLES_BY_DEFAULT: usize = 1;
const REGIONS_BY_DEFAULT: usize = 2;

/// Which route this process runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Route {
    /// One pass over the reads writes the psp and the census together.
    DuringTheWalk,
    /// The walk writes the psp alone; a second pass reads it back and builds the census.
    AfterTheWalk,
}

impl Route {
    fn of(word: &str) -> Option<Self> {
        match word {
            "during-the-walk" => Some(Self::DuringTheWalk),
            "after-the-walk" => Some(Self::AfterTheWalk),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::DuringTheWalk => "during-the-walk",
            Self::AfterTheWalk => "after-the-walk",
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [route, fasta, catalog, bed, crams] = args.as_slice() else {
        eprintln!(
            "usage: ng_census_route_cost <during-the-walk|after-the-walk> <reference.fa> \
             <catalog.parquet> <regions.bed> <cram-or-dir>\n\
             runs one route to a census and reports what it cost.\n\
             NG_SAMPLES=n takes the first n files of a directory (default \
             {SAMPLES_BY_DEFAULT}); NG_REGIONS=n the first n BED intervals (default \
             {REGIONS_BY_DEFAULT}); NG_WORK=dir is where the files land."
        );
        return ExitCode::from(2);
    };
    let Some(route) = Route::of(route) else {
        eprintln!("the route must be during-the-walk or after-the-walk, and it is '{route}'");
        return ExitCode::from(2);
    };
    match run(
        route,
        Path::new(fasta),
        Path::new(catalog),
        Path::new(bed),
        Path::new(crams),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            let mut cause = error.source();
            while let Some(next) = cause {
                eprintln!("  caused by: {next}");
                cause = next.source();
            }
            ExitCode::FAILURE
        }
    }
}

/// How many of something this run takes: the environment's answer, or `fallback`.
fn how_many(name: &str, fallback: usize) -> Result<usize, Box<dyn Error>> {
    let Ok(value) = std::env::var(name) else {
        return Ok(fallback);
    };
    let count: usize = value
        .parse()
        .map_err(|_| format!("{name} must be a count, and it is '{value}'"))?;
    if count == 0 {
        return Err(format!("{name} must be at least 1, and it is 0").into());
    }
    Ok(count)
}

fn run(
    route: Route,
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    crams: &Path,
) -> Result<(), Box<dyn Error>> {
    let setting_up = Instant::now();

    // **Read with an observer**, because the census selection has to know where the reference is
    // sequence at all: a position inside a run of `N` has no base to compare a read against.
    let mut callable = UnambiguousRuns::default();
    let with_checksums = Arc::new(read_reference_observing_or_creating_fai(
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
        &mut callable,
    )?);
    let unambiguous = callable.into_selectable()?;
    let contigs = with_checksums.contig_list();
    let reference = OpenReference::new(Arc::clone(&with_checksums));

    let scratch = tempfile::tempdir()?;
    let trimmed = first_regions_of(
        bed,
        how_many("NG_REGIONS", REGIONS_BY_DEFAULT)?,
        scratch.path(),
    )?;
    let ground = GroundRequest {
        reference: fasta,
        catalog: Some(catalog_path),
        regions: Some(&trimmed),
        routing: RepeatRouting {
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        },
    };
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();
    let segmentation = Arc::new(run_ground::segments_over(
        &ground,
        &analysed,
        &with_checksums,
    )?);

    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let plan = CensusPlan::of_run(
        CensusSelection::SHIPPED,
        &catalog,
        &analysed,
        &unambiguous,
        &with_checksums,
        &segmentation.inputs().repeat_tract_criteria,
    )?;
    drop(catalog);

    let mut paths: Vec<PathBuf> = if crams.is_dir() {
        std::fs::read_dir(crams)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .filter(|path| {
                path.extension()
                    .is_some_and(|kind| kind == "cram" || kind == "bam")
            })
            .collect()
    } else {
        vec![crams.to_path_buf()]
    };
    paths.sort();
    paths.truncate(how_many("NG_SAMPLES", SAMPLES_BY_DEFAULT)?);
    if paths.is_empty() {
        return Err(format!("no .cram or .bam under {}", crams.display()).into());
    }

    let work_dir = PathBuf::from(
        std::env::var("NG_WORK")
            .unwrap_or_else(|_| format!("tmp/ng_census_route_cost/{}", route.name())),
    );
    std::fs::create_dir_all(&work_dir)?;

    println!("# route: {}", route.name());
    println!("# analysed bases: {analysed_bases}");
    println!("# segments: {}", segmentation.segments().len());
    println!("# census positions: {}", plan.loci.generic().len());
    println!("# samples: {}", paths.len());
    println!("# setup: {:.2} s", setting_up.elapsed().as_secs_f64());

    // **The clock starts after the setup**, because the setup is the same work whichever route
    // runs and including it would dilute the difference the measurement is about.
    let working = Instant::now();
    let mut psp_bytes = 0_u64;
    let mut census_bytes = 0_u64;
    // The psps whose trailer is copied out for the wrapping script, once the clock has stopped.
    let mut to_copy_out: Vec<(PathBuf, PathBuf)> = Vec::new();
    for path in &paths {
        let alignments = [path.clone()];
        let gatherer = SampleObservationGatherer::open(
            SampleWalkInputs {
                alignments: &alignments,
                reference: &reference,
                read_filters: ReadFilterConfig::default(),
                locus_generator_settings: PileupGeneratorConfig::default(),
                build_index_if_missing: false,
            },
            Arc::clone(&segmentation),
            provenance(),
            // **The one difference between the two routes at the walk**: the route that builds
            // its census afterwards walks without a plan, so nothing is accumulated as the
            // loci go past.
            match route {
                Route::DuringTheWalk => Some(&plan),
                Route::AfterTheWalk => None,
            },
        )?;
        let sample = gatherer.sample_name().to_string();
        let psp_path = work_dir.join(format!("{sample}.psp"));
        let census_path = work_dir.join(format!("{sample}.census"));

        // **The census goes into the psp's trailer on the first route and is built from the
        // stored psp on the second**, which is what the two routes are: one walk that accumulates
        // it, and one that reads the psp back afterwards. **This is all of the work either route
        // does, and it is what the clock is around.**
        let (stats, census_size) = match route {
            Route::DuringTheWalk => {
                let (stats, _) = gatherer.write_psp(&psp_path)?;
                // The footer, not the payload: how big the census is, without reading it.
                let inside = PspReader::open(&psp_path)?.footer().trailer_bytes;
                to_copy_out.push((psp_path.clone(), census_path.clone()));
                (stats, inside)
            }
            Route::AfterTheWalk => {
                let (stats, _) = gatherer.write_psp(&psp_path)?;
                let produced = census_from_psp(&psp_path, &plan, &segmentation)?;
                let mut file = std::fs::File::create(&census_path)?;
                write_census(&produced.evidence, &mut file)?;
                (stats, std::fs::metadata(&census_path)?.len())
            }
        };
        println!(
            "{sample}: {} records, psp {} bytes, census {census_size} bytes",
            stats.records, stats.bytes,
        );
        psp_bytes += stats.bytes;
        census_bytes += census_size;
    }
    let seconds = working.elapsed().as_secs_f64();

    // **After the clock, on purpose.** The wrapping script compares the two routes' censuses, so
    // both have to leave a `<sample>.census` behind — but on this route the census is already in
    // the psp, and reading the trailer back and writing it out again is a quarter of a megabyte a
    // sample each way, pure double handling. Charged to the route, it would bias the one
    // measurement this harness exists to make, against the route that is cheaper.
    //
    // **The other route's write stays inside the clock**, because its real counterpart makes one
    // of the same size: `regenerate-census` puts the census it built back into the psp's trailer
    // (`replace_trailer`), so one write of a census is work that route pays for either way.
    //
    // **Neither side is encoded with a pileup identity** — the psp-header digest and record count
    // (`psp_census_pair.md` §3). Neither real route writes one now that both censuses live in a
    // trailer, so writing one here would make the two sides differ on a field neither is wrong
    // about.
    for (psp, census) in &to_copy_out {
        let trailer = PspReader::open(psp)?.trailer()?;
        std::fs::write(census, &trailer)
            .map_err(|source| format!("writing {}: {source}", census.display()))?;
    }

    println!(
        "route={} samples={} psp_bytes={psp_bytes} census_bytes={census_bytes} seconds={seconds:.2}",
        route.name(),
        paths.len(),
    );
    Ok(())
}

/// The first `count` intervals of `bed`, written to a file of their own.
fn first_regions_of(bed: &Path, count: usize, scratch: &Path) -> Result<PathBuf, Box<dyn Error>> {
    let whole = std::fs::read_to_string(bed)?;
    let kept: Vec<&str> = whole
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .take(count)
        .collect();
    if kept.is_empty() {
        return Err(format!("{} holds no intervals to analyse", bed.display()).into());
    }
    let trimmed = scratch.join("analysed.bed");
    std::fs::write(&trimmed, format!("{}\n", kept.join("\n")))?;
    Ok(trimmed)
}

/// What this harness records about itself in every psp it writes.
///
/// **The command line is a constant and not this process's own**, which the comparison the
/// wrapper makes depends on. A psp's header carries its provenance, so recording the route word
/// here would give the two routes psps that differ in their headers — and the wrapper's `cmp`
/// would then fail on the harness rather than on anything about the routes. Until plan step E2 a
/// census also named its psp by a digest of that header, so the two routes' *censuses* differed
/// too, which is how this was measured before the constant went in.
fn provenance() -> WriterProvenance {
    WriterProvenance {
        tool: "ng".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        subcommand: "ng_census_route_cost".to_string(),
        input_alignments: Vec::new(),
        input_reference: String::new(),
        command_line: "ng_census_route_cost".to_string(),
        parameters: std::collections::BTreeMap::from([(
            "depth-cap".to_string(),
            ParameterValue::Integer(300),
        )]),
        created: "2026-09-05T00:00:00Z".parse().expect("a datetime"),
    }
}
