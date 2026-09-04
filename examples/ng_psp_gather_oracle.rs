//! **Does a psp gathered from real reads hold exactly what the walk streamed?**
//!
//! Plan step B2's real-data oracle (`doc/devel/ng/impl_plan/run_driver_psp_mode.md`): one
//! sample's real alignment file, walked through the real repeat catalog, written to a psp by
//! [`SampleObservationGatherer::write_psp`] — and then the file read back and compared
//! **record for record, field for field** against the same sample walked again in memory.
//! The fixture-scale version of this comparison lives in `src/ng/run/gatherer.rs`'s tests;
//! what this harness adds is real sequenced DNA, ground with repeat tracts, gaps and contig
//! ends in it, and record values the synthetic fixtures never produce.
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example ng_psp_gather_oracle -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <cram-or-dir>
//! ```
//!
//! `NG_SAMPLES=n` gathers the first `n` alignment files of a directory (default 1);
//! `NG_REGIONS=n` analyses the first `n` BED intervals (default 2); `NG_WORK=dir` says where
//! the psp files land (default `tmp/ng_psp_gather_oracle`, inside the project); `NG_TWICE=1`
//! also gathers each sample a second time and compares the two files as bytes — plan step
//! B3's identity oracle on real reads (spec §12.1).
//!
//! # What each comparison proves
//!
//! - **the file read back == the walk in memory** — the store loses and changes nothing on
//!   real records: every coordinate, every allele sequence, every support count, every
//!   chain-id list, every tract motif and flank. This is the plan's north star for the walk
//!   stage. The in-memory side is a second gatherer over the same inputs; that a gatherer
//!   *is* the bare direct-mode chain is pinned separately, at fixture scale, by
//!   `the_gatherer_yields_what_the_direct_walk_yields`.
//! - **the header read back == the header the gatherer fixed** (plus the one compression
//!   parameter the store adds at create) — what a calling run will later check is what the
//!   walk wrote.
//!
//! The harness checks every sample and exits non-zero if any sample's file and walk
//! disagree, naming the record index and the first differing field group.
//!
//! # What it measured, 2026-09-03
//!
//! One tomato accession (`SRR7279481.p1.bench.cram`; its read group names the sample
//! `SRS3394712`) over the first two intervals of `benchmarks/tomato1/regions.bed` — 200 kb
//! of SL4.0 at about three reads a position — in the container, release build:
//!
//! - **183,807 records** (4 blocks, 948,689 bytes on disk), gathered in 0.32 s; the
//!   comparison walk another 0.20 s; setup (reference + catalog + segments) 2.75 s;
//! - **1,217 of the walk's records are repeat tracts**, so the tract path is on this oracle,
//!   not only the generic one;
//! - 2,713 of the 2,991 segments were walked by a filled generator (the rest are repeat
//!   bundles and satellites, refused by kind as unbuilt/out of scope);
//! - **the file IS the walk**: all 183,807 records equal field for field, header equal;
//! - under `NG_TWICE=1`, **the second gather is byte-identical** — all 948,689 bytes,
//!   timestamp included, because the timestamp is this harness's own fixed value.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::{LocusKind, SampleLocusObservations};
use pop_var_caller::ng::psp::WriterProvenance;
use pop_var_caller::ng::psp::{
    Header, ParameterValue, PspReader, WriteStats, ZSTD_COMPRESSION_LEVEL,
    ZSTD_COMPRESSION_LEVEL_KEY,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::GenomeRegions;
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::run::{
    RunError, SampleObservationGatherer, SampleWalkInputs, Segmentation,
};
use pop_var_caller::regions::ContigBounds;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::pileup::PileupGeneratorConfig;

const SAMPLES_BY_DEFAULT: usize = 1;
const REGIONS_BY_DEFAULT: usize = 2;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [fasta, catalog, bed, crams] = args.as_slice() else {
        eprintln!(
            "usage: ng_psp_gather_oracle <reference.fa> <catalog.parquet> <regions.bed> \
             <cram-or-dir>\n\
             gathers each sample to a psp and compares the file, record for record, against \
             the same sample walked in memory.\n\
             NG_SAMPLES=n takes the first n files of a directory (default {SAMPLES_BY_DEFAULT}); \
             NG_REGIONS=n the first n BED intervals (default {REGIONS_BY_DEFAULT}); \
             NG_WORK=dir is where the psp files land (default tmp/ng_psp_gather_oracle); \
             NG_TWICE=1 also gathers each sample a second time and compares the two files \
             as bytes (plan step B3's identity oracle)."
        );
        return ExitCode::from(2);
    };
    match run(
        Path::new(fasta),
        Path::new(catalog),
        Path::new(bed),
        Path::new(crams),
    ) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
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
///
/// **Fallible rather than panicking**, so a mistyped variable arrives as this harness's
/// ordinary exit-1 error naming the variable and what it held, not as an exit-101 panic
/// naming neither.
fn how_many(name: &str, fallback: usize) -> Result<usize, Box<dyn std::error::Error>> {
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
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    crams: &Path,
) -> Result<bool, Box<dyn std::error::Error>> {
    let setting_up = Instant::now();
    // The reference with its FASTA read to the end: the catalog check compares checksums,
    // and a reference whose background verification has not been joined carries none.
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let with_checksums = match verify {
        Some(handle) => handle.join()?,
        None => Arc::clone(&info),
    };
    let contigs: ContigList = info.contig_list();
    let reference = OpenReference::new(info);

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: entry.length as u32,
        })
        .collect();
    let scratch = tempfile::tempdir()?;
    let trimmed = first_regions_of(
        bed,
        how_many("NG_REGIONS", REGIONS_BY_DEFAULT)?,
        scratch.path(),
    )?;
    let analysed = GenomeRegions::from_bed_path(&trimmed, &bounds)?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();

    let criteria = StrRepeatCriteria::default();
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog.genome_segments(&criteria, ReadScope::Regions(&spans))?;
    let segmentation = Arc::new(Segmentation::build(
        segments,
        analysed,
        catalog.header().clone(),
        criteria,
        catalog_path.to_path_buf(),
    )?);

    // **Every entry or none.** An entry that fails to read is an error, not a file to skip:
    // this harness's whole claim is that exit 0 means every sample was gathered and checked,
    // and a silently dropped entry would be a sample that escaped it.
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
        std::env::var("NG_WORK").unwrap_or_else(|_| "tmp/ng_psp_gather_oracle".to_string()),
    );
    std::fs::create_dir_all(&work_dir)?;

    println!("# reference: {}", fasta.display());
    println!("# analysed intervals: {}", spans.len());
    println!("# analysed bases: {analysed_bases}");
    println!("# segments: {}", segmentation.segments().len());
    println!("# samples: {}", paths.len());
    println!("# setup: {:.2} s", setting_up.elapsed().as_secs_f64());

    let mut all_agree = true;
    for path in &paths {
        all_agree &= gather_and_compare(path, &reference, &segmentation, &work_dir)?;
    }
    Ok(all_agree)
}

/// Gather `path`'s sample to a psp and compare the file against the same walk in memory.
/// `Ok(true)` when they agree everywhere.
fn gather_and_compare(
    path: &Path,
    reference: &OpenReference,
    segmentation: &Arc<Segmentation>,
    work_dir: &Path,
) -> Result<bool, Box<dyn std::error::Error>> {
    let alignments = [path.to_path_buf()];
    let open_gatherer = || {
        SampleObservationGatherer::open(
            SampleWalkInputs {
                alignments: &alignments,
                reference,
                read_filters: ReadFilterConfig::default(),
                locus_generator_settings: PileupGeneratorConfig::default(),
                build_index_if_missing: false,
            },
            Arc::clone(segmentation),
            provenance(),
            None,
        )
    };

    let gathering = Instant::now();
    let gatherer = open_gatherer()?;
    let sample = gatherer.sample_name().to_string();
    // Taken before `write_psp` consumes the gatherer, so the comparison costs one open
    // rather than two — on real data an open re-reads the alignment headers.
    let header_the_gatherer_fixed = gatherer.header().clone();
    let psp_path = work_dir.join(format!("{sample}.psp"));
    let (stats, counts) = gatherer.write_psp(&psp_path, None)?;
    let gather_seconds = gathering.elapsed().as_secs_f64();

    let walking = Instant::now();
    let walked: Vec<SampleLocusObservations> = open_gatherer()?.collect::<Result<_, _>>()?;
    let walk_seconds = walking.elapsed().as_secs_f64();

    let tract_records = walked
        .iter()
        .filter(|observation| matches!(observation.kind, LocusKind::Ssr(_)))
        .count();
    println!(
        "\n== {sample}: {} records ({} blocks, {} bytes) in {gather_seconds:.2} s; \
         comparison walk {walk_seconds:.2} s; {} of the walk's records are repeat tracts; \
         regions walked {} of {}",
        stats.records,
        stats.blocks,
        stats.bytes,
        tract_records,
        counts.regions_handled,
        counts.regions_in,
    );
    // **The tract path is part of what this oracle claims to cover**, so a walk without a
    // single tract record is a refusal rather than a quiet pass: the ground moved, or
    // routing dropped the tracts from walk and file alike.
    if tract_records == 0 {
        println!(
            "!! no repeat-tract records on this walk — the tract path is NOT on this oracle; \
             widen NG_REGIONS or choose tract-bearing intervals"
        );
        return Ok(false);
    }

    let mut reader = PspReader::open(&psp_path)?;
    if !header_matches(&mut reader, header_the_gatherer_fixed) {
        return Ok(false);
    }
    if !gathered_twice_is_byte_identical(&open_gatherer, work_dir, &sample, &psp_path, &stats)? {
        return Ok(false);
    }
    if !file_matches_walk(&mut reader, &walked, &sample)? {
        return Ok(false);
    }
    Ok(true)
}

/// The header read back is the one the gatherer fixed, plus the compression level the store
/// records at `create`. `false` after printing which field differs.
fn header_matches(reader: &mut PspReader, expected: Header) -> bool {
    let mut expected = expected;
    expected.writer.parameters.insert(
        ZSTD_COMPRESSION_LEVEL_KEY.to_string(),
        ParameterValue::Integer(i64::from(ZSTD_COMPRESSION_LEVEL)),
    );
    let found = reader.header();
    if found == &expected {
        return true;
    }
    println!("!! the header read back is not the header the gatherer fixed:");
    // Named field by field rather than as two whole dumps, so the line says what moved.
    for (what, differs) in [
        ("sample", found.sample != expected.sample),
        ("reference", found.reference != expected.reference),
        ("contigs", found.contigs != expected.contigs),
        ("read groups", found.read_groups != expected.read_groups),
        (
            "observation reach ceiling",
            found.observation_reach_ceiling_bp != expected.observation_reach_ceiling_bp,
        ),
        ("writer provenance", found.writer != expected.writer),
        (
            "segmentation inputs",
            found.segmentation_inputs != expected.segmentation_inputs,
        ),
        ("manifest", found.manifest != expected.manifest),
        (
            "format version",
            found.format_version != expected.format_version,
        ),
    ] {
        if differs {
            println!("   {what} differs");
        }
    }
    false
}

/// Plan step B3's oracle on real reads: gather the same sample a second time and compare
/// whole files as bytes. Skipped unless `NG_TWICE=1`.
///
/// The provenance — timestamp included — is this harness's own fixed value, so identity must
/// hold over **every** byte here, not "all but the timestamp" as spec §12.1 allows a run
/// that stamps the clock.
fn gathered_twice_is_byte_identical(
    open_gatherer: &dyn Fn() -> Result<SampleObservationGatherer, RunError>,
    work_dir: &Path,
    sample: &str,
    first_psp: &Path,
    first_stats: &WriteStats,
) -> Result<bool, Box<dyn std::error::Error>> {
    if std::env::var("NG_TWICE").is_ok_and(|value| value == "1") {
        let again_path = work_dir.join(format!("{sample}.again.psp"));
        let (again_stats, _) = open_gatherer()?.write_psp(&again_path, None)?;
        let first_bytes = std::fs::read(first_psp)?;
        let again_bytes = std::fs::read(&again_path)?;
        if first_bytes != again_bytes {
            println!(
                "!! two gathers differ: {} vs {} bytes ({} vs {} records)",
                first_bytes.len(),
                again_bytes.len(),
                first_stats.records,
                again_stats.records,
            );
            return Ok(false);
        }
        println!(
            "== {sample}: gathered twice, byte-identical ({} bytes)",
            first_bytes.len(),
        );
    }
    Ok(true)
}

/// B2's oracle proper: every record in the file equals the record the walk streamed at that
/// position, and the two hold the same number of records.
fn file_matches_walk(
    reader: &mut PspReader,
    walked: &[SampleLocusObservations],
    sample: &str,
) -> Result<bool, Box<dyn std::error::Error>> {
    let mut records_read_back = 0usize;
    for (index, streamed) in reader.records()?.enumerate() {
        // PANIC-FREE: `records()` takes no predicate, and a body is `None` only where a
        // predicate skipped it (`psp::block::StreamedRecord`).
        let record = streamed?.record.expect("records() builds every body");
        let Some(walked_record) = walked.get(index) else {
            println!("!! the file holds more records than the walk: index {index}");
            return Ok(false);
        };
        if &record != walked_record {
            println!("!! record {index} differs between the file and the walk:");
            report_first_difference(&record, walked_record);
            return Ok(false);
        }
        records_read_back += 1;
    }
    if records_read_back != walked.len() {
        println!(
            "!! the file holds {records_read_back} records where the walk streamed {}",
            walked.len(),
        );
        return Ok(false);
    }
    println!("== {sample}: the file IS the walk — {records_read_back} records equal, header equal");
    Ok(true)
}

/// Say which part of a record disagrees, so a failure is a lead rather than a riddle.
fn report_first_difference(
    from_file: &SampleLocusObservations,
    from_walk: &SampleLocusObservations,
) {
    // No `..`: a field added to the record must be dispositioned here, or a difference in
    // it would be reported as "the observations differ" and send the reader hunting.
    let SampleLocusObservations {
        region: _,
        reference_bases: _,
        observations: _,
        reads_without_observation: _,
        reads_discarded_by_cap: _,
        kind: _,
    } = from_file;
    if from_file.region != from_walk.region {
        println!(
            "   region: file {:?} vs walk {:?}",
            from_file.region, from_walk.region
        );
    } else if from_file.kind != from_walk.kind {
        println!("   locus kind differs at {:?}", from_file.region);
    } else if from_file.reference_bases != from_walk.reference_bases {
        println!("   reference bases differ at {:?}", from_file.region);
    } else if from_file.reads_without_observation != from_walk.reads_without_observation
        || from_file.reads_discarded_by_cap != from_walk.reads_discarded_by_cap
    {
        println!(
            "   the two no-observation counts differ at {:?}",
            from_file.region
        );
    } else {
        println!(
            "   the observations differ at {:?} ({} in the file, {} in the walk)",
            from_file.region,
            from_file.observations.len(),
            from_walk.observations.len(),
        );
    }
}

/// Fixed provenance, so two runs of this harness produce comparable files: the harness is
/// the caller, and the caller supplies what only it knows (the gatherer overwrites the
/// input names and records the filters itself).
fn provenance() -> WriterProvenance {
    WriterProvenance {
        tool: "ng".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        subcommand: "generate-psps".to_string(),
        input_alignments: Vec::new(),
        input_reference: String::new(),
        command_line: std::env::args().collect::<Vec<_>>().join(" "),
        parameters: std::collections::BTreeMap::new(),
        created: "2026-09-03T00:00:00Z".parse().expect("a datetime"),
    }
}

/// The first `how_many` intervals of `bed`, written to a trimmed copy under `scratch` —
/// truncating rather than sampling keeps a re-run comparable with the one before it.
fn first_regions_of(
    bed: &Path,
    count: usize,
    scratch: &Path,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
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
