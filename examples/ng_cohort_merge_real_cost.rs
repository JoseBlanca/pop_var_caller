//! What the cohort merge costs on **observations a real cohort's reads produced**, and what
//! producing them cost beside it.
//!
//! **Every number the merge has ever been measured on came from a fabricated fixture** — one
//! record per sample at the same positions, one base long, one observation each
//! (`examples/ng_cohort_merge_parallel_cost.rs`). That fixture cannot answer the two questions
//! the module's remaining decisions turn on, because both depend on what real observations
//! look like:
//!
//! - **How wide should a building region be?** The fabricated sweep says 200 bases on ground
//!   with a record every four and 1,000–2,000 on ground with a record every hundred. Those are
//!   an order of magnitude apart and nothing fabricated says which one real data resembles.
//! - **How often is a building region empty?** A builder over a region no observation begins in
//!   returns without opening its walk, which is worth about a third of the merge on the
//!   fabricated ground — but there every sample carries a record at the *same* positions, so a
//!   region is empty for the whole cohort at once. Real samples share most of their positions
//!   and not all of them, and with a thousand samples the union is far denser than any one
//!   sample's.
//!
//! And one question nothing in this module can answer at all: **is the merge worth optimising
//! further?** It is one stage of a run whose other stages have never been timed against it.
//! This walks the generic locus generator over real reads, which is the stage immediately
//! upstream, and prints what each cost.
//!
//! ```text
//! ng_cohort_merge_real_cost <reference.fa> <cram-dir> <regions.bed>
//! ```
//!
//! `NG_REAL_SAMPLES=n` walks only the first `n` CRAMs of the directory, in name order;
//! `NG_REAL_REGIONS=n` only the first `n` intervals of the BED. Both default to everything,
//! and both matter: the observations of every sample over every interval are held at once,
//! which is what the merge consumes and is the memory this probe's peak is made of.
//!
//! **The regions go to the generator as `Generic`, with no repeat catalog.** A run reads the
//! catalog beside its reference to route repeat tracts to the STR generator; there is none
//! beside the tomato reference, so the tracts here are walked as ordinary ground and produce
//! generic records where a run would produce microsatellite ones. What that changes is the
//! *shape* of a few records, not the density of them, which is what the merge's cost turns on.
//!
//! **The generator's time is one sample after another on one thread.** A run would walk the
//! cohort's samples in parallel, so the ratio printed here is the pessimistic one for the
//! merge — divide the generator's share by the threads a real run would give it.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, SampleLocusObservations, SampleLocusObservationsIterator,
    UnhandledReason,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegion};
use pop_var_caller::ng::run::cohort_merge::observation_cache::{
    ObservationCache, building_regions_of,
};
use pop_var_caller::ng::run::cohort_merge::parallel::merge_cohort_in_parallel;
use pop_var_caller::ng::run::cohort_merge::serial::{
    merge_cohort_serially, merge_cohort_through_cache,
};
use pop_var_caller::ng::run::cohort_merge::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltObs,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::reference_check_from_env;

/// This probe's sources cannot fail: the observations are already in memory.
#[derive(Debug)]
struct Never;

/// How many repeats each merge time is the median of.
const REPEATS: usize = 5;

/// The widths swept, in reference bases. 200 is the default; the rest bracket it either way.
const WIDTHS: [u32; 5] = [20, 100, 200, 500, 1_000];

/// One sample's observations over the analysed ground, and what walking its reads cost.
struct WalkedSample {
    observations: Vec<SampleLocusObservations>,
    seconds: f64,
}

/// The BED's intervals as the analysed regions of a merge, in the reference's contig order.
///
/// **One-based inclusive, where a BED is zero-based half-open** — `GenomeRegion` is the
/// caller's own coordinate system (`typed_regions.md` §1.1) and a probe that handed over BED
/// coordinates unchanged would measure a merge over ground shifted by a base.
fn analysed_regions_of(
    bed: &Path,
    contig_index: impl Fn(&str) -> Option<u32>,
    limit: Option<usize>,
) -> Result<Vec<GenomeRegion>, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(bed)?;
    let mut regions = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some(contig), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
        else {
            return Err(format!("a BED line has fewer than three fields: {line}").into());
        };
        let Some(contig) = contig_index(contig) else {
            return Err(format!("the reference has no contig named {contig}").into());
        };
        regions.push(GenomeRegion {
            contig: ContigId(contig),
            start: Position(start.parse::<u64>()? + 1),
            end: Position(end.parse::<u64>()?),
        });
    }
    regions.sort_by_key(|region| (region.contig.0, region.start.0));
    if let Some(limit) = limit {
        regions.truncate(limit);
    }
    Ok(regions)
}

/// Walk one sample's reads over `analysed` and keep every observation, in coordinate order.
///
/// This is [`ng_generic_walk_probe`](../ng_generic_walk_probe)'s pipeline — the same reference
/// reading, the same `SampleReads`, the same `LeftAlignPreparer`, the same
/// `PileupGeneratorConfig` — with the one difference that it *retains* what that probe drops,
/// because the retained observations are what this one is here to merge.
fn walk_one_sample(
    fasta: &Path,
    cram: &Path,
    analysed: &[GenomeRegion],
    cache: &Arc<ReferenceInfoCache>,
) -> Result<WalkedSample, Box<dyn std::error::Error>> {
    let check = reference_check_from_env()?;
    let (info, _verify) =
        read_reference_verifying_or_creating_fai(cache, fasta.to_path_buf(), check)?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));

    let reference = OpenReference::new(info);
    let reads = SampleReads::open_only_sample(
        &[cram.to_path_buf()],
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;

    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new takes Arc and this accessor is file-backed and single-threaded, as in ng_generic_walk_probe"
    )]
    let shared = Arc::new(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
    let make_reference = {
        let fasta = fasta.to_path_buf();
        let contigs = contigs.clone();
        let index = index.clone();
        move || WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
    };
    let generator = PileupGenerator::new(
        shared,
        make_reference,
        preparer,
        PileupGeneratorConfig::default(),
    )?;
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    let regions: Vec<Result<TypedRegion, pop_var_caller::ng::repeat_catalog::RepeatCatalogError>> =
        analysed
            .iter()
            .map(|region| {
                Ok(TypedRegion {
                    region: *region,
                    kind: RegionKind::Generic,
                })
            })
            .collect();

    let mut observations = Vec::new();
    let mut stream = SampleLocusObservationsIterator::new(regions.into_iter(), reads, generators);
    let started = Instant::now();
    for locus in &mut stream {
        observations.push(locus?);
    }
    Ok(WalkedSample {
        seconds: started.elapsed().as_secs_f64(),
        observations,
    })
}

/// The median, fastest and slowest of [`REPEATS`] runs, in milliseconds, with `prepare` outside
/// the clock — the discipline `ng_cohort_merge_parallel_cost` had to learn.
fn timed<T>(mut prepare: impl FnMut() -> T, mut one_merge: impl FnMut(T)) -> (f64, f64, f64) {
    one_merge(prepare());
    let mut each: Vec<f64> = (0..REPEATS)
        .map(|_| {
            let prepared = prepare();
            let started = Instant::now();
            one_merge(prepared);
            started.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    each.sort_by(f64::total_cmp);
    (each[each.len() / 2], each[0], each[each.len() - 1])
}

/// One source per sample over `cohort`, each a fresh copy the cache then owns.
fn sources_over(
    cohort: &[Vec<SampleLocusObservations>],
) -> Vec<std::vec::IntoIter<Result<SampleLocusObservations, Never>>> {
    cohort
        .iter()
        .map(|sample| {
            sample
                .iter()
                .cloned()
                .map(Ok)
                .collect::<Vec<_>>()
                .into_iter()
        })
        .collect()
}

/// How many of the building regions at `width` hold an observation beginning in them, and how
/// many do not — the hit rate of the skip in `build_region`, which the fabricated fixture makes
/// as favourable as it can possibly be.
fn regions_with_a_locus_start(
    analysed: &[GenomeRegion],
    cohort: &[Vec<SampleLocusObservations>],
    width: CohortLocusBuilderRegionsLen,
) -> (u64, u64) {
    // One pass per sample over its own records, marking the region each start falls in — rather
    // than asking each region about each sample, which is the cost the skip exists to avoid and
    // would make this probe slower than the merge it measures.
    let mut occupied: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();
    for sample in cohort {
        for record in sample {
            let start = record.region.start.min(record.region.end);
            occupied.insert((record.region.contig.0, start.get()));
        }
    }
    let (mut with, mut without) = (0u64, 0u64);
    for analysed_region in analysed {
        for building_region in building_regions_of(*analysed_region, width) {
            let holds = (building_region.start.get()..=building_region.end.get())
                .any(|base| occupied.contains(&(building_region.contig.0, base)));
            if holds { with += 1 } else { without += 1 }
        }
    }
    (with, without)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, fasta, crams, bed] = args.as_slice() else {
        eprintln!(
            "usage: ng_cohort_merge_real_cost <reference.fa> <cram-dir> <regions.bed>\n\
             NG_REAL_SAMPLES=n walks the first n CRAMs; NG_REAL_REGIONS=n the first n BED \
             intervals."
        );
        return ExitCode::from(2);
    };
    match run(Path::new(fasta), Path::new(crams), Path::new(bed)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            let mut source = error.source();
            while let Some(cause) = source {
                eprintln!("  caused by: {cause}");
                source = cause.source();
            }
            ExitCode::FAILURE
        }
    }
}

fn run(fasta: &Path, crams: &Path, bed: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let limit_of = |name: &str| -> Option<usize> {
        std::env::var(name)
            .ok()
            .map(|value| value.parse().expect("a count"))
    };

    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, _verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        reference_check_from_env()?,
    )?;
    let contigs = info.contig_list();
    let analysed = analysed_regions_of(
        bed,
        |name| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|at| at as u32)
        },
        limit_of("NG_REAL_REGIONS"),
    )?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();

    let mut cram_paths: Vec<PathBuf> = std::fs::read_dir(crams)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .is_some_and(|kind| kind == "cram" || kind == "bam")
        })
        .collect();
    cram_paths.sort();
    if let Some(limit) = limit_of("NG_REAL_SAMPLES") {
        cram_paths.truncate(limit);
    }
    if cram_paths.is_empty() {
        return Err(format!("no .cram or .bam under {}", crams.display()).into());
    }

    println!("# analysed intervals: {}", analysed.len());
    println!("# analysed bases: {analysed_bases}");
    println!("# samples: {}", cram_paths.len());
    println!("# threads available: {}", rayon::current_num_threads());

    let mut cohort: Vec<Vec<SampleLocusObservations>> = Vec::with_capacity(cram_paths.len());
    let mut generator_seconds = 0.0;
    for cram in &cram_paths {
        let walked = walk_one_sample(fasta, cram, &analysed, &cache)?;
        generator_seconds += walked.seconds;
        cohort.push(walked.observations);
    }

    let records: usize = cohort.iter().map(Vec::len).sum();
    let observations: usize = cohort
        .iter()
        .flatten()
        .map(|record| record.observations.len())
        .sum();
    let record_bases: u64 = cohort
        .iter()
        .flatten()
        .map(|record| record.region.len())
        .sum();
    println!("# generator seconds (one sample after another, one thread): {generator_seconds:.2}");
    println!("# records: {records}");
    println!(
        "# records per sample: {:.0}",
        records as f64 / cohort.len() as f64
    );
    println!(
        "# bases between one sample's records: {:.1}",
        analysed_bases as f64 / (records as f64 / cohort.len() as f64)
    );
    println!(
        "# observations per record: {:.2}",
        observations as f64 / records as f64
    );
    println!(
        "# bases per record: {:.2}",
        record_bases as f64 / records as f64
    );

    println!("\nwidth_bases, regions_with_a_start, regions_without");
    for bases in WIDTHS {
        let width =
            CohortLocusBuilderRegionsLen(std::num::NonZeroU32::new(bases).expect("non-zero"));
        let (with, without) = regions_with_a_locus_start(&analysed, &cohort, width);
        println!(
            "{bases}, {with}, {without}   # {:.0}% hold nothing",
            100.0 * without as f64 / (with + without) as f64
        );
    }

    let slices: Vec<&[SampleLocusObservations]> = cohort.iter().map(Vec::as_slice).collect();
    let (median, fastest, slowest) = timed(
        || (),
        |()| {
            std::hint::black_box(&merge_cohort_serially(
                &analysed,
                &slices,
                MaxCohortLocusSpan::DEFAULT,
                MinAltObs::DEFAULT,
            ));
        },
    );
    println!("\ndriver, region_bases, median_ms, min_ms, max_ms");
    println!("oracle, -, {median:.2}, {fastest:.2}, {slowest:.2}");

    for bases in WIDTHS {
        let width =
            CohortLocusBuilderRegionsLen(std::num::NonZeroU32::new(bases).expect("non-zero"));
        let (median, fastest, slowest) = timed(
            || ObservationCache::over(sources_over(&cohort)),
            |mut cache| {
                std::hint::black_box(
                    &merge_cohort_through_cache(
                        &analysed,
                        &mut cache,
                        width,
                        MaxCohortLocusSpan::DEFAULT,
                        MinAltObs::DEFAULT,
                    )
                    .expect("the probe's sources cannot fail"),
                );
            },
        );
        println!("one reader per sample, {bases}, {median:.2}, {fastest:.2}, {slowest:.2}");

        for threads in [1usize, 4, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("a pool of the asked-for size");
            let in_flight = CohortLocusBuilderRegionsInFlight(
                std::num::NonZeroUsize::new(threads).expect("non-zero"),
            );
            let (median, fastest, slowest) = pool.install(|| {
                timed(
                    || ObservationCache::over(sources_over(&cohort)),
                    |mut cache| {
                        std::hint::black_box(
                            &merge_cohort_in_parallel(
                                &analysed,
                                &mut cache,
                                width,
                                in_flight,
                                MaxCohortLocusSpan::DEFAULT,
                                MinAltObs::DEFAULT,
                            )
                            .expect("the probe's sources cannot fail"),
                        );
                    },
                )
            });
            println!("pool of {threads}, {bases}, {median:.2}, {fastest:.2}, {slowest:.2}");
        }
    }

    Ok(())
}
