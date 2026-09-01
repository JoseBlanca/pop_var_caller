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

// **The allocator this run measures, named in its own output.** The merge's Linux profile is
// 39% allocator, so which allocator is installed is part of what every number here means.
// `--features alloc-mimalloc` swaps it and `--features dhat-heap` counts it; `ALLOCATOR` below
// is printed so a run cannot claim one allocator while carrying another.
//
// Only one `#[global_allocator]` may exist and dhat takes the slot when both features are on,
// which is why the mimalloc arm excludes it rather than colliding with it.
#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Which global allocator this build installed.
#[cfg(feature = "dhat-heap")]
const ALLOCATOR: &str = "dhat (counting)";
/// Which global allocator this build installed.
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
const ALLOCATOR: &str = "mimalloc";
/// Which global allocator this build installed.
#[cfg(not(any(feature = "alloc-mimalloc", feature = "dhat-heap")))]
const ALLOCATOR: &str = "system";

/// How many heap blocks this process has allocated so far, counting from its start.
///
/// **Zero without the `dhat-heap` feature**, because nothing is counting: the project forbids
/// `unsafe`, so the only counting allocator available is dhat's. A run that wants the count
/// asks for it at the build.
fn allocated_blocks() -> u64 {
    #[cfg(feature = "dhat-heap")]
    {
        dhat::HeapStats::get().total_blocks
    }
    #[cfg(not(feature = "dhat-heap"))]
    {
        0
    }
}

/// How many heap blocks are live right now. Zero without `dhat-heap`, as above.
///
/// **This is what prices the merge's deallocation**, which is most of what it asks the
/// allocator for and which no count of allocations can see: the records a cached driver evicts
/// were allocated before the clock started, so they raise nothing in `allocated_blocks` and
/// lower this.
fn live_blocks() -> u64 {
    #[cfg(feature = "dhat-heap")]
    {
        dhat::HeapStats::get().curr_blocks as u64
    }
    #[cfg(not(feature = "dhat-heap"))]
    {
        0
    }
}

/// How many heap bytes this process has allocated so far. Zero without `dhat-heap`, as above.
fn allocated_bytes() -> u64 {
    #[cfg(feature = "dhat-heap")]
    {
        dhat::HeapStats::get().total_bytes
    }
    #[cfg(not(feature = "dhat-heap"))]
    {
        0
    }
}

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
use pop_var_caller::ng::run::cohort_merge::build::RegionOutcome;
use pop_var_caller::ng::run::cohort_merge::close::{LocusCloser, Verdict};
use pop_var_caller::ng::run::cohort_merge::observation_cache::{
    ObservationCache, ObservationSource, building_regions_of,
};
use pop_var_caller::ng::run::cohort_merge::parallel::merge_cohort_in_parallel;
use pop_var_caller::ng::run::cohort_merge::serial::{
    merge_cohort_serially, merge_cohort_through_cache,
};
use pop_var_caller::ng::run::cohort_merge::timing as merge_timing;
use pop_var_caller::ng::run::cohort_merge::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltObs,
    MinAltReadShare, MinAltReads,
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

/// One merge's answer as text, one line per built locus and one per refused span.
///
/// **This is the module's own definition of "the same answer"** — `cohort_merge`'s test fixtures
/// render an outcome exactly this way, through each locus's `Debug`, so a field the type gains is
/// compared here without this probe being told about it.
fn rendered(outcome: &RegionOutcome) -> Vec<String> {
    let RegionOutcome {
        cohort_observations,
        failed_locus_spans,
    } = outcome;
    cohort_observations
        .iter()
        .map(|observed| format!("{observed:?}"))
        .chain(
            failed_locus_spans
                .iter()
                .map(|span| format!("failed {span}")),
        )
        .collect()
}

/// Refuse any difference between the oracle's answer and a driver's, naming the first line that
/// differs.
///
/// **The oracle is `merge_cohort_serially`**, which holds every sample's observations at once and
/// does not divide the ground at all, so it shares no scheduling machinery with what it checks.
fn refuse_any_difference(what: &str, oracle: &[String], actual: &[String]) -> Result<(), String> {
    if let Some(first) = oracle
        .iter()
        .zip(actual)
        .position(|(one, other)| one != other)
    {
        return Err(format!(
            "{what} does not give the oracle's answer, first at locus {first}:\n  oracle: {}\n  \
             {what}: {}",
            oracle[first], actual[first],
        ));
    }
    if oracle.len() != actual.len() {
        return Err(format!(
            "{what} built {} loci where the oracle built {}",
            actual.len(),
            oracle.len(),
        ));
    }
    Ok(())
}

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

/// The most memory this process has ever held, as the kernel counted it.
///
/// **The allocator question needs this beside the clock.** Swapping the allocator is the
/// largest single lever measured on this merge, and the project's reason for existing is
/// trading memory for sample-count scaling — so a speed-up that costs resident memory is a
/// different decision from one that does not, and the probe should not make the reader go and
/// find out separately.
///
/// Read from `/proc/self/status`, which only Linux has; the container is where the release
/// numbers are taken, so that is where it answers. Elsewhere it says so rather than guessing.
fn peak_resident() -> String {
    match std::fs::read_to_string("/proc/self/status") {
        Ok(status) => status
            .lines()
            .find_map(|line| line.strip_prefix("VmHWM:"))
            .map_or_else(
                || "unknown (no VmHWM in /proc/self/status)".to_string(),
                |value| value.trim().to_string(),
            ),
        Err(_) => "not measured (no /proc — this is not Linux)".to_string(),
    }
}

/// One sample's records, produced one at a time as the cache asks for them.
///
/// **Both ways of producing them do it inside the clock, and that is deliberate.** The cache
/// consumes its readers, so every round needs the sample's records made again; producing them
/// up front and handing over a finished vector charges the merge for none of it, and the
/// question here is precisely what producing a record costs. So both variants below mint on
/// demand, and the only difference between them is where the memory comes from.
enum ProbeSource<'a> {
    /// A fresh record every time, which is what the generator does today.
    Minting {
        template: &'a [SampleLocusObservations],
        at: usize,
    },
    /// The record the merge handed back, filled again — what a generator that leased its
    /// records would do, standing in for one so the saving can be measured before it is built.
    Leasing {
        template: &'a [SampleLocusObservations],
        at: usize,
    },
    /// **Records made before the clock started, handed over by moving them.**
    ///
    /// Neither source above can say what the merge's *freeing* costs, because both build the
    /// record inside the merge's clock — one clones a template, one refills the returned
    /// buffer — and on this ground those two cost about the same, so their arms tie and the
    /// free is invisible. A real run pays neither: the generator fills the record once,
    /// upstream. These two arms do that. They differ in one thing only:
    ///
    /// - `Handing` lets the merge drop the record it hands back, which is what happens today;
    /// - `Hoarding` keeps it instead, so nothing is freed inside the clock.
    ///
    /// **The difference between them is the merge's free, and nothing else.** What `Hoarding`
    /// costs is memory — it holds every record the merge gives back — which is why it is a
    /// measuring device and not a design.
    Handing {
        made: std::vec::IntoIter<SampleLocusObservations>,
    },
    /// [`Handing`](Self::Handing) keeping what the merge hands back rather than dropping it.
    Hoarding {
        made: std::vec::IntoIter<SampleLocusObservations>,
        kept: Vec<SampleLocusObservations>,
    },
}

impl ObservationSource for ProbeSource<'_> {
    type Error = Never;

    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, Never>> {
        // The two arms that hand over a record made earlier answer first: they neither
        // clone nor refill, so the only allocator work left inside the clock is the merge's.
        match self {
            ProbeSource::Handing { made } => {
                // Dropped here, which is where the merge's free happens for a record the
                // cache offered back.
                drop(spare);
                return made.next().map(Ok);
            }
            ProbeSource::Hoarding { made, kept } => {
                if let Some(spare) = spare {
                    kept.push(spare);
                }
                return made.next().map(Ok);
            }
            _ => {}
        }
        let (template, at, leasing) = match self {
            ProbeSource::Minting { template, at } => (*template, at, false),
            ProbeSource::Leasing { template, at } => (*template, at, true),
            ProbeSource::Handing { .. } | ProbeSource::Hoarding { .. } => {
                unreachable!("the two arms that hand over a made record returned above")
            }
        };
        let next = template.get(*at)?;
        *at += 1;
        Some(Ok(match (leasing, spare) {
            (true, Some(mut spare)) => {
                refill(&mut spare, next);
                spare
            }
            (_, spare) => {
                drop(spare);
                next.clone()
            }
        }))
    }
}

/// Fill `into` with what `from` holds, keeping every buffer of `into` that is the right size.
///
/// On generic ground every record covers one base and carries about one sequence, so nothing
/// here reallocates after the first window — which is the whole point of the exercise.
fn refill(into: &mut SampleLocusObservations, from: &SampleLocusObservations) {
    into.region = from.region;
    if into.reference_bases.len() == from.reference_bases.len() {
        into.reference_bases.copy_from_slice(&from.reference_bases);
    } else {
        into.reference_bases = from.reference_bases.clone();
    }
    into.reads_without_observation = from.reads_without_observation;
    into.reads_discarded_by_cap = from.reads_discarded_by_cap;
    into.kind = from.kind.clone();

    let reused = into.observations.len().min(from.observations.len());
    into.observations.truncate(reused);
    for (slot, source) in into.observations.iter_mut().zip(&from.observations) {
        if slot.bases.len() == source.bases.len() {
            slot.bases.copy_from_slice(&source.bases);
        } else {
            slot.bases = source.bases.clone();
        }
        slot.read_witness = source.read_witness.clone();
        slot.read_group = source.read_group;
        slot.num_obs = source.num_obs;
        slot.num_fwd = source.num_fwd;
        slot.q_sum = source.q_sum;
        slot.mapq_sum = source.mapq_sum;
        slot.mapq_sum_sq = source.mapq_sum_sq;
        slot.placed_left = source.placed_left;
        slot.chain_ids.clear();
        slot.chain_ids.extend_from_slice(&source.chain_ids);
    }
    for extra in &from.observations[reused..] {
        into.observations.push(extra.clone());
    }
}

/// One source per sample over `cohort`, minting a record a draw or refilling the one the
/// merge handed back.
///
/// **The choice is the caller's rather than the environment's**, so that one process can run
/// both and the two are compared on the same machine in the same minute. `NG_REAL_LEASE` is
/// what sets it, and in the one-driver arm it takes a list — `0,1` runs both.
// **`to_vec().into_iter()` is deliberate in two arms below, and clippy's `iter().cloned()` would
// break what they measure.** That rewrite makes each record's copy happen *inside* the merge's
// clock, one per draw, which is the cost these arms exist to keep out of it — see the comment on
// the arms themselves.
#[allow(
    clippy::unnecessary_to_owned,
    reason = "the copy is made up front on purpose; cloning per draw is the thing being measured"
)]
fn sources_over(cohort: &[Vec<SampleLocusObservations>], supply: Supply) -> Vec<ProbeSource<'_>> {
    cohort
        .iter()
        .map(|sample| match supply {
            Supply::Minted => ProbeSource::Minting {
                template: sample,
                at: 0,
            },
            Supply::Leased => ProbeSource::Leasing {
                template: sample,
                at: 0,
            },
            // **The copy is made here, outside the merge's clock**, which is the whole point of
            // these two arms: a real run's generator makes the record once and the merge is not
            // charged for it.
            Supply::Handed => ProbeSource::Handing {
                made: sample.to_vec().into_iter(),
            },
            Supply::Hoarded => ProbeSource::Hoarding {
                made: sample.to_vec().into_iter(),
                kept: Vec::new(),
            },
        })
        .collect()
}

/// Where the merge's records come from, and what happens to the ones it gives back.
///
/// **The last two are a measuring device, not designs.** They exist to price the merge's
/// freeing on its own: both hand over a record made before the clock started, and they differ
/// only in whether the merge's returned record is dropped or kept.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Supply {
    /// A fresh record every draw — what the generator does today.
    Minted,
    /// The returned record filled again.
    Leased,
    /// A record made earlier, handed over; the returned one dropped.
    Handed,
    /// A record made earlier, handed over; the returned one kept.
    Hoarded,
}

impl Supply {
    /// Read one setting from the `NG_REAL_LEASE` list.
    fn from_name(name: &str) -> Self {
        match name.trim() {
            "0" | "mint" | "minted" => Self::Minted,
            "1" | "lease" | "leased" => Self::Leased,
            "hand" | "handed" => Self::Handed,
            "hoard" | "hoarded" => Self::Hoarded,
            other => panic!("NG_REAL_LEASE takes mint, lease, hand or hoard, not {other}"),
        }
    }

    /// What the run's table calls it.
    fn name(self) -> &'static str {
        match self {
            Self::Minted => "minted",
            Self::Leased => "leased",
            Self::Handed => "handed",
            Self::Hoarded => "hoarded",
        }
    }
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
    // The counting allocator's bookkeeping, which `allocated_blocks` reads. Dropping it at the
    // end of `main` writes `dhat-heap.json` beside the run.
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

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
    println!("# global allocator: {ALLOCATOR}");

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

    // **The keep rule, and the two knobs that sweep it** (`cohort_merge.md` §4.3).
    // `NG_REAL_MIN_ALT=n` sets the floor and `NG_REAL_ALT_SHARE=f` the share, so what the
    // built loci cost can be read off by raising either until nothing builds.
    let min_alt_reads = MinAltReads {
        floor: match std::env::var("NG_REAL_MIN_ALT").ok() {
            Some(raw) => MinAltObs(
                std::num::NonZeroU32::new(raw.parse().expect("NG_REAL_MIN_ALT is a number"))
                    .expect("NG_REAL_MIN_ALT is non-zero"),
            ),
            None => MinAltObs::DEFAULT,
        },
        share: match std::env::var("NG_REAL_ALT_SHARE").ok() {
            Some(raw) => MinAltReadShare::new(raw.parse().expect("NG_REAL_ALT_SHARE is a number"))
                .expect("NG_REAL_ALT_SHARE is a fraction of one"),
            None => MinAltReadShare::DEFAULT,
        },
    };

    // **What the merge discards before assembling anything, and the depth the share is
    // taken of.** The keep rule asks each sample about its own compared reads, so the
    // mean below is the number that decides which half of the rule is doing the work:
    // under about a hundred the floor decides, above it the share does.
    {
        let all: Vec<&[SampleLocusObservations]> = cohort.iter().map(Vec::as_slice).collect();
        let (mut built, mut quiet, mut failed) = (0usize, 0usize, 0usize);
        let (mut compared_total, mut contributing_members) = (0u64, 0u64);
        for locus in LocusCloser::over(&all, MaxCohortLocusSpan::DEFAULT, min_alt_reads) {
            match locus.verdict {
                Verdict::Build => built += 1,
                Verdict::TooQuiet => quiet += 1,
                Verdict::Failed => failed += 1,
            }
            for member in &locus.members {
                compared_total += member
                    .observations
                    .iter()
                    .map(|obs| u64::from(obs.reads_compared_with_reference()))
                    .sum::<u64>();
                contributing_members += 1;
            }
        }
        let closed = built + quiet + failed;
        println!(
            "# keep rule: floor {} reads or {:.1}% of a sample's compared reads",
            min_alt_reads.floor.get(),
            100.0 * min_alt_reads.share.get(),
        );
        println!("# cohort loci closed: {closed}");
        println!(
            "# of those built: {built}, too quiet: {quiet} ({:.1}%), failed on width: {failed}",
            100.0 * quiet as f64 / closed.max(1) as f64
        );
        println!(
            "# mean compared reads per sample at a closed locus: {:.1}, so the rule asks {} \
             of a sample at that depth",
            compared_total as f64 / contributing_members.max(1) as f64,
            min_alt_reads.required_of(
                (compared_total / contributing_members.max(1)).min(u64::from(u32::MAX)) as u32
            ),
        );
    }
    // **One driver by itself, for a profiler to look at.** A sampling profiler attributes
    // whatever the process was doing when it interrupted it, so a run that walks the reads and
    // then sweeps five widths across three pool sizes gives a tree in which the merge is a
    // minority of the samples and each driver a minority of that. `NG_REAL_ONLY=oracle|cache|
    // parallel` runs one driver and nothing else, `NG_REAL_ROUNDS` times (default 20), at
    // `NG_REAL_WIDTH` bases (default 200) on `NG_REAL_THREADS` threads (default 8).
    //
    // **Rebuilding the sources between rounds is unavoidable and is not hidden**: the cache
    // consumes its readers, so every round after the first needs a fresh copy of every
    // sample's observations. That copying is real time in the process and shows in the profile
    // under `sources_over`, which is a subtree of its own — read the merge's share under
    // `merge_cohort_*` and leave the copying where it is. The oracle needs no copy at all.
    if let Ok(driver) = std::env::var("NG_REAL_ONLY") {
        let rounds = limit_of("NG_REAL_ROUNDS").unwrap_or(20).max(1);
        let bases = u32::try_from(limit_of("NG_REAL_WIDTH").unwrap_or(200)).expect("a width");
        // **A list, and every entry is run in this one process.** Comparing thread counts
        // across container starts compares two machines: this box has drifted by a factor of
        // two between runs of one unchanged binary, which is larger than the effect the sweep
        // is looking for. `NG_REAL_THREADS=1,2,4,8` puts them in one sitting.
        let thread_counts: Vec<usize> = match std::env::var("NG_REAL_THREADS") {
            Ok(list) => list
                .split(',')
                .map(|count| {
                    count
                        .trim()
                        .parse::<usize>()
                        .expect("a thread count")
                        .max(1)
                })
                .collect(),
            Err(_) => vec![8],
        };
        // **Both settings of the record supply in one process**, for the reason the thread
        // list is one: `0` mints a record a draw, which is what the generator does today, and
        // `1` refills the record the merge handed back. `NG_REAL_LEASE=0,1` runs both.
        let lease_settings: Vec<Supply> = match std::env::var("NG_REAL_LEASE") {
            Ok(list) => list.split(',').map(Supply::from_name).collect(),
            Err(_) => vec![Supply::Minted],
        };
        let width =
            CohortLocusBuilderRegionsLen(std::num::NonZeroU32::new(bases).expect("non-zero"));
        let slices: Vec<&[SampleLocusObservations]> = cohort.iter().map(Vec::as_slice).collect();

        // **The oracle check, and it runs before any clock.** The parallel driver divides the
        // ground and hands pieces of it to builders; the oracle divides nothing and holds the
        // whole cohort at once, so it is the only comparison that can catch a fault the dividing
        // itself introduces. The suite makes it on fixtures — this makes it on observations real
        // reads produced. `NG_REAL_ORACLE_CHECK=0` turns it off for a run that only wants a clock.
        if std::env::var("NG_REAL_ORACLE_CHECK").as_deref() != Ok("0") {
            let oracle = rendered(&merge_cohort_serially(
                &analysed,
                &slices,
                MaxCohortLocusSpan::DEFAULT,
                min_alt_reads,
            ));
            for regions in [1usize, 4, 16] {
                let in_flight = CohortLocusBuilderRegionsInFlight(
                    std::num::NonZeroUsize::new(regions).expect("non-zero"),
                );
                // Minted, which is what a run does today; the check is about the schedule.
                let mut cache = ObservationCache::over(sources_over(&cohort, Supply::Minted));
                let merged = merge_cohort_in_parallel(
                    &analysed,
                    &mut cache,
                    width,
                    in_flight,
                    MaxCohortLocusSpan::DEFAULT,
                    min_alt_reads,
                )
                .expect("the probe's sources cannot fail");
                refuse_any_difference(
                    &format!("the parallel driver at {bases} bases, {regions} in flight"),
                    &oracle,
                    &rendered(&merged),
                )?;
            }
            println!(
                "# oracle check: the parallel driver gives the oracle's answer at {bases} bases \
                 and 1, 4 and 16 regions in flight, over {} loci",
                oracle.len(),
            );
        }

        println!("# profile-start: {driver}");
        println!(
            "\ndriver, rounds, width_bases, threads, record_supply, loci_per_round, \
             median_ms, min_ms, max_ms, seconds_all_rounds_including_copies"
        );

        for threads in thread_counts {
            for supply in lease_settings.iter().copied() {
                let in_flight = CohortLocusBuilderRegionsInFlight(
                    std::num::NonZeroUsize::new(
                        limit_of("NG_REAL_IN_FLIGHT").unwrap_or(threads).max(1),
                    )
                    .expect("non-zero"),
                );
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("a pool of the asked-for size");

                // **The cohort's observations are copied between rounds and that is not the merge.**
                // A cache consumes its readers, so every round after the first needs a fresh copy of
                // every sample's records; timing it with the merge would charge the two cached drivers
                // for work the oracle is never asked to do. The copy happens before the inner clock
                // starts, so `merge_seconds` below is the merge and nothing else.
                // **Where the merge's own wall time went**, summed over the rounds below and printed
                // after them. Every counter is zero unless the build asked for `--features
                // merge-timing`, in which case the merge itself is what timed each part
                // (`pop_var_caller::ng::run::cohort_merge::timing`) — no sampling, no attribution.
                merge_timing::reset();
                let started = Instant::now();
                // **Every round's own time, not a running sum.** One descheduled round moves a mean
                // and leaves nothing to show it did, and this machine's swing between two runs of one
                // unchanged binary has been measured at 13–26% — larger than most of the effects this
                // probe is used to judge. The three other probes in this module already print the
                // median with its two extremes; this is the one that did not.
                let mut each_round: Vec<f64> = Vec::with_capacity(rounds);
                let mut loci = 0usize;
                // **What the merge allocates, which is the same on every run of the same code on the
                // same input where the milliseconds are not.** Counted only under `--features
                // dhat-heap`; zero otherwise, and printed only when it is non-zero.
                let mut blocks = 0u64;
                let mut bytes = 0u64;
                let mut freed = 0u64;
                for _ in 0..rounds {
                    match driver.as_str() {
                        "oracle" => {
                            let blocks_at = allocated_blocks();
                            let bytes_at = allocated_bytes();
                            let live_at = live_blocks();
                            let at = Instant::now();
                            let merged = merge_cohort_serially(
                                &analysed,
                                &slices,
                                MaxCohortLocusSpan::DEFAULT,
                                min_alt_reads,
                            );
                            each_round.push(at.elapsed().as_secs_f64() * 1e3);
                            let allocated_here = allocated_blocks() - blocks_at;
                            blocks += allocated_here;
                            bytes += allocated_bytes() - bytes_at;
                            freed += (allocated_here as i64 + live_at as i64 - live_blocks() as i64)
                                .max(0) as u64;
                            loci += merged.cohort_observations.len();
                        }
                        // **The oracle asked to own what it merges, so the two sides free the same
                        // records inside the same clock.** The plain `oracle` above borrows the
                        // cohort and so frees nothing while it is timed, while both cached drivers
                        // own their copy and free every record of it through eviction. Comparing
                        // those two charges the cache for a deallocation the oracle's caller pays
                        // later and out of shot. This arm gives the oracle its own copy — built
                        // before the clock, exactly as `sources_over` is — and drops it inside.
                        "oracle_owned" => {
                            let owned: Vec<Vec<SampleLocusObservations>> = cohort.clone();
                            let blocks_at = allocated_blocks();
                            let bytes_at = allocated_bytes();
                            let live_at = live_blocks();
                            let at = Instant::now();
                            let merged = {
                                let slices: Vec<&[SampleLocusObservations]> =
                                    owned.iter().map(Vec::as_slice).collect();
                                merge_cohort_serially(
                                    &analysed,
                                    &slices,
                                    MaxCohortLocusSpan::DEFAULT,
                                    min_alt_reads,
                                )
                            };
                            let built = merged.cohort_observations.len();
                            drop(owned);
                            each_round.push(at.elapsed().as_secs_f64() * 1e3);
                            let allocated_here = allocated_blocks() - blocks_at;
                            blocks += allocated_here;
                            bytes += allocated_bytes() - bytes_at;
                            freed += (allocated_here as i64 + live_at as i64 - live_blocks() as i64)
                                .max(0) as u64;
                            loci += built;
                            drop(merged);
                        }
                        "cache" => {
                            let mut cache = ObservationCache::over(sources_over(&cohort, supply));
                            let blocks_at = allocated_blocks();
                            let bytes_at = allocated_bytes();
                            let live_at = live_blocks();
                            let at = Instant::now();
                            let merged = merge_cohort_through_cache(
                                &analysed,
                                &mut cache,
                                width,
                                MaxCohortLocusSpan::DEFAULT,
                                min_alt_reads,
                            )
                            .expect("the probe's sources cannot fail");
                            each_round.push(at.elapsed().as_secs_f64() * 1e3);
                            let allocated_here = allocated_blocks() - blocks_at;
                            blocks += allocated_here;
                            bytes += allocated_bytes() - bytes_at;
                            freed += (allocated_here as i64 + live_at as i64 - live_blocks() as i64)
                                .max(0) as u64;
                            loci += merged.cohort_observations.len();
                        }
                        "parallel" => {
                            let mut cache = ObservationCache::over(sources_over(&cohort, supply));
                            let blocks_at = allocated_blocks();
                            let bytes_at = allocated_bytes();
                            let live_at = live_blocks();
                            let at = Instant::now();
                            let merged = pool
                                .install(|| {
                                    merge_cohort_in_parallel(
                                        &analysed,
                                        &mut cache,
                                        width,
                                        in_flight,
                                        MaxCohortLocusSpan::DEFAULT,
                                        min_alt_reads,
                                    )
                                })
                                .expect("the probe's sources cannot fail");
                            each_round.push(at.elapsed().as_secs_f64() * 1e3);
                            let allocated_here = allocated_blocks() - blocks_at;
                            blocks += allocated_here;
                            bytes += allocated_bytes() - bytes_at;
                            freed += (allocated_here as i64 + live_at as i64 - live_blocks() as i64)
                                .max(0) as u64;
                            loci += merged.cohort_observations.len();
                        }
                        other => {
                            return Err(format!(
                        "NG_REAL_ONLY must be oracle, oracle_owned, cache or parallel, not {other}"
                    )
                    .into());
                        }
                    }
                }
                let seconds = started.elapsed().as_secs_f64();
                if blocks > 0 {
                    println!(
                        "# heap blocks allocated inside the merge, per round: {}",
                        blocks / rounds as u64
                    );
                    println!(
                        "# heap bytes allocated inside the merge, per round: {}",
                        bytes / rounds as u64
                    );
                    println!(
                        "# heap blocks freed inside the merge, per round: {}",
                        freed / rounds as u64
                    );
                }
                each_round.sort_by(f64::total_cmp);
                println!(
                    "{driver}, {rounds}, {bases}, {threads}, {}, {}, {:.2}, {:.2}, {:.2}, \
                 {seconds:.2}",
                    supply.name(),
                    loci / rounds,
                    each_round[each_round.len() / 2],
                    each_round[0],
                    each_round[each_round.len() - 1],
                );
                // Only the parallel driver has a pool; the others run whatever the calling thread
                // gives them, so the "spread perfectly" line must divide by one for those.
                let breakdown =
                    merge_timing::report(if driver == "parallel" { threads } else { 1 });
                if breakdown.merge_wall_ms > 0.0 {
                    println!(
                        "# the merge's own stopwatches on {threads} threads with records {}, \
                 summed over all {rounds} rounds",
                        supply.name(),
                    );
                    println!("{breakdown}");
                }
                println!(
                    "# peak resident after {threads} threads: {}",
                    peak_resident()
                );
            }
        }
        return Ok(());
    }

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
                min_alt_reads,
            ));
        },
    );
    println!("\ndriver, region_bases, median_ms, min_ms, max_ms");
    println!("oracle, -, {median:.2}, {fastest:.2}, {slowest:.2}");

    // The width sweep takes one setting of the record supply, not a list: it is already five
    // widths across four drivers, and the setting the sweep is about is the width.
    let supply =
        std::env::var("NG_REAL_LEASE").map_or(Supply::Minted, |value| Supply::from_name(&value));

    for bases in WIDTHS {
        let width =
            CohortLocusBuilderRegionsLen(std::num::NonZeroU32::new(bases).expect("non-zero"));
        let (median, fastest, slowest) = timed(
            || ObservationCache::over(sources_over(&cohort, supply)),
            |mut cache| {
                std::hint::black_box(
                    &merge_cohort_through_cache(
                        &analysed,
                        &mut cache,
                        width,
                        MaxCohortLocusSpan::DEFAULT,
                        min_alt_reads,
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
                    || ObservationCache::over(sources_over(&cohort, supply)),
                    |mut cache| {
                        std::hint::black_box(
                            &merge_cohort_in_parallel(
                                &analysed,
                                &mut cache,
                                width,
                                in_flight,
                                MaxCohortLocusSpan::DEFAULT,
                                min_alt_reads,
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
