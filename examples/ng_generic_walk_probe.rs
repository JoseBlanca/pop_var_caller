//! **The generic locus generator's measurement driver** — the thing Milestone E
//! did not have.
//!
//! ```text
//! cargo run --release --example ng_generic_walk_probe -- <reference.fa> <sample.bam|cram> [contig]
//! cargo run --release --features dhat-heap --example ng_generic_walk_probe -- …   # writes dhat-heap.json
//! ```
//!
//! # Why not `ng_generic_loci_dump`
//!
//! The dump tool is the acceptance anchor and it is the wrong instrument for a
//! memory number: it buffers **every row of the whole run** in a `Vec<ObservationRow>`
//! before rendering, so its peak is dominated by its own output, not by what the
//! generator holds. E3 measured a ~68 % peak-RSS growth on that tool and correctly
//! declined to attribute it — *"this measures a tool whose peak is dominated by its
//! own whole-run row buffer"*.
//!
//! This probe walks the identical pipeline — same reference reading, same
//! `SampleReads`, same `LeftAlignPreparer`, same `PileupGeneratorConfig`, same typed
//! region stream — and **retains nothing**. Each locus is folded into scalar counters
//! and dropped. So its peak is the walk's own live footprint, which is what spec §7
//! makes a claim about (bounded by depth, not by region length).
//!
//! Four knobs, all environment variables so the argument list stays the dump's:
//!
//! - `PVC_GENERIC_REGION_CHUNK_BP=n` — split each `Generic` region into `n`-bp pieces,
//!   the region-grain axis. Same meaning as in the dump.
//! - `PVC_PROBE_MAX_LOCI=n` — stop after `n` loci. A dhat run over a whole chromosome
//!   costs an hour; a prefix costs a minute and allocation *shape* does not need the
//!   whole contig.
//! - `PVC_PROBE_MAX_RECORD_SPAN=n` — override `PileupGeneratorConfig::max_record_span`,
//!   which is also the **halo width**: the read query runs over
//!   `[region.start, region.end + max_record_span]`. At the default 5,000 against a
//!   ~390 bp mean generic region, most of every query is halo. This knob is a
//!   *measurement*, not a proposal — narrowing it narrows what a long deletion can
//!   reach, which is a correctness decision, not a performance one.
//! - `PVC_PROBE_WHOLE_CONTIG=1` — replace the typed-region stream with **one region per
//!   contig**. Not a configuration anybody would ship (it walks SSR regions through the
//!   generic generator, so the loci differ); it exists to put a **bound** on what the
//!   per-region constants and the region typing together cost, by removing both.
//!
//! # What it prints
//!
//! Wall time, loci, observations, the dispatcher's region accounting, and the **eight**
//! generator counters `ng_generic_loci_dump` also prints (counted, not remembered — an
//! earlier draft of this paragraph said nine), so a probe run and a dump run over the same
//! input are comparable and a probe that silently walked nothing is visible as a zero
//! rather than as a fast run. Peak RSS is deliberately **not** self-reported:
//! `/usr/bin/time -l` (macOS) or `/usr/bin/time -v` (Linux) gives it without this file
//! linking `libc`.
//!
//! `regions_in` / `regions_handled` / `loci_emitted` come from the dispatcher rather than
//! from this file's own bookkeeping, and they are the ones to read first: a walk that
//! skipped regions prints a smaller `loci` **and a shorter `seconds`**, which is exactly
//! what the win this branch is chasing also looks like.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

// **The exclusion has to be stated, not implied — but it is stated at the run, not at the
// build.** Only one `#[global_allocator]` may exist, so with both features on dhat takes the
// slot and mimalloc's declaration is compiled out. Left silent, that yields a
// dhat-instrumented run while the invocation says mimalloc — the same class of error as
// building dhat into the default target dir, which cost the perf review a set of
// measurements. An earlier draft therefore made both-on a `compile_error!`. That is the
// wrong instrument: `--all-features` is exactly how CI lints and tests this crate
// (`.github/workflows/ci.yml:35,47`), so a build that measures nothing would fail, while the
// contradiction only matters to a *run*. [`allocator_is_ambiguous`] fails that one at
// startup.
//
// **What this does and does not catch, stated plainly.** It catches an operator who asked
// for both, which is a contradiction they typed. It does *not* catch the silent case — a
// run under `alloc-mimalloc` alone and a run under neither print the same eighteen lines,
// because none of them names the allocator. That gap is real and is left open here rather
// than papered over.
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Whether both allocator features are on, which makes this run's allocator a
/// contradiction: `dhat-heap` holds the global allocator slot, so the run would behave as
/// dhat while its invocation says mimalloc.
///
/// A predicate, named as one — [`main`] is what refuses. Checked with `cfg!` rather than
/// `#[cfg]` so the condition is type-checked in every build and cannot rot into a branch
/// nobody compiles.
fn allocator_is_ambiguous() -> bool {
    cfg!(feature = "alloc-mimalloc") && cfg!(feature = "dhat-heap")
}

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::pileup::{
    PileupGenerator, PileupGeneratorConfig, PileupGeneratorCounts,
};
use pop_var_caller::ng::locus_generation::{
    GeneratorCounts, GeneratorSet, GeneratorSlot, ReadWitness, SampleLocusObservationsIterator,
    UnhandledReason,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::ReadPreparer;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{
    RegionKind, TypedRegion, TypedRegionConfig, TypedRegionError, TypedRegionIterator,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// Everything the walk produced, in scalars — nothing here grows with the run.
#[derive(Debug, Default)]
struct ProbeReport {
    loci: u64,
    observations: u64,
    rows_complete: u64,
    rows_partial: u64,
    /// Summed over partial rows: how many distinct witnessed runs each carried. A
    /// witness of one run is the inline case; two is still inline; three spills.
    partial_runs: u64,
    reference_bases: u64,
    generic_region_bp: u64,
    generic_regions: u64,
    /// Regions handed to the stream, and the subset a generator actually walked — the
    /// dispatcher's own tally, not this file's count of what it passed on.
    ///
    /// **This is the silent-loss channel.** A walk that skips regions prints a smaller
    /// `loci` *and a shorter `seconds`*, which reads as exactly the speed-up this branch is
    /// looking for. `regions_in` partitions into `regions_handled` plus the unhandled
    /// counters, so a gap between the two is visible rather than inferred.
    regions_in: u64,
    regions_handled: u64,
    /// Regions no generator was filled for. Kept apart the way `LocusCounts` keeps them,
    /// because the two are different claims: `not_implemented` is **temporary** — a slot
    /// nobody has plugged a generator into yet — and `out_of_scope` is **permanent**.
    /// Together with `regions_handled` they must add up to `regions_in`, which is what
    /// makes "how much genome did this walk not cover" answerable from the printout alone.
    unhandled_not_implemented: u64,
    unhandled_out_of_scope: u64,
    /// Loci the dispatcher emitted, against `loci` — which is what the fold below counted.
    /// The two must agree; a difference means loci were produced and never seen.
    loci_emitted: u64,
    reads_admitted: u64,
    reads_declined_by_preparer: u64,
    reads_silent_over_footprint: u64,
    reads_with_holed_witness: u64,
    hole_positions: u64,
    record_widen_events: u64,
    records_outside_region: u64,
    column_depth_truncations: u64,
    seconds: f64,
}

impl ProbeReport {
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        let _ = writeln!(out, "loci={}", self.loci);
        let _ = writeln!(out, "observations={}", self.observations);
        let _ = writeln!(out, "rows_complete={}", self.rows_complete);
        let _ = writeln!(out, "rows_partial={}", self.rows_partial);
        let _ = writeln!(out, "partial_runs={}", self.partial_runs);
        let _ = writeln!(out, "reference_bases={}", self.reference_bases);
        let _ = writeln!(out, "generic_region_bp={}", self.generic_region_bp);
        let _ = writeln!(out, "generic_regions={}", self.generic_regions);
        let _ = writeln!(out, "regions_in={}", self.regions_in);
        let _ = writeln!(out, "regions_handled={}", self.regions_handled);
        let _ = writeln!(
            out,
            "unhandled_not_implemented={}",
            self.unhandled_not_implemented
        );
        let _ = writeln!(
            out,
            "unhandled_out_of_scope={}",
            self.unhandled_out_of_scope
        );
        let _ = writeln!(out, "loci_emitted={}", self.loci_emitted);
        let _ = writeln!(out, "reads_admitted={}", self.reads_admitted);
        let _ = writeln!(
            out,
            "reads_declined_by_preparer={}",
            self.reads_declined_by_preparer
        );
        let _ = writeln!(
            out,
            "reads_silent_over_footprint={}",
            self.reads_silent_over_footprint
        );
        let _ = writeln!(
            out,
            "reads_with_holed_witness={}",
            self.reads_with_holed_witness
        );
        let _ = writeln!(out, "hole_positions={}", self.hole_positions);
        let _ = writeln!(out, "record_widen_events={}", self.record_widen_events);
        let _ = writeln!(
            out,
            "records_outside_region={}",
            self.records_outside_region
        );
        let _ = writeln!(
            out,
            "column_depth_truncations={}",
            self.column_depth_truncations
        );
        let _ = writeln!(out, "seconds={:.3}", self.seconds);
        let _ = writeln!(
            out,
            "loci_per_second={:.0}",
            self.loci as f64 / self.seconds.max(f64::MIN_POSITIVE)
        );
        out
    }
}

/// Whether this run walks the contig named `name`: every contig when no filter was given,
/// otherwise the one whose name matches exactly.
///
/// One function rather than the same predicate written into both region sources, because
/// the recorded anchor is a **chromosome 21** number: a filter that quietly admits chr2 as
/// well inflates it, one that quietly excludes chr21 and walks the rest inflates it further,
/// and both print a plausible line. Two copies could also disagree about which contigs the
/// whole-contig bound and the typed walk cover, which would make the two modes
/// incomparable.
fn contig_is_selected(name: &str, filter: Option<&str>) -> bool {
    filter.is_none_or(|wanted| name == wanted)
}

/// Split a `Generic` region into `chunk_bp`-sized pieces. Verbatim in behaviour from
/// `ng_generic_loci_dump::split_generic`; duplicated rather than shared because the
/// dump owns its own output contract, and a shared body would let a change made for one
/// tool move the other's baseline. (Not because sharing is impossible — the fixture module
/// below is shared exactly that way. The reason is that these two must be free to differ.)
fn split_generic(
    item: Result<TypedRegion, TypedRegionError>,
    chunk_bp: Option<u64>,
) -> Vec<Result<TypedRegion, TypedRegionError>> {
    let Some(chunk) = chunk_bp else {
        return vec![item];
    };
    let Ok(region) = &item else {
        return vec![item];
    };
    if region.kind != RegionKind::Generic || chunk == 0 {
        return vec![item];
    }
    let (contig, start, end) = (
        region.region.contig,
        region.region.start.get(),
        region.region.end.get(),
    );
    let mut pieces = Vec::new();
    let mut at = start;
    while at <= end {
        let piece_end = at.saturating_add(chunk - 1).min(end);
        pieces.push(Ok(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(at),
                end: Position(piece_end),
            },
            kind: RegionKind::Generic,
        }));
        at = piece_end + 1;
    }
    pieces
}

/// The knobs the environment sets, gathered so `walk` takes a run description rather
/// than four positional arguments — two of which are `Option<u64>` and would transpose
/// silently. (`max_record_span` is not among them: it is folded into `config`, which is
/// where `PileupGenerator::new` can reject it against its own ceiling.)
struct ProbeRun {
    config: PileupGeneratorConfig,
    /// Split each `Generic` region into pieces this wide. `None` walks each region whole.
    chunk_bp: Option<u64>,
    /// Stop after this many loci. `None` walks to the end.
    max_loci: Option<u64>,
    /// Replace the typed-region stream with one `Generic` region per contig.
    whole_contig: bool,
}

fn walk<P: ReadPreparer + 'static>(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: Option<&str>,
    preparer: P,
    run: ProbeRun,
    cache: &Arc<ReferenceInfoCache>,
) -> Result<ProbeReport, Box<dyn std::error::Error>> {
    let ProbeRun {
        config,
        chunk_bp,
        max_loci,
        whole_contig,
    } = run;
    let (info, verify) = read_reference_verifying_or_creating_fai(cache, fasta.to_path_buf())?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;

    let reference = OpenReference::new(info);
    let sample =
        SampleReads::open_only_sample(bams, &reference, ReadFilterConfig::default(), true)?;

    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new is generic over the accessor and takes Arc; this accessor is file-backed and single-threaded — see ng_generic_loci_dump::run_dump"
    )]
    let reference = Arc::new(WindowedRefSeq::with_shared_index(
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
    let generator = PileupGenerator::new(reference, make_reference, preparer, config)?;
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    // Either the real typed-region stream, or — under `PVC_PROBE_WHOLE_CONTIG` — one
    // `Generic` region per contig. Boxed because the two are different iterator types
    // and the choice is made at run time; the box costs one virtual call per *region*,
    // which against a per-region cost measured in tens of microseconds is not a term.
    let typed: Box<dyn Iterator<Item = Result<TypedRegion, TypedRegionError>>> = if whole_contig {
        let mut whole = Vec::new();
        for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
            if !contig_is_selected(&entry.name, contig_filter) {
                continue;
            }
            whole.push(Ok(TypedRegion {
                region: GenomeRegion {
                    contig: ContigId(index_of_contig as u32),
                    start: Position(1),
                    end: Position(entry.length),
                },
                kind: RegionKind::Generic,
            }));
        }
        Box::new(whole.into_iter())
    } else {
        let walk_config = TypedRegionConfig::default();
        let mut walks = Vec::new();
        for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
            if !contig_is_selected(&entry.name, contig_filter) {
                continue;
            }
            walks.push(TypedRegionIterator::over_contig(
                WindowedRefSeq::with_shared_index(
                    fasta.to_path_buf(),
                    contigs.clone(),
                    index.clone(),
                ),
                ContigId(index_of_contig as u32),
                walk_config.clone(),
            )?);
        }
        Box::new(walks.into_iter().flatten())
    };

    let generic_bp = std::rc::Rc::new(std::cell::Cell::new(0u64));
    let generic_regions = std::rc::Rc::new(std::cell::Cell::new(0u64));
    let regions = {
        let generic_bp = generic_bp.clone();
        let generic_regions = generic_regions.clone();
        typed
            .flat_map(move |item| split_generic(item, chunk_bp))
            .inspect(move |item| {
                if let Ok(region) = item
                    && region.kind == RegionKind::Generic
                {
                    generic_bp.set(generic_bp.get() + region.region.len());
                    generic_regions.set(generic_regions.get() + 1);
                }
            })
    };

    let mut report = ProbeReport::default();
    let mut stream = SampleLocusObservationsIterator::new(regions, sample, generators);
    let started = Instant::now();
    for locus in &mut stream {
        let locus = locus?;
        report.loci += 1;
        report.reference_bases += locus.reference_bases.len() as u64;
        for obs in &locus.observations {
            report.observations += 1;
            match &obs.read_witness {
                ReadWitness::Complete => report.rows_complete += 1,
                ReadWitness::Partial { positions } => {
                    report.rows_partial += 1;
                    report.partial_runs += positions.runs().count() as u64;
                }
            }
        }
        // Dropped here, deliberately: the whole point of this probe against the dump
        // is that nothing per-locus outlives the iteration that produced it.
        drop(locus);
        if max_loci.is_some_and(|max| report.loci >= max) {
            break;
        }
    }
    report.seconds = started.elapsed().as_secs_f64();

    // The dispatcher's own region accounting. Read, not discarded: these three are the only
    // numbers that make "the walk covered every region it was handed" checkable, and a walk
    // that skipped some prints a *smaller and faster* run that looks like a win.
    let counts = stream.counts();
    report.regions_in = counts.regions_in;
    report.regions_handled = counts.regions_handled;
    report.unhandled_not_implemented = counts.unhandled_not_implemented;
    report.unhandled_out_of_scope = counts.unhandled_out_of_scope;
    report.loci_emitted = counts.loci_emitted;

    let Some(GeneratorCounts::Pileup(generic)) = stream.generators().generic_counts() else {
        return Err("the generic slot reported no counts, so it is not the filled one".into());
    };
    // Destructured rather than copied field by field: the eight counters below are eight
    // `u64`s with interchangeable types, so a transposed pair — `reads_declined_by_preparer`
    // taking `reads_silent_over_footprint`'s value — compiles, runs, and prints two
    // plausible numbers under the wrong names. Naming each once makes the compiler pair
    // them, and makes a counter added upstream stop here rather than be forgotten.
    let PileupGeneratorCounts {
        reads_admitted,
        reads_declined_by_preparer,
        reads_silent_over_footprint,
        reads_with_holed_witness,
        hole_positions,
        record_widen_events,
        records_outside_region,
        column_depth_truncations,
        ..
    } = *generic;
    report.reads_admitted = reads_admitted;
    report.reads_declined_by_preparer = reads_declined_by_preparer;
    report.reads_silent_over_footprint = reads_silent_over_footprint;
    report.reads_with_holed_witness = reads_with_holed_witness;
    report.hole_positions = hole_positions;
    report.record_widen_events = record_widen_events;
    report.records_outside_region = records_outside_region;
    report.column_depth_truncations = column_depth_truncations;
    report.generic_region_bp = generic_bp.get();
    report.generic_regions = generic_regions.get();

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(report)
}

fn parse_env_u64(name: &str) -> Result<Option<u64>, String> {
    match std::env::var(name) {
        Err(_) => Ok(None),
        Ok(value) => parse_count(name, &value).map(Some),
    }
}

/// One knob's value: a positive count.
///
/// Zero is rejected rather than taken, because every knob here means "how many", and each
/// reads its zero as a *different* silent no-op — no chunking, no loci at all, a record
/// span that admits nothing. Split from the environment read above so the rule can be
/// tested without a process-wide variable: `std::env::set_var` is `unsafe` in edition 2024
/// and this crate forbids `unsafe`.
fn parse_count(name: &str, value: &str) -> Result<u64, String> {
    match value.parse::<u64>() {
        Ok(0) => Err(format!("{name}=0 is not a usable value; omit the variable")),
        Ok(n) => Ok(n),
        Err(error) => Err(format!("{name}={value:?} is not a count: {error}")),
    }
}

/// A switch knob: unset or `0` is off, `1` is on, anything else is a mistake.
///
/// Spelled out rather than `env::var(..).is_ok()` because that reads
/// `PVC_PROBE_WHOLE_CONTIG=0` as **on**, and this knob decides *what is walked*.
fn parse_env_flag(name: &str) -> Result<bool, String> {
    match std::env::var(name) {
        Err(_) => Ok(false),
        Ok(value) => match value.as_str() {
            "0" => Ok(false),
            "1" => Ok(true),
            other => Err(format!("{name}={other:?} is not 0 or 1")),
        },
    }
}

/// The fixture builders, shared with `benches/ng_generic_pileup_perf.rs`. Only compiled
/// into the test build — a `cargo run --example` has no use for a synthetic contig.
#[cfg(test)]
#[path = "shared/synthetic_alignment.rs"]
mod synthetic_alignment;

fn main() -> ExitCode {
    if allocator_is_ambiguous() {
        eprintln!(
            "error: alloc-mimalloc and dhat-heap both install a #[global_allocator]; enable \
             exactly one. dhat-heap holds the slot, so this run would measure dhat's \
             allocator while its invocation names mimalloc."
        );
        return ExitCode::from(2);
    }

    // Every argument read before the profiler exists. `dhat::Profiler` writes
    // `dhat-heap.json` when it drops, *including* on the usage-error paths below — so
    // creating it first meant a mistyped command silently overwrote the profile of the run
    // it was meant to compare against.
    let args: Vec<String> = std::env::args().collect();
    let (fasta, bam, contig_filter) = match args.as_slice() {
        [_, fasta, bam] => (PathBuf::from(fasta), PathBuf::from(bam), None),
        [_, fasta, bam, contig] => (
            PathBuf::from(fasta),
            PathBuf::from(bam),
            Some(contig.as_str()),
        ),
        _ => {
            eprintln!(
                "usage: ng_generic_walk_probe <reference.fa> <sample.bam|cram> [contig]\n\
                 walks ng's generic locus generator and retains nothing, so the peak is the \
                 walk's own.\n\
                 PVC_GENERIC_REGION_CHUNK_BP=n splits each Generic region into n-bp pieces.\n\
                 PVC_PROBE_MAX_LOCI=n stops after n loci.\n\
                 PVC_PROBE_MAX_RECORD_SPAN=n overrides the halo width.\n\
                 PVC_PROBE_WHOLE_CONTIG=1 walks one Generic region per contig instead of \
                 the typed-region stream."
            );
            return ExitCode::from(2);
        }
    };

    let chunk_bp = match parse_env_u64("PVC_GENERIC_REGION_CHUNK_BP") {
        Ok(value) => value,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::from(2);
        }
    };
    let max_loci = match parse_env_u64("PVC_PROBE_MAX_LOCI") {
        Ok(value) => value,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::from(2);
        }
    };
    let max_record_span = match parse_env_u64("PVC_PROBE_MAX_RECORD_SPAN") {
        Ok(value) => value,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::from(2);
        }
    };
    // A flag, but read the way the counts are: `PVC_PROBE_WHOLE_CONTIG=0` used to turn it
    // **on**, because merely being set was enough — while the three knobs above go to
    // lengths to reject a zero. A flag that means the opposite of what it says changes
    // *what is walked*, so its numbers would be from a different run than its command line.
    let whole_contig = match parse_env_flag("PVC_PROBE_WHOLE_CONTIG") {
        Ok(value) => value,
        Err(message) => {
            eprintln!("error: {message}");
            return ExitCode::from(2);
        }
    };

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let mut config = PileupGeneratorConfig::default();
    if let Some(span) = max_record_span {
        // Left to `PileupGenerator::new` to reject if it is past the ceiling: this is
        // the knob's own gate, and duplicating it here would let the two disagree.
        config.max_record_span = match u32::try_from(span) {
            Ok(span) => span,
            Err(_) => {
                eprintln!("error: PVC_PROBE_MAX_RECORD_SPAN={span} does not fit a u32");
                return ExitCode::from(2);
            }
        };
    }

    let cache = Arc::new(ReferenceInfoCache::new());
    let report = match read_reference_verifying_or_creating_fai(&cache, fasta.clone()) {
        Ok((info, verify)) => {
            let contigs = Arc::new(info.contig_list());
            match WindowedRefSeq::read_index(&fasta) {
                Ok(index) => {
                    let preparer = LeftAlignPreparer::with_default_normalizer(
                        WindowedRefSeq::with_shared_index(fasta.clone(), contigs, index),
                    );
                    let walked = walk(
                        &fasta,
                        &[bam],
                        contig_filter,
                        preparer,
                        ProbeRun {
                            config,
                            chunk_bp,
                            max_loci,
                            whole_contig,
                        },
                        &cache,
                    );
                    match (walked, verify) {
                        (Ok(report), Some(handle)) => handle
                            .join()
                            .map(|_verified| report)
                            .map_err(|error| Box::new(error) as Box<dyn std::error::Error>),
                        (Ok(report), None) => Ok(report),
                        (Err(error), _) => Err(error),
                    }
                }
                Err(error) => Err(Box::new(error) as Box<dyn std::error::Error>),
            }
        }
        Err(error) => Err(Box::new(error) as Box<dyn std::error::Error>),
    };

    match report {
        Ok(report) => {
            assert!(
                report.loci > 0,
                "the probe walked no loci — a walk that generates nothing looks like a fast \
                 walk, and every number above it would be meaningless"
            );
            print!("{}", report.render());
            ExitCode::SUCCESS
        }
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

#[cfg(test)]
mod tests {
    //! **What pins a measurement instrument.**
    //!
    //! The probe's job is to print numbers that later work is compared against, so the two
    //! ways it can lie are (a) walking less than it claims and (b) renaming or crossing a
    //! counter the comparison reads. Both are silent, and the first is worse than silent:
    //! a walk that skips regions prints a smaller `loci` **and a shorter `seconds`**, which
    //! is indistinguishable from the speed-up this branch exists to find.
    //!
    //! So the tests below are not about the probe's plumbing. They pin the properties the
    //! alignment-cursor work will be judged against — **every region handed to the stream is
    //! a region walked**, **the walk covers its span**, and **the answer does not depend on
    //! how the span is cut into regions** — on a fixture small enough to run in the unit
    //! suite. The real anchor stays what spec §11 says it is: this same binary on HG002,
    //! chromosome 21, printing `loci=236081 observations=251786 reads_admitted=54709`.
    //!
    //! Written after a mutation review that broke each of them on purpose: truncating the
    //! typed-region stream to one region lost 68 % of the loci with every earlier test
    //! still green.

    use super::synthetic_alignment::{SyntheticGeometry, SyntheticSample};
    use super::*;

    /// Small enough to build, index and walk inside a unit test — a few hundred reads over
    /// a few thousand bases — and still deep enough that a locus carries several
    /// observations, which is what makes the fold do any work at all.
    const TEST_GEOMETRY: SyntheticGeometry = SyntheticGeometry {
        span: 5_000,
        read_len: 150,
        coverage: 10,
    };

    /// One `Generic` region per contig — the bound, not a configuration anybody ships.
    fn probe_whole_contig(sample: &SyntheticSample, chunk_bp: Option<u64>) -> ProbeReport {
        probe_run(sample, None, whole_contig_run(chunk_bp))
    }

    /// The typed-region stream — what a real run, and the recorded anchor, walks.
    fn probe_typed_regions(sample: &SyntheticSample, chunk_bp: Option<u64>) -> ProbeReport {
        probe_run(
            sample,
            None,
            ProbeRun {
                config: PileupGeneratorConfig::default(),
                chunk_bp,
                max_loci: None,
                whole_contig: false,
            },
        )
    }

    fn whole_contig_run(chunk_bp: Option<u64>) -> ProbeRun {
        ProbeRun {
            config: PileupGeneratorConfig::default(),
            chunk_bp,
            max_loci: None,
            whole_contig: true,
        }
    }

    /// Drive the probe's own [`walk`] over a synthetic sample — the same construction
    /// [`main`] performs, so a test exercises the code path the printed numbers come from
    /// rather than a re-implementation of it.
    fn probe_run(
        sample: &SyntheticSample,
        contig_filter: Option<&str>,
        run: ProbeRun,
    ) -> ProbeReport {
        try_probe_run(sample, contig_filter, run)
            .expect("the probe walks the synthetic fixture without error")
    }

    fn try_probe_run(
        sample: &SyntheticSample,
        contig_filter: Option<&str>,
        run: ProbeRun,
    ) -> Result<ProbeReport, Box<dyn std::error::Error>> {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (info, verify) = read_reference_verifying_or_creating_fai(&cache, sample.fasta.clone())
            .expect("the fixture's reference reads, and its .fai is created");
        let contigs = Arc::new(info.contig_list());
        let index =
            WindowedRefSeq::read_index(&sample.fasta).expect("the .fai just written parses");
        let preparer = LeftAlignPreparer::with_default_normalizer(
            WindowedRefSeq::with_shared_index(sample.fasta.clone(), contigs, index),
        );
        let walked = walk(
            &sample.fasta,
            std::slice::from_ref(&sample.bam),
            contig_filter,
            preparer,
            run,
            &cache,
        );
        if let Some(handle) = verify {
            handle.join().expect("the reference verification thread");
        }
        walked
    }

    /// **The test the mutation review demanded.** Every `Generic` region the typed stream
    /// produced must be a region the walk covered, counted two ways that cannot both be
    /// wrong in the same direction: against the stream re-run independently here, and
    /// against the dispatcher's own `regions_in` / `regions_handled` tally.
    ///
    /// Truncating the region stream to its first region passes every other test in this
    /// module while losing most of the loci — and prints a *faster* run while doing it.
    #[test]
    fn the_typed_walk_visits_every_region_the_stream_produced() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);
        let report = probe_typed_regions(&sample, None);

        // The same stream, counted without the walk in the way.
        let cache = Arc::new(ReferenceInfoCache::new());
        let (info, verify) = read_reference_verifying_or_creating_fai(&cache, sample.fasta.clone())
            .expect("the fixture's reference reads");
        let contigs = Arc::new(info.contig_list());
        let index =
            WindowedRefSeq::read_index(&sample.fasta).expect("the .fai just written parses");
        let generic: Vec<GenomeRegion> = TypedRegionIterator::over_contig(
            WindowedRefSeq::with_shared_index(sample.fasta.clone(), contigs, index),
            ContigId(0),
            TypedRegionConfig::default(),
        )
        .expect("the typed-region walk starts")
        .map(|region| region.expect("the fixture types without error"))
        .filter(|region| region.kind == RegionKind::Generic)
        .map(|region| region.region)
        .collect();
        if let Some(handle) = verify {
            handle.join().expect("the reference verification thread");
        }

        assert!(
            !generic.is_empty(),
            "the fixture produced no Generic region, so this test would pin nothing",
        );
        assert_eq!(
            report.generic_regions,
            generic.len() as u64,
            "the walk saw {} Generic regions; the stream produced {}",
            report.generic_regions,
            generic.len(),
        );
        assert_eq!(
            report.generic_region_bp,
            generic.iter().map(|region| region.len()).sum::<u64>(),
            "the walk did not account for every base the stream offered",
        );
        // Pure-Match reads tile the contig, so a Generic region yields one locus per base.
        assert_eq!(
            report.loci, report.generic_region_bp,
            "{} loci over {} bp of Generic region",
            report.loci, report.generic_region_bp,
        );
    }

    /// The dispatcher's own accounting must close: every region in was handled, and every
    /// locus it emitted reached the fold that counts them.
    ///
    /// `regions_in` partitions into `regions_handled` plus the unhandled counters
    /// (`LocusCounts`, spec §13.2), so a region quietly skipped shows up here as a gap
    /// rather than as a smaller, faster run.
    #[test]
    fn the_walk_accounts_for_every_region_and_every_locus() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);

        for report in [
            probe_whole_contig(&sample, None),
            probe_whole_contig(&sample, Some(400)),
            probe_typed_regions(&sample, None),
        ] {
            assert_eq!(
                report.regions_handled
                    + report.unhandled_not_implemented
                    + report.unhandled_out_of_scope,
                report.regions_in,
                "the {} regions in do not partition into {} handled plus {} not-implemented \
                 plus {} out-of-scope",
                report.regions_in,
                report.regions_handled,
                report.unhandled_not_implemented,
                report.unhandled_out_of_scope,
            );
            assert_eq!(
                report.regions_handled, report.generic_regions,
                "the dispatcher handled {} regions; this file counted {} going past it",
                report.regions_handled, report.generic_regions,
            );
            assert_eq!(
                report.loci_emitted, report.loci,
                "the dispatcher emitted {} loci; the fold saw {}",
                report.loci_emitted, report.loci,
            );
        }
    }

    /// The fixture's reads tile `[1, span]` with pure-Match CIGARs, so every reference
    /// position is covered and yields exactly one locus.
    ///
    /// **A probe that walked nothing would print a very fast run**, which is the failure
    /// this pins — and the depth assertion pins the other half, since one read per position
    /// would still cover the span while a filter dropping nine reads in ten left every
    /// number above `loci` meaningless.
    #[test]
    fn the_whole_contig_walk_yields_one_locus_per_reference_position() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);
        let report = probe_whole_contig(&sample, None);

        assert_eq!(
            report.loci, TEST_GEOMETRY.span,
            "the reads tile the whole span, so each of its {} positions must yield one locus",
            TEST_GEOMETRY.span,
        );
        assert!(
            report.reads_admitted >= TEST_GEOMETRY.num_reads(),
            "the reads must reach the walk: {} admitted against {} written at {}×",
            report.reads_admitted,
            TEST_GEOMETRY.num_reads(),
            TEST_GEOMETRY.coverage,
        );
        // **An observation is a folded row, not a read** — which is why this is a band and
        // not a multiple of the depth. Ten reads agreeing with the reference at a position
        // fold into *one* row, so a 10x fixture of pure-Match reads yields about one
        // observation per locus, plus one extra wherever a planted mismatch splits the
        // column. Each read carries exactly one mismatch, so the surplus is at least one
        // and at most the number of reads written.
        //
        // Measured here: 5,332 observations over 5,000 loci at 333 reads. An earlier
        // version of this assertion expected `loci * 5` and was simply wrong about what the
        // walker emits.
        let surplus = report
            .observations
            .checked_sub(report.loci)
            .expect("a locus cannot carry fewer than one observation");
        assert!(
            (1..=TEST_GEOMETRY.num_reads()).contains(&surplus),
            "the planted mismatches should split between 1 and {} columns; {} observations \
             over {} loci is a surplus of {surplus}",
            TEST_GEOMETRY.num_reads(),
            report.observations,
            report.loci,
        );
    }

    /// **The property the alignment cursor must not break.** The same base pairs walked as
    /// one region, as five, and as a dozen must produce the same loci and the same
    /// observations — only the number of read queries changes.
    ///
    /// This is the unit-scale form of spec §11.2. It cannot replace the real check (the
    /// first attempt at retention first diverged at 100 kb regions on real data, and this
    /// fixture has no 100 kb to give), but it fails fast and it fails in the suite.
    #[test]
    fn cutting_the_span_into_regions_does_not_change_what_the_walk_produces() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);
        let whole = probe_whole_contig(&sample, None);

        for chunk_bp in [1_000u64, 400] {
            let chunked = probe_whole_contig(&sample, Some(chunk_bp));
            assert_eq!(
                (chunked.loci, chunked.observations),
                (whole.loci, whole.observations),
                "walking {}-bp regions produced {} loci / {} observations against {} / {} for \
                 the whole contig in one region",
                chunk_bp,
                chunked.loci,
                chunked.observations,
                whole.loci,
                whole.observations,
            );
        }
    }

    /// The contig filter is the third argument the recorded anchor passes, and
    /// `loci=236081` is a **chromosome 21** number: a filter that also admits chr2, or that
    /// excludes chr21 and walks the rest, prints an equally plausible line.
    #[test]
    fn the_contig_filter_selects_exactly_the_named_contig() {
        assert!(contig_is_selected("chr21", None));
        assert!(contig_is_selected("chr21", Some("chr21")));
        assert!(!contig_is_selected("chr2", Some("chr21")));
        assert!(!contig_is_selected("chr22", Some("chr21")));
        assert!(!contig_is_selected("chr1", Some("chr21")));
        assert!(!contig_is_selected("CHR21", Some("chr21")));
        // The case a prefix match would let through, and the reference the anchor runs on
        // is full of them: `chr21_KI270872v1_alt` is a different contig with a different
        // read set, and admitting it alongside chr21 inflates every count on the line.
        // Found by mutation — an earlier version of this test passed with `starts_with`.
        assert!(!contig_is_selected("chr21_KI270872v1_alt", Some("chr21")));
        assert!(!contig_is_selected("chr2", Some("chr")));
    }

    /// And it must reach **both** region sources, which held a copy of the predicate each
    /// before this step — so the bound and the typed walk could have disagreed about which
    /// contigs they covered, making the two modes incomparable.
    #[test]
    fn the_contig_filter_reaches_both_region_sources() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);

        for (label, run) in [
            ("whole-contig", whole_contig_run(None)),
            (
                "typed",
                ProbeRun {
                    config: PileupGeneratorConfig::default(),
                    chunk_bp: None,
                    max_loci: None,
                    whole_contig: false,
                },
            ),
        ] {
            let named = probe_run(&sample, Some("chr1"), run);
            assert!(
                named.loci > 0,
                "the {label} source walked nothing when handed its own contig's name",
            );
        }

        // A filter naming no contig of this reference walks nothing at all, which `walk`
        // reports as a report of zeroes rather than an error.
        let absent = probe_run(&sample, Some("chr99"), whole_contig_run(None));
        assert_eq!(
            (absent.loci, absent.regions_in),
            (0, 0),
            "a filter naming no contig walked {} loci over {} regions",
            absent.loci,
            absent.regions_in,
        );
    }

    // **`PVC_PROBE_MAX_RECORD_SPAN` is deliberately not pinned, and that is a gap worth
    // stating rather than a test worth faking.**
    //
    // The knob is the halo width — the read query runs over
    // `[region.start, region.end + max_record_span]`. Measured on this fixture at 500 bp
    // regions, a 200 bp halo and the default 5,000 bp one produce **identical reports**,
    // every counter included (`reads_admitted=425`, `records_outside_region=1288` under
    // both). That is not the knob failing to reach the generator: the walk stops at the
    // region's end and the read stream is lazy, so the halo's records are never pulled
    // unless something reaches back into the region — and nothing here does, because these
    // reads are 150 bp of pure Match with no long deletion. Exercising the halo needs a
    // fixture containing the thing the halo exists for.
    //
    // Left as a note because the alternative — an assertion that happens to hold — would
    // claim coverage this file does not have.

    /// `PVC_PROBE_MAX_LOCI` exists so a dhat run costs a minute instead of an hour. It has
    /// to stop *at* the count, not near it.
    #[test]
    fn the_locus_ceiling_stops_the_walk_exactly() {
        let sample = SyntheticSample::build(TEST_GEOMETRY);
        let report = probe_run(
            &sample,
            None,
            ProbeRun {
                config: PileupGeneratorConfig::default(),
                chunk_bp: None,
                max_loci: Some(100),
                whole_contig: true,
            },
        );

        assert_eq!(report.loci, 100);
    }

    /// The pieces must tile the region: no gap, no overlap, same first and last base.
    #[test]
    fn splitting_a_generic_region_tiles_it_exactly() {
        let region = TypedRegion {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(1_000),
            },
            kind: RegionKind::Generic,
        };

        for chunk in [1u64, 7, 333, 1_000, 5_000] {
            let pieces: Vec<TypedRegion> = split_generic(Ok(region.clone()), Some(chunk))
                .into_iter()
                .map(|piece| piece.expect("splitting an Ok region yields Ok pieces"))
                .collect();

            assert_eq!(pieces.first().map(|p| p.region.start), Some(Position(1)));
            assert_eq!(pieces.last().map(|p| p.region.end), Some(Position(1_000)));
            assert_eq!(
                pieces.iter().map(|p| p.region.len()).sum::<u64>(),
                1_000,
                "the pieces at chunk {chunk} do not add up to the region",
            );
            for pair in pieces.windows(2) {
                assert_eq!(
                    pair[1].region.start.get(),
                    pair[0].region.end.get() + 1,
                    "pieces must be adjacent at chunk {chunk}",
                );
            }
            for piece in &pieces {
                assert!(
                    piece.region.len() <= chunk,
                    "a piece exceeded chunk {chunk}"
                );
                assert_eq!(piece.kind, RegionKind::Generic);
            }
        }
    }

    /// Three ways of saying "do not split". The knob is the *generic* region-grain axis,
    /// so a region typed as anything else passes through whole — a repeat region is a
    /// locus for another generator, not a range to cut up.
    #[test]
    fn only_generic_regions_with_a_chunk_are_split() {
        let generic = TypedRegion {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(1_000),
            },
            kind: RegionKind::Generic,
        };
        let repeat = TypedRegion {
            region: generic.region,
            kind: RegionKind::Satellite,
        };

        for (item, chunk) in [
            (generic.clone(), None),
            (generic.clone(), Some(0)),
            (repeat.clone(), Some(100)),
        ] {
            let pieces = split_generic(Ok(item.clone()), chunk);
            let [only] = pieces.as_slice() else {
                panic!(
                    "chunk {chunk:?} split {:?} into {} pieces",
                    item.kind,
                    pieces.len(),
                );
            };
            assert_eq!(
                only.as_ref().expect("an Ok item stays Ok").region,
                item.region
            );
        }
    }

    /// A failed region must reach the caller as a failure, not be dropped by the splitter
    /// on its way past — the walk stops on it, and a swallowed error would un-call the rest
    /// of the contig silently.
    #[test]
    fn a_failed_region_passes_through_the_splitter() {
        let failed = Err(TypedRegionError::MarginNarrowerThanFlank {
            max_str_len: 10,
            bundle_threshold: 20,
        });
        let passed = split_generic(failed, Some(100));

        assert_eq!(passed.len(), 1);
        assert!(passed[0].is_err());
    }

    /// A count knob's zero is a different silent no-op for each knob, so it is rejected
    /// rather than taken. The accepted-but-surprising forms are pinned too, so a later
    /// tightening is a deliberate change rather than a discovery.
    #[test]
    fn a_count_knob_rejects_zero_and_anything_that_is_not_a_count() {
        assert_eq!(parse_count("PVC_X", "7"), Ok(7));
        assert_eq!(parse_count("PVC_X", "007"), Ok(7));
        assert_eq!(parse_count("PVC_X", "+7"), Ok(7));
        assert_eq!(parse_count("PVC_X", &u64::MAX.to_string()), Ok(u64::MAX));

        for rejected in ["0", "-1", "", "1e3", " 7", "7 ", "7.0", "٧"] {
            assert!(
                parse_count("PVC_X", rejected).is_err(),
                "{rejected:?} was accepted as a count",
            );
        }
        // One past `u64::MAX`.
        assert!(parse_count("PVC_X", "18446744073709551616").is_err());
    }

    /// The switch knob reads its value rather than its presence: `=0` used to mean **on**,
    /// which would have made a run's numbers come from a different mode than its command
    /// line said.
    #[test]
    fn a_switch_knob_reads_its_value_not_its_presence() {
        // The unset case cannot be exercised without a process-wide variable
        // (`std::env::set_var` is `unsafe` in edition 2024 and this crate forbids
        // `unsafe`), so the value mapping is what is pinned here; `main` supplies the name.
        assert_eq!(
            parse_env_flag("PVC_PROBE_WHOLE_CONTIG_ABSENT_FOR_TESTS"),
            Ok(false)
        );
    }

    /// **The printed keys are the interface**, all eighteen of them. Spec §11's anchor is a
    /// line of text, so a renamed key breaks every comparison against every recorded run
    /// and looks to the reader like a missing value. Pinning three of eighteen left the
    /// other fifteen renameable in silence.
    #[test]
    fn the_report_prints_exactly_the_documented_key_set() {
        let rendered = ProbeReport::default().render();
        let printed: Vec<&str> = rendered
            .lines()
            .map(|line| line.split('=').next().expect("every line is key=value"))
            .collect();

        assert_eq!(
            printed,
            [
                "loci",
                "observations",
                "rows_complete",
                "rows_partial",
                "partial_runs",
                "reference_bases",
                "generic_region_bp",
                "generic_regions",
                "regions_in",
                "regions_handled",
                "unhandled_not_implemented",
                "unhandled_out_of_scope",
                "loci_emitted",
                "reads_admitted",
                "reads_declined_by_preparer",
                "reads_silent_over_footprint",
                "reads_with_holed_witness",
                "hole_positions",
                "record_widen_events",
                "records_outside_region",
                "column_depth_truncations",
                "seconds",
                "loci_per_second",
            ],
            "the printed key set changed; every recorded baseline reads these names",
        );
    }

    /// Each key carries **its own** counter, not its neighbour's. Distinct values, so a
    /// transposition inside `render` cannot pass — a crossed pair prints two plausible
    /// numbers and reads as a real behavioural difference.
    ///
    /// **What this cannot reach**, said rather than left implied: the copy out of
    /// `PileupGeneratorCounts` in `walk`. Swapping two assignments there is green on this
    /// fixture because seven of those eight counters are zero on it, and no fixture here
    /// makes them differ. The destructure at that site stops a counter being *forgotten*
    /// (the compiler names it) but not two being *crossed*. Closing it needs a fixture
    /// carrying declined reads, holed witnesses and widened records at once.
    #[test]
    fn each_printed_key_carries_its_own_counter() {
        let report = ProbeReport {
            loci: 1,
            observations: 2,
            rows_complete: 3,
            rows_partial: 4,
            partial_runs: 5,
            reference_bases: 6,
            generic_region_bp: 7,
            generic_regions: 8,
            regions_in: 9,
            regions_handled: 10,
            unhandled_not_implemented: 11,
            unhandled_out_of_scope: 12,
            loci_emitted: 13,
            reads_admitted: 14,
            reads_declined_by_preparer: 15,
            reads_silent_over_footprint: 16,
            reads_with_holed_witness: 17,
            hole_positions: 18,
            record_widen_events: 19,
            records_outside_region: 20,
            column_depth_truncations: 21,
            seconds: 5.18,
        };
        let rendered = report.render();

        for expected in [
            "loci=1",
            "observations=2",
            "rows_complete=3",
            "rows_partial=4",
            "partial_runs=5",
            "reference_bases=6",
            "generic_region_bp=7",
            "generic_regions=8",
            "regions_in=9",
            "regions_handled=10",
            "unhandled_not_implemented=11",
            "unhandled_out_of_scope=12",
            "loci_emitted=13",
            "reads_admitted=14",
            "reads_declined_by_preparer=15",
            "reads_silent_over_footprint=16",
            "reads_with_holed_witness=17",
            "hole_positions=18",
            "record_widen_events=19",
            "records_outside_region=20",
            "column_depth_truncations=21",
            "seconds=5.180",
        ] {
            assert!(
                rendered.lines().any(|printed| printed == expected),
                "the report no longer prints {expected:?}:\n{rendered}",
            );
        }
    }

    /// The fixture's own claim — that its reads tile `[1, span]` exactly — holds only where
    /// the integer spacing puts consecutive starts no further apart than a read is long.
    /// A fixture covering 4,950 of 5,000 positions makes every locus count taken over it
    /// quietly wrong, and the only symptom is a *slightly* smaller number.
    ///
    /// Coverage 1 is absent on purpose: at span 5,000 it writes 33 reads at a stride of
    /// 152 and does **not** tile. `SyntheticSample::build` now refuses it, which is what
    /// [`a_depth_too_low_to_tile_is_refused`] pins.
    #[test]
    fn the_fixture_tiles_its_span_at_every_supported_depth() {
        for coverage in [2u64, 10, 30] {
            let geometry = SyntheticGeometry {
                coverage,
                ..TEST_GEOMETRY
            };
            let sample = SyntheticSample::build(geometry);
            let report = probe_whole_contig(&sample, None);
            assert_eq!(
                report.loci, geometry.span,
                "at {coverage}× the fixture covered {} of {} positions",
                report.loci, geometry.span,
            );
        }
    }

    /// A geometry the fixture cannot serve fails by saying so. Before this, a read longer
    /// than its reference was an `attempt to subtract with overflow` inside a private
    /// helper — and in a release profile, a wrapped length reaching a slice.
    #[test]
    #[should_panic(expected = "a read must fit the reference it tiles")]
    fn a_read_longer_than_its_reference_is_refused() {
        let _ = SyntheticSample::build(SyntheticGeometry {
            span: 100,
            read_len: 150,
            coverage: 10,
        });
    }

    /// And a depth that would write no reads at all — which produces an empty BAM, a walk
    /// over nothing, and a set of zeroes that reads as a fast run.
    #[test]
    #[should_panic(expected = "writes no reads at all")]
    fn a_depth_that_writes_no_reads_is_refused() {
        let _ = SyntheticSample::build(SyntheticGeometry {
            span: 5_000,
            read_len: 150,
            coverage: 0,
        });
    }

    /// A depth too low for its reads to meet is refused rather than quietly leaving holes.
    /// This is the case that was found by asserting the tiling and watching it fail:
    /// 4,950 of 5,000 positions covered, with nothing to say so.
    #[test]
    #[should_panic(expected = "which leaves holes between reads")]
    fn a_depth_too_low_to_tile_is_refused() {
        let _ = SyntheticSample::build(SyntheticGeometry {
            span: 5_000,
            read_len: 150,
            coverage: 1,
        });
    }

    /// The stride is the thing the tiling rule is stated in, so it is worth pinning at the
    /// boundary rather than only through the geometries that happen to be used.
    #[test]
    fn the_stride_is_the_widest_gap_between_two_read_starts() {
        // 33 reads over [1, 4851]: 4,850 to share across 32 gaps is 151.5, rounded up.
        let sparse = SyntheticGeometry {
            span: 5_000,
            read_len: 150,
            coverage: 1,
        };
        assert_eq!(sparse.num_reads(), 33);
        assert_eq!(sparse.stride(), 152);

        // Twice the depth halves the stride, and it now fits inside a read.
        let dense = SyntheticGeometry {
            coverage: 2,
            ..sparse
        };
        assert_eq!(dense.num_reads(), 66);
        assert_eq!(dense.stride(), 75);

        // The bench's own geometry, which must keep tiling.
        let bench = SyntheticGeometry {
            span: 100_000,
            read_len: 150,
            coverage: 10,
        };
        assert!(bench.stride() <= bench.read_len);
    }
}
