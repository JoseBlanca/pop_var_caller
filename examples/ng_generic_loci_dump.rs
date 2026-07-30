//! **The generic locus generator's dump tool** — `locus_generation_pileup.md` §13, the
//! acceptance anchor. Runs the real pipeline over a reference + a sample's reads and prints,
//! per covered reference position, what that sample's reads showed there: the observed
//! sequences, how much of the locus each read witnessed, which read group it came from, and
//! its support. This generator emits no calls, so "done" must not decay into "compiles" — the
//! output is inspectable and, on committed fixtures, **asserted**.
//!
//! ```text
//! ng_generic_loci_dump <reference.fa> <sample.bam|cram> [contig]
//! ```
//!
//! The reference's sibling `<reference.fa>.fai` is used (created if absent); reads go through
//! ng's real ingestion (`SampleReads`, the filtered reads a caller sees) and ng's real read
//! preparation (`LeftAlignPreparer`). An optional third argument restricts the walk to one
//! contig by name. `PVC_GENERIC_REGION_CHUNK_BP=n` splits each `Generic` region into `n`-bp
//! pieces before generating — the knob the region-boundary check uses, and a real one: in the
//! pipeline the regions are handed in, and their size is the caller's.
//!
//! Output: a `#`-prefixed `key=value` counts header, a bare TSV column line, then one
//! tab-separated row per observed sequence.
//!
//! # What the columns say that a `PileupRecord` could not
//!
//! - **`read_coverage`** is `complete` or `observed:<offset>+<positions>`. A partial witness
//!   is a *censored* observation — the sequence is at least this long — and production has no
//!   way to say it: it fills the positions the read did not witness from the reference and
//!   folds the result as if the read had seen them (spec §4, §6). This column is the whole
//!   reason the generic generator exists.
//! - **`read_group`** joins the row identity, so an allele carried by two lanes is two rows.
//! - **`ref_bases`** is the whole footprint, and **`end`** its last position: a one-base SNP
//!   locus and a widened deletion locus are told apart by the region, not by guessing.
//!
//! # The counts header, and why every line of it is load-bearing
//!
//! Every read the walk saw has to be somewhere. A read is **admitted** and then either
//! contributes at some position (and so appears in a locus's rows, or loses its row to a
//! non-contiguous witness, or is dropped by a depth cap) or contributes **nowhere** — silent
//! over every footprint it touched, every base `N` or adaptor-masked. That last class is
//! invisible to both per-locus counters and is `reads_silent_over_footprint`. A read the
//! **preparer** declines never reaches the walk at all, and is `reads_declined_by_preparer`.
//!
//! **The header prints them on two lines because they are in two units, and this doc used to
//! promise a sum that cannot exist.** It said printing all of them "makes 'every fetched read
//! is accounted for' a thing a reader can check rather than a claim (spec §13)". It does not:
//! `reads_admitted`, `reads_declined_by_preparer` and `reads_silent_over_footprint` are once
//! per read for the whole run, while `reads_without_observation` and `reads_discarded_by_cap`
//! are **per-locus** quantities summed over loci — a read is counted once per locus it was
//! affected at, and a read can be in a row at one locus and unobserved at another. On the
//! four-read accounting fixture with both depth caps at 1 the header reads `reads_admitted=4`
//! beside `reads_discarded_by_cap=55`. The per-locus pair therefore carries a `locus_sum_`
//! prefix on its own line.
//!
//! **What is actually checkable, and is checked:** `reads_complete + reads_observed` equals the
//! summed `reads` column over every row, asserted below. The identity spec §13 asks for — every
//! *fetched* read accounted for exactly once — needs a run-level "appeared in some row" flag
//! beside `ever_contributed`, which the walk does not carry; until it does, the honest statement
//! is the one above rather than a sum a reader would fail to reproduce.

use std::cell::Cell;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::rc::Rc;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorCounts, GeneratorSet, GeneratorSlot, LocusKind, ReadCoverage, SampleLocusObservations,
    SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::ReadPreparer;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{
    RegionKind, TypedRegion, TypedRegionConfig, TypedRegionError, TypedRegionIterator,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// One TSV row: an observed sequence at a locus, with its support.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ObservationRow {
    contig: String,
    /// 1-based locus coordinates (inclusive) — the record's whole footprint.
    start: u64,
    end: u64,
    ref_bases: Vec<u8>,
    /// The locus's **complete**-observation depth, shown on every row of the locus so a
    /// partial row is read against the depth that actually pinned it.
    depth: u32,
    read_coverage: String,
    read_group: u32,
    observed: Vec<u8>,
    reads: u32,
    chain_ids: Vec<u64>,
}

/// The whole dump: the run-level counts plus the rows. Separated from rendering so the tests
/// assert on numbers rather than by parsing text.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct DumpReport {
    /// Loci emitted — one per covered position, after the region clamp.
    generic_loci: u64,
    /// Reference bases inside the `Generic` regions the walk was given. The denominator for
    /// "uncovered positions yield no locus": `generic_loci` must be **strictly less** on any
    /// fixture with an uncovered stretch.
    generic_region_bp: u64,
    /// `Generic` regions handed to the generator (after any chunking).
    generic_regions: u64,
    /// Rows by coverage class, and the reads behind them. Both, because "rows" is about the
    /// emitted shape and "reads" about the evidence, and the read-group split moves one
    /// without moving the other.
    rows_complete: u64,
    rows_observed: u64,
    reads_complete: u64,
    reads_observed: u64,
    /// Summed over loci: reads that covered a locus and produced nothing (A5), and reads a
    /// depth cap discarded (B3).
    reads_without_observation: u64,
    reads_discarded_by_cap: u64,
    /// The generator's own run-level counts (`PileupGeneratorCounts`), read through the boxed
    /// generator — so this tool is also the check that they are reachable at all.
    reads_admitted: u64,
    reads_declined_by_preparer: u64,
    reads_silent_over_footprint: u64,
    record_widen_events: u64,
    records_outside_region: u64,
    column_depth_truncations: u64,
    /// The dispatcher's shared tally.
    regions_in: u64,
    regions_handled: u64,
    loci_emitted: u64,
    rows: Vec<ObservationRow>,
}

impl DumpReport {
    /// Fold one locus into the report, emitting a row per observation.
    ///
    /// # The one invariant asserted here rather than in a test
    ///
    /// **No row claims a position the locus does not have, and no row claims zero** (spec §13):
    /// every `Observed` row satisfies `offset_in_locus + positions_covered <= footprint` and
    /// `positions_covered > 0`. Asserted on every locus of every run, including the tool's real
    /// ones, because a run over new data is where it would first be violated.
    ///
    /// **This is a bound on positions, and it does not read `bases` — deliberately.** The
    /// clause used to be worded as "the consistency check between `bases` and `read_coverage`",
    /// and that check cannot be written: `positions_covered` is *derived from* the read's
    /// events, so checking it against them is tautological, and §8's trap is that no inequality
    /// relates `bases.len()` to the footprint in general — an insertion adds bases without
    /// positions, a deletion positions without bases, so a row's byte count may exceed or fall
    /// short of its position count either way. Spec §13 now asks for what is checkable
    /// (corrected 2026-07-30); the code here did not change.
    ///
    /// What ties the bases to the coverage is construction, not an assertion: both come out of
    /// one `apply_events_into` call, which returns the bases and the witnessed extent together
    /// and returns `None` rather than a run it cannot describe honestly (spec §6).
    fn push_locus(&mut self, locus: &SampleLocusObservations, contig: &str) {
        self.generic_loci += 1;
        let footprint = locus.region.len();
        self.reads_without_observation += u64::from(locus.reads_without_observation);
        self.reads_discarded_by_cap += u64::from(locus.reads_discarded_by_cap);
        let depth: u32 = locus.complete_observations().map(|obs| obs.num_obs).sum();
        assert_eq!(
            locus.kind,
            LocusKind::Generic,
            "the generic generator emitted a locus of another kind at {:?}",
            locus.region,
        );
        for obs in &locus.observed_sequences {
            match obs.read_coverage {
                ReadCoverage::Complete => {
                    self.rows_complete += 1;
                    self.reads_complete += u64::from(obs.num_obs);
                }
                ReadCoverage::Observed {
                    offset_in_locus,
                    positions_covered,
                } => {
                    self.rows_observed += 1;
                    self.reads_observed += u64::from(obs.num_obs);
                    let reach = u64::from(offset_in_locus) + u64::from(positions_covered);
                    assert!(
                        reach <= footprint,
                        "a row at {:?} claims positions {}..{} of a {footprint}-position \
                         locus — a witness cannot reach past the footprint it is measured \
                         against",
                        locus.region,
                        offset_in_locus,
                        reach,
                    );
                    assert!(
                        positions_covered > 0,
                        "a row at {:?} witnessed zero positions, which is not an \
                         observation at all — it belongs in reads_without_observation",
                        locus.region,
                    );
                }
            }
            self.rows.push(ObservationRow {
                contig: contig.to_string(),
                start: locus.region.start.get(),
                end: locus.region.end.get(),
                ref_bases: locus.reference_bases.to_vec(),
                depth,
                read_coverage: coverage_label(obs.read_coverage),
                read_group: obs.read_group.0,
                observed: obs.bases.to_vec(),
                reads: obs.num_obs,
                chain_ids: obs.chain_ids.clone(),
            });
        }
    }

    /// The dump text: the `#` header lines, the TSV column line, then the rows.
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        // Writing into a `String` via `fmt::Write` never fails, so results are ignored.
        let _ = writeln!(
            out,
            "# generic_loci={} generic_regions={} generic_region_bp={} \
             records_outside_region={}",
            self.generic_loci,
            self.generic_regions,
            self.generic_region_bp,
            self.records_outside_region,
        );
        // **Two lines, because these counters are in two units that do not add up.** They used
        // to print on one, under a module doc telling the reader every admitted read is
        // somewhere — and no arithmetic relates them: the three below are once per read for the
        // run, while `reads_without_observation` and `reads_discarded_by_cap` are *per-locus*
        // quantities summed over loci, so one read counts once per locus it was affected at.
        // Measured on the accounting fixture with both depth caps forced to 1:
        // `reads_admitted=4` beside `reads_discarded_by_cap=55`. The `locus_sum_` prefix and
        // the split are what stop a reader trying to balance them.
        let _ = writeln!(
            out,
            "# reads_admitted={} reads_declined_by_preparer={} \
             reads_silent_over_footprint={}",
            self.reads_admitted, self.reads_declined_by_preparer, self.reads_silent_over_footprint,
        );
        let _ = writeln!(
            out,
            "# locus_sum_reads_without_observation={} locus_sum_reads_discarded_by_cap={}",
            self.reads_without_observation, self.reads_discarded_by_cap,
        );
        let _ = writeln!(
            out,
            "# rows_complete={} rows_observed={} reads_complete={} reads_observed={}",
            self.rows_complete, self.rows_observed, self.reads_complete, self.reads_observed,
        );
        let _ = writeln!(
            out,
            "# record_widen_events={} column_depth_truncations={} regions_in={} \
             regions_handled={} loci_emitted={}",
            self.record_widen_events,
            self.column_depth_truncations,
            self.regions_in,
            self.regions_handled,
            self.loci_emitted,
        );
        out.push_str(
            "contig\tstart\tend\tref_bases\tdepth\tread_coverage\tread_group\tobserved\treads\tchain_ids\n",
        );
        for row in &self.rows {
            let ids: Vec<String> = row.chain_ids.iter().map(|id| id.to_string()).collect();
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.contig,
                row.start,
                row.end,
                String::from_utf8_lossy(&row.ref_bases),
                row.depth,
                row.read_coverage,
                row.read_group,
                String::from_utf8_lossy(&row.observed),
                row.reads,
                ids.join(","),
            );
        }
        out
    }
}

/// The tag a coverage carries in the `read_coverage` column.
///
/// Destructured rather than matched with `_`, so a future `ReadCoverage` variant is a compile
/// error here instead of a silently mislabelled row. The **run** is printed rather than a
/// side label: on the generic path a read can be blind in the *middle* of a footprint, where
/// "partial:left" would be a lie, and the offset is what a consumer needs to place the
/// evidence anyway.
fn coverage_label(coverage: ReadCoverage) -> String {
    match coverage {
        ReadCoverage::Complete => "complete".to_string(),
        ReadCoverage::Observed {
            offset_in_locus,
            positions_covered,
        } => format!("observed:{offset_in_locus}+{positions_covered}"),
    }
}

/// Split a typed region into `chunk_bp`-sized pieces when it is `Generic`, leaving every other
/// kind alone.
///
/// Lazy by construction (it returns a small `Vec` per input item, consumed through
/// `flat_map`), because a chromosome's typed regions must not be collected — the whole point
/// of the region stream is that nothing region-shaped is resident (spec §7).
///
/// Adjacent pieces tile the original exactly: the last piece ends where the region does, so a
/// record whose footprint crosses a piece boundary is the halo case, which is what the
/// boundary check exercises (spec §2).
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
        // `saturating_add`, because `chunk` comes from the environment: a plain `+` panics on
        // `PVC_GENERIC_REGION_CHUNK_BP=18446744073709551615` in debug and wraps in release,
        // where wrapping would make `piece_end` land *below* `at` and the loop never terminate.
        // Saturating gives the honest answer for a chunk wider than the contig — one piece.
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

/// Run the whole pipeline over `fasta` + `bams`, optionally restricted to `contig_filter` (a
/// contig name), preparing reads with `preparer`, and build the [`DumpReport`].
///
/// # One `.fai`, parsed once
///
/// The generator needs **two** reference accessors — one for the walk's REF fetches and one
/// per file for the read query's mismatch-fraction filter — and the factory is called at every
/// region. A fresh `WindowedRefSeq::new` per region parses the whole `.fai` on its first fetch:
/// measured at **189 µs on a GRCh38-shaped reference against a ~120 µs per-region walk
/// constant**, so the factory would more than double the cost of a region. `read_index` plus
/// [`WindowedRefSeq::with_shared_index`] takes that to **19 µs**, which is the `open(2)` each
/// file's own cursor costs and the floor unless handles are pooled. This tool is that fix's
/// first caller.
fn run_dump<P: ReadPreparer + 'static>(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: Option<&str>,
    preparer: P,
    config: PileupGeneratorConfig,
    chunk_bp: Option<u64>,
    // **The caller's cache, not one of our own.** `main` has already read this FASTA to build
    // the preparer's accessor; taking its cache makes this call a hit instead of a second read
    // and a second verification of the same file.
    cache: &Arc<ReferenceInfoCache>,
) -> Result<DumpReport, Box<dyn std::error::Error>> {
    let (info, verify) = read_reference_verifying_or_creating_fai(cache, fasta.to_path_buf())?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;

    let sample = SampleReads::open_only_sample(bams, &info, ReadFilterConfig::default(), true)?;

    // **An `Arc` around a type that is neither `Send` nor `Sync`, and clippy is right to
    // notice.** `WindowedRefSeq` holds an open reader behind a `RefCell`, so it is `!Sync`,
    // and `Arc<T>` is only `Send` when `T: Send + Sync`.
    //
    // **The `Arc` stays, and an `Rc` would not say anything truer.** That was the carried
    // finding and it was wrong on the decisive point: `Rc` and `Arc` are *equally* unable to
    // cross a thread here, so swapping them buys no honesty — it only moves which pointer
    // names the same single-threaded fact. What makes the signature right is that
    // `PileupGenerator` is generic over the accessor: for an `InMemoryRefSeq`, which *is*
    // `Send + Sync`, the `Arc` is meaningful, and it is `ref_seq.rs`'s recorded choice of
    // `Arc` over `Rc` for the per-sample fan-out to come. This call site is the one where the
    // accessor happens to be file-backed, so this is where the lint is allowed rather than
    // where the signature changes. (Owner, 2026-07-30.)
    //
    // **What the fan-out will actually change is not this pointer, and the spec already says
    // so.** `locus_generation.md` §9: *"a shared `WindowedRefSeq` is `Send` but not `Sync`, so
    // each worker needs its own accessor"* — so the planned shape is one accessor per worker,
    // never one shared handle, and no pointer type makes a shared one work. That is also what
    // `with_shared_index` is for: `WindowedRefSeq` already holds `Arc<ContigList>` and
    // `Arc<fai::Index>`, sharing the immutable index while keeping the reader per accessor,
    // because k accessors over one index must still hold k cursors. Separately,
    // `PileupGenerator` is already `!Send` through its `Rc<RefCell<ReadPreparation>>`, so the
    // handle on the reference is not what stands between this and a worker.
    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new is generic over the accessor and takes Arc, which is meaningful for the Send+Sync in-memory ones; this accessor is file-backed and single-threaded — see the comment above"
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

    // **One region stream for the whole run, not one per contig.** The locus iterator owns the
    // stream, the reads and the generator set; handing it a fresh one per contig would mean
    // taking those back out between contigs, and the generator in particular must not be
    // rebuilt — its chain-id allocator is run-lifetime, so restarting it would give two
    // fragments on two contigs the same identity (spec §8).
    //
    // The per-contig typed-region walks are built up front and chained. Each is a path, a
    // contig table and a shared index — no I/O until its first fetch — so the ones a contig
    // filter excludes are never built at all, and the ones it keeps cost nothing until walked.
    let walk_config = TypedRegionConfig::default();
    let mut walks = Vec::new();
    for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
        if contig_filter.is_some_and(|name| entry.name != name) {
            continue;
        }
        walks.push(TypedRegionIterator::over_contig(
            WindowedRefSeq::with_shared_index(fasta.to_path_buf(), contigs.clone(), index.clone()),
            ContigId(index_of_contig as u32),
            walk_config.clone(),
        )?);
    }

    // Counted inside the stream, because the stream is moved into the locus iterator and
    // nothing downstream sees the regions: the dump wants to say how many *generic* base pairs
    // it walked, which is the denominator for "uncovered positions yield no locus".
    let generic_bp = Rc::new(Cell::new(0u64));
    let generic_regions = Rc::new(Cell::new(0u64));
    let regions = {
        let generic_bp = generic_bp.clone();
        let generic_regions = generic_regions.clone();
        walks
            .into_iter()
            .flatten()
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

    let mut report = DumpReport::default();
    let mut stream = SampleLocusObservationsIterator::new(regions, sample, generators);
    for locus in &mut stream {
        let locus = locus?;
        let contig = contigs
            .entries
            .get(locus.region.contig.0 as usize)
            .map(|entry| entry.name.as_str())
            .unwrap_or("?");
        report.push_locus(&locus, contig);
    }
    let counts = stream.counts();
    report.regions_in = counts.regions_in;
    report.regions_handled = counts.regions_handled;
    report.loci_emitted = counts.loci_emitted;
    // **Read through the boxed generator**, which is the surface `GeneratorSet::generic_counts`
    // exists for: a `GeneratorSlot` erases the type, so before it these nine counters had no
    // reader outside the generator's own module — and a walk that emitted nothing for a covered
    // region was therefore totally silent, with the truncations that explained it counted and
    // unreachable (Milestone C review). This tool is that accessor's first non-test caller.
    let Some(GeneratorCounts::Pileup(generic)) = stream.generators().generic_counts() else {
        return Err("the generic slot reported no counts, so it is not the filled one".into());
    };
    report.reads_admitted = generic.reads_admitted;
    report.reads_declined_by_preparer = generic.reads_declined_by_preparer;
    report.reads_silent_over_footprint = generic.reads_silent_over_footprint;
    report.record_widen_events = generic.record_widen_events;
    report.records_outside_region = generic.records_outside_region;
    report.column_depth_truncations = generic.column_depth_truncations;
    report.generic_region_bp = generic_bp.get();
    report.generic_regions = generic_regions.get();

    // The background `.fai` verification (only when a pre-existing `.fai` was used) is joined
    // at the end, so a stale index is a failure rather than a silently-wrong walk.
    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(report)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_generic_loci_dump <reference.fa> <sample.bam|cram> [contig]\n\
             dumps, per covered reference position, the sequences one sample's reads showed \
             and how much of each locus each read witnessed.\n\
             PVC_GENERIC_REGION_CHUNK_BP=n splits each Generic region into n-bp pieces."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter = args.get(3).map(String::as_str);
    // **Parsed strictly, because the old `.ok().and_then(|v| v.parse().ok())` turned off the
    // check this knob exists for.** `0`, `-1`, `abc` and an out-of-range value all became
    // `None`, which is "do not chunk" — so a typo silently produced a single-region run that
    // still exits 0 and still prints a plausible report, and the region-boundary comparison the
    // knob is for was simply not made. A knob whose misuse is indistinguishable from its absence
    // is worse than no knob.
    let chunk_bp = match std::env::var("PVC_GENERIC_REGION_CHUNK_BP") {
        Err(_) => None,
        Ok(value) => match value.parse::<u64>() {
            Ok(0) => {
                eprintln!(
                    "error: PVC_GENERIC_REGION_CHUNK_BP=0 — a zero-width piece would not \
                     advance; omit the variable to walk each Generic region whole"
                );
                return ExitCode::from(2);
            }
            Ok(bp) => Some(bp),
            Err(error) => {
                eprintln!(
                    "error: PVC_GENERIC_REGION_CHUNK_BP={value:?} is not a base count: {error}"
                );
                return ExitCode::from(2);
            }
        },
    };

    // ng's real read preparation: left-alignment over the canonical (uppercased) reference.
    // Its own reference accessor shares the parsed index with the generator's, for the reason
    // `run_dump` gives.
    // **One cache, shared with `run_dump`, and the handle is joined rather than dropped.**
    // Two defects lived here. `run_dump` builds its own `ReferenceInfoCache`, so the FASTA was
    // read and verified **twice** from two caches — passing this one in makes the second call a
    // cache hit, which is what a caller-held single-flight cache is for. And the
    // `#[must_use] VerificationHandle` was bound to `_verify` and dropped, which printed
    // *"verification handle dropped without join(); a possibly-pending reference-verification
    // error went unobserved"* on every successful run — visibly, at the top of D3's own
    // measurement output — and meant a stale `.fai` produced a silently wrong walk.
    //
    // Joined **after** the run, not before: verification hashes the whole FASTA on a background
    // thread, and overlapping the walk with it is the entire point of the handle. Blocking here
    // would trade a real defect for a slow tool.
    let cache = Arc::new(ReferenceInfoCache::new());
    let report = match read_reference_verifying_or_creating_fai(&cache, fasta.clone()) {
        Ok((info, verify)) => {
            let contigs = Arc::new(info.contig_list());
            match WindowedRefSeq::read_index(&fasta) {
                Ok(index) => {
                    let preparer = LeftAlignPreparer::with_default_normalizer(
                        WindowedRefSeq::with_shared_index(fasta.clone(), contigs, index),
                    );
                    let walked = run_dump(
                        &fasta,
                        &[bam],
                        contig_filter,
                        preparer,
                        PileupGeneratorConfig::default(),
                        chunk_bp,
                        &cache,
                    );
                    match (walked, verify) {
                        (Ok(report), Some(handle)) => handle
                            .join()
                            .map(|_verified| report)
                            .map_err(|error| Box::new(error) as Box<dyn std::error::Error>),
                        (Ok(report), None) => Ok(report),
                        // The walk's own error wins: it happened first and describes the run.
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
            print!("{}", report.render());
            ExitCode::SUCCESS
        }
        Err(error) => {
            // **The whole `#[source]` chain, because the top of it names nothing actionable.**
            // A missing BAM printed only `error: reading the input files' read groups failed` —
            // the path and the `No such file or directory` were both in the chain and both
            // dropped, so the one message a user of this tool is most likely to see told them
            // neither what was missing nor where it was looked for.
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
    use super::*;
    use std::collections::BTreeMap;
    use std::fs::File;

    use pop_var_caller::ng::read::{AlignedRead, PreparedRead, ReadPrepError};
    use tempfile::TempDir;

    // ---------------------------------------------------------------------
    // The reference, and the preparer the fixtures need
    // ---------------------------------------------------------------------

    /// A 120 bp aperiodic contig — three 40 bp blocks with no tandem repeat long enough for
    /// region typing to carve an `SsrSegment` out of, so the whole contig types as `Generic`
    /// and the generic generator is what the walk reaches. (They are the *flank* constants of
    /// `ng_ssr_loci_dump`'s fixture, chosen there for the same aperiodicity.)
    const BLOCK_A: &[u8] = b"GATCTTGCAAGCTGGAATCCGTTACGATCGGATCAAGCTT";
    const BLOCK_B: &[u8] = b"CAGTTGCACGATCCTAAGGCTTGACCATGGATCCAAGTTG";
    const BLOCK_C: &[u8] = b"GGTTCAAGATCCGGATCTTGCAATCGGATCAAGCTTGACT";

    fn contig() -> Vec<u8> {
        [BLOCK_A, BLOCK_B, BLOCK_C].concat()
    }

    /// The reference base at 1-based `pos`.
    fn ref_base(pos: u64) -> u8 {
        contig()[pos as usize - 1]
    }

    /// The reference bases over 1-based `start..=end`.
    fn ref_span(start: u64, end: u64) -> Vec<u8> {
        contig()[start as usize - 1..end as usize].to_vec()
    }

    /// Write `chr1` as a one-line FASTA (its `.fai` is created by `run_dump`).
    fn write_fasta(path: &Path, bases: &[u8]) {
        use std::io::Write as _;
        let mut file = File::create(path).unwrap();
        writeln!(file, ">chr1").unwrap();
        file.write_all(bases).unwrap();
        writeln!(file).unwrap();
    }

    /// **The preparer the fixtures drive**: ng's real passthrough, plus a per-read adaptor
    /// boundary applied by qname.
    ///
    /// # Why the boundary is placed here and not by the real preparer
    ///
    /// Two of spec §12's six fixtures need a read whose bases are **silenced over part of a
    /// footprint**, and the G1 adaptor filter is the mechanism the walk has for that
    /// (`cigar_cursor::base_in_adaptor`). In a real run the boundary is derived from mate
    /// geometry — a fragment shorter than the read, so the 3′ end reads into the adaptor — and
    /// deriving it *is* production's code, tested there. Building a BAM whose mate geometry
    /// happens to land the boundary at the position a fixture needs would make the fixture a
    /// test of that derivation instead of a test of the walk, and would move the boundary
    /// whenever the derivation changed.
    ///
    /// So the boundary is stated — **on the way in**, on the `AlignedRead`, which is exactly
    /// where mate geometry puts it in a real run. Everything after that is ng's real
    /// preparation: this wraps [`LeftAlignPreparer`] and delegates. The property under test is
    /// the one spec §4 changed: **what the fold does with a position the read did not
    /// witness.**
    struct MaskingPreparer {
        inner: LeftAlignPreparer<WindowedRefSeq>,
        /// qname prefix → the boundary to set. Forward strand, so every base at
        /// `ref_pos >= boundary` is masked (`cigar_cursor::base_in_adaptor`).
        boundaries: Vec<(&'static str, u32)>,
    }

    impl ReadPreparer for MaskingPreparer {
        type Scratch = <LeftAlignPreparer<WindowedRefSeq> as ReadPreparer>::Scratch;

        fn prepare_read(
            &self,
            mut read: AlignedRead,
            scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            const REVERSE: u16 = 0x10;
            for (prefix, boundary) in &self.boundaries {
                if read.qname.starts_with(prefix.as_bytes()) {
                    assert_eq!(
                        read.flag & REVERSE,
                        0,
                        "the fixture boundary is a forward-strand one; on a reverse read \
                         `ref_pos <= boundary` would silence the other half of the read and \
                         the fixture would be testing the opposite of what it says"
                    );
                    read.adaptor_boundary = Some(*boundary);
                }
            }
            self.inner.prepare_read(read, scratch)
        }
    }

    /// A preparer over the fixture's reference, masking the named reads from the given
    /// positions.
    fn masking(fasta: &Path, boundaries: &[(&'static str, u32)]) -> MaskingPreparer {
        let cache = Arc::new(ReferenceInfoCache::new());
        let (info, _verify) =
            read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf()).unwrap();
        let contigs = Arc::new(info.contig_list());
        let index = WindowedRefSeq::read_index(fasta).unwrap();
        MaskingPreparer {
            inner: LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
                fasta.to_path_buf(),
                contigs,
                index,
            )),
            boundaries: boundaries.to_vec(),
        }
    }

    // ---------------------------------------------------------------------
    // Reads and BAMs
    // ---------------------------------------------------------------------

    /// A read at 1-based `start` with `cigar`, whose sequence is `seq`.
    ///
    /// `seq` is given explicitly rather than sliced from the contig, because three of the
    /// fixtures need a base that is *not* the reference one (an `N`) or a length the CIGAR
    /// dictates (a deletion consumes reference without consuming read).
    fn read(
        qname: &str,
        start: usize,
        cigar: &[(noodles_sam::alignment::record::cigar::op::Kind, usize)],
        seq: Vec<u8>,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_core::Position;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::{Flags, MappingQuality};
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
        let quals = vec![40u8; seq.len()];
        noodles_sam::alignment::RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(
                cigar
                    .iter()
                    .map(|(kind, len)| Op::new(*kind, *len))
                    .collect(),
            )
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(quals))
            .build()
    }

    /// A reference-matching read spanning `start..start + len - 1`.
    fn matching_read(qname: &str, start: usize, len: usize) -> noodles_sam::alignment::RecordBuf {
        use noodles_sam::alignment::record::cigar::op::Kind;
        let bases = contig()[start - 1..start - 1 + len].to_vec();
        read(qname, start, &[(Kind::Match, len)], bases)
    }

    /// Stamp a record with an `RG` tag, so a multi-read-group file can say which group each
    /// record belongs to.
    fn in_read_group(
        mut record: noodles_sam::alignment::RecordBuf,
        read_group: &str,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;
        record.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(read_group.as_bytes().to_vec().into()),
        );
        record
    }

    /// Write a coordinate-sorted single-contig BAM declaring `read_groups`, all naming one
    /// sample (a `SampleReads` open serves exactly one).
    fn write_bam(
        path: &Path,
        contig_len: usize,
        reads: &[noodles_sam::alignment::RecordBuf],
        read_groups: &[&str],
    ) {
        use bstr::BString;
        use noodles_bam as bam;
        use noodles_sam as sam;
        use sam::alignment::io::Write as _;
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::Version;
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
        use std::num::NonZero;

        let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
        let sq = Map::<ReferenceSequence>::new(NonZero::new(contig_len).unwrap());
        let mut builder = sam::Header::builder()
            .set_header(hd)
            .add_reference_sequence(b"chr1".to_vec(), sq);
        for name in read_groups {
            let mut rg = Map::<ReadGroup>::default();
            rg.other_fields_mut()
                .insert(SAMPLE, BString::from("sample0"));
            builder = builder.add_read_group(name.as_bytes().to_vec(), rg);
        }
        let header = builder.build();

        let mut writer = bam::io::Writer::new(File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        for record in reads {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    /// Write a fixture (FASTA + BAM) into a fresh temporary directory.
    fn fixture(
        reads: &[noodles_sam::alignment::RecordBuf],
        read_groups: &[&str],
    ) -> (TempDir, PathBuf, PathBuf) {
        let bases = contig();
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("ref.fa");
        let bam = dir.path().join("sample.bam");
        write_fasta(&fasta, &bases);
        write_bam(&bam, bases.len(), reads, read_groups);
        (dir, fasta, bam)
    }

    /// Dump a fixture with the given preparer.
    fn dump(
        fasta: &Path,
        bam: &Path,
        boundaries: &[(&'static str, u32)],
        chunk_bp: Option<u64>,
    ) -> DumpReport {
        let preparer = masking(fasta, boundaries);
        // A cache of its own here, unlike `main`'s: these fixtures have no earlier read to
        // share with, so a fresh one is the honest thing rather than a hit.
        let cache = Arc::new(ReferenceInfoCache::new());
        run_dump(
            fasta,
            std::slice::from_ref(&bam.to_path_buf()),
            None,
            preparer,
            PileupGeneratorConfig::default(),
            chunk_bp,
            &cache,
        )
        .expect("the dump succeeds")
    }

    /// The rows at 1-based `start`, in emission order.
    fn rows_at(report: &DumpReport, start: u64) -> Vec<&ObservationRow> {
        report
            .rows
            .iter()
            .filter(|row| row.start == start)
            .collect()
    }

    /// The row at `start` whose observed bases are `observed` — panicking if there is not
    /// exactly one, because "the row I meant" is the claim every assertion below rests on.
    fn row_at<'a>(report: &'a DumpReport, start: u64, observed: &[u8]) -> &'a ObservationRow {
        let matching: Vec<&ObservationRow> = rows_at(report, start)
            .into_iter()
            .filter(|row| row.observed == observed)
            .collect();
        assert_eq!(
            matching.len(),
            1,
            "expected exactly one row at {start} observing {:?}, found {:?}",
            String::from_utf8_lossy(observed),
            rows_at(report, start)
                .iter()
                .map(|row| (
                    String::from_utf8_lossy(&row.observed).to_string(),
                    row.read_coverage.clone(),
                    row.read_group,
                ))
                .collect::<Vec<_>>(),
        );
        matching[0]
    }

    // ---------------------------------------------------------------------
    // The whole-fixture properties (spec §13)
    // ---------------------------------------------------------------------

    /// **A covered position yields exactly one locus; an uncovered one yields none.**
    ///
    /// The fixture leaves the contig's first ten and last ten positions uncovered, so
    /// `generic_loci` must be strictly less than `generic_region_bp` — a walk that emitted a
    /// locus per *region* position rather than per covered position would fail here, and so
    /// would one that emitted two loci for a position two reads cover.
    #[test]
    fn every_covered_position_yields_one_locus_and_no_uncovered_one_does() {
        let (_dir, fasta, bam) = fixture(
            &[matching_read("a", 11, 40), matching_read("b", 21, 40)],
            &["rg0"],
        );
        let report = dump(&fasta, &bam, &[], None);

        let anchors: Vec<u64> = {
            let mut anchors: Vec<u64> = report.rows.iter().map(|row| row.start).collect();
            anchors.dedup();
            anchors
        };
        assert_eq!(
            anchors,
            (11..=60).collect::<Vec<u64>>(),
            "the union of the two reads' spans, every position once"
        );
        assert_eq!(report.generic_loci, 50, "one locus per covered position");
        assert_eq!(
            report.generic_loci, report.loci_emitted,
            "the dispatcher's tally and the rows must describe the same loci"
        );
        assert!(
            report.generic_region_bp > report.generic_loci,
            "the fixture leaves positions uncovered ({} bp of generic region against {} \
             loci), and they must yield nothing",
            report.generic_region_bp,
            report.generic_loci,
        );
    }

    /// **Every read the walk saw is accounted for** (spec §13), on a fixture where each class
    /// is populated by hand — which is the only honest oracle for counters production does not
    /// keep (spec §3, class 3).
    ///
    /// Four reads, one per class: one that contributes normally, one silenced over every
    /// position it touches, one blind in the middle of a footprint (so it folds and then loses
    /// its row), and — through the depth cap — one discarded.
    #[test]
    fn every_read_the_walk_saw_is_accounted_for() {
        use noodles_sam::alignment::record::cigar::op::Kind;
        // `blind` matches 11..40 except for an `N` at position 22, which the cursor emits no
        // event for — and 22 is *inside* the record `deleter` opens at 20, so the witness
        // there is non-contiguous and the read loses its row.
        // The `N` goes at position **22**, inside the record `deleter` opens at 20 (20..=24) —
        // a hole outside every footprint costs the read nothing and the counter would stay
        // zero, which is how the first draft of this fixture passed for the wrong reason.
        let mut blind_seq = contig()[10..40].to_vec();
        blind_seq[11] = b'N';
        let reads = vec![
            matching_read("speaks", 11, 30),
            // Silenced from its own first position: admitted, never a contributor.
            matching_read("silent", 11, 30),
            read("blind", 11, &[(Kind::Match, 30)], blind_seq),
            // A deletion anchored at 20 covering 21..24, so the record at 20 spans 20..=24 and
            // `blind`'s hole at 25 is *outside* it — the hole matters at the loci 25 anchors.
            read(
                "deleter",
                11,
                &[(Kind::Match, 10), (Kind::Deletion, 4), (Kind::Match, 20)],
                [contig()[10..20].to_vec(), contig()[24..44].to_vec()].concat(),
            ),
        ];
        let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
        let report = dump(&fasta, &bam, &[("silent", 11)], None);

        assert_eq!(
            report.reads_admitted, 4,
            "all four reads reached the walk — silencing happens inside it"
        );
        assert_eq!(
            report.reads_declined_by_preparer, 0,
            "this preparer declines nothing; the counter exists so that fact is stated"
        );
        assert_eq!(
            report.reads_silent_over_footprint, 1,
            "`silent` contributed at no position, and neither per-locus counter can see it"
        );
        assert!(
            report.reads_without_observation > 0,
            "`blind`'s hole must cost it its row at the loci whose footprint spans it"
        );
        // Every row's reads, plus the reads counted out, must be the observations the walk
        // made. `reads_admitted` is not that number — a read contributes at many positions —
        // so the identity is per class, which is why all five counters are printed.
        assert_eq!(
            report.reads_complete + report.reads_observed,
            report
                .rows
                .iter()
                .map(|row| u64::from(row.reads))
                .sum::<u64>(),
            "the per-class read totals are the rows' reads, split by coverage"
        );
    }

    /// **No observation carries a chain id for a read that agreed with the reference across
    /// everything it witnessed** — spec §6's invariant, stated per read because row splitting
    /// means "the REF row" is no longer a unique row (§13).
    ///
    /// # It runs on the masking fixture, because that is where the two rules differ
    ///
    /// Production's rule is positional — drop the ids of `alleles[0]` — and ng's is per read:
    /// did *this read* agree with the reference across everything it witnessed. They coincide
    /// exactly when every row is a whole-footprint witness, so a fixture of complete reads
    /// cannot tell them apart, and the first draft of this test (a matching read and a deletion
    /// read) **passed with the positional rule substituted**.
    ///
    /// The case that separates them is a **partial** row whose bases are a prefix of the
    /// footprint's reference bytes: `masked` witnessed positions 20 and 21 and matched the
    /// reference at both, so it agreed with the reference across everything it saw and carries
    /// no id — while its row is not `alleles[0]` (which holds all five reference bases), so the
    /// positional rule would tag it. All three cases are asserted here: that row, a whole-
    /// footprint reference match, and a read that genuinely departed.
    #[test]
    fn only_the_rows_that_departed_from_the_reference_carry_a_chain_id() {
        let (_dir, fasta, bam) = fixture(&masking_fixture_reads(), &["rg0"]);
        let report = dump(&fasta, &bam, &[("masked", 22)], None);

        // 1. A plain matched position: one row, the reference base, no id.
        let plain = rows_at(&report, 11);
        assert_eq!(plain.len(), 1, "position 11 is a clean pileup column");
        assert_eq!(plain[0].observed, vec![ref_base(11)]);
        assert!(
            plain[0].chain_ids.is_empty(),
            "a whole-footprint reference match needs no chain id — the reference is the \
             default haplotype and a default needs no tag"
        );

        // 2. The partial row that matched the reference over everything it witnessed. **This is
        //    the one production's positional rule gets wrong.**
        let partial = rows_at(&report, 20)
            .into_iter()
            .find(|row| row.read_coverage == "observed:0+2")
            .expect("the masked read's partial row");
        assert_eq!(
            partial.observed,
            ref_span(20, 21),
            "the masked read matched the reference at both positions it witnessed"
        );
        assert!(
            partial.chain_ids.is_empty(),
            "a read that agreed with the reference across everything it witnessed carries no \
             chain id, even though its row is not the locus's reference row — carrying ids {:?}",
            partial.chain_ids,
        );

        // 3. And a row that genuinely departed keeps its id, or nothing could be phased.
        let departed = rows_at(&report, 20)
            .into_iter()
            .find(|row| row.observed == vec![ref_base(20)])
            .expect("the opener's row: one base, having deleted the other four positions");
        assert_eq!(
            departed.read_coverage, "complete",
            "the opener witnessed all five positions — it deleted four of them"
        );
        assert!(
            !departed.chain_ids.is_empty(),
            "a read carrying a deletion departed from the reference and must keep its id"
        );
    }

    /// **Output is byte-identical across repeated runs** (spec §13), which for a per-read
    /// re-derivation off an `AHashMap` is not free — the row order and every `q_sum` depend on
    /// the order the reads are taken in.
    ///
    /// Within one process `ahash` seeds once, so this cannot see a *seed* dependency; that is
    /// `parity::ng_emits_the_same_bytes_in_a_second_process`'s job, in two processes. What this
    /// catches is the cheaper failure: any dependence on state left over from the previous run.
    #[test]
    fn the_dump_is_byte_identical_across_runs() {
        let (_dir, fasta, bam) = fixture(
            &[matching_read("a", 11, 30), matching_read("b", 21, 30)],
            &["rg0"],
        );
        let first = dump(&fasta, &bam, &[], None).render();
        let second = dump(&fasta, &bam, &[], None).render();
        assert_eq!(first, second);
        assert!(
            first.starts_with("# generic_loci="),
            "the counts header comes first, `#`-prefixed (spec §13):\n{first}"
        );
    }

    // ---------------------------------------------------------------------
    // The six fixtures spec §12 owes, because no inherited test exercises the defect
    // ---------------------------------------------------------------------

    /// The reads shared by the masking fixtures: a deletion that opens a **five-position**
    /// record at 20, a read whose own deletion covers the tail of that footprint, and a read
    /// masked from 22.
    ///
    /// `deleter2` and `masked` are the pair that matters: both observe the same two bases —
    /// positions 20 and 21 — and they observe them for **different reasons**. `deleter2`
    /// witnessed all five positions and deleted three of them; `masked` witnessed two and is
    /// silent about the rest. Production cannot tell them apart: it folds `masked` with the
    /// reference bases at 22, 23 and 24 appended, which lands it in `deleter2`'s bucket or the
    /// REF one depending on the bases. ng emits **two rows with the same bases** and different
    /// `read_coverage`, which is what proves the field is computed rather than defaulted
    /// (spec §13).
    /// **Every read here carries at least 30 sequenced bases, and that is a constraint of the
    /// real filter, not of the walk.** `ReadFilterConfig`'s `DEFAULT_MIN_READ_LENGTH` is 30, so
    /// a shorter read never reaches the generator at all — it is dropped at ingestion, which
    /// looks from the dump exactly like a read the walk chose to ignore. The first draft of
    /// these fixtures had 26-base indel reads and reported `reads_admitted` 3 of 4.
    fn masking_fixture_reads() -> Vec<noodles_sam::alignment::RecordBuf> {
        use noodles_sam::alignment::record::cigar::op::Kind;
        vec![
            // Opens the record at 20: matched 11..20, deletes 21..24, matched 25..44.
            read(
                "opener",
                11,
                &[(Kind::Match, 10), (Kind::Deletion, 4), (Kind::Match, 20)],
                [contig()[10..20].to_vec(), contig()[24..44].to_vec()].concat(),
            ),
            // Witnesses all five positions, emitting bases only at 20 and 21: matched 11..21,
            // deletes 22..24, matched 25..43.
            read(
                "deleter2",
                11,
                &[(Kind::Match, 11), (Kind::Deletion, 3), (Kind::Match, 19)],
                [contig()[10..21].to_vec(), contig()[24..43].to_vec()].concat(),
            ),
            // Reference-matching across the whole span; the preparer masks it from 22.
            matching_read("masked", 11, 30),
        ]
    }

    /// **Fixture 1 — a read adaptor-masked over part of a multi-base record is `Observed`, not
    /// `Complete`** (spec §12.1), *and* **fixture 1b — an `Observed` row is separate from a
    /// `Complete` one with the same bases** (spec §13).
    ///
    /// A span-derived implementation passes every other test in the suite and fails this one:
    /// `masked`'s alignment spans the whole footprint, so a coverage read off the span says
    /// `Complete`. Only the *events* know that positions 22, 23 and 24 were silenced.
    #[test]
    fn a_read_masked_over_part_of_a_footprint_is_observed_not_complete() {
        let (_dir, fasta, bam) = fixture(&masking_fixture_reads(), &["rg0"]);
        let report = dump(&fasta, &bam, &[("masked", 22)], None);

        let locus_ref = ref_span(20, 24);
        assert_eq!(locus_ref.len(), 5, "the record at 20 spans 20..=24");
        let two_bases = ref_span(20, 21);

        // The two rows that carry the *same* bases for different reasons.
        let rows: Vec<&ObservationRow> = rows_at(&report, 20)
            .into_iter()
            .filter(|row| row.observed == two_bases)
            .collect();
        let coverages: Vec<&str> = rows.iter().map(|row| row.read_coverage.as_str()).collect();
        assert_eq!(
            coverages,
            vec!["complete", "observed:0+2"],
            "the same two bases must appear twice at locus 20 — once from a read that \
             witnessed all five positions and deleted three, once from a read that witnessed \
             two and was silenced over the rest. Rows: {:?}",
            rows_at(&report, 20)
                .iter()
                .map(|row| (
                    String::from_utf8_lossy(&row.observed).to_string(),
                    row.read_coverage.clone()
                ))
                .collect::<Vec<_>>(),
        );
        // And the masked read's bases are what it saw, not the footprint's. Selected by its
        // *coverage*, because the whole point of the pair above is that the bases do not
        // identify it.
        let masked_row = rows
            .iter()
            .find(|row| row.read_coverage == "observed:0+2")
            .expect("the masked read's row");
        assert_ne!(
            masked_row.observed, locus_ref,
            "the masked read must not carry the five reference bases of the footprint — three \
             of them it was silenced over"
        );
        assert_eq!(
            report.end_of(20),
            24,
            "the locus's end is the footprint's last position, which only ng's type can say"
        );
    }

    /// **Fixture 2 — an interior `N`, and a ref-skip: no observation, and counted**
    /// (spec §12.2).
    ///
    /// The witness is *non-contiguous*, which is the one case that yields no row at all: a
    /// sequence with a hole in it is not an observation of anything, and A5's answer is to drop
    /// the row and record the read. Both spellings are checked because they reach it by
    /// different routes — an `N` base emits no `Match` event, a ref-skip consumes reference
    /// without emitting one — and a fill-preserving implementation reports both as reference
    /// matches with no hole at all.
    #[test]
    fn a_read_blind_inside_a_footprint_yields_no_observation_and_is_counted() {
        use noodles_sam::alignment::record::cigar::op::Kind;
        for (label, blind) in [
            ("an interior N", {
                let mut seq = contig()[10..40].to_vec();
                seq[11] = b'N'; // position 22
                read("blind", 11, &[(Kind::Match, 30)], seq)
            }),
            ("a ref-skip", {
                // Matched 11..21, skips 22..23, matched 24..42 — the skip witnesses nothing.
                read(
                    "blind",
                    11,
                    &[(Kind::Match, 11), (Kind::Skip, 2), (Kind::Match, 19)],
                    [contig()[10..21].to_vec(), contig()[23..42].to_vec()].concat(),
                )
            }),
        ] {
            // The masking fixture's two indel reads (which open and shape the record at 20)
            // without its masked read, plus the blind one. Rebuilt rather than filtered, so
            // the list is what it says it is.
            let mut reads = masking_fixture_reads();
            reads.truncate(2);
            reads.push(blind);
            let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
            let report = dump(&fasta, &bam, &[], None);

            assert!(
                report.reads_without_observation > 0,
                "{label}: the blind read's witness is non-contiguous over the record at 20, \
                 so it must yield no row and be counted"
            );
            // The hole is at 22, inside the record at 20..=24, so no row may claim to have
            // witnessed the whole footprint on this read's behalf.
            let blind_bases = [ref_span(20, 21), ref_span(23, 24)].concat();
            assert!(
                rows_at(&report, 20)
                    .iter()
                    .all(|row| row.observed != blind_bases),
                "{label}: a row spliced the two sides of the hole together, which claims a \
                 sequence the read never witnessed as one run"
            );
        }
    }

    /// **Fixture 3 — a record widened past an already-expired read: its bases must not have
    /// grown** (spec §12.3, §13).
    ///
    /// This is the check production **cannot pass by construction**: its `widen` appends the
    /// new reference bases to every bucket, including one belonging to a read that left the
    /// active set before the widen happened and whose cursor is gone.
    ///
    /// The geometry has to be exact, and getting it exact means reckoning with ng's **real read
    /// preparation**, which left-aligns every indel before the walk sees it. The first draft of
    /// this fixture placed a deletion at 20..24 and the normalizer shifted it to 19..23,
    /// because the base before it equals the last base it deletes — so the record opened one
    /// position earlier than the test asserted and nothing at 20 existed at all. Each deletion
    /// below is therefore chosen so that `ref[anchor] != ref[anchor + len]`, which is exactly
    /// the condition under which a left-shift is not available.
    ///
    /// What the three reads do:
    ///
    /// - `opener` deletes 21..25 after matching to 20. `ref[20] = C` and `ref[25] = C`, so this
    ///   one **does** shift, to anchor 19 deleting 20..24 — deliberately, because it makes the
    ///   record 19..=24 and gives the widen something six positions wide to grow.
    /// - `short` is aligned 11..23 and gone. It folds with five of those six positions and
    ///   **expires at 24**. Its 17-base soft clip is what lets it be that short: the real
    ///   filter's 30-base minimum counts *sequenced* bases, and a read aligned over 30 could
    ///   not end inside this footprint at all.
    /// - `widener` deletes 25..34, anchored at **24** — the record's last position, so it is
    ///   inside the footprint and grows it to 19..=34. `ref[24] = A`, `ref[34] = C`: no shift.
    ///   The walker is at 24 when this happens, and `short` ended at 23, so it is **already out
    ///   of the active set** and nothing can re-fold it.
    #[test]
    fn a_locus_widened_past_an_expired_read_shows_it_observed_at_its_own_length() {
        use noodles_sam::alignment::record::cigar::op::Kind;
        let reads = vec![
            read(
                "opener",
                11,
                &[(Kind::Match, 10), (Kind::Deletion, 5), (Kind::Match, 20)],
                [contig()[10..20].to_vec(), contig()[25..45].to_vec()].concat(),
            ),
            read(
                "short",
                11,
                &[(Kind::SoftClip, 17), (Kind::Match, 13)],
                [vec![b'A'; 17], contig()[10..23].to_vec()].concat(),
            ),
            // Matched 11..24, deletes 25..34, matched 35..50.
            read(
                "widener",
                11,
                &[(Kind::Match, 14), (Kind::Deletion, 10), (Kind::Match, 16)],
                [contig()[10..24].to_vec(), contig()[34..50].to_vec()].concat(),
            ),
        ];
        let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
        let report = dump(&fasta, &bam, &[], None);

        assert!(
            report.record_widen_events > 0,
            "no record widened, so this fixture is not exercising the widen at all"
        );
        assert_eq!(
            report.end_of(19),
            34,
            "the record at 19 must have grown to the widener's deletion end, or the fixture \
             is not the one this test needs. Rows: {:?}",
            report
                .rows
                .iter()
                .map(|row| (row.start, row.end, row.read_coverage.clone()))
                .collect::<Vec<_>>(),
        );
        let short_bases = ref_span(19, 23);
        let row = row_at(&report, 19, &short_bases);
        assert_eq!(
            row.read_coverage, "observed:0+5",
            "`short` witnessed positions 19..23 of a sixteen-position footprint"
        );
        assert_eq!(
            row.observed.len(),
            5,
            "and its bases are the five it sequenced — production appends the reference bases \
             for 24..34 here, which is eleven bases it never saw"
        );
        assert!(
            rows_at(&report, 19)
                .iter()
                .all(|row| row.observed != ref_span(19, 34)),
            "some row carries the whole sixteen-base reference footprint, which no read in \
             this fixture witnessed as a run of matches"
        );
    }

    /// **Fixture 4 — an indel whose own anchor base was masked: the residual is one base, and
    /// no more** (spec §12.4).
    ///
    /// A2's rule is that nothing is written into an observation its read did not witness, and
    /// it has **exactly one recorded exception**: an indel needs an anchor base, and when the
    /// `Match` that would have emitted it is silenced the fold borrows it from `ref_seq`. That
    /// residual is deliberate and bounded, and this pins the bound — the row is that one base
    /// and nothing else, where the read's remaining matches (25..40, all masked) contribute
    /// nothing.
    #[test]
    fn an_indels_masked_anchor_borrows_exactly_one_reference_base() {
        use noodles_sam::alignment::record::cigar::op::Kind;
        let reads = vec![
            read(
                "anchor_masked",
                11,
                &[(Kind::Match, 10), (Kind::Deletion, 4), (Kind::Match, 20)],
                [contig()[10..20].to_vec(), contig()[24..44].to_vec()].concat(),
            ),
            matching_read("whole", 11, 30),
        ];
        let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
        // Masked from 20 — its own deletion's anchor — so every base it would have emitted
        // over this record is silenced and only the indel event survives.
        let report = dump(&fasta, &bam, &[("anchor_masked", 20)], None);

        let anchor_only = vec![ref_base(20)];
        let row = row_at(&report, 20, &anchor_only);
        assert_eq!(
            row.observed.len(),
            1,
            "the residual is the indel's anchor base and nothing else"
        );
        assert_eq!(
            row.observed[0],
            ref_base(20),
            "and it is the reference base at the anchor, which is where it was borrowed from"
        );
    }

    /// **Fixture 5 — two read groups on one allele: two rows summing to the single-group
    /// total** (spec §12.5, §13).
    ///
    /// Both halves of §13's claim are checked, because they fail in opposite directions. At
    /// **two** groups an allele both carry must be two rows, and their `num_obs` must sum to
    /// what one group reported — a split that lost or duplicated a read would pass a row count
    /// and fail this. At **one** group the row count must be *identical* to a run that ignores
    /// the field, which is what "free at one read group" has to mean in practice.
    #[test]
    fn two_read_groups_split_one_allele_into_two_rows_that_sum_back() {
        let base_reads = vec![
            matching_read("r0", 11, 30),
            matching_read("r1", 11, 30),
            matching_read("r2", 11, 30),
        ];
        let (_dir_one, fasta_one, bam_one) = fixture(&base_reads, &["rg0"]);
        let one_group = dump(&fasta_one, &bam_one, &[], None);

        let split_reads: Vec<_> = base_reads
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, record)| in_read_group(record, if index == 0 { "rg1" } else { "rg0" }))
            .collect();
        let (_dir_two, fasta_two, bam_two) = fixture(&split_reads, &["rg0", "rg1"]);
        let two_groups = dump(&fasta_two, &bam_two, &[], None);

        assert_eq!(
            one_group.generic_loci, two_groups.generic_loci,
            "the group is part of a row's identity and of nothing else — it cannot change \
             which loci exist"
        );
        assert_eq!(
            one_group.rows.len() * 2,
            two_groups.rows.len(),
            "every locus's single row must have become two, one per group"
        );
        assert_eq!(
            one_group.reads_complete, two_groups.reads_complete,
            "and the reads behind them are the same reads"
        );

        // Per locus: the split rows sum to the single-group row, field for field on `reads`.
        let mut per_locus: BTreeMap<u64, u32> = BTreeMap::new();
        for row in &two_groups.rows {
            *per_locus.entry(row.start).or_default() += row.reads;
        }
        for row in &one_group.rows {
            assert_eq!(
                per_locus.get(&row.start).copied(),
                Some(row.reads),
                "at locus {}, the two groups' rows must sum to the single group's {} reads",
                row.start,
                row.reads,
            );
        }
        let groups: std::collections::BTreeSet<u32> =
            two_groups.rows.iter().map(|row| row.read_group).collect();
        assert_eq!(
            groups.len(),
            2,
            "both groups must appear, or the fixture is not split at all"
        );
    }

    /// **A chunk wider than the contig is one piece, not a panic or a hang.**
    ///
    /// `PVC_GENERIC_REGION_CHUNK_BP` reaches `chunk_typed_regions` straight from the
    /// environment, where `at + chunk - 1` panicked with "attempt to add with overflow" in
    /// debug and **wrapped** in release — and a wrapped `piece_end` lands below `at`, so the
    /// `while at <= end` loop would never terminate. `saturating_add` gives the honest answer
    /// instead: one piece covering the region.
    ///
    /// Asserted as **equality with the unchunked walk**, because that is the property a caller
    /// of this knob relies on — the output is invariant to the region grain — and it is what a
    /// silently-wrapped bound would break.
    ///
    /// Mutation: `at.saturating_add(chunk - 1)` → `at + chunk - 1` panics here in debug.
    #[test]
    fn a_chunk_wider_than_the_contig_is_one_piece() {
        let reads = masking_fixture_reads();
        let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
        let whole = dump(&fasta, &bam, &[("masked", 22)], None);
        let huge = dump(&fasta, &bam, &[("masked", 22)], Some(u64::MAX));
        assert_eq!(
            huge.generic_regions, whole.generic_regions,
            "a chunk wider than the contig must not split it",
        );
        assert_eq!(
            huge.render(),
            whole.render(),
            "the dump is invariant to the region grain, and a chunk wider than the contig is \
             the degenerate case of that",
        );
    }

    /// **Fixture 6 — a deletion at a region boundary: its support matches a single-region
    /// walk** (spec §12.6, §2).
    ///
    /// The halo is why this holds: a region's query runs to `end + max_record_span`, so a
    /// record anchored inside the region still sees the reads that fold into it from beyond the
    /// boundary. Without it the record at 20 would lose everything past position 20, and its
    /// support would depend on where a caller happened to cut the regions — the one thing a
    /// region-parallel caller cannot tolerate.
    ///
    /// Checked as an equality against the *same fixture* walked as one region, which is the
    /// only comparison that isolates the boundary: the loci must concatenate coordinate-sorted
    /// and duplicate-free, with no gap at the join (spec §13).
    /// # The read that starts *past* the boundary is the whole test
    ///
    /// A first draft chunked the masking fixture, whose reads all start at 11 — inside the first
    /// piece — and **passed with the halo removed entirely**, because a read overlapping the
    /// unwidened query span is fetched either way. The halo only matters for a read that begins
    /// **after** the region's end and still folds into a record anchored inside it, so this
    /// fixture adds one: `beyond` starts at 21, one base past the 20 bp boundary, and
    /// contributes positions 21..24 of the record anchored at 20.
    #[test]
    fn a_deletion_across_a_region_boundary_keeps_the_support_a_single_region_walk_gives() {
        let mut reads = masking_fixture_reads();
        reads.push(matching_read("beyond", 21, 30));
        let (_dir, fasta, bam) = fixture(&reads, &["rg0"]);
        let whole = dump(&fasta, &bam, &[("masked", 22)], None);
        assert!(
            whole
                .rows
                .iter()
                .any(|row| row.start == 20 && row.read_coverage == "observed:1+4"),
            "`beyond` must fold into the record at 20 from outside it — without that row this \
             fixture cannot tell a halo from no halo. Rows at 20: {:?}",
            rows_at(&whole, 20)
                .iter()
                .map(|row| row.read_coverage.clone())
                .collect::<Vec<_>>(),
        );
        // 20 bp pieces put a boundary at 20|21, inside the footprint of the record anchored at
        // 20 — and another at 40|41, inside no record, so both cases are on the walk.
        let chunked = dump(&fasta, &bam, &[("masked", 22)], Some(20));

        assert!(
            chunked.generic_regions > whole.generic_regions,
            "the chunked run must actually have more regions ({} against {}), or this test \
             compares one configuration with itself",
            chunked.generic_regions,
            whole.generic_regions,
        );
        assert_eq!(
            chunked.rows, whole.rows,
            "the rows must not depend on where the regions were cut"
        );
        assert_eq!(
            chunked.generic_loci, whole.generic_loci,
            "and neither must the loci: no duplicate at a join, no gap either"
        );
        assert!(
            chunked.records_outside_region > 0,
            "no record was dropped by the region clamp, so the halo was never walked past a \
             boundary and this fixture proves nothing about it"
        );
        let anchors: Vec<u64> = {
            let mut anchors: Vec<u64> = chunked.rows.iter().map(|row| row.start).collect();
            anchors.dedup();
            anchors
        };
        let mut sorted = anchors.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            anchors, sorted,
            "the concatenated stream must be coordinate-sorted and duplicate-free"
        );
    }
}

#[cfg(test)]
impl DumpReport {
    /// The `end` of the locus anchored at `start` — the footprint's last position, which is
    /// the field a `PileupRecord` has no counterpart for.
    fn end_of(&self, start: u64) -> u64 {
        let ends: std::collections::BTreeSet<u64> = self
            .rows
            .iter()
            .filter(|row| row.start == start)
            .map(|row| row.end)
            .collect();
        assert_eq!(
            ends.len(),
            1,
            "the rows of the locus at {start} disagree about its end: {ends:?}"
        );
        *ends.iter().next().expect("one end")
    }
}
