//! **A repeatable harness for ng's generic locus generator** — the thing plan 3
//! finished without: both of its measurements came from throwaway probes that were
//! deleted after use, so every later question started from scratch (arch
//! `locus_generation_pileup.md`, *Test & bench shape*).
//!
//! This file is **infrastructure, not a verdict**. It answers "how do I measure ng's
//! walker again next week", and it deliberately tunes nothing: performance is parked
//! on the owner's instruction until correctness is done.
//!
//! # What it drives, and why not the walker directly
//!
//! `pileup_walker_scaling` benches **production's** walker by calling
//! `pileup::walker::run` with a `Vec<PreparedRead>` it builds in memory. The obvious
//! mirror — call ng's copy of `run` the same way — is not available and should not
//! be: ng's walker vocabulary is `pub(crate)`, and `run` was demoted from the public
//! API precisely because a caller reaching it bypasses
//! `PileupGeneratorConfig::check()` and its record-span ceiling.
//!
//! So this bench enters where a real caller does: [`PileupGenerator`], over a
//! [`SampleReads`] opened on a real (synthetic, temporary) BAM, with ng's real
//! [`LeftAlignPreparer`]. That makes the measurement **wider** than production's
//! walker bench — it includes BAM decode, read filtering and preparation — and the
//! two numbers are therefore not comparable to each other. They are comparable to
//! *themselves* across commits, which is what a regression harness is for.
//!
//! # The two axes
//!
//! - **Coverage** at a fixed span. Spec §7's claim is that what the generator holds
//!   is bounded by depth, not by region length, so depth is the axis that claim
//!   lives on.
//! - **Region grain** at fixed coverage and span: the same base pairs walked as one
//!   region, as ten, as a hundred. The halo means reads near a boundary are prepared
//!   twice, and the grain is the caller's choice in the real pipeline, so a harness
//!   that cannot vary it cannot answer questions about it. **This is the axis D1
//!   flattened**: a region used to re-open a read query and re-fetch the reference, and
//!   now it repositions a cursor that stayed open and hands back reads it already
//!   decoded, so the slope here is what the change is for.
//!
//! Both fixtures are deterministic: the contig comes from a fixed-seed LCG and the
//! reads are laid down by an integer rule, so two runs on one commit compare.
//!
//! ```text
//! cargo bench --bench ng_generic_pileup_perf
//! ```

use std::hint::black_box;
use std::sync::Arc;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// The fixture builders, shared with `examples/ng_generic_walk_probe.rs` so a bench
/// number and that probe's assertions describe the same reads.
#[path = "../examples/shared/synthetic_alignment.rs"]
mod synthetic_alignment;
use synthetic_alignment::{SyntheticGeometry, SyntheticSample};

/// Reference span of every fixture, in base pairs. Big enough that per-region
/// constants do not dominate — since D1 those are a cursor reposition rather than a
/// fresh query and a first reference fetch — small enough that writing the BAM stays
/// a setup cost measured in tenths of a second.
const SPAN: u64 = 100_000;

/// Read length. Illumina-shaped, and the length production's own walker bench sweeps
/// around.
const READ_LEN: u64 = 150;

/// The fixture's shape at one depth. Span and read length are this bench's; only
/// coverage moves along its first axis.
fn geometry(coverage: u64) -> SyntheticGeometry {
    SyntheticGeometry {
        span: SPAN,
        read_len: READ_LEN,
        coverage,
    }
}

// ---------------------------------------------------------------------------
// Driving the generator
// ---------------------------------------------------------------------------

/// The generator plus the reads it is driven against.
///
/// They are returned together because `next_locus` borrows the `SampleReads`, and a
/// bench body that rebuilt either between iterations would be timing the setup.
struct Driver {
    generator: PileupGenerator<WindowedRefSeq, LeftAlignPreparer<WindowedRefSeq>>,
    reads: SampleReads,
}

/// Open the fixture and build a generator over it — the same construction
/// `ng_generic_loci_dump` uses, including the one thing that is not obvious: the
/// `.fai` is parsed once and shared, because a fresh `WindowedRefSeq::new` re-parses
/// it on its first fetch. Since D1 the factory is called once per file per
/// *chromosome* rather than once per region, so this matters far less than it did —
/// it is kept because it is the shape the dump tool uses and this bench exists to
/// measure that shape.
fn driver(fixture: &SyntheticSample) -> Driver {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) =
        read_reference_verifying_or_creating_fai(&cache, fixture.fasta.clone()).unwrap();
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&fixture.fasta).unwrap();

    // One reference for every file this run opens (the CRAM repository-sharing shape `main`
    // moved to while this branch was in flight).
    let reference = OpenReference::new(info);
    let reads = SampleReads::open_only_sample(
        std::slice::from_ref(&fixture.bam),
        &reference,
        ReadFilterConfig::default(),
        true,
    )
    .unwrap();

    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fixture.fasta.clone(),
        contigs.clone(),
        index.clone(),
    ));

    // Single-threaded and file-backed, exactly as in the dump tool — see its
    // `run_dump` for why the `Arc` is the right signature even here.
    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new is generic over the accessor and takes Arc; this bench's accessor is file-backed and single-threaded"
    )]
    let reference = Arc::new(WindowedRefSeq::with_shared_index(
        fixture.fasta.clone(),
        contigs.clone(),
        index.clone(),
    ));
    // Not boxed here any more: `PileupGenerator::new` boxes the factory itself, which is
    // what keeps it off the type this `Driver` names (arch §3.6).
    let make_reference = {
        let fasta = fixture.fasta.clone();
        move || WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
    };
    let generator = PileupGenerator::new(
        reference,
        make_reference,
        preparer,
        PileupGeneratorConfig::default(),
    )
    .unwrap();

    // The reference verification runs on its own thread and its error is only
    // reachable through the handle; a bench that dropped it would print the
    // library's "went unobserved" warning on every run. `None` when the `.fai`
    // was created rather than read — there is nothing to verify against.
    if let Some(verify) = verify {
        verify.join().unwrap();
    }

    Driver { generator, reads }
}

/// Split `[1, SPAN]` into pieces of `grain` base pairs — the region list a caller
/// would hand in.
fn regions(grain: u64) -> Vec<GenomeRegion> {
    let mut out = Vec::new();
    let mut at = 1u64;
    while at <= SPAN {
        let end = at.saturating_add(grain - 1).min(SPAN);
        out.push(GenomeRegion {
            contig: ContigId(0),
            start: Position(at),
            end: Position(end),
        });
        at = end + 1;
    }
    out
}

/// The timed body: walk every region and count the loci.
///
/// Counting rather than discarding, and returned rather than dropped, so the
/// optimiser cannot delete the walk — and so a bench run that silently generated
/// nothing shows up as a zero rather than as a very fast result.
fn walk(driver: &mut Driver, regions: &[GenomeRegion]) -> u64 {
    let mut loci = 0u64;
    for region in regions {
        driver.generator.begin_segment(*region);
        while driver
            .generator
            .next_locus(&driver.reads)
            .expect("the generator walks the fixture without error")
            .is_some()
        {
            loci += 1;
        }
    }
    loci
}

/// One untimed walk before each timed set, asserting the fixture produces the loci
/// **and** the depth it claims.
///
/// **A bench measuring a walk that generates nothing looks like a fast walk**, and
/// nothing else here would notice: `criterion` reports a time, not an answer. Two
/// things are pinned, because either can fail silently and they fail differently:
///
/// - **Every position yields a locus.** The reads tile `[1, SPAN]` with pure-Match
///   CIGARs, so each reference position is covered and yields exactly one locus, and
///   the regions tile disjointly — splitting the span moves loci between regions
///   without changing the total. This caught the first draft, whose read spacing left
///   the last fifteen positions uncovered.
/// - **The reads reach the walk.** The locus count alone cannot see this: one read
///   per position would still cover the span, so a read filter dropping nine reads in
///   ten would leave the coverage *axis* meaning nothing while every locus was still
///   emitted. `reads_admitted` is read as a delta around this walk, and the floor is
///   the number written, not an equality: a walk split into regions admits a read
///   again in every region whose halo reaches it.
fn check_the_walk_covers_the_span(
    driver: &mut Driver,
    regions: &[GenomeRegion],
    geometry: SyntheticGeometry,
) {
    let admitted_before = driver.generator.counts().reads_admitted;
    let loci = walk(driver, regions);
    let admitted = driver.generator.counts().reads_admitted - admitted_before;
    assert_eq!(
        loci,
        geometry.span,
        "the fixture's reads tile the whole span, so every one of its {} positions \
         must yield exactly one locus across the {} region(s) walked",
        geometry.span,
        regions.len(),
    );
    assert!(
        admitted >= geometry.num_reads(),
        "the depth axis is only meaningful if the reads reach the walk: {admitted} admitted \
         against {} written at {}×",
        geometry.num_reads(),
        geometry.coverage,
    );
}

// ---------------------------------------------------------------------------
// The benches
// ---------------------------------------------------------------------------

/// Cost against depth, at one region for the whole span.
fn bench_coverage(c: &mut Criterion) {
    let mut group = c.benchmark_group("ng_generic_pileup_coverage");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    let whole = regions(SPAN);
    for &coverage in &[10u64, 30, 100] {
        let geometry = geometry(coverage);
        let fixture = SyntheticSample::build(geometry);
        let mut driver = driver(&fixture);
        check_the_walk_covers_the_span(&mut driver, &whole, geometry);
        group.bench_with_input(BenchmarkId::from_parameter(coverage), &coverage, |b, _| {
            b.iter(|| black_box(walk(&mut driver, &whole)))
        });
    }
    group.finish();
}

/// Cost against region grain, at one depth. The same base pairs and the same reads
/// every time — only the number of `begin_segment` calls changes.
fn bench_region_grain(c: &mut Criterion) {
    let mut group = c.benchmark_group("ng_generic_pileup_region_grain");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    const GRAIN_COVERAGE: u64 = 30;
    let geometry = geometry(GRAIN_COVERAGE);
    let fixture = SyntheticSample::build(geometry);
    let mut driver = driver(&fixture);
    for &grain in &[SPAN, 10_000, 1_000] {
        let regions = regions(grain);
        check_the_walk_covers_the_span(&mut driver, &regions, geometry);
        group.bench_with_input(BenchmarkId::from_parameter(grain), &grain, |b, _| {
            b.iter(|| black_box(walk(&mut driver, &regions)));
        });
    }
    group.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(3))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_coverage, bench_region_grain
}

criterion_main!(benches);
