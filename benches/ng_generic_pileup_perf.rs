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
//!   region, as ten, as a hundred. Each region re-opens a read query and re-fetches
//!   reference; the halo means reads near a boundary are prepared twice. The grain
//!   is the caller's choice in the real pipeline, so a harness that cannot vary it
//!   cannot answer questions about it.
//!
//! Both fixtures are deterministic: the contig comes from a fixed-seed LCG and the
//! reads are laid down by an integer rule, so two runs on one commit compare.
//!
//! ```text
//! cargo bench --bench ng_generic_pileup_perf
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use tempfile::TempDir;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// Reference span of every fixture, in base pairs. Big enough that per-region
/// constants (opening the query, the first reference fetch) do not dominate, small
/// enough that writing the BAM stays a setup cost measured in tenths of a second.
const SPAN: u64 = 100_000;

/// Read length. Illumina-shaped, and the length production's own walker bench sweeps
/// around.
const READ_LEN: u64 = 150;

// ---------------------------------------------------------------------------
// The fixture: a contig, its reads, and the files they live in
// ---------------------------------------------------------------------------

/// `len` deterministic bases from a fixed-seed LCG.
///
/// Not a real sequence and not meant to be: this bench drives the generator
/// directly, so nothing here types regions or looks for repeats. What the bases must
/// do is vary — a homopolymer contig would give the left-aligner and the fold an
/// unrepresentatively easy time — and be the same on every run, which a seeded
/// generator gives and `rand` would not.
fn synthetic_contig(len: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            BASES[(state >> 33) as usize % BASES.len()]
        })
        .collect()
}

/// Reads written for one coverage — the fixture's own definition of depth, and the
/// floor the admission check holds the walk to.
fn num_reads(coverage: u64) -> u64 {
    SPAN * coverage / READ_LEN
}

/// The base that is not this one — used to plant a mismatch.
fn other_base(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// Reads tiling `[1, SPAN]` at approximately `coverage`×, each carrying **one**
/// mismatch.
///
/// The mismatch is what makes the fold do the work the generator exists for: a run
/// of reads that match the reference everywhere produces one row per locus and never
/// exercises a second bucket, a split row or a chain id. One per read, at a position
/// that moves with the read index, keeps every locus's row count small and realistic
/// rather than turning the fixture into a variant caller stress test.
///
/// Starts ascend, which is both what a coordinate-sorted BAM needs and what the
/// walker requires of its input.
fn synthetic_reads(contig: &[u8], coverage: u64) -> Vec<noodles_sam::alignment::RecordBuf> {
    use noodles_core::Position as CorePosition;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    let last_start = SPAN - READ_LEN + 1;
    let num_reads = num_reads(coverage);
    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        // Even spacing over the startable positions; several reads share a start
        // once `num_reads` exceeds them, which is what real depth looks like.
        //
        // The divisor is `num_reads - 1`, not `num_reads`, so the last read starts
        // at `last_start` and its final base is `SPAN`. With `num_reads` the walk
        // stopped fifteen positions short of the span at 10× — invisible to a timing
        // number, which is why the coverage check downstream asserts on the loci.
        let start = 1 + (i * (last_start - 1) / (num_reads - 1).max(1));
        let mut seq = contig[start as usize - 1..(start + READ_LEN - 1) as usize].to_vec();
        let at = (i % READ_LEN) as usize;
        seq[at] = other_base(seq[at]);
        reads.push(
            noodles_sam::alignment::RecordBuf::builder()
                .set_name(format!("r{i}").as_bytes())
                .set_reference_sequence_id(0)
                .set_flags(Flags::empty())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_alignment_start(CorePosition::try_from(start as usize).unwrap())
                .set_cigar(
                    [Op::new(Kind::Match, READ_LEN as usize)]
                        .into_iter()
                        .collect(),
                )
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![40u8; READ_LEN as usize]))
                .build(),
        );
    }
    reads
}

/// Write `chr1` as a one-line FASTA. The `.fai` is created by the reference read.
fn write_fasta(path: &Path, bases: &[u8]) {
    use std::io::Write as _;
    let mut file = std::fs::File::create(path).unwrap();
    writeln!(file, ">chr1").unwrap();
    file.write_all(bases).unwrap();
    writeln!(file).unwrap();
}

/// Write a coordinate-sorted single-contig, single-read-group BAM.
fn write_bam(path: &Path, contig_len: usize, reads: &[noodles_sam::alignment::RecordBuf]) {
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
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from("sample0"));
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(b"chr1".to_vec(), sq)
        .add_read_group(b"rg0".to_vec(), rg)
        .build();

    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for record in reads {
        writer.write_alignment_record(&header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

/// A reference and a sample's reads on disk, held together so the temporary
/// directory outlives them.
struct Fixture {
    _dir: TempDir,
    fasta: PathBuf,
    bam: PathBuf,
}

/// Build the files for one coverage. Untimed by construction — every bench body
/// takes this already built.
fn fixture(coverage: u64) -> Fixture {
    let bases = synthetic_contig(SPAN as usize);
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let bam = dir.path().join("sample.bam");
    write_fasta(&fasta, &bases);
    write_bam(&bam, bases.len(), &synthetic_reads(&bases, coverage));
    Fixture {
        _dir: dir,
        fasta,
        bam,
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
    generator: PileupGenerator<
        WindowedRefSeq,
        Box<dyn FnMut() -> WindowedRefSeq>,
        LeftAlignPreparer<WindowedRefSeq>,
    >,
    reads: SampleReads,
}

/// Open the fixture and build a generator over it — the same construction
/// `ng_generic_loci_dump` uses, including the one thing that is not obvious: the
/// `.fai` is parsed once and shared, because a fresh `WindowedRefSeq::new` per region
/// re-parses it and more than doubles a region's cost.
fn driver(fixture: &Fixture) -> Driver {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) =
        read_reference_verifying_or_creating_fai(&cache, fixture.fasta.clone()).unwrap();
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&fixture.fasta).unwrap();

    let reads = SampleReads::open_only_sample(
        std::slice::from_ref(&fixture.bam),
        &info,
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
    let make_reference: Box<dyn FnMut() -> WindowedRefSeq> = {
        let fasta = fixture.fasta.clone();
        Box::new(move || {
            WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
        })
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
fn check_the_walk_covers_the_span(driver: &mut Driver, regions: &[GenomeRegion], coverage: u64) {
    let admitted_before = driver.generator.counts().reads_admitted;
    let loci = walk(driver, regions);
    let admitted = driver.generator.counts().reads_admitted - admitted_before;
    assert_eq!(
        loci,
        SPAN,
        "the fixture's reads tile the whole span, so every one of its {SPAN} positions \
         must yield exactly one locus across the {} region(s) walked",
        regions.len(),
    );
    assert!(
        admitted >= num_reads(coverage),
        "the depth axis is only meaningful if the reads reach the walk: {admitted} admitted \
         against {} written at {coverage}×",
        num_reads(coverage),
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
        let fixture = fixture(coverage);
        let mut driver = driver(&fixture);
        check_the_walk_covers_the_span(&mut driver, &whole, coverage);
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
    let fixture = fixture(GRAIN_COVERAGE);
    let mut driver = driver(&fixture);
    for &grain in &[SPAN, 10_000, 1_000] {
        let regions = regions(grain);
        check_the_walk_covers_the_span(&mut driver, &regions, GRAIN_COVERAGE);
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
