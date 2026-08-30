//! ng's psp store: what a walk costs, and what a write costs.
//!
//! **The store had no benchmark at all before this file.** `benches/psp_reader_perf.rs` and
//! `benches/psp_writer_perf.rs` measure *production's* `src/psp/`, from the May 2026 review; ng's
//! `src/ng/psp/` had only milestone harnesses under `examples/`, which measure ratios (the head
//! skip, H5) and resident memory (H4) rather than throughput. A profile needs something
//! reproducible under it, and this is it.
//!
//! # What is timed
//!
//! Four workloads, each over a store written once outside the timed region:
//!
//! - `ng_psp/walk_full/*` — [`PspReader::records`], every body built. The shape a caller that
//!   wants all the evidence reads in.
//! - `ng_psp/walk_heads/*` — the same walk with a predicate that declines every record, so no body
//!   is ever built and only the head is decoded: the position offset, the span, the
//!   non-reference read count, the body length and the chain ids' live-set changes. This is the
//!   floor a walk cannot go below, because the bytes still arrive and the changes are still parsed
//!   (spec `psp_file_format.md` §4.3).
//! - `ng_psp/walk_one_in_100/*` — one record in a hundred built. The cohort's first pass
//!   (`cohort_merge.md`: roughly one position in a hundred varies at the measured corner).
//! - `ng_psp/write/*` — [`PspWriter::push`] over the same records plus the `finish` that writes the
//!   index and the footer. Compression is most of it; that is the point of measuring it.
//!
//! # The two corpora, and why the depth axis is here
//!
//! `CLAUDE.md` requires every measurement to say which corner of the range it was taken on, and
//! for this store the axis that matters is read depth: the chain ids' live-set changes ride in the
//! **head**, so they are paid by a walk that skips every body, and they grow with depth — spec
//! `psp_record_encoding.md` §6 measures 0.432 bytes a position at 11.4 reads and 6.42 at 293. So:
//!
//! - `shallow` — 10 reads a position, which is about the tomato panel's 10.3.
//! - `deep` — 280 reads a position, which is about the HG002 corpus's 280.0.
//!
//! **These are synthesised, not real records**, so that the bench is reproducible from the
//! repository alone and does not need a corpus under `tmp/`. What they reproduce is the *shape*:
//! one record a position, a live set of about `depth` reads sliding forward as 150-base reads
//! start and end, one observation carrying nearly all of them, and a second observation at one
//! record in a hundred. What they do **not** reproduce is a real sample's variety of witnesses,
//! read groups and locus kinds; for record-for-record fidelity against a real corpus, see
//! `examples/ng_psp_parity.rs`.
//!
//! ⚠ **The heads here are ng-shaped, not production-shaped.** The stores under `tmp/` that the
//! milestone harnesses walk are built from a production `.psp`, which names about 3.4 % of the
//! reads ng will name (the owner's ruling of 2026-08-17), so their heads carry a small fraction of
//! the chain-id changes ng's will. This corpus names one identifier per read, which is what ng
//! will write — so a profile taken here weights the head the way the shipped run will, and a
//! walk timed here is **not** comparable with one timed on a `tmp/` store.

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc ...).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::collections::BTreeMap;
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use pop_var_caller::ng::psp::{
    ContigIdentity, FORMAT_VERSION, Header, Manifest, PspReader, PspWriter, ReferenceIdentity,
    WriterProvenance,
};
use pop_var_caller::ng::types::{
    Bp, ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError,
};
use pop_var_caller::pileup_record::ChainId;

/// How many reference positions a read covers. **150, the Illumina read this project's corpora
/// are** — it is what sets how fast the live set turns over: `depth / 150` reads start at each
/// position, and that arrival rate is what the head has to encode.
const READ_LENGTH_POSITIONS: u64 = 150;

/// One record in this many carries a second observation, so the residual derivation has something
/// to subtract. **A hundred**, matching `cohort_merge.md`'s measured corner: roughly one position
/// in a hundred varies.
const SECOND_OBSERVATION_ONE_IN: u64 = 100;

/// The trailer these stores carry. Opaque to the format (spec §3.4); a marker, not data.
const THE_TRAILER: &[u8] = b"ng_psp_perf";

/// One corner of the depth range, and how much of it to build.
///
/// **The record count and the block grid move together**, so that both corpora cut the same
/// number of blocks: a walk that never crosses a boundary never pays for the live set being
/// restated whole, which is a real per-block cost and would be missing from the deep arm
/// otherwise. The deep corpus is shorter because it holds about `depth` identifiers a record —
/// 280 of them at eight bytes is what decides how much memory the *bench* needs, not the store.
struct CorpusShape {
    name: &'static str,
    records: u64,
    reads_a_position: u32,
    genomic_block_size_bp: u64,
}

/// The tomato panel's corner: 10.3 reads a record, measured (`ng_psp_h5_2026-08-30.md`).
const SHALLOW: CorpusShape = CorpusShape {
    name: "shallow_10_reads",
    records: 300_000,
    reads_a_position: 10,
    genomic_block_size_bp: 100_000,
};

/// The HG002 corner: 280.0 reads a record, measured (`ng_psp_h5_2026-08-30.md`).
const DEEP: CorpusShape = CorpusShape {
    name: "deep_280_reads",
    records: 60_000,
    reads_a_position: 280,
    genomic_block_size_bp: 20_000,
};

/// The write workloads use fewer records than the walk ones: zstd at level 9 is most of a write,
/// and a criterion sample that compresses 300,000 records takes long enough to make the bench
/// unusable in a review loop.
const RECORDS_WRITTEN: u64 = 60_000;

// ---------------------------------------------------------------------
// The corpus
// ---------------------------------------------------------------------

/// The identifiers live at `position`: a sliding window of about `depth` reads.
///
/// `depth / READ_LENGTH_POSITIONS` reads start at each position and one ends, so the set turns
/// over completely every 150 positions — which is what puts arrivals and departures in every
/// record's head at 280 reads and in one record in fifteen at 10.
fn live_at(position: u64, depth: u64) -> std::ops::RangeInclusive<u64> {
    let started = position * depth / READ_LENGTH_POSITIONS;
    let oldest_still_live = started.saturating_sub(depth.saturating_sub(1));
    oldest_still_live..=started
}

/// One record: one reference position, one observation holding nearly every live read, and — at
/// one position in a hundred — a second observation holding the rest.
///
/// **`num_obs` equals the identifier count** so that the writer's residual derivation is allowed:
/// it refuses to derive a list whose length does not sit between half the read count and the read
/// count (`record.rs::check_a_read_list_against_its_read_count`), and a record it refuses stores
/// every list instead, which is not the shape a real one has.
fn a_record(position: u64, depth: u64) -> SampleLocusObservations {
    let live = live_at(position, depth);
    let ids: Vec<ChainId> = live.collect();
    let split = if position.is_multiple_of(SECOND_OBSERVATION_ONE_IN) && ids.len() >= 4 {
        ids.len() - ids.len() / 4
    } else {
        ids.len()
    };
    let (reference_reads, alternative_reads) = ids.split_at(split);

    let mut observations = Vec::with_capacity(2);
    observations.push(an_observation(b"A", reference_reads, 0));
    if !alternative_reads.is_empty() {
        observations.push(an_observation(b"C", alternative_reads, 1));
    }
    SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(0),
            start: Position(position),
            end: Position(position),
        },
        reference_bases: Box::from(&b"A"[..]),
        observations,
        reads_without_observation: (position % 7) as u32,
        reads_discarded_by_cap: (position % 3) as u32,
        kind: LocusKind::Generic,
    }
}

/// One observation, its support scaled off its read count so that no field is a constant a codec
/// could invent.
fn an_observation(bases: &[u8], reads: &[ChainId], group: u32) -> SequenceObservation {
    let num_obs = reads.len() as u32;
    SequenceObservation {
        bases: Box::from(bases),
        read_witness: ReadWitness::Complete,
        read_group: ReadGroupId(group),
        num_obs,
        num_fwd: num_obs / 2,
        q_sum: SummedLogError::from_steps(-i64::from(num_obs) * 4_096),
        mapq_sum: num_obs * 60,
        mapq_sum_sq: u64::from(num_obs) * 3_600,
        placed_left: num_obs / 3,
        chain_ids: reads.to_vec(),
    }
}

/// Every record of a corpus, built once.
///
/// **Materialised rather than generated inside the timed region**, because the write workload
/// times `push` and building a record is not part of what a writer costs.
fn a_corpus(shape: &CorpusShape, records: u64) -> Vec<SampleLocusObservations> {
    let depth = u64::from(shape.reads_a_position);
    (0..records).map(|at| a_record(at, depth)).collect()
}

/// The header these stores carry: one contig long enough for the corpus, and the grid this shape
/// cuts blocks on.
fn a_header(shape: &CorpusShape, records: u64) -> Header {
    Header {
        format_version: FORMAT_VERSION,
        sample: "bench".to_string(),
        reference: ReferenceIdentity {
            name: "bench.fa".to_string(),
            md5: None,
        },
        contigs: vec![ContigIdentity {
            name: "chr1".to_string(),
            length: records + 1_024,
            md5: None,
        }],
        writer: WriterProvenance {
            tool: "ng_psp_perf".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "bench".to_string(),
            input_alignments: vec!["bench.bam".to_string()],
            input_reference: "bench.fa".to_string(),
            command_line: "ng_psp_perf".to_string(),
            parameters: BTreeMap::new(),
            created: "2026-08-30T00:00:00Z"
                .parse()
                .expect("a valid RFC 3339 stamp"),
        },
        manifest: Manifest {
            genomic_block_size_bp: Bp(shape.genomic_block_size_bp),
            ..Manifest::as_this_build_writes_it()
        },
    }
}

/// Write a store, and report how many blocks it cut — quoted in the bench's own output so a
/// reading can say what block structure it was taken over.
fn write_a_store(path: &Path, header: Header, records: &[SampleLocusObservations]) -> u64 {
    let mut writer = PspWriter::create(path, header).expect("create the store");
    for record in records {
        writer.push(record).expect("push a record");
    }
    writer.finish(THE_TRAILER).expect("finish the store").blocks
}

// ---------------------------------------------------------------------
// The workloads
// ---------------------------------------------------------------------

/// Walk every record, building every body.
fn walk_full(path: &Path) -> u64 {
    let mut reader = PspReader::open(path).expect("open the store");
    let mut built = 0u64;
    for record in reader.records().expect("start the walk") {
        let record = record.expect("a record");
        built += u64::from(record.record.is_some());
        black_box(&record);
    }
    built
}

/// Walk every record, building one body in `keep_one_in` — `u64::MAX` for none of them.
fn walk_building_one_in(path: &Path, keep_one_in: u64) -> u64 {
    let mut reader = PspReader::open(path).expect("open the store");
    let mut seen = 0u64;
    let mut built = 0u64;
    let walk = reader
        .records()
        .expect("start the walk")
        .building_only_where(|_head| {
            let wanted = keep_one_in != u64::MAX && seen.is_multiple_of(keep_one_in);
            seen += 1;
            wanted
        });
    for record in walk {
        let record = record.expect("a record");
        built += u64::from(record.record.is_some());
        black_box(&record);
    }
    built
}

fn ng_psp_walks(c: &mut Criterion) {
    let scratch = tempfile::tempdir().expect("a scratch directory");
    for shape in [&SHALLOW, &DEEP] {
        let records = a_corpus(shape, shape.records);
        let path = scratch.path().join(format!("{}.ngpsp", shape.name));
        let blocks = write_a_store(&path, a_header(shape, shape.records), &records);
        // Freed before anything is timed: the corpus is the writer's input, and holding 280
        // identifiers a record beside the reader would price the allocator against a heap the
        // shipped run never has.
        drop(records);
        assert!(blocks >= 3, "{} cut {blocks} blocks", shape.name);

        let mut group = c.benchmark_group("ng_psp");
        group.throughput(Throughput::Elements(shape.records));
        group.measurement_time(Duration::from_secs(10));

        group.bench_function(format!("walk_full/{}", shape.name), |b| {
            b.iter(|| black_box(walk_full(&path)));
        });
        group.bench_function(format!("walk_heads/{}", shape.name), |b| {
            b.iter(|| black_box(walk_building_one_in(&path, u64::MAX)));
        });
        group.bench_function(format!("walk_one_in_100/{}", shape.name), |b| {
            b.iter(|| black_box(walk_building_one_in(&path, SECOND_OBSERVATION_ONE_IN)));
        });

        group.finish();
    }
}

fn ng_psp_writes(c: &mut Criterion) {
    let scratch = tempfile::tempdir().expect("a scratch directory");
    for shape in [&SHALLOW, &DEEP] {
        let records = a_corpus(shape, RECORDS_WRITTEN);
        let path = scratch.path().join(format!("write_{}.ngpsp", shape.name));

        let mut group = c.benchmark_group("ng_psp");
        group.throughput(Throughput::Elements(RECORDS_WRITTEN));
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(10));
        group.bench_function(format!("write/{}", shape.name), |b| {
            b.iter(|| {
                black_box(write_a_store(
                    &path,
                    a_header(shape, RECORDS_WRITTEN),
                    &records,
                ))
            });
        });
        group.finish();
    }
}

criterion_group!(benches, ng_psp_walks, ng_psp_writes);
criterion_main!(benches);
