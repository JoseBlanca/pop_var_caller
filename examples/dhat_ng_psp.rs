//! Count the heap allocations ng's psp store makes, per record, on both sides.
//!
//! **The counterpart the store did not have.** `examples/dhat_psp_reader.rs` and
//! `dhat_psp_writer.rs` profile *production's* `src/psp/`; `src/ng/psp/` had no allocation
//! oracle at all, so an allocations review could only gate on wall time — which moves 20–25 %
//! run to run — instead of on a count that is identical every run.
//!
//! Run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_ng_psp \
//!     --no-default-features --features dhat-heap
//! ```
//!
//! `--no-default-features` is what takes mimalloc out of the way; the crate's mimalloc
//! declaration is already gated on `not(feature = "dhat-heap")`, so `--all-features` builds too,
//! but a run wants only one global allocator's bookkeeping.
//!
//! # What it prints
//!
//! One line a phase, each a delta over the phase alone — the corpus is built and the store
//! written *outside* the region the numbers cover, so what is reported is the walk's own
//! allocations, or the write's.
//!
//! ```text
//! phase                 records   blocks   bytes   blocks/record  bytes/record
//! ```
//!
//! `blocks` is dhat's word for allocations. **It is the number to gate an allocation fix on**:
//! the same code path allocates the same number of times on every run.
//!
//! The corpus is the one `benches/ng_psp_perf.rs` synthesises, at the same two depths and with
//! every read named, which is what ng will write — **paired reads, so an identifier goes absent
//! and comes back**, the shape spec `psp_record_encoding.md` §6 measures at 83 % of identifiers on
//! the human sample and 91 % on tomato. It is deliberately smaller than the bench's (a tenth of the
//! records) because a count a record does not need a long run to be exact.
//!
//! **The corpus builder is duplicated from the bench rather than shared**, because an example
//! cannot import a bench target. The two must stay the same shape or the counts here stop
//! describing what the bench times; if you change one, change the other.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::collections::BTreeMap;
use std::path::Path;

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

const MATE_LENGTH_POSITIONS: u64 = 150;
const INNER_GAP_POSITIONS: u64 = 200;
const FRAGMENT_SPAN_POSITIONS: u64 = 2 * MATE_LENGTH_POSITIONS + INNER_GAP_POSITIONS;
const SECOND_OBSERVATION_ONE_IN: u64 = 100;
const THE_TRAILER: &[u8] = b"dhat_ng_psp";

struct CorpusShape {
    name: &'static str,
    records: u64,
    reads_a_position: u32,
    genomic_block_size_bp: u64,
}

const SHALLOW: CorpusShape = CorpusShape {
    name: "shallow_10_reads",
    records: 30_000,
    reads_a_position: 10,
    genomic_block_size_bp: 10_000,
};

const DEEP: CorpusShape = CorpusShape {
    name: "deep_280_reads",
    records: 6_000,
    reads_a_position: 280,
    genomic_block_size_bp: 2_000,
};

/// The identifiers live at `position`, ascending — two stretches, because a chain id names a read
/// pair. Mirrors `benches/ng_psp_perf.rs::live_at`; see this file's header.
fn live_at(position: u64, depth: u64, into: &mut Vec<ChainId>) {
    into.clear();
    let first_starting_after = |p: u64| p * depth / (2 * MATE_LENGTH_POSITIONS);
    let second_mate_from = first_starting_after(position.saturating_sub(FRAGMENT_SPAN_POSITIONS));
    let second_mate_to =
        first_starting_after(position.saturating_sub(MATE_LENGTH_POSITIONS + INNER_GAP_POSITIONS));
    let first_mate_from = first_starting_after(position.saturating_sub(MATE_LENGTH_POSITIONS));
    let first_mate_to = first_starting_after(position);
    into.extend(second_mate_from..second_mate_to);
    into.extend(first_mate_from..first_mate_to);
}

fn a_record(position: u64, depth: u64, live: &mut Vec<ChainId>) -> SampleLocusObservations {
    live_at(position, depth, live);
    let split = if position.is_multiple_of(SECOND_OBSERVATION_ONE_IN) && live.len() >= 4 {
        live.len() - live.len() / 4
    } else {
        live.len()
    };
    let (reference_reads, alternative_reads) = live.split_at(split);

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
        reference_bases: b"A".to_vec(),
        observations,
        reads_without_observation: (position % 7) as u32,
        reads_discarded_by_cap: (position % 3) as u32,
        kind: LocusKind::Generic,
    }
}

fn an_observation(bases: &[u8], reads: &[ChainId], group: u32) -> SequenceObservation {
    let num_obs = reads.len() as u32;
    SequenceObservation {
        bases: bases.to_vec(),
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

fn a_corpus(shape: &CorpusShape) -> Vec<SampleLocusObservations> {
    let depth = u64::from(shape.reads_a_position);
    let mut live = Vec::new();
    (0..shape.records)
        .map(|at| a_record(at, depth, &mut live))
        .collect()
}

fn a_header(shape: &CorpusShape) -> Header {
    Header {
        format_version: FORMAT_VERSION,
        sample: "dhat".to_string(),
        reference: ReferenceIdentity {
            name: "dhat.fa".to_string(),
            md5: None,
        },
        contigs: vec![ContigIdentity {
            name: "chr1".to_string(),
            length: shape.records + 1_024,
            md5: None,
        }],
        writer: WriterProvenance {
            tool: "dhat_ng_psp".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "dhat".to_string(),
            input_alignments: vec!["dhat.bam".to_string()],
            input_reference: "dhat.fa".to_string(),
            command_line: "dhat_ng_psp".to_string(),
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

fn write_a_store(path: &Path, header: Header, records: &[SampleLocusObservations]) -> u64 {
    let mut writer = PspWriter::create(path, header).expect("create the store");
    for record in records {
        writer.push(record).expect("push a record");
    }
    writer.finish(THE_TRAILER).expect("finish the store").blocks
}

/// Open a store and drop it — `rounds` times. The header, the block index and the footer.
fn open_only(path: &Path, rounds: u64) -> u64 {
    let mut opened = 0u64;
    for _ in 0..rounds {
        let reader = PspReader::open(path).expect("open the store");
        opened += 1;
        std::hint::black_box(&reader);
    }
    opened
}

/// Start a walk on an already-open reader, take one record, drop the walk — `rounds` times.
///
/// **What a second walk over the same open sample costs**, which is the number that decides
/// whether a caller reading region by region should be handed a reader that keeps its buffers.
fn start_a_walk_again(path: &Path, rounds: u64) -> u64 {
    let mut reader = PspReader::open(path).expect("open the store");
    let mut taken = 0u64;
    for _ in 0..rounds {
        let mut walk = reader.records().expect("start the walk");
        if let Some(Ok(record)) = walk.next() {
            taken += 1;
            std::hint::black_box(&record);
        }
    }
    taken
}

/// Open a store, start a walk, take one record and drop it — `rounds` times.
///
/// **What a caller pays to *begin* looking at a sample**, which is a different number from what
/// it pays a record: every `records()` builds a fresh `BlockStream`, and a cohort walking region
/// by region builds one per region per sample.
fn open_and_start_a_walk(path: &Path, rounds: u64) -> u64 {
    let mut taken = 0u64;
    for _ in 0..rounds {
        let mut reader = PspReader::open(path).expect("open the store");
        let mut walk = reader.records().expect("start the walk");
        if let Some(Ok(record)) = walk.next() {
            taken += 1;
            std::hint::black_box(&record);
        }
    }
    taken
}

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
        std::hint::black_box(&record);
    }
    built
}

/// What one phase cost, as a delta over the counters either side of it.
#[cfg(feature = "dhat-heap")]
fn report(phase: &str, records: u64, before: &dhat::HeapStats, after: &dhat::HeapStats) {
    let blocks = after.total_blocks - before.total_blocks;
    let bytes = after.total_bytes - before.total_bytes;
    println!(
        "{phase:<34} records={records:>7}  blocks={blocks:>9}  bytes={bytes:>11}  \
         blocks/record={:>8.3}  bytes/record={:>9.1}",
        blocks as f64 / records as f64,
        bytes as f64 / records as f64,
    );
}

#[cfg(not(feature = "dhat-heap"))]
fn report(phase: &str, records: u64, _before: &(), _after: &()) {
    println!("{phase:<34} records={records:>7}  (build with --features dhat-heap for counts)");
}

#[cfg(feature = "dhat-heap")]
fn stats() -> dhat::HeapStats {
    dhat::HeapStats::get()
}

#[cfg(not(feature = "dhat-heap"))]
fn stats() {}

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::builder().testing().build();

    let scratch = tempfile::tempdir().expect("a scratch directory");

    for shape in [&SHALLOW, &DEEP] {
        let records = a_corpus(shape);
        let path = scratch.path().join(format!("{}.ngpsp", shape.name));

        let before = stats();
        let blocks = write_a_store(&path, a_header(shape), &records);
        let after = stats();
        report(
            &format!("write/{}", shape.name),
            shape.records,
            &before,
            &after,
        );
        assert!(blocks >= 2, "{} cut {blocks} blocks", shape.name);

        // Freed before the walk, so the walk's counters are the reader's alone.
        drop(records);

        let before = stats();
        let built = walk_building_one_in(&path, 1);
        let after = stats();
        assert_eq!(built, shape.records);
        report(
            &format!("walk_full/{}", shape.name),
            shape.records,
            &before,
            &after,
        );

        let before = stats();
        let built = walk_building_one_in(&path, u64::MAX);
        let after = stats();
        assert_eq!(built, 0);
        report(
            &format!("walk_heads/{}", shape.name),
            shape.records,
            &before,
            &after,
        );

        let before = stats();
        let opened = open_and_start_a_walk(&path, 200);
        let after = stats();
        assert_eq!(opened, 200);
        report(
            &format!("open_and_start/{}", shape.name),
            200,
            &before,
            &after,
        );

        let before = stats();
        let opened = open_only(&path, 200);
        let after = stats();
        assert_eq!(opened, 200);
        report(&format!("open_only/{}", shape.name), 200, &before, &after);

        let before = stats();
        let started = start_a_walk_again(&path, 200);
        let after = stats();
        assert_eq!(started, 200);
        report(
            &format!("walk_start_again/{}", shape.name),
            200,
            &before,
            &after,
        );

        let before = stats();
        let built = walk_building_one_in(&path, SECOND_OBSERVATION_ONE_IN);
        let after = stats();
        assert!(built > 0);
        report(
            &format!("walk_one_in_100/{}", shape.name),
            shape.records,
            &before,
            &after,
        );
    }
}
