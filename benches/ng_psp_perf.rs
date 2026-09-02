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
//! one record a position, a live set of about `depth` reads sliding forward as paired 150-base
//! reads start and end, one observation carrying nearly all of them, and a second observation at
//! one record in a hundred. **The reads are paired and their two mates are 200 positions apart**,
//! so an identifier goes absent and comes back — which spec `psp_record_encoding.md` §6 measures
//! at 83 % of identifiers on the human sample and 91 % on tomato, and which is the case the
//! reader's live-set merge cannot take its cheap path on. What they do **not** reproduce is a real
//! sample's variety of witnesses, read groups and locus kinds; for record-for-record fidelity
//! against a real corpus, see `examples/ng_psp_parity.rs`.
//!
//! # Two things these numbers are not
//!
//! **Every walk arm re-opens one store that criterion has already walked, so the file is in the
//! page cache from the second iteration on.** These are decode costs, not end-to-end read costs;
//! a figure from here may not be quoted as what a cohort pays to read a store off a disk.
//!
//! **A `--release` profile of these arms cannot say whether the head or the body holds the time.**
//! `lto = "fat"` with `codegen-units = 1` inlines `read_record_head`, `decode_record_body` and
//! `BlockStream::next_record_where` into `RecordIter::next`, which then shows up as one large
//! entry. Take the breakdown under `--profile profiling` instead (release codegen, no fat LTO,
//! full debug info) and the shares under `--release` for what the shipped binary does — the two
//! builds are 1.48× apart on this walk, so never compare a timing from one with a timing from the
//! other.
//!
//! ⚠ **The heads here are ng-shaped, not production-shaped.** The stores under `tmp/` that the
//! milestone harnesses walk are built from a production `.psp`, which names about 3.4 % of the
//! reads ng will name (the owner's ruling of 2026-08-17), so their heads carry a small fraction of
//! the chain-id changes ng's will. This corpus names one identifier per read, which is what ng
//! will write — so a profile taken here weights the head the way the shipped run will, and a
//! walk timed here is **not** comparable with one timed on a `tmp/` store.

// **The allocator every shipped binary links, and it has to be asked for in each one.**
// `alloc-mimalloc` is a *default* feature, so this is opt-out (`--no-default-features`) and not
// opt-in — but a `#[global_allocator]` is resolved per **binary**, so a bench or example without
// this declaration links the system allocator whatever the feature says. Measured on
// `examples/ng_psp_skip_value.rs`, two builds of one source differing only in this declaration,
// over 7,687,686 tomato records: the full walk took 2.594 s under the system allocator against
// 1.837 s under mimalloc, while the body-skipping walk moved from 0.649 s to 0.628 s — so the
// ratio that harness prints came out 3.99 instead of 2.93.
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
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

/// How many reference positions one mate covers. **150, the Illumina read this project's corpora
/// are.**
const MATE_LENGTH_POSITIONS: u64 = 150;

/// How many positions separate a pair's two mates — the unsequenced middle of the fragment.
///
/// **200, so a fragment spans 500 positions**, which is the ordinary shape of a paired library.
const INNER_GAP_POSITIONS: u64 = 200;

/// A whole fragment: first mate, gap, second mate.
const FRAGMENT_SPAN_POSITIONS: u64 = 2 * MATE_LENGTH_POSITIONS + INNER_GAP_POSITIONS;

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

/// How many records a write workload pushes, for a shape that cuts a block every
/// `genomic_block_size_bp` records.
///
/// **Three and a half blocks' worth, and it has to be a function of the shape.** Blocks are cut on
/// the genomic grid alone, so at one record a reference position a shape cuts a block every
/// `genomic_block_size_bp` records. An earlier version of this file pushed a flat 60,000 records
/// in both arms; against the shallow shape's 100 kb grid that is **one block**, and
/// [`PspWriter::push`] only compresses and writes when a block closes — so the arm named `write`
/// never once ran the writer's compress-and-put path, and did all of its compression inside
/// `finish`. Its number was 15.02 ms against the deep arm's 125.2 ms, and the two were not
/// measuring the same program.
fn records_written(shape: &CorpusShape) -> u64 {
    shape.genomic_block_size_bp * 7 / 2
}

// ---------------------------------------------------------------------
// The corpus
// ---------------------------------------------------------------------

/// The identifiers live at `position`, ascending — **two stretches, because a chain id names a
/// read pair and a pair covers the reference twice with a gap between**.
///
/// **The two stretches are the whole point, and an earlier version of this file got it wrong.** It
/// returned one contiguous range, so every identifier that arrived sorted *above* every identifier
/// already live and the reader's merge could always append. That is the easy case, and it is not
/// the one the data has: spec `psp_record_encoding.md` §6 measures **83 % of identifiers on the
/// human sample and 91 % on tomato covering two stretches**. A benchmark that never makes an
/// identifier come back exercises none of the merging that costs, and reports a saving on the
/// merge that the real shape would not give. Caught by the data-layout review, 2026-08-30.
///
/// The model: identifiers start at a steady rate, each covers `MATE_LENGTH_POSITIONS`, goes absent
/// for `INNER_GAP_POSITIONS`, and covers `MATE_LENGTH_POSITIONS` again. So at any position two
/// bands of identifiers are live — the pairs showing their second mate, and the younger pairs
/// showing their first — and every second-mate arrival sorts **below** every identifier in the
/// younger band. At 280 reads a position that is 0.93 returning identifiers a record, against the
/// 0.8 implied by the spec's 83–91 %; at 10 reads a position it is 0.033.
///
/// A pair covers `2 * MATE_LENGTH_POSITIONS` reference positions in all, so holding `depth` reads
/// live at every position needs `depth / (2 * MATE_LENGTH_POSITIONS)` pairs starting at each.
fn live_at(position: u64, depth: u64, into: &mut Vec<ChainId>) {
    into.clear();
    // `first_starting_after(p)` is the lowest identifier whose pair starts after position `p`,
    // scaled so that `depth / (2 * MATE_LENGTH_POSITIONS)` of them start at every position. The
    // arithmetic is done on the numerator to keep it in integers.
    let first_starting_after = |p: u64| p * depth / (2 * MATE_LENGTH_POSITIONS);
    // The older band: pairs showing their **second** mate here.
    let second_mate_from = first_starting_after(position.saturating_sub(FRAGMENT_SPAN_POSITIONS));
    let second_mate_to =
        first_starting_after(position.saturating_sub(MATE_LENGTH_POSITIONS + INNER_GAP_POSITIONS));
    // The younger band: pairs showing their **first** mate here.
    let first_mate_from = first_starting_after(position.saturating_sub(MATE_LENGTH_POSITIONS));
    let first_mate_to = first_starting_after(position);
    into.extend(second_mate_from..second_mate_to);
    into.extend(first_mate_from..first_mate_to);
}

/// One record: one reference position, one observation holding nearly every live read, and — at
/// one position in a hundred — a second observation holding the rest.
///
/// **`num_obs` equals the identifier count** so that the writer's residual derivation is allowed:
/// it refuses to derive a list whose length does not sit between half the read count and the read
/// count (`record.rs::check_a_read_list_against_its_read_count`), and a record it refuses stores
/// every list instead, which is not the shape a real one has.
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

/// One observation, its support scaled off its read count so that no field is a constant a codec
/// could invent.
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

/// Every record of a corpus, built once.
///
/// **Materialised rather than generated inside the timed region**, because the write workload
/// times `push` and building a record is not part of what a writer costs.
fn a_corpus(shape: &CorpusShape, records: u64) -> Vec<SampleLocusObservations> {
    let depth = u64::from(shape.reads_a_position);
    // One scratch buffer for the whole corpus rather than one allocation a record: the corpus is
    // built outside every timed region, but a 300,000-record build is still worth not making slow.
    let mut live = Vec::new();
    let built: Vec<SampleLocusObservations> = (0..records)
        .map(|at| a_record(at, depth, &mut live))
        .collect();
    check_the_corpus_is_the_shape_it_claims(shape, &built);
    built
}

/// **The corpus's depth and its returning identifiers are load-bearing claims, so they are checked
/// here rather than in a comment.** Both were wrong once: an earlier `live_at` returned one
/// contiguous range, which held the depth right and made every arrival sort above the whole live
/// set, so the reader's merge never once took the branch the real data makes it take.
fn check_the_corpus_is_the_shape_it_claims(shape: &CorpusShape, built: &[SampleLocusObservations]) {
    let depth = u64::from(shape.reads_a_position);
    let mut live_ids = 0u64;
    let mut returning = 0u64;
    // The identifiers live at the record before this one, to see which arrivals sort below them.
    let mut was_live: Vec<ChainId> = Vec::new();
    let mut now_live: Vec<ChainId> = Vec::new();
    // The first records of the corpus are still filling their first fragment span, so the live set
    // has not reached its steady size; they are built but not counted.
    let warm = (FRAGMENT_SPAN_POSITIONS as usize).min(built.len() / 2);
    for (at, record) in built.iter().enumerate() {
        now_live.clear();
        for observation in &record.observations {
            now_live.extend_from_slice(&observation.chain_ids);
        }
        if at >= warm {
            live_ids += now_live.len() as u64;
            // **An arrival the reader's merge cannot simply append is one that sorts below an
            // identifier already live** — that is, below the highest of them, since both lists are
            // ascending. `binary_search` rather than `contains`, because a linear membership test
            // over 280 identifiers a record for 60,000 records is 4.7 billion comparisons at bench
            // setup.
            let highest_live = was_live.last().copied();
            returning += now_live
                .iter()
                .filter(|id| was_live.binary_search(id).is_err())
                .filter(|id| highest_live.is_some_and(|highest| **id < highest))
                .count() as u64;
        }
        std::mem::swap(&mut was_live, &mut now_live);
    }
    let counted = (built.len() - warm) as u64;
    let mean_live = live_ids as f64 / counted as f64;
    assert!(
        (mean_live - depth as f64).abs() < depth as f64 * 0.1,
        "{} holds {mean_live:.1} identifiers a record where it claims {depth}",
        shape.name
    );
    assert!(
        returning > 0,
        "{} has no identifier that goes absent and comes back, so the reader's merge only ever \
         appends — which is the case spec psp_record_encoding.md §6 says 83 % to 91 % of \
         identifiers are not",
        shape.name
    );
    eprintln!(
        "corpus {}: {mean_live:.1} identifiers a record, {:.3} of them a record arriving below \
         one already live",
        shape.name,
        returning as f64 / counted as f64
    );
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

        // **Each arm's name is a claim about how many bodies were built, and each arm checks
        // its own.** A predicate that silently stopped being consulted would leave `walk_heads`
        // and `walk_full` timing the same work with nothing in the output to say so.
        let every_record = shape.records;
        let one_in_a_hundred = shape.records.div_ceil(SECOND_OBSERVATION_ONE_IN);
        group.bench_function(format!("walk_full/{}", shape.name), |b| {
            b.iter(|| {
                let built = walk_full(&path);
                assert_eq!(built, every_record, "a full walk builds every body");
                black_box(built)
            });
        });
        group.bench_function(format!("walk_heads/{}", shape.name), |b| {
            b.iter(|| {
                let built = walk_building_one_in(&path, u64::MAX);
                assert_eq!(built, 0, "a head-only walk builds no body at all");
                black_box(built)
            });
        });
        group.bench_function(format!("walk_one_in_100/{}", shape.name), |b| {
            b.iter(|| {
                let built = walk_building_one_in(&path, SECOND_OBSERVATION_ONE_IN);
                assert_eq!(
                    built, one_in_a_hundred,
                    "one body in a hundred, and no more"
                );
                black_box(built)
            });
        });

        group.finish();
    }
}

fn ng_psp_writes(c: &mut Criterion) {
    let scratch = tempfile::tempdir().expect("a scratch directory");
    for shape in [&SHALLOW, &DEEP] {
        let pushed = records_written(shape);
        let records = a_corpus(shape, pushed);
        let path = scratch.path().join(format!("write_{}.ngpsp", shape.name));

        let mut group = c.benchmark_group("ng_psp");
        group.throughput(Throughput::Elements(pushed));
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(10));
        group.bench_function(format!("write/{}", shape.name), |b| {
            b.iter(|| {
                let blocks = write_a_store(&path, a_header(shape, pushed), &records);
                // **The arm's name is a claim about what ran, so it is checked here.** A write
                // that cuts one block compresses once, inside `finish`, and never reaches the
                // branch of `PspWriter::push` that closes a block — which is what this arm exists
                // to time. See `records_written`.
                assert!(
                    blocks >= 3,
                    "{} wrote {blocks} block(s); a write arm that cuts fewer than three never \
                     pays for a block boundary",
                    shape.name
                );
                black_box(blocks)
            });
        });
        group.finish();
    }
}

criterion_group!(benches, ng_psp_walks, ng_psp_writes);
criterion_main!(benches);
