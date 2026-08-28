//! **Does `src/ng/psp/` give back the records it was handed, and does an independent encoder
//! of the same shape agree with it?**
//!
//! Milestone H1's oracle. Every byte figure and every claim of fidelity the psp specs carry was
//! taken on the measuring prototype (`examples/psp_row_stream_roundtrip.rs`) over production's
//! `PileupRecord`. This harness is where the module that shipped is checked instead — against a
//! real corpus rather than a forty-record fixture, and against a codec that is not it.
//!
//! ```text
//! cargo run --release --example ng_psp_parity -- <a production .psp> [--limit N] [--work DIR]
//! ```
//!
//! # The three streams, and what each pair proves
//!
//! One production `.psp` is the source. Every record it holds is turned into an ng record; both
//! stores are written from that same sequence, and all three are then walked in lockstep.
//!
//! | pair | what a disagreement means |
//! |---|---|
//! | the ng record pushed → the ng record read back | `src/ng/psp/` lost or changed a field |
//! | the ng record read back → the prototype's record read back | the two encoders disagree, and only one of them is the code under test |
//! | the prototype's record → the source record | the prototype moved a number (its own quantisation) |
//!
//! **The second pair is the one this step exists for.** A round-trip through one codec proves
//! self-consistency: an encoder that writes a field wrongly and a decoder that reads it back the
//! same wrong way agree with each other perfectly. The prototype was written before this module
//! and shares no code with it.
//!
//! # The strictness, field by field
//!
//! **A blanket tolerance would pass while a chain-id list was being corrupted**, which is why
//! there is no blanket tolerance here:
//!
//! - **Exactly**, on both pairs: the coordinate, the reference bases, every allele sequence,
//!   every read witness, every read group, every count (`num_obs`, `num_fwd`, `mapq_sum`,
//!   `mapq_sum_sq`, `placed_left`, the two counts of reads that showed nothing), the locus kind,
//!   and **every chain-id list**.
//! - **The summed log-error inside its own step.** It is an integer count of
//!   1/[`SummedLogError::STEPS_PER_NAT`] of a natural log since Milestone B3, so it round-trips
//!   through ng's store *exactly*; what it is not is exactly the `f64` production stored, and
//!   the distance to that is reported against half a step, which is the most rounding to nearest
//!   can cost.
//! - **And the prototype is asked for the same step**, rather than for its default. Its scale is
//!   a number in its file, so this harness sets it to ng's — and then the two stores must agree
//!   on the summed log-error *exactly* as well, which is a far stronger check than a tolerance
//!   wide enough for both.
//!
//! **Chain-id lists are compared as sets**, because that is what the format stores: ng's encoder
//! writes each list as ascending gaps, so a list arrives back sorted and deduplicated whatever
//! order it went in. The harness normalises the source's list the same way and **reports how
//! many lists the normalisation actually changed**, so the reader can see whether that
//! concession was doing any work.
//!
//! # What the corpus is, and what this harness adds to it
//!
//! ng cannot yet produce a psp of its own, so the records come from a production `.psp` — and
//! production's record is missing three things ng's carries (arch §2). Two of them,
//! **`read_witness` and `read_group`, are synthesised here rather than left at their empty
//! values**, and that is deliberate: a corpus where every witness is `Complete` and every read
//! group is 0 cannot tell an encoder that stores those fields from one that drops them. They are
//! derived from each record's own coordinate, so the corpus is a function of the input alone and
//! two runs give the same file.
//!
//! **Every precondition the comparison rests on is asserted before anything is proven** — that
//! the corpus contains records with chain ids, with more than one observation, with a witness
//! that is not `Complete`, with more than one read group, and enough records to cut more than
//! one block. A fixture that cannot separate the thing a test names is the failure this milestone
//! keeps finding.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation, WitnessedLocusPositions,
};
use pop_var_caller::ng::psp::{
    ContigIdentity, DEFAULT_LOOK_BACK_WINDOW_LOG, FORMAT_VERSION, Header, Manifest, ParameterValue,
    PspReader, PspWriter, ReferenceIdentity, WriterProvenance, ZSTD_COMPRESSION_LEVEL,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};
use pop_var_caller::pileup_record::{ChainId, PileupRecord};
use pop_var_caller::psp::PspReader as ProductionPspReader;

// **The prototype, included whole and unchanged.** An example cannot depend on another example,
// and copying its encoder in would make the oracle a fork of itself — the one property it has
// that this module's own tests do not is that it was written by someone who was not writing this
// module. `#[path]` keeps one copy.
//
// **The two `allow`s are the prototype's own and are deliberately not fixed there.** Including it
// puts it in this target's lint selection for the first time — examples are not in the gate
// otherwise — and an oracle edited to please a lint in the code it is checking is a little less
// of an oracle. Only its visibility was widened.
#[allow(dead_code, unused_assignments, clippy::too_many_arguments)]
#[path = "psp_row_stream_roundtrip.rs"]
mod prototype;

/// How many bytes of a block the prototype's cut is allowed before it closes one early. Large
/// enough that the genomic grid is what cuts every block, which is what ng's writer does.
const PROTOTYPE_BLOCK_BYTES: usize = 1 << 20;

/// The genomic grid both stores cut on, in base pairs.
///
/// **Smaller than the settled 100 kb**, so that a corpus of a few hundred thousand records cuts
/// hundreds of blocks rather than a handful: every running difference in the format resets at a
/// block boundary, and a parity run that crossed two of them would be proving very little about
/// the reset.
const GRID_BP: u32 = 10_000;

// ---------------------------------------------------------------------
// The corpus: a production record as an ng record
// ---------------------------------------------------------------------

/// One production record as an ng record, or `None` for a record the ng encoder has no span for.
///
/// **Not the same mapping as `ng_psp_head_encoding.rs`'s**, and the difference is the point.
/// That harness compares two encodings of the same four head fields, so it wants the records to
/// carry as little as possible beside them and leaves the witness `Complete` and the read group
/// 0 throughout. Here those two fields are under test, and a corpus that holds one value of each
/// cannot fail if they are dropped.
fn as_an_ng_record(record: &PileupRecord) -> Option<SampleLocusObservations> {
    let reference_bases: Box<[u8]> = record.alleles.first()?.seq.clone().into_boxed_slice();
    let span = u64::try_from(reference_bases.len()).ok()?;
    if span == 0 {
        return None;
    }
    let start = u64::from(record.pos);
    let observations = record
        .alleles
        .iter()
        .enumerate()
        .map(|(at, allele)| SequenceObservation {
            bases: allele.seq.clone().into_boxed_slice(),
            read_witness: a_witness_for(start, at, span),
            read_group: ReadGroupId(((start + at as u64) % 4) as u32),
            num_obs: allele.support.num_obs,
            num_fwd: allele.support.fwd,
            q_sum: SummedLogError::from_nats(allele.support.q_sum),
            mapq_sum: allele.support.mapq_sum,
            mapq_sum_sq: allele.support.mapq_sum_sq,
            placed_left: allele.support.placed_left,
            chain_ids: allele.chain_ids.clone(),
        })
        .collect();
    Some(SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(record.chrom_id),
            start: Position(start),
            end: Position(start + span - 1),
        },
        reference_bases,
        observations,
        // **Synthesised, like the witness and the read group, and for the same reason.**
        // Production has no equivalent of either count (arch §2), so left at zero they make the
        // corpus unable to fail: an encoder writing a constant 0 for both passes every
        // comparison. Measured — with both at zero, a mutant that wrote 0 for
        // `reads_discarded_by_cap` survived a 3,000-record parity run.
        reads_without_observation: (start % 7) as u32,
        reads_discarded_by_cap: (start % 3) as u32,
        kind: LocusKind::Generic,
    })
}

/// A witness for one observation, chosen from the record's own coordinate so that the corpus is
/// a function of its input.
///
/// Three shapes, and each is a different thing for the codec to get wrong: `Complete`, which is
/// a run count of zero on the wire; a `Partial` that covers the whole locus, which must **not**
/// come back as `Complete`; and a `Partial` with a hole in it, which is the shape a spliced read
/// across a widened record has and the one a two-number witness could not express.
fn a_witness_for(start: u64, at: usize, span: u64) -> ReadWitness {
    let locus_len = u16::try_from(span.min(u64::from(u16::MAX))).expect("clamped to u16");
    if locus_len == 0 {
        return ReadWitness::Complete;
    }
    let runs: Vec<(u16, u16)> = match (start + at as u64) % 3 {
        0 => return ReadWitness::Complete,
        1 => vec![(0, locus_len)],
        _ if locus_len >= 3 => vec![(0, 1), (locus_len - 1, locus_len)],
        _ => vec![(0, 1)],
    };
    match WitnessedLocusPositions::from_half_open_runs(runs) {
        Some(positions) => ReadWitness::Partial { positions },
        None => ReadWitness::Complete,
    }
}

/// A list of chain ids as the format stores it: ascending and without duplicates.
fn as_a_read_set(ids: &[ChainId]) -> Vec<ChainId> {
    let mut set = ids.to_vec();
    set.sort_unstable();
    set.dedup();
    set
}

// ---------------------------------------------------------------------
// Writing the two stores
// ---------------------------------------------------------------------

/// ng's header, built from the production file's own so that the two stores describe the same
/// reference and the same contig numbering.
fn an_ng_header(source: &Path, production: &ProductionPspReader<BufReader<File>>) -> Header {
    let parsed = production.header();
    Header {
        format_version: FORMAT_VERSION,
        sample: parsed.sample.clone(),
        reference: ReferenceIdentity {
            name: parsed.reference.clone(),
            md5: None,
        },
        contigs: parsed
            .chromosomes
            .iter()
            .map(|chromosome| ContigIdentity {
                name: chromosome.name.clone(),
                length: u64::from(chromosome.length),
                md5: as_an_md5(&chromosome.md5),
            })
            .collect(),
        writer: WriterProvenance {
            tool: "ng_psp_parity".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "parity".to_string(),
            // A basename, because the header refuses a path with a directory component in it —
            // a psp is meant to be movable and a recorded absolute path is not.
            input_alignments: vec![
                source
                    .file_name()
                    .map(|name| name.to_string_lossy().into_owned())
                    .unwrap_or_else(|| "unknown.psp".to_string()),
            ],
            input_reference: basename_of(&parsed.reference),
            command_line: std::env::args().collect::<Vec<_>>().join(" "),
            parameters: BTreeMap::from([(
                "genomic-block-size-bp".to_string(),
                ParameterValue::Integer(i64::from(GRID_BP)),
            )]),
            created: "2026-08-28T00:00:00Z"
                .parse()
                .expect("a valid RFC 3339 stamp"),
        },
        manifest: Manifest {
            genomic_block_size_bp: pop_var_caller::ng::types::Bp(u64::from(GRID_BP)),
            ..Manifest::as_this_build_writes_it()
        },
    }
}

/// The last path component of `name`, since the header records a basename and not a path.
fn basename_of(name: &str) -> String {
    Path::new(name)
        .file_name()
        .map(|it| it.to_string_lossy().into_owned())
        .unwrap_or_else(|| name.to_string())
}

/// 32 hex characters as sixteen bytes, or `None` for anything else — an older production file
/// leaves the field empty.
fn as_an_md5(hex: &str) -> Option<[u8; 16]> {
    if hex.len() != 32 {
        return None;
    }
    let mut bytes = [0u8; 16];
    for (at, byte) in bytes.iter_mut().enumerate() {
        *byte = u8::from_str_radix(hex.get(at * 2..at * 2 + 2)?, 16).ok()?;
    }
    Some(bytes)
}

/// Write the ng store, and say how many source records had no ng record and how many blocks it
/// cut.
fn write_the_ng_store(source: &Path, out: &Path, limit: usize) -> (u64, u64, u64) {
    let mut production = open_production(source);
    let header = an_ng_header(source, &production);
    let mut writer = PspWriter::create(out, header).expect("create the ng store");
    let (mut pushed, mut skipped) = (0u64, 0u64);
    {
        let records = production.records();
        for record in records {
            let record = record.expect("a production record");
            match as_an_ng_record(&record) {
                Some(record) => {
                    writer.push(&record).expect("push");
                    pushed += 1;
                }
                None => skipped += 1,
            }
            if (pushed + skipped) as usize >= limit {
                break;
            }
        }
    }
    let stats = writer
        .finish(b"ng_psp_parity")
        .expect("finish the ng store");
    (pushed, skipped, stats.blocks)
}

fn open_production(path: &Path) -> ProductionPspReader<BufReader<File>> {
    ProductionPspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(path).expect("open the production psp"),
    ))
    .expect("read the production psp")
}

// ---------------------------------------------------------------------
// The lockstep comparison
// ---------------------------------------------------------------------

/// What the corpus turned out to contain. **Asserted before any parity claim is made**: each of
/// these is a shape the comparison would otherwise be unable to fail on.
#[derive(Default)]
struct CorpusShape {
    records: u64,
    observations: u64,
    with_chain_ids: u64,
    chain_ids: u64,
    lists_reordered_by_normalising: u64,
    with_several_observations: u64,
    complete_witnesses: u64,
    partial_witnesses: u64,
    witnesses_with_a_hole: u64,
    read_groups: [u64; 4],
    records_with_reads_that_showed_nothing: u64,
    records_with_reads_the_cap_discarded: u64,
    blocks: u64,
}

/// The worst distance seen between a stored summed log-error and the `f64` it came from.
#[derive(Default)]
struct WorstDrift {
    ng_against_the_source: f64,
    prototype_against_the_source: f64,
}

fn main() {
    let mut args = std::env::args().skip(1);
    let source = PathBuf::from(
        args.next()
            .expect("usage: ng_psp_parity <a production .psp> [--limit N] [--work DIR]"),
    );
    let mut limit = usize::MAX;
    let mut work = PathBuf::from("tmp/ng_psp_parity");
    let mut label = source
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .unwrap_or_else(|| "corpus".to_string());
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--limit" => limit = args.next().expect("N").parse().expect("N"),
            "--work" => work = PathBuf::from(args.next().expect("DIR")),
            "--label" => label = args.next().expect("L"),
            other => panic!("unknown flag {other}"),
        }
    }
    std::fs::create_dir_all(&work).expect("create the work directory");
    let ng_store = work.join(format!("{label}.ngpsp"));
    let prototype_store = work.join(format!("{label}.ngs"));

    let (pushed, skipped, blocks) = write_the_ng_store(&source, &ng_store, limit);
    assert!(pushed > 0, "the corpus produced no ng records at all");

    // **The prototype is asked for ng's own step**, so the two stores must then agree on the
    // summed log-error exactly rather than inside a tolerance wide enough to hold both.
    let scales = prototype::Scales {
        q_sum: SummedLogError::STEPS_PER_NAT as f64,
        ..prototype::Scales::default()
    };
    prototype::encode_streaming(
        source.to_str().expect("a utf-8 path"),
        prototype_store.to_str().expect("a utf-8 path"),
        PROTOTYPE_BLOCK_BYTES,
        GRID_BP,
        false,
        false,
        true,
        ZSTD_COMPRESSION_LEVEL,
        u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG),
        scales,
    );

    let (shape, worst) = compare(&source, &ng_store, &prototype_store, limit, blocks);
    report(&label, &shape, &worst, pushed, skipped);
}

/// Walk the source, the ng store and the prototype's store together, and fail on the first
/// record that disagrees.
fn compare(
    source: &Path,
    ng_store: &Path,
    prototype_store: &Path,
    limit: usize,
    blocks: u64,
) -> (CorpusShape, WorstDrift) {
    let mut production = open_production(source);
    let mut ng = PspReader::open(ng_store).expect("open the ng store");
    let mut proto = prototype::StreamingStore::open(prototype_store);
    let mut shape = CorpusShape {
        blocks,
        ..CorpusShape::default()
    };
    let mut worst = WorstDrift::default();

    let mut ng_records = ng.records().expect("walk the ng store");
    let source_records = production.records();
    for record in source_records.take(limit) {
        let from_production = record.expect("a production record");
        // The prototype writes every source record, so its stream advances even for the ones
        // that have no ng record — otherwise the two walks would slide apart by the count of
        // skipped records and every later comparison would be against the wrong record.
        let from_prototype = proto.next().expect("the prototype's store ended early");
        let Some(want) = as_an_ng_record(&from_production) else {
            continue;
        };
        let got = ng_records
            .next()
            .expect("the ng store ended early")
            .expect("an ng record")
            .record
            .expect("a walk that declines no body builds every one");
        let at = shape.records;
        check_the_ng_round_trip(at, &want, &got, &mut shape);
        check_against_the_prototype(at, &got, &from_prototype);
        note_the_drift(&from_production, &got, &from_prototype, &mut worst);
        shape.records += 1;
    }
    assert!(
        ng_records.next().is_none(),
        "the ng store holds records the source does not"
    );
    (shape, worst)
}

/// Every field of the record that was pushed, against the record that came back. **Exactly**,
/// field by field — the summed log-error included, since it is an integer count of steps by the
/// time the store sees it.
fn check_the_ng_round_trip(
    at: u64,
    want: &SampleLocusObservations,
    got: &SampleLocusObservations,
    shape: &mut CorpusShape,
) {
    assert_eq!(got.region, want.region, "record {at}: region");
    assert_eq!(
        got.reference_bases, want.reference_bases,
        "record {at}: reference bases"
    );
    assert_eq!(got.kind, want.kind, "record {at}: locus kind");
    assert_eq!(
        got.reads_without_observation, want.reads_without_observation,
        "record {at}: reads without observation"
    );
    assert_eq!(
        got.reads_discarded_by_cap, want.reads_discarded_by_cap,
        "record {at}: reads discarded by the cap"
    );
    if want.reads_without_observation > 0 {
        shape.records_with_reads_that_showed_nothing += 1;
    }
    if want.reads_discarded_by_cap > 0 {
        shape.records_with_reads_the_cap_discarded += 1;
    }
    assert_eq!(
        got.observations.len(),
        want.observations.len(),
        "record {at}: observation count"
    );
    if want.observations.len() > 1 {
        shape.with_several_observations += 1;
    }
    let mut any_chain_ids = false;
    for (which, (g, w)) in got.observations.iter().zip(&want.observations).enumerate() {
        assert_eq!(g.bases, w.bases, "record {at}, observation {which}: bases");
        assert_eq!(
            g.read_witness, w.read_witness,
            "record {at}, observation {which}: read witness"
        );
        assert_eq!(
            g.read_group, w.read_group,
            "record {at}, observation {which}: read group"
        );
        assert_eq!(
            g.num_obs, w.num_obs,
            "record {at}, observation {which}: num_obs"
        );
        assert_eq!(
            g.num_fwd, w.num_fwd,
            "record {at}, observation {which}: num_fwd"
        );
        assert_eq!(
            g.q_sum.steps(),
            w.q_sum.steps(),
            "record {at}, observation {which}: summed log-error"
        );
        assert_eq!(
            g.mapq_sum, w.mapq_sum,
            "record {at}, observation {which}: mapq_sum"
        );
        assert_eq!(
            g.mapq_sum_sq, w.mapq_sum_sq,
            "record {at}, observation {which}: mapq_sum_sq"
        );
        assert_eq!(
            g.placed_left, w.placed_left,
            "record {at}, observation {which}: placed_left"
        );
        let wanted = as_a_read_set(&w.chain_ids);
        if wanted != w.chain_ids {
            shape.lists_reordered_by_normalising += 1;
        }
        assert_eq!(
            g.chain_ids, wanted,
            "record {at}, observation {which}: chain ids"
        );
        shape.observations += 1;
        shape.chain_ids += wanted.len() as u64;
        any_chain_ids |= !wanted.is_empty();
        match &w.read_witness {
            ReadWitness::Complete => shape.complete_witnesses += 1,
            ReadWitness::Partial { positions } => {
                shape.partial_witnesses += 1;
                if positions.runs().len() > 1 {
                    shape.witnesses_with_a_hole += 1;
                }
            }
        }
        let group = w.read_group.get() as usize;
        if group < shape.read_groups.len() {
            shape.read_groups[group] += 1;
        }
    }
    if any_chain_ids {
        shape.with_chain_ids += 1;
    }
}

/// The ng store's record against the prototype's, on every field the prototype's record has.
///
/// **The three ng carries that production does not** — the witness, the read group and the two
/// counts of reads that showed nothing — have no counterpart here, and are held by
/// [`check_the_ng_round_trip`] instead.
fn check_against_the_prototype(at: u64, got: &SampleLocusObservations, proto: &PileupRecord) {
    assert_eq!(
        got.region.contig.0, proto.chrom_id,
        "record {at}: contig, against the prototype"
    );
    assert_eq!(
        got.region.start.0,
        u64::from(proto.pos),
        "record {at}: position, against the prototype"
    );
    assert_eq!(
        got.observations.len(),
        proto.alleles.len(),
        "record {at}: allele count, against the prototype"
    );
    for (which, (g, p)) in got.observations.iter().zip(&proto.alleles).enumerate() {
        assert_eq!(
            g.bases.as_ref(),
            p.seq.as_slice(),
            "record {at}, allele {which}: sequence, against the prototype"
        );
        assert_eq!(
            g.chain_ids,
            as_a_read_set(&p.chain_ids),
            "record {at}, allele {which}: chain ids, against the prototype"
        );
        assert_eq!(
            g.num_obs, p.support.num_obs,
            "record {at}, allele {which}: num_obs, against the prototype"
        );
        assert_eq!(
            g.num_fwd, p.support.fwd,
            "record {at}, allele {which}: num_fwd, against the prototype"
        );
        assert_eq!(
            g.mapq_sum, p.support.mapq_sum,
            "record {at}, allele {which}: mapq_sum, against the prototype"
        );
        assert_eq!(
            g.mapq_sum_sq, p.support.mapq_sum_sq,
            "record {at}, allele {which}: mapq_sum_sq, against the prototype"
        );
        // Exact, because the prototype was handed ng's own step: both stores round the same
        // `f64` to the same integer count of 1/4,096 of a natural log, so a difference here is
        // one of them storing a different number, not one of them rounding differently.
        assert_eq!(
            g.q_sum.steps(),
            (p.support.q_sum * SummedLogError::STEPS_PER_NAT as f64).round() as i64,
            "record {at}, allele {which}: summed log-error, against the prototype"
        );
    }
}

/// How far each store's summed log-error ended up from the `f64` production held. Reported, not
/// asserted per record — the assertion is against half a step, once, at the end.
fn note_the_drift(
    source: &PileupRecord,
    ng: &SampleLocusObservations,
    proto: &PileupRecord,
    worst: &mut WorstDrift,
) {
    for (which, allele) in source.alleles.iter().enumerate() {
        if let Some(observation) = ng.observations.get(which) {
            let drift = (observation.q_sum.nats() - allele.support.q_sum).abs();
            worst.ng_against_the_source = worst.ng_against_the_source.max(drift);
        }
        if let Some(other) = proto.alleles.get(which) {
            let drift = (other.support.q_sum - allele.support.q_sum).abs();
            worst.prototype_against_the_source = worst.prototype_against_the_source.max(drift);
        }
    }
}

// ---------------------------------------------------------------------
// The report, and the preconditions it rests on
// ---------------------------------------------------------------------

fn report(label: &str, shape: &CorpusShape, worst: &WorstDrift, pushed: u64, skipped: u64) {
    let half_a_step = 0.5 / SummedLogError::STEPS_PER_NAT as f64;
    assert_the_corpus_can_fail(shape);
    assert!(
        worst.ng_against_the_source <= half_a_step,
        "the ng store's summed log-error is {} from the source, and rounding to the nearest \
         1/{} of a natural log cannot cost more than {half_a_step}",
        worst.ng_against_the_source,
        SummedLogError::STEPS_PER_NAT
    );

    println!("phase\tng-psp-parity");
    println!("corpus\t{label}");
    println!("records-compared\t{}", shape.records);
    println!("records-pushed\t{pushed}");
    println!("records-with-no-ng-record\t{skipped}");
    println!("blocks\t{}", shape.blocks);
    println!("observations\t{}", shape.observations);
    println!(
        "records-with-several-observations\t{}",
        shape.with_several_observations
    );
    println!("records-naming-reads\t{}", shape.with_chain_ids);
    println!("chain-ids-compared\t{}", shape.chain_ids);
    println!(
        "chain-id-lists-the-normalising-changed\t{}",
        shape.lists_reordered_by_normalising
    );
    println!("witnesses-complete\t{}", shape.complete_witnesses);
    println!("witnesses-partial\t{}", shape.partial_witnesses);
    println!("witnesses-with-a-hole\t{}", shape.witnesses_with_a_hole);
    for (group, count) in shape.read_groups.iter().enumerate() {
        println!("observations-in-read-group-{group}\t{count}");
    }
    println!(
        "records-with-reads-that-showed-nothing\t{}",
        shape.records_with_reads_that_showed_nothing
    );
    println!(
        "records-with-reads-the-cap-discarded\t{}",
        shape.records_with_reads_the_cap_discarded
    );
    println!(
        "worst-summed-log-error-drift-ng\t{:.9}\thalf-a-step\t{half_a_step:.9}",
        worst.ng_against_the_source
    );
    println!(
        "worst-summed-log-error-drift-prototype\t{:.9}",
        worst.prototype_against_the_source
    );
}

/// **Every shape the comparison would otherwise be unable to fail on.** A corpus where no record
/// names a read proves nothing about the chain-id column; one where every witness is `Complete`
/// proves nothing about the witness; one that cuts a single block proves nothing about the
/// running differences that reset at a boundary.
fn assert_the_corpus_can_fail(shape: &CorpusShape) {
    assert!(shape.records > 0, "no records were compared");
    assert!(
        shape.blocks > 1,
        "the corpus cut {} block; the running differences reset at a boundary and a \
         single-block corpus never crosses one",
        shape.blocks
    );
    assert!(
        shape.with_chain_ids > 0,
        "no record in the corpus names a read, so the chain-id comparison cannot fail"
    );
    assert!(
        shape.with_several_observations > 0,
        "every record in the corpus has one observation, so the residual derivation is never \
         exercised"
    );
    assert!(
        shape.complete_witnesses > 0 && shape.partial_witnesses > 0,
        "the corpus holds {} complete and {} partial witnesses; it needs both, or an encoder \
         that writes one value for every witness passes",
        shape.complete_witnesses,
        shape.partial_witnesses
    );
    assert!(
        shape.witnesses_with_a_hole > 0,
        "no witness in the corpus has more than one run, so a codec that keeps only the \
         outermost edges passes"
    );
    assert!(
        shape.read_groups.iter().filter(|count| **count > 0).count() > 1,
        "every observation in the corpus is in one read group, so a dropped read group passes"
    );
    assert!(
        shape.records_with_reads_that_showed_nothing > 0
            && shape.records_with_reads_the_cap_discarded > 0,
        "the corpus has {} records with reads that showed nothing and {} with reads the cap \
         discarded; with either at zero an encoder writing a constant 0 for that count passes",
        shape.records_with_reads_that_showed_nothing,
        shape.records_with_reads_the_cap_discarded
    );
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

/// **What these hold is that the comparison can fail**, which is the one thing a parity harness
/// cannot demonstrate by passing. Each takes a record that agrees, changes one field, and
/// requires the check to refuse it — so a comparison that quietly stopped looking at a field
/// fails here rather than reporting a clean run over eight million records.
///
/// They also hold the two concessions the comparison makes — chain-id lists compared as sets,
/// and the summed log-error compared against the source inside its own step — since a concession
/// nothing tests is indistinguishable from no check at all.
#[cfg(test)]
mod tests {
    use super::*;
    use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats};

    /// A record that names reads, holds two observations, and carries a witness with a hole in
    /// it — the three shapes the comparison would otherwise be unable to fail on.
    fn a_record() -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(1),
                start: Position(100),
                end: Position(103),
            },
            reference_bases: b"ACGT".to_vec().into_boxed_slice(),
            observations: vec![
                SequenceObservation {
                    bases: b"ACGT".to_vec().into_boxed_slice(),
                    read_witness: ReadWitness::Complete,
                    read_group: ReadGroupId(0),
                    num_obs: 4,
                    num_fwd: 2,
                    q_sum: SummedLogError::from_steps(-8_192),
                    mapq_sum: 240,
                    mapq_sum_sq: 14_400,
                    placed_left: 1,
                    chain_ids: vec![7, 9, 11],
                },
                SequenceObservation {
                    bases: b"ACCT".to_vec().into_boxed_slice(),
                    read_witness: ReadWitness::Partial {
                        positions: WitnessedLocusPositions::from_half_open_runs([(0, 1), (3, 4)])
                            .expect("two runs with a gap between them"),
                    },
                    read_group: ReadGroupId(1),
                    num_obs: 2,
                    num_fwd: 1,
                    q_sum: SummedLogError::from_steps(-4_096),
                    mapq_sum: 120,
                    mapq_sum_sq: 7_200,
                    placed_left: 0,
                    chain_ids: vec![13, 17],
                },
            ],
            reads_without_observation: 3,
            reads_discarded_by_cap: 1,
            kind: LocusKind::Generic,
        }
    }

    /// The same record as production would have held it: no witness, no read group, and the
    /// summed log-error still a float.
    fn the_same_record_as_production_holds_it() -> PileupRecord {
        let record = a_record();
        PileupRecord::new(
            record.region.contig.0,
            record.region.start.0 as u32,
            record
                .observations
                .iter()
                .map(|observation| {
                    AlleleObservation::new(
                        observation.bases.to_vec(),
                        AlleleSupportStats::new(
                            observation.num_obs,
                            observation.q_sum.nats(),
                            observation.num_fwd,
                            observation.placed_left,
                            0,
                            observation.mapq_sum,
                            observation.mapq_sum_sq,
                        ),
                        observation.chain_ids.clone(),
                    )
                })
                .collect(),
        )
    }

    /// A corpus whose every precondition is met, so that a test which breaks one of them is
    /// testing that one.
    fn a_corpus_that_can_fail() -> CorpusShape {
        CorpusShape {
            records: 10,
            observations: 20,
            with_chain_ids: 5,
            chain_ids: 30,
            lists_reordered_by_normalising: 0,
            with_several_observations: 5,
            complete_witnesses: 10,
            partial_witnesses: 10,
            witnesses_with_a_hole: 2,
            read_groups: [5, 5, 5, 5],
            records_with_reads_that_showed_nothing: 8,
            records_with_reads_the_cap_discarded: 6,
            blocks: 3,
        }
    }

    fn round_trip(want: &SampleLocusObservations, got: &SampleLocusObservations) {
        check_the_ng_round_trip(0, want, got, &mut CorpusShape::default());
    }

    #[test]
    fn a_record_that_came_back_unchanged_is_accepted_on_both_pairs() {
        let record = a_record();
        round_trip(&record, &record);
        check_against_the_prototype(0, &record, &the_same_record_as_production_holds_it());
    }

    #[test]
    #[should_panic(expected = "observation 1: read witness")]
    fn a_witness_that_came_back_complete_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.observations[1].read_witness = ReadWitness::Complete;
        round_trip(&want, &got);
    }

    #[test]
    #[should_panic(expected = "observation 1: read witness")]
    fn a_witness_that_lost_its_hole_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.observations[1].read_witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 4)])
                .expect("one run over the whole locus"),
        };
        round_trip(&want, &got);
    }

    #[test]
    #[should_panic(expected = "observation 1: read group")]
    fn a_read_group_that_came_back_as_zero_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.observations[1].read_group = ReadGroupId(0);
        round_trip(&want, &got);
    }

    #[test]
    #[should_panic(expected = "observation 0: chain ids")]
    fn a_chain_id_list_missing_one_read_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.observations[0].chain_ids.pop();
        round_trip(&want, &got);
    }

    #[test]
    #[should_panic(expected = "observation 0: summed log-error")]
    fn a_summed_log_error_one_step_out_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.observations[0].q_sum = SummedLogError::from_steps(-8_191);
        round_trip(&want, &got);
    }

    #[test]
    #[should_panic(expected = "reads without observation")]
    fn a_lost_count_of_reads_that_showed_nothing_is_refused() {
        let want = a_record();
        let mut got = want.clone();
        got.reads_without_observation = 0;
        round_trip(&want, &got);
    }

    /// **The one concession the round-trip check makes.** ng's encoder writes a chain-id list as
    /// ascending gaps, so a list arrives back sorted and deduplicated whatever order it went in;
    /// the check normalises the pushed list the same way rather than demanding an order the
    /// format does not keep. What it must still refuse is a list with a *different* set in it,
    /// which the test above holds.
    #[test]
    fn a_chain_id_list_that_came_back_sorted_is_accepted_and_counted() {
        let mut want = a_record();
        want.observations[0].chain_ids = vec![11, 7, 9, 7];
        let mut got = want.clone();
        got.observations[0].chain_ids = vec![7, 9, 11];
        let mut shape = CorpusShape::default();
        check_the_ng_round_trip(0, &want, &got, &mut shape);
        assert_eq!(
            shape.lists_reordered_by_normalising, 1,
            "the harness has to report that the normalising did work, or a reader cannot tell \
             this concession from an exact comparison"
        );
    }

    #[test]
    #[should_panic(expected = "allele 1: summed log-error, against the prototype")]
    fn a_summed_log_error_the_prototype_disagrees_with_is_refused() {
        let got = a_record();
        let mut proto = the_same_record_as_production_holds_it();
        proto.alleles[1].support.q_sum += 1.0 / SummedLogError::STEPS_PER_NAT as f64;
        check_against_the_prototype(0, &got, &proto);
    }

    #[test]
    #[should_panic(expected = "allele 0: sequence, against the prototype")]
    fn an_allele_sequence_the_prototype_disagrees_with_is_refused() {
        let got = a_record();
        let mut proto = the_same_record_as_production_holds_it();
        proto.alleles[0].seq = b"ACGA".to_vec();
        check_against_the_prototype(0, &got, &proto);
    }

    #[test]
    #[should_panic(expected = "allele 0: chain ids, against the prototype")]
    fn a_chain_id_list_the_prototype_disagrees_with_is_refused() {
        let got = a_record();
        let mut proto = the_same_record_as_production_holds_it();
        proto.alleles[0].chain_ids.push(19);
        check_against_the_prototype(0, &got, &proto);
    }

    /// **A `Partial` that covers every position of its locus is not `Complete`**, and the
    /// difference is what the reader's likelihood gates on: `Complete` means the read reached
    /// both borders. The corpus mints one of these deliberately, so a codec that folded the two
    /// together would be caught.
    #[test]
    fn a_witness_covering_the_whole_locus_stays_partial() {
        let witness = a_witness_for(1, 0, 4);
        match witness {
            ReadWitness::Partial { positions } => {
                assert_eq!(
                    positions.runs().collect::<Vec<_>>(),
                    vec![(0, 4)],
                    "one run flush with both borders"
                );
            }
            ReadWitness::Complete => panic!("the whole-locus witness collapsed to Complete"),
        }
    }

    #[test]
    fn the_three_witness_shapes_all_appear_and_a_hole_needs_three_positions() {
        assert_eq!(a_witness_for(0, 0, 4), ReadWitness::Complete);
        let with_a_hole = a_witness_for(2, 0, 4);
        let ReadWitness::Partial { positions } = &with_a_hole else {
            panic!("the third shape is a partial witness");
        };
        assert_eq!(
            positions.runs().collect::<Vec<_>>(),
            vec![(0, 1), (3, 4)],
            "a run at each border with a gap between them"
        );
        // Below three positions there is no room for a gap, so the same arm gives one run —
        // which is why the corpus's count of witnesses with a hole is small on a corpus whose
        // records are mostly one position long.
        let ReadWitness::Partial { positions } = a_witness_for(2, 0, 2) else {
            panic!("still a partial witness");
        };
        assert_eq!(positions.runs().len(), 1);
    }

    #[test]
    fn a_corpus_meeting_every_precondition_is_accepted() {
        assert_the_corpus_can_fail(&a_corpus_that_can_fail());
    }

    #[test]
    #[should_panic(expected = "never crosses one")]
    fn a_corpus_of_one_block_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            blocks: 1,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "names a read")]
    fn a_corpus_naming_no_reads_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            with_chain_ids: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "residual derivation is never exercised")]
    fn a_corpus_of_single_observation_records_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            with_several_observations: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "it needs both")]
    fn a_corpus_whose_witnesses_are_all_complete_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            partial_witnesses: 0,
            witnesses_with_a_hole: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "outermost edges")]
    fn a_corpus_with_no_hole_in_any_witness_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            witnesses_with_a_hole: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "one read group")]
    fn a_corpus_of_one_read_group_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            read_groups: [20, 0, 0, 0],
            ..a_corpus_that_can_fail()
        });
    }

    /// **This one is not hypothetical.** With both counts left at production's zero, a mutant
    /// that wrote a constant 0 for `reads_discarded_by_cap` survived a 3,000-record parity run
    /// over hg002 — the only one of six injected defects that did.
    #[test]
    #[should_panic(expected = "reads the cap")]
    fn a_corpus_where_no_read_was_discarded_by_the_cap_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            records_with_reads_the_cap_discarded: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    fn an_md5_is_sixteen_bytes_of_32_hex_characters_and_nothing_else() {
        assert_eq!(
            as_an_md5("000102030405060708090a0b0c0d0e0f"),
            Some([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        );
        assert_eq!(
            as_an_md5(""),
            None,
            "an older production file leaves it empty"
        );
        assert_eq!(as_an_md5("0102"), None, "too short");
        assert_eq!(
            as_an_md5("0g0102030405060708090a0b0c0d0e0f"),
            None,
            "32 characters, not all of them hex"
        );
    }

    #[test]
    fn a_recorded_input_is_a_basename_because_the_header_refuses_a_path() {
        assert_eq!(basename_of("tmp/d2/hg002_chr21.psp"), "hg002_chr21.psp");
        assert_eq!(basename_of("hg002_chr21.psp"), "hg002_chr21.psp");
    }
}
