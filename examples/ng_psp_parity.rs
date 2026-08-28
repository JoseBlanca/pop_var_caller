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
//! - **Exactly, on both pairs**: the coordinate and the extent, the reference bases, every
//!   allele sequence, every support count (`num_obs`, `num_fwd`, `mapq_sum`, `mapq_sum_sq`,
//!   `placed_left`), and **every chain-id list**.
//! - **Exactly, on the round-trip pair alone**, because production's record has no counterpart
//!   for them: the read witness, the read group, the two counts of reads that showed nothing,
//!   and the locus kind — which is `Generic` on every record of both corpora, so no parity run
//!   has yet moved the `Ssr` arm of that codec.
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
//! production's record has no counterpart for **four** of ng's fields (arch §2): the read
//! witness, the read group, and the two counts of reads that showed nothing. **All four are
//! synthesised here rather than left at their empty values**, and that is deliberate: a corpus
//! where every witness is `Complete`, every read group is 0 and both counts are 0 cannot tell an
//! encoder that stores those fields from one that drops them — measured, an encoder writing a
//! constant 0 for `reads_discarded_by_cap` passed a 3,000-record run before they were
//! synthesised. They are derived from each record's own coordinate, so the corpus is a function
//! of the input alone and two runs give the same file.
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
    LocusKind, LocusLen, ReadWitness, SampleLocusObservations, SequenceObservation,
    WitnessedLocusPositions,
};
use pop_var_caller::ng::psp::{
    ContigIdentity, DEFAULT_LOOK_BACK_WINDOW_LOG, FORMAT_VERSION, Header, Manifest, ParameterValue,
    PspReader, PspWriter, ReferenceIdentity, WriterProvenance, ZSTD_COMPRESSION_LEVEL,
};
use pop_var_caller::ng::types::{
    Bp, ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError,
};
use pop_var_caller::pileup_record::{ChainId, PileupRecord};
use pop_var_caller::psp::PspReader as ProductionPspReader;

// **The prototype, included whole and unchanged.** An example cannot depend on another example,
// and copying its encoder in would make the oracle a fork of itself — the one property it has
// that this module's own tests do not is that it was written by someone who was not writing this
// module. `#[path]` keeps one copy.
//
// **`dead_code` is this harness's doing, not the prototype's.** The include pulls in all 2,047
// lines and this file calls four of them, `fn main` among the rest. The prototype's own two lint
// findings are allowed on the items they belong to, in that file, where they cover its own
// target as well — CI runs `clippy --all-targets`, which builds it either way.
//
// Remove this the day the prototype is retired; nothing here should outlive it.
#[allow(dead_code)]
#[path = "psp_row_stream_roundtrip.rs"]
mod prototype;

/// How many bytes of a block the prototype's cut is allowed before it closes one early. Large
/// enough that the genomic grid is what cuts every block, which is what ng's writer does.
const PROTOTYPE_BLOCK_BYTES: usize = 1 << 20;

/// How many read groups the corpus cycles through.
///
/// **More than one, and the tally array is sized from this**: a corpus whose every observation
/// sits in one read group cannot fail on an encoder that drops the field, and a modulus widened
/// without the array would silently stop counting the groups past its end.
const READ_GROUPS_IN_THE_CORPUS: usize = 4;

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
fn as_an_ng_record_with_synthesised_fields(
    record: &PileupRecord,
) -> Option<SampleLocusObservations> {
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
            read_group: ReadGroupId(
                ((start + at as u64) % READ_GROUPS_IN_THE_CORPUS as u64) as u32,
            ),
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
        // **Two moduli, coprime**, so the two counts do not move together: an encoder writing
        // one of them into the other's slot is still caught.
        reads_without_observation: (start % 7) as u32,
        reads_discarded_by_cap: (start % 3) as u32,
        kind: LocusKind::Generic,
    })
}

/// A witness for one observation, chosen from the record's own coordinate so that the corpus is
/// a function of its input.
///
/// **Four shapes, and each is a different thing for the codec to get wrong**: `Complete`, which
/// is a run count of zero on the wire; a `Partial` covering the whole locus, which must **not**
/// come back as `Complete`; a `Partial` with a hole in it, which is what a spliced read across a
/// widened record looks like; and a single run flush with **neither** border.
///
/// ⚠ **The fourth shape was added by the H1 review, and the reason is exact.** The first three
/// are all reproduced by a witness stored as one prefix length and one suffix length — the hole
/// among them, since its two runs sit on the two borders — so a corpus holding only those cannot
/// fail on a codec that keeps the outermost edges and throws the interior away. A run flush with
/// neither border is the shape that cannot be spelled that way, and
/// [`CorpusShape::interior_witnesses`] counts it so the precondition can require one.
fn a_witness_for(start: u64, at: usize, span: u64) -> ReadWitness {
    let locus_len = u16::try_from(span.min(u64::from(u16::MAX))).expect("clamped to u16");
    if locus_len == 0 {
        return ReadWitness::Complete;
    }
    let runs: Vec<(u16, u16)> = match (start + at as u64) % 4 {
        0 => return ReadWitness::Complete,
        1 => vec![(0, locus_len)],
        2 if locus_len >= 3 => vec![(0, 1), (locus_len - 1, locus_len)],
        _ if locus_len >= 3 => vec![(1, locus_len - 1)],
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
            genomic_block_size_bp: Bp(u64::from(GRID_BP)),
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

/// What writing the ng store produced.
///
/// **Three counts of the same shape, so they are named rather than returned positionally**:
/// `blocks_cut` is what the single-block precondition gates on, and a transposition with
/// `records_pushed` would turn that guard into a tautology on any real corpus.
struct WhatTheNgStoreWriteProduced {
    records_pushed: u64,
    records_with_no_ng_record: u64,
    blocks_cut: u64,
    /// The header the writer was handed, kept so the reader's copy can be compared against it.
    header_written: Header,
}

/// Write the ng store from the source corpus.
///
/// **`take(limit)` and not a trailing break**, so this pass and the comparison pass consume the
/// same source records: the two used to implement the limit differently and disagreed by one
/// record at `--limit 0`, which then surfaced as *the ng store holds records the source does
/// not* — a data-integrity message for an argument-parsing difference.
fn write_the_ng_store(source: &Path, out: &Path, limit: usize) -> WhatTheNgStoreWriteProduced {
    let mut production = open_production(source);
    let header = an_ng_header(source, &production);
    let mut writer = PspWriter::create(out, header.clone()).expect("create the ng store");
    let (mut pushed, mut skipped) = (0u64, 0u64);
    for record in production.records().take(limit) {
        let record = record.expect("a production record, while writing the ng store");
        match as_an_ng_record_with_synthesised_fields(&record) {
            Some(record) => {
                writer.push(&record).unwrap_or_else(|why| {
                    panic!(
                        "pushing the record at {}:{}: {why}",
                        record.region.contig.0, record.region.start.0
                    )
                });
                pushed += 1;
            }
            None => skipped += 1,
        }
    }
    let stats = writer.finish(THE_TRAILER).expect("finish the ng store");
    WhatTheNgStoreWriteProduced {
        records_pushed: pushed,
        records_with_no_ng_record: skipped,
        blocks_cut: stats.blocks,
        header_written: header,
    }
}

/// What the ng store's trailer holds. **Opaque to the format** (spec §3.4); here it is a marker
/// the comparison reads back, because nothing else in the harness would notice a trailer that
/// came back wrong.
const THE_TRAILER: &[u8] = b"ng_psp_parity";

/// The source of records, opened as production's own reader.
///
/// **The two failures here name the file and say what to do about it.** The harness's whole
/// input is a file someone else wrote, possibly months ago, so the likeliest failure is a psp
/// written before the column set this build requires — which arrives as a serde error with the
/// entire TOML header inlined and no path in it.
fn open_production(path: &Path) -> ProductionPspReader<BufReader<File>> {
    let file = File::open(path)
        .unwrap_or_else(|why| panic!("cannot open {} as the source corpus: {why}", path.display()));
    ProductionPspReader::new(BufReader::with_capacity(1 << 20, file)).unwrap_or_else(|why| {
        panic!(
            "{} is not a production .psp this build can read — one written before the current \
             column set fails here and has to be regenerated: {why}",
            path.display()
        )
    })
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
    /// Witnesses holding a run flush with **neither** border — the shape a witness stored as one
    /// prefix length and one suffix length cannot express. See [`a_witness_for`].
    interior_witnesses: u64,
    read_groups: [u64; READ_GROUPS_IN_THE_CORPUS],
    records_with_reads_that_showed_nothing: u64,
    records_with_reads_the_cap_discarded: u64,
    blocks: u64,
}

impl CorpusShape {
    /// Tally one observation's witness, by the shape it takes on the wire.
    fn note_the_witness(&mut self, witness: &ReadWitness, locus_len: LocusLen) {
        match witness {
            ReadWitness::Complete => self.complete_witnesses += 1,
            ReadWitness::Partial { positions } => {
                self.partial_witnesses += 1;
                if positions.runs().len() > 1 {
                    self.witnesses_with_a_hole += 1;
                }
                if !positions.is_flush_left() || !positions.is_flush_right(locus_len) {
                    self.interior_witnesses += 1;
                }
            }
        }
    }

    /// Tally one observation's read group. **Panics rather than dropping one out of range**: a
    /// modulus widened past the array would otherwise weaken the read-group precondition without
    /// saying so.
    fn note_the_read_group(&mut self, group: ReadGroupId) {
        let group = group.get() as usize;
        let groups = self.read_groups.len();
        *self.read_groups.get_mut(group).unwrap_or_else(|| {
            panic!("the corpus minted read group {group} and the tally holds {groups}")
        }) += 1;
    }
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
            "--limit" => {
                let value = args
                    .next()
                    .expect("--limit needs a record count, e.g. --limit 300000");
                limit = value
                    .parse()
                    .unwrap_or_else(|why| panic!("--limit {value:?} is not a record count: {why}"));
            }
            "--work" => {
                work = PathBuf::from(args.next().expect("--work needs a directory to write into"))
            }
            "--label" => {
                label = args.next().expect(
                    "--label needs a name; it names both stores and the report's corpus line",
                )
            }
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(
        limit > 0,
        "--limit needs at least one record; 0 compares nothing"
    );
    std::fs::create_dir_all(&work).expect("create the work directory");
    let ng_store = work.join(format!("{label}.ngpsp"));
    let prototype_store = work.join(format!("{label}.ngs"));

    let written = write_the_ng_store(&source, &ng_store, limit);
    assert!(
        written.records_pushed > 0,
        "the corpus produced no ng records at all"
    );

    // **The prototype is asked for ng's own step**, so the two stores must then agree on the
    // summed log-error exactly rather than inside a tolerance wide enough to hold both.
    //
    // **Every field spelled out, no `..default()`.** This is the load-bearing setting of the
    // whole cross-arm comparison, and a field added to `Scales` would otherwise be quantised at
    // the prototype's default here while ng used its own — narrowing the equality to `q_sum`
    // with no compile error. The prototype's own two literals spell every field for the same
    // reason; measured by the H1 review, adding a fourth field flagged those two and not this.
    let prototype_defaults = prototype::Scales::default();
    let scales = prototype::Scales {
        gc: prototype_defaults.gc,
        coverage: prototype_defaults.coverage,
        q_sum: SummedLogError::STEPS_PER_NAT as f64,
    };
    // The prototype's three format switches, named here because its signature takes them
    // positionally: the store carries every field rather than the light subset, does not
    // length-prefix its records, and does write a record head — which is the shape ng writes.
    let light_only = false;
    let length_prefix = false;
    let record_head = true;
    prototype::encode_streaming(
        source.to_str().expect("a utf-8 path"),
        prototype_store.to_str().expect("a utf-8 path"),
        PROTOTYPE_BLOCK_BYTES,
        GRID_BP,
        light_only,
        length_prefix,
        record_head,
        ZSTD_COMPRESSION_LEVEL,
        u32::from(DEFAULT_LOOK_BACK_WINDOW_LOG),
        scales,
    );

    let (shape, worst) =
        walk_the_three_streams_in_lockstep(&source, &ng_store, &prototype_store, limit, &written);
    assert_the_run_proves_something(&shape, &worst);
    print_the_report(&label, &work, &shape, &worst, &written);
}

/// Walk the source, the ng store and the prototype's store together, and fail on the first
/// record that disagrees.
fn walk_the_three_streams_in_lockstep(
    source: &Path,
    ng_store: &Path,
    prototype_store: &Path,
    limit: usize,
    written: &WhatTheNgStoreWriteProduced,
) -> (CorpusShape, WorstDrift) {
    let mut production = open_production(source);
    let mut ng = PspReader::open(ng_store).expect("open the ng store");
    let mut proto = prototype::StreamingStore::open(prototype_store);

    check_everything_the_store_holds_that_is_not_a_record(&mut ng, written);

    // **The block count comes from the file's own index, not from the writer's tally.** The
    // single-block precondition exists to guarantee the corpus crosses a boundary, and a writer
    // that miscounted — or an index that recorded fewer blocks than were cut — would satisfy it
    // on a file that crosses none. The writer's number is kept as a cross-check.
    let blocks_in_the_file = ng.block_index().len() as u64;
    assert_eq!(
        blocks_in_the_file, written.blocks_cut,
        "the writer says it cut {} blocks and the file's index holds {blocks_in_the_file}",
        written.blocks_cut
    );
    let mut shape = CorpusShape {
        blocks: blocks_in_the_file,
        ..CorpusShape::default()
    };
    let mut worst = WorstDrift::default();

    let mut ng_records = ng.records().expect("walk the ng store");
    let source_records = production.records();
    for record in source_records.take(limit) {
        let from_production = record.expect("a production record, while comparing");
        // The prototype writes every source record, so its stream advances even for the ones
        // that have no ng record — otherwise the two walks would slide apart by the count of
        // skipped records and every later comparison would be against the wrong record.
        let from_prototype = proto.next().expect("the prototype's store ended early");
        let Some(pushed) = as_an_ng_record_with_synthesised_fields(&from_production) else {
            continue;
        };
        let read_back = ng_records
            .next()
            .expect("the ng store ended early")
            .expect("an ng record")
            .record
            .expect("a walk that declines no body builds every one");
        let at = shape.records;
        check_the_ng_round_trip(at, &pushed, &read_back, &mut shape);
        check_against_the_prototype(at, &read_back, &from_prototype);
        note_the_drift(&from_production, &read_back, &from_prototype, &mut worst);
        shape.records += 1;
    }
    assert!(
        ng_records.next().is_none(),
        "the ng store holds records the source does not"
    );
    // **Both stores are checked for a surplus, not just ng's.** Under `--limit` the prototype
    // legitimately holds the rest of the file, so this only applies to a whole run.
    if limit == usize::MAX {
        assert!(
            proto.next().is_none(),
            "the prototype's store holds records the source does not"
        );
    }
    (shape, worst)
}

/// The header, the trailer and the block index — **everything the store holds that is not a
/// record**.
///
/// ⚠ **Nothing read any of this until the H1 review**, and the gap was measured rather than
/// argued: with defects planted in `src/ng/psp/`, a contig length written one too long, a block
/// index whose every first position was a base too far, and a trailer truncated by one byte all
/// passed a 100,000-record run unremarked. The cohort merge reaches for all three — the contig
/// list to agree on a reference, the index to seek, the trailer for the per-sample summary — so
/// an oracle that only walks records is not measuring the store a run will use.
fn check_everything_the_store_holds_that_is_not_a_record(
    ng: &mut PspReader,
    written: &WhatTheNgStoreWriteProduced,
) {
    // **Destructured with no `..`**, so a field added to `Header` has to be classified here
    // rather than dropping out of this check.
    let Header {
        format_version,
        sample,
        reference,
        contigs,
        writer,
        manifest,
    } = ng.header();
    let handed = &written.header_written;
    assert_eq!(
        *format_version, handed.format_version,
        "the format version read back is not the one written"
    );
    assert_eq!(*sample, handed.sample, "the sample name read back");
    assert_eq!(*reference, handed.reference, "the reference read back");
    assert_eq!(
        *contigs, handed.contigs,
        "the contig list read back — every name, length and md5"
    );
    assert_eq!(*manifest, handed.manifest, "the manifest read back");
    // **The provenance is compared in two halves, because `create` adds to it.** It records the
    // compression level the blocks were written at (`writer.rs`, F3), so the header in the file
    // is deliberately not the header it was handed — everything else must survive unchanged, and
    // every parameter handed in must still be there.
    assert_eq!(
        writer.tool, handed.writer.tool,
        "the writer provenance read back: tool"
    );
    assert_eq!(
        writer.version, handed.writer.version,
        "the writer provenance read back: version"
    );
    assert_eq!(
        writer.subcommand, handed.writer.subcommand,
        "the writer provenance read back: subcommand"
    );
    assert_eq!(
        writer.input_alignments, handed.writer.input_alignments,
        "the writer provenance read back: input alignments"
    );
    assert_eq!(
        writer.input_reference, handed.writer.input_reference,
        "the writer provenance read back: input reference"
    );
    assert_eq!(
        writer.command_line, handed.writer.command_line,
        "the writer provenance read back: command line"
    );
    assert_eq!(
        writer.created, handed.writer.created,
        "the writer provenance read back: creation stamp"
    );
    for (name, value) in &handed.writer.parameters {
        assert_eq!(
            writer.parameters.get(name),
            Some(value),
            "the parameter {name:?} the writer was handed is not in the file's header"
        );
    }
    assert!(
        writer.parameters.contains_key("zstd-compression-level"),
        "`create` records the level its blocks were written at, and the header read back does \
         not carry it"
    );
    assert_eq!(
        ng.trailer().expect("read the trailer back"),
        THE_TRAILER,
        "the trailer read back is not the one `finish` was given"
    );
    // **The index is checked two ways, and only the second has teeth against a wrong
    // coordinate.** Entering by ordinal must give a suffix of the whole walk — that holds the
    // offsets. But `records_from` searches the index on the *coordinate*, so an index whose
    // every entry named a position one base too far would pass the suffix check untouched:
    // measured, that defect passed a 20,000-record run until the second assertion below existed.
    let entries: Vec<_> = ng.block_index().to_vec();
    let every_position: Vec<(u32, u64)> = ng
        .records()
        .expect("walk the ng store for its coordinates")
        .map(|record| {
            let region = record.expect("an ng record").head.region;
            (region.contig.0, region.start.0)
        })
        .collect();
    for (ordinal, entry) in entries.iter().enumerate() {
        let from_here: Vec<(u32, u64)> = ng
            .records_from_block(ordinal)
            .expect("enter the store at a block the index names")
            .map(|record| {
                let region = record.expect("an ng record").head.region;
                (region.contig.0, region.start.0)
            })
            .collect();
        assert!(
            every_position.ends_with(&from_here),
            "block {ordinal}, which the index puts at {}:{}, does not yield a suffix of the \
             whole walk",
            entry.first_position.contig.0,
            entry.first_position.position.0
        );
        let (contig, position) = *from_here
            .first()
            .expect("a block the index names holds at least one record");
        assert_eq!(
            (
                entry.first_position.contig.0,
                entry.first_position.position.0
            ),
            (contig, position),
            "the index says block {ordinal} begins at {}:{}, and its first record is at \
             {contig}:{position}",
            entry.first_position.contig.0,
            entry.first_position.position.0
        );
    }
}

/// Every field of the record that was pushed, against the record that came back. **Exactly**,
/// field by field — the summed log-error included, since it is an integer count of steps by the
/// time the store sees it.
///
/// **Both records are destructured with no `..`, and that is load-bearing.** A field added to
/// [`SampleLocusObservations`] or [`SequenceObservation`] is then a compile error here rather
/// than a field this oracle silently stops comparing — the same discipline
/// `encode_record_body` keeps in `src/ng/psp/record.rs`, and for the same reason. The H1
/// review found the drift had already happened once on the other arm.
fn check_the_ng_round_trip(
    record_ordinal: u64,
    pushed: &SampleLocusObservations,
    read_back: &SampleLocusObservations,
    shape: &mut CorpusShape,
) {
    let at = record_ordinal;
    let SampleLocusObservations {
        region: pushed_region,
        reference_bases: pushed_reference_bases,
        observations: pushed_observations,
        reads_without_observation: pushed_reads_without_observation,
        reads_discarded_by_cap: pushed_reads_discarded_by_cap,
        kind: pushed_kind,
    } = pushed;
    let SampleLocusObservations {
        region: read_back_region,
        reference_bases: read_back_reference_bases,
        observations: read_back_observations,
        reads_without_observation: read_back_reads_without_observation,
        reads_discarded_by_cap: read_back_reads_discarded_by_cap,
        kind: read_back_kind,
    } = read_back;

    assert_eq!(read_back_region, pushed_region, "record {at}: region");
    assert_eq!(
        read_back_reference_bases, pushed_reference_bases,
        "record {at}: reference bases"
    );
    assert_eq!(read_back_kind, pushed_kind, "record {at}: locus kind");
    assert_eq!(
        read_back_reads_without_observation, pushed_reads_without_observation,
        "record {at}: reads without observation"
    );
    assert_eq!(
        read_back_reads_discarded_by_cap, pushed_reads_discarded_by_cap,
        "record {at}: reads discarded by the cap"
    );
    if *pushed_reads_without_observation > 0 {
        shape.records_with_reads_that_showed_nothing += 1;
    }
    if *pushed_reads_discarded_by_cap > 0 {
        shape.records_with_reads_the_cap_discarded += 1;
    }
    assert_eq!(
        read_back_observations.len(),
        pushed_observations.len(),
        "record {at}: observation count"
    );
    if pushed_observations.len() > 1 {
        shape.with_several_observations += 1;
    }

    let mut any_chain_ids = false;
    for (which, (read_back_observation, pushed_observation)) in read_back_observations
        .iter()
        .zip(pushed_observations)
        .enumerate()
    {
        let SequenceObservation {
            bases: pushed_bases,
            read_witness: pushed_read_witness,
            read_group: pushed_read_group,
            num_obs: pushed_num_obs,
            num_fwd: pushed_num_fwd,
            q_sum: pushed_q_sum,
            mapq_sum: pushed_mapq_sum,
            mapq_sum_sq: pushed_mapq_sum_sq,
            placed_left: pushed_placed_left,
            chain_ids: pushed_chain_ids,
        } = pushed_observation;
        let SequenceObservation {
            bases: read_back_bases,
            read_witness: read_back_read_witness,
            read_group: read_back_read_group,
            num_obs: read_back_num_obs,
            num_fwd: read_back_num_fwd,
            q_sum: read_back_q_sum,
            mapq_sum: read_back_mapq_sum,
            mapq_sum_sq: read_back_mapq_sum_sq,
            placed_left: read_back_placed_left,
            chain_ids: read_back_chain_ids,
        } = read_back_observation;

        assert_eq!(
            read_back_bases, pushed_bases,
            "record {at}, observation {which}: bases"
        );
        assert_eq!(
            read_back_read_witness, pushed_read_witness,
            "record {at}, observation {which}: read witness"
        );
        assert_eq!(
            read_back_read_group, pushed_read_group,
            "record {at}, observation {which}: read group"
        );
        assert_eq!(
            read_back_num_obs, pushed_num_obs,
            "record {at}, observation {which}: num_obs"
        );
        assert_eq!(
            read_back_num_fwd, pushed_num_fwd,
            "record {at}, observation {which}: num_fwd"
        );
        assert_eq!(
            read_back_q_sum.steps(),
            pushed_q_sum.steps(),
            "record {at}, observation {which}: summed log-error"
        );
        assert_eq!(
            read_back_mapq_sum, pushed_mapq_sum,
            "record {at}, observation {which}: mapq_sum"
        );
        assert_eq!(
            read_back_mapq_sum_sq, pushed_mapq_sum_sq,
            "record {at}, observation {which}: mapq_sum_sq"
        );
        assert_eq!(
            read_back_placed_left, pushed_placed_left,
            "record {at}, observation {which}: placed_left"
        );
        let wanted = as_a_read_set(pushed_chain_ids);
        if &wanted != pushed_chain_ids {
            shape.lists_reordered_by_normalising += 1;
        }
        assert_eq!(
            read_back_chain_ids, &wanted,
            "record {at}, observation {which}: chain ids"
        );

        shape.observations += 1;
        shape.chain_ids += wanted.len() as u64;
        any_chain_ids |= !wanted.is_empty();
        shape.note_the_witness(pushed_read_witness, LocusLen::of_region(*pushed_region));
        shape.note_the_read_group(*pushed_read_group);
    }
    if any_chain_ids {
        shape.with_chain_ids += 1;
    }
}

/// The ng store's record against the prototype's, on **every field the prototype's record has**.
///
/// **Four fields of ng's record have no counterpart here** and are held by
/// [`check_the_ng_round_trip`] alone: the read witness, the read group, and the two counts of
/// reads that showed nothing. Everything else production carries is compared, and the
/// observation is destructured with no `..` so that stays true — the four are named as `_` here
/// rather than skipped silently.
///
/// ⚠ **`placed_left`, the reference bases and the record's extent were missing from this arm
/// until the H1 review**, and the omission was not visible: a defect planted in this harness's
/// own reading of `placed_left`, and one making every region a base too long, both passed a
/// 74,623-record run with a clean report. That is exactly the class of defect this second arm
/// exists to catch, since the round-trip arm carries such a defect on both sides.
fn check_against_the_prototype(
    record_ordinal: u64,
    read_back: &SampleLocusObservations,
    from_prototype: &PileupRecord,
) {
    let at = record_ordinal;
    assert_eq!(
        read_back.region.contig.0, from_prototype.chrom_id,
        "record {at}: contig, against the prototype"
    );
    assert_eq!(
        read_back.region.start.0,
        u64::from(from_prototype.pos),
        "record {at}: position, against the prototype"
    );
    // The prototype has no extent of its own — production anchors a record at one position — so
    // the thing to compare is the extent this harness *derived*, against the reference allele it
    // derived it from. Without this, a mapping that widened every record by a base passes.
    assert_eq!(
        read_back.region.len(),
        from_prototype
            .alleles
            .first()
            .map_or(0, |allele| allele.seq.len() as u64),
        "record {at}: reference span, against the prototype"
    );
    assert_eq!(
        read_back.reference_bases.as_ref(),
        from_prototype
            .alleles
            .first()
            .map_or(&[][..], |allele| allele.seq.as_slice()),
        "record {at}: reference bases, against the prototype"
    );
    assert_eq!(
        read_back.observations.len(),
        from_prototype.alleles.len(),
        "record {at}: allele count, against the prototype"
    );
    for (which, (observation, allele)) in read_back
        .observations
        .iter()
        .zip(&from_prototype.alleles)
        .enumerate()
    {
        let SequenceObservation {
            bases,
            // The four production has no equivalent of. Named, not skipped with `..`, so that a
            // field added to the type has to be classified here rather than silently dropping
            // out of this arm.
            read_witness: _,
            read_group: _,
            num_obs,
            num_fwd,
            q_sum,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
            chain_ids,
        } = observation;
        assert_eq!(
            bases.as_ref(),
            allele.seq.as_slice(),
            "record {at}, allele {which}: sequence, against the prototype"
        );
        assert_eq!(
            chain_ids,
            &as_a_read_set(&allele.chain_ids),
            "record {at}, allele {which}: chain ids, against the prototype"
        );
        assert_eq!(
            *num_obs, allele.support.num_obs,
            "record {at}, allele {which}: num_obs, against the prototype"
        );
        assert_eq!(
            *num_fwd, allele.support.fwd,
            "record {at}, allele {which}: num_fwd, against the prototype"
        );
        assert_eq!(
            *mapq_sum, allele.support.mapq_sum,
            "record {at}, allele {which}: mapq_sum, against the prototype"
        );
        assert_eq!(
            *mapq_sum_sq, allele.support.mapq_sum_sq,
            "record {at}, allele {which}: mapq_sum_sq, against the prototype"
        );
        assert_eq!(
            *placed_left, allele.support.placed_left,
            "record {at}, allele {which}: placed_left, against the prototype"
        );
        // Exact, because the prototype was handed ng's own step: both stores round the same
        // `f64` to the same integer count of 1/4,096 of a natural log, so a difference here is
        // one of them storing a different number, not one of them rounding differently.
        assert_eq!(
            q_sum.steps(),
            (allele.support.q_sum * SummedLogError::STEPS_PER_NAT as f64).round() as i64,
            "record {at}, allele {which}: summed log-error, against the prototype"
        );
    }
}

/// How far each store's summed log-error ended up from the `f64` production held. Reported, not
/// asserted per record — **the ng store's** distance is asserted against half a step once, at the
/// end; the prototype's is printed beside it and not asserted.
///
/// **The two `PileupRecord`s are named rather than passed as a positional pair**: they are the
/// same type, and transposing them would report a prototype drift of exactly zero for every
/// record — a *passing* number, so nothing would flag it.
fn note_the_drift(
    from_production: &PileupRecord,
    from_the_ng_store: &SampleLocusObservations,
    from_the_prototype: &PileupRecord,
    worst: &mut WorstDrift,
) {
    for (which, allele) in from_production.alleles.iter().enumerate() {
        if let Some(observation) = from_the_ng_store.observations.get(which) {
            let drift = (observation.q_sum.nats() - allele.support.q_sum).abs();
            worst.ng_against_the_source = worst.ng_against_the_source.max(drift);
        }
        if let Some(other) = from_the_prototype.alleles.get(which) {
            let drift = (other.support.q_sum - allele.support.q_sum).abs();
            worst.prototype_against_the_source = worst.prototype_against_the_source.max(drift);
        }
    }
}

// ---------------------------------------------------------------------
// The report, and the preconditions it rests on
// ---------------------------------------------------------------------

/// Half of one step of the summed log-error — the most that rounding to the nearest step can
/// cost, and therefore the bound the ng store's drift from the source `f64` must sit inside.
fn half_a_step() -> f64 {
    0.5 / SummedLogError::STEPS_PER_NAT as f64
}

/// **The two claims a passing run makes**, kept out of the printing so that neither can be lost
/// behind a later `--quiet`: that the corpus holds every shape the comparison needs in order to
/// be *able* to fail, and that ng's summed log-error never moved further from production's `f64`
/// than rounding to the nearest step can cost.
fn assert_the_run_proves_something(shape: &CorpusShape, worst: &WorstDrift) {
    assert_the_corpus_can_fail(shape);
    assert!(
        worst.ng_against_the_source <= half_a_step(),
        "the ng store's summed log-error is {} from the source, and rounding to the nearest \
         1/{} of a natural log cannot cost more than {}",
        worst.ng_against_the_source,
        SummedLogError::STEPS_PER_NAT,
        half_a_step()
    );
}

/// What the run measured, one `key<TAB>value` line each.
///
/// **The prototype prints its own `encode-streaming` block above this one**, with `records` and
/// `blocks` lines of its own that mean different things — its counts are over the whole source
/// file, these are over the compared prefix. Every line here is prefixed `ng-` where the two
/// could be confused.
fn print_the_report(
    label: &str,
    work: &Path,
    shape: &CorpusShape,
    worst: &WorstDrift,
    written: &WhatTheNgStoreWriteProduced,
) {
    let half_a_step = half_a_step();
    println!("phase\tng-psp-parity");
    println!("corpus\t{label}");
    println!("work-dir\t{}", work.display());
    println!("genomic-grid-bp\t{GRID_BP}");
    println!("records-compared\t{}", shape.records);
    println!("records-pushed\t{}", written.records_pushed);
    println!(
        "records-with-no-ng-record\t{}",
        written.records_with_no_ng_record
    );
    println!("ng-blocks\t{}", shape.blocks);
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
    println!(
        "witnesses-flush-with-neither-border\t{}",
        shape.interior_witnesses
    );
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
    // **This equals the line above by construction**, because the prototype was handed ng's own
    // step, so both stores round the same `f64` the same way. It is printed anyway: the day the
    // two differ, the equality the cross-arm comparison rests on has stopped holding, and this
    // is the line that says so.
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
        "every record in the corpus holds one observation, so the format's cheapest chain-id \
         path — one observation's read list derived as the live set minus the others' — is \
         never taken"
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
        "no witness in the corpus has more than one run, so a codec that stores one run passes"
    );
    assert!(
        shape.interior_witnesses > 0,
        "no witness in the corpus has a run flush with neither border, so a codec that keeps \
         one prefix length and one suffix length reproduces every witness here and passes"
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
            interior_witnesses: 2,
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

    /// **Allele 1, not allele 0.** The reference bases are compared against allele 0's sequence
    /// before the loop, so changing that one is refused as a reference-bases difference and this
    /// test would pass without the per-allele sequence comparison existing at all.
    #[test]
    #[should_panic(expected = "allele 1: sequence, against the prototype")]
    fn an_allele_sequence_the_prototype_disagrees_with_is_refused() {
        let got = a_record();
        let mut proto = the_same_record_as_production_holds_it();
        proto.alleles[1].seq = b"ACGA".to_vec();
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
    #[should_panic(expected = "cheapest chain-id path")]
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
    #[should_panic(expected = "has more than one run")]
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

    // -----------------------------------------------------------------
    // Every comparison, proved able to fail
    // -----------------------------------------------------------------

    /// Run `check` on a perturbed record and give back the message it refused with, or `None` if
    /// it accepted the change.
    fn what_it_refuses_with(check: impl FnOnce() + std::panic::UnwindSafe) -> Option<String> {
        std::panic::catch_unwind(check).err().map(|payload| {
            payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|it| (*it).to_string()))
                .unwrap_or_else(|| "a panic carrying no message".to_string())
        })
    }

    /// **One row per comparison in the round-trip arm, each proved able to fail.**
    ///
    /// ⚠ **This exists because the individual `should_panic` tests above did not cover the arm.**
    /// The H1 review neutered each comparison in turn — rewriting `assert_eq!(a.x, b.x)` as
    /// `assert_eq!(a.x, a.x)` — and found that **eighteen of the twenty-six comparisons in the
    /// two arms survived with all 23 tests green**, among them the region, the reference bases,
    /// the locus kind, `num_obs`, `mapq_sum` and `placed_left`. A comparison nothing can fail is
    /// a comparison that could be deleted, after which the harness reports a clean run over eight
    /// million records and proves nothing about that field.
    ///
    /// The table is the guard *and* the inventory: a comparison added without a row here is a
    /// comparison nothing holds.
    #[test]
    fn every_comparison_in_the_round_trip_arm_can_fail() {
        type Perturb = fn(&mut SampleLocusObservations);
        let rows: &[(&str, Perturb)] = &[
            ("record 0: region", |r| {
                r.region.end = Position(r.region.end.0 + 1)
            }),
            ("record 0: reference bases", |r| {
                r.reference_bases = b"AAAA".to_vec().into_boxed_slice()
            }),
            ("record 0: locus kind", |r| r.kind = LocusKind::SsrBundle),
            ("record 0: reads without observation", |r| {
                r.reads_without_observation += 1
            }),
            ("record 0: reads discarded by the cap", |r| {
                r.reads_discarded_by_cap += 1
            }),
            ("record 0: observation count", |r| {
                r.observations.pop();
            }),
            ("observation 0: bases", |r| {
                r.observations[0].bases = b"TTTT".to_vec().into_boxed_slice()
            }),
            ("observation 1: read witness", |r| {
                r.observations[1].read_witness = ReadWitness::Complete
            }),
            ("observation 1: read group", |r| {
                r.observations[1].read_group = ReadGroupId(0)
            }),
            ("observation 0: num_obs", |r| r.observations[0].num_obs += 1),
            ("observation 0: num_fwd", |r| r.observations[0].num_fwd += 1),
            ("observation 0: summed log-error", |r| {
                r.observations[0].q_sum = SummedLogError::from_steps(-1)
            }),
            ("observation 0: mapq_sum", |r| {
                r.observations[0].mapq_sum += 1
            }),
            ("observation 0: mapq_sum_sq", |r| {
                r.observations[0].mapq_sum_sq += 1
            }),
            ("observation 0: placed_left", |r| {
                r.observations[0].placed_left += 1
            }),
            ("observation 0: chain ids", |r| {
                r.observations[0].chain_ids.pop();
            }),
        ];
        for (names, perturb) in rows {
            let pushed = a_record();
            let mut read_back = pushed.clone();
            perturb(&mut read_back);
            assert_ne!(
                read_back, pushed,
                "the row for {names:?} did not change the record, so it proves nothing"
            );
            let refused = what_it_refuses_with(|| {
                check_the_ng_round_trip(0, &pushed, &read_back, &mut CorpusShape::default())
            });
            let refused = refused.unwrap_or_else(|| {
                panic!("the round-trip arm accepted a record whose {names} was changed")
            });
            assert!(
                refused.contains(names),
                "the round-trip arm refused a changed {names}, but named it {refused:?}"
            );
        }
    }

    /// The same, for the arm that compares the two stores against each other.
    ///
    /// ⚠ **Three of these rows are the H1 review's Blocker**: `placed_left`, the reference bases
    /// and the reference span were not compared here at all, and the omission was invisible — a
    /// defect planted in this harness's own reading of `placed_left`, and one making every region
    /// a base too long, each passed a 74,623-record run with a clean report. The round-trip arm
    /// cannot catch either, because it carries the defect on both sides.
    #[test]
    fn every_comparison_in_the_prototype_arm_can_fail() {
        type Perturb = fn(&mut SampleLocusObservations, &mut PileupRecord);
        let rows: &[(&str, Perturb)] = &[
            ("record 0: contig", |_, p| p.chrom_id += 1),
            ("record 0: position", |_, p| p.pos += 1),
            ("record 0: reference span", |r, _| {
                r.region.end = Position(r.region.end.0 + 1)
            }),
            ("record 0: reference bases", |r, _| {
                r.reference_bases = b"AAAA".to_vec().into_boxed_slice()
            }),
            ("record 0: allele count", |r, _| {
                r.observations.pop();
            }),
            ("allele 1: sequence", |r, _| {
                r.observations[1].bases = b"AAAA".to_vec().into_boxed_slice()
            }),
            ("allele 0: chain ids", |_, p| {
                p.alleles[0].chain_ids.push(99)
            }),
            ("allele 0: num_obs", |_, p| {
                p.alleles[0].support.num_obs += 1
            }),
            ("allele 0: num_fwd", |_, p| p.alleles[0].support.fwd += 1),
            ("allele 0: mapq_sum", |_, p| {
                p.alleles[0].support.mapq_sum += 1
            }),
            ("allele 0: mapq_sum_sq", |_, p| {
                p.alleles[0].support.mapq_sum_sq += 1
            }),
            ("allele 0: placed_left", |_, p| {
                p.alleles[0].support.placed_left += 1
            }),
            ("allele 0: summed log-error", |_, p| {
                p.alleles[0].support.q_sum += 1.0 / SummedLogError::STEPS_PER_NAT as f64
            }),
        ];
        for (names, perturb) in rows {
            let mut read_back = a_record();
            let mut from_prototype = the_same_record_as_production_holds_it();
            perturb(&mut read_back, &mut from_prototype);
            let refused = what_it_refuses_with(|| {
                check_against_the_prototype(0, &read_back, &from_prototype)
            });
            let refused = refused.unwrap_or_else(|| {
                panic!("the prototype arm accepted a record whose {names} was changed")
            });
            assert!(
                refused.contains(names) && refused.contains("against the prototype"),
                "the prototype arm refused a changed {names}, but named it {refused:?}"
            );
        }
    }

    /// **A record that agrees is accepted by both tables' machinery**, so a `check` that panicked
    /// unconditionally would not be mistaken for one that discriminates.
    #[test]
    fn the_two_arms_accept_a_record_that_agrees() {
        let record = a_record();
        assert!(
            what_it_refuses_with(|| check_the_ng_round_trip(
                0,
                &record,
                &record,
                &mut CorpusShape::default()
            ))
            .is_none()
        );
        assert!(
            what_it_refuses_with(|| check_against_the_prototype(
                0,
                &record,
                &the_same_record_as_production_holds_it()
            ))
            .is_none()
        );
    }

    // -----------------------------------------------------------------
    // The fourth witness shape, and the precondition that requires it
    // -----------------------------------------------------------------

    /// **A run flush with neither border** — the shape a witness stored as one prefix length and
    /// one suffix length cannot express, and the reason `a_witness_for` has a fourth arm.
    #[test]
    fn the_fourth_witness_shape_is_flush_with_neither_border() {
        let ReadWitness::Partial { positions } = a_witness_for(3, 0, 5) else {
            panic!("the fourth shape is a partial witness");
        };
        assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(1, 4)]);
        assert!(
            !positions.is_flush_left() && !positions.is_flush_right(LocusLen::from_positions(5)),
            "the point of this shape is that it touches neither border"
        );
    }

    /// The three shapes the corpus minted before the H1 review are **all** reproduced by a
    /// witness kept as one prefix length and one suffix length — the hole among them, since its
    /// two runs sit on the two borders. That is why the fourth shape had to be added.
    #[test]
    fn the_first_three_witness_shapes_are_all_flush_with_a_border() {
        for witness in [
            a_witness_for(0, 0, 5),
            a_witness_for(1, 0, 5),
            a_witness_for(2, 0, 5),
        ] {
            let ReadWitness::Partial { positions } = witness else {
                continue; // `Complete` carries no runs at all
            };
            assert!(
                positions.is_flush_left() || positions.is_flush_right(LocusLen::from_positions(5)),
                "shape {positions:?} would have needed the fourth arm"
            );
        }
    }

    #[test]
    #[should_panic(expected = "no records were compared")]
    fn a_corpus_of_no_records_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            records: 0,
            ..a_corpus_that_can_fail()
        });
    }

    /// The other half of the witness precondition. A corpus of *only* partial witnesses is
    /// exactly as unable to fail as one of only complete ones, and only the first direction was
    /// tested.
    #[test]
    #[should_panic(expected = "it needs both")]
    fn a_corpus_whose_witnesses_are_all_partial_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            complete_witnesses: 0,
            ..a_corpus_that_can_fail()
        });
    }

    #[test]
    #[should_panic(expected = "flush with neither border")]
    fn a_corpus_with_no_interior_witness_is_refused() {
        assert_the_corpus_can_fail(&CorpusShape {
            interior_witnesses: 0,
            ..a_corpus_that_can_fail()
        });
    }

    /// **A read group past the tally's end is a panic, not a silent drop.** The corpus mints
    /// groups with `% READ_GROUPS_IN_THE_CORPUS`, so widening that modulus without widening the
    /// array would otherwise stop counting the new groups — and the read-group precondition
    /// would keep passing while measuring less than it claims.
    #[test]
    #[should_panic(expected = "the corpus minted read group 9")]
    fn a_read_group_past_the_tally_is_refused() {
        CorpusShape::default().note_the_read_group(ReadGroupId(9));
    }
}
