//! **Are a record head's four fields cheaper as variable-length integers or as fixed-width
//! ones — after the block has been compressed?**
//!
//! The container spec leaves the width to the manifest and records that a fixed width is
//! quicker to read and "costs less than it looks after compression", with the two never having
//! been compared *in place* (spec `psp_file_format.md` §4.3). The architecture document asks for
//! the comparison at implementation time rather than by argument (arch §7). It could not be
//! taken before a compressor existed; a compressor exists at Milestone D2, and this is it.
//!
//! ```text
//! cargo run --release --example ng_psp_head_encoding -- <a production .psp> [--label tomato]
//! ```
//!
//! # What the corpus is, and what it is not
//!
//! ng cannot yet produce a psp, so the records come from a **production `.psp`** and are turned
//! into ng records. That is faithful for three of the four head fields and approximate for the
//! fourth:
//!
//! - **position-offset** — exact. The coordinates are the real ones.
//! - **reference-span** — exact for these records, and *narrower than ng's will be*: production
//!   anchors a record at one position and derives its span from the reference allele, so almost
//!   every span is 1. ng widens a record across a deletion.
//! - **non-reference-reads** — exact. It is derived from the same allele support.
//! - **record-body-byte-count** — **approximate, and the approximation runs in both
//!   directions.** ng's body carries a read witness and a read group that production has no
//!   equivalent of, which makes these bodies smaller than ng's will be; and the chain ids are
//!   not encoded until Milestone E, which makes them smaller again — at depth, by a lot.
//!
//! **So the head's share of the file here is an overstatement**, because the bodies are
//! smaller than ng's finished ones. What the comparison is about — which of two encodings of
//! the same four values is cheaper — is not affected by that, since both arms see the same
//! bodies.

use std::fs::File;
use std::io::BufReader;

use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use pop_var_caller::ng::psp::{
    BlockCompressor, BlockHead, DEFAULT_GENOMIC_BLOCK_SIZE_BP, DEFAULT_LOOK_BACK_WINDOW_LOG,
    encode_record_body,
};
use pop_var_caller::ng::types::{
    Bp, ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError,
};
use pop_var_caller::pileup_record::PileupRecord;
use pop_var_caller::psp::PspReader;

/// How a run writes the four fields of a record's head.
#[derive(Clone, Copy)]
enum HeadEncoding {
    /// What the format writes today: LEB128, a byte for a small value.
    Varint,
    /// Every field at a declared fixed width, little-endian. The widths have to hold every
    /// value the *format* allows, not every value this corpus happens to contain.
    Fixed {
        offset: u8,
        span: u8,
        reads: u8,
        body: u8,
    },
}

impl HeadEncoding {
    fn name(self) -> String {
        match self {
            Self::Varint => "varint".to_string(),
            Self::Fixed {
                offset,
                span,
                reads,
                body,
            } => format!("fixed {offset}/{span}/{reads}/{body}"),
        }
    }
}

/// Write one head, or say which field did not fit its declared width.
fn put_head(
    out: &mut Vec<u8>,
    encoding: HeadEncoding,
    offset: u64,
    span: u64,
    reads: u64,
    body: u64,
) -> Result<(), &'static str> {
    match encoding {
        HeadEncoding::Varint => {
            // LEB128 by hand, so this harness needs nothing from a private module: seven
            // bits a byte, the top bit set while more follow.
            for mut value in [offset, span, reads, body] {
                while value >= 0x80 {
                    out.push((value as u8) | 0x80);
                    value >>= 7;
                }
                out.push(value as u8);
            }
            Ok(())
        }
        HeadEncoding::Fixed {
            offset: offset_bytes,
            span: span_bytes,
            reads: reads_bytes,
            body: body_bytes,
        } => {
            for (value, width, field) in [
                (offset, offset_bytes, "position-offset"),
                (span, span_bytes, "reference-span"),
                (reads, reads_bytes, "non-reference-reads"),
                (body, body_bytes, "record-body-byte-count"),
            ] {
                let room = if width >= 8 {
                    u64::MAX
                } else {
                    (1u64 << (8 * u32::from(width))) - 1
                };
                if value > room {
                    return Err(field);
                }
                out.extend_from_slice(&value.to_le_bytes()[..width as usize]);
            }
            Ok(())
        }
    }
}

/// One production record as an ng record. `None` for a record with no reference allele, which
/// has no span and which the record encoder refuses.
fn as_an_ng_record(record: &PileupRecord) -> Option<SampleLocusObservations> {
    let reference_bases: Vec<u8> = record.alleles.first()?.seq.clone();
    let span = u64::try_from(reference_bases.len()).ok()?;
    if span == 0 {
        return None;
    }
    let start = u64::from(record.pos);
    let observations = record
        .alleles
        .iter()
        .map(|allele| SequenceObservation {
            bases: allele.seq.clone(),
            // Production has no witness and no read group; every read is treated as having
            // spanned the locus, which is what makes the head's read count comparable.
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: allele.support.num_obs,
            num_fwd: allele.support.fwd,
            q_sum: SummedLogError::from_nats(allele.support.q_sum),
            mapq_sum: allele.support.mapq_sum,
            mapq_sum_sq: allele.support.mapq_sum_sq,
            placed_left: allele.support.placed_left,
            // Dropped by the encoder until Milestone E; carried so the day it stops being
            // dropped this harness measures the fuller record without being touched.
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
        reads_without_observation: 0,
        reads_discarded_by_cap: 0,
        kind: LocusKind::Generic,
    })
}

/// What one arm of the comparison produced.
#[derive(Default)]
struct Tally {
    records: u64,
    blocks: u64,
    head_bytes: u64,
    body_bytes: u64,
    payload_bytes: u64,
    compressed_bytes: u64,
}

/// Cut `records` on the grid, write each head the way `encoding` says, and compress every
/// block.
///
/// **The cut is reproduced here rather than driven through `BlockBuilder`**, because the
/// fixed-width arm needs a different head and the builder writes the format's own. Both arms go
/// through this one function, so the only difference between them is the head — and the varint
/// arm's payloads are checked against `BlockBuilder`'s, byte for byte, by `check_against_the_shipped_cut`.
fn run(
    records: &[SampleLocusObservations],
    encoding: HeadEncoding,
    grid: Bp,
    window_log: u8,
) -> Result<(Tally, Vec<Vec<u8>>), String> {
    let mut compressor = BlockCompressor::new(window_log).map_err(|e| e.to_string())?;
    let mut tally = Tally::default();
    let mut payloads: Vec<Vec<u8>> = Vec::new();

    let mut body = Vec::new();
    let mut heads_and_bodies: Vec<u8> = Vec::new();
    let mut open: Option<(ContigId, u64, Position, u64)> = None; // contig, cell, first, count
    let mut measured_from = Position(0);

    let close = |open: &mut Option<(ContigId, u64, Position, u64)>,
                 records_bytes: &mut Vec<u8>,
                 tally: &mut Tally,
                 payloads: &mut Vec<Vec<u8>>,
                 compressor: &mut BlockCompressor|
     -> Result<(), String> {
        let Some((contig, _, first, count)) = open.take() else {
            return Ok(());
        };
        let mut payload = Vec::new();
        BlockHead {
            contig,
            first_position: first,
            record_count: std::num::NonZeroU64::new(count).expect("a block holds a record"),
        }
        .encode(&mut payload);
        payload.extend_from_slice(records_bytes);
        records_bytes.clear();
        tally.payload_bytes += payload.len() as u64;
        tally.compressed_bytes += compressor
            .compress(&payload)
            .map_err(|e| e.to_string())?
            .len() as u64;
        tally.blocks += 1;
        payloads.push(payload);
        Ok(())
    };

    for record in records {
        let cell = record.region.start.get() / grid.get();
        let cuts = match open {
            Some((contig, open_cell, _, _)) => contig != record.region.contig || open_cell != cell,
            None => true,
        };
        if cuts {
            close(
                &mut open,
                &mut heads_and_bodies,
                &mut tally,
                &mut payloads,
                &mut compressor,
            )?;
            open = Some((record.region.contig, cell, record.region.start, 0));
            measured_from = record.region.start;
        }

        body.clear();
        encode_record_body(record, &mut body);
        let (non_reference_reads, _) = record.non_reference_and_compared_reads();
        let offset = record.region.start.get() - measured_from.get();
        let span = record.region.end.get() - record.region.start.get() + 1;

        let before = heads_and_bodies.len();
        put_head(
            &mut heads_and_bodies,
            encoding,
            offset,
            span,
            u64::from(non_reference_reads),
            body.len() as u64,
        )
        .map_err(|field| format!("{field} does not fit its declared width"))?;
        tally.head_bytes += (heads_and_bodies.len() - before) as u64;
        heads_and_bodies.extend_from_slice(&body);
        tally.body_bytes += body.len() as u64;

        measured_from = record.region.start;
        tally.records += 1;
        if let Some((_, _, _, count)) = open.as_mut() {
            *count += 1;
        }
    }
    close(
        &mut open,
        &mut heads_and_bodies,
        &mut tally,
        &mut payloads,
        &mut compressor,
    )?;
    Ok((tally, payloads))
}

/// **The check that keeps this harness honest.** The varint arm reproduces the shipped cut by
/// hand, so it has to produce the shipped bytes: every block payload byte-identical to what
/// `BlockBuilder` writes over the same records. Without it, a harness that quietly built
/// something else would report a difference between two encodings neither of which is ours.
fn check_against_the_shipped_cut(
    records: &[SampleLocusObservations],
    grid: Bp,
    ours: &[Vec<u8>],
) -> Result<(), String> {
    use pop_var_caller::ng::psp::BlockBuilder;
    let mut builder = BlockBuilder::new(grid, None).map_err(|e| e.to_string())?;
    let mut shipped = Vec::new();
    for record in records {
        if let Some(closed) = builder.push(record).map_err(|e| e.to_string())? {
            shipped.push(closed.to_vec());
        }
    }
    if let Some(last) = builder.finish() {
        shipped.push(last);
    }
    if shipped.len() != ours.len() {
        return Err(format!(
            "the harness cut {} blocks where the shipped builder cuts {}",
            ours.len(),
            shipped.len()
        ));
    }
    for (index, (a, b)) in ours.iter().zip(&shipped).enumerate() {
        if a != b {
            return Err(format!(
                "block {index}: the harness wrote {} bytes and the shipped builder {}",
                a.len(),
                b.len()
            ));
        }
    }
    Ok(())
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap_or_else(|| {
        eprintln!(
            "usage: ng_psp_head_encoding <production.psp> [--label NAME] [--grid BP] \
             [--window-log N] [--limit N]"
        );
        std::process::exit(2);
    });
    let mut label = "sample".to_string();
    let mut grid = DEFAULT_GENOMIC_BLOCK_SIZE_BP;
    let mut window_log = DEFAULT_LOOK_BACK_WINDOW_LOG;
    let mut limit = u64::MAX;
    while let Some(flag) = args.next() {
        let value = args.next().unwrap_or_else(|| {
            eprintln!("{flag} needs a value");
            std::process::exit(2);
        });
        match flag.as_str() {
            "--label" => label = value,
            "--grid" => grid = Bp(value.parse().expect("a number of base pairs")),
            "--window-log" => window_log = value.parse().expect("an exponent"),
            "--limit" => limit = value.parse().expect("a record count"),
            other => {
                eprintln!("unknown flag {other}");
                std::process::exit(2);
            }
        }
    }

    let read_at = std::time::Instant::now();
    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(&path).expect("the production psp opens"),
    ))
    .expect("it is a production psp");
    let mut records = Vec::new();
    let mut skipped = 0u64;
    {
        let iter = reader.records();
        for record in iter {
            let record = record.expect("a production record");
            match as_an_ng_record(&record) {
                Some(ours) => records.push(ours),
                None => skipped += 1,
            }
            if records.len() as u64 >= limit {
                break;
            }
        }
    }
    // **The corpus describes itself, rather than being described from memory.** The depth
    // label on a sample is exactly the kind of number this project keeps getting wrong: the
    // specs call this tomato accession "three reads a position" and it measures ten.
    let reads: u64 = records
        .iter()
        .map(|record| {
            record
                .observations
                .iter()
                .map(|observation| u64::from(observation.num_obs))
                .sum::<u64>()
        })
        .sum();
    // Which records could not be written with a two-byte position offset, and how far past it
    // the worst one is — the fact behind a narrow fixed-width arm that cannot encode a sample.
    let mut over_a_two_byte_offset = 0u64;
    let mut widest_offset = 0u64;
    let mut measured_from: Option<(ContigId, u64, Position)> = None;
    for record in &records {
        let cell = record.region.start.get() / grid.get();
        let offset = match measured_from {
            Some((contig, open_cell, previous))
                if contig == record.region.contig && open_cell == cell =>
            {
                record.region.start.get() - previous.get()
            }
            _ => 0,
        };
        if offset > u64::from(u16::MAX) {
            over_a_two_byte_offset += 1;
        }
        widest_offset = widest_offset.max(offset);
        measured_from = Some((record.region.contig, cell, record.region.start));
    }

    println!("sample\t{label}");
    println!("psp\t{path}");
    println!("records\t{}", records.len());
    println!("records-skipped-no-reference-allele\t{skipped}");
    println!(
        "mean-reads-a-position\t{:.2}",
        reads as f64 / records.len() as f64
    );
    println!("genomic-block-size-bp\t{}", grid.get());
    println!("look-back-window-bytes\t{}", 1u64 << window_log);
    println!("widest-within-block-position-offset\t{widest_offset}");
    println!("records-whose-offset-passes-65535\t{over_a_two_byte_offset}");
    println!("read-seconds\t{:.1}", read_at.elapsed().as_secs_f64());

    // Every field at the width the *format* allows: a position offset is bounded by the grid,
    // a span and a read count and a body length are each what a `u32` holds.
    let widths = [
        HeadEncoding::Varint,
        HeadEncoding::Fixed {
            offset: 4,
            span: 4,
            reads: 4,
            body: 4,
        },
        // And the narrow widths this corpus happens to fit in, which a manifest could declare
        // only for a writer that knew its data in advance.
        HeadEncoding::Fixed {
            offset: 4,
            span: 2,
            reads: 2,
            body: 4,
        },
        HeadEncoding::Fixed {
            offset: 2,
            span: 1,
            reads: 2,
            body: 2,
        },
    ];

    println!(
        "\nencoding\tblocks\thead-bytes-a-record\tuncompressed-bytes-a-record\tcompressed-bytes-a-record\tagainst-varint"
    );
    let mut varint_compressed = 0f64;
    for encoding in widths {
        match run(&records, encoding, grid, window_log) {
            Err(why) => println!("{}\t—\t—\t—\t—\t{why}", encoding.name()),
            Ok((tally, payloads)) => {
                if matches!(encoding, HeadEncoding::Varint) {
                    check_against_the_shipped_cut(&records, grid, &payloads)
                        .expect("the harness writes the bytes the shipped builder writes");
                }
                let per = |n: u64| n as f64 / tally.records as f64;
                let compressed = per(tally.compressed_bytes);
                if matches!(encoding, HeadEncoding::Varint) {
                    varint_compressed = compressed;
                }
                println!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:+.2}%",
                    encoding.name(),
                    tally.blocks,
                    per(tally.head_bytes),
                    per(tally.payload_bytes),
                    compressed,
                    100.0 * (compressed - varint_compressed) / varint_compressed
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Read the four head fields back out of a fixed-width head, so the arm that produces the
    /// comparison's fixed rows has an oracle of its own.
    ///
    /// **The varint arm is checked against the shipped `BlockBuilder` and the fixed arm was
    /// checked against nothing** — and the fixed arm is the one whose bytes the conclusion
    /// turns on. A field order or a width silently wrong there would report a difference
    /// between varint and *something nobody wrote*.
    fn take_fixed_head(bytes: &[u8], widths: [u8; 4]) -> ([u64; 4], usize) {
        let mut at = 0usize;
        let mut values = [0u64; 4];
        for (index, width) in widths.into_iter().enumerate() {
            let mut whole = [0u8; 8];
            whole[..width as usize].copy_from_slice(&bytes[at..at + width as usize]);
            values[index] = u64::from_le_bytes(whole);
            at += width as usize;
        }
        (values, at)
    }

    #[test]
    fn a_fixed_width_head_reads_back_field_for_field() {
        let widths = [4u8, 2, 2, 4];
        let encoding = HeadEncoding::Fixed {
            offset: widths[0],
            span: widths[1],
            reads: widths[2],
            body: widths[3],
        };
        // Values that need every byte of their declared width, so a width read one short or
        // one long is visible rather than absorbed by leading zeros.
        let written = [3_000_000_007u64, 40_001, 60_002, 4_000_000_009];

        let mut bytes = Vec::new();
        put_head(
            &mut bytes, encoding, written[0], written[1], written[2], written[3],
        )
        .expect("every value fits its declared width");

        let (read, used) = take_fixed_head(&bytes, widths);
        assert_eq!(read, written, "field for field");
        assert_eq!(used, bytes.len(), "and nothing else was written");
        assert_eq!(used, 12, "four plus two plus two plus four");
    }

    /// **The width bound is exact.** A value of exactly what a width holds is written; one more
    /// is refused and names its own field. An off-by-one there would silently truncate a
    /// position offset — the field that decides whether the narrow arm can encode a sample at
    /// all.
    #[test]
    fn a_value_one_past_its_declared_width_is_refused_and_names_its_field() {
        for (width, largest) in [(1u8, 255u64), (2, 65_535), (4, 4_294_967_295)] {
            let encoding = HeadEncoding::Fixed {
                offset: width,
                span: 8,
                reads: 8,
                body: 8,
            };
            let mut bytes = Vec::new();
            put_head(&mut bytes, encoding, largest, 0, 0, 0).unwrap_or_else(|field| {
                panic!("{largest} must fit {width} bytes; {field} refused")
            });
            assert_eq!(bytes.len(), width as usize + 24);

            let mut bytes = Vec::new();
            assert_eq!(
                put_head(&mut bytes, encoding, largest + 1, 0, 0, 0),
                Err("position-offset"),
                "{} must not fit {width} bytes",
                largest + 1
            );
        }
    }

    /// An eight-byte width holds every value there is, so the bound must not overflow computing
    /// its own ceiling.
    #[test]
    fn an_eight_byte_width_holds_every_value_there_is() {
        let encoding = HeadEncoding::Fixed {
            offset: 8,
            span: 8,
            reads: 8,
            body: 8,
        };
        let mut bytes = Vec::new();
        put_head(&mut bytes, encoding, u64::MAX, u64::MAX, u64::MAX, u64::MAX)
            .expect("eight bytes hold a u64");
        assert_eq!(bytes.len(), 32);
        assert_eq!(take_fixed_head(&bytes, [8, 8, 8, 8]).0, [u64::MAX; 4]);
    }

    /// The varint arm writes LEB128, checked against an independent decode rather than against
    /// itself. Its agreement with the shipped format is checked at run time, against
    /// `BlockBuilder`'s own bytes.
    #[test]
    fn the_varint_arm_writes_leb128() {
        let mut bytes = Vec::new();
        put_head(&mut bytes, HeadEncoding::Varint, 0, 1, 127, 128).expect("varints never refuse");
        assert_eq!(bytes, vec![0x00, 0x01, 0x7f, 0x80, 0x01]);
    }
}
