//! What does naming every read at every position cost on disk, and does storing
//! only the changes make it cheaper?
//!
//! ng's per-sample store names, at every covered reference position, every read
//! folded there — the owner's ruling of 2026-08-17, because a read that covered a
//! position and agreed with the reference has to be distinguishable from one that
//! never reached it. Production stores about 3.4 % of those names, so no
//! measurement taken on a `.psp` says anything about what ng's column will cost.
//!
//! This program takes that measurement from the alignments themselves. It walks a
//! sample's reads, gives each read — or each read *pair*, mates collapsed — one
//! identifier allocated in order, and works out which identifiers are live at each
//! reference position. Then it writes that same information three ways and
//! compresses each with the same settings:
//!
//! - **whole list, raw** — a count then one 8-byte identifier per live read at
//!   every position, which is what production's column format does;
//! - **whole list, deltas** — the same list, but each identifier stored as its
//!   distance from the one before, as a variable-length integer;
//! - **changes only** — per position, which identifiers started covering it and
//!   which stopped, the arrivals as deltas and the departures as their place in
//!   the live set. Every frame restates the whole live set first, so a reader can
//!   start at any frame without reading the ones before it.
//!
//! It also counts the thing that decides whether the third form needs a re-entry
//! rule: how often one identifier goes live, stops, and goes live again, which
//! happens when a pair's two mates have a gap between them.
//!
//! ```text
//! ng_chain_id_column_cost <reference.fa> <sample.cram|bam> [contig] [--frame-kib K]
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;

use pop_var_caller::bam::alignment_input::build_fasta_repository;

// ---------------------------------------------------------------------
// varint
// ---------------------------------------------------------------------

fn put_varint(buf: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            buf.push(byte);
            return;
        }
        buf.push(byte | 0x80);
    }
}

// ---------------------------------------------------------------------
// One encoding of the column, cut into independent frames
// ---------------------------------------------------------------------

struct Arm {
    name: &'static str,
    raw: Vec<u8>,
    frame_bytes: usize,
    raw_total: u64,
    compressed_total: u64,
    /// Set when the arm has to restate its whole state at a frame boundary.
    restates_live_set: bool,
    samples: Vec<Vec<u8>>,
}

impl Arm {
    fn new(name: &'static str, frame_bytes: usize, restates_live_set: bool) -> Self {
        Self {
            name,
            raw: Vec::with_capacity(frame_bytes * 2),
            frame_bytes,
            raw_total: 0,
            compressed_total: 0,
            restates_live_set,
            samples: Vec::new(),
        }
    }

    /// True when this position filled the frame, so the caller knows to restate
    /// the live set at the start of the next one.
    fn end_position(&mut self, comp: &mut zstd::bulk::Compressor<'static>, force: bool) -> bool {
        if !force && self.raw.len() < self.frame_bytes {
            return false;
        }
        if self.raw.is_empty() {
            return false;
        }
        self.flush(comp);
        true
    }

    fn flush(&mut self, comp: &mut zstd::bulk::Compressor<'static>) {
        if self.raw.is_empty() {
            return;
        }
        self.raw_total += self.raw.len() as u64;
        self.compressed_total += comp.compress(&self.raw).expect("compress").len() as u64;
        if self.samples.len() < 128 {
            self.samples.push(self.raw.clone());
        }
        self.raw.clear();
    }
}

// ---------------------------------------------------------------------
// Reading the alignments into (start, end, identifier) intervals
// ---------------------------------------------------------------------

struct Interval {
    start: u64,
    end: u64,
    id: u64,
}

/// The reference stretch a read covers, as the aligned blocks of its CIGAR.
/// Deletions count as covered — the read is evidence there. Skips (`N`) do not.
fn covered_spans(
    start: u64,
    cigar: impl Iterator<Item = (Kind, usize)>,
    out: &mut Vec<(u64, u64)>,
) {
    out.clear();
    let mut pos = start;
    let mut open: Option<(u64, u64)> = None;
    for (kind, len) in cigar {
        let len = len as u64;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                match &mut open {
                    Some((_, end)) => *end = pos + len,
                    None => open = Some((pos, pos + len)),
                }
                pos += len;
            }
            Kind::Skip => {
                if let Some(span) = open.take() {
                    out.push(span);
                }
                pos += len;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if let Some(span) = open.take() {
        out.push(span);
    }
}

struct Collected {
    intervals: Vec<Interval>,
    n_reads: u64,
    n_ids: u64,
    /// Identifiers whose coverage is more than one stretch — a pair whose mates
    /// have a gap between them.
    n_ids_with_a_gap: u64,
}

fn collect(reference: &Path, alignment: &Path, contig: Option<&str>) -> Collected {
    let is_cram = alignment.extension().and_then(|e| e.to_str()) == Some("cram");
    let mut intervals: Vec<Interval> = Vec::new();
    let mut by_name: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut spans_per_id: HashMap<u64, u32> = HashMap::new();
    let mut next_id = 0u64;
    let mut n_reads = 0u64;
    let mut spans = Vec::new();

    let mut take = |name: Option<&[u8]>,
                    flags: u16,
                    ref_name: Option<&str>,
                    start: u64,
                    spans: &[(u64, u64)],
                    intervals: &mut Vec<Interval>| {
        let _ = start;
        // Unmapped, secondary, supplementary and duplicate reads are not folded.
        if flags & 0x4 != 0 || flags & 0x100 != 0 || flags & 0x800 != 0 || flags & 0x400 != 0 {
            return;
        }
        if let (Some(want), Some(got)) = (contig, ref_name)
            && want != got
        {
            return;
        }
        let key = name.unwrap_or(b"").to_vec();
        let id = *by_name.entry(key).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
        n_reads += 1;
        for &(s, e) in spans {
            intervals.push(Interval {
                start: s,
                end: e,
                id,
            });
            *spans_per_id.entry(id).or_insert(0) += 1;
        }
    };

    if is_cram {
        let repository = build_fasta_repository(reference).expect("reference repository");
        let mut reader = open_cram(alignment, &repository);
        let header = read_cram_header(&mut reader);
        let names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
            .collect();
        for result in reader.records(&header) {
            let record = result.expect("cram record");
            let Some(start) = record.alignment_start() else {
                continue;
            };
            let ref_name = record
                .reference_sequence_id()
                .and_then(|i| names.get(i).map(String::as_str));
            let start = usize::from(start) as u64;
            covered_spans(
                start,
                record
                    .cigar()
                    .as_ref()
                    .iter()
                    .map(|op| (op.kind(), op.len())),
                &mut spans,
            );
            let flags = u16::from(record.flags());
            take(
                record.name().map(|n| n.as_ref()),
                flags,
                ref_name,
                start,
                &spans,
                &mut intervals,
            );
        }
    } else {
        let mut reader = bam::io::reader::Builder
            .build_from_path(alignment)
            .expect("open bam");
        let header = reader.read_header().expect("bam header");
        let names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
            .collect();
        let mut record = sam::alignment::RecordBuf::default();
        while reader
            .read_record_buf(&header, &mut record)
            .expect("bam record")
            != 0
        {
            let Some(start) = record.alignment_start() else {
                continue;
            };
            let ref_name = record
                .reference_sequence_id()
                .and_then(|i| names.get(i).map(String::as_str));
            let start = usize::from(start) as u64;
            covered_spans(
                start,
                record
                    .cigar()
                    .as_ref()
                    .iter()
                    .map(|op| (op.kind(), op.len())),
                &mut spans,
            );
            let flags = u16::from(record.flags());
            take(
                record.name().map(|n| n.as_ref()),
                flags,
                ref_name,
                start,
                &spans,
                &mut intervals,
            );
        }
    }

    intervals.sort_by_key(|i| (i.start, i.id));
    // An identifier's coverage is discontinuous only when two of its stretches
    // do not touch — mates that overlap are one stretch, not two.
    let mut per_id: HashMap<u64, Vec<(u64, u64)>> = HashMap::new();
    for iv in &intervals {
        per_id.entry(iv.id).or_default().push((iv.start, iv.end));
    }
    let mut n_ids_with_a_gap = 0u64;
    for spans in per_id.values_mut() {
        spans.sort_unstable();
        let mut runs = 1;
        let mut end = spans[0].1;
        for &(s, e) in &spans[1..] {
            if s > end {
                runs += 1;
            }
            end = end.max(e);
        }
        if runs > 1 {
            n_ids_with_a_gap += 1;
        }
    }
    let _ = &spans_per_id;
    Collected {
        intervals,
        n_reads,
        n_ids: next_id,
        n_ids_with_a_gap,
    }
}

fn open_cram(path: &Path, repository: &fasta::Repository) -> cram::io::Reader<File> {
    cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .expect("open cram")
}

fn read_cram_header(reader: &mut cram::io::Reader<File>) -> sam::Header {
    reader.read_file_definition().expect("cram file definition");
    reader.read_file_header().expect("cram header")
}

// ---------------------------------------------------------------------
// The sweep
// ---------------------------------------------------------------------

fn main() {
    let mut args = std::env::args().skip(1);
    let reference = PathBuf::from(
        args.next()
            .expect("usage: <reference.fa> <sample.cram|bam> [contig] [--frame-kib K] [--level L]"),
    );
    let alignment = PathBuf::from(args.next().expect("alignment path"));
    let mut contig: Option<String> = None;
    let mut frame_kib = 32usize;
    let mut level = 9i32;
    // When set, every arm's frame is cut after this many covered positions —
    // the shape a shared file has, where the record stream decides the cut and
    // the chain-id stream must restate whatever it carries at the same points.
    let mut frame_positions: Option<u64> = None;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--frame-kib" => frame_kib = args.next().expect("K").parse().expect("K"),
            "--level" => level = args.next().expect("L").parse().expect("L"),
            "--frame-positions" => {
                frame_positions = Some(args.next().expect("N").parse().expect("N"))
            }
            other => contig = Some(other.to_string()),
        }
    }

    eprintln!("reading alignments…");
    let collected = collect(&reference, &alignment, contig.as_deref());
    assert!(!collected.intervals.is_empty(), "no reads collected");

    let frame_bytes = frame_kib * 1024;
    let mut comp = zstd::bulk::Compressor::new(level).expect("compressor");
    comp.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
        .expect("content size");

    let mut whole_raw = Arm::new("whole list, raw 8-byte ids", frame_bytes, false);
    let mut whole_delta = Arm::new("whole list, deltas", frame_bytes, false);
    let mut changes = Arm::new("changes only", frame_bytes, true);

    // The sweep: positions in order, a sorted live set, and the two events.
    let mut live: Vec<u64> = Vec::new();
    // Ends of the live intervals, parallel to `live` by identifier.
    let mut ends: HashMap<u64, u64> = HashMap::new();
    let mut next_interval = 0usize;
    let mut arrivals: Vec<u64> = Vec::new();
    let mut departures: Vec<u32> = Vec::new();

    let first_pos = collected.intervals[0].start;
    let last_pos = collected
        .intervals
        .iter()
        .map(|i| i.end)
        .max()
        .expect("end");
    let mut n_positions = 0u64;
    let mut n_id_mentions = 0u64;
    let mut restate_positions = 0u64;
    let mut changes_needs_restate = true;

    eprintln!("sweeping {} positions…", last_pos - first_pos);
    for pos in first_pos..last_pos {
        arrivals.clear();
        while next_interval < collected.intervals.len()
            && collected.intervals[next_interval].start == pos
        {
            let iv = &collected.intervals[next_interval];
            arrivals.push(iv.id);
            // A read can be named twice at one position only if two intervals of
            // the same identifier meet, which the sort makes adjacent.
            let end = ends.entry(iv.id).or_insert(0);
            *end = (*end).max(iv.end);
            next_interval += 1;
        }
        departures.clear();
        let mut i = 0;
        while i < live.len() {
            if ends[&live[i]] <= pos {
                departures.push(i as u32);
                live.remove(i);
            } else {
                i += 1;
            }
        }
        for &id in &arrivals {
            match live.binary_search(&id) {
                Ok(_) => {}
                Err(at) => live.insert(at, id),
            }
        }
        if live.is_empty() {
            continue;
        }
        n_positions += 1;
        n_id_mentions += live.len() as u64;

        // Whole list, raw — production's column shape.
        put_varint(&mut whole_raw.raw, live.len() as u64);
        for &id in &live {
            whole_raw.raw.extend_from_slice(&id.to_le_bytes());
        }

        // Whole list, deltas.
        put_varint(&mut whole_delta.raw, live.len() as u64);
        let mut prev = 0u64;
        for &id in &live {
            put_varint(&mut whole_delta.raw, id - prev);
            prev = id;
        }

        // Changes only, restating the live set whenever a frame has just closed.
        if changes_needs_restate {
            restate_positions += 1;
            put_varint(&mut changes.raw, live.len() as u64);
            let mut prev = 0u64;
            for &id in &live {
                put_varint(&mut changes.raw, id - prev);
                prev = id;
            }
            changes_needs_restate = false;
        } else {
            put_varint(&mut changes.raw, departures.len() as u64);
            for &at in &departures {
                put_varint(&mut changes.raw, at as u64);
            }
            put_varint(&mut changes.raw, arrivals.len() as u64);
            let mut prev = 0u64;
            for &id in &arrivals {
                put_varint(&mut changes.raw, id - prev);
                prev = id;
            }
        }

        let force = frame_positions.is_some_and(|n| n_positions.is_multiple_of(n));
        whole_raw.end_position(&mut comp, force);
        whole_delta.end_position(&mut comp, force);
        if changes.end_position(&mut comp, force) {
            changes_needs_restate = true;
        }
    }
    whole_raw.flush(&mut comp);
    whole_delta.flush(&mut comp);
    changes.flush(&mut comp);

    println!("# alignment\t{}", alignment.display());
    if let Some(c) = &contig {
        println!("# contig\t{c}");
    }
    println!("# reads-folded\t{}", collected.n_reads);
    println!("# read-identifiers\t{}", collected.n_ids);
    println!(
        "# identifiers-covering-more-than-one-stretch\t{}\t{:.4}",
        collected.n_ids_with_a_gap,
        collected.n_ids_with_a_gap as f64 / collected.n_ids as f64
    );
    println!("# covered-positions\t{n_positions}");
    println!("# name-mentions\t{n_id_mentions}");
    println!(
        "# mean-live-reads-per-position\t{:.2}",
        n_id_mentions as f64 / n_positions as f64
    );
    println!("# frame-kib\t{frame_kib}");
    match frame_positions {
        Some(n) => println!("# frame-cut-every-positions\t{n}"),
        None => println!("# frame-cut-every-positions\tby-bytes-only"),
    }
    println!("# zstd-level\t{level}");
    println!(
        "# live-set-restatements\t{restate_positions}\t{:.4} of positions",
        restate_positions as f64 / n_positions as f64
    );
    println!("encoding\traw_bytes\tcompressed_bytes\tbytes_per_position\tbytes_per_name");
    for arm in [&whole_raw, &whole_delta, &changes] {
        println!(
            "{}\t{}\t{}\t{:.3}\t{:.4}",
            arm.name,
            arm.raw_total,
            arm.compressed_total,
            arm.compressed_total as f64 / n_positions as f64,
            arm.compressed_total as f64 / n_id_mentions as f64
        );
    }

    // The same three, primed with a dictionary. **Trained on the odd-numbered
    // frames and measured on the even ones** — a dictionary trained on the very
    // bytes it is then scored against would report a saving no reader ever gets.
    println!("encoding+dictionary\traw_bytes\tcompressed_bytes\tbytes_per_position_scaled");
    for arm in [&whole_raw, &whole_delta, &changes] {
        if arm.samples.len() < 16 {
            continue;
        }
        let train: Vec<&Vec<u8>> = arm.samples.iter().skip(1).step_by(2).collect();
        let test: Vec<&Vec<u8>> = arm.samples.iter().step_by(2).collect();
        let Ok(dict) = zstd::dict::from_samples(&train, 112 * 1024) else {
            continue;
        };
        let mut c = zstd::bulk::Compressor::with_dictionary(level, &dict).expect("dict compressor");
        c.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
            .expect("content size");
        let raw: u64 = test.iter().map(|s| s.len() as u64).sum();
        let total: u64 = test
            .iter()
            .map(|s| c.compress(s).expect("compress").len() as u64)
            .sum();
        // Scaled onto the same per-position axis by the arm's own raw-byte share.
        let per_position =
            (total as f64 / raw as f64) * (arm.raw_total as f64 / n_positions as f64);
        println!("{}\t{}\t{}\t{:.3}", arm.name, raw, total, per_position);
        let _ = arm.restates_live_set;
    }
}
