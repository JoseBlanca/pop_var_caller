//! What does a per-sample record store cost, in bytes on disk, for a given
//! amount of memory a reader must hold?
//!
//! The production `.psp` answers that question with one knob. It groups records
//! into blocks of a target uncompressed size, transposes each block into one
//! buffer per field, and compresses every field buffer as its own zstd frame.
//! The block is both the unit a reader must decode before it can hand out the
//! first record **and** the furthest back the compressor may look for a repeat.
//! Shrinking it to save memory therefore also costs ratio.
//!
//! Those two roles can be separated. zstd can compress an unbounded stream with
//! a sliding window of a chosen size, so the compressor's reach (the window) and
//! the reader's materialisation unit (a batch of records, or a single record)
//! become independent knobs. This probe measures the whole grid:
//!
//! - **field order** — records laid end to end (`row`), or transposed into one
//!   buffer per field (`col`);
//! - **field width** — fixed-width little-endian integers as the production
//!   columns store them (`fixed`), or LEB128 varints with deltas (`varint`);
//! - **batch** — how many bytes of records go into one materialisation unit;
//! - **framing** — an independent zstd frame per batch (what production does,
//!   with or without a trained dictionary), or one continuous stream whose
//!   window is set independently of the batch, with and without a flush at each
//!   batch boundary. A flush is what makes a stream resumable: it ends a zstd
//!   block so a reader can start there. Its cost is measured, not assumed.
//!
//! Each row of the output also carries `reader_kib`: the bytes a reader has to
//! hold to yield records, counted as the zstd window plus whatever the layout
//! forces it to materialise at once (a whole batch of decoded columns for `col`,
//! one record for `row`). That is the axis to read the table against — a layout
//! is only better if it is cheaper at the same `reader_kib`.
//!
//! The `col`/`frames` configurations also report each field's own compressed
//! bytes, which says where the file actually is before any layout is chosen.
//!
//! Run:
//!
//! ```text
//! cargo run --release --example psp_record_stream_compression -- \
//!     tmp/em_bench/giab/psp/HG002.30x.psp --max-records 3000000
//! ```

use std::fs::File;
use std::io::{BufReader, Write};
use std::time::Instant;

use pop_var_caller::pileup_record::PileupRecord;
use pop_var_caller::psp::PspReader;

// ---------------------------------------------------------------------
// Field identity — one slot per production column, so the transposed layout
// writes exactly the buffers the `.psp` writes today.
// ---------------------------------------------------------------------

const N_SLOTS: usize = 14;

const SLOT_NAMES: [&str; N_SLOTS] = [
    "delta-pos",
    "n-alleles",
    "allele-seq-len",
    "allele-seq",
    "windowed-gc",
    "windowed-coverage",
    "allele-obs-count",
    "allele-q-sum-log",
    "allele-fwd-count",
    "allele-placed-left",
    "allele-placed-start",
    "allele-mapq-sum",
    "allele-mapq-sum-sq",
    "allele-chain-ids",
];

const S_DELTA_POS: usize = 0;
const S_N_ALLELES: usize = 1;
const S_SEQ_LEN: usize = 2;
const S_SEQ: usize = 3;
const S_GC: usize = 4;
const S_COV: usize = 5;
const S_OBS: usize = 6;
const S_QSUM: usize = 7;
const S_FWD: usize = 8;
const S_PLACED_LEFT: usize = 9;
const S_PLACED_START: usize = 10;
const S_MAPQ_SUM: usize = 11;
const S_MAPQ_SUM_SQ: usize = 12;
const S_CHAIN: usize = 13;

// ---------------------------------------------------------------------
// Primitive encodings
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

fn put_svarint(buf: &mut Vec<u8>, v: i64) {
    put_varint(buf, ((v << 1) ^ (v >> 63)) as u64);
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Width {
    Fixed,
    Varint,
    /// Varints, plus the three floating-point fields replaced by fixed-point
    /// integers: the window's GC fraction to 1 part in 10,000, its mean coverage
    /// to 1/16 of a read and stored as a difference from the position before,
    /// and each allele's summed log-error to 1/256 of a natural-log unit. All
    /// three are quantities whose low mantissa bits are noise, and noise is what
    /// a compressor cannot shrink.
    Quant,
}

impl Width {
    fn name(self) -> &'static str {
        match self {
            Width::Fixed => "fixed",
            Width::Varint => "varint",
            Width::Quant => "quant",
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Order {
    Row,
    Col,
}

impl Order {
    fn name(self) -> &'static str {
        match self {
            Order::Row => "row",
            Order::Col => "col",
        }
    }
}

// ---------------------------------------------------------------------
// Batch builder: identical bytes, routed either into one buffer or into one
// buffer per field.
// ---------------------------------------------------------------------

struct BatchBuilder {
    order: Order,
    width: Width,
    row: Vec<u8>,
    slots: Vec<Vec<u8>>,
    prev_pos: u32,
    prev_cov_q: i64,
    qsum_scale: f64,
    qsum_predict: bool,
    qsum_mean_q: i64,
    last_chain: u64,
    n_records: u64,
    slot_totals: [u64; N_SLOTS],
}

impl BatchBuilder {
    fn new(order: Order, width: Width) -> Self {
        Self {
            order,
            width,
            row: Vec::with_capacity(1 << 20),
            slots: (0..N_SLOTS).map(|_| Vec::with_capacity(1 << 16)).collect(),
            prev_pos: 0,
            prev_cov_q: 0,
            qsum_scale: std::env::var("QSUM_SCALE")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(256.0),
            qsum_predict: std::env::var("QSUM_PREDICT").is_ok(),
            qsum_mean_q: 0,
            last_chain: 0,
            n_records: 0,
            slot_totals: [0; N_SLOTS],
        }
    }

    #[inline]
    fn buf(&mut self, slot: usize) -> &mut Vec<u8> {
        match self.order {
            Order::Row => &mut self.row,
            Order::Col => &mut self.slots[slot],
        }
    }

    #[inline]
    fn slot_len(&self, slot: usize) -> usize {
        match self.order {
            Order::Row => self.row.len(),
            Order::Col => self.slots[slot].len(),
        }
    }

    fn size(&self) -> usize {
        match self.order {
            Order::Row => self.row.len(),
            Order::Col => self.slots.iter().map(Vec::len).sum(),
        }
    }

    #[inline]
    fn put_u32(&mut self, slot: usize, v: u32) {
        let width = self.width;
        let before = self.slot_len(slot);
        let buf = self.buf(slot);
        match width {
            Width::Fixed => buf.extend_from_slice(&v.to_le_bytes()),
            Width::Varint | Width::Quant => put_varint(buf, v as u64),
        }
        self.slot_totals[slot] += (self.slot_len(slot) - before) as u64;
    }

    fn push(&mut self, rec: &PileupRecord) {
        let delta = if self.n_records == 0 {
            0
        } else {
            rec.pos.saturating_sub(self.prev_pos) as u64
        };
        {
            // `delta-pos` is a varint in the production column too, so both
            // widths agree here.
            let before = self.slot_len(S_DELTA_POS);
            let buf = self.buf(S_DELTA_POS);
            put_varint(buf, delta);
            self.slot_totals[S_DELTA_POS] += (self.slot_len(S_DELTA_POS) - before) as u64;
        }
        self.prev_pos = rec.pos;

        {
            let before = self.slot_len(S_N_ALLELES);
            let buf = self.buf(S_N_ALLELES);
            put_varint(buf, rec.alleles.len() as u64);
            self.slot_totals[S_N_ALLELES] += (self.slot_len(S_N_ALLELES) - before) as u64;
        }
        if self.width == Width::Quant {
            // A missing window is a real value here (an `N` reference position
            // has none), so 0 is reserved for it and every present value is
            // shifted by one.
            let before = self.slot_len(S_GC);
            let code = if rec.windowed_gc.is_nan() {
                0
            } else {
                1 + (rec.windowed_gc.clamp(0.0, 1.0) * 10_000.0).round() as u64
            };
            let buf = self.buf(S_GC);
            put_varint(buf, code);
            self.slot_totals[S_GC] += (self.slot_len(S_GC) - before) as u64;

            let before = self.slot_len(S_COV);
            let (code, next_prev) = if rec.windowed_coverage.is_nan() {
                (0u64, self.prev_cov_q)
            } else {
                let q = (rec.windowed_coverage.max(0.0) * 16.0).round() as i64;
                let d = q - self.prev_cov_q;
                (1 + (((d << 1) ^ (d >> 63)) as u64), q)
            };
            let buf = self.buf(S_COV);
            put_varint(buf, code);
            self.prev_cov_q = next_prev;
            self.slot_totals[S_COV] += (self.slot_len(S_COV) - before) as u64;
        } else {
            let gc = rec.windowed_gc.to_le_bytes();
            self.buf(S_GC).extend_from_slice(&gc);
            self.slot_totals[S_GC] += 4;
            let cov = rec.windowed_coverage.to_le_bytes();
            self.buf(S_COV).extend_from_slice(&cov);
            self.slot_totals[S_COV] += 4;
        }

        for allele in &rec.alleles {
            {
                let before = self.slot_len(S_SEQ_LEN);
                let buf = self.buf(S_SEQ_LEN);
                put_varint(buf, allele.seq.len() as u64);
                self.slot_totals[S_SEQ_LEN] += (self.slot_len(S_SEQ_LEN) - before) as u64;
            }
            self.buf(S_SEQ).extend_from_slice(&allele.seq);
            self.slot_totals[S_SEQ] += allele.seq.len() as u64;

            let s = &allele.support;
            self.put_u32(S_OBS, s.num_obs);
            if self.width == Width::Quant {
                let before = self.slot_len(S_QSUM);
                let q = (s.q_sum * self.qsum_scale).round() as i64;
                // The summed log-error is close to the read count times a
                // per-read average that barely moves within a sample, so what is
                // stored is the difference from that prediction.
                let v = if self.qsum_predict {
                    q - (s.num_obs as i64) * self.qsum_mean_q
                } else {
                    q
                };
                if self.qsum_predict && s.num_obs > 0 {
                    let per_read = q / (s.num_obs as i64);
                    self.qsum_mean_q = (self.qsum_mean_q * 15 + per_read) / 16;
                }
                let buf = self.buf(S_QSUM);
                put_svarint(buf, v);
                self.slot_totals[S_QSUM] += (self.slot_len(S_QSUM) - before) as u64;
            } else {
                let bits = s.q_sum.to_le_bytes();
                self.buf(S_QSUM).extend_from_slice(&bits);
                self.slot_totals[S_QSUM] += 8;
            }
            self.put_u32(S_FWD, s.fwd);
            self.put_u32(S_PLACED_LEFT, s.placed_left);
            self.put_u32(S_PLACED_START, s.placed_start);
            self.put_u32(S_MAPQ_SUM, s.mapq_sum);
            {
                let width = self.width;
                let before = self.slot_len(S_MAPQ_SUM_SQ);
                let buf = self.buf(S_MAPQ_SUM_SQ);
                match width {
                    Width::Fixed => buf.extend_from_slice(&s.mapq_sum_sq.to_le_bytes()),
                    Width::Varint | Width::Quant => put_varint(buf, s.mapq_sum_sq),
                }
                self.slot_totals[S_MAPQ_SUM_SQ] += (self.slot_len(S_MAPQ_SUM_SQ) - before) as u64;
            }
            {
                let width = self.width;
                let mut last = self.last_chain;
                let before = self.slot_len(S_CHAIN);
                {
                    let buf = self.buf(S_CHAIN);
                    put_varint(buf, allele.chain_ids.len() as u64);
                    for &id in &allele.chain_ids {
                        match width {
                            Width::Fixed => buf.extend_from_slice(&id.to_le_bytes()),
                            Width::Varint | Width::Quant => {
                                put_svarint(buf, (id as i64).wrapping_sub(last as i64));
                                last = id;
                            }
                        }
                    }
                }
                self.last_chain = last;
                self.slot_totals[S_CHAIN] += (self.slot_len(S_CHAIN) - before) as u64;
            }
        }
        self.n_records += 1;
    }

    /// Hand out the batch and reset the delta bases: every batch has to stand on
    /// its own, whether the framing is an independent frame or a flush point in
    /// a stream.
    fn take(&mut self, out: &mut Vec<Vec<u8>>) {
        out.clear();
        match self.order {
            Order::Row => {
                out.push(std::mem::take(&mut self.row));
                self.row = Vec::with_capacity(1 << 20);
            }
            Order::Col => {
                for s in &mut self.slots {
                    out.push(std::mem::take(s));
                    *s = Vec::with_capacity(1 << 16);
                }
            }
        }
        self.prev_pos = 0;
        self.prev_cov_q = 0;
        self.last_chain = 0;
        self.n_records = 0;
    }
}

// ---------------------------------------------------------------------
// Framings
// ---------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
enum Framing {
    /// One zstd frame per batch. `per_field` splits it into one frame per field,
    /// which is what production writes; `dict` primes each frame with a
    /// dictionary trained on this layout's own batches.
    Frames { per_field: bool, dict: bool },
    /// One zstd stream over the whole file, window capped at `window`.
    Stream { window: usize, flush: bool },
}

/// A sink that counts bytes and discards them: the probe measures size, and
/// writing a hundred candidate files would measure the disk.
struct CountingSink(u64);

impl Write for CountingSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0 += buf.len() as u64;
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

struct Engine {
    order: Order,
    width: Width,
    batch: usize,
    framing: Framing,
    compressor: Option<zstd::bulk::Compressor<'static>>,
    stream: Option<zstd::stream::write::Encoder<'static, CountingSink>>,
    out_bytes: u64,
    /// Compressed bytes attributable to each field — only filled by the
    /// per-field framing, where each field is its own frame.
    slot_bytes: [u64; N_SLOTS],
    encode_nanos: u128,
}

impl Engine {
    fn new(
        order: Order,
        width: Width,
        batch: usize,
        framing: Framing,
        level: i32,
        dict: Option<&[u8]>,
    ) -> Self {
        let (compressor, stream) = match framing {
            Framing::Stream { window, .. } => {
                let mut enc =
                    zstd::stream::write::Encoder::new(CountingSink(0), level).expect("encoder");
                enc.set_parameter(zstd::zstd_safe::CParameter::WindowLog(
                    window.trailing_zeros(),
                ))
                .expect("window log");
                enc.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
                    .expect("content size");
                (None, Some(enc))
            }
            Framing::Frames {
                dict: want_dict, ..
            } => {
                let mut c = match (want_dict, dict) {
                    (true, Some(d)) => {
                        zstd::bulk::Compressor::with_dictionary(level, d).expect("dict compressor")
                    }
                    _ => zstd::bulk::Compressor::new(level).expect("compressor"),
                };
                c.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
                    .expect("content size");
                (Some(c), None)
            }
        };
        Self {
            order,
            width,
            batch,
            framing,
            compressor,
            stream,
            out_bytes: 0,
            slot_bytes: [0; N_SLOTS],
            encode_nanos: 0,
        }
    }

    fn push_batch(&mut self, batch: &[Vec<u8>], scratch: &mut Vec<u8>) {
        let t0 = Instant::now();
        match self.framing {
            Framing::Frames {
                per_field: true, ..
            } => {
                let c = self.compressor.as_mut().expect("compressor");
                for (i, part) in batch.iter().enumerate() {
                    if part.is_empty() {
                        continue;
                    }
                    // Each block carries a length per field in its manifest.
                    let n = 3 + c.compress(part).expect("compress").len() as u64;
                    self.out_bytes += n;
                    self.slot_bytes[i] += n;
                }
            }
            Framing::Frames {
                per_field: false, ..
            } => {
                pack(batch, scratch);
                let c = self.compressor.as_mut().expect("compressor");
                self.out_bytes += c.compress(scratch).expect("compress").len() as u64;
            }
            Framing::Stream { flush, .. } => {
                pack(batch, scratch);
                let enc = self.stream.as_mut().expect("stream");
                enc.write_all(scratch).expect("write");
                if flush {
                    enc.flush().expect("flush");
                }
            }
        }
        self.encode_nanos += t0.elapsed().as_nanos();
    }

    fn finish(&mut self) {
        if let Some(enc) = self.stream.take() {
            self.out_bytes += enc.finish().expect("finish").0;
        }
    }

    /// Bytes a reader has to hold to yield records: the compressor's window,
    /// plus whatever the layout forces it to materialise at once. A transposed
    /// batch has to be decoded whole before its first record exists; a row batch
    /// is consumed a record at a time.
    fn reader_bytes(&self) -> usize {
        let window = match self.framing {
            Framing::Stream { window, .. } => window,
            Framing::Frames { .. } => self.batch,
        };
        let materialised = match self.order {
            Order::Col => self.batch,
            Order::Row => 4096,
        };
        window + materialised
    }

    fn framing_name(&self) -> String {
        match self.framing {
            Framing::Frames { per_field, dict } => format!(
                "frames/{}{}",
                if per_field { "per-field" } else { "one" },
                if dict { "+dict" } else { "" }
            ),
            Framing::Stream { window, flush } => format!(
                "stream/{}{}",
                window / 1024,
                if flush { "/flush" } else { "/noflush" }
            ),
        }
    }
}

fn pack(batch: &[Vec<u8>], scratch: &mut Vec<u8>) {
    scratch.clear();
    if batch.len() > 1 {
        for part in batch {
            put_varint(scratch, part.len() as u64);
        }
    }
    for part in batch {
        scratch.extend_from_slice(part);
    }
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

const DICT_MAX_BATCH: usize = 256 * 1024;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: <file.psp> [--max-records N] [--level L] [--batches K,K,..]");
    let mut max_records = u64::MAX;
    let mut level = 9i32;
    let mut batches_kib: Vec<usize> = vec![8, 32, 128, 512, 2048];
    let mut big_window_kib = 8192usize;
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--max-records" => max_records = args.next().expect("N").parse().expect("N"),
            "--level" => level = args.next().expect("L").parse().expect("L"),
            "--big-window-kib" => {
                big_window_kib = args.next().expect("K").parse().expect("K");
            }
            "--batches" => {
                batches_kib = args
                    .next()
                    .expect("K,K")
                    .split(',')
                    .map(|s| s.parse().expect("K"))
                    .collect()
            }
            other => panic!("unknown flag {other}"),
        }
    }

    let file_bytes = std::fs::metadata(&path).expect("stat").len();

    let combos: Vec<(Order, Width, usize)> = batches_kib
        .iter()
        .flat_map(|&kib| {
            [Order::Row, Order::Col].into_iter().flat_map(move |o| {
                [Width::Fixed, Width::Varint, Width::Quant]
                    .into_iter()
                    .map(move |w| (o, w, kib * 1024))
            })
        })
        .collect();

    eprintln!("training dictionaries…");
    let dicts: Vec<Option<Vec<u8>>> = combos
        .iter()
        .map(|&(o, w, b)| {
            if b <= DICT_MAX_BATCH {
                train_dictionary(&path, o, w, b)
            } else {
                None
            }
        })
        .collect();

    // One builder and its engines per (order, width, batch): the batch cut
    // depends on all three.
    let mut builders: Vec<BatchBuilder> = combos
        .iter()
        .map(|&(o, w, _)| BatchBuilder::new(o, w))
        .collect();
    let mut engine_sets: Vec<Vec<Engine>> = combos
        .iter()
        .zip(&dicts)
        .map(|(&(o, w, b), dict)| {
            let mut v = vec![
                Engine::new(
                    o,
                    w,
                    b,
                    Framing::Frames {
                        per_field: false,
                        dict: false,
                    },
                    level,
                    None,
                ),
                Engine::new(
                    o,
                    w,
                    b,
                    Framing::Stream {
                        window: b,
                        flush: true,
                    },
                    level,
                    None,
                ),
                Engine::new(
                    o,
                    w,
                    b,
                    Framing::Stream {
                        window: big_window_kib * 1024,
                        flush: true,
                    },
                    level,
                    None,
                ),
                Engine::new(
                    o,
                    w,
                    b,
                    Framing::Stream {
                        window: big_window_kib * 1024,
                        flush: false,
                    },
                    level,
                    None,
                ),
            ];
            if o == Order::Col {
                v.push(Engine::new(
                    o,
                    w,
                    b,
                    Framing::Frames {
                        per_field: true,
                        dict: false,
                    },
                    level,
                    None,
                ));
            }
            if let Some(d) = dict {
                v.push(Engine::new(
                    o,
                    w,
                    b,
                    Framing::Frames {
                        per_field: false,
                        dict: true,
                    },
                    level,
                    Some(d),
                ));
            }
            v
        })
        .collect();

    let file = File::open(&path).expect("open");
    let mut reader = PspReader::new(BufReader::with_capacity(1 << 20, file)).expect("psp header");
    let sample = reader.header().sample.clone();
    let n_blocks = reader.block_index().len();

    let mut batch: Vec<Vec<u8>> = Vec::new();
    let mut scratch: Vec<u8> = Vec::new();
    let mut n_records = 0u64;
    let mut n_alleles = 0u64;
    let mut n_chain_ids = 0u64;
    let mut n_ref_chain_ids = 0u64;
    let mut cur_chrom = u32::MAX;

    eprintln!("encoding…");
    let t_start = Instant::now();
    {
        let records = reader.records();
        for rec in records {
            if n_records >= max_records {
                break;
            }
            let rec = rec.expect("record");
            n_records += 1;
            n_alleles += rec.alleles.len() as u64;
            n_chain_ids += rec
                .alleles
                .iter()
                .map(|a| a.chain_ids.len() as u64)
                .sum::<u64>();
            n_ref_chain_ids += rec.alleles[0].chain_ids.len() as u64;
            let chrom_changed = rec.chrom_id != cur_chrom;
            cur_chrom = rec.chrom_id;

            for i in 0..builders.len() {
                let target = combos[i].2;
                // A batch never crosses a chromosome, matching the production
                // block rule, so cut before the first record of a new one.
                if chrom_changed && builders[i].n_records > 0 {
                    builders[i].take(&mut batch);
                    for e in &mut engine_sets[i] {
                        e.push_batch(&batch, &mut scratch);
                    }
                }
                builders[i].push(&rec);
                if builders[i].size() >= target {
                    builders[i].take(&mut batch);
                    for e in &mut engine_sets[i] {
                        e.push_batch(&batch, &mut scratch);
                    }
                }
            }
        }
    }
    for i in 0..builders.len() {
        if builders[i].n_records > 0 {
            builders[i].take(&mut batch);
            for e in &mut engine_sets[i] {
                e.push_batch(&batch, &mut scratch);
            }
        }
        for e in &mut engine_sets[i] {
            e.finish();
        }
    }
    let elapsed = t_start.elapsed();

    println!("# file\t{path}");
    println!("# sample\t{sample}");
    println!("# on-disk-bytes\t{file_bytes}");
    println!("# blocks\t{n_blocks}");
    println!("# records\t{n_records}");
    println!("# alleles\t{n_alleles}");
    println!("# chain-ids\t{n_chain_ids}");
    println!("# chain-ids-on-ref-allele\t{n_ref_chain_ids}");
    println!("# zstd-level\t{level}");
    println!("# wall-seconds\t{:.1}", elapsed.as_secs_f64());

    // Uncompressed size per field, from any one builder of each width.
    for width in [Width::Fixed, Width::Varint, Width::Quant] {
        for (i, b) in builders.iter().enumerate() {
            if b.width == width && combos[i].0 == Order::Col && combos[i].2 == batches_kib[0] * 1024
            {
                let total: u64 = b.slot_totals.iter().sum();
                println!(
                    "# uncompressed\t{}\t{}\t{:.3}",
                    width.name(),
                    total,
                    total as f64 / n_records as f64
                );
                for (s, name) in SLOT_NAMES.iter().enumerate() {
                    println!(
                        "#   raw\t{}\t{}\t{}\t{:.3}",
                        width.name(),
                        name,
                        b.slot_totals[s],
                        b.slot_totals[s] as f64 / n_records as f64
                    );
                }
            }
        }
    }

    // Compressed bytes per field, from the per-field framings — where the file
    // actually is, before any layout choice.
    for (i, set) in engine_sets.iter().enumerate() {
        for e in set {
            if matches!(
                e.framing,
                Framing::Frames {
                    per_field: true,
                    ..
                }
            ) {
                println!(
                    "# per-field-compressed\t{}\t{}",
                    e.width.name(),
                    combos[i].2 / 1024
                );
                for (s, name) in SLOT_NAMES.iter().enumerate() {
                    println!(
                        "#   comp\t{}\t{}\t{}\t{}\t{:.3}",
                        e.width.name(),
                        combos[i].2 / 1024,
                        name,
                        e.slot_bytes[s],
                        e.slot_bytes[s] as f64 / n_records as f64
                    );
                }
            }
        }
    }

    println!("order\twidth\tbatch_kib\tframing\treader_kib\tbytes\tbytes_per_record\tenc_MBps");
    let mut rows: Vec<(f64, String)> = Vec::new();
    for (i, set) in engine_sets.iter().enumerate() {
        for e in set {
            let bpr = e.out_bytes as f64 / n_records as f64;
            let mbps = if e.encode_nanos > 0 {
                (e.out_bytes as f64) / (e.encode_nanos as f64 / 1e9) / 1e6
            } else {
                f64::NAN
            };
            rows.push((
                bpr,
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.1}",
                    e.order.name(),
                    e.width.name(),
                    combos[i].2 / 1024,
                    e.framing_name(),
                    e.reader_bytes() / 1024,
                    e.out_bytes,
                    bpr,
                    mbps
                ),
            ));
        }
    }
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (_, r) in rows {
        println!("{r}");
    }
}

/// Train a zstd dictionary on this layout's own batches — training on one
/// layout and applying it to another would measure the mismatch.
fn train_dictionary(path: &str, order: Order, width: Width, batch: usize) -> Option<Vec<u8>> {
    let file = File::open(path).ok()?;
    let mut reader = PspReader::new(BufReader::with_capacity(1 << 20, file)).ok()?;
    let mut builder = BatchBuilder::new(order, width);
    let want = (8 * 1024 * 1024 / batch).clamp(16, 256);
    let mut samples: Vec<Vec<u8>> = Vec::new();
    let mut parts = Vec::new();
    let mut scratch = Vec::new();
    {
        let records = reader.records();
        for rec in records {
            let Ok(rec) = rec else { break };
            builder.push(&rec);
            if builder.size() >= batch {
                builder.take(&mut parts);
                pack(&parts, &mut scratch);
                samples.push(scratch.clone());
                if samples.len() >= want {
                    break;
                }
            }
        }
    }
    if samples.len() < 8 {
        return None;
    }
    zstd::dict::from_samples(&samples, 112 * 1024).ok()
}
