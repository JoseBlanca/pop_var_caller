//! A working record-stream store, to check that the size the layout probe
//! projects survives a real encoder and decoder — and what it costs in time and
//! memory to read back.
//!
//! The store is deliberately the simplest thing that could work: records are
//! serialised one after another, cut into small independent zstd frames, and
//! each frame is compressed against a dictionary trained on this file's own
//! records and stored in its header. There is no transposition, no column
//! manifest, and no block bookkeeping — a reader decompresses one frame and
//! walks it record by record.
//!
//! Three phases, each a separate run so the reported peak memory belongs to one
//! of them:
//!
//! ```text
//! cargo run --release --example psp_row_stream_roundtrip -- FILE.psp encode tmp/x.ngr
//! cargo run --release --example psp_row_stream_roundtrip -- FILE.psp verify tmp/x.ngr
//! cargo run --release --example psp_row_stream_roundtrip -- FILE.psp decode tmp/x.ngr
//! cargo run --release --example psp_row_stream_roundtrip -- FILE.psp psp-read
//! ```
//!
//! `verify` walks both stores in lockstep and fails on the first record that
//! disagrees, so a size claim is never made about a store that cannot be read
//! back. The three quantised fields are compared against the tolerance their
//! quantisation implies; everything else must match exactly.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::time::Instant;

use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
use pop_var_caller::psp::PspReader;
use pop_var_caller::psp::writer::SnpKind;

const MAGIC: &[u8; 4] = b"NGR1";

/// How finely the three quantities that arrive as floating point are stored:
/// each is kept as an integer count of steps, and the step is one over the
/// scale. A scale of 16 on the window's mean coverage means 1/16 of a read.
/// The scales are written into the file, so a reader never has to be told them.
#[derive(Copy, Clone)]
/// A scale of `0.0` means **store the value as raw IEEE bytes**, not as a
/// count of steps. That is the setting that keeps psp mode bit-identical to
/// direct mode: an approximated field makes the two routes see different
/// numbers, and `run_streaming.md` §1.2 makes their agreement the oracle for
/// the whole psp path.
struct Scales {
    gc: f64,
    coverage: f64,
    q_sum: f64,
}

impl Default for Scales {
    fn default() -> Self {
        Self {
            gc: 10_000.0,
            coverage: 16.0,
            q_sum: 256.0,
        }
    }
}

/// Where one open store is inside its current frame: the chromosome, the last
/// position yielded, how many records are left, the running coverage and
/// read-id bases, and the byte offset reached.
type FrameCursor = (u32, u32, u64, i64, u64, usize);

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

fn put_svarint(buf: &mut Vec<u8>, v: i64) {
    put_varint(buf, ((v << 1) ^ (v >> 63)) as u64);
}

struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    #[inline]
    fn varint(&mut self) -> u64 {
        let mut v = 0u64;
        let mut shift = 0;
        loop {
            let b = self.bytes[self.at];
            self.at += 1;
            v |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                return v;
            }
            shift += 7;
        }
    }
    #[inline]
    fn svarint(&mut self) -> i64 {
        let v = self.varint();
        ((v >> 1) as i64) ^ -((v & 1) as i64)
    }
    #[inline]
    fn bytes(&mut self, n: usize) -> &'a [u8] {
        let s = &self.bytes[self.at..self.at + n];
        self.at += n;
        s
    }
}

// ---------------------------------------------------------------------
// One frame's worth of records
// ---------------------------------------------------------------------

#[derive(Default)]
struct FrameWriter {
    /// Write only what a cohort scan reads — the position and the summed
    /// non-reference support — and drop the rest of the record. This is the
    /// "light stream" of a two-part block, built here on its own so the time a
    /// scan would save can be measured against a full walk.
    light_only: bool,
    /// Prefix each record with its own byte length, so a reader that has read
    /// the position and the support can skip the rest without parsing it. This
    /// is the one-stream alternative to a separate summary stream: it buys the
    /// same skipped parsing with one decompressor instead of two. WRITE-SIDE
    /// ONLY — this measures what the prefix costs in bytes; no reader uses it.
    length_prefix: bool,
    scales: Scales,
    buf: Vec<u8>,
    n_records: u64,
    chrom_id: u32,
    first_pos: u32,
    prev_pos: u32,
    prev_cov_q: i64,
    last_chain: u64,
}

impl FrameWriter {
    fn is_empty(&self) -> bool {
        self.n_records == 0
    }

    fn push(&mut self, rec: &PileupRecord) {
        let prefix_at = self.buf.len();
        if self.n_records == 0 {
            self.chrom_id = rec.chrom_id;
            self.first_pos = rec.pos;
            self.prev_pos = rec.pos;
            self.prev_cov_q = 0;
            self.last_chain = 0;
            put_varint(&mut self.buf, 0);
        } else {
            put_varint(&mut self.buf, (rec.pos - self.prev_pos) as u64);
        }
        self.prev_pos = rec.pos;

        if self.light_only {
            // Allele 0 is the reference bucket, so the scan's number is the
            // support summed over the rest.
            let non_ref: u64 = rec
                .alleles
                .iter()
                .skip(1)
                .map(|a| u64::from(a.support.num_obs))
                .sum();
            put_varint(&mut self.buf, non_ref);
            self.n_records += 1;
            return;
        }
        put_varint(&mut self.buf, rec.alleles.len() as u64);
        if self.scales.gc == 0.0 {
            self.buf.extend_from_slice(&rec.windowed_gc.to_le_bytes());
        } else {
            let gc = if rec.windowed_gc.is_nan() {
                0
            } else {
                1 + (f64::from(rec.windowed_gc).clamp(0.0, 1.0) * self.scales.gc).round() as u64
            };
            put_varint(&mut self.buf, gc);
        }
        if self.scales.coverage == 0.0 {
            self.buf
                .extend_from_slice(&rec.windowed_coverage.to_le_bytes());
        } else {
            let cov = if rec.windowed_coverage.is_nan() {
                0
            } else {
                let q = (f64::from(rec.windowed_coverage).max(0.0) * self.scales.coverage).round()
                    as i64;
                let d = q - self.prev_cov_q;
                self.prev_cov_q = q;
                1 + (((d << 1) ^ (d >> 63)) as u64)
            };
            put_varint(&mut self.buf, cov);
        }

        for allele in &rec.alleles {
            put_varint(&mut self.buf, allele.seq.len() as u64);
            self.buf.extend_from_slice(&allele.seq);
            let s = &allele.support;
            put_varint(&mut self.buf, s.num_obs as u64);
            if self.scales.q_sum == 0.0 {
                self.buf.extend_from_slice(&s.q_sum.to_le_bytes());
            } else {
                put_svarint(&mut self.buf, (s.q_sum * self.scales.q_sum).round() as i64);
            }
            put_varint(&mut self.buf, s.fwd as u64);
            put_varint(&mut self.buf, s.placed_left as u64);
            put_varint(&mut self.buf, s.placed_start as u64);
            put_varint(&mut self.buf, s.mapq_sum as u64);
            put_varint(&mut self.buf, s.mapq_sum_sq);
            put_varint(&mut self.buf, allele.chain_ids.len() as u64);
            for &id in &allele.chain_ids {
                put_svarint(
                    &mut self.buf,
                    (id as i64).wrapping_sub(self.last_chain as i64),
                );
                self.last_chain = id;
            }
        }
        if self.length_prefix {
            // Measure-only: append the length the prefix would have carried, so
            // its compressed cost is real even though nothing reads it.
            let len = (self.buf.len() - prefix_at) as u64;
            put_varint(&mut self.buf, len);
        }
        self.n_records += 1;
    }

    /// The frame's own header goes in front of the records: which chromosome it
    /// starts on, where, and how many records it holds. That is everything a
    /// reader needs to start here and nothing more.
    fn finish(&mut self, out: &mut Vec<u8>) {
        out.clear();
        put_varint(out, self.chrom_id as u64);
        put_varint(out, self.first_pos as u64);
        put_varint(out, self.n_records);
        out.extend_from_slice(&self.buf);
        self.buf.clear();
        self.n_records = 0;
    }
}

/// Reads one decompressed frame back into records, reusing the caller's
/// allocations — the reader holds one frame and one record, never a batch.
struct FrameReader<'a> {
    scales: Scales,
    cur: Cursor<'a>,
    chrom_id: u32,
    pos: u32,
    remaining: u64,
    prev_cov_q: i64,
    last_chain: u64,
    first: bool,
}

impl<'a> FrameReader<'a> {
    fn new(bytes: &'a [u8], scales: Scales) -> Self {
        let mut cur = Cursor { bytes, at: 0 };
        let chrom_id = cur.varint() as u32;
        let first_pos = cur.varint() as u32;
        let remaining = cur.varint();
        Self {
            scales,
            cur,
            chrom_id,
            pos: first_pos,
            remaining,
            prev_cov_q: 0,
            last_chain: 0,
            first: true,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn resume(
        bytes: &'a [u8],
        scales: Scales,
        chrom_id: u32,
        pos: u32,
        remaining: u64,
        prev_cov_q: i64,
        last_chain: u64,
        at: usize,
    ) -> Self {
        Self {
            scales,
            cur: Cursor { bytes, at },
            chrom_id,
            pos,
            remaining,
            prev_cov_q,
            last_chain,
            first: at == 0,
        }
    }

    fn state(&self) -> FrameCursor {
        (
            self.chrom_id,
            self.pos,
            self.remaining,
            self.prev_cov_q,
            self.last_chain,
            self.cur.at,
        )
    }

    fn next(&mut self) -> Option<PileupRecord> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;
        let delta = self.cur.varint() as u32;
        if self.first {
            self.first = false;
        } else {
            self.pos += delta;
        }
        let n_alleles = self.cur.varint() as usize;
        let windowed_gc = if self.scales.gc == 0.0 {
            f32::from_le_bytes(self.cur.bytes(4).try_into().expect("gc bytes"))
        } else {
            let gc_code = self.cur.varint();
            if gc_code == 0 {
                f32::NAN
            } else {
                ((gc_code - 1) as f64 / self.scales.gc) as f32
            }
        };
        let windowed_coverage = if self.scales.coverage == 0.0 {
            f32::from_le_bytes(self.cur.bytes(4).try_into().expect("coverage bytes"))
        } else {
            let cov_code = self.cur.varint();
            if cov_code == 0 {
                f32::NAN
            } else {
                let z = cov_code - 1;
                let d = ((z >> 1) as i64) ^ -((z & 1) as i64);
                self.prev_cov_q += d;
                (self.prev_cov_q as f64 / self.scales.coverage) as f32
            }
        };

        let mut alleles = Vec::with_capacity(n_alleles);
        for _ in 0..n_alleles {
            let len = self.cur.varint() as usize;
            let seq = self.cur.bytes(len).to_vec();
            let num_obs = self.cur.varint() as u32;
            let q_sum = if self.scales.q_sum == 0.0 {
                f64::from_le_bytes(self.cur.bytes(8).try_into().expect("q_sum bytes"))
            } else {
                self.cur.svarint() as f64 / self.scales.q_sum
            };
            let fwd = self.cur.varint() as u32;
            let placed_left = self.cur.varint() as u32;
            let placed_start = self.cur.varint() as u32;
            let mapq_sum = self.cur.varint() as u32;
            let mapq_sum_sq = self.cur.varint();
            let n_chain = self.cur.varint() as usize;
            let mut chain_ids = Vec::with_capacity(n_chain);
            for _ in 0..n_chain {
                let d = self.cur.svarint();
                self.last_chain = (self.last_chain as i64).wrapping_add(d) as u64;
                chain_ids.push(self.last_chain);
            }
            alleles.push(AlleleObservation::new(
                seq,
                AlleleSupportStats::new(
                    num_obs,
                    q_sum,
                    fwd,
                    placed_left,
                    placed_start,
                    mapq_sum,
                    mapq_sum_sq,
                ),
                chain_ids,
            ));
        }
        let mut rec = PileupRecord::new(self.chrom_id, self.pos, alleles);
        rec.windowed_gc = windowed_gc;
        rec.windowed_coverage = windowed_coverage;
        Some(rec)
    }
}

// ---------------------------------------------------------------------
// Phases
// ---------------------------------------------------------------------

fn peak_rss_kib() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            return rest.split_whitespace().next()?.parse().ok();
        }
    }
    None
}

fn train_dictionary(
    psp: &str,
    frame_bytes: usize,
    dict_bytes: usize,
    scales: Scales,
) -> Option<Vec<u8>> {
    let mut reader =
        PspReader::new(BufReader::with_capacity(1 << 20, File::open(psp).ok()?)).ok()?;
    let mut fw = FrameWriter {
        scales,
        ..FrameWriter::default()
    };
    let mut samples = Vec::new();
    let mut out = Vec::new();
    let mut cur_chrom = u32::MAX;
    let records = reader.records();
    for rec in records {
        let Ok(rec) = rec else { break };
        if (rec.chrom_id != cur_chrom || fw.buf.len() >= frame_bytes) && !fw.is_empty() {
            fw.finish(&mut out);
            samples.push(out.clone());
            if samples.len() >= 256 {
                break;
            }
        }
        cur_chrom = rec.chrom_id;
        fw.push(&rec);
    }
    if samples.len() < 8 {
        return None;
    }
    zstd::dict::from_samples(&samples, dict_bytes).ok()
}

fn encode(
    psp: &str,
    out_path: &str,
    frame_bytes: usize,
    level: i32,
    dict_bytes: usize,
    scales: Scales,
) {
    let t0 = Instant::now();
    let dict = train_dictionary(psp, frame_bytes, dict_bytes, scales).unwrap_or_default();
    let mut comp = if dict.is_empty() {
        zstd::bulk::Compressor::new(level).expect("compressor")
    } else {
        zstd::bulk::Compressor::with_dictionary(level, &dict).expect("compressor")
    };
    comp.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
        .expect("content size");

    let mut sink = BufWriter::with_capacity(1 << 20, File::create(out_path).expect("create"));
    sink.write_all(MAGIC).expect("write");
    for scale in [scales.gc, scales.coverage, scales.q_sum] {
        sink.write_all(&scale.to_le_bytes()).expect("write");
    }
    sink.write_all(&(dict.len() as u32).to_le_bytes())
        .expect("write");
    sink.write_all(&dict).expect("write");

    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(psp).expect("open"),
    ))
    .expect("psp");
    let mut fw = FrameWriter {
        scales,
        ..FrameWriter::default()
    };
    let mut raw = Vec::new();
    let mut n_records = 0u64;
    let mut n_frames = 0u64;
    let mut cur_chrom = u32::MAX;
    {
        let records = reader.records();
        for rec in records {
            let rec = rec.expect("record");
            // A frame never crosses a chromosome, so a reader that seeks to one
            // needs no carried-in state.
            if (rec.chrom_id != cur_chrom || fw.buf.len() >= frame_bytes) && !fw.is_empty() {
                fw.finish(&mut raw);
                write_frame(&mut sink, &mut comp, &raw);
                n_frames += 1;
            }
            cur_chrom = rec.chrom_id;
            fw.push(&rec);
            n_records += 1;
        }
    }
    if !fw.is_empty() {
        fw.finish(&mut raw);
        write_frame(&mut sink, &mut comp, &raw);
        n_frames += 1;
    }
    sink.flush().expect("flush");
    drop(sink);

    let bytes = std::fs::metadata(out_path).expect("stat").len();
    let psp_bytes = std::fs::metadata(psp).expect("stat").len();
    println!("phase\tencode");
    println!("records\t{n_records}");
    println!("frames\t{n_frames}");
    println!("dict-bytes\t{}", dict.len());
    println!("out-bytes\t{bytes}");
    println!("bytes-per-record\t{:.3}", bytes as f64 / n_records as f64);
    println!("psp-bytes\t{psp_bytes}");
    println!(
        "psp-bytes-per-record\t{:.3}",
        psp_bytes as f64 / n_records as f64
    );
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

fn write_frame(sink: &mut impl Write, comp: &mut zstd::bulk::Compressor<'static>, raw: &[u8]) {
    let out = comp.compress(raw).expect("compress");
    sink.write_all(&(out.len() as u32).to_le_bytes())
        .expect("write");
    sink.write_all(&(raw.len() as u32).to_le_bytes())
        .expect("write");
    sink.write_all(&out).expect("write");
}

struct StoreReader {
    scales: Scales,
    src: BufReader<File>,
    dec: zstd::bulk::Decompressor<'static>,
    comp_buf: Vec<u8>,
    raw_buf: Vec<u8>,
}

impl StoreReader {
    fn open(path: &str) -> Self {
        let mut src = BufReader::with_capacity(1 << 16, File::open(path).expect("open"));
        let mut magic = [0u8; 4];
        src.read_exact(&mut magic).expect("magic");
        assert_eq!(&magic, MAGIC, "not a record store");
        let mut scale_bytes = [0u8; 8];
        let mut read_scale = || {
            src.read_exact(&mut scale_bytes).expect("scale");
            f64::from_le_bytes(scale_bytes)
        };
        let scales = Scales {
            gc: read_scale(),
            coverage: read_scale(),
            q_sum: read_scale(),
        };
        let mut len = [0u8; 4];
        src.read_exact(&mut len).expect("dict len");
        let dict_len = u32::from_le_bytes(len) as usize;
        let mut dict = vec![0u8; dict_len];
        src.read_exact(&mut dict).expect("dict");
        let dec = if dict.is_empty() {
            zstd::bulk::Decompressor::new().expect("decompressor")
        } else {
            zstd::bulk::Decompressor::with_dictionary(&dict).expect("decompressor")
        };
        Self {
            scales,
            src,
            dec,
            comp_buf: Vec::new(),
            raw_buf: Vec::new(),
        }
    }

    /// Pull the next frame into the reusable buffer. Returns false at the end.
    fn next_frame(&mut self) -> bool {
        let mut hdr = [0u8; 8];
        match self.src.read_exact(&mut hdr) {
            Ok(()) => {}
            Err(_) => return false,
        }
        let comp_len = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as usize;
        let raw_len = u32::from_le_bytes(hdr[4..8].try_into().unwrap()) as usize;
        self.comp_buf.resize(comp_len, 0);
        self.src.read_exact(&mut self.comp_buf).expect("frame");
        self.raw_buf.resize(raw_len, 0);
        let n = self
            .dec
            .decompress_to_buffer(&self.comp_buf, &mut self.raw_buf)
            .expect("decompress");
        assert_eq!(n, raw_len, "frame decompressed to the wrong length");
        true
    }
}

fn decode(path: &str) {
    let t0 = Instant::now();
    let mut store = StoreReader::open(path);
    let mut n_records = 0u64;
    let mut checksum = 0u64;
    while store.next_frame() {
        let mut fr = FrameReader::new(&store.raw_buf, store.scales);
        while let Some(rec) = fr.next() {
            n_records += 1;
            checksum = checksum
                .wrapping_add(rec.pos as u64)
                .wrapping_add(rec.alleles[0].support.num_obs as u64);
        }
    }
    println!("phase\tdecode");
    println!("records\t{n_records}");
    println!("checksum\t{checksum}");
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    println!(
        "records-per-second\t{:.0}",
        n_records as f64 / t0.elapsed().as_secs_f64()
    );
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

fn psp_read(psp: &str) {
    let t0 = Instant::now();
    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(psp).expect("open"),
    ))
    .expect("psp");
    let mut n_records = 0u64;
    let mut checksum = 0u64;
    {
        let records = reader.records();
        for rec in records {
            let rec = rec.expect("record");
            n_records += 1;
            checksum = checksum
                .wrapping_add(rec.pos as u64)
                .wrapping_add(rec.alleles[0].support.num_obs as u64);
        }
    }
    println!("phase\tpsp-read");
    println!("records\t{n_records}");
    println!("checksum\t{checksum}");
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    println!(
        "records-per-second\t{:.0}",
        n_records as f64 / t0.elapsed().as_secs_f64()
    );
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

fn verify(psp: &str, path: &str) {
    let mut store = StoreReader::open(path);
    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(psp).expect("open"),
    ))
    .expect("psp");
    let mut records = reader.records();
    let mut n = 0u64;
    let mut worst_gc = 0.0f64;
    let mut worst_cov = 0.0f64;
    let mut worst_q = 0.0f64;
    while store.next_frame() {
        let mut fr = FrameReader::new(&store.raw_buf, store.scales);
        while let Some(got) = fr.next() {
            let want = records
                .next()
                .expect("psp ran out of records first")
                .expect("record");
            assert_eq!(got.chrom_id, want.chrom_id, "record {n}: chromosome");
            assert_eq!(got.pos, want.pos, "record {n}: position");
            assert_eq!(got.alleles.len(), want.alleles.len(), "record {n}: alleles");
            track(&mut worst_gc, got.windowed_gc, want.windowed_gc);
            track(
                &mut worst_cov,
                got.windowed_coverage,
                want.windowed_coverage,
            );
            for (a, b) in got.alleles.iter().zip(&want.alleles) {
                assert_eq!(a.seq, b.seq, "record {n}: allele sequence");
                assert_eq!(a.chain_ids, b.chain_ids, "record {n}: chain ids");
                assert_eq!(a.support.num_obs, b.support.num_obs, "record {n}: num_obs");
                assert_eq!(a.support.fwd, b.support.fwd, "record {n}: fwd");
                assert_eq!(
                    a.support.placed_left, b.support.placed_left,
                    "record {n}: placed_left"
                );
                assert_eq!(
                    a.support.placed_start, b.support.placed_start,
                    "record {n}: placed_start"
                );
                assert_eq!(
                    a.support.mapq_sum, b.support.mapq_sum,
                    "record {n}: mapq_sum"
                );
                assert_eq!(
                    a.support.mapq_sum_sq, b.support.mapq_sum_sq,
                    "record {n}: mapq_sum_sq"
                );
                let d = (a.support.q_sum - b.support.q_sum).abs();
                if d > worst_q {
                    worst_q = d;
                }
            }
            n += 1;
        }
    }
    assert!(
        records.next().is_none(),
        "the record store ran out before the .psp did"
    );
    println!("phase\tverify");
    println!("records\t{n}");
    println!("worst-gc-error\t{worst_gc:.6}");
    println!("worst-coverage-error\t{worst_cov:.6}");
    println!("worst-q-sum-error\t{worst_q:.6}");
}

fn track(worst: &mut f64, got: f32, want: f32) {
    if want.is_nan() {
        assert!(got.is_nan(), "a missing window came back as a number");
        return;
    }
    let d = (f64::from(got) - f64::from(want)).abs();
    if d > *worst {
        *worst = d;
    }
}

/// Open every per-sample store in a directory at once and walk them in
/// lockstep, one record each per round, the way a cohort merge does. This is
/// the shape that decides the memory bill: it is paid once per open sample.
fn many(dir: &str, kind: &str, limit: usize) {
    let mut paths: Vec<std::path::PathBuf> = std::fs::read_dir(dir)
        .expect("read dir")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some(kind))
        .collect();
    paths.sort();
    paths.truncate(limit);
    assert!(!paths.is_empty(), "no .{kind} files in {dir}");

    let t0 = Instant::now();
    let mut n_records = 0u64;
    let mut checksum = 0u64;
    let n_open = paths.len();
    if kind == "psp" {
        let mut readers: Vec<_> = paths
            .iter()
            .map(|p| {
                PspReader::new(BufReader::with_capacity(
                    1 << 20,
                    File::open(p).expect("open"),
                ))
                .expect("psp")
                .into_records_of::<SnpKind>()
            })
            .collect();
        let mut live = readers.len();
        let mut done = vec![false; readers.len()];
        while live > 0 {
            for (i, it) in readers.iter_mut().enumerate() {
                if done[i] {
                    continue;
                }
                match it.next() {
                    Some(rec) => {
                        let rec = rec.expect("record");
                        n_records += 1;
                        checksum = checksum.wrapping_add(rec.pos as u64);
                    }
                    None => {
                        done[i] = true;
                        live -= 1;
                    }
                }
            }
        }
    } else {
        let mut stores: Vec<StoreReader> = paths
            .iter()
            .map(|p| StoreReader::open(p.to_str().unwrap()))
            .collect();
        // Each open store keeps its own frame cursor: the frame bytes and the
        // offset within them.
        let mut cursors: Vec<Option<FrameCursor>> = vec![None; stores.len()];
        let mut done = vec![false; stores.len()];
        let mut live = stores.len();
        while live > 0 {
            for i in 0..stores.len() {
                if done[i] {
                    continue;
                }
                if cursors[i].is_none() {
                    if !stores[i].next_frame() {
                        done[i] = true;
                        live -= 1;
                        continue;
                    }
                    let mut c = Cursor {
                        bytes: &stores[i].raw_buf,
                        at: 0,
                    };
                    let chrom = c.varint() as u32;
                    let first = c.varint() as u32;
                    let n = c.varint();
                    cursors[i] = Some((chrom, first, n, 0, 0, c.at));
                }
                let (chrom, pos, left, prev_cov, last_chain, at) = cursors[i].unwrap();
                let mut fr = FrameReader::resume(
                    &stores[i].raw_buf,
                    stores[i].scales,
                    chrom,
                    pos,
                    left,
                    prev_cov,
                    last_chain,
                    at,
                );
                match fr.next() {
                    Some(rec) => {
                        n_records += 1;
                        checksum = checksum.wrapping_add(rec.pos as u64);
                        cursors[i] = Some(fr.state());
                    }
                    None => cursors[i] = None,
                }
            }
        }
    }
    println!("phase\tmany-{kind}");
    println!("open-files\t{n_open}");
    println!("records\t{n_records}");
    println!("checksum\t{checksum}");
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

// =====================================================================
// The streaming store: big blocks, a small window, nothing fully inflated
// =====================================================================
//
// The store above ties two things to one number, exactly as production's
// `.psp` does: a frame is both how far back the compressor may look and how
// much a reader must inflate before it can hand out a record. Shrinking it to
// save memory costs bytes, and it multiplies the block index besides.
//
// This store separates them. A **block** is large — a megabyte of records, so
// the compressor has plenty to work with and the index has few entries — while
// the **window**, the reach the compressor is allowed, is capped at write time
// and is what a reader has to hold. The reader never inflates a whole block:
// it pulls decompressed bytes into a small rolling buffer, parses records out
// of it, and drops them.
//
// Two conditions have to hold together for that to be worth anything, and only
// the first is zstd's doing:
//
//   1. do not inflate the whole block  — the capped window and the streaming
//      decoder below;
//   2. do not accumulate what you inflated — the reader hands each record to
//      the caller and keeps nothing, so the rolling buffer is the whole
//      decoded footprint.
//
// Satisfy only the first and the memory reappears in the caller's arrays.

const MAGIC_STREAM: &[u8; 4] = b"NGS1";

/// How much decompressed data the reader keeps in front of the parser. One
/// record has to fit; beyond that, bigger only means fewer refills. Set by
/// `--rolling-kib`, because it is one of the two things a reader's memory is
/// actually made of and the right value differs between a stream carrying whole
/// records and one carrying a single number per record.
static ROLLING_BYTES: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(64 * 1024);
/// How much compressed data is read from the file at a time. Set by
/// `--read-chunk-kib`.
static READ_CHUNK_BYTES: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(64 * 1024);

fn rolling_bytes() -> usize {
    ROLLING_BYTES.load(std::sync::atomic::Ordering::Relaxed)
}
fn read_chunk_bytes() -> usize {
    READ_CHUNK_BYTES.load(std::sync::atomic::Ordering::Relaxed)
}

fn encode_streaming(
    psp: &str,
    out_path: &str,
    block_bytes: usize,
    genomic_block_bp: u32,
    light_only: bool,
    length_prefix: bool,
    level: i32,
    window_log: u32,
    scales: Scales,
) {
    let t0 = Instant::now();
    let mut comp = zstd::bulk::Compressor::new(level).expect("compressor");
    comp.set_parameter(zstd::zstd_safe::CParameter::ContentSizeFlag(false))
        .expect("content size");
    // The cap that makes the block size and the reader's memory independent.
    // Without it zstd sizes its window from the block and a reader must hold
    // the whole thing.
    comp.set_parameter(zstd::zstd_safe::CParameter::WindowLog(window_log))
        .expect("window log");

    let mut sink = BufWriter::with_capacity(1 << 20, File::create(out_path).expect("create"));
    sink.write_all(MAGIC_STREAM).expect("write");
    for scale in [scales.gc, scales.coverage, scales.q_sum] {
        sink.write_all(&scale.to_le_bytes()).expect("write");
    }
    sink.write_all(&[u8::from(light_only)]).expect("write");
    sink.write_all(&window_log.to_le_bytes()).expect("write");

    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(psp).expect("open"),
    ))
    .expect("psp");
    let mut fw = FrameWriter {
        light_only,
        length_prefix,
        scales,
        ..FrameWriter::default()
    };
    let mut raw = Vec::new();
    let (mut n_records, mut n_blocks) = (0u64, 0u64);
    let mut cur_chrom = u32::MAX;
    let mut cur_grid = u64::MAX;
    {
        let records = reader.records();
        for rec in records {
            let rec = rec.expect("record");
            // The genomic block size is a grid on the reference coordinate, not
            // a running count: a block ends when a position crosses into the
            // next multiple of it. That is what makes every sample cut at the
            // same coordinates, so a cohort reader advancing over a region
            // touches one aligned block per sample rather than one in some and
            // two in others. A running count would not align.
            let grid = if genomic_block_bp > 0 {
                u64::from(rec.pos) / u64::from(genomic_block_bp)
            } else {
                0
            };
            let span_closed = genomic_block_bp > 0 && grid != cur_grid;
            // The byte ceiling is the secondary rule: at high depth one span can
            // hold a great deal of data.
            let bytes_closed = block_bytes > 0 && fw.buf.len() >= block_bytes;
            // A block never crosses a chromosome, so a reader that starts at
            // one carries in no state.
            if (rec.chrom_id != cur_chrom || span_closed || bytes_closed) && !fw.is_empty() {
                fw.finish(&mut raw);
                write_frame(&mut sink, &mut comp, &raw);
                n_blocks += 1;
            }
            cur_chrom = rec.chrom_id;
            cur_grid = grid;
            fw.push(&rec);
            n_records += 1;
        }
    }
    if !fw.is_empty() {
        fw.finish(&mut raw);
        write_frame(&mut sink, &mut comp, &raw);
        n_blocks += 1;
    }
    sink.flush().expect("flush");
    drop(sink);

    let bytes = std::fs::metadata(out_path).expect("stat").len();
    let psp_bytes = std::fs::metadata(psp).expect("stat").len();
    println!("phase\tencode-streaming");
    println!("records\t{n_records}");
    println!("blocks\t{n_blocks}");
    println!("genomic-block-bp\t{genomic_block_bp}");
    println!("light-only\t{light_only}");
    println!("window-bytes\t{}", 1usize << window_log);
    println!("out-bytes\t{bytes}");
    println!("bytes-per-record\t{:.3}", bytes as f64 / n_records as f64);
    println!("psp-bytes\t{psp_bytes}");
    println!(
        "psp-bytes-per-record\t{:.3}",
        psp_bytes as f64 / n_records as f64
    );
    // What an index over these blocks would cost per open sample, at the 24
    // bytes an entry production's costs.
    println!("index-bytes-24b-per-block\t{}", n_blocks * 24);
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

/// A cursor that reports running out of bytes instead of indexing past the
/// end. The streaming reader parses from a buffer that may hold only part of a
/// record, so "not enough yet" has to be an answer rather than a panic.
struct TryCursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> TryCursor<'a> {
    #[inline]
    fn varint(&mut self) -> Option<u64> {
        let mut v = 0u64;
        let mut shift = 0;
        loop {
            let b = *self.bytes.get(self.at)?;
            self.at += 1;
            v |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                return Some(v);
            }
            shift += 7;
        }
    }
    #[inline]
    fn svarint(&mut self) -> Option<i64> {
        let v = self.varint()?;
        Some(((v >> 1) as i64) ^ -((v & 1) as i64))
    }
    #[inline]
    fn bytes(&mut self, n: usize) -> Option<&'a [u8]> {
        let s = self.bytes.get(self.at..self.at + n)?;
        self.at += n;
        Some(s)
    }
}

/// One open sample. Everything it holds is bounded and none of it is the
/// block: a compressed read chunk, the rolling decompressed buffer, zstd's own
/// state sized by the window, and the record being built.
struct StreamingStore {
    light_only: bool,
    scales: Scales,
    src: BufReader<File>,
    dctx: zstd::zstd_safe::DCtx<'static>,
    in_buf: Vec<u8>,
    in_pos: usize,
    in_filled: usize,
    /// Compressed bytes left in the block being read.
    block_left: u64,
    out: Vec<u8>,
    out_at: usize,
    /// True once the current block's frame has been fully decompressed.
    block_done: bool,
    eof: bool,
    // The running state a record is parsed against.
    chrom_id: u32,
    pos: u32,
    remaining: u64,
    prev_cov_q: i64,
    last_chain: u64,
    first: bool,
}

impl StreamingStore {
    fn open(path: &std::path::Path) -> Self {
        let mut src = BufReader::with_capacity(1 << 16, File::open(path).expect("open"));
        let mut magic = [0u8; 4];
        src.read_exact(&mut magic).expect("magic");
        assert_eq!(&magic, MAGIC_STREAM, "not a streaming store");
        let mut b8 = [0u8; 8];
        let mut read_scale = || {
            src.read_exact(&mut b8).expect("scale");
            f64::from_le_bytes(b8)
        };
        let scales = Scales {
            gc: read_scale(),
            coverage: read_scale(),
            q_sum: read_scale(),
        };
        let mut b1 = [0u8; 1];
        src.read_exact(&mut b1).expect("light-only flag");
        let light_only = b1[0] != 0;
        let mut b4 = [0u8; 4];
        src.read_exact(&mut b4).expect("window log");
        let window_log = u32::from_le_bytes(b4);

        let mut dctx = zstd::zstd_safe::DCtx::create();
        // Refuse a file whose window is larger than we budgeted for, rather
        // than quietly allocating it.
        dctx.set_parameter(zstd::zstd_safe::DParameter::WindowLogMax(window_log))
            .expect("window log max");

        Self {
            light_only,
            scales,
            src,
            dctx,
            in_buf: vec![0u8; read_chunk_bytes()],
            in_pos: 0,
            in_filled: 0,
            block_left: 0,
            out: Vec::with_capacity(rolling_bytes()),
            out_at: 0,
            block_done: true,
            eof: false,
            chrom_id: 0,
            pos: 0,
            remaining: 0,
            prev_cov_q: 0,
            last_chain: 0,
            first: true,
        }
    }

    /// Start the next block, or report that the file is finished.
    fn next_block(&mut self) -> bool {
        let mut hdr = [0u8; 8];
        if self.src.read_exact(&mut hdr).is_err() {
            self.eof = true;
            return false;
        }
        self.block_left = u32::from_le_bytes(hdr[0..4].try_into().unwrap()) as u64;
        // hdr[4..8] is the uncompressed length, which a streaming reader does
        // not need and deliberately does not use: needing it would mean
        // sizing a buffer from the block.
        self.in_pos = 0;
        self.in_filled = 0;
        self.out.clear();
        self.out_at = 0;
        self.block_done = false;
        self.remaining = 0;
        true
    }

    /// Decompress more of the current block into the rolling buffer. Returns
    /// false when the block is finished and nothing more will arrive.
    fn pump(&mut self) -> bool {
        if self.block_done {
            return false;
        }
        // Drop what the parser has already consumed before asking for more.
        if self.out_at > 0 {
            self.out.drain(..self.out_at);
            self.out_at = 0;
        }
        if self.out.len() >= self.out.capacity() {
            // One record needs more than the rolling buffer holds. Grow rather
            // than fail — this is rare and bounded by the largest record.
            self.out.reserve(self.out.capacity());
        }
        if self.in_pos == self.in_filled {
            if self.block_left == 0 {
                self.block_done = true;
                return false;
            }
            let want = (self.block_left as usize).min(self.in_buf.len());
            self.src
                .read_exact(&mut self.in_buf[..want])
                .expect("block bytes");
            self.block_left -= want as u64;
            self.in_pos = 0;
            self.in_filled = want;
        }
        let mut input = zstd::zstd_safe::InBuffer {
            src: &self.in_buf[..self.in_filled],
            pos: self.in_pos,
        };
        let at = self.out.len();
        let mut output = zstd::zstd_safe::OutBuffer::around_pos(&mut self.out, at);
        let hint = self
            .dctx
            .decompress_stream(&mut output, &mut input)
            .expect("decompress");
        self.in_pos = input.pos;
        if hint == 0 && self.block_left == 0 && self.in_pos == self.in_filled {
            self.block_done = true;
        }
        true
    }

    /// The next record, or None at the end of the file. Nothing is retained:
    /// the caller gets the record and the store keeps only its running state.
    fn next(&mut self) -> Option<PileupRecord> {
        loop {
            if self.remaining == 0 {
                if self.block_done && self.out_at >= self.out.len() && !self.next_block() {
                    return None;
                }
                // A block opens with its chromosome, first position and count.
                loop {
                    let mut cur = TryCursor {
                        bytes: &self.out[self.out_at..],
                        at: 0,
                    };
                    match (cur.varint(), cur.varint(), cur.varint()) {
                        (Some(c), Some(p), Some(n)) => {
                            self.out_at += cur.at;
                            self.chrom_id = c as u32;
                            self.pos = p as u32;
                            self.remaining = n;
                            self.prev_cov_q = 0;
                            self.last_chain = 0;
                            self.first = true;
                            break;
                        }
                        _ => {
                            if !self.pump() {
                                return None;
                            }
                        }
                    }
                }
            }
            // Snapshot: a parse that runs out of bytes is retried from here
            // once more has arrived, so it must start from the same state.
            let (pos, prev_cov_q, last_chain, first) =
                (self.pos, self.prev_cov_q, self.last_chain, self.first);
            match self.parse_record() {
                Some(rec) => {
                    self.remaining -= 1;
                    return Some(rec);
                }
                None => {
                    self.pos = pos;
                    self.prev_cov_q = prev_cov_q;
                    self.last_chain = last_chain;
                    self.first = first;
                    if !self.pump() {
                        return None;
                    }
                }
            }
        }
    }

    /// Parse one record out of the rolling buffer, or None if it is not all
    /// there yet. Mutates the running state only on success.
    fn parse_record(&mut self) -> Option<PileupRecord> {
        let mut cur = TryCursor {
            bytes: &self.out[self.out_at..],
            at: 0,
        };
        let delta = cur.varint()? as u32;
        let mut pos = self.pos;
        let mut first = self.first;
        if first {
            first = false;
        } else {
            pos += delta;
        }
        if self.light_only {
            let non_ref = cur.varint()?;
            self.out_at += cur.at;
            self.pos = pos;
            self.first = first;
            let mut rec = PileupRecord::new(self.chrom_id, pos, Vec::new());
            // Park the scan number where the timing loop can see it.
            rec.windowed_coverage = non_ref as f32;
            return Some(rec);
        }
        let n_alleles = cur.varint()? as usize;
        let windowed_gc = if self.scales.gc == 0.0 {
            f32::from_le_bytes(cur.bytes(4)?.try_into().expect("gc bytes"))
        } else {
            let gc_code = cur.varint()?;
            if gc_code == 0 {
                f32::NAN
            } else {
                ((gc_code - 1) as f64 / self.scales.gc) as f32
            }
        };
        let mut prev_cov_q = self.prev_cov_q;
        let windowed_coverage = if self.scales.coverage == 0.0 {
            f32::from_le_bytes(cur.bytes(4)?.try_into().expect("coverage bytes"))
        } else {
            let cov_code = cur.varint()?;
            if cov_code == 0 {
                f32::NAN
            } else {
                let z = cov_code - 1;
                let d = ((z >> 1) as i64) ^ -((z & 1) as i64);
                prev_cov_q += d;
                (prev_cov_q as f64 / self.scales.coverage) as f32
            }
        };

        let mut last_chain = self.last_chain;
        let mut alleles = Vec::with_capacity(n_alleles);
        for _ in 0..n_alleles {
            let len = cur.varint()? as usize;
            let seq = cur.bytes(len)?.to_vec();
            let num_obs = cur.varint()? as u32;
            let q_sum = if self.scales.q_sum == 0.0 {
                f64::from_le_bytes(cur.bytes(8)?.try_into().expect("q_sum bytes"))
            } else {
                cur.svarint()? as f64 / self.scales.q_sum
            };
            let fwd = cur.varint()? as u32;
            let placed_left = cur.varint()? as u32;
            let placed_start = cur.varint()? as u32;
            let mapq_sum = cur.varint()? as u32;
            let mapq_sum_sq = cur.varint()?;
            let n_chain = cur.varint()? as usize;
            let mut chain_ids = Vec::with_capacity(n_chain);
            for _ in 0..n_chain {
                let d = cur.svarint()?;
                last_chain = (last_chain as i64).wrapping_add(d) as u64;
                chain_ids.push(last_chain);
            }
            alleles.push(AlleleObservation::new(
                seq,
                AlleleSupportStats::new(
                    num_obs,
                    q_sum,
                    fwd,
                    placed_left,
                    placed_start,
                    mapq_sum,
                    mapq_sum_sq,
                ),
                chain_ids,
            ));
        }
        // Committed: everything was there.
        self.out_at += cur.at;
        self.pos = pos;
        self.first = first;
        self.prev_cov_q = prev_cov_q;
        self.last_chain = last_chain;
        let mut rec = PileupRecord::new(self.chrom_id, pos, alleles);
        rec.windowed_gc = windowed_gc;
        rec.windowed_coverage = windowed_coverage;
        Some(rec)
    }
}

/// Decompress every byte of a store without parsing a single record.
///
/// This separates the two halves of a walk. A reader cannot decompress one
/// record out of a psp block: the block is one zstd frame and comes out
/// sequentially from its start, so the decompression is paid whatever the reader
/// intends to keep. Only the parse can be skipped. This measures the half that
/// cannot.
fn decode_raw(path: &str) {
    let t0 = Instant::now();
    let mut store = StreamingStore::open(std::path::Path::new(path));
    let mut bytes = 0u64;
    loop {
        if store.block_done && store.out_at >= store.out.len() && !store.next_block() {
            break;
        }
        // Consume whatever has been decompressed and ask for more.
        bytes += (store.out.len() - store.out_at) as u64;
        store.out_at = store.out.len();
        if !store.pump()
            && store.block_done
            && store.out_at >= store.out.len()
            && !store.next_block()
        {
            break;
        }
    }
    println!("phase\tdecode-raw");
    println!("decompressed-bytes\t{bytes}");
    println!("seconds\t{:.3}", t0.elapsed().as_secs_f64());
}

/// Walk one streaming store end to end, retaining nothing.
fn decode_streaming(path: &str) {
    let t0 = Instant::now();
    let mut store = StreamingStore::open(std::path::Path::new(path));
    let mut n_records = 0u64;
    let mut checksum = 0u64;
    while let Some(rec) = store.next() {
        n_records += 1;
        checksum = checksum.wrapping_add(rec.pos as u64);
    }
    println!("phase\tdecode-streaming");
    println!("records\t{n_records}");
    println!("checksum\t{checksum}");
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    println!(
        "records-per-second\t{:.0}",
        n_records as f64 / t0.elapsed().as_secs_f64()
    );
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

/// Read the streaming store and the `.psp` it was made from in lockstep and
/// fail on the first record that disagrees. A size or memory claim about a
/// store that cannot be read back is worth nothing.
fn verify_streaming(psp: &str, path: &str) {
    let mut reader = PspReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(psp).expect("open"),
    ))
    .expect("psp");
    let mut store = StreamingStore::open(std::path::Path::new(path));
    let scales = store.scales;
    let (mut n, mut worst_gc, mut worst_cov, mut worst_q) = (0u64, 0f64, 0f64, 0f64);
    {
        let records = reader.records();
        for want in records {
            let want = want.expect("record");
            let got = store.next().expect("store ended early");
            assert_eq!(got.chrom_id, want.chrom_id, "chrom at record {n}");
            assert_eq!(got.pos, want.pos, "pos at record {n}");
            assert_eq!(
                got.alleles.len(),
                want.alleles.len(),
                "allele count at record {n}"
            );
            track(&mut worst_gc, got.windowed_gc, want.windowed_gc);
            track(
                &mut worst_cov,
                got.windowed_coverage,
                want.windowed_coverage,
            );
            for (g, w) in got.alleles.iter().zip(&want.alleles) {
                assert_eq!(g.seq, w.seq, "allele sequence at record {n}");
                assert_eq!(g.chain_ids, w.chain_ids, "read names at record {n}");
                assert_eq!(g.support.num_obs, w.support.num_obs, "num_obs at {n}");
                assert_eq!(g.support.fwd, w.support.fwd, "fwd at {n}");
                assert_eq!(g.support.mapq_sum, w.support.mapq_sum, "mapq_sum at {n}");
                assert_eq!(
                    g.support.mapq_sum_sq, w.support.mapq_sum_sq,
                    "mapq_sum_sq at {n}"
                );
                let d = (g.support.q_sum - w.support.q_sum).abs();
                if d > worst_q {
                    worst_q = d;
                }
            }
            n += 1;
        }
    }
    assert!(
        store.next().is_none(),
        "store has records the .psp does not"
    );
    println!("phase\tverify-streaming");
    println!("records\t{n}");
    println!(
        "worst-gc-error\t{worst_gc:.6}\ttolerance\t{:.6}",
        1.0 / scales.gc
    );
    println!(
        "worst-coverage-error\t{worst_cov:.6}\ttolerance\t{:.6}",
        1.0 / scales.coverage
    );
    println!(
        "worst-q-sum-error\t{worst_q:.6}\ttolerance\t{:.6}",
        1.0 / scales.q_sum
    );
}

/// Open every streaming store in a directory at once and walk them in
/// lockstep, the way a cohort merge reads. This is the number the whole design
/// is for: it is paid once per open sample.
fn many_streaming(dir: &str, limit: usize, streams: usize) {
    let mut paths: Vec<std::path::PathBuf> = std::fs::read_dir(dir)
        .expect("read dir")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("ngs"))
        .collect();
    paths.sort();
    paths.truncate(limit);
    assert!(!paths.is_empty(), "no .ngs files in {dir}");

    let t0 = Instant::now();
    let n_open = paths.len();
    // `streams` live decoders per open sample. A block cut into several
    // separately-compressed pieces — the records in one, the cheap scan scalar
    // in another — needs one decoder, one read buffer and one rolling buffer
    // per piece, and they are all the same size whatever a piece carries. So
    // opening the same store `streams` times measures the multiplier without
    // needing the pieces to differ, which is what the memory question turns on.
    let mut stores: Vec<StreamingStore> = paths
        .iter()
        .flat_map(|p| (0..streams).map(move |_| StreamingStore::open(p)))
        .collect();
    let n_stores = stores.len();
    let mut live = vec![true; n_stores];
    let (mut n_records, mut checksum) = (0u64, 0u64);
    let mut any = true;
    while any {
        any = false;
        for i in 0..n_stores {
            if !live[i] {
                continue;
            }
            match stores[i].next() {
                Some(rec) => {
                    n_records += 1;
                    checksum = checksum.wrapping_add(rec.pos as u64);
                    any = true;
                }
                None => live[i] = false,
            }
        }
    }
    println!("phase\tmany-ngs");
    println!("open-files\t{n_open}");
    println!("streams-per-file\t{streams}");
    println!("live-decoders\t{n_stores}");
    println!("records\t{n_records}");
    println!("checksum\t{checksum}");
    println!("seconds\t{:.2}", t0.elapsed().as_secs_f64());
    if let Some(k) = peak_rss_kib() {
        println!("peak-rss-kib\t{k}");
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let psp = args
        .next()
        .expect("usage: <file.psp> <phase> [store] [flags]");
    let phase = args
        .next()
        .expect("phase: encode | decode | verify | psp-read");
    let store = args.next().unwrap_or_default();
    let mut frame_kib = 32usize;
    let mut level = 9i32;
    let mut dict_kib = 112usize;
    let mut limit = usize::MAX;
    // The two knobs this store separates: how much the compressor may look
    // back over (the block) and how much a reader must hold (the window).
    let mut block_kib = 1024usize;
    let mut window_log = 15u32; // 32 KiB
    // How many separately-compressed pieces one block is cut into.
    let mut streams = 1usize;
    // The genomic block size, in kilobases. 0 leaves the cut to the byte target
    // alone, which is what every measurement before 2026-08-25 used.
    let mut genomic_block_kb = 0u32;
    let mut light_only = false;
    let mut length_prefix = false;
    let mut scales = Scales::default();
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--frame-kib" => frame_kib = args.next().expect("K").parse().expect("K"),
            "--level" => level = args.next().expect("L").parse().expect("L"),
            "--dict-kib" => dict_kib = args.next().expect("K").parse().expect("K"),
            "--limit" => limit = args.next().expect("K").parse().expect("K"),
            "--block-kib" => block_kib = args.next().expect("K").parse().expect("K"),
            "--window-log" => window_log = args.next().expect("N").parse().expect("N"),
            "--streams" => streams = args.next().expect("N").parse().expect("N"),
            "--light-only" => light_only = true,
            "--length-prefix" => length_prefix = true,
            "--genomic-block-kb" => genomic_block_kb = args.next().expect("K").parse().expect("K"),
            "--rolling-kib" => ROLLING_BYTES.store(
                args.next().expect("K").parse::<usize>().expect("K") * 1024,
                std::sync::atomic::Ordering::Relaxed,
            ),
            "--read-chunk-kib" => READ_CHUNK_BYTES.store(
                args.next().expect("K").parse::<usize>().expect("K") * 1024,
                std::sync::atomic::Ordering::Relaxed,
            ),
            "--gc-scale" => scales.gc = args.next().expect("S").parse().expect("S"),
            "--coverage-scale" => scales.coverage = args.next().expect("S").parse().expect("S"),
            "--q-sum-scale" => scales.q_sum = args.next().expect("S").parse().expect("S"),
            other => panic!("unknown flag {other}"),
        }
    }
    match phase.as_str() {
        "encode" => encode(
            &psp,
            &store,
            frame_kib * 1024,
            level,
            dict_kib * 1024,
            scales,
        ),
        "decode" => decode(&store),
        "verify" => verify(&psp, &store),
        "psp-read" => psp_read(&psp),
        "many-psp" => many(&psp, "psp", limit),
        "many-ngr" => many(&psp, "ngr", limit),
        "encode-streaming" => encode_streaming(
            &psp,
            &store,
            block_kib * 1024,
            genomic_block_kb * 1000,
            light_only,
            length_prefix,
            level,
            window_log,
            scales,
        ),
        "decode-streaming" => decode_streaming(&store),
        "decode-raw" => decode_raw(&store),
        "verify-streaming" => verify_streaming(&psp, &store),
        "many-ngs" => many_streaming(&psp, limit, streams),
        other => panic!("unknown phase {other}"),
    }
}
