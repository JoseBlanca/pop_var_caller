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

        put_varint(&mut self.buf, rec.alleles.len() as u64);
        let gc = if rec.windowed_gc.is_nan() {
            0
        } else {
            1 + (f64::from(rec.windowed_gc).clamp(0.0, 1.0) * self.scales.gc).round() as u64
        };
        put_varint(&mut self.buf, gc);
        let cov = if rec.windowed_coverage.is_nan() {
            0
        } else {
            let q =
                (f64::from(rec.windowed_coverage).max(0.0) * self.scales.coverage).round() as i64;
            let d = q - self.prev_cov_q;
            self.prev_cov_q = q;
            1 + (((d << 1) ^ (d >> 63)) as u64)
        };
        put_varint(&mut self.buf, cov);

        for allele in &rec.alleles {
            put_varint(&mut self.buf, allele.seq.len() as u64);
            self.buf.extend_from_slice(&allele.seq);
            let s = &allele.support;
            put_varint(&mut self.buf, s.num_obs as u64);
            put_svarint(&mut self.buf, (s.q_sum * self.scales.q_sum).round() as i64);
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
        let gc_code = self.cur.varint();
        let windowed_gc = if gc_code == 0 {
            f32::NAN
        } else {
            ((gc_code - 1) as f64 / self.scales.gc) as f32
        };
        let cov_code = self.cur.varint();
        let windowed_coverage = if cov_code == 0 {
            f32::NAN
        } else {
            let z = cov_code - 1;
            let d = ((z >> 1) as i64) ^ -((z & 1) as i64);
            self.prev_cov_q += d;
            (self.prev_cov_q as f64 / self.scales.coverage) as f32
        };

        let mut alleles = Vec::with_capacity(n_alleles);
        for _ in 0..n_alleles {
            let len = self.cur.varint() as usize;
            let seq = self.cur.bytes(len).to_vec();
            let num_obs = self.cur.varint() as u32;
            let q_sum = self.cur.svarint() as f64 / self.scales.q_sum;
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
    let mut scales = Scales::default();
    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--frame-kib" => frame_kib = args.next().expect("K").parse().expect("K"),
            "--level" => level = args.next().expect("L").parse().expect("L"),
            "--dict-kib" => dict_kib = args.next().expect("K").parse().expect("K"),
            "--limit" => limit = args.next().expect("K").parse().expect("K"),
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
        other => panic!("unknown phase {other}"),
    }
}
