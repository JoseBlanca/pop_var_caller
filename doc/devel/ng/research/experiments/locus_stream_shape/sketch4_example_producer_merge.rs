//! **Sketch 4 — a block-filling generator on one side of the cohort merge.**
//!
//! Throwaway experiment for `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md`.
//! Deleted after the decision; commits no interface.
//!
//! Sketch 2 measured the merge two ways. Reading `.psp` files, folding one cheap column
//! across samples before materialising anything was 2.00× fewer instructions. Reading ng's
//! generator directly it was worth nothing — because the generator's only emitter returns a
//! fully owned record, so the fold's light column had to be derived from records that were
//! already built. Sketch 1 then taught the walk to fill blocks. Nobody put the two
//! together. This does.
//!
//! ## The four states, and the two floors that separate the merge from the producer
//!
//! | mode | producer | merge |
//! |---|---|---|
//! | `rec-drop` | records (shipped) | none — every locus dropped. Floor for A and B. |
//! | `rec-records` (**A**) | records | one owned record of lookahead per sample, O(N) head scan |
//! | `rec-fold` (**B**) | records | per-sample owned queues, light column derived per record, folded |
//! | `blk-drop` | blocks, no summary column | none. Floor for C. |
//! | `blk-drop-sum` | blocks, summary column on | none. Floor for D, and `blk-drop-sum − blk-drop` prices the column. |
//! | `blk-records` (**C**) | blocks, no summary | each locus refilled into a per-sample scratch record, then the same O(N) head scan |
//! | `blk-fold` (**D**) | blocks, summary on | the block's own summary column folded across samples; only variable positions ever read the heavy arrays |
//!
//! A and B are re-measurements of sketch 2's ng-fed arms, not citations: different
//! worktree, different binaries, and the two producers have to be compared inside one
//! experiment.
//!
//! ## Usage
//!
//! ```text
//! sketch4_producer_merge <mode> <full|digest> <ref.fa> <regions.bed> <n_regions> <a.cram> ...
//! ```

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::SampleLocusObservations;
use pop_var_caller::ng::locus_generation::block::{LocusBlock, LocusSink};
use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

const MIN_ALT_OBS: u32 = 2;
const MAX_ALLELES: usize = 4;
const EM_ITERS: usize = 5;
/// Loci arm B pulls into each sample's owned buffer per round. Sketch 2's value.
const FOLD_BATCH: usize = 4096;
/// The block budget. Sketch 1 measured instructions flat from 4 KiB to 16 MiB and peak
/// resident memory growing by about twice the budget, so this is the middle of a free range.
const BLOCK_BYTES_DEFAULT: usize = 256 * 1024;
/// Ceiling on the fold's dense window, in positions. The position key packs the contig id
/// into its high 32 bits, so an unguarded window would ask for terabytes at a contig change.
const WINDOW_CAP: u64 = 1 << 22;

type Gen = PileupGenerator<WindowedRefSeq, LeftAlignPreparer<WindowedRefSeq>>;

// =====================================================================
// Merged evidence — sketch 2's shape, unchanged, so the digests connect
// =====================================================================

#[derive(Default)]
struct AlleleEvidence {
    seq: Vec<u8>,
    chain_ids: Vec<u64>,
    obs: u32,
    q_sum: f64,
    fwd: u32,
    placed_left: u32,
    mapq_sum: u32,
    mapq_sum_sq: u64,
}

#[derive(Default)]
struct SampleEvidence {
    sample: u32,
    alleles: Vec<AlleleEvidence>,
    n_alleles: usize,
}

#[derive(Default)]
struct MergedPosition {
    key: u64,
    samples: Vec<SampleEvidence>,
    n: usize,
}

impl MergedPosition {
    fn begin(&mut self, key: u64) {
        self.key = key;
        self.n = 0;
    }
    fn open_sample(&mut self, sample: u32) -> &mut SampleEvidence {
        if self.n == self.samples.len() {
            self.samples.push(SampleEvidence::default());
        }
        let se = &mut self.samples[self.n];
        self.n += 1;
        se.sample = sample;
        se.n_alleles = 0;
        se
    }
}

impl SampleEvidence {
    fn open_allele(&mut self) -> &mut AlleleEvidence {
        if self.n_alleles == self.alleles.len() {
            self.alleles.push(AlleleEvidence::default());
        }
        let a = &mut self.alleles[self.n_alleles];
        self.n_alleles += 1;
        a
    }
}

/// Copy one owned locus's observations into the merged evidence. Returns bytes copied.
fn copy_locus(mp: &mut MergedPosition, sample: u32, locus: &SampleLocusObservations) -> u64 {
    let se = mp.open_sample(sample);
    let mut bytes = 0u64;
    for o in &locus.observations {
        let dst = se.open_allele();
        dst.seq.clear();
        dst.seq.extend_from_slice(&o.bases);
        dst.chain_ids.clear();
        dst.chain_ids.extend_from_slice(&o.chain_ids);
        dst.obs = o.num_obs;
        dst.q_sum = o.q_sum;
        dst.fwd = o.num_fwd;
        dst.placed_left = o.placed_left;
        dst.mapq_sum = o.mapq_sum;
        dst.mapq_sum_sq = o.mapq_sum_sq;
        bytes += dst.seq.len() as u64 + (dst.chain_ids.len() as u64) * 8 + 32;
    }
    bytes
}

/// The same, straight out of a block's columns — no record in between.
fn copy_locus_from_block(
    mp: &mut MergedPosition,
    sample: u32,
    block: &LocusBlock,
    i: usize,
) -> u64 {
    let se = mp.open_sample(sample);
    let mut bytes = 0u64;
    for j in block.observations_of(i) {
        let dst = se.open_allele();
        dst.seq.clear();
        dst.seq.extend_from_slice(block.bases_of(j));
        dst.chain_ids.clear();
        dst.chain_ids.extend_from_slice(block.chain_ids_of(j));
        dst.obs = block.obs_num_obs[j];
        dst.q_sum = block.obs_q_sum[j];
        dst.fwd = block.obs_num_fwd[j];
        dst.placed_left = block.obs_placed_left[j];
        dst.mapq_sum = block.obs_mapq_sum[j];
        dst.mapq_sum_sq = block.obs_mapq_sum_sq[j];
        bytes += dst.seq.len() as u64 + (dst.chain_ids.len() as u64) * 8 + 32;
    }
    bytes
}

/// The same, out of a scratch record the consumer refilled.
fn copy_locus_from_scratch(mp: &mut MergedPosition, sample: u32, s: &LocusScratch) -> u64 {
    let se = mp.open_sample(sample);
    let mut bytes = 0u64;
    for o in &s.observations[..s.used] {
        let dst = se.open_allele();
        dst.seq.clear();
        dst.seq.extend_from_slice(&o.bases);
        dst.chain_ids.clear();
        dst.chain_ids.extend_from_slice(&o.chain_ids);
        dst.obs = o.num_obs;
        dst.q_sum = o.q_sum;
        dst.fwd = o.num_fwd;
        dst.placed_left = o.placed_left;
        dst.mapq_sum = o.mapq_sum;
        dst.mapq_sum_sq = o.mapq_sum_sq;
        bytes += dst.seq.len() as u64 + (dst.chain_ids.len() as u64) * 8 + 32;
    }
    bytes
}

/// The light column ng's record path does not provide: summed observations whose bases
/// differ from the locus's reference bases. Costs a byte comparison per observation, paid
/// by the consumer because the producer handed over a record and nothing else.
#[inline]
fn nonref_obs_of(locus: &SampleLocusObservations) -> u32 {
    let mut s = 0u32;
    for o in &locus.observations {
        if o.bases.as_ref() != locus.reference_bases.as_ref() {
            s = s.saturating_add(o.num_obs);
        }
    }
    s
}

#[inline]
fn nonref_obs_of_scratch(s: &LocusScratch) -> u32 {
    let mut t = 0u32;
    for o in &s.observations[..s.used] {
        if o.bases.as_slice() != s.reference_bases.as_slice() {
            t = t.saturating_add(o.num_obs);
        }
    }
    t
}

#[inline]
fn key_of(locus: &SampleLocusObservations) -> u64 {
    (u64::from(locus.region.contig.0) << 32) | locus.region.start.0
}

#[inline]
fn key_in_block(block: &LocusBlock, i: usize) -> u64 {
    (u64::from(block.contig[i]) << 32) | u64::from(block.start[i])
}

// =====================================================================
// The consumer's own record buffer — sketch 1's arm C, one per sample
// =====================================================================

#[derive(Default)]
struct ScratchObservation {
    bases: Vec<u8>,
    complete: bool,
    runs: Vec<(u16, u16)>,
    read_group: u32,
    num_obs: u32,
    num_fwd: u32,
    q_sum: f64,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    placed_left: u32,
    chain_ids: Vec<u64>,
}

/// One locus, record-shaped, in storage the consumer owns and refills. This is the shape
/// that survived in production, rather than the borrow into the block that did not.
#[derive(Default)]
struct LocusScratch {
    contig: u32,
    start: u32,
    end: u32,
    reference_bases: Vec<u8>,
    reads_without_observation: u32,
    reads_discarded_by_cap: u32,
    observations: Vec<ScratchObservation>,
    used: usize,
    bytes_copied: u64,
}

impl LocusScratch {
    /// Refill from locus `i` of `block`. Everything is copied; nothing borrows the block.
    fn refill_from(&mut self, block: &LocusBlock, i: usize) {
        self.contig = block.contig[i];
        self.start = block.start[i];
        self.end = block.end[i];
        self.reads_without_observation = block.reads_without_observation[i];
        self.reads_discarded_by_cap = block.reads_discarded_by_cap[i];
        self.reference_bases.clear();
        self.reference_bases
            .extend_from_slice(block.reference_bases_of(i));
        let mut copied = block.reference_bases_of(i).len();
        self.used = 0;
        for j in block.observations_of(i) {
            if self.used == self.observations.len() {
                self.observations.push(ScratchObservation::default());
            }
            let slot = &mut self.observations[self.used];
            slot.bases.clear();
            slot.bases.extend_from_slice(block.bases_of(j));
            slot.complete = block.obs_witness_kind[j] == 0;
            slot.runs.clear();
            slot.runs.extend_from_slice(block.runs_of(j));
            slot.read_group = block.obs_read_group[j];
            slot.num_obs = block.obs_num_obs[j];
            slot.num_fwd = block.obs_num_fwd[j];
            slot.q_sum = block.obs_q_sum[j];
            slot.mapq_sum = block.obs_mapq_sum[j];
            slot.mapq_sum_sq = block.obs_mapq_sum_sq[j];
            slot.placed_left = block.obs_placed_left[j];
            slot.chain_ids.clear();
            slot.chain_ids.extend_from_slice(block.chain_ids_of(j));
            copied += slot.bases.len() + slot.runs.len() * 4 + slot.chain_ids.len() * 8 + 41;
            self.used += 1;
        }
        self.bytes_copied += copied as u64;
    }

    #[inline]
    fn key(&self) -> u64 {
        (u64::from(self.contig) << 32) | u64::from(self.start)
    }
}

// =====================================================================
// Counters + sink
// =====================================================================

#[derive(Default, Debug)]
struct Counters {
    loci_produced: u64,
    positions_seen: u64,
    positions_kept: u64,
    merge_objects: u64,
    merge_bytes: u64,
    /// Bytes the producer owns in the records it built (heap payload only) — record path.
    producer_bytes: u64,
    /// Bytes the producer wrote into blocks — block path.
    block_bytes: u64,
    /// Loci ever materialised as a record-shaped object by anyone.
    loci_materialised: u64,
    /// Blocks handed over.
    blocks: u64,
    block_payload_high_water: u64,
    /// Bytes arm C copied out of blocks into its own record-shaped buffers.
    scratch_bytes: u64,
    /// Emitted observations, over every locus — what the summary column's byte comparison
    /// is charged per.
    observations: u64,
    /// Reads at a position, `Σ num_obs`. Divided by loci this is the fixture's depth.
    visits: u64,
}

fn locus_owned_bytes(l: &SampleLocusObservations) -> u64 {
    let mut b = l.reference_bases.len() as u64;
    b += l.observations.len() as u64 * 96;
    for o in &l.observations {
        b += o.bases.len() as u64 + (o.chain_ids.len() as u64) * 8;
    }
    b
}

#[derive(Clone, Copy, PartialEq)]
enum Consumer {
    Digest,
    Full,
}

struct Sink {
    consumer: Consumer,
    digest: u64,
    table: Vec<Vec<u8>>,
    table_len: usize,
    counts: Vec<u32>,
    errs: Vec<f64>,
    gl: Vec<f64>,
    post: Vec<f64>,
    freqs: Vec<f64>,
    loglik_acc: f64,
}

#[inline]
fn fnv(h: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *h ^= u64::from(b);
        *h = h.wrapping_mul(0x100_0000_01b3);
    }
}

impl Sink {
    fn new(consumer: Consumer) -> Self {
        Self {
            consumer,
            digest: 0xcbf2_9ce4_8422_2325,
            table: Vec::new(),
            table_len: 0,
            counts: Vec::new(),
            errs: Vec::new(),
            gl: Vec::new(),
            post: Vec::new(),
            freqs: Vec::new(),
            loglik_acc: 0.0,
        }
    }

    fn accept(&mut self, mp: &MergedPosition) {
        let mut h = self.digest;
        fnv(&mut h, &mp.key.to_le_bytes());
        fnv(&mut h, &(mp.n as u32).to_le_bytes());
        for se in &mp.samples[..mp.n] {
            fnv(&mut h, &se.sample.to_le_bytes());
            fnv(&mut h, &(se.n_alleles as u32).to_le_bytes());
            for a in &se.alleles[..se.n_alleles] {
                fnv(&mut h, &a.seq);
                fnv(&mut h, &a.obs.to_le_bytes());
                fnv(&mut h, &a.q_sum.to_bits().to_le_bytes());
                fnv(&mut h, &a.fwd.to_le_bytes());
                fnv(&mut h, &a.placed_left.to_le_bytes());
                fnv(&mut h, &a.mapq_sum.to_le_bytes());
                fnv(&mut h, &a.mapq_sum_sq.to_le_bytes());
                for &c in &a.chain_ids {
                    fnv(&mut h, &c.to_le_bytes());
                }
            }
        }
        self.digest = h;
        if self.consumer == Consumer::Full {
            let ll = self.call(mp);
            self.loglik_acc += ll;
            let mut h = self.digest;
            fnv(&mut h, &ll.to_bits().to_le_bytes());
            self.digest = h;
        }
    }

    fn call(&mut self, mp: &MergedPosition) -> f64 {
        self.table_len = 0;
        for se in &mp.samples[..mp.n] {
            for a in &se.alleles[..se.n_alleles] {
                if a.obs == 0 {
                    continue;
                }
                let known = self.table[..self.table_len]
                    .iter()
                    .any(|t| t.as_slice() == a.seq.as_slice());
                if !known && self.table_len < MAX_ALLELES {
                    if self.table_len == self.table.len() {
                        self.table.push(Vec::new());
                    }
                    let slot = &mut self.table[self.table_len];
                    slot.clear();
                    slot.extend_from_slice(&a.seq);
                    self.table_len += 1;
                }
            }
        }
        let k = self.table_len.max(1);
        let n_gt = k * (k + 1) / 2;
        let s_n = mp.n;

        self.counts.clear();
        self.counts.resize(s_n * k, 0);
        self.errs.clear();
        self.errs.resize(s_n * k, 0.0);
        for (si, se) in mp.samples[..s_n].iter().enumerate() {
            for a in &se.alleles[..se.n_alleles] {
                if a.obs == 0 {
                    continue;
                }
                let Some(idx) = self.table[..self.table_len]
                    .iter()
                    .position(|t| t.as_slice() == a.seq.as_slice())
                else {
                    continue;
                };
                self.counts[si * k + idx] += a.obs;
                self.errs[si * k + idx] = (a.q_sum / f64::from(a.obs)).exp().clamp(1e-9, 0.5);
            }
        }

        self.gl.clear();
        self.gl.resize(s_n * n_gt, 0.0);
        for si in 0..s_n {
            let mut g = 0usize;
            for i in 0..k {
                for j in i..k {
                    let mut ll = 0.0f64;
                    for a in 0..k {
                        let n = self.counts[si * k + a];
                        if n == 0 {
                            continue;
                        }
                        let e = self.errs[si * k + a];
                        let dose = f64::from(u8::from(a == i)) + f64::from(u8::from(a == j));
                        let p = (dose / 2.0) * (1.0 - e) + (1.0 - dose / 2.0) * (e / 3.0);
                        ll += f64::from(n) * p.max(1e-300).ln();
                    }
                    self.gl[si * n_gt + g] = ll;
                    g += 1;
                }
            }
        }

        self.freqs.clear();
        self.freqs.resize(k, 1.0 / k as f64);
        self.post.clear();
        self.post.resize(n_gt, 0.0);
        let mut loglik = 0.0f64;
        for _ in 0..EM_ITERS {
            let mut new_f = [0.0f64; MAX_ALLELES];
            loglik = 0.0;
            for si in 0..s_n {
                let mut norm = 0.0f64;
                let mut g = 0usize;
                for i in 0..k {
                    for j in i..k {
                        let prior = if i == j {
                            self.freqs[i] * self.freqs[i]
                        } else {
                            2.0 * self.freqs[i] * self.freqs[j]
                        };
                        let w = prior * self.gl[si * n_gt + g].exp();
                        self.post[g] = w;
                        norm += w;
                        g += 1;
                    }
                }
                if norm <= 0.0 {
                    continue;
                }
                loglik += norm.ln();
                let mut g = 0usize;
                for i in 0..k {
                    for j in i..k {
                        let p = self.post[g] / norm;
                        new_f[i] += p;
                        new_f[j] += p;
                        g += 1;
                    }
                }
            }
            let tot: f64 = new_f[..k].iter().sum();
            if tot > 0.0 {
                for i in 0..k {
                    self.freqs[i] = new_f[i] / tot;
                }
            }
        }
        loglik
    }
}

// =====================================================================
// Per-sample generator
// =====================================================================

struct NgSample {
    generator: Gen,
    reads: SampleReads,
    head: Option<SampleLocusObservations>,
}

impl NgSample {
    fn pull(&mut self, c: &mut Counters) -> Result<(), Box<dyn std::error::Error>> {
        self.head = self.generator.next_locus(&self.reads)?;
        if let Some(l) = &self.head {
            c.loci_produced += 1;
            c.loci_materialised += 1;
            c.producer_bytes += locus_owned_bytes(l);
        }
        Ok(())
    }

    #[inline]
    fn key(&self) -> u64 {
        self.head.as_ref().map_or(u64::MAX, key_of)
    }
}

/// Where a sample is inside the block its generator most recently filled.
#[derive(Default)]
struct BlockCursor {
    at: usize,
    exhausted: bool,
    /// Only arm C uses this; arm D reads the block's columns directly.
    has_head: bool,
}

// =====================================================================
// The record producer's three modes — A, B, and their floor
// =====================================================================

fn mode_rec_drop(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    for region in regions {
        for s in samples.iter_mut() {
            s.generator.begin_segment(*region);
            while let Some(l) = s.generator.next_locus(&s.reads)? {
                c.loci_produced += 1;
                c.loci_materialised += 1;
                c.producer_bytes += locus_owned_bytes(&l);
                drop(l);
            }
        }
    }
    Ok(())
}

/// Not a timed mode. It walks with the record producer and counts what the fixture is made
/// of, so the report can state the depth rather than infer it. Kept out of every timed
/// path so no timed binary carries its counting.
fn mode_stats(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    for region in regions {
        for s in samples.iter_mut() {
            s.generator.begin_segment(*region);
            while let Some(l) = s.generator.next_locus(&s.reads)? {
                c.loci_produced += 1;
                c.observations += l.observations.len() as u64;
                for o in &l.observations {
                    c.visits += u64::from(o.num_obs);
                }
            }
        }
    }
    Ok(())
}

fn mode_rec_records(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    sink: &mut Sink,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = samples.len();
    let mut mp = MergedPosition::default();
    for region in regions {
        for s in samples.iter_mut() {
            s.generator.begin_segment(*region);
            s.head = None;
        }
        for si in 0..n {
            samples[si].pull(c)?;
        }
        loop {
            let mut min_key = u64::MAX;
            for s in samples.iter() {
                min_key = min_key.min(s.key());
            }
            if min_key == u64::MAX {
                break;
            }
            c.positions_seen += 1;

            let mut max_nonref = 0u32;
            for s in samples.iter() {
                if s.key() == min_key {
                    // PANIC-FREE: a finite key means the head is present.
                    max_nonref = max_nonref.max(nonref_obs_of(
                        s.head.as_ref().expect("a finite key means a present head"),
                    ));
                }
            }
            let kept = max_nonref >= MIN_ALT_OBS;
            if kept {
                c.positions_kept += 1;
                mp.begin(min_key);
            }

            for si in 0..n {
                if samples[si].key() != min_key {
                    continue;
                }
                if kept {
                    // PANIC-FREE: as above.
                    let locus = samples[si]
                        .head
                        .as_ref()
                        .expect("a finite key means a present head");
                    c.merge_bytes += copy_locus(&mut mp, si as u32, locus);
                    c.merge_objects += 1;
                }
                samples[si].pull(c)?;
            }
            if kept {
                sink.accept(&mp);
            }
        }
    }
    Ok(())
}

fn mode_rec_fold(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    sink: &mut Sink,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = samples.len();
    let mut buf: Vec<Vec<SampleLocusObservations>> = (0..n).map(|_| Vec::new()).collect();
    let mut light_key: Vec<Vec<u64>> = (0..n).map(|_| Vec::new()).collect();
    let mut light_nonref: Vec<Vec<u32>> = (0..n).map(|_| Vec::new()).collect();
    let mut cursor = vec![0usize; n];
    let mut exhausted = vec![false; n];
    let mut mat_cursor = vec![0usize; n];
    let mut window_key: Vec<u64> = Vec::new();
    let mut window_nonref: Vec<u32> = Vec::new();
    let mut mp = MergedPosition::default();

    for region in regions {
        for s in samples.iter_mut() {
            s.generator.begin_segment(*region);
        }
        for s in 0..n {
            buf[s].clear();
            light_key[s].clear();
            light_nonref[s].clear();
            cursor[s] = 0;
            exhausted[s] = false;
        }
        loop {
            for s in 0..n {
                while !exhausted[s] && buf[s].len() - cursor[s] < FOLD_BATCH {
                    match samples[s].generator.next_locus(&samples[s].reads)? {
                        Some(l) => {
                            c.loci_produced += 1;
                            c.loci_materialised += 1;
                            c.producer_bytes += locus_owned_bytes(&l);
                            light_key[s].push(key_of(&l));
                            light_nonref[s].push(nonref_obs_of(&l));
                            buf[s].push(l);
                        }
                        None => exhausted[s] = true,
                    }
                }
            }

            let mut watermark = u64::MAX;
            let mut any = false;
            for s in 0..n {
                if cursor[s] < buf[s].len() {
                    any = true;
                }
                if !exhausted[s] {
                    if let Some(&k) = light_key[s].last() {
                        watermark = watermark.min(k);
                    } else {
                        watermark = 0;
                    }
                }
            }
            if !any {
                break;
            }
            if watermark == u64::MAX {
                watermark = light_key
                    .iter()
                    .filter_map(|v| v.last().copied())
                    .max()
                    .unwrap_or(0);
            }

            let win_lo = (0..n)
                .filter(|&s| cursor[s] < light_key[s].len())
                .map(|s| light_key[s][cursor[s]])
                .min()
                .unwrap_or(watermark);
            let win_len = (watermark.saturating_sub(win_lo) + 1).min(WINDOW_CAP) as usize;
            let watermark = win_lo + win_len as u64 - 1;
            window_nonref.clear();
            window_nonref.resize(win_len, u32::MAX); // MAX = "no sample here"
            for s in 0..n {
                let mut i = cursor[s];
                while i < light_key[s].len() && light_key[s][i] <= watermark {
                    let off = (light_key[s][i] - win_lo) as usize;
                    let v = light_nonref[s][i];
                    let slot = &mut window_nonref[off];
                    *slot = if *slot == u32::MAX { v } else { (*slot).max(v) };
                    i += 1;
                }
            }
            window_key.clear();
            for (off, &v) in window_nonref.iter().enumerate() {
                if v != u32::MAX {
                    window_key.push(win_lo + off as u64);
                }
            }
            c.positions_seen += window_key.len() as u64;

            for s in 0..n {
                mat_cursor[s] = cursor[s];
            }
            for &key in &window_key {
                if window_nonref[(key - win_lo) as usize] < MIN_ALT_OBS {
                    continue;
                }
                c.positions_kept += 1;
                mp.begin(key);
                for s in 0..n {
                    let mut i = mat_cursor[s];
                    while i < light_key[s].len() && light_key[s][i] < key {
                        i += 1;
                    }
                    mat_cursor[s] = i;
                    if i >= light_key[s].len() || light_key[s][i] != key {
                        continue;
                    }
                    c.merge_bytes += copy_locus(&mut mp, s as u32, &buf[s][i]);
                    c.merge_objects += 1;
                }
                sink.accept(&mp);
            }

            for s in 0..n {
                let mut i = cursor[s];
                while i < light_key[s].len() && light_key[s][i] <= watermark {
                    i += 1;
                }
                buf[s].drain(..i);
                light_key[s].drain(..i);
                light_nonref[s].drain(..i);
                cursor[s] = 0;
            }
        }
    }
    Ok(())
}

// =====================================================================
// The block producer's three modes — C, D, and their floors
// =====================================================================

/// Tally a filled block into the counters. Called once per handover, before it is read.
fn account_block(block: &LocusBlock, c: &mut Counters) {
    c.blocks += 1;
    c.loci_produced += block.len() as u64;
    c.block_bytes += block.payload_bytes() as u64;
    c.block_payload_high_water = c.block_payload_high_water.max(block.payload_bytes() as u64);
}

fn mode_blk_drop(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    for region in regions {
        for s in samples.iter_mut() {
            s.generator.begin_segment(*region);
            loop {
                let more = s.generator.next_block(&s.reads)?;
                // PANIC-FREE: the sink is installed for every block mode.
                let sink = s
                    .generator
                    .block_sink_mut()
                    .expect("a block mode installs the sink");
                account_block(sink.block(), c);
                sink.block_mut().clear();
                if !more {
                    break;
                }
            }
        }
    }
    Ok(())
}

/// Advance sample `s` to its next locus, refilling the block when the cursor runs off its
/// end. The scratch buffer lives outside `samples`, so the block can be borrowed while it
/// is written.
fn blk_pull_scratch(
    sample: &mut NgSample,
    cursor: &mut BlockCursor,
    scratch: &mut LocusScratch,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    loop {
        // PANIC-FREE: the sink is installed for every block mode.
        let sink = sample
            .generator
            .block_sink_mut()
            .expect("a block mode installs the sink");
        if cursor.at < sink.block().len() {
            scratch.refill_from(sink.block(), cursor.at);
            cursor.at += 1;
            cursor.has_head = true;
            c.loci_materialised += 1;
            return Ok(());
        }
        if cursor.exhausted {
            cursor.has_head = false;
            return Ok(());
        }
        sink.block_mut().clear();
        let more = sample.generator.next_block(&sample.reads)?;
        // PANIC-FREE: as above.
        let sink = sample
            .generator
            .block_sink_mut()
            .expect("a block mode installs the sink");
        account_block(sink.block(), c);
        cursor.at = 0;
        cursor.exhausted = !more;
    }
}

fn mode_blk_records(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    sink: &mut Sink,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = samples.len();
    let mut cursors: Vec<BlockCursor> = (0..n).map(|_| BlockCursor::default()).collect();
    let mut scratch: Vec<LocusScratch> = (0..n).map(|_| LocusScratch::default()).collect();
    let mut mp = MergedPosition::default();

    for region in regions {
        for (si, s) in samples.iter_mut().enumerate() {
            s.generator.begin_segment(*region);
            // PANIC-FREE: the sink is installed for every block mode.
            s.generator
                .block_sink_mut()
                .expect("a block mode installs the sink")
                .block_mut()
                .clear();
            cursors[si] = BlockCursor::default();
        }
        for si in 0..n {
            blk_pull_scratch(&mut samples[si], &mut cursors[si], &mut scratch[si], c)?;
        }
        loop {
            let mut min_key = u64::MAX;
            for si in 0..n {
                if cursors[si].has_head {
                    min_key = min_key.min(scratch[si].key());
                }
            }
            if min_key == u64::MAX {
                break;
            }
            c.positions_seen += 1;

            let mut max_nonref = 0u32;
            for si in 0..n {
                if cursors[si].has_head && scratch[si].key() == min_key {
                    max_nonref = max_nonref.max(nonref_obs_of_scratch(&scratch[si]));
                }
            }
            let kept = max_nonref >= MIN_ALT_OBS;
            if kept {
                c.positions_kept += 1;
                mp.begin(min_key);
            }

            for si in 0..n {
                if !cursors[si].has_head || scratch[si].key() != min_key {
                    continue;
                }
                if kept {
                    c.merge_bytes += copy_locus_from_scratch(&mut mp, si as u32, &scratch[si]);
                    c.merge_objects += 1;
                }
                blk_pull_scratch(&mut samples[si], &mut cursors[si], &mut scratch[si], c)?;
            }
            if kept {
                sink.accept(&mp);
            }
        }
    }
    c.scratch_bytes = scratch.iter().map(|s| s.bytes_copied).sum();
    Ok(())
}

/// **Arm E.** The record-shaped merge again, but allowed to read the block's summary
/// column for its keep rule, so a locus is only refilled into a scratch record when the
/// cohort has decided to keep its position.
///
/// It exists to separate two things arm D confounds: *folding a cheap column across
/// samples* and *the block carrying a cheap column at all*. If E lands on D, the fold's
/// dense window buys nothing and what pays is the summary column.
///
/// It also needs no round structure. It borrows one sample's block at a time, so no
/// moment exists where every sample's block must be held at once.
fn mode_blk_records_sum(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    sink: &mut Sink,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = samples.len();
    let mut cursors: Vec<BlockCursor> = (0..n).map(|_| BlockCursor::default()).collect();
    let mut scratch: Vec<LocusScratch> = (0..n).map(|_| LocusScratch::default()).collect();
    let mut head_key = vec![u64::MAX; n];
    let mut head_nonref = vec![0u32; n];
    let mut mp = MergedPosition::default();

    for region in regions {
        for (si, s) in samples.iter_mut().enumerate() {
            s.generator.begin_segment(*region);
            // PANIC-FREE: the sink is installed for every block mode.
            s.generator
                .block_sink_mut()
                .expect("a block mode installs the sink")
                .block_mut()
                .clear();
            cursors[si] = BlockCursor::default();
        }
        for si in 0..n {
            blk_seat_head(
                &mut samples[si],
                &mut cursors[si],
                &mut head_key[si],
                &mut head_nonref[si],
                c,
            )?;
        }
        loop {
            let mut min_key = u64::MAX;
            for &k in head_key.iter() {
                min_key = min_key.min(k);
            }
            if min_key == u64::MAX {
                break;
            }
            c.positions_seen += 1;

            let mut max_nonref = 0u32;
            for si in 0..n {
                if head_key[si] == min_key {
                    max_nonref = max_nonref.max(head_nonref[si]);
                }
            }
            let kept = max_nonref >= MIN_ALT_OBS;
            if kept {
                c.positions_kept += 1;
                mp.begin(min_key);
            }

            for si in 0..n {
                if head_key[si] != min_key {
                    continue;
                }
                if kept {
                    // Only now does a record exist.
                    let block = samples[si]
                        .generator
                        .block_sink()
                        .expect("a block mode installs the sink")
                        .block();
                    scratch[si].refill_from(block, cursors[si].at);
                    c.loci_materialised += 1;
                    c.merge_bytes += copy_locus_from_scratch(&mut mp, si as u32, &scratch[si]);
                    c.merge_objects += 1;
                }
                cursors[si].at += 1;
                blk_seat_head(
                    &mut samples[si],
                    &mut cursors[si],
                    &mut head_key[si],
                    &mut head_nonref[si],
                    c,
                )?;
            }
            if kept {
                sink.accept(&mp);
            }
        }
    }
    c.scratch_bytes = scratch.iter().map(|s| s.bytes_copied).sum();
    Ok(())
}

/// Seat arm E's head on the locus at the cursor, refilling the block if the cursor ran off
/// its end. Reads three per-locus arrays and nothing ragged.
fn blk_seat_head(
    sample: &mut NgSample,
    cursor: &mut BlockCursor,
    key: &mut u64,
    nonref: &mut u32,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    loop {
        // PANIC-FREE: the sink is installed for every block mode.
        let sink = sample
            .generator
            .block_sink_mut()
            .expect("a block mode installs the sink");
        if cursor.at < sink.block().len() {
            *key = key_in_block(sink.block(), cursor.at);
            *nonref = sink.block().nonref_obs_of(cursor.at);
            cursor.has_head = true;
            return Ok(());
        }
        if cursor.exhausted {
            *key = u64::MAX;
            cursor.has_head = false;
            return Ok(());
        }
        sink.block_mut().clear();
        let more = sample.generator.next_block(&sample.reads)?;
        // PANIC-FREE: as above.
        let sink = sample
            .generator
            .block_sink_mut()
            .expect("a block mode installs the sink");
        account_block(sink.block(), c);
        cursor.at = 0;
        cursor.exhausted = !more;
    }
}

fn mode_blk_fold(
    samples: &mut [NgSample],
    regions: &[GenomeRegion],
    sink: &mut Sink,
    c: &mut Counters,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = samples.len();
    let mut cursors: Vec<BlockCursor> = (0..n).map(|_| BlockCursor::default()).collect();
    let mut window_nonref: Vec<u32> = Vec::new();
    let mut window_key: Vec<u64> = Vec::new();
    let mut mat_cursor = vec![0usize; n];
    let mut mp = MergedPosition::default();

    for region in regions {
        for (si, s) in samples.iter_mut().enumerate() {
            s.generator.begin_segment(*region);
            // PANIC-FREE: the sink is installed for every block mode.
            s.generator
                .block_sink_mut()
                .expect("a block mode installs the sink")
                .block_mut()
                .clear();
            cursors[si] = BlockCursor::default();
        }
        loop {
            // ---- mutable phase: every sample that has run off its block gets a new one.
            for si in 0..n {
                loop {
                    // PANIC-FREE: the sink is installed for every block mode.
                    let sink_ref = samples[si]
                        .generator
                        .block_sink_mut()
                        .expect("a block mode installs the sink");
                    if cursors[si].exhausted || cursors[si].at < sink_ref.block().len() {
                        break;
                    }
                    sink_ref.block_mut().clear();
                    let more = samples[si].generator.next_block(&samples[si].reads)?;
                    // PANIC-FREE: as above.
                    let sink_ref = samples[si]
                        .generator
                        .block_sink_mut()
                        .expect("a block mode installs the sink");
                    account_block(sink_ref.block(), c);
                    cursors[si].at = 0;
                    cursors[si].exhausted = !more;
                }
            }

            // ---- shared phase: every sample's block is borrowed at once. N distinct
            // owners can all be borrowed shared, which is why no reference count and no
            // per-sample owned copy is needed. Nothing may advance until the views drop.
            let blocks: Vec<&LocusBlock> = samples
                .iter()
                .map(|s| {
                    // PANIC-FREE: the sink is installed for every block mode.
                    s.generator
                        .block_sink()
                        .expect("a block mode installs the sink")
                        .block()
                })
                .collect();

            let mut any = false;
            let mut watermark = u64::MAX;
            for si in 0..n {
                let len = blocks[si].len();
                if cursors[si].at < len {
                    any = true;
                }
                if !cursors[si].exhausted && len > 0 {
                    watermark = watermark.min(key_in_block(blocks[si], len - 1));
                }
            }
            if !any {
                break;
            }
            if watermark == u64::MAX {
                watermark = (0..n)
                    .filter(|&si| blocks[si].len() > 0)
                    .map(|si| key_in_block(blocks[si], blocks[si].len() - 1))
                    .max()
                    .unwrap_or(0);
            }

            let win_lo = (0..n)
                .filter(|&si| cursors[si].at < blocks[si].len())
                .map(|si| key_in_block(blocks[si], cursors[si].at))
                .min()
                .unwrap_or(watermark);
            let win_len = (watermark.saturating_sub(win_lo) + 1).min(WINDOW_CAP) as usize;
            let watermark = win_lo + win_len as u64 - 1;

            // The fold itself: three arrays per sample — contig, start, and the block's own
            // summary column. `obs_bases` is not touched, and neither is anything ragged.
            window_nonref.clear();
            window_nonref.resize(win_len, u32::MAX); // MAX = "no sample here"
            for si in 0..n {
                let block = blocks[si];
                let mut i = cursors[si].at;
                while i < block.len() {
                    let key = key_in_block(block, i);
                    if key > watermark {
                        break;
                    }
                    let off = (key - win_lo) as usize;
                    let v = block.nonref_obs_of(i);
                    let slot = &mut window_nonref[off];
                    *slot = if *slot == u32::MAX { v } else { (*slot).max(v) };
                    i += 1;
                }
            }
            window_key.clear();
            for (off, &v) in window_nonref.iter().enumerate() {
                if v != u32::MAX {
                    window_key.push(win_lo + off as u64);
                }
            }
            c.positions_seen += window_key.len() as u64;

            // Only the variable positions ever read the heavy arrays.
            for si in 0..n {
                mat_cursor[si] = cursors[si].at;
            }
            for &key in &window_key {
                if window_nonref[(key - win_lo) as usize] < MIN_ALT_OBS {
                    continue;
                }
                c.positions_kept += 1;
                mp.begin(key);
                for si in 0..n {
                    let block = blocks[si];
                    let mut i = mat_cursor[si];
                    while i < block.len() && key_in_block(block, i) < key {
                        i += 1;
                    }
                    mat_cursor[si] = i;
                    if i >= block.len() || key_in_block(block, i) != key {
                        continue;
                    }
                    c.merge_bytes += copy_locus_from_block(&mut mp, si as u32, block, i);
                    c.merge_objects += 1;
                }
                sink.accept(&mp);
            }

            // Advance the cursors past the window, then let the views drop.
            for si in 0..n {
                let block = blocks[si];
                let mut i = cursors[si].at;
                while i < block.len() && key_in_block(block, i) <= watermark {
                    i += 1;
                }
                cursors[si].at = i;
            }
            drop(blocks);
        }
    }
    Ok(())
}

// =====================================================================
// Driver
// =====================================================================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 6 {
        eprintln!(
            "usage: sketch4_producer_merge \
             <rec-drop|rec-records|rec-fold|blk-drop|blk-drop-sum|blk-records|blk-fold> \
             <full|digest> <ref.fa> <regions.bed> <n_regions> <a.cram> ..."
        );
        std::process::exit(2);
    }
    let mode = args[0].clone();
    let consumer = match args[1].as_str() {
        "full" => Consumer::Full,
        "digest" => Consumer::Digest,
        other => panic!("unknown consumer {other}"),
    };
    let fasta = PathBuf::from(&args[2]);
    let bed = PathBuf::from(&args[3]);
    let n_regions: usize = args[4].parse()?;
    let crams: Vec<PathBuf> = args[5..].iter().map(PathBuf::from).collect();
    let block_bytes: usize = std::env::var("PVC_SKETCH4_BLOCK_KB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .map_or(BLOCK_BYTES_DEFAULT, |kb| kb * 1024);

    let uses_block = mode.starts_with("blk-");
    let wants_summary = matches!(
        mode.as_str(),
        "blk-drop-sum" | "blk-fold" | "blk-records-sum"
    );

    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.clone(),
        ReferenceCheck::TrustIndexWithoutChecking,
    )?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&fasta)?;
    if let Some(v) = verify {
        v.join().expect("reference verification thread");
    }

    let mut regions: Vec<GenomeRegion> = Vec::new();
    for line in BufReader::new(File::open(&bed)?).lines() {
        let line = line?;
        let mut f = line.split('\t');
        let (Some(name), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else {
            continue;
        };
        let Some(cid) = contigs.entries.iter().position(|c| c.name == name) else {
            continue;
        };
        regions.push(GenomeRegion {
            contig: ContigId(cid as u32),
            start: Position(s.parse::<u64>()? + 1),
            end: Position(e.parse::<u64>()?),
        });
        if regions.len() >= n_regions {
            break;
        }
    }

    let reference = OpenReference::new(info);
    let mut samples: Vec<NgSample> = Vec::with_capacity(crams.len());
    for cram in &crams {
        let reads = SampleReads::open_only_sample(
            std::slice::from_ref(cram),
            &reference,
            ReadFilterConfig::default(),
            true,
        )?;
        let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
            fasta.clone(),
            contigs.clone(),
            index.clone(),
        ));
        #[allow(
            clippy::arc_with_non_send_sync,
            reason = "the file-backed accessor is single-threaded here, as in benches/ng_generic_pileup_perf.rs"
        )]
        let refseq = Arc::new(WindowedRefSeq::with_shared_index(
            fasta.clone(),
            contigs.clone(),
            index.clone(),
        ));
        let make_reference = {
            let fasta = fasta.clone();
            let contigs = contigs.clone();
            let index = index.clone();
            move || WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
        };
        let mut generator = PileupGenerator::new(
            refseq,
            make_reference,
            preparer,
            PileupGeneratorConfig::default(),
        )?;
        if uses_block {
            generator.install_block_sink(LocusSink::columns_with_summary(
                block_bytes,
                wants_summary,
            ));
        }
        samples.push(NgSample {
            generator,
            reads,
            head: None,
        });
    }

    let mut sink = Sink::new(consumer);
    let mut c = Counters::default();
    match mode.as_str() {
        "stats" => mode_stats(&mut samples, &regions, &mut c)?,
        "rec-drop" => mode_rec_drop(&mut samples, &regions, &mut c)?,
        "rec-records" => mode_rec_records(&mut samples, &regions, &mut sink, &mut c)?,
        "rec-fold" => mode_rec_fold(&mut samples, &regions, &mut sink, &mut c)?,
        "blk-drop" | "blk-drop-sum" => mode_blk_drop(&mut samples, &regions, &mut c)?,
        "blk-records" => mode_blk_records(&mut samples, &regions, &mut sink, &mut c)?,
        "blk-records-sum" => mode_blk_records_sum(&mut samples, &regions, &mut sink, &mut c)?,
        "blk-fold" => mode_blk_fold(&mut samples, &regions, &mut sink, &mut c)?,
        other => panic!("unknown mode {other}"),
    }

    println!("mode={mode} consumer={}", args[1]);
    println!("samples={} regions={}", crams.len(), regions.len());
    println!("block_bytes={block_bytes} summary={wants_summary}");
    println!("digest=0x{:016x}", sink.digest);
    println!("loglik_acc={:.6}", sink.loglik_acc);
    println!("loci_produced={}", c.loci_produced);
    println!("loci_materialised={}", c.loci_materialised);
    println!("positions_seen={}", c.positions_seen);
    println!("positions_kept={}", c.positions_kept);
    println!("producer_bytes={}", c.producer_bytes);
    println!("block_bytes_written={}", c.block_bytes);
    println!("scratch_bytes={}", c.scratch_bytes);
    println!("observations={}", c.observations);
    println!("visits={}", c.visits);
    println!("blocks={}", c.blocks);
    println!(
        "block_payload_high_water={}",
        c.block_payload_high_water
    );
    println!("merge_objects={}", c.merge_objects);
    println!("merge_bytes={}", c.merge_bytes);
    Ok(())
}
