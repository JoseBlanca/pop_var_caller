//! **Refinement's memory and output on a drawn cohort** — what the window budget counts, what the
//! allocator holds, and a hash of the lines written, so two builds can be compared byte for byte.
//!
//! ```text
//! refine_mem_probe <samples> <contig_kb> <records_per_kb> <depth> <core_kb> <margin_kb> [density]
//! ```
//!
//! Draws a multi-allelic cohort (2 to 6 alleles in tomato25's shares: 57, 20, 12, 6 and 5 records
//! in 100), each sample copying stretches of 24 founders, 90 in 100 samples fully inbred, with reads
//! as 300-base fragments at `depth` fragments covering a base. Pushes every line and its read
//! alleles through `NeighbourEvidenceStage` with windows of a fixed `core_kb` and `margin_kb`, and
//! margins thinned to `density` records a kb if given. After each window it prints the bytes the
//! budget counted for it; at the end, the hash of every line written, in order.
//!
//! Environment:
//!
//! - `PROBE_BUDGET_MB` and `PROBE_FLOOR_KB`: a window budget, and the margin and core floor it may
//!   shrink a window to; without them the budget never binds, so the windows do not depend on the
//!   accounting being measured.
//! - `PROBE_OUT`: a file the written lines go to as well.
//! - `PROBE_STATS=1`: instead of running the stage, the average bytes one record holds a sample —
//!   its line, its numbers, its read alleles, packed and in byte planes.
//!
//! Built with `--no-default-features --features dhat-heap`, dhat counts the bytes asked for and
//! prints the peak at the end; with the default mimalloc, read the peak resident memory from
//! `/usr/bin/time -l` (macOS) or `-v` (Linux). `RAYON_NUM_THREADS` sets the threads.
#![allow(
    clippy::disallowed_methods,
    reason = "a research probe outside the portable-float guarantee (clippy.toml)"
)]

use std::hash::{Hash, Hasher};
use std::time::Instant;

use pop_var_caller::neighbour_evidence::line::site_of;
use pop_var_caller::neighbour_evidence::stage::{
    NeighbourEvidenceStage, NeighbourEvidenceSwitch, StageStop, decide,
};
use pop_var_caller::neighbour_evidence::windows::MarginThinning;
use pop_var_caller::read_phase::fragments::{
    FragmentAllele, PackedRecordFragments, RecordFragments,
};
use pop_var_caller::types::{ContigId, GenomePosition, Ploidy, Position};
use pop_var_caller::vcf::RecordPlace;

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

const MB: f64 = 1e6;

struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut mixed = self.0;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        mixed ^ (mixed >> 31)
    }
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }
    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }
    fn poisson(&mut self, mean: f64) -> usize {
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit {
            product *= self.uniform();
            count += 1;
        }
        count
    }
    fn pick(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut target = self.uniform() * total;
        for (index, &weight) in weights.iter().enumerate() {
            if target < weight {
                return index;
            }
            target -= weight;
        }
        weights.len() - 1
    }
}

const FOUNDERS: usize = 24;
const FRAGMENT_BASES: u64 = 300;
const SEGMENT_BASES: f64 = 200_000.0;
const INBREEDING: f64 = 0.9;
const READ_ERROR: f64 = 0.01;
const KEPT: usize = 32;
const ALTS: [&str; 5] = ["C", "G", "T", "AC", "AG"];

#[derive(Default)]
struct Stats {
    records: usize,
    line_plain: usize,
    line_compressed: usize,
    record_heap: usize,
    reads_heap: usize,
    read_alleles: usize,
    packed: usize,
    packed_numbers: usize,
    packed_planes: usize,
    offsets_fit_a_byte: usize,
    phred_fits_a_byte: usize,
}

struct Fragment {
    id: u32,
    end: u64,
    copy: usize,
}

struct SampleState {
    founders: [usize; 2],
    inbred: bool,
    fragments: Vec<Fragment>,
    next_fragment: u32,
}

fn genotype_index(a: usize, b: usize) -> usize {
    let (a, b) = (a.min(b), a.max(b));
    b * (b + 1) / 2 + a
}

#[allow(clippy::too_many_lines)]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 7 {
        eprintln!(
            "usage: refine_mem_probe <samples> <contig_kb> <records_per_kb> <depth> <core_kb> <margin_kb> [density]"
        );
        std::process::exit(2);
    }
    let samples: usize = args[1].parse().expect("samples");
    let contig_bases: u64 = args[2].parse::<u64>().expect("contig kb") * 1_000;
    let per_kb: f64 = args[3].parse().expect("records a kb");
    let depth: f64 = args[4].parse().expect("depth");
    let core: u64 = args[5].parse::<u64>().expect("core kb") * 1_000;
    let margin: u64 = args[6].parse::<u64>().expect("margin kb") * 1_000;
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();
    let spacing = (1_000.0 / per_kb).round().max(1.0) as u64;

    let mut decision = decide(
        NeighbourEvidenceSwitch::On,
        &vec![INBREEDING; samples],
        Ploidy::try_new(2).expect("two"),
    );
    let env_number = |name: &str| {
        std::env::var(name).ok().map(|text| {
            text.parse::<u64>()
                .unwrap_or_else(|_| panic!("{name} is a number"))
        })
    };
    let floor = env_number("PROBE_FLOOR_KB").map_or(margin, |kb| kb * 1_000);
    decision.sizes.core_bases = core;
    decision.sizes.margin_bases = margin;
    decision.sizes.margin_floor_bases = floor.min(margin);
    decision.sizes.core_floor_bases = floor.min(core).min(margin);
    decision.sizes.budget_bytes =
        env_number("PROBE_BUDGET_MB").map_or(u64::MAX, |mb| mb * 1_000_000);
    if let Some(density) = args.get(7) {
        decision = decision.with_margin_thinning(MarginThinning::parse(density).expect("density"));
    }
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![contig_bases + 1_000]);

    let mut draw = Draw(0x5eed);
    let mut states: Vec<SampleState> = (0..samples)
        .map(|_| {
            let inbred = draw.uniform() < INBREEDING;
            let first = draw.below(FOUNDERS);
            SampleState {
                founders: [first, if inbred { first } else { draw.below(FOUNDERS) }],
                inbred,
                fragments: Vec::new(),
                next_fragment: 0,
            }
        })
        .collect();
    let switch = 1.0 - (-(spacing as f64) / SEGMENT_BASES).exp();

    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    let (mut written, mut lines_in_window, mut windows) = (0_u64, 0_u64, 0_usize);
    let mut most_counted = 0_u64;
    let started = Instant::now();
    let mut window_started = Instant::now();
    // Prints the last window planned and returns the bytes the budget counted for it.
    let report = |windows: usize, stage: &NeighbourEvidenceStage, lines: u64, seconds: f64| {
        let Some(plan) = stage.last_window() else {
            return 0;
        };
        eprintln!(
            "window {windows}: {lines} lines; counted {:.1} MB ({:.1} MB outside it) over {} \
             records ({} core, core {} kb, margins {} kb); {seconds:.1} s",
            plan.bytes as f64 / MB,
            plan.outside_bytes as f64 / MB,
            plan.records(),
            plan.core_end - plan.core_start,
            plan.core_bases / 1_000,
            plan.margin_bases / 1_000,
        );
        plan.bytes
    };

    let stats_only = std::env::var("PROBE_STATS").is_ok();
    let mut stats = Stats::default();
    let mut out = std::env::var("PROBE_OUT")
        .ok()
        .map(|path| std::io::BufWriter::new(std::fs::File::create(path).expect("PROBE_OUT")));
    let mut rows_over_255 = 0_usize;
    let mut position = 1_u64;
    let mut alleles_seen = [0_usize; 7];
    while position <= contig_bases {
        // The record: its alleles and their frequencies.
        let allele_count = 2 + draw.pick(&[57.0, 20.0, 12.0, 6.0, 5.0]);
        alleles_seen[allele_count] += 1;
        let reference = 0.3 + 0.6 * draw.uniform();
        let weights: Vec<f64> = (1..allele_count).map(|_| -draw.uniform().ln()).collect();
        let total: f64 = weights.iter().sum();
        let mut frequencies = vec![reference];
        frequencies.extend(
            weights
                .iter()
                .map(|weight| (1.0 - reference) * weight / total),
        );
        let founder_alleles: Vec<usize> = (0..FOUNDERS).map(|_| draw.pick(&frequencies)).collect();
        let genotypes = allele_count * (allele_count + 1) / 2;

        let mut columns = Vec::with_capacity(samples);
        let mut record_max_phred = 0_u32;
        let mut entries = Vec::new();
        for (sample, state) in states.iter_mut().enumerate() {
            for copy in 0..2 {
                if draw.uniform() < switch {
                    state.founders[copy] = draw.below(FOUNDERS);
                }
            }
            if state.inbred {
                state.founders[1] = state.founders[0];
            }
            state.fragments.retain(|fragment| fragment.end >= position);
            for _ in 0..draw.poisson(depth * spacing as f64 / FRAGMENT_BASES as f64) {
                state.fragments.push(Fragment {
                    id: state.next_fragment,
                    end: position + FRAGMENT_BASES,
                    copy: draw.below(2),
                });
                state.next_fragment += 1;
            }
            let mut reads = Vec::with_capacity(state.fragments.len());
            for fragment in &state.fragments {
                let mut allele = founder_alleles[state.founders[fragment.copy]];
                if draw.uniform() < READ_ERROR {
                    allele = (allele + 1 + draw.below(allele_count - 1)) % allele_count;
                }
                let phred = 25 + draw.below(16) as u8;
                reads.push((fragment.id, allele, phred));
            }
            for &(fragment, allele, error_phred) in reads.iter().take(KEPT) {
                entries.push(FragmentAllele {
                    sample: sample as u32,
                    fragment,
                    allele: allele as u16,
                    error_phred,
                });
            }
            // PL from every read; GP from PL and the inbred prior, three decimals.
            let mut log_likelihoods = vec![0.0_f64; genotypes];
            for b in 0..allele_count {
                for a in 0..=b {
                    let mut total = 0.0;
                    for &(_, allele, _) in &reads {
                        let given = |x: usize| {
                            if allele == x {
                                1.0 - READ_ERROR
                            } else {
                                READ_ERROR / (allele_count - 1) as f64
                            }
                        };
                        total += (0.5 * given(a) + 0.5 * given(b)).ln();
                    }
                    log_likelihoods[genotype_index(a, b)] = total;
                }
            }
            let best = log_likelihoods.iter().copied().fold(f64::MIN, f64::max);
            let phred: Vec<u32> = log_likelihoods
                .iter()
                .map(|&value| (-10.0 * (value - best) / std::f64::consts::LN_10).round() as u32)
                .collect();
            record_max_phred = record_max_phred.max(phred.iter().copied().max().unwrap_or(0));
            rows_over_255 += usize::from(phred.iter().any(|&value| value > 255));
            let mut posterior = vec![0.0_f64; genotypes];
            for b in 0..allele_count {
                for a in 0..=b {
                    let prior = if a == b {
                        INBREEDING * frequencies[a]
                            + (1.0 - INBREEDING) * frequencies[a] * frequencies[a]
                    } else {
                        2.0 * (1.0 - INBREEDING) * frequencies[a] * frequencies[b]
                    };
                    let index = genotype_index(a, b);
                    posterior[index] = 10_f64.powf(-f64::from(phred[index]) / 10.0) * prior;
                }
            }
            let total: f64 = posterior.iter().sum();
            let posterior: Vec<f64> = posterior
                .iter()
                .map(|value| (value / total * 1000.0).round() / 1000.0)
                .collect();
            let call = (0..genotypes).fold(0, |best, index| {
                if posterior[index] > posterior[best] {
                    index
                } else {
                    best
                }
            });
            let genotype = if reads.is_empty() {
                "./.".to_string()
            } else {
                let mut found = (0, 0);
                for b in 0..allele_count {
                    for a in 0..=b {
                        if genotype_index(a, b) == call {
                            found = (a, b);
                        }
                    }
                }
                format!("{}/{}", found.0, found.1)
            };
            let join = |values: Vec<String>| values.join(",");
            columns.push(format!(
                "{genotype}:{}:{}",
                join(phred.iter().map(ToString::to_string).collect()),
                join(posterior.iter().map(ToString::to_string).collect())
            ));
        }
        let line = format!(
            "chr1\t{position}\t.\tA\t{}\t50\tPASS\tAF={};AC=0;AN=0\tGT:PL:GP\t{}",
            ALTS[..allele_count - 1].join(","),
            frequencies[1..]
                .iter()
                .map(|frequency| format!("{frequency:.4}"))
                .collect::<Vec<_>>()
                .join(","),
            columns.join("\t")
        );
        let place = RecordPlace {
            at: GenomePosition {
                contig: ContigId(0),
                position: Position(position),
            },
            is_repeat_tract: false,
        };
        if stats_only {
            // What one record holds, a sample at a time: its line compressed alone at level 1, its
            // numbers and reads, packed as a skipped far-margin record is, and its read alleles.
            let fragments = RecordFragments::from_entries(entries);
            let mut max_offset = 0;
            let mut first: Option<(u32, u32)> = None;
            for entry in fragments.entries() {
                match first {
                    Some((sample, fragment)) if sample == entry.sample => {
                        max_offset = max_offset.max(entry.fragment - fragment);
                    }
                    _ => first = Some((entry.sample, entry.fragment)),
                }
            }
            stats.read_alleles += fragments.len();
            stats.offsets_fit_a_byte += usize::from(max_offset <= 255);
            let packed = PackedRecordFragments::from(fragments);
            stats.reads_heap += packed.heap_bytes();
            let mut record = site_of(line.as_bytes(), samples).expect("a site");
            stats.record_heap += record.heap_bytes();
            let numbers = record.take_numbers();
            let mut plain = Vec::new();
            numbers.write_to(&mut plain);
            let numbers_alone = zstd::bulk::compress(&plain, 1).expect("zstd").len();
            packed.write_to(&mut plain);
            stats.packed += zstd::bulk::compress(&plain, 1).expect("zstd").len();
            stats.packed_numbers += numbers_alone;
            // The same bytes, each array's values split into byte planes (the k-th byte of every
            // value together), which zstd finds more alike.
            stats.packed_planes += zstd::bulk::compress(&planes(&plain), 1)
                .expect("zstd")
                .len();
            stats.line_plain += line.len();
            stats.line_compressed += zstd::bulk::compress(line.as_bytes(), 1)
                .expect("zstd")
                .len();
            stats.records += 1;
            stats.phred_fits_a_byte += usize::from(record_max_phred <= 255);
            position += spacing;
            continue;
        }
        let outcome = stage.push_writing(
            place,
            line.into_bytes(),
            RecordFragments::from_entries(entries),
            &mut |_, line: Vec<u8>| -> Result<(), std::convert::Infallible> {
                line.hash(&mut hasher);
                if let Some(out) = out.as_mut() {
                    use std::io::Write;
                    out.write_all(&line).expect("write");
                    out.write_all(b"\n").expect("write");
                }
                lines_in_window += 1;
                Ok(())
            },
        );
        match outcome {
            Ok(()) => {}
            Err(StageStop::OverBudget(over)) => panic!("over budget: {over:?}"),
            Err(StageStop::Write(never)) => match never {},
            Err(StageStop::Stopped) => unreachable!("the probe stops at the first stop"),
        }
        if lines_in_window > 0 {
            windows += 1;
            most_counted = most_counted.max(report(
                windows,
                &stage,
                lines_in_window,
                window_started.elapsed().as_secs_f64(),
            ));
            written += lines_in_window;
            lines_in_window = 0;
            window_started = Instant::now();
        }
        position += spacing;
    }
    if stats_only {
        let per = |bytes: usize| bytes as f64 / (stats.records * samples) as f64;
        eprintln!(
            "bytes a sample a record over {} records: line {:.2} plain, {:.2} compressed; record heap {:.2}; read alleles {:.2} heap, {:.2} entries; packed numbers+reads {:.2} (numbers alone {:.2}; in byte planes {:.2}); records whose offsets fit a byte {:.4}; records whose PL fit a byte {:.4}; sample rows over 255 {:.6}",
            stats.records,
            per(stats.line_plain),
            per(stats.line_compressed),
            per(stats.record_heap),
            per(stats.reads_heap),
            per(stats.read_alleles),
            per(stats.packed),
            per(stats.packed_numbers),
            per(stats.packed_planes),
            stats.offsets_fit_a_byte as f64 / stats.records as f64,
            stats.phred_fits_a_byte as f64 / stats.records as f64,
            rows_over_255 as f64 / (stats.records * samples) as f64
        );
        return;
    }
    match stage.finish_writing(
        &mut |_, line: Vec<u8>| -> Result<(), std::convert::Infallible> {
            line.hash(&mut hasher);
            if let Some(out) = out.as_mut() {
                use std::io::Write;
                out.write_all(&line).expect("write");
                out.write_all(b"\n").expect("write");
            }
            lines_in_window += 1;
            Ok(())
        },
    ) {
        Ok(()) => {}
        Err(StageStop::OverBudget(over)) => panic!("over budget: {over:?}"),
        Err(StageStop::Write(never)) => match never {},
        Err(StageStop::Stopped) => unreachable!("the probe stops at the first stop"),
    }
    most_counted = most_counted.max(report(
        windows + 1,
        &stage,
        lines_in_window,
        window_started.elapsed().as_secs_f64(),
    ));
    written += lines_in_window;
    #[cfg(feature = "dhat-heap")]
    eprintln!(
        "dhat: peak {:.1} MB requested, {:.1} GB asked for in all",
        dhat::HeapStats::get().max_bytes as f64 / MB,
        dhat::HeapStats::get().total_bytes as f64 / 1e9
    );
    eprintln!(
        "alleles: {:?}; {written} lines; most counted {:.1} MB; {:.1} s; output hash {:016x}",
        &alleles_seen[2..],
        most_counted as f64 / MB,
        started.elapsed().as_secs_f64(),
        hasher.finish()
    );
}

/// **A packed record's bytes with every array split into byte planes**: what `RecordNumbers` and
/// then `PackedRecordFragments` write (three length-prefixed arrays of `f32`, `u16`, `f32`; then the
/// sample table of three `u32`s a sample, a form byte, a count and three-byte entries), each
/// array's k-th bytes together. The same bytes, reordered.
fn planes(plain: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(plain.len());
    let mut at = 0;
    let u32_at = |at: &mut usize| {
        let value = u32::from_le_bytes(plain[*at..*at + 4].try_into().expect("4")) as usize;
        *at += 4;
        value
    };
    let split = |out: &mut Vec<u8>, bytes: &[u8], width: usize| {
        for plane in 0..width {
            out.extend(bytes.iter().skip(plane).step_by(width));
        }
    };
    for width in [4, 2, 4] {
        let count = u32_at(&mut at);
        split(&mut out, &plain[at..at + count * width], width);
        at += count * width;
    }
    let samples = u32_at(&mut at);
    split(&mut out, &plain[at..at + 12 * samples], 12);
    at += 12 * samples;
    let form = plain[at];
    at += 1;
    let count = u32_at(&mut at);
    let width = [3, 4, 7][usize::from(form)];
    split(&mut out, &plain[at..at + count * width], width);
    out
}
