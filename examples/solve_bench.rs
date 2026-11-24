//! **Times refinement's solve on a drawn cohort** — a measurement tool for the speed review of
//! 2026-11-24 (`doc/devel/reports/reviews/refine_solve_speed_opportunities_2026-11-24.md`).
//!
//! ```text
//! solve_bench <samples> <records> <spacing> <depth> <window_kb> <seed> [out.vcf]
//! ```
//!
//! Draws a cohort of mostly inbred samples whose copies are long stretches of a few founder
//! haplotypes, draws paired fragments of 350 bases at `depth` reads a position, writes the `PL`,
//! `GP` and `AF` calling would from those reads, and pushes every line with its read alleles
//! through `NeighbourEvidenceStage` at windows of at most `window_kb` and the default margin
//! thinning (5 records a kb). The output lines go to `out.vcf` when given; the wall time of the
//! drawing and of the stage go to stderr. Threads: `RAYON_NUM_THREADS`. Margin thinning:
//! `DENSITY`, as `--refine-margin-density` takes it (`0` walks every record).
//!
//! The allocator is mimalloc, as `pvc`'s. Built with `--no-default-features --features dhat-heap`,
//! it is dhat, and the peak of the bytes asked for is printed at the end.
#![allow(
    clippy::disallowed_methods,
    reason = "a measurement tool outside the portable-float guarantee (clippy.toml)"
)]

use std::io::Write;
use std::time::Instant;

use pop_var_caller::neighbour_evidence::stage::{
    NeighbourEvidenceStage, NeighbourEvidenceSwitch, decide,
};
use pop_var_caller::neighbour_evidence::windows::MarginThinning;
use pop_var_caller::read_phase::fragments::{FragmentAllele, RecordFragments};
use pop_var_caller::types::{ContigId, GenomePosition, Ploidy, Position};
use pop_var_caller::vcf::RecordPlace;

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
    fn exponential(&mut self, mean: f64) -> f64 {
        -mean * (1.0 - self.uniform()).ln()
    }
}

const FOUNDERS: usize = 12;
const SEGMENT_BASES: f64 = 150_000.0;
const FRAGMENT_BASES: u64 = 350;
const READ_ERROR: f64 = 0.01;
const ERROR_PHRED: u8 = 20;
const MUTATION: f64 = 0.003;

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 7 {
        eprintln!(
            "usage: solve_bench <samples> <records> <spacing> <depth> <window_kb> <seed> [out.vcf]"
        );
        std::process::exit(2);
    }
    let samples: usize = args[1].parse().unwrap();
    let records: usize = args[2].parse().unwrap();
    let spacing: u64 = args[3].parse().unwrap();
    let depth: f64 = args[4].parse().unwrap();
    let window_bases: u64 = args[5].parse::<u64>().unwrap() * 1_000;
    let seed: u64 = args[6].parse().unwrap();
    let out_path = args.get(7);

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();
    print_windows();
    let started = Instant::now();
    let mut draw = Draw(seed);
    let positions: Vec<u64> = (0..records as u64).map(|i| 1 + i * spacing).collect();
    let frequencies: Vec<f64> = (0..records).map(|_| 0.05 + 0.45 * draw.uniform()).collect();
    let founder_alleles: Vec<Vec<u8>> = (0..FOUNDERS)
        .map(|_| {
            frequencies
                .iter()
                .map(|&f| u8::from(draw.uniform() < f))
                .collect()
        })
        .collect();
    // Each sample: inbred (one mosaic for both copies) with chance 0.85, else two mosaics.
    let inbreeding: Vec<f64> = (0..samples)
        .map(|_| if draw.uniform() < 0.85 { 0.95 } else { 0.3 })
        .collect();
    let switch = 1.0 - (-(spacing as f64) / SEGMENT_BASES).exp();
    // copies[sample][copy][record] -> allele
    let mut copies: Vec<[Vec<u8>; 2]> = Vec::with_capacity(samples);
    for &sample_inbreeding in &inbreeding {
        let mosaic = |draw: &mut Draw| {
            let mut founder = draw.below(FOUNDERS);
            (0..records)
                .map(|record| {
                    if record > 0 && draw.uniform() < switch {
                        founder = draw.below(FOUNDERS);
                    }
                    let allele = founder_alleles[founder][record];
                    if draw.uniform() < MUTATION {
                        1 - allele
                    } else {
                        allele
                    }
                })
                .collect::<Vec<u8>>()
        };
        let first = mosaic(&mut draw);
        let second = if sample_inbreeding > 0.5 {
            first.clone()
        } else {
            mosaic(&mut draw)
        };
        copies.push([first, second]);
    }
    // Fragments: starts a Poisson process at depth / FRAGMENT_BASES a base.
    let length = positions[records - 1] + FRAGMENT_BASES;
    let mut fragments_at: Vec<Vec<FragmentAllele>> = vec![Vec::new(); records];
    let mut alt_reads = vec![vec![0u16; samples]; records];
    let mut all_reads = vec![vec![0u16; samples]; records];
    for sample in 0..samples {
        let mut start = draw.exponential(FRAGMENT_BASES as f64 / depth);
        let mut fragment = 0u32;
        while (start as u64) < length {
            let begin = start as u64;
            let copy = draw.below(2);
            let first = positions.partition_point(|&p| p < begin);
            let end = positions.partition_point(|&p| p < begin + FRAGMENT_BASES);
            for record in first..end {
                let mut allele = copies[sample][copy][record];
                if draw.uniform() < READ_ERROR {
                    allele = 1 - allele;
                }
                all_reads[record][sample] += 1;
                alt_reads[record][sample] += u16::from(allele);
                if end - first > 1 {
                    fragments_at[record].push(FragmentAllele {
                        sample: sample as u32,
                        fragment,
                        allele: u16::from(allele),
                        error_phred: ERROR_PHRED,
                    });
                }
            }
            fragment += 1;
            start += draw.exponential(FRAGMENT_BASES as f64 / depth);
        }
    }
    let mut lines = Vec::with_capacity(records);
    for record in 0..records {
        let frequency = frequencies[record];
        let mut columns = Vec::with_capacity(samples);
        for sample in 0..samples {
            let (reads, alternative) = (
                all_reads[record][sample] as f64,
                alt_reads[record][sample] as f64,
            );
            let log_likelihoods: Vec<f64> = [0.0, 0.5, 1.0]
                .iter()
                .map(|&share: &f64| {
                    let p = share * (1.0 - READ_ERROR) + (1.0 - share) * READ_ERROR;
                    alternative * p.ln() + (reads - alternative) * (1.0 - p).ln()
                })
                .collect();
            let best = log_likelihoods.iter().copied().fold(f64::MIN, f64::max);
            let phred: Vec<u32> = log_likelihoods
                .iter()
                .map(|&l| (-10.0 * (l - best) / std::f64::consts::LN_10).round() as u32)
                .collect();
            let f = inbreeding[sample];
            let prior = [
                (1.0 - f) * (1.0 - frequency) * (1.0 - frequency) + f * (1.0 - frequency),
                (1.0 - f) * 2.0 * frequency * (1.0 - frequency),
                (1.0 - f) * frequency * frequency + f * frequency,
            ];
            let unscaled: Vec<f64> = phred
                .iter()
                .zip(prior)
                .map(|(&p, prior)| 10f64.powf(-f64::from(p) / 10.0) * prior)
                .collect();
            let total: f64 = unscaled.iter().sum();
            let posteriors: Vec<f64> = unscaled
                .iter()
                .map(|v| (v / total * 1000.0).round() / 1000.0)
                .collect();
            let call = (0..3).fold(0, |b, g| if posteriors[g] > posteriors[b] { g } else { b });
            let genotype = if reads == 0.0 {
                "./."
            } else {
                ["0/0", "0/1", "1/1"][call]
            };
            columns.push(format!(
                "{genotype}:{},{},{}:{},{},{}",
                phred[0], phred[1], phred[2], posteriors[0], posteriors[1], posteriors[2]
            ));
        }
        lines.push(format!(
            "chr1\t{}\t.\tA\tC\t50\tPASS\tAF={frequency:.4};AC=0;AN=0\tGT:PL:GP\t{}",
            positions[record],
            columns.join("\t")
        ));
    }
    let read_alleles: usize = fragments_at.iter().map(Vec::len).sum();
    eprintln!(
        "drawn in {:.1} s: {samples} samples, {records} records, {:.2} read alleles in fragments a sample a record",
        started.elapsed().as_secs_f64(),
        read_alleles as f64 / (samples * records) as f64
    );

    let decision = decide(
        NeighbourEvidenceSwitch::On,
        &inbreeding,
        Ploidy::try_new(2).unwrap(),
    )
    .with_window_max(window_bases)
    .with_window_floor(window_bases.min(500_000))
    .with_margin_thinning(
        MarginThinning::parse(&std::env::var("DENSITY").unwrap_or("5".into())).unwrap(),
    );
    let mut stage = NeighbourEvidenceStage::new(&decision, vec![length + 1]);
    let solving = Instant::now();
    let mut out: Vec<(RecordPlace, Vec<u8>)> = Vec::new();
    if let Ok(path) = std::env::var("DUMP_IN") {
        std::fs::write(path, lines.join("\n")).unwrap();
    }
    for (record, (line, fragments)) in lines.into_iter().zip(fragments_at).enumerate() {
        let place = RecordPlace {
            at: GenomePosition {
                contig: ContigId(0),
                position: Position(positions[record]),
            },
            is_repeat_tract: false,
        };
        stage
            .push(
                place,
                line.into_bytes(),
                RecordFragments::from_entries(fragments),
                &mut out,
            )
            .unwrap();
    }
    stage.finish(&mut out).unwrap();
    eprintln!(
        "stage: {:.2} s, {} lines out",
        solving.elapsed().as_secs_f64(),
        out.len()
    );
    #[cfg(feature = "dhat-heap")]
    eprintln!(
        "dhat: peak {:.1} MB requested",
        dhat::HeapStats::get().max_bytes as f64 / 1e6
    );
    if let Some(path) = out_path {
        let mut file = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        for (_, line) in &out {
            file.write_all(line).unwrap();
            file.write_all(b"\n").unwrap();
        }
    }
}

/// Prints every `tracing` event's message to stderr: the stage's one line a window.
struct PrintEvents;

struct Message(String);

impl tracing::field::Visit for Message {
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        if field.name() == "message" {
            self.0 = format!("{value:?}");
        }
    }
}

impl<S: tracing::Subscriber> tracing_subscriber::Layer<S> for PrintEvents {
    fn on_event(&self, event: &tracing::Event<'_>, _: tracing_subscriber::layer::Context<'_, S>) {
        let mut message = Message(String::new());
        event.record(&mut message);
        if message.0.starts_with("refinement window") {
            eprintln!("{}", message.0);
        }
    }
}

fn print_windows() {
    use tracing_subscriber::layer::SubscriberExt;
    let subscriber = tracing_subscriber::registry().with(PrintEvents);
    tracing::subscriber::set_global_default(subscriber).unwrap();
}
