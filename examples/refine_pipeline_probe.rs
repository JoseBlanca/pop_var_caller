//! **The refinement pipeline's wall time and output on a drawn cohort, written through the real
//! `VcfWriter`** — the probe behind
//! `doc/devel/reports/implementations/refine_pipeline_build_2026-11-24.md`.
//!
//! ```text
//! refine_pipeline_probe gen <out.txt> <samples> <records> <spacing> <depth> <seed>
//! refine_pipeline_probe run <in.txt> <samples> <threads> <density> <out.vcf[.gz]>
//! ```
//!
//! `gen` draws an inbred cohort (eight founders, 300-kb stretches, inbreeding 0.95) with 150-base
//! reads, `depth` a sample a record on average, and writes each record as two text lines: the VCF
//! line calling would write (`GT:GQ:DP:AD:PL:GP`) and its fragments (`sample,fragment,allele,phred`
//! separated by `;`), at most 32 read alleles a sample a record, as capture keeps. Records are
//! biallelic; `refine_mem_probe` draws multi-allelic ones.
//!
//! `run` loads every record into memory before the clock starts, then a thread of its own hands
//! them in batches of 256 to this thread, which writes each through `VcfWriter` with the pass on
//! (`density` as `--refine-margin-density`, 0 for none), on a pool of `threads`. It prints the wall
//! time, the process's CPU time, and the CPU time of this thread — the writing thread every line
//! passes through, whose work the pool cannot share. Hash the output to compare two builds: a
//! `.vcf.gz` is cut into blocks by its bytes alone, so equal text gives equal files.
//!
//! Runs with mimalloc, the allocator `pvc` declares: under the system allocator the rewrite cost
//! about three times as much. Built with `--no-default-features --features dhat-heap`, it prints
//! also the peak of the bytes asked of the allocator while writing: a measure that does not move
//! with the allocator's own caching, as the resident memory does. What was allocated before the
//! clock started — the loaded records, and so each line's text until the stage compresses it — is
//! not counted. dhat writes `dhat-heap.json` in the working directory.
#![allow(
    clippy::disallowed_methods,
    reason = "a research probe outside the portable-float guarantee (clippy.toml)"
)]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc;
use std::time::Instant;

use pop_var_caller::neighbour_evidence::stage::{NeighbourEvidenceSwitch, decide};
use pop_var_caller::neighbour_evidence::windows::MarginThinning;
use pop_var_caller::read_phase::fragments::{FragmentAllele, RecordFragments};
use pop_var_caller::types::{ContigId, GenomePosition, Ploidy, Position};
use pop_var_caller::vcf::{HeaderContig, RecordPlace, VcfHeaderMetadata, VcfWriter};

type Failure = Box<dyn std::error::Error + Send + Sync>;

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> Result<(), Failure> {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("gen") if args.len() == 8 => generate(&args[2..]),
        Some("run") if args.len() == 7 => run(&args[2..]),
        _ => Err("usage: see the file's doc comment".into()),
    }
}
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
        // Normal approximation past 30, so deep cohorts draw quickly.
        if mean > 30.0 {
            let u = self.uniform().max(1e-12);
            let v = self.uniform();
            let normal = (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos();
            return (mean + normal * mean.sqrt()).round().max(0.0) as usize;
        }
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit {
            product *= self.uniform();
            count += 1;
        }
        count
    }
}

const READ_LENGTH: u64 = 150;
const FOUNDERS: usize = 8;
const SEGMENT_BASES: f64 = 300_000.0;
const READ_ERROR: f64 = 0.01;
const KEPT: usize = 32;

fn generate(args: &[String]) -> Result<(), Failure> {
    let out = &args[0];
    let samples: usize = args[1].parse()?;
    let records: usize = args[2].parse()?;
    let spacing: u64 = args[3].parse()?;
    let depth: f64 = args[4].parse()?;
    let mut draw = Draw(args[5].parse()?);
    let mut file = BufWriter::with_capacity(1 << 22, std::fs::File::create(out)?);
    let switch = 1.0 - (-(spacing as f64) / SEGMENT_BASES).exp();
    let mut founder_of: Vec<usize> = (0..samples).map(|_| draw.below(FOUNDERS)).collect();
    // Each sample's reads still covering the current record: (fragment, last position, allele copy).
    let mut active: Vec<Vec<(u32, u64)>> = vec![Vec::new(); samples];
    let mut next_fragment = vec![0_u32; samples];
    let starts_per_record = depth * spacing.min(READ_LENGTH) as f64 / READ_LENGTH as f64;
    let mut line = String::new();
    let mut fragments = String::new();
    for record in 0..records {
        let position = 1 + record as u64 * spacing;
        let frequency = 0.05 + 0.45 * draw.uniform();
        let founder_alleles: Vec<usize> = (0..FOUNDERS)
            .map(|_| usize::from(draw.uniform() < frequency))
            .collect();
        line.clear();
        fragments.clear();
        let mut columns = String::new();
        let mut alternative_copies = 0;
        let mut called_copies = 0;
        for sample in 0..samples {
            if record > 0 && draw.uniform() < switch {
                founder_of[sample] = draw.below(FOUNDERS);
            }
            let allele = founder_alleles[founder_of[sample]];
            active[sample].retain(|&(_, last)| last >= position);
            // A record further than a read from the last one starts every read afresh.
            let new = if spacing >= READ_LENGTH {
                draw.poisson(depth)
            } else {
                draw.poisson(starts_per_record)
            };
            for _ in 0..new {
                active[sample].push((next_fragment[sample], position + READ_LENGTH - 1));
                next_fragment[sample] += 1;
            }
            let reads = active[sample].len();
            let mut alternative_reads = 0;
            let mut kept = 0;
            for &(fragment, last) in &active[sample] {
                let mut read_allele = allele;
                if draw.uniform() < READ_ERROR {
                    read_allele = 1 - read_allele;
                }
                alternative_reads += read_allele;
                // A read reaching another record is a fragment.
                if last >= position + spacing && kept < KEPT {
                    kept += 1;
                    if !fragments.is_empty() {
                        fragments.push(';');
                    }
                    fragments.push_str(&format!("{sample},{fragment},{read_allele},30"));
                }
            }
            let log_likelihoods: Vec<f64> = [0.0, 0.5, 1.0]
                .iter()
                .map(|&share: &f64| {
                    let alternative = share * (1.0 - READ_ERROR) + (1.0 - share) * READ_ERROR;
                    alternative_reads as f64 * alternative.ln()
                        + (reads - alternative_reads) as f64 * (1.0 - alternative).ln()
                })
                .collect();
            let best = log_likelihoods.iter().copied().fold(f64::MIN, f64::max);
            let phred: Vec<u32> = log_likelihoods
                .iter()
                .map(|&value| (-10.0 * (value - best) / std::f64::consts::LN_10).round() as u32)
                .collect();
            let inbreeding = 0.95;
            let prior = [
                (1.0 - inbreeding) * (1.0 - frequency) * (1.0 - frequency)
                    + inbreeding * (1.0 - frequency),
                (1.0 - inbreeding) * 2.0 * frequency * (1.0 - frequency),
                (1.0 - inbreeding) * frequency * frequency + inbreeding * frequency,
            ];
            let unscaled: Vec<f64> = phred
                .iter()
                .zip(prior)
                .map(|(&phred, prior)| 10f64.powf(-f64::from(phred) / 10.0) * prior)
                .collect();
            let total: f64 = unscaled.iter().sum();
            let posteriors: Vec<f64> = unscaled
                .iter()
                .map(|value| (value / total * 1000.0).round() / 1000.0)
                .collect();
            let call = (0..3).fold(0, |best, genotype| {
                if posteriors[genotype] > posteriors[best] {
                    genotype
                } else {
                    best
                }
            });
            let quality = ((1.0 - posteriors[call]).max(1e-10).log10() * -10.0).round() as u32;
            columns.push('\t');
            if reads == 0 {
                columns.push_str("./.:.:0:0,0:0,0,0:");
                columns.push_str(&format!(
                    "{},{},{}",
                    posteriors[0], posteriors[1], posteriors[2]
                ));
            } else {
                called_copies += 2;
                alternative_copies += call;
                columns.push_str(&format!(
                    "{}:{}:{}:{},{}:{},{},{}:{},{},{}",
                    ["0/0", "0/1", "1/1"][call],
                    quality.min(99),
                    reads,
                    reads - alternative_reads,
                    alternative_reads,
                    phred[0],
                    phred[1],
                    phred[2],
                    posteriors[0],
                    posteriors[1],
                    posteriors[2]
                ));
            }
        }
        line.push_str(&format!(
            "chr1\t{position}\t.\tA\tC\t50\tPASS\tAF={frequency:.4};AC={alternative_copies};AN={called_copies}\tGT:GQ:DP:AD:PL:GP{columns}"
        ));
        writeln!(file, "{line}")?;
        writeln!(file, "{fragments}")?;
    }
    file.flush()?;
    Ok(())
}

/// One record as the reading thread hands it over.
type Record = (u64, Vec<u8>, RecordFragments);

fn parse_fragments(text: &str) -> Result<RecordFragments, Failure> {
    if text.is_empty() {
        return Ok(RecordFragments::default());
    }
    let entries = text
        .split(';')
        .map(|entry| -> Result<FragmentAllele, Failure> {
            let fields: Vec<u32> = entry.split(',').map(str::parse).collect::<Result<_, _>>()?;
            let [sample, fragment, allele, error_phred] = fields[..] else {
                return Err(format!("`{entry}` is not four numbers").into());
            };
            Ok(FragmentAllele {
                sample,
                fragment,
                allele: u16::try_from(allele)?,
                error_phred: u8::try_from(error_phred)?,
            })
        })
        .collect::<Result<_, _>>()?;
    Ok(RecordFragments::from_entries(entries))
}

fn cpu_seconds(clock: rustix::time::ClockId) -> f64 {
    let time = rustix::time::clock_gettime(clock);
    time.tv_sec as f64 + time.tv_nsec as f64 / 1e9
}

fn run(args: &[String]) -> Result<(), Failure> {
    let input = args[0].clone();
    let samples: usize = args[1].parse()?;
    let threads: usize = args[2].parse()?;
    let thinning = MarginThinning::parse(&args[3])?;
    let out = Path::new(&args[4]);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    let names: Vec<String> = (0..samples).map(|sample| format!("s{sample}")).collect();
    let decision = decide(
        NeighbourEvidenceSwitch::On,
        &vec![0.95; samples],
        Ploidy::try_new(2)?,
    )
    .with_margin_thinning(thinning);
    // The contig is long enough for any cohort drawn here.
    let metadata = VcfHeaderMetadata::try_new(
        vec![HeaderContig {
            name: "chr1".to_owned(),
            length: 100_000_000,
            md5: None,
        }],
        names,
        "probe".to_owned(),
        "reference.fa".to_owned(),
        "probe.log".to_owned(),
        "defaults".to_owned(),
    )?
    .the_neighbour_evidence_decision(decision);
    let mut writer = VcfWriter::create(out, metadata, Ploidy::try_new(2)?)?;

    // Every record read into memory before the clock starts, so reading the input is not
    // measured: the thread that hands them over only sends batches already made.
    let mut loaded: Vec<Vec<Record>> = Vec::new();
    {
        let mut lines = BufReader::with_capacity(1 << 22, std::fs::File::open(input)?).lines();
        let mut batch = Vec::with_capacity(256);
        while let Some(line) = lines.next() {
            let line = line?;
            let fragments = parse_fragments(&lines.next().ok_or("a record without fragments")??)?;
            let position: u64 = line.split('\t').nth(1).ok_or("no POS")?.parse()?;
            batch.push((position, line.into_bytes(), fragments));
            if batch.len() == 256 {
                loaded.push(std::mem::take(&mut batch));
                batch.reserve(256);
            }
        }
        if !batch.is_empty() {
            loaded.push(batch);
        }
    }
    let (batches, receive) = mpsc::sync_channel::<Vec<Record>>(1);
    let handing = std::thread::spawn(move || {
        for batch in loaded {
            if batches.send(batch).is_err() {
                break;
            }
        }
    });

    #[cfg(feature = "dhat-heap")]
    let profiler = dhat::Profiler::new_heap();
    let began = Instant::now();
    let process_began = cpu_seconds(rustix::time::ClockId::ProcessCPUTime);
    let thread_began = cpu_seconds(rustix::time::ClockId::ThreadCPUTime);
    let mut records = 0_u64;
    while let Ok(batch) = receive.recv() {
        for (position, line, fragments) in batch {
            let place = RecordPlace {
                at: GenomePosition {
                    contig: ContigId(0),
                    position: Position(position),
                },
                is_repeat_tract: false,
            };
            writer.write_line(place, line, fragments)?;
            records += 1;
        }
    }
    let finished = writer.finish()?;
    let wall = began.elapsed().as_secs_f64();
    let process = cpu_seconds(rustix::time::ClockId::ProcessCPUTime) - process_began;
    let writing_thread = cpu_seconds(rustix::time::ClockId::ThreadCPUTime) - thread_began;
    handing.join().map_err(|_| "the handing thread panicked")?;
    #[cfg(feature = "dhat-heap")]
    {
        let heap = dhat::HeapStats::get();
        eprintln!(
            "dhat: peak {:.1} MB asked of the allocator while writing",
            heap.max_bytes as f64 / 1e6
        );
        drop(profiler);
    }

    eprintln!(
        "records {records}, written {}, threads {threads}: wall {wall:.2} s, cpu {process:.2} s \
         ({:.2} cores busy), writing thread cpu {writing_thread:.2} s",
        finished.records_written,
        process / wall
    );
    Ok(())
}
