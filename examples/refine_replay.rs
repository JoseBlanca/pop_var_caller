//! **Replays genotype refinement over a slice of a run's spill** — a measurement tool for the
//! refinement stage's memory and its output, at a real cohort's scale, without calling again.
//!
//! ```text
//! refine_replay replay    <slice> <entries> <header.vcf> <parameters.toml> <out.vcf> <window_bytes> [floor]
//! refine_replay breakdown <slice> <entries> <header.vcf> <out.tsv>
//! ```
//!
//! `replay` pushes every entry's line and read alleles through the stage the writer runs
//! (`NeighbourEvidenceStage`), as the filter's pass three does, and writes what it hands back, one
//! line each, to `out.vcf` (no header): two builds compared on the same slice must write the same
//! bytes. The decision is the run's: refinement on, each sample's inbreeding from the parameters'
//! `[inbreeding]` table in the header's sample order. Peak memory and the wall time go to stderr.
//!
//! `breakdown` solves nothing: one row a record of what the window budget counts for it —
//! `position line line_compressed is_site record_heap solve_held read_alleles_heap read_alleles
//! core alleles` — to be summed over any window.
//!
//! `header.vcf` is the run's header, through `#CHROM`; `entries` is what `spill_slice` printed.
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::ExitCode;
use std::time::Instant;

use pop_var_caller::neighbour_evidence::line::site_of;
use pop_var_caller::neighbour_evidence::refine::record_window_bytes;
use pop_var_caller::neighbour_evidence::stage::{
    NeighbourEvidenceStage, NeighbourEvidenceSwitch, StageStop, decide,
};
use pop_var_caller::read_phase::fragments::PackedRecordFragments;
use pop_var_caller::run::paralog_filter::spill::SpillReaderWithFragments;
use pop_var_caller::types::Ploidy;
use pop_var_caller::vcf::RecordPlace;

type Failure = Box<dyn std::error::Error>;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let outcome = match args.get(1).map(String::as_str) {
        Some("replay") if args.len() == 8 || args.len() == 9 => replay(&args[2..]),
        Some("breakdown") if args.len() == 6 => breakdown(&args[2..]),
        _ => {
            eprintln!(
                "usage: refine_replay replay <slice> <entries> <header.vcf> <parameters.toml> <out.vcf> <window_bytes> [floor]\n       refine_replay breakdown <slice> <entries> <header.vcf> <out.tsv>"
            );
            return ExitCode::from(2);
        }
    };
    match outcome {
        Ok(code) => code,
        Err(why) => {
            eprintln!("refine_replay: {why}");
            ExitCode::FAILURE
        }
    }
}

/// The header's sample names and each contig's length, in order.
fn read_header(path: &str) -> Result<(Vec<String>, Vec<u64>), Failure> {
    let mut lengths = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if let Some(contig) = line.strip_prefix("##contig=<") {
            let length = contig
                .split(',')
                .find_map(|field| field.strip_prefix("length="))
                .ok_or("a contig line without its length")?;
            lengths.push(length.trim_end_matches('>').parse()?);
        } else if line.starts_with("#CHROM") {
            let samples = line.split('\t').skip(9).map(str::to_string).collect();
            return Ok((samples, lengths));
        }
    }
    Err("no #CHROM line in the header".into())
}

/// Each sample's inbreeding coefficient from the parameters, in `samples`' order.
fn inbreeding(parameters: &str, samples: &[String]) -> Result<Vec<f64>, Failure> {
    let table: toml::Table = std::fs::read_to_string(parameters)?.parse()?;
    let rows = table
        .get("inbreeding")
        .and_then(|section| section.get("by_sample"))
        .and_then(toml::Value::as_array)
        .ok_or("no [inbreeding] by_sample in the parameters")?;
    let mut by_name = HashMap::new();
    for row in rows {
        let name = row.get("sample").and_then(toml::Value::as_str);
        let value = row
            .get("inbreeding_coefficient")
            .and_then(|coefficient| coefficient.get("value"))
            .and_then(toml::Value::as_float);
        if let (Some(name), Some(value)) = (name, value) {
            by_name.insert(name.to_string(), value);
        }
    }
    samples
        .iter()
        .map(|sample| {
            by_name
                .get(sample)
                .copied()
                .ok_or_else(|| format!("no inbreeding for sample {sample}").into())
        })
        .collect()
}

fn spill(path: &str, entries: &str) -> Result<SpillReaderWithFragments<BufReader<File>>, Failure> {
    Ok(SpillReaderWithFragments::new(
        BufReader::with_capacity(16 << 20, File::open(path)?),
        entries.parse()?,
    ))
}

/// This process's peak resident memory, in GB, from the kernel.
fn peak_gb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find_map(|line| line.strip_prefix("VmHWM:"))
                .and_then(|kb| kb.trim().trim_end_matches("kB").trim().parse::<f64>().ok())
        })
        .map_or(f64::NAN, |kb| kb * 1024.0 / 1e9)
}

/// The replay's output, and the lines of the window being written.
struct Out {
    file: BufWriter<File>,
    lines_now: u64,
    first: Option<u64>,
    last: u64,
}

impl Out {
    fn write(&mut self, place: RecordPlace, line: &[u8]) -> std::io::Result<()> {
        self.file.write_all(line)?;
        self.file.write_all(b"\n")?;
        self.first.get_or_insert(place.at.position.0);
        self.last = place.at.position.0;
        self.lines_now += 1;
        Ok(())
    }
}

/// Starts the kernel's count of this process's peak memory again from what it holds now.
fn reset_peak() {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
}

fn replay(args: &[String]) -> Result<ExitCode, Failure> {
    let (samples, lengths) = read_header(&args[2])?;
    let inbreeding = inbreeding(&args[3], &samples)?;
    let window_bytes: u64 = args[5].parse()?;
    let mut decision = decide(
        NeighbourEvidenceSwitch::On,
        &inbreeding,
        Ploidy::try_new(2).map_err(|why| format!("{why:?}"))?,
    )
    .with_window_bytes(window_bytes);
    if let Some(floor) = args.get(6) {
        decision = decision.with_window_floor(floor.parse()?);
    }
    // `REPLAY_FIXED=<core>,<margin>`: every window that size, whatever the budget counts — so two
    // builds that count differently still solve the same windows, and their outputs compare.
    if let Ok(fixed) = std::env::var("REPLAY_FIXED") {
        let (core, margin) = fixed
            .split_once(',')
            .ok_or("REPLAY_FIXED is <core>,<margin>")?;
        let (core, margin): (u64, u64) = (core.parse()?, margin.parse()?);
        decision.sizes.core_bases = core;
        decision.sizes.margin_bases = margin;
        decision.sizes.margin_floor_bases = margin;
        decision.sizes.budget_bytes = u64::MAX;
        eprintln!("every window: core {core}, margins {margin}");
    }
    eprintln!("decision: {}", decision.why);
    let mut stage = NeighbourEvidenceStage::new(&decision, lengths);
    let mut out = Out {
        file: BufWriter::with_capacity(16 << 20, File::create(&args[4])?),
        lines_now: 0,
        first: None,
        last: 0,
    };
    let started = Instant::now();
    let (mut pushed, mut written) = (0u64, 0u64);
    let mut window_started = Instant::now();

    for item in spill(&args[0], &args[1])? {
        let (entry, fragments) = item?;
        let place = RecordPlace::from(&entry);
        match stage.push_writing(place, entry.line, fragments, &mut |place, line| {
            out.write(place, &line)
        }) {
            Ok(()) => {}
            Err(StageStop::OverBudget(over)) => {
                eprintln!("over budget after {pushed} entries: {over:?}");
                eprintln!(
                    "peak {:.1} GB, {:.0} s",
                    peak_gb(),
                    started.elapsed().as_secs_f64()
                );
                return Ok(ExitCode::from(3));
            }
            Err(StageStop::Write(why)) => return Err(why.into()),
        }
        pushed += 1;
        if out.lines_now > 0 {
            // A window's core was written: what it held at its peak, and how long since the last.
            eprintln!(
                "window core {}..{}: {} lines, peak {:.1} GB, {:.0} s since the last, {pushed} pushed",
                out.first.unwrap_or(0),
                out.last,
                out.lines_now,
                peak_gb(),
                window_started.elapsed().as_secs_f64()
            );
            written += out.lines_now;
            out.lines_now = 0;
            out.first = None;
            reset_peak();
            window_started = Instant::now();
        }
    }
    match stage.finish_writing(&mut |place, line| out.write(place, &line)) {
        Ok(()) => {}
        Err(StageStop::OverBudget(over)) => {
            eprintln!("over budget at the end: {over:?}");
            return Ok(ExitCode::from(3));
        }
        Err(StageStop::Write(why)) => return Err(why.into()),
    }
    written += out.lines_now;
    out.file.flush()?;
    eprintln!(
        "done: {pushed} pushed, {written} written, peak {:.1} GB, {:.0} s",
        peak_gb(),
        started.elapsed().as_secs_f64()
    );
    Ok(ExitCode::SUCCESS)
}

fn breakdown(args: &[String]) -> Result<ExitCode, Failure> {
    let (samples, _) = read_header(&args[2])?;
    let samples = samples.len();
    let mut out = BufWriter::new(File::create(&args[3])?);
    writeln!(
        out,
        "position\tline\tline_compressed\tis_site\trecord_heap\tsolve_held\tread_alleles_heap\tread_alleles\tcore\talleles"
    )?;
    let mut compressor = zstd::bulk::Compressor::new(1)?;
    for item in spill(&args[0], &args[1])? {
        let (entry, fragments) = item?;
        let compressed = compressor.compress(&entry.line)?.len();
        let read_alleles = fragments.len();
        let (is_site, heap, held, alleles_heap, core, alleles) = match site_of(&entry.line, samples)
        {
            Some(record) => {
                let (held, core) = record_window_bytes(&record);
                let heap = record.heap_bytes() as u64;
                let packed = PackedRecordFragments::from(fragments);
                (
                    1,
                    heap,
                    held - heap,
                    packed.heap_bytes() as u64,
                    core,
                    record.allele_count,
                )
            }
            None => (0, 0, 0, 0, 0, 0),
        };
        writeln!(
            out,
            "{}\t{}\t{compressed}\t{is_site}\t{heap}\t{held}\t{alleles_heap}\t{read_alleles}\t{core}\t{alleles}",
            entry.position.0,
            entry.line.len()
        )?;
    }
    out.flush()?;
    Ok(ExitCode::SUCCESS)
}
