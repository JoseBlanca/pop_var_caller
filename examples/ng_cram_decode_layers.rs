//! **Where does the time go when ng reads a CRAM?** One file, one contiguous run of
//! containers, timed layer by layer.
//!
//! The calling path's own profile can only say "a third of the serial run is inside the CRAM
//! reader". It cannot say which *part* of the reader, because the parts are nested: pulling
//! bytes off disk, inflating the compressed blocks, decoding the records out of those blocks,
//! and copying each record into ng's own buffers all appear under one frame. This runs the
//! same containers four times, each pass doing one layer more than the last, so each layer's
//! own cost is the difference between two passes.
//!
//! ```text
//! cargo run --release --example ng_cram_decode_layers -- \
//!     --reference <reference.fa> --cram <sample.cram> \
//!     --contig <name> [--start <1-based>] [--containers <n>] [--repeats <n>]
//! ```
//!
//! ## The four passes
//!
//! | pass | what it does |
//! |---|---|
//! | `read` | seek, and read each container's bytes into memory. No decompression. |
//! | `blocks` | the above, plus parse the compression header and inflate every block. |
//! | `records` | the above, plus decode every record out of the inflated blocks. This is where noodles verifies the reference digest. |
//! | `convert` | the above, plus what ng itself does with a record: resolve its read group, clone it into a `RecordBuf`, and flatten it into flat buffers. |
//!
//! Each pass reads the same containers, so the differences are layer costs and the last pass
//! is the whole cost of ng's CRAM reading. The first pass warms the page cache; `--repeats`
//! runs the set several times and reports the fastest of each, because on a shared machine the
//! minimum is the number least polluted by other load.
//!
//! ## The digest is measured, not inferred
//!
//! For every slice, noodles hashes the reference bases the slice spans and compares against
//! the digest in the slice header (CRAM v3 §8.5). That work is inside the `records` pass and
//! cannot be switched off through any published setting. So this harness also computes,
//! separately, the same digest over the same spans, and reports it — the figure is what a
//! build that skipped the check would stop paying.

use std::fs::File;
use std::io::{self, SeekFrom};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use md5::{Digest, Md5};
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use pop_var_caller::bam::alignment_input::build_fasta_repository;

fn main() {
    let args = match Args::parse() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            eprintln!();
            eprintln!("{USAGE}");
            std::process::exit(2);
        }
    };

    if let Err(error) = run(&args) {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

const USAGE: &str = "\
usage: ng_cram_decode_layers --reference <fa> --cram <cram> --contig <name>
                             [--start <1-based>] [--containers <n>] [--repeats <n>]";

struct Args {
    reference: PathBuf,
    cram: PathBuf,
    contig: String,
    start: u64,
    containers: usize,
    repeats: usize,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut reference = None;
        let mut cram = None;
        let mut contig = None;
        let mut start = 1u64;
        let mut containers = 50usize;
        let mut repeats = 3usize;

        let mut argv = std::env::args().skip(1);
        while let Some(flag) = argv.next() {
            let mut value = || argv.next().ok_or_else(|| format!("{flag} needs a value"));
            match flag.as_str() {
                "--reference" => reference = Some(PathBuf::from(value()?)),
                "--cram" => cram = Some(PathBuf::from(value()?)),
                "--contig" => contig = Some(value()?),
                "--start" => start = value()?.parse().map_err(|_| "--start is a number")?,
                "--containers" => {
                    containers = value()?.parse().map_err(|_| "--containers is a number")?
                }
                "--repeats" => repeats = value()?.parse().map_err(|_| "--repeats is a number")?,
                other => return Err(format!("unknown argument {other}")),
            }
        }

        Ok(Self {
            reference: reference.ok_or("--reference is required")?,
            cram: cram.ok_or("--cram is required")?,
            contig: contig.ok_or("--contig is required")?,
            start,
            containers,
            repeats,
        })
    }
}

fn run(args: &Args) -> io::Result<()> {
    let repository = build_fasta_repository(&args.reference)
        .map_err(|error| io::Error::other(error.to_string()))?;

    let mut reader = cram::io::Reader::new(File::open(&args.cram)?);
    reader.read_file_definition()?;
    let header = reader.read_file_header()?;

    let contig_id = header
        .reference_sequences()
        .get_index_of(args.contig.as_bytes())
        .ok_or_else(|| io::Error::other(format!("{} is not in the header", args.contig)))?;

    let offsets = container_offsets(&index_path(&args.cram), contig_id, args.start)?;
    let offsets = &offsets[..offsets.len().min(args.containers)];
    if offsets.is_empty() {
        return Err(io::Error::other(
            "the index has no container at or after that position",
        ));
    }

    println!(
        "file      {}\ncontig    {} (id {contig_id}) from {}\ncontainers {} \
         (index offsets, de-duplicated)\nrepeats   {}\n",
        args.cram.display(),
        args.contig,
        args.start,
        offsets.len(),
        args.repeats,
    );

    let spans = spans_from_index(&index_path(&args.cram), contig_id, args.start, offsets)?;

    // Warm the page cache and collect the shape of the work, so the timings below are of
    // decoding rather than of the first disk read.
    let mut shape = walk_shape(&mut reader, &header, &repository, offsets)?;
    shape.reference_bases = spans.iter().map(|(start, end)| end - start + 1).sum();
    println!(
        "work      {} containers, {} slices, {} records, {:.1} MiB of container bytes\n\
         reference {:.1} Mb of bases spanned by those slices\n",
        shape.containers,
        shape.slices,
        shape.records,
        shape.container_bytes as f64 / (1024.0 * 1024.0),
        shape.reference_bases as f64 / 1e6,
    );

    let mut read = Duration::MAX;
    let mut blocks = Duration::MAX;
    let mut records = Duration::MAX;
    let mut convert = Duration::MAX;
    let mut digest = Duration::MAX;

    for _ in 0..args.repeats {
        read = read.min(time_pass(
            &mut reader,
            &header,
            &repository,
            offsets,
            Pass::Read,
        )?);
        blocks = blocks.min(time_pass(
            &mut reader,
            &header,
            &repository,
            offsets,
            Pass::Blocks,
        )?);
        records = records.min(time_pass(
            &mut reader,
            &header,
            &repository,
            offsets,
            Pass::Records,
        )?);
        convert = convert.min(time_pass(
            &mut reader,
            &header,
            &repository,
            offsets,
            Pass::Convert,
        )?);
        digest = digest.min(time_digest(&repository, &args.contig, &spans)?);
    }

    let total = convert.as_secs_f64();
    let row = |name: &str, cumulative: Duration, own: f64| {
        println!(
            "{name:<10} {:>8.3} s cumulative   {:>8.3} s own   {:>5.1} % of the whole",
            cumulative.as_secs_f64(),
            own,
            100.0 * own / total,
        );
    };

    println!("Fastest of {} runs of each pass:\n", args.repeats);
    row("read", read, read.as_secs_f64());
    row(
        "blocks",
        blocks,
        (blocks.as_secs_f64() - read.as_secs_f64()).max(0.0),
    );
    row(
        "records",
        records,
        (records.as_secs_f64() - blocks.as_secs_f64()).max(0.0),
    );
    row(
        "convert",
        convert,
        (convert.as_secs_f64() - records.as_secs_f64()).max(0.0),
    );

    println!(
        "\nOf the `records` layer, hashing the reference spans is {:.3} s \
         ({:.1} % of the whole), measured by digesting the same spans here.",
        digest.as_secs_f64(),
        100.0 * digest.as_secs_f64() / total,
    );

    // Which compression method the block-decode layer spent its time in. The counters live in
    // the vendored copy of noodles-cram and are exact in blocks and bytes; the nanoseconds are
    // one clock read per block.
    #[cfg(feature = "cram-perf-counters")]
    print_codec_costs(&mut reader, &header, &repository, offsets)?;

    println!(
        "\nEvery decoded field of all {} records hashes to {}.",
        shape.records,
        checksum_pass(&mut reader, &header, &repository, offsets, Pass::Convert)?,
    );

    println!(
        "\nThroughput at the last layer: {:.2} M records/s, {:.1} MiB/s of container bytes.",
        shape.records as f64 / total / 1e6,
        shape.container_bytes as f64 / total / (1024.0 * 1024.0),
    );

    Ok(())
}

#[derive(Clone, Copy, PartialEq)]
enum Pass {
    Read,
    Blocks,
    Records,
    Convert,
}

fn time_pass(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
    pass: Pass,
) -> io::Result<Duration> {
    time_pass_checked(reader, header, repository, offsets, pass, false).map(|(elapsed, _)| elapsed)
}

/// Every field of every decoded record, hashed — the oracle a change to the decoder has to
/// leave unmoved. Sequence and CIGAR are reconstructed from the reference, so a wrong
/// reference or a wrong feature decode moves this and nothing else in the harness would.
fn checksum_pass(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
    pass: Pass,
) -> io::Result<String> {
    let (_, digest) = time_pass_checked(reader, header, repository, offsets, pass, true)?;
    Ok(digest)
}

fn time_pass_checked(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
    pass: Pass,
    checksum: bool,
) -> io::Result<(Duration, String)> {
    let mut digest_of_records = Md5::new();
    let mut container = cram::io::reader::Container::default();
    let mut record_buf = sam::alignment::RecordBuf::default();
    let mut sink = Sink::default();

    let started = Instant::now();
    for &offset in offsets {
        reader.seek(SeekFrom::Start(offset))?;
        if reader.read_container(&mut container)? == 0 {
            break;
        }
        if pass == Pass::Read {
            continue;
        }

        let compression_header = container.compression_header()?;
        for slice in container.slices() {
            let slice = slice?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            if pass == Pass::Blocks {
                sink.absorb_usize(core_data_src.len() + external_data_srcs.len());
                continue;
            }

            let records = slice.records(
                repository.clone(),
                header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;
            if pass == Pass::Records {
                sink.absorb_usize(records.len());
                continue;
            }

            for record in &records {
                if checksum {
                    hash_record(&mut digest_of_records, header, record)?;
                }
                // What ng does per record: ask which read group it is in (the question that
                // decides whether it is this sample's at all), then copy it out of the
                // borrowed block into owned buffers.
                sink.absorb_usize(read_group_of(record)?.map_or(0, |id| id + 1));
                record_buf.try_clone_from_alignment_record(header, record)?;
                sink.absorb_record(&record_buf);
            }
        }
    }
    let elapsed = started.elapsed();

    // Keep every pass's result observable so nothing above is optimised away.
    std::hint::black_box(&sink);
    let digest: [u8; 16] = digest_of_records.finalize().into();
    let rendered = digest.iter().map(|byte| format!("{byte:02x}")).collect();
    Ok((elapsed, rendered))
}

/// Fold one decoded record's every field into the running digest.
fn hash_record(
    digest: &mut Md5,
    header: &sam::Header,
    record: &cram::Record<'_>,
) -> io::Result<()> {
    use sam::alignment::Record as _;

    digest.update(record.flags()?.bits().to_le_bytes());
    digest.update(
        record
            .alignment_start()
            .transpose()?
            .map_or(0u64, |position| usize::from(position) as u64)
            .to_le_bytes(),
    );
    digest.update(
        record
            .mapping_quality()
            .transpose()?
            .map_or(255u8, u8::from)
            .to_le_bytes(),
    );
    digest.update(record.name().unwrap_or_default());
    for base in record.sequence().iter() {
        digest.update([base]);
    }
    for score in record.quality_scores().iter() {
        digest.update([score?]);
    }
    for op in record.cigar().iter() {
        let op = op?;
        digest.update([op.kind() as u8]);
        digest.update((op.len() as u32).to_le_bytes());
    }
    let _ = header;
    Ok(())
}

/// The reference digest noodles computes inside `records`, computed here instead, so its cost
/// is a number rather than a share read off a profile.
fn time_digest(
    repository: &fasta::Repository,
    contig: &str,
    spans: &[(u64, u64)],
) -> io::Result<Duration> {
    let started = Instant::now();
    let mut sink = Sink::default();
    for (start, end) in spans {
        let sequence = repository
            .get(contig.as_bytes())
            .transpose()?
            .expect("the contig is in the reference");
        let bases: &[u8] = (*sequence).as_ref();
        let bases = &bases[(*start as usize - 1)..(*end as usize)];
        let digest = normalized_sequence_digest(bases);
        sink.absorb_usize(digest[0] as usize);
    }
    let elapsed = started.elapsed();
    std::hint::black_box(&sink);
    Ok(elapsed)
}

struct Shape {
    containers: usize,
    slices: usize,
    records: usize,
    container_bytes: u64,
    reference_bases: u64,
}

fn walk_shape(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
) -> io::Result<Shape> {
    let mut shape = Shape {
        containers: 0,
        slices: 0,
        records: 0,
        container_bytes: 0,
        reference_bases: 0,
    };
    let mut container = cram::io::reader::Container::default();

    for &offset in offsets {
        reader.seek(SeekFrom::Start(offset))?;
        let read = reader.read_container(&mut container)?;
        if read == 0 {
            break;
        }
        shape.containers += 1;
        shape.container_bytes += read as u64;

        let compression_header = container.compression_header()?;
        for slice in container.slices() {
            let slice = slice?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let records = slice.records(
                repository.clone(),
                header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            );
            // A slice whose records need external reference bases cannot be decoded against an
            // empty repository; count what it holds and move on. The record count comes from
            // the successful arm where it can.
            shape.slices += 1;
            if let Ok(records) = records {
                shape.records += records.len();
            }
        }
    }

    Ok(shape)
}

fn read_group_of(record: &cram::Record<'_>) -> io::Result<Option<usize>> {
    use sam::alignment::Record as _;
    use sam::alignment::record::data::field::{Tag, Value};

    match record.data().get(&Tag::READ_GROUP).transpose()? {
        Some(Value::String(name)) => Ok(Some(name.len())),
        Some(_) => Ok(None),
        None => Ok(None),
    }
}

/// A place to put every value a pass produces, so the optimiser cannot delete the work.
#[derive(Default)]
struct Sink {
    total: u64,
}

impl Sink {
    fn absorb_usize(&mut self, value: usize) {
        self.total = self.total.wrapping_add(value as u64);
    }

    fn absorb_record(&mut self, record: &sam::alignment::RecordBuf) {
        self.total = self
            .total
            .wrapping_add(record.sequence().len() as u64)
            .wrapping_add(record.cigar().as_ref().len() as u64);
    }
}

fn index_path(cram: &Path) -> PathBuf {
    let mut path = cram.as_os_str().to_owned();
    path.push(".crai");
    PathBuf::from(path)
}

/// Container offsets from the `.crai`, for one contig, from `start` onward, in file order and
/// with the repeats removed — a container holding several slices appears once per slice.
fn container_offsets(index: &Path, contig_id: usize, start: u64) -> io::Result<Vec<u64>> {
    let mut offsets = Vec::new();
    for entry in cram::crai::fs::read(index)? {
        if !covers(&entry, contig_id, start) {
            continue;
        }
        if offsets.last() != Some(&entry.offset()) {
            offsets.push(entry.offset());
        }
    }
    Ok(offsets)
}

/// Whether an index entry belongs to this contig and ends at or after `start`.
fn covers(entry: &cram::crai::Record, contig_id: usize, start: u64) -> bool {
    if entry.reference_sequence_id() != Some(contig_id) {
        return false;
    }
    let Some(alignment_start) = entry.alignment_start() else {
        return false;
    };
    usize::from(alignment_start) as u64 + entry.alignment_span() as u64 >= start
}

/// The reference span of every slice in the walked containers — the same
/// `alignment_start..=alignment_end` noodles digests, read out of the index.
fn spans_from_index(
    index: &Path,
    contig_id: usize,
    start: u64,
    containers: &[u64],
) -> io::Result<Vec<(u64, u64)>> {
    let mut spans = Vec::new();
    for entry in cram::crai::fs::read(index)? {
        if !covers(&entry, contig_id, start) || !containers.contains(&entry.offset()) {
            continue;
        }
        let alignment_start = usize::from(entry.alignment_start().unwrap()) as u64;
        spans.push((
            alignment_start,
            alignment_start + entry.alignment_span() as u64 - 1,
        ));
    }
    Ok(spans)
}

/// The digest noodles compares a slice's reference span against: MD5 of the bases with
/// lowercase folded up and anything outside the printable range dropped (CRAM v3 §8.5).
fn normalized_sequence_digest(bases: &[u8]) -> [u8; 16] {
    let mut hasher = Md5::new();
    if bases
        .iter()
        .all(|&base| base.is_ascii_uppercase() || !base.is_ascii_alphabetic())
    {
        hasher.update(bases);
    } else {
        let mut chunk = [0u8; 512];
        for window in bases.chunks(512) {
            let mut kept = 0;
            for &base in window {
                if base.is_ascii_graphic() {
                    chunk[kept] = base.to_ascii_uppercase();
                    kept += 1;
                }
            }
            hasher.update(&chunk[..kept]);
        }
    }
    hasher.finalize().into()
}

/// What each compression method cost, from the vendored copy's counters.
#[cfg(feature = "cram-perf-counters")]
fn print_codec_costs(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
) -> io::Result<()> {
    cram::perf::reset_counters();
    let _ = time_pass(reader, header, repository, offsets, Pass::Blocks)?;
    let counters = cram::perf::read_counters();
    println!("\nBlock decompression, by compression method (one pass):\n");
    println!(
        "{:<16} {:>8} {:>12} {:>12} {:>9} {:>8}",
        "method", "blocks", "compressed", "inflated", "seconds", "share"
    );
    let decompression: f64 = counters
        .iter()
        .map(|cost| cost.nanos as f64 / 1e9)
        .sum::<f64>()
        .max(f64::MIN_POSITIVE);
    for (name, cost) in cram::perf::METHOD_NAMES.iter().zip(&counters) {
        if cost.blocks == 0 {
            continue;
        }
        println!(
            "{name:<16} {:>8} {:>10.1} MB {:>10.1} MB {:>9.3} {:>7.1} %",
            cost.blocks,
            cost.compressed_bytes as f64 / 1e6,
            cost.uncompressed_bytes as f64 / 1e6,
            cost.nanos as f64 / 1e9,
            100.0 * (cost.nanos as f64 / 1e9) / decompression,
        );
    }

    Ok(())
}
