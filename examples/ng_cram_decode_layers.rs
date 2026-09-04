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
    /// Run one pass alone and report nothing but its wall time — so the process holds only
    /// what that pass needs, and `/usr/bin/time -l` measures that pass's memory rather than
    /// the union of every pass's.
    only: Option<Pass>,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut reference = None;
        let mut cram = None;
        let mut contig = None;
        let mut start = 1u64;
        let mut containers = 50usize;
        let mut repeats = 3usize;
        let mut only = None;

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
                "--only" => {
                    only = Some(match value()?.as_str() {
                        "contig" => Pass::DiscardTags,
                        "window" => Pass::Window,
                        other => return Err(format!("--only takes contig or window, not {other}")),
                    })
                }
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
            only,
        })
    }
}

fn run(args: &Args) -> io::Result<()> {
    let repository = build_fasta_repository(&args.reference)
        .map_err(|error| io::Error::other(error.to_string()))?;

    let mut reference: IndexedFasta =
        fasta::io::indexed_reader::Builder::default().build_from_path(&args.reference)?;

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

    let index_entries = cram::crai::fs::read(&index_path(&args.cram))?.len();
    println!(
        "index     {index_entries} entries, {:.1} MB held for the whole file at {} bytes each",
        (index_entries * std::mem::size_of::<cram::crai::Record>()) as f64 / 1e6,
        std::mem::size_of::<cram::crai::Record>(),
    );
    println!(
        "file      {}\ncontig    {} (id {contig_id}) from {}\ncontainers {} \
         (index offsets, de-duplicated)\nrepeats   {}\n",
        args.cram.display(),
        args.contig,
        args.start,
        offsets.len(),
        args.repeats,
    );

    // One pass alone, so the process holds only what that pass needs and an external memory
    // measurement is of that pass rather than of every pass together.
    if let Some(pass) = args.only {
        let mut elapsed = Duration::MAX;
        for _ in 0..args.repeats {
            elapsed = elapsed.min(time_pass(
                &mut reader,
                &header,
                &repository,
                &mut reference,
                offsets,
                pass,
            )?);
        }
        println!(
            "{} containers, {}: {:.3} s",
            offsets.len(),
            match pass {
                Pass::Window => "decoded against a per-slice window of the reference",
                _ => "decoded against the whole contig, held resident",
            },
            elapsed.as_secs_f64(),
        );
        return Ok(());
    }

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
    let mut direct = Duration::MAX;
    let mut discard_tags = Duration::MAX;
    let mut digest = Duration::MAX;

    for _ in 0..args.repeats {
        read = read.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::Read,
        )?);
        blocks = blocks.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::Blocks,
        )?);
        records = records.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::Records,
        )?);
        convert = convert.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::Convert,
        )?);
        direct = direct.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::Direct,
        )?);
        discard_tags = discard_tags.min(time_pass(
            &mut reader,
            &header,
            &repository,
            &mut reference,
            offsets,
            Pass::DiscardTags,
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
    print_codec_costs(&mut reader, &header, &repository, &mut reference, offsets)?;

    println!(
        "\nReading each field straight off the CRAM record instead of through `RecordBuf`: \
         {:.3} s for the whole walk against {:.3} s, so the last layer costs {:.3} s rather \
         than {:.3} s.",
        direct.as_secs_f64(),
        convert.as_secs_f64(),
        (direct.as_secs_f64() - records.as_secs_f64()).max(0.0),
        (convert.as_secs_f64() - records.as_secs_f64()).max(0.0),
    );

    println!(
        "\nDecoding the auxiliary tags but not keeping them, which is every tag ng reads: \
         {:.3} s for the whole walk against {:.3} s.",
        discard_tags.as_secs_f64(),
        direct.as_secs_f64(),
    );

    let through_record_buf = checksum_pass(
        &mut reader,
        &header,
        &repository,
        &mut reference,
        offsets,
        Pass::Convert,
    )?;
    let read_directly = checksum_pass(
        &mut reader,
        &header,
        &repository,
        &mut reference,
        offsets,
        Pass::Direct,
    )?;
    let without_tags = checksum_pass(
        &mut reader,
        &header,
        &repository,
        &mut reference,
        offsets,
        Pass::DiscardTags,
    )?;
    let over_a_window = checksum_pass(
        &mut reader,
        &header,
        &repository,
        &mut reference,
        offsets,
        Pass::Window,
    )?;
    println!(
        "\nEvery decoded field of all {} records hashes to {through_record_buf} through \
         `RecordBuf` and to {read_directly} read directly — {}.",
        shape.records,
        if through_record_buf == read_directly {
            "the same records"
        } else {
            "THESE DIFFER"
        },
    );
    println!(
        "With the tags decoded and discarded they hash to {without_tags} — {}.",
        if without_tags == read_directly {
            "still the same records"
        } else {
            "THESE DIFFER"
        },
    );
    println!(
        "Decoded against a per-slice window of the reference they hash to {over_a_window} — {}.",
        if over_a_window == read_directly {
            "still the same records"
        } else {
            "THESE DIFFER"
        },
    );

    println!(
        "\nThroughput at the last layer: {:.2} M records/s, {:.1} MiB/s of container bytes.",
        shape.records as f64 / total / 1e6,
        shape.container_bytes as f64 / total / (1024.0 * 1024.0),
    );

    let slice = cram::perf::read_slice_counters();
    if slice.records > 0 {
        println!(
            "\nTurning {} decoded slices' blocks into records:\n\
             \x20 fetching and digesting the reference   {:.3} s\n\
             \x20 allocating the record vector           {:.3} s\n\
             \x20 reading the records                    {:.3} s, of which\n\
             \x20   the auxiliary tags                   {:.3} s\n\
             \x20 linking each record to its mate        {:.3} s",
            slice.records,
            slice.reference_nanos as f64 / 1e9,
            slice.allocate_nanos as f64 / 1e9,
            slice.read_nanos as f64 / 1e9,
            slice.tag_nanos as f64 / 1e9,
            slice.resolve_mates_nanos as f64 / 1e9,
        );
    }

    Ok(())
}

#[derive(Clone, Copy, PartialEq)]
enum Pass {
    Read,
    Blocks,
    Records,
    Convert,
    /// What `Convert` does, but reading each field straight off the CRAM record into flat
    /// buffers instead of going through `RecordBuf` and the eight boxed accessors behind it.
    Direct,
    /// What `Direct` does, and the auxiliary tags decoded but not kept.
    DiscardTags,
    /// What `DiscardTags` does, decoding against a window of the reference fetched for each
    /// slice rather than against a whole contig held resident.
    Window,
}

/// One window of the reference, read out of the indexed FASTA — which is what a caller
/// holding no contig would do, and the reason this harness does not take it from the
/// repository: taking it from there would load the contig this is meant to avoid.
fn fetch_window(
    reference: &mut IndexedFasta,
    contig: &[u8],
    start: noodles_core::Position,
    end: noodles_core::Position,
) -> io::Result<fasta::Record> {
    let region =
        noodles_core::Region::new(contig, noodles_core::region::Interval::from(start..=end));
    reference.query(&region)
}

fn time_pass(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    reference: &mut IndexedFasta,
    offsets: &[u64],
    pass: Pass,
) -> io::Result<Duration> {
    time_pass_checked(reader, header, repository, reference, offsets, pass, false)
        .map(|(elapsed, _)| elapsed)
}

/// The FASTA, open and indexed, for the windowed pass.
type IndexedFasta = fasta::io::IndexedReader<fasta::io::BufReader<File>>;

/// Every field of every decoded record, hashed — the oracle a change to the decoder has to
/// leave unmoved. Sequence and CIGAR are reconstructed from the reference, so a wrong
/// reference or a wrong feature decode moves this and nothing else in the harness would.
fn checksum_pass(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    reference: &mut IndexedFasta,
    offsets: &[u64],
    pass: Pass,
) -> io::Result<String> {
    let (_, digest) =
        time_pass_checked(reader, header, repository, reference, offsets, pass, true)?;
    Ok(digest)
}

fn time_pass_checked(
    reader: &mut cram::io::Reader<File>,
    header: &sam::Header,
    repository: &fasta::Repository,
    reference: &mut IndexedFasta,
    offsets: &[u64],
    pass: Pass,
    checksum: bool,
) -> io::Result<(Duration, String)> {
    let mut digest_of_records = Md5::new();
    let mut container = cram::io::reader::Container::default();
    let mut record_buf = sam::alignment::RecordBuf::default();
    let mut sink = Sink::default();
    // The buffers the direct pass reuses across every record of the walk, standing in for the
    // flat buffers ng's decoded container keeps.
    let mut bases: Vec<u8> = Vec::new();
    let mut quality_scores: Vec<u8> = Vec::new();
    let mut cigar: Vec<sam::alignment::record::cigar::Op> = Vec::new();
    let mut window: Vec<u8> = Vec::new();
    let mut payload: Vec<u8> = Vec::new();
    let mut cigar_ops: Vec<sam::alignment::record::cigar::Op> = Vec::new();

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

            let records = if pass == Pass::Window {
                // What a caller walking the genome in order would do: ask the slice which
                // bases it needs, fetch exactly those, and hand them over. Fetched here from
                // the same repository for want of a second reference reader in this harness —
                // the point being measured is what the decode needs resident, not how the
                // window is read.
                let (contig, start, end) = slice
                    .reference_span()
                    .expect("a mapped slice names one reference sequence");
                let name = header
                    .reference_sequences()
                    .get_index(contig)
                    .map(|(name, _)| name.clone())
                    .expect("the slice's reference sequence is in the header");
                window.clear();
                let fetched = fetch_window(reference, &name, start, end)?;
                window.extend_from_slice(fetched.sequence().as_ref());
                slice.records_over_window(
                    &window,
                    start,
                    repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?
            } else if pass == Pass::DiscardTags {
                slice.records_discarding_tags(
                    repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?
            } else {
                slice.records(
                    repository.clone(),
                    header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?
            };
            if pass == Pass::Records {
                sink.absorb_usize(records.len());
                continue;
            }

            if matches!(pass, Pass::Direct | Pass::DiscardTags | Pass::Window) {
                for record in &records {
                    sink.absorb_usize(record.read_group_index().map_or(0, |index| index + 1));
                    if checksum {
                        hash_record_directly(&mut digest_of_records, record)?;
                    }
                    record.write_bases_into(&mut bases);
                    record.write_quality_scores_into(&mut quality_scores);
                    record.write_cigar_into(&mut cigar);
                    payload.extend_from_slice(record.name_bytes().unwrap_or_default());
                    payload.extend_from_slice(&bases);
                    payload.extend_from_slice(&quality_scores);
                    cigar_ops.extend_from_slice(&cigar);
                    sink.absorb_usize(payload.len() + cigar_ops.len());
                    payload.clear();
                    cigar_ops.clear();
                }
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

/// The same fold as [`hash_record`], over the fields read straight off the CRAM record. The
/// two must agree: that is what says the direct path decodes what the boxed one does.
fn hash_record_directly(digest: &mut Md5, record: &cram::Record<'_>) -> io::Result<()> {
    use sam::alignment::Record as _;

    let mut bases = Vec::new();
    let mut quality_scores = Vec::new();
    let mut cigar = Vec::new();

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
    digest.update(record.name_bytes().unwrap_or_default());
    record.write_bases_into(&mut bases);
    for base in &bases {
        digest.update([*base]);
    }
    record.write_quality_scores_into(&mut quality_scores);
    for score in &quality_scores {
        digest.update([*score]);
    }
    record.write_cigar_into(&mut cigar);
    for op in &cigar {
        digest.update([op.kind() as u8]);
        digest.update((op.len() as u32).to_le_bytes());
    }
    Ok(())
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
    reference: &mut IndexedFasta,
    offsets: &[u64],
) -> io::Result<()> {
    cram::perf::reset_counters();
    cram::perf::reset_rans_counters();
    cram::perf::reset_slice_counters();
    let _ = time_pass(reader, header, repository, reference, offsets, Pass::Blocks)?;
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

    let rans = cram::perf::read_rans_counters();
    if rans.order_0_blocks + rans.order_1_blocks > 0 {
        println!(
            "\nOf the rANS blocks: {} are order-0 over {:.1} MB and {} are order-1 over {:.1} MB.\n\
             Building the decode tables costs {:.3} s and decoding the symbols {:.3} s.",
            rans.order_0_blocks,
            rans.order_0_bytes as f64 / 1e6,
            rans.order_1_blocks,
            rans.order_1_bytes as f64 / 1e6,
            rans.table_nanos as f64 / 1e9,
            rans.decode_nanos as f64 / 1e9,
        );
    }

    Ok(())
}
