//! **Copies the entries of a hidden-duplication spill that fall in one region into a smaller
//! spill** — a measurement tool, so refinement can be replayed on a run's real records without
//! calling them again.
//!
//! ```text
//! spill_slice <spill> <out> <contig_index> <start> <end>
//! ```
//!
//! `contig_index` is the contig's place in the run's header, from 0; `start` and `end` are written
//! positions, the end excluded. Each entry is copied byte for byte, both frames, so the output is
//! a spill [`SpillReaderWithFragments`](pop_var_caller::run::paralog_filter::spill) reads. Only the
//! start of each body is inflated, for its contig and position; the walk stops at the first entry
//! past the region. Prints the number of entries copied, which the reader needs.
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

/// Reads one LEB128 varint, appending its bytes to `raw`; `None` at a clean end of file.
fn varint(source: &mut impl Read, raw: &mut Vec<u8>) -> std::io::Result<Option<u64>> {
    let mut value = 0u64;
    let mut shift = 0;
    loop {
        let mut byte = [0u8; 1];
        match source.read(&mut byte)? {
            0 if shift == 0 => return Ok(None),
            0 => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "a varint cut short",
                ));
            }
            _ => {}
        }
        raw.push(byte[0]);
        value |= u64::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
    }
}

/// A varint from the front of `bytes`, and how many bytes it took.
fn varint_in(bytes: &[u8]) -> Option<(u64, usize)> {
    let mut value = 0u64;
    for (index, byte) in bytes.iter().enumerate().take(10) {
        value |= u64::from(byte & 0x7f) << (7 * index);
        if byte & 0x80 == 0 {
            return Some((value, index + 1));
        }
    }
    None
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 6 {
        eprintln!("usage: spill_slice <spill> <out> <contig_index> <start> <end>");
        return ExitCode::from(2);
    }
    let (contig, start, end): (u64, u64, u64) =
        match (args[3].parse(), args[4].parse(), args[5].parse()) {
            (Ok(contig), Ok(start), Ok(end)) => (contig, start, end),
            _ => {
                eprintln!("contig_index, start and end are whole numbers");
                return ExitCode::from(2);
            }
        };
    match slice(&args[1], &args[2], contig, start, end) {
        Ok(copied) => {
            println!("{copied}");
            ExitCode::SUCCESS
        }
        Err(why) => {
            eprintln!("spill_slice: {why}");
            ExitCode::FAILURE
        }
    }
}

fn slice(spill: &str, out: &str, contig: u64, start: u64, end: u64) -> std::io::Result<u64> {
    let mut source = BufReader::with_capacity(16 << 20, File::open(spill)?);
    let mut sink = BufWriter::with_capacity(16 << 20, File::create(out)?);
    let mut raw = Vec::new();
    let mut frame = Vec::new();
    let mut head = Vec::with_capacity(32);
    let (mut read, mut copied, mut bytes_read) = (0u64, 0u64, 0u64);
    let started = std::time::Instant::now();
    loop {
        raw.clear();
        let Some(body_bytes) = varint(&mut source, &mut raw)? else {
            break;
        };
        frame.resize(body_bytes as usize, 0);
        source.read_exact(&mut frame)?;
        raw.extend_from_slice(&frame);
        // The contig and the position open the body: inflate only that much.
        head.clear();
        zstd::stream::read::Decoder::with_buffer(&frame[..])?
            .take(24)
            .read_to_end(&mut head)?;
        let (entry_contig, used) = varint_in(&head).ok_or_else(|| bad("the body's contig"))?;
        let (position, _) = varint_in(&head[used..]).ok_or_else(|| bad("the body's position"))?;
        let fragments_bytes =
            varint(&mut source, &mut raw)?.ok_or_else(|| bad("an entry cut short"))?;
        frame.resize(fragments_bytes as usize, 0);
        source.read_exact(&mut frame)?;
        raw.extend_from_slice(&frame);
        read += 1;
        bytes_read += raw.len() as u64;
        if entry_contig > contig || (entry_contig == contig && position >= end) {
            break;
        }
        if entry_contig == contig && position >= start {
            sink.write_all(&raw)?;
            copied += 1;
        }
        if read % 200_000 == 0 {
            eprintln!(
                "{read} entries read, {:.1} GB, at contig {entry_contig} position {position}, {copied} copied, {:.0} s",
                bytes_read as f64 / 1e9,
                started.elapsed().as_secs_f64()
            );
        }
    }
    sink.flush()?;
    eprintln!(
        "done: {read} entries read ({:.1} GB), {copied} copied, {:.0} s",
        bytes_read as f64 / 1e9,
        started.elapsed().as_secs_f64()
    );
    Ok(copied)
}

fn bad(what: &str) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("could not read {what}"),
    )
}
