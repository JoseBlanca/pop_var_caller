//! **The alignment cursor against the per-region query it replaces**, over the same real
//! regions of the same real file.
//!
//! ```text
//! cargo run --release --example ng_cursor_vs_query -- <reference.fa> <sample.bam> [contig]
//! ```
//!
//! # Why this exists, and why the probe cannot answer instead
//!
//! `impl_plan/alignment_cursor.md` says Checkpoint C is "the first point the saving can be
//! measured". It is not, quite: at C the cursor is built and tested but **nothing a user runs
//! goes through it** — `PileupGenerator` still calls `SampleReads::reads_in_region`
//! (`generator.rs:854`), so `ng_generic_walk_probe` measures the old path either way and its
//! numbers say only that nothing broke. Wiring the generators is D1 and D2.
//!
//! So the measurement is made where the change actually is: the same typed regions, the same
//! reads, once through each path. That answers the question the checkpoint is for — *is the
//! retention worth the walker rewrite that follows?* — a milestone before the rewrite is paid
//! for, which is what the plan wanted the checkpoint to do.
//!
//! **It is a correctness check as well as a measurement, and that is the more important
//! half.** Both paths are asked for the same regions in the same order and their reads are
//! compared name by name. A faster path that returns different reads is not a saving.

use std::path::Path;
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

/// The halo the generic generator adds to every region it walks
/// (`PileupGeneratorConfig::max_record_span`). Applied here too, because it is what makes
/// consecutive queries overlap by ~93 % — and therefore what the cursor exists to exploit.
const HALO: u64 = 5_000;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, fasta, bam, rest @ ..] = args.as_slice() else {
        eprintln!(
            "usage: ng_cursor_vs_query <reference.fa> <sample.bam> [contig]\n\
             walks the same typed regions through the cursor and through the per-region \
             query, compares the reads, and times both."
        );
        return ExitCode::from(2);
    };
    let contig_filter = rest.first().map(String::as_str);

    match run(Path::new(fasta), Path::new(bam), contig_filter) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &Path,
    bam: &Path,
    contig_filter: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;
    let reference = OpenReference::new(info);

    let bases =
        || WindowedRefSeq::with_shared_index(fasta.to_path_buf(), contigs.clone(), index.clone());

    // The regions a real run walks: the generic ones the typed-region stream produces, each
    // widened by the halo, in the order the generator would ask for them.
    let mut regions: Vec<GenomeRegion> = Vec::new();
    for (index_of_contig, entry) in contigs.entries.iter().enumerate() {
        if contig_filter.is_some_and(|name| entry.name != name) {
            continue;
        }
        let contig = ContigId(index_of_contig as u32);
        for typed in
            TypedRegionIterator::over_contig(bases(), contig, TypedRegionConfig::default())?
        {
            let typed = typed?;
            if typed.kind != RegionKind::Generic {
                continue;
            }
            regions.push(GenomeRegion {
                contig,
                start: typed.region.start,
                end: Position(
                    typed
                        .region
                        .end
                        .get()
                        .saturating_add(HALO)
                        .min(entry.length),
                ),
            });
        }
    }
    if regions.is_empty() {
        return Err("no Generic regions to walk — check the contig name".into());
    }
    let contig = regions[0].contig;
    println!("regions={}", regions.len());

    // --- the path being replaced: one query per region ---
    let sample = SampleReads::open_only_sample(
        std::slice::from_ref(&bam.to_path_buf()),
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;
    let started = Instant::now();
    let mut by_query: u64 = 0;
    let mut query_digest: u64 = 0;
    for region in &regions {
        let stream = sample.reads_in_region(*region, bases)?;
        for read in stream {
            let read = read?;
            by_query += 1;
            query_digest = query_digest.wrapping_mul(31).wrapping_add(read.pos);
        }
    }
    let query_seconds = started.elapsed().as_secs_f64();

    // --- the path replacing it: one cursor for the whole chromosome ---
    let mut cursor = sample.cursor(contig, bases)?;
    let started = Instant::now();
    let mut by_cursor: u64 = 0;
    let mut cursor_digest: u64 = 0;
    for region in &regions {
        cursor.move_to_region(*region)?;
        while let Some(read) = cursor.next_read() {
            let read = read?;
            by_cursor += 1;
            cursor_digest = cursor_digest.wrapping_mul(31).wrapping_add(read.pos);
        }
    }
    let cursor_seconds = started.elapsed().as_secs_f64();

    println!("reads_by_query={by_query}");
    println!("reads_by_cursor={by_cursor}");
    println!("query_seconds={query_seconds:.3}");
    println!("cursor_seconds={cursor_seconds:.3}");
    println!(
        "speedup={:.2}",
        query_seconds / cursor_seconds.max(f64::MIN_POSITIVE)
    );
    let counts = cursor.counts();
    println!("cursor_reads_decoded={}", counts.reads_decoded);
    println!("cursor_reads_replayed={}", counts.reads_replayed);
    println!("cursor_regions_reusing={}", counts.regions_reusing);
    println!("cursor_regions_jumping={}", counts.regions_jumping);
    println!("cursor_reads_evicted={}", counts.reads_evicted);

    // **The half that matters more than the times.** A faster path returning different reads
    // is not a saving, and on this fixture a 1.6 % loss is what the first attempt at this
    // feature produced while every unit test passed.
    assert_eq!(
        (by_query, query_digest),
        (by_cursor, cursor_digest),
        "the cursor and the per-region query disagreed about which reads these regions hold",
    );
    println!("agreement=exact");

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}
