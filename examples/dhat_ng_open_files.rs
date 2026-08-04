//! **What does one *open* alignment file cost in memory?**
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_ng_open_files --features dhat-heap -- \
//!     <reference.fa> <regions.bed> <n> <sample.cram> [sample ...]
//! ```
//!
//! ## The question
//!
//! After the reference-sharing and one-contig fixes
//! (`doc/devel/reports/implementations/ng_cram_reference_sharing_2026-07-29.md`),
//! a cohort's peak memory is ~166 MiB fixed plus **~14.6 MiB for every file
//! opened**. The fixed part is accounted for — one chromosome of reference
//! bases. The per-file part is not: it was attributed by inspection to "the
//! reader pool, the index, the `.crai` table and the decoded container", and
//! that guess is already partly wrong, since these `.crai` files are ~1.6 KB
//! each. This measures it instead.
//!
//! ## Why this shape
//!
//! The cost is a property of an **open, queried file**, not of the STR walk
//! above it, so this runs the input layer alone: open k files as one cohort,
//! then stream every BED span through every file. That is the same
//! open-and-query traffic the cohort dump generates, with the locus generator
//! and the output writer removed — so whatever scales with k here scales with k
//! there, and nothing else is in the way.
//!
//! Run it at two file counts and diff. A term that is per-file appears k times;
//! the fixed terms cancel. `dhat-heap.json` holds the full attribution; compare
//! two runs with `tmp/dhat_diff.py`.
//!
//! ## Read `t-gmax` with care — it cost this investigation an hour
//!
//! DHAT's `At t-gmax` is live heap **at one instant**: the global peak. Below
//! ~16 files that instant falls during reference loading, so `t-gmax` came out
//! identical at 1 and 8 files and appeared to prove that nothing is retained per
//! file. It proved no such thing — the per-file retention simply had not yet
//! grown past the reference spike, and above 16 files the peak moves and it
//! shows up. Comparing a peak across configurations is only valid when the peak
//! falls in the same phase.
//!
//! The reliable measurement is the explicit `HeapStats::get()` probe at the end
//! of the walk, with the files still open: that is *this* many files' retained
//! heap, at a known point, with nothing to interpret.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

// Opt-in mimalloc global allocator (--features alloc-mimalloc). Kept because
// "is this just the allocator failing to return freed pages?" is the first
// alternative explanation to rule out, and swapping allocators rules it out in
// one run — here it did not move the per-file slope (12.0 -> 10.8 MiB), which
// is what sent the investigation back to looking for genuinely live memory.
// Mutually exclusive with dhat, which needs the global allocator slot itself.
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::GenomeRegions;
use pop_var_caller::ng::types::ContigId;
use pop_var_caller::regions::ContigBounds;

fn main() -> ExitCode {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        eprintln!(
            "usage: dhat_ng_open_files <reference.fa> <regions.bed> <n> <sample.cram> [sample ...]\n\
             opens the first <n> of the given files as one cohort and streams every BED span \
             through each, so peak heap can be split into fixed and per-file terms."
        );
        return ExitCode::from(2);
    }

    let fasta = PathBuf::from(&args[0]);
    let bed = PathBuf::from(&args[1]);
    let wanted: usize = args[2].parse().expect("<n> is a file count");
    let paths: Vec<PathBuf> = args[3..].iter().take(wanted).map(PathBuf::from).collect();
    eprintln!("opening {} file(s)", paths.len());

    match run(&fasta, &bed, &paths) {
        Ok((reads, spans)) => {
            eprintln!("{reads} reads over {spans} spans");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &std::path::Path,
    bed: &std::path::Path,
    paths: &[PathBuf],
) -> Result<(u64, u64), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs: ContigList = info.contig_list();
    let reference = OpenReference::new(info);

    // One `SampleReads` per sample, exactly as the cohort dump builds them — so
    // the k files are open simultaneously, which is the condition being
    // measured.
    let read_groups = build_read_groups(paths)?;
    let samples: Vec<SampleReads> = read_groups
        .read_groups_per_sample()
        .iter()
        .map(|entry| {
            SampleReads::open(
                entry,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                true,
            )
        })
        .collect::<Result<_, _>>()?;

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|e| ContigBounds {
            name: &e.name,
            length: e.length as u32,
        })
        .collect();
    let spans = GenomeRegions::from_bed_path(bed, &bounds)?;

    // A reference accessor per query, as the real caller does: `RawRefSeq` impls
    // are stateful readers, so one shared accessor would need a lock the moment
    // queries run concurrently.
    #[expect(
        clippy::arc_with_non_send_sync,
        reason = "RawRefSeq is implemented for Arc only; this harness is single-threaded"
    )]
    let make_reference = || Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));

    // `PVC_OPEN_ONLY` stops here, with every file open and nothing decoded.
    // That splits the per-file cost in two: what an *open* file costs, and what
    // *decoding through* one costs — which are different bugs with different
    // fixes, and the profile cannot tell them apart on its own.
    if std::env::var_os("PVC_OPEN_ONLY").is_some() {
        eprintln!(
            "open-only: {} sample(s) open, nothing decoded",
            samples.len()
        );
        if let Some(handle) = verify {
            handle.join()?;
        }
        return Ok((0, 0));
    }

    // `PVC_SPAN_LIMIT` caps the spans walked, so "one span per file" can be compared against a
    // full walk: the first span is what opens each sample's cursor and fills what it holds.
    let span_limit: usize = std::env::var("PVC_SPAN_LIMIT")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(usize::MAX);

    // **One cursor per sample per chromosome, which is what this harness now measures.** Before
    // the alignment cursor this opened a reader per span from a per-file pool, and what an open
    // file *held* was that pool plus a decoded CRAM container. A cursor holds its own reader,
    // its reference accessor and the reads it has kept for the next span — so the question the
    // harness asks is unchanged ("what are k open, queried files holding?") while the answer's
    // composition is not, and a number from before the cursor is not comparable to one after.
    //
    // Minted lazily and rebuilt at a chromosome change, exactly as a generator does: a cursor
    // covers one chromosome and nothing in it survives a change of chromosome.
    let mut cursors: Vec<Option<(ContigId, _)>> = samples.iter().map(|_| None).collect();

    let mut total_reads = 0u64;
    let mut total_spans = 0u64;
    for region in spans.iter().take(span_limit) {
        total_spans += 1;
        for (sample, slot) in samples.iter().zip(cursors.iter_mut()) {
            if slot.as_ref().is_none_or(|(on, _)| *on != region.contig) {
                *slot = Some((region.contig, sample.cursor(region.contig, make_reference)?));
            }
            let cursor = &mut slot.as_mut().expect("just ensured").1;
            cursor.move_to_region(region)?;
            while let Some(read) = cursor.next_read() {
                let _ = read?;
                total_reads += 1;
            }
        }
    }

    // Live heap **while the files are still open**, which the global peak
    // cannot tell us: `t-gmax` is one instant, and if it falls during reference
    // loading it is blind to what the per-file reader pools retain afterwards.
    // Reported here, before `samples` is dropped, so the number is exactly "what
    // k open, queried files are holding".
    #[cfg(feature = "dhat-heap")]
    eprintln!(
        "live heap with {} file(s) open and queried: {} bytes",
        samples.len(),
        dhat::HeapStats::get().curr_bytes
    );

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok((total_reads, total_spans))
}
