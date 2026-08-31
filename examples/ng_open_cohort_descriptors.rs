//! **What does opening a real cohort cost in file descriptors, and does a 63-sample cohort open
//! at all?**
//!
//! Run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example ng_open_cohort_descriptors -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <sample.cram> [sample ...]
//! ```
//!
//! ## The two questions
//!
//! **First, that the object opens a cohort of the committed size.** `AlignedFilesVariantCaller`
//! is tested against fixtures of at most three samples, and the range this caller commits to is
//! one sample to several thousand (`doc/devel/ng/spec/design_principles.md` §0). Sixty-three
//! tomato accessions is the largest cohort on this machine; opening it is the first evidence
//! that the construction path — every file opened, six refusals evaluated, every sample's contig
//! checksums compared — survives past a handful of samples.
//!
//! **Second, what a run actually spends per alignment file.** `src/ng/run/callers.rs`'s
//! `DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS` is 2, and until this ran it was spec §7.1a's estimate
//! ("a CRAM and its index are two descriptors each") rather than anything counted. So this counts
//! `/proc/self/fd` at three points and reports the per-file slope between them:
//!
//! - after the reference is open and verified, before any sample is;
//! - after `AlignedFilesVariantCaller::open` — every sample's files open, nothing decoded;
//! - after one cursor per sample is made and moved onto a region — the shape a run holds while
//!   it walks.
//!
//! **What it found, on the 63 tomato accessions over `benchmarks/tomato1/regions.bed`:** 3
//! descriptors before, **4** with all 63 files open, **130** with a cursor on each. So the
//! constant is right at 2 a file and the spec's reason for it is wrong twice over — 63 open files
//! cost one descriptor between them, because the index is parsed into memory and an open
//! `AlignmentFile` keeps no handle; what costs 2 is a **cursor**, one for the file's reader and
//! one for the per-file reference accessor, which opens the FASTA.
//!
//! The third point is the one the constant has to cover, because a run holds a cursor per file
//! for the whole walk (spec §5.1). It is also the number Milestone E can change: several callers
//! in flight is a question nobody has counted, and this is the probe to re-run there.
//!
//! **`/proc/self/fd` is Linux**, so this reports the count only there; on other platforms it
//! still opens the cohort and says the count was unavailable. That is the ordinary case here —
//! every build in this project runs in the container.
//!
//! ## The catalog is given by path, not found beside the reference
//!
//! `RepeatCatalog::open_beside_reference` is the convention, and it cannot be used here: the
//! benchmark reference lives under `$HOME/genomes`, which the dev container mounts read-only, so
//! no catalog can be written beside it. The catalog for that reference is inside the project
//! tree instead, and this takes its path.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::calling::allele_candidates::CandidateSelectionConfig;
use pop_var_caller::ng::calling::inference::CallingLoopConfig;
use pop_var_caller::ng::calling::parameters_file::DeclaredInbreeding;
use pop_var_caller::ng::calling::run_parameters::RunParameters;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::GenomeRegions;
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::run::{
    AlignedFilesVariantCaller, AlignmentInputs, MergeParameters, Segmentation,
};
use pop_var_caller::ng::types::Ploidy;
use pop_var_caller::regions::ContigBounds;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        eprintln!(
            "usage: ng_open_cohort_descriptors <reference.fa> <catalog.parquet> <regions.bed> \
             <sample.cram> [sample ...]\n\
             opens the given files as one cohort through AlignedFilesVariantCaller and reports \
             how many descriptors that costs per alignment file."
        );
        return ExitCode::from(2);
    }

    let fasta = PathBuf::from(&args[0]);
    let catalog = PathBuf::from(&args[1]);
    let bed = PathBuf::from(&args[2]);
    let paths: Vec<PathBuf> = args[3..].iter().map(PathBuf::from).collect();

    match run(&fasta, &catalog, &bed, &paths) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

/// How many descriptors this process holds, or `None` where the kernel does not say.
///
/// **Counted rather than estimated**, which is the whole point of this probe. The listing
/// itself holds one open directory descriptor while it runs; that is subtracted, so the number
/// is what the process holds when nothing is looking.
fn open_descriptors() -> Option<usize> {
    let entries = std::fs::read_dir("/proc/self/fd").ok()?;
    let counted = entries.filter(Result::is_ok).count();
    Some(counted.saturating_sub(1))
}

fn report(stage: &str, count: Option<usize>) {
    match count {
        Some(count) => println!("{stage}: {count} open descriptors"),
        None => println!("{stage}: /proc/self/fd is not readable on this platform"),
    }
}

fn run(
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    paths: &[PathBuf],
) -> Result<(), Box<dyn std::error::Error>> {
    // The reference, **with its FASTA read to the end**: the assembly check and the catalog
    // check both compare checksums, and a reference whose background verification has not been
    // joined carries none (`callers.rs`, `AlignmentInputs::reference_with_checksums`).
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let with_checksums = match verify {
        Some(handle) => handle.join()?,
        None => Arc::clone(&info),
    };
    let contigs: ContigList = info.contig_list();
    let reference = OpenReference::new(info);

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: entry.length as u32,
        })
        .collect();
    let analysed = GenomeRegions::from_bed_path(bed, &bounds)?;

    let criteria = StrRepeatCriteria::default();
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog.genome_segments(&criteria, ReadScope::Regions(&spans))?;
    let segmentation = Segmentation::build(
        segments,
        analysed,
        catalog.header().clone(),
        criteria,
        catalog_path.to_path_buf(),
    )?;
    println!(
        "segmentation: {} segments over {} analysed region(s)",
        segmentation.segments().len(),
        segmentation.analysed_regions().len(),
    );

    let read_groups = build_read_groups(paths)?;
    let parameters = RunParameters::of_defaults(
        &read_groups,
        Ploidy::try_new(2).expect("a diploid"),
        &DeclaredInbreeding::nothing_said(),
    );

    let before_opening = open_descriptors();
    report("reference open, no sample open", before_opening);

    let caller = AlignedFilesVariantCaller::open(
        AlignmentInputs {
            read_groups: &read_groups,
            reference: &reference,
            read_filters: ReadFilterConfig::default(),
            build_index_if_missing: false,
            reference_with_checksums: &with_checksums,
        },
        segmentation,
        parameters,
        CallingLoopConfig::DEFAULT
            .validate()
            .expect("the shipped calling-loop settings are runnable"),
        CandidateSelectionConfig::DEFAULT,
        MergeParameters::DEFAULT,
    )?;

    let alignment_files: usize = caller.samples().map(|sample| sample.file_count()).sum();
    println!(
        "cohort: {} samples over {} alignment files; assembly check: {:?}",
        caller.sample_count(),
        alignment_files,
        caller.assembly_check(),
    );

    let after_opening = open_descriptors();
    report("every sample open, nothing decoded", after_opening);

    // **One cursor per sample, on the first analysed region's contig, moved onto that region** —
    // the shape a run holds while it walks (spec §5.1). A cursor holds one reader and one
    // reference accessor *per file*, so this is where a per-file descriptor would appear if
    // there is one.
    let first = caller
        .segmentation()
        .analysed_regions()
        .first()
        .copied()
        .ok_or("the BED gave no regions to walk")?;
    let shared_contigs = Arc::new(contigs.clone());
    let shared_index = WindowedRefSeq::read_index(fasta)?;
    let make_reference = || {
        WindowedRefSeq::with_shared_index(
            fasta.to_path_buf(),
            Arc::clone(&shared_contigs),
            Arc::clone(&shared_index),
        )
    };
    let mut cursors = Vec::with_capacity(caller.sample_count());
    for sample in caller.samples() {
        let mut cursor = sample.cursor(first.contig, make_reference)?;
        cursor.move_to_region(first)?;
        while let Some(read) = cursor.next_read() {
            read?;
        }
        cursors.push(cursor);
    }

    let while_walking = open_descriptors();
    report(
        "one cursor per sample, walked over one region",
        while_walking,
    );

    // **And the two accessors a locus generator holds per sample**, which the cursor count above
    // does not include and which a run holds for its whole walk. `PileupGenerator` keeps one for
    // the walk's own REF fetches and the read preparer keeps a second; both are
    // `WindowedRefSeq`s, and each opens a reader on the FASTA at its first fetch and keeps it.
    // They are per *sample*, not per file, so they are a second term in the arithmetic and not a
    // correction to the first.
    let mut held_by_generators = Vec::with_capacity(caller.sample_count() * 2);
    for _ in 0..caller.sample_count() {
        for _ in 0..2 {
            let accessor = WindowedRefSeq::with_shared_index(
                fasta.to_path_buf(),
                Arc::clone(&shared_contigs),
                Arc::clone(&shared_index),
            );
            // Fetch, because the reader is opened lazily: an accessor nobody has asked for holds
            // nothing, and counting before the first fetch would report zero.
            accessor
                .fetch(first.contig, first.start.get(), 1)
                .map_err(|error| format!("the reference fetches: {error}"))?;
            held_by_generators.push(accessor);
        }
    }

    let while_generating = open_descriptors();
    report(
        "and the two accessors a generator holds per sample",
        while_generating,
    );

    if let (Some(before), Some(open), Some(walking), Some(generating)) = (
        before_opening,
        after_opening,
        while_walking,
        while_generating,
    ) {
        let files = alignment_files as f64;
        let samples = caller.sample_count() as f64;
        println!(
            "per alignment file: {:.2} descriptors when open, {:.2} with a cursor on it",
            (open.saturating_sub(before)) as f64 / files,
            (walking.saturating_sub(before)) as f64 / files,
        );
        println!(
            "per sample, on top of that: {:.2} descriptors for the generator's own accessors",
            (generating.saturating_sub(walking)) as f64 / samples,
        );
        println!(
            "a walking run therefore holds {} descriptors for {alignment_files} files over {} samples",
            generating.saturating_sub(before),
            caller.sample_count(),
        );
    }

    drop(cursors);
    Ok(())
}
