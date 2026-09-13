//! **How wide are the records a census reads, and what does its per-position rule do to them?**
//!
//! A calling run's locus is variable-width: one sample's deletion chains the positions it covers
//! into one cohort locus, and every sample is then genotyped over that whole stretch
//! ([`cohort_merge.md`](../doc/devel/ng/spec/cohort_merge.md) §3.2). A census works the other way
//! round — its sites are single reference positions chosen before a read is seen — so a record
//! wider than one base is taken apart by
//! [`CensusWriter::add_generic`](../src/parameter_estimation/joint/census.rs), which records
//! that record's depth at every position it covers and its **alleles at none of them**, because
//! the observations describe the whole span and not any one base of it.
//!
//! This counts what that rule reaches.
//!
//! ```text
//! cargo run --release --example ng_census_locus_spans -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <a.psp> [--keep-one-in K]
//! ```
//!
//! Everything here is read from record heads, so no body is decoded: a head carries the region a
//! record covers and how many of its reads disagreed with the reference, which is both halves of
//! the question.

use std::error::Error;
use std::path::Path;
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::cli::run_ground::{self, GroundRequest, RepeatRouting};
use pop_var_caller::parameter_estimation::joint::loci::UnambiguousRuns;
use pop_var_caller::psp::PspReader;
use pop_var_caller::reference_info::{ReferenceCheck, read_reference_observing_or_creating_fai};
use pop_var_caller::region_typing::DEFAULT_MAX_STR_LEN;
use pop_var_caller::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use pop_var_caller::repeat_catalog::RepeatCatalog;
use pop_var_caller::run::{CensusPlan, CensusSelection};
use pop_var_caller::types::{GenomePosition, GenomeRegion};

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut keep_one_in: u64 = 1;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--keep-one-in" => match args.next().and_then(|v| v.parse().ok()) {
                Some(k) if k >= 1 => keep_one_in = k,
                _ => {
                    eprintln!("--keep-one-in wants a count of at least 1");
                    return ExitCode::from(2);
                }
            },
            other => positional.push(other.to_string()),
        }
    }
    let [fasta, catalog, bed, psp] = positional.as_slice() else {
        eprintln!(
            "usage: ng_census_locus_spans <reference.fa> <catalog.parquet> <regions.bed> \
             <a.psp> [--keep-one-in K]"
        );
        return ExitCode::from(2);
    };
    match run(
        Path::new(fasta),
        Path::new(catalog),
        Path::new(bed),
        Path::new(psp),
        keep_one_in,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            let mut cause = error.source();
            while let Some(next) = cause {
                eprintln!("  caused by: {next}");
                cause = next.source();
            }
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    psp_path: &Path,
    keep_one_in: u64,
) -> Result<(), Box<dyn Error>> {
    let mut callable = UnambiguousRuns::default();
    let with_checksums = Arc::new(read_reference_observing_or_creating_fai(
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
        &mut callable,
    )?);
    let unambiguous = callable.into_selectable()?;
    let contigs = with_checksums.contig_list();
    let ground = GroundRequest {
        reference: fasta,
        catalog: Some(catalog_path),
        regions: Some(bed),
        routing: RepeatRouting {
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        },
    };
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let plan = CensusPlan::of_run(
        CensusSelection {
            seed: CensusSelection::SHIPPED.seed,
            generic_target: (analysed_bases / keep_one_in).max(1),
            ssr_cap: CensusSelection::SHIPPED.ssr_cap,
        },
        &catalog,
        &analysed,
        &unambiguous,
        &with_checksums,
        &segmentation.inputs().repeat_tract_criteria,
    )?;
    drop(catalog);
    let kept = plan.loci.generic();

    let mut reader = PspReader::open(psp_path)?;
    let mut records = 0_u64;
    let mut wide_records = 0_u64;
    let mut widest = 0_u64;
    let mut wide_with_a_disagreeing_read = 0_u64;
    // Kept census positions, by the width of the record covering them.
    let mut kept_under_one_base = 0_u64;
    let mut kept_under_a_wide_record = 0_u64;
    let mut kept_under_a_wide_record_with_a_disagreeing_read = 0_u64;
    // The same counts for one-base records, so the two are read against each other.
    let mut kept_under_one_base_with_a_disagreeing_read = 0_u64;

    for streamed in reader.records_where(|_| false)? {
        let head = streamed?.head;
        records += 1;
        let span = head.region.len();
        let disagreeing = head.non_reference_reads > 0;
        let covered = kept_positions_in(kept, head.region);
        if span > 1 {
            wide_records += 1;
            widest = widest.max(span);
            if disagreeing {
                wide_with_a_disagreeing_read += 1;
            }
            kept_under_a_wide_record += covered;
            if disagreeing {
                kept_under_a_wide_record_with_a_disagreeing_read += covered;
            }
        } else {
            kept_under_one_base += covered;
            if disagreeing {
                kept_under_one_base_with_a_disagreeing_read += covered;
            }
        }
    }

    println!("# psp: {}", psp_path.display());
    println!("# analysed bases: {analysed_bases}");
    println!("# census positions kept: {}", kept.len());
    println!("records: {records}");
    println!(
        "records wider than one base: {wide_records} ({:.2}% of records), widest {widest} bases",
        100.0 * wide_records as f64 / records.max(1) as f64
    );
    println!(
        "  of those, carrying at least one read that disagrees with the reference: {wide_with_a_disagreeing_read}"
    );
    println!("kept positions under a one-base record: {kept_under_one_base}");
    println!(
        "  of those, in a record with a disagreeing read: {kept_under_one_base_with_a_disagreeing_read}"
    );
    println!(
        "kept positions under a record wider than one base: {kept_under_a_wide_record} \
         ({:.2}% of kept positions the walk reached)",
        100.0 * kept_under_a_wide_record as f64
            / (kept_under_one_base + kept_under_a_wide_record).max(1) as f64
    );
    println!(
        "  of those, in a record with a disagreeing read — depth recorded, alleles dropped: \
         {kept_under_a_wide_record_with_a_disagreeing_read}"
    );
    println!(
        "so of every kept position the census sees carrying non-reference evidence, \
         {} of {} have that evidence dropped",
        kept_under_a_wide_record_with_a_disagreeing_read,
        kept_under_one_base_with_a_disagreeing_read
            + kept_under_a_wide_record_with_a_disagreeing_read,
    );
    Ok(())
}

/// How many kept census positions fall inside `region`.
fn kept_positions_in(kept: &[GenomePosition], region: GenomeRegion) -> u64 {
    let first = kept.partition_point(|position| {
        (position.contig.get(), position.position.get()) < (region.contig.get(), region.start.get())
    });
    let mut count = 0_u64;
    for position in &kept[first..] {
        if position.contig != region.contig || position.position.get() > region.end.get() {
            break;
        }
        count += 1;
    }
    count
}
