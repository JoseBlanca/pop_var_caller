//! **The reference runs of every STR bundle a calling run types** — where in each bundle's hull
//! stutter will be priced, and in which motif (`str_bundle_loci.md` §6.1).
//!
//! ```text
//! ng_str_bundle_runs <reference.fa> [regions.bed]
//! ```
//!
//! Types the ground the way `generate-psps` does, from the repeat catalog beside the FASTA at
//! the calling floors, and for each `ssr_bundle` region finds its runs from the hull's reference
//! bases with [`find_reference_runs`] — the function the bundle generator calls.
//!
//! Output: a `#` header, a bare column line, then one row per bundle — `chrom start end members
//! runs`, the hull 1-based inclusive, `members` as `motif` list, `runs` as `motif:offset:length`
//! joined by `,` (offsets from the hull's first base), `.` for a bundle with no run. A summary
//! goes to standard error: bundles by their number of runs, and runs by period.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use pop_var_caller::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::region_typing::reference_runs::find_reference_runs;
use pop_var_caller::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::regions::ContigBounds;
use pop_var_caller::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::types::GenomeRegion;

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::reference_check_from_env;

fn run(fasta: &Path, bed: Option<&PathBuf>) -> Result<(), Box<dyn std::error::Error>> {
    let cache = std::sync::Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        reference_check_from_env()?,
    )?;
    let contigs = info.contig_list();
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    let reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());

    let mut bounds: Vec<ContigBounds<'_>> = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        bounds.push(ContigBounds {
            name: &entry.name,
            length: u32::try_from(entry.length)
                .map_err(|_| format!("contig {} is longer than 4 Gb", entry.name))?,
        });
    }
    let ground: Vec<GenomeRegion> = match bed {
        Some(path) => GenomeRegions::from_bed_path(path, &bounds)?
            .iter()
            .collect(),
        None => GenomeRegions::whole_contigs(&bounds).iter().collect(),
    };
    let config = TypedRegionConfig::default();
    let criteria = StrRepeatCriteria::from(&config);
    let min_copies = &criteria.classification.min_copies;

    println!("## tool: ng_str_bundle_runs");
    println!("## reference: {}", fasta.display());
    println!("#chrom\tstart\tend\tmembers\truns");

    let mut bundles_by_run_count: BTreeMap<usize, u64> = BTreeMap::new();
    let mut runs_by_period: BTreeMap<usize, u64> = BTreeMap::new();
    let mut hull = Vec::new();
    let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(&ground))?;
    for region in walk.by_ref() {
        let region = region?;
        let RegionKind::SsrBundle { members } = &region.kind else {
            continue;
        };
        let span = region.region;
        reference.fetch_into(span.contig, span.start.get(), span.len(), &mut hull)?;
        let runs = find_reference_runs(&hull, min_copies);
        *bundles_by_run_count.entry(runs.len()).or_default() += 1;
        for run in &runs {
            *runs_by_period.entry(run.period()).or_default() += 1;
        }
        let members_cell: Vec<String> = members
            .iter()
            .map(|member| String::from_utf8_lossy(member.motif().as_bytes()).into_owned())
            .collect();
        let runs_cell: Vec<String> = runs
            .iter()
            .map(|run| {
                format!(
                    "{}:{}:{}",
                    String::from_utf8_lossy(run.motif().as_bytes()),
                    run.offset(),
                    run.length()
                )
            })
            .collect();
        println!(
            "{}\t{}\t{}\t{}\t{}",
            contigs.entries[span.contig.get() as usize].name,
            span.start.get(),
            span.end.get(),
            members_cell.join(","),
            if runs_cell.is_empty() {
                ".".to_string()
            } else {
                runs_cell.join(",")
            }
        );
    }

    let bundles: u64 = bundles_by_run_count.values().sum();
    eprintln!("bundles: {bundles}");
    for (count, bundles) in &bundles_by_run_count {
        eprintln!("  with {count} runs: {bundles}");
    }
    for (period, runs) in &runs_by_period {
        eprintln!("  runs of period {period}: {runs}");
    }
    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: ng_str_bundle_runs <reference.fa> [regions.bed]");
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bed = args.get(2).map(PathBuf::from);
    match run(&fasta, bed.as_ref()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
