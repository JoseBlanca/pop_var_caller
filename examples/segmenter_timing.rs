//! Experiment for the 2026-10-09 segmenter performance review: time each part of typing a
//! whole genome from its stored repeat catalog, with no output, so the typing's own cost is
//! separated from printing and from checking the reference.
//!
//! ```text
//! segmenter_timing <reference.fa> [catalog|calling] [phases] [regions.bed]
//! ```
//!
//! `phases` is a comma list of `segments`, `loci`, `rows`, `census` (default
//! `segments,loci,rows`). `census` runs what a psp-producing run does to plan its census:
//! segments + loci again (`select_kept_loci`'s two passes).

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use pop_var_caller::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::regions::ContigBounds;
use pop_var_caller::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::types::GenomeRegion;

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::reference_check_from_env;

fn run(
    fasta: &Path,
    calling: bool,
    phases: &str,
    bed: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    let t = Instant::now();
    let reference_check = reference_check_from_env()?;
    let cache = std::sync::Arc::new(ReferenceInfoCache::new());
    let (info, verify) =
        read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf(), reference_check)?;
    let contigs = info.contig_list();
    eprintln!("reference_open_s\t{:.3}", t.elapsed().as_secs_f64());

    let t = Instant::now();
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    eprintln!("catalog_open_s\t{:.3}", t.elapsed().as_secs_f64());

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|e| ContigBounds {
            name: &e.name,
            length: u32::try_from(e.length).expect("contig < 4 Gb"),
        })
        .collect();
    let ground: Vec<GenomeRegion> = match bed {
        Some(path) => GenomeRegions::from_bed_path(Path::new(path), &bounds)?.iter().collect(),
        None => GenomeRegions::whole_contigs(&bounds).iter().collect(),
    };
    eprintln!("spans\t{}", ground.len());

    let criteria = if calling {
        StrRepeatCriteria::from(&TypedRegionConfig::default())
    } else {
        StrRepeatCriteria::default()
    };

    for phase in phases.split(',') {
        let t = Instant::now();
        match phase {
            "segments" => {
                let (mut n, mut generic, mut loci) = (0u64, 0u64, 0u64);
                for region in catalog.genome_segments(&criteria, ReadScope::Regions(&ground))? {
                    let region = region?;
                    n += 1;
                    match region.kind {
                        RegionKind::Generic => generic += 1,
                        RegionKind::SsrSegment(_) => loci += 1,
                        _ => {}
                    }
                }
                eprintln!(
                    "segments_s\t{:.3}\tregions={n}\tgeneric={generic}\tloci={loci}",
                    t.elapsed().as_secs_f64()
                );
            }
            "segments_whole" => {
                let mut n = 0u64;
                for region in catalog.genome_segments(&criteria, ReadScope::WholeReference)? {
                    region?;
                    n += 1;
                }
                eprintln!(
                    "segments_whole_s\t{:.3}\tregions={n}",
                    t.elapsed().as_secs_f64()
                );
            }
            "loci" => {
                let mut n = 0u64;
                for locus in catalog.str_loci(&criteria, ReadScope::Regions(&ground))? {
                    locus?;
                    n += 1;
                }
                eprintln!("loci_s\t{:.3}\tloci={n}", t.elapsed().as_secs_f64());
            }
            "rows" => {
                let rows = catalog.repeats(ReadScope::WholeReference)?;
                eprintln!(
                    "rows_s\t{:.3}\trows={}",
                    t.elapsed().as_secs_f64(),
                    rows.len()
                );
            }
            "census" => {
                let mut n = 0u64;
                for region in catalog.genome_segments(&criteria, ReadScope::Regions(&ground))? {
                    region?;
                    n += 1;
                }
                let (_, sample) = catalog.sample_loci_per_stratum(
                    &criteria,
                    ReadScope::Regions(&ground),
                    200,
                    7,
                )?;
                let _ = sample;
                eprintln!("census_s\t{:.3}\tregions={n}", t.elapsed().as_secs_f64());
            }
            other => return Err(format!("unknown phase {other}").into()),
        }
    }

    if let Some(handle) = verify {
        let t = Instant::now();
        handle.join()?;
        eprintln!("reference_verify_join_s\t{:.3}", t.elapsed().as_secs_f64());
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: segmenter_timing <reference.fa> [catalog|calling] [phases]");
        return ExitCode::FAILURE;
    }
    let calling = args.get(2).map(|s| s == "calling").unwrap_or(false);
    let phases = args.get(3).map(String::as_str).unwrap_or("segments,loci,rows");
    match run(Path::new(&args[1]), calling, phases, args.get(4).map(String::as_str)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
