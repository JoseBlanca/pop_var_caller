//! **What a calling run decides each stretch of ground is** — the same typed regions
//! `call-from-alignments` routes on, printed so they can be intersected with anything else.
//!
//! ```text
//! ng_typed_region_dump <reference.fa> [regions.bed] [catalog|calling]
//! ```
//!
//! A run cuts the ground into four kinds and sends each to its own locus generator: an
//! ordinary stretch (`generic`), one microsatellite tract (`ssr_locus`), a cluster of tracts
//! too close to separate (`ssr_bundle`), and a tandem array too long to type as callable
//! (`satellite`). **Only the first has a generator wired in today**, so everything printed
//! here as one of the other three is ground the caller builds no locus over and can say
//! nothing about. The run report gives that as a percentage of bases; this gives the spans,
//! which is what it takes to ask *which known variants fall in them*.
//!
//! The kinds come from the reference's repeat catalog, opened beside the FASTA — the same
//! file the run reads, so this answers the run's question rather than a similar one. Build it
//! first with `pop_var_caller_exp repeat-catalog --reference <reference.fa>`.
//!
//! **The third argument is which floors decide what an STR tract is**, and the two differ by
//! a lot. `catalog` (the default) is what `call-from-alignments` uses today: the floors the
//! catalog file was *stored* at, deliberately permissive so one file can serve any later
//! policy — a 5-copy mononucleotide is a tract. `calling` is ng's short-read calling policy,
//! whose copy floors were measured as the point where a repeat starts to stutter, and below
//! which `MinCopies::default`'s own documentation says the generic SNP/indel caller should
//! handle the tract. Running both and diffing says how much ground the run currently sends
//! down the STR path that the calling policy would leave on the generic one.
//!
//! Without a BED the whole genome is typed, which on a large reference is a large output;
//! with one, only the stretches overlapping it. A tract that straddles a BED edge is printed
//! whole, because half a tract is not a tract — that is the run's rule too, and the reason a
//! run's "bases spoken for" can exceed the bases asked for.
//!
//! Output: a `#` comment header, a bare TSV column line, then one row per typed region —
//! `chrom  start  end  kind  motif  period  copies  purity`, coordinates 1-based inclusive,
//! with `.` wherever a kind carries no such detail.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::GenomeRegion;
use pop_var_caller::regions::ContigBounds;

/// The shared "should this run check the reference?" rule — see the module's own docs for
/// why the default is to check even in a tool that is re-run constantly.
#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::reference_check_from_env;

fn run(
    fasta: &Path,
    bed: Option<&PathBuf>,
    floors: Floors,
) -> Result<(), Box<dyn std::error::Error>> {
    let reference_check = reference_check_from_env()?;
    let cache = std::sync::Arc::new(ReferenceInfoCache::new());
    let (info, verify) =
        read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf(), reference_check)?;
    let contigs = info.contig_list();
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;

    // The ground to type, resolved the way the run resolves it: BED spans against the
    // contig table, so an unknown contig name is a refusal rather than an empty answer;
    // whole contigs when no BED is given.
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

    let criteria = match floors {
        Floors::Catalog => StrRepeatCriteria::default(),
        Floors::Calling => StrRepeatCriteria::from(&TypedRegionConfig::default()),
    };
    println!("## tool: ng_typed_region_dump");
    println!("## reference: {}", fasta.display());
    println!("## floors: {}", floors.name());
    let floors_by_period: Vec<u32> = (criteria.classification.periods.min()
        ..=criteria.classification.periods.max())
        .map(|period| criteria.classification.min_copies.for_period(period))
        .collect();
    println!(
        "## periods: {}..={}  min_copies: {:?}  max_str_len_bp: {}",
        criteria.classification.periods.min(),
        criteria.classification.periods.max(),
        floors_by_period,
        criteria.max_str_len_bp.get(),
    );
    match bed {
        Some(path) => println!("## regions: {}", path.display()),
        None => println!("## regions: whole genome"),
    }
    println!("#chrom\tstart\tend\tkind\tmotif\tperiod\tcopies\tpurity");

    let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(&ground))?;
    for region in walk.by_ref() {
        let region = region?;
        let name = &contigs.entries[region.region.contig.get() as usize].name;
        let (kind, motif, period, copies, purity) = match &region.kind {
            RegionKind::Generic => (
                "generic",
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
            ),
            RegionKind::Satellite => (
                "satellite",
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
            ),
            RegionKind::SsrBundle { tracts } => (
                "ssr_bundle",
                ".".to_string(),
                ".".to_string(),
                format!("{}", tracts.len()),
                ".".to_string(),
            ),
            RegionKind::SsrSegment(segment) => (
                "ssr_locus",
                String::from_utf8_lossy(segment.motif().as_bytes()).into_owned(),
                format!("{}", segment.period()),
                // Copies as the tract measures them: whole tract over one unit.
                format!(
                    "{:.1}",
                    segment.tract_len() as f64 / segment.period() as f64
                ),
                format!("{:.2}", segment.purity_fraction()),
            ),
        };
        println!(
            "{name}\t{}\t{}\t{kind}\t{motif}\t{period}\t{copies}\t{purity}",
            region.region.start.get(),
            region.region.end.get(),
        );
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

/// Which copy-number floors decide what counts as an STR tract.
#[derive(Clone, Copy)]
enum Floors {
    /// The floors the catalog file was stored at — what `call-from-alignments` routes on
    /// today.
    Catalog,
    /// ng's short-read calling policy: the copy numbers at which a repeat starts to stutter.
    Calling,
}

impl Floors {
    fn name(self) -> &'static str {
        match self {
            Floors::Catalog => "catalog (what call-from-alignments uses)",
            Floors::Calling => "calling (SsrSegmentCriteria::default)",
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "usage: ng_typed_region_dump <reference.fa> [regions.bed] [catalog|calling]\n\
             prints the typed regions a calling run routes on, one per line."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bed = args.get(2).map(PathBuf::from);
    let floors = match args.get(3).map(String::as_str) {
        None | Some("catalog") => Floors::Catalog,
        Some("calling") => Floors::Calling,
        Some(other) => {
            eprintln!("unknown floors '{other}': expected 'catalog' or 'calling'");
            return ExitCode::from(2);
        }
    };
    match run(&fasta, bed.as_ref(), floors) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
