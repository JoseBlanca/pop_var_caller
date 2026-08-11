//! How much of a chromosome a windowed catalog read actually touches.
//!
//! Reads the same catalog twice — once for a whole contig, once for a 1 Mb stretch of it —
//! and reports the wall clock and the row count of each. The point is the ratio: if the
//! file's page index is doing its job, asking for a hundredth of a chromosome costs far
//! less than asking for all of it.

use std::path::PathBuf;
use std::time::Instant;

use pop_var_caller::ng::reference_info::{ReferenceSource, read_reference_info};
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind};
use pop_var_caller::ng::repeat_catalog::{RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::regions::ContigBounds;

fn main() {
    let mut args = std::env::args().skip(1);
    let reference = PathBuf::from(
        args.next()
            .expect("usage: <reference.fa> <catalog.parquet> <bed>"),
    );
    let catalog_path = PathBuf::from(args.next().expect("a catalog"));
    let bed = PathBuf::from(args.next().expect("a bed"));

    let info = read_reference_info(ReferenceSource::Fasta {
        fasta: reference,
        fai: None,
    })
    .expect("the reference reads");
    let catalog =
        RepeatCatalog::open_checking_against_reference(&catalog_path, &info).expect("opens");
    let criteria = StrRepeatCriteria::default();

    let bounds: Vec<ContigBounds> = info
        .contigs
        .iter()
        .map(|c| ContigBounds {
            name: &c.name,
            length: c.length as u32,
        })
        .collect();
    let spans = GenomeRegions::from_bed_path(&bed, &bounds).expect("a readable bed");
    let contig = spans.iter().next().expect("a span").contig;

    let whole = Instant::now();
    let mut whole_regions = 0usize;
    let mut whole_loci = 0usize;
    for region in catalog
        .genome_segments(&criteria, Some(contig))
        .expect("servable")
    {
        let region = region.expect("a region");
        whole_regions += 1;
        if matches!(region.kind, RegionKind::SsrSegment(_)) {
            whole_loci += 1;
        }
    }
    let whole_seconds = whole.elapsed().as_secs_f64();

    let windowed = Instant::now();
    let mut window_regions = 0usize;
    let mut window_loci = 0usize;
    for region in catalog
        .genome_segments_in(&criteria, &spans)
        .expect("servable")
    {
        let region = region.expect("a region");
        window_regions += 1;
        if matches!(region.kind, RegionKind::SsrSegment(_)) {
            window_loci += 1;
        }
    }
    let window_seconds = windowed.elapsed().as_secs_f64();

    println!("whole_contig_seconds\t{whole_seconds:.3}");
    println!("whole_contig_regions\t{whole_regions}");
    println!("whole_contig_loci\t{whole_loci}");
    println!("window_seconds\t{window_seconds:.3}");
    println!("window_regions\t{window_regions}");
    println!("window_loci\t{window_loci}");
    println!("speedup\t{:.1}x", whole_seconds / window_seconds.max(1e-9));
}
