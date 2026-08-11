//! **STR delimiter divergence — the per-read view.** Step 1 of the "which aligner is right?"
//! investigation: dump the individual reads where the two tract delimiters *disagree*, with each
//! aligner's measurement and the read's own sequence, so a human (or a slow oracle) can adjudicate
//! who read the repeat correctly.
//!
//! ```text
//! ng_ssr_divergent_reads <reference.fa> <sample.bam|cram> [contig ...]
//! ```
//!
//! Both delimiters run on the **same** reads (one region-typing walk, the reservoir kept-order is
//! identical across aligners), via [`SsrGenerator::delimit_segment_reads`] — the per-read view
//! behind the tally. A read is *divergent* when the two aligners give it a different measurement:
//! different tract bases, a different coverage class (complete / partial-left / partial-right), or
//! one anchors it and the other does not. Each divergent read becomes one TSV row.
//!
//! Output: a `#` counts header, a bare TSV column line, then one row per divergent read — the locus
//! (coords, motif, reference tract), the read name, each aligner's `class:bases`, and the read
//! sequence. A per-category tally is written to stderr.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::locus_generation::LocusGenerator;
use pop_var_caller::ng::locus_generation::ReadWitness;
use pop_var_caller::ng::locus_generation::ssr::{
    RepeatDelimiter, SegmentDelimitations, SsrGenerator, SsrGeneratorConfig,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig};
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::{Bp, ContigId};

#[path = "shared/catalog_regions.rs"]
mod catalog_regions;

/// A read's measurement under one aligner, normalised for comparison: the coverage class and the
/// measured tract bases (the partial *reach* integer is dropped — it is not the tract measurement).
/// `None` = the read anchored no border / was gated out.
type Measurement = Option<(&'static str, Vec<u8>)>;

fn measurement(obs: &Option<(ReadWitness, Vec<u8>)>) -> Measurement {
    obs.as_ref().map(|(cov, bases)| {
        // The side is a derivation since the reshape: a run flush with the left border is a
        // prefix. On this path a partial always anchors one border, so "not flush left" is
        // exactly "flush right" and the locus length is not needed to tell them apart — an
        // interior run would be mislabelled here, but the STR generator cannot mint one.
        //
        // Destructured rather than guarded on `_`, so a future variant is a compile error.
        let class = match cov {
            ReadWitness::Complete => "complete",
            run @ ReadWitness::Partial { .. } if run.is_flush_left() => "partialL",
            ReadWitness::Partial { .. } => "partialR",
        };
        (class, bases.clone())
    })
}

/// Render a measurement as `class:bases` (bases empty for a partial of length 0), or `none`.
fn render(m: &Measurement) -> String {
    match m {
        Some((class, bases)) => format!("{class}:{}", String::from_utf8_lossy(bases)),
        None => "none".to_string(),
    }
}

/// The category of a divergence — a coarse triage of *why* the two aligners disagree.
fn category(flat: &Measurement, unit: &Measurement) -> &'static str {
    match (flat, unit) {
        (None, Some(_)) | (Some(_), None) => "anchor_vs_none",
        (Some((cf, bf)), Some((cu, bu))) => {
            if cf != cu {
                "class_differ"
            } else if bf.len() != bu.len() {
                "length_differ"
            } else {
                "same_length_diff_bases" // an interior substitution / interruption call
            }
        }
        (None, None) => "agree", // unreachable — filtered out before here
    }
}

fn make_generator<A: RepeatDelimiter>(
    fasta: &Path,
    contigs: &ContigList,
    aligner: A,
    config: SsrGeneratorConfig,
    bundle_threshold: Bp,
) -> Result<SsrGenerator<WindowedRefSeq, A>, Box<dyn std::error::Error>> {
    Ok(SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        aligner,
        config,
        bundle_threshold,
    )?)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_divergent_reads <reference.fa> <sample.bam|cram> [contig ...]\n\
             dumps the individual reads where algorithm 3 (flat-gap) and algorithm 4 (unit-slip) \
             give a different tract measurement, with each aligner's call and the read sequence."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter: Vec<String> = args[3..].to_vec();

    match run(&fasta, &[bam], &contig_filter) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: &[String],
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs = info.contig_list();
    // **The typed regions come from the catalog beside the reference**, checked against what
    // the pass just reported. No catalog, no run: the error names the command that writes one.
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    // One reference for every file this run opens — and so one copy of the bases.
    let reference = OpenReference::new(info);
    let sample =
        SampleReads::open_only_sample(bams, &reference, ReadFilterConfig::default(), true)?;

    let walk_config = TypedRegionConfig::default();
    let criteria = StrRepeatCriteria::from(&walk_config);
    let bundle_threshold = Bp(walk_config.criteria.bundle_threshold);
    let config = SsrGeneratorConfig::default();
    let emission = PerQualityEmission::new();

    let mut gen_flat = make_generator(
        fasta,
        &contigs,
        SsrFlatGapAligner::new(emission),
        config.clone(),
        bundle_threshold,
    )?;
    let mut gen_unit = make_generator(
        fasta,
        &contigs,
        SsrUnitSlipAligner::new(emission),
        config,
        bundle_threshold,
    )?;

    let mut rows = String::new();
    let mut categories: BTreeMap<&'static str, u64> = BTreeMap::new();
    let mut divergent_reads: u64 = 0;
    let mut divergent_loci: u64 = 0;
    let mut total_reads: u64 = 0;
    let mut order_mismatches: u64 = 0;

    for (index, entry) in contigs.entries.iter().enumerate() {
        if !contig_filter.is_empty() && !contig_filter.iter().any(|name| name == &entry.name) {
            continue;
        }
        let this_contig = [catalog_regions::whole_contig(
            ContigId(index as u32),
            entry.length,
        )];
        let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(&this_contig))?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            gen_flat.begin_segment(region.region);
            gen_unit.begin_segment(region.region);
            let flat = gen_flat.delimit_segment_reads(segment, &sample)?;
            let unit = gen_unit.delimit_segment_reads(segment, &sample)?;

            let locus_divergent = fold_locus(
                segment,
                &flat,
                &unit,
                &mut rows,
                &mut categories,
                &mut divergent_reads,
                &mut total_reads,
                &mut order_mismatches,
            );
            if locus_divergent {
                divergent_loci += 1;
            }
        }
    }

    if let Some(handle) = verify {
        handle.join()?;
    }

    print!(
        "# divergent_reads={divergent_reads} divergent_loci={divergent_loci} total_reads={total_reads}"
    );
    for (cat, n) in &categories {
        print!(" {cat}={n}");
    }
    println!();
    println!(
        "contig\tstart\tend\tperiod\tref_len\tmotif\tref_tract\tleft_flank\tright_flank\tqname\tcategory\tflat_gap\tunit_slip\tread_seq"
    );
    print!("{rows}");

    eprintln!("\n--- divergence summary (flat_gap vs unit_slip) ---");
    eprintln!("reads compared:   {total_reads}");
    eprintln!(
        "divergent reads:  {divergent_reads} ({:.2}%) over {divergent_loci} loci",
        100.0 * divergent_reads as f64 / total_reads.max(1) as f64
    );
    for (cat, n) in &categories {
        eprintln!("  {cat:<24} {n}");
    }
    if order_mismatches > 0 {
        eprintln!(
            "WARNING: {order_mismatches} reads where the two aligners' kept-order qnames disagreed \
             (skipped) — the index-zip assumption did not hold."
        );
    }
    Ok(())
}

/// Compare the two aligners' per-read measurements for one segment, emit a row per divergent read,
/// and update the tallies. Returns whether this locus had any divergence.
#[allow(clippy::too_many_arguments)]
fn fold_locus(
    segment: &SsrSegment,
    flat: &SegmentDelimitations,
    unit: &SegmentDelimitations,
    rows: &mut String,
    categories: &mut BTreeMap<&'static str, u64>,
    divergent_reads: &mut u64,
    total_reads: &mut u64,
    order_mismatches: &mut u64,
) -> bool {
    use std::fmt::Write as _;
    let ref_tract = String::from_utf8_lossy(&flat.reference_tract);
    let left_flank = String::from_utf8_lossy(&flat.left_flank);
    let right_flank = String::from_utf8_lossy(&flat.right_flank);
    let motif = String::from_utf8_lossy(segment.motif().as_bytes()).to_string();
    let period = motif.len();
    let ref_len = flat.reference_tract.len();

    let mut any = false;
    // Same config → identical reservoir kept order, so the two vectors align by index.
    for (f, u) in flat.reads.iter().zip(unit.reads.iter()) {
        if f.qname != u.qname {
            *order_mismatches += 1;
            continue;
        }
        *total_reads += 1;
        let mf = measurement(&f.observation);
        let mu = measurement(&u.observation);
        if mf == mu {
            continue;
        }
        any = true;
        *divergent_reads += 1;
        *categories.entry(category(&mf, &mu)).or_insert(0) += 1;
        let _ = writeln!(
            rows,
            "{contig}\t{start}\t{end}\t{period}\t{ref_len}\t{motif}\t{ref_tract}\t{left_flank}\t{right_flank}\t{qname}\t{cat}\t{flat}\t{unit}\t{seq}",
            contig = segment.chrom(),
            start = segment.start(),
            end = segment.end(),
            qname = String::from_utf8_lossy(&f.qname),
            cat = category(&mf, &mu),
            flat = render(&mf),
            unit = render(&mu),
            seq = String::from_utf8_lossy(&f.read_seq),
        );
    }
    any
}
