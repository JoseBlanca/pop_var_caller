//! **Real-data validation of the candidate delimiters.** The synthetic scorecard rewards fixing the
//! fabricated-complete failure, but it does not measure whether a candidate *preserves* the genuine
//! partial observations the ng design keeps (censored lower bounds that carry per-position depth).
//! This checks that on real reads: it runs unit_slip (the baseline) alongside each candidate and
//! tallies the per-read class transition, binned by reference tract length.
//!
//! The failure mode to catch: a candidate that converts real one-anchored partials into `none`
//! (destroying evidence) rather than only fixing over-claimed completes. A good candidate demotes
//! completes → partials concentrated at long tracts, and leaves genuine partials as partials.
//!
//! ```text
//! ng_ssr_anchor_firm_validate <reference.fa> <sample.bam|cram> [contig ...]
//! ```

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_anchor_firm::SsrAnchorFirmAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::alignment::ssr_noise_robust::SsrNoiseRobustAligner;
use pop_var_caller::ng::alignment::ssr_robust_indel::SsrRobustIndelAligner;
use pop_var_caller::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner;
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
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId};

fn class(obs: &Option<(ReadWitness, Vec<u8>)>) -> &'static str {
    match obs {
        None => "none",
        Some((ReadWitness::Complete, _)) => "complete",
        Some((ReadWitness::Partial { .. }, _)) => "partial",
    }
}

fn band(len: usize) -> &'static str {
    match len {
        0..=14 => "<15",
        15..=24 => "15-24",
        25..=39 => "25-39",
        40..=59 => "40-59",
        _ => "60+",
    }
}

#[derive(Default)]
struct Tally {
    transitions: BTreeMap<(&'static str, &'static str), u64>,
    complete_by_band: BTreeMap<&'static str, u64>,
    demote_by_band: BTreeMap<&'static str, u64>,
    partial_lost: u64, // unit=partial -> candidate=none (destroyed evidence)
    partial_kept: u64, // unit=partial -> candidate=partial
    promotions: u64, // candidate=complete while unit != complete (should be ~0 for a demote-only fix)
}

impl Tally {
    fn fold(&mut self, unit: &SegmentDelimitations, cand: &SegmentDelimitations) {
        let b = band(unit.reference_tract.len());
        for (ur, cr) in unit.reads.iter().zip(cand.reads.iter()) {
            let (uc, cc) = (class(&ur.observation), class(&cr.observation));
            *self.transitions.entry((uc, cc)).or_default() += 1;
            if uc == "complete" {
                *self.complete_by_band.entry(b).or_default() += 1;
                if cc == "partial" {
                    *self.demote_by_band.entry(b).or_default() += 1;
                }
            }
            if uc == "partial" {
                match cc {
                    "none" => self.partial_lost += 1,
                    "partial" => self.partial_kept += 1,
                    _ => {}
                }
            }
            if cc == "complete" && uc != "complete" {
                self.promotions += 1;
            }
        }
    }

    fn report(&self, name: &str) {
        println!("\n================ {name} vs unit_slip ================");
        println!("class transition unit_slip -> {name}:");
        for ((u, c), n) in &self.transitions {
            let tag = if u == c {
                ""
            } else if *u == "complete" && *c == "partial" {
                "   <- intended demotion"
            } else if *u == "partial" && *c == "none" {
                "   <- PARTIAL DESTROYED"
            } else {
                "   <- change"
            };
            println!("  {u:>8} -> {c:<8} {n:>7}{tag}");
        }
        println!(
            "promotions (candidate=complete, unit!=complete): {}",
            self.promotions
        );
        let ptot = self.partial_kept + self.partial_lost;
        if ptot > 0 {
            println!(
                "genuine partials preserved: {}/{} ({:.1}% kept, {:.1}% destroyed→none)",
                self.partial_kept,
                ptot,
                100.0 * self.partial_kept as f64 / ptot as f64,
                100.0 * self.partial_lost as f64 / ptot as f64,
            );
        }
        println!("complete→partial demotion rate by reference tract length:");
        for b in ["<15", "15-24", "25-39", "40-59", "60+"] {
            let c = *self.complete_by_band.get(b).unwrap_or(&0);
            let d = *self.demote_by_band.get(b).unwrap_or(&0);
            if c > 0 {
                println!(
                    "  {b:>7} completes={c:>7} demoted={d:>6} ({:.1}%)",
                    100.0 * d as f64 / c as f64
                );
            }
        }
    }
}

fn make_gen<A: RepeatDelimiter>(
    fasta: &Path,
    contigs: &ContigList,
    aligner: A,
    bundle: Bp,
) -> Result<
    SsrGenerator<WindowedRefSeq, impl FnMut() -> WindowedRefSeq, A>,
    Box<dyn std::error::Error>,
> {
    Ok(SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        aligner,
        SsrGeneratorConfig::default(),
        bundle,
    )?)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_anchor_firm_validate <reference.fa> <sample.bam|cram> [contig ...]"
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter: Vec<String> = args[3..].to_vec();
    match run(&fasta, &bam, &contig_filter) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    fasta: &Path,
    bam: &Path,
    contig_filter: &[String],
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();
    // One reference for every file this run opens — and so one copy of the bases.
    let reference = OpenReference::new(info);
    let sample = SampleReads::open_only_sample(
        std::slice::from_ref(&bam.to_path_buf()),
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;

    let walk_config = TypedRegionConfig::default();
    let bundle = Bp(walk_config.criteria.bundle_threshold);
    let e = PerQualityEmission::new();
    let mut gen_unit = make_gen(fasta, &contigs, SsrUnitSlipAligner::new(e), bundle)?;
    let mut gen_af = make_gen(fasta, &contigs, SsrAnchorFirmAligner::new(e), bundle)?;
    let mut gen_nr = make_gen(fasta, &contigs, SsrNoiseRobustAligner::new(e), bundle)?;
    let mut gen_ri = make_gen(fasta, &contigs, SsrRobustIndelAligner::new(e), bundle)?;
    let mut gen_ur = make_gen(fasta, &contigs, SsrUnitRobustAligner::new(e), bundle)?;

    let mut af = Tally::default();
    let mut nr = Tally::default();
    let mut ri = Tally::default();
    let mut ur = Tally::default();
    let mut reads = 0u64;

    for (index, entry) in contigs.entries.iter().enumerate() {
        if !contig_filter.is_empty() && !contig_filter.iter().any(|n| n == &entry.name) {
            continue;
        }
        let walk_ref = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
        let mut walk = TypedRegionIterator::over_contig(
            walk_ref,
            ContigId(index as u32),
            walk_config.clone(),
        )?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            gen_unit.begin_segment(region.region);
            gen_af.begin_segment(region.region);
            gen_nr.begin_segment(region.region);
            gen_ri.begin_segment(region.region);
            gen_ur.begin_segment(region.region);
            let u = gen_unit.delimit_segment_reads(segment, &sample)?;
            let a = gen_af.delimit_segment_reads(segment, &sample)?;
            let n = gen_nr.delimit_segment_reads(segment, &sample)?;
            let r = gen_ri.delimit_segment_reads(segment, &sample)?;
            let ur_d = gen_ur.delimit_segment_reads(segment, &sample)?;
            reads += u.reads.len() as u64;
            af.fold(&u, &a);
            nr.fold(&u, &n);
            ri.fold(&u, &r);
            ur.fold(&u, &ur_d);
        }
    }
    if let Some(h) = verify {
        h.join()?;
    }

    println!("reads compared: {reads}");
    af.report("anchor_firm");
    nr.report("noise_robust");
    ri.report("robust_indel");
    ur.report("unit_robust");
    Ok(())
}
