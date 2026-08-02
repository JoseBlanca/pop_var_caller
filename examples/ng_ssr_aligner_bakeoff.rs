//! **STR delimiter bake-off dump** — run every tract delimiter over the same reference + sample and
//! emit a tidy, per-(locus, aligner, observation) table for the period×length "shape of the data"
//! dashboard (`ng_proposal.md` §step-3: *how much stuttering exists for different repeat sizes and
//! lengths, and is there a difference between the aligners*).
//!
//! **⚠ Since 2026-07-28 a row is one `(bases, read_witness, read_group)` CELL, not one allele.**
//! `SequenceObservation` gained the read group as part of its identity, so on a sample declaring
//! several `@RG`s one allele becomes several rows — and **this dump has no read-group column**, so
//! those rows are indistinguishable in the output and the per-row counters count cells rather than
//! alleles. Single-read-group samples are unaffected, which is every fixture here so far. Adding
//! the column is an open question at Checkpoint B: it would change an artifact the marimo
//! dashboards parse, so it is not done silently.
//!
//! ```text
//! ng_ssr_aligner_bakeoff <reference.fa> <sample.bam|cram> [contig ...]
//! ```
//!
//! It walks region typing **once** and feeds each `SsrSegment` to three generators — algorithm 3
//! (`SsrFlatGapAligner`, the flat-gap production-parity port), algorithm 4 (`SsrUnitSlipAligner`,
//! the former default) and algorithm 4u (`SsrUnitRobustAligner`, algorithm 4 hardened with a narrow
//! junction guard and an evidence-based anchor test — **ng's current default**). Because the walk is
//! shared, all three aligners see the *identical* locus set, so their rows join exactly on
//! `(contig, start, end)` — the join a bake-off needs and separate `ng_ssr_loci_dump` runs cannot
//! guarantee. One or more trailing contig names restrict the walk; none = the whole reference.
//!
//! Output: a `#`-prefixed run-level counts header (one line per aligner), a bare TSV column line,
//! then one row per observation. Each covered locus contributes, per aligner, one row per distinct
//! observed sequence (tagged `complete` / `partial:left` / `partial:right`) plus, when non-zero, a
//! synthetic `no_border` row (reads that reached the aligner and anchored nothing) and a `capped`
//! row (reads the depth cap discarded). Zero-coverage loci are omitted — they carry no reads and are
//! aligner-independent. The dashboard derives period = `len(motif)`, tract length = `len(ref_tract)`,
//! and stutter = `obs_len − ref_len`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner;
#[cfg(test)]
use pop_var_caller::ng::locus_generation::WitnessedLocusPositions;
use pop_var_caller::ng::locus_generation::ssr::{
    RepeatDelimiter, SsrGenerator, SsrGeneratorConfig,
};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, LocusLen, ReadWitness, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId, GenomeRegion};

/// The side derivation, shared with the other two STR dumps so the three cannot drift apart
/// again (D4). Each tool keeps its own strings — see `witness_label`.
#[path = "shared/witness_side.rs"]
mod witness_side;
use witness_side::{WitnessSide, witness_side};

/// The three delimiters, and the tag each carries in the `aligner` column.
const FLAT_GAP: &str = "flat_gap"; // algorithm 3 — the production-parity flat-gap port
const UNIT_SLIP: &str = "unit_slip"; // algorithm 4 — the former default
const UNIT_ROBUST: &str = "unit_robust"; // algorithm 4u — algorithm 4 hardened; ng's default

/// One TSV row: one distinct observation (or a synthetic no-border/capped tally) at one locus,
/// under one aligner.
#[derive(Debug, Clone)]
struct Row {
    aligner: &'static str,
    contig: String,
    /// 1-based tract coordinates (inclusive).
    start: u64,
    end: u64,
    motif: Vec<u8>,
    ref_tract: Vec<u8>,
    /// Complete-observation depth (sum of the complete rows' read counts), repeated on every row of
    /// the locus so a partial or no-border row is read against the depth that actually pinned it.
    depth: u32,
    coverage: &'static str,
    observed: Vec<u8>,
    reads: u32,
}

/// Run-level totals for one aligner — the accounting header, mirroring `ng_ssr_loci_dump`.
#[derive(Debug, Clone, Copy, Default)]
struct AlignerCounts {
    zero_coverage: u64,
    reads_fetched: u64,
    reads_capped: u64,
    reads_without_observation: u64,
    obs_complete: u64,
    obs_partial: u64,
}

/// The whole dump: the shared locus count, per-aligner run totals, and the rows.
#[derive(Debug, Clone, Default)]
struct BakeoffReport {
    /// Loci walked (one per `SsrSegment`); shared by every aligner.
    ssr_loci: u64,
    /// One tally per delimiter, in the order the header lines are written.
    counts: Vec<(&'static str, AlignerCounts)>,
    rows: Vec<Row>,
}

impl BakeoffReport {
    /// Render the dump: one `#` header line per aligner, the TSV column line, then the rows.
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        for (tag, c) in &self.counts {
            let _ = writeln!(
                out,
                "# aligner={tag} ssr_loci={} zero_coverage={} reads_fetched={} reads_capped={} \
                 reads_without_observation={} obs_complete={} obs_partial={}",
                self.ssr_loci,
                c.zero_coverage,
                c.reads_fetched,
                c.reads_capped,
                c.reads_without_observation,
                c.obs_complete,
                c.obs_partial,
            );
        }
        out.push_str(
            "aligner\tcontig\tstart\tend\tmotif\tref_tract\tdepth\tcoverage\tobserved\treads\n",
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.aligner,
                row.contig,
                row.start,
                row.end,
                String::from_utf8_lossy(&row.motif),
                String::from_utf8_lossy(&row.ref_tract),
                row.depth,
                row.coverage,
                String::from_utf8_lossy(&row.observed),
                row.reads,
            );
        }
        out
    }
}

/// Fold one aligner's locus into `rows`: an observation row per distinct sequence, plus synthetic
/// `no_border` / `capped` rows when those tallies are non-zero. Zero-coverage loci (no reads reached
/// them at all) emit nothing and only bump the counter.
fn push_locus(
    rows: &mut Vec<Row>,
    aligner: &'static str,
    counts: &mut AlignerCounts,
    locus: &SampleLocusObservations,
    segment: &SsrSegment,
) {
    let has_reads = !locus.observations.is_empty()
        || locus.reads_without_observation > 0
        || locus.reads_discarded_by_cap > 0;
    if !has_reads {
        counts.zero_coverage += 1;
        return;
    }
    let depth: u32 = locus.complete_observations().map(|obs| obs.num_obs).sum();
    let motif = match &locus.kind {
        LocusKind::Ssr(detail) => detail.motif.as_bytes().to_vec(),
        _ => Vec::new(),
    };
    let mut push = |coverage: &'static str, observed: Vec<u8>, reads: u32| {
        rows.push(Row {
            aligner,
            contig: segment.chrom().to_string(),
            start: locus.region.start.get(),
            end: locus.region.end.get(),
            motif: motif.clone(),
            ref_tract: locus.reference_bases.to_vec(),
            depth,
            coverage,
            observed,
            reads,
        });
    };
    for obs in &locus.observations {
        push(
            witness_label(&obs.read_witness, locus.locus_len()),
            obs.bases.to_vec(),
            obs.num_obs,
        );
    }
    if locus.reads_without_observation > 0 {
        push("no_border", Vec::new(), locus.reads_without_observation);
    }
    if locus.reads_discarded_by_cap > 0 {
        push("capped", Vec::new(), locus.reads_discarded_by_cap);
    }
}

/// The tag a witness carries in the `coverage` column.
fn witness_label(witness: &ReadWitness, locus_len: LocusLen) -> &'static str {
    // The derivation is shared (`shared/witness_side.rs`); the spelling is this tool's. **These
    // two strings moved at D4**: `partial_left` / `partial_right` became `partial:left` /
    // `partial:right`, because this function already said `partial:interior` beside them, so a
    // consumer grepping `partial:` got this tool's interiors and none of its sides.
    match witness_side(witness, locus_len) {
        WitnessSide::Complete => "complete",
        WitnessSide::Left => "partial:left",
        WitnessSide::Right => "partial:right",
        WitnessSide::BothBorders => "partial:both",
        WitnessSide::Interior => "partial:interior",
    }
}

/// Build a generator over the given aligner, sharing the reference + read-query factory shape the
/// dump tool uses.
fn make_generator<A: RepeatDelimiter>(
    fasta: &Path,
    contigs: &ContigList,
    aligner: A,
    config: SsrGeneratorConfig,
    bundle_threshold: Bp,
) -> Result<SsrGenerator<WindowedRefSeq, A>, Box<dyn std::error::Error>> {
    let generator = SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        aligner,
        config,
        bundle_threshold,
    )?;
    Ok(generator)
}

/// Run one delimiter over one `SsrSegment` and append its rows — the per-aligner half of the walk,
/// written once so adding a delimiter is one call, not another copy of the loop.
fn dump_segment<A>(
    generator: &mut SsrGenerator<WindowedRefSeq, A>,
    aligner: &'static str,
    counts: &mut AlignerCounts,
    rows: &mut Vec<Row>,
    region: GenomeRegion,
    segment: &SsrSegment,
    sample: &SampleReads,
) -> Result<(), Box<dyn std::error::Error>>
where
    A: RepeatDelimiter,
{
    generator.begin_segment(region);
    while let Some(locus) = generator.next_locus(segment, sample)? {
        push_locus(rows, aligner, counts, &locus, segment);
    }
    Ok(())
}

/// Walk region typing once over `fasta` + `bams` (optionally restricted to `contig_filter`, a set of
/// contig names) and dump every `SsrSegment` through **every** delimiter, building the paired
/// [`BakeoffReport`].
fn run_bakeoff(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: &[String],
) -> Result<BakeoffReport, Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();
    // One reference for every file this run opens — and so one copy of the bases.
    let reference = OpenReference::new(info);

    let sample =
        SampleReads::open_only_sample(bams, &reference, ReadFilterConfig::default(), true)?;

    let walk_config = TypedRegionConfig::default();
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
        config.clone(),
        bundle_threshold,
    )?;
    let mut gen_robust = make_generator(
        fasta,
        &contigs,
        SsrUnitRobustAligner::new(emission),
        config,
        bundle_threshold,
    )?;

    let mut report = BakeoffReport::default();
    let (mut flat, mut unit, mut robust) =
        <(AlignerCounts, AlignerCounts, AlignerCounts)>::default();
    for (index, entry) in contigs.entries.iter().enumerate() {
        if !contig_filter.is_empty() && !contig_filter.iter().any(|name| name == &entry.name) {
            continue;
        }
        let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
        let mut walk = TypedRegionIterator::over_contig(
            walk_reference,
            ContigId(index as u32),
            walk_config.clone(),
        )?;
        for region in walk.by_ref() {
            let region = region?;
            if let RegionKind::SsrSegment(segment) = &region.kind {
                report.ssr_loci += 1;
                let rows = &mut report.rows;
                let region = region.region;
                dump_segment(
                    &mut gen_flat,
                    FLAT_GAP,
                    &mut flat,
                    rows,
                    region,
                    segment,
                    &sample,
                )?;
                dump_segment(
                    &mut gen_unit,
                    UNIT_SLIP,
                    &mut unit,
                    rows,
                    region,
                    segment,
                    &sample,
                )?;
                dump_segment(
                    &mut gen_robust,
                    UNIT_ROBUST,
                    &mut robust,
                    rows,
                    region,
                    segment,
                    &sample,
                )?;
            }
        }
    }

    if let Some(handle) = verify {
        handle.join()?;
    }

    // Run-level per-aligner totals come from the generators' own counters, so the header numbers are
    // the authoritative accounting identity, not a re-tally of the emitted rows.
    for (tag, counts, gen_counts) in [
        (FLAT_GAP, &mut flat, gen_flat.counts()),
        (UNIT_SLIP, &mut unit, gen_unit.counts()),
        (UNIT_ROBUST, &mut robust, gen_robust.counts()),
    ] {
        counts.reads_fetched = gen_counts.reads_fetched;
        counts.reads_capped = gen_counts.reads_discarded_by_cap;
        counts.obs_complete = gen_counts.observations_complete;
        counts.obs_partial = gen_counts.observations_partial;
        // The sum lives on the counts type, behind an exhaustive destructure. This tool did
        // sum the reasons itself and was the one left out when C0 added a fourth, under-
        // reporting by 6,704 reads of ~9,265 on tomato chr01 (Milestone C review, F5).
        counts.reads_without_observation = gen_counts.reads_without_observation();
        report.counts.push((tag, *counts));
    }

    Ok(report)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_aligner_bakeoff <reference.fa> <sample.bam|cram> [contig ...]\n\
             dumps, per microsatellite tract, what ALL THREE STR delimiters (algorithm 3 flat-gap, \
             algorithm 4 unit-slip and algorithm 4u unit-robust, the default) observed — the \
             period×length bake-off input. Trailing contig names restrict the walk; none = the \
             whole reference."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter: Vec<String> = args[3..].to_vec();

    match run_bakeoff(&fasta, &[bam], &contig_filter) {
        Ok(report) => {
            print!("{}", report.render());
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The four strings this tool prints in its `coverage` column**, pinned — see the identical
    /// test in `ng_ssr_cohort_stutter` for why (D4's rename, and the Milestone D reliability
    /// review finding that a binary with no tests notices neither a rename nor a mutation
    /// labelling every partial `complete`).
    ///
    /// This tool's column is the sharper case: `ng_ssr_aligner_bakeoff_dashboard.py` maps it into
    /// outcome classes, and an unmapped label becomes a `NaN` that every downstream count drops.
    /// The notebook now asserts on an unknown label; this asserts the labels it is given.
    #[test]
    fn the_coverage_column_spells_the_four_cases_the_dashboard_maps() {
        let len = LocusLen::from_positions(10);
        let partial = |runs: &[(u16, u16)]| {
            witness_label(
                &ReadWitness::Partial {
                    positions: WitnessedLocusPositions::from_half_open_runs(runs.iter().copied())
                        .expect("a non-empty set of runs"),
                },
                len,
            )
        };
        assert_eq!(witness_label(&ReadWitness::Complete, len), "complete");
        assert_eq!(partial(&[(0, 4)]), "partial:left");
        assert_eq!(partial(&[(6, 10)]), "partial:right");
        assert_eq!(partial(&[(3, 7)]), "partial:interior");
        assert_eq!(
            partial(&[(0, 10)]),
            "partial:both",
            "a repeat read that ran out covers the tract end to end without measuring it",
        );
        assert_eq!(
            partial(&[(0, 3), (7, 10)]),
            "partial:both",
            "so does a read blind in the middle — neither is a measurement",
        );
    }
}
