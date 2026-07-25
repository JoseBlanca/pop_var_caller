//! **STR delimiter bake-off dump** — run BOTH tract delimiters over the same reference + sample and
//! emit a tidy, per-(locus, aligner, observation) table for the period×length "shape of the data"
//! dashboard (`ng_proposal.md` §step-3: *how much stuttering exists for different repeat sizes and
//! lengths, and is there a difference between the aligners*).
//!
//! ```text
//! ng_ssr_aligner_bakeoff <reference.fa> <sample.bam|cram> [contig ...]
//! ```
//!
//! It walks region typing **once** and feeds each `SsrSegment` to two generators — algorithm 3
//! (`SsrFlatGapAligner`, the flat-gap production-parity port) and algorithm 4 (`SsrUnitSlipAligner`,
//! the recommended unit-slip delimiter). Because the walk is shared, the two aligners see the
//! *identical* locus set, so their rows join exactly on `(contig, start, end)` — the join a bake-off
//! needs and two separate `ng_ssr_loci_dump` runs cannot guarantee. One or more trailing contig
//! names restrict the walk; none = the whole reference.
//!
//! Output: a `#`-prefixed run-level counts header (one line per aligner), a bare TSV column line,
//! then one row per observation. Each covered locus contributes, per aligner, one row per distinct
//! observed sequence (tagged `complete` / `partial_left` / `partial_right`) plus, when non-zero, a
//! synthetic `no_border` row (reads that reached the aligner and anchored nothing) and a `capped`
//! row (reads the depth cap discarded). Zero-coverage loci are omitted — they carry no reads and are
//! aligner-independent. The dashboard derives period = `len(motif)`, tract length = `len(ref_tract)`,
//! and stutter = `obs_len − ref_len`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::locus_generation::ssr::{RepeatDelimiter, SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, ReadCoverage, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId};

/// The two delimiters, and the tag each carries in the `aligner` column.
const FLAT_GAP: &str = "flat_gap"; // algorithm 3 — the production-parity flat-gap port
const UNIT_SLIP: &str = "unit_slip"; // algorithm 4 — the recommended unit-slip delimiter

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
    /// Loci walked (one per `SsrSegment`); shared by both aligners.
    ssr_loci: u64,
    flat: AlignerCounts,
    unit: AlignerCounts,
    rows: Vec<Row>,
}

impl BakeoffReport {
    /// Fold one aligner's locus into the report: emit an observation row per distinct sequence, and
    /// synthetic `no_border` / `capped` rows when those tallies are non-zero. Zero-coverage loci
    /// (no reads reached them at all) emit nothing.
    fn push_locus(
        &mut self,
        aligner: &'static str,
        counts: &mut AlignerCounts,
        locus: &SampleLocusObservations,
        segment: &SsrSegment,
    ) {
        let has_reads = !locus.observed_sequences.is_empty()
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
            self.rows.push(Row {
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
        for obs in &locus.observed_sequences {
            push(
                coverage_label(obs.read_coverage),
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

    /// Render the dump: two `#` header lines (one per aligner), the TSV column line, then the rows.
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        let header = |out: &mut String, tag: &str, c: &AlignerCounts| {
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
        };
        header(&mut out, FLAT_GAP, &self.flat);
        header(&mut out, UNIT_SLIP, &self.unit);
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

/// The tag a read-coverage carries in the `coverage` column.
fn coverage_label(coverage: ReadCoverage) -> &'static str {
    match coverage {
        ReadCoverage::Complete => "complete",
        ReadCoverage::PartialLeft(_) => "partial_left",
        ReadCoverage::PartialRight(_) => "partial_right",
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
) -> Result<
    SsrGenerator<WindowedRefSeq, impl FnMut() -> WindowedRefSeq, A>,
    Box<dyn std::error::Error>,
> {
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

/// Walk region typing once over `fasta` + `bams` (optionally restricted to `contig_filter`, a set of
/// contig names) and dump every `SsrSegment` through **both** delimiters, building the paired
/// [`BakeoffReport`].
fn run_bakeoff(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: &[String],
) -> Result<BakeoffReport, Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();

    let sample = SampleReads::open(bams, &info, ReadFilterConfig::default(), true)?;

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
        config,
        bundle_threshold,
    )?;

    let mut report = BakeoffReport::default();
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
                gen_flat.begin_segment(region.region);
                while let Some(locus) = gen_flat.next_locus(segment, &sample)? {
                    let mut counts = report.flat;
                    report.push_locus(FLAT_GAP, &mut counts, &locus, segment);
                    report.flat = counts;
                }
                gen_unit.begin_segment(region.region);
                while let Some(locus) = gen_unit.next_locus(segment, &sample)? {
                    let mut counts = report.unit;
                    report.push_locus(UNIT_SLIP, &mut counts, &locus, segment);
                    report.unit = counts;
                }
            }
        }
    }

    if let Some(handle) = verify {
        handle.join()?;
    }

    // Run-level per-aligner totals come from the generators' own counters, so the header numbers are
    // the authoritative accounting identity, not a re-tally of the emitted rows.
    for (counts, gen_counts) in [
        (&mut report.flat, gen_flat.counts()),
        (&mut report.unit, gen_unit.counts()),
    ] {
        counts.reads_fetched = gen_counts.reads_fetched;
        counts.reads_capped = gen_counts.reads_discarded_by_cap;
        counts.obs_complete = gen_counts.observations_complete;
        counts.obs_partial = gen_counts.observations_partial;
        counts.reads_without_observation = gen_counts.no_border_anchored
            + gen_counts.low_quality
            + gen_counts.window_truncated;
    }

    Ok(report)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_aligner_bakeoff <reference.fa> <sample.bam|cram> [contig ...]\n\
             dumps, per microsatellite tract, what BOTH STR delimiters (algorithm 3 flat-gap and \
             algorithm 4 unit-slip) observed — the period×length bake-off input. Trailing contig \
             names restrict the walk; none = the whole reference."
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
