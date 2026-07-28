//! **Per-sample STR stutter dump** — walk region typing once and delimit every sample's reads at
//! the *same* microsatellite tracts, emitting one tidy row per (sample, locus, observation).
//!
//! **⚠ Since 2026-07-28 a row is one `(bases, read_coverage, read_group)` CELL, not one allele.**
//! `ObservedSequence` gained the read group as part of its identity, so on a sample declaring
//! several `@RG`s one allele becomes several rows — and **this dump has no read-group column**, so
//! those rows are indistinguishable in the output and the per-row counters below count cells
//! rather than alleles. Single-read-group samples are unaffected, which is every fixture here so
//! far. Adding the column is an open question at Checkpoint B: it would change an artifact the
//! marimo dashboards parse, so it is not done silently.
//!
//! ```text
//! ng_ssr_cohort_stutter [--contigs a,b] [--regions r.bed] <reference.fa> <sample.cram> [sample ...]
//! ```
//!
//! `--contigs` and `--regions` restrict the walk. They are flags, not environment variables,
//! because the dev container wrapper forwards only a fixed set of variables — an env knob would be
//! silently ignored inside it and walk the whole reference.
//!
//! **`--regions` is what makes a cohort run finish.** Typing a tract is cheap but not free, and a
//! walk types every tract in the reference whether or not a read reaches it. When the inputs are
//! slices cut to a target BED — as the ssr_tomato1 cohort CRAMs are — the overwhelming majority of
//! typed tracts have no reads in any sample, and walking them is pure waste: over 90 Mb of tomato
//! chromosome 1 that waste is ~8 minutes per sample, against ~1.4 Mb of actual targets. Pass the
//! same BED the inputs were sliced to.
//!
//! The bake-off dump (`ng_ssr_aligner_bakeoff`) varies the *delimiter* over one sample; this one
//! fixes the delimiter at ng's default (`SsrUnitRobustAligner`, algorithm 4u) and varies the
//! *sample*. That is the shape a cohort question needs: how much a locus stutters is a property of
//! the sample's library chemistry (PCR amplification above all), so the comparison worth making is
//! sample-to-sample **at identical loci**. Sharing one region-typing walk guarantees the identical
//! locus set, so rows join exactly on `(contig, start, end)`.
//!
//! Rows stream to stdout as they are produced rather than accumulating — a 50-sample cohort emits
//! millions of rows, and buffering them all would cost gigabytes for no benefit.
//!
//! Output: a `#`-prefixed run-level line per sample, a bare TSV column line, then the rows. Each
//! covered locus contributes, per sample, one row per distinct observed sequence (tagged
//! `complete` / `partial_left` / `partial_right`) plus, when non-zero, a synthetic `no_border` row
//! (reads that reached the aligner and anchored nothing) and a `capped` row (reads the depth cap
//! discarded). A locus no read of that sample reaches emits nothing for that sample.
//!
//! The dashboard derives period = `len(motif)`, tract length = `len(ref_tract)`, and two stutter
//! measures: `obs_len − ref_len` (off-reference, comparable across samples but confounded by
//! genuine non-reference alleles) and `obs_len − modal obs_len of that sample at that locus`
//! (off-mode — the stutter-specific one, since a sample's own modal allele is its best available
//! stand-in for its true genotype).

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, LocusLen, ReadCoverage, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegionConfig, TypedRegionIterator,
};
use pop_var_caller::ng::types::{Bp, ContigId};
use pop_var_caller::regions::ContigBounds;

/// Run-level totals for one sample — the accounting identity `reads_fetched = complete + partial +
/// no_border + capped`, tallied from the rows because the generator's own counters are shared by
/// every sample it serves.
#[derive(Debug, Clone, Default)]
struct SampleCounts {
    name: String,
    covered_loci: u64,
    obs_complete: u64,
    obs_partial: u64,
    reads_complete: u64,
    reads_partial: u64,
    reads_no_border: u64,
    reads_capped: u64,
}

/// Stream one sample's locus into `out`: an observation row per distinct sequence, plus synthetic
/// `no_border` / `capped` rows when those tallies are non-zero. A locus this sample has no reads at
/// emits nothing — it is that sample's zero-coverage, and the row set stays sparse.
fn write_locus<W: Write>(
    out: &mut W,
    counts: &mut SampleCounts,
    locus: &SampleLocusObservations,
    segment: &SsrSegment,
) -> std::io::Result<()> {
    let has_reads = !locus.observed_sequences.is_empty()
        || locus.reads_without_observation > 0
        || locus.reads_discarded_by_cap > 0;
    if !has_reads {
        return Ok(());
    }
    counts.covered_loci += 1;
    let depth: u32 = locus.complete_observations().map(|obs| obs.num_obs).sum();
    let motif = match &locus.kind {
        LocusKind::Ssr(detail) => detail.motif.as_bytes().to_vec(),
        _ => Vec::new(),
    };
    let chrom = segment.chrom();
    let start = locus.region.start.get();
    let end = locus.region.end.get();
    let ref_tract = String::from_utf8_lossy(&locus.reference_bases);
    let motif_str = String::from_utf8_lossy(&motif);

    let row = |out: &mut W, coverage: &str, observed: &str, reads: u32| {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            counts.name, chrom, start, end, motif_str, ref_tract, depth, coverage, observed, reads,
        )
    };

    for obs in &locus.observed_sequences {
        let label = coverage_label(obs.read_coverage, locus.locus_len());
        match obs.read_coverage {
            ReadCoverage::Complete => {
                counts.obs_complete += 1;
                counts.reads_complete += u64::from(obs.num_obs);
            }
            _ => {
                counts.obs_partial += 1;
                counts.reads_partial += u64::from(obs.num_obs);
            }
        }
        row(
            out,
            label,
            &String::from_utf8_lossy(&obs.bases),
            obs.num_obs,
        )?;
    }
    if locus.reads_without_observation > 0 {
        counts.reads_no_border += u64::from(locus.reads_without_observation);
        row(out, "no_border", "", locus.reads_without_observation)?;
    }
    if locus.reads_discarded_by_cap > 0 {
        counts.reads_capped += u64::from(locus.reads_discarded_by_cap);
        row(out, "capped", "", locus.reads_discarded_by_cap)?;
    }
    Ok(())
}

/// The tag a read-coverage carries in the `coverage` column.
fn coverage_label(coverage: ReadCoverage, locus_len: LocusLen) -> &'static str {
    // Since the reshape the side is a **derivation**, not a variant: a run flush with the left
    // border is a prefix constraint, one flush with the right border a suffix. A run flush with
    // neither is interior — the STR path cannot mint one (it anchors a border or yields nothing),
    // so it never appears here, but naming it keeps the label honest for the generic path.
    // Destructured rather than guarded on `_`, so a future `ReadCoverage` variant is a
    // compile error here. The guard form is what this migration used and it is exactly what
    // let the compiler stop forcing these sites to be revisited.
    match coverage {
        ReadCoverage::Complete => "complete",
        run @ ReadCoverage::Observed { .. } => {
            match (run.is_flush_left(), run.is_flush_right(locus_len)) {
                (true, _) => "partial_left",
                (false, true) => "partial_right",
                (false, false) => "partial:interior",
            }
        }
    }
}

/// Walk region typing once over `fasta` (optionally restricted to the contigs named in
/// `PVC_CONTIGS`) and stream every `SsrSegment` through **every** sample, with ng's default
/// delimiter. Returns the per-sample tallies for the header, which is written last to stderr-free
/// stdout by the caller — hence the two-pass shape: rows first, counts reported by the caller.
fn run_cohort(
    fasta: &Path,
    crams: &[PathBuf],
    contig_filter: &[String],
    regions_bed: Option<&Path>,
) -> Result<Vec<SampleCounts>, Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs: ContigList = info.contig_list();

    // Group the inputs by the sample their read groups name, rather than assuming one file is one
    // sample. Several files of a single sample is the normal case for a library sequenced across
    // lanes — the madrid_herb1 herbarium specimen is eight. Getting this wrong would not merely
    // mislabel: it would split one library into k pseudo-samples whose stutter is correlated by
    // construction, inflating any between-sample statistic.
    //
    // The read-group pre-pass does the grouping now. It reads every header once, mints an
    // identifier per `@RG`, and hands back the read groups grouped by sample — which is what this
    // loop used to reconstruct by opening each file on its own just to read its header.
    let read_groups = build_read_groups(crams)?;
    for entry in read_groups.read_groups_per_sample() {
        let files: usize = {
            let mut paths: Vec<&Path> = entry
                .read_groups
                .iter()
                .map(|id| &*read_groups.get(*id).file)
                .collect();
            paths.sort();
            paths.dedup();
            paths.len()
        };
        if files > 1 {
            eprintln!("  {}: {files} files merged into one sample", entry.sample);
        }
    }
    let samples: Vec<SampleReads> = read_groups
        .read_groups_per_sample()
        .iter()
        .map(|entry| {
            SampleReads::open(
                entry,
                &read_groups,
                &info,
                ReadFilterConfig::default(),
                true,
            )
        })
        .collect::<Result<_, _>>()?;
    let mut counts: Vec<SampleCounts> = samples
        .iter()
        .map(|s| SampleCounts {
            name: s.sample_name().to_string(),
            ..Default::default()
        })
        .collect();

    let walk_config = TypedRegionConfig::default();
    let bundle_threshold = Bp(walk_config.criteria.bundle_threshold);
    // **One reference reader for the whole walk, shared** — the margin fetch and the per-query
    // read filter both hold the same `Arc`. Building a fresh `WindowedRefSeq` per query (which is
    // what `FnMut() -> R` invited) meant re-reading the whole `.fai` and re-`open`ing the FASTA
    // before serving one ~150-base window: 14% of a cohort run, ~564k `open(2)`s per chromosome.
    // The shared reader establishes its window once and slides.
    // `Arc` rather than `Rc` even though this walk is single-threaded: `RawRefSeq`
    // is implemented for `Arc<T>` and nothing else (`ref_seq.rs`), so the generator
    // cannot take an `Rc`. Clippy's usual remedy therefore does not apply here.
    #[expect(
        clippy::arc_with_non_send_sync,
        reason = "RawRefSeq is implemented for Arc only; this walk is single-threaded"
    )]
    let reference = Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));
    let mut generator = SsrGenerator::with_default_aligner(
        Arc::clone(&reference),
        {
            let reference = Arc::clone(&reference);
            move || Arc::clone(&reference)
        },
        SsrGeneratorConfig::default(),
        bundle_threshold,
    )?;

    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    writeln!(
        out,
        "sample\tcontig\tstart\tend\tmotif\tref_tract\tdepth\tcoverage\tobserved\treads"
    )?;

    // The BED path walks only the targeted spans; without one, every contig end to end. Both feed
    // the same per-segment body, so restricting the walk cannot change what a covered locus is —
    // only how much untargeted sequence is typed on the way to it.
    let bed_spans = regions_bed
        .map(|bed| {
            let bounds: Vec<ContigBounds<'_>> = contigs
                .entries
                .iter()
                .map(|e| ContigBounds {
                    name: &e.name,
                    length: e.length as u32,
                })
                .collect();
            GenomeRegions::from_bed_path(bed, &bounds)
        })
        .transpose()?;

    // `over_regions` takes the whole BED; a contig restriction is applied to the segments it
    // yields rather than to the span list, because narrowing a `GenomeRegions` is not part of its
    // surface and a BED walk is cheap enough that typing the unwanted contigs' spans is not worth
    // a new constructor.
    let wanted_contig = |contig: ContigId| {
        contig_filter.is_empty()
            || contigs
                .entries
                .get(contig.0 as usize)
                .is_some_and(|e| contig_filter.iter().any(|n| n == &e.name))
    };

    let mut walks: Vec<(String, TypedRegionIterator<WindowedRefSeq>)> = Vec::new();
    match bed_spans {
        Some(spans) => {
            eprintln!("  {} BED spans", spans.iter().count());
            let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
            walks.push((
                "BED".to_string(),
                TypedRegionIterator::over_regions(walk_reference, spans, walk_config.clone())?,
            ));
        }
        None => {
            for (index, entry) in contigs.entries.iter().enumerate() {
                if !wanted_contig(ContigId(index as u32)) {
                    continue;
                }
                let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
                walks.push((
                    entry.name.clone(),
                    TypedRegionIterator::over_contig(
                        walk_reference,
                        ContigId(index as u32),
                        walk_config.clone(),
                    )?,
                ));
            }
        }
    }

    for (label, mut walk) in walks {
        // A cohort walk runs for a long time; say what is in hand so it can be watched.
        eprintln!("  walking {label}");
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            if !wanted_contig(region.region.contig) {
                continue;
            }
            for (sample, counts) in samples.iter().zip(counts.iter_mut()) {
                generator.begin_segment(region.region);
                while let Some(locus) = generator.next_locus(segment, sample)? {
                    write_locus(&mut out, counts, &locus, segment)?;
                }
            }
        }
        // Flush per walk so a long run shows progress in the output file, not just at the end.
        out.flush()?;
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(counts)
}

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut contig_filter: Vec<String> = Vec::new();
    let mut regions_bed: Option<PathBuf> = None;
    let mut rest = std::env::args().skip(1);
    while let Some(arg) = rest.next() {
        match arg.as_str() {
            "--regions" => match rest.next() {
                Some(path) => regions_bed = Some(PathBuf::from(path)),
                None => {
                    eprintln!("error: --regions needs a BED path");
                    return ExitCode::from(2);
                }
            },
            "--contigs" => match rest.next() {
                Some(list) => {
                    contig_filter = list.split(',').map(|s| s.trim().to_string()).collect()
                }
                None => {
                    eprintln!("error: --contigs needs a comma-separated list");
                    return ExitCode::from(2);
                }
            },
            _ => positional.push(arg),
        }
    }
    if positional.len() < 2 {
        eprintln!(
            "usage: ng_ssr_cohort_stutter [--contigs a,b] [--regions r.bed] <reference.fa> \
             <sample.bam|cram> [sample ...]\n\
             dumps, per microsatellite tract and per sample, what ng's default STR delimiter \
             observed — the per-sample stutter input. Pass --regions the BED the inputs were \
             sliced to; without it the walk types the whole reference."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&positional[0]);
    let crams: Vec<PathBuf> = positional[1..].iter().map(PathBuf::from).collect();

    eprintln!(
        "walking {} sample(s){}",
        crams.len(),
        if contig_filter.is_empty() {
            String::new()
        } else {
            format!(" over {}", contig_filter.join(", "))
        }
    );

    match run_cohort(&fasta, &crams, &contig_filter, regions_bed.as_deref()) {
        Ok(counts) => {
            // The accounting goes to stderr so the TSV on stdout stays a clean single table for a
            // cohort this size — a per-sample header block would be 50 comment lines.
            eprintln!(
                "sample\tcovered_loci\tobs_complete\tobs_partial\treads_complete\treads_partial\t\
                 reads_no_border\treads_capped"
            );
            for c in &counts {
                eprintln!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    c.name,
                    c.covered_loci,
                    c.obs_complete,
                    c.obs_partial,
                    c.reads_complete,
                    c.reads_partial,
                    c.reads_no_border,
                    c.reads_capped,
                );
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
