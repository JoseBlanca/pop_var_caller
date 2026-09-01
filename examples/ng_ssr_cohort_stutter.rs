//! **Per-sample STR stutter dump** — walk region typing once and delimit every sample's reads at
//! the *same* microsatellite tracts, emitting one tidy row per (sample, locus, observation).
//!
//! A row is one `(bases, read_witness, read_group)` **cell**, not one allele:
//! `SequenceObservation` carries the read group as part of its identity, so on a sample declaring
//! several `@RG`s one allele becomes several rows. The `read_group` column below is what keeps
//! those rows distinguishable — without it they would be indistinguishable in the output and the
//! per-row counters would silently count cells while reading as alleles.
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
//! Output: the run's **read-group table** as `#rg` lines, a bare TSV column line, then the rows.
//! Each covered locus contributes, per sample, one row per distinct (sequence, witness, read
//! group) — the witness tagged `complete` / `partial:left` / `partial:right` — plus, when non-zero,
//! a synthetic `no_border` row (reads that reached the aligner and anchored nothing) and a `capped`
//! row (reads the depth cap discarded). A locus no read of that sample reaches emits nothing for
//! that sample.
//!
//! **The read group is on every observation row**, because that is the grain chemistry actually
//! varies at. An allele seen from two groups is two rows, so a per-group model gets the allele ×
//! group cross *with its quality moments* rather than a count that has already been merged. Rows
//! carry only the numeric id; the `#rg` table maps it to sample, library, experiment, platform and
//! file, so a consumer picks its own grain — read group, library, experiment or sample — instead of
//! this step guessing one. That matters here: one BioSample can hold several libraries (and in this
//! archive one holds sixteen), so folding to the sample would destroy exactly the contrast a
//! chemistry question is asking about.
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
#[cfg(test)]
use pop_var_caller::ng::locus_generation::WitnessedLocusPositions;
use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, LocusLen, ReadWitness, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::{NameOrigin, build_read_groups};
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::GenomeRegion;
use pop_var_caller::ng::types::{Bp, ContigId};
use pop_var_caller::regions::ContigBounds;

#[path = "shared/catalog_regions.rs"]
mod catalog_regions;
/// The side derivation, shared with the other two STR dumps so the three cannot drift apart
/// again (D4). Each tool keeps its own strings — see `witness_label`.
#[path = "shared/witness_side.rs"]
mod witness_side;
use witness_side::{WitnessSide, witness_side};

/// Run-level totals for one sample — the accounting identity `reads_fetched = complete + partial +
/// no_border + capped`, tallied from the rows.
///
/// Tallied from the rows rather than read off the generator, which is a habit from when one
/// generator served every sample and its counters could not be attributed. Each sample now has its
/// own generator, so its counters *would* be per sample — but the rows are what the dashboards
/// parse, so counting them keeps the header and the body answering from one source.
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
    let has_reads = !locus.observations.is_empty()
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

    // `read_group` is the numeric id from the run's table, not a name: the names live once in the
    // `#rg` header rather than repeated across millions of rows. Empty for the synthetic
    // `no_border` / `capped` tallies, which are per-locus counters and belong to no one group.
    let row = |out: &mut W, rg: &str, coverage: &str, observed: &str, reads: u32| {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            counts.name,
            rg,
            chrom,
            start,
            end,
            motif_str,
            ref_tract,
            depth,
            coverage,
            observed,
            reads,
        )
    };

    for obs in &locus.observations {
        let label = witness_label(&obs.read_witness, locus.locus_len());
        match obs.read_witness {
            ReadWitness::Complete => {
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
            &obs.read_group.get().to_string(),
            label,
            &String::from_utf8_lossy(&obs.bases),
            obs.num_obs,
        )?;
    }
    if locus.reads_without_observation > 0 {
        counts.reads_no_border += u64::from(locus.reads_without_observation);
        row(out, "", "no_border", "", locus.reads_without_observation)?;
    }
    if locus.reads_discarded_by_cap > 0 {
        counts.reads_capped += u64::from(locus.reads_discarded_by_cap);
        row(out, "", "capped", "", locus.reads_discarded_by_cap)?;
    }
    Ok(())
}

/// Whether a grouping name is the file's own or one this run invented — the reader has to be able
/// to tell, because a synthesized library is a guess and a declared one is evidence.
fn origin_label(origin: NameOrigin) -> &'static str {
    match origin {
        NameOrigin::Declared => "declared",
        NameOrigin::Synthesized => "synthesized",
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
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs: ContigList = info.contig_list();
    // **The typed regions come from the catalog beside the reference**, checked against what
    // the pass just reported. No catalog, no run — the error names the command that writes one.
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    // **One reference for the whole cohort, and so one copy of its bases.** A
    // `fasta::Repository` memoises whole contigs and never evicts, so the
    // per-file repository this replaces cost ~752 MiB of resident tomato
    // genome per open CRAM — 51 samples asked for 38 GiB against a 16 GB cap
    // and were OOM-killed at ~80 s. Handing every `SampleReads::open` the same
    // `OpenReference` makes that one genome, once.
    let reference = OpenReference::new(info);

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
                &reference,
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
    let reference = Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));
    // **One generator per sample, and that is a requirement rather than a tidiness.** A
    // generator opens a reader for one sample's files and keeps it positioned for a whole
    // chromosome, so a generator shared between samples would answer every sample out of the
    // first one's files — this tool's whole question is whether samples differ, and it would
    // have reported one sample N times, with no error and rows of exactly the right shape.
    // `SsrGenerator` now refuses a sample it was not opened for, so the mistake is a loud one;
    // this is the shape that does not make it.
    //
    // The *reference* is still shared across all of them: it is read-only sliding-window access
    // to the same FASTA, and giving each sample its own would re-read the whole `.fai` and
    // re-`open(2)` the file per sample.
    let mut generators = samples
        .iter()
        .map(|_| {
            SsrGenerator::with_default_aligner(
                Arc::clone(&reference),
                {
                    let reference = Arc::clone(&reference);
                    move || Arc::clone(&reference)
                },
                SsrGeneratorConfig::default(),
                bundle_threshold,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;

    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());

    // The read-group table, once, as `#`-prefixed lines above the data. The rows carry only the
    // numeric id, so this is what makes them resolvable — and it is what lets a consumer fold to
    // whatever grain it wants. `library` and `experiment` each carry the origin of their name,
    // because a grouping this module synthesized and one the file declared are not equally
    // trustworthy and a chemistry report has to be able to say which it used.
    for (id, group) in read_groups.iter() {
        writeln!(
            out,
            "#rg\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            id.get(),
            group.id,
            group.sample,
            group.library.value,
            origin_label(group.library.origin),
            group.experiment.value,
            origin_label(group.experiment.origin),
            group.platform.as_deref().unwrap_or(""),
            group.file.display(),
        )?;
    }
    // `read_group` is minted per run and means nothing across runs; `(file, rg_id)` is the stable
    // identity, because the SAM specification makes `@RG ID` unique within its file. Anything that
    // merges the output of two runs — a batched cohort, say — has to renumber on that pair.
    writeln!(
        out,
        "#rg_columns\tread_group\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\t\
         experiment_origin\tplatform\tfile"
    )?;
    writeln!(
        out,
        "sample\tread_group\tcontig\tstart\tend\tmotif\tref_tract\tdepth\tcoverage\tobserved\treads"
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

    // What to ask the catalog for, labelled so a long run can be watched. The spans are
    // collected first and the reader built per batch, because a reader borrows the catalog.
    let mut batches: Vec<(String, Vec<GenomeRegion>)> = Vec::new();
    match bed_spans {
        Some(spans) => {
            let wanted: Vec<GenomeRegion> = spans.iter().collect();
            eprintln!("  {} BED spans", wanted.len());
            batches.push(("BED".to_string(), wanted));
        }
        None => {
            for (index, entry) in contigs.entries.iter().enumerate() {
                if !wanted_contig(ContigId(index as u32)) {
                    continue;
                }
                batches.push((
                    entry.name.clone(),
                    vec![catalog_regions::whole_contig(
                        ContigId(index as u32),
                        entry.length,
                    )],
                ));
            }
        }
    }

    let criteria = StrRepeatCriteria::from(&walk_config);
    for (label, spans) in batches {
        // A cohort run takes a long time; say what is in hand so it can be watched.
        eprintln!("  walking {label}");
        let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(&spans))?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            if !wanted_contig(region.region.contig) {
                continue;
            }
            for ((sample, counts), generator) in samples
                .iter()
                .zip(counts.iter_mut())
                .zip(generators.iter_mut())
            {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// **The four strings this tool prints in its `coverage` column**, pinned.
    ///
    /// D4 moved two of them — `partial_left` / `partial_right` became `partial:left` /
    /// `partial:right`, so the three STR dumps stop disagreeing about a separator — and the
    /// Milestone D reliability review then found that nothing in the tree noticed either the
    /// rename or a mutation labelling *every* partial `complete`, because this binary had no
    /// tests at all. It has one now, and it is not decoration: the dashboards key on this column
    /// (`ng_ssr_cohort_stutter_dashboard.py:172` selects `coverage == "complete"` to decide which
    /// reads carry an exact length), so a partial mislabelled `complete` feeds a censored lower
    /// bound into a stutter distribution silently.
    ///
    /// The derivation itself lives in `shared/witness_side.rs` and is exercised by
    /// `ng_ssr_loci_dump`'s fixtures; what is this tool's own, and only this tool's, is the
    /// spelling.
    #[test]
    fn the_coverage_column_spells_the_four_cases_the_dashboards_read() {
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
