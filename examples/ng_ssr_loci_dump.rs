//! **The STR locus generator's dump tool** — `locus_generation_ssr.md` Milestone E, the
//! acceptance anchor. Runs the real pipeline over a reference + a sample's reads and prints, per
//! microsatellite tract, what that sample's reads showed: the observed tract sequences and their
//! support. This generator emits no calls, so "done" must not decay into "compiles" — the output
//! is inspectable and, on a committed fixture, asserted (spec §9).
//!
//! ```text
//! ng_ssr_loci_dump <reference.fa> <sample.bam|cram> [contig]
//! ```
//!
//! The reference's sibling `<reference.fa>.fai` is used (created if absent); reads go through ng's
//! real ingestion (`SampleReads`, the filtered reads a caller sees). An optional third argument
//! restricts the walk to one contig by name. The pipeline is: region typing over the reference →
//! one `SsrSegment` per tract → the `SsrGenerator` turns each into one locus → its observations
//! become rows.
//!
//! Output: a `#`-prefixed `key=value` counts header, a bare TSV column line, then one tab-separated
//! row per observed sequence. `depth` sums the reads behind **complete** observations only (the ones
//! that pinned the tract); partial reads appear on their own rows, tagged `partial:left` /
//! `partial:right`, because conflating a lower bound with an exact length is the mistake spec §3
//! exists to prevent.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::alignment::PerQualityEmission;
use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner;
use pop_var_caller::ng::locus_generation::ssr::{
    RepeatDelimiter, SsrGenerator, SsrGeneratorConfig,
};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, LocusLen, ReadWitness, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegionConfig, TypedRegionIterator};
use pop_var_caller::ng::types::{Bp, ContigId};

/// One TSV row: an observed tract sequence at a locus, with its support.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ObservationRow {
    contig: String,
    /// 1-based tract coordinates (inclusive).
    start: u64,
    end: u64,
    motif: Vec<u8>,
    ref_tract: Vec<u8>,
    /// The locus's complete-observation depth (sum of the complete rows' read counts) — shown on
    /// every row of the locus, so a partial row is read against the depth that actually pinned it.
    depth: u32,
    read_witness: &'static str,
    observed: Vec<u8>,
    reads: u32,
}

/// The whole dump: the run-level counts plus the rows. Separated from rendering so the tests can
/// assert on the numbers directly, not by parsing text.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct DumpReport {
    /// Loci emitted — one per `SsrSegment`, including uncovered ones.
    ssr_loci: u64,
    /// Loci no read reached at all (present-but-empty).
    zero_coverage: u64,
    /// Reads the depth cap discarded, across all loci.
    reads_capped: u64,
    /// Reads that reached the aligner and yielded nothing (by any reason).
    reads_without_observation: u64,
    /// Complete / partial observation totals (reads that became each kind).
    obs_complete: u64,
    obs_partial: u64,
    /// Reads fetched over all loci — **per fetch event, not per distinct read** (a read
    /// overlapping two loci' query spans is fetched, and counted, at each). This is the left side
    /// of the accounting identity (`= obs_complete + obs_partial + reads_capped +
    /// reads_without_observation`), which holds because every fetched read folds into exactly one
    /// of those run-level counters.
    reads_fetched: u64,
    /// The typed-region walk's own `ssr_loci` — `ssr_loci` must equal it (one locus per segment).
    walk_ssr_loci: u64,
    rows: Vec<ObservationRow>,
}

impl DumpReport {
    /// Fold one locus into the report: count it, note zero-coverage, and emit a row per observation.
    fn push_locus(&mut self, locus: &SampleLocusObservations, segment: &SsrSegment) {
        self.ssr_loci += 1;
        // Zero coverage is "no read reached the tract" — distinct from "reads reached it and said
        // nothing" (that is `reads_without_observation`) and from "the cap dropped them all".
        if locus.observations.is_empty()
            && locus.reads_without_observation == 0
            && locus.reads_discarded_by_cap == 0
        {
            self.zero_coverage += 1;
        }
        let depth: u32 = locus.complete_observations().map(|obs| obs.num_obs).sum();
        // The locus length in positions — what turns a run's offset into "flush left" or
        // "flush right", now that the side is derived rather than tagged.
        let locus_len = locus.locus_len();
        let motif = match &locus.kind {
            LocusKind::Ssr(detail) => detail.motif.as_bytes().to_vec(),
            _ => Vec::new(),
        };
        for obs in &locus.observations {
            self.rows.push(ObservationRow {
                contig: segment.chrom().to_string(),
                start: locus.region.start.get(),
                end: locus.region.end.get(),
                motif: motif.clone(),
                ref_tract: locus.reference_bases.to_vec(),
                depth,
                read_witness: witness_label(&obs.read_witness, locus_len),
                observed: obs.bases.to_vec(),
                reads: obs.num_obs,
            });
        }
    }

    /// The dump text: the two `#` header lines, the TSV column line, then the rows.
    fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        // Writing into a `String` via `fmt::Write` never fails, so the `write!` results are ignored.
        let _ = writeln!(
            out,
            "# ssr_loci={} zero_coverage={} reads_capped={} reads_without_observation={}",
            self.ssr_loci, self.zero_coverage, self.reads_capped, self.reads_without_observation
        );
        let _ = writeln!(
            out,
            "# obs_complete={} obs_partial={}",
            self.obs_complete, self.obs_partial
        );
        out.push_str(
            "contig\tstart\tend\tmotif\tref_tract\tdepth\tread_witness\tobserved\treads\n",
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.contig,
                row.start,
                row.end,
                String::from_utf8_lossy(&row.motif),
                String::from_utf8_lossy(&row.ref_tract),
                row.depth,
                row.read_witness,
                String::from_utf8_lossy(&row.observed),
                row.reads,
            );
        }
        out
    }
}

/// The tag a witness carries in the `read_witness` column.
fn witness_label(witness: &ReadWitness, locus_len: LocusLen) -> &'static str {
    // Since the reshape the side is a **derivation**, not a variant: a run flush with the left
    // border is a prefix constraint, one flush with the right border a suffix. A run flush with
    // neither is interior — the STR path cannot mint one (it anchors a border or yields nothing),
    // so it never appears here, but naming it keeps the label honest for the generic path.
    // Destructured rather than guarded on `_`, so a future `ReadWitness` variant is a
    // compile error here. The guard form is what this migration used and it is exactly what
    // let the compiler stop forcing these sites to be revisited.
    match witness {
        ReadWitness::Complete => "complete",
        run @ ReadWitness::Partial { .. } => {
            match (run.is_flush_left(), run.is_flush_right(locus_len)) {
                (true, _) => "partial:left",
                (false, true) => "partial:right",
                (false, false) => "partial:interior",
            }
        }
    }
}

/// Run the whole pipeline over `fasta` + `bams`, optionally restricted to `contig_filter` (a contig
/// name), delimiting reads with `aligner` (algorithm 3 or 4 — [`RepeatDelimiter`], monomorphised so
/// the per-read `align` is a direct call), building the [`DumpReport`]. `gen_config` is the STR
/// generator's config (the tests vary its cap). The reference's `.fai` is created if absent; the BAM
/// index likewise (`SampleReads`).
fn run_dump<A: RepeatDelimiter>(
    fasta: &Path,
    bams: &[PathBuf],
    contig_filter: Option<&str>,
    aligner: A,
    gen_config: SsrGeneratorConfig,
) -> Result<DumpReport, Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(&cache, fasta.to_path_buf())?;
    let contigs = info.contig_list();

    let sample = SampleReads::open_only_sample(bams, &info, ReadFilterConfig::default(), true)?;

    let walk_config = TypedRegionConfig::default();
    // The generator holds its own reference (margin fetch) and a factory (the per-file read query),
    // both windowed over the same FASTA — the reference seam the STR generator's doc calls the Arc
    // gap. Both are cheap: a path plus the contig table, nothing resident until a fetch.
    let mut generator = SsrGenerator::new(
        WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()),
        {
            let fasta = fasta.to_path_buf();
            let contigs = contigs.clone();
            move || WindowedRefSeq::new(fasta.clone(), contigs.clone())
        },
        aligner,
        gen_config,
        Bp(walk_config.criteria.bundle_threshold),
    )?;

    let mut report = DumpReport::default();
    for (index, entry) in contigs.entries.iter().enumerate() {
        if contig_filter.is_some_and(|name| entry.name != name) {
            continue;
        }
        // A fresh windowed reference per contig for the walk (it takes the reference by value).
        let walk_reference = WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone());
        let mut walk = TypedRegionIterator::over_contig(
            walk_reference,
            ContigId(index as u32),
            walk_config.clone(),
        )?;
        for region in walk.by_ref() {
            let region = region?;
            if let RegionKind::SsrSegment(segment) = &region.kind {
                generator.begin_segment(region.region);
                while let Some(locus) = generator.next_locus(segment, &sample)? {
                    report.push_locus(&locus, segment);
                }
            }
        }
        report.walk_ssr_loci += walk.counts().ssr_loci;
    }

    // The background `.fai` verification (only when a pre-existing `.fai` was used) is joined at the
    // end, so a stale index is a failure rather than a silently-wrong walk.
    if let Some(handle) = verify {
        handle.join()?;
    }

    let counts = generator.counts();
    report.reads_fetched = counts.reads_fetched;
    report.reads_capped = counts.reads_discarded_by_cap;
    report.obs_complete = counts.observations_complete;
    report.obs_partial = counts.observations_partial;
    // Every reason, named — `outside_tract` is the largest of the four on real data, so
    // summing three of them would report a fraction of the reads that yielded nothing.
    report.reads_without_observation = counts.no_border_anchored
        + counts.low_quality
        + counts.window_truncated
        + counts.outside_tract;
    Ok(report)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_ssr_loci_dump <reference.fa> <sample.bam|cram> [contig]\n\
             dumps, per microsatellite tract, the observed tract sequences one sample's reads showed.\n\
             delimiter: PVC_SSR_DELIMITER=unit-robust (default, algorithm 4u — the bake-off winner) | \
             unit-slip (algorithm 4) | flat-gap (algorithm 3, the production-parity port)."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let bam = PathBuf::from(&args[2]);
    let contig_filter = args.get(3).map(String::as_str);
    let config = SsrGeneratorConfig::default();

    // The delimiter is chosen once here and monomorphised into `run_dump` — the per-read `align` in
    // the walk is a static call, never a `dyn` one. Default: algorithm 4u (the recommended
    // unit-robust aligner, the delimiter bake-off winner); `unit-slip` is algorithm 4, `flat-gap` is
    // algorithm 3 (the byte-parity port) — set it for a side-by-side comparison.
    let delimiter = std::env::var("PVC_SSR_DELIMITER").unwrap_or_default();
    let emission = PerQualityEmission::new();
    let report = match delimiter.as_str() {
        "flat-gap" => run_dump(
            &fasta,
            &[bam],
            contig_filter,
            SsrFlatGapAligner::new(emission),
            config,
        ),
        "unit-slip" => run_dump(
            &fasta,
            &[bam],
            contig_filter,
            SsrUnitSlipAligner::new(emission),
            config,
        ),
        _ => run_dump(
            &fasta,
            &[bam],
            contig_filter,
            SsrUnitRobustAligner::new(emission),
            config,
        ),
    };

    match report {
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
    use std::fs::File;
    use tempfile::TempDir;

    // chr1 carries TWO clean (AC)×12 tracts, 40 bp apart: LEFT + tract1 + MID + tract2 + RIGHT.
    // Each tract (period 2, 12 copies, purity 1.0) is past the copy-number floor; every flank
    // (LEFT/MID/RIGHT) exceeds the 30 bp bundle threshold, and the tracts are 40 bp apart (>30), so
    // neither is bundled and both sit ≥30 bp from the contig ends — region typing emits both as
    // clean `SsrSegment`s. The reads only reach tract1 (all end ≤ position 74, clear of tract2's
    // query span 75..158), so tract2 is a genuine **zero-coverage** locus. The flanks are aperiodic
    // and their joins with the AC tracts ("…CTT|AC…", "…AC|CAG…", "…AC|GG…") do not extend the run.
    const LEFT: &[u8] = b"GATCTTGCAAGCTGGAATCCGTTACGATCGGATCAAGCTT";
    const TRACT: &[u8] = b"ACACACACACACACACACACACAC"; // (AC)×12 = 24 bp; tract1 = 41..64, tract2 = 105..128
    const MID: &[u8] = b"CAGTTGCACGATCCTAAGGCTTGACCATGGATCCAAGTTG";
    const RIGHT: &[u8] = b"GGTTCAAGATCCGGATCTTGCAATCGGATCAAGCTTGACT";

    fn contig() -> Vec<u8> {
        assert_eq!(LEFT.len(), 40);
        assert_eq!(TRACT.len(), 24);
        assert_eq!(MID.len(), 40);
        assert_eq!(RIGHT.len(), 40);
        [LEFT, TRACT, MID, TRACT, RIGHT].concat()
    }

    /// A reference-allele read at 1-based `start` spanning `len` bases — its sequence is the contig
    /// slice it maps to, so it carries zero mismatches and clears the mismatch filter.
    fn read(
        contig: &[u8],
        qname: &str,
        start: usize,
        len: usize,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_core::Position;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::cigar::op::Kind;
        use noodles_sam::alignment::record::{Flags, MappingQuality};
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
        let seq = contig[start - 1..start - 1 + len].to_vec();
        noodles_sam::alignment::RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![40u8; len]))
            .build()
    }

    /// Write `chr1` as a one-line FASTA (its `.fai` is created by `run_dump`).
    fn write_fasta(path: &Path, contig: &[u8]) {
        use std::io::Write as _;
        let mut file = File::create(path).unwrap();
        writeln!(file, ">chr1").unwrap();
        file.write_all(contig).unwrap();
        writeln!(file).unwrap();
    }

    /// Stamp a record with an `RG` tag, so a multi-read-group file can say which group each
    /// record belongs to (a file declaring several resolves **per record**).
    fn in_read_group(
        mut record: noodles_sam::alignment::RecordBuf,
        read_group: &str,
    ) -> noodles_sam::alignment::RecordBuf {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;
        record.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(read_group.as_bytes().to_vec().into()),
        );
        record
    }

    /// Write a coordinate-sorted single-contig BAM (`@RG SM`, no `@SQ M5`) holding `reads`,
    /// declaring one read group.
    fn write_bam(path: &Path, contig_len: usize, reads: &[noodles_sam::alignment::RecordBuf]) {
        write_bam_with_read_groups(path, contig_len, reads, &["rg0"]);
    }

    /// The same, declaring `read_groups` — all naming **one sample**, since a `SampleReads` open
    /// serves exactly one. With several declared, each record's own `RG` names which.
    fn write_bam_with_read_groups(
        path: &Path,
        contig_len: usize,
        reads: &[noodles_sam::alignment::RecordBuf],
        read_groups: &[&str],
    ) {
        use bstr::BString;
        use noodles_bam as bam;
        use noodles_sam as sam;
        use sam::alignment::io::Write as _;
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::Version;
        use sam::header::record::value::map::header::tag::SORT_ORDER;
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
        use std::num::NonZero;

        let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
        hd.other_fields_mut()
            .insert(SORT_ORDER, BString::from("coordinate"));
        let sq = Map::<ReferenceSequence>::new(NonZero::new(contig_len).unwrap());
        let mut builder = sam::Header::builder()
            .set_header(hd)
            .add_reference_sequence(b"chr1".to_vec(), sq);
        for name in read_groups {
            let mut rg = Map::<ReadGroup>::default();
            rg.other_fields_mut()
                .insert(SAMPLE, BString::from("sample0"));
            builder = builder.add_read_group(name.as_bytes().to_vec(), rg);
        }
        let header = builder.build();

        let mut writer = bam::io::Writer::new(File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        for record in reads {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    /// The committed fixture: the reference and a BAM whose reads produce, at **tract1**, four
    /// complete observations, one left-partial and one right-partial — while **tract2** is left
    /// uncovered (zero coverage). Every read ends at or before position 74, clear of tract2's query
    /// span (75..158). Reads are in coordinate order.
    fn fixture() -> (TempDir, PathBuf, PathBuf) {
        let contig = contig();
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("ref.fa");
        let bam = dir.path().join("sample.bam");
        write_fasta(&fasta, &contig);

        let reads = vec![
            // Four identical reads spanning tract1 + both its flanks (21..74) → complete "AC"×12.
            read(&contig, "c0", 21, 54),
            read(&contig, "c1", 21, 54),
            read(&contig, "c2", 21, 54),
            read(&contig, "c3", 21, 54),
            // Anchors tract1's left flank, runs off inside the tract (21..60) → partial:left,
            // reaching 20 of the tract's 24 positions.
            read(&contig, "pl", 21, 40),
            // The same, stopping earlier (21..52) → partial:left reaching 12. **Its reach is what
            // makes this fixture able to catch a left/right mix-up at all**: with only the
            // symmetric `pl`/`pr` pair below, the two partial rows carry identical bases and
            // counts, so swapping the sides merely exchanges their labels and the dump does not
            // move. A partial whose reach no other partial shares breaks that symmetry.
            read(&contig, "pl2", 21, 32),
            // Begins inside tract1, anchors its right flank (45..74) → partial:right.
            read(&contig, "pr", 45, 30),
        ];
        write_bam(&bam, contig.len(), &reads);
        (dir, fasta, bam)
    }

    /// The same seven reads, split across **two read groups** of one sample: the four completes
    /// go 3-to-`rg0` / 1-to-`rg1`, the two left partials to `rg0` and the right partial to `rg1`.
    ///
    /// Deliberately the *same reads* as [`fixture`], so the two dumps are comparable row for row
    /// and "the split sums back" is a statement about the same evidence rather than about two
    /// different fixtures.
    fn two_read_group_fixture() -> (TempDir, PathBuf, PathBuf) {
        let contig = contig();
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("ref.fa");
        let bam = dir.path().join("sample.bam");
        write_fasta(&fasta, &contig);

        let reads = vec![
            in_read_group(read(&contig, "c0", 21, 54), "rg0"),
            in_read_group(read(&contig, "c1", 21, 54), "rg0"),
            in_read_group(read(&contig, "c2", 21, 54), "rg0"),
            in_read_group(read(&contig, "c3", 21, 54), "rg1"),
            in_read_group(read(&contig, "pl", 21, 40), "rg0"),
            in_read_group(read(&contig, "pl2", 21, 32), "rg0"),
            in_read_group(read(&contig, "pr", 45, 30), "rg1"),
        ];
        write_bam_with_read_groups(&bam, contig.len(), &reads, &["rg0", "rg1"]);
        (dir, fasta, bam)
    }

    /// Dump with the default (recommended) delimiter, algorithm 4.
    fn dump(fasta: &Path, bam: &Path, config: SsrGeneratorConfig) -> DumpReport {
        dump_with(
            fasta,
            bam,
            SsrUnitSlipAligner::new(PerQualityEmission::new()),
            config,
        )
    }

    /// Dump with an explicit delimiter — the seam a bake-off drives.
    fn dump_with<A: RepeatDelimiter>(
        fasta: &Path,
        bam: &Path,
        aligner: A,
        config: SsrGeneratorConfig,
    ) -> DumpReport {
        run_dump(
            fasta,
            std::slice::from_ref(&bam.to_path_buf()),
            None,
            aligner,
            config,
        )
        .expect("the fixture dumps")
    }

    /// The two delimiters are **selectable** and agree on the clean fixture: algorithm 4 (the
    /// unit-slip default) and algorithm 3 (the flat-gap port) produce byte-identical output on
    /// reference-allele reads, so the choice is a real, comparable knob — not a behaviour fork on
    /// clean data (they diverge only on the harder alleles a bake-off exists to weigh).
    #[test]
    fn the_two_delimiters_are_selectable_and_agree_on_the_clean_fixture() {
        let (_dir, fasta, bam) = fixture();
        let unit_slip = dump_with(
            &fasta,
            &bam,
            SsrUnitSlipAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig::default(),
        );
        let flat_gap = dump_with(
            &fasta,
            &bam,
            SsrFlatGapAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig::default(),
        );
        assert_eq!(
            unit_slip.render(),
            flat_gap.render(),
            "algorithms 3 and 4 agree on clean reference-allele reads"
        );
    }

    /// Both AC tracts are emitted as loci (one per `SsrSegment`, including the uncovered tract2),
    /// and every fetched read is accounted for — a complete observation, a partial one, a cap
    /// discard, or a no-observation (spec §9.1, §9.2).
    #[test]
    fn every_segment_is_one_locus_and_every_read_is_accounted() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());

        assert_eq!(report.walk_ssr_loci, 2, "both AC tracts are detected");
        assert_eq!(
            report.ssr_loci, report.walk_ssr_loci,
            "one emitted locus per typed SsrSegment (spec §9.1)"
        );
        assert_eq!(
            report.zero_coverage, 1,
            "tract2 is present but uncovered (spec §9.1) — 'looked and saw nothing' ≠ 'never looked'"
        );
        assert_eq!(
            report.reads_fetched,
            report.obs_complete
                + report.obs_partial
                + report.reads_capped
                + report.reads_without_observation,
            "every fetched read is accounted for (spec §9.2)"
        );
        // tract1's seven reads land as four complete + three partial observations; tract2 gets
        // none.
        assert_eq!(report.obs_complete, 4);
        assert_eq!(report.obs_partial, 3);
        assert_eq!(report.reads_without_observation, 0);
        assert_eq!(report.reads_capped, 0);
    }

    /// The rendered text is exactly the spec §9 shape: two `#` header lines, the TSV column line,
    /// then a tab-separated row per observation. This pins the **format** — column order, the tabs,
    /// the header keys — which the structured assertions (reading `report.rows`) cannot catch.
    #[test]
    fn render_emits_the_spec_9_header_and_tsv_rows() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let text = report.render();
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(
            lines[0],
            "# ssr_loci=2 zero_coverage=1 reads_capped=0 reads_without_observation=0"
        );
        assert_eq!(lines[1], "# obs_complete=4 obs_partial=3");
        assert_eq!(
            lines[2],
            "contig\tstart\tend\tmotif\tref_tract\tdepth\tread_witness\tobserved\treads"
        );
        // The complete row, built from the report's own coordinates so the assertion pins the TSV
        // format (tabs, column order, values) without hard-coding region typing's tract bounds.
        let complete = report
            .rows
            .iter()
            .find(|row| row.read_witness == "complete")
            .expect("a complete row");
        let expected = format!(
            "chr1\t{}\t{}\tAC\tACACACACACACACACACACACAC\t4\tcomplete\tACACACACACACACACACACACAC\t4",
            complete.start, complete.end
        );
        assert!(
            text.lines().any(|line| line == expected),
            "the complete row should render as:\n{expected}\ngot:\n{text}"
        );
    }

    /// Partial observations exist — which proves the relevance gate admitted the partially-covering
    /// reads (a spanning-only gate would have dropped them, spec §9.4). Two left and one right, on
    /// their own rows, tagged and distinct from the complete rows.
    #[test]
    fn partial_observations_are_present_and_tagged() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());

        assert!(
            report.obs_partial >= 1,
            "the relevance gate admitted partials"
        );
        let left = report
            .rows
            .iter()
            .find(|row| row.read_witness == "partial:left")
            .expect("a left partial row");
        let right = report
            .rows
            .iter()
            .find(|row| row.read_witness == "partial:right")
            .expect("a right partial row");
        // A partial is a lower bound — shorter than the complete tract it sits under.
        assert!(left.observed.len() < left.ref_tract.len());
        assert!(right.observed.len() < right.ref_tract.len());
        // Depth on a partial row is the locus's *complete* depth (the four spanning reads).
        assert_eq!(left.depth, 4);
        assert_eq!(right.depth, 4);
    }

    /// The complete row is the reference tract, at depth 4, and the row for it carries the tract as
    /// both `ref_tract` and `observed`.
    #[test]
    fn the_complete_observation_is_the_reference_tract() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let complete = report
            .rows
            .iter()
            .find(|row| row.read_witness == "complete")
            .expect("a complete row");
        assert_eq!(complete.observed, TRACT);
        assert_eq!(complete.ref_tract, TRACT);
        assert_eq!(complete.motif, b"AC");
        assert_eq!(complete.reads, 4);
        assert_eq!(complete.depth, 4);
    }

    /// The output is byte-identical across repeated runs (spec §9.5).
    #[test]
    fn output_is_deterministic_across_runs() {
        let (_dir, fasta, bam) = fixture();
        let first = dump(&fasta, &bam, SsrGeneratorConfig::default()).render();
        let second = dump(&fasta, &bam, SsrGeneratorConfig::default()).render();
        assert_eq!(first, second);
    }

    /// Raising the cap above the deepest locus (6 reads) leaves the output unchanged; a cap *below*
    /// it changes the output — that is what a cap does, not a determinism failure (spec §9.5).
    #[test]
    fn a_cap_above_the_depth_is_invisible_and_a_cap_below_it_bites() {
        let (_dir, fasta, bam) = fixture();
        let uncapped = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: None,
            },
        )
        .render();
        let cap_above = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: Some(100),
            },
        )
        .render();
        assert_eq!(uncapped, cap_above, "a cap above the depth changes nothing");

        let cap_below = dump(
            &fasta,
            &bam,
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: Some(2),
            },
        );
        assert_ne!(
            uncapped,
            cap_below.render(),
            "a cap below the depth changes the output"
        );
        assert!(cap_below.reads_capped > 0, "the cap discarded reads");
    }

    /// **The fixture can tell a left partial from a right one** — which it could not until `pl2`
    /// was added, and the reason is worth keeping.
    ///
    /// With only the symmetric `pl`/`pr` pair, both partial rows carried identical `bases` and
    /// identical `reads`, differing solely in the label; sorting put the left-flush run first
    /// either way, so **swapping left for right at the mint site left the dump byte-identical**
    /// (checked, during B1, by applying the swap). The acceptance anchor was blind to the one
    /// property the `ReadWitness` reshape most affects.
    ///
    /// `pl2` breaks the symmetry by reaching a *different* number of tract positions, so the two
    /// left-flush rows carry bases no right-flush row shares and a swap moves the output. This
    /// test states the property the fixture must keep: **no two partial rows are
    /// interchangeable.**
    #[test]
    fn the_fixtures_partials_are_asymmetric_and_so_can_catch_a_side_swap() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let partials: Vec<&ObservationRow> = report
            .rows
            .iter()
            .filter(|row| row.read_witness.starts_with("partial"))
            .collect();

        assert_eq!(partials.len(), 3, "two left partials and one right");
        let left: Vec<&&ObservationRow> = partials
            .iter()
            .filter(|row| row.read_witness == "partial:left")
            .collect();
        let right: Vec<&&ObservationRow> = partials
            .iter()
            .filter(|row| row.read_witness == "partial:right")
            .collect();
        assert_eq!(left.len(), 2);
        assert_eq!(right.len(), 1);

        // The discriminating property: swapping the sides would have to move some row's bases
        // between the two labels, and that is only invisible if every left row has a right
        // counterpart with the same bases. Assert no such pairing exists.
        assert!(
            left.iter()
                .any(|l| right.iter().all(|r| r.observed != l.observed)),
            "at least one partial's bases must be unique to its side, or a side swap is invisible"
        );
    }

    /// **B2's oracle: the read group splits rows, and the split sums back.**
    ///
    /// The same seven reads as [`fixture`], dealt across two read groups of one sample. The
    /// complete allele was seen by three `rg0` reads and one `rg1` read, so it becomes **two
    /// rows** carrying 3 and 1 — where the single-group dump has one row of 4. Summing the group
    /// axis away must recover the single-group dump exactly, which is what makes the finer grain
    /// safe rather than merely conservative (spec §6).
    ///
    /// The row *count* rising is the half that proves the field is computed rather than
    /// defaulted: with `read_group` constant these reads would merge back into one row.
    #[test]
    fn two_read_groups_split_the_rows_and_the_counts_sum_back() {
        let (_single_dir, single_fasta, single_bam) = fixture();
        let single = dump(&single_fasta, &single_bam, SsrGeneratorConfig::default());

        let (_split_dir, split_fasta, split_bam) = two_read_group_fixture();
        let split = dump(&split_fasta, &split_bam, SsrGeneratorConfig::default());

        // The run-level totals are per read, not per row, so they cannot move: the same six
        // reads produced the same four complete and two partial observations.
        assert_eq!(single.obs_complete, split.obs_complete);
        assert_eq!(single.obs_partial, split.obs_partial);
        assert_eq!(single.reads_fetched, split.reads_fetched);

        assert_eq!(
            single.rows.len(),
            4,
            "one complete row and the three partials, ungrouped"
        );
        assert_eq!(
            split.rows.len(),
            5,
            "the complete allele splits in two; the partials were one read each already"
        );

        // Fold the group axis away and the two dumps must agree row for row.
        let collapse = |rows: &[ObservationRow]| {
            let mut totals: Vec<(Vec<u8>, &'static str, u32)> = Vec::new();
            for row in rows {
                match totals.iter_mut().find(|(bases, witness, _)| {
                    bases == &row.observed && *witness == row.read_witness
                }) {
                    Some((_, _, reads)) => *reads += row.reads,
                    None => totals.push((row.observed.clone(), row.read_witness, row.reads)),
                }
            }
            totals
        };
        assert_eq!(
            collapse(&split.rows),
            collapse(&single.rows),
            "collapsing the group axis recovers the single-group dump exactly"
        );
    }

    /// **A one-read-group fixture is unchanged by the field**, which is what "free at one read
    /// group" has to mean in practice: the row count is identical to a run that ignored the group
    /// entirely, and the rendered text does not move.
    ///
    /// Guards the direction the split could go wrong cheaply — splitting rows on a sample that
    /// has nothing to split.
    #[test]
    fn a_single_read_group_fixture_is_unchanged_by_the_group_axis() {
        let (_dir, fasta, bam) = fixture();
        let report = dump(&fasta, &bam, SsrGeneratorConfig::default());
        let distinct_cells: std::collections::HashSet<(Vec<u8>, &str)> = report
            .rows
            .iter()
            .map(|row| (row.observed.clone(), row.read_witness))
            .collect();
        assert_eq!(
            report.rows.len(),
            distinct_cells.len(),
            "with one read group no (bases, coverage) cell may appear twice"
        );
    }
}
