//! Opening a catalog, checking it describes the reference this run is working on, and
//! reading its rows back.
//!
//! **The FASTA is the source of truth for every digest** (spec §4.3). The header does not
//! state what the reference is; it states what the builder claims it was, and the claim is
//! checked here against MD5s this run computed from the FASTA. Nothing is trusted because
//! it was written down.

use std::fs::File;
use std::path::{Path, PathBuf};

use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::ng::reference_info::{ContigInfo, ReferenceInfo};
use crate::ng::region_typing::segment_criteria::SsrSegment;
use crate::ng::region_typing::{GenomeRegions, TypedRegion};
use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::ng::repeat_catalog::parquet_file::{HEADER_KEY, decode_header, row_from_batch};
use crate::ng::repeat_catalog::segments::{
    loci_of_contig, segments_of_contig, segments_of_contig_in,
};
use crate::ng::repeat_catalog::strata::{StratumCounts, StratumSample, StratumSampler};
use crate::ng::repeat_catalog::{FoundRepeat, RepeatCatalogError, RepeatCatalogHeader};
use crate::ng::tandem_repeat::ScanParams;
use crate::ng::types::ContigId;
use crate::ng::types::{Bp, GenomeRegion};

/// An opened catalog: its header, already checked against this run's reference.
///
/// Holding one is the evidence that the file describes the reference in hand — which is why
/// [`open_checking_against_reference`](Self::open_checking_against_reference) is the only
/// way to get one, and why no row-reading method re-checks.
#[derive(Debug)]
pub struct RepeatCatalog {
    path: PathBuf,
    header: RepeatCatalogHeader,
}

impl RepeatCatalog {
    /// Open the catalog at `path` and check it against `reference`.
    ///
    /// Two failures, and they are kept apart because a caller can act on one and not the
    /// other (spec §5.3): **there is no catalog** — which names the command that writes
    /// one — or **the catalog is not this reference's**, which names the contig whose
    /// digest differs. Neither is a warning and neither triggers a rebuild: a short answer
    /// here is a genome-wide count taken over the wrong genome.
    ///
    /// Reads the footer only; no row is touched.
    pub fn open_checking_against_reference(
        path: &Path,
        reference: &ReferenceInfo,
    ) -> Result<Self, RepeatCatalogError> {
        if !path.exists() {
            return Err(RepeatCatalogError::NotFound {
                path: path.to_path_buf(),
            });
        }
        let header = read_header(path)?;
        check_against_reference(&header, reference)?;
        Ok(Self {
            path: path.to_path_buf(),
            header,
        })
    }

    /// What the file says about itself: the contig table, the criteria it was built under,
    /// the scoring weights and the tool version (spec §3.4).
    pub fn header(&self) -> &RepeatCatalogHeader {
        &self.header
    }

    /// The contig table — name, length and `@SQ M5` per contig, in reference order.
    pub fn contigs(&self) -> &[ContigInfo] {
        &self.header.contigs
    }

    /// Every tandem repeat the file holds for `contig`, or for the whole reference when
    /// `contig` is `None` — **exactly as stored, with no criteria applied**.
    ///
    /// Rows arrive in the file's order: contigs in reference order, then by start, then by
    /// period, then by end. The iterator streams a row group at a time, so peak memory is
    /// one contig's row group and not the genome.
    ///
    /// Restricting to a contig reads only that contig's row group, which is what one row
    /// group per contig buys (spec §3.5).
    pub fn repeats_in_region(
        &self,
        contig: Option<ContigId>,
    ) -> Result<
        impl Iterator<Item = Result<FoundRepeat, RepeatCatalogError>> + use<>,
        RepeatCatalogError,
    > {
        let file = File::open(&self.path).map_err(|source| RepeatCatalogError::Write {
            path: self.path.clone(),
            source,
        })?;
        let builder =
            ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| self.unreadable(e))?;

        // One row group per contig, in reference order, so a contig's rows are the row
        // group whose rows carry its index. Groups with no rows are never written.
        let wanted_groups: Vec<usize> = match contig {
            None => (0..builder.metadata().num_row_groups()).collect(),
            Some(id) => {
                let mut groups = Vec::new();
                for group in 0..builder.metadata().num_row_groups() {
                    if row_group_contig(&builder, group) == Some(id) {
                        groups.push(group);
                    }
                }
                groups
            }
        };

        let reader = builder
            .with_row_groups(wanted_groups)
            .build()
            .map_err(|e| self.unreadable(e))?;

        let path = self.path.clone();
        Ok(reader.flat_map(move |batch| {
            let path = path.clone();
            match batch {
                Err(e) => vec![Err(RepeatCatalogError::Unreadable {
                    path,
                    source: Box::new(e),
                })],
                Ok(batch) => (0..batch.num_rows())
                    .map(|i| {
                        row_from_batch(&batch, i).ok_or_else(|| RepeatCatalogError::Unreadable {
                            path: path.clone(),
                            source: format!("row {i} does not match the catalog schema").into(),
                        })
                    })
                    .collect::<Vec<_>>(),
            }
        }))
    }

    /// The genome's segments in coordinate order — the STR segments and the generic
    /// stretches between them — derived from the file with no FASTA open (spec §5.1).
    ///
    /// **The criteria are checked once, here, before any row is read.** A refusal is a fact
    /// about the file and the policy, not about a row, so it arrives at the call rather
    /// than part-way through a loop the caller has already started.
    ///
    /// The iterator then streams **one contig at a time**: a contig's rows and its segments
    /// are in memory, the genome's are not.
    pub fn genome_segments(
        &self,
        criteria: &StrRepeatCriteria,
        contig: Option<ContigId>,
    ) -> Result<GenomeSegments<'_>, RepeatCatalogError> {
        self.check_serves(criteria)?;
        Ok(GenomeSegments {
            catalog: self,
            criteria: criteria.clone(),
            contigs: self.contig_walk(contig).into_iter(),
            buffered: Vec::new().into_iter(),
        })
    }

    /// The genome's segments **inside a region set**, in coordinate order.
    ///
    /// The whole of each touched contig is typed and only then clipped, because a tract
    /// outside the requested span still bundles with one inside it and a satellite outside
    /// still swallows a locus inside — which is exactly what step 3's walk does with its
    /// scan span and its requested span.
    ///
    /// Contigs the region set does not name are not read at all.
    pub fn genome_segments_in(
        &self,
        criteria: &StrRepeatCriteria,
        wanted: &GenomeRegions,
    ) -> Result<Vec<TypedRegion>, RepeatCatalogError> {
        self.check_serves(criteria)?;

        let spans: Vec<GenomeRegion> = wanted.iter().collect();
        let mut out = Vec::new();
        for (contig, name, length) in self.contig_walk(None) {
            if !spans.iter().any(|s| s.contig == contig) {
                continue;
            }
            let rows = rows_of_contig(self, contig)?;
            out.extend(segments_of_contig_in(
                &name, contig, length, &rows, criteria, &spans,
            ));
        }
        Ok(out)
    }

    /// The contigs a read covers — all of them, or the one asked for — with the names and
    /// lengths the header carries.
    fn contig_walk(&self, contig: Option<ContigId>) -> Vec<(ContigId, String, Bp)> {
        self.header
            .contigs
            .iter()
            .enumerate()
            .map(|(i, info)| (ContigId(i as u32), info.name.clone(), Bp(info.length)))
            .filter(|(id, _, _)| contig.is_none_or(|wanted| *id == wanted))
            .collect()
    }

    /// The surviving STR loci alone, without the generic stretches between them.
    ///
    /// The same admission `genome_segments` runs — a locus here is a locus there — but the
    /// spans a caller would immediately discard are never built (spec §5.3). Streams one
    /// contig at a time, and refuses eagerly.
    pub fn str_loci(
        &self,
        criteria: &StrRepeatCriteria,
        contig: Option<ContigId>,
    ) -> Result<StrLoci<'_>, RepeatCatalogError> {
        self.check_serves(criteria)?;
        Ok(StrLoci {
            catalog: self,
            criteria: criteria.clone(),
            contigs: self.contig_walk(contig).into_iter(),
            buffered: Vec::new().into_iter(),
        })
    }

    /// How many loci the catalog holds in each `(period, repeat count)` stratum.
    ///
    /// The count is on the **trimmed** span — the count `classify` applies its copy floor to
    /// and the one a stratum is keyed by (spec §3.1). This is what the sampling consumer
    /// reweights by, so a short count here propagates into a diversity estimate rather than
    /// into a crash (spec §10.3).
    pub fn count_loci_per_stratum(
        &self,
        criteria: &StrRepeatCriteria,
        contig: Option<ContigId>,
    ) -> Result<StratumCounts, RepeatCatalogError> {
        let (counts, _) = self.sample_loci_per_stratum(criteria, contig, 0, 0)?;
        Ok(counts)
    }

    /// Up to `cap` loci from each stratum — the ones whose `hash(contig, start, seed)` is
    /// lowest — **and the counts from the same pass**.
    ///
    /// The pre-pass wants both, and the tally is a counter beside the heaps rather than a
    /// second traversal (spec §5.3). One forward pass, bounded state, and a result that does
    /// not depend on the order loci arrive in.
    ///
    /// The seed is the caller's: two runs at the same seed must keep the identical loci, and
    /// that is a property of the run rather than of the file.
    pub fn sample_loci_per_stratum(
        &self,
        criteria: &StrRepeatCriteria,
        contig: Option<ContigId>,
        cap: usize,
        seed: u64,
    ) -> Result<(StratumCounts, StratumSample), RepeatCatalogError> {
        let mut sampler = StratumSampler::new(cap, seed);
        for locus in self.str_loci(criteria, contig)? {
            let locus = locus?;
            let period = locus.motif().period() as u8;
            let repeat_count = (locus.end() - locus.start() + 1) / u64::from(period);
            sampler.offer(period, repeat_count, locus);
        }
        Ok(sampler.finish())
    }

    /// Whether this catalog can answer a reader asking with `criteria`, and whether it was
    /// scored the same way (spec §4.1, §4.2).
    ///
    /// The weights are checked for **equality**, not permissiveness: a different match
    /// reward or mismatch penalty is a different set of tracts, not a subset of this one,
    /// and no filter over the file reproduces it.
    fn check_serves(&self, criteria: &StrRepeatCriteria) -> Result<(), RepeatCatalogError> {
        self.header.built_under.serves(criteria)?;
        Ok(())
    }

    /// The scoring weights the file was built with, for a caller that scans as well as
    /// reads and must not mix the two.
    pub fn check_scored_with(&self, scan: &ScanParams) -> Result<(), RepeatCatalogError> {
        if self.header.scan.match_reward != scan.match_reward
            || self.header.scan.mismatch_penalty != scan.mismatch_penalty
        {
            return Err(RepeatCatalogError::ScoringWeightsDiffer {
                built_match_reward: self.header.scan.match_reward,
                built_mismatch_penalty: self.header.scan.mismatch_penalty,
                wanted_match_reward: scan.match_reward,
                wanted_mismatch_penalty: scan.mismatch_penalty,
            });
        }
        Ok(())
    }

    fn unreadable(
        &self,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> RepeatCatalogError {
        RepeatCatalogError::Unreadable {
            path: self.path.clone(),
            source: Box::new(source),
        }
    }
}

/// The genome's segments, one contig at a time.
///
/// Holding a contig's rows and its segments at once is the memory this costs; the genome's
/// are never all in hand, which is what "streams" means here (spec §5.3).
pub struct GenomeSegments<'a> {
    catalog: &'a RepeatCatalog,
    criteria: StrRepeatCriteria,
    contigs: std::vec::IntoIter<(ContigId, String, Bp)>,
    buffered: std::vec::IntoIter<TypedRegion>,
}

impl Iterator for GenomeSegments<'_> {
    type Item = Result<TypedRegion, RepeatCatalogError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(region) = self.buffered.next() {
                return Some(Ok(region));
            }
            let (contig, name, length) = self.contigs.next()?;
            // One contig's rows, whole, because that is what classification needs: bundling
            // reads a tract's neighbours, the satellite cap reads their merged coverage, and
            // the generic stretches are the gaps between the results. A row at a time would
            // answer none of them — and a contig is the unit the file is grouped by anyway.
            let rows = match rows_of_contig(self.catalog, contig) {
                Ok(rows) => rows,
                Err(e) => return Some(Err(e)),
            };
            self.buffered =
                segments_of_contig(&name, contig, length, &rows, &self.criteria).into_iter();
        }
    }
}

/// The STR loci alone, one contig at a time.
///
/// The same contig walk as [`GenomeSegments`], stopping one step earlier: `loci_of_contig`
/// resolves the repeat features and takes the loci, so the generic stretches between them —
/// a whole contig's worth of spans a caller would discard — are never built (spec §5.3).
pub struct StrLoci<'a> {
    catalog: &'a RepeatCatalog,
    criteria: StrRepeatCriteria,
    contigs: std::vec::IntoIter<(ContigId, String, Bp)>,
    buffered: std::vec::IntoIter<SsrSegment>,
}

impl Iterator for StrLoci<'_> {
    type Item = Result<SsrSegment, RepeatCatalogError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(locus) = self.buffered.next() {
                return Some(Ok(locus));
            }
            let (contig, name, length) = self.contigs.next()?;
            let rows = match rows_of_contig(self.catalog, contig) {
                Ok(rows) => rows,
                Err(e) => return Some(Err(e)),
            };
            self.buffered =
                loci_of_contig(&name, contig, length, &rows, &self.criteria).into_iter();
        }
    }
}

/// One contig's rows, in file order.
fn rows_of_contig(
    catalog: &RepeatCatalog,
    contig: ContigId,
) -> Result<Vec<FoundRepeat>, RepeatCatalogError> {
    catalog.repeats_in_region(Some(contig))?.collect()
}

/// Which contig a row group holds, from its first row's `contig` column.
fn row_group_contig(
    builder: &ParquetRecordBatchReaderBuilder<File>,
    group: usize,
) -> Option<ContigId> {
    use parquet::file::statistics::Statistics;

    let column = builder.metadata().row_group(group).column(0);
    match column.statistics()? {
        Statistics::Int32(s) => s.min_opt().map(|v| ContigId(*v as u32)),
        Statistics::Int64(s) => s.min_opt().map(|v| ContigId(*v as u32)),
        _ => None,
    }
}

/// Read and decode the footer's header block.
fn read_header(path: &Path) -> Result<RepeatCatalogHeader, RepeatCatalogError> {
    let unreadable =
        |source: Box<dyn std::error::Error + Send + Sync>| RepeatCatalogError::Unreadable {
            path: path.to_path_buf(),
            source,
        };

    let file = File::open(path).map_err(|e| unreadable(Box::new(e)))?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| unreadable(Box::new(e)))?;
    let text = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|entry| entry.key == HEADER_KEY))
        .and_then(|entry| entry.value.clone())
        .ok_or_else(|| unreadable("not a repeat catalog: no header in the footer".into()))?;

    decode_header(&text, path)
}

/// Compare the header's contig table against the digests this run computed from the FASTA.
///
/// The order matters as much as the contents: rows address contigs by index, so the same
/// sequences in a different order describe different coordinates. A length table alone
/// would not notice that, which is why the per-contig MD5s are stored (spec §3.4).
fn check_against_reference(
    header: &RepeatCatalogHeader,
    reference: &ReferenceInfo,
) -> Result<(), RepeatCatalogError> {
    if header.contigs.len() != reference.contigs.len() {
        return Err(RepeatCatalogError::ContigTableMismatch {
            detail: format!(
                "the catalog has {} contigs, the reference {}",
                header.contigs.len(),
                reference.contigs.len()
            ),
        });
    }

    for (stored, actual) in header.contigs.iter().zip(&reference.contigs) {
        if stored.name != actual.name {
            return Err(RepeatCatalogError::ContigTableMismatch {
                detail: format!(
                    "contig {} of the catalog is `{}`, of the reference `{}` \
                     (same contigs in a different order describe different coordinates)",
                    stored.name, stored.name, actual.name
                ),
            });
        }
        if stored.length != actual.length {
            return Err(RepeatCatalogError::ContigTableMismatch {
                detail: format!(
                    "contig {} is {} bases in the catalog and {} in the reference",
                    stored.name, stored.length, actual.length
                ),
            });
        }
        // A `.fai`-only read of the reference carries no digests, and a digest that was
        // never computed cannot disagree — the length and order checks above are what such
        // a run gets. A FASTA read always carries them.
        if let (Some(stored_md5), Some(actual_md5)) = (stored.md5, actual.md5)
            && stored_md5 != actual_md5
        {
            return Err(RepeatCatalogError::ContigDigestMismatch {
                contig: stored.name.clone(),
                found: hex(&actual_md5),
                expected: hex(&stored_md5),
            });
        }
    }

    if let Some(actual) = reference.md5
        && actual != header.reference_md5
    {
        return Err(RepeatCatalogError::ContigTableMismatch {
            detail: format!(
                "every contig matches but the whole-reference digest does not: \
                 {} here, {} in the catalog",
                hex(&actual),
                hex(&header.reference_md5)
            ),
        });
    }

    Ok(())
}

fn hex(digest: &[u8; 16]) -> String {
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::repeat_catalog::parquet_file::{RepeatCatalogWriter, temp_path_for};
    use crate::ng::repeat_catalog::{StrRepeatCriteria, TractSpan};
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::{Motif, Position};

    fn contig(name: &str, length: u64, md5: u8) -> ContigInfo {
        ContigInfo {
            name: name.to_string(),
            length,
            offset: 6,
            line_bases: 60,
            line_width: 61,
            md5: Some([md5; 16]),
        }
    }

    fn reference() -> ReferenceInfo {
        ReferenceInfo {
            md5: Some([42u8; 16]),
            contigs: vec![contig("chr1", 1_000, 7), contig("chr2", 500, 9)],
            fasta_path: None,
        }
    }

    fn header_for(reference: &ReferenceInfo) -> RepeatCatalogHeader {
        RepeatCatalogHeader {
            contigs: reference.contigs.clone(),
            reference_md5: reference.md5.expect("a digest"),
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test-0.1".to_string(),
        }
    }

    fn row(contig: u32, start: u64, period: u8) -> FoundRepeat {
        let end = start + 17;
        FoundRepeat {
            contig: ContigId(contig),
            detected: TractSpan {
                start: Position(start),
                end: Position(end),
            },
            trimmed: Some(TractSpan {
                start: Position(start),
                end: Position(end),
            }),
            period,
            score: 30,
            motif: Motif::new(b"CAG").expect("valid"),
            purity: Some(1.0),
        }
    }

    /// Write a catalog whose header describes `reference`, with two rows on chr1 and one on
    /// chr2.
    fn write_catalog(path: &Path, reference: &ReferenceInfo) {
        let mut writer = RepeatCatalogWriter::create(path).expect("writer");
        writer
            .write_contig(&[row(0, 100, 3), row(0, 400, 3)])
            .expect("chr1");
        writer.write_contig(&[row(1, 50, 3)]).expect("chr2");
        writer.finish(&header_for(reference)).expect("footer");
    }

    #[test]
    fn a_catalog_of_this_reference_opens_and_its_rows_come_back() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        let reference = reference();
        write_catalog(&path, &reference);

        let catalog =
            RepeatCatalog::open_checking_against_reference(&path, &reference).expect("opens");
        assert_eq!(catalog.contigs().len(), 2);
        assert_eq!(catalog.header().tool_version, "test-0.1");

        let rows: Vec<FoundRepeat> = catalog
            .repeats_in_region(None)
            .expect("rows")
            .map(|r| r.expect("a row"))
            .collect();
        assert_eq!(rows, vec![row(0, 100, 3), row(0, 400, 3), row(1, 50, 3)]);
    }

    /// One row group per contig means a contig's rows are read without touching the others.
    #[test]
    fn one_contig_can_be_read_on_its_own() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        let reference = reference();
        write_catalog(&path, &reference);

        let catalog =
            RepeatCatalog::open_checking_against_reference(&path, &reference).expect("opens");
        let rows: Vec<FoundRepeat> = catalog
            .repeats_in_region(Some(ContigId(1)))
            .expect("rows")
            .map(|r| r.expect("a row"))
            .collect();
        assert_eq!(rows, vec![row(1, 50, 3)]);
    }

    /// The two failures a caller reacts to differently must not be the same error.
    #[test]
    fn a_missing_catalog_names_the_command_that_builds_one() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("absent.parquet");

        let err = RepeatCatalog::open_checking_against_reference(&path, &reference())
            .expect_err("no file");
        assert!(matches!(err, RepeatCatalogError::NotFound { .. }));
        assert!(
            format!("{err}").contains("repeat-catalog"),
            "the message must say how to make one: {err}"
        );
    }

    #[test]
    fn a_catalog_of_another_reference_names_the_contig_that_differs() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        write_catalog(&path, &reference());

        // The same contigs, but chr2's bases changed.
        let mut edited = reference();
        edited.contigs[1].md5 = Some([99u8; 16]);

        let err = RepeatCatalog::open_checking_against_reference(&path, &edited)
            .expect_err("a different reference");
        match &err {
            RepeatCatalogError::ContigDigestMismatch { contig, .. } => {
                assert_eq!(contig, "chr2", "the message names which contig moved");
            }
            other => panic!("expected a digest mismatch, got {other}"),
        }
        assert!(format!("{err}").contains("chr2"));
    }

    /// Two contigs of equal length swapped: the lengths still match, so only the digests
    /// and the order catch it — and rows address contigs by index, so it matters.
    #[test]
    fn reordered_contigs_of_equal_length_are_caught() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");

        let same_length = ReferenceInfo {
            md5: Some([42u8; 16]),
            contigs: vec![contig("chrA", 1_000, 7), contig("chrB", 1_000, 9)],
            fasta_path: None,
        };
        write_catalog(&path, &same_length);

        let swapped = ReferenceInfo {
            contigs: vec![contig("chrB", 1_000, 9), contig("chrA", 1_000, 7)],
            ..same_length
        };
        let err = RepeatCatalog::open_checking_against_reference(&path, &swapped)
            .expect_err("the order changed");
        assert!(matches!(
            err,
            RepeatCatalogError::ContigTableMismatch { .. }
        ));
    }

    /// Every contig agreeing while the whole-reference digest does not means the reference
    /// is not the concatenation the catalog was built from.
    #[test]
    fn a_matching_contig_table_with_a_different_reference_digest_is_refused() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        write_catalog(&path, &reference());

        let elsewhere = ReferenceInfo {
            md5: Some([1u8; 16]),
            ..reference()
        };
        assert!(matches!(
            RepeatCatalog::open_checking_against_reference(&path, &elsewhere),
            Err(RepeatCatalogError::ContigTableMismatch { .. })
        ));
    }

    /// A build that died mid-write leaves a footerless `.tmp`; pointed at it, the reader
    /// refuses rather than reading a short catalog (spec §6).
    #[test]
    fn a_truncated_catalog_will_not_open() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");

        let mut writer = RepeatCatalogWriter::create(&path).expect("writer");
        writer.write_contig(&[row(0, 100, 3)]).expect("chr1");
        drop(writer);

        let leftover = temp_path_for(&path);
        assert!(matches!(
            RepeatCatalog::open_checking_against_reference(&leftover, &reference()),
            Err(RepeatCatalogError::Unreadable { .. })
        ));
    }

    /// A parquet file that is not a catalog has no header, and that is a refusal rather
    /// than an empty read.
    #[test]
    fn a_parquet_file_that_is_not_a_catalog_is_refused() {
        use arrow_array::{RecordBatch, UInt64Array};
        use arrow_schema::{DataType, Field, Schema};
        use parquet::arrow::ArrowWriter;
        use std::sync::Arc;

        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("something_else.parquet");
        let schema = Arc::new(Schema::new(vec![Field::new("x", DataType::UInt64, false)]));
        let batch = RecordBatch::try_new(
            Arc::clone(&schema),
            vec![Arc::new(UInt64Array::from(vec![1u64, 2, 3]))],
        )
        .expect("batch");
        let mut writer =
            ArrowWriter::try_new(File::create(&path).expect("create"), schema, None).expect("w");
        writer.write(&batch).expect("write");
        writer.close().expect("close");

        let err = RepeatCatalog::open_checking_against_reference(&path, &reference())
            .expect_err("not a catalog");
        assert!(matches!(err, RepeatCatalogError::Unreadable { .. }));
    }
}
