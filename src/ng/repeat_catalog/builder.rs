//! Building a catalog while the reference streams past: accumulate a contig, scan it whole,
//! write its rows.
//!
//! **Whole contigs, not windows** (spec §2.3). A scanner fed buffer-sized chunks needs a
//! margin carried across each one, a rule for which side a straddling detection belongs to,
//! and a cap on the repeat length it can promise to catch whole. Scanning a contig in one
//! slice removes all three: a satellite of any size comes out as one row. What it costs is
//! that contig resident while it is scanned — 90 MB for tomato's largest chromosome, 250 MB
//! for human chromosome 1.

use std::path::{Path, PathBuf};

use crate::ng::reference_info::{ContigInfo, ReferenceBasesObserver, ReferenceInfo};
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::parquet_file::RepeatCatalogWriter;
use crate::ng::repeat_catalog::row::{RowRejection, row_for_interval};
use crate::ng::repeat_catalog::{
    FoundRepeat, RepeatCatalogError, RepeatCatalogHeader, RowsByPeriod,
};
use crate::ng::tandem_repeat::{ScanParams, find_tandem_repeats};
use crate::ng::types::{Bp, ContigId};

/// What a build wrote, and what it did not.
///
/// The rejections are here because a count of what was kept, on its own, cannot say whether
/// the floors did what they were meant to (spec §4.4).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BuildTally {
    /// Rows written, per period.
    pub rows: RowsByPeriod,
    /// Detections dropped for holding fewer copies than the catalog's floor, measured on
    /// the **detected** span.
    pub below_copy_floor: u64,
    /// Detections dropped for sitting closer to a contig's end than the flank floor.
    pub too_close_to_contig_end: u64,
    /// Detections the scanner emitted that could not form a motif at their own period.
    /// Zero in practice; counted rather than assumed.
    pub malformed: u64,
}

impl BuildTally {
    fn charge(&mut self, rejection: RowRejection) {
        match rejection {
            RowRejection::CopyFloor => self.below_copy_floor += 1,
            RowRejection::TooCloseToContigEnd => self.too_close_to_contig_end += 1,
            RowRejection::NoMotif | RowRejection::PeriodOutOfRange => self.malformed += 1,
        }
    }
}

/// Builds a catalog from the reference pass's bases.
///
/// Attach one to [`crate::ng::reference_info::read_reference_info_observing`] and call
/// [`finish`](Self::finish) with the `ReferenceInfo` that pass returned. The header is
/// written from that value, so the digests a later reader checks are the ones this run
/// computed.
///
/// **Errors are held, not thrown.** The observer seam is infallible so that
/// `reference_info` stays a leaf (spec §2.2); a write that fails puts the builder into a
/// state where it ignores the rest of the reference, and `finish` returns the error.
pub struct RepeatCatalogBuilder {
    criteria: StrRepeatCriteria,
    scan: ScanParams,
    writer: Option<RepeatCatalogWriter>,
    path: PathBuf,

    /// The contig being accumulated: its index, its bases, and the rows found in it.
    contig: Option<ContigId>,
    bases: Vec<u8>,
    rows: Vec<FoundRepeat>,

    tally: BuildTally,
    failure: Option<RepeatCatalogError>,
}

impl RepeatCatalogBuilder {
    /// Open a builder writing to `path`.
    ///
    /// `criteria`'s period range, copy floors and minimum flank are what it applies; the
    /// rest are recorded in the header for provenance and applied by readers (spec §4.2).
    pub fn create(
        path: &Path,
        criteria: StrRepeatCriteria,
        scan: ScanParams,
    ) -> Result<Self, RepeatCatalogError> {
        Ok(Self {
            criteria,
            scan,
            writer: Some(RepeatCatalogWriter::create(path)?),
            path: path.to_path_buf(),
            contig: None,
            bases: Vec::new(),
            rows: Vec::new(),
            tally: BuildTally::default(),
            failure: None,
        })
    }

    /// Write the header and move the file into place.
    ///
    /// `reference` is what the pass this builder rode returned — its contig table and
    /// digests are the header's (spec §3.4). Fails if any write during the pass did, or if
    /// the pass saw contigs this builder did not.
    pub fn finish(
        mut self,
        reference: &ReferenceInfo,
        tool_version: &str,
    ) -> Result<BuildTally, RepeatCatalogError> {
        if let Some(failure) = self.failure.take() {
            return Err(failure);
        }
        let writer = self
            .writer
            .take()
            .expect("the writer is taken only here, and `failure` guards the other path");

        let header = RepeatCatalogHeader {
            contigs: reference.contigs.clone(),
            reference_md5: reference
                .md5
                .ok_or_else(|| RepeatCatalogError::Unreadable {
                    path: self.path.clone(),
                    source: "the reference was read without computing its digest, so the catalog \
                         would have nothing to be checked against"
                        .into(),
                })?,
            built_under: self.criteria.clone(),
            scan: self.scan,
            tool_version: tool_version.to_string(),
        };

        let rows = writer.finish(&header)?;
        debug_assert_eq!(rows, self.tally.rows, "the writer and the tally disagree");
        Ok(self.tally)
    }

    /// Scan the accumulated contig and write its rows as one row group.
    fn scan_and_write(&mut self, info: &ContigInfo) {
        let Some(contig) = self.contig else { return };
        let Some(writer) = self.writer.as_mut() else {
            return;
        };

        // One call over the whole contig: no window, so no margin and no length at which a
        // tract stops being caught whole (spec §2.3).
        let intervals = find_tandem_repeats(
            &self.bases,
            self.criteria.classification.periods,
            &self.scan,
        );

        self.rows.clear();
        self.rows.reserve(intervals.len());
        for interval in &intervals {
            match row_for_interval(
                contig,
                Bp(info.length),
                &self.bases,
                interval,
                &self.criteria,
            ) {
                Ok(row) => self.rows.push(row),
                Err(rejection) => self.tally.charge(rejection),
            }
        }

        // The file's order (spec §3.1). `find_tandem_repeats` emits period by period, so
        // rows arrive grouped by period rather than by position; the sort is what makes
        // the file's order a property of the format instead of of the detector's loop.
        self.rows
            .sort_by_key(|row| (row.detected.start, row.period, row.detected.end));
        for row in &self.rows {
            self.tally.rows.count(row.period);
        }

        if let Err(e) = writer.write_contig(&self.rows) {
            self.failure = Some(e);
            self.writer = None;
        }
    }
}

impl ReferenceBasesObserver for RepeatCatalogBuilder {
    fn contig_started(&mut self, _name: &str, index: usize) {
        self.contig = Some(ContigId(index as u32));
        // Load / use / clear: the buffer is reused across contigs rather than reallocated.
        self.bases.clear();
    }

    fn bases(&mut self, upper: &[u8]) {
        if self.failure.is_some() {
            return;
        }
        self.bases.extend_from_slice(upper);
    }

    fn contig_finished(&mut self, info: &ContigInfo) {
        if self.failure.is_some() {
            return;
        }
        debug_assert_eq!(
            self.bases.len() as u64,
            info.length,
            "the accumulated bases are the contig the pass just finished"
        );
        self.scan_and_write(info);
        // Nothing accumulates across contigs (spec §6).
        self.bases.clear();
        self.bases.shrink_to_fit();
        self.rows.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::RepeatCatalog;
    use crate::ng::types::Position;
    use std::io::Write;

    /// Write a FASTA of `(name, sequence)` contigs, 60 bases to the line.
    fn write_fasta(dir: &Path, contigs: &[(&str, String)]) -> PathBuf {
        let path = dir.join("ref.fa");
        let mut file = std::fs::File::create(&path).expect("create");
        for (name, seq) in contigs {
            writeln!(file, ">{name}").expect("header");
            for chunk in seq.as_bytes().chunks(60) {
                file.write_all(chunk).expect("bases");
                file.write_all(b"\n").expect("newline");
            }
        }
        path
    }

    /// Filler with no tandem structure of its own at periods 1..=6.
    fn filler(len: usize) -> String {
        const CYCLE: &[u8] = b"ACGTTGCAAGGTCCAT";
        (0..len).map(|i| CYCLE[i % CYCLE.len()] as char).collect()
    }

    fn build(dir: &Path, contigs: &[(&str, String)]) -> (PathBuf, ReferenceInfo, BuildTally) {
        let fasta = write_fasta(dir, contigs);
        let catalog_path = dir.join("ref.fa.repeats.parquet");
        let mut builder = RepeatCatalogBuilder::create(
            &catalog_path,
            StrRepeatCriteria::default(),
            ScanParams::default(),
        )
        .expect("builder");
        let reference = read_reference_info_observing(
            ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            },
            &mut builder,
        )
        .expect("the pass runs");
        let tally = builder.finish(&reference, "test-0.1").expect("finish");
        (catalog_path, reference, tally)
    }

    fn rows_of(path: &Path, reference: &ReferenceInfo) -> Vec<FoundRepeat> {
        let catalog =
            RepeatCatalog::open_checking_against_reference(path, reference).expect("opens");
        catalog
            .repeats_in_region(None)
            .expect("rows")
            .map(|r| r.expect("a row"))
            .collect()
    }

    /// The end-to-end shape: a reference with one planted tract yields one row, in the
    /// right place, with the right motif — and the catalog it wrote opens against the very
    /// reference the pass computed.
    #[test]
    fn a_planted_tract_becomes_a_row_at_its_own_coordinates() {
        let dir = tempfile::tempdir().expect("tmp");
        let seq = format!("{}{}{}", filler(100), "CAG".repeat(8), filler(100));
        let (path, reference, tally) = build(dir.path(), &[("chr1", seq)]);

        let rows = rows_of(&path, &reference);
        let planted: Vec<&FoundRepeat> = rows.iter().filter(|r| r.period == 3).collect();
        assert_eq!(planted.len(), 1, "one trimer tract: {rows:#?}");
        let row = planted[0];
        assert_eq!(row.motif.as_bytes(), b"CAG");
        assert_eq!(
            row.detected.start,
            Position(101),
            "1-based, after 100 filler bases"
        );
        assert_eq!(row.detected.end, Position(124), "24 bases of CAG");
        assert_eq!(row.trimmed, Some(row.detected));
        assert_eq!(row.stratum(), Some((3, 8)));
        assert_eq!(tally.rows.for_period(3), 1);
    }

    /// Spec §2.3's whole point: a tract far longer than any satellite cap is **one row**,
    /// not several, because a contig is scanned in one slice.
    #[test]
    fn a_two_kilobase_tract_comes_out_as_one_row() {
        let dir = tempfile::tempdir().expect("tmp");
        let seq = format!("{}{}{}", filler(100), "AT".repeat(1_000), filler(100));
        let (path, reference, _) = build(dir.path(), &[("chr1", seq)]);

        let rows = rows_of(&path, &reference);
        let long: Vec<&FoundRepeat> = rows
            .iter()
            .filter(|r| r.detected.len_bp() > 1_000)
            .collect();
        assert_eq!(long.len(), 1, "one 2 kb tract, whole: {long:#?}");
        // At least the planted 2,000 bases: the detector reaches into the flanking filler
        // when its last bases happen to continue the tiling, which is its business and not
        // this test's. What is this test's is that the tract is **one** row of that length
        // rather than several capped ones.
        assert!(
            long[0].detected.len_bp() >= 2_000,
            "the whole tract, in one row: {:?}",
            long[0].detected
        );
        assert_eq!(long[0].period, 2);
    }

    /// The flank floor at both contig ends, through the whole pipeline rather than in the
    /// row builder alone: a tract that starts at base 1 is not in the file.
    #[test]
    fn a_tract_at_a_contigs_very_start_is_not_in_the_file() {
        let dir = tempfile::tempdir().expect("tmp");
        let seq = format!("{}{}", "CAG".repeat(8), filler(200));
        let (path, reference, tally) = build(dir.path(), &[("chr1", seq)]);

        let rows = rows_of(&path, &reference);
        assert!(
            rows.iter().all(|r| r.detected.start.get() > 1),
            "a tract abutting base 1 has no flank: {rows:#?}"
        );
        assert!(
            tally.too_close_to_contig_end >= 1,
            "and it is counted, not silently gone: {tally:?}"
        );
    }

    /// Rows are written contig by contig, in reference order, and each contig is its own
    /// row group — which is what lets a reader ask for one contig.
    #[test]
    fn rows_come_out_in_reference_order_one_group_per_contig() {
        let dir = tempfile::tempdir().expect("tmp");
        let contigs = [
            (
                "chr1",
                format!("{}{}{}", filler(100), "CAG".repeat(8), filler(100)),
            ),
            (
                "chr2",
                format!("{}{}{}", filler(100), "AT".repeat(12), filler(100)),
            ),
        ];
        let (path, reference, _) = build(dir.path(), &contigs);

        let rows = rows_of(&path, &reference);
        let contig_order: Vec<u32> = rows.iter().map(|r| r.contig.get()).collect();
        let mut sorted = contig_order.clone();
        sorted.sort_unstable();
        assert_eq!(contig_order, sorted, "contigs in reference order");

        let catalog =
            RepeatCatalog::open_checking_against_reference(&path, &reference).expect("opens");
        let only_chr2: Vec<FoundRepeat> = catalog
            .repeats_in_region(Some(ContigId(1)))
            .expect("rows")
            .map(|r| r.expect("a row"))
            .collect();
        assert!(!only_chr2.is_empty());
        assert!(only_chr2.iter().all(|r| r.contig == ContigId(1)));
    }

    /// Within a contig, rows are ordered by start, then period, then end — the order the
    /// format promises, which the detector's period-by-period loop does not produce on its
    /// own.
    #[test]
    fn rows_within_a_contig_are_ordered_by_start_then_period() {
        let dir = tempfile::tempdir().expect("tmp");
        // Two tracts of different periods, the trimer first.
        let seq = format!(
            "{}{}{}{}{}",
            filler(100),
            "CAG".repeat(8),
            filler(100),
            "AT".repeat(12),
            filler(100)
        );
        let (path, reference, _) = build(dir.path(), &[("chr1", seq)]);

        let rows = rows_of(&path, &reference);
        let keys: Vec<(u64, u8, u64)> = rows
            .iter()
            .map(|r| (r.detected.start.get(), r.period, r.detected.end.get()))
            .collect();
        let mut sorted = keys.clone();
        sorted.sort_unstable();
        assert_eq!(keys, sorted, "rows: {rows:#?}");
    }

    /// A build over a reference read without digests has nothing for a reader to check
    /// against, and says so instead of writing a catalog nobody can validate.
    #[test]
    fn a_reference_without_a_digest_cannot_be_catalogued() {
        let dir = tempfile::tempdir().expect("tmp");
        let path = dir.path().join("catalog.parquet");
        let builder = RepeatCatalogBuilder::create(
            &path,
            StrRepeatCriteria::default(),
            ScanParams::default(),
        )
        .expect("builder");

        let digestless = ReferenceInfo {
            md5: None,
            contigs: Vec::new(),
            fasta_path: None,
        };
        assert!(matches!(
            builder.finish(&digestless, "test-0.1"),
            Err(RepeatCatalogError::Unreadable { .. })
        ));
    }
}
