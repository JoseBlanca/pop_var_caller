//! Building a catalog while the reference streams past: accumulate a contig, scan it whole,
//! write its rows.
//!
//! **Whole contigs, not windows** (spec §2.3). A scanner fed buffer-sized chunks needs a
//! margin carried across each one, a rule for which side a straddling detection belongs to,
//! and a cap on the repeat length it can promise to catch whole. Scanning a contig in one
//! slice removes all three: a satellite of any size comes out as one row.
//!
//! **What it costs is set by the largest contig, at about 16 bytes a base** — measured at
//! 1.49 GB on tomato (largest contig 90 Mb) and 3.89 GB on GRCh38 (248 Mb), times the thread
//! count. Most of that is the scan's own working set rather than the bases, which are a
//! sixteenth of it (`doc/devel/reports/implementations/ng-repeat-catalog-E_2026-08-10.md`).

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

    /// The contig being accumulated: its index and its bases.
    contig: Option<ContigId>,
    bases: Vec<u8>,
    /// Contigs finished but not yet scanned — at most [`threads`](Self::with_threads) of
    /// them, which is what bounds the memory a parallel build uses (spec §2.4).
    pending: Vec<PendingContig>,
    threads: usize,

    tally: BuildTally,
    failure: Option<RepeatCatalogError>,
    /// How many contigs the pass has handed over, complete. Counted rather than inferred
    /// from [`longest_tract`](Self::longest_tract), which a contig holding no repeat at all
    /// never extends — so the two are different questions and only this one answers "did
    /// this builder see the reference `finish` is being given?".
    contigs_seen: usize,
    /// The longest tract written per contig, in contig order — the header carries it so a
    /// later reader can window its reads exactly (`RepeatCatalogHeader::reach_into_window`).
    longest_tract: Vec<u64>,
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
            pending: Vec::new(),
            threads: 1,
            tally: BuildTally::default(),
            failure: None,
            contigs_seen: 0,
            longest_tract: Vec::new(),
        })
    }

    /// Write the header and move the file into place.
    ///
    /// `reference` is what the pass this builder rode returned — its contig table and
    /// digests are the header's (spec §3.4). Fails if any write during the pass did, if the
    /// reference was read without computing its digests, or if the pass saw more contigs
    /// than this builder scanned.
    ///
    /// **The tool version is stamped here, not passed in.** A caller that could choose it
    /// could stamp a catalog with a version that never wrote it, and the reader's refusal
    /// (`check_written_by_this_build`) would then be checking a claim rather than a fact.
    pub fn finish(mut self, reference: &ReferenceInfo) -> Result<BuildTally, RepeatCatalogError> {
        // Whatever is still held back is scanned and written before the footer, in
        // reference order like everything else.
        self.flush_pending();
        if let Some(failure) = self.failure.take() {
            return Err(failure);
        }
        let writer = self
            .writer
            .take()
            .expect("the writer is taken only here, and `failure` guards the other path");

        // **The contig table and the rows must describe the same reference.** Rows address
        // contigs by index into this table, so a builder that saw fewer contigs than the
        // pass reports writes a file whose header promises contigs no row group holds — and
        // every read of those contigs comes back empty and `Ok`. It happens when a builder
        // is attached to one pass and finished with another pass's `ReferenceInfo`.
        if self.contigs_seen != reference.contigs.len() {
            return Err(RepeatCatalogError::ContigTableMismatch {
                detail: format!(
                    "the builder saw {} contigs and the reference reports {}; the catalog \
                     would name contigs it holds no rows for",
                    self.contigs_seen,
                    reference.contigs.len()
                ),
            });
        }

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
            tool_version: crate::ng::repeat_catalog::reader::running_tool_version().to_string(),
            // One number per contig, and it is what lets a reader ask for a stretch instead
            // of a whole contig: no row outside a window can reach further into it than the
            // longest tract there is (§*A windowed read*, `reach_into_window`).
            longest_tract_bp: {
                let mut per_contig = vec![0u64; reference.contigs.len()];
                for (index, longest) in self.longest_tract.iter().enumerate() {
                    if let Some(slot) = per_contig.get_mut(index) {
                        *slot = *longest;
                    }
                }
                per_contig
            },
        };

        let rows = writer.finish(&header)?;
        debug_assert_eq!(rows, self.tally.rows, "the writer and the tally disagree");
        Ok(self.tally)
    }

    /// Scan up to `threads` contigs at once instead of one at a time.
    ///
    /// **A speed knob and nothing else** (spec §2.4): contigs are written in the
    /// reference's order however many are scanned together, so the file is byte-identical
    /// at every thread count. What it costs is `threads` contigs resident rather than one.
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads.max(1);
        self
    }

    /// Scan every contig held back, and write their rows in **reference order**.
    ///
    /// Held-back contigs are scanned together when there is more than one thread; the
    /// results are collected in order before any is written, so completion order cannot
    /// reach the file.
    fn flush_pending(&mut self) {
        if self.pending.is_empty() || self.writer.is_none() {
            self.pending.clear();
            return;
        }

        let criteria = &self.criteria;
        let scan = &self.scan;
        let scanned: Vec<ScannedContig> = if self.threads == 1 {
            self.pending
                .iter()
                .map(|contig| scan_contig(contig, criteria, scan))
                .collect()
        } else {
            use rayon::prelude::*;
            self.pending
                .par_iter()
                .map(|contig| scan_contig(contig, criteria, scan))
                .collect()
        };

        for result in scanned {
            if let Some(contig) = result.rows.first().map(|row| row.contig.get() as usize) {
                if self.longest_tract.len() <= contig {
                    self.longest_tract.resize(contig + 1, 0);
                }
                let longest = result
                    .rows
                    .iter()
                    .map(|row| row.detected.len_bp())
                    .max()
                    .unwrap_or(0);
                self.longest_tract[contig] = self.longest_tract[contig].max(longest);
            }
            for rejection in result.rejections {
                self.tally.charge(rejection);
            }
            for row in &result.rows {
                self.tally.rows.count(row.period);
            }
            let Some(writer) = self.writer.as_mut() else {
                break;
            };
            if let Err(e) = writer.write_contig(&result.rows) {
                self.failure = Some(e);
                self.writer = None;
                break;
            }
        }

        // Nothing accumulates across flushes (spec §6).
        self.pending.clear();
        self.pending.shrink_to_fit();
    }
}

/// One contig's bases, held until it is scanned.
struct PendingContig {
    contig: ContigId,
    length: Bp,
    bases: Vec<u8>,
}

/// What scanning one contig produced. Rejections travel with the rows so that a parallel
/// scan's tally is assembled in reference order rather than in completion order.
struct ScannedContig {
    rows: Vec<FoundRepeat>,
    rejections: Vec<RowRejection>,
}

/// Detect every tandem repeat in one contig and turn each into a row.
///
/// Free-standing rather than a method because it is what runs on a worker: it borrows the
/// criteria and the weights, and touches nothing the builder mutates.
fn scan_contig(
    pending: &PendingContig,
    criteria: &StrRepeatCriteria,
    scan: &ScanParams,
) -> ScannedContig {
    // One call over the whole contig: no window, so no margin and no length at which a
    // tract stops being caught whole (spec §2.3).
    let intervals = find_tandem_repeats(&pending.bases, criteria.classification.periods, scan);

    let mut rows = Vec::with_capacity(intervals.len());
    let mut rejections = Vec::new();
    for interval in &intervals {
        match row_for_interval(
            pending.contig,
            pending.length,
            &pending.bases,
            interval,
            criteria,
        ) {
            Ok(row) => rows.push(row),
            Err(rejection) => rejections.push(rejection),
        }
    }

    // The file's order (spec §3.1). `find_tandem_repeats` emits period by period, so rows
    // arrive grouped by period rather than by position; the sort is what makes the file's
    // order a property of the format instead of of the detector's loop.
    rows.sort_by_key(|row| (row.detected.start, row.period, row.detected.end));
    ScannedContig { rows, rejections }
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
        let Some(contig) = self.contig else { return };
        self.contigs_seen += 1;

        self.pending.push(PendingContig {
            contig,
            length: Bp(info.length),
            // The bases move to the pending contig; the next one starts from an empty
            // buffer rather than reusing this one, since a worker may still be reading it.
            bases: std::mem::take(&mut self.bases),
        });
        if self.pending.len() >= self.threads {
            self.flush_pending();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::{ReadScope, RepeatCatalog};
    use crate::ng::types::GenomeRegion;
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
        build_with_threads(dir, contigs, 1, "ref.fa.repeats.parquet")
    }

    fn build_with_threads(
        dir: &Path,
        contigs: &[(&str, String)],
        threads: usize,
        catalog_name: &str,
    ) -> (PathBuf, ReferenceInfo, BuildTally) {
        let fasta = write_fasta(dir, contigs);
        let catalog_path = dir.join(catalog_name);
        let mut builder = RepeatCatalogBuilder::create(
            &catalog_path,
            StrRepeatCriteria::default(),
            ScanParams::default(),
        )
        .expect("builder")
        .with_threads(threads);
        let reference = read_reference_info_observing(
            ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            },
            &mut builder,
        )
        .expect("the pass runs");
        let tally = builder.finish(&reference).expect("finish");
        (catalog_path, reference, tally)
    }

    fn rows_of(path: &Path, reference: &ReferenceInfo) -> Vec<FoundRepeat> {
        let catalog =
            RepeatCatalog::open_checking_against_reference(path, reference).expect("opens");
        catalog.repeats(ReadScope::WholeReference).expect("rows")
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
            .repeats(ReadScope::Regions(&[GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(u64::MAX),
            }]))
            .expect("rows");
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

    /// Spec §10.5: `--threads` is a speed knob and must not reach the file. The contigs
    /// differ in size on purpose, so workers finish out of order and the writer has to put
    /// them back.
    #[test]
    fn the_thread_count_does_not_change_the_bytes() {
        let dir = tempfile::tempdir().expect("tmp");
        let contigs = [
            (
                "chr_big",
                format!(
                    "{}{}{}{}{}",
                    filler(4_000),
                    "CAG".repeat(20),
                    filler(4_000),
                    "AT".repeat(30),
                    filler(4_000)
                ),
            ),
            (
                "chr_small",
                format!("{}{}{}", filler(80), "TTA".repeat(9), filler(80)),
            ),
            (
                "chr_middling",
                format!("{}{}{}", filler(600), "GAAT".repeat(15), filler(600)),
            ),
        ];

        let mut bytes_at = Vec::new();
        let mut tallies = Vec::new();
        for threads in [1usize, 2, 4] {
            let (path, _reference, tally) = build_with_threads(
                dir.path(),
                &contigs,
                threads,
                &format!("t{threads}.parquet"),
            );
            bytes_at.push(std::fs::read(&path).expect("read back"));
            tallies.push(tally);
        }

        assert!(!bytes_at[0].is_empty());
        assert_eq!(bytes_at[0], bytes_at[1], "1 thread vs 2");
        assert_eq!(bytes_at[0], bytes_at[2], "1 thread vs 4");
        assert_eq!(tallies[0], tallies[1], "the tally is order-independent too");
        assert_eq!(tallies[0], tallies[2]);
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
            builder.finish(&digestless),
            Err(RepeatCatalogError::Unreadable { .. })
        ));
    }
}
