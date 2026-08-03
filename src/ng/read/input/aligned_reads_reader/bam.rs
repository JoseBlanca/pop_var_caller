//! Records from a BAM: the index says where to start, and reading carries on from there.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_sam as sam;

use crate::bam::index_preflight::AlignmentIndex;
use crate::ng::read::aligned_read::NoodlesRawAlignedRead;
use crate::ng::read::input::AlignmentFileError;
use crate::ng::types::GenomeRegion;

/// A BAM reader that stays where it is between regions.
///
/// **What it yields is undecoded** — a [`NoodlesRawAlignedRead`], which is a SAM flag and a
/// mapping quality readable without unpacking anything else. The name no longer says so, so
/// this does: nothing here builds an [`AlignedRead`](crate::ng::read::AlignedRead), and the
/// conversion happens above, only for the reads that clear the flag/MAPQ filter.
///
/// # The one thing this must get right
///
/// [`begin_region`](Self::begin_region) **positions; it never bounds.** After it, reading on
/// yields every record from that point to the end of the chromosome — not just the ones inside
/// the region it was handed. That is the contract in `aligned_reads_reader/mod.rs`, and it is
/// load-bearing rather than tidy: the cursor above serves a forward region by *not
/// repositioning at all*, so a reader that had quietly stopped at the previous region's end
/// would lose every record past it, silently, for every region after the first.
///
/// **This is exactly the trap noodles' `query()` sets.** That returns an iterator bounded by
/// the interval it was given, which satisfies every other word of the contract. So the index
/// is asked here for `region.start` **to the end of the contig**, and the chunks it answers
/// with are scanned in order. On a coordinate-sorted file that costs almost nothing extra: a
/// contig's records are contiguous, so the bins covering a suffix merge into a short run of
/// chunks rather than one per region.
pub(crate) struct BamAlignedReadsReader {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    /// Parsed once at open and shared, never re-read per region.
    header: Arc<sam::Header>,
    /// Also parsed once at open. Every variant is an `Arc` inside, so holding one here is a
    /// share rather than a copy of the index.
    index: AlignmentIndex,
    path: Arc<Path>,
    /// The chunks left to scan, from the last [`begin_region`](Self::begin_region).
    chunks: std::vec::IntoIter<Chunk>,
    /// Where the chunk being scanned ends, or `None` to take the next one.
    current_chunk_end: Option<bgzf::VirtualPosition>,
    /// Latched when the chunks run out or a read fails. Cleared by the next
    /// [`begin_region`](Self::begin_region) — and **only** by it, so a reader that has run off
    /// the end of the contig stays quiet until it is repositioned rather than restarting.
    done: bool,
}

/// Hand-written because noodles' reader is not `Debug`, and so the output says what
/// identifies this reader rather than dumping a parsed index.
impl std::fmt::Debug for BamAlignedReadsReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BamAlignedReadsReader")
            .field("path", &self.path)
            .field("chunks_left", &self.chunks.len())
            .field("done", &self.done)
            .finish_non_exhaustive()
    }
}

impl BamAlignedReadsReader {
    pub(crate) fn new(
        reader: bam::io::Reader<bgzf::io::Reader<File>>,
        header: Arc<sam::Header>,
        index: AlignmentIndex,
        path: Arc<Path>,
    ) -> Self {
        Self {
            reader,
            header,
            index,
            path,
            chunks: Vec::new().into_iter(),
            current_chunk_end: None,
            // Nothing has been asked for yet, so there is nothing to read. `begin_region` is
            // what starts this reader, exactly as it restarts it.
            done: true,
        }
    }

    pub(crate) fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Point the reader at `region.start`, and leave it able to read to the end of the contig.
    ///
    /// The index query is an in-memory lookup on the index parsed at open — never a re-parse,
    /// and never a file open.
    pub(crate) fn begin_region(&mut self, region: GenomeRegion) -> io::Result<()> {
        let chunks = self.chunks_from(region).map_err(io::Error::other)?;
        self.chunks = chunks.into_iter();
        self.current_chunk_end = None;
        self.done = false;
        Ok(())
    }

    /// The chunks holding every record from `region.start` to the end of its contig.
    fn chunks_from(&self, region: GenomeRegion) -> Result<Vec<Chunk>, AlignmentFileError> {
        let invalid_region = || AlignmentFileError::Region { region };

        // **There is no `region.is_empty()` check here, and its absence is the point.** One
        // was carried over from the per-region planner, where it mattered: that source
        // *bounded* its scan by `region.end`, so an inverted region asked the index a
        // nonsense interval and got an answer. This reader positions and never bounds —
        // `region.end` is not read below, only the contig's length is — so an inverted region
        // is a perfectly well-defined position for it and the check guarded nothing. It now
        // lives one layer up, where `end` is meaningful, and covers both formats and both
        // paths (`CursorError::InvalidRegion`).
        let reference_sequence_id =
            usize::try_from(region.contig.get()).map_err(|_| invalid_region())?;
        let (_, reference_sequence) = self
            .header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .ok_or_else(invalid_region)?;

        // **To the end of the contig, not to the end of the region** — see the type's doc.
        let start = noodles_core::Position::new(
            usize::try_from(region.start.get()).map_err(|_| invalid_region())?,
        )
        .ok_or_else(invalid_region)?;
        let end = noodles_core::Position::new(usize::from(reference_sequence.length()))
            .ok_or_else(invalid_region)?;
        if start > end {
            // A region beginning past the contig's last base. There is nothing to read, and
            // saying so is better than asking the index about an inverted interval.
            return Ok(Vec::new());
        }
        let interval: noodles_core::region::Interval = (start..=end).into();

        match &self.index {
            AlignmentIndex::BamCsi(index) => index.query(reference_sequence_id, interval),
            AlignmentIndex::BamBai(index) => index.query(reference_sequence_id, interval),
            // `open` pairs a `.crai` only with a CRAM, which reads through the CRAM arm.
            // Reaching here means the reader was built wrong — this module's bug, not the
            // caller's — so it must not masquerade as bad input.
            AlignmentIndex::Crai(_) => unreachable!(
                "a .crai index reached the BAM aligned-reads reader; \
                 open pairs .crai only with CRAM"
            ),
        }
        .map_err(|source| AlignmentFileError::IndexQuery {
            path: self.path.to_path_buf(),
            region,
            source,
        })
    }

    /// Test-only: the chunks a region resolves to, so the *positions, never bounds* contract
    /// can be checked directly rather than inferred from what a walk happened to return.
    #[cfg(test)]
    pub(crate) fn chunks_for(&self, region: GenomeRegion) -> Vec<Chunk> {
        self.chunks_from(region)
            .expect("the fixture region resolves")
    }

    /// Test-only: what the index answers for the region's **own** interval — the bounded
    /// query [`chunks_from`](Self::chunks_from) deliberately does not make.
    ///
    /// Exists so a test can show the two are different questions on a fixture whose index has
    /// the resolution to tell them apart. Without it, a chunk-equality assertion passes on any
    /// fixture too small to matter, which is exactly how this contract would go unchecked.
    #[cfg(test)]
    pub(crate) fn bounded_chunks_for(&self, region: GenomeRegion) -> Vec<Chunk> {
        let reference_sequence_id = region.contig.get() as usize;
        let start = noodles_core::Position::new(region.start.get() as usize).expect("1-based");
        let end = noodles_core::Position::new(region.end.get() as usize).expect("1-based");
        let interval: noodles_core::region::Interval = (start..=end).into();
        match &self.index {
            AlignmentIndex::BamCsi(index) => index.query(reference_sequence_id, interval),
            AlignmentIndex::BamBai(index) => index.query(reference_sequence_id, interval),
            AlignmentIndex::Crai(_) => unreachable!("BAM fixture"),
        }
        .expect("the fixture index answers")
    }

    /// The next record, raw: no contig test, no overlap test, no read-group resolution. Those
    /// belong to the layer above, written once for every format.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawAlignedRead) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        // Cleared and never set here: records come out raw, and a stale group left in the
        // reused buffer would attribute this record to the previous one's read group.
        buf.read_group = None;

        loop {
            // Land inside a chunk, seeking to its start if we have just taken it.
            let chunk_end = match self.current_chunk_end {
                Some(end) => end,
                None => {
                    let Some(chunk) = self.chunks.next() else {
                        self.done = true;
                        return Ok(false);
                    };
                    self.reader
                        .get_mut()
                        .seek(chunk.start())
                        .inspect_err(|_| self.done = true)?;
                    self.current_chunk_end = Some(chunk.end());
                    chunk.end()
                }
            };

            // The previous read may have carried us past this chunk's end.
            if self.reader.get_ref().virtual_position() >= chunk_end {
                self.current_chunk_end = None;
                continue;
            }

            let bytes_read = self
                .reader
                .read_record_buf(&self.header, &mut buf.record)
                .inspect_err(|_| self.done = true)?;
            if bytes_read == 0 {
                // End of file inside a chunk — try the next one.
                self.current_chunk_end = None;
                continue;
            }
            return Ok(true);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::index_preflight::load_alignment_index;
    use crate::ng::read::input::test_fixtures::{bam_header, indexed_bam, read_named_with_length};
    use crate::ng::types::{ContigId, Position};
    use noodles_sam::alignment::RecordBuf;

    const CONTIG_LENGTH: usize = 200_000;

    /// A contig long enough that the index has **resolution**: BAI's finest bins are 16 kb, so
    /// a fixture smaller than that resolves every region to the same single chunk and cannot
    /// tell a reader that positions from one that bounds. This one spans a dozen bins and
    /// several BGZF blocks.
    fn big_fixture() -> (tempfile::TempDir, std::path::PathBuf, sam::Header) {
        let header = bam_header(&[("chr1", CONTIG_LENGTH, None)]);
        let records: Vec<RecordBuf> = (0..2_000)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 90, 30))
            .collect();
        let (dir, path) = indexed_bam(&header, &records);
        (dir, path, header)
    }

    fn reader_over(path: &std::path::Path, header: &sam::Header) -> BamAlignedReadsReader {
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        reader.read_header().expect("read header");
        BamAlignedReadsReader::new(
            reader,
            Arc::new(header.clone()),
            load_alignment_index(path).expect("the fixture is indexed"),
            Arc::from(path),
        )
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **`begin_region` positions; it never bounds** — the contract in
    /// `aligned_reads_reader/mod.rs`, checked where it can actually fail.
    ///
    /// The cursor above serves a forward region by *not repositioning at all*, so a reader
    /// that stopped at the previous region's end would lose every record past it, silently,
    /// for every region after the first. noodles' `query()` returns exactly such a bounded
    /// iterator, which satisfies every other word of the contract.
    ///
    /// Checked on the chunks rather than on a walk, because a walk over a small fixture cannot
    /// tell the two apart: below 16 kb every region resolves to one chunk. Here the answer
    /// must not depend on `region.end` **at all**.
    #[test]
    fn begin_region_resolves_the_same_chunks_whatever_the_region_ends_at() {
        let (_dir, path, header) = big_fixture();
        let reader = reader_over(&path, &header);

        let from_the_start = reader.chunks_for(region(1, 100));

        // **The fixture must be able to tell the two apart**, or the assertion below passes
        // on any file too small for the index to have resolution — which is precisely how
        // this contract would go unchecked.
        assert_ne!(
            reader.bounded_chunks_for(region(1, 100)),
            from_the_start,
            "this fixture's index cannot distinguish a bounded query from a positioning one, \
             so nothing below is being tested",
        );

        for end in [100u64, 20_000, 100_000, CONTIG_LENGTH as u64] {
            assert_eq!(
                reader.chunks_for(region(1, end)),
                from_the_start,
                "the chunks changed with the region's end — the reader is bounding, not \
                 positioning",
            );
        }
    }

    /// And the same thing end to end: after being positioned at a narrow region, reading on
    /// must keep yielding records far past that region's end.
    #[test]
    fn reading_on_yields_records_far_past_the_region_it_was_positioned_for() {
        let (_dir, path, header) = big_fixture();
        let mut reader = reader_over(&path, &header);

        reader.begin_region(region(1, 100)).expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        let mut last_start = 0u64;
        let mut records = 0u64;
        while reader.read_next(&mut buf).expect("the fixture reads") {
            records += 1;
            if let Some(start) = buf.record.alignment_start() {
                last_start = usize::from(start) as u64;
            }
        }

        assert_eq!(
            records, 2_000,
            "reading on must reach every record on the contig"
        );
        assert!(
            last_start > 100_000,
            "reading on stopped at {last_start}, nowhere near the end of the contig",
        );
    }

    /// A reader that has run off the end stays quiet until it is repositioned, rather than
    /// restarting — which would hand the whole contig over a second time.
    #[test]
    fn a_reader_that_ran_out_stays_quiet_until_it_is_repositioned() {
        let (_dir, path, header) = big_fixture();
        let mut reader = reader_over(&path, &header);
        reader
            .begin_region(region(190_000, 200_000))
            .expect("positions");

        let mut buf = NoodlesRawAlignedRead::default();
        while reader.read_next(&mut buf).expect("the fixture reads") {}
        assert!(!reader.read_next(&mut buf).expect("still quiet"));

        reader.begin_region(region(1, 100)).expect("repositions");
        assert!(
            reader.read_next(&mut buf).expect("reads again"),
            "a repositioned reader must start again",
        );
    }

    /// A record comes out raw: no read group stamped, because resolving it belongs to the
    /// layer above and a stale one in the reused buffer would attribute this record to the
    /// previous one's group.
    #[test]
    fn a_record_comes_out_of_the_bam_arm_raw() {
        let (_dir, path, header) = big_fixture();
        let mut reader = reader_over(&path, &header);
        reader.begin_region(region(1, 100)).expect("positions");

        // Both fields spelled out: the spread would fill exactly one and would
        // silently absorb a third if the buffer ever grows one — on the very
        // literal whose subject is the field being spread past.
        let mut buf = NoodlesRawAlignedRead {
            record: RecordBuf::default(),
            read_group: Some(crate::ng::types::ReadGroupId(7)),
        };
        assert!(reader.read_next(&mut buf).expect("reads"));
        assert!(
            buf.read_group.is_none(),
            "the BAM arm must not stamp a read group"
        );
    }
}
