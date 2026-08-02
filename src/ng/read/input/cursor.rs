//! **The alignment cursor — a reader that stays where it is.** So far only its errors:
//! [`CursorError`]. The cursor itself lands in Milestone B, and the record readers beneath it
//! in the next step; the rest of this module doc is the problem being built for.
//!
//! Reading a sorted alignment file today opens a fresh query per region: the index is
//! consulted, the reader seeks, and the block it lands in is decompressed and decoded — and
//! then the next region, 390 bases along, does all of it again. Measured on chromosome 21
//! of a tandem-repeat-targeted HG002 run, **82 % of the seeks land in the block the reader
//! already holds**, and the same 35,228 records — that is the probe's **whole-contig** mode;
//! the typed-region walk counts 34,633, and spec §11.5 requires the mode to be named because
//! both figures are in circulation — are decoded 1,067,729 times. Consecutive
//! queries overlap by about 93 %, because the caller widens every region by 5,000 bases
//! against regions averaging 390.
//!
//! A cursor is the reader kept between regions instead: positioned in one chromosome of one
//! file, holding the reads it has already decoded *and filtered*, and handing them back when
//! the next region can use them.
//!
//! Design: `doc/devel/ng/spec/alignment_cursor.md` (what and why) and
//! `doc/devel/ng/arch/alignment_cursor.md` (types and interfaces). Build order in
//! `doc/devel/ng/impl_plan/alignment_cursor.md`.
//!
//! # What is here so far
//!
//! The cursor, over a scripted list of records — [`RecordReader::InMemory`]. The BAM arm
//! lands at Milestone C and CRAM at E; the sample-level merge that callers actually hold at
//! C4.
//!
//! # Two hazards this step found, for the step that adds the forget rule
//!
//! Both were reached by experiment during review, and both are invisible until kept reads
//! start surviving a reposition — which is exactly what B2 does.
//!
//! - **The overlap test must not be written twice.** A read whose CIGAR consumes no reference
//!   — all soft-clip — clears every step-1 filter and reaches this layer. noodles maps its
//!   zero span to *no* span, so `alignment_end()` reports the one-base footprint
//!   `start..=start` and [`RegionRecords`] **accepts** it. A second, hand-written test above
//!   the filter that treats a zero span as "touches nothing" would **reject** the same read —
//!   so it would be yielded when read fresh and dropped when replayed, which is a read lost
//!   with nothing failing. Whatever tests a kept read for overlap must apply the rule
//!   `RegionRecords` applies, not a second one that looks equivalent.
//! - **Reuse needs more than a comparison.** `RegionRecords::move_to` repositions the reader
//!   unconditionally, and the filter arm below does not check `kept` before pushing. So
//!   simply *not clearing* `kept` does not give reuse — it gives the same read twice. Spec §4's
//!   "partly held — hand over the kept reads, then carry on reading, no jump" has no code path
//!   yet; it needs the layers to agree on where reading resumes, not just a rule about what to
//!   keep.
//!
//! # The forget rule, which is the whole point and one comparison
//!
//! **Reuse what is held only when the new region begins at or after the last one served;
//! otherwise drop everything and jump.** Eviction is its mirror image: drop a kept read once
//! it ends before the current region begins, because every later region begins at or after
//! this one.
//!
//! That single test is sufficient, and the argument is short enough to check. Split the
//! records this region needs by where they begin. One beginning at or before the *last*
//! region's end reaches forward into this region, and this region begins at or after the last
//! one did — so it overlapped the last region too, was read then, and is held now. One
//! beginning after the last region's end was where the scan stopped, so the reader is sitting
//! on it and it is read forward with no jump. Every record is in one of those two groups.
//!
//! **There is nothing to tune and no index is consulted.** An earlier design derived a byte
//! cut-off from the index and was unsound three ways — `min_offset` collapses to byte 0 past
//! the last populated window on *both* index kinds, and a byte range is not a record set
//! (spec §6).

use std::collections::VecDeque;
use std::path::Path;
use std::sync::Arc;

use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::filtering::{ReadFilter, ReadFilterConfig, ReadFilterError, ReadGroupCounts};
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::record_reader::RecordReader;
use crate::ng::read::input::region_records::{RegionRecords, read_end, read_overlaps};
use crate::ng::ref_seq::{RawRefSeq, RefSeqError};
use crate::ng::types::{ContigId, GenomeRegion};

/// What can go wrong once a cursor exists.
///
/// Two conditions, and they are not peers. [`ReadRecord`](Self::ReadRecord) is the file
/// failing under a read that was asked for correctly; [`WrongChromosome`](Self::WrongChromosome) is a **caller bug** — a
/// guard, not a step in normal control flow. Correct code compares against the cursor's own
/// chromosome first and never sees it.
///
/// There is deliberately **no ordering variant**. Within its chromosome a cursor answers any
/// region in any order and the answer is always right; only the *speed* depends on how close
/// the region is to the last one (spec §4). Requiring regions to move forward was weighed
/// and rejected: a backward jump costs a seek and a block, which is what every request costs
/// today, so the restriction would protect against almost nothing while putting error
/// handling at every call site. Requiring a new cursor per chromosome was kept, because
/// there the cost prevented is chromosome-sized — on CRAM, re-reading hundreds of megabytes
/// of reference bases — and because nothing in a cursor survives a chromosome change anyway.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CursorError {
    /// A region on a chromosome this cursor does not cover. Make a cursor for that
    /// chromosome; this one is unharmed and still good for its own.
    ///
    /// **The contigs are reported as the numbers they are, and the file is named so the
    /// numbers mean something.** A cursor holds no name table, so it cannot say `chr21` — but
    /// an index is only interpretable against a particular table, and the path says which:
    /// [`AlignmentFile::contigs`](super::open_bam::AlignmentFile::contigs) on that file turns
    /// both numbers into names. Without it, a run holding one cursor per chromosome per
    /// generator per worker — 32 of them at one file per sample, 320 at ten — reports two bare
    /// integers and no way to tell which of the 320 produced them.
    #[error(
        "cursor on '{}' covers contig {} but the region is on contig {}",
        path.display(),
        cursor_contig.get(),
        requested_contig.get()
    )]
    WrongChromosome {
        path: Arc<Path>,
        /// The chromosome this cursor was made for.
        cursor_contig: ContigId,
        /// The chromosome the region the caller asked for lies on.
        requested_contig: ContigId,
    },

    /// The same read reached this sample from two of its files.
    ///
    /// **A sample's files must not overlap.** Two copies of one read are two votes at a
    /// locus: the depth doubles, the allele counts double, and every number derived from them
    /// is wrong in a way that looks like real evidence. The reachable route is a file listed
    /// twice under two names — a copy, or a symlink — which the path-based check at open
    /// cannot see.
    #[error(
        "the read '{}' appears in two of this sample's files (cursor {} and cursor {}) at \
         contig {} position {}",
        String::from_utf8_lossy(qname),
        first_file,
        second_file,
        contig.get(),
        position
    )]
    DuplicateReadAcrossFiles {
        qname: Vec<u8>,
        contig: ContigId,
        position: u64,
        first_file: usize,
        second_file: usize,
    },

    /// A read came out before the one before it.
    ///
    /// **The guarantee every layer above depends on and none of them re-checks.** The walker
    /// treats out-of-order input as a hard error, the merge's argmin is only sound over sorted
    /// inputs, and the sorted early stop below is only sound if the file really is sorted —
    /// the open gate proved the file *claims* to be, which is not the same thing. The
    /// per-region query this cursor replaces wrapped every stream in this check; a cursor that
    /// dropped it would be trusting a claim rather than the data.
    #[error(
        "alignment file '{}' yielded a read at position {} after one at {}, within one \
         region: the file is not coordinate-sorted",
        path.display(),
        position,
        after
    )]
    OutOfOrderRead {
        path: Arc<Path>,
        position: u64,
        after: u64,
    },

    /// A region was asked for after a read had already failed.
    ///
    /// **The alternative is worse than an error.** A cursor whose file failed cannot serve
    /// any later region, but it *can* still answer from the reads it happens to be holding —
    /// so without this the caller gets a plausible, silently short answer for every remaining
    /// region instead of being told the run is over. The failure was reported once, when it
    /// happened; this says the cursor has not recovered from it.
    #[error("cursor on '{}' cannot serve more regions: a read already failed", path.display())]
    AfterFailure { path: Arc<Path> },

    /// Reading the next record from the file failed.
    #[error("reading alignment file '{}' failed", path.display())]
    ReadRecord {
        path: Arc<Path>,
        #[source]
        source: std::io::Error,
    },
}

/// A reader positioned in one chromosome of one file, holding the reads it has recently
/// decoded and filtered so a nearby region can be answered without unpacking again.
///
/// **Not `Sync`, and that is a design statement rather than an omission.** An open file
/// position belongs to one consumer. Parallelism comes from more cursors — one per worker,
/// sharing nothing — never from sharing one (spec §3).
///
/// # What it promises
///
/// **Ask for any region on its chromosome, in any order, and the answer is right.** Whether
/// it is *fast* depends on how close the region is to the last one. A region on another
/// chromosome is refused, and refused before anything is touched, so the cursor is left
/// exactly as it was and is still good for its own (spec §10).
///
/// Nothing is unpacked ahead of demand: a caller that pulls one read and moves elsewhere has
/// unpacked at most one block, and abandoning a region costs nothing to unwind because there
/// is no stream object to give back.
pub struct AlignmentCursor<R: RawRefSeq> {
    /// The whole chain below, owned: the filter holds [`RegionRecords`], which holds the
    /// [`RecordReader`]. Not a cycle — the filter's source is the layer *below* this cursor,
    /// not the cursor itself.
    ///
    /// The reference accessor the mismatch filter needs sits inside, taken once here rather
    /// than rebuilt per query (perf review L2).
    filter: ReadFilter<RegionRecords, R>,
    /// **Our** reads — decoded and filtered — in the order they came off the file.
    ///
    /// Held above the filter, so serving one again skips both decode and filtering: that is
    /// the whole saving, since a read is returned by about a dozen consecutive regions and
    /// would otherwise be turned into an `AlignedRead` a dozen times (spec §5).
    kept: VecDeque<AlignedRead>,
    /// How far into `kept` the region being served has looked — the index of the next kept
    /// read to consider. **Not a drain**: a read this region has already been given may still
    /// be owed to the next one, so it stays where it is and only this marker moves.
    examined: usize,
    /// Where the last region served began. **This one number is the entire forget rule**
    /// (spec §6). `None` before the first region.
    last_region_start: Option<u64>,
    /// The region being served, or `None` before the first [`move_to_region`](Self::move_to_region).
    region: Option<GenomeRegion>,
    contig: ContigId,
    /// The position of the last read handed to the caller **for the region being served**,
    /// or `None` before the first. The order guard, per region — see [`next_read`](Self::next_read).
    last_emitted: Option<u64>,
    /// Named in errors, so two bare contig numbers mean something in a run holding hundreds
    /// of cursors over many files.
    path: Arc<Path>,
    counts: CursorCounts,
}

/// What a cursor did, as opposed to what it returned.
///
/// **The saving this whole design exists for is invisible in the reads.** A cursor that keeps
/// nothing and one that keeps everything hand back exactly the same reads for exactly the
/// same regions — that is the correctness requirement — so the only way to tell whether it is
/// working is to count what it *avoided*. On the fixture the perf review measured, the same
/// 35,228 records were decoded 1,067,729 times; the number that has to come down is
/// `reads_decoded`, and nothing else in the output moves when it does.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CursorCounts {
    /// Reads that came up through the filter — decoded and step-1 filtered exactly once each,
    /// if the rule is working. **This is the number the feature is about.**
    pub reads_decoded: u64,
    /// Reads handed to a caller from what was already held, skipping decode and filtering
    /// altogether. The saving, counted from the other side.
    pub reads_replayed: u64,
    /// Regions that reused what was held, and regions that dropped it and repositioned. They
    /// sum to the regions asked for, so a walk that quietly stopped reusing shows up as a
    /// ratio rather than as a slower run nobody attributes.
    pub regions_reusing: u64,
    pub regions_jumping: u64,
    /// Kept reads dropped because they ended before a region began. Without this, "the kept
    /// set is bounded" is an argument rather than an observation.
    pub reads_evicted: u64,
}

/// Summing two cursors' tallies — a sample's k files, or a run's chromosomes.
///
/// **It is an impl rather than five `+=` lines at each site because there are four such
/// sites**, and review found that dropping four of the five fields at one of them left the
/// whole suite green: `regions_reusing` and `regions_jumping` are the only observable that
/// says whether the cursor is being *kept*, and a fold that quietly loses them makes the
/// feature undetectable again on every chromosome but the live one. A sixth field added to
/// this struct is now folded everywhere by construction instead of in three places out of
/// four.
impl std::ops::AddAssign for CursorCounts {
    fn add_assign(&mut self, other: Self) {
        // Exhaustive destructure, no `..`: a new field is a compile error here until it is
        // folded, which is the property this impl exists to make free.
        let Self {
            reads_decoded,
            reads_replayed,
            regions_reusing,
            regions_jumping,
            reads_evicted,
        } = other;
        self.reads_decoded += reads_decoded;
        self.reads_replayed += reads_replayed;
        self.regions_reusing += regions_reusing;
        self.regions_jumping += regions_jumping;
        self.reads_evicted += reads_evicted;
    }
}

impl<R: RawRefSeq> AlignmentCursor<R> {
    /// A cursor over a scripted list of records, with no file behind it.
    ///
    /// The oracle's constructor: what a region *should* return is answerable by scanning the
    /// same list by hand, which is what makes this the thing the forget rule is judged
    /// against before an indexed file can hide a defect in it.
    pub(crate) fn over_records(
        reader: RecordReader,
        contig: ContigId,
        resolution: ReadGroupResolution,
        reference: R,
        config: ReadFilterConfig,
        path: Arc<Path>,
    ) -> Result<Self, RefSeqError> {
        let records = RegionRecords::new(reader, contig, resolution);
        Ok(Self {
            filter: ReadFilter::new(records, reference, config)?,
            kept: VecDeque::new(),
            examined: 0,
            last_region_start: None,
            last_emitted: None,
            region: None,
            contig,
            path,
            counts: CursorCounts::default(),
        })
    }

    /// The chromosome this cursor covers.
    ///
    /// A caller compares against this and mints a new cursor at a chromosome boundary; the
    /// error below exists as a guard against a bug, not as a step in normal control flow.
    pub fn contig(&self) -> ContigId {
        self.contig
    }

    /// Point the cursor at `region`.
    ///
    /// **The chromosome is checked first, before any state moves**, so a refusal leaves the
    /// cursor exactly as it was — which is what makes "unharmed and still good for its own"
    /// true by construction rather than by care (spec §10).
    pub fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), CursorError> {
        if region.contig != self.contig {
            return Err(CursorError::WrongChromosome {
                path: Arc::clone(&self.path),
                cursor_contig: self.contig,
                requested_contig: region.contig,
            });
        }

        // A cursor whose file has failed serves nothing further. Checked here, after the
        // chromosome test and before anything moves, so a dead cursor says so rather than
        // answering later regions out of whatever it is still holding.
        if self.filter.has_failed() {
            return Err(CursorError::AfterFailure {
                path: Arc::clone(&self.path),
            });
        }

        // **The forget rule, and it is one comparison** (spec §6). Reuse what is held only
        // when the new region begins at or after the last one served; otherwise drop it all
        // and jump.
        //
        // Why that single test is enough, in the two cases every needed record falls into.
        // Take a record this region needs. Either it begins at or before the *last* region's
        // end — and since this region begins at or after the last one did, such a record
        // overlapped the last region too, so it was read then and is held now. Or it begins
        // after the last region's end — and the scan stopped at the first of those, so the
        // reader is sitting on it and it is read forward, with no jump.
        //
        // The rule consults no index and has nothing to tune. An earlier design derived a
        // byte cut-off from the index instead and was unsound three ways: `min_offset`
        // collapses to byte 0 past the last populated window on **both** index kinds, and a
        // byte range is not a record set (spec §6).
        let reuse = self
            .last_region_start
            .is_some_and(|last| region.start.get() >= last);

        if reuse {
            // Eviction is the mirror image of the rule: a read that ends before this region
            // begins cannot touch this region or any later one, because every later region
            // begins at or after this one. Dropped from the front, which is the oldest end,
            // so this stays a walk rather than a scan.
            while self
                .kept
                .front()
                .is_some_and(|read| read_end(read) < region.start.get())
            {
                self.kept.pop_front();
                self.counts.reads_evicted += 1;
            }
            self.counts.regions_reusing += 1;
        } else {
            self.counts.reads_evicted += self.kept.len() as u64;
            self.kept.clear();
            self.counts.regions_jumping += 1;
        }
        // Every kept read is offered to the new region, including ones the last region was
        // already given: consecutive regions overlap, and a read touching both is owed to
        // both.
        self.examined = 0;
        self.region = Some(region);
        self.last_region_start = Some(region.start.get());
        // Per region, because a new region rewinds through what is kept: positions ascend
        // *within* a region, never across one.
        self.last_emitted = None;

        // Undo the *clean* stop the previous region ended with. Without this the first region
        // a cursor drains silences it for the whole chromosome: a region boundary reaches the
        // filter as an ordinary end of input, because that is the only end a `RecordSource`
        // can report. A filter stopped by a **fatal error** is not restarted, and must not be.
        self.filter.restart_after_end_of_input();

        if reuse {
            self.filter.source_mut().continue_into(region);
            Ok(())
        } else {
            self.filter
                .source_mut()
                .jump_to(region)
                .map_err(|source| CursorError::ReadRecord {
                    path: Arc::clone(&self.path),
                    source,
                })
        }
    }

    /// The next read of the current region, or `None` at the end of it.
    ///
    /// Shaped as [`Iterator::next`] so `Iterator` stays available to build on later. A cursor
    /// that has not been pointed anywhere yields `None` rather than guessing at a region.
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>> {
        // Not pointed anywhere yields nothing rather than guessing at a region. The layer
        // below would answer the same, but a type is better for stating its own contract than
        // for leaning on a neighbour's.
        let region = self.region?;

        // **What is already held**, in the order it came off the file — so a read this region
        // shares with the last one skips decode and filtering entirely. That is the saving:
        // a read is returned by about a dozen consecutive regions and would otherwise be
        // turned into an `AlignedRead` a dozen times.
        while let Some(read) = self.kept.get(self.examined) {
            // Held reads are in position order, so once one begins past this region's end,
            // neither it nor any later one can reach back into it, and neither can anything
            // still in the file. The answer is complete.
            //
            // **A short-circuit, and — traced carefully — nothing more.** Walking on instead
            // of stopping would skip every one of these reads anyway, since overlap needs
            // `pos <= region.end`; it would then reach the layer below, which either hands
            // back the record its own early stop is already holding (and stops again) or
            // finds the reader at the end. Same answer, same reads decoded, same counters.
            //
            // So mutating `return None` into `continue` fails no test *today*, and the
            // reason is worth stating exactly: kept reads never decrease in position, so
            // reaching the layer below would only make it early-stop for the same reason.
            //
            // **It is not, however, unobservable in principle**, and an earlier version of
            // this comment overclaimed that. Walking on consumes records off the reader that
            // stopping leaves alone — and once the reader is a real file rather than a
            // scripted list, that read can fail, turning a clean `None` into an error and a
            // permanently dead filter. So this is a short-circuit that also narrows what can
            // go wrong, not merely a saving.
            if read.pos > region.end.get() {
                return None;
            }
            self.examined += 1;
            if read_overlaps(read, region) {
                let read = read.clone();
                self.counts.reads_replayed += 1;
                return Some(self.emit(read));
            }
        }

        // Then read on, from wherever the reader is. Everything the filter yields already
        // overlaps the region — `RegionRecords` narrowed it below — so what is kept here is
        // exactly what a later region may be able to reuse.
        match self.filter.next() {
            None => None,
            Some(Err(error)) => Some(Err(self.read_failure(error))),
            Some(Ok(read)) => {
                self.counts.reads_decoded += 1;
                self.kept.push_back(read.clone());
                self.examined = self.kept.len();
                Some(self.emit(read))
            }
        }
    }

    #[allow(
        dead_code,
        reason = "the memory bound this reports is what Milestone D measures and what the \
                  deferred per-chromosome reference registry is triggered by; until then it \
                  is read only by this module's tests"
    )]
    /// How many reads this cursor is holding. The memory bound made observable, and what the
    /// eviction rule is judged on.
    pub(crate) fn kept_reads(&self) -> usize {
        self.kept.len()
    }

    /// What this cursor did — see [`CursorCounts`], and read `reads_decoded` first.
    pub fn counts(&self) -> CursorCounts {
        self.counts
    }

    /// Step-1's per-read-group tally, for as much of the chromosome as this cursor has read.
    ///
    /// **It needs no field here and no hand-over**, which is the point arch §2.3 makes: the
    /// filter has been keeping a running tally all along, and now that it lives as long as the
    /// cursor rather than as long as one region, reading it at any moment gives a whole-cursor
    /// total. The per-query sources this replaces had to fold their counts back into the file
    /// as each stream ended, or the drops they recorded vanished with them.
    pub fn read_group_counts(&self) -> Vec<ReadGroupCounts> {
        self.filter.counts()
    }

    /// Hand a read to the caller, checking it does not go backwards.
    ///
    /// The order guard the per-region query kept in `OrderVerified` and this cursor would
    /// otherwise have dropped. Scoped **to the region**: a new region rewinds through what is
    /// kept, so positions ascend within a region and not across one.
    fn emit(&mut self, read: AlignedRead) -> Result<AlignedRead, CursorError> {
        if let Some(last) = self.last_emitted
            && read.pos < last
        {
            return Err(CursorError::OutOfOrderRead {
                path: Arc::clone(&self.path),
                position: read.pos,
                after: last,
            });
        }
        self.last_emitted = Some(read.pos);
        Ok(read)
    }

    /// A filter failure, named with the file it came from.
    ///
    /// `ReadFilterError` says what went wrong and not *where*: the filter does not know which
    /// file it is reading, and in a run holding hundreds of cursors that is the first thing
    /// an operator needs.
    fn read_failure(&self, error: ReadFilterError) -> CursorError {
        CursorError::ReadRecord {
            path: Arc::clone(&self.path),
            source: std::io::Error::other(error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::record_reader::InMemoryRecordReader;
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, bam_header, fixture_read_group, matching_contigs, read_named_with_length,
    };
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::Position;
    use noodles_sam::alignment::RecordBuf;

    /// An all-`A` reference over the fixture contigs, so every fixture read matches perfectly
    /// and nothing is dropped for mismatching — this milestone is about *which* reads come
    /// back, not about filtering.
    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A cursor over a scripted list, on contig 0.
    fn cursor_over(records: Vec<RecordBuf>) -> AlignmentCursor<InMemoryRefSeq> {
        AlignmentCursor::over_records(
            RecordReader::InMemory(InMemoryRecordReader::new(
                bam_header(&matching_contigs()),
                records,
            )),
            ContigId(0),
            fixture_read_group(),
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        )
        .expect("the fixture header resolves against the in-memory reference")
    }

    /// A 30-base read on contig 0 — long enough to clear the default minimum read length,
    /// which a shorter one is silently dropped by, and all `A` so it matches the reference
    /// everywhere.
    const READ_LENGTH: u64 = 30;

    fn read_at(qname: &str, start: usize) -> RecordBuf {
        read_named_with_length(qname, 0, start, READ_LENGTH as usize)
    }

    /// Every read a region yields, by name.
    fn reads_of(cursor: &mut AlignmentCursor<InMemoryRefSeq>, region: GenomeRegion) -> Vec<String> {
        cursor
            .move_to_region(region)
            .expect("the region is on this cursor's chromosome");
        let mut names = Vec::new();
        while let Some(read) = cursor.next_read() {
            let read = read.expect("the scripted reads decode");
            names.push(String::from_utf8_lossy(&read.qname).into_owned());
        }
        names
    }

    /// **What a region should return, answered by hand.** A linear scan of the same script,
    /// with the same overlap test — the oracle every assertion below is stated against, so a
    /// cursor and the answer it is compared with cannot share a mistake.
    fn by_linear_scan(script: &[RecordBuf], region: GenomeRegion) -> Vec<String> {
        script
            .iter()
            .filter(|record| {
                let (Some(first), Some(last)) = (record.alignment_start(), record.alignment_end())
                else {
                    return false;
                };
                record.reference_sequence_id() == Some(region.contig.0 as usize)
                    && usize::from(first) as u64 <= region.end.get()
                    && usize::from(last) as u64 >= region.start.get()
            })
            .map(|record| String::from_utf8_lossy(record.name().expect("named")).into_owned())
            .collect()
    }

    /// The script every walk test below runs on: five 30-base reads starting every 15 bases
    /// across contig 0, which is 100 long. They overlap each other two deep — the shape a
    /// real region walk meets — and the last ends at 90, inside the contig, because a read
    /// running off the end is a reference fetch failure rather than an interesting case here.
    fn script() -> Vec<RecordBuf> {
        (0..5)
            .map(|i| read_at(&format!("r{i}"), 1 + i * 15))
            .collect()
    }

    // -----------------------------------------------------------------
    // The walk (B1)
    // -----------------------------------------------------------------

    /// **The oracle, one region at a time.** Every region's answer is compared against a
    /// linear scan of the same script, which is the whole reason the in-memory arm exists:
    /// a cursor that agrees with a hand scan on a list nobody indexed is a cursor whose
    /// *bookkeeping* is right, before any file can hide a defect in it.
    #[test]
    fn a_region_returns_exactly_what_a_linear_scan_of_the_same_script_returns() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        for region in [
            region(1, 200),
            region(1, 10),
            region(15, 25),
            region(41, 50),
            region(150, 200),
            region(1, 1),
        ] {
            assert_eq!(
                reads_of(&mut cursor, region),
                by_linear_scan(&script, region),
                "region {}..={}",
                region.start.get(),
                region.end.get(),
            );
        }
    }

    /// **The case a single-query test cannot reach, and the one this whole milestone exists
    /// for.** A run of regions through *one* cursor must give what the same regions give
    /// through fresh ones. The previous attempt at this feature passed 1,471 unit tests while
    /// losing 3,830 loci precisely because every read-path test drove a single query.
    #[test]
    fn a_run_of_regions_through_one_cursor_matches_the_same_regions_one_at_a_time() {
        let script = script();
        let regions = [
            region(1, 20),
            region(21, 40),
            region(41, 60),
            region(61, 80),
            region(81, 100),
        ];

        let mut shared = cursor_over(script.clone());
        for region in regions {
            let through_one = reads_of(&mut shared, region);
            let through_a_fresh_one = reads_of(&mut cursor_over(script.clone()), region);
            assert_eq!(
                through_one,
                by_linear_scan(&script, region),
                "region {}..={} through a shared cursor",
                region.start.get(),
                region.end.get(),
            );
            assert_eq!(through_one, through_a_fresh_one);
        }
    }

    /// Regions **backwards**, which the design allows and which is where the first attempt
    /// lost reads: its rule assumed moving forward along the chromosome meant moving forward
    /// through the file, and the index makes that false.
    #[test]
    fn regions_walked_backwards_return_what_a_linear_scan_returns() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        for region in [region(81, 100), region(41, 60), region(1, 20)] {
            assert_eq!(
                reads_of(&mut cursor, region),
                by_linear_scan(&script, region),
                "region {}..={} walked backwards",
                region.start.get(),
                region.end.get(),
            );
        }
    }

    /// The same region twice must answer the same both times — a cursor whose state drifted
    /// with use would show up here and nowhere else.
    #[test]
    fn asking_for_the_same_region_twice_answers_the_same_twice() {
        let script = script();
        let mut cursor = cursor_over(script.clone());
        let asked = region(21, 60);

        let first = reads_of(&mut cursor, asked);
        let second = reads_of(&mut cursor, asked);
        assert_eq!(first, by_linear_scan(&script, asked));
        assert_eq!(first, second);
    }

    /// **Spec §7's abandoned region.** A caller that stops pulling and moves elsewhere leaves
    /// the cursor mid-region — next to where the first attempt lost reads — and the next
    /// region must still be answered in full.
    #[test]
    fn a_region_abandoned_half_way_does_not_disturb_the_next_one() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        let _first = cursor.next_read().expect("a read").expect("it decodes");
        // …and the caller walks away without draining.

        let next = region(41, 60);
        assert_eq!(reads_of(&mut cursor, next), by_linear_scan(&script, next));
    }

    /// A region no read touches is an empty answer, not an error — and it must not leave the
    /// cursor unable to answer the next one.
    #[test]
    fn a_region_with_no_reads_is_empty_and_harmless() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        assert!(reads_of(&mut cursor, region(500, 600)).is_empty());
        let after = region(1, 20);
        assert_eq!(reads_of(&mut cursor, after), by_linear_scan(&script, after));
    }

    /// **A cursor whose file has failed must say so, not answer short.**
    ///
    /// After a fatal read the filter is finished for good — but the cursor is still holding
    /// reads, so it can keep producing *plausible* answers for every later region: not empty,
    /// which might be noticed, but **truncated**, which would not be. The failure was
    /// reported once, when it happened; a later region has to be told the run is over.
    ///
    /// Found in review by driving a failing read and then asking for a region past it: the
    /// cursor answered `Ok` with `[]` where a linear scan had two reads.
    #[test]
    fn a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short() {
        let mut script = script();
        // A record that decodes but whose *filtering* cannot be completed: its footprint runs
        // off the end of the contig, so the reference fetch fails, which is fatal.
        script.push(read_named_with_length(
            "overruns",
            0,
            95,
            READ_LENGTH as usize,
        ));
        let mut cursor = cursor_over(script);

        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        let mut saw_failure = false;
        while let Some(read) = cursor.next_read() {
            if read.is_err() {
                saw_failure = true;
                break;
            }
        }
        assert!(saw_failure, "the fixture must actually reach a fatal read");

        let after = cursor.move_to_region(region(1, 100));
        assert!(
            matches!(after, Err(CursorError::AfterFailure { .. })),
            "a cursor that met a fatal error answered a later region instead of refusing: \
             {after:?}",
        );
    }

    /// A cursor that has not been pointed anywhere yields nothing rather than guessing at a
    /// region — which would make the first region's reads depend on call order.
    #[test]
    fn a_cursor_pointed_nowhere_yields_nothing() {
        let mut cursor = cursor_over(script());

        assert!(cursor.next_read().is_none());
        assert_eq!(cursor.kept_reads(), 0);
    }

    /// **The refusal, and the promise that comes with it** (spec §10, owner 2026-08-02): a
    /// region on another chromosome is refused — never silently ignored, never answered from
    /// this chromosome's reads — and the cursor is left able to serve its own.
    #[test]
    fn a_region_on_another_chromosome_is_refused_and_the_cursor_survives() {
        let script = script();
        let mut cursor = cursor_over(script.clone());
        let served = region(1, 40);
        let before = reads_of(&mut cursor, served);

        let foreign = GenomeRegion {
            contig: ContigId(1),
            start: Position(1),
            end: Position(40),
        };
        let error = cursor
            .move_to_region(foreign)
            .expect_err("a region on contig 1 is not this cursor's business");
        assert!(matches!(
            error,
            CursorError::WrongChromosome {
                cursor_contig: ContigId(0),
                requested_contig: ContigId(1),
                ..
            }
        ));

        assert_eq!(cursor.contig(), ContigId(0));
        assert_eq!(
            reads_of(&mut cursor, served),
            before,
            "the cursor stopped serving its own chromosome after refusing another's",
        );
    }

    /// **The obligation the refusal rests on**: the check runs before anything is touched, so
    /// a refused region cannot disturb the region already in progress. Pinned by refusing
    /// *mid-walk* and then continuing the walk.
    #[test]
    fn a_refused_region_does_not_disturb_the_walk_in_progress() {
        let script = script();
        let mut cursor = cursor_over(script.clone());
        let served = region(1, 100);
        let whole = by_linear_scan(&script, served);

        cursor.move_to_region(served).expect("on this chromosome");
        let first = cursor.next_read().expect("a read").expect("it decodes");
        let held = cursor.kept_reads();

        let foreign = GenomeRegion {
            contig: ContigId(1),
            start: Position(1),
            end: Position(40),
        };
        assert!(cursor.move_to_region(foreign).is_err());
        assert_eq!(
            cursor.kept_reads(),
            held,
            "the refusal touched the cursor's state — the check must run before anything moves",
        );

        // The walk carries on from where it was, as though the refusal had not happened.
        let mut names = vec![String::from_utf8_lossy(&first.qname).into_owned()];
        while let Some(read) = cursor.next_read() {
            names.push(String::from_utf8_lossy(&read.expect("it decodes").qname).into_owned());
        }
        assert_eq!(names, whole);
    }

    /// Reads are **kept** as a region is walked, and a region beginning at or after the last
    /// one keeps them.
    #[test]
    fn reads_are_kept_while_a_region_is_walked_and_survive_a_forward_move() {
        let mut cursor = cursor_over(script());

        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            0,
            "nothing is read before it is asked for"
        );
        while cursor.next_read().is_some() {}
        assert_eq!(cursor.kept_reads(), 5, "every read of the region is held");

        // Forward, and no read has ended before base 1, so everything is still reachable.
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            5,
            "a region beginning at or after the last one reuses what is held",
        );
    }

    /// **A backward region drops everything**, because the rule's argument does not hold for
    /// it: a record it needs may begin before the last region's start, and the reader has
    /// already gone past. That is the case the first attempt got wrong — it assumed moving
    /// forward along the chromosome meant moving forward through the file, and a bin index
    /// makes that false.
    #[test]
    fn a_backward_region_drops_everything_and_jumps() {
        let mut cursor = cursor_over(script());

        let _ = reads_of(&mut cursor, region(46, 100));
        assert!(cursor.kept_reads() > 0);

        cursor
            .move_to_region(region(1, 30))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            0,
            "a region beginning before the last one served cannot reuse what is held",
        );
    }

    /// **Eviction, which is the rule's mirror image**: a read that ends before this region
    /// begins cannot touch it, or any later region, because every later region begins at or
    /// after this one.
    #[test]
    fn a_forward_region_evicts_only_the_reads_that_ended_before_it() {
        let mut cursor = cursor_over(script());

        // Reads start at 1, 16, 31, 46, 61 and are 30 long, so they end at 30, 45, 60, 75, 90.
        let _ = reads_of(&mut cursor, region(1, 100));
        assert_eq!(cursor.kept_reads(), 5);

        cursor
            .move_to_region(region(61, 100))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            2,
            "the three reads ending at 30, 45 and 60 are all before 61",
        );

        cursor
            .move_to_region(region(90, 100))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            1,
            "only the read ending at 90 survives"
        );

        cursor
            .move_to_region(region(91, 100))
            .expect("on this chromosome");
        assert_eq!(cursor.kept_reads(), 0, "and then nothing does");
    }

    /// **The oracle the plan asks for, over the five shapes it names** (B2): ascending,
    /// backward, overlapping, adjacent and far-apart regions, driven through *one* cursor and
    /// compared against a linear scan of the same script at every step.
    ///
    /// This is the test the first attempt at this feature did not have. It passed 1,471 unit
    /// tests while losing 3,830 of 236,081 loci, because every read-path test drove a single
    /// query.
    #[test]
    fn a_run_of_every_region_shape_matches_a_linear_scan_at_every_step() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        let regions = [
            // ascending and adjacent
            region(1, 20),
            region(21, 40),
            region(41, 60),
            // overlapping, forward
            region(50, 80),
            region(55, 85),
            // far apart, forward
            region(95, 100),
            // backward
            region(1, 30),
            // backward again, further
            region(1, 5),
            // forward from there, overlapping
            region(3, 40),
            // exactly the same region twice
            region(3, 40),
            // a region no read touches, then back into the reads
            region(97, 100),
            region(16, 45),
        ];

        for region in regions {
            assert_eq!(
                reads_of(&mut cursor, region),
                by_linear_scan(&script, region),
                "region {}..={} in a run through one cursor",
                region.start.get(),
                region.end.get(),
            );
        }
    }

    /// **The record the early stop consumed must survive into the next region.** The stop
    /// fires *on* a record — the first beginning past the region's end — and that record has
    /// already been taken from the reader. A region that continues from here, rather than
    /// jumping, has nothing else that will ever produce it.
    #[test]
    fn the_read_the_early_stop_consumed_reaches_the_next_region() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        // Ends at 35, so the stop fires on the read starting at 46.
        let first = region(1, 35);
        assert_eq!(reads_of(&mut cursor, first), by_linear_scan(&script, first));

        // Forward, so this reuses — and the read at 46 is the one the stop took.
        let second = region(36, 70);
        assert_eq!(
            reads_of(&mut cursor, second),
            by_linear_scan(&script, second),
            "the read the previous region's early stop consumed was lost",
        );
    }

    /// A read owed to two consecutive regions is given to **both** — kept reads are examined
    /// afresh by every region, not drained by the first that sees them.
    #[test]
    fn a_read_touching_two_regions_is_given_to_both() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        // The read starting at 31 ends at 60, so it touches both of these.
        let left = region(31, 45);
        let right = region(46, 60);
        assert!(by_linear_scan(&script, left).contains(&"r2".to_string()));
        assert!(by_linear_scan(&script, right).contains(&"r2".to_string()));

        assert_eq!(reads_of(&mut cursor, left), by_linear_scan(&script, left));
        assert_eq!(reads_of(&mut cursor, right), by_linear_scan(&script, right));
    }

    /// Reuse must not hand the same read to one region twice — which is what simply *not
    /// clearing* the kept set would do, because the reader carries on and the filter arm does
    /// not check what is already held.
    #[test]
    fn no_read_is_given_to_one_region_twice() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        for region in [region(1, 40), region(20, 70), region(30, 100)] {
            let names = reads_of(&mut cursor, region);
            let mut unique = names.clone();
            unique.sort();
            unique.dedup();
            assert_eq!(
                names.len(),
                unique.len(),
                "a read was yielded twice in one region",
            );
        }
    }

    /// **Random scripts, random runs of regions, one cursor.** The hand-written tables above
    /// name the shapes a reader can think of; this looks for the ones nobody thought of.
    ///
    /// It earned its place immediately during review of the previous step: an equivalent
    /// property test killed an `end >= region.start` boundary mutation that a table of six
    /// regions did not reach.
    ///
    /// Every region's answer is compared against a linear scan of the same script, so a
    /// cursor and the answer it is checked against cannot share a mistake — the scan does not
    /// call any of the code under test.
    #[test]
    fn any_run_of_regions_through_one_cursor_matches_a_linear_scan() {
        use proptest::prelude::*;

        // Reads on a 100-base contig, long enough to clear the default minimum length.
        let read_starts = prop::collection::vec(1usize..=70, 0..8);
        let regions = prop::collection::vec(
            (1u64..=100, 0u64..=40).prop_map(|(start, width)| (start, (start + width).min(100))),
            1..10,
        );

        proptest!(|(starts in read_starts, asked in regions)| {
            // A coordinate-sorted file, which is what every layer below assumes.
            let mut starts = starts;
            starts.sort_unstable();
            let script: Vec<RecordBuf> = starts
                .iter()
                .enumerate()
                .map(|(i, start)| read_at(&format!("r{i}"), *start))
                .collect();

            let mut cursor = cursor_over(script.clone());
            for (start, end) in asked {
                let asked = region(start, end);
                prop_assert_eq!(
                    reads_of(&mut cursor, asked),
                    by_linear_scan(&script, asked),
                    "region {}..={}",
                    start,
                    end
                );
            }
        });
    }

    /// Nothing is unpacked ahead of demand: pulling one read must not walk the script.
    #[test]
    fn nothing_is_read_ahead_of_demand() {
        let mut cursor = cursor_over(script());
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        assert!(cursor.next_read().is_some());
        assert_eq!(
            cursor.kept_reads(),
            1,
            "one read asked for, one read decoded",
        );
    }

    // -----------------------------------------------------------------
    // What the cursor did, as opposed to what it returned (B3)
    // -----------------------------------------------------------------

    /// **The claim the whole feature rests on: a forward walk decodes each read once.**
    ///
    /// Every region overlaps the last, which is what a real walk looks like — the caller
    /// widens each region by 5,000 bases against regions averaging 390 — so almost every read
    /// is owed to several regions. Decoding is what must not repeat.
    #[test]
    fn a_forward_walk_decodes_each_read_once() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        // Overlapping, ascending: 1..40, 20..60, 40..80, 60..100.
        let mut served = 0;
        for (start, end) in [(1, 40), (20, 60), (40, 80), (60, 100)] {
            served += reads_of(&mut cursor, region(start, end)).len();
        }

        let counts = cursor.counts();
        assert_eq!(
            counts.reads_decoded,
            script.len() as u64,
            "a forward walk decoded {} reads out of a script of {}",
            counts.reads_decoded,
            script.len(),
        );
        assert_eq!(
            counts.regions_jumping, 1,
            "only the first region has nothing to reuse",
        );
        assert_eq!(counts.regions_reusing, 3);
        assert!(
            counts.reads_replayed > 0,
            "reads owed to several regions must be served from what is held",
        );
        assert_eq!(
            counts.reads_decoded + counts.reads_replayed,
            served as u64,
            "every read handed to a caller was either decoded or replayed, and not both",
        );
    }

    /// The counters against the same walk with the rule *not* applying: regions that go
    /// backwards cannot reuse, so each one decodes again. The contrast is the measurement.
    #[test]
    fn a_backward_walk_cannot_reuse_and_decodes_again() {
        let mut forward = cursor_over(script());
        for (start, end) in [(1, 40), (20, 60), (40, 80)] {
            let _ = reads_of(&mut forward, region(start, end));
        }

        let mut backward = cursor_over(script());
        for (start, end) in [(40, 80), (20, 60), (1, 40)] {
            let _ = reads_of(&mut backward, region(start, end));
        }

        assert_eq!(backward.counts().regions_reusing, 0);
        assert_eq!(backward.counts().reads_replayed, 0);
        assert!(
            backward.counts().reads_decoded > forward.counts().reads_decoded,
            "the backward walk decoded {} against the forward walk's {}",
            backward.counts().reads_decoded,
            forward.counts().reads_decoded,
        );
    }

    /// **A region whose answer is entirely held must reach the file not at all** — the
    /// saving stated as the thing it is, rather than inferred from a faster run.
    ///
    /// (This is *not* what pins the kept walk's early `return None`. That short-circuit has
    /// no observable effect — see its comment — and no test can be written that it would
    /// fail; this one passes with or without it.)
    #[test]
    fn a_region_answered_from_what_is_held_reads_nothing() {
        let mut cursor = cursor_over(script());

        // Read the whole contig, so everything is held.
        let _ = reads_of(&mut cursor, region(1, 100));
        let after_the_first_pass = cursor.counts().reads_decoded;
        assert_eq!(after_the_first_pass, 5);

        // A narrower region inside it: every read it needs is held, and the reads beyond it
        // are held too — so the walk can answer without touching the file.
        let _ = reads_of(&mut cursor, region(16, 40));
        assert_eq!(
            cursor.counts().reads_decoded,
            after_the_first_pass,
            "a region whose answer was entirely held still read from the file",
        );
    }

    /// Eviction is counted, so "the kept set is bounded" is an observation rather than an
    /// argument — and the two ways a read leaves the set both show up.
    #[test]
    fn every_read_that_leaves_the_kept_set_is_counted() {
        let mut cursor = cursor_over(script());

        let _ = reads_of(&mut cursor, region(1, 100));
        assert_eq!(cursor.kept_reads(), 5);
        assert_eq!(cursor.counts().reads_evicted, 0);

        // Forward past three of them: evicted one at a time.
        cursor
            .move_to_region(region(61, 100))
            .expect("on this chromosome");
        assert_eq!(cursor.counts().reads_evicted, 3);
        assert_eq!(cursor.kept_reads(), 2);

        // Backward: the rest go together, and they are counted the same way.
        cursor
            .move_to_region(region(1, 10))
            .expect("on this chromosome");
        assert_eq!(
            cursor.counts().reads_evicted,
            5,
            "a jump drops what it holds, and dropping is dropping however it happens",
        );
        assert_eq!(cursor.kept_reads(), 0);
    }

    /// **The tally survives across regions with no field and no hand-over** (arch §2.3), which
    /// is what a filter living as long as the cursor buys: `read_group_counts` is a
    /// whole-chromosome total rather than whichever region happened to be last.
    ///
    /// The per-query sources this replaces had to fold their counts back into the file as each
    /// stream ended, or the drops they recorded vanished with the stream.
    #[test]
    fn the_step_one_tally_accumulates_across_regions() {
        let mut cursor = cursor_over(script());

        let _ = reads_of(&mut cursor, region(1, 40));
        let after_one = cursor.read_group_counts();
        assert!(
            !after_one.is_empty(),
            "the walk met at least one read group"
        );
        let kept_after_one: u64 = after_one.iter().map(|(_, counts)| counts.kept).sum();
        assert!(kept_after_one > 0);

        // A backward region, so nothing is replayed and every read is filtered again.
        let _ = reads_of(&mut cursor, region(1, 100));
        let kept_after_two: u64 = cursor
            .read_group_counts()
            .iter()
            .map(|(_, counts)| counts.kept)
            .sum();
        assert!(
            kept_after_two > kept_after_one,
            "the tally did not accumulate: {kept_after_one} then {kept_after_two}",
        );
    }

    /// A cursor that has done nothing says so, rather than reporting a number nobody produced.
    #[test]
    fn a_fresh_cursor_has_counted_nothing() {
        let cursor = cursor_over(script());

        assert_eq!(cursor.counts(), CursorCounts::default());
    }

    /// The message names both contigs, names them apart, **and names the file** — two bare
    /// integers are only meaningful against a particular contig table, and a run holds up to
    /// 320 cursors over many files.
    #[test]
    fn the_wrong_chromosome_message_names_the_file_and_both_contigs() {
        let error = CursorError::WrongChromosome {
            path: Arc::from(Path::new("/data/sample.bam")),
            cursor_contig: ContigId(20),
            requested_contig: ContigId(7),
        };

        assert_eq!(
            error.to_string(),
            "cursor on '/data/sample.bam' covers contig 20 but the region is on contig 7",
        );
    }

    /// The path is rendered, not debug-printed: `Path` has no `Display`, so the naive
    /// `{path}` does not compile and the naive `{path:?}` prints quotes and escapes into an
    /// operator-facing message.
    #[test]
    fn the_read_failure_message_renders_the_path_and_keeps_the_cause() {
        let error = CursorError::ReadRecord {
            path: Arc::from(Path::new("/data/sample.bam")),
            source: std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated block"),
        };

        assert_eq!(
            error.to_string(),
            "reading alignment file '/data/sample.bam' failed",
        );
        // The cause survives as a `#[source]`, so the renderer that walks the chain reaches
        // it — without it, "reading … failed" would be the whole story.
        let source = std::error::Error::source(&error).expect("the io error is the source");
        assert_eq!(source.to_string(), "truncated block");
    }
}
