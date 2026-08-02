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
//! **What is deliberately absent is the forget rule.** This cursor keeps reads, and throws
//! all of them away on every `move_to_region`: correct, and no faster than the code it will
//! replace. Deciding which kept reads *may be reused* is the next step, on its own, because
//! it is the one part of this design that can lose reads without anything failing — a rule
//! that drops a read it should have kept produces a wrong genotype, not a crash. Building the
//! machinery first means the rule arrives as a diff against something already known to be
//! correct (spec §6, and the plan's principles).

#![allow(
    dead_code,
    reason = "Milestone B builds the cursor against a scripted list; the callers that hold one \
              arrive with the BAM arm at Milestone C and the generators at D. Until then every \
              item here is reached from its own tests. Remove at C3, where the cursor is wired \
              to a real file — if it is still needed then, nothing is reading through the \
              cursor that was built for it."
)]

use std::collections::VecDeque;
use std::path::Path;
use std::sync::Arc;

use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::filtering::{ReadFilter, ReadFilterConfig, ReadFilterError};
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::record_reader::RecordReader;
use crate::ng::read::input::region_records::RegionRecords;
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
    /// The region being served, or `None` before the first [`move_to_region`](Self::move_to_region).
    region: Option<GenomeRegion>,
    contig: ContigId,
    /// Named in errors, so two bare contig numbers mean something in a run holding hundreds
    /// of cursors over many files.
    path: Arc<Path>,
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
            region: None,
            contig,
            path,
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

        // **Everything is dropped, every time — for now.** Deciding which of these may be
        // reused is the next step, deliberately on its own. Until it lands this cursor is
        // correct and no faster than a fresh query, which is the right order for the one
        // rule in this design that can lose reads without anything failing.
        self.kept.clear();
        self.region = Some(region);

        // Undo the *clean* stop the previous region ended with. Without this the first region
        // a cursor drains silences it for the whole chromosome: a region boundary reaches the
        // filter as an ordinary end of input, because that is the only end a `RecordSource`
        // can report. A filter stopped by a **fatal error** is not restarted, and must not be.
        self.filter.restart_after_end_of_input();

        self.filter
            .source_mut()
            .move_to(region)
            .map_err(|source| CursorError::ReadRecord {
                path: Arc::clone(&self.path),
                source,
            })
    }

    /// The next read of the current region, or `None` at the end of it.
    ///
    /// Shaped as [`Iterator::next`] so `Iterator` stays available to build on later. A cursor
    /// that has not been pointed anywhere yields `None` rather than guessing at a region.
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>> {
        // Not pointed anywhere yields nothing rather than guessing at a region. The layer
        // below would answer the same, but a type is better for stating its own contract than
        // for leaning on a neighbour's.
        let _region = self.region?;

        // Everything the filter yields already overlaps the region — `RegionRecords` narrowed
        // it below — so what is kept here is exactly what a later region may be able to
        // reuse. **Nothing reads from `kept` yet**: until the forget rule lets kept reads
        // survive a reposition there is never anything in it to replay, so the walk that
        // would do the replaying lands with the rule that makes it reachable (B2). Writing it
        // here would have been code no test could execute — which is what a first draft of
        // this step did, and a reviewer found it by putting a `panic!` in the loop body and
        // watching all thirteen tests pass.
        match self.filter.next() {
            None => None,
            Some(Err(error)) => Some(Err(self.read_failure(error))),
            Some(Ok(read)) => {
                self.kept.push_back(read.clone());
                Some(Ok(read))
            }
        }
    }

    /// How many reads this cursor is holding. The memory bound made observable, and what the
    /// eviction rule will be judged on.
    pub(crate) fn kept_reads(&self) -> usize {
        self.kept.len()
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

    /// Reads are **kept**, which is the mechanism the next step turns into a saving. Until
    /// the forget rule lands they are all dropped at every reposition, and that is stated
    /// here so the step that changes it has something to change.
    #[test]
    fn reads_are_kept_while_a_region_is_walked_and_dropped_when_it_moves() {
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

        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        assert_eq!(
            cursor.kept_reads(),
            0,
            "B1 keeps nothing across a reposition; B2 is what changes this",
        );
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
