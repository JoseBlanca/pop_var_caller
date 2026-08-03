//! **The alignment cursor — a reader that stays where it is.** So far only its errors:
//! [`CursorError`]. The cursor itself lands in Milestone B, and the aligned-reads readers beneath it
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
//! The cursor, over a scripted list of records — [`AlignedReadsReader::InMemory`]. The BAM arm
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
//!   `start..=start` and [`RegionRawAlignedReads`] **accepts** it. A second, hand-written test
//!   above the filter that treats a zero span as "touches nothing" would **reject** the same read —
//!   so it would be yielded when read fresh and dropped when replayed, which is a read lost
//!   with nothing failing. Whatever tests a kept read for overlap must apply the rule
//!   `RegionRawAlignedReads` applies, not a second one that looks equivalent.
//! - **Reuse needs more than a comparison.** `RegionRawAlignedReads::move_to` repositions the
//!   reader unconditionally, and the filter arm below does not check `kept` before pushing. So
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
use std::io;
use std::path::Path;
use std::sync::Arc;

use crate::ng::read::aligned_read::{AlignedRead, NoodlesRawAlignedRead, RawAlignedRead};
use crate::ng::read::filtering::{
    FilterVerdict, ReadFilterConfig, ReadFilterCounts, verdict_on_aligned_read, verdict_on_raw_read,
};
use crate::ng::read::input::aligned_reads_reader::AlignedReadsReader;
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::region_raw_aligned_reads::{
    RegionRawAlignedReads, RegionReadError, read_end, read_overlaps,
};
use crate::ng::ref_seq::{EvictableRefSeq, RawRefSeq, RefSeqError};
use crate::ng::types::{ContigId, GenomeRegion, ReadGroupId};

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

    /// The region is not a region: it ends before it begins, or it starts at base 0.
    ///
    /// **Refused rather than answered, and refused in one place for both formats.** Neither
    /// shape is dangerous on its own — an inverted region overlaps nothing, and a `0` start
    /// behaves like `1` — but "answered emptily" and "refused" are different things to a
    /// caller, and before this the answer depended on which arm and which path it reached. The
    /// per-region query rejected both in its planners; when Milestone F deleted those, the
    /// only survivor was a copy inside the **BAM** reader, so a BAM refused on a jump, a CRAM
    /// never checked, and a *forward* region reached neither because it is served without
    /// repositioning at all.
    ///
    /// So the check moved here, to the one entry point above both arms and both paths. It is
    /// **not** the reader's to make: an aligned-reads reader positions and never bounds, so
    /// `region.end` does not reach it and an inverted region is a perfectly well-defined
    /// position for it. `end` is only meaningful where the overlap test and the early stop
    /// are, which is above the reader.
    ///
    /// Coordinates are 1-based and inclusive, so `start == end` is an ordinary one-base region
    /// and not this error.
    #[error(
        "cursor on '{}' was asked for an invalid region: contig {} [{}, {}]",
        path.display(),
        region.contig.get(),
        region.start.get(),
        region.end.get()
    )]
    InvalidRegion {
        path: Arc<Path>,
        region: GenomeRegion,
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

/// A fatal, run-level failure of step 1. **Three conditions, one per piece of the work**: the
/// read off the file, the conversion, and the second filter's reference fetch.
///
/// It is yielded **in the cursor's item stream** — a fatal condition makes
/// [`AlignmentCursor::next_read`](crate::ng::read::input::cursor::AlignmentCursor::next_read) return
/// `Some(Err(..))` once and then `None` — so a caller cannot mistake it for a clean end of
/// input: `let read = read?;` propagates it. It is never folded into a per-read drop or a
/// silent end of input.
///
/// **Raised by the cursor, not here.** This module states the keep-or-drop rules; the loop that
/// meets these failures lives in `read/input/cursor.rs` (spec §5).
#[derive(Debug, thiserror::Error)]
pub(crate) enum ReadFilterError {
    /// The reader failed to hand over the next record — a truncated file, a bad block.
    #[error("reading the next alignment record failed")]
    Source(#[source] io::Error),
    /// A record's read group could not be resolved against its file's `@RG` table — an absent
    /// tag in a file declaring several groups, or a tag naming a group the file does not.
    ///
    /// **Its own variant since 2026-08-03, and it used to be `Source`'s.** Both failures leave
    /// `RegionRawAlignedReads::read_next`, so while that returned an `io::Result` the cursor
    /// could not tell them apart and charged both here to *"reading the next alignment record
    /// failed"*. An operator meeting that message goes looking for a truncated file, when what
    /// is wrong is the `@RG` header — a different fault, in a different file, wanting a
    /// different fix.
    #[error("resolving a record's read group failed")]
    ReadGroup(#[source] io::Error),
    /// A record that cleared the first filter failed to convert.
    ///
    /// **No input can reach this, and that is measured rather than assumed** (C1, 2026-08-03;
    /// confirmed independently by three reviewers). The conversion refuses exactly three things
    /// — a record with no reference sequence id, one with no alignment start, and a buffer with
    /// no read group stamped — and [`RegionRawAlignedReads::read_next`] guarantees all three
    /// before it yields: it drops anything not on this contig, `overlaps` is false without both
    /// an alignment start and an end, and the read group is resolved and stamped on the record
    /// actually handed over. An earlier version of this doc named "the unmapped flag clear yet
    /// no position" as the cause, which is one of the shapes the layer below discards first.
    ///
    /// So the variant is **defence in depth against the narrowing regressing**, not a response
    /// to any input — and it is untestable through the chain that makes it unreachable. Its two
    /// remaining constructions went with the test doubles at C2/C3, so rewriting this arm as a
    /// silent `continue` now survives the whole suite. Kept, and recorded, because "this cannot
    /// happen" is a claim that rots: if the narrowing ever stops guaranteeing one of the three,
    /// this note is what says a test became possible.
    #[error("decoding an alignment record failed")]
    Decode(#[source] io::Error),
    /// Filter #8's reference fetch failed — corrupt input, or a read reaching past a contig's
    /// end.
    ///
    /// **Narrow since B1.** `AlignmentFile::cursor` compares the two contig tables *and* proves
    /// the accessor can serve this cursor's own contig before building anything, so what is
    /// left to fail here is a read whose footprint runs off the end of the contig.
    #[error("reference access failed during filtering")]
    Reference(#[source] RefSeqError),
}

/// One read group's step-1 tally, keyed by the read group it belongs to.
///
/// `None` keys the records whose reader never stamped a read group — which
/// [`RawAlignedRead::decode`] refuses, so they are counted apart rather than charged to an
/// arbitrary group.
pub type ReadGroupCounts = (Option<ReadGroupId>, ReadFilterCounts);

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
    /// The chain below, owned: [`RegionRawAlignedReads`] holds the [`AlignedReadsReader`].
    reads: RegionRawAlignedReads,
    /// The single raw read buffer reused across the whole walk, so the pass allocates one
    /// record rather than one per read.
    buffer: NoodlesRawAlignedRead,
    /// Raw reference bases for the second filter's mismatch check.
    ///
    /// **Held here because only that filter needs it and it is the cursor that keeps it alive
    /// for the chromosome** (spec §9 Q1). Taken once rather than rebuilt per query.
    reference: R,
    /// Reused scratch the second filter's reference fetch reads into — touched only when
    /// filter #8 runs, so a walk with it disabled allocates nothing here.
    ref_buf: Vec<u8>,
    config: ReadFilterConfig,
    /// One tally per read group met, in first-seen order.
    ///
    /// A `Vec` scanned linearly rather than a map: a file declares a handful of read groups, and
    /// the scan is a few integer compares on a path that has just decoded a record. Per read
    /// group rather than per *file* because a drop rate is a read group's property — one bad
    /// library shows up as one read group with an anomalous MAPQ or mismatch rate, and a
    /// per-file tally over a file holding several would average that away (spec §7).
    ///
    /// **Named `tally` rather than `counts`, and not for variety.** The cursor already has a
    /// `counts` field — [`CursorCounts`], what the cursor *did* — and a public
    /// [`read_group_counts`](Self::read_group_counts) method whose value differs from this
    /// field's: the method stamps the `other_sample` rider onto the first entry. A field and a
    /// method spelled the same and returning different things, both reachable in one function,
    /// is a reading a caller has to know to distrust.
    read_group_tally: Vec<ReadGroupCounts>,
    /// What the layer below had already skipped as another sample's when the current tally
    /// window opened — `0` until [`reset_read_group_counts`](Self::reset_read_group_counts) is
    /// first called.
    ///
    /// **A window that reset the drops but not this would not start empty**, which is the one
    /// thing `reset_read_group_counts` promises. The skipped-record count lives on the narrowing
    /// and is cumulative for the life of the cursor — the narrowing has no notion of a window —
    /// so the window is expressed here, as the offset to subtract.
    other_sample_at_window_start: u64,
    /// **Whether this cursor has met a fatal error.** Replaces the three-way `FilterState` the
    /// filter used to keep: that state existed only because a separate filter could not tell why
    /// its source stopped, and reached "end of input" at the end of *every* region. The cursor
    /// is the thing that *causes* region ends, so it never has to ask — which is what collapses
    /// three states to one flag (spec §5).
    ///
    /// **Set by both routes into a stopped cursor**, not only by a failed read: a failed
    /// reposition leaves the reader at an unknown position, so no later region can be served
    /// from it either.
    failed: bool,
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
    ///
    /// **Infallible, and that is what changed here.** It used to return `Result` for one
    /// reason: it built a `ReadFilter` through a constructor that fetched a zero-length window
    /// on every contig of the header to prove each one resolved in `reference`. Nothing else
    /// about constructing a cursor can fail. The caller now proves the same thing — and more —
    /// by comparing the two contig tables before it ever gets here
    /// ([`AlignmentFile::cursor`](crate::ng::read::input::AlignmentFile::cursor)), and since
    /// C2 there is no separate filter to build at all.
    ///
    /// **The precondition is the caller's, and it is not a formality.**
    /// `with_validated_contigs` assumes the file's `@SQ` list has been proved *equal* to the
    /// accessor's table, order included. A permuted list would resolve on every fetch and make
    /// filter #8 compare each read against the wrong contig's bases, silently — which is why
    /// the check upstream is an equality and not a resolvability test.
    pub(crate) fn over_records(
        reader: AlignedReadsReader,
        contig: ContigId,
        resolution: ReadGroupResolution,
        reference: R,
        config: ReadFilterConfig,
        path: Arc<Path>,
    ) -> Self {
        Self {
            reads: RegionRawAlignedReads::new(reader, contig, resolution),
            buffer: NoodlesRawAlignedRead::default(),
            reference,
            ref_buf: Vec::new(),
            config,
            read_group_tally: Vec::new(),
            other_sample_at_window_start: 0,
            failed: false,
            kept: VecDeque::new(),
            examined: 0,
            last_region_start: None,
            last_emitted: None,
            region: None,
            contig,
            path,
            counts: CursorCounts::default(),
        }
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
    /// **Every refusal happens before any state moves**, so a rejected region leaves the
    /// cursor exactly as it was — which is what makes "unharmed and still good for its own"
    /// true by construction rather than by care (spec §10). The chromosome is checked first,
    /// then the region's shape, then whether this cursor is still alive.
    pub fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), CursorError> {
        if region.contig != self.contig {
            return Err(CursorError::WrongChromosome {
                path: Arc::clone(&self.path),
                cursor_contig: self.contig,
                requested_contig: region.contig,
            });
        }

        // **The one place a malformed region is refused, for both formats and both paths.**
        // See `CursorError::InvalidRegion`: below here the arms diverge — a reader positions
        // at `region.start` and never sees `region.end` at all — so this is the last point at
        // which the two shapes mean anything, and the only one at which a *forward* region is
        // looked at before being served from what is held.
        if region.is_empty() || region.start.get() == 0 {
            return Err(CursorError::InvalidRegion {
                path: Arc::clone(&self.path),
                region,
            });
        }

        // A cursor whose file has failed serves nothing further. Checked here, after the
        // chromosome test and before anything moves, so a dead cursor says so rather than
        // answering later regions out of whatever it is still holding.
        if self.failed {
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

        // **The reposition comes first, because it is the only thing here that can fail.**
        // Everything below it — the eviction, the counters, the region state — is committed
        // only once the reader is really where this region needs it. See the block after it.
        if !reuse && let Err(source) = self.reads.jump_to(region) {
            self.failed = true;
            return Err(CursorError::ReadRecord {
                path: Arc::clone(&self.path),
                source,
            });
        }

        if reuse {
            self.reads.continue_into(region);
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

        // **Everything from the reposition down is committed only because the jump succeeded**,
        // and that ordering is the point: a cursor left pointing at a region it never reached
        // would serve that region from wherever the failed seek abandoned the reader — not
        // empty, which might be noticed, but *wrong*, which would not be. A failed jump now
        // leaves the kept set, the counters and the region state exactly as they were.
        //
        // A failed reposition also **stops the cursor for good**: the reader's position is
        // unknown afterwards, so no later region can be served from it, the reuse path included.
        // Found by mutation at C1 — swallowing the failure entirely left all 2,845 tests green
        // — and widened to the `failed` flag on the owner's ruling.
        //
        // **The two halves are independent, and only one of them used to be pinned.** The flag
        // masks the ordering, so reverting the ordering left the whole suite green until
        // `a_failed_reposition_leaves_the_cursor_untouched` drove a seek that fails *after* a
        // region has been served and checked what survived.
        //
        // Every kept read is offered to the new region, including ones the last region was
        // already given: consecutive regions overlap, and a read touching both is owed to
        // both.
        self.examined = 0;
        self.region = Some(region);
        self.last_region_start = Some(region.start.get());
        // Per region, because a new region rewinds through what is kept: positions ascend
        // *within* a region, never across one.
        self.last_emitted = None;

        // **Nothing restarts anything here any more**, and its absence is the point. A filter
        // that lived apart from the cursor reached "end of input" at the end of *every* region
        // — a region boundary was the only end its source could report — so a cursor had to undo
        // that on each move or the first region silenced it for the whole chromosome. The cursor
        // is the thing that causes region ends, so it never has to ask: the loop below simply
        // reads on, and `failed` is the only stop it keeps (spec §5).
        Ok(())
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

        // Then read on, from wherever the reader is. Everything below already overlaps the
        // region — `RegionRawAlignedReads` narrowed it — so what is kept here is exactly what a
        // later region may be able to reuse.
        match self.next_filtered_read() {
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

    /// **Step 1, one read at a time** — the loop that used to be `ReadFilter::next`.
    ///
    /// Reads the next raw aligned read into the one reused buffer, rejects it on its flag and
    /// mapping quality, converts only what survives, rejects *that* on its length, CIGAR and
    /// mismatch fraction, and charges every drop to the read group it came from. Returns the
    /// first read that passes.
    ///
    /// **The order is the whole design** (spec §2): the six flag/MAPQ filters read values the
    /// raw read already carries, so a read they drop never pays for a conversion; the other
    /// three read fields only the conversion produces. Moving either filter across the
    /// conversion changes no output and would quietly undo that.
    ///
    /// Fused on `failed`: a fatal condition is yielded once and then this returns `None`
    /// forever, and `move_to_region` refuses every later region.
    fn next_filtered_read(&mut self) -> Option<Result<AlignedRead, ReadFilterError>> {
        if self.failed {
            return None;
        }
        loop {
            match self.reads.read_next(&mut self.buffer) {
                Ok(true) => {}
                // The region is done. **Not the end of the file**, and the cursor knows the
                // difference because it is the thing that set the region — which is why no
                // three-way state is needed here (spec §5).
                Ok(false) => return None,
                // The two faults the layer below can meet are charged apart: a file that cannot
                // be read, and a record whose read group cannot be resolved.
                Err(RegionReadError::Read(error)) => {
                    return self.fail(ReadFilterError::Source(error));
                }
                Err(RegionReadError::ReadGroup(error)) => {
                    return self.fail(ReadFilterError::ReadGroup(error));
                }
            }

            // The first filter — flag and mapping quality, before any conversion. Exhaustive
            // so a new `FilterVerdict` variant cannot silently fall through to the conversion.
            match verdict_on_raw_read(self.buffer.flag(), self.buffer.mapq(), &self.config) {
                FilterVerdict::Keep => {}
                FilterVerdict::Drop(reason) => {
                    self.tally_for_buffered_read().record_drop(reason);
                    continue;
                }
            }

            // The conversion, for the survivors only.
            let read = match self.buffer.decode() {
                Ok(read) => read,
                Err(error) => return self.fail(ReadFilterError::Decode(error)),
            };

            // The second filter — length, CIGAR, mismatch fraction, in that order.
            match verdict_on_aligned_read(&read, &self.reference, &self.config, &mut self.ref_buf) {
                Ok(FilterVerdict::Keep) => {
                    self.tally_for_buffered_read().kept += 1;
                    return Some(Ok(read));
                }
                Ok(FilterVerdict::Drop(reason)) => {
                    self.tally_for_buffered_read().record_drop(reason);
                    continue;
                }
                Err(error) => return self.fail(ReadFilterError::Reference(error)),
            }
        }
    }

    /// Stop the cursor and yield the failure that stopped it. The three fatal arms share it.
    fn fail(&mut self, error: ReadFilterError) -> Option<Result<AlignedRead, ReadFilterError>> {
        self.failed = true;
        Some(Err(error))
    }

    /// The tally of the read group the **buffered** read belongs to, created on first sight.
    ///
    /// Keyed on the buffer rather than passed in, because the first filter drops a read before
    /// anything owns it — the buffer is the only thing that knows, and it knows because the
    /// region narrowing resolved the read group before handing it over.
    fn tally_for_buffered_read(&mut self) -> &mut ReadFilterCounts {
        let read_group = self.buffer.read_group();
        let index = match self
            .read_group_tally
            .iter()
            .position(|(id, _)| *id == read_group)
        {
            Some(index) => index,
            None => {
                self.read_group_tally
                    .push((read_group, ReadFilterCounts::default()));
                self.read_group_tally.len() - 1
            }
        };
        &mut self.read_group_tally[index].1
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

    /// Release the reference bases this cursor's read filter has gone past.
    ///
    /// **A cursor keeps its reference reader as long as it keeps its file**, and the filter
    /// reads that reference once per surviving read to check its mismatch fraction. Reading
    /// is all it does — the window only ever grows — so on a densely covered chromosome the
    /// reader ends up holding one byte for every base walked: about 250 MB on human
    /// chromosome 1, against a walk that otherwise peaks around 25 MB.
    ///
    /// The caller says when, because only the caller knows what it will ask for next. Passing
    /// too high a position costs a re-read and never an answer — eviction is a hint.
    pub fn evict_reference_before(&self, pos: u64)
    where
        R: EvictableRefSeq,
    {
        self.reference.evict_before(pos);
    }

    /// How many reference bases this cursor's reader is holding — the bound made observable.
    pub fn resident_reference_bases(&self) -> usize
    where
        R: EvictableRefSeq,
    {
        self.reference.resident_bases()
    }

    /// Step-1's per-read-group tally, for as much of the chromosome as this cursor has read.
    ///
    /// **One entry per read group met, in first-seen order, and never summed here.** A drop rate
    /// is a read group's property: one bad library is one read group with an anomalous rate, and
    /// adding them up erases exactly that. A caller wanting a total adds them itself.
    ///
    /// **Cumulative** — a whole-cursor total rather than whichever region happened to be last.
    /// The per-query sources this replaces had to fold their counts back into the file as each
    /// stream ended, or the drops they recorded vanished with the stream.
    ///
    /// The skipped-record count belongs to no single read group — a foreign record is never
    /// yielded to the filters, and its own group is not one of these — so it is reported once,
    /// against the first entry.
    pub fn read_group_counts(&self) -> Vec<ReadGroupCounts> {
        let mut counts = self.read_group_tally.clone();
        match counts.first_mut() {
            Some((_, first)) => first.other_sample = self.other_sample_this_window(),
            // The tally is empty — either every record was another sample's, or the window has
            // just been reset. **The `None` key here is a fabrication**, and it is spelled the
            // same as the genuine one `ReadGroupCounts` documents: a read whose reader never
            // stamped a group. A caller cannot tell the two apart, which is a wart worth knowing
            // about rather than one worth a second key nobody asked for.
            None => counts.push((
                None,
                ReadFilterCounts {
                    other_sample: self.other_sample_this_window(),
                    ..ReadFilterCounts::default()
                },
            )),
        }
        counts
    }

    /// Records skipped as another sample's **since the current tally window opened**.
    ///
    /// The narrowing counts them for the life of the cursor and has no notion of a window, so
    /// the subtraction happens here.
    ///
    /// **The subtraction cannot underflow**, and an earlier version of this comment justified
    /// the saturation with the wrong direction. The underlying count only ever grows — every
    /// contribution is a `+=` — and the baseline is sampled *from* that same number, so the
    /// baseline can never exceed it. `saturating_sub` is free belt-and-braces, not a guard
    /// against a reachable case.
    ///
    /// **The hazard that is real runs the other way, and it is the CRAM arm's**, pre-dating this
    /// window: `CramAlignedReadsReader` adds a container's foreign-record count each time it
    /// decodes one, so a container decoded twice — which a backward reposition can cause —
    /// contributes twice. The field's own doc already calls the number container-granular and
    /// says it can run ahead of where a walk has reached. A window therefore *starts* at zero
    /// honestly and may *over*-report afterwards on CRAM, exactly as the unwindowed number does.
    /// Recorded rather than fixed: it is the reader's accounting, not the window's.
    fn other_sample_this_window(&self) -> u64 {
        self.reads
            .other_sample_records()
            .saturating_sub(self.other_sample_at_window_start)
    }

    /// Start a fresh tally window: [`read_group_counts`](Self::read_group_counts) reports only
    /// what happens from here.
    ///
    /// **The caller chooses the window, and the cursor never chooses one for itself** — not per
    /// region, and not on a reposition (spec §7). Regions overlap by about 93 % and a read is
    /// filtered once, when it is first read off the file, so a per-region tally would record
    /// where the reader happened to be when it met a bad read rather than how bad the region is;
    /// the numbers would not sum to the chromosome's total and would not be comparable between
    /// regions.
    ///
    /// **It resets the tally and nothing else.** The reads being held, the walk's own counters
    /// ([`counts`](Self::counts)), the region being served and whether the cursor has failed are
    /// all untouched — this is a question about a number, not a way to rewind a cursor.
    ///
    /// Named for the tally it resets rather than `reset_counts`, because the cursor keeps two
    /// unrelated tallies: this one, and [`CursorCounts`] — what the cursor *did*.
    pub fn reset_read_group_counts(&mut self) {
        self.read_group_tally.clear();
        // Without this the window would start with every foreign record the cursor has ever
        // stepped over already in it.
        self.other_sample_at_window_start = self.reads.other_sample_records();
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
    use crate::ng::read::input::aligned_reads_reader::InMemoryAlignedReadsReader;
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, bam_header, fixture_read_group, fixture_reference_bases, matching_contigs,
        only_tally, read_named_with_length,
    };
    use crate::ng::ref_seq::{InMemoryRefSeq, RefSeqError};
    use crate::ng::types::Position;
    use noodles_sam::alignment::RecordBuf;

    // For the drop-tally fixture below, which builds records by hand rather than through
    // `read_named_with_length` because each one has to fail a *named* filter.
    use crate::bam::alignment_input::{
        FLAG_DUPLICATE, FLAG_PAIRED, FLAG_QC_FAIL, FLAG_SECONDARY, FLAG_SUPPLEMENTARY,
        FLAG_UNMAPPED,
    };
    use crate::ng::read::filtering::ReadFilterCounts;
    use noodles_core::Position as RecordPosition;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    /// An all-`A` reference over the fixture contigs, so every fixture read matches perfectly
    /// and nothing is dropped for mismatching — this milestone is about *which* reads come
    /// back, not about filtering.
    fn reference_bases() -> InMemoryRefSeq {
        fixture_reference_bases()
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
        cursor_over_reader(InMemoryAlignedReadsReader::new(
            bam_header(&matching_contigs()),
            records,
        ))
    }

    /// A cursor over a scripted list whose reader **fails at read `read_index`** instead of
    /// handing back a record — the truncated file, driven from the bottom of the real chain.
    fn cursor_whose_reader_fails_at_read(
        records: Vec<RecordBuf>,
        read_index: usize,
    ) -> AlignmentCursor<InMemoryRefSeq> {
        cursor_over_reader(
            InMemoryAlignedReadsReader::new(bam_header(&matching_contigs()), records)
                .with_failure_at_read(read_index),
        )
    }

    /// A cursor on contig 0 over an already-built reader — the one place the six-argument
    /// constructor is spelled out, so a scripted variant cannot drift from the plain one.
    fn cursor_over_reader(reader: InMemoryAlignedReadsReader) -> AlignmentCursor<InMemoryRefSeq> {
        AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(reader),
            ContigId(0),
            fixture_read_group(),
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        )
    }

    /// The step-1 failure a cursor error raised **while reading** is carrying, or a panic naming
    /// what came instead.
    ///
    /// Scoped to the failures `AlignmentCursor::read_failure` builds, which is what `next_read`
    /// yields. `CursorError::ReadRecord` has a **second** construction site —
    /// `move_to_region`'s failed reposition (`cursor.rs`, the `jump_to` arm) — whose `source` is
    /// the reader's raw `io::Error` with no step-1 error inside it, so this would panic on one
    /// of those rather than answer. That is deliberate: a reposition failure is not a step-1
    /// failure and has no variant to discriminate.
    ///
    /// Asserting by type rather than on the rendered message is the point: a message substring
    /// would still pass if `Source` and `Reference` were swapped at the call site, which is
    /// exactly the wiring these tests exist to pin.
    fn step_one_failure(error: &CursorError) -> &ReadFilterError {
        match error {
            CursorError::ReadRecord { source, .. } => source
                .get_ref()
                .and_then(|cause| cause.downcast_ref::<ReadFilterError>())
                .expect("a failure yielded by next_read carries the step-1 error that caused it"),
            other => panic!("expected CursorError::ReadRecord, got {other:?}"),
        }
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

    // -----------------------------------------------------------------
    // Which fatal condition a fault is charged to (C1)
    // -----------------------------------------------------------------
    //
    // `filtering.rs` pinned these against two test doubles — `ErroringSource`, whose every read
    // failed, and `FakeSource`, a `RecordSource` implemented in the test module. Both stood
    // where `RegionRawAlignedReads` and `AlignedReadsReader` stand in a real run, so neither
    // said anything about the two layers in between.
    //
    // **That gap was already closed from the other end, on real inputs, and these tests are not
    // what closes it.** `open_bam.rs`'s `t10_a_truncated_file_fails_once_and_then_refuses_
    // later_regions` truncates an indexed BAM mid-walk, and
    // `a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short` above drives
    // a reference fetch off the end of a contig. Both carry a fault up the whole chain and both
    // stop the cursor; making the narrowing swallow a reader failure, or making a reference
    // failure a silent drop, fails them.
    //
    // What neither can see is **which** `ReadFilterError` the fault was charged to — both match
    // `Err(_)`. Swapping `Source` and `Reference` at their two call sites in `ReadFilter::next`
    // leaves the whole suite green. These two tests pin the charge, and each is the only test in
    // the tree that fails when its own variant is swapped. That is their whole job, and it is
    // why a *scripted* fault is worth the mechanism: the script chose the kind, so the test can
    // assert it. C3 deleted the doubles, and this is what had to exist first.

    /// **A failure reading off the file is fatal, charged to `Source`, and the cursor never
    /// recovers from it.**
    ///
    /// The fault is scripted into the reader at the *second* read, so the walk is working
    /// before it breaks: a chain that never delivered a read at all would satisfy "the second
    /// call is an error" without proving anything travelled through it.
    #[test]
    fn a_failure_reading_off_the_file_is_fatal_through_the_whole_chain() {
        let mut cursor = cursor_whose_reader_fails_at_read(script(), 1);
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        let first = cursor
            .next_read()
            .expect("a read")
            .expect("the first record is clean");
        assert_eq!(first.qname, b"r0", "the reads before the fault flow");

        let error = cursor
            .next_read()
            .expect("the fault is yielded, not swallowed into a clean end of input")
            .expect_err("the scripted read failure");
        assert!(
            matches!(step_one_failure(&error), ReadFilterError::Source(_)),
            "a read failure must be charged to the source, got {error:?}",
        );

        // Yielded once, and then the cursor is finished: a fused walk, and every later region
        // refused rather than answered short out of what is still held.
        assert!(cursor.next_read().is_none(), "the walk stays stopped");
        assert!(matches!(
            cursor.move_to_region(region(1, 100)),
            Err(CursorError::AfterFailure { .. })
        ));
    }

    /// **An unresolvable read group is fatal, and charged to `ReadGroup` — not to `Source`.**
    ///
    /// The fourth fatal condition, and until 2026-08-03 it wore the third one's name: both
    /// failures leave `RegionRawAlignedReads::read_next`, so while that returned an `io::Result`
    /// the cursor could not tell them apart and rendered this one as *"reading the next alignment
    /// record failed"*. An operator meeting that goes looking for a truncated file, when what is
    /// wrong is the `@RG` header.
    ///
    /// The fixture is a file declaring two read groups and a record carrying neither.
    #[test]
    fn an_unresolvable_read_group_is_fatal_and_charged_to_its_own_condition() {
        use crate::ng::read::input::read_groups::RecordOwner;
        use crate::ng::types::ReadGroupId;

        let resolution = ReadGroupResolution::PerRecord(
            vec![
                (
                    "rg1".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(0)),
                ),
                (
                    "rg2".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(1)),
                ),
            ]
            .into_iter()
            .collect(),
        );
        let untagged = record_with_seq(
            "untagged",
            10,
            60,
            Flags::from(FLAG_PAIRED),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        let mut cursor = AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                vec![untagged],
            )),
            ContigId(1),
            resolution,
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        );
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");

        let error = cursor
            .next_read()
            .expect("the unresolvable read group is yielded, not swallowed")
            .expect_err("a record with no RG tag in a two-group file cannot be attributed");
        assert!(
            matches!(step_one_failure(&error), ReadFilterError::ReadGroup(_)),
            "an unresolvable read group must not be charged to the file failing to read: \
             {error:?}",
        );

        assert!(cursor.next_read().is_none(), "the walk stays stopped");
        assert!(matches!(
            cursor.move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(50),
            }),
            Err(CursorError::AfterFailure { .. })
        ));
    }

    /// **A reference fetch that fails mid-walk is fatal too, and charged to `Reference`.**
    ///
    /// Filter #8 fetches the bases a read covers; a read whose footprint runs off the end of
    /// the contig has none to fetch. Under the fatal error model that is corrupt input to stop
    /// on, not a read to quietly drop — a validly-aligned read cannot cover positions the
    /// contig does not have. Needs no scripted fault: the record is ordinary and it is the
    /// *reference* that cannot answer.
    ///
    /// **Distinct from `a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_
    /// short` in exactly one assertion**, over the same `overruns` fixture — and that assertion
    /// is the reason it exists. The older test matches `Err(_)`, so charging this failure to
    /// `Source` instead of `Reference` passes it. Deleting either as "the duplicate" loses one
    /// of the two properties.
    #[test]
    fn a_reference_fetch_failure_mid_walk_is_fatal_through_the_whole_chain() {
        // Contig 0 is 100 bases; a 30-base read at 95 reaches 124.
        let mut cursor = cursor_over(vec![
            read_at("clean", 1),
            read_named_with_length("overruns", 0, 95, READ_LENGTH as usize),
        ]);
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        let first = cursor
            .next_read()
            .expect("a read")
            .expect("the first record is clean");
        assert_eq!(first.qname, b"clean");

        let error = cursor
            .next_read()
            .expect("the fetch failure is yielded")
            .expect_err("the read runs off the end of the contig");
        assert!(
            matches!(
                step_one_failure(&error),
                ReadFilterError::Reference(RefSeqError::OutOfBounds { .. }),
            ),
            "an out-of-bounds fetch must be charged to the reference, got {error:?}",
        );

        assert!(cursor.next_read().is_none(), "the walk stays stopped");
        assert!(matches!(
            cursor.move_to_region(region(1, 100)),
            Err(CursorError::AfterFailure { .. })
        ));
    }

    /// **A reposition that fails is refused, not answered** — the third fatal route, and the one
    /// no test reached.
    ///
    /// A reader can break in two places, and only one of them is a read. On a BAM,
    /// `begin_region` runs an index query, so a corrupt index fails the **move** and no read is
    /// ever attempted. Swallowing that serves the region from wherever the reader happened to be
    /// left — not empty, which might be noticed, but *wrong*, which would not be. It is the same
    /// condition `CursorError::AfterFailure` exists to make loud one layer later.
    ///
    /// Found by mutation during C1's review: replacing the whole arm with
    /// `let _ = self.filter.source_mut().jump_to(region); Ok(())` passed all 2,845 tests. Every
    /// other test that touches that line asserts the move *succeeds*.
    ///
    /// # And it stops the cursor, which is the half C2 added
    ///
    /// Writing this test at C1 showed the defect was bigger than a missing assertion:
    /// `move_to_region` set `region` and `last_region_start` **before** the fallible `jump_to`,
    /// so a failed reposition left the cursor pointed at a region it never reached, serving from
    /// wherever the seek abandoned the reader — and because `last_region_start` had moved too,
    /// the *next* forward region took the reuse path and read on without jumping at all. The
    /// `AfterFailure` guard missed it because it asked the filter, and a failed reposition never
    /// reached the filter.
    ///
    /// C2 repairs both halves on the owner's ruling: the reposition happens before any of the
    /// region's state is committed, and it sets the cursor's own `failed` flag — which covers
    /// **both** routes into a stopped cursor, not only the one that replaced `FilterState`.
    ///
    /// Latent in production: both real arms' `begin_region` are effectively infallible (the BAM
    /// arm queries an in-memory index, the CRAM arm resets state), which is why a scripted fault
    /// is what reaches it.
    #[test]
    fn a_reposition_that_fails_is_refused_and_stops_the_cursor() {
        let mut cursor = cursor_over_reader(
            InMemoryAlignedReadsReader::new(bam_header(&matching_contigs()), script())
                .with_failing_seek_at(0),
        );

        // The first region always jumps — there is nothing held to continue from.
        let moved = cursor.move_to_region(region(1, 100));
        assert!(
            matches!(moved, Err(CursorError::ReadRecord { .. })),
            "a failed reposition must be reported, not swallowed: {moved:?}",
        );
        assert!(
            cursor.next_read().is_none(),
            "a region that was never reached was served anyway",
        );

        // And no later region is answered out of a reader whose position is unknown — including
        // a *forward* one, which is the case the reuse path would otherwise serve without
        // repositioning at all.
        assert!(
            matches!(
                cursor.move_to_region(region(50, 100)),
                Err(CursorError::AfterFailure { .. })
            ),
            "a cursor whose reposition failed answered a later region",
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

    // -----------------------------------------------------------------
    // The order guard (F5) — inherited from the per-region query's T4a..T4d
    // -----------------------------------------------------------------
    //
    // The per-region query wrapped every stream in an `OrderVerified` iterator and pinned it
    // with four tests. That wrapper is gone and its job moved into `emit`, which had **no
    // test at all** — the guard was written, reviewed and shipped without one, and a guard
    // that never fires looks exactly like a guard that works. These are those four rules
    // restated against the cursor, dropping only the one that no longer exists here: the
    // contig-order half of the old key, which a cursor cannot reach because
    // `move_to_region` refuses a foreign chromosome before anything moves
    // (`a_region_on_another_chromosome_is_refused_and_the_cursor_survives`).

    /// **T4a — a read going backwards within one region is fatal**, and the message says
    /// where.
    ///
    /// The open gate proved the file *claims* `SO:coordinate`; this is the file that claims it
    /// and lies, which the header check structurally cannot see. Mutation-verified: deleting
    /// the comparison in [`AlignmentCursor::emit`] makes this test fail.
    #[test]
    fn a_read_going_backwards_within_a_region_is_a_fatal_error() {
        // Scripted out of order on purpose. The in-memory reader hands its list over as
        // given (`an_out_of_order_script_is_not_quietly_sorted`), so the fault reaches the
        // guard rather than being tidied away below it.
        let mut cursor = cursor_over(vec![read_at("a", 1), read_at("c", 61), read_at("b", 31)]);

        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        let mut names = Vec::new();
        let mut error = None;
        while let Some(item) = cursor.next_read() {
            match item {
                Ok(read) => names.push(String::from_utf8_lossy(&read.qname).into_owned()),
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }

        assert_eq!(names, vec!["a", "c"], "the reads before the break flow");
        match error.expect("the regression must be fatal") {
            CursorError::OutOfOrderRead {
                path,
                position,
                after,
            } => {
                assert_eq!(position, 31);
                assert_eq!(after, 61);
                // Asserted because `path` is a constructor argument: a wiring mistake would
                // produce a correct-looking error naming the wrong file, in a run holding
                // hundreds of cursors.
                assert_eq!(path.as_ref(), Path::new("/fixture/sample.bam"));
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
    }

    /// **T4c — equal positions are ordinary, not a fault.** Several reads may start at the
    /// same base, and a guard written with `<=` would reject every pile-up in every real file.
    #[test]
    fn reads_sharing_a_start_position_are_not_out_of_order() {
        let script = vec![
            read_at("a", 1),
            read_at("b", 1),
            read_at("c", 1),
            read_at("d", 2),
        ];
        let mut cursor = cursor_over(script.clone());

        assert_eq!(
            reads_of(&mut cursor, region(1, 100)),
            by_linear_scan(&script, region(1, 100)),
        );
    }

    /// **T4d — the guard is scoped to the region, not to the cursor.** A caller is entitled to
    /// ask for a later region and then an earlier one; the second is a new walk over what is
    /// kept, not a regression. The old query got this for free by building a fresh guard per
    /// query — a cursor has to reset `last_emitted` on every move, and forgetting to would
    /// turn every backward jump into a fatal error.
    ///
    /// The last assertion is what makes this able to fail in the other direction: a guard
    /// that had been silently disarmed would also pass the first two.
    #[test]
    fn an_earlier_region_after_a_later_one_is_not_out_of_order() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        let later = region(61, 100);
        assert_eq!(reads_of(&mut cursor, later), by_linear_scan(&script, later),);
        let earlier = region(1, 20);
        assert_eq!(
            reads_of(&mut cursor, earlier),
            by_linear_scan(&script, earlier),
            "an earlier region is a new walk, not a regression",
        );

        // And the guard is still armed after the backward move.
        let mut planted = cursor_over(vec![read_at("a", 1), read_at("c", 61), read_at("b", 31)]);
        let _ = reads_of(&mut planted, region(61, 100));
        planted
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        let mut saw_error = false;
        while let Some(item) = planted.next_read() {
            if item.is_err() {
                saw_error = true;
                break;
            }
        }
        assert!(
            saw_error,
            "the guard must still catch a regression inside the region it was reset for",
        );
    }

    /// **An inverted or zero-start region is refused, and the cursor survives it.**
    ///
    /// Both shapes were rejected by the per-region query's planners. Milestone F deleted those
    /// and left the rule with one home — a copy inside the *BAM* reader — so a BAM refused on
    /// a jump, a CRAM never checked, and a **forward** region reached neither, because a
    /// forward region is served without repositioning at all. Disabling the survivor left the
    /// whole suite green.
    ///
    /// It now lives in `move_to_region`, above both arms and both paths. Asserted through the
    /// reuse path as well as the jump path, which is the case the old placement structurally
    /// could not reach.
    #[test]
    fn a_malformed_region_is_refused_and_the_cursor_survives_it() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        // Establish a served region first, so the *next* move is a candidate for reuse rather
        // than a jump — the path the reader-level check could never see.
        let served = region(1, 40);
        assert_eq!(
            reads_of(&mut cursor, served),
            by_linear_scan(&script, served)
        );

        for malformed in [
            // Ends before it begins, and forward of the last region served, so the forget rule
            // would reuse rather than reposition.
            region(80, 70),
            // The same, backwards, so the jump path is covered too.
            region(20, 10),
            // Base 0 does not exist in 1-based inclusive coordinates.
            GenomeRegion {
                contig: ContigId(0),
                start: Position(0),
                end: Position(50),
            },
        ] {
            match cursor.move_to_region(malformed) {
                Err(CursorError::InvalidRegion { path, region }) => {
                    assert_eq!(region, malformed);
                    assert_eq!(path.as_ref(), Path::new("/fixture/sample.bam"));
                }
                other => panic!("expected InvalidRegion for {malformed:?}, got {other:?}"),
            }
        }

        // **Nothing moved.** The cursor is still serving the region it was given, and a fresh
        // valid region is answered in full — a refusal that had disturbed the walk would show
        // up here and nowhere else.
        let next = region(41, 100);
        assert_eq!(reads_of(&mut cursor, next), by_linear_scan(&script, next));
    }

    /// `start == end` is an ordinary one-base region, not a malformed one — coordinates are
    /// 1-based and inclusive, and a check written with `<=` would reject every single-base
    /// region a caller asks for.
    #[test]
    fn a_one_base_region_is_not_malformed() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        let one_base = region(20, 20);
        assert_eq!(
            reads_of(&mut cursor, one_base),
            by_linear_scan(&script, one_base),
        );
    }

    // -----------------------------------------------------------------
    // The step-1 drop tally, by hand count — moved here from `filtering.rs`
    // -----------------------------------------------------------------
    //
    // `filtering.rs` asserted this twice, over a real BAM and a real CRAM, through the
    // whole-file `BamRecordSource`/`CramRecordSource` it owned. Those sources are gone: a
    // filter module has no business opening files. What the two tests actually pinned is the
    // *filter's accounting*, and it is pinned here, where the chain that accounting belongs to
    // composes — `AlignedReadsReader` → `RegionRawAlignedReads` → `ReadFilter` →
    // `AlignmentCursor`.
    //
    // The fixture sits on **contig 1**, which the fixture table makes 200 bases long; its
    // records reach base 149 and would not fit on contig 0's 100.

    fn record_with_seq(qname: &str, start: usize, mapq: u8, flags: Flags, seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname)
            .set_reference_sequence_id(1usize)
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq in range"))
            .set_alignment_start(RecordPosition::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, seq.len())].into_iter().collect())
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .build()
    }

    /// A 30-base read on contig 1 carrying an `RG` tag, so a `PerRecord` resolution has something
    /// to resolve. The tally tests need this; `Sole` short-circuits before the tag is read.
    fn tagged_read(qname: &str, start: usize, flags: u16, read_group: &str) -> RecordBuf {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut record = record_with_seq(
            qname,
            start,
            60,
            Flags::from(flags),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        record.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(read_group.as_bytes().to_vec().into()),
        );
        record
    }

    fn boundary_deletion_record(start: usize, seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name("bad")
            .set_reference_sequence_id(1usize)
            .set_flags(Flags::from(FLAG_PAIRED))
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_alignment_start(RecordPosition::try_from(start).unwrap())
            .set_cigar(
                [Op::new(Kind::Deletion, 2), Op::new(Kind::Match, seq.len())]
                    .into_iter()
                    .collect(),
            )
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .build()
    }

    /// The fixture: two kept reads plus one read per *mapped* drop reason (#1–#4,
    /// #6–#9), all on a single all-`A` contig. The #5 unmapped drop is covered
    /// separately (`read_filter_charges_an_unmapped_read_end_to_end`): a realistic
    /// unmapped read has MAPQ 0 and is charged to #2 first, and a fake MAPQ does
    /// not survive a CRAM round-trip. Returns the records and the
    /// `ReadFilterCounts` a correct pass must produce; the counts are asserted, not
    /// read order.
    fn drop_fixture() -> (Vec<RecordBuf>, ReadFilterCounts) {
        let clean = |name: &str, start: usize| {
            record_with_seq(
                name,
                start,
                60,
                Flags::from(FLAG_PAIRED),
                b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            )
        };
        let records = vec![
            clean("kept1", 10), // kept
            clean("kept2", 20), // kept
            record_with_seq(
                "dup",
                30,
                60,
                Flags::from(FLAG_PAIRED | FLAG_DUPLICATE),
                b"AAAAAAAAAA",
            ), // #1
            record_with_seq("lowmapq", 40, 5, Flags::from(FLAG_PAIRED), b"AAAAAAAAAA"), // #2 (mapq 5 < 20)
            record_with_seq(
                "supp",
                50,
                60,
                Flags::from(FLAG_PAIRED | FLAG_SUPPLEMENTARY),
                b"AAAAAAAAAA",
            ), // #3
            record_with_seq(
                "sec",
                60,
                60,
                Flags::from(FLAG_PAIRED | FLAG_SECONDARY),
                b"AAAAAAAAAA",
            ), // #4
            record_with_seq(
                "qcfail",
                70,
                60,
                Flags::from(FLAG_PAIRED | FLAG_QC_FAIL),
                b"AAAAAAAAAA",
            ), // #6
            record_with_seq("tooshort", 80, 60, Flags::from(FLAG_PAIRED), b"AAAAA"), // #7 (len 5 < 30)
            // #8: 5 non-reference bases out of 30 = 16.7% > 10%.
            record_with_seq(
                "highmismatch",
                90,
                60,
                Flags::from(FLAG_PAIRED),
                b"CCCCCAAAAAAAAAAAAAAAAAAAAAAAAA",
            ),
            // #9 bad CIGAR — AND high-mismatch (5 `C`s). Because it fails both #9
            // and #8, the exact counts below discriminate the ng #9-before-#8
            // order: it must land in `bad_cigar` (not `high_mismatch_fraction`).
            boundary_deletion_record(120, b"CCCCCAAAAAAAAAAAAAAAAAAAAAAAAA"), // #9
        ];
        let expected = ReadFilterCounts {
            kept: 2,
            duplicate: 1,
            low_mapq: 1,
            supplementary: 1,
            secondary: 1,
            unmapped: 0,
            qc_fail: 1,
            too_short: 1,
            high_mismatch_fraction: 1,
            bad_cigar: 1,
            other_sample: 0,
        };
        (records, expected)
    }

    /// A cursor over a scripted list on **contig 1**, whose 200 bases the drop fixture needs.
    /// `reference_bases()` is all-`A` over both fixture contigs, which is what makes the
    /// mismatch counts below the fixture's own property rather than the reference's.
    fn cursor_over_contig_one(records: Vec<RecordBuf>) -> AlignmentCursor<InMemoryRefSeq> {
        AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                records,
            )),
            ContigId(1),
            fixture_read_group(),
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        )
    }

    /// **Every drop reason, hand-counted, through the whole chain.**
    ///
    /// Two kept reads and one read per mapped drop reason. The assertion is on the *counts*,
    /// never on read order, so it survives any reordering of the walk that does not change
    /// what is dropped and why.
    ///
    /// **The last record is the load-bearing one.** It has a leading deletion (#9, bad CIGAR)
    /// *and* five mismatching bases out of thirty (#8, above the 10 % ceiling), so it fails
    /// both — and the exact tally below is what pins ng's **#9-before-#8** order. Charge it to
    /// `high_mismatch_fraction` instead and this fails while every "the right reads survive"
    /// test stays green.
    #[test]
    fn a_walk_charges_every_drop_reason_by_hand_count() {
        let (records, expected) = drop_fixture();
        let mut cursor = cursor_over_contig_one(records);

        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");

        let mut kept = Vec::new();
        while let Some(read) = cursor.next_read() {
            kept.push(read.expect("the fixture decodes"));
        }

        assert_eq!(kept.len(), 2, "exactly the two clean reads survive");
        assert_eq!(only_tally(&cursor.read_group_counts()), expected);
    }

    // -----------------------------------------------------------------
    // The tally's remaining properties, re-homed from `filtering.rs` (C2)
    // -----------------------------------------------------------------
    //
    // Ten tests died with `ReadFilter`, and C2's first account of where their properties went
    // was wrong on two rows. The review corrected both by experiment; this is the version that
    // survived it.
    //
    // - **Four have no successor, by design.** They drove `source_mut`,
    //   `restart_after_end_of_input` and the three-way `FilterState` — the machinery spec §5
    //   collapses to one flag, because the cursor causes region ends and never has to ask why
    //   reading stopped.
    // - **Two covered the fatal stop**, and the cursor's own fatal-path tests cover it — but
    //   only for `Source` and `Reference`. One of them drove the **`Decode`** arm, which no
    //   test now reaches: C1 established no input can, so the arm is unreachable rather than
    //   untested, which is a better disposition than the one first written but not the same one.
    // - **One was the tally surviving a reposition**, and it was *not* covered by
    //   `the_step_one_tally_accumulates_across_regions` as first claimed: that test's second
    //   region begins at or after its first, so it reuses and never jumps. Its successor is
    //   `the_tally_survives_a_reposition_that_drops_everything_and_jumps`.
    // - **Three are re-homed below.**
    //
    // The two tally tests after them are not re-homings: they close mutations that survived the
    // whole suite, which is how the first account's gaps were found.

    /// **The #5 counter the drop fixture omits**: an unmapped read that clears the mapping-quality
    /// filter is charged to `unmapped`, through the whole chain.
    ///
    /// The fixture leaves it out because a *realistic* unmapped read has MAPQ 0 and is charged to
    /// #2 first, and because a fake mapping quality does not survive a CRAM round trip. Given one
    /// anyway, the flag must decide it.
    ///
    /// It needs a placed unmapped read — a contig and a position — because one without them has
    /// no footprint and the region narrowing drops it below the filters, uncounted. That is the
    /// same measurement C1 recorded against the plan's D2.
    #[test]
    fn an_unmapped_read_that_clears_the_mapping_quality_filter_is_charged_to_unmapped() {
        let unmapped = record_with_seq(
            "unmapped",
            10,
            60,
            Flags::from(FLAG_PAIRED | FLAG_UNMAPPED),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        let clean = record_with_seq(
            "clean",
            50,
            60,
            Flags::from(FLAG_PAIRED),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        let mut cursor = cursor_over_contig_one(vec![unmapped, clean]);
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");

        let mut kept = Vec::new();
        while let Some(read) = cursor.next_read() {
            kept.push(read.expect("the fixture decodes"));
        }

        assert_eq!(kept.len(), 1, "only the clean read survives");
        let tally = only_tally(&cursor.read_group_counts());
        assert_eq!(tally.unmapped, 1);
        assert_eq!(tally.kept, 1);
    }

    /// A cursor over an empty script yields nothing and counts nothing — and still answers with a
    /// tally rather than an empty vector, because a caller folding several cursors together needs
    /// an entry to fold.
    #[test]
    fn a_cursor_over_an_empty_script_yields_nothing_and_counts_nothing() {
        let mut cursor = cursor_over(Vec::new());
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        assert!(cursor.next_read().is_none());
        assert_eq!(
            only_tally(&cursor.read_group_counts()),
            ReadFilterCounts::default()
        );
    }

    /// **The tally is running, not final**: it is readable part-way through a walk and already
    /// reflects the drops the walk has met.
    ///
    /// Stated because the opposite is just as plausible a design — a tally folded together only
    /// when the walk ends — and because a caller reading it early would then get zeros that look
    /// like a clean file.
    #[test]
    fn the_tally_is_readable_before_the_walk_is_finished() {
        let duplicate = record_with_seq(
            "dup",
            10,
            60,
            Flags::from(FLAG_PAIRED | FLAG_DUPLICATE),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        let clean = record_with_seq(
            "clean",
            50,
            60,
            Flags::from(FLAG_PAIRED),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        );
        let mut cursor = cursor_over_contig_one(vec![duplicate, clean]);
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");

        // One read pulled — which is the *second* record, the first having been dropped.
        assert!(matches!(cursor.next_read(), Some(Ok(_))));
        let part_way = only_tally(&cursor.read_group_counts());
        assert_eq!(part_way.duplicate, 1, "the drop is already counted");
        assert_eq!(part_way.kept, 1);

        assert!(cursor.next_read().is_none(), "and then the script is done");
        assert_eq!(only_tally(&cursor.read_group_counts()).kept, 1);
    }

    /// **A failed reposition leaves the cursor untouched** — the kept set, the counters and the
    /// region state all as they were.
    ///
    /// The owner's ruling has two halves: the reposition happens before any of the new region's
    /// state is committed, *and* it sets `failed`. **Only the flag was pinned.** The flag masks
    /// the ordering — nothing can observe the stale state once the cursor refuses everything —
    /// so hoisting the commits back above the jump left all 2,855 tests green.
    ///
    /// Reaching it needs a seek that fails **after a region has been served**, which the
    /// all-or-nothing fault could not express: a reader whose *first* seek fails has served
    /// nothing, so "left exactly as it was" and "was never anywhere" are the same observation.
    /// That is why `with_failing_seek_at` is positional.
    #[test]
    fn a_failed_reposition_leaves_the_cursor_untouched() {
        let mut cursor = cursor_over_reader(
            InMemoryAlignedReadsReader::new(bam_header(&matching_contigs()), script())
                // The first region's jump succeeds; the second one's fails.
                .with_failing_seek_at(1),
        );

        let served = reads_of(&mut cursor, region(46, 100));
        assert!(!served.is_empty(), "the first region must really be served");
        let kept_before = cursor.kept_reads();
        let counts_before = cursor.counts();
        let tally_before = cursor.read_group_counts();
        assert!(kept_before > 0, "and it must leave reads held");

        // Backwards, so the forget rule jumps — and the jump fails.
        let failed = cursor.move_to_region(region(1, 20));
        assert!(
            matches!(failed, Err(CursorError::ReadRecord { .. })),
            "the failing seek must be reported: {failed:?}",
        );

        assert_eq!(
            cursor.kept_reads(),
            kept_before,
            "a failed reposition dropped the reads it was holding",
        );
        assert_eq!(
            cursor.counts(),
            counts_before,
            "a failed reposition moved the counters for a region that never happened",
        );
        assert_eq!(
            cursor.read_group_counts(),
            tally_before,
            "a failed reposition disturbed the step-1 tally",
        );
    }

    /// **The tally survives a reposition that drops everything and jumps.**
    ///
    /// The half of the deleted `repositioning_the_source_does_not_reset_the_running_tally` that
    /// `the_step_one_tally_accumulates_across_regions` does **not** reach, and C2's accounting
    /// first mis-filed as covered by it. That test's second region begins at or after its first,
    /// so it takes the *reuse* path and never repositions — instrumented, it reports
    /// `jumping=1 reusing=1 replayed=3`. Nothing in the tree drove the tally across a jump.
    ///
    /// The deleted test said in as many words why non-erasure needs its own assertion:
    /// accumulation alone is also satisfied by a tally that stopped counting. Measured — adding
    /// `self.read_group_tally.clear()` to the jump branch of the forget rule left all 2,855
    /// tests green.
    ///
    /// This is the surface the plan singles C2 out for, and **C4 is about to add
    /// `reset_counts`** — the first legitimate caller of "clear the tally" — into a file where
    /// clearing it on the wrong edge is invisible.
    #[test]
    fn the_tally_survives_a_reposition_that_drops_everything_and_jumps() {
        let mut cursor = cursor_over(script());

        let _ = reads_of(&mut cursor, region(46, 100));
        let after_forward: u64 = cursor.read_group_counts().iter().map(|(_, c)| c.kept).sum();
        assert!(after_forward > 0, "the first region kept something");
        assert_eq!(cursor.counts().regions_jumping, 1, "the first region jumps");

        // Backwards, which is the path that drops the kept set and repositions.
        let _ = reads_of(&mut cursor, region(1, 20));
        assert_eq!(
            cursor.counts().regions_jumping,
            2,
            "a backward region must jump, or this test pins nothing",
        );
        let after_jump: u64 = cursor.read_group_counts().iter().map(|(_, c)| c.kept).sum();
        assert!(
            after_jump > after_forward,
            "the reposition reset the running tally: {after_forward} then {after_jump}",
        );
    }

    /// **Two read groups are tallied apart, not summed** — the property spec §7 exists for, and
    /// the one nothing pinned.
    ///
    /// A drop rate is a read group's property: one bad library shows up as one read group with an
    /// anomalous mapping-quality or mismatch rate, and adding them together erases exactly that
    /// signal. The failure is silent in the strongest sense — it changes no output, no dump and
    /// no read, only a number nobody is looking at.
    ///
    /// Found by mutation while moving the tally onto the cursor: keying every read onto the first
    /// entry, so a file's libraries merge into one, **left all 2,853 tests green**. The existing
    /// multi-read-group test collects `(qname, read_group)` off the *reads* and never looks at the
    /// tally.
    ///
    /// Each group is given a *different* drop reason, so a fold that merged them would have to
    /// lose one of the two counters to pass.
    #[test]
    fn two_read_groups_are_tallied_apart_rather_than_summed() {
        use crate::ng::read::input::read_groups::RecordOwner;
        use crate::ng::types::ReadGroupId;

        let resolution = ReadGroupResolution::PerRecord(
            vec![
                (
                    "rg1".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(0)),
                ),
                (
                    "rg2".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(1)),
                ),
            ]
            .into_iter()
            .collect(),
        );

        let mut cursor = AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                vec![
                    // rg1 loses one read to the duplicate filter, rg2 one to QC-fail.
                    tagged_read("dup", 10, FLAG_PAIRED | FLAG_DUPLICATE, "rg1"),
                    tagged_read("kept1", 20, FLAG_PAIRED, "rg1"),
                    tagged_read("qcfail", 30, FLAG_PAIRED | FLAG_QC_FAIL, "rg2"),
                    tagged_read("kept2", 40, FLAG_PAIRED, "rg2"),
                ],
            )),
            ContigId(1),
            resolution,
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        );
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");
        while let Some(read) = cursor.next_read() {
            read.expect("the fixture decodes");
        }

        let counts = cursor.read_group_counts();
        assert_eq!(
            counts.len(),
            2,
            "the two libraries were folded into one tally: {counts:?}",
        );

        let of = |id: ReadGroupId| {
            counts
                .iter()
                .find(|(group, _)| *group == Some(id))
                .map(|(_, tally)| tally.clone())
                .unwrap_or_else(|| panic!("no tally for {id:?} in {counts:?}"))
        };
        let first = of(ReadGroupId(0));
        assert_eq!((first.duplicate, first.qc_fail, first.kept), (1, 0, 1));
        let second = of(ReadGroupId(1));
        assert_eq!((second.duplicate, second.qc_fail, second.kept), (0, 1, 1));
    }

    /// **The other-sample count rides on the first entry, and is not a drop.**
    ///
    /// A record belonging to another sample says nothing about how *this* sample's read groups
    /// behaved, so it is deliberately outside every other counter — charging it as a drop would
    /// make a shared file look like a low-quality one. It belongs to no single read group either,
    /// so it is reported once, against the **first**.
    ///
    /// **It takes two of our read groups and a foreign record to state that**, which is why no
    /// earlier fixture could: with one read group the first entry *is* the last, so moving the
    /// rider to the last entry is undetectable. Measured — that mutation survived all 2,855 tests
    /// until this fixture existed.
    #[test]
    fn the_other_sample_count_rides_on_the_first_entry_and_is_not_a_drop() {
        use crate::ng::read::input::read_groups::RecordOwner;
        use crate::ng::types::ReadGroupId;

        let resolution = ReadGroupResolution::PerRecord(
            vec![
                (
                    "mine1".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(0)),
                ),
                (
                    "mine2".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(1)),
                ),
                (
                    "theirs".to_string().into_boxed_str(),
                    RecordOwner::OtherSample,
                ),
            ]
            .into_iter()
            .collect(),
        );

        let mut cursor = AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                vec![
                    tagged_read("ours-a", 10, FLAG_PAIRED, "mine1"),
                    tagged_read("foreign", 20, FLAG_PAIRED, "theirs"),
                    tagged_read("ours-b", 30, FLAG_PAIRED, "mine2"),
                ],
            )),
            ContigId(1),
            resolution,
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        );
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(1),
                start: Position(1),
                end: Position(FIXTURE_CONTIGS[1].1 as u64),
            })
            .expect("on this chromosome");
        while let Some(read) = cursor.next_read() {
            read.expect("the fixture decodes");
        }

        let counts = cursor.read_group_counts();
        assert_eq!(
            counts.len(),
            2,
            "two of ours, and the foreign one is not a read group of this sample"
        );

        let (_, first) = &counts[0];
        assert_eq!(
            first.other_sample, 1,
            "the foreign record must be reported against the first entry: {counts:?}",
        );
        let (_, second) = &counts[1];
        assert_eq!(
            second.other_sample, 0,
            "…and against that entry only, or a caller summing them double-counts it",
        );

        // **Not a drop, and that is the point of the separate counter.** Every drop counter of
        // both entries is zero: the foreign record was skipped, not rejected.
        for (group, tally) in &counts {
            assert_eq!(
                (tally.kept, tally.duplicate, tally.low_mapq, tally.unmapped),
                (1, 0, 0, 0),
                "the foreign record was charged as a drop against {group:?}",
            );
        }
    }

    // -----------------------------------------------------------------
    // Choosing the tally window (C4)
    // -----------------------------------------------------------------

    /// **A fresh window starts empty, and nothing else on the cursor moves.**
    ///
    /// The tally is cumulative until the caller says otherwise, so the only way to scope it is
    /// to say so — and the reset must be a question about a *number*, not a way to rewind a
    /// cursor. Both halves are asserted: the tally is empty afterwards, and the reads being
    /// held, the walk's own counters, and the cursor's ability to carry on are not.
    #[test]
    fn resetting_the_tally_starts_a_fresh_window_and_moves_nothing_else() {
        let script = script();
        let mut cursor = cursor_over(script.clone());

        // A forward region first, so the *next* one can be backward and therefore re-read.
        let first = reads_of(&mut cursor, region(46, 100));
        assert!(!first.is_empty(), "the first window saw reads");
        let kept: u64 = cursor
            .read_group_counts()
            .iter()
            .map(|(_, counts)| counts.kept)
            .sum();
        assert!(kept > 0, "…and tallied them");

        let held_before = cursor.kept_reads();
        let walk_before = cursor.counts();

        cursor.reset_read_group_counts();

        // The window is empty — every counter of every entry, not merely `kept`.
        for (group, counts) in cursor.read_group_counts() {
            assert_eq!(
                counts,
                ReadFilterCounts::default(),
                "the fresh window already had {group:?}'s drops in it",
            );
        }

        // …and nothing else moved.
        assert_eq!(
            cursor.kept_reads(),
            held_before,
            "resetting the tally dropped the reads the cursor was holding",
        );
        assert_eq!(
            cursor.counts(),
            walk_before,
            "resetting the tally disturbed the walk's own counters",
        );

        // The cursor still works, and the new window fills — from a **backward** region, which
        // drops what is held and re-reads. A forward one would replay the kept reads instead,
        // and a replayed read is not filtered again, so it is not tallied again either
        // (spec §7) — which is why this cannot simply repeat the region it just served.
        let after_reset = reads_of(&mut cursor, region(1, 20));
        assert_eq!(
            after_reset,
            by_linear_scan(&script, region(1, 20)),
            "the cursor stopped serving its own chromosome",
        );
        let kept_after: u64 = cursor
            .read_group_counts()
            .iter()
            .map(|(_, counts)| counts.kept)
            .sum();
        assert!(
            kept_after > 0,
            "the fresh window never filled, so the reset broke the tally rather than scoping it",
        );
        assert!(
            kept_after <= kept,
            "the fresh window carried the first one's count forward: {kept} then {kept_after}",
        );
    }

    /// **A reset mid-region leaves the region being served untouched.**
    ///
    /// The two tests above reset at a region *boundary* and then call `reads_of`, which begins
    /// with `move_to_region` and re-establishes everything a mutation could have disturbed — so
    /// neither can see the state that governs a region **in flight**. Measured: setting
    /// `region = None` in the reset truncates the region to nothing, and setting `examined = 0`
    /// serves a read **twice**; both survived all 2,860 tests.
    ///
    /// The duplicate is the nastier one, because the order guard cannot catch it: it compares
    /// with `<`, so a read replayed at its own position sails through.
    #[test]
    fn resetting_the_tally_mid_region_leaves_the_region_being_served_untouched() {
        let script = script();

        let mut undisturbed = cursor_over(script.clone());
        let expected = reads_of(&mut undisturbed, region(1, 100));

        let mut cursor = cursor_over(script.clone());
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        let first = cursor.next_read().expect("a first read").expect("decodes");
        let mut served = vec![String::from_utf8_lossy(&first.qname).into_owned()];

        cursor.reset_read_group_counts();

        while let Some(read) = cursor.next_read() {
            let read = read.expect("the reset must not stop the region being served");
            served.push(String::from_utf8_lossy(&read.qname).into_owned());
        }
        assert_eq!(
            served, expected,
            "resetting the tally changed which reads the region being served hands back",
        );
    }

    /// **A reset mid-region leaves the order guard armed**, for the region already in flight.
    ///
    /// `last_emitted` is per region, so clearing it in the reset would disarm the guard for the
    /// rest of the region a caller is walking — and a file that lies about being sorted would
    /// then be served silently. Survived all 2,860 before this existed.
    #[test]
    fn resetting_the_tally_leaves_the_order_guard_armed_for_the_region_being_served() {
        let mut cursor = cursor_over(vec![read_at("a", 1), read_at("c", 61), read_at("b", 31)]);
        cursor
            .move_to_region(region(1, 100))
            .expect("on this chromosome");

        let mut names = Vec::new();
        for _ in 0..2 {
            let read = cursor.next_read().expect("a read").expect("decodes");
            names.push(String::from_utf8_lossy(&read.qname).into_owned());
        }
        assert_eq!(names, vec!["a", "c"], "the reads before the break flow");

        cursor.reset_read_group_counts();

        match cursor
            .next_read()
            .expect("the backwards read")
            .expect_err("resetting the tally disarmed the order guard mid-region")
        {
            CursorError::OutOfOrderRead {
                position, after, ..
            } => {
                assert_eq!(position, 31);
                assert_eq!(after, 61);
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
    }

    /// **A reset does not make the next forward region reposition.**
    ///
    /// `last_region_start` *is* the forget rule — one number. Clearing it in the reset would make
    /// every later region jump, dropping and re-reading what it was holding. Invisible in the
    /// reads and in the acceptance dumps, which is why it survived all 2,860: only the cursor's
    /// own counters say whether the rule is still working.
    #[test]
    fn resetting_the_tally_does_not_make_the_next_forward_region_reposition() {
        let script = script();
        let mut cursor = cursor_over(script.clone());
        let _ = reads_of(&mut cursor, region(1, 50));
        let before = cursor.counts();

        cursor.reset_read_group_counts();

        let _ = reads_of(&mut cursor, region(20, 80));
        let after = cursor.counts();
        assert_eq!(
            after.regions_jumping, before.regions_jumping,
            "resetting the tally forgot where the last region began, so the next forward region \
             dropped everything it was holding and re-read it",
        );
        assert_eq!(after.regions_reusing, before.regions_reusing + 1);
    }

    /// **The window is applied on the arm a real walk takes**, not only when the tally is empty.
    ///
    /// `read_group_counts` folds the other-sample rider in two places — onto the first entry when
    /// the walk met read groups of its own, and into a fabricated entry when it met none. The
    /// test below exercises the *empty* arm, because its fixture consumes the whole contig before
    /// the reset. The `first_mut` arm — the one every real cohort walk takes — was reverted to
    /// its pre-C4 body and **all 2,860 tests stayed green**: the single line the window field
    /// exists to protect was the untested one.
    #[test]
    fn read_group_counts_scopes_the_other_sample_rider_when_the_new_window_met_its_own_reads() {
        use crate::ng::read::input::read_groups::RecordOwner;
        use crate::ng::types::ReadGroupId;

        let resolution = ReadGroupResolution::PerRecord(
            vec![
                (
                    "mine".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(0)),
                ),
                (
                    "theirs".to_string().into_boxed_str(),
                    RecordOwner::OtherSample,
                ),
            ]
            .into_iter()
            .collect(),
        );
        let mut cursor = AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                vec![
                    tagged_read("mine_early", 10, FLAG_PAIRED, "mine"),
                    tagged_read("theirs_early", 20, FLAG_PAIRED, "theirs"),
                    tagged_read("mine_late", 100, FLAG_PAIRED, "mine"),
                    tagged_read("theirs_late", 110, FLAG_PAIRED, "theirs"),
                ],
            )),
            ContigId(1),
            resolution,
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        );
        let at = |start: u64, end: u64| GenomeRegion {
            contig: ContigId(1),
            start: Position(start),
            end: Position(end),
        };

        cursor.move_to_region(at(100, 200)).expect("on contig 1");
        while let Some(read) = cursor.next_read() {
            read.expect("the fixture decodes");
        }
        assert_eq!(
            only_tally(&cursor.read_group_counts()).other_sample,
            1,
            "the first window must meet a foreign record, or this test pins nothing",
        );

        cursor.reset_read_group_counts();

        // A *backward* region, so it is re-read rather than replayed — and it meets one of our
        // own reads, so the tally is non-empty and the rider lands on the `first_mut` arm.
        cursor.move_to_region(at(1, 60)).expect("on contig 1");
        while let Some(read) = cursor.next_read() {
            read.expect("the fixture decodes");
        }
        let tally = cursor.read_group_counts();
        assert!(
            tally.first().expect("an entry").0.is_some(),
            "the second window must meet one of our own reads, or it tests the empty-tally arm",
        );
        assert_eq!(
            only_tally(&tally).other_sample,
            1,
            "the non-empty-tally arm reported the cursor's whole history, not this window's",
        );
    }

    /// **A window served entirely from replayed reads tallies nothing**, and that is worth its
    /// own test rather than a sentence in someone else's comment.
    ///
    /// A read is filtered once, when first read off the file — never again when replayed (spec
    /// §7). So a window can serve five reads and count none of them, and a caller reading these
    /// numbers as "what this window served" is reading them wrong.
    ///
    /// Without this, a later "fix" that tallied replayed reads would double-count every read the
    /// cursor keeps, and only the acceptance dumps' *counts* would move — the reads themselves
    /// would not.
    #[test]
    fn read_group_counts_stays_empty_when_a_window_is_served_only_from_replayed_reads() {
        let script = script();
        let mut cursor = cursor_over(script.clone());
        let first = reads_of(&mut cursor, region(1, 100));
        assert!(!first.is_empty(), "the first walk served reads");

        cursor.reset_read_group_counts();
        let replayed_before = cursor.counts().reads_replayed;

        let again = reads_of(&mut cursor, region(1, 100));
        assert_eq!(again, first, "the same region serves the same reads");
        assert_eq!(
            cursor.counts().reads_replayed,
            replayed_before + again.len() as u64,
            "every read of the second walk came from what was held",
        );

        assert_eq!(
            cursor.read_group_counts(),
            vec![(None, ReadFilterCounts::default())],
            "a replayed read was tallied a second time, or the empty window grew an entry",
        );
    }

    /// **The other-sample rider is part of the window too**, which is the half a `clear()` on the
    /// tally alone would miss.
    ///
    /// That count lives on the layer below and is cumulative for the life of the cursor — the
    /// narrowing has no notion of a window — so a reset that only emptied the tally would open a
    /// window with every foreign record the cursor had ever stepped over already in it.
    #[test]
    fn resetting_the_tally_also_starts_the_other_sample_count_from_zero() {
        use crate::ng::read::input::read_groups::RecordOwner;
        use crate::ng::types::ReadGroupId;

        let resolution = ReadGroupResolution::PerRecord(
            vec![
                (
                    "mine".to_string().into_boxed_str(),
                    RecordOwner::Mine(ReadGroupId(0)),
                ),
                (
                    "theirs".to_string().into_boxed_str(),
                    RecordOwner::OtherSample,
                ),
            ]
            .into_iter()
            .collect(),
        );
        let mut cursor = AlignmentCursor::over_records(
            AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
                bam_header(&matching_contigs()),
                vec![
                    tagged_read("ours", 10, FLAG_PAIRED, "mine"),
                    tagged_read("foreign", 20, FLAG_PAIRED, "theirs"),
                ],
            )),
            ContigId(1),
            resolution,
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new("/fixture/sample.bam")),
        );
        let whole = GenomeRegion {
            contig: ContigId(1),
            start: Position(1),
            end: Position(FIXTURE_CONTIGS[1].1 as u64),
        };
        cursor.move_to_region(whole).expect("on this chromosome");
        while let Some(read) = cursor.next_read() {
            read.expect("the fixture decodes");
        }
        assert_eq!(
            only_tally(&cursor.read_group_counts()).other_sample,
            1,
            "the fixture must actually meet a foreign record, or this test pins nothing",
        );

        cursor.reset_read_group_counts();

        assert_eq!(
            only_tally(&cursor.read_group_counts()).other_sample,
            0,
            "the fresh window opened with the cursor's whole other-sample history in it",
        );
    }

    /// The message names the file and both positions — two bare offsets are meaningless
    /// against a run reading hundreds of files at once.
    #[test]
    fn the_out_of_order_message_names_the_file_and_both_positions() {
        let error = CursorError::OutOfOrderRead {
            path: Arc::from(Path::new("/data/sample.bam")),
            position: 150,
            after: 200,
        };

        assert_eq!(
            error.to_string(),
            "alignment file '/data/sample.bam' yielded a read at position 150 after one at \
             200, within one region: the file is not coordinate-sorted",
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

    /// The colliding read's name is what lets a user grep the inputs and confirm the
    /// diagnosis, and the cursor indices have to read as *files* rather than as another
    /// coordinate pair sitting beside `contig 0 position 99`.
    ///
    /// Moved from `IngestError::DuplicateReadAcrossFiles` at Milestone F, which the
    /// per-region merge raised and which is now the same condition under one name.
    #[test]
    fn the_duplicate_read_message_names_the_read_and_both_files() {
        let error = CursorError::DuplicateReadAcrossFiles {
            qname: b"read-1".to_vec(),
            contig: ContigId(0),
            position: 99,
            first_file: 0,
            second_file: 1,
        };

        assert_eq!(
            error.to_string(),
            "the read 'read-1' appears in two of this sample's files (cursor 0 and cursor 1) \
             at contig 0 position 99",
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
