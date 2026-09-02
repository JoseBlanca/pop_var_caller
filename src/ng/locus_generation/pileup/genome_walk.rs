//! The walker loop. Single-threaded; drives all the building
//! blocks (active set, chain-id allocator, open-record table) through
//! the closure rule and yields finalised `PileupRecord`s through
//! a pull-shaped `Iterator`.
//!
//! The walker is a state machine: each call to `Iterator::next()`
//! advances it until at least one record is ready (or end-of-input
//! is reached), then yields one record at a time. A single walker
//! tick may emit 0, 1, or many records; the iterator buffers them
//! in a small `VecDeque` and drains across successive `next()`
//! calls.
//!
//! **No longer a verbatim copy — released from `copy_fidelity.rs` at A0 (plan 3),
//! and it has diverged since.** Copied from `src/pileup/walker/driver.rs`, then:
//!
//! - **A0:** the reference accessor is bound by ng's [`RefSeq`] rather than
//!   production's `MultiChromRefFetcher`, and the field and parameters carrying
//!   it are named `reference` rather than `ref_fetcher` — there is no fetcher
//!   any more.
//! - **C2 (plan 3):** `stop_after`, the region walk's right bound.
//! - **C3 (plan 3):** `adopting_chain_ids` / `into_chain_ids`, because ng lends
//!   the allocator rather than owning it.
//! - **D1 (the alignment cursor):** [`LookAhead`] in place of `Peekable`,
//!   [`RegionReadSource`], [`PileupWalker::move_to_region`] and
//!   [`WalkerState::begin_region`] — a walker now lives for a chromosome and is
//!   pointed at each of its regions in turn, where production builds one per
//!   chromosome and ng used to build one per region.
//!
//! The per-position walk itself — admit, process, expire, close, advance — is
//! still the transcription.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};

use crate::pileup_record::ChainId;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::GenomeRegion;

use super::active_read_set::ActiveReads;
use super::chain_id_allocator::{ChainIdAllocator, ChainIdAllocatorCounters};
use super::decompose::ReadEvent;
use super::errors::WalkerError;
use super::fast_column::FastColumnScratch;
use super::open_record::{
    OpenPileupRecord, OpenPileupRecordTable, ReadContribution, process_position,
};
use super::read_sampling;
use super::{PreparedRead, ReadLengthError, WalkerConfig};
use crate::ng::ref_seq::RefSeq;

/// Construct a [`PileupWalker`] over a coordinate-sorted stream of
/// prepared reads. The walker is an `Iterator<Item = Result<PileupRecord,
/// WalkerError>>`; callers drive it by repeatedly calling `next()`
/// (or by collecting / for-looping). After iteration ends, the
/// run's cumulative counters are available via
/// [`PileupWalker::summary`].
///
/// Coordinate-order invariant: every read pulled from `reads` must
/// have `(chrom_id, alignment_start)` non-decreasing relative to
/// the previous one. A regression is a hard error — stale or
/// malformed input shouldn't pass silently.
pub fn run<R, F>(reads: R, reference: F, config: &WalkerConfig) -> PileupWalker<R::IntoIter, F>
where
    R: IntoIterator<Item = PreparedRead>,
    F: RefSeq,
{
    PileupWalker::new(reads.into_iter(), reference, config)
}

/// A read source a walk can be pointed at one region after another — **ng's, D1.**
///
/// The walker's ordinary source is a plain `Iterator<Item = PreparedRead>` and stays one:
/// production's walk, the stage-1 differential and every unit test hand it a list, and none
/// of them has a region to be pointed at. This is the extra thing a *long-lived* walker
/// needs, so it is a separate trait bounding one method rather than a bound on the walker.
///
/// Fallible, with the error left to the implementor: repositioning reaches a file, and the
/// walker consumes an infallible item type, so an error met here has nowhere to go except
/// back to the caller that asked for the region.
pub trait RegionReadSource: Iterator<Item = PreparedRead> {
    /// What repositioning can fail with.
    type Error;

    /// Point the source at `region`. Every subsequent read belongs to it.
    ///
    /// # An implementor must **replay**, and this is not a detail
    ///
    /// A read this source has already handed out **must be offered again** to any later
    /// region that overlaps it. Consecutive regions overlap by design — each is asked for a
    /// halo past its end while the next is asked from its own start — so a source that
    /// consumed each read once would be short by every read that straddles a boundary.
    ///
    /// The walker relies on this in a way that is invisible from the outside: it holds one
    /// read of look-ahead, and [`move_to_region`](PileupWalker::move_to_region) *throws that
    /// look-ahead away* rather than carrying a previous region's read into the next region's
    /// walk. That read is not lost only because this method will offer it again.
    ///
    /// A source that does not replay therefore drops one read per region boundary, and it
    /// drops it **silently**: no error, no counter, just a genotype computed from less
    /// evidence than the file holds. The first attempt at this feature lost 3,830 of 236,081
    /// loci while all 1,471 unit tests passed (`spec/alignment_cursor.md` §6, §11).
    ///
    /// [`AlignmentCursor`](crate::ng::read::input::cursor::AlignmentCursor) honours this: it
    /// keeps every read it decodes and evicts one only once it ends before the current
    /// region begins.
    fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), Self::Error>;
}

/// A read source with one read of look-ahead, and a way to throw that look-ahead away.
///
/// **`Peekable` cannot do the second, and the second is what D1 needs.** The walk decides
/// where to advance to by looking at the next read without taking it. While a walker is
/// rebuilt per region that look-ahead is discarded along with the walker, so `Peekable` was
/// enough; a walker that stays alive across regions has to be able to drop it deliberately,
/// or it carries the previous region's peeked read into the next one's walk.
///
/// Behaviourally identical to `Peekable` today — every existing walk exercises it, and the
/// stage-1 differential and both dumps agree byte for byte. `forget_lookahead` is the one
/// thing it adds, and it is unused until the walker is pointed at a second region.
///
/// **Throwing the look-ahead away will lose no read, but only because of what sits
/// underneath**: a cursor keeps every read it hands out, so the next region is offered it
/// again. Against a source that does not keep its reads, discarding it would discard a read —
/// which is why this type is not offered as a general utility.
struct LookAhead<I: Iterator<Item = PreparedRead>> {
    inner: I,
    peeked: Option<PreparedRead>,
}

impl<I: Iterator<Item = PreparedRead>> LookAhead<I> {
    fn new(inner: I) -> Self {
        Self {
            inner,
            peeked: None,
        }
    }

    fn peek(&mut self) -> Option<&PreparedRead> {
        if self.peeked.is_none() {
            self.peeked = self.inner.next();
        }
        self.peeked.as_ref()
    }

    fn next(&mut self) -> Option<PreparedRead> {
        match self.peeked.take() {
            Some(read) => Some(read),
            None => self.inner.next(),
        }
    }

    /// Drop the look-ahead, so the next `peek` asks the source again.
    fn forget_lookahead(&mut self) {
        self.peeked = None;
    }
}

impl<I: RegionReadSource> LookAhead<I> {
    /// Point the source at `region`, **throwing the look-ahead away first**.
    ///
    /// The order is the whole reason this type exists. A read peeked but not taken was
    /// pulled for the region being left; carried into the next one it would be admitted
    /// first, ahead of reads that begin before it, and the walk would reject its own input
    /// as out of order — or worse, not.
    ///
    /// **Discarding it loses no read, and only because of what sits underneath**: a cursor
    /// keeps every read it hands out and offers it again to the next region that can use it
    /// (`spec/alignment_cursor.md` §6). Against a source that does not, this would discard a
    /// read, which is why `LookAhead` is not offered as a general utility.
    fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), I::Error> {
        self.forget_lookahead();
        self.inner.move_to_region(region)
    }
}

/// Pull-shaped walker over a coordinate-sorted stream of prepared
/// reads. See [`run`] for the convenience constructor.
pub struct PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeq,
{
    reads: LookAhead<I>,
    reference: F,
    state: WalkerState,
    /// Records produced by walker ticks but not yet consumed by
    /// `Iterator::next`. A single tick may emit 0–many records
    /// (e.g. a wide deletion at an earlier anchor unblocks several
    /// narrower records simultaneously); they're appended here in
    /// emission order and drained via `pop_front`.
    pending: VecDeque<SampleLocusObservations>,
    /// `true` once end-of-input has been flushed *or* a `next()`
    /// call has returned an error. Both terminal states stop the
    /// iterator from doing further work.
    done: bool,
    /// The last position worth walking, or `None` for production's
    /// unbounded walk — **ng's addition, C2 (plan 3).**
    ///
    /// A region walk queries a halo of `max_record_span` past the region's
    /// end, so a record anchored inside the region still sees the reads that
    /// fold into it from beyond the boundary (spec §2). Querying the halo is
    /// not enough: the walk has no right bound of its own, so it would walk
    /// all 5,000 halo positions at full depth, finalise every record in them
    /// and throw them away at the region clamp — a tax that, with regions
    /// tiling the genome, can exceed the region interiors.
    ///
    /// So the walk stops once it is past this position **and nothing anchored
    /// at or before it is still open**. The second half is what makes the stop
    /// safe rather than merely early: a record anchored inside the region can
    /// have a footprint running far past the boundary, and it is not finished
    /// until the walker has passed all of it.
    stop_after: Option<u32>,
}

impl<I, F> PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeq,
{
    pub fn new(reads: I, reference: F, config: &WalkerConfig) -> Self {
        let mut reads = LookAhead::new(reads);
        let mut state = WalkerState::new(*config);
        // Initial chromosome anchor: the first peeked read sets
        // `chrom_id`, `walker_pos = 1`, and
        // `last_admitted_chrom_id`. Subsequent re-anchors happen
        // inside the chromosome-boundary block of `fill_pending`,
        // right after `flush_chromosome_into`.
        if let Some(first) = reads.peek() {
            state.enter_chrom(first.chrom_id);
        }
        Self {
            reads,
            reference,
            state,
            pending: VecDeque::new(),
            done: false,
            stop_after: None,
        }
    }

    /// Walk with a chain-id allocator handed in from outside, replacing the one
    /// [`new`](Self::new) built — **ng's addition, C3 (plan 3).**
    ///
    /// Production builds a walker per chromosome, so its allocator can be the
    /// walker's own. ng walks one *region* at a time, and a fresh allocator per
    /// region would give two fragments of two regions the same chain id, which a
    /// later phasing step would chain together (spec §8). So the allocator lives
    /// on the generator and is lent to each walk; [`into_chain_ids`](Self::into_chain_ids)
    /// is how it comes back.
    ///
    /// Must be called before the walk starts — the walker is lazy, so "before
    /// the first `next()`" is all that means, and the assert says so in the
    /// build that can check it. Swapped in **mid-walk** it discards the ids
    /// already issued for the reads still active, and the allocations they
    /// represent go missing from the summary: one fragment, two identities,
    /// which is the corruption a run-lifetime allocator exists to prevent
    /// (review).
    ///
    /// `#[must_use]`, because a consuming builder called as a bare statement
    /// compiles, discards the walker and adopts nothing.
    #[must_use]
    pub fn adopting_chain_ids(mut self, chain_ids: ChainIdAllocator) -> Self {
        debug_assert!(
            self.state.summary.reads_admitted == 0,
            "adopting_chain_ids after {} reads: the ids already issued for the active reads \
             would be discarded",
            self.state.summary.reads_admitted,
        );
        self.state.chain_ids = chain_ids;
        self
    }

    /// Take the chain-id allocator back out at the end of a region's walk, so
    /// the next region continues the same `next_id` sequence.
    ///
    /// Consuming rather than swapping: a swap needs a placeholder allocator, and
    /// a placeholder that starts at zero is exactly the state this exists to
    /// avoid ever being in.
    pub fn into_chain_ids(self) -> ChainIdAllocator {
        self.state.chain_ids
    }

    /// Offer a finished record back, so the next column is filled into it — **G1.**
    ///
    /// Refused past the pool's bound; see [`RecordPool`](super::record_pool::RecordPool).
    pub fn recycle(&mut self, record: crate::ng::locus_generation::SampleLocusObservations) {
        self.state.record_pool.put(record);
    }

    /// The chain-id allocator's counters as they stand now — **ng's, D1.**
    ///
    /// The baseline a region's own contribution is measured against. A walker that lives
    /// for a chromosome never lets go of the allocator between regions, so the caller can
    /// no longer read the counters off the allocator it was about to lend; it reads them
    /// through here instead, at the moment the region opens. See
    /// `PileupGeneratorCounts::fold_region_walk` for what they are for and what goes wrong
    /// without them.
    pub fn chain_id_counters(&self) -> ChainIdAllocatorCounters {
        self.state.chain_ids.counters()
    }

    /// The read source, for a caller that needs to ask it something — **ng's, D1.**
    ///
    /// A walker that lives for a chromosome swallows its source for that long, and the
    /// source is where the cursor's tallies live. Without this the one thing that says
    /// whether the cursor is *working* — reads decoded against reads replayed
    /// (`spec/alignment_cursor.md` §11.5) — is unreachable from above, and the feature can
    /// be switched off with every test still green. It was, in review.
    ///
    /// Shared, not mutable: repositioning goes through
    /// [`move_to_region`](Self::move_to_region), which has a reset to run alongside it.
    pub fn reads(&self) -> &I {
        &self.reads.inner
    }

    /// Whether the walk has passed its right bound with nothing left that
    /// could still belong to the region.
    ///
    /// **Both halves matter.** `walker_pos > stop` alone would cut a record
    /// anchored inside the region short of its own footprint — the long
    /// deletions the halo exists for are exactly the records still open here.
    /// The open-record table is keyed by anchor, so "anything anchored at or
    /// before `stop`" is its first key.
    fn reached_stop(&self) -> bool {
        let Some(stop) = self.stop_after else {
            return false;
        };
        // **A held ordinary-column locus counts as an open record here**, and must: it
        // stands where a record the general path would still be holding stands, and a stop
        // rule blind to it cuts the region one position short — see `WalkerState::sealed`.
        let first_outstanding_anchor = match (
            self.state.sealed.as_ref().map(|l| l.region.start.0 as u32),
            self.state.open_records.first_open_anchor(),
        ) {
            (Some(held), Some(open)) => Some(held.min(open)),
            (Some(held), None) => Some(held),
            (None, Some(open)) => Some(open),
            (None, None) => None,
        };
        self.state.walker_pos > stop && first_outstanding_anchor.is_none_or(|anchor| anchor > stop)
    }

    /// Cumulative counters for the run so far. Safe to call
    /// mid-stream; the final summary is the value observed after
    /// `next()` has returned `None`.
    pub fn summary(&self) -> RunSummary {
        self.state.summary()
    }

    /// Drive the walker until at least one record is ready in
    /// `pending`, or until end-of-input. End-of-input also flushes
    /// the remaining chromosome and sets `done = true`.
    fn fill_pending(&mut self) -> Result<(), WalkerError> {
        loop {
            // Terminal condition (ng's, C2): the walk is past its right bound
            // and nothing anchored at or before it is still open. Flushed the
            // same way end-of-input is — the records still open are anchored
            // past the bound, so the region clamp discards them, but they are
            // emitted rather than dropped on the floor so the clamp can count
            // them.
            if self.reached_stop() {
                self.state.flush_chromosome_into(&mut self.pending)?;
                self.done = true;
                return Ok(());
            }

            // Terminal condition: no more reads to pull and the
            // active set is empty. Stopping at "reads empty" alone
            // would leak open records whose anchors sit ahead of
            // the walker but inside the active set's coverage.
            if self.reads.peek().is_none() && self.state.active_reads.is_empty() {
                self.state.flush_chromosome_into(&mut self.pending)?;
                self.done = true;
                return Ok(());
            }

            // Chromosome boundary: the next peeked read sits on a
            // new chromosome. Finalise everything still in flight
            // from the previous chromosome and re-anchor.
            //
            // The forward-direction check has to run *before* the
            // flush so a backward chromosome change errors out
            // without first emitting the previous chromosome's
            // records into `pending`.
            //
            // We pull the chrom_id and qname out of the peek
            // borrow into locals so the subsequent calls into
            // `self.state` aren't blocked by the peek borrow.
            let chrom_transition: Option<u32> = {
                let peeked = self.reads.peek();
                match (peeked, self.state.last_admitted_chrom_id) {
                    (Some(p), Some(prev)) if prev != p.chrom_id => {
                        if p.chrom_id < prev {
                            let prev_pos = self
                                .state
                                .last_admitted_locus
                                .map(|l| l.pos)
                                .unwrap_or(self.state.walker_pos);
                            return Err(WalkerError::OutOfOrder {
                                qname: p.qname.to_string(),
                                prev_chrom_id: prev,
                                prev_pos,
                                chrom_id: p.chrom_id,
                                pos: p.alignment_start,
                            });
                        }
                        Some(p.chrom_id)
                    }
                    _ => None,
                }
            };
            if let Some(new_chrom) = chrom_transition {
                self.state.flush_chromosome_into(&mut self.pending)?;
                self.state.enter_chrom(new_chrom);
            }

            // Pull every read with alignment_start ≤ walker_pos
            // (only on the current chromosome; reads on later
            // chromosomes wait for the chromosome flush above).
            while let Some(peeked_read) = self.reads.peek() {
                if peeked_read.chrom_id != self.state.chrom_id {
                    break;
                }
                if peeked_read.alignment_start > self.state.walker_pos {
                    break;
                }
                // PANIC-FREE: `peek()` returned Some on the loop
                // condition above, and `self.reads` has not been
                // advanced between then and here.
                let r = self.reads.next().expect("peek matched");
                self.state.admit_read(r)?;
            }

            // Process events at walker_pos, expire passed reads,
            // and close aged records. The expire-before-close
            // ordering also keeps the active-read count accurate
            // when an emitted record's footprint coincides with a
            // read's `alignment_end`.
            self.state
                .process_position(&self.reference, &mut self.pending)?;
            self.state.expire_passed_reads()?;
            self.state.close_aged_records_into(&mut self.pending);

            // Advance walker_pos to the next interesting position:
            // walker_pos+1 if any active read still has events at
            // or beyond it, otherwise jump to the next read's
            // alignment_start (skip uncovered span).
            self.state.advance(self.reads.peek())?;

            if !self.pending.is_empty() {
                return Ok(());
            }
        }
    }
}

impl<I, F> PileupWalker<I, F>
where
    I: RegionReadSource,
    F: RefSeq,
{
    /// Point this walker at `region`, and stop the walk once it is past `stop_after` with
    /// nothing anchored at or before that position still open — **ng's, D1.**
    ///
    /// This is what replaces building a walker per region. Everything scoped to the region
    /// being left is thrown away and everything scoped to the *run* is kept, which is the
    /// whole of the difficulty: see [`WalkerState::begin_region`] for the field-by-field
    /// decision and for the one field that must survive.
    ///
    /// **`stop_after` is passed rather than derived from `region`.** The walk's right bound
    /// and the span the source is pointed at are not the same coordinate — the caller asks
    /// the source for a halo past the region so a record anchored inside it still sees the
    /// reads that fold into it (`spec/alignment_cursor.md` §2), and then stops the walk at
    /// the region's own end so the halo is not walked at full depth and thrown away. Which
    /// span is which is the caller's to know; this type only needs the two numbers.
    ///
    /// **What the order here does and does not mean.** Both the reposition and the reset
    /// must happen before the re-anchoring peek at the end, and they do. Which of the two
    /// runs first is *not* load-bearing — an earlier version of this comment claimed the
    /// reposition had to lead "because the reset re-anchors, and that means peeking", which
    /// is untrue: the peek is after both. It leads because a source that cannot be
    /// repositioned leaves the walker untouched, which is the smaller mess.
    pub fn move_to_region(
        &mut self,
        region: GenomeRegion,
        stop_after: u32,
    ) -> Result<(), I::Error> {
        self.reads.move_to_region(region)?;
        self.state.begin_region();
        // Records produced by the region being left and never collected. The per-region
        // walker took them to the grave; so does this.
        self.pending.clear();
        self.done = false;
        self.stop_after = Some(stop_after);
        // The same anchor `new` takes, for the same reason: the first read's chromosome is
        // where the walk starts. A region with no reads leaves the state `begin_region`
        // set, which is what a freshly built walker over an empty source holds too.
        if let Some(first) = self.reads.peek() {
            self.state.enter_chrom(first.chrom_id);
        }
        Ok(())
    }
}

impl<I, F> Iterator for PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeq,
{
    type Item = Result<SampleLocusObservations, WalkerError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(record) = self.pending.pop_front() {
            return Some(Ok(record));
        }
        if self.done {
            return None;
        }
        match self.fill_pending() {
            Ok(()) => self.pending.pop_front().map(Ok),
            Err(e) => {
                // Terminal-on-first-error: stop yielding after this.
                // The previous push-based shape stopped emission as
                // soon as `run` returned `Err`; preserve that.
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Cumulative counters reported back by `run` so callers can log a
/// per-sample summary.
#[derive(Debug, Default, Clone, Copy)]
#[non_exhaustive]
pub struct RunSummary {
    pub reads_admitted: u64,
    pub records_emitted: u64,
    pub record_widen_events: u64,
    pub mate_overlap_positions: u64,
    pub chain_allocations: u64,
    pub active_reads_high_water: u32,
    pub mate_lookup_evictions: u64,
    /// Number of columns where the contributor list was truncated
    /// because depth exceeded the applicable per-column cap (see
    /// `WalkerConfig::max_snp_column_depth` /
    /// `max_indel_column_depth`). A non-zero value flags
    /// pathologically deep regions; QC pipelines may want to look
    /// at those samples / regions specifically.
    pub column_depth_truncations: u64,
    /// **ng's — production's `RunSummary` has no counterpart.**
    ///
    /// Reads the walk refused to admit because `max_active_reads` were already open.
    /// Non-zero means some region was deeper than the walk will hold, and the evidence
    /// there is a subsample of the reads that were available.
    ///
    /// **Not comparable with `column_depth_truncations`, which counts *positions*.** This
    /// counts *reads*, and the two caps act on different quantities: the column cap limits
    /// how many reads are used at one position, this one limits how many are held open at
    /// once. Where a region trips this cap the column cap becomes unreachable, because a
    /// position can no longer gather more contributors than the walk holds reads.
    ///
    /// Production cannot state it, so `parity.rs` binds it by name and drops it from the
    /// counter comparison.
    pub reads_shed_at_admission: u64,
    /// **ng's, added by D2 — production's `RunSummary` has no counterpart.**
    ///
    /// Reads that were admitted and left the active set having never been a contributor
    /// at any position: every base `N` or adaptor-masked, so the fold never heard of
    /// them and *neither per-locus counter can see them* — they produced no observation,
    /// but they also never reached the path that records `reads_without_observation`
    /// (spec §6). Read off the active set, which is where a read leaves.
    ///
    /// Because production cannot state it, `parity.rs` binds this field by name and
    /// drops it from the counter comparison; the exhaustive destructure there is what
    /// forces that decision to be made rather than defaulted.
    pub reads_silent_over_footprint: u64,
    /// **ng's, and production has no counterpart** — folded reads whose witness is more
    /// than one run, i.e. reads blind in the *middle* of a record's footprint. A spliced
    /// read across a record widened over its intron is the case that produces one, and
    /// recording it instead of discarding it is what this whole representation exists for.
    ///
    /// **Why it is on the walk's own summary and not only on the parity census.** The
    /// census counts the same thing, but it lives behind `#[cfg(test)]` and only measures
    /// loci where *production's* walker also produced a record — so it can never answer
    /// "how often does this fire on a real spliced BAM", which is the open measurement the
    /// change was made for (spec §8). Here, any BAM the generator is pointed at reports it,
    /// and the dump tool prints it (owner, 2026-07-31).
    ///
    /// Expected to read **zero on DNA-seq**, structurally: a ref-skip emits no event, so an
    /// intron cannot widen a record on its own, and modern Illumina puts `N`s at read ends
    /// where they cannot make a hole.
    pub reads_with_holed_witness: u64,
    /// The positions inside those holes, summed over the reads — what a holed read was
    /// blind over, which is the quantity that says how much evidence the old drop threw
    /// away rather than merely how many reads it threw away.
    pub hole_positions: u64,
    /// **ng's, added by the depth-cap change of 2026-08-05.**
    ///
    /// Reads the walk gave back after admitting them, because the hold ceiling was full
    /// and the arriving read had a smaller sampling key than one already held. The read
    /// counted here is the one *given up*, not the one that arrived.
    ///
    /// Read beside [`reads_shed_at_admission`](Self::reads_shed_at_admission): together
    /// the two are every read the ceiling removed from the evidence, and the split says
    /// which rule removed it. A run where both are zero is a run whose ceiling never
    /// shaped the output.
    pub reads_evicted_at_ceiling: u64,
    /// **ng's, added by the depth-cap change — and the owner's test of success.**
    ///
    /// Positions the walk folded **fewer reads than the per-position cap allows, while
    /// reads covering that position had been given up by the hold ceiling**. Every one
    /// of these is a position whose depth is lower than the BAM's for no reason the data
    /// justifies: not "too deep to fold", simply "the walk was not holding it".
    ///
    /// Counted against the ceiling's two losses only — refusals and evictions — because
    /// those are the ones with no defence. A position truncated by the per-position cap
    /// is not short: it folded exactly the cap.
    pub positions_short_of_cap: u64,
    /// The reads missing from those positions, summed — how *much* was lost, where
    /// [`positions_short_of_cap`](Self::positions_short_of_cap) says how often. Each
    /// position contributes the smaller of the gap to the cap and the number of
    /// ceiling-lost reads actually covering it, so it can never overstate.
    pub short_of_cap_deficit: u64,
    /// **ng's** — the deepest column the walk ever assembled, counted **before** the
    /// per-position cap and **after** the mate-overlap collapse. The quantity the caps
    /// are set against, so a run says whether its caps were near being reached rather
    /// than only whether they fired.
    pub column_depth_high_water: u32,
    /// **ng's** — the largest the chain-id allocator's first-mates-awaiting-a-partner map
    /// ever got, carried up from `ChainIdAllocatorCounters` so the headroom under
    /// `MAX_PENDING_MATES` is reportable.
    pub pending_mates_high_water: u32,
    /// **ng's — the four counters that say whether a capped position's evidence is a fair
    /// sample.** All four are touched only at a position the cap acts on, so they cost
    /// nothing on a run that never caps.
    ///
    /// A capped position keeps `cap` of the reads covering it. Whether that is a *sample*
    /// or a *selection* cannot be read off the kept reads alone — it needs the population
    /// they were drawn from. So both are counted, and both are counted twice: once
    /// altogether, and once for the reads that begin **left of the position**.
    ///
    /// Read as two fractions. `…_placed_left / …` over the seen reads is what the position
    /// actually looks like; the same over the kept reads is what the walk decided it looks
    /// like. **Equal fractions mean the cap sampled; a larger kept fraction means it
    /// preferred reads that reach the position from the left**, which is what keeping a
    /// prefix of an arrival-ordered list does, and what tilts `placed_left` and witness
    /// extent in everything downstream.
    pub capped_column_reads_seen: u64,
    pub capped_column_reads_seen_placed_left: u64,
    pub capped_column_reads_kept: u64,
    pub capped_column_reads_kept_placed_left: u64,
}

impl RunSummary {
    /// Fold the chain-id allocator's counters into this summary at
    /// run-end. The walker tracks reads/records/widens/overlap
    /// itself; the allocator owns chain-id bookkeeping.
    fn merge_chain_id_counters(mut self, counters: ChainIdAllocatorCounters) -> Self {
        self.chain_allocations = counters.chain_allocations;
        self.active_reads_high_water = counters.active_reads_high_water;
        self.mate_lookup_evictions = counters.mate_lookup_evictions;
        self.pending_mates_high_water = counters.pending_mates_high_water;
        self
    }

    // **Production's `merge` is deliberately not copied here — deleted at the
    // Milestone C review.** It totals one region's summary into another, which is
    // right for production (a fresh walker, and so a fresh chain-id allocator,
    // per region) and **wrong for ng**: ng shares one allocator across regions and
    // `reset()` preserves its counters, so summing region summaries
    // triangular-sums `chain_allocations` and `mate_lookup_evictions` (spec §8).
    // It had no ng caller, and a `pub fn merge` sitting on the type that
    // `PileupGeneratorCounts::fold_region_walk` exists to fold *correctly* is a
    // trap name-completion offers first. The exhaustive-destructure idiom it
    // carried lives on in `fold_region_walk` and in `parity.rs`'s summary
    // comparison.
}

/// A genomic locus: a position on a specific chromosome. The
/// `Ord` derive gives lexicographic ordering — chromosomes
/// compared first, position within a chromosome second — which is
/// exactly the walker's coordinate-sort invariant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Locus {
    chrom_id: u32,
    /// 1-based reference position.
    pos: u32,
}

struct WalkerState {
    chrom_id: u32,
    walker_pos: u32,
    /// `None` until the first read is admitted. Tracks the
    /// chromosome the walker has been processing so the
    /// flush-on-chromosome-change logic knows whether to flush
    /// before re-anchoring.
    last_admitted_chrom_id: Option<u32>,
    /// Last admitted locus for the coordinate-order invariant.
    last_admitted_locus: Option<Locus>,
    active_reads: ActiveReads,
    chain_ids: ChainIdAllocator,
    open_records: OpenPileupRecordTable,
    summary: RunSummary,
    config: WalkerConfig,
    /// Reusable per-step buffer for the contributors list. Hoisted
    /// here so the per-walker-step `Vec<ReadContribution>` is
    /// allocated once and reused via `clear()` between steps. L6
    /// in `ia/reviews/perf_pileup_2026-05-10.md`.
    contributors_buf: Vec<ReadContribution>,
    /// Reusable per-step buffer for the read ids the column-depth cap removed. Hoisted for
    /// the same reason as `contributors_buf`: this runs once per covered reference base,
    /// and a truncated column is rare, so the buffer is usually empty and never grows.
    truncated_read_ids_buf: Vec<u32>,
    /// **ng's, added by the depth-cap change** — the alignment ends of the reads the hold
    /// ceiling gave up, smallest first, so the walk can say at any position **how many
    /// reads covering it the ceiling took away**. That is the whole of
    /// `positions_short_of_cap`, and there is no other way to know it: a read the ceiling
    /// drops is never decomposed, so nothing downstream can be asked which positions it
    /// would have reached.
    ///
    /// Costs nothing when the ceiling never binds, which is the case this ships for: the
    /// heap is empty, and the per-position work is one `is_empty` on an empty `Vec`.
    ceiling_losses_by_end: BinaryHeap<Reverse<u32>>,
    /// **ng's, added by the depth-cap change** — scratch for the per-position sample:
    /// `(sampling key, index into the contributor list)` for every contributor, and then
    /// the indices the sample keeps. Hoisted for the reason every buffer here is: a deep
    /// enough region caps at many consecutive positions, and each would otherwise
    /// allocate twice.
    sample_keys_buf: Vec<(u64, u32)>,
    sample_kept_buf: Vec<u32>,
    /// EXPERIMENT E2: reusable scratch for `resolve_mate_overlap_at_pos`, replacing
    /// the per-column `AHashMap<ChainId, Vec<usize>>` and its two companion `Vec`s.
    mate_overlap_buf: MateOverlapScratch,
    /// Reusable per-step buffer for `close_aged_records`'s drained
    /// records. Paired with `OpenPileupRecordTable::closing_keys_buf`;
    /// together they remove the two per-walker-step `Vec` allocations
    /// `drain_aged` paid in round-1. H2 in
    /// `ia/reviews/perf_pileup_2026-05-12.md`.
    drained_buf: Vec<OpenPileupRecord>,
    /// Reusable buffers for the ordinary-column path — see [`fast_column`](super::fast_column).
    fast_column_buf: FastColumnScratch,
    /// **Retired locus records, waiting to be filled again — G1.** Scratch like the
    /// buffers above, and kept across regions for the same reason: the merge hands a
    /// record back once it has evicted it, and the next region's first locus is filled
    /// into it rather than allocated.
    ///
    /// **It is not carried across chromosomes, unlike the chain-id allocator**, and the
    /// difference is that the allocator's *contents* are load-bearing — two fragments of
    /// one read must not be given the same id — while a pool holds nothing but buffers.
    /// A walker is minted per chromosome, so a cold pool costs the allocations of one
    /// locus, twelve times on this reference. Carrying it would be machinery bought for
    /// twelve allocations a run.
    record_pool: super::record_pool::RecordPool,
    /// **Measurement knob only.** `PVC_FAST_COLUMN=0` sends every column down the general
    /// path, so the two can be A/B-ed inside one binary rather than across two builds — the
    /// only way to alternate runs on a host several other measurements are sharing. Read
    /// once, here, so the per-column cost of carrying it is one predictable branch.
    fast_column_enabled: bool,
    /// **The record the ordinary-column path did not have to open** — the locus it built,
    /// held for exactly as long as the general path would have held the record.
    ///
    /// The fast lane finishes a column's locus at the position it walks; the general path
    /// leaves a one-base record open and drains it one step later. **That one step is
    /// observable and had to be reproduced**, in two ways a fully-consumed walk cannot show
    /// and two tests do:
    ///
    /// - A walk that **aborts** — a reference fetch past the contig end — loses whatever is
    ///   still open. Emitting a step early hands the consumer one locus more before the
    ///   error (`parity::both_walkers_report_the_same_error_on_the_same_malformed_input`).
    /// - A walk that is **abandoned** part-way stops where its consumer stopped pulling, and
    ///   a locus offered a step early moves that point back by one position — so fewer reads
    ///   are admitted and fewer chain ids allocated
    ///   (`generator::tests::an_abandoned_walk_does_not_leak_its_active_reads_into_the_next_region`).
    ///
    /// **One slot is always enough.** The fast lane fires only when no record covers the
    /// base it is on, which means every record still open ends at or before it — so by the
    /// time the next position is reached the table is empty, and a locus held here is the
    /// only thing outstanding. That is also why emitting it first is always coordinate
    /// order.
    sealed: Option<SampleLocusObservations>,
}

impl WalkerState {
    fn new(config: WalkerConfig) -> Self {
        Self {
            chrom_id: 0,
            walker_pos: 1,
            last_admitted_chrom_id: None,
            last_admitted_locus: None,
            active_reads: ActiveReads::new(),
            // **The allocator is handed a cap it can never reach, on purpose.**
            //
            // `max_active_reads` is still the ceiling — `admit_read` enforces it, by
            // refusing the read rather than by failing the walk. Passing the same number
            // here as well would leave the allocator's two responses to a deep region live
            // underneath a walk that no longer wants either: a hard
            // `ActiveReadsExhausted` that can no longer be reached, and a one-shot warning
            // at three-quarters of the cap telling the user *"the run will fail"*, which
            // by then is untrue. Neither can be edited — `chain_id_allocator.rs` is locked
            // byte-identical to production's — so they are put out of reach instead.
            //
            // What is lost is a backstop on a bug in the shed itself. It was worth little:
            // its response to an over-full active set was to abort the run, which is the
            // behaviour being removed here, and the same over-full set is visible after
            // the fact in `active_reads_high_water`.
            chain_ids: ChainIdAllocator::with_caps(u32::MAX, config.mate_lookup_window),
            open_records: OpenPileupRecordTable::with_cap(config.max_record_span),
            summary: RunSummary::default(),
            config,
            contributors_buf: Vec::new(),
            truncated_read_ids_buf: Vec::new(),
            ceiling_losses_by_end: BinaryHeap::new(),
            sample_keys_buf: Vec::new(),
            sample_kept_buf: Vec::new(),
            mate_overlap_buf: MateOverlapScratch::default(),
            drained_buf: Vec::new(),
            fast_column_buf: FastColumnScratch::default(),
            record_pool: super::record_pool::RecordPool::new(),
            fast_column_enabled: std::env::var_os("PVC_FAST_COLUMN")
                .is_none_or(|value| value != "0"),
            sealed: None,
        }
    }

    /// Anchor the walker to a chromosome. Called twice in the
    /// run lifecycle: once before the loop with the first peeked
    /// read's chrom, and again inside the boundary block right
    /// after `flush_chromosome` when the next peeked read sits on
    /// a new chrom. Walker_pos resets per chromosome — the walker
    /// only emits within a chromosome, so position numbering
    /// restarts from 1.
    ///
    /// `last_admitted_locus` is deliberately preserved across
    /// chromosome boundaries (Mi14 in
    /// `ia/reviews/pileup_2026-05-11.md`): the per-read tuple
    /// comparison in `admit_read` correctly admits a forward
    /// chrom change (`(new_chrom, _) > (old_chrom, _)` holds
    /// whenever `new_chrom > old_chrom`), and keeping the locus
    /// sticky lets the outer chrom regression's error message
    /// report the *actual* last admitted (chrom, pos) instead of
    /// falling back to a misleading `walker_pos`.
    fn enter_chrom(&mut self, chrom_id: u32) {
        self.chrom_id = chrom_id;
        self.walker_pos = 1;
        self.last_admitted_chrom_id = Some(chrom_id);
    }

    /// Put this state back where a freshly built walker's would be, **except the chain-id
    /// allocator** — ng's, D1.
    ///
    /// # The trap, stated first because it is invisible in the output
    ///
    /// ng shares one chain-id allocator across every region of a chromosome, so two
    /// fragments of two regions never carry the same id. `ChainIdAllocator::reset` exists
    /// for that: it drops `pending_mates` and `active_count` and **preserves `next_id` and
    /// the three counters**. `PileupGeneratorCounts::fold_region_walk` then folds two of
    /// those counters as *deltas* against the value they held when the region opened.
    ///
    /// A reset that replaced the allocator — `WalkerState::new(config)` is one keystroke
    /// away and compiles — would zero those counters, collapse both deltas, and leave
    /// `active_reads_high_water` looking right because it is a max. That is what makes the
    /// corruption look selective enough to rationalise, and it has happened here before.
    ///
    /// The mirror-image trap sits one field away: `ActiveReads::reset` *preserves*
    /// `silent_exits` as a run total, and `fold_region_walk` sums that one **per region**.
    /// So the active set is put back with [`ActiveReads::begin_region`], which zeroes it,
    /// and not with `reset`.
    ///
    /// The destructure is exhaustive on purpose: a field added to this struct is a compile
    /// error here until someone decides which side of that line it falls on.
    fn begin_region(&mut self) {
        let Self {
            chrom_id,
            walker_pos,
            last_admitted_chrom_id,
            last_admitted_locus,
            active_reads,
            chain_ids,
            open_records,
            summary,
            // Not region-scoped: the knobs the walker was built with.
            config: _,
            // Scratch, deliberately untouched. Each is cleared at the point of use, and
            // keeping its capacity across a region boundary is the entire reason these are
            // fields rather than locals.
            contributors_buf: _,
            truncated_read_ids_buf: _,
            sample_keys_buf: _,
            sample_kept_buf: _,
            mate_overlap_buf: _,
            drained_buf: _,
            fast_column_buf: _,
            // Scratch, and kept for the same reason as the buffers above: a record handed
            // back during region N is exactly the right shape for region N+1's first locus.
            record_pool: _,
            // A knob the walker was built with, like `config`.
            fast_column_enabled: _,
            // Region-scoped, and **dropped rather than emitted** — for the reason the open
            // records a few lines below are dropped: it belongs to the region being left,
            // whose output nobody collects.
            sealed,
            // **Not scratch, and region-scoped.** It holds reads the *previous* region's
            // ceiling gave up; carried across, they would be counted as covering
            // positions in a region they may not even be on, and
            // `positions_short_of_cap` would blame this region for the last one's losses.
            ceiling_losses_by_end,
        } = self;
        *sealed = None;
        ceiling_losses_by_end.clear();

        // What `new` starts with. `enter_chrom` overwrites the first three as soon as the
        // new region's first read is peeked; they are set here so a region with no reads at
        // all is still in a defined state rather than the previous region's.
        *chrom_id = 0;
        *walker_pos = 1;
        *last_admitted_chrom_id = None;
        // **The one that would fire on real data rather than in a test.** The source is
        // asked for a halo past the region, so the last read admitted for region N can
        // begin far past region N+1's start. Carried across, the coordinate-order check in
        // `admit_read` would reject the next region's first read as going backwards.
        *last_admitted_locus = None;

        active_reads.begin_region();

        // Records still open belong to the region being left. Drained rather than
        // finalised: finalising would emit records nobody asked for and tally their
        // witnesses into a summary that is about to be discarded. The per-region walker
        // dropped them, and so does this.
        for _ in open_records.drain_all() {}
        open_records.reset();

        // Cleared, never replaced — see the trap above.
        chain_ids.reset();

        // The walk's own counters are per region, because `fold_region_walk` sums them.
        *summary = RunSummary::default();
    }

    fn admit_read(&mut self, read: PreparedRead) -> Result<(), WalkerError> {
        // Order-invariant check.
        let read_locus = Locus {
            chrom_id: read.chrom_id,
            pos: read.alignment_start,
        };
        if let Some(prev) = self.last_admitted_locus
            && read_locus < prev
        {
            return Err(WalkerError::OutOfOrder {
                qname: read.qname.to_string(),
                prev_chrom_id: prev.chrom_id,
                prev_pos: prev.pos,
                chrom_id: read_locus.chrom_id,
                pos: read_locus.pos,
            });
        }
        // Zero-ref-span check.
        if read.alignment_end < read.alignment_start {
            return Err(WalkerError::ZeroRefSpan {
                qname: read.qname.to_string(),
                chrom_id: read.chrom_id,
                pos: read.alignment_start,
            });
        }
        // Length invariants the cursor relies on. See
        // `PreparedRead::length` for the rationale.
        read.length()
            .map_err(|e| malformed_read_from_length_err(e, &read))?;

        // **ng's — the ceiling on what the walk holds, and the one place it can be
        // enforced.** The per-position caps act after the reads are already open, so they
        // bound what the walk *uses* and not what it *costs*; this bounds the two
        // structures that actually fill up on a deep region — the active set and the
        // allocator's map of first mates waiting for a partner — because a read that
        // never gets in enters neither.
        //
        // **The ceiling stays** (owner, 2026-08-05: *"still with a high enough cap,
        // otherwise we could run out of memory"*). What changed is its default, from
        // 4,096 to `DEFAULT_MAX_ACTIVE_READS`, and **how it chooses**.
        //
        // # First-come-first-served was the fault, not the ceiling
        //
        // The old rule refused whichever read arrived when the set was full. Reads arrive
        // sorted by alignment start, so that rule kept the leftmost-starting reads of a
        // deep region and threw away everything after — and a read it threw away
        // contributed at **no** position, which is what left positions covered by fewer
        // reads than the BAM had for them. The owner's words: *"what it is wrong is to
        // leave positions with less coverage because we have discarded reads that cover
        // it."*
        //
        // # What it does instead
        //
        // When the set is full, the arriving read and the reads already held are compared
        // on the **same deterministic function of the read** a capped position uses —
        // `read_sampling::sampling_key`, a hash of the query name. Smallest keys survive.
        // If the arrival beats the worst held read, that read is given back and the
        // arrival takes its place; otherwise the arrival is refused. Either way the set
        // ends up holding an unbiased subsample of the reads over the region instead of
        // its leftmost-starting prefix.
        //
        // **Nothing is paid for this on the ordinary path.** No key is stored on a read,
        // no index is kept in key order: the keys are computed by a scan, inside this
        // branch, which is entered only when the set is already full. At the ceiling this
        // ships, that is no position on any real fixture measured.
        //
        // **Eviction is an early expiry and adds no case to the fold.** A read given back
        // here may already have folded into open records; so may a read whose end the
        // walker passes while a record it folded into is still open, and
        // `refold_live_reads` has always skipped a folded read that is no longer live.
        // What it does mean is that the read's witness stops being updated as the record
        // widens — again, exactly what an ordinary expiry means.
        //
        // **What is still knowingly given up.** A refused read carries no `read_id` and is
        // never decomposed, so nothing downstream can say which loci it would have
        // reached; unlike the per-position cap this cannot feed a per-record
        // `reads_discarded_by_cap`. `reads_shed_at_admission`, `reads_evicted_at_ceiling`
        // and `positions_short_of_cap` are the whole report. And a read the ceiling
        // removes leaves its mate unpartnered, so that pair's two views of an overlap are
        // not collapsed.
        if self.active_reads.len() >= self.config.max_active_reads as usize {
            let arriving_key = read_sampling::sampling_key(&read);
            match self.active_reads.worst_sampling_key() {
                // Strictly worse, so an exact tie refuses the arrival — a rule that needs
                // no tie-break beyond the one already inside `worst_sampling_key`, and
                // that is reached about once in 2^64 draws.
                Some((worst_key, worst_read_id)) if worst_key > arriving_key => {
                    let end_of_evicted = self
                        .active_reads
                        .get_by_read_id(worst_read_id)
                        .map(|active| active.read.alignment_end);
                    self.active_reads.evict_by_read_id(
                        worst_read_id,
                        &mut self.chain_ids,
                        self.walker_pos,
                    )?;
                    self.summary.reads_evicted_at_ceiling += 1;
                    if let Some(end) = end_of_evicted {
                        self.ceiling_losses_by_end.push(Reverse(end));
                    }
                }
                _ => {
                    self.summary.reads_shed_at_admission += 1;
                    self.ceiling_losses_by_end.push(Reverse(read.alignment_end));
                    return Ok(());
                }
            }
        }

        self.last_admitted_locus = Some(read_locus);
        self.active_reads.admit(read, &mut self.chain_ids)?;
        self.summary.reads_admitted += 1;
        Ok(())
    }

    fn process_position<F: RefSeq>(
        &mut self,
        reference: &F,
        out: &mut VecDeque<SampleLocusObservations>,
    ) -> Result<(), WalkerError> {
        let walker_pos = self.walker_pos;
        // **Asked once per column, and read by both paths below.** See
        // [`ActiveReads::may_have_mate_overlap_at`]: `false` means no pair of this set's
        // reads has both alignments still on the reference here, so no two contributors can
        // share a chain id. The general path skips its whole mate-overlap step on that
        // answer; the ordinary-column path skips the sort it would otherwise run to reach
        // the same conclusion. It is hoisted above both because it needs the active set
        // mutably — pruning is how the heap stays O(1) to consult — and both paths below
        // borrow it immutably.
        let may_have_mate_overlap = self.active_reads.may_have_mate_overlap_at(walker_pos);
        // **The ordinary column, answered in scalars.** Roughly eight columns in ten are one
        // covered base at which every read simply matches; for those the general machine
        // below is asked to express a handful of numbers. The predicate is decided before any
        // per-read work and handing back costs nothing — see `fast_column`.
        //
        // **The one condition the fast lane cannot answer for is a position the hold ceiling
        // took reads from** (2026-08-05, merging the ordinary-column path with the depth-cap
        // change). The fast lane emits its locus and returns, above the block below that
        // prunes `ceiling_losses_by_end` and raises `positions_short_of_cap` — so a fast
        // column would leave the owner's test of success unable to fire, and would let the
        // heap grow across a run of columns that never prune it. While the heap holds
        // anything the general path takes the column, which is where the accounting lives.
        //
        // It costs one `is_empty` per column. The heap receives a push only when the set is
        // at the ceiling, and on every fixture measured at the shipping ceiling of 32,768 it
        // never receives one — so this is a branch that is false for the whole of a normal
        // run, and the fast lane's coverage is unchanged wherever the ceiling does not bind.
        let attempt = if self.fast_column_enabled && self.ceiling_losses_by_end.is_empty() {
            super::fast_column::try_ordinary_column(
                walker_pos,
                self.chrom_id,
                &self.active_reads,
                &self.open_records,
                reference,
                &mut self.fast_column_buf,
                &mut self.record_pool,
                self.config.max_snp_column_depth as usize,
                may_have_mate_overlap,
            )?
        } else {
            super::fast_column::FastColumn::Fallback
        };
        match attempt {
            super::fast_column::FastColumn::Emitted(locus) => {
                // **The previous position's locus goes out here, and only after this one was
                // built without error.** It is one step old, which is where the general path
                // would have drained its record, and the fetch that could still fail has
                // already succeeded — so an aborting walk loses exactly what it lost before.
                // The table is empty whenever anything is held (see `sealed`), so there is
                // nothing this could be pushed ahead of.
                self.emit_held_locus_into(out);
                self.sealed = Some(locus);
                return Ok(());
            }
            super::fast_column::FastColumn::Fallback => {}
        }

        // Step 1: query each active read's cursor for events
        // anchored at walker_pos. Reads with no event here are
        // silent (deletion interior or N-skip), so they are not
        // added as contributors at all.
        // Hoisted buffer; cleared per step. L6.
        self.contributors_buf.clear();
        let contributors = &mut self.contributors_buf;

        for (active_index, active_read) in self.active_reads.iter().enumerate() {
            let events_at_pos = active_read.cursor.events_at(walker_pos, &active_read.read);

            if events_at_pos.is_empty() {
                continue;
            }

            // Fallback to BQ=0 when this contributor has only indel events at
            // walker_pos (no Match). BQ=0 → ln(P_err)=0 in `phred_to_ln_perr`,
            // so the contributor's Match-side q_sum contribution is zero; the
            // indel BQ itself flows through events_at_pos and is folded
            // separately. Not a recovered error — there is no Match BQ here.
            let bq_at_walker = events_at_pos
                .iter()
                .find_map(|e| match e {
                    ReadEvent::Match { bq_baq, .. } => Some(*bq_baq),
                    ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. } => None,
                })
                .unwrap_or(0);

            // **ng's, added by D2.** Set here — before the mate-overlap collapse and
            // before the depth cap — because both of those remove reads the walk plainly
            // saw, and each has its own counter. What this flag exists to find is the read
            // that reaches *no* contributor list anywhere: every base `N` or adaptor-masked,
            // admitted and expired without the fold ever hearing of it, and so invisible to
            // both per-locus counters (spec §6).
            active_read.ever_contributed.set(true);
            contributors.push(ReadContribution {
                read_id: active_read.read_id,
                active_index: active_index as u32,
                chain_id: active_read.chain_id,
                events_at_pos,
                bq_baq_at_walker_pos: bq_at_walker,
                alignment_start: active_read.read.alignment_start,
                mate_role: active_read.read.mate_role,
                bq_zero_in_window: false,
                bq_override_at_walker_pos: None,
            });
        }

        // Step 2: resolve mate overlap on events at this walker
        // position. For each pair of contributors whose reads
        // share a chain_id, compare BQ; the lower-BQ side
        // has its `bq_baq_at_walker_pos` zeroed (still one
        // observation, contributing ln(1)=0 log-likelihood mass)
        // and is flagged so any window event the fold pulls from
        // its cursor also gets BQ-zeroed.
        //
        // **Skipped outright wherever no pair is present, which is most columns at the
        // depths this walk is for.** Reconciled columns are 1,664 in 10,000 on a 130×
        // tomato whole-genome sample and 1,911 in 10,000 on a 30× human one; the step's own
        // no-pair exit still costs a depth-sized tuple build and sort at every one of the
        // other eight in ten. The active set can rule the column out in O(1) instead,
        // because a shared chain id is a mate pair and a mate pair is known at admission
        // (`ActiveReads::may_have_mate_overlap_at`). At 300× a pair is present at most
        // columns and the skip simply stops firing.
        if may_have_mate_overlap {
            resolve_mate_overlap_at_pos(
                contributors,
                &mut self.summary,
                &mut self.mate_overlap_buf,
            );
        } else {
            debug_assert!(
                !column_shares_a_chain_id(contributors),
                "the mate-overlap skip fired at {}:{} on a column where two contributors \
                 share a chain id — the reconciliation was silently lost",
                self.chrom_id,
                walker_pos,
            );
        }

        // Step 2b: per-column depth cap. Adopted from samtools'
        // mpileup (see `WalkerConfig` doc-comment). Apply *after*
        // mate-overlap so the cap counts genuine post-collapse
        // observations, not per-mate. Detect "indel column" from
        // the post-collapse contributor events — any Insertion or
        // Deletion at this anchor flips the column to the tighter
        // indel cap.
        //
        // **Which reads survive is a fact about the reads** (ng's, 2026-08-05). It used
        // to be `contributors.truncate(cap)` — keep the first `cap` in whatever order the
        // active set happened to hold them. That had two faults, and the second is the
        // one that matters: the order is a permutation produced by `swap_remove`, so a
        // change to how the set stores reads moved 88,351 of 341,094 emitted rows on a
        // 300× walk; and wherever arrival order *was* preserved, "first `cap`" meant
        // "leftmost-starting `cap`", which tilts `placed_left` and witness extent.
        // `read_sampling` replaces it with the `cap` smallest sampling keys — see that
        // module for what the rule does and does not guarantee.
        //
        // **The ids the cap removes are kept** (B3). The locus type reports reads a cap
        // discarded *per record*, and this is the only moment they are knowable: the drop
        // happens before step 3, so a dropped read opens and widens nothing and is
        // invisible to every record it would have reached. They are candidates rather than
        // a count — see `OpenPileupRecord::reads_discarded_by_cap`.
        let cap = column_depth_cap(contributors, &self.config);
        let depth = contributors.len();
        if depth as u32 > self.summary.column_depth_high_water {
            self.summary.column_depth_high_water = depth as u32;
        }
        self.truncated_read_ids_buf.clear();
        if depth > cap {
            // The fairness census, taken before and after the cap acts — see the four
            // fields on `RunSummary`. A read is "placed left" when it began before this
            // position, i.e. it reaches the position from the left rather than starting
            // at it.
            self.summary.capped_column_reads_seen += depth as u64;
            self.summary.capped_column_reads_seen_placed_left += contributors
                .iter()
                .filter(|contrib| contrib.alignment_start < walker_pos)
                .count() as u64;
            sample_to_cap(
                contributors,
                cap,
                &self.active_reads,
                &mut self.sample_keys_buf,
                &mut self.sample_kept_buf,
                &mut self.truncated_read_ids_buf,
            );
            self.summary.capped_column_reads_kept += contributors.len() as u64;
            self.summary.capped_column_reads_kept_placed_left += contributors
                .iter()
                .filter(|contrib| contrib.alignment_start < walker_pos)
                .count() as u64;
            self.summary.column_depth_truncations += 1;
        }

        // **The owner's test of success, counted rather than argued.** The heap holds the
        // alignment ends of the reads the hold ceiling gave up; anything ending before the
        // walker no longer covers this position and is dropped. What is left is exactly
        // the reads that cover *here* and are not in the fold because the walk was not
        // holding them — so a position below the cap folded fewer reads than the cap
        // allows for a reason the data does not justify.
        //
        // A position *at or above* the cap is not short, whatever the ceiling did: it
        // folded the cap's worth, which is all it was ever going to fold. The pruning
        // still runs there, so the heap cannot grow across a run of capped positions.
        //
        // The whole block is behind an `is_empty` on a heap that never received a push
        // unless the ceiling bound, so on a run where it never binds this is one branch
        // per position and nothing else.
        if !self.ceiling_losses_by_end.is_empty() {
            while let Some(&Reverse(end)) = self.ceiling_losses_by_end.peek() {
                if end < walker_pos {
                    self.ceiling_losses_by_end.pop();
                } else {
                    break;
                }
            }
            let lost_here = self.ceiling_losses_by_end.len();
            if lost_here > 0 && depth < cap {
                self.summary.positions_short_of_cap += 1;
                self.summary.short_of_cap_deficit += (cap - depth).min(lost_here) as u64;
            }
        }

        // **Measurement only** — how often the column the walk is about to fold is the
        // ordinary one. Off unless `PVC_COLUMN_CENSUS=1`; see `column_census`.
        if super::column_census::enabled() && !contributors.is_empty() {
            use super::column_census as census;
            census::add(&census::COLUMNS, 1);
            census::add(&census::CONTRIBUTORS, contributors.len() as u64);

            let record_already_open = self
                .open_records
                .find_overlapping(walker_pos, walker_pos.saturating_add(1))
                .is_some();
            let indel_event = contributors.iter().any(|c| {
                c.events_at_pos
                    .iter()
                    .any(|e| !matches!(e, ReadEvent::Match { .. }))
                    || c.events_at_pos.len() != 1
            });
            let read_has_deletion = contributors.iter().any(|c| {
                self.active_reads
                    .get_by_read_id(c.read_id)
                    .is_some_and(|a| !a.cursor.spans_only_its_anchors())
            });
            let read_has_indel = contributors.iter().any(|c| {
                self.active_reads
                    .get_by_read_id(c.read_id)
                    .is_some_and(|a| !a.cursor.matches_only())
            });
            let mate_overlap = contributors
                .iter()
                .any(|c| c.bq_zero_in_window || c.bq_override_at_walker_pos.is_some());
            let depth_cap = !self.truncated_read_ids_buf.is_empty();
            let first_group = self
                .active_reads
                .get_by_read_id(contributors[0].read_id)
                .map(|a| a.read.read_group);
            let multi_read_group = contributors.iter().any(|c| {
                self.active_reads
                    .get_by_read_id(c.read_id)
                    .map(|a| a.read.read_group)
                    != first_group
            });

            if record_already_open {
                census::add(&census::REJECT_RECORD_ALREADY_OPEN, 1);
            }
            if indel_event {
                census::add(&census::REJECT_INDEL_EVENT, 1);
            }
            if read_has_deletion {
                census::add(&census::REJECT_READ_HAS_DELETION, 1);
            }
            if mate_overlap {
                census::add(&census::REJECT_MATE_OVERLAP, 1);
            }
            if depth_cap {
                census::add(&census::REJECT_DEPTH_CAP, 1);
            }
            if multi_read_group {
                census::add(&census::REJECT_MULTI_READ_GROUP, 1);
            }
            if read_has_indel {
                census::add(&census::REJECT_READ_HAS_INDEL, 1);
            }
            if !(record_already_open
                || indel_event
                || read_has_deletion
                || mate_overlap
                || depth_cap
                || multi_read_group)
            {
                census::add(&census::COLUMNS_ORDINARY, 1);
                census::add(&census::CONTRIBUTORS_ORDINARY, contributors.len() as u64);
            }
            // The predicate a cheap per-read probe can actually decide: no contributor's
            // read carries an indel op anywhere, so `events_at` can only be one `Match`.
            let simple = !(record_already_open || read_has_indel || depth_cap || multi_read_group);
            if simple && !mate_overlap {
                census::add(&census::COLUMNS_SIMPLE, 1);
                census::add(&census::CONTRIBUTORS_SIMPLE, contributors.len() as u64);
            }
            if simple {
                census::add(&census::COLUMNS_SIMPLE_WITH_MATE, 1);
                census::add(
                    &census::CONTRIBUTORS_SIMPLE_WITH_MATE,
                    contributors.len() as u64,
                );
            }
        }

        // Step 3–6: fold contributors into the records affected
        // at this walker_pos. The fold queries each contributor's
        // cursor through `&self.active`; the returned outcome
        // counts only the records that *widened* during this
        // call — fresh opens and re-finds against an
        // already-large-enough footprint do not count.
        let outcome = process_position(
            &mut self.open_records,
            walker_pos,
            self.chrom_id,
            contributors,
            &self.truncated_read_ids_buf,
            &self.active_reads,
            reference,
        )?;
        self.summary.record_widen_events += outcome.widen_count;

        Ok(())
    }

    /// Finalise any open records whose footprint is fully behind
    /// the walker and append them to `out` in emission order.
    /// Returns without touching `out` if there are no aged records
    /// to drain.
    fn close_aged_records_into(&mut self, out: &mut VecDeque<SampleLocusObservations>) {
        // The ordinary column's locus, if one is a step old. First, because the table is
        // empty whenever one is held — see `sealed`.
        self.emit_held_locus_into(out);
        self.open_records
            .drain_aged_into(self.walker_pos, &mut self.drained_buf);
        if self.drained_buf.is_empty() {
            return;
        }
        // `finalise()` consumes each `OpenPileupRecord` by value, so
        // we drain the hoisted buffer rather than `into_iter()`ing
        // it; the backing `Vec` stays allocated and reusable.
        for open in self.drained_buf.drain(..) {
            // The witness tally is resolved at `finalise`. **Resolved there** because a
            // witness is a read's extent measured against the record's *final* footprint,
            // which only `finalise` knows and no later caller can reconstruct — the reads may
            // have expired. Two of its counts the locus already carries
            // (`reads_without_observation`, `reads_discarded_by_cap`) and the complete/partial
            // split is read only by tests; the **holed** pair is kept, because nothing else in
            // a non-test build can state how often a read was blind in the middle of a record.
            let (record, witness, storage) = open.finalise_recycling();
            self.open_records.recycle(storage);
            self.summary.reads_with_holed_witness += u64::from(witness.reads_with_holed_witness);
            self.summary.hole_positions += u64::from(witness.hole_positions);
            out.push_back(record);
            self.summary.records_emitted += 1;
        }
    }

    /// Hand over an ordinary column's locus once the walker has passed it — the moment the
    /// general path drains the one-base record it would have left open. See `sealed`.
    fn emit_held_locus_into(&mut self, out: &mut VecDeque<SampleLocusObservations>) {
        if self
            .sealed
            .as_ref()
            .is_some_and(|locus| locus.region.end.0 < u64::from(self.walker_pos))
        {
            // PANIC-FREE: the predicate above only holds on `Some`.
            out.push_back(self.sealed.take().expect("just tested as Some"));
            self.summary.records_emitted += 1;
        }
    }

    fn expire_passed_reads(&mut self) -> Result<(), WalkerError> {
        // Reads whose alignment_end < walker_pos can no longer
        // contribute. Their last event (if any) was processed at
        // their alignment_end position; once walker advances past
        // that they're done.
        debug_assert!(
            self.walker_pos > 0,
            "walker_pos starts at 1 and never decreases below 1",
        );
        self.active_reads
            .expire_passed(self.walker_pos, &mut self.chain_ids)
    }

    fn advance(&mut self, next_pulled: Option<&PreparedRead>) -> Result<(), WalkerError> {
        // Default: one position forward, so any active read's REF
        // contribution gets folded at every position it sits on.
        // `checked_add` guards against the (extreme) case of a
        // chromosome longer than `u32::MAX` bp — silent wrap to 0
        // would corrupt all subsequent record positions. Realistic
        // mammal/plant genomes don't approach this, but a few
        // salamander/lungfish genomes do (≥ 4 Gbp).
        let mut next_pos = self
            .walker_pos
            .checked_add(1)
            .ok_or_else(|| WalkerError::Internal {
                detail: format!("walker_pos overflowed u32 at {}", self.walker_pos),
                qname: String::new(),
                chrom_id: self.chrom_id,
                pos: self.walker_pos,
            })?;

        // If the active set is empty and the next pulled read
        // starts past the walker, skip the uncovered span.
        if self.active_reads.is_empty()
            && let Some(peeked_read) = next_pulled
            && peeked_read.chrom_id == self.chrom_id
            && peeked_read.alignment_start > self.walker_pos
        {
            next_pos = peeked_read.alignment_start;
        }

        self.walker_pos = next_pos;
        Ok(())
    }

    /// Finalise everything still in flight at end-of-chromosome
    /// (or end-of-input), appending the records to `out` in
    /// emission order.
    fn flush_chromosome_into(
        &mut self,
        out: &mut VecDeque<SampleLocusObservations>,
    ) -> Result<(), WalkerError> {
        // The held ordinary-column locus goes with the records it stood among, and first:
        // its anchor is the smallest of them. Unconditional here, as `drain_all` is — a
        // chromosome flush closes everything still in flight regardless of the walker.
        if let Some(locus) = self.sealed.take() {
            out.push_back(locus);
            self.summary.records_emitted += 1;
        }
        // Drain remaining open records (anything that was still
        // open at end-of-chromosome is by definition ready to
        // close — there are no future reads on this chromosome).
        for open in self.open_records.drain_all() {
            // Same as `close_aged_records_into` — see the note there. The storage is
            // dropped rather than handed back: this runs once per chromosome, against
            // 86 million times for the loop that does hand it back, so the pool is
            // better left holding whatever the walk was using mid-chromosome.
            let (record, witness, _storage) = open.finalise_recycling();
            self.summary.reads_with_holed_witness += u64::from(witness.reads_with_holed_witness);
            self.summary.hole_positions += u64::from(witness.hole_positions);
            out.push_back(record);
            self.summary.records_emitted += 1;
        }
        // Release any active-set reads so the active-count
        // bookkeeping is correct at the chromosome boundary. With
        // unique chain ids there are no per-record lifecycle marks
        // to stamp — the released ids are simply unused going
        // forward, and the next chromosome's first read mints a
        // fresh id that has never appeared in the file.
        self.active_reads
            .flush_all(&mut self.chain_ids, self.walker_pos)?;
        // Reset chromosome-scoped state. `self.open_records.reset()`
        // keeps the perf-hoisted `allele_seq_buf` capacity across
        // the chromosome boundary (Mi11). `chain_ids.reset()` clears
        // the active-read count and pending-mates map but
        // preserves the file-scoped `next_id` counter so chain
        // ids stay unique across chromosomes.
        self.chain_ids.reset();
        self.active_reads.reset();
        self.open_records.reset();
        // **ng's** — reads do not span chromosomes, and `walker_pos` restarts at 1, so an
        // end left here would never be popped and would make every position of the next
        // chromosome look as if a read were missing from it.
        self.ceiling_losses_by_end.clear();
        Ok(())
    }

    fn summary(&self) -> RunSummary {
        let mut summary = self
            .summary
            .merge_chain_id_counters(self.chain_ids.counters());
        // Read off the active set at every ask rather than accumulated as the walk goes:
        // the set is where a read *leaves*, and asking it means the number cannot drift
        // from the exits that produced it. Reads still active have not left, so a summary
        // taken mid-walk reports the reads that have — which is what every other counter
        // here does too.
        summary.reads_silent_over_footprint = self.active_reads.silent_exits();
        summary
    }
}

/// Wrap a [`ReadLengthError`] (which carries only raw lengths) into
/// a [`WalkerError::MalformedRead`] with the offending read's locus
/// context attached.
fn malformed_read_from_length_err(err: ReadLengthError, read: &PreparedRead) -> WalkerError {
    let reason = match err {
        ReadLengthError::SeqBqMismatch {
            seq_len,
            bq_baq_len,
        } => format!("seq.len ({seq_len}) != bq_baq.len ({bq_baq_len})"),
        ReadLengthError::CigarSeqMismatch {
            cigar_consumed,
            seq_len,
        } => format!("CIGAR consumes {cigar_consumed} read bases but seq.len = {seq_len}"),
    };
    WalkerError::MalformedRead {
        reason,
        qname: read.qname.to_string(),
        chrom_id: read.chrom_id,
        pos: read.alignment_start,
    }
}

/// The scratch space `resolve_mate_overlap_at_pos` needs for one column, hoisted so it is
/// reused instead of rebuilt.
///
/// It used to build a hash map keyed by chain id plus one `Vec` per distinct chain id — lists
/// almost always of length one — and two more `Vec`s besides. Grouping by sorting a single
/// `(chain_id, index)` vector gives the same runs for a `clear()` and a sort per column
/// rather than `1 + N + 2` allocations.
#[derive(Debug, Default)]
struct MateOverlapScratch {
    /// `(chain_id, contributor index)`, sorted so equal chain ids form one run.
    by_chain_id: Vec<(ChainId, usize)>,
    to_remove: Vec<usize>,
    bq_updates: Vec<(usize, u8, bool)>,
}

/// Does any pair of contributors at this column share a chain id?
///
/// **The pin on the O(1) skip, and it only exists in debug builds.** The skip's claim is
/// that a column the active set rules out cannot contain two contributors with one chain
/// id; a wrong skip loses a mate reconciliation, changes the emitted bytes, and shows up
/// nowhere else — `mate_overlap_positions` would simply be smaller. This is the all-pairs
/// scan the skip replaces, asserted every time the skip fires under `cargo test` and under
/// any debug-built dump, at the cost the old code paid on every column of every build.
///
/// Not `#[cfg(debug_assertions)]`: `debug_assert!` expands to a `cfg!`-guarded `assert!`,
/// so the call is type-checked in every build and elided from the release one.
fn column_shares_a_chain_id(contributors: &[ReadContribution]) -> bool {
    contributors.iter().enumerate().any(|(i, a)| {
        contributors[i + 1..]
            .iter()
            .any(|b| a.chain_id == b.chain_id)
    })
}

/// Resolve mate-overlap at the current walker position.
///
/// Two regimes, distinguished by whether either side carries an
/// indel anchored at this position:
///
/// - **Match-only overlap.** Both mates have only `Match` events
///   at this anchor. The lower-BQ side has its event BQs zeroed
///   in the local fold (so its `q_sum` contribution becomes
///   `ln(1) = 0`); both still count as observations and both
///   contribute the shared chain id. Tie-break: first mate
///   keeps its BQ.
///
/// - **Indel overlap.** Either both mates report an indel at the
///   same anchor, or one reports an indel and the other a clean
///   Match (mates disagree on indel presence). The pair collapses
///   to a single observation: the loser is removed from the
///   contributor list at this walker step, so it contributes
///   nothing to the anchor record. Tie-break on `bq_baq_at_walker_pos`
///   (Match BQ where present, indel `bq_proxy` mapped through 0
///   when the loser carries no Match here); ties go to the first
///   mate.
// `&mut Vec<_>` is intentional: the function `swap_remove`s
// contributors on indel-overlap, which requires the owning `Vec`,
// not a slice.
#[allow(clippy::ptr_arg)]
fn resolve_mate_overlap_at_pos(
    contributors: &mut Vec<ReadContribution>,
    summary: &mut RunSummary,
    scratch: &mut MateOverlapScratch,
) {
    // Build a small index: chain_id → list of contributor
    // indices. Anything with a list length >= 2 is a candidate.
    // ahash::AHashMap matches the rest of the module — std HashMap's
    // RandomState would make iteration non-deterministic between runs
    // and is also slower for this hot path. Mi4 in
    // `ia/reviews/pileup_2026-05-09.md`.
    //
    // EXPERIMENT E2: built as a sorted `(chain_id, index)` list in a reused buffer
    // rather than an `AHashMap<ChainId, Vec<usize>>` rebuilt per column. Groups are
    // the runs of equal chain id; within a run the indices stay ascending, which is
    // the order `values()` produced, and the two loops below are order-independent
    // across groups (each contributor belongs to exactly one group, and `to_remove`
    // is sorted and deduped before it is applied).
    let MateOverlapScratch {
        by_chain_id,
        to_remove,
        bq_updates,
    } = scratch;
    by_chain_id.clear();
    by_chain_id.extend(
        contributors
            .iter()
            .enumerate()
            .map(|(i, c)| (c.chain_id, i)),
    );
    by_chain_id.sort_unstable();

    // **The no-pair exit, read off the sort instead of hunting for it.** Mate overlap
    // needs two contributors at this position sharing a chain id, and in paired-end
    // data most columns have none — so this exit is the common case and its cost is
    // what matters.
    //
    // It used to be an all-pairs scan over the unsorted contributors, breaking as soon
    // as a shared id turned up. That breaks early only when a pair *exists*; the
    // common case ran the full n(n−1)/2 comparisons, which makes the cheap path the
    // quadratic one. Its comment priced it at "typical n ≤ ~30 contributors per
    // column", and at that depth it was the right call. A whole-genome sample at ~130×
    // puts n near 130, where the scan is ~19× more work per column than at 30 and the
    // function became the largest single site in the walk (13.8 % of a tomato
    // `SL4.0ch01` profile, against 4.2 % on a 30× fixture).
    //
    // Sorting first and looking for an adjacent equal pair answers the same question in
    // n log n, and the sort is not new work: every column that *does* have a pair
    // needed it anyway, three lines below. Equal chain ids are adjacent after a sort on
    // the `(chain_id, index)` tuple, so "some id repeats" and "some neighbour repeats"
    // are the same statement.
    if !by_chain_id.windows(2).any(|w| w[0].0 == w[1].0) {
        return;
    }

    // Indices to discard outright (indel-overlap losers).
    to_remove.clear();
    // (idx, new_bq_at_walker_pos, zero_in_window) — applied to
    // each contributor of a match-only overlap pair (S7). Agree-
    // case keeper gets the summed BQ (capped at 200, zero_in_window
    // = false); disagree-case winner gets `0.8 * bq` truncated
    // (zero_in_window = false). Other / loser gets new_bq=0 with
    // zero_in_window=true. The fold honours `bq_zero_in_window`
    // (zeros every window event from this contributor's cursor)
    // and `bq_override_at_walker_pos` (rewrites walker_pos events'
    // BQ on top of the cursor pull).
    bq_updates.clear();

    let mut run_start = 0usize;
    while run_start < by_chain_id.len() {
        let chain = by_chain_id[run_start].0;
        let mut run_end = run_start + 1;
        while run_end < by_chain_id.len() && by_chain_id[run_end].0 == chain {
            run_end += 1;
        }
        let group = &by_chain_id[run_start..run_end];
        run_start = run_end;
        if group.len() < 2 {
            continue;
        }
        // Spec invariant: only mate pairs share a chain id, so at
        // most two contributors per chain. Assert here so a future
        // change that admits a third reader of the same chain
        // (e.g. supplementary alignments slipping past upstream
        // filters) surfaces in tests instead of in production.
        debug_assert!(
            group.len() <= 2,
            "more than two contributors share chain_id {:?}",
            group,
        );
        // All-pairs comparison so a future relaxation of the
        // invariant doesn't silently miss the (i, j>i+1) cases
        // that `indices.windows(2)` skips.
        for i in 0..group.len() {
            for j in (i + 1)..group.len() {
                let (a, b) = (group[i].1, group[j].1);
                summary.mate_overlap_positions += 1;
                let any_indel_here = pair_has_indel(&contributors[a], &contributors[b]);
                if any_indel_here {
                    // Indel on either side at this walker_pos:
                    // collapse to a single observation by removing
                    // the loser entirely. Tie-break: BQ first,
                    // then first-of-pair, then alignment_start.
                    let loser_idx = pick_overlap_loser(contributors, a, b);
                    to_remove.push(loser_idx);
                } else {
                    // Match-only mate overlap (S7): apply
                    // samtools-style BQ math.
                    //
                    // PANIC-FREE: inside the !any_indel_here branch, every
                    // event at walker_pos on either side is a Match by
                    // definition of `pair_has_indel`, so `match_base_at_pos`
                    // returns Some.
                    let base_a = match_base_at_pos(&contributors[a])
                        .expect("match-only overlap: each side has a Match event at walker_pos");
                    let base_b = match_base_at_pos(&contributors[b])
                        .expect("match-only overlap: each side has a Match event at walker_pos");
                    if base_a == base_b {
                        // Agree case: sum BQs (cap 200), keeper
                        // takes the sum, other is zeroed.
                        let combined_bq = sum_bq_capped_at_200(
                            contributors[a].bq_baq_at_walker_pos,
                            contributors[b].bq_baq_at_walker_pos,
                        );
                        let keeper_idx = pick_agree_keeper(contributors, a, b);
                        let other_idx = if keeper_idx == a { b } else { a };
                        bq_updates.push((keeper_idx, combined_bq, false));
                        bq_updates.push((other_idx, 0, true));
                    } else {
                        // Disagree case: higher-BQ side keeps its
                        // BQ scaled by 0.8 (samtools' "we trust
                        // this less" haircut); loser zeroed.
                        let winner_idx = pick_disagree_winner(contributors, a, b);
                        let loser_idx = if winner_idx == a { b } else { a };
                        let scaled_bq =
                            scale_bq_by_0_8(contributors[winner_idx].bq_baq_at_walker_pos);
                        bq_updates.push((winner_idx, scaled_bq, false));
                        bq_updates.push((loser_idx, 0, true));
                    }
                }
            }
        }
    }

    // Apply bq updates in place. The fold honours both
    // `bq_zero_in_window` (zeroing every window event from this
    // contributor's cursor) and `bq_override_at_walker_pos`
    // (rewriting walker_pos events' BQ on top of the cursor
    // pull). Update the local contribution's `bq_baq_at_walker_pos`
    // and `events_at_pos` for consistency with the override.
    for (idx, new_bq, zero_in_window) in bq_updates.drain(..) {
        contributors[idx].bq_baq_at_walker_pos = new_bq;
        for ev in contributors[idx].events_at_pos.iter_mut() {
            set_match_event_bq(ev, new_bq);
        }
        if zero_in_window {
            contributors[idx].bq_zero_in_window = true;
        } else {
            contributors[idx].bq_override_at_walker_pos = Some(new_bq);
        }
    }

    // Drop indel-overlap losers from the contributor list, **in place**: the list arrives
    // in ascending `read_id` from the active set, and everything downstream of here — the
    // fold order, and with it the fold table's shape — is built on it still being so.
    // `swap_remove` would put the last contributor in the loser's place and undo that for
    // the whole column. The shift costs a memmove of the tail, and only on a column that
    // has a mate overlap *with an indel on one side*; `swap_remove` charged nothing there
    // and charged the ordering everywhere.
    //
    // Removal is applied back to front so earlier indices stay valid either way.
    to_remove.sort_unstable();
    to_remove.dedup();
    for idx in to_remove.drain(..).rev() {
        contributors.remove(idx);
    }
}

/// Loser-selection for the indel-overlap case. BQ first, then
/// first-of-pair, then `alignment_start`. Matches the pre-S7
/// semantics (which combined match-only and indel paths behind
/// the same loser-selection).
fn pick_overlap_loser(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let bq_a = contributors[a].bq_baq_at_walker_pos;
    let bq_b = contributors[b].bq_baq_at_walker_pos;
    match bq_a.cmp(&bq_b) {
        std::cmp::Ordering::Less => a,
        std::cmp::Ordering::Greater => b,
        std::cmp::Ordering::Equal => {
            let a_first = contributors[a].mate_role.is_first_of_pair();
            let b_first = contributors[b].mate_role.is_first_of_pair();
            match (a_first, b_first) {
                (true, false) => b,
                (false, true) => a,
                _ => {
                    if contributors[a].alignment_start <= contributors[b].alignment_start {
                        b
                    } else {
                        a
                    }
                }
            }
        }
    }
}

/// Keeper-selection for the agree case (S7). The choice is
/// statistically irrelevant — the surviving side carries the
/// summed BQ regardless — but must be deterministic. samtools
/// uses a qname hash; we mirror our existing tie-break logic
/// (first-of-pair, then `alignment_start`).
fn pick_agree_keeper(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let a_first = contributors[a].mate_role.is_first_of_pair();
    let b_first = contributors[b].mate_role.is_first_of_pair();
    match (a_first, b_first) {
        (true, false) => a,
        (false, true) => b,
        _ => {
            if contributors[a].alignment_start <= contributors[b].alignment_start {
                a
            } else {
                b
            }
        }
    }
}

/// Winner-selection for the disagree case (S7): higher BQ wins;
/// ties fall back to first-of-pair, then `alignment_start`
/// (samtools uses a qname hash on ties).
fn pick_disagree_winner(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let bq_a = contributors[a].bq_baq_at_walker_pos;
    let bq_b = contributors[b].bq_baq_at_walker_pos;
    match bq_a.cmp(&bq_b) {
        std::cmp::Ordering::Greater => a,
        std::cmp::Ordering::Less => b,
        std::cmp::Ordering::Equal => pick_agree_keeper(contributors, a, b),
    }
}

/// Extract the base from the `Match` event in `events_at_pos`.
/// In a match-only mate-overlap (no indel anchored at walker_pos
/// on either side), each contributor has exactly one Match event
/// at walker_pos.
fn match_base_at_pos(c: &ReadContribution) -> Option<u8> {
    c.events_at_pos.iter().find_map(|e| match e {
        ReadEvent::Match { base, .. } => Some(*base),
        ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. } => None,
    })
}

/// `min(a + b, 200)` in u8 space without overflow. Cap from
/// samtools (`tweak_overlap_quality` in
/// [`htslib/sam.c:5919-5921`](../../../htslib/sam.c#L5919-L5921)) —
/// quality values above ~Q200 are effectively meaningless.
fn sum_bq_capped_at_200(a: u8, b: u8) -> u8 {
    let sum = (a as u16) + (b as u16);
    sum.min(200) as u8
}

/// `(bq * 0.8)` truncated to u8, matching samtools' C `0.8 *
/// uint8_t` cast at
/// [`htslib/sam.c:5927`](../../../htslib/sam.c#L5927) (truncation,
/// not rounding).
fn scale_bq_by_0_8(bq: u8) -> u8 {
    (bq as f64 * 0.8) as u8
}

/// In-place BQ rewrite on a `Match` event. No-op on indel
/// events — the S7 BQ math only applies to match-only overlaps.
fn set_match_event_bq(ev: &mut ReadEvent, bq: u8) {
    if let ReadEvent::Match { bq_baq, .. } = ev {
        *bq_baq = bq;
    }
}

/// True iff at least one of the two contributors has an
/// Insertion or Deletion anchored at the current walker_pos.
fn pair_has_indel(a: &ReadContribution, b: &ReadContribution) -> bool {
    let has_indel = |c: &ReadContribution| {
        c.events_at_pos
            .iter()
            .any(|e| matches!(e, ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. },))
    };
    has_indel(a) || has_indel(b)
}

/// **ng's** — cut `contributors` down to `cap` by keeping the `cap` reads with the
/// smallest [`sampling_key`](read_sampling::sampling_key), and put the read ids of
/// everything dropped into `dropped`.
///
/// # Why it is a select and not a sort
///
/// `select_nth_unstable` partitions in linear time, which is what a rule applied at
/// every position of a deep region should cost. It leaves the kept side unordered,
/// so the kept indices are sorted afterwards — `cap` of them, at most 8,000 — and the
/// contributors compacted in index order. **Keeping the original relative order is not
/// cosmetic**: the fold's tie-breaks read the contributor list in the order the active
/// set holds it, so preserving that order is what confines this change to *which* reads
/// are kept and stops it from also changing how the kept ones are folded.
///
/// # The key comes from the active set, not from the contributor
///
/// A `ReadContribution` does not carry its read's query name, and putting the key on it
/// would cost eight bytes on every contributor at every position for a decision taken at
/// 909 positions in 251,792. So the key is looked up here, through the same secondary
/// index the fold uses, and only on the positions that cap.
///
/// A contributor whose read is somehow not in the active set — a state no path reaches,
/// since the contributor list is built from that very set two dozen lines earlier — is
/// given `u64::MAX`, i.e. dropped first. That is the answer that keeps the invariant
/// "the kept set is a function of the reads" true even if the impossible happens.
fn sample_to_cap(
    contributors: &mut Vec<ReadContribution>,
    cap: usize,
    active_reads: &ActiveReads,
    keys: &mut Vec<(u64, u32)>,
    kept: &mut Vec<u32>,
    dropped: &mut Vec<u32>,
) {
    keys.clear();
    keys.extend(contributors.iter().enumerate().map(|(index, contrib)| {
        let key = active_reads
            .get_by_read_id(contrib.read_id)
            .map_or(u64::MAX, |active| read_sampling::sampling_key(&active.read));
        (key, index as u32)
    }));
    // `(key, index)`: the index breaks a key tie deterministically, and it is the
    // contributor's own position in a list built by walking the active set in order — so
    // a tie is broken the same way a `sort_unstable` on keys alone could not promise.
    keys.select_nth_unstable(cap);
    dropped.extend(
        keys[cap..]
            .iter()
            .map(|&(_, index)| contributors[index as usize].read_id),
    );
    kept.clear();
    kept.extend(keys[..cap].iter().map(|&(_, index)| index));
    kept.sort_unstable();
    // Compaction in place. `kept` is ascending and its `rank`-th entry is at least
    // `rank`, so every swap moves an element forward into a slot already dealt with;
    // no kept contributor can be displaced past the prefix before its turn comes.
    for (rank, &index) in kept.iter().enumerate() {
        contributors.swap(rank, index as usize);
    }
    contributors.truncate(cap);
}

/// Per-column depth cap. Returns the lower indel cap if any
/// contributor reports an Insertion or Deletion at this anchor;
/// otherwise the SNP/REF cap.
fn column_depth_cap(contributors: &[ReadContribution], config: &WalkerConfig) -> usize {
    let any_indel = contributors.iter().any(|c| {
        c.events_at_pos
            .iter()
            .any(|e| matches!(e, ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. },))
    });
    if any_indel {
        config.max_indel_column_depth as usize
    } else {
        config.max_snp_column_depth as usize
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::cigar_cursor::EventsAt;
    use super::*;
    use crate::ng::locus_generation::pileup::tests::{Locus, MockFasta, snp_read};
    use crate::ng::types::{ContigId, Position};

    /// **The cap's discarded read ids reach the record** — the plumbing B3 adds, end to
    /// end through `run` rather than through `process_position`.
    ///
    /// Worth its own test because the per-record fixtures drive the fold directly and hand
    /// it the truncated ids themselves, so they cannot see the walk failing to *collect*
    /// them. Deleting the collection leaves every one of those fixtures green.
    ///
    /// Six reads at one position with a cap of two: four are truncated, none of them folds
    /// anywhere, and the record they would have reached says so.
    #[test]
    fn reads_the_column_cap_removed_are_reported_on_the_record() {
        use crate::ng::locus_generation::pileup::tests::{MockFasta, snp_read};

        let config = WalkerConfig {
            max_snp_column_depth: 2,
            ..WalkerConfig::default()
        };
        let reads: Vec<_> = (0..6)
            .map(|index| snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]))
            .collect();
        let loci: Vec<_> = run(reads, &MockFasta::new("ACG"), &config)
            .map(|item| item.expect("the walk succeeds"))
            .collect();

        assert!(!loci.is_empty(), "the fixture must emit loci");
        // **The count, per locus — not `> 0` summed.** Summing and asserting non-zero
        // survives an off-by-one in the slice the walk collects (`contributors[cap + 1..]`
        // reports three where there are four and still passes), which leaves the very
        // number B3 exists to get right unpinned on every record that left the walker.
        for locus in &loci {
            assert_eq!(
                locus.reads_discarded_by_cap, 4,
                "six reads at a cap of two truncates four, and none of them folded \
                 anywhere: {locus:?}"
            );
        }
        for locus in &loci {
            let folded: u32 = locus
                .observations
                .iter()
                .map(|observation| observation.num_obs)
                .sum();
            assert_eq!(
                folded, 2,
                "each column folds exactly the cap's worth, or the fixture is not capping"
            );
        }
    }

    // ------------------------------------------------------------------
    // The depth cap that acts at the door
    // ------------------------------------------------------------------

    /// **A region deeper than the walk will hold is walked, not failed.**
    ///
    /// This is the whole change, at its smallest: six reads at one position against a
    /// ceiling of two. Before it, the seventh line of this test was
    /// `WalkerError::ActiveReadsExhausted` and there were no loci at all — the walk of a
    /// real 100×-coverage tomato sample died this way 33 Mb into its first chromosome,
    /// where 4,143 reads pass the filters at a single base.
    ///
    /// Every read is accounted for on purpose. A read the ceiling removes must be absent
    /// from the fold *and* present in one of the two ceiling counters: dropping it from
    /// both would leave a walk that silently loses reads and reports a clean run, which
    /// is the failure these counters exist to make impossible.
    ///
    /// **Which two the ceiling keeps changed on 2026-08-05** — see
    /// [`the_ceiling_keeps_the_smallest_sampling_keys`] for that half. Here the claim is
    /// only that the walk survives and that the arithmetic closes.
    #[test]
    fn reads_past_the_active_read_cap_are_shed_and_the_walk_survives() {
        let config = WalkerConfig {
            max_active_reads: 2,
            ..WalkerConfig::default()
        };
        let reads: Vec<_> = (0..6)
            .map(|index| snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]))
            .collect();
        let reference = MockFasta::new("ACG");
        let mut walker = run(reads, &reference, &config);
        let loci: Vec<_> = walker
            .by_ref()
            .map(|item| item.expect("the walk must survive a region deeper than its cap"))
            .collect();
        let summary = walker.summary();

        assert_eq!(
            summary.reads_admitted + summary.reads_shed_at_admission,
            6,
            "every read the walk saw is either admitted or refused, and none is lost: \
             {summary:?}"
        );
        assert_eq!(
            summary.reads_evicted_at_ceiling,
            summary.reads_admitted - 2,
            "a read admitted beyond the ceiling's worth is one the ceiling gave back, so \
             the evictions are exactly the admissions past the two the set can hold"
        );
        assert!(!loci.is_empty(), "the region still produces loci");
        for locus in &loci {
            let folded: u32 = locus
                .observations
                .iter()
                .map(|observation| observation.num_obs)
                .sum();
            assert_eq!(
                folded, 2,
                "the evidence is the reads the set was holding when the column was folded, \
                 and only those — every one of these reads arrives before the first \
                 position is processed, so an evicted read folds nowhere"
            );
        }
    }

    /// **The ceiling keeps a fair subsample, not the reads that arrived first.**
    ///
    /// Six reads at one position against a ceiling of two. Which two survive is decided
    /// by [`read_sampling::sampling_key`] — a hash of the query name — so this test
    /// computes the same two names independently and asserts the walk kept them. Under
    /// the old rule it was always `r0` and `r1`, whatever the reads were, because reads
    /// arrive sorted by alignment start and the first two through the door won.
    ///
    /// The fixture is chosen so the expected pair is **not** `r0`/`r1`: an implementation
    /// that quietly went back to first-come would fail here rather than pass by
    /// coincidence, and the assertion says which pair it expected.
    #[test]
    fn the_ceiling_keeps_the_smallest_sampling_keys() {
        let config = WalkerConfig {
            max_active_reads: 2,
            ..WalkerConfig::default()
        };
        // **Each read gets a power-of-two MAPQ**, so the `mapq_sum` of the surviving
        // column names the surviving *set* and not merely its size — the one field on an
        // emitted observation that can distinguish six otherwise identical reads.
        let reads: Vec<_> = (0..6)
            .map(|index| {
                let mut read = snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]);
                read.mapq = 1 << index;
                read
            })
            .collect();

        // The two smallest keys, worked out from the reads alone.
        let mut ranked: Vec<_> = reads
            .iter()
            .map(|read| {
                (
                    read_sampling::sampling_key(read),
                    read.qname.to_string(),
                    read.mapq,
                )
            })
            .collect();
        ranked.sort();
        let expected_names: Vec<&str> = ranked[..2]
            .iter()
            .map(|(_, name, _)| name.as_str())
            .collect();
        let expected_mapq_sum: u64 = ranked[..2]
            .iter()
            .map(|(_, _, mapq)| u64::from(*mapq))
            .sum();
        assert!(
            !(expected_names.contains(&"r0") && expected_names.contains(&"r1")),
            "this fixture is only a test of the rule if the fair answer differs from the \
             first-come answer (r0, r1); it no longer does, so change the read names"
        );

        let reference = MockFasta::new("ACG");
        let mut walker = run(reads, &reference, &config);
        let loci: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let mapq_sum: u64 = loci[0]
            .observations
            .iter()
            .map(|observation| u64::from(observation.mapq_sum))
            .sum();

        assert_eq!(
            mapq_sum, expected_mapq_sum,
            "the ceiling kept a pair whose MAPQs sum to {mapq_sum}; the two smallest \
             sampling keys are {expected_names:?}, summing to {expected_mapq_sum}"
        );
    }

    /// **Nothing is shed below the cap** — the guard against a cap that fires early and
    /// silently subsamples ordinary data, which no output comparison would catch until
    /// the numbers had already moved.
    #[test]
    fn a_region_within_the_cap_sheds_nothing() {
        let reads: Vec<_> = (0..6)
            .map(|index| snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]))
            .collect();
        let reference = MockFasta::new("ACG");
        let mut walker = run(reads, &reference, &WalkerConfig::default());
        let _: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let summary = walker.summary();

        assert_eq!(summary.reads_admitted, 6);
        assert_eq!(summary.reads_shed_at_admission, 0);
    }

    // ------------------------------------------------------------------
    // The two properties the depth-cap change exists for
    // ------------------------------------------------------------------

    /// **No position is left short of the cap while reads covering it exist** — the
    /// owner's own test of the change, as a test.
    ///
    /// Forty reads over one base against a hold ceiling of eight and a SNP cap of eight.
    /// Every position must fold exactly eight, and `positions_short_of_cap` must be zero:
    /// the ceiling and the cap are the same number here, so the ceiling can shed as much
    /// as it likes and every position still reaches the cap.
    ///
    /// Then the same reads against a ceiling of four. Now the ceiling is *below* the cap,
    /// which is the configuration that leaves positions with less coverage than the input
    /// had for them, and the counter must say so — because a counter that reads zero
    /// whatever the configuration is not measuring anything.
    #[test]
    fn no_position_is_short_of_the_cap_while_reads_covering_it_exist() {
        let reference = MockFasta::new("ACG");
        let reads = || {
            (0..40)
                .map(|index| snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]))
                .collect::<Vec<_>>()
        };

        let at_the_cap = WalkerConfig {
            max_active_reads: 8,
            max_snp_column_depth: 8,
            ..WalkerConfig::default()
        };
        let mut walker = run(reads(), &reference, &at_the_cap);
        let loci: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let summary = walker.summary();
        assert_eq!(
            summary.positions_short_of_cap, 0,
            "the ceiling is at the cap, so no position can be short: {summary:?}"
        );
        assert_eq!(summary.short_of_cap_deficit, 0);
        for locus in &loci {
            let folded: u32 = locus
                .observations
                .iter()
                .map(|observation| observation.num_obs)
                .sum();
            assert_eq!(folded, 8, "every position folds the cap's worth");
        }

        let below_the_cap = WalkerConfig {
            max_active_reads: 4,
            ..at_the_cap
        };
        let mut walker = run(reads(), &reference, &below_the_cap);
        let _: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let summary = walker.summary();
        assert!(
            summary.positions_short_of_cap > 0,
            "with the ceiling below the cap every position is short, and a counter that \
             cannot say so is measuring nothing: {summary:?}"
        );
        assert!(
            summary.short_of_cap_deficit >= summary.positions_short_of_cap,
            "each short position is at least one read short: {summary:?}"
        );
    }

    /// **The kept set does not depend on the order the active set holds reads in.**
    ///
    /// The same reads at the same positions with the same names, offered to the walk in
    /// two different arrival orders within a position, must produce the same loci. Under
    /// the old `truncate` rule they did not: the cap kept a prefix of the active set's
    /// storage order, so re-ordering the input re-ordered the evidence.
    ///
    /// **Permuting the input is the strongest permutation available from outside**, and
    /// it is a genuine one: `admit_read` requires non-decreasing `(chrom, start)` and
    /// nothing more, so reads sharing a start may be offered in any order, and that order
    /// is exactly what the active set's iteration order was. Reversing it is the
    /// permutation the in-flight ordered-active-set work would otherwise impose.
    #[test]
    fn the_kept_set_does_not_depend_on_the_order_the_set_holds_reads_in() {
        let config = WalkerConfig {
            max_snp_column_depth: 5,
            ..WalkerConfig::default()
        };
        let reference = MockFasta::new("ACG");
        let bases: [&[u8]; 4] = [b"ACG", b"CCG", b"AGG", b"ACT"];
        let forwards: Vec<_> = (0..20)
            .map(|index| {
                let mut read = snp_read(
                    &format!("q{index}"),
                    1,
                    bases[index % bases.len()],
                    &[30; 3],
                );
                read.mapq = 20 + index as u8;
                read
            })
            .collect();
        let backwards: Vec<_> = forwards.iter().rev().cloned().collect();

        let mut walker = run(forwards, &reference, &config);
        let one: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let mut walker = run(backwards, &reference, &config);
        let other: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();

        assert!(!one.is_empty(), "the fixture must emit loci");
        // `q_sum` is a float sum whose order follows the arrival order, so it is compared
        // through the same tolerance the parity harness uses rather than exactly. Every
        // other field — the alleles, their counts, the strand and MAPQ moments — is a
        // property of the *set* that folded and must match bit for bit.
        assert_eq!(one.len(), other.len(), "the same loci, either way round");
        for (a, b) in one.iter().zip(other.iter()) {
            let key = |locus: &SampleLocusObservations| {
                let mut rows: Vec<_> = locus
                    .observations
                    .iter()
                    .map(|observation| {
                        (
                            observation.bases.to_vec(),
                            observation.num_obs,
                            observation.num_fwd,
                            observation.mapq_sum,
                            observation.mapq_sum_sq,
                            observation.placed_left,
                        )
                    })
                    .collect();
                rows.sort();
                (locus.region, locus.reference_bases.to_vec(), rows)
            };
            assert_eq!(
                key(a),
                key(b),
                "the same reads offered in two orders kept two different sets"
            );
        }
    }

    /// **The ceiling is on reads held open at once, not on reads seen.** A slot freed by
    /// a read the walker has passed is available to the next one, so a deep pile costs
    /// the reads inside it and nothing after it.
    ///
    /// Three reads at position 1 against a ceiling of two: one is refused. The fourth
    /// read starts at 100, long after the first three have expired, and is admitted —
    /// four reads seen, three admitted, one shed. A cap that counted admissions rather
    /// than residents would shed this last read too, and the depth of a whole
    /// chromosome would collapse after its first crowded base.
    #[test]
    fn a_slot_freed_by_an_expired_read_admits_the_next_one() {
        let config = WalkerConfig {
            max_active_reads: 2,
            ..WalkerConfig::default()
        };
        let reference = MockFasta::new(&("ACG".to_string() + &"T".repeat(96) + "ACG"));
        let mut reads: Vec<_> = (0..3)
            .map(|index| snp_read(&format!("r{index}"), 1, b"ACG", &[30; 3]))
            .collect();
        reads.push(snp_read("late", 100, b"ACG", &[30; 3]));

        let mut walker = run(reads, &reference, &config);
        let loci: Vec<_> = walker.by_ref().map(|item| item.expect("walks")).collect();
        let summary = walker.summary();

        assert_eq!(
            summary.reads_admitted + summary.reads_shed_at_admission,
            4,
            "four reads seen, each either admitted or refused"
        );
        // **The point of the test, and the one part the sampling rule cannot move.** The
        // late read is alone over positions 100..=102, so a locus there with one folded
        // read is proof it was let in — and it can only have been let in if the ceiling
        // counts residents rather than admissions.
        let late_locus = loci
            .iter()
            .find(|locus| locus.region.start.0 == 100)
            .expect("the late read must produce a locus of its own at 100");
        let folded: u32 = late_locus
            .observations
            .iter()
            .map(|observation| observation.num_obs)
            .sum();
        assert_eq!(
            folded, 1,
            "the late read starts long after the crowded base has cleared, so its slot is \
             free; a ceiling that counted admissions would have refused it and the depth \
             of a whole chromosome would collapse after its first crowded base"
        );
    }

    // ------------------------------------------------------------------
    // D1 — the walker pointed at one region after another
    // ------------------------------------------------------------------

    /// A [`RegionReadSource`] over a fixed list: for each region, the reads overlapping it,
    /// in position order. **A cursor with no file behind it**, and the one property that
    /// makes it a fair stand-in is that reads are *replayed* rather than consumed — the
    /// real cursor keeps every read it hands out and offers it to the next region that can
    /// use it (`spec/alignment_cursor.md` §6), which is what makes throwing the walker's
    /// look-ahead away safe.
    struct ScriptedRegionSource {
        reads: Vec<PreparedRead>,
        region: Option<GenomeRegion>,
        served: usize,
    }

    impl ScriptedRegionSource {
        fn new(reads: Vec<PreparedRead>) -> Self {
            Self {
                reads,
                region: None,
                served: 0,
            }
        }
    }

    impl Iterator for ScriptedRegionSource {
        type Item = PreparedRead;

        fn next(&mut self) -> Option<PreparedRead> {
            let region = self.region?;
            while let Some(read) = self.reads.get(self.served) {
                self.served += 1;
                let on_contig = read.chrom_id == region.contig.get();
                let overlaps = u64::from(read.alignment_start) <= region.end.get()
                    && u64::from(read.alignment_end) >= region.start.get();
                if on_contig && overlaps {
                    return Some(read.clone());
                }
            }
            None
        }
    }

    impl RegionReadSource for ScriptedRegionSource {
        type Error = std::convert::Infallible;

        fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), Self::Error> {
            self.region = Some(region);
            self.served = 0;
            Ok(())
        }
    }

    fn walk_region(
        walker: &mut PileupWalker<ScriptedRegionSource, &MockFasta>,
        region: GenomeRegion,
        stop_after: u32,
    ) -> Vec<Locus> {
        walker
            .move_to_region(region, stop_after)
            .expect("the scripted source cannot fail");
        walker
            .by_ref()
            .map(|item| item.expect("the walk succeeds"))
            .collect()
    }

    fn scripted_region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// The fixture the D1 tests share, over a 100-base all-`A` reference.
    ///
    /// **The lengths are not uniform, and that is the whole design.** Consecutive regions
    /// overlap, because each is asked for a halo past its end while the next one is asked
    /// from its own start — so a region serves reads an earlier one already admitted. `r8`
    /// runs 8..=27 and `r12` runs 12..=21: a walk of `1..=14` admits both and ends with
    /// `last_admitted_locus` at **12**, and the next region, from 15, is served `r8` at
    /// **8** first. That is the coordinate-order check in `admit_read` firing on ordinary
    /// forward progress, and it is what `begin_region` has to clear.
    ///
    /// **One read is silent at every position**, its adaptor boundary at its own first base
    /// — the shape `reads_silent_over_footprint` exists for. Without it the active set's
    /// `silent_exits` reads zero everywhere, and the difference between
    /// `ActiveReads::begin_region` (which zeroes that tally, as a per-region fold needs) and
    /// `ActiveReads::reset` (which preserves it, as a run total) is invisible.
    fn scripted_reads() -> Vec<PreparedRead> {
        let mut silent = snp_read("silent", 16, &[b'A'; 10], &[30; 10]);
        silent.adaptor_boundary = Some(silent.alignment_start);
        vec![
            snp_read("r1", 1, &[b'A'; 10], &[30; 10]),
            snp_read("r5", 5, &[b'A'; 10], &[30; 10]),
            snp_read("r8", 8, &[b'A'; 20], &[30; 20]),
            snp_read("r12", 12, &[b'A'; 10], &[30; 10]),
            silent,
            snp_read("r40", 40, &[b'A'; 10], &[30; 10]),
        ]
    }

    /// **A walker pointed at a second region must be indistinguishable from a fresh one
    /// pointed at the same region.** D1's whole contract, and the only oracle that covers
    /// every field `WalkerState::begin_region` has to decide about at once: a summary left
    /// un-cleared, an active read left behind, a stale look-ahead or a stale
    /// `last_admitted_locus` all show up here as a difference against the fresh walker.
    #[test]
    fn a_reused_walker_answers_a_region_exactly_as_a_fresh_one_does() {
        let reference = MockFasta::new(&"A".repeat(100));
        let config = WalkerConfig::default();
        // Halo-widened spans, as the generator asks for: the source is pointed 30 past each
        // region's end while the walk stops at the end itself.
        let regions = [(1u64, 14u64), (15, 30), (31, 60)];

        let mut reused = run(
            ScriptedRegionSource::new(scripted_reads()),
            &reference,
            &config,
        );

        for (start, end) in regions {
            let query = scripted_region(start, end + 30);
            let stop_after = end as u32;

            let from_reused = walk_region(&mut reused, query, stop_after);
            let reused_summary = reused.summary();

            let mut fresh = run(
                ScriptedRegionSource::new(scripted_reads()),
                &reference,
                &config,
            );
            let from_fresh = walk_region(&mut fresh, query, stop_after);

            // **Up to what the chain ids are called** — a reused walker carries its
            // allocator's counter forward from the previous region, so its ids start higher
            // than a fresh walker's while every other byte is the same (the owner's ruling
            // of 2026-08-17 put an id on every observation, which is what made the
            // difference visible).
            super::super::assert_same_evidence_up_to_chain_renaming(
                &from_reused,
                &from_fresh,
                &format!(
                    "region {start}..={end}: a reused walker emitted different loci from a \
                     fresh one"
                ),
            );
            // The two counters `fold_region_walk` sums region by region. A summary carried
            // across the boundary reads high here, and the caller would triangular-sum it.
            assert_eq!(
                (
                    reused_summary.reads_admitted,
                    reused_summary.records_emitted,
                    reused_summary.reads_silent_over_footprint,
                ),
                (
                    fresh.summary().reads_admitted,
                    fresh.summary().records_emitted,
                    fresh.summary().reads_silent_over_footprint,
                ),
                "region {start}..={end}: the reused walker's per-region counters differ \
                 from a fresh walker's",
            );
        }
    }

    /// **The chain-id allocator is the one thing that must *not* be restarted.**
    ///
    /// It is the run's, lent to the walker for a chromosome, and
    /// `PileupGeneratorCounts::fold_region_walk` folds two of its counters as deltas
    /// against the value they held when the region opened. A `begin_region` that replaced
    /// it — `WalkerState::new(config)` is one keystroke away and compiles — would zero
    /// them, and the deltas would collapse to nothing while `active_reads_high_water`
    /// survived as a max, which is what would make the corruption look selective.
    #[test]
    fn the_chain_id_allocators_counters_survive_a_region_boundary() {
        let reference = MockFasta::new(&"A".repeat(100));
        let config = WalkerConfig::default();
        let mut walker = run(
            ScriptedRegionSource::new(scripted_reads()),
            &reference,
            &config,
        );

        walk_region(&mut walker, scripted_region(1, 44), 14);
        let after_first = walker.chain_id_counters().chain_allocations;
        walk_region(&mut walker, scripted_region(31, 90), 60);
        let after_second = walker.chain_id_counters().chain_allocations;

        assert!(
            after_first > 0,
            "the first region must allocate something, or this test cannot fail"
        );
        assert!(
            after_second > after_first,
            "the second region allocated ids too, so the run-to-date total must have \
             grown: {after_first} then {after_second}",
        );
    }

    /// **A region walk admits reads far past the next region's start, and the next region
    /// must still be walkable.** The generator points the source at a halo past each
    /// region's end, so this is not an edge case — it is every region boundary on real
    /// data. Carried across, `last_admitted_locus` makes the next region's first read look
    /// like a read going backwards and the walk fails with `OutOfOrder`.
    #[test]
    fn a_region_that_admitted_reads_past_the_next_regions_start_still_walks_it() {
        let reference = MockFasta::new(&"A".repeat(100));
        let config = WalkerConfig::default();
        let mut walker = run(
            ScriptedRegionSource::new(scripted_reads()),
            &reference,
            &config,
        );

        // Segment `1..=14`, asked for its halo out to 44. It admits `r1`, `r5`, `r8` and
        // `r12`, so it ends with `last_admitted_locus` at 12.
        let first = walk_region(&mut walker, scripted_region(1, 44), 14);
        assert!(!first.is_empty(), "the first region must emit something");
        assert_eq!(
            walker.summary().reads_admitted,
            4,
            "the first region must admit through `r12` at 12, or the next region's `r8` at \
             8 is not a step backwards and this test cannot fail",
        );

        // Segment `15..=30`. `r8` runs 8..=27, so it overlaps and is served **first** — at
        // position 8, four positions before the last read the previous region admitted.
        walker
            .move_to_region(scripted_region(15, 60), 30)
            .expect("the scripted source cannot fail");
        let second: Vec<_> = walker.by_ref().collect();

        assert!(
            second.iter().all(|item| item.is_ok()),
            "the second region must walk without the first one's reads making its own look \
             out of order: {:?}",
            second.iter().find(|item| item.is_err()),
        );
        assert!(!second.is_empty(), "the second region must emit something");
    }

    /// **A region abandoned half-walked must leak no record into the next region** — the
    /// `open_records` half of `begin_region`, which the review found deletable with the
    /// whole suite green.
    ///
    /// The failure it hides is not a stale counter, which is what makes it worth its own
    /// test: a record left open is *finalised by the next region's walk*, the moment
    /// `close_aged_records_into` passes its footprint. So the next region leads its output
    /// with a locus at a coordinate nobody asked about — and with reads folded into it under
    /// `read_id`s that `ActiveReads::begin_region` has since restarted at zero.
    ///
    /// Abandoning is not exotic: `begin_segment` on a half-drained walk does exactly this,
    /// and the generator's own tests cover that path.
    #[test]
    fn a_region_abandoned_half_walked_leaks_no_record_into_the_next_one() {
        let reference = MockFasta::new(&"A".repeat(100));
        let config = WalkerConfig::default();
        let mut walker = run(
            ScriptedRegionSource::new(scripted_reads()),
            &reference,
            &config,
        );

        // Take **one** locus of the first region and walk away, leaving its later positions
        // open in the table.
        walker
            .move_to_region(scripted_region(1, 44), 14)
            .expect("the scripted source cannot fail");
        let abandoned = walker
            .next()
            .expect("the first region emits at least one locus")
            .expect("the walk succeeds");
        assert_eq!(
            abandoned.region.start,
            Position(1),
            "the fixture's first locus anchors at 1, so anything the next region leads with \
             below its own start came from here",
        );

        // A region well clear of the first, so nothing it emits can legitimately anchor
        // before position 40.
        let second = walk_region(&mut walker, scripted_region(40, 90), 60);

        assert!(!second.is_empty(), "the second region must emit something");
        for locus in &second {
            assert!(
                locus.region.start >= Position(40),
                "the second region emitted a locus at {:?}, which is the abandoned \
                 region's record finalised by this region's walk",
                locus.region.start,
            );
        }
    }

    fn contribution(
        bq: u8,
        is_first_mate: bool,
        alignment_start: u32,
        events: EventsAt,
    ) -> ReadContribution {
        ReadContribution {
            read_id: 0,
            active_index: 0,
            chain_id: 0,
            events_at_pos: events,
            bq_baq_at_walker_pos: bq,
            alignment_start,
            mate_role: if is_first_mate {
                super::super::MateRole::FirstOfPair
            } else {
                super::super::MateRole::SecondOfPair
            },
            bq_zero_in_window: false,
            bq_override_at_walker_pos: None,
        }
    }

    fn match_evs(pos: u32, base: u8, bq: u8) -> EventsAt {
        let mut v = EventsAt::new();
        v.push(ReadEvent::Match {
            ref_pos: pos,
            base,
            bq_baq: bq,
        });
        v
    }

    fn indel_ins_evs(anchor: u32, bq: u8) -> EventsAt {
        let mut v = EventsAt::new();
        v.push(ReadEvent::Match {
            ref_pos: anchor,
            base: b'A',
            bq_baq: bq,
        });
        v.push(ReadEvent::Insertion {
            anchor_ref_pos: anchor,
            seq: b"A".to_vec(),
            bq_proxy: bq,
        });
        v
    }

    // --- M19: pick_* tertiary tie-break tests --------------------
    //
    // All three pick functions fall through to comparing
    // `alignment_start` when both contributors' first-of-pair bits
    // agree. A flipped comparison here would silently change
    // tie-break determinism — exactly the kind of bug that surfaces
    // as "different VCF on the same input".

    #[test]
    fn pick_agree_keeper_breaks_remaining_tie_by_earlier_alignment_start() {
        // Both first-mate, BQ tie (no BQ check in this function).
        // Earlier alignment_start (50) wins over later (100).
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_agree_keeper(&c, 0, 1), 1);
        // Swap order: now index 0 has the earlier alignment_start.
        let c = vec![
            contribution(30, true, 50, match_evs(1, b'A', 30)),
            contribution(30, true, 100, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_agree_keeper(&c, 0, 1), 0);
    }

    #[test]
    fn pick_overlap_loser_breaks_bq_and_first_mate_tie_by_alignment_start() {
        // BQ tie + first-mate tie → the loser is the one with the
        // larger alignment_start.
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        // a=0 (start 100), b=1 (start 50) → loser is a (later start).
        assert_eq!(pick_overlap_loser(&c, 0, 1), 0);
    }

    #[test]
    fn pick_disagree_winner_on_bq_tie_delegates_to_pick_agree_keeper() {
        // BQ tie + first-mate tie → falls back to alignment_start.
        // The winner under the agree-keeper rule is the earlier
        // alignment_start.
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_disagree_winner(&c, 0, 1), 1);
    }

    // --- Boundary tests for the samtools-C parity helpers --------

    #[test]
    fn sum_bq_capped_at_200_caps_exactly_at_200() {
        assert_eq!(sum_bq_capped_at_200(0, 0), 0);
        assert_eq!(sum_bq_capped_at_200(100, 100), 200);
        assert_eq!(sum_bq_capped_at_200(150, 100), 200);
        assert_eq!(sum_bq_capped_at_200(255, 255), 200);
        assert_eq!(sum_bq_capped_at_200(99, 100), 199);
    }

    #[test]
    fn scale_bq_by_0_8_truncates_not_rounds() {
        assert_eq!(scale_bq_by_0_8(0), 0);
        assert_eq!(scale_bq_by_0_8(5), 4); // 4.0 exact
        assert_eq!(scale_bq_by_0_8(7), 5); // 5.6 → trunc 5 (round would give 6)
        assert_eq!(scale_bq_by_0_8(30), 24); // 24.0
    }

    // --- column_depth_cap: any-indel rule ------------------------

    #[test]
    fn column_depth_cap_returns_indel_cap_when_only_some_contributors_have_indel() {
        // Mixed SNP + one indel contributor at the same anchor.
        // The "any" rule must flip the column to the indel cap.
        let cfg = WalkerConfig {
            max_snp_column_depth: 8000,
            max_indel_column_depth: 2,
            ..WalkerConfig::default()
        };
        let v = vec![
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, indel_ins_evs(1, 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
        ];
        assert_eq!(column_depth_cap(&v, &cfg), 2);
    }
}
