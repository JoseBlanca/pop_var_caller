//! The generic locus generator itself — the knobs it takes, the run-level counts
//! it keeps, and the state it holds across segments.
//!
//! **ng's own file, not a copy.** Everything beside it in this module was
//! transcribed from `src/pileup/walker/`; this is the part production has no
//! counterpart for, because production drives its walker from a per-chromosome
//! pipeline stage rather than from a region-at-a-time generator
//! (`doc/devel/ng/arch/locus_generation_pileup.md` §1.1, §2.1).
//!
//! C1 lands the types and the one invariant that cannot be inherited: the
//! ceiling on [`PileupGeneratorConfig::max_record_span`]. The walk over a region
//! ([`begin_segment`](super::super::LocusGenerator::begin_segment) and the halo)
//! is C2's, and the allocator's lifetime across segments is C3's.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use crate::ng::locus_generation::{LocusGenerationError, SampleLocusObservations};
use crate::ng::read::input::{SampleReads, SampleRegionReads};
use crate::ng::read::{PreparedRead, ReadPrepError, ReadPreparer};
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::types::{GenomeRegion, Position};

use super::chain_id_allocator::{ChainIdAllocator, ChainIdAllocatorCounters};
use super::genome_walk::{PileupWalker, RunSummary};
use super::{DEFAULT_MAX_ACTIVE_READS, WalkerConfig};

/// The widest `max_record_span` this generator accepts: 65,535 reference
/// positions, the widest footprint a [`ReadCoverage`] run can describe.
///
/// [`ReadCoverage::Observed`] carries `offset_in_locus` and `positions_covered`
/// as `u16`s, minted through [`LocusLen::from_positions`], which **saturates**
/// rather than failing. A record footprint wider than this therefore makes a
/// partial witness report a *truncated* `positions_covered` — a wrong number, no
/// error, at exactly the long-deletion loci this generator exists to get right.
///
/// [`ReadCoverage`]: crate::ng::locus_generation::ReadCoverage
/// [`ReadCoverage::Observed`]: crate::ng::locus_generation::ReadCoverage::Observed
/// [`LocusLen::from_positions`]: crate::ng::locus_generation::LocusLen::from_positions
pub const MAX_RECORD_SPAN_CEILING: u32 = u16::MAX as u32;

/// **Production's default `max_record_span` is inside ng's ceiling**, which is
/// why [`PileupGeneratorConfig::default`] can inherit the knob by name without
/// shipping a configuration [`PileupGenerator::new`] would reject.
///
/// A compile-time check rather than a test, because it compares two constants:
/// should production ever raise its default past 65,535, this breaks the build
/// and the divergence becomes a decision rather than a runtime surprise.
const _: () = assert!(
    crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN <= MAX_RECORD_SPAN_CEILING,
    "production's default max_record_span no longer fits a ReadCoverage run: ng must either \
     widen the run or stop inheriting the default",
);

/// This generator's knobs — owned, taken at construction, and **production's
/// values, inherited and never measured by ng** (spec §7).
///
/// Raw `u32` rather than `Bp`, because the copied walker speaks production's
/// integer widths and a port must not change types under itself (spec §3). The
/// five are exactly production's [`WalkerConfig`] fields, reached from
/// production's own `pub const`s **by name** in [`Default`], so there is one
/// source of truth until ng deliberately diverges and the divergence shows up as
/// a diff.
///
/// # One knob is not simply production's
///
/// [`max_record_span`](Self::max_record_span) is capped at
/// [`MAX_RECORD_SPAN_CEILING`] and [`check`](Self::check) rejects more.
/// Production's `--max-record-span` is an unbounded `u32`; inheriting the knob by
/// name would inherit a silent truncation ng's coverage runs cannot survive. The
/// cap is not a constraint in practice — production's own default of 5,000 is
/// already unreachable with Illumina reads, and the ceiling is thirteen times
/// that again (owner, 2026-07-29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupGeneratorConfig {
    /// Reads folded at a position with no indel anchored there. Production: 8,000.
    pub max_snp_column_depth: u32,
    /// Reads folded at a position where any read has an indel. Production: 250.
    pub max_indel_column_depth: u32,
    /// Widest record footprint before the walk fails. Production: 5,000; ng
    /// additionally rejects anything above [`MAX_RECORD_SPAN_CEILING`].
    pub max_record_span: u32,
    /// How far a first mate stays available for pairing. Production: 10,000.
    pub mate_lookup_window: u32,
    /// Active-read ceiling. Production: 4,096.
    pub max_active_reads: u32,
}

impl Default for PileupGeneratorConfig {
    /// Production's five constants, **by name** — never by literal, so a retune
    /// on production's side reaches ng as a diff rather than as drift.
    fn default() -> Self {
        Self {
            max_snp_column_depth: crate::pileup::walker::DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: crate::pileup::walker::DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW,
            // ng's copy of the constant, which `walker_vocabulary_tests` pins equal
            // to production's — the one `DEFAULT_*` the verbatim copy forked, because
            // it is declared inside `chain_id_allocator.rs`.
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        }
    }
}

impl PileupGeneratorConfig {
    /// Reject a configuration a [`ReadCoverage`] run could not describe.
    ///
    /// Called by [`PileupGenerator::new`], so a bad knob never reaches a locus.
    /// `coverage_of` carries a `debug_assert` stating the same envelope; this is
    /// what makes it provable rather than hopeful, in release builds too.
    ///
    /// [`ReadCoverage`]: crate::ng::locus_generation::ReadCoverage
    pub fn check(&self) -> Result<(), PileupGeneratorConfigError> {
        if self.max_record_span > MAX_RECORD_SPAN_CEILING {
            return Err(PileupGeneratorConfigError::RecordSpanExceedsCoverageRun {
                max_record_span: self.max_record_span,
                ceiling: MAX_RECORD_SPAN_CEILING,
            });
        }
        Ok(())
    }

    /// The same five knobs in the shape the copied walker reads them.
    ///
    /// Written as an exhaustive struct literal rather than with `..default()`:
    /// a knob production adds to [`WalkerConfig`] is then a compile error here,
    /// which is a decision to make rather than a default to inherit silently.
    pub fn to_walker_config(&self) -> WalkerConfig {
        WalkerConfig {
            max_snp_column_depth: self.max_snp_column_depth,
            max_indel_column_depth: self.max_indel_column_depth,
            max_record_span: self.max_record_span,
            mate_lookup_window: self.mate_lookup_window,
            max_active_reads: self.max_active_reads,
        }
    }
}

/// A configuration this generator cannot walk under. `#[non_exhaustive]`; raised
/// at construction, never mid-walk.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PileupGeneratorConfigError {
    /// `max_record_span` is wider than a coverage run can describe, so a partial
    /// witness of a maximally-wide record would report a saturated
    /// `positions_covered` — silently, since [`LocusLen::from_positions`] clamps
    /// rather than failing.
    ///
    /// [`LocusLen::from_positions`]: crate::ng::locus_generation::LocusLen::from_positions
    #[error(
        "max_record_span ({max_record_span}) exceeds the widest footprint a coverage run can \
         describe ({ceiling}): a partially-witnessed record that wide would report a truncated \
         positions_covered rather than an error"
    )]
    RecordSpanExceedsCoverageRun { max_record_span: u32, ceiling: u32 },
}

/// Run-level counts for this generator, kept alongside the shared
/// [`LocusCounts`](super::super::LocusCounts) the dispatcher owns.
///
/// The first seven mirror production's [`RunSummary`](super::RunSummary) field for
/// field (spec §7); the last two are ng's. **`records_emitted` is deliberately not
/// mirrored**: what a caller wants is loci *kept*, which is the walk's emissions
/// minus [`records_outside_region`](Self::records_outside_region), and the kept
/// count is already `LocusCounts::loci_emitted`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PileupGeneratorCounts {
    /// Reads the walk admitted, across every segment.
    pub reads_admitted: u64,
    /// Times an open record's footprint grew under a longer deletion.
    pub record_widen_events: u64,
    /// Positions where two overlapping mates were reconciled.
    pub mate_overlap_positions: u64,
    /// Chain ids allocated — fragments the walk gave an identity.
    pub chain_allocations: u64,
    /// The largest number of concurrently-active reads seen. A **max**, not a
    /// sum: regions are walked one at a time, so the run's peak is the largest
    /// single region's peak.
    pub active_reads_high_water: u32,
    /// First mates evicted from `pending_mates` unpaired, having waited longer
    /// than `mate_lookup_window`.
    pub mate_lookup_evictions: u64,
    /// Columns where the contributor list was truncated by the applicable
    /// per-column depth cap.
    pub column_depth_truncations: u64,
    /// Reads silent over a whole record footprint — every base `N` or
    /// adaptor-masked, so never a contributor at any position and invisible to
    /// the per-locus `reads_without_observation` tally (spec §6).
    ///
    /// **Nothing increments this yet.** It needs a per-active-read "ever
    /// contributed" flag in the walk, which no step through C4 adds; the dump
    /// tool's read-accounting assertion (spec §13, plan D2) is what forces it.
    /// Until then the field reads zero, which is not the same claim as "no read
    /// was silent".
    pub reads_silent_over_footprint: u64,
    /// Records the region clamp dropped, because their anchor fell outside the
    /// region being walked. Observably zero-sum across neighbouring regions,
    /// which is how the gap-free tiling argument stays checkable (spec §7).
    pub records_outside_region: u64,
}

impl PileupGeneratorCounts {
    /// Fold one region's walk into the run's totals.
    ///
    /// **Not [`RunSummary::merge`], and the difference is the trap spec §8
    /// names.** Production builds a fresh walker — and so a fresh allocator —
    /// per region, so `summary()` *assigns* the allocator's three counters and
    /// summing region summaries is right. ng shares one allocator across
    /// regions and `reset()` deliberately preserves those counters, so every
    /// region's summary reports the **run-to-date** total: adding them up
    /// triangular-sums `chain_allocations` and `mate_lookup_evictions`, roughly
    /// by the region count, in a plausible-looking `u64`. They are therefore
    /// folded as **deltas** against `chain_ids_at_open`, the allocator's
    /// counters as they stood when this walk took it.
    ///
    /// `active_reads_high_water` is a max and survives either treatment — which
    /// is what would make the corruption look selective enough to rationalise.
    fn fold_region_walk(
        &mut self,
        summary: &RunSummary,
        chain_ids_at_open: ChainIdAllocatorCounters,
    ) {
        // Exhaustive destructure (no `..`), production's own idiom for this
        // struct: a new `RunSummary` field is a compile error here until it is
        // explicitly folded, and the fold is not uniform — two fields are
        // deltas, one is a max, one is deliberately dropped.
        let RunSummary {
            reads_admitted,
            records_emitted,
            record_widen_events,
            mate_overlap_positions,
            chain_allocations,
            active_reads_high_water,
            mate_lookup_evictions,
            column_depth_truncations,
        } = *summary;
        debug_assert!(
            chain_allocations >= chain_ids_at_open.chain_allocations
                && mate_lookup_evictions >= chain_ids_at_open.mate_lookup_evictions,
            "the allocator's counters only grow; a baseline above the current value means \
             the walk was handed a different allocator than the one it gave back",
        );
        self.reads_admitted += reads_admitted;
        self.record_widen_events += record_widen_events;
        self.mate_overlap_positions += mate_overlap_positions;
        self.chain_allocations +=
            chain_allocations.saturating_sub(chain_ids_at_open.chain_allocations);
        self.active_reads_high_water = self.active_reads_high_water.max(active_reads_high_water);
        self.mate_lookup_evictions +=
            mate_lookup_evictions.saturating_sub(chain_ids_at_open.mate_lookup_evictions);
        self.column_depth_truncations += column_depth_truncations;
        // `records_emitted` is deliberately not mirrored: what a caller wants is
        // loci *kept*, which is this minus `records_outside_region` — see the
        // type's own doc.
        let _ = records_emitted;
    }
}

/// What the generator lends the per-segment read stream: the preparer, its
/// scratch, and the slot the stream latches a fatal error into.
///
/// **Shared because the stream lives inside the walk and the walk lives inside
/// the generator.** The walk owns its read stream (arch §2.2) and the stream
/// has to reach the preparer for every read; a borrow of a sibling field is not
/// expressible, so the three travel together behind one handle the generator
/// keeps and each walk clones. Keeping the *scratch* here rather than in the
/// stream is the point of the arrangement — a stream-owned scratch would be
/// reallocated at every region boundary, which is what a generator-level
/// scratch field exists to prevent (arch §1.1).
///
/// One `RefCell` borrow per read, held only across `prepare_read`; the walk is
/// single-threaded (arch §9) and nothing else touches the cell while a read is
/// being prepared.
struct ReadPreparation<P: ReadPreparer> {
    preparer: P,
    scratch: P::Scratch,
    /// The first fatal error the stream hit, kept because an
    /// `Iterator<Item = PreparedRead>` has nowhere to report one.
    ///
    /// The walker consumes an infallible item type, so a failed query or a
    /// failed preparation has to shed its error and report end-of-stream; the
    /// generator takes it back once the walker drains. Production hit the same
    /// seam and solved it the same way (`ErrorSheddingAdapter` in
    /// `pop_var_caller/cli/error_bridge.rs`) — ng keeps its own because it
    /// shed*s* two error types into one and because a locus generator has no
    /// business importing from the CLI.
    latched_error: Option<LocusGenerationError>,
}

/// One region's read stream, prepared: the query's `AlignedRead`s canonicalised
/// into the `PreparedRead`s the walk consumes, with fatal errors shed into the
/// shared [`ReadPreparation`] and reported as end-of-stream.
///
/// **Lazy, and it must stay lazy.** The query is a pull iterator that collects
/// nothing, so what is resident is the reads overlapping the walker's current
/// position — bounded by depth, not by the region's length. Collecting here to
/// make an ownership problem go away would turn a depth-shaped footprint into a
/// region-shaped one, on regions that run to hundreds of kilobases (spec §7).
struct PreparedRegionReads<R: RawRefSeq, P: ReadPreparer> {
    reads: SampleRegionReads<R>,
    preparation: Rc<RefCell<ReadPreparation<P>>>,
    /// Set once the query is exhausted or an error was latched. Both are
    /// terminal: a stream that shed an error must not resume, or the walk would
    /// carry on over a hole it was never told about.
    done: bool,
}

impl<R: RawRefSeq, P: ReadPreparer> Iterator for PreparedRegionReads<R, P> {
    type Item = PreparedRead;

    fn next(&mut self) -> Option<PreparedRead> {
        while !self.done {
            let read = match self.reads.next() {
                None => {
                    self.done = true;
                    return None;
                }
                Some(Ok(read)) => read,
                Some(Err(source)) => return self.shed(source.into()),
            };
            let mut preparation = self.preparation.borrow_mut();
            // Split the borrow: `prepare_read` takes `&self` and `&mut
            // Self::Scratch`, which are two fields of the same cell.
            let ReadPreparation {
                preparer, scratch, ..
            } = &mut *preparation;
            match preparer.prepare_read(read, scratch) {
                Ok(Some(prepared)) => return Some(prepared),
                // "No usable observation" — the preparer declined this read and
                // the run continues. **No v1 preparer returns it** (the only
                // step that could decline was BAQ, deferred), so there is no
                // tally for it yet; when a declining preparer lands, its count
                // belongs beside `reads_silent_over_footprint`.
                Ok(None) => continue,
                Err(source) => {
                    drop(preparation);
                    // Matched exhaustively rather than through a catch-all:
                    // `ReadPrepError` is `#[non_exhaustive]`, which binds other
                    // crates but not this one, so a preparation failure ng
                    // cannot yet describe is a compile error here instead of a
                    // reference error it is not.
                    let error = match source {
                        ReadPrepError::Reference(source) => LocusGenerationError::Reference(source),
                    };
                    return self.shed(error);
                }
            }
        }
        None
    }
}

impl<R: RawRefSeq, P: ReadPreparer> PreparedRegionReads<R, P> {
    /// Latch a fatal error and report end-of-stream. The **first** error wins:
    /// a later one describes a walk that was already over.
    fn shed(&mut self, error: LocusGenerationError) -> Option<PreparedRead> {
        self.done = true;
        let mut preparation = self.preparation.borrow_mut();
        if preparation.latched_error.is_none() {
            preparation.latched_error = Some(error);
        }
        None
    }
}

/// The walk over one region: **owns** its read stream, so nothing borrows the
/// `SampleReads` it was built from and nothing has to be materialised (arch
/// §2.2). One per region, built on the first `next_locus` and drained across the
/// calls that follow.
struct RegionWalk<R: RawRefSeq, P: ReadPreparer> {
    walker: PileupWalker<PreparedRegionReads<R, P>, Arc<R>>,
    /// The region records are clamped to — the walk emits records anchored just
    /// outside it, because the query returns every read *overlapping* the
    /// region and the halo widens that further still.
    region: GenomeRegion,
    /// The chain-id allocator's counters as they stood when this walk took it —
    /// the baseline its own contribution is measured against (C3).
    ///
    /// Snapshotted here rather than at `begin_segment`, which the plan's wording
    /// suggests, because `begin_segment` opens nothing: between it and the
    /// query, nothing touches the allocator, so the two moments hold the same
    /// numbers and this one cannot drift from the walk it belongs to.
    chain_ids_at_open: ChainIdAllocatorCounters,
}

/// ng's generic locus generator: a streaming pileup walk over one `Generic`
/// region, emitting one
/// [`SampleLocusObservations`](super::super::SampleLocusObservations) per covered
/// position.
///
/// # What it holds, and why it holds it
///
/// A generator owns its accessors as fields (`locus_generation.md` §2). **Two
/// reference accessors, not one:** the preparer carries its own (read
/// preparation's rule), and `reference` serves the walk's REF fetches. Neither
/// is rebuilt per segment — a fresh accessor per region throws away the sliding
/// buffer at every boundary and re-pays a `.fai` parse plus two `open(2)`s,
/// which is the ~564k-opens trap the STR side already paid for (spec §8). That
/// is why it is held behind an `Arc`: the walk must *own* a `RefSeq`, and an
/// `Arc` clone is the owned handle that does not rebuild the reader
/// ([`ref_seq`](crate::ng::ref_seq)'s shared-handle impls exist for this).
///
/// `make_reference` is a third accessor and deliberately a **factory**: the read
/// query gives each of a sample's k files its own raw accessor, because they are
/// stateful readers and sharing one would make k streams share a file cursor
/// ([`SampleReads::reads_in_region`]).
///
/// The chain-id allocator likewise lives here rather than inside the walk: ng
/// walks one region where production walks a chromosome, and a fresh allocator
/// per segment would give two fragments of different regions the same id (spec
/// §8). What that costs — `reset()` between segments, and counters that must be
/// folded as deltas because `reset()` preserves them — is C3's.
pub struct PileupGenerator<R: RawRefSeq, MakeReference, P: ReadPreparer>
where
    MakeReference: FnMut() -> R,
{
    /// The reference the walk fetches REF bases from. Built once, for the run,
    /// and handed to each walk as a shared handle.
    reference: Arc<R>,
    /// Builds a fresh read-query accessor for each file of the sample, which is
    /// what [`SampleReads::reads_in_region`]'s mismatch-fraction filter reads.
    make_reference: MakeReference,
    /// The preparer and its scratch, lent to each region's read stream.
    preparation: Rc<RefCell<ReadPreparation<P>>>,
    /// Lives across segments so `next_id` never repeats — **`None` exactly while
    /// a walk holds it**, which is the invariant `open_walk` and `end_walk`
    /// keep between them.
    ///
    /// An `Option` rather than a swap with a fresh allocator: a placeholder
    /// starting at zero is the state this whole arrangement exists to avoid, and
    /// it would fail silently. Absent, it fails loudly.
    chain_ids: Option<ChainIdAllocator>,
    config: PileupGeneratorConfig,
    counts: PileupGeneratorCounts,
    /// The region [`begin_segment`](Self::begin_segment) was given, before any
    /// query has been opened for it.
    current_region: Option<GenomeRegion>,
    /// The walk over `current_region`, opened by the first `next_locus`.
    walk: Option<RegionWalk<R, P>>,
}

impl<R: RawRefSeq, MakeReference, P: ReadPreparer> PileupGenerator<R, MakeReference, P>
where
    MakeReference: FnMut() -> R,
{
    /// Build a generator over `reference` (the walk's REF fetches),
    /// `make_reference` (the read query's per-file accessor factory) and
    /// `preparer` (per-read canonicalisation), with `config` checked before
    /// anything is held.
    ///
    /// Fails only on a configuration a coverage run could not describe — see
    /// [`PileupGeneratorConfig::check`].
    pub fn new(
        reference: Arc<R>,
        make_reference: MakeReference,
        preparer: P,
        config: PileupGeneratorConfig,
    ) -> Result<Self, PileupGeneratorConfigError> {
        config.check()?;
        Ok(Self {
            reference,
            make_reference,
            preparation: Rc::new(RefCell::new(ReadPreparation {
                preparer,
                scratch: P::Scratch::default(),
                latched_error: None,
            })),
            chain_ids: Some(ChainIdAllocator::with_caps(
                config.max_active_reads,
                config.mate_lookup_window,
            )),
            config,
            counts: PileupGeneratorCounts::default(),
            current_region: None,
            walk: None,
        })
    }

    /// The knobs this generator was built with.
    pub fn config(&self) -> &PileupGeneratorConfig {
        &self.config
    }

    /// The run-level counts accumulated across every segment walked so far.
    pub fn counts(&self) -> &PileupGeneratorCounts {
        &self.counts
    }

    /// Start a region: record it and **open nothing**.
    ///
    /// It cannot fail, and opening a read query can — so the query is opened by
    /// the first [`next_locus`](Self::next_locus), which is where an
    /// [`IngestError`](crate::ng::read::input::IngestError) surfaces (arch
    /// §2.1). Any unfinished walk of the previous region is ended here: the
    /// contract is one segment at a time, and a half-drained walk of a region
    /// nobody is asking about any more has nothing left to contribute — but it
    /// still holds the chain-id allocator, and **an abandoned walk is the one
    /// case where `active_count` does not fall back to zero on its own** (spec
    /// §8), so it is ended rather than dropped.
    pub fn begin_segment(&mut self, region: GenomeRegion) {
        self.end_walk();
        self.current_region = Some(region);
    }

    /// End the walk in flight, if any: fold its counters into the run's totals,
    /// take the chain-id allocator back, and `reset()` it for the next region.
    ///
    /// **The only place `walk` is cleared**, which is what keeps "the allocator
    /// is absent exactly while a walk holds it" true without anyone having to
    /// remember it. `reset()` clears `pending_mates` and `active_count` while
    /// preserving `next_id` — carried across regions blindly, a pending first
    /// mate from one region pairs with a read in another (the eviction test
    /// compares raw positions) and `active_count` climbs toward
    /// `ActiveReadsExhausted` (spec §8).
    fn end_walk(&mut self) {
        let Some(walk) = self.walk.take() else {
            return;
        };
        self.counts
            .fold_region_walk(&walk.walker.summary(), walk.chain_ids_at_open);
        let mut chain_ids = walk.walker.into_chain_ids();
        chain_ids.reset();
        self.chain_ids = Some(chain_ids);
    }

    /// The next locus of the region begun, or `None` once the walk drains.
    ///
    /// Records whose **anchor** falls outside the region are dropped and
    /// tallied in [`records_outside_region`](PileupGeneratorCounts::records_outside_region)
    /// — the rule that makes neighbouring regions tile without duplicates or
    /// holes, since typed regions tile the genome gap-free and disjointly
    /// (spec §2).
    pub fn next_locus(
        &mut self,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        let Some(region) = self.current_region else {
            return Ok(None);
        };
        if self.walk.is_none() {
            self.walk = Some(self.open_walk(region, reads)?);
        }
        loop {
            // PANIC-FREE: opened immediately above, and every path that clears
            // the walk returns rather than looping.
            let walk = self.walk.as_mut().expect("the walk was just opened");
            let clamp = walk.region;
            match walk.walker.next() {
                Some(Ok(locus)) => {
                    if clamp.contains(locus.region.start) {
                        return Ok(Some(locus));
                    }
                    self.counts.records_outside_region += 1;
                }
                Some(Err(source)) => {
                    // Fatal and terminal: the walker yields nothing after an
                    // error, so the walk ends here and gives back the allocator
                    // it was lent.
                    self.end_walk();
                    return Err(LocusGenerationError::Walker(source));
                }
                None => {
                    // The walk is over. A read stream that shed a fatal error
                    // also reported end-of-stream, so the walker drains
                    // normally and the error is only visible here — checked
                    // **before** `Ok(None)`, or a broken query would read as an
                    // empty region.
                    self.end_walk();
                    if let Some(error) = self.preparation.borrow_mut().latched_error.take() {
                        return Err(error);
                    }
                    return Ok(None);
                }
            }
        }
    }

    /// Open the query for `region` and build the walk over it.
    ///
    /// **The query is `[region.start, region.end + max_record_span]` — the
    /// halo.** A record anchored inside the region can have a footprint
    /// reaching `max_record_span` past its end (a long deletion does exactly
    /// that), and reads that fold into it may lie *entirely beyond* the region.
    /// A query for "reads overlapping the region" never returns them, and the
    /// record is then emitted by the right region with part of its support
    /// missing — **and no counter notices**, because the record itself is not
    /// lost (spec §2). The extra reads are walked and their records clamped
    /// away unless anchored in the region.
    ///
    /// The walk is **stopped** at the region's end rather than run to the end of
    /// the halo: see [`PileupWalker::stopping_after`].
    fn open_walk(
        &mut self,
        region: GenomeRegion,
        reads: &SampleReads,
    ) -> Result<RegionWalk<R, P>, LocusGenerationError> {
        let query = GenomeRegion {
            contig: region.contig,
            start: region.start,
            end: Position(
                region
                    .end
                    .get()
                    .saturating_add(u64::from(self.config.max_record_span)),
            ),
        };
        // Reborrowed into a local: the factory is a field, and `reads_in_region`
        // would otherwise hold a borrow of `self` across the call that stores
        // the walk back into `self`.
        let make_reference = &mut self.make_reference;
        let stream = reads.reads_in_region(query, make_reference)?;
        let prepared = PreparedRegionReads {
            reads: stream,
            preparation: Rc::clone(&self.preparation),
            done: false,
        };
        // Saturating, and in the safe direction: the walker's positions are
        // `u32`, so a region ending past 4 Gb cannot be walked at all, and a
        // stop bound clamped to `u32::MAX` makes the walk run long rather than
        // stop early.
        let stop_after = region.end.get().min(u64::from(u32::MAX)) as u32;
        // Taken **after** the query, which is the only fallible step: taken
        // before, an unqueryable region would carry the run's chain-id sequence
        // off with the error.
        let chain_ids = self
            .chain_ids
            .take()
            .expect("the allocator is absent only while a walk holds it, and no walk is open");
        let chain_ids_at_open = chain_ids.counters();
        let walker = super::genome_walk::run(
            prepared,
            Arc::clone(&self.reference),
            &self.config.to_walker_config(),
        )
        .stopping_after(stop_after)
        .adopting_chain_ids(chain_ids);
        Ok(RegionWalk {
            walker,
            region,
            chain_ids_at_open,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use noodles_sam::alignment::RecordBuf;

    use crate::ng::locus_generation::LocusLen;
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, fixture_reference, header, indexed_bam, matching_contigs,
        read_named_with_length,
    };
    use crate::ng::read::{AlignedRead, PreparedRead, ReadPrepError};
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::ContigId;
    use crate::pileup_record::ChainId;

    /// Production's `--no-baq` build, re-attached to ng's read type — the read
    /// preparation these tests want out of the way, so what they measure is the
    /// walk. Mirrors `ng::read`'s own `prepare_via_passthrough`.
    struct PassthroughPreparer;

    impl ReadPreparer for PassthroughPreparer {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: AlignedRead,
            _scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            let read_group = read.read_group;
            let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
            Ok(Some(PreparedRead::from_production(
                crate::pileup::per_sample::baq_engine::prepare_passthrough(
                    read.into_mapped_read(),
                    chrom_id,
                ),
                read_group,
            )))
        }
    }

    /// A preparer that fails on every read — the fatal-error path, which the
    /// stream has to shed and the generator has to report.
    struct FailsToPrepare;

    impl ReadPreparer for FailsToPrepare {
        type Scratch = ();
        fn prepare_read(
            &self,
            _read: AlignedRead,
            _scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            Err(ReadPrepError::Reference(
                crate::ng::ref_seq::RefSeqError::UnknownContig(ContigId(7)),
            ))
        }
    }

    /// The reference the walk fetches REF bases from: the fixture contigs, all
    /// `A`, which is what `build_fasta` writes and what the fixture reads carry —
    /// so every locus is REF-only and the tests are about *which* loci exist and
    /// what supports them.
    fn fixture_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_named_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(name, length)| ((*name).to_string(), vec![b'A'; *length]))
                .collect(),
        )
    }

    fn a_generator_with<P: ReadPreparer>(
        config: PileupGeneratorConfig,
        preparer: P,
    ) -> Result<
        PileupGenerator<InMemoryRefSeq, impl FnMut() -> InMemoryRefSeq, P>,
        PileupGeneratorConfigError,
    > {
        PileupGenerator::new(Arc::new(fixture_bases()), fixture_bases, preparer, config)
    }

    fn a_generator(
        config: PileupGeneratorConfig,
    ) -> Result<
        PileupGenerator<InMemoryRefSeq, impl FnMut() -> InMemoryRefSeq, PassthroughPreparer>,
        PileupGeneratorConfigError,
    > {
        a_generator_with(config, PassthroughPreparer)
    }

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// Open a `SampleReads` over one indexed BAM of `records`, against the
    /// fixture reference. The two `TempDir`s must outlive the sample.
    fn sample_reads_with(
        records: &[RecordBuf],
    ) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::ng::read::filtering::ReadFilterConfig;
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// Drain a whole segment: every locus `next_locus` yields for `region`.
    fn loci_of<R: RawRefSeq, MF: FnMut() -> R, P: ReadPreparer>(
        generator: &mut PileupGenerator<R, MF, P>,
        region: GenomeRegion,
        reads: &SampleReads,
    ) -> Vec<SampleLocusObservations> {
        generator.begin_segment(region);
        let mut loci = Vec::new();
        while let Some(locus) = generator.next_locus(reads).expect("the walk succeeds") {
            loci.push(locus);
        }
        loci
    }

    /// The 1-based anchor of each locus, in emission order.
    fn anchors(loci: &[SampleLocusObservations]) -> Vec<u64> {
        loci.iter().map(|locus| locus.region.start.get()).collect()
    }

    /// Total observations supporting a locus, across every row.
    fn total_obs(locus: &SampleLocusObservations) -> u32 {
        locus
            .observed_sequences
            .iter()
            .map(|observation| observation.num_obs)
            .sum()
    }

    /// A 30 bp read carrying exactly one mismatch, at `mismatch_offset` bases
    /// into it.
    ///
    /// The fixture reference is all `A`, and a read that agrees with the
    /// reference across everything it witnessed carries **no chain id** — so a
    /// test about chain ids needs a read that disagrees somewhere. One base in
    /// thirty is 3.3 %, comfortably under the 10 % the read filter drops at.
    fn read_with_one_mismatch(
        qname: &str,
        reference_sequence_id: usize,
        start: usize,
        mismatch_offset: usize,
    ) -> RecordBuf {
        use noodles_sam::alignment::record_buf::Sequence;

        const LENGTH: usize = 30;
        let mut bases = vec![b'A'; LENGTH];
        bases[mismatch_offset] = b'C';
        let mut record = read_named_with_length(qname, reference_sequence_id, start, LENGTH);
        *record.sequence_mut() = Sequence::from(bases);
        record
    }

    /// A read whose CIGAR is spelled out — the fixture helpers only build plain
    /// matches, and the halo exists for records a deletion widens.
    fn read_with_cigar(
        qname: &str,
        reference_sequence_id: usize,
        start: usize,
        ops: &[(noodles_sam::alignment::record::cigar::op::Kind, usize)],
    ) -> RecordBuf {
        use noodles_core::Position as RecordPosition;
        use noodles_sam::alignment::record::Flags;
        use noodles_sam::alignment::record::MappingQuality;
        use noodles_sam::alignment::record::cigar::op::{Kind, Op};
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

        let read_bases: usize = ops
            .iter()
            .filter(|(kind, _)| matches!(kind, Kind::Match | Kind::Insertion | Kind::SoftClip))
            .map(|(_, len)| len)
            .sum();
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(reference_sequence_id)
            .set_flags(Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
            .set_alignment_start(RecordPosition::try_from(start).expect("start is 1-based"))
            .set_cigar(ops.iter().map(|(kind, len)| Op::new(*kind, *len)).collect())
            .set_sequence(Sequence::from(vec![b'A'; read_bases]))
            .set_quality_scores(QualityScores::from(vec![30u8; read_bases]))
            .build()
    }

    /// **The defaults are production's, read from production's own constants.**
    ///
    /// Asserted against the `pub const`s rather than against literals: a literal
    /// here would let production retune a knob and ng silently keep the old
    /// value, which is the drift the "by name, not by literal" rule exists to
    /// prevent (arch §1.1).
    #[test]
    fn the_default_knobs_are_productions_five_constants() {
        let config = PileupGeneratorConfig::default();
        assert_eq!(
            config.max_snp_column_depth,
            crate::pileup::walker::DEFAULT_MAX_SNP_COLUMN_DEPTH
        );
        assert_eq!(
            config.max_indel_column_depth,
            crate::pileup::walker::DEFAULT_MAX_INDEL_COLUMN_DEPTH
        );
        assert_eq!(
            config.max_record_span,
            crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN
        );
        assert_eq!(
            config.mate_lookup_window,
            crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW
        );
        assert_eq!(
            config.max_active_reads,
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS
        );
    }

    /// The whole config, not just the one knob C1 constrains: production's
    /// defaults must all be walkable, or the `Default` impl ships a generator
    /// that cannot be constructed.
    #[test]
    fn the_default_config_is_accepted() {
        assert!(PileupGeneratorConfig::default().check().is_ok());
        assert!(a_generator(PileupGeneratorConfig::default()).is_ok());
    }

    /// **The ceiling is exactly where the narrowing starts lying.** This is the
    /// reason the knob is capped, asserted rather than described: at the ceiling
    /// `LocusLen` still reports the span it was given; one position past it, it
    /// reports a smaller number and no error.
    #[test]
    fn a_span_one_past_the_ceiling_is_where_a_coverage_run_starts_lying() {
        let at_ceiling = u64::from(MAX_RECORD_SPAN_CEILING);
        assert_eq!(
            u64::from(LocusLen::from_positions(at_ceiling).get()),
            at_ceiling,
            "a run at the ceiling must still describe itself exactly"
        );
        assert!(
            u64::from(LocusLen::from_positions(at_ceiling + 1).get()) < at_ceiling + 1,
            "one position past the ceiling the run saturates — silently, which is what \
             `check` exists to make unreachable"
        );
    }

    /// The boundary is inclusive: a span **at** the ceiling is describable, so it
    /// is accepted.
    #[test]
    fn a_max_record_span_at_the_ceiling_is_accepted() {
        let config = PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING,
            ..PileupGeneratorConfig::default()
        };
        assert!(config.check().is_ok());
        assert!(a_generator(config).is_ok());
    }

    /// One position past it is rejected, and the error names both numbers —
    /// production's own `--max-record-span` is an unbounded `u32`, so this is the
    /// knob a user is most likely to raise past what ng can describe.
    #[test]
    fn a_max_record_span_past_the_ceiling_is_rejected() {
        let config = PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING + 1,
            ..PileupGeneratorConfig::default()
        };
        let error = config
            .check()
            .expect_err("a span wider than a coverage run must be rejected");
        let PileupGeneratorConfigError::RecordSpanExceedsCoverageRun {
            max_record_span,
            ceiling,
        } = error;
        assert_eq!(max_record_span, MAX_RECORD_SPAN_CEILING + 1);
        assert_eq!(ceiling, MAX_RECORD_SPAN_CEILING);
    }

    /// **The rejection reaches the constructor**, which is the only door the
    /// generator has: a `check` nobody calls would leave the invariant a comment.
    #[test]
    fn the_constructor_rejects_a_span_no_coverage_run_could_describe() {
        let error = a_generator(PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING + 1,
            ..PileupGeneratorConfig::default()
        })
        .err()
        .expect("the constructor must reject what `check` rejects");
        assert!(
            error.to_string().contains("65535"),
            "the message must name the ceiling the user has to come back under, got: {error}"
        );
    }

    /// **Every knob reaches the walker, and reaches the right field.** Five
    /// deliberately distinct, non-default values: with the defaults, or with two
    /// knobs sharing a value, a transposed pair of fields in `to_walker_config`
    /// would pass.
    #[test]
    fn every_knob_reaches_the_walker_config_it_names() {
        let config = PileupGeneratorConfig {
            max_snp_column_depth: 11,
            max_indel_column_depth: 22,
            max_record_span: 33,
            mate_lookup_window: 44,
            max_active_reads: 55,
        };
        let walker_config = config.to_walker_config();
        assert_eq!(walker_config.max_snp_column_depth, 11);
        assert_eq!(walker_config.max_indel_column_depth, 22);
        assert_eq!(walker_config.max_record_span, 33);
        assert_eq!(walker_config.mate_lookup_window, 44);
        assert_eq!(walker_config.max_active_reads, 55);
    }

    /// A fresh generator has counted nothing — the baseline every later delta is
    /// folded onto (C3).
    #[test]
    fn a_fresh_generator_has_counted_nothing() {
        let generator = a_generator(PileupGeneratorConfig::default()).expect("default config");
        assert_eq!(*generator.counts(), PileupGeneratorCounts::default());
        assert_eq!(*generator.config(), PileupGeneratorConfig::default());
    }

    // --- C2: the region walk ------------------------------------------------

    /// A locus per covered position, and none where nothing was read — the shape
    /// of the whole output, before any boundary is involved.
    #[test]
    fn the_walk_emits_a_locus_at_every_covered_position() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let loci = loci_of(&mut generator, region(0, 1, 100), &reads);

        assert_eq!(
            anchors(&loci),
            (10..40).collect::<Vec<u64>>(),
            "one locus per position the read covered, and nothing outside it"
        );
    }

    /// **Records are dropped on their anchor, and the drop is counted.** The
    /// query returns every read *overlapping* the region, so a read starting
    /// before it produces records the region does not own.
    #[test]
    fn a_record_anchored_before_the_region_is_dropped_and_counted() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let loci = loci_of(&mut generator, region(0, 20, 100), &reads);

        assert_eq!(
            anchors(&loci),
            (20..40).collect::<Vec<u64>>(),
            "positions 10..=19 are covered but anchored outside the region"
        );
        assert_eq!(
            generator.counts().records_outside_region,
            10,
            "the ten dropped records are tallied, not silently discarded"
        );
    }

    /// **Neighbouring regions tile: no duplicates, no holes, coordinate order
    /// preserved** — the acceptance property for the clamp (plan's verification
    /// table, spec §2). Split at 50, and every locus of the whole-region walk
    /// appears exactly once, in the same order, with the same support.
    #[test]
    fn two_adjacent_regions_concatenate_into_the_single_region_walk() {
        let records: Vec<RecordBuf> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 5 + i * 20, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        let mut whole = a_generator(PileupGeneratorConfig::default()).expect("config");
        let expected = loci_of(&mut whole, region(0, 1, 100), &reads);

        let mut split = a_generator(PileupGeneratorConfig::default()).expect("config");
        let mut joined = loci_of(&mut split, region(0, 1, 50), &reads);
        joined.extend(loci_of(&mut split, region(0, 51, 100), &reads));

        assert!(!expected.is_empty(), "the fixture must produce loci at all");
        assert_eq!(
            joined, expected,
            "two adjacent regions must emit exactly what one region does"
        );
        assert_eq!(
            whole.counts().records_outside_region,
            0,
            "nothing is anchored outside a region covering every read"
        );
    }

    /// **The halo: a record anchored inside the region keeps the support lying
    /// beyond the boundary.** `widener`'s deletion grows a record anchored at 99
    /// out to 139; `beyond` starts at 120, entirely past the region's end, and
    /// folds into that record. A query for "reads overlapping the region" never
    /// returns `beyond`, and the record would then be emitted — by the right
    /// region, with a counter reading zero — carrying half its evidence (spec §2).
    #[test]
    fn the_halo_keeps_the_support_that_lies_past_the_region_end() {
        use noodles_sam::alignment::record::cigar::op::Kind;

        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_with_cigar(
                "widener",
                1,
                95,
                &[(Kind::Match, 5), (Kind::Deletion, 40), (Kind::Match, 30)],
            ),
            read_named_with_length("beyond", 1, 120, 30),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 60,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        let loci = loci_of(&mut generator, region(1, 1, 100), &reads);

        let widened = loci
            .iter()
            .find(|locus| locus.region.start == Position(99))
            .expect("the deletion opens a record at its anchor, inside the region");
        assert_eq!(
            widened.region.end,
            Position(139),
            "the record spans the whole deletion"
        );
        assert_eq!(
            total_obs(widened),
            2,
            "both the widener and the read lying entirely beyond the region support it"
        );
    }

    /// **The halo is stopped, not walked.** Nothing is anchored inside the
    /// region, so the walk has no reason to look at the reads filling the halo —
    /// and every record it finalised in there would be built at full depth and
    /// then thrown away by the clamp. `records_outside_region` is what that waste
    /// would show up as.
    #[test]
    fn the_walk_stops_at_the_region_end_instead_of_walking_the_whole_halo() {
        let records: Vec<RecordBuf> = (0..4)
            .map(|i| read_named_with_length(&format!("beyond{i}"), 1, 110 + i * 20, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 90,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        let loci = loci_of(&mut generator, region(1, 1, 100), &reads);

        assert!(loci.is_empty(), "no read covers the region itself");
        assert_eq!(
            generator.counts().records_outside_region,
            0,
            "a walk that ran the halo out would finalise ~80 records and clamp them all away"
        );
    }

    /// **A record anchored inside the region keeps growing past the boundary.**
    /// The stop rule has two halves and this is the second one: at position 101
    /// the record anchored at 99 is still open, so the walk must carry on rather
    /// than finalise it where it stands.
    ///
    /// The fixture is built so that only *widening* distinguishes the two halves.
    /// A record's footprint is fixed when it opens, so an early flush does not
    /// shorten it — the observable difference is the widen that has not happened
    /// yet: `later` is anchored at 99, inside the region and inside a plain query,
    /// and its deletion anchors at **101**, one position past the bound. Stop at
    /// 101 and the record is emitted spanning 99..=109; wait for it and it spans
    /// 99..=131.
    #[test]
    fn a_record_anchored_inside_the_region_still_widens_past_the_boundary() {
        use noodles_sam::alignment::record::cigar::op::Kind;

        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_with_cigar(
                "opener",
                1,
                95,
                &[(Kind::Match, 5), (Kind::Deletion, 10), (Kind::Match, 25)],
            ),
            read_with_cigar(
                "later",
                1,
                99,
                &[(Kind::Match, 3), (Kind::Deletion, 30), (Kind::Match, 27)],
            ),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 60,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        let loci = loci_of(&mut generator, region(1, 1, 100), &reads);

        let widened = loci
            .iter()
            .find(|locus| locus.region.start == Position(99))
            .expect("the first deletion opens a record at 99, inside the region");
        assert_eq!(
            widened.region.end,
            Position(131),
            "the walk stopped before the second deletion widened the record"
        );
    }

    /// A region whose halo runs past the end of the contig is still walkable —
    /// the last region of every contig is this case, so it must not be an error.
    #[test]
    fn a_region_ending_at_the_contig_end_still_queries_its_halo() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 61, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let loci = loci_of(&mut generator, region(0, 51, 100), &reads);

        assert_eq!(anchors(&loci), (61..91).collect::<Vec<u64>>());
    }

    /// **`begin_segment` cannot fail; the first `next_locus` is where the query
    /// does.** A contig the header does not have is the cheapest way to make the
    /// open fail, and the point is *where* the failure surfaces (arch §2.1).
    #[test]
    fn the_query_opens_at_the_first_locus_not_at_begin_segment() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        generator.begin_segment(region(9, 1, 100));

        let error = generator
            .next_locus(&reads)
            .expect_err("a region on a contig the file does not have cannot be queried");
        assert!(
            matches!(error, LocusGenerationError::Reads(_)),
            "the read query's failure reaches the caller as a read failure, got {error:?}"
        );
    }

    /// **A stream that failed is not an empty region.** The walker consumes an
    /// infallible item type, so a preparation failure has to be shed and
    /// reported as end-of-stream; if the generator read that as "the walk
    /// drained", a broken run would look like a region with no coverage.
    #[test]
    fn a_read_preparation_failure_is_reported_rather_than_read_as_an_empty_region() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), FailsToPrepare).expect("config");

        generator.begin_segment(region(0, 1, 100));

        let error = generator
            .next_locus(&reads)
            .expect_err("the preparer failed on the only read");
        assert!(
            matches!(error, LocusGenerationError::Reference(_)),
            "a failed preparation reaches the caller as the reference failure it is, got {error:?}"
        );
    }

    /// A failure inside the walk reaches the caller as `Walker` — the variant
    /// that exists because none of the other three describes a walk over inputs
    /// it already accepted. A deletion wider than `max_record_span` is the
    /// cheapest one to provoke.
    #[test]
    fn a_walk_failure_reaches_the_caller_as_a_walker_error() {
        use noodles_sam::alignment::record::cigar::op::Kind;

        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[read_with_cigar(
            "too_wide",
            1,
            95,
            &[(Kind::Match, 5), (Kind::Deletion, 40), (Kind::Match, 30)],
        )]);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 20,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        generator.begin_segment(region(1, 1, 100));

        let mut result = generator.next_locus(&reads);
        while matches!(result, Ok(Some(_))) {
            result = generator.next_locus(&reads);
        }
        let error = result.expect_err("a 41-position record exceeds a 20-position cap");
        assert!(
            matches!(error, LocusGenerationError::Walker(_)),
            "got {error:?}"
        );
    }

    /// Asking for a locus before any segment has begun is a normal `None`, not a
    /// panic and not an error: the pairing of `begin_segment` and `next_locus` is
    /// a contract nothing in the types enforces.
    #[test]
    fn next_locus_before_any_segment_yields_none() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        assert!(
            generator
                .next_locus(&reads)
                .expect("no segment is not a failure")
                .is_none()
        );
    }

    // --- C3: the allocator across segments ----------------------------------

    /// **The chain-id counters are folded as deltas, not summed.** Four reads
    /// across two regions mint four chain ids. The allocator is shared and
    /// `reset()` preserves its counters, so the second region's `RunSummary`
    /// reports the **run-to-date** four; adding the two summaries gives six —
    /// the triangular sum spec §8 names, in a number nothing else contradicts.
    #[test]
    fn the_allocators_counters_are_folded_as_deltas_across_regions() {
        let records: Vec<RecordBuf> = [5, 40, 90, 130]
            .iter()
            .map(|start| read_named_with_length(&format!("r{start}"), 1, *start, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 5,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        loci_of(&mut generator, region(1, 1, 80), &reads);
        loci_of(&mut generator, region(1, 81, 190), &reads);

        assert_eq!(
            generator.counts().reads_admitted,
            4,
            "two reads are anchored in each region, and no read spans the join"
        );
        assert_eq!(
            generator.counts().chain_allocations,
            4,
            "one chain id per solo read — summing the summaries would report six"
        );
    }

    /// **`reset()` between regions, and the abandoned walk is why it matters.**
    /// `active_count` falls back to zero on its own only when a walk drains; a
    /// walk abandoned mid-region leaves its reads counted as active for ever.
    /// With a cap of two, the second region's first admission is where that
    /// shows — as `ActiveReadsExhausted`, a failure with nothing to say about
    /// the region it fires in.
    #[test]
    fn an_abandoned_walk_does_not_leak_its_active_reads_into_the_next_region() {
        let records: Vec<RecordBuf> = [5, 6, 90, 91]
            .iter()
            .map(|start| read_named_with_length(&format!("r{start}"), 1, *start, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_active_reads: 2,
            max_record_span: 5,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        // Abandon the first region with both its reads still active.
        generator.begin_segment(region(1, 1, 80));
        generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .expect("the reads cover the region");

        let second = loci_of(&mut generator, region(1, 81, 190), &reads);

        assert!(
            !second.is_empty(),
            "the second region's reads must still be admissible"
        );
        assert_eq!(
            generator.counts().chain_allocations,
            4,
            "the abandoned walk's allocations are folded too, and its ids are not reissued"
        );
    }

    /// **Two regions never issue the same chain id.** This is what the shared
    /// allocator is *for*: an id names the fragment a read came from, and a
    /// later phasing step chains reads that share one — so two fragments in two
    /// regions holding the same id would be chained into a haplotype neither
    /// supports (spec §8).
    ///
    /// Asserted on the ids themselves rather than on the allocation *count*,
    /// which is blind to the failure: a fresh allocator per region mints one id
    /// per read either way, and both reads would come back as `ChainId(0)`.
    #[test]
    fn a_second_region_does_not_reissue_the_first_regions_chain_ids() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_with_one_mismatch("first", 1, 5, 10),
            read_with_one_mismatch("second", 1, 90, 10),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 5,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        let ids_of = |loci: &[SampleLocusObservations]| -> Vec<ChainId> {
            let mut ids: Vec<ChainId> = loci
                .iter()
                .flat_map(|locus| locus.observed_sequences.iter())
                .flat_map(|observation| observation.chain_ids.iter().copied())
                .collect();
            ids.sort_unstable();
            ids.dedup();
            ids
        };

        let first = ids_of(&loci_of(&mut generator, region(1, 1, 80), &reads));
        let second = ids_of(&loci_of(&mut generator, region(1, 81, 190), &reads));

        assert_eq!(first.len(), 1, "the mismatching read carries one chain id");
        assert_eq!(second.len(), 1);
        assert_ne!(
            first[0], second[0],
            "a fresh allocator per region would hand both fragments the same id"
        );
    }

    /// An abandoned walk gives the allocator back, so the next region can open
    /// at all. Without it the generator has no allocator to lend and cannot
    /// build a second walk.
    #[test]
    fn an_abandoned_walk_returns_the_allocator() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        generator.begin_segment(region(0, 1, 100));
        generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .expect("the read covers the region");
        generator.begin_segment(region(0, 1, 100));

        assert!(
            generator
                .next_locus(&reads)
                .expect("the second walk opens")
                .is_some()
        );
    }

    /// **Beginning a segment abandons the walk in flight.** One segment at a
    /// time is the contract; a half-drained walk of a region nobody is asking
    /// about must not leak into the next one's stream.
    #[test]
    fn beginning_a_segment_abandons_the_walk_in_flight() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        generator.begin_segment(region(0, 1, 100));
        let first = generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .expect("the read covers the region");
        assert_eq!(first.region.start, Position(10));

        generator.begin_segment(region(0, 1, 100));
        let restarted = generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .expect("the second segment starts over");
        assert_eq!(
            restarted.region.start,
            Position(10),
            "the new segment starts at its own first covered position"
        );
    }
}
