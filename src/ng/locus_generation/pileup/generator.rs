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

use crate::ng::locus_generation::{
    GeneratorCounts, LocusGenerationError, LocusGenerator, SampleLocusObservations,
};
use crate::ng::read::input::cursor::CursorCounts;
use crate::ng::read::input::sample_cursor::SampleCursor;
use crate::ng::read::input::{SampleIdentity, SampleReads};
use crate::ng::read::{PreparedRead, ReadPrepError, ReadPreparer};
use crate::ng::ref_seq::{ContigTable, EvictableRefSeq, RawRefSeq};
use crate::ng::types::{ContigId, GenomeRegion, Position};

use super::chain_id_allocator::{ChainIdAllocator, ChainIdAllocatorCounters};
use super::genome_walk::{PileupWalker, RegionReadSource, RunSummary};
use super::{DEFAULT_MAX_ACTIVE_READS, WalkerConfig};

/// The widest `max_record_span` this generator accepts: 65,535 reference
/// positions, the widest footprint a [`ReadWitness`] run can describe.
///
/// [`ReadWitness::Partial`] carries `offset_in_locus` and `positions_covered`
/// as `u16`s, minted through [`LocusLen::from_positions`], which **saturates**
/// rather than failing. A record footprint wider than this therefore makes a
/// partial witness report a *truncated* `positions_covered` — a wrong number, no
/// error, at exactly the long-deletion loci this generator exists to get right.
///
/// [`ReadWitness`]: crate::ng::locus_generation::ReadWitness
/// [`ReadWitness::Partial`]: crate::ng::locus_generation::ReadWitness::Partial
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
    "production's default max_record_span no longer fits a ReadWitness run: ng must either \
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
/// name would inherit a silent truncation ng's witness runs cannot survive. The
/// cap is not a constraint in practice — production's own default of 5,000 is
/// already unreachable with Illumina reads, and the ceiling is thirteen times
/// that again (owner, 2026-07-29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupGeneratorConfig {
    // The values are named, never spelled: `Default` reaches production's
    // constants by name so a retune arrives as a diff, and a literal in the
    // prose here would go stale silently while every test stayed green (review).
    /// Reads folded at a position with no indel anchored there. Defaults to
    /// [`DEFAULT_MAX_SNP_COLUMN_DEPTH`](crate::pileup::walker::DEFAULT_MAX_SNP_COLUMN_DEPTH).
    pub max_snp_column_depth: u32,
    /// Reads folded at a position where any read has an indel. Defaults to
    /// [`DEFAULT_MAX_INDEL_COLUMN_DEPTH`](crate::pileup::walker::DEFAULT_MAX_INDEL_COLUMN_DEPTH).
    pub max_indel_column_depth: u32,
    /// Widest record footprint before the walk fails. Defaults to
    /// [`DEFAULT_MAX_RECORD_SPAN`](crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN);
    /// ng additionally rejects anything above [`MAX_RECORD_SPAN_CEILING`].
    pub max_record_span: u32,
    /// How far a first mate stays available for pairing. Defaults to
    /// [`DEFAULT_MATE_LOOKUP_WINDOW`](crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW).
    pub mate_lookup_window: u32,
    /// Active-read ceiling: how many reads the walk will hold open at once. Defaults
    /// to [`DEFAULT_MAX_ACTIVE_READS`], **32,768 since 2026-08-05** (it was 4,096).
    ///
    /// **A cap, not a limit that fails.** Reaching it makes the walk give up a read —
    /// either the arriving one (`reads_shed_at_admission`) or the held read furthest
    /// from being kept (`reads_evicted_at_ceiling`) — and carry on; production's walker
    /// aborts instead. This is the only place a cap can bound what the walk *holds*,
    /// which is what bounds its memory: the per-column caps below act after the reads
    /// are already open.
    ///
    /// **What the ceiling is for, stated so it is not confused with the caps below.** It
    /// exists to keep a pathological pile-up from exhausting memory, and for nothing
    /// else. It is *not* a way to sample a deep position — that is
    /// `max_snp_column_depth`'s job, and a ceiling low enough to do it leaves positions
    /// covered by fewer reads than the input had for them, which is what 4,096 was doing
    /// (19,725 reads refused on one ~130× tomato chromosome). The number to check is
    /// `positions_short_of_cap`: zero means the ceiling shaped no position's evidence.
    ///
    /// **Where 32,768 comes from.** The deepest column measured on real data is 12,792
    /// reads (~130× tomato `SL4.0ch01` near 33,037,565, against a local typical of
    /// 86–133); the ceiling sits above the *active set*, which is larger than any one
    /// column, and the measured peaks are in `depth_cap.md`. The cost is bounded and
    /// linear: a held read is its sequence, its base qualities, its CIGAR and its cursor,
    /// so the ceiling's worst case is what a run must be able to afford.
    ///
    /// **The pending-mate ceiling was raised with it.** Every open read whose mate lies
    /// further along occupies an entry in the allocator's pending-mate map, whose own
    /// defensive cap was 10,000 — unreachable while the walk held 4,096 reads, and a hard
    /// `PendingMatesExhausted` as soon as it held more. That constant is now 1,000,000,
    /// and `pending_mates_high_water` reports how close a run came.
    ///
    /// **It no longer subsumes `max_snp_column_depth`.** At 4,096 a position could not
    /// gather 8,000 contributors, so the SNP cap could never fire; at 32,768 it can, and
    /// on the deepest real position it does.
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
    /// Reject a configuration a [`ReadWitness`] run could not describe.
    ///
    /// Called by [`PileupGenerator::new`], so a bad knob never reaches a locus.
    /// `witness_of` carries a `debug_assert` stating the same envelope; this is
    /// what makes it provable rather than hopeful, in release builds too.
    ///
    /// [`ReadWitness`]: crate::ng::locus_generation::ReadWitness
    pub fn check(&self) -> Result<(), PileupGeneratorConfigError> {
        if self.max_record_span > MAX_RECORD_SPAN_CEILING {
            return Err(PileupGeneratorConfigError::RecordSpanExceedsWitnessRun {
                max_record_span: self.max_record_span,
                ceiling: MAX_RECORD_SPAN_CEILING,
            });
        }
        // **Every knob has a floor as well as the one that has a ceiling** — the
        // review measured what a zero does, and all five are bad in one of three
        // ways: `max_active_reads` and `max_record_span` abort the walk on its
        // first read with an error naming a region rather than a knob;
        // `mate_lookup_window` silently fails to collapse a mate pair, so one
        // fragment gets two chain ids; and `max_snp_column_depth` walks
        // **successfully** and returns *zero loci for a covered region*, with the
        // truncations that explain it counted into a struct the dispatcher
        // cannot read. The last is the one this module's "no silent caps"
        // principle rules out outright.
        for (knob, value) in [
            ("max_snp_column_depth", self.max_snp_column_depth),
            ("max_indel_column_depth", self.max_indel_column_depth),
            ("max_record_span", self.max_record_span),
            ("mate_lookup_window", self.mate_lookup_window),
            ("max_active_reads", self.max_active_reads),
        ] {
            if value == 0 {
                return Err(PileupGeneratorConfigError::KnobIsZero { knob });
            }
        }
        Ok(())
    }

    /// The same five knobs in the shape the copied walker reads them.
    ///
    /// Written as an exhaustive struct literal rather than with `..default()`:
    /// a knob production adds to [`WalkerConfig`] is then a compile error here,
    /// which is a decision to make rather than a default to inherit silently.
    ///
    /// `pub(super)`, not `pub`: public, it was a **bypass** — the ceiling on
    /// `max_record_span` is enforced by [`PileupGenerator::new`], and a caller
    /// could reach a `WalkerConfig` from an unchecked config and hand it
    /// straight to this module's re-exported `run` (review).
    pub(super) fn to_walker_config(self) -> WalkerConfig {
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
    /// `max_record_span` is wider than a witness run can describe, so a partial
    /// witness of a maximally-wide record would report a saturated
    /// `positions_covered` — silently, since [`LocusLen::from_positions`] clamps
    /// rather than failing.
    ///
    /// [`LocusLen::from_positions`]: crate::ng::locus_generation::LocusLen::from_positions
    #[error(
        "max_record_span ({max_record_span}) exceeds the widest footprint a witness run can \
         describe ({ceiling}): a partially-witnessed record that wide would report a truncated \
         positions_covered rather than an error"
    )]
    RecordSpanExceedsWitnessRun { max_record_span: u32, ceiling: u32 },
    /// A knob set to zero. Each of the five means something the walk cannot do
    /// at all — fold no reads at a column, allow no active reads, keep no mate
    /// pending, allow no record span — and the walk's answer to each is either a
    /// run-fatal error naming a region rather than the knob, or a plausible
    /// wrong result. The knob is named so the message points at what to change.
    #[error("{knob} must be at least 1, got 0")]
    KnobIsZero { knob: &'static str },
}

/// Run-level counts for this generator, kept alongside the shared
/// [`LocusCounts`](super::super::LocusCounts) the dispatcher owns.
///
/// **Ten fields of three kinds**, and the kinds matter because only the first two come
/// off a walk at all:
///
/// - **Seven mirror production's
///   [`RunSummary`](crate::pileup::walker::RunSummary) field for field** (spec §7) —
///   everything on it bar `records_emitted`.
/// - **One mirrors ng's copy of `RunSummary` and has no production counterpart**:
///   `reads_silent_over_footprint`, the ninth field D2 added, which `parity.rs`'s counter
///   comparison drops by name for exactly that reason.
/// - **Two are the generator's own, off the walk entirely**:
///   `reads_declined_by_preparer` (a read the preparer refused never reaches a walk) and
///   `records_outside_region` (the region clamp is the generator's, not the walker's).
///
/// The link is to **production's** type deliberately: this module's `RunSummary`
/// re-export is `#[cfg(test)]`, and a `pub` item's doc must not link into a private one.
///
/// **`records_emitted` is deliberately not mirrored**: what a caller wants is loci
/// *kept*, which is the walk's emissions minus
/// [`records_outside_region`](Self::records_outside_region), and the kept count is
/// already `LocusCounts::loci_emitted`.
///
/// **These counters are the only thing this generator reports, and that is a
/// decision** (owner-facing, 2026-07-30, arch §3 *"the generator does not log"*):
/// nothing here writes to stderr, and the crate has no logging framework to write
/// to instead. A caller that wants to know what a walk saw reads this struct — the
/// dump tool prints every field of it — and a caller that wants a *message* is
/// asking for a policy a library cannot set for it. The one warning in this module
/// is `chain_id_allocator`'s high-water line, inherited verbatim from production.
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
    /// Reads the walk refused at admission because `max_active_reads` were already
    /// open — a depth cap acting on what the walk *holds*, where
    /// `column_depth_truncations` above counts positions where a cap acted on what
    /// it *uses*.
    ///
    /// Non-zero says the reads at some locus are a subsample of what the input
    /// offered. It cannot say *which* locus: a refused read has no `read_id` and its
    /// alignment is never decomposed, so nothing knows where it would have landed.
    /// The region is findable from `active_reads_high_water` sitting at the cap.
    pub reads_shed_at_admission: u64,
    /// Reads the ceiling gave back **after** admitting them, to make room for a read
    /// with a smaller sampling key. Beside `reads_shed_at_admission` above, the two are
    /// every read the ceiling removed from the evidence; the split says which rule
    /// removed it.
    pub reads_evicted_at_ceiling: u64,
    /// **The number that says whether the walk is throwing coverage away.** Positions the
    /// walk folded fewer reads than the per-position cap allows, while reads covering
    /// that position had been given up by the hold ceiling.
    ///
    /// A zero is the claim "no position in this run is short of the cap while reads
    /// covering it exist". It can only ever hold *up to the ceiling*: raise the coverage
    /// past what the walk holds and this starts counting again, which is the honest
    /// statement of a bounded-memory walk.
    pub positions_short_of_cap: u64,
    /// The reads missing from those positions, summed. `positions_short_of_cap` says how
    /// often; this says how much.
    pub short_of_cap_deficit: u64,
    /// The deepest column the walk assembled, before the per-position cap and after the
    /// mate-overlap collapse. A **max**, like `active_reads_high_water`.
    pub column_depth_high_water: u32,
    /// The largest the allocator's first-mates-awaiting-a-partner map ever got — the
    /// headroom under `MAX_PENDING_MATES`, measured. A **max**.
    pub pending_mates_high_water: u32,
    /// **Whether a capped position's evidence is a fair sample of the reads over it.**
    /// The reads seen at every position the cap acted on, and the reads it kept, each
    /// counted altogether and again for the reads that began left of the position. Two
    /// fractions to compare: equal means the cap sampled, a larger kept fraction means it
    /// preferred reads reaching the position from the left. See `RunSummary`'s copies.
    pub capped_column_reads_seen: u64,
    pub capped_column_reads_seen_placed_left: u64,
    pub capped_column_reads_kept: u64,
    pub capped_column_reads_kept_placed_left: u64,
    /// Reads silent over a whole record footprint — every base `N` or
    /// adaptor-masked, so never a contributor at any position and invisible to
    /// the per-locus `reads_without_observation` tally (spec §6).
    ///
    /// **Incremented since D2**, which is what spec §13's read-accounting
    /// assertion forced: the walk's active set carries a per-read "ever
    /// contributed" flag ([`ActiveRead::ever_contributed`](super::super::pileup))
    /// and tallies the reads that leave without it. Set before the mate-overlap
    /// collapse and before the depth cap, so a read either of those removed is
    /// counted where it belongs and not here as well.
    pub reads_silent_over_footprint: u64,
    /// Reads the **preparer** declined — `Ok(None)` from
    /// [`ReadPreparer::prepare_read`](crate::ng::read::ReadPreparer::prepare_read),
    /// "no usable observation", which ends nothing and is not an error.
    ///
    /// **Added by D2, and it reads zero for a reason that is not a bug:** no v1
    /// preparer declines anything (the only step that could was BAQ, deferred to
    /// spec §10). It exists because the read accounting has to balance —
    /// `reads_admitted` counts what reached the walk, and a declined read never
    /// does — and because "no preparer declines today" is a fact about today
    /// that a counter can state and a missing field cannot. The counter is
    /// exercised by a test preparer that declines.
    pub reads_declined_by_preparer: u64,
    /// Records the region clamp dropped, because their anchor fell outside the
    /// region being walked. Observably zero-sum across neighbouring regions,
    /// which is how the gap-free tiling argument stays checkable (spec §7).
    pub records_outside_region: u64,
    /// Reads that witnessed a locus in **more than one run** — blind in the
    /// middle of it — summed over the run. The case the witnessed-set
    /// representation exists for: a spliced read crossing a record widened over
    /// its intron used to be discarded whole, and is now recorded with both of
    /// its pieces.
    ///
    /// **Expected to read zero on DNA-seq**, structurally: a ref-skip emits no
    /// event, so an intron cannot widen a record by itself, and modern Illumina
    /// puts `N`s at read ends where they cannot make a hole. The number worth
    /// having is RNA-seq's, and this is the counter that can state it for any
    /// BAM — the divergence census counts the same thing but is `#[cfg(test)]`
    /// and only measures where production's walker also produced a record
    /// (owner, 2026-07-31).
    pub reads_with_holed_witness: u64,
    /// The locus positions inside those holes, summed over the reads — how much
    /// evidence the old discard threw away, rather than how many reads it threw
    /// away.
    pub hole_positions: u64,
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
            // Deliberately not mirrored: what a caller wants is loci *kept*,
            // which is this minus `records_outside_region`, and the kept count is
            // already `LocusCounts::loci_emitted`. Bound to `_` in the pattern
            // rather than discarded in the body, so it reads as a decision at the
            // one place the exhaustive destructure forces one.
            records_emitted: _,
            record_widen_events,
            mate_overlap_positions,
            chain_allocations,
            active_reads_high_water,
            mate_lookup_evictions,
            column_depth_truncations,
            reads_shed_at_admission,
            reads_silent_over_footprint,
            reads_with_holed_witness,
            hole_positions,
            reads_evicted_at_ceiling,
            positions_short_of_cap,
            short_of_cap_deficit,
            column_depth_high_water,
            pending_mates_high_water,
            capped_column_reads_seen,
            capped_column_reads_seen_placed_left,
            capped_column_reads_kept,
            capped_column_reads_kept_placed_left,
        } = *summary;
        debug_assert!(
            chain_allocations >= chain_ids_at_open.chain_allocations
                && mate_lookup_evictions >= chain_ids_at_open.mate_lookup_evictions,
            "the allocator's counters only grow; a baseline above the current value means \
             this walk is being folded against a baseline taken from a different allocator \
             — since D1 the walker keeps the allocator for a whole chromosome, so the two \
             can only disagree if a chromosome boundary was crossed between them",
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
        // A plain sum, for the plainest of the reasons on this page: the walk increments
        // this itself, and `begin_region` zeroes the whole summary — so each region
        // reports its own refusals and nothing is carried forward to be re-added.
        self.reads_shed_at_admission += reads_shed_at_admission;
        // A plain sum — and **the reason it is safe changed at D1**, so do not read the
        // fold without reading this. It used to be that each region's walk owned its own
        // active set. It does not: one walker, and so one `ActiveReads`, now lives for a
        // whole chromosome. The sum is correct because `ActiveReads::begin_region` *zeroes*
        // `silent_exits` at every region start, where `ActiveReads::reset` preserves it as a
        // run total. Swap the one call for the other and this line triangular-sums by the
        // region count, in a plausible-looking `u64`, exactly as the delta folds above would
        // if the allocator were replaced instead of reset.
        self.reads_silent_over_footprint += reads_silent_over_footprint;
        // Plain sums for the same reason: each region's walk finalises its own records, so a
        // holed read is counted by exactly the walk whose record it folded into.
        self.reads_with_holed_witness += reads_with_holed_witness;
        self.hole_positions += hole_positions;
        // Plain sums for the same reason as `reads_shed_at_admission` above: the walk
        // increments them itself and `begin_region` zeroes the whole summary.
        self.reads_evicted_at_ceiling += reads_evicted_at_ceiling;
        self.positions_short_of_cap += positions_short_of_cap;
        self.short_of_cap_deficit += short_of_cap_deficit;
        // Maxes, like `active_reads_high_water` — a peak over regions is the largest
        // single region's peak, not their sum. `pending_mates_high_water` comes off the
        // allocator, which the walker keeps for a whole chromosome, so a *delta* would be
        // meaningless for it where it is right for the two counters above it.
        self.column_depth_high_water = self.column_depth_high_water.max(column_depth_high_water);
        self.pending_mates_high_water = self.pending_mates_high_water.max(pending_mates_high_water);
        // Plain sums: each region's walk counts the columns it capped.
        self.capped_column_reads_seen += capped_column_reads_seen;
        self.capped_column_reads_seen_placed_left += capped_column_reads_seen_placed_left;
        self.capped_column_reads_kept += capped_column_reads_kept;
        self.capped_column_reads_kept_placed_left += capped_column_reads_kept_placed_left;
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
/// single-threaded (`locus_generation.md` §9) and nothing else touches the cell while a read is
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
    /// Reads the preparer declined since the last fold — see
    /// [`PileupGeneratorCounts::reads_declined_by_preparer`].
    ///
    /// It lives here, beside the latched error, for the same reason that does: the
    /// stream is inside the walk and cannot reach the generator's counts, and this cell
    /// is the one thing both hold. Taken by `end_walk`, so a region's tally is folded
    /// with that region's walk rather than leaking into the next one.
    declined: u64,
}

/// The sample's reads, prepared: the cursor's `AlignedRead`s canonicalised into the
/// `PreparedRead`s the walk consumes, with fatal errors shed into the shared
/// [`ReadPreparation`] and reported as end-of-stream.
///
/// **One of these for a chromosome, not for a region — D2.** It holds a
/// [`SampleCursor`], which is made once per chromosome and stays positioned, so the
/// stream is pointed at each region in turn rather than rebuilt for it. That is the
/// whole feature: the reads two consecutive regions share are decoded and filtered once
/// (`spec/alignment_cursor.md` §1).
///
/// **Lazy, and it must stay lazy.** The cursor is a pull source that collects nothing,
/// so what is resident is the reads overlapping the walker's current position — bounded
/// by depth, not by the region's length — plus what the cursor keeps for the next region,
/// which is bounded by the halo. Collecting here to make an ownership problem go away
/// would turn a depth-shaped footprint into a region-shaped one, on regions that run to
/// hundreds of kilobases (spec §7).
struct PreparedSampleReads<R: RawRefSeq, P: ReadPreparer> {
    reads: SampleCursor<R>,
    /// The region this stream's failures are attributed to — **the segment, not
    /// the halo-widened span it queries**. The span is an implementation detail
    /// of the halo; the segment is the unit a caller can act on.
    region: GenomeRegion,
    /// How far past the segment's end the cursor is pointed: `max_record_span`, the
    /// generator's own knob, carried here because this is where the widening now happens
    /// — see [`query`](Self::query).
    halo: u32,
    preparation: Rc<RefCell<ReadPreparation<P>>>,
    /// Set once the region is exhausted or an error was latched. Both are
    /// terminal **for the region**: a stream that shed an error must not resume, or the
    /// walk would carry on over a hole it was never told about. Cleared by
    /// [`move_to_region`](RegionReadSource::move_to_region) — a shed error has by then
    /// been taken by `end_walk`, which is what stops it outliving its region.
    done: bool,
}

impl<R: RawRefSeq, P: ReadPreparer> Iterator for PreparedSampleReads<R, P> {
    type Item = PreparedRead;

    fn next(&mut self) -> Option<PreparedRead> {
        while !self.done {
            let read = match self.reads.next_read() {
                None => {
                    self.done = true;
                    return None;
                }
                Some(Ok(read)) => read,
                Some(Err(source)) => {
                    let region = self.region;
                    return self.shed(LocusGenerationError::Reads {
                        region,
                        source: source.into(),
                    });
                }
            };
            let mut preparation = self.preparation.borrow_mut();
            // Split the borrow: `prepare_read` takes `&self` and `&mut
            // Self::Scratch`, which are two fields of the same cell.
            let ReadPreparation {
                preparer,
                scratch,
                declined,
                ..
            } = &mut *preparation;
            match preparer.prepare_read(read, scratch) {
                Ok(Some(prepared)) => return Some(prepared),
                // "No usable observation" — the preparer declined this read and the run
                // continues. **No v1 preparer returns it** (the only step that could
                // decline was BAQ, deferred), so this counter reads zero on every real
                // run today; it is counted anyway, because a read that never reached the
                // walk is missing from `reads_admitted` and the accounting has to say
                // where it went (D2).
                Ok(None) => {
                    *declined += 1;
                    continue;
                }
                Err(source) => {
                    drop(preparation);
                    // Matched exhaustively rather than through a catch-all:
                    // `ReadPrepError` is `#[non_exhaustive]`, which binds other
                    // crates but not this one, so a preparation failure ng
                    // cannot yet describe is a compile error here instead of a
                    // reference error it is not.
                    let region = self.region;
                    let error = match source {
                        ReadPrepError::Reference(source) => {
                            LocusGenerationError::Reference { region, source }
                        }
                    };
                    return self.shed(error);
                }
            }
        }
        None
    }
}

impl<R: RawRefSeq, P: ReadPreparer> RegionReadSource for PreparedSampleReads<R, P> {
    type Error = LocusGenerationError;

    /// Point the cursor at `segment`'s **halo-widened span**, and start the region.
    ///
    /// Returned rather than shed, unlike the failures `next` meets. A stream error
    /// happens with the walk already running and has nowhere to go but the latch; this
    /// one happens before the region has begun, so it goes back the way it came — which
    /// is where the per-region query's `OpenReadQuery` used to surface, and by the same
    /// `next_locus` call.
    fn move_to_region(&mut self, segment: GenomeRegion) -> Result<(), LocusGenerationError> {
        let query = self.query(segment);
        // **Pointed first, adopted after.** `AlignmentCursor::move_to_region` can fail
        // *late* — after it has already dropped what it kept and set its own region — so a
        // failure leaves the cursor somewhere unknown. Setting these first would leave the
        // stream claiming to be mid-region on a region it was never pointed at, and marked
        // not-done, which is a readable stream over a cursor in an unknown place. Every
        // such failure is terminal for the run today, so this costs nothing and stops the
        // property resting on that.
        self.reads
            .move_to_region(query)
            .map_err(|source| LocusGenerationError::OpenReadQuery {
                region: segment,
                source: source.into(),
            })?;
        self.region = segment;
        self.done = false;
        Ok(())
    }
}

impl<R: RawRefSeq, P: ReadPreparer> PreparedSampleReads<R, P> {
    /// The span the cursor is pointed at for `segment`: **the segment plus the halo.**
    ///
    /// A record anchored inside the segment can have a footprint reaching
    /// `max_record_span` past its end — a long deletion does exactly that — and reads
    /// that fold into it may lie *entirely beyond* the segment. Asking only for "reads
    /// overlapping the segment" never returns them, and the record is then emitted by the
    /// right region with part of its support missing, **and no counter notices**, because
    /// the record itself is not lost (spec §2). The extra reads are walked and their
    /// records clamped away unless anchored in the segment.
    ///
    /// **The widening lives here rather than at the call site now, and that is D2's one
    /// re-layering.** It moved with the query: `move_to_region` above is the seam the
    /// walker forwards through, and it takes the segment, because the segment is what a
    /// failure is attributed to. Deciding *which* span to ask for is still the caller's
    /// business and the caller is still this generator — the read layer is not consulted
    /// (spec §2).
    ///
    /// Widened only on the right, which is what makes consecutive queries nested forward
    /// and the cursor's forget rule as cheap at the call site as inside it (spec §13).
    fn query(&self, segment: GenomeRegion) -> GenomeRegion {
        GenomeRegion {
            contig: segment.contig,
            start: segment.start,
            end: Position(segment.end.get().saturating_add(u64::from(self.halo))),
        }
    }

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

/// The walk over one chromosome — **D2, and the thing that replaces a walker per region.**
///
/// It **owns** its read stream, which owns the sample's cursor, so nothing borrows the
/// `SampleReads` it was built from and nothing has to be materialised (arch §2.2). One per
/// chromosome, minted by the first `next_locus` that names a new one, and pointed at each
/// of that chromosome's regions in turn.
struct ChromosomeWalk<R: RawRefSeq, P: ReadPreparer> {
    /// The chromosome every region this walker is pointed at must lie on. Compared before
    /// each region, so the cursor's own `WrongChromosome` guard is never reached in
    /// correct operation (`spec/alignment_cursor.md` §4).
    contig: ContigId,
    walker: PileupWalker<PreparedSampleReads<R, P>, Arc<R>>,
}

/// The region walk in flight: what the walker is currently pointed at, and what its
/// contribution to the run's counters is measured against.
///
/// Not the walker itself any more — that lives for the chromosome ([`ChromosomeWalk`]).
/// This is `Some` exactly while a region has been opened and not yet ended, which is what
/// tells `end_walk` whether there is anything to fold.
struct RegionWalk {
    /// The region records are clamped to — the walk emits records anchored just
    /// outside it, because the query returns every read *overlapping* the
    /// region and the halo widens that further still.
    region: GenomeRegion,
    /// The chain-id allocator's counters as they stood when this region opened — the
    /// baseline its own contribution is measured against (C3).
    ///
    /// **Read through the walker now, not off an allocator being lent** (D1). A walker
    /// that lives for a chromosome never hands the allocator back between regions, so
    /// this is taken via [`PileupWalker::chain_id_counters`] the moment the region opens.
    /// It is still the same moment: nothing touches the allocator between `begin_segment`
    /// and the region opening.
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
/// `make_reference` is a third accessor and deliberately a **factory**: each of a
/// sample's k files gets its own raw accessor, because they are stateful readers and
/// sharing one would make k cursors share a file position and one sliding window
/// ([`SampleReads::cursor`]).
///
/// The chain-id allocator likewise lives here rather than inside the walk: ng
/// walks one region where production walks a chromosome, and a fresh allocator
/// per segment would give two fragments of different regions the same id (spec
/// §8). What that costs — `reset()` between segments, and counters that must be
/// folded as deltas because `reset()` preserves them — is C3's.
pub struct PileupGenerator<R: RawRefSeq + EvictableRefSeq + ContigTable, P: ReadPreparer> {
    /// The reference the walk fetches REF bases from. Built once, for the run,
    /// and handed to each walk as a shared handle.
    reference: Arc<R>,
    /// Builds a fresh read accessor for each file of the sample, which is what the
    /// cursor's mismatch-fraction filter reads.
    ///
    /// **Called once per file per *chromosome* now, not once per region — that is
    /// perf-review finding L2 closed** (D2). It used to be called at every
    /// `begin_segment`'s first `next_locus`, and the cost is not in building the accessor
    /// — [`WindowedRefSeq::new`](crate::ng::ref_seq::WindowedRefSeq) is a path and a contig
    /// table, no I/O — but in its **first fetch**, which reaches
    /// `RawChromReader::for_contig` and parses the whole `.fai` (~2,580 records on GRCh38)
    /// before opening the FASTA. That was 9,820 FASTA opens on chromosome 21; it is now one
    /// per file per chromosome.
    ///
    /// **Boxed, which is what drops `PileupGenerator`'s third type parameter** (arch §3.6).
    /// The indirection costs one virtual call per file per chromosome — the factory is not
    /// on any hot path any more, which is the point. It stays a *factory* rather than
    /// becoming a shared handle because a sample's k cursors interleave over the same
    /// coordinates: one accessor between them would be one file position and one sliding
    /// window serving k readers, and its eviction would be driven by whichever asked last.
    ///
    /// **`+ Send` costs nothing *here* and is kept for the fan-out.** A bare `Box<dyn FnMut>`
    /// makes this generator `!Send` whatever `R` and `P` are. Nothing is blocked today — the
    /// `Rc<RefCell<ReadPreparation>>` below already blocks it — but `locus_generation.md` §9
    /// gives each worker its own generator, so this would become a second blocker to remove
    /// for no gain.
    ///
    /// **The STR generator deliberately does *not* carry the bound**, and the difference is a
    /// caller rather than a preference: `ng_ssr_cohort_stutter` hands every file a clone of one
    /// `Arc<WindowedRefSeq>`, which is `!Send` because `WindowedRefSeq` holds a `RefCell` and is
    /// therefore `!Sync`. Every caller of *this* constructor builds a fresh accessor per call,
    /// capturing only a `PathBuf` and two `Arc`s.
    make_reference: Box<dyn FnMut() -> R + Send>,
    /// The preparer and its scratch, lent to each region's read stream.
    preparation: Rc<RefCell<ReadPreparation<P>>>,
    /// Lives across chromosomes so `next_id` never repeats — **`None` exactly while a
    /// [`ChromosomeWalk`] holds it**, which is the invariant `enter_chromosome` keeps.
    ///
    /// **The window it is absent for widened at D2**, from a region to a chromosome: a
    /// walker built per region borrowed it per region, and one built per chromosome borrows
    /// it for the chromosome. `reset()` between regions still happens — inside the walker,
    /// at [`WalkerState::begin_region`], where the counters it preserves are documented
    /// along with what breaks if it is replaced instead.
    ///
    /// An `Option` rather than a swap with a fresh allocator: a placeholder
    /// starting at zero is the state this whole arrangement exists to avoid, and
    /// it would fail silently. Absent, it fails loudly.
    chain_ids: Option<ChainIdAllocator>,
    config: PileupGeneratorConfig,
    counts: PileupGeneratorCounts,
    /// The region [`begin_segment`](Self::begin_segment) was given, before the walker has
    /// been pointed at it.
    current_region: Option<GenomeRegion>,
    /// The walker and the sample's cursor, for as long as the run stays on one chromosome.
    /// Minted by the first `next_locus` on each, and rebuilt at every boundary — nothing in
    /// a cursor survives a chromosome change (`spec/alignment_cursor.md` §4).
    chromosome: Option<ChromosomeWalk<R, P>>,
    /// The region walk in flight — `Some` between the `next_locus` that opened a region and
    /// the `end_walk` that closes it.
    walk: Option<RegionWalk>,
    /// What the cursors of *already-retired* chromosomes did. A cursor's tallies die with
    /// it at the boundary, so they are taken there; the live one is asked directly. See
    /// [`cursor_counts`](Self::cursor_counts).
    retired_cursor_counts: CursorCounts,
    /// The sample this generator opened its first cursor for. **One sample per generator** —
    /// a second one is refused rather than answered out of the first one's files. See
    /// [`open_walk`](Self::open_walk).
    ///
    /// `None` until the first cursor opens, so a generator is unclaimed until it reads
    /// something.
    sample: Option<SampleIdentity>,
    /// A fatal error that has happened and not yet been reported.
    ///
    /// The stream's errors are **shed** (it hands the walker an infallible item
    /// type), so they surface a call later than they happen; this is where they
    /// wait, moved out of the walk's cell by `end_walk` so an error can never
    /// outlive the region that raised it. Never `Some` while
    /// [`failed`](Self::failed) is `false`.
    pending_failure: Option<LocusGenerationError>,
    /// Set the moment a fatal error happens, and never cleared.
    ///
    /// **Every error this generator can raise is terminal for the run** (spec
    /// §7), so after one it emits nothing more: `next_locus` reports the
    /// error(s) once each and then answers `Ok(None)` for ever, whatever
    /// `begin_segment` is called with. Without the latch it re-opened the query
    /// and re-walked the region — the same reads admitted twice, and one
    /// fragment handed two chain ids, which is the corruption the run-lifetime
    /// allocator exists to prevent, arrived at from the other direction (review).
    failed: bool,
}

impl<R: RawRefSeq + EvictableRefSeq + ContigTable, P: ReadPreparer> PileupGenerator<R, P> {
    /// Build a generator over `reference` (the walk's REF fetches),
    /// `make_reference` (the cursor's per-file accessor factory) and
    /// `preparer` (per-read canonicalisation), with `config` checked before
    /// anything is held.
    ///
    /// **`make_reference` is boxed here**, which is what keeps it off this type's
    /// parameter list (arch §3.6). The `'static` that comes with the box is the one real
    /// constraint: the factory is called for the whole life of the generator, at every
    /// chromosome boundary, so it must not borrow anything shorter-lived — which every
    /// caller already satisfies by capturing owned paths and shared indexes. `Send` is required
    /// for the reason on the field, and every caller satisfies that too.
    ///
    /// Fails only on a configuration a witness run could not describe — see
    /// [`PileupGeneratorConfig::check`].
    pub fn new(
        reference: Arc<R>,
        make_reference: impl FnMut() -> R + Send + 'static,
        preparer: P,
        config: PileupGeneratorConfig,
    ) -> Result<Self, PileupGeneratorConfigError> {
        config.check()?;
        Ok(Self {
            reference,
            make_reference: Box::new(make_reference),
            preparation: Rc::new(RefCell::new(ReadPreparation {
                preparer,
                scratch: P::Scratch::default(),
                latched_error: None,
                declined: 0,
            })),
            chain_ids: Some(ChainIdAllocator::with_caps(
                config.max_active_reads,
                config.mate_lookup_window,
            )),
            config,
            counts: PileupGeneratorCounts::default(),
            current_region: None,
            chromosome: None,
            walk: None,
            retired_cursor_counts: CursorCounts::default(),
            sample: None,
            pending_failure: None,
            failed: false,
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

    /// **What the cursors did**, summed over every chromosome walked so far — reads decoded
    /// against reads replayed from memory, and regions reused against regions jumped.
    ///
    /// # Why this exists, and why it is not on [`counts`](Self::counts)
    ///
    /// **The saving this whole design is for is invisible in the output.** A generator that
    /// mints a cursor per region and one that keeps a cursor per chromosome emit exactly the
    /// same loci — that is the correctness requirement — so the only way to tell them apart
    /// is to count what was *avoided* (`spec/alignment_cursor.md` §11.5).
    ///
    /// Review proved the point: replacing the chromosome comparison in
    /// [`open_walk`](Self::open_walk) with `if true`, so that every region built a fresh
    /// cursor, left **all 1,557 tests green**. `CursorCounts` existed, and nothing above the
    /// cursor could read it.
    ///
    /// It is kept off `PileupGeneratorCounts` deliberately: that struct is printed verbatim
    /// by the dump tools, whose output is asserted byte-identical against committed
    /// baselines, and this is a fact about the *reader* rather than about the loci.
    ///
    /// Cursors retired at a chromosome boundary are added up as they go; the live one is
    /// asked directly, so the answer is current at any moment.
    pub fn cursor_counts(&self) -> CursorCounts {
        let mut total = self.retired_cursor_counts;
        if let Some(chromosome) = &self.chromosome {
            total += chromosome.walker.reads().reads.counts();
        }
        total
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

    /// End the region walk in flight, if any: fold its counters into the run's totals and
    /// take its shed error out of the shared preparation cell.
    ///
    /// **The only place `walk` is cleared**, which is what keeps "a region is open exactly
    /// while `walk` is `Some`" true without anyone having to remember it.
    ///
    /// **What it no longer does — D2.** It used to take the chain-id allocator back and
    /// `reset()` it, because the walker was dropped here. The walker now lives for the
    /// chromosome and keeps the allocator, so the reset moved to
    /// [`WalkerState::begin_region`], which runs at the *next* region's start. The reset
    /// itself is unchanged and still necessary: carried across regions un-reset, a pending
    /// first mate from one region pairs with a read in another (the eviction test compares
    /// raw positions) and `active_count` climbs toward `ActiveReadsExhausted` (spec §8).
    ///
    /// **The fold has to happen here and not there**, before the walker is pointed at
    /// anything else: it reads the walk's summary, which `begin_region` is about to clear.
    ///
    /// **It also takes the walk's shed error with it**, which is the third piece
    /// of per-region state and the one that went missing (review): the error
    /// slot lives on the preparation cell, which the generator holds for its
    /// whole life, so an error shed by a region **abandoned** before it drained
    /// stayed there and became the *next* region's failure — reported after that
    /// region had emitted every one of its loci, and against a region that never
    /// saw the failing read. Moved here, an error cannot outlive the walk that
    /// raised it.
    fn end_walk(&mut self) {
        // **The region ends with its walk.** `begin_segment` sets it again
        // immediately; every other caller of this wants the segment over. Left
        // set, a `next_locus` after the walk drained saw `walk.is_none()`, opened
        // the region again and **re-emitted the whole of it** — with every read
        // admitted a second time and its fragment handed a second chain id
        // (review).
        self.current_region = None;
        let Some(walk) = self.walk.take() else {
            return;
        };
        // PANIC-FREE: a region walk is only ever opened through `open_walk`, which puts a
        // chromosome walk in place first and never clears one while `walk` is `Some`.
        let summary = self
            .chromosome
            .as_ref()
            .expect("a region walk is open, so its chromosome's walker is here")
            .walker
            .summary();
        self.counts
            .fold_region_walk(&summary, walk.chain_ids_at_open);
        let mut preparation = self.preparation.borrow_mut();
        // Taken, not read: the cell outlives every walk, so a tally left in it would be
        // folded again at the next region's end — the same shape as the shed error below.
        self.counts.reads_declined_by_preparer += std::mem::take(&mut preparation.declined);
        if let Some(error) = preparation.latched_error.take() {
            self.failed = true;
            // First error wins: a later one describes a run that is already over.
            self.pending_failure.get_or_insert(error);
        }
    }

    /// Mark the run failed and hand the error straight back to the caller.
    ///
    /// Every error this generator raises is terminal (spec §7), so there is no
    /// path that records one and keeps walking.
    fn fail(&mut self, error: LocusGenerationError) -> LocusGenerationError {
        self.failed = true;
        error
    }

    /// The next locus of the region begun, or `None` once the walk drains.
    ///
    /// Records whose **anchor** falls outside the region are dropped and
    /// tallied in [`records_outside_region`](PileupGeneratorCounts::records_outside_region)
    /// — the rule that makes neighbouring regions tile without duplicates or
    /// holes, since typed regions tile the genome gap-free and disjointly
    /// (spec §2).
    ///
    /// # After a fatal error, nothing more comes out
    ///
    /// Every error here is terminal for the run (spec §7). Each is reported
    /// once — a shed stream error surfaces one call after it happened, so there
    /// can be two — and then this answers `Ok(None)` for ever, whatever
    /// `begin_segment` is called with afterwards. Without that latch a caller
    /// that logged the error and asked again (a shape the return type positively
    /// invites) got the whole region re-walked: the same reads admitted twice,
    /// and one fragment carrying two chain ids.
    pub fn next_locus(
        &mut self,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        if let Some(error) = self.pending_failure.take() {
            return Err(error);
        }
        if self.failed {
            return Ok(None);
        }
        let Some(region) = self.current_region else {
            return Ok(None);
        };
        if self.walk.is_none()
            && let Err(error) = self.open_walk(region, reads)
        {
            return Err(self.fail(error));
        }
        loop {
            // PANIC-FREE: opened immediately above, and every path that clears
            // the walk returns rather than looping.
            let clamp = self.walk.as_ref().expect("the walk was just opened").region;
            let walker = &mut self
                .chromosome
                .as_mut()
                .expect("the walk was just opened, so its chromosome's walker is here")
                .walker;
            match walker.next() {
                Some(Ok(locus)) => {
                    if clamp.contains(locus.region.start) {
                        return Ok(Some(locus));
                    }
                    self.counts.records_outside_region += 1;
                }
                Some(Err(source)) => {
                    // Fatal and terminal: the walker yields nothing after an
                    // error, so the walk ends here and gives back the allocator
                    // it was lent. The walk's own failure is the proximate one
                    // and is reported now; a stream error shed earlier in the
                    // same walk is left pending rather than lost, and the next
                    // call reports that one before fusing.
                    self.end_walk();
                    return Err(self.fail(LocusGenerationError::Walker {
                        region: clamp,
                        source,
                    }));
                }
                None => {
                    // The walk is over. A read stream that shed a fatal error
                    // also reported end-of-stream, so the walker drains
                    // normally and the error is only visible after `end_walk`
                    // has moved it out of the walk — checked **before**
                    // `Ok(None)`, or a broken query would read as an empty
                    // region.
                    self.end_walk();
                    if let Some(error) = self.pending_failure.take() {
                        return Err(error);
                    }
                    return Ok(None);
                }
            }
        }
    }

    /// Point the walk at `region`, minting the chromosome's walker and cursor first if this
    /// is the first region on it.
    ///
    /// The walk is **stopped** at the region's own end rather than run to the end of the
    /// halo the cursor is pointed at — see [`PileupWalker::move_to_region`], and
    /// [`PreparedSampleReads::query`] for the halo itself.
    ///
    /// # One sample per generator, and it is checked here
    ///
    /// A cursor is opened for **one sample's files** and kept for a whole chromosome. Handing
    /// this generator a second sample would answer that sample out of the first one's files —
    /// with no error, and output of exactly the right shape. So the sample is remembered and a
    /// foreign one is refused ([`ForeignSample`](LocusGenerationError::ForeignSample)).
    ///
    /// **No caller does this today, and the STR generator's did.** The check is here because
    /// the two generators are handed a `&SampleReads` per call by the same trait, so the
    /// mistake is equally available on both — it was simply made on the other one first.
    fn open_walk(
        &mut self,
        region: GenomeRegion,
        reads: &SampleReads,
    ) -> Result<(), LocusGenerationError> {
        let sample = reads.identity();
        if self.sample.as_ref().is_some_and(|opened| *opened != sample) {
            return Err(LocusGenerationError::ForeignSample { region });
        }
        if self
            .chromosome
            .as_ref()
            .is_none_or(|chromosome| chromosome.contig != region.contig)
        {
            self.enter_chromosome(region, reads)?;
            // Adopted only once a cursor exists, so a failed open leaves the generator
            // unclaimed and a retry with the same sample still works.
            self.sample = Some(sample);
        }
        self.release_reference_behind(region);
        // Saturating, and in the safe direction: the walker's positions are
        // `u32`, so a region ending past 4 Gb cannot be walked at all, and a
        // stop bound clamped to `u32::MAX` makes the walk run long rather than
        // stop early.
        let stop_after = region.end.get().min(u64::from(u32::MAX)) as u32;
        // PANIC-FREE: minted immediately above when it was absent or on another chromosome.
        let walker = &mut self
            .chromosome
            .as_mut()
            .expect("the chromosome walk was just ensured")
            .walker;
        walker.move_to_region(region, stop_after)?;
        // **After the move, not before**, so a region that could not be reached leaves no
        // half-open walk behind for `end_walk` to fold.
        let chain_ids_at_open = walker.chain_id_counters();
        self.walk = Some(RegionWalk {
            region,
            chain_ids_at_open,
        });
        Ok(())
    }

    /// Let every reference reader drop the bases this walk has gone past.
    ///
    /// # The problem, which is not the cursor's and predates it
    ///
    /// A reference reader here is a **sliding window over the FASTA**, not the whole
    /// chromosome. It extends forward as it is asked for more, and it only ever shrinks when
    /// someone tells it to. Nobody told it to. So a run walking a chromosome forward ended up
    /// holding **one byte in memory for every base it had walked** — measured directly at
    /// 1,999,995 bytes resident after a 2 Mb walk, against 150 bytes with releasing turned on.
    /// On human chromosome 1 that is around 250 MB, against a walk that otherwise peaks near
    /// 25 MB.
    ///
    /// It never showed up, because every measurement behind this generator used a
    /// tandem-repeat-targeted file. Its coverage is sparse enough that consecutive reads are
    /// far apart, and the reader jumps rather than extending — which is exactly the case that
    /// cannot grow.
    ///
    /// **Three readers, because they cannot be one.** The walk asks *what base is here*, the
    /// read preparer asks for the window under each read it left-aligns, and each input file's
    /// read filter asks for the window under each read it mismatch-checks. They are separate
    /// because each is a stateful reader with its own file position, and sharing one would
    /// give several consumers one position between them.
    ///
    /// # The bound, and why it is generous
    ///
    /// Everything any of them will read from here on lies at or after this region's start,
    /// less a margin: a read overlapping the region may begin before it, and a record's
    /// footprint may reach back to the read that opened it. `max_record_span` is a bound on
    /// both, so this releases everything before `region.start - max_record_span`.
    ///
    /// **Being wrong here costs a re-read, never an answer.** Releasing is a hint: a base
    /// asked for after it was released is simply read again. That is what lets the margin be
    /// generous rather than exact.
    fn release_reference_behind(&mut self, region: GenomeRegion) {
        let keep_from = region
            .start
            .get()
            .saturating_sub(u64::from(self.config.max_record_span));
        self.reference.evict_before(keep_from);
        self.preparation
            .borrow()
            .preparer
            .evict_reference_before(keep_from);
        if let Some(chromosome) = &self.chromosome {
            chromosome
                .walker
                .reads()
                .reads
                .evict_reference_before(keep_from);
        }
    }

    /// How many reference bases this generator's readers are holding between them — the walk's,
    /// the preparer's, and the cursor's.
    ///
    /// **The bound made observable.** Without it "the window is bounded" is an argument about
    /// where [`release_reference_behind`](Self::release_reference_behind) is called from, and a
    /// call site quietly dropped costs memory without changing a single answer — which is how
    /// it went unnoticed in the first place.
    pub fn resident_reference_bases(&self) -> usize {
        self.reference.resident_bases()
            + self
                .preparation
                .borrow()
                .preparer
                .resident_reference_bases()
            + self.chromosome.as_ref().map_or(0, |chromosome| {
                chromosome.walker.reads().reads.resident_reference_bases()
            })
    }

    /// Retire the walker on the chromosome being left and build the one for `region`'s.
    ///
    /// **This is where the cursor is made — once per chromosome, never per region**, which
    /// is the whole of D2 (`spec/alignment_cursor.md` §1). A cursor keeps the reads it has
    /// decoded and hands them back to the next region that can use them; nothing in one
    /// survives a chromosome change, so the boundary is where a new one is minted rather
    /// than a place the old one is repositioned (§4).
    ///
    /// The chain-id allocator crosses the boundary with the run, not with the walker:
    /// production's own walk preserves `next_id` across chromosomes so two fragments on two
    /// chromosomes never share an identity, and taking it out of the retiring walker and
    /// into the new one is how that is kept.
    fn enter_chromosome(
        &mut self,
        region: GenomeRegion,
        reads: &SampleReads,
    ) -> Result<(), LocusGenerationError> {
        if let Some(retiring) = self.chromosome.take() {
            // Taken before the walker is consumed: a cursor's tallies die with it, and they
            // are what says whether the cursor did anything (see `cursor_counts`).
            self.retired_cursor_counts += retiring.walker.reads().reads.counts();

            let mut chain_ids = retiring.walker.into_chain_ids();
            // **Redundant with `WalkerState::begin_region`, and kept anyway.** The next
            // region's `move_to_region` resets the allocator one call later, so on the
            // ordinary path this changes nothing — an earlier comment here claimed it was
            // what stopped a pending first mate pairing across a chromosome, which
            // `begin_region` already does. It earns its place on the failure path: if the
            // `cursor` call below fails, the allocator stays in `self.chain_ids` with this
            // chromosome's pending mates in it, and no `begin_region` is coming.
            chain_ids.reset();
            self.chain_ids = Some(chain_ids);
        }
        // Reborrowed into a local: the factory is a field, and `cursor` would otherwise
        // hold a borrow of `self` across the call that stores the walker back into `self`.
        let make_reference = &mut self.make_reference;
        let cursor = reads
            .cursor(region.contig, make_reference)
            .map_err(|source| LocusGenerationError::OpenReadQuery { region, source })?;
        let prepared = PreparedSampleReads {
            reads: cursor,
            region,
            halo: self.config.max_record_span,
            preparation: Rc::clone(&self.preparation),
            // Not pointed at a region yet, so there is nothing to read. `move_to_region`
            // clears this; without it the walker's construction peek would ask a cursor
            // that has no region and read the answer as end of input.
            done: true,
        };
        // Taken **after** the cursor, which is the only fallible step: taken before, an
        // unopenable chromosome would carry the run's chain-id sequence off with the error.
        let chain_ids = self.chain_ids.take().expect(
            "the allocator is absent only while a chromosome walk holds it, and the one \
             that did has just been retired",
        );
        let walker = super::genome_walk::run(
            prepared,
            Arc::clone(&self.reference),
            &self.config.to_walker_config(),
        )
        .adopting_chain_ids(chain_ids);
        self.chromosome = Some(ChromosomeWalk {
            contig: region.contig,
            walker,
        });
        Ok(())
    }
}

/// The generic path's slot in the dispatcher, filled (C4).
///
/// The segment payload is `()` because `RegionKind::Generic` carries none — a
/// generic region is its geometry and nothing else, and the region reaches the
/// generator through [`begin_segment`](LocusGenerator::begin_segment) alone.
///
/// Both methods delegate to the inherent ones of the same name, which is what
/// the tests in this module drive: an inherent method wins name resolution
/// against a trait method, so `generator.next_locus(&reads)` is the two-argument
/// one below and never a mis-resolved trait call.
impl<R: RawRefSeq + EvictableRefSeq + ContigTable, P: ReadPreparer> LocusGenerator<()>
    for PileupGenerator<R, P>
{
    fn begin_segment(&mut self, region: GenomeRegion) {
        PileupGenerator::begin_segment(self, region);
    }

    fn next_locus(
        &mut self,
        _segment: &(),
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        PileupGenerator::next_locus(self, reads)
    }

    /// **The only way these ten counters are reachable once the generator is boxed.**
    /// `GeneratorSlot` erases the type, so before this they had no reader that was not
    /// a test — and a walk that emitted nothing for a covered region counted the
    /// truncations explaining it into a struct nobody could see (Milestone C review).
    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        Some(GeneratorCounts::Pileup(PileupGenerator::counts(self)))
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

    /// A passthrough preparer that counts the reads pulled through it — the
    /// only seam that sees the query being consumed, one read at a time.
    #[derive(Default)]
    struct CountingPreparer {
        prepared: std::cell::Cell<u32>,
    }

    impl ReadPreparer for CountingPreparer {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: AlignedRead,
            scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            self.prepared.set(self.prepared.get() + 1);
            PassthroughPreparer.prepare_read(read, scratch)
        }
    }

    /// A preparer that **declines** the reads whose qname starts with `Self::0`, and passes
    /// every other through.
    ///
    /// `Ok(None)` is the "no usable observation" answer the trait documents and **no v1
    /// preparer gives** — the only step that could was BAQ, deferred (spec §10). So
    /// `reads_declined_by_preparer` cannot be exercised by any real preparer, and a counter
    /// nothing can move is indistinguishable from one that is wired to nothing. This is what
    /// moves it.
    struct DeclinesRead(&'static str);

    impl ReadPreparer for DeclinesRead {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: AlignedRead,
            scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            if read.qname.starts_with(self.0.as_bytes()) {
                return Ok(None);
            }
            PassthroughPreparer.prepare_read(read, scratch)
        }
    }

    /// A preparer that silences **every base** of the reads whose qname starts with
    /// `Self::0`, by placing their adaptor boundary at their own first position.
    ///
    /// That is the shape of the read `reads_silent_over_footprint` exists for: admitted,
    /// walked past, and never a contributor anywhere, because the G1 filter answers "in the
    /// adaptor" at every position (`cigar_cursor::base_in_adaptor` — forward strand, so
    /// `ref_pos >= boundary` is every base). An all-`N` read would do the same thing to the
    /// walk but not reach it: against an all-`A` reference every base mismatches, and the
    /// real read filter drops it before admission.
    struct SilencesRead(&'static str);

    impl ReadPreparer for SilencesRead {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: AlignedRead,
            scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            let silence = read.qname.starts_with(self.0.as_bytes());
            let prepared = PassthroughPreparer.prepare_read(read, scratch)?;
            Ok(prepared.map(|mut read| {
                if silence {
                    debug_assert!(
                        !read.is_reverse_strand,
                        "the boundary is a forward-strand one; a reverse read would be \
                         silenced by `ref_pos <= boundary` instead and this fixture would \
                         silence nothing",
                    );
                    read.adaptor_boundary = Some(read.alignment_start);
                }
                read
            }))
        }
    }

    /// A preparer that fails on **one named read** and passes every other.
    ///
    /// The selectivity is the point: a region that fails on its first read never
    /// gets to be a region with an *unreported* error, because the walk drains
    /// and reports immediately. To abandon a region mid-drain with its failure
    /// still pending, the failure has to come from a read the walker reaches
    /// while earlier reads are still producing loci.
    struct FailsToPrepareRead(&'static str);

    impl ReadPreparer for FailsToPrepareRead {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: AlignedRead,
            scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            if read.qname == self.0.as_bytes() {
                return Err(ReadPrepError::Reference(
                    crate::ng::ref_seq::RefSeqError::UnknownContig(ContigId(7)),
                ));
            }
            PassthroughPreparer.prepare_read(read, scratch)
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
    ) -> Result<PileupGenerator<InMemoryRefSeq, P>, PileupGeneratorConfigError> {
        PileupGenerator::new(Arc::new(fixture_bases()), fixture_bases, preparer, config)
    }

    fn a_generator(
        config: PileupGeneratorConfig,
    ) -> Result<PileupGenerator<InMemoryRefSeq, PassthroughPreparer>, PileupGeneratorConfigError>
    {
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
    fn loci_of<R: RawRefSeq + EvictableRefSeq + ContigTable, P: ReadPreparer>(
        generator: &mut PileupGenerator<R, P>,
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

    /// Total supporting reads at a locus, summed over its observations.
    fn total_obs(locus: &SampleLocusObservations) -> u32 {
        locus
            .observations
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

    /// **Four of the five defaults are production's, read from production's own
    /// constants.**
    ///
    /// Asserted against the `pub const`s rather than against literals: a literal
    /// here would let production retune a knob and ng silently keep the old
    /// value, which is the drift the "by name, not by literal" rule exists to
    /// prevent (arch §1.1).
    ///
    /// **The fifth is ng's own from 2026-08-05**, and the assertion says so rather than
    /// being deleted: `max_active_reads` is 32,768 where production's constant is 4,096,
    /// because production's value was refusing reads at the door on ordinary
    /// whole-genome data. A test that merely stopped checking the knob would let ng drift
    /// back without anybody noticing; this one fails if either number moves.
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
            config.max_active_reads, 32_768,
            "ng's own, not production's"
        );
        assert_eq!(
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS,
            4096,
            "production's, unchanged — ng is frozen out of editing it, and this states \
             what ng diverged *from*"
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
    fn a_span_one_past_the_ceiling_is_where_a_witness_run_starts_lying() {
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
            .expect_err("a span wider than a witness run must be rejected");
        let PileupGeneratorConfigError::RecordSpanExceedsWitnessRun {
            max_record_span,
            ceiling,
        } = error
        else {
            panic!("wrong variant: {error:?}");
        };
        assert_eq!(max_record_span, MAX_RECORD_SPAN_CEILING + 1);
        assert_eq!(ceiling, MAX_RECORD_SPAN_CEILING);
    }

    /// **Every knob has a floor, and the review measured why.** A zero is
    /// accepted by arithmetic and then means something the walk cannot do: two
    /// of the five abort the walk on its first read with an error naming a
    /// region rather than the knob, one silently fails to collapse a mate pair,
    /// and `max_snp_column_depth` walks **successfully** and returns zero loci
    /// for a covered region.
    #[test]
    fn a_knob_set_to_zero_is_rejected_at_construction() {
        let sound = PileupGeneratorConfig::default();
        let zeroed = [
            (
                "max_snp_column_depth",
                PileupGeneratorConfig {
                    max_snp_column_depth: 0,
                    ..sound
                },
            ),
            (
                "max_indel_column_depth",
                PileupGeneratorConfig {
                    max_indel_column_depth: 0,
                    ..sound
                },
            ),
            (
                "max_record_span",
                PileupGeneratorConfig {
                    max_record_span: 0,
                    ..sound
                },
            ),
            (
                "mate_lookup_window",
                PileupGeneratorConfig {
                    mate_lookup_window: 0,
                    ..sound
                },
            ),
            (
                "max_active_reads",
                PileupGeneratorConfig {
                    max_active_reads: 0,
                    ..sound
                },
            ),
        ];
        for (name, config) in zeroed {
            let error = config
                .check()
                .expect_err(&format!("{name} = 0 must be rejected"));
            let PileupGeneratorConfigError::KnobIsZero { knob } = error else {
                panic!("{name} = 0 was rejected for the wrong reason: {error:?}");
            };
            assert_eq!(knob, name, "the message must name the knob to change");
            assert!(
                a_generator(config).is_err(),
                "{name} = 0 must be rejected at construction, not mid-walk"
            );
        }
    }

    /// **The rejection reaches the constructor**, which is the only door the
    /// generator has: a `check` nobody calls would leave the invariant a comment.
    #[test]
    fn the_constructor_rejects_a_span_no_witness_run_could_describe() {
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
        assert!(
            split.counts().records_outside_region > 0,
            "and the split run really does exercise the clamp — asserting only the \
             whole-region zero would pass with the clamp deleted"
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
            matches!(error, LocusGenerationError::OpenReadQuery { .. }),
            "an *open* that failed is its own variant, distinct from a stream that broke \
             mid-region, got {error:?}"
        );
        assert_eq!(
            error.region(),
            Some(region(9, 1, 100)),
            "and it names the region it was opening for"
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
            matches!(error, LocusGenerationError::Reference { .. }),
            "a failed preparation reaches the caller as the reference failure it is, got {error:?}"
        );
        assert_eq!(
            error.region(),
            Some(region(0, 1, 100)),
            "attributed to the region being walked, not to the halo-widened span queried"
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
            matches!(error, LocusGenerationError::Walker { .. }),
            "got {error:?}"
        );
        assert_eq!(error.region(), Some(region(1, 1, 100)));
        assert!(
            error.to_string().contains("contig 1:1-100"),
            "the message names the region, which is what a log line has to go on: {error}"
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

    // --- C4: wired into the dispatcher --------------------------------------

    /// **End to end through the public surface**: a filled `Generic` slot, a
    /// typed-region stream, and the loci come out of
    /// `SampleLocusObservationsIterator` with the shared tally agreeing.
    ///
    /// The satellite region is in the fixture because the tally has to keep
    /// telling the two kinds of nothing apart (spec §5): it is `OutOfScope`
    /// permanently, while the STR and bundle slots stay `NotImplemented` — and
    /// the generic slot, which was one of those until this step, is now neither.
    #[test]
    fn a_filled_generic_slot_mints_loci_through_the_public_iterator() {
        use crate::ng::locus_generation::{
            GeneratorSet, GeneratorSlot, SampleLocusObservationsIterator, UnhandledReason,
        };
        use crate::ng::region_typing::{RegionKind, TypedRegion};

        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let generator = a_generator(PileupGeneratorConfig::default()).expect("config");
        let generators = GeneratorSet::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Generator(Box::new(generator)),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );
        let typed = vec![
            Ok(TypedRegion {
                region: region(0, 1, 50),
                kind: RegionKind::Generic,
            }),
            Ok(TypedRegion {
                region: region(0, 51, 100),
                kind: RegionKind::Satellite,
            }),
        ];

        let mut stream = SampleLocusObservationsIterator::new(typed.into_iter(), reads, generators);
        let mut loci = Vec::new();
        for locus in &mut stream {
            loci.push(locus.expect("the walk succeeds"));
        }

        assert_eq!(
            anchors(&loci),
            (10..40).collect::<Vec<u64>>(),
            "the generic region's covered positions, and only those"
        );
        let counts = stream.counts();
        assert_eq!(counts.regions_in, 2);
        assert_eq!(counts.regions_handled, 1, "the generic region is handled");
        assert_eq!(counts.loci_emitted, 30);
        assert_eq!(
            counts.unhandled_out_of_scope, 1,
            "the satellite stays permanently out of scope"
        );
        assert_eq!(
            counts.unhandled_not_implemented, 0,
            "nothing is unimplemented in this fixture — the generic slot is filled"
        );
    }

    /// **The generator's own counters are reachable through the boxed generator, and
    /// they are the running ones.**
    ///
    /// `GeneratorSlot` erases the type, so until `LocusGenerator::counts` existed these
    /// nine had no reader at all outside this file — which is what made a walk that
    /// emitted nothing for a covered region *totally* silent: the truncations that
    /// explained it were counted and then unreachable (Milestone C review).
    ///
    /// Asserted on a **non-zero** value, and through the public iterator rather than
    /// the concrete type: a surface returning a default-constructed struct, or one
    /// wired to a counter nothing increments, passes an `is_some()` and fails this.
    #[test]
    fn the_generators_own_counts_are_reachable_through_the_public_iterator() {
        use crate::ng::locus_generation::{
            GeneratorCounts, GeneratorSet, GeneratorSlot, SampleLocusObservationsIterator,
            UnhandledReason,
        };
        use crate::ng::region_typing::{RegionKind, TypedRegion};

        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let generators = GeneratorSet::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Generator(Box::new(
                a_generator(PileupGeneratorConfig::default()).expect("config"),
            )),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );
        let typed = vec![Ok(TypedRegion {
            region: region(0, 1, 100),
            kind: RegionKind::Generic,
        })];

        let mut stream = SampleLocusObservationsIterator::new(typed.into_iter(), reads, generators);
        let loci = (&mut stream)
            .collect::<Result<Vec<_>, _>>()
            .expect("the walk succeeds");
        assert_eq!(loci.len(), 30);

        let Some(GeneratorCounts::Pileup(counts)) = stream.generators().generic_counts() else {
            panic!("the generic slot's counts must be reachable once it is filled");
        };
        assert_eq!(
            counts.reads_admitted, 1,
            "the running tally the walk actually kept, not a fresh one"
        );
        assert!(
            stream.generators().ssr_counts().is_none(),
            "an unfilled slot counts nothing"
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
                .flat_map(|locus| locus.observations.iter())
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

    // --- D2: the cursor per chromosome, and what crosses the boundary ---------

    /// **A run that changes chromosome mints a new cursor and keeps walking.**
    ///
    /// A cursor covers one chromosome and refuses a region on any other
    /// (`spec/alignment_cursor.md` §4), so the generator has to notice the boundary and
    /// build a new one. Nothing before D2 could see this: a walker built per region was
    /// handed a fresh query each time and never had a chromosome of its own to be wrong
    /// about.
    ///
    /// **Loci from both, and chain ids that do not repeat across the boundary.** The
    /// allocator is the run's, not the walker's, so retiring one chromosome's walker must
    /// hand it on rather than let it go — otherwise two fragments on two chromosomes carry
    /// the same identity, which is the corruption the shared allocator exists to prevent.
    #[test]
    fn a_region_on_a_new_chromosome_mints_a_cursor_and_keeps_the_chain_ids_going() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_with_one_mismatch("on-chr1", 0, 10, 10),
            read_with_one_mismatch("on-chr2", 1, 10, 10),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let ids_of = |loci: &[SampleLocusObservations]| -> Vec<ChainId> {
            let mut ids: Vec<ChainId> = loci
                .iter()
                .flat_map(|locus| locus.observations.iter())
                .flat_map(|observation| observation.chain_ids.iter().copied())
                .collect();
            ids.sort_unstable();
            ids.dedup();
            ids
        };

        let first = loci_of(&mut generator, region(0, 1, 100), &reads);
        let second = loci_of(&mut generator, region(1, 1, 200), &reads);

        assert!(!first.is_empty(), "chr1's read must yield loci");
        assert!(
            !second.is_empty(),
            "chr2's read must yield loci too — a cursor stuck on chr1 would answer with \
             nothing, or refuse"
        );
        assert_eq!(
            generator.counts().reads_admitted,
            2,
            "one read on each chromosome, each admitted by its own walk"
        );

        let (first_ids, second_ids) = (ids_of(&first), ids_of(&second));
        assert_eq!(first_ids.len(), 1);
        assert_eq!(second_ids.len(), 1);
        assert_ne!(
            first_ids[0], second_ids[0],
            "a chromosome boundary that dropped the allocator instead of handing it on \
             would give both fragments the same id"
        );
    }

    /// **Going back to a chromosome already walked is a new cursor, not a rewind.**
    ///
    /// Nothing in a cursor survives a chromosome change, so the boundary is where one is
    /// minted rather than a place an old one is repositioned (`spec/alignment_cursor.md`
    /// §4). The failure this pins is the cheap-looking optimisation: keeping a cursor per
    /// chromosome in a map, or comparing against the wrong one, and answering the third
    /// region out of the first's kept reads.
    #[test]
    fn returning_to_a_chromosome_already_walked_yields_its_loci_again() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("on-chr1", 0, 10, 30),
            read_named_with_length("on-chr2", 1, 10, 30),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let first = anchors(&loci_of(&mut generator, region(0, 1, 100), &reads));
        let between = anchors(&loci_of(&mut generator, region(1, 1, 200), &reads));
        let again = anchors(&loci_of(&mut generator, region(0, 1, 100), &reads));

        assert!(!first.is_empty(), "chr1 must yield loci");
        assert!(
            !between.is_empty(),
            "chr2 must yield loci too — discarded, a chromosome that silently returned \
             nothing would leave the round trip below passing",
        );
        assert_eq!(
            first, again,
            "the same region on the same chromosome must yield the same loci, whatever was \
             walked in between"
        );
    }

    /// **`reads_silent_over_footprint` is summed region by region, so each region's walk
    /// must report only its own.**
    ///
    /// The active set keeps this tally, and `ActiveReads::reset` deliberately *preserves*
    /// it as a run total — which was the same thing while every region got a fresh set. A
    /// walker that lives for a chromosome shares one set across every region on it, so the
    /// per-region reset has to zero the tally or `fold_region_walk`'s plain sum
    /// triangular-sums it.
    ///
    /// **Three regions, one silent read each, and the number that fails is 6 against 3.**
    /// Two regions would give 3 against 2, which is close enough to read as an off-by-one
    /// somewhere else.
    #[test]
    fn each_regions_silent_reads_are_folded_once_and_not_again_at_the_next_region() {
        let records: Vec<RecordBuf> = [10u32, 50, 90]
            .iter()
            .flat_map(|start| {
                [
                    read_named_with_length(&format!("silent{start}"), 1, *start as usize, 30),
                    read_named_with_length(&format!("speaks{start}"), 1, *start as usize, 30),
                ]
            })
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator_with(
            PileupGeneratorConfig {
                max_record_span: 5,
                ..PileupGeneratorConfig::default()
            },
            SilencesRead("silent"),
        )
        .expect("config");

        loci_of(&mut generator, region(1, 1, 40), &reads);
        loci_of(&mut generator, region(1, 41, 80), &reads);
        loci_of(&mut generator, region(1, 81, 120), &reads);

        assert_eq!(
            generator.counts().reads_admitted,
            6,
            "two reads anchored in each of the three regions, and no read spans a join"
        );
        assert_eq!(
            generator.counts().reads_silent_over_footprint,
            3,
            "one silent read per region — a tally carried across the boundaries would \
             report six"
        );
    }

    /// **The cursor is kept across regions, and this is the test that can tell.**
    ///
    /// Review replaced the chromosome comparison in `open_walk` with `if true`, so that
    /// every region minted a fresh cursor — the per-region query D2 exists to remove — and
    /// **all 1,557 tests passed.** They had to: the loci are identical either way, which is
    /// the correctness requirement. The only thing that moves is what the reader *avoided*,
    /// and until `cursor_counts` nothing above the cursor could see it
    /// (`spec/alignment_cursor.md` §11.5).
    ///
    /// **Asserted on `regions_reusing`, not on `reads_decoded`.** Decode counts depend on
    /// how the fixture's reads fall across the regions; "the second and third regions did
    /// not jump" is the property, and it is what a cursor per region cannot produce — a
    /// fresh cursor has no last region to compare against, so its first move always jumps.
    #[test]
    fn the_cursor_is_kept_across_regions_rather_than_minted_for_each() {
        // The reads at 35 and 75 are what makes replay happen at all: each runs 30 bases,
        // so it overlaps two consecutive regions and is owed to both. Without one straddling
        // every join a cursor can reuse honestly and still replay nothing.
        let records: Vec<RecordBuf> = [10u32, 35, 50, 75, 90]
            .iter()
            .map(|start| read_named_with_length(&format!("r{start}"), 1, *start as usize, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 5,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        loci_of(&mut generator, region(1, 1, 40), &reads);
        loci_of(&mut generator, region(1, 41, 80), &reads);
        loci_of(&mut generator, region(1, 81, 120), &reads);

        let counts = generator.cursor_counts();
        assert_eq!(
            (counts.regions_jumping, counts.regions_reusing),
            (1, 2),
            "three ascending regions through one cursor jump once — at the first, which \
             has nothing to reuse — and reuse twice: {counts:?}",
        );
        assert!(
            counts.reads_replayed > 0,
            "and the point of reusing is that reads come back without being decoded again: \
             {counts:?}",
        );
    }

    /// **A cursor's tallies survive the chromosome it belonged to.**
    ///
    /// They live on the cursor, and a cursor dies at every boundary, so a run over several
    /// chromosomes would report only the last one's unless they are taken as it is retired.
    /// That would understate the saving by the number of chromosomes — silently, and in the
    /// flattering direction.
    #[test]
    fn cursor_tallies_are_taken_from_a_chromosome_before_it_is_retired() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("on-chr1", 0, 10, 30),
            read_named_with_length("on-chr2", 1, 10, 30),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        loci_of(&mut generator, region(0, 1, 100), &reads);
        let after_first = generator.cursor_counts();
        loci_of(&mut generator, region(1, 1, 200), &reads);
        let after_second = generator.cursor_counts();

        assert!(
            after_first.reads_decoded > 0,
            "the first chromosome must decode something, or this test cannot fail: \
             {after_first:?}",
        );
        assert!(
            after_second.reads_decoded > after_first.reads_decoded,
            "the second chromosome's reads must add to the first's, not replace them: \
             {after_first:?} then {after_second:?}",
        );
    }

    /// **A second sample through one generator is refused, not answered.**
    ///
    /// A generator opens a reader for one sample's files and keeps it for a whole chromosome,
    /// but the `LocusGenerator` trait hands it a `&SampleReads` afresh on every call. Without
    /// this check it would answer the second sample out of the **first** sample's files, with
    /// no error and output of exactly the right shape.
    ///
    /// No caller of *this* generator loops samples today. The STR generator's did, and review
    /// caught it there; the two are handed their reads by the same trait, so the mistake is
    /// equally available here.
    ///
    /// Sample B's read is 60 bases from the region, so any locus it produces is A's.
    #[test]
    fn a_second_sample_through_one_generator_is_refused() {
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 10, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 1, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        assert!(
            !loci_of(&mut generator, region(0, 1, 40), &sample_a).is_empty(),
            "sample A's read must yield loci, or this test cannot fail",
        );

        generator.begin_segment(region(0, 1, 40));
        match generator.next_locus(&sample_b) {
            Err(LocusGenerationError::ForeignSample { region: named }) => {
                assert_eq!(
                    named,
                    region(0, 1, 40),
                    "the refusal names the region asked for"
                );
            }
            other => panic!("expected a refusal naming the foreign sample, got {other:?}"),
        }
    }

    /// **A generator per sample is the shape that works** — the other half of the test above,
    /// because "refuses everything" would pass that one on its own.
    #[test]
    fn a_generator_per_sample_gives_each_sample_its_own_reads() {
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 10, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 1, 10, 30)]);

        let anchors_for = |reads: &SampleReads| {
            let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");
            anchors(&loci_of(&mut generator, region(0, 1, 40), reads))
        };

        assert!(
            !anchors_for(&sample_a).is_empty(),
            "sample A's read is on this contig and must yield loci",
        );
        assert!(
            anchors_for(&sample_b).is_empty(),
            "sample B's only read is on another contig and must yield nothing — which is the \
             whole point",
        );
    }

    /// A reference that records what it was told to release, and delegates everything else.
    ///
    /// The fixture reference is in memory, so it holds no window and `resident_reference_bases`
    /// reads zero however badly the release is wired. Asking *what was released, and when* is
    /// falsifiable where asking *how much is held* is not.
    struct ReleaseSpy {
        inner: InMemoryRefSeq,
        released: Arc<std::sync::Mutex<Vec<u64>>>,
    }

    impl EvictableRefSeq for ReleaseSpy {
        fn evict_before(&self, pos: u64) {
            self.released.lock().expect("no panic holds this").push(pos);
        }
    }

    /// Delegated, so the spy answers the table it actually reads from. A spy that
    /// invented its own would be refused by `AlignmentFile::cursor`'s contig-table
    /// check before the walk it exists to observe ever started.
    impl crate::ng::ref_seq::ContigTable for ReleaseSpy {
        fn contigs(&self) -> &crate::fasta::ContigList {
            self.inner.contigs()
        }
    }

    impl RawRefSeq for ReleaseSpy {
        fn fetch_raw_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), crate::ng::ref_seq::RefSeqError> {
            self.inner.fetch_raw_into(contig, start_1based, length, dst)
        }
    }

    impl crate::ng::ref_seq::RefSeq for ReleaseSpy {
        fn fetch_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), crate::ng::ref_seq::RefSeqError> {
            self.inner.fetch_into(contig, start_1based, length, dst)
        }
    }

    /// **Every region tells the reference readers what they may release, and the mark moves
    /// forward with the walk.**
    ///
    /// A reference reader is a sliding window over the FASTA. It extends forward as it is asked
    /// for more, and it only shrinks when it is told to — and nothing was telling it. A run
    /// walking a chromosome forward therefore held one byte in memory for every base it had
    /// walked: measured directly at 1,999,995 bytes resident after a 2 Mb walk, against 150
    /// bytes once the release is wired (`ref_seq`'s own tests own that measurement).
    ///
    /// **It went unseen because of the fixture.** Every measurement behind this generator used
    /// a tandem-repeat-targeted file, whose coverage is sparse enough that the reader jumps
    /// between reads rather than extending — the one case that cannot grow.
    ///
    /// This asserts the two things the generator is responsible for: that it releases at all,
    /// and that the mark advances. What is *safe* to release is the margin below.
    #[test]
    fn every_region_releases_the_reference_behind_it() {
        let records: Vec<RecordBuf> = (0..5)
            .map(|index| read_named_with_length(&format!("r{index}"), 1, 30 + index * 20, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        // Two recorders, because there are two readers to release: the walk's own, and the
        // one each input file's read filter holds inside the cursor.
        let released = Arc::new(std::sync::Mutex::new(Vec::new()));
        let released_by_cursor = Arc::new(std::sync::Mutex::new(Vec::new()));
        let spy = ReleaseSpy {
            inner: fixture_bases(),
            released: Arc::clone(&released),
        };
        let make_reference = {
            let released_by_cursor = Arc::clone(&released_by_cursor);
            move || ReleaseSpy {
                inner: fixture_bases(),
                released: Arc::clone(&released_by_cursor),
            }
        };
        let mut generator = PileupGenerator::new(
            Arc::new(spy),
            make_reference,
            PassthroughPreparer,
            PileupGeneratorConfig {
                max_record_span: 5,
                ..PileupGeneratorConfig::default()
            },
        )
        .expect("config");

        for index in 0..5u64 {
            let start = 30 + index * 20;
            loci_of(&mut generator, region(1, start, start + 19), &reads);
        }

        let marks = released.lock().expect("no panic holds this").clone();
        assert_eq!(
            marks.len(),
            5,
            "one release per region, so a walk of any length holds one region's worth: {marks:?}",
        );
        assert!(
            marks.windows(2).all(|pair| pair[1] > pair[0]),
            "the mark must move forward with the walk, or it releases the same nothing every \
             time: {marks:?}",
        );
        // Region 1 starts at 30 and the margin is `max_record_span` = 5.
        assert_eq!(
            marks[0], 25,
            "the mark is the region's start less one `max_record_span`, which bounds both how \
             far back a read overlapping the region can begin and how far a record's footprint \
             can reach: {marks:?}",
        );

        // The cursor's reader is released too, and it is the one that was growing hardest —
        // the read filter reads the reference once per surviving read, so it touches every
        // covered base of the chromosome.
        let by_cursor = released_by_cursor
            .lock()
            .expect("no panic holds this")
            .clone();
        assert_eq!(
            by_cursor, marks,
            "the cursor's reader is released at the same mark and as often as the walk's — \
             the release runs after the cursor is minted, so even the first region reaches it: \
             {by_cursor:?}",
        );
    }

    /// **A failure names the region it happened in, not the first region of the
    /// chromosome.**
    ///
    /// The stream's `region` field is what every shed error carries, and it is set when the
    /// chromosome's stream is built — so without reassigning it at each
    /// `move_to_region` every failure on a chromosome is blamed on whichever region opened
    /// it. Review proved the omission passes the whole suite: the only shipped test for a
    /// preparation failure walks a single region, where the two regions are the same one.
    #[test]
    fn a_failure_names_the_region_it_happened_in_not_the_chromosomes_first() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("fine", 1, 10, 30),
            read_named_with_length("boom", 1, 90, 30),
        ]);
        let mut generator = a_generator_with(
            PileupGeneratorConfig {
                max_record_span: 5,
                ..PileupGeneratorConfig::default()
            },
            FailsToPrepareRead("boom"),
        )
        .expect("config");

        // A clean first region, so the chromosome's stream is opened against it.
        let first = region(1, 1, 40);
        assert!(
            !loci_of(&mut generator, first, &reads).is_empty(),
            "the first region must walk cleanly, or the failure below has nothing to be \
             misattributed to",
        );

        let second = region(1, 81, 120);
        generator.begin_segment(second);
        let mut error = None;
        loop {
            match generator.next_locus(&reads) {
                Ok(Some(_)) => continue,
                Ok(None) => break,
                Err(reported) => {
                    error = Some(reported);
                    break;
                }
            }
        }

        let error = error.expect("the second region's read fails preparation");
        assert_eq!(
            error.region(),
            Some(second),
            "the failure is in {second}, not in the region that happened to open the \
             chromosome's stream: {error}",
        );
    }

    // --- the review's fixes: properties that were load-bearing and unpinned --

    /// **The halo is exactly `max_record_span` wide.** The shipped halo test put
    /// its far read 20 positions into a 60-position halo, so a halo of *half*
    /// the width passed it (review). Here the far read starts at 131 — one past
    /// half of a 60-position halo over a region ending at 100 — and folds into a
    /// record anchored at the region's last position whose deletion runs the
    /// full span.
    #[test]
    fn the_halo_reaches_a_read_at_its_far_end() {
        use noodles_sam::alignment::record::cigar::op::Kind;

        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_with_cigar(
                "opener",
                1,
                96,
                &[(Kind::Match, 5), (Kind::Deletion, 59), (Kind::Match, 30)],
            ),
            read_named_with_length("far", 1, 131, 30),
        ]);
        let mut generator = a_generator(PileupGeneratorConfig {
            max_record_span: 60,
            ..PileupGeneratorConfig::default()
        })
        .expect("config");

        let loci = loci_of(&mut generator, region(1, 1, 100), &reads);

        let widened = loci
            .iter()
            .find(|locus| locus.region.start == Position(100))
            .expect("the deletion opens a record at 100, the region's last position");
        assert_eq!(widened.region.end, Position(159));
        assert_eq!(
            total_obs(widened),
            2,
            "a halo narrower than max_record_span never returns the read at 131"
        );
    }

    /// **The locus anchored exactly at the region's last position is emitted.**
    /// The stop rule is `walker_pos > region.end`; as `>=` it stops one position
    /// early and leaves a **one-base hole at every region boundary** (review).
    ///
    /// **The fixture has to leave nothing open at the boundary, or the rule's
    /// other half covers for the mutation** — which is how it survived the first
    /// version of this test, and every other test in the file. Coverage starting
    /// *at* `region.end` is that case: the walker jumps the uncovered span, so
    /// when it arrives at 100 no record is open and only the comparison decides
    /// whether position 100 is walked at all.
    #[test]
    fn the_locus_at_the_regions_last_position_is_not_lost_to_the_stop() {
        // chr2 (200 bp), so the read may run past the region's end without
        // running past the contig's.
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 1, 100, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let loci = loci_of(&mut generator, region(1, 1, 100), &reads);

        assert_eq!(
            anchors(&loci),
            vec![100],
            "the read covers 100..=129 and the region owns exactly its first position"
        );
    }

    /// **All three of `fold_region_walk`'s rules, on the two counters no BAM
    /// fixture reaches.** `mate_lookup_evictions` needs paired reads at eviction
    /// distance and `active_reads_high_water` needs a max rather than a sum;
    /// dropping the delta on the first and summing the second both survived the
    /// end-to-end tests (review), so the fold is exercised directly.
    #[test]
    fn the_fold_deltas_the_allocators_counters_and_maxes_the_high_water() {
        let mut counts = PileupGeneratorCounts::default();
        let baseline = ChainIdAllocatorCounters {
            chain_allocations: 10,
            active_reads_high_water: 7,
            mate_lookup_evictions: 4,
            pending_mates_high_water: 3,
        };
        let summary = RunSummary {
            reads_admitted: 3,
            records_emitted: 99,
            record_widen_events: 1,
            mate_overlap_positions: 2,
            // Run-to-date values, as `summary()` reports them: the allocator is
            // shared and `reset()` preserves its counters.
            chain_allocations: 14,
            active_reads_high_water: 5,
            mate_lookup_evictions: 6,
            column_depth_truncations: 8,
            // ng's, and a plain sum for the same reason as the line below: the walk
            // increments it and `begin_region` zeroes the summary it lives on.
            reads_shed_at_admission: 9,
            // ng's, and a plain sum — each region's walk owns its own active set (D2).
            reads_silent_over_footprint: 2,
            // Not folded by this test's subject; named so the destructure stays exhaustive.
            reads_with_holed_witness: 0,
            hole_positions: 0,
            // ng's, and plain sums for the same reason as `reads_shed_at_admission`.
            reads_evicted_at_ceiling: 4,
            positions_short_of_cap: 7,
            short_of_cap_deficit: 11,
            // Maxes, like `active_reads_high_water` beside them.
            column_depth_high_water: 12,
            pending_mates_high_water: 3,
            // Plain sums; named so the destructure stays exhaustive.
            capped_column_reads_seen: 0,
            capped_column_reads_seen_placed_left: 0,
            capped_column_reads_kept: 0,
            capped_column_reads_kept_placed_left: 0,
        };

        counts.fold_region_walk(&summary, baseline);
        // A second region, whose allocator counters have grown further still.
        counts.fold_region_walk(
            &RunSummary {
                chain_allocations: 20,
                active_reads_high_water: 9,
                mate_lookup_evictions: 9,
                column_depth_high_water: 6,
                pending_mates_high_water: 8,
                ..summary
            },
            ChainIdAllocatorCounters {
                chain_allocations: 14,
                active_reads_high_water: 5,
                mate_lookup_evictions: 6,
                pending_mates_high_water: 3,
            },
        );

        assert_eq!(
            counts.chain_allocations,
            (14 - 10) + (20 - 14),
            "deltas, not the run-to-date totals — summing them is the triangular sum"
        );
        assert_eq!(
            counts.mate_lookup_evictions,
            (6 - 4) + (9 - 6),
            "the other counter the same trap applies to"
        );
        assert_eq!(
            counts.active_reads_high_water, 9,
            "a peak is the largest region's, never the sum"
        );
        assert_eq!(counts.reads_admitted, 6, "per-walk counters add");
        assert_eq!(counts.column_depth_truncations, 16);
        assert_eq!(
            counts.reads_shed_at_admission, 18,
            "refusals add per walk, like every other counter the walk owns itself"
        );
        assert_eq!(
            counts.reads_silent_over_footprint, 4,
            "each walk owns its active set, so the silent exits add like the walker's own \
             counters and not like the shared allocator's"
        );
        assert_eq!(counts.reads_evicted_at_ceiling, 8, "evictions add per walk");
        assert_eq!(counts.positions_short_of_cap, 14, "positions add per walk");
        assert_eq!(counts.short_of_cap_deficit, 22, "the deficit adds per walk");
        assert_eq!(
            counts.column_depth_high_water, 12,
            "the deepest column is a peak: 12 in the first walk, 6 in the second"
        );
        assert_eq!(
            counts.pending_mates_high_water, 8,
            "the pending-mate peak is a peak too — and it comes off the shared allocator, \
             so a *delta* would be the wrong fold for it"
        );
    }

    /// **A read that contributes nowhere is counted, not lost** — `reads_silent_over_footprint`
    /// (D2), the counter spec §13's read accounting cannot do without.
    ///
    /// A read every base of which the G1 adaptor filter silences is admitted, walked past and
    /// expired without ever appearing in a contributor list. **Both per-locus counters are
    /// blind to it:** it produced no observation, but it never reached the fold that records
    /// `reads_without_observation` either (spec §6), so before this counter the only honest
    /// thing to say about such a read was nothing at all.
    ///
    /// The neighbouring read is what makes the test discriminating: it proves the walk still
    /// emits loci here, so "no observations" cannot be mistaken for "no coverage".
    #[test]
    fn a_read_silent_at_every_position_is_counted_rather_than_lost() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("silent", 0, 10, 30),
            read_named_with_length("speaks", 0, 10, 30),
        ]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), SilencesRead("silent"))
                .expect("the default config is valid");

        let loci = loci_of(&mut generator, region(0, 1, 100), &reads);

        assert_eq!(
            generator.counts().reads_admitted,
            2,
            "both reads reached the walk — silencing happens inside it, not at admission"
        );
        assert_eq!(
            generator.counts().reads_silent_over_footprint,
            1,
            "the silenced read contributed at no position and must be counted as such"
        );
        assert!(!loci.is_empty(), "the speaking read still yields loci");
        for locus in &loci {
            assert_eq!(
                total_obs(locus),
                1,
                "only the speaking read supports {:?} — the silent one contributes nothing",
                locus.region,
            );
            assert_eq!(
                (
                    locus.reads_without_observation,
                    locus.reads_discarded_by_cap
                ),
                (0, 0),
                "and neither per-locus counter can see the silent read, which is why the \
                 run-level one has to exist",
            );
        }
    }

    /// **A read still active when the walk stops is judged on what it contributed, not on how it
    /// left** — the `ever_contributed` guard in `ActiveReads::flush_all` (D2).
    ///
    /// # Why this needs its own test, and how the review found that out
    ///
    /// `reads_silent_over_footprint` is fed by the active set's **two** exits, and only one of
    /// them was pinned. `expire_passed` is the ordinary one — a read whose `alignment_end` the
    /// walker has passed — and
    /// [`a_read_silent_at_every_position_is_counted_rather_than_lost`] covers it. `flush_all` is
    /// the other, and on the generic path it is not an edge case at all: **a region walk stops at
    /// `region.end` while the reads reaching into the halo are still active**, so every bounded
    /// walk ends by flushing reads that never expired.
    ///
    /// Milestone D's review deleted the guard — making every read that leaves through `flush_all`
    /// count as silent — and **the whole 2,724-test suite passed**. The counter would then have
    /// over-reported on every region of every real run, which is the failure mode this whole
    /// milestone exists to make impossible: a number nobody can see being wrong.
    ///
    /// **One silent read and *two* contributing ones, and the asymmetry is the point.** With one
    /// of each, the correct guard and a guard *inverted* to count the contributors both yield a
    /// total of 1, so the test could not tell them apart — which is this branch's recurring
    /// defect, and it was in the first draft of this very test. Two contributors make the two
    /// answers 1 and 2.
    #[test]
    fn a_read_still_active_when_the_walk_stops_is_counted_by_what_it_contributed() {
        // The reads span 10..=39; the region stops at 20, so the walk ends with all three of
        // them still in the active set — they leave through `flush_all`, never through
        // `expire_passed`.
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("silent", 0, 10, 30),
            read_named_with_length("speaks", 0, 10, 30),
            read_named_with_length("speaks_too", 0, 10, 30),
        ]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), SilencesRead("silent"))
                .expect("the default config is valid");

        let loci = loci_of(&mut generator, region(0, 1, 20), &reads);

        assert_eq!(
            anchors(&loci),
            (10..=20).collect::<Vec<u64>>(),
            "the walk must stop at the region's end with the reads still active — if it ran to \
             their ends instead they would expire, and this test would be covering the other \
             exit"
        );
        assert_eq!(
            generator.counts().reads_silent_over_footprint,
            1,
            "exactly the one read that contributed nowhere: the two speaking reads reached the \
             fold at eleven positions each and must not be counted merely because the walk \
             stopped under them, and `silent` must be counted even though it never expired. A \
             guard that counted the contributors instead would say 2 here"
        );
    }

    /// **A read the preparer declines is counted** — `reads_declined_by_preparer` (D2).
    ///
    /// It never reaches the walk, so `reads_admitted` cannot account for it and no per-locus
    /// counter ever hears of it. No v1 preparer declines anything, so this counter reads zero
    /// on every real run; [`DeclinesRead`] is what makes it move, because a counter nothing
    /// can move is indistinguishable from one wired to nothing.
    #[test]
    fn a_read_the_preparer_declines_is_counted_and_never_admitted() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("declined", 0, 10, 30),
            read_named_with_length("kept", 0, 10, 30),
        ]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), DeclinesRead("declined"))
                .expect("the default config is valid");

        let loci = loci_of(&mut generator, region(0, 1, 100), &reads);

        assert_eq!(
            generator.counts().reads_declined_by_preparer,
            1,
            "the declined read has to be accounted for somewhere, and this is where"
        );
        assert_eq!(
            generator.counts().reads_admitted,
            1,
            "a declined read never reaches the walk, so it is not admitted"
        );
        assert!(!loci.is_empty(), "the kept read still yields loci");
        for locus in &loci {
            assert_eq!(
                total_obs(locus),
                1,
                "only the kept read supports {:?}",
                locus.region,
            );
        }
    }

    /// **Each region folds its own declines once** — the `std::mem::take` in
    /// [`end_walk`](PileupGenerator::end_walk), which had no test.
    ///
    /// `ReadPreparation` outlives every walk, so a tally left in its cell is folded
    /// again at the next region's end. `end_walk` takes rather than reads, and says so
    /// in a comment — and **the whole 2,725-test suite stayed green with the `take`
    /// replaced by a plain read**, because the test above walks *one* region and one
    /// region cannot see a fold that repeats.
    ///
    /// That is the third counter on this branch whose triangular-summing guard nothing
    /// could see, after the allocator's `reset`/`summary` pair (Checkpoint C) and
    /// `flush_all`'s `ever_contributed` (Milestone D's review). In all three the guard
    /// was right and the comment above it named the trap; what was missing every time
    /// was a **second region**. Any future counter folded at a region boundary wants a
    /// test of this shape, not of the one above.
    ///
    /// **One region per contig, so the expected number depends on nothing else.** The
    /// halo (`region.end + max_record_span`) covers the whole of a 100 bp fixture
    /// contig, so a same-contig pair would make the count a function of the halo
    /// reaching the second region's reads. Two contigs make each region's declines its
    /// own: 2 on `chr1`, 1 on `chr2`.
    ///
    /// Mutation: `std::mem::take(&mut preparation.declined)` → `preparation.declined`
    /// reports **5** instead of 3 — `chr1`'s 2 folded twice.
    #[test]
    fn each_regions_declined_reads_are_folded_once_and_not_again_at_the_next_region() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            // chr1, coordinate-ordered: two declined and one kept.
            read_named_with_length("declined_1a", 0, 10, 30),
            read_named_with_length("kept_1", 0, 10, 30),
            read_named_with_length("declined_1b", 0, 40, 30),
            // chr2: one declined and one kept.
            read_named_with_length("declined_2", 1, 10, 30),
            read_named_with_length("kept_2", 1, 10, 30),
        ]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), DeclinesRead("declined"))
                .expect("the default config is valid");

        let chr1_loci = loci_of(&mut generator, region(0, 1, 100), &reads);
        let chr2_loci = loci_of(&mut generator, region(1, 1, 200), &reads);

        assert_eq!(
            generator.counts().reads_declined_by_preparer,
            3,
            "chr1 declines two and chr2 one; 5 means chr1's tally was folded again when \
             chr2's walk ended"
        );
        assert_eq!(
            generator.counts().reads_admitted,
            2,
            "one kept read per contig reached the walk — a zero here would make the \
             count above vacuous"
        );
        assert!(
            !chr1_loci.is_empty() && !chr2_loci.is_empty(),
            "both regions yielded loci, so both walks ran and both ends folded"
        );
    }

    /// **The query is consumed as the walk advances, never up front.** Spec §7
    /// makes this the property a port can quietly destroy: collecting the query
    /// into a `Vec` to make an ownership problem go away turns a depth-shaped
    /// footprint into a region-shaped one, and a `Generic` region runs to
    /// hundreds of kilobases. The review ablated the stream to a `collect()` and
    /// **the entire library suite stayed green**, parity included — so nothing
    /// was watching.
    ///
    /// The preparer is the seam that sees each read arrive, so counting there
    /// counts pulls. One read is in hand and one is peeked, so the first locus
    /// must have cost two of the five, not five.
    #[test]
    fn the_query_is_pulled_as_the_walk_advances_not_collected_up_front() {
        let records: Vec<RecordBuf> = [10, 50, 90, 130, 170]
            .iter()
            .map(|start| read_named_with_length(&format!("r{start}"), 1, *start, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = a_generator_with(
            PileupGeneratorConfig::default(),
            CountingPreparer::default(),
        )
        .expect("config");

        generator.begin_segment(region(1, 1, 200));
        generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .expect("the first read covers the region");

        let after_first = generator.preparation.borrow().preparer.prepared.get();
        assert!(
            after_first <= 2,
            "the first locus needs the read under the walker and the one peeked ahead, \
             not the region's whole read set; {after_first} of 5 were pulled"
        );

        while generator
            .next_locus(&reads)
            .expect("the walk succeeds")
            .is_some()
        {}
        assert_eq!(
            generator.preparation.borrow().preparer.prepared.get(),
            5,
            "and by the end every read has been pulled exactly once"
        );
    }

    /// A failed open must not strand the chain-id allocator. Since a fatal error
    /// now fuses the generator, nothing re-opens a walk to trip over it — so the
    /// ordering is pinned where it is decided rather than by its consequence
    /// (review).
    #[test]
    fn a_failed_query_leaves_the_allocator_with_the_generator() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        generator.begin_segment(region(9, 1, 100));
        assert!(generator.next_locus(&reads).is_err());

        assert!(
            generator.chain_ids.is_some(),
            "the allocator is taken only after the one fallible step of opening a walk"
        );
    }

    // --- the review's fixes: the generator ends, and stays ended -------------

    /// **A drained segment stays drained.** `next_locus` returning `None` is the
    /// caller's signal to move on; asking once more used to re-open the query
    /// and re-emit every locus of the region, with every read admitted a second
    /// time and its fragment handed a second chain id.
    #[test]
    fn a_drained_segment_yields_nothing_when_asked_again() {
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r", 0, 10, 30)]);
        let mut generator = a_generator(PileupGeneratorConfig::default()).expect("config");

        let loci = loci_of(&mut generator, region(0, 1, 100), &reads);
        assert_eq!(loci.len(), 30);
        let admitted = generator.counts().reads_admitted;

        assert!(
            generator
                .next_locus(&reads)
                .expect("a drained segment is not an error")
                .is_none(),
            "the region was already drained"
        );
        assert_eq!(
            generator.counts().reads_admitted,
            admitted,
            "nothing was walked a second time"
        );
    }

    /// **A failed run stays failed.** Every error here is terminal (spec §7);
    /// asking again used to re-open the region, re-emit its loci and re-issue
    /// its chain ids — the corruption the run-lifetime allocator exists to
    /// prevent, reached from the other direction.
    #[test]
    fn a_walk_that_failed_does_not_restart() {
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
        assert!(matches!(result, Err(LocusGenerationError::Walker { .. })));
        let admitted = generator.counts().reads_admitted;
        assert!(
            admitted > 0,
            "the erroring region's counters are folded too — `end_walk` runs on the error \
             path, and nothing begins another region after a fatal error to fold them later"
        );

        assert!(
            generator
                .next_locus(&reads)
                .expect("a failed run reports its error once, then yields nothing")
                .is_none()
        );
        assert_eq!(generator.counts().reads_admitted, admitted);
    }

    /// **A shed error belongs to the region that shed it.** The stream's errors
    /// surface a call after they happen, and the slot they wait in used to live
    /// on the generator rather than on the walk — so a region **abandoned**
    /// before it drained carried its failure into the next region, which
    /// reported it after emitting every one of its own loci.
    ///
    /// **The first version of this test could not fail for its stated reason**,
    /// and adding the region to the error is what exposed it: it called
    /// `begin_segment` twice in a row, and `begin_segment` opens nothing — so
    /// the "abandoned" region never ran, never shed anything, and the error it
    /// caught came from the *second* region's own read. It failed under the fix
    /// being reverted, which is why it looked sound, but what it pinned was "a
    /// shed error is reported at all", which another test already covers.
    ///
    /// The shape it needs: the failing read must be one the walker reaches while
    /// an earlier read is still producing loci, so the caller can stop mid-drain
    /// with the failure latched and unreported. `bad` on chr2 is that read;
    /// `good` gives the region loci to emit first.
    #[test]
    fn an_abandoned_regions_shed_error_is_not_charged_to_the_next_region() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("chr1_read", 0, 10, 30),
            read_named_with_length("good", 1, 10, 30),
            read_named_with_length("bad", 1, 50, 30),
        ]);
        let mut generator =
            a_generator_with(PileupGeneratorConfig::default(), FailsToPrepareRead("bad"))
                .expect("config");

        // chr2 emits at least one locus and sheds on `bad` — then is abandoned
        // with that failure still unreported.
        generator.begin_segment(region(1, 1, 100));
        let first = generator
            .next_locus(&reads)
            .expect("the walk succeeds up to the failing read")
            .expect("`good` covers the region");
        assert_eq!(first.region.start, Position(10));

        generator.begin_segment(region(0, 1, 100));

        let error = generator
            .next_locus(&reads)
            .expect_err("the abandoned region's failure is reported before the new region walks");
        assert!(matches!(error, LocusGenerationError::Reference { .. }));
        assert_eq!(
            error.region(),
            Some(region(1, 1, 100)),
            "the error belongs to chr2, whose read failed — not to the chr1 region current \
             when it surfaced"
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
