//! The observation cache — one forward reader per sample, and the window the builders read.
//!
//! Upstream produces a sample's observations in one forward pass; the builders want them in
//! several places at once (`doc/devel/ng/spec/cohort_merge.md` §6.4). Nothing can serve that
//! by seeking, and giving every builder its own reader would mean as many readers per sample
//! as there are builders. So there is **one reader per sample for the whole run**, advancing
//! forward only, and this is what sits between it and the builders: it draws the readers
//! forward until a region's ground is covered, hands the covering observations out for the
//! length of a call, and drops what nothing can reach any more.
//!
//! **A *region* here is a building region** — the stretch of genome assigned to one builder,
//! `cohort_locus_builder_regions_len` bases (spec §6.1). The run's own intervals are called
//! *analysed regions*, as in [`merge_cohort_serially`](super::serial::merge_cohort_serially),
//! and are much longer.
//!
//! **Why the merge needs this rather than merely liking it.** A builder closes loci from the
//! beginning of whatever observations it is handed and discards those that opened before its
//! own ground ([`build_region`](super::build::build_region)), so handing every builder the
//! whole analysed stretch costs each of them the whole prefix — about **3.3 µs per prefix
//! base at 63 samples**, measured in a release build by the C1 review
//! (`doc/devel/reports/reviews/ng_cohort_merge_c1_2026-08-17.md`). The same effect end to
//! end, measured on **one sample** by the C2 review
//! (`doc/devel/reports/reviews/ng_cohort_merge_c2_2026-08-17.md`): 20,000 observations cost
//! **5.4 ms merged as one analysed region and 184 ms as a thousand**. Short building regions
//! are only affordable when each builder is handed a window over its own ground, and that
//! window is what this file produces.
//!
//! **Builders read it and never write it** (spec §6.4, goal 1):
//! [`cover`](ObservationCache::cover) and [`evict_before`](ObservationCache::evict_before)
//! take `&mut self` and belong to the organiser, while
//! [`with_observations`](ObservationCache::with_observations) takes `&self` and is a builder's
//! only way in.
//!
//! **The organiser is the other half of this file** — [`Organiser`], which takes the builders'
//! outcomes and releases their loci in genome order. Ordered release has landed; resolving the
//! overlaps between neighbouring regions, and drawing the cache forward for the builders, are
//! the next two steps (`doc/devel/ng/impl_plan/cohort_merge.md`, E2 and E3).
//!
//! **Whether the two belong in one file is open, and the argument for it is weaker than it
//! first looked.** The case was that the organiser would become the cache's only writer, at
//! which point [`cover`](ObservationCache::cover) and
//! [`evict_before`](ObservationCache::evict_before) could turn from `pub` into private to this
//! file — a split enforced by the compiler rather than by convention. E1's review found the
//! premise false: `super::serial::merge_cohort_through_cache` calls both today, from a sibling
//! module, and nothing in the plan removes that driver. So the reachable narrowing is
//! `pub(super)`, which a file of its own would get just as well. Against that, this file is
//! about 1,950 lines covering two types that share no field and no function. The split is
//! recorded for the owner at Checkpoint E rather than taken here, because it would move D1 and
//! D2's code and the architecture's file tree names `organise.rs`.

use std::collections::{BTreeMap, VecDeque};
use std::iter::Fuse;

use super::CohortLocusBuilderRegionsLen;
use super::build::{CohortObservation, RegionOutcome};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion, Position};

/// Every sample's observations over the ground currently assigned to builders (spec §6.4).
///
/// **One forward reader per sample, and a window of what it has drawn.** The cache never
/// seeks and never goes back: the only ways the window moves are [`cover`](Self::cover),
/// which draws it forward, and [`evict_before`](Self::evict_before), which drops from its
/// left edge.
///
/// **This is the module's dominant memory** (spec §8), which is why building regions are
/// short: the ground it spans is `builders × cohort_locus_builder_regions_len` — 320 bases at
/// 16 builders on the default 20-base regions — plus the tail of observations reaching past
/// it.
///
/// **A source is any iterator of one sample's observations in coordinate order, and its error
/// type is its own.** The architecture writes `cover(…) -> Result<(), RunError>` and names
/// `run_streaming.md` arch §2's `ObservationSource` as what the cache is handed; neither
/// exists yet — both belong to the run's own document, which is out of this plan's scope — so
/// the cache is generic over the source and passes its failures through untouched. When
/// `ObservationSource` lands, its `observations_in` iterator is exactly this shape and
/// `RunError` is exactly this `E`.
///
/// **`E` must identify the sample it came from.** This cache knows which sample's reader
/// failed and passes the error on without adding to it, so a source whose error does not name
/// its own sample produces a failure an operator cannot act on at a cohort of thousands —
/// arch §5's two variants both carry `sample` for exactly that reason.
pub struct ObservationCache<S> {
    /// One per sample, in the run's sample order — the order every consumer indexes by.
    samples: Vec<SampleWindow<S>>,
    /// How far a **successful** [`cover`](Self::cover) has drawn, genome-wide.
    ///
    /// A cover that failed does not move it, which is what lets
    /// [`with_observations`](Self::with_observations) refuse ground no reader reached rather
    /// than hand out a window that is short — and short is a locus closed over the wrong
    /// ground, which is a wrong answer rather than a failure.
    covered_to: Option<GenomePosition>,
}

/// One sample's reader and the observations drawn from it that have not been evicted.
struct SampleWindow<S> {
    /// The forward reader. Never seeks, never rewinds.
    ///
    /// **Fused**, because an iterator is not required to be: a source that yielded `Some`
    /// after a `None` would be drawn in behind the window's own right edge, and so silently
    /// out of coordinate order. [`Fuse`] is the standard library's guard for exactly that,
    /// and using it means there is no flag of our own to set on the right branch — in
    /// particular, a failure is `Some(Err(_))` and so leaves the source live, which is what
    /// lets a cover be made again.
    source: Fuse<S>,
    /// What has been drawn and not yet evicted, in coordinate order.
    ///
    /// **A `Vec` and not a `VecDeque`**, because a builder is handed a contiguous slice of it
    /// ([`ObservationCache::with_observations`]) and a deque's two halves are not one.
    /// Eviction pays a move of what survives, which is the window — short by construction.
    held_observations: Vec<SampleLocusObservations>,
    /// Where the last observation drawn from `source` began — the ordering check's memory.
    last_drawn: Option<GenomePosition>,
}

impl<S: Iterator> ObservationCache<S> {
    /// A cache over one source per sample, in the run's sample order.
    ///
    /// Zero samples is not an error here: refusing a zero-sample *run* happens where the run
    /// is configured (spec §7.2), and a cache over an empty cohort covers nothing and hands
    /// nothing out.
    pub fn over(sources: Vec<S>) -> Self {
        Self {
            samples: sources
                .into_iter()
                .map(|source| SampleWindow {
                    source: source.fuse(),
                    held_observations: Vec::new(),
                    last_drawn: None,
                })
                .collect(),
            covered_to: None,
        }
    }
}

impl<S> ObservationCache<S> {
    /// Every sample's held observations from the first one that reaches into `span`, for the
    /// length of the call — a builder's only way in, and read-only by construction.
    ///
    /// The slices are handed out in the run's sample order, which is what
    /// [`build_region`](super::build::build_region) indexes by.
    ///
    /// **`span`'s two ends do different jobs, and neither is the obvious one.** Its start
    /// selects the window's left edge. Its end does **not** bound the right — trimming there
    /// would cut a locus that opens inside `span` and reaches past it, precisely the deletion
    /// the ownership rule exists to keep whole (spec §6.1) — but neither is it ignored: it is
    /// what this checks against the ground [`cover`](Self::cover) actually reached, so a
    /// builder asking for ground no reader has drawn is refused rather than handed a window
    /// one observation short.
    ///
    /// **Observations that end before `span` are skipped rather than evicted**, so a builder
    /// pays no prefix for ground the organiser has not released yet. Skipping them changes no
    /// locus the caller can own: a locus whose first position is at or after `span.start` has
    /// every member at or after it too, and a locus that chains back past `span.start` does so
    /// through an observation that reaches into `span` — which is kept, so the chain is seen
    /// and the locus still opens before `span`, and is still skipped by the ownership rule.
    pub fn with_observations<R>(
        &self,
        span: GenomeRegion,
        f: impl FnOnce(&[&[SampleLocusObservations]]) -> R,
    ) -> R {
        // `min`/`max` for the reason `cover` and `SampleLocusObservations::reach` use them:
        // `GenomeRegion` has public fields and no constructor ordering them, and these two
        // methods must read an inverted region the same way, or one draws ground the other
        // then trims away.
        let left_edge = GenomePosition {
            contig: span.contig,
            position: span.start.min(span.end),
        };
        let last_base = GenomePosition {
            contig: span.contig,
            position: span.start.max(span.end),
        };

        // **A release check, for the same reason the coordinate-order one is**: a window over
        // ground no cover reached is short, and short is a locus closed over the wrong ground
        // rather than a failure. What it says is that the ground was *drawn* — what is still
        // held is `evict_before`'s caller's business.
        assert!(
            self.covered_to
                .is_some_and(|covered_to| last_base <= covered_to),
            "the window was asked for {span} but cover has only reached {:?} — a cover that \
             failed, or was never made, leaves a window short of what a locus there reaches",
            self.covered_to,
        );

        // One slice reference per sample, one allocation per call. The closer copies them
        // again into its own array (`LocusCloser::over`); doing it twice costs nothing beside
        // the walk, and the alternative — a scratch buffer on the cache — would need interior
        // mutability on a `&self` method that several builders will hold at once.
        let observations_per_sample: Vec<&[SampleLocusObservations]> = self
            .samples
            .iter()
            .map(|sample| {
                let held = &sample.held_observations;
                &held[first_reaching_index(held, left_edge)..]
            })
            .collect();
        f(&observations_per_sample)
    }

    /// How many observations are held, summed across samples — the size of the window this
    /// cache is the memory of (spec §8).
    ///
    /// **It exists so that "eviction keeps up" can be asserted rather than assumed.** Nothing
    /// in a merge's output shows whether [`evict_before`](Self::evict_before) was ever called:
    /// a driver that never evicted would produce exactly the right answer and hold the whole
    /// stretch while doing it. This is what a test — and later the run's memory report — reads
    /// to tell those two apart.
    pub fn held_observations_len(&self) -> usize {
        self.samples
            .iter()
            .map(|sample| sample.held_observations.len())
            .sum()
    }

    /// Drop everything that ends before `position`. Called once the organiser has released
    /// every locus that could have started there.
    ///
    /// **An observation that began before `position` and reaches into it stays**, which is
    /// what makes eviction safe rather than merely cheap: that observation is what chains a
    /// locus across the evicted point, and dropping it would let a locus that really opens
    /// before `position` appear to open at it — and so be claimed and built by a builder that
    /// never saw its first bases.
    ///
    /// **It is a prefix drain, and on legal input that is a distinction without a
    /// difference.** A sample's records are disjoint and ascending — `build_region` asserts
    /// it — so reach is monotone across one sample's window and the first survivor is the last
    /// non-survivor's successor: nothing behind it could have been kept, nothing after it
    /// could have been dropped. The prefix form is chosen for its cost, which is proportional
    /// to what it drops rather than to the window, and not for a difference in what it keeps.
    pub fn evict_before(&mut self, position: GenomePosition) {
        for sample in &mut self.samples {
            let first_survivor = first_reaching_index(&sample.held_observations, position);
            sample.held_observations.drain(..first_survivor);
        }
    }
}

impl<S, E> ObservationCache<S>
where
    S: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    /// Draw every sample forward until `region` is covered, and far enough past it to hold
    /// what a locus starting inside it can reach (spec §6.4).
    ///
    /// **The second half is the whole difficulty, and it is a fixpoint across samples.** A
    /// locus opening inside `region` closes only when the next observation begins beyond its
    /// reach, and each observation the reach pulls in may push it further — so the chain's
    /// reach starts at `region`'s last base and grows with every observation that begins at or
    /// before it, which is the same `<=` the closer chains on (spec §4.1). One sample's
    /// deletion is what makes another sample's later observation part of the locus, so the
    /// samples are swept repeatedly until a whole sweep moves the reach no further. A single
    /// sweep stops short whenever the widening sample is swept after the sample it widens
    /// onto, and three samples can need three sweeps
    /// (`a_chain_that_needs_a_third_sweep_is_drawn_whole`): the sweep count is a property of
    /// the data, not of the code.
    ///
    /// **How far past `region` this can go is not bounded here.** The reach grows only through
    /// observations that chain into it, so what limits it is the widest observation the
    /// generator can mint — the reach ceiling, `max_record_span` (spec §1.3) — times the
    /// length of the chain. On ground where observations overlap wall to wall (spec §7.1) one
    /// cover can draw a whole segment. Nothing here reads or checks that ceiling: it is the
    /// generator's, and the psp header field that would carry it is deferred (spec §13).
    ///
    /// **What a cover costs.** Each sweep visits every sample and re-reads its held window, so
    /// a cover is `sweeps × (samples + held)`. Two sweeps is the ordinary count; the worst
    /// case is `samples` sweeps, when a chain of overlapping observations runs through the
    /// cohort in decreasing sample order. Measured by the D1 review on a synthetic cohort in a
    /// release build, over 20-base regions: one extra sweep costs 3.1 µs at 3,000 samples
    /// against 2.87 ms for the whole cover, and that worst case costs 28 ms for one 11-base
    /// region at 3,000 samples. The `held` term stays short **only while the organiser evicts
    /// at the pace it releases ground**, which this module cannot enforce — `evict_before` is
    /// the organiser's call (milestone E). At 1,000 samples the same walk costs 616 µs a cover
    /// with 4 observations held per sample and 1,028 µs with 200.
    ///
    /// **The window overshoots by at most one observation per sample**, and that is what a
    /// forward reader costs: the only way to know whether the next observation begins beyond
    /// the chain's reach is to draw it, and once drawn it is held rather than thrown away.
    ///
    /// Failures from a source end the cover and are passed through unchanged. The window keeps
    /// whatever was drawn before the failure and the cover can be made again — **which asks of
    /// a source that it may be polled after yielding `Err`**, something `Iterator` does not
    /// grant on its own. A source that cannot honour it must yield `None` after failing, which
    /// this reads as a spent sample rather than as a retry.
    pub fn cover(&mut self, region: GenomeRegion) -> Result<(), E> {
        // `max` for the same reason `SampleLocusObservations::reach` uses it: `GenomeRegion`
        // has public fields and no constructor enforcing `start <= end`, and an inverted
        // region must not put the chain's reach before the ground it is meant to cover.
        let mut chain_reach = GenomePosition {
            contig: region.contig,
            position: region.end.max(region.start),
        };

        // The fixpoint: sweep until a whole sweep moves nothing.
        while self.sweep(&mut chain_reach)? {}

        self.covered_to = Some(
            self.covered_to
                .map_or(chain_reach, |reached| reached.max(chain_reach)),
        );
        Ok(())
    }

    /// One sweep of every sample against `chain_reach`. Answers whether any of them moved it.
    fn sweep(&mut self, chain_reach: &mut GenomePosition) -> Result<bool, E> {
        let mut reach_grew = false;
        for sample in &mut self.samples {
            // Not `reach_grew |= sample.draw_to(…)?`: that reads as though the call could be
            // skipped, and a later `||` in its place would skip it.
            if sample.draw_to(chain_reach)? {
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }
}

impl<S, E> SampleWindow<S>
where
    S: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    /// Widen `chain_reach` with this sample's observations while they begin at or before it,
    /// drawing from the source as the held ones run out. Answers whether the reach moved.
    ///
    /// **It starts from the window's first observation every time, rather than remembering
    /// where the last sweep stopped.** Re-reading one that is already inside the reach cannot
    /// move it — the reach only ever grows — so the scan is idempotent, and keeping no mark
    /// means there is no mark for eviction to correct or for a later region to inherit: a
    /// cover is a function of the window, the source and the region it was asked for. A mark
    /// carried across covers is not merely stale but can point past the end of a window that
    /// lost two entries at once, which `a_survivor_of_an_eviction_still_widens_the_next_reach`
    /// and `two_evicted_at_once_leave_the_window_sound` exist to catch.
    fn draw_to(&mut self, chain_reach: &mut GenomePosition) -> Result<bool, E> {
        let mut reach_grew = false;
        let mut considered = 0;
        loop {
            if considered == self.held_observations.len() {
                // Nothing left in this sample: it cannot widen the reach again.
                let Some(observation) = self.draw_next()? else {
                    break;
                };
                self.held_observations.push(observation);
            }
            let observation = &self.held_observations[considered];
            if observation.start_position() > *chain_reach {
                // Held, but beyond the reach of any locus this cover can see. It stays in the
                // window, and the next cover reconsiders it against a later reach.
                break;
            }
            considered += 1;
            let reach = observation.reach_position();
            if reach > *chain_reach {
                *chain_reach = reach;
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }

    /// The next observation from the source, or `None` once it is spent.
    fn draw_next(&mut self) -> Result<Option<SampleLocusObservations>, E> {
        let Some(next) = self.source.next().transpose()? else {
            return Ok(None);
        };

        // **A release check, not a `debug_assert!`** — the release profile is the one this
        // repo runs, and a source that goes backwards is silent otherwise: the draw loop stops
        // at the first observation beyond the reach, so an out-of-order one would end a cover
        // early and hand a builder a locus cut short, which is a wrong answer rather than a
        // failure. The closer's own ordering check (`LocusCloser::over`) compares only within
        // one locus and would not see it. One comparison per observation drawn.
        //
        // **When observations are decoded from a psp file this must become a `RunError`**
        // beside `ObservationExceedsReachCeiling` (arch §5), like the producer-guarantee
        // checks in `build.rs` — a source out of coordinate order is then a fact about the
        // file rather than a bug in this crate, and this is the first such check the psp path
        // reaches.
        let start = next.start_position();
        if let Some(previous) = self.last_drawn {
            assert!(
                start >= previous,
                "this sample's source is not in coordinate order: {} follows {previous:?}",
                next.region,
            );
        }
        self.last_drawn = Some(start);

        Ok(Some(next))
    }
}

/// One analysed region divided into the regions single builders own — `building_region_width`
/// bases each, adjacent, in genome order, the last one clamped to the analysed region's own
/// last base.
///
/// **It lives here, beside the cache, because it is the organiser's geometry**: milestone E
/// hands these regions out to builders, and a second derivation of the clamp there is exactly
/// the defect this function's tests exist to catch. Its first caller is the serial driver that
/// reads through the cache (`super::serial::merge_cohort_through_cache`).
///
/// **It is a function of its own so that the division can be checked**, which no merge's output
/// can do: dividing the ground changes nothing a caller sees — that is the whole claim of the
/// cached driver — so a driver that ignored the width and built each analysed region as one
/// would give the right answer and lose everything the division buys.
///
/// **The clamp is load-bearing.** A locus belongs to the builder whose region holds its first
/// position, so a last building region running past the analysed ground would claim loci the
/// run never analysed.
///
/// An analysed region whose ends are the wrong way round is read as the ground it names, the
/// same defence `SampleLocusObservations::reach` makes — belt-and-braces, since the drivers
/// refuse such a region before they get here
/// (`super::serial::refuse_malformed_analysed_regions`).
pub fn building_regions_of(
    analysed_region: GenomeRegion,
    building_region_width: CohortLocusBuilderRegionsLen,
) -> impl Iterator<Item = GenomeRegion> {
    let first_base = analysed_region.start.min(analysed_region.end);
    let last_base = analysed_region.end.max(analysed_region.start);
    let bases_per_region = u64::from(building_region_width.get());

    let building_region_from = move |start: Position| GenomeRegion {
        contig: analysed_region.contig,
        start,
        // Saturating, then clamped: the arithmetic must not wrap at the top of the coordinate
        // space, and the region must not reach past the ground the run analysed. The `- 1`
        // cannot underflow — the width wraps a `NonZeroU32`.
        end: Position(
            start
                .0
                .saturating_add(bases_per_region - 1)
                .min(last_base.0),
        ),
    };

    std::iter::successors(Some(building_region_from(first_base)), move |previous| {
        // `checked_add`, because a region ending on the last base of the coordinate space has
        // no successor and the addition would wrap to zero — and with the saturating arithmetic
        // above, a wrap would divide the genome for ever rather than panic.
        let next_first_base = previous.end.0.checked_add(1)?;
        (next_first_base <= last_base.0).then(|| building_region_from(Position(next_first_base)))
    })
}

/// The first held observation that reaches `position` or beyond — the window's left edge for
/// a call at `position`, and the window's length when every one of them ends before it.
fn first_reaching_index(
    held_observations: &[SampleLocusObservations],
    position: GenomePosition,
) -> usize {
    held_observations
        .iter()
        .position(|observation| observation.reach_position() >= position)
        .unwrap_or(held_observations.len())
}

/// Where one building region falls in the order the run hands them out — 0 for the run's
/// first, 1 for its second, counting on across every analysed region of the run.
///
/// **It is the run's numbering, not a coordinate**, and that is the point: the organiser
/// releases on a gapless run of indexes, so it can tell "region 7 has not arrived yet" from
/// "region 7 found nothing", which two genome positions cannot say. The regions themselves
/// are [`building_regions_of`]'s, taken in the order that iterator yields them.
#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct RegionIndex(pub u64);

impl std::fmt::Display for RegionIndex {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}", self.0)
    }
}

/// A run that ended owing output — the release-level guard production keeps as
/// `WriterError::MissingChunks` (`var_calling/vcf_writer.rs:152-158`), and what arch §5 will
/// fold into `RunError` when that type exists.
///
/// **The whole of E1 is that a gap must be an error and not a truncation** (spec §6.3): the
/// output would simply stop early, and the failed-locus total would be short by everything the
/// lost regions refused — the one number that says the width bound is charging more than
/// expected (spec §3.3). Both are silent, and a partly-lost run answers exactly like a complete
/// one.
///
/// **An enum, because the two ways to end short can happen together and can happen apart, and
/// no combination of them means "nothing was lost".** Arch §5 sketches one struct with one
/// count, which production can afford because it emits from inside its own submit and has no
/// second step to forget; here taking the released loci is the caller's own call, so a run can
/// also end with loci released and never taken.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RunEndedShort {
    /// Regions the run handed out whose loci never reached the caller: the first index that
    /// stalled the drain, and how many regions from it onwards were never released. **The
    /// count covers both halves of a stall** — the region that never delivered, and every
    /// region behind it that delivered and could not be let out.
    #[error(
        "region {first_stalled} never delivered its result, so {regions} region(s) were \
         handed out and never released"
    )]
    RegionsNeverReleased {
        first_stalled: RegionIndex,
        regions: u64,
    },

    /// Loci released in order and never taken from the organiser. Nothing stalled; the caller
    /// stopped draining.
    #[error("{loci} released locus/loci were never taken from the organiser")]
    LociNeverDrained { loci: u64 },

    /// Both at once, which one count could not express.
    #[error(
        "region {first_stalled} never delivered its result, so {regions} region(s) were \
         handed out and never released, and {loci} released locus/loci were never taken"
    )]
    RegionsNeverReleasedAndLociNeverDrained {
        first_stalled: RegionIndex,
        regions: u64,
        loci: u64,
    },
}

/// Holds the builders' outcomes and releases their loci in genome order (spec §6.1, §6.3).
///
/// **Builders finish out of order and the consumer must see genome order**, so the outcomes
/// are keyed by region index and drained while the head equals the next expected index — the
/// reorder map production's VCF writer drains on `next_expected`, carried whole
/// (`var_calling/vcf_writer.rs:168-176`).
///
/// **Every region delivers exactly one outcome, including the empty ones** (spec §6.3). The
/// drain advances only along an unbroken run of indexes, so a region that found nothing must
/// still submit — an empty [`RegionOutcome`] — or every region behind it is held for ever,
/// and what those regions found is lost with nothing to say so. That is what
/// [`RunEndedShort`] refuses at the end of the run.
///
/// **Waiting for the predecessor is not merely about order.** A region's loci can only be
/// confirmed once the region before it has been resolved, because that is what says whether a
/// locus owned earlier already covers this region's ground (spec §6.3). That resolution is the
/// next step's — this one releases what it is given, in order, and drops nothing — and
/// [`release_regions_in_turn`](Self::release_regions_in_turn) is where it will go (the plan's
/// E2).
///
/// **It does not yet hold the observation cache**, which arch §4 gives it. The cache is drawn
/// forward and evicted by whoever hands regions out, and that is the parallel arrangement's
/// shape to settle (the plan's E3); until then the two live side by side in this file rather
/// than one inside the other.
#[derive(Debug, Default)]
pub struct Organiser {
    /// The next region index that may be released. Everything below it has been released
    /// already; nothing at or above it has.
    next_expected_region: RegionIndex,
    /// Outcomes that arrived before their turn, keyed by region index.
    held_outcomes: BTreeMap<RegionIndex, RegionOutcome>,
    /// Loci already released, in genome order, waiting for the caller to take them. They are
    /// held rather than handed straight on because releasing and taking are two different
    /// moments: a region becomes releasable when its predecessor arrives, which is not when
    /// the caller next asks.
    released_loci: VecDeque<CohortObservation>,
    /// Summed over the regions **released so far**, not over those submitted — the number
    /// spec §3.3 requires to reach the run summary. Counting at release rather than at
    /// submission is what lets the next step drop a failed locus that an earlier one displaced
    /// without the total having counted it already.
    failed_locus_count: u64,
}

impl Organiser {
    /// An organiser expecting region index 0 and holding nothing.
    ///
    /// Written out field by field rather than deferring to `Default`, so that a field E2 or E3
    /// adds has to be given a starting value here instead of taking whatever its own `Default`
    /// happens to be.
    #[must_use]
    pub fn new() -> Self {
        Self {
            next_expected_region: RegionIndex(0),
            held_outcomes: BTreeMap::new(),
            released_loci: VecDeque::new(),
            failed_locus_count: 0,
        }
    }

    /// Take one region's outcome, and release everything the arrival makes releasable.
    ///
    /// **Releasing here rather than in [`drain_ready`](Self::drain_ready) is production's
    /// shape** (`var_calling/vcf_writer.rs:246`): a region becomes releasable the moment the
    /// index in front of it arrives, which has nothing to do with when the caller next asks
    /// for loci. Keeping the two apart is also what makes
    /// [`failed_locus_count`](Self::failed_locus_count) meaningful before anything has been
    /// drained.
    ///
    /// **Panics** when `index` was submitted already or has been released. Both are bugs in
    /// whoever hands the regions out rather than facts about the data, and both are caught
    /// mid-flight, where the release order is already wrong and nothing coherent can follow —
    /// which is what separates them from [`RunEndedShort`], a bug too, but one caught at the
    /// end of a run where what was lost can still be named and reported.
    pub fn submit(&mut self, index: RegionIndex, outcome: RegionOutcome) {
        assert!(
            index >= self.next_expected_region,
            "the outcome of region {index} arrived after it was released (next expected {}), \
             so its loci would be released a second time, out of order",
            self.next_expected_region,
        );
        // **Checked before the insert, not inside it.** Writing the map inside the assertion's
        // condition would put the module's only insertion in an expression a later edit to
        // `debug_assert!` would stop evaluating — and this crate's release profile leaves debug
        // assertions off, so every region's outcome would be dropped in the shipped binary and
        // in no test. Checking first also keeps the first outcome, rather than replacing it
        // and then panicking.
        assert!(
            !self.held_outcomes.contains_key(&index),
            "region {index} delivered a second outcome, and a region owns its loci exactly once",
        );
        self.held_outcomes.insert(index, outcome);
        self.release_regions_in_turn();
    }

    /// Everything released and not yet taken, in genome order.
    ///
    /// **What is left when the iterator is dropped early stays here**, to be handed out by the
    /// next call: this takes one locus at a time from the front rather than emptying the buffer
    /// into the iterator, so a caller that stops halfway loses nothing. Loci still here at the
    /// end of a run are what [`finish`](Self::finish) refuses.
    #[must_use = "the loci stay in the organiser until the iterator is consumed"]
    pub fn drain_ready(&mut self) -> impl Iterator<Item = CohortObservation> + '_ {
        let released_loci = &mut self.released_loci;
        std::iter::from_fn(move || released_loci.pop_front())
    }

    /// The failed loci of every region released so far, summed (spec §3.3).
    #[must_use]
    pub fn failed_locus_count(&self) -> u64 {
        self.failed_locus_count
    }

    /// Nothing outstanding: every region the run handed out has been released, and every
    /// released locus taken (arch §4).
    ///
    /// **`regions_handed_out` is how many building regions the run dealt out**, their indexes
    /// being `0..regions_handed_out`. Without it the organiser cannot see a gap at the *tail*
    /// of a run — indexes that never submitted, with no later index behind them to hold — and
    /// a run that lost its last regions would look exactly like one that finished.
    #[must_use]
    pub fn is_finished(&self, regions_handed_out: u64) -> bool {
        self.next_expected_region.0 >= regions_handed_out
            && self.held_outcomes.is_empty()
            && self.released_loci.is_empty()
    }

    /// End the run: the failed-locus total, or a refusal naming what would have been lost.
    ///
    /// Consuming, because there is nothing to ask an organiser afterwards — production's
    /// `VcfWriter::finish` has the same shape and the same reason
    /// (`var_calling/vcf_writer.rs:256`). `regions_handed_out` is
    /// [`is_finished`](Self::is_finished)'s, and this returns `Ok` on exactly the runs that
    /// method calls finished.
    ///
    /// **Panics** when an outcome was submitted for an index the run says it never handed out,
    /// which is the same class of hand-out bug [`submit`](Self::submit) refuses.
    pub fn finish(self, regions_handed_out: u64) -> Result<u64, RunEndedShort> {
        // Destructured rather than field-accessed, so that anything the organiser gains at E2
        // or E3 has to be answered for here — drained by the end of the run, or deliberately
        // not — instead of being left behind in silence.
        let Self {
            next_expected_region,
            held_outcomes,
            released_loci,
            failed_locus_count,
        } = self;

        if let Some((last_held, _)) = held_outcomes.last_key_value() {
            assert!(
                last_held.0 < regions_handed_out,
                "region {last_held} delivered an outcome though the run handed out only \
                 {regions_handed_out} region(s)",
            );
        }

        // Every index from the cursor to the last one handed out: the region that stalled the
        // drain, and every region behind it, whether it delivered or not. `saturating_sub`
        // because a caller that under-reports its own hand-out is the panic above, not a wrap.
        let regions = regions_handed_out.saturating_sub(next_expected_region.0);
        let loci = released_loci.len() as u64;

        match (regions, loci) {
            (0, 0) => Ok(failed_locus_count),
            (0, loci) => Err(RunEndedShort::LociNeverDrained { loci }),
            (regions, 0) => Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: next_expected_region,
                regions,
            }),
            (regions, loci) => Err(RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
                first_stalled: next_expected_region,
                regions,
                loci,
            }),
        }
    }

    /// Release the unbroken run of regions now at the head, oldest first.
    fn release_regions_in_turn(&mut self) {
        while let Some(outcome) = self.held_outcomes.remove(&self.next_expected_region) {
            // Destructured for the reason `build.rs` gives where it consumes a comparable
            // type: anything `RegionOutcome` gains has to be answered for here — carried or
            // dropped deliberately — rather than vanishing.
            let RegionOutcome {
                cohort_observations,
                failed_locus_spans,
            } = outcome;
            self.released_loci.extend(cohort_observations);
            // Saturating where the cursor below is checked, and the difference is deliberate:
            // a saturated total is a truer answer than a wrap for a count, while a cursor that
            // stopped advancing would release the same region for ever.
            self.failed_locus_count = self
                .failed_locus_count
                .saturating_add(failed_locus_spans.len() as u64);
            // PANIC-FREE: the cursor rises once per released region, and a run cannot hand out
            // more than `u64::MAX` building regions.
            self.next_expected_region = RegionIndex(
                self.next_expected_region
                    .0
                    .checked_add(1)
                    .expect("a run cannot hand out more than u64::MAX building regions"),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{SourceFailed, position_on, region, region_on};
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, Position, ReadGroupId};
    use std::cell::Cell;
    use std::rc::Rc;

    /// One sample's record over `region`. What it showed is irrelevant to the cache — which
    /// reads nothing but where an observation begins and how far it reaches — so every fixture
    /// observation carries the same single sighting. It is a well-formed one, because the
    /// differential test at the end of this file hands these to a real builder.
    fn observation_over(region: GenomeRegion) -> SampleLocusObservations {
        let width = usize::try_from(region.end.0.saturating_sub(region.start.0) + 1)
            .expect("a fixture region fits in memory");
        SampleLocusObservations {
            region,
            reference_bases: vec![b'A'; width].into_boxed_slice(),
            observations: vec![SequenceObservation {
                bases: Box::from(&b"C"[..]),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 3,
                q_sum: -6.0,
                mapq_sum: 180,
                mapq_sum_sq: 10_800,
                placed_left: 0,
                chain_ids: vec![1, 2, 3],
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// One sample's reader, counting what was drawn out of it. The count is what says whether
    /// the cache stopped drawing where it claims to.
    struct CountingSource {
        remaining: std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>>,
        drawn: Rc<Cell<usize>>,
    }

    impl Iterator for CountingSource {
        type Item = Result<SampleLocusObservations, SourceFailed>;

        fn next(&mut self) -> Option<Self::Item> {
            let next = self.remaining.next();
            if next.is_some() {
                self.drawn.set(self.drawn.get() + 1);
            }
            next
        }
    }

    fn width(bases: u32) -> CohortLocusBuilderRegionsLen {
        CohortLocusBuilderRegionsLen(
            std::num::NonZeroU32::new(bases).expect("a fixture width is non-zero"),
        )
    }

    /// A source over the given regions, with the counter the test reads afterwards.
    fn source_over(regions: &[GenomeRegion]) -> (CountingSource, Rc<Cell<usize>>) {
        let drawn = Rc::new(Cell::new(0));
        let items: Vec<_> = regions
            .iter()
            .map(|region| Ok(observation_over(*region)))
            .collect();
        (
            CountingSource {
                remaining: items.into_iter(),
                drawn: Rc::clone(&drawn),
            },
            drawn,
        )
    }

    /// The regions each sample's window hands out over `span`, in sample order.
    fn handed_out<S>(cache: &ObservationCache<S>, span: GenomeRegion) -> Vec<Vec<GenomeRegion>> {
        cache.with_observations(span, |per_sample| {
            per_sample
                .iter()
                .map(|observations| {
                    observations
                        .iter()
                        .map(|observation| observation.region)
                        .collect()
                })
                .collect()
        })
    }

    /// The plain case: covering a region draws the observations over it, and they come back.
    #[test]
    fn covering_a_region_draws_the_observations_over_it() {
        let (source, _drawn) = source_over(&[region(45, 45), region(48, 48)]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(45, 45), region(48, 48)]],
        );
    }

    /// **A window covering a region also holds an observation that started before it and
    /// reaches in** — the plan's own D1 case. The deletion opening at 30 is what a locus over
    /// this ground chains back to, so a cache that handed out only what begins inside the
    /// region would let a builder claim a locus whose first bases it never saw.
    #[test]
    fn an_observation_that_began_earlier_and_reaches_in_is_handed_out() {
        let (source, _drawn) = source_over(&[region(30, 44), region(45, 45)]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(30, 44), region(45, 45)]],
            "30–44 ends inside the region, so it is part of what a locus there can chain to",
        );
    }

    /// And one that ends before the span is not handed out, though the forward reader had to
    /// pass over it. That is the prefix a builder does not pay for.
    ///
    /// **The two assertions sit either side of the one base that decides it.** 20–39 ends one
    /// base before a span opening at 40 and is trimmed away; asked again over a span opening
    /// on 39 — its own last base — it comes back, because a locus opening there chains into
    /// it. An edge drawn one base further right would lose the locus's true first position.
    #[test]
    fn an_observation_that_ends_before_the_span_is_drawn_but_not_handed_out() {
        let (source, drawn) = source_over(&[region(20, 39), region(41, 45), region(46, 46)]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(41, 45), region(46, 46)]],
            "20–39 ends one base before the span and can reach no locus inside it",
        );
        assert_eq!(
            handed_out(&cache, region(39, 50)),
            vec![vec![region(20, 39), region(41, 45), region(46, 46)]],
            "and it comes back the moment the span opens on its own last base",
        );
        assert_eq!(
            drawn.get(),
            3,
            "it was still read: a forward reader passes over it to reach 41",
        );
    }

    /// **The chain's reach follows a widening in a later sample.** Sample 1's deletion, from
    /// 48 to 70, is what makes sample 0's observations part of a locus opening inside the
    /// region — and only a second sweep, after that deletion has moved the reach, draws them.
    ///
    /// **The second of sample 0's records is what makes the sweep visible.** One observation
    /// beyond the reach is held whatever the cover does, because the draw that discovers it
    /// keeps it; what a single sweep loses is the observation *after* it. Measured by the D1
    /// review: with only the first, this test passes under a single-sweep cover.
    #[test]
    fn the_chain_reach_follows_a_widening_in_a_later_sample() {
        let (far, far_drawn) = source_over(&[region(60, 60), region(65, 65)]);
        let (widening, _widening_drawn) = source_over(&[region(48, 70)]);
        let mut cache = ObservationCache::over(vec![far, widening]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(far_drawn.get(), 2, "a single sweep would draw only 60");
        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![vec![region(60, 60), region(65, 65)], vec![region(48, 70)]],
            "the deletion reaching to 70 is what draws 65 into the window; 60 was already \
             kept by the draw that discovered it",
        );
    }

    /// **The fixpoint needs a third sweep here, and a cover that stopped after two would hand
    /// out a window one observation short.** Sample 2's deletion carries the reach to 70; only
    /// then is sample 1's 60–120 inside it, which carries it to 120; only then is sample 0's
    /// 100–300 inside it, which carries it to 300 — and 260 is in the window only because of
    /// that third widening.
    ///
    /// Found by the D1 review: capping the sweeps at two passed every other test in this file,
    /// and disagreed with a whole-stretch builder on 410 of 600 random layouts.
    #[test]
    fn a_chain_that_needs_a_third_sweep_is_drawn_whole() {
        let (far, _far_drawn) = source_over(&[region(100, 300)]);
        let (middle, _middle_drawn) =
            source_over(&[region(60, 120), region(250, 250), region(260, 260)]);
        let (opening, _opening_drawn) = source_over(&[region(50, 70)]);
        let mut cache = ObservationCache::over(vec![far, middle, opening]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![
                vec![region(100, 300)],
                vec![region(60, 120), region(250, 250), region(260, 260)],
                vec![region(50, 70)],
            ],
            "260 is in the window only because the third sweep carried the reach to 300",
        );
    }

    /// **An observation beginning on the reach's own base is inside it, and what it reaches is
    /// drawn too.** This is the same `<=` the closer chains on (spec §4.1), and it is the
    /// boundary that decides whether the window is complete: the third sample's second
    /// observation is in the window only because sample 0's record — which begins exactly
    /// where sample 1's deletion ends — carried the reach out to 200.
    ///
    /// **The observation at the boundary alone would not show it.** The one draw that
    /// discovers an observation is beyond the reach leaves it in the window either way, so a
    /// cache that stopped at the boundary would still hand out sample 0's record; what it
    /// would not hold is what that record reaches, which is why the late sample carries two.
    #[test]
    fn an_observation_beginning_on_the_reach_is_drawn_and_widens_it() {
        let (widened, widened_drawn) = source_over(&[region(70, 200), region(500, 500)]);
        let (widening, _widening_drawn) = source_over(&[region(48, 70)]);
        let (late, late_drawn) = source_over(&[region(150, 150), region(160, 160)]);
        let mut cache = ObservationCache::over(vec![widened, widening, late]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![
                vec![region(70, 200), region(500, 500)],
                vec![region(48, 70)],
                vec![region(150, 150), region(160, 160)],
            ],
            "70 begins on the reach the deletion left at 70, so its own reach to 200 is what \
             puts 160 in the window",
        );
        assert_eq!(widened_drawn.get(), 2, "500 is the one draw past the chain");
        assert_eq!(
            late_drawn.get(),
            2,
            "both of the late sample's observations were drawn, not merely the first",
        );
    }

    /// **The window overshoots by exactly one observation per sample, and stops.** The
    /// observation at 60 is beyond the reach and had to be drawn to find that out; the one at
    /// 80 was never read at all.
    #[test]
    fn drawing_stops_one_observation_past_the_reach() {
        let (source, drawn) = source_over(&[region(45, 45), region(60, 60), region(80, 80)]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            2,
            "45 is inside the region and 60 is the one draw that says the chain has ended",
        );
        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(45, 45), region(60, 60)]],
            "the overshoot is held rather than thrown away, and a later cover reconsiders it",
        );
    }

    /// A second cover carries on from where the first stopped, and re-reads nothing.
    #[test]
    fn a_second_cover_carries_on_without_re_reading() {
        let (source, drawn) = source_over(&[region(45, 45), region(60, 60), region(80, 80)]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "the second cover drew only 80 — 60 was already held from the first",
        );
        assert_eq!(
            handed_out(&cache, region(51, 70)),
            vec![vec![region(60, 60), region(80, 80)]],
            "and 45 is behind the span, so the second builder is not handed it",
        );
    }

    /// **Eviction drops nothing a live region can still reach** — the plan's second D1 case.
    /// The observation at 20–21 cannot touch the ground from 40 on and goes.
    ///
    /// **The second observation ends exactly on the evicted base**, which is the case the rule
    /// turns on: 25–40 still covers position 40, so a locus opening there chains back to 25,
    /// and a cache that dropped it would let that locus appear to open at 40 — claimed and
    /// built by a builder that never saw its first fifteen bases.
    ///
    /// The window is asked over ground reaching back to base 1, so a survivor that eviction
    /// had merely hidden would show here.
    #[test]
    fn eviction_keeps_what_reaches_past_the_evicted_point() {
        let (source, _drawn) = source_over(&[
            region(20, 21),
            region(25, 40),
            region(41, 44),
            region(45, 45),
        ]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(0, 40));

        assert_eq!(
            handed_out(&cache, region(1, 50)),
            vec![vec![region(25, 40), region(41, 44), region(45, 45)]],
            "20–21 is gone; 25–40 ends on the evicted base and stays",
        );
    }

    /// **Every sample's window is evicted, not merely the first.** Both samples here have a
    /// record behind the evicted point and one ahead of it, so a drop that reached only sample
    /// 0 would leave sample 1's 22–23 in the window — and the cache's whole memory bound is
    /// that this loop runs to the end (spec §8). Nothing in the output would show it.
    #[test]
    fn eviction_drops_from_every_sample_not_only_the_first() {
        let (first, _first_drawn) = source_over(&[region(20, 21), region(45, 45)]);
        let (second, _second_drawn) = source_over(&[region(22, 23), region(46, 46)]);
        let mut cache = ObservationCache::over(vec![first, second]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(0, 40));

        assert_eq!(
            handed_out(&cache, region(1, 50)),
            vec![vec![region(45, 45)], vec![region(46, 46)]],
            "both samples lost the record that ends before 40",
        );
    }

    /// **A contig boundary is past every position on the contig before it.** The held record
    /// at 900 on contig 0 is a much higher *number* than the position evicted at, 10 on contig
    /// 1, and it still goes: nothing on contig 1 can chain back to it.
    #[test]
    fn eviction_at_a_later_contig_drops_the_previous_contigs_window() {
        let (source, _drawn) = source_over(&[region_on(0, 900, 900), region_on(1, 45, 45)]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region_on(1, 40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(1, 10));

        assert_eq!(
            handed_out(&cache, region_on(0, 1, 999)),
            vec![vec![region_on(1, 45, 45)]],
            "contig 0's record is behind the evicted point though 900 is ahead of 10",
        );
    }

    /// Eviction leaves the window usable: the next cover draws forward from where it stood,
    /// and what survived is still handed out.
    #[test]
    fn a_cover_after_an_eviction_draws_forward_and_keeps_the_survivor() {
        let (source, drawn) = source_over(&[region(20, 21), region(30, 44), region(60, 60)]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(drawn.get(), 3, "everything the fixture holds, and no more");
        assert_eq!(
            handed_out(&cache, region(1, 70)),
            vec![vec![region(30, 44), region(60, 60)]],
            "the survivor of the eviction is still there, and the new ground was drawn",
        );
    }

    /// **The window's own observations are re-read at every cover, so an eviction leaves no
    /// mark to correct.** After the eviction, 30–80 is still what a locus in 51–70 chains
    /// through, so it must widen the reach again — a cover that remembered where the last
    /// sweep stopped would skip it, stop at 70, and draw one observation too far.
    #[test]
    fn a_survivor_of_an_eviction_still_widens_the_next_reach() {
        let (source, drawn) = source_over(&[
            region(20, 21),
            region(30, 80),
            region(90, 90),
            region(200, 200),
        ]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "30–80 reaches to 80, and 90 is the draw that ends the chain there",
        );
        assert_eq!(
            handed_out(&cache, region(1, 80)),
            vec![vec![region(30, 80), region(90, 90)]],
        );
    }

    /// Two observations evicted at once, with more still held — the case where an index
    /// remembered across covers would be not merely stale but out of range.
    #[test]
    fn two_evicted_at_once_leave_the_window_sound() {
        let (source, _drawn) = source_over(&[
            region(10, 11),
            region(20, 21),
            region(30, 80),
            region(100, 150),
            region(160, 160),
            region(400, 400),
        ]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(1, 80)),
            vec![vec![region(30, 80), region(100, 150)]],
            "both of the evicted pair are gone and the window behind them is intact",
        );
    }

    /// **A region whose ends are the wrong way round covers the ground it names**, rather than
    /// a reach before its own first base. `GenomeRegion` has public fields and no constructor
    /// ordering them, and a reach at the lower end would end the cover before the ground was
    /// drawn — a window short of what a locus there reaches, which is a wrong answer and not a
    /// failure. The same defence `SampleLocusObservations::reach` makes, and the second
    /// assertion is what says `with_observations` reads an inverted span the same way.
    #[test]
    fn a_region_given_end_first_still_covers_its_ground() {
        let (source, drawn) = source_over(&[region(45, 45), region(60, 60)]);
        let mut cache = ObservationCache::over(vec![source]);
        let inverted = GenomeRegion {
            contig: ContigId(0),
            start: Position(50),
            end: Position(40),
        };

        cache.cover(inverted).expect("the fixture source holds");

        assert_eq!(drawn.get(), 2, "45 was drawn, and 60 ended the chain");
        assert_eq!(
            handed_out(&cache, inverted),
            vec![vec![region(45, 45), region(60, 60)]],
            "asked back over the same inverted region, the window is the ground it names",
        );
    }

    /// A contig boundary is beyond every reach on the contig before it, whatever the positions
    /// say: covering a region on contig 1 draws past contig 0's observations and stops at
    /// contig 2's.
    #[test]
    fn a_cover_stops_at_the_next_contig() {
        let (source, drawn) = source_over(&[
            region_on(0, 90, 90),
            region_on(1, 45, 45),
            region_on(2, 10, 10),
            region_on(2, 20, 20),
        ]);
        let mut cache = ObservationCache::over(vec![source]);

        cache
            .cover(region_on(1, 40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "contig 0's observation, contig 1's, and the one draw that ends the chain",
        );
        assert_eq!(
            handed_out(&cache, region_on(1, 40, 50)),
            vec![vec![region_on(1, 45, 45), region_on(2, 10, 10)]],
            "contig 0's observation is behind the span even though position 90 is ahead of 40",
        );
    }

    /// A source's failure ends the cover and comes back unchanged.
    #[test]
    fn a_failing_source_ends_the_cover() {
        let drawn = Rc::new(Cell::new(0));
        let source = CountingSource {
            remaining: vec![
                Ok(observation_over(region(45, 45))),
                Err(SourceFailed("the block would not decode")),
            ]
            .into_iter(),
            drawn: Rc::clone(&drawn),
        };
        let mut cache = ObservationCache::over(vec![source]);

        let outcome = cache.cover(region(40, 50));

        assert_eq!(outcome, Err(SourceFailed("the block would not decode")));
        assert_eq!(drawn.get(), 2, "it failed on its second observation");
    }

    /// **A window over ground no cover reached is refused**, which is the state a failed cover
    /// leaves behind. Handing one out instead would give a builder a window short of what a
    /// locus there reaches — a locus closed over the wrong ground, and no failure anywhere.
    #[test]
    #[should_panic(expected = "cover has only reached")]
    fn a_window_over_ground_no_cover_reached_is_refused() {
        let (source, _drawn) = source_over(&[region(45, 45), region(120, 120)]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        let _ = handed_out(&cache, region(100, 130));
    }

    /// And a window is refused for reaching **past** covered ground, not merely for starting
    /// past it: the span below opens at 45, which was drawn, and ends at 130, which was not.
    /// That is the shape a builder handed too wide a region takes, and the loci it would lose
    /// are the ones in the ground the reader never reached.
    #[test]
    #[should_panic(expected = "cover has only reached")]
    fn a_window_reaching_past_the_covered_ground_is_refused() {
        let (source, _drawn) = source_over(&[region(45, 45), region(120, 120)]);
        let mut cache = ObservationCache::over(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        let _ = handed_out(&cache, region(45, 130));
    }

    /// **A cover that failed can be made again**, and the sample whose source failed is drawn
    /// on from where it stopped — which is why nothing latches on `Err`. The reach is rebuilt
    /// by re-reading the held window, so the retry still reaches 70 without re-reading the
    /// source, and what was drawn before the failure is still there.
    #[test]
    fn a_cover_can_be_made_again_after_a_failure() {
        let (steady, _steady_drawn) = source_over(&[region(45, 45)]);
        let recovering = CountingSource {
            remaining: vec![
                Ok(observation_over(region(48, 70))),
                Err(SourceFailed("the block would not decode")),
                Ok(observation_over(region(60, 60))),
            ]
            .into_iter(),
            drawn: Rc::new(Cell::new(0)),
        };
        let mut cache = ObservationCache::over(vec![steady, recovering]);

        let first = cache.cover(region(40, 50));
        assert_eq!(first, Err(SourceFailed("the block would not decode")));

        cache.cover(region(40, 50)).expect("the retry holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(45, 45)], vec![region(48, 70), region(60, 60)]],
            "the failed sample was drawn on, and the re-read window still reached 70",
        );
    }

    /// A spent source is not drawn from again. An iterator is free to yield `Some` after a
    /// `None`, and a cache that kept asking would draw observations behind the window's own
    /// reach, and so silently out of order.
    #[test]
    fn a_spent_source_is_not_drawn_from_again() {
        struct Resurrecting {
            polls: Rc<Cell<usize>>,
        }

        impl Iterator for Resurrecting {
            type Item = Result<SampleLocusObservations, SourceFailed>;

            fn next(&mut self) -> Option<Self::Item> {
                self.polls.set(self.polls.get() + 1);
                match self.polls.get() {
                    1 => Some(Ok(observation_over(region(45, 45)))),
                    2 => None,
                    _ => Some(Ok(observation_over(region(10, 10)))),
                }
            }
        }

        let polls = Rc::new(Cell::new(0));
        let mut cache = ObservationCache::over(vec![Resurrecting {
            polls: Rc::clone(&polls),
        }]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache
            .cover(region(51, 60))
            .expect("the fixture source holds");

        assert_eq!(
            polls.get(),
            2,
            "the second cover did not poll a spent source"
        );
        assert_eq!(
            handed_out(&cache, region(1, 60)),
            vec![vec![region(45, 45)]]
        );
    }

    /// A source whose observations go backwards is refused rather than quietly ending the
    /// cover where the walk cannot see it.
    #[test]
    #[should_panic(expected = "not in coordinate order")]
    fn a_source_that_goes_backwards_is_refused() {
        let (source, _drawn) = source_over(&[region(45, 45), region(30, 30)]);
        let mut cache = ObservationCache::over(vec![source]);

        let _ = cache.cover(region(40, 50));
    }

    /// And one that goes back to an **earlier contig** is refused too, though its positions
    /// rise. `GenomePosition`'s order is contig first, so this is the same check — but a
    /// source merging per-contig files is where it actually fires, and a comparison made only
    /// within a contig would let it through.
    #[test]
    #[should_panic(expected = "not in coordinate order")]
    fn a_source_that_goes_back_a_contig_is_refused() {
        let (source, _drawn) = source_over(&[region_on(1, 45, 45), region_on(0, 90, 90)]);
        let mut cache = ObservationCache::over(vec![source]);

        let _ = cache.cover(region_on(1, 40, 50));
    }

    /// A cohort of no samples covers nothing and hands out nothing, rather than failing:
    /// refusing a zero-sample run happens where the run is configured (spec §7.2).
    #[test]
    fn a_cache_over_no_samples_covers_nothing() {
        let mut cache: ObservationCache<CountingSource> = ObservationCache::over(Vec::new());

        cache.cover(region(40, 50)).expect("nothing can fail");

        assert!(handed_out(&cache, region(40, 50)).is_empty());
    }

    /// A sample that has run out does not stop the others: the cohort's ground is still
    /// covered, and the spent sample simply hands out nothing.
    #[test]
    fn a_sample_that_has_run_out_does_not_stop_the_cover() {
        let (spent, _spent_drawn) = source_over(&[region(10, 10)]);
        let (covering, _covering_drawn) = source_over(&[region(45, 45)]);
        let mut cache = ObservationCache::over(vec![spent, covering]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![Vec::new(), vec![region(45, 45)]],
        );
    }

    /// A sample with no records at all is not a special case either — an ordinary thing at 63
    /// samples, where one accession may have nothing on a contig.
    #[test]
    fn covering_a_region_with_an_empty_source_hands_out_nothing() {
        let (empty, empty_drawn) = source_over(&[]);
        let (covering, _covering_drawn) = source_over(&[region(45, 45)]);
        let mut cache = ObservationCache::over(vec![empty, covering]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        assert_eq!(
            empty_drawn.get(),
            0,
            "an empty source yields nothing to count",
        );
        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![Vec::new(), vec![region(45, 45)]],
        );
    }

    /// **The held count is the whole cohort's, not the first sample's**, and it is zero before
    /// anything is drawn. A count over one sample would under-report the cache by the cohort
    /// size, which is what this accessor exists to measure (spec §8).
    #[test]
    fn the_held_count_sums_every_samples_window() {
        let (first, _first_drawn) = source_over(&[region(45, 45)]);
        let (second, _second_drawn) = source_over(&[region(46, 46)]);
        let mut cache = ObservationCache::over(vec![first, second]);

        assert_eq!(
            cache.held_observations_len(),
            0,
            "nothing has been drawn yet"
        );
        cache
            .cover(region(40, 50))
            .expect("the fixture sources hold");

        assert_eq!(
            cache.held_observations_len(),
            2,
            "one from each sample, not one in all",
        );
    }

    /// And a cache over no samples holds nothing rather than failing.
    #[test]
    fn the_held_count_is_zero_over_no_samples() {
        let cache: ObservationCache<CountingSource> = ObservationCache::over(Vec::new());

        assert_eq!(cache.held_observations_len(), 0);
    }

    /// **The building regions tile the analysed ground exactly** — adjacent, in genome order,
    /// none of them empty, and the last one clamped to the analysed region's own last base.
    ///
    /// **This is what a merge's output cannot check.** Dividing the ground changes nothing a
    /// caller sees, so a driver that ignored the width and built each analysed region as one
    /// would pass every byte-identity test in `serial.rs` while losing the whole point of the
    /// cache. Fifty bases divided five ways, including two widths that do not divide it.
    #[test]
    fn the_building_regions_tile_the_analysed_ground_exactly() {
        for bases in [1, 7, 20, 50, 600] {
            let divided: Vec<_> = building_regions_of(region(1, 50), width(bases)).collect();

            assert_eq!(divided[0].start, Position(1), "at {bases} bases");
            assert_eq!(
                divided.last().expect("at least one region").end,
                Position(50),
                "the last region is clamped to the analysed ground, at {bases} bases",
            );
            for pair in divided.windows(2) {
                assert_eq!(
                    pair[1].start.0,
                    pair[0].end.0 + 1,
                    "the regions are adjacent, at {bases} bases",
                );
            }
            for one in &divided {
                assert!(one.start <= one.end, "no empty region, at {bases} bases");
                assert!(
                    one.end.0 - one.start.0 + 1 <= u64::from(bases),
                    "no region wider than the width asked for, at {bases} bases",
                );
            }
            assert_eq!(
                divided.len(),
                usize::try_from(50u64.div_ceil(u64::from(bases))).expect("a small count"),
                "as many regions as the width divides fifty bases into",
            );
        }
    }

    /// A building region ending on the last base of the coordinate space has no successor, and
    /// the division stops rather than wrapping to zero and dividing the genome for ever.
    ///
    /// **Bounded, because the failure this guards against is unbounded.** With the wrap in
    /// place the iterator never ends, and an unbounded `collect` would take the whole test
    /// binary down with it rather than print a difference.
    #[test]
    fn the_division_stops_at_the_coordinate_ceiling() {
        let at_the_ceiling = GenomeRegion {
            contig: ContigId(0),
            start: Position(u64::MAX - 1),
            end: Position(u64::MAX),
        };

        let divided: Vec<_> = building_regions_of(at_the_ceiling, width(20))
            .take(3)
            .collect();

        assert_eq!(
            divided,
            vec![at_the_ceiling],
            "one region, and no second lap"
        );
    }

    /// **An analysed region whose ends are the wrong way round is read as the ground it
    /// names** — belt-and-braces behind the drivers' own refusal
    /// (`super::serial::refuse_malformed_analysed_regions`), and the reason the division cannot
    /// yield one empty inverted region instead.
    #[test]
    fn the_division_reads_an_inverted_region_as_the_ground_it_names() {
        let inverted = GenomeRegion {
            contig: ContigId(0),
            start: Position(50),
            end: Position(1),
        };

        assert_eq!(
            building_regions_of(inverted, width(20)).collect::<Vec<_>>(),
            vec![region(1, 20), region(21, 40), region(41, 50)],
        );
    }

    /// **The window is sufficient**, which is the property this whole file exists for and the
    /// one no single fixture can state: a builder fed from the cache closes exactly the loci it
    /// would close given every observation in the genome
    /// ([`build_region`](super::super::build::build_region)'s input contract).
    ///
    /// Randomised over 200 layouts from a seeded generator, so a disagreement is reproducible
    /// from the seed it prints. Records are disjoint and ascending within a sample, which is
    /// what the generator mints and what `build_region` asserts; one observation in ten is a
    /// deletion up to 150 bases wide, so chains across samples are common; the building regions
    /// have gaps, so the cache is asked to jump forward over ground nobody builds.
    ///
    /// Written by the D1 review, which measured it: it agrees on 600 of 600 layouts as the code
    /// stands, and disagrees on 410 of 600 under a cover capped at two sweeps.
    ///
    /// **The evict/cover/hand-out loop is written out here rather than calling the serial
    /// driver**, deliberately: this test leaves *gaps* between its builder regions, and the
    /// driver tiles the analysed ground without gaps. Unifying the two would quietly drop the
    /// gaps, which are what make the cache jump forward over ground nobody builds.
    #[test]
    fn a_builder_fed_from_the_cache_closes_the_loci_a_whole_stretch_would() {
        use super::super::build::build_region;
        use super::super::{MaxCohortLocusSpan, MinAltObs};

        /// A seeded linear congruential generator — no dependency, and the seed is in the
        /// failure message.
        struct Seeded(u64);
        impl Seeded {
            fn next(&mut self, bound: u64) -> u64 {
                self.0 = self
                    .0
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (self.0 >> 33) % bound
            }
        }

        let max_span = MaxCohortLocusSpan::DEFAULT;
        let keep = MinAltObs::DEFAULT;
        let ground_end = 400u64;
        let mut disagreements = Vec::new();

        for seed in 0..200u64 {
            let mut draw = Seeded(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xDEAD_BEEF);
            let samples = 1 + draw.next(8) as usize;
            let contigs = 1 + draw.next(2) as u32;

            let mut layouts: Vec<Vec<GenomeRegion>> = Vec::new();
            for _ in 0..samples {
                let mut regions = Vec::new();
                for contig in 0..contigs {
                    let mut at_base = 1 + draw.next(10);
                    while at_base <= ground_end {
                        let width = match draw.next(10) {
                            0 => 1 + draw.next(150),
                            1 | 2 => 1 + draw.next(20),
                            _ => 1 + draw.next(4),
                        };
                        let end = at_base + width - 1;
                        regions.push(region_on(contig, at_base, end));
                        at_base = end + 1 + draw.next(6);
                    }
                }
                layouts.push(regions);
            }

            let width = 1 + draw.next(12);
            let mut builder_regions = Vec::new();
            for contig in 0..contigs {
                let mut at_base = 1u64;
                while at_base <= ground_end {
                    if draw.next(4) != 0 {
                        builder_regions.push(region_on(contig, at_base, at_base + width - 1));
                    }
                    at_base += width;
                }
            }

            // The oracle: every builder is handed the whole stretch.
            let whole: Vec<Vec<SampleLocusObservations>> = layouts
                .iter()
                .map(|regions| regions.iter().map(|at| observation_over(*at)).collect())
                .collect();
            let whole_slices: Vec<&[SampleLocusObservations]> =
                whole.iter().map(Vec::as_slice).collect();
            let oracle: Vec<String> = builder_regions
                .iter()
                .map(|at| format!("{:?}", build_region(*at, &whole_slices, max_span, keep)))
                .collect();

            // Under test: every builder is handed the cache's window.
            let sources: Vec<CountingSource> = layouts
                .iter()
                .map(|regions| source_over(regions).0)
                .collect();
            let mut cache = ObservationCache::over(sources);
            let mut through_cache = Vec::new();
            for at in &builder_regions {
                cache.evict_before(GenomePosition {
                    contig: at.contig,
                    position: at.start,
                });
                cache.cover(*at).expect("the fixture sources hold");
                through_cache.push(cache.with_observations(*at, |windows| {
                    format!("{:?}", build_region(*at, windows, max_span, keep))
                }));
            }

            if oracle != through_cache {
                let first = oracle
                    .iter()
                    .zip(&through_cache)
                    .position(|(whole, windowed)| whole != windowed)
                    .expect("they differ somewhere");
                disagreements.push(format!(
                    "seed {seed}, region {}\n  whole stretch: {}\n  through cache: {}",
                    builder_regions[first], oracle[first], through_cache[first],
                ));
            }
        }

        assert!(
            disagreements.is_empty(),
            "{} of 200 layouts disagreed with the whole-stretch builder; the first:\n{}",
            disagreements.len(),
            disagreements.first().expect("not empty"),
        );
    }

    /// A cohort locus over `region`. The organiser reads nothing but where a locus sits, so
    /// the allele table and the per-sample support are left empty: what these tests check is
    /// which loci come out and in what order, and a fabricated one carries that in its region.
    fn locus_at(region: GenomeRegion) -> CohortObservation {
        CohortObservation {
            region,
            alleles: Vec::new(),
            per_sample: Vec::new(),
        }
    }

    /// One region's outcome: loci over `locus_regions`, and failures over `failed_regions`.
    fn outcome_of(
        locus_regions: &[GenomeRegion],
        failed_regions: &[GenomeRegion],
    ) -> RegionOutcome {
        RegionOutcome {
            cohort_observations: locus_regions.iter().copied().map(locus_at).collect(),
            failed_locus_spans: failed_regions.to_vec(),
        }
    }

    /// The regions of everything the organiser will now part with, in the order it parts
    /// with them.
    fn drained_regions(organiser: &mut Organiser) -> Vec<GenomeRegion> {
        organiser
            .drain_ready()
            .map(|observation| observation.region)
            .collect()
    }

    /// The plain case: regions arriving in their own order release as they arrive.
    #[test]
    fn a_region_arriving_in_its_turn_releases_at_once() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(1, 3)]);

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// The reorder buffer's whole job: a region cannot pass the one in front of it.
    #[test]
    fn a_region_that_arrives_early_waits_for_the_one_before_it() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        assert!(
            drained_regions(&mut organiser).is_empty(),
            "region 1 was released before region 0 had arrived, so its loci were never \
             offered the chance to be displaced by an earlier owner's",
        );

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(21, 24)],
        );
    }

    /// A whole convoy held behind one missing index comes out at once, in index order, when
    /// that index lands — and the order is the regions' own, not the order they arrived in.
    #[test]
    fn a_convoy_held_behind_a_gap_releases_in_region_order_when_the_gap_closes() {
        let mut organiser = Organiser::new();

        for (index, locus) in [
            (3, region(61, 63)),
            (1, region(21, 24)),
            (2, region(41, 41)),
        ] {
            organiser.submit(RegionIndex(index), outcome_of(&[locus], &[]));
        }
        assert!(drained_regions(&mut organiser).is_empty());

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(21, 24), region(41, 41), region(61, 63)],
        );
    }

    /// **Exactly one outcome per region, empty ones included** (spec §6.3). The empty region
    /// carries no loci and no failures, and delivering it is the only thing that lets the
    /// region behind it out.
    #[test]
    fn an_empty_region_still_lets_the_region_behind_it_out() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(RegionIndex(0), outcome_of(&[], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// The other half of the rule above: a region that never submits holds back every region
    /// after it, for ever. This is what a builder dropping its result costs.
    ///
    /// **The three held regions carry one, two and three loci**, so the refusal's count of
    /// *regions* (3) differs from the loci they hold (6): a count that added up the held loci
    /// instead would read 6 and be caught here. The reviewer found that mutation alive against
    /// an earlier fixture of three one-locus regions, where the two rules give the same number.
    #[test]
    fn a_region_that_never_submits_holds_back_every_region_after_it() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(
            RegionIndex(2),
            outcome_of(&[region(31, 34), region(36, 36)], &[]),
        );
        organiser.submit(
            RegionIndex(3),
            outcome_of(&[region(41, 44), region(46, 46), region(48, 48)], &[]),
        );

        assert!(drained_regions(&mut organiser).is_empty());
        assert!(!organiser.is_finished(4));
        assert_eq!(
            organiser.finish(4),
            Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: RegionIndex(0),
                regions: 4,
            }),
        );
    }

    /// **A gap at the tail of the run is a gap too**, and it is the one the organiser cannot
    /// see from its own state: nothing is held, because the regions that went missing have no
    /// later index behind them. Five regions dealt out, two heard from.
    #[test]
    fn finish_refuses_a_run_whose_last_regions_never_submitted() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(1, 60)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        drained_regions(&mut organiser);

        assert!(
            !organiser.is_finished(5),
            "three of the five regions handed out never delivered an outcome",
        );
        assert_eq!(
            organiser.finish(5),
            Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: RegionIndex(2),
                regions: 3,
            }),
        );
    }

    /// Within one region the builder's order is the genome's, and the organiser keeps it.
    #[test]
    fn the_loci_of_one_region_come_out_in_the_order_the_builder_gave_them() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7), region(11, 18)], &[]),
        );

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(7, 7), region(11, 18)],
        );
    }

    /// The count spec §3.3 needs is the sum over the regions, and it reaches it whatever
    /// order the regions arrive in.
    #[test]
    fn the_failed_locus_count_sums_every_released_region() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[], &[region(21, 90)]));
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[], &[region(1, 60), region(3, 70)]),
        );
        organiser.submit(RegionIndex(2), outcome_of(&[], &[]));

        assert_eq!(organiser.failed_locus_count(), 3);
    }

    /// **The count is of what has been released, not of what has been submitted**, so that
    /// the next step can drop a failed locus an earlier one displaced without the total
    /// having counted it already. Region 1 is in the organiser's hands and uncounted until
    /// region 0 lets it out.
    #[test]
    fn the_failed_locus_count_ignores_a_region_still_held_behind_a_gap() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(1),
            outcome_of(&[], &[region(21, 90), region(25, 95)]),
        );
        assert_eq!(organiser.failed_locus_count(), 0);

        organiser.submit(RegionIndex(0), outcome_of(&[], &[region(1, 60)]));
        assert_eq!(organiser.failed_locus_count(), 3);
    }

    /// **Releasing happens when the region arrives, not when the caller asks.** Nothing has
    /// been drained here, and the failed loci of both regions are already counted — which is
    /// what says the release ran inside `submit`.
    #[test]
    fn a_region_is_released_on_arrival_rather_than_at_the_next_drain() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(1, 60)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[], &[region(21, 90)]));

        assert_eq!(organiser.failed_locus_count(), 2);
    }

    /// A drain the caller stops halfway leaves the rest where it was: the loci are taken one
    /// at a time from the front, not emptied into the iterator and dropped with it.
    #[test]
    fn a_drain_stopped_halfway_leaves_the_rest_for_the_next_call() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7), region(11, 18)], &[]),
        );

        let first = organiser
            .drain_ready()
            .next()
            .expect("the region's first locus");
        assert_eq!(first.region, region(1, 3));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(7, 7), region(11, 18)],
        );
    }

    /// Nothing outstanding, which is what a run asserts at its end — and `finish` returns `Ok`
    /// on exactly those runs.
    #[test]
    fn an_organiser_is_finished_once_every_region_is_released_and_drained() {
        let mut organiser = Organiser::new();
        assert!(Organiser::new().is_finished(0));

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert!(
            !organiser.is_finished(1),
            "a locus released and not yet taken is a locus the run still owes its output",
        );

        drained_regions(&mut organiser);
        assert!(organiser.is_finished(1));
        assert_eq!(organiser.finish(1), Ok(0));
    }

    /// Loci released and never taken truncate the output exactly as a missing region does,
    /// so `finish` refuses them too — and says so without inventing a gap.
    #[test]
    fn finish_refuses_a_run_with_loci_released_and_never_drained() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7)], &[]),
        );

        assert_eq!(
            organiser.finish(1),
            Err(RunEndedShort::LociNeverDrained { loci: 2 }),
        );
    }

    /// **Both ways of ending short at once**, each named by its own count: region 1 stalled
    /// the drain and region 2 is held behind it — two regions handed out and never released,
    /// though they hold three loci between them — while region 0's two loci were released and
    /// never taken. A refusal that could carry one count at a time would drop half of this.
    #[test]
    fn finish_names_both_counts_when_a_stall_and_undrained_loci_coincide() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7)], &[]),
        );
        organiser.submit(
            RegionIndex(2),
            outcome_of(&[region(41, 41), region(45, 45), region(51, 51)], &[]),
        );

        assert_eq!(
            organiser.finish(3),
            Err(RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
                first_stalled: RegionIndex(1),
                regions: 2,
                loci: 2,
            }),
        );
    }

    /// The refusal reads back the way a person will read it, each count against its own noun.
    /// The two counts are deliberately different, so swapping them in the message cannot pass.
    #[test]
    fn the_refusal_names_each_count_against_its_own_noun() {
        let both = RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
            first_stalled: RegionIndex(4),
            regions: 3,
            loci: 7,
        };
        assert_eq!(
            both.to_string(),
            "region 4 never delivered its result, so 3 region(s) were handed out and never \
             released, and 7 released locus/loci were never taken",
        );

        assert_eq!(
            RunEndedShort::LociNeverDrained { loci: 7 }.to_string(),
            "7 released locus/loci were never taken from the organiser",
            "a run that lost nothing to a stall must not have one described to it",
        );
    }

    /// The end of a healthy run: the failed-locus total comes back, and nothing is owed.
    #[test]
    fn finish_returns_the_failed_locus_total_of_a_run_that_ended_clean() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(1, 60)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[], &[]));
        drained_regions(&mut organiser);

        assert_eq!(organiser.finish(2), Ok(1));
    }

    /// A run over no regions at all is finished and owes nothing.
    #[test]
    fn an_organiser_that_was_given_nothing_finishes_clean() {
        assert_eq!(Organiser::new().finish(0), Ok(0));
    }

    /// One outcome per region, whether or not its loci have gone out yet.
    #[test]
    #[should_panic(expected = "delivered a second outcome")]
    fn a_region_submitted_twice_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
    }

    /// The duplicate is refused **and the first outcome is the one kept** — the check runs
    /// before the map is written, so a second delivery cannot displace what the region
    /// already owns on its way to the panic.
    #[test]
    fn a_refused_second_outcome_does_not_displace_the_first() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));

        let refused = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            organiser.submit(RegionIndex(1), outcome_of(&[region(91, 94)], &[]));
        }));
        assert!(
            refused.is_err(),
            "a second outcome for region 1 must be refused"
        );

        organiser.submit(RegionIndex(0), outcome_of(&[], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// An index that has already been released cannot be submitted again either — the drain
    /// has passed it, so its loci would come out behind loci that are already downstream.
    #[test]
    #[should_panic(expected = "arrived after it was released")]
    fn a_region_submitted_after_it_was_released_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
    }

    /// An outcome for a region the run says it never handed out is the same class of hand-out
    /// bug as a duplicate, and is refused where it is first visible.
    #[test]
    #[should_panic(expected = "though the run handed out only")]
    fn an_outcome_for_a_region_never_handed_out_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(7), outcome_of(&[region(1, 3)], &[]));
        let _ = organiser.finish(3);
    }
}
