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
//! **The organiser is [`super::organise`], not here** — the two shared a file until the owner
//! split them at Checkpoint E. The argument for keeping them together had been that the
//! organiser would become the cache's only writer, at which point
//! [`cover`](ObservationCache::cover) and [`evict_before`](ObservationCache::evict_before)
//! could turn private to that file. E1's review found the premise false: the cached serial
//! driver calls both from a sibling module and nothing in the plan removes it. So the
//! reachable narrowing is `pub(super)` — taken here — which a file of its own gets equally,
//! and against it stood one file of two thousand-odd lines covering two types that share no
//! field and no function.

use super::CohortLocusBuilderRegionsLen;
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion, Position};

/// One sample's observations in coordinate order, **and the place a record goes when the
/// merge has finished with it**.
///
/// **The second half is why this is a trait and not an `Iterator` bound.** The merge frees
/// far more than it allocates — 6.4 million blocks a round against 216 thousand on the
/// tomato panel, counted — because the records it walks were allocated by the stage upstream
/// and are released as it passes them. Measured by making the merge leak instead of free,
/// that is **63% of the eight-thread merge** and 23% of the single-threaded one. No
/// scheduling change reaches it: the work is real, it is just work nobody needs done.
///
/// So a source is offered its own spent records back. A source that mints them — the
/// generator, or a psp reader decoding into buffers — can fill the record it is handed
/// instead of allocating a new one, and then neither side allocates per position after the
/// first window. A source that cannot simply ignores the offer, which is what the blanket
/// implementation below does for every plain iterator, so nothing that exists today has to
/// change.
///
/// **The spare is an offer and not an obligation**, and that is the whole of the contract: a
/// source may fill it, drop it, or keep it for later, and the record it returns need have
/// nothing to do with it. Making it an obligation would mean a source could not decide per
/// record whether reuse is possible — which a decoder must, since a record whose buffers are
/// the wrong size is cheaper to allocate than to reshape.
pub trait ObservationSource {
    /// What a failed read is. The cache adds nothing to it and passes it through, so it must
    /// name the sample it came from (arch §5).
    type Error;

    /// The next observation, or `None` once this sample is spent.
    ///
    /// `spare` is a record the merge will not read again, offered for reuse. It is `None`
    /// when the cache has none to hand back.
    ///
    /// **Once this answers `None` it is never called again**, which the cache guarantees with
    /// a flag of its own — the guard [`Fuse`](std::iter::Fuse) used to give, and the reason it
    /// matters is unchanged: a source that yielded `Some` after a `None` would be drawn in
    /// behind the window's own right edge and so silently out of coordinate order. A
    /// *failure* is `Some(Err(_))` and leaves the source live, which is what lets a cover be
    /// made again.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, Self::Error>>;
}

/// Every iterator of one sample's observations is a source that does not reuse.
///
/// This is what keeps the trait from being a migration: the fixtures, the probes and the
/// direct path all hand the cache a plain iterator and go on working, paying exactly the
/// allocator traffic they paid before. Reuse is something a source opts into by implementing
/// the trait itself.
impl<I, E> ObservationSource for I
where
    I: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    type Error = E;

    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, E>> {
        drop(spare);
        self.next()
    }
}

/// Every sample's observations over the ground currently assigned to builders (spec §6.4).
///
/// **One forward reader per sample, and a window of what it has drawn.** The cache never
/// seeks and never goes back: the only ways the window moves are [`cover`](Self::cover),
/// which draws it forward, and [`evict_before`](Self::evict_before), which drops from its
/// left edge.
///
/// **This is the module's dominant memory** (spec §8), and it is what a building region's
/// width is paid for in: the ground it spans is `builders × cohort_locus_builder_regions_len`
/// — 3,200 bases at 16 builders on the default 200-base regions — plus the tail of
/// observations reaching past it.
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
/// **The window and the reader are separate fields, and that is what lets a cover run beside
/// the builders.** A builder reads [`released`](Self::released); a cover advances
/// [`readers`](Self::readers). Held in one `Vec<SampleWindow>` as they used to be, the two
/// could not be borrowed apart, so the organiser had to finish covering before any builder
/// started — the round barrier that
/// [`merge_cohort_in_parallel`](super::parallel::merge_cohort_in_parallel) pays. Split, they
/// can be handed to the two halves of a `rayon::join`.
pub struct ObservationCache<S> {
    /// What builders may read, one window per sample in the run's sample order — the order
    /// every consumer indexes by.
    ///
    /// **Immutable while a round is being built.** A cover running beside that round appends
    /// to `readers[sample].drawn` instead, and [`promote`](Self::promote) moves those in when
    /// the round is over.
    released: Vec<Vec<SampleLocusObservations>>,
    /// Drawn by a cover that ran beside a round, not yet released, one per sample.
    ///
    /// **Empty everywhere except between a [`cover_beside`](Self::cover_beside) and its
    /// [`promote`](Self::promote).** The two covers that do not overlap a round draw straight
    /// into `released`, so they never touch this and pay nothing for its existence.
    ///
    /// It is a field of the cache rather than of a reader so that a cover can hold it mutably
    /// while the builders hold `released` — the whole point of the split.
    drawn: Vec<Vec<SampleLocusObservations>>,
    /// One forward reader per sample, in the same order.
    readers: Vec<SampleReader<S>>,
    /// How far a **successful** [`cover`](Self::cover) has drawn, genome-wide.
    ///
    /// A cover that failed does not move it, which is what lets
    /// [`with_observations`](Self::with_observations) refuse ground no reader reached rather
    /// than hand out a window that is short — and short is a locus closed over the wrong
    /// ground, which is a wrong answer rather than a failure.
    covered_to: Option<GenomePosition>,
    /// How far the ground in `released` reaches — what a builder's window is checked against.
    ///
    /// Equal to `covered_to` except while a cover is running beside a round of builders, when
    /// it is one round behind: what those builders may read is what the *previous* cover drew,
    /// and reading against the newer mark would let one of them ask for ground whose records
    /// are still in `drawn`.
    released_to: Option<GenomePosition>,
}

/// One sample's forward reader, and what it has drawn that its window has not yet been given.
struct SampleReader<S> {
    /// The forward reader. Never seeks, never rewinds.
    source: S,
    /// Whether the source has already answered `None`.
    ///
    /// **This is the guard [`Fuse`](std::iter::Fuse) used to be**, kept because the reason
    /// for it is unchanged and moved here because a source is now a trait rather than an
    /// iterator: one that yielded `Some` after a `None` would be drawn in behind the window's
    /// own right edge, and so silently out of coordinate order. Set only on `None`, so a
    /// *failure* — `Some(Err(_))` — leaves the source live and a cover can be made again.
    spent: bool,
    /// Records this sample's window has finished with, offered back to its source.
    ///
    /// **Capped at what the sample currently holds**, which is the round's own ground and is
    /// the term spec §8 already prices — so recycling cannot make the cache's memory a
    /// different shape from the one that section bounds. What does not fit is freed as
    /// before.
    spare: Vec<SampleLocusObservations>,
    /// Where the last observation drawn from `source` began — the ordering check's memory.
    last_drawn: Option<GenomePosition>,
}

impl<S> ObservationCache<S> {
    /// A cache over one source per sample, in the run's sample order.
    ///
    /// Zero samples is not an error here: refusing a zero-sample *run* happens where the run
    /// is configured (spec §7.2), and a cache over an empty cohort covers nothing and hands
    /// nothing out.
    pub fn over(sources: Vec<S>) -> Self {
        Self {
            released: sources.iter().map(|_| Vec::new()).collect(),
            drawn: sources.iter().map(|_| Vec::new()).collect(),
            readers: sources
                .into_iter()
                .map(|source| SampleReader {
                    source,
                    spent: false,
                    spare: Vec::new(),
                    last_drawn: None,
                })
                .collect(),
            covered_to: None,
            released_to: None,
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
        self.released_windows().with_observations(span, f)
    }

    /// The windows builders may read, borrowed apart from the readers.
    ///
    /// **This is the whole of what a builder is allowed to see**, and giving it a name is what
    /// lets [`cover_beside`](Self::cover_beside) hand it to one half of a `rayon::join` while
    /// the other half advances the readers.
    pub(super) fn released_windows(&self) -> ReleasedWindows<'_> {
        ReleasedWindows {
            per_sample: &self.released,
            covered_to: self.released_to,
        }
    }

    /// How many observations are held, summed across samples — including any a cover has drawn
    /// but not yet released.
    pub fn held_observations_len(&self) -> usize {
        self.released.iter().map(Vec::len).sum::<usize>()
            + self.drawn.iter().map(Vec::len).sum::<usize>()
    }
}

/// Every sample's released observations — what a round's builders read, and nothing else.
///
/// **Borrowed for the length of a round and immutable while it runs**, which is what makes it
/// safe to build one round while the next round's ground is being drawn: the drawing appends to
/// the readers' own buffers and cannot move a record a builder is looking at.
pub(super) struct ReleasedWindows<'a> {
    /// One window per sample, in the run's sample order.
    per_sample: &'a [Vec<SampleLocusObservations>],
    /// How far these windows were covered — the mark a builder's span is checked against.
    covered_to: Option<GenomePosition>,
}

impl ReleasedWindows<'_> {
    /// Every sample's held observations from the first one that reaches into `span`, for the
    /// length of the call — a builder's only way in, and read-only by construction.
    ///
    /// The rules are [`ObservationCache::with_observations`]'s, which is this method's only
    /// other caller and where they are stated.
    pub(super) fn with_observations<R>(
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
        let windowing = super::timing::Stopwatch::start();
        let observations_per_sample: Vec<&[SampleLocusObservations]> = self
            .per_sample
            .iter()
            .map(|held| &held[first_reaching_index(held, left_edge)..])
            .collect();
        windowing.add_to(&super::timing::WINDOW_NANOS);
        f(&observations_per_sample)
    }
}

impl<S> ObservationCache<S> {
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
    pub(super) fn evict_before(&mut self, position: GenomePosition) {
        for (held, reader) in self.released.iter_mut().zip(&mut self.readers) {
            let first_survivor = first_reaching_index(held, position);
            let room = held.len();
            for record in held.drain(..first_survivor) {
                if reader.spare.len() < room {
                    reader.spare.push(record);
                }
            }
        }
    }

    /// [`evict_before`](Self::evict_before), moving what it drops into `graveyard` instead of
    /// freeing it here.
    ///
    /// **Same records dropped, freed somewhere else.** What eviction costs is not deciding what
    /// to drop — that is one walk of the window's prefix — but returning every record's
    /// buffers to the allocator, and that runs on whichever thread calls it. Moving the records
    /// out first makes the free a job of its own that a caller can put beside work rather than
    /// in front of it; the caller owns `graveyard` and decides when it dies.
    ///
    /// **The records moved out are unreachable by construction**, which is what makes this safe
    /// to hand to another thread: they are exactly the ones
    /// [`evict_before`](Self::evict_before) would have dropped, so nothing the cache still
    /// holds and no window a builder can be given refers to them.
    ///
    /// A caller that never empties `graveyard` holds every record the run ever evicted, which
    /// is the whole cohort — so it is a buffer to drain each round, not one to accumulate.
    pub(super) fn evict_before_into(
        &mut self,
        position: GenomePosition,
        graveyard: &mut Vec<SampleLocusObservations>,
    ) {
        for (held, reader) in self.released.iter_mut().zip(&mut self.readers) {
            let first_survivor = first_reaching_index(held, position);
            let room = held.len();
            for record in held.drain(..first_survivor) {
                // **The source's offer comes first, the graveyard second.** A record the
                // source will refill costs nothing to keep and saves the allocation the next
                // draw would make; only what the sample cannot hold goes on to be freed, and
                // that freeing is what the caller's partner thread is for.
                if reader.spare.len() < room {
                    reader.spare.push(record);
                } else {
                    graveyard.push(record);
                }
            }
        }
    }

    /// [`evict_before_into`](Self::evict_before_into) with the samples evicted at the same time
    /// as each other.
    ///
    /// **Eviction is per sample and shares nothing between them**, so this is the same walk on a
    /// pool: each sample decides its own first survivor, drains its own prefix into its own
    /// spare list, and pushes what will not fit into its own graveyard. The graveyards come back
    /// as one so the caller still has a single buffer to free beside the next round's builders.
    ///
    /// **Why it is worth a method rather than left serial**: it runs between rounds, when every
    /// worker but one is idle, and on the tomato benchmark's 63 accessions it was 13.6% of the
    /// merge — the largest part that no thread but the organiser's ever ran.
    ///
    /// The evicted records are freed by whichever worker drained them, which is not the worker
    /// that drew them; what that costs is the same cross-thread free the merge already pays
    /// everywhere else.
    pub(super) fn evict_before_in_parallel(
        &mut self,
        position: GenomePosition,
        graveyard: &mut Vec<SampleLocusObservations>,
    ) where
        S: Send,
    {
        use rayon::prelude::*;

        let dead: Vec<Vec<SampleLocusObservations>> = self
            .released
            .par_iter_mut()
            .zip(self.readers.par_iter_mut())
            .map(|(held, reader)| {
                let first_survivor = first_reaching_index(held, position);
                let room = held.len();
                let mut dropped = Vec::new();
                for record in held.drain(..first_survivor) {
                    if reader.spare.len() < room {
                        reader.spare.push(record);
                    } else {
                        dropped.push(record);
                    }
                }
                dropped
            })
            .collect();
        for mut sample in dead {
            graveyard.append(&mut sample);
        }
    }

    /// Give the builders everything the last cover drew — the round boundary.
    ///
    /// **Called only by the driver that covers beside a round** ([`cover_beside`](Self::cover_beside));
    /// the two covers that do not overlap a round draw straight into the released windows and
    /// have nothing to promote.
    ///
    /// **What it costs is a move of each drawn record**, not of what they point at: the
    /// records' own buffers stay where they are and only the `SampleLocusObservations` structs
    /// are copied to the end of the window. That is the price of overlapping, and it is paid
    /// once per record per merge.
    pub(super) fn promote(&mut self) {
        for (held, drawn) in self.released.iter_mut().zip(&mut self.drawn) {
            held.append(drawn);
        }
        self.released_to = self.covered_to;
    }
}

impl<S, E> ObservationCache<S>
where
    S: ObservationSource<Error = E>,
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
    pub(super) fn cover(&mut self, region: GenomeRegion) -> Result<(), E> {
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
        self.released_to = self.covered_to;
        Ok(())
    }

    /// [`cover`](Self::cover), with each sweep's samples swept concurrently.
    ///
    /// **Same fixpoint, a different schedule.** The serial sweep is Gauss-Seidel — sample `j`
    /// sees the reach sample `i < j` widened inside the same sweep — and this is Jacobi: every
    /// sample is drawn against the reach the last sweep ended on, and the sweep's answer is the
    /// widest of theirs. Drawing is monotone in the reach it is given and the reach only ever
    /// grows, so both schedules climb to the same least fixpoint above the region's last base,
    /// and the held window is the same because every sample's last draw is against that
    /// fixpoint. What differs is the sweep count: a chain that runs through the cohort in
    /// decreasing sample order costs one sweep per link either way, but a chain the serial form
    /// follows within one sweep costs this one a sweep per link.
    pub(super) fn cover_in_parallel(&mut self, region: GenomeRegion) -> Result<(), E>
    where
        S: Send,
        E: Send,
    {
        // **Straight into the released windows**, because no builder is running: this cover
        // scans and appends to the same buffer, exactly as the serial sweep does, and so has
        // nothing to promote afterwards.
        let ObservationCache {
            released, readers, ..
        } = self;
        let chain_reach = draw_every_sample(readers, None, released, region)?;

        self.covered_to = Some(
            self.covered_to
                .map_or(chain_reach, |reached| reached.max(chain_reach)),
        );
        self.released_to = self.covered_to;
        Ok(())
    }

    /// Build a round while the next round's ground is drawn — the two halves of the merge that
    /// used to have to take turns.
    ///
    /// **The builders read what the last cover released and the cover writes somewhere else**,
    /// so neither can move a record the other is looking at: `build` is handed
    /// [`ReleasedWindows`], and the drawing appends to each reader's own buffer. Nothing the
    /// cover draws becomes visible until [`promote`](Self::promote), which the caller runs when
    /// the round is over.
    ///
    /// **Both halves run whatever the other does**, including when the cover fails: the round's
    /// builders have everything they need already, so a source that dies part-way leaves a
    /// finished round and an error, and the caller decides. That is the same contract
    /// [`cover`](Self::cover) has — the window keeps what was drawn before the failure — with
    /// the round's output added.
    ///
    /// **What it costs is a second round of observations resident**: the released windows hold
    /// the round being built and the readers hold the round being drawn (spec §8's term,
    /// doubled).
    pub(super) fn cover_beside<R>(
        &mut self,
        region: GenomeRegion,
        build: impl FnOnce(&ReleasedWindows<'_>) -> R + Send,
    ) -> (R, Result<(), E>)
    where
        S: Send,
        E: Send,
        R: Send,
    {
        let ObservationCache {
            released,
            drawn,
            readers,
            covered_to,
            released_to,
        } = self;
        let released = &*released;
        let windows = ReleasedWindows {
            per_sample: released,
            covered_to: *released_to,
        };
        let (built, covered) = rayon::join(
            move || build(&windows),
            || draw_every_sample(readers, Some(released), drawn, region),
        );
        match covered {
            Ok(chain_reach) => {
                *covered_to =
                    Some(covered_to.map_or(chain_reach, |reached| reached.max(chain_reach)));
                (built, Ok(()))
            }
            Err(failed) => (built, Err(failed)),
        }
    }

    /// One sweep of every sample against `chain_reach`. Answers whether any of them moved it.
    fn sweep(&mut self, chain_reach: &mut GenomePosition) -> Result<bool, E> {
        let mut reach_grew = false;
        super::timing::COVER_SWEEPS.add(1);
        for (held, reader) in self.released.iter_mut().zip(&mut self.readers) {
            // Not `reach_grew |= reader.draw_to(…)?`: that reads as though the call could be
            // skipped, and a later `||` in its place would skip it.
            //
            // **The serial cover draws straight into the released window**, which it can
            // because no builder is running: it passes that window as both the ground to scan
            // and the buffer to append to, and so has nothing to promote afterwards.
            let drawing = super::timing::Stopwatch::start();
            let grew = reader.draw_to(&[], held, chain_reach)?;
            drawing.add_to(&super::timing::COVER_BUSY_NANOS);
            if grew {
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }
}

/// How many of a sample's held records begin at or before `reach` — the prefix a cover admits.
///
/// **A bisection over two slices that are one window.** `front` is what builders may read and
/// `back` what a cover has drawn since; together they are the sample's records in coordinate
/// order, so their starts ascend across the join and the answer is one `partition_point` on
/// whichever half holds the boundary.
fn starts_at_or_before(
    front: &[SampleLocusObservations],
    back: &[SampleLocusObservations],
    reach: GenomePosition,
) -> usize {
    let in_front = front.partition_point(|observation| observation.start_position() <= reach);
    if in_front < front.len() {
        return in_front;
    }
    front.len() + back.partition_point(|observation| observation.start_position() <= reach)
}

/// The record at `at` in the window `front` then `back`.
fn observation_at<'a>(
    front: &'a [SampleLocusObservations],
    back: &'a [SampleLocusObservations],
    at: usize,
) -> &'a SampleLocusObservations {
    front.get(at).unwrap_or_else(|| &back[at - front.len()])
}

/// Draw every sample forward until `region` is covered and the chain from its last base closes,
/// with the samples drawn at the same time as each other — [`ObservationCache::cover`]'s
/// fixpoint on a pool.
///
/// **Jacobi, where the serial sweep is Gauss-Seidel**, and the argument that both reach the same
/// least fixpoint is [`ObservationCache::cover_in_parallel`]'s.
///
/// Each sample scans its released window and then what it has already drawn, and appends what it
/// draws to the second. **Whether the released window is also where the drawing lands is the
/// caller's choice**, and it is what separates a cover that overlaps a round of builders from one
/// that does not.
fn draw_every_sample<S, E>(
    readers: &mut [SampleReader<S>],
    fronts: Option<&[Vec<SampleLocusObservations>]>,
    backs: &mut [Vec<SampleLocusObservations>],
    region: GenomeRegion,
) -> Result<GenomePosition, E>
where
    S: ObservationSource<Error = E> + Send,
    E: Send,
{
    use rayon::prelude::*;

    let mut chain_reach = GenomePosition {
        contig: region.contig,
        position: region.end.max(region.start),
    };
    loop {
        let snapshot = chain_reach;
        let widest = readers
            .par_iter_mut()
            .zip(backs.par_iter_mut())
            .enumerate()
            .map(|(sample, (reader, back))| {
                let front = fronts.map_or(&[][..], |per_sample| &per_sample[sample][..]);
                let drawing = super::timing::Stopwatch::start();
                let mut reach = snapshot;
                reader.draw_to(front, back, &mut reach)?;
                drawing.add_to(&super::timing::COVER_BUSY_NANOS);
                Ok(reach)
            })
            .try_reduce(|| snapshot, |left, right| Ok(left.max(right)))?;
        super::timing::COVER_SWEEPS.add(1);
        if widest == chain_reach {
            break;
        }
        chain_reach = widest;
    }
    Ok(chain_reach)
}

impl<S, E> SampleReader<S>
where
    S: ObservationSource<Error = E>,
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
    fn draw_to(
        &mut self,
        front: &[SampleLocusObservations],
        back: &mut Vec<SampleLocusObservations>,
        chain_reach: &mut GenomePosition,
    ) -> Result<bool, E> {
        let mut reach_grew = false;
        super::timing::HELD_WHEN_DRAWING.add((front.len() + back.len()) as u64);
        super::timing::COVER_DRAWS.add(1);

        // **The records already in the window are folded in by bisection, not by walking
        // them.** A sample's records are disjoint and ascending, so their starts ascend and
        // their reaches are monotone — the same two facts `first_reaching_index` rests on. So
        // the ones a given reach admits are a prefix, and the furthest any of them reaches is
        // the reach of the last one in it: finding that prefix is a search, and the whole
        // prefix need never be touched.
        //
        // **Widening can admit more, so it repeats until the prefix stops growing**, and each
        // turn strictly lengthens it. What the walk it replaces cost was one comparison per
        // record in the window per cover, and the window is a whole round of ground.
        let mut considered = 0;
        loop {
            let admitted = starts_at_or_before(front, back, *chain_reach);
            if admitted == considered {
                break;
            }
            considered = admitted;
            let furthest = observation_at(front, back, admitted - 1).reach_position();
            if furthest <= *chain_reach {
                break;
            }
            *chain_reach = furthest;
            reach_grew = true;
        }

        // **Then the records that are not in the window yet, one at a time**, because where the
        // next one begins is not known until it is drawn. This can only start once every held
        // record has been admitted: a record still in the window that the reach does not admit
        // begins after it, and everything drawn from here on begins later still.
        while considered == front.len() + back.len() {
            let Some(observation) = self.draw_next()? else {
                break;
            };
            let starts_within = observation.start_position() <= *chain_reach;
            let reach = observation.reach_position();
            // Held either way. One beyond the reach stays in the window, and the next cover
            // reconsiders it against a later reach.
            back.push(observation);
            if !starts_within {
                break;
            }
            considered += 1;
            if reach > *chain_reach {
                *chain_reach = reach;
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }

    /// The next observation from the source, or `None` once it is spent.
    fn draw_next(&mut self) -> Result<Option<SampleLocusObservations>, E> {
        if self.spent {
            return Ok(None);
        }
        // **The offer, and it is the whole of the recycling.** A source that mints records
        // fills this one instead of allocating; a source that cannot drops it and allocates,
        // which is what every plain iterator does.
        let spare = self.spare.pop();
        let Some(next) = self.source.next_observation(spare).transpose()? else {
            self.spent = true;
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
    // **A bisection, and it is what makes the window's ordering load-bearing.** A sample's
    // records are disjoint and ascending — `build_region` refuses a sample whose are not — so
    // reach is monotone across the window and "does this one reach `position`" is false over a
    // prefix and true over the rest, which is the shape `partition_point` needs. The scan this
    // replaced would have given the same answer on a window that was not ordered; this one
    // gives a wrong answer instead of a slow one, which is why the precondition is stated here
    // rather than left to `evict_before`'s doc.
    //
    // It is worth nothing to the cached serial driver, whose window starts at the left edge
    // because it evicts immediately before every cover, and about a tenth of the parallel
    // driver's merge, where eviction opens a whole round and a region late in that round would
    // otherwise walk past every earlier region's records first.
    held_observations.partition_point(|observation| observation.reach_position() < position)
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
                // The width of an inclusive region, which is what the `+ 1` is for.
                let width = one.end.0 - one.start.0 + 1;
                assert!(
                    width <= u64::from(bases),
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
        use super::super::{MaxCohortLocusSpan, MinAltReads};

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
        let keep = MinAltReads::DEFAULT;
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
}
