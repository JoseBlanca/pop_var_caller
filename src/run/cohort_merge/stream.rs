//! **The round driver without its barrier**: the next rounds' ground is read while earlier rounds'
//! regions are still being built, and a thread that finishes a region takes the next one, whichever
//! round it belongs to (`spec/cohort_merge_overlapping_rounds.md`).
//!
//! [`super::rounds::merge_and_call_in_rounds`] reads a round's ground, builds its regions on the
//! pool, and waits for the slowest before reading the next: while it reads, no region is built, and
//! while the slowest region finishes, the other threads wait. Here each round's records are sealed
//! once read ([`ObservationCache::seal_round`]), so its regions read the seal and not the cache, and
//! the cache is free to read on.
//!
//! **Three threads besides the pool, each for a reason that was measured.**
//!
//! - *A reader*, which reads one round's ground after another, on a pool of its own
//!   ([`ReadPool`]). As a task on the regions' pool, a read waited behind region builds at every
//!   one of its parallel steps and took 4 to 6 s a round where the round driver's took 0.6.
//! - *A dispatcher*, which hands regions to the pool as their round is sealed, a few more than
//!   there are threads at a time, so that work is always waiting.
//! - *The calling thread*, which hands finished regions over in genome order through `take` — the
//!   record pipeline, which runs for seconds at a time and, done on the dispatcher's thread, left
//!   the pool without regions while it ran.
//!
//! **It holds up to [`ROUNDS_HELD`] rounds' records at once**, and the caller cuts its ground into
//! rounds that much smaller (`crate::cli::calling_run::round_shape_for`), so the run holds about the
//! ground one round of the round driver held. Many small rounds rather than two large ones because
//! the next read can start only once the oldest round is built: with two rounds the pool ran out of
//! work at every read. Measured on the tomato25 cohort's stored files (2,176 samples, about three
//! reads a position, 48 threads, 24,000 bases of ground), over 1 Mb of chromosome 0:
//!
//! | rounds held × ground each | calling | peak resident |
//! |---|---|---|
//! | the round driver, 1 × 24,000 | 3 m 17 s | 25.8 GB |
//! | 2 × 12,000 | 2 m 53 s | 37.3 GB |
//! | 4 × 6,000 | 2 m 38 s | 31.3 GB |
//! | 6 × 4,000 | 2 m 24 s | 32.4 GB |
//! | 8 × 3,000 | 2 m 28 s | 31.4 GB |
//! | 12 × 2,000 | 2 m 34 s | 30.4 GB |
//!
//! The VCF is byte-identical in every row. Over 300 kb the six rounds hold what the round driver
//! held (21.3 GB against 21.2, 43 s against 57–59); over 1 Mb, about a quarter more: what each
//! sealed round carries over from the one before, and the reads ahead of their builders.
//!
//! **What comes out is the round driver's, in the same order.** Each region is built from the
//! records the round driver would have built it from, and finished regions are handed to `take` in
//! genome order through the same frontier rule ([`super::rounds::claim`]). A failure is reported at
//! the first round, in genome order, that has one — a read that failed, or a region whose build did
//! — after every region before it was handed over.

use std::collections::VecDeque;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};

use super::build::{CohortObservation, LocusReads, build_region_handing_over_windowed};
use super::observation_cache::{
    ObservationCache, ObservationSource, ReferenceUnreadable, RoundView, building_regions_of,
};
use super::organise::first_base_of;
use super::rounds::{CalledRegion, claim};
use super::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan,
    MinAltReads, refuse_malformed_analysed_regions,
};
use crate::locus_generation::SampleLocusObservations;
use crate::types::{GenomePosition, GenomeRegion, Position};

/// **How many rounds' records are held at once, at most** — rounds with regions still to build,
/// besides the one being read. The module's table is where six comes from.
pub const ROUNDS_HELD: usize = 6;

/// **How many regions wait for a thread, at most, beyond one a thread.**
const REGIONS_QUEUED_AHEAD: usize = 4;

/// What a task on the pool tells the calling thread.
enum Done<T, E> {
    /// A round's ground was read and sealed, or the read failed.
    Covered {
        round: u64,
        sealed: Result<Arc<RoundView<E>>, E>,
    },
    /// One of a round's regions was built and called, or its build failed.
    Built {
        round: u64,
        at: usize,
        outcome: Result<CalledRegion<T>, E>,
    },
}

/// One round, from being planned to being handed over.
struct Round<T, E> {
    number: u64,
    regions: Vec<GenomeRegion>,
    /// Each region's outcome once built, in region order.
    outcomes: Vec<Option<Result<CalledRegion<T>, E>>>,
    /// Regions not built yet; the round's seal is dropped when this reaches zero.
    unbuilt: usize,
    /// Regions handed to the pool so far, in region order.
    spawned: usize,
    /// The round's seal, while regions are still to be built from it.
    sealed: Option<Arc<RoundView<E>>>,
    /// Whether the round's ground has been read.
    covered: bool,
    /// Why the round's ground could not be read, if it could not.
    cover_failed: Option<E>,
}

/// Merge and call the cohort over `analysed` with no barrier between rounds — the arguments, the
/// result and the order of what reaches `take` are [`super::rounds::merge_and_call_in_rounds`]'s.
///
/// # Errors
///
/// The first failure in genome order: a read of a round's ground, or a region's build.
#[expect(
    clippy::too_many_arguments,
    reason = "the round driver's arguments, for the reason given there"
)]
pub fn merge_and_call_streaming<S, E, Sc, T, MK, CL>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    regions_in_flight: CohortLocusBuilderRegionsInFlight,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
    make_scratch: MK,
    call_locus: CL,
    take: &mut impl FnMut(T),
    refused: &mut Vec<GenomeRegion>,
) -> Result<u64, E>
where
    S: ObservationSource<Error = E> + Sync + Send,
    E: Send + 'static + From<ReferenceUnreadable>,
    T: Send,
    Sc: Send,
    MK: Fn() -> Sc + Sync + Send,
    CL: Fn(&mut Sc, CohortObservation, &LocusReads<'_, '_>) -> T + Sync + Send,
{
    refuse_malformed_analysed_regions(analysed);

    // The rounds, cut as the round driver cuts them: up to `regions_in_flight` regions, never
    // across an analysed region.
    let mut planned = analysed.iter().flat_map(|analysed_region| {
        let mut building = building_regions_of(*analysed_region, cohort_locus_builder_regions_len);
        std::iter::from_fn(move || {
            let round: Vec<GenomeRegion> =
                building.by_ref().take(regions_in_flight.get()).collect();
            (!round.is_empty()).then_some(round)
        })
    });

    let mut owned_through: Option<GenomePosition> = None;
    let mut displaced: u64 = 0;
    let mut hand_over = |outcome: CalledRegion<T>| {
        let CalledRegion { calls, failed } = outcome;
        let mut calls = calls.into_iter().peekable();
        for span in failed {
            let span_first = first_base_of(span);
            while let Some((locus, made)) =
                calls.next_if(|(locus, _)| first_base_of(*locus) <= span_first)
            {
                if claim(&mut owned_through, &mut displaced, locus) {
                    take(made);
                }
            }
            if claim(&mut owned_through, &mut displaced, span) {
                refused.push(span);
            }
        }
        for (locus, made) in calls {
            if claim(&mut owned_through, &mut displaced, locus) {
                take(made);
            }
        }
    };

    // **Regions are handed out on a thread of their own, and handed over on this one.** Handing a
    // round over runs the record pipeline behind `take`, and while one thread did both, no region
    // went out during it and the pool ran dry. A finished region comes here, in genome order,
    // through a channel a few rounds deep, which holds the dispatcher back if handing over falls
    // behind.
    let (to_hand_over, arrived) =
        mpsc::sync_channel::<CalledRegion<T>>(ROUNDS_HELD.saturating_mul(regions_in_flight.get()));
    let (told, hear) = mpsc::channel::<Done<T, E>>();
    // **A round's ground is read on a thread of its own, not as a task on the pool.** Its reading
    // is parallel over samples, and a pool thread waiting for its share of that work takes up
    // whatever task is queued — a whole region — so a read run as a task took 4 to 6 s a round
    // where the round driver's took 0.6. A thread outside the pool only waits.
    let (ask_to_read, to_read) = mpsc::channel::<(u64, GenomeRegion, GenomeRegion)>();
    let first_failure = std::thread::scope(|threads| {
        let reader_told = told.clone();
        threads.spawn(move || {
            let read_pool = ReadPool::beside(rayon::current_num_threads());
            for (round, first, last) in to_read {
                let sealed = read_pool.run(|| read_and_seal(&mut *cache, first, last));
                if reader_told.send(Done::Covered { round, sealed }).is_err() {
                    break;
                }
            }
        });
        let dispatcher = threads.spawn(move || {
            // **Workers' scratch, reused across regions** rather than made per region: a region is small
            // and a scratch is not.
            let scratches: Mutex<Vec<Sc>> = Mutex::new(Vec::new());
            let scratches = &scratches;
            let make_scratch = &make_scratch;
            let call_locus = &call_locus;
            let mut rounds: VecDeque<Round<T, E>> = VecDeque::new();
            let mut next_number: u64 = 0;
            let mut reading = false;
            let mut planning_done = false;
            let mut first_failure: Option<E> = None;
            let mut tasks_out: usize = 0;
            let mut regions_out: usize = 0;
            // **Regions go to the pool only a little ahead of the threads**, so that a read started
            // behind them waits for about one region rather than for a whole round's.
            let regions_out_at_most = rayon::current_num_threads() + REGIONS_QUEUED_AHEAD;

            rayon::in_place_scope(|scope| {
                loop {
                    // **Read the next round when nothing is being read, the run has not failed, and
                    // fewer than `ROUNDS_HELD` rounds still have regions to build.**
                    let held = rounds
                        .iter()
                        .filter(|round| round.covered && round.unbuilt > 0)
                        .count();
                    if !reading && first_failure.is_none() && !planning_done && held < ROUNDS_HELD {
                        if let Some(regions) = planned.next() {
                            let number = next_number;
                            next_number += 1;
                            let (first, last) = (regions[0], regions[regions.len() - 1]);
                            rounds.push_back(Round {
                                number,
                                outcomes: (0..regions.len()).map(|_| None).collect(),
                                unbuilt: regions.len(),
                                spawned: 0,
                                regions,
                                sealed: None,
                                covered: false,
                                cover_failed: None,
                            });
                            reading = true;
                            tasks_out += 1;
                            let _ = ask_to_read.send((number, first, last));
                        } else {
                            planning_done = true;
                        }
                    }

                    // **Hand regions out, oldest round first.**
                    if first_failure.is_none() {
                        for state in rounds.iter_mut() {
                            let Some(sealed) = state.sealed.as_ref() else {
                                continue;
                            };
                            while state.spawned < state.regions.len()
                                && regions_out < regions_out_at_most
                            {
                                let at = state.spawned;
                                let region = state.regions[at];
                                let round = state.number;
                                let sealed = Arc::clone(sealed);
                                let told = told.clone();
                                state.spawned += 1;
                                regions_out += 1;
                                tasks_out += 1;
                                scope.spawn(move |_| {
                                    let mut scratch = scratches
                                        .lock()
                                        .expect("no worker panicked holding the scratches")
                                        .pop()
                                        .unwrap_or_else(make_scratch);
                                    let outcome = build_one_region(
                                        &sealed,
                                        region,
                                        max_cohort_locus_span,
                                        min_alt_reads,
                                        &mut scratch,
                                        call_locus,
                                    );
                                    scratches
                                        .lock()
                                        .expect("no worker panicked holding the scratches")
                                        .push(scratch);
                                    drop(sealed);
                                    let _ = told.send(Done::Built { round, at, outcome });
                                });
                            }
                        }
                    }

                    // **Hand over every finished round at the front, in order.**
                    while first_failure.is_none() {
                        let Some(front) = rounds.front_mut() else {
                            break;
                        };
                        if let Some(failed) = front.cover_failed.take() {
                            first_failure = Some(failed);
                            break;
                        }
                        if !front.covered || front.unbuilt > 0 {
                            break;
                        }
                        let finished = rounds.pop_front().expect("a front, just seen");
                        for outcome in finished.outcomes {
                            match outcome.expect("every region of a finished round has an outcome")
                            {
                                Ok(called) => {
                                    // The calling thread receives until this thread is done.
                                    let _ = to_hand_over.send(called);
                                }
                                Err(failed) => {
                                    first_failure = Some(failed);
                                    break;
                                }
                            }
                        }
                    }
                    if first_failure.is_some() {
                        // Nothing new is started; what is out finishes and is dropped.
                        rounds.clear();
                    }

                    if tasks_out == 0
                        && (first_failure.is_some() || (planning_done && rounds.is_empty()))
                    {
                        break;
                    }

                    let Ok(done) = hear.recv() else {
                        unreachable!("the calling thread holds a sender, so the channel stays open")
                    };
                    tasks_out -= 1;
                    match done {
                        Done::Covered { round, sealed } => {
                            reading = false;
                            // A round dropped after a failure: its seal is not wanted.
                            let Some(state) = rounds.iter_mut().find(|state| state.number == round)
                            else {
                                continue;
                            };
                            state.covered = true;
                            match sealed {
                                Err(failed) => {
                                    state.cover_failed = Some(failed);
                                    state.unbuilt = 0;
                                }
                                Ok(sealed) => state.sealed = Some(sealed),
                            }
                        }
                        Done::Built { round, at, outcome } => {
                            regions_out -= 1;
                            let Some(state) = rounds.iter_mut().find(|state| state.number == round)
                            else {
                                continue;
                            };
                            state.outcomes[at] = Some(outcome);
                            state.unbuilt -= 1;
                            if state.unbuilt == 0 {
                                state.sealed = None;
                            }
                        }
                    }
                }
            });
            first_failure
        });
        for called in arrived {
            hand_over(called);
        }
        dispatcher
            .join()
            .unwrap_or_else(|panicked| std::panic::resume_unwind(panicked))
    });

    match first_failure {
        Some(failed) => Err(failed),
        None => Ok(displaced),
    }
}

/// **Threads of their own for reading a round's ground**, beside the pool the regions are built
/// on — for the reason the module gives.
///
/// **Its size**: [`read_threads`] of the run's workers, **in addition to them**, for the reason the
/// hand-over's pool is (`crate::run::psp_caller`): they are busy only while a round is read.
struct ReadPool(Option<rayon::ThreadPool>);

impl ReadPool {
    /// The pool beside `workers` region threads, or none — reads then use the regions' pool —
    /// when [`read_threads`] gives none or the threads cannot be started.
    fn beside(workers: usize) -> Self {
        let threads = read_threads(workers);
        if threads == 0 {
            return Self(None);
        }
        Self(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .thread_name(|at| format!("pvc-read-{at}"))
                .build()
                .ok(),
        )
    }

    /// Run `work` on the pool, or where it was called when there is none.
    fn run<R: Send>(&self, work: impl FnOnce() -> R + Send) -> R {
        match &self.0 {
            Some(pool) => pool.install(work),
            None => work(),
        }
    }
}

/// **How many read threads a run of `workers` region threads gets**: one for every three, at most
/// sixteen, and none below twelve workers, where a round's regions are few and a read waits little
/// behind them. Sixteen at 48 workers is the point measured: over 300 kb at 2,176 samples, with two
/// rounds held, 16 read threads called in 52 s and 8 in 57 s; the other points are not measured.
fn read_threads(workers: usize) -> usize {
    if workers < 12 {
        return 0;
    }
    (workers / 3).clamp(1, 16)
}

/// Evict behind a round, read its ground, and seal what every sample holds.
fn read_and_seal<S, E>(
    cache: &mut ObservationCache<S>,
    first: GenomeRegion,
    last: GenomeRegion,
) -> Result<Arc<RoundView<E>>, E>
where
    S: ObservationSource<Error = E> + Sync + Send,
    E: Send + 'static + From<ReferenceUnreadable>,
{
    let mut dead: Vec<SampleLocusObservations> = Vec::new();
    cache.evict_before_in_parallel(
        GenomePosition {
            contig: first.contig,
            position: first.start,
        },
        &mut dead,
    );
    drop(dead);
    cache.cover_in_parallel(GenomeRegion {
        contig: first.contig,
        start: first.start,
        end: last.end,
    })?;
    Ok(Arc::new(cache.seal_round(GenomePosition {
        contig: last.contig,
        position: Position(last.end.0 + 1),
    })))
}

/// Build and call one region from its round's seal, as the round driver's worker does.
fn build_one_region<E, Sc, T, CL>(
    sealed: &RoundView<E>,
    region: GenomeRegion,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
    scratch: &mut Sc,
    call_locus: &CL,
) -> Result<CalledRegion<T>, E>
where
    CL: Fn(&mut Sc, CohortObservation, &LocusReads<'_, '_>) -> T,
{
    let mut calls: Vec<(GenomeRegion, T)> = Vec::new();
    let mut failed: Vec<GenomeRegion> = Vec::new();
    // The window's own prefix, added back to every member index, for the reason the round
    // driver's worker gives.
    let mut starts: Vec<usize> = Vec::new();
    sealed.window_starts_into(region, &mut starts);
    let starts = &starts;
    sealed.with_observations(region, |window| {
        build_region_handing_over_windowed(
            region,
            window,
            max_cohort_locus_span,
            min_alt_reads,
            &mut |built: CohortObservation, reads: &LocusReads<'_, '_>| {
                let span = built.region;
                calls.push((span, call_locus(scratch, built, reads)));
            },
            &mut failed,
            &|sample, index| sealed.build_at(sample, starts[sample] + index),
            sealed.dropped_loci(),
        )
    })?;
    Ok(CalledRegion { calls, failed })
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{
        KeepingSource, SourceFailed, in_flight, member, region, source_of,
        three_samples_over_six_hundred_bases, width,
    };
    use super::super::serial::merge_cohort_handing_each_locus_over;
    use super::*;

    /// What a driver made of a fixture: each locus rendered, and the spans the width bound
    /// refused — the round driver's tests' shape, for their reason.
    #[derive(Debug, PartialEq, Eq)]
    struct WhatCameOut {
        loci: Vec<String>,
        refused: Vec<GenomeRegion>,
    }

    /// One region at a time, which every driver must reproduce.
    fn one_at_a_time(
        analysed: &[GenomeRegion],
        per_sample: &[&[SampleLocusObservations]],
        region_width: CohortLocusBuilderRegionsLen,
    ) -> WhatCameOut {
        let mut cache = ObservationCache::over_fixture(
            per_sample.iter().map(|sample| source_of(sample)).collect(),
        );
        let mut loci = Vec::new();
        let mut refused = Vec::new();
        merge_cohort_handing_each_locus_over(
            analysed,
            &mut cache,
            region_width,
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
            &mut |built, _| loci.push(format!("{built:?}")),
            &mut refused,
        )
        .expect("the fixture sources hold");
        WhatCameOut { loci, refused }
    }

    /// The streaming driver over the same fixture, through `cache`.
    fn streamed<S>(
        analysed: &[GenomeRegion],
        mut cache: ObservationCache<S>,
        region_width: CohortLocusBuilderRegionsLen,
        regions: usize,
    ) -> WhatCameOut
    where
        S: ObservationSource<Error = super::super::fixtures::SourceFailed> + Send + Sync,
    {
        let mut loci = Vec::new();
        let mut refused = Vec::new();
        let displaced = merge_and_call_streaming(
            analysed,
            &mut cache,
            region_width,
            in_flight(regions),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
            || (),
            |(), built, _| format!("{built:?}"),
            &mut |rendered| loci.push(rendered),
            &mut refused,
        )
        .expect("the fixture sources hold");
        assert_eq!(displaced, 0, "a fixture merge displaced a locus");
        WhatCameOut { loci, refused }
    }

    fn keeping(per_sample: &[&[SampleLocusObservations]]) -> ObservationCache<KeepingSource> {
        ObservationCache::over_fixture(
            per_sample
                .iter()
                .map(|sample| KeepingSource::over(sample))
                .collect(),
        )
    }

    /// **The claim of the whole arrangement**: streamed, the merge calls what one region at a
    /// time calls, at every width and round size, over sources that keep their evidence (the
    /// psp path's shape, built from each round's seal) and over sources that build as they draw.
    /// The widths straddle the fixture's 600 bases and put a boundary inside its 305–330
    /// deletion; a round of one region holds the fewest records a round can.
    #[test]
    fn the_stream_calls_what_one_region_at_a_time_calls() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 600)];
        for bases in [20, 100, 305, 610] {
            let expected = one_at_a_time(&analysed, &per_sample, width(bases));
            for regions in [1, 2, 3, 7, 16] {
                assert_eq!(
                    expected,
                    streamed(&analysed, keeping(&per_sample), width(bases), regions),
                    "kept evidence: {regions} regions of {bases} bases",
                );
                let building = ObservationCache::over_fixture(
                    per_sample.iter().map(|sample| source_of(sample)).collect(),
                );
                assert_eq!(
                    expected,
                    streamed(&analysed, building, width(bases), regions),
                    "built records: {regions} regions of {bases} bases",
                );
            }
        }
    }

    /// **Each region's members come from that region's own bodies**, the round driver's
    /// regression for the window's prefix, through seals: one record a region, six regions, and
    /// rounds of one, two and six, so later rounds are built from seals whose first records were
    /// carried over from the round before.
    #[test]
    fn a_streamed_region_builds_its_members_from_its_own_bodies() {
        let sample: Vec<SampleLocusObservations> = (0..6)
            .map(|at| {
                let base = 100 * at + 50;
                member(
                    region(base, base),
                    b"A",
                    if at % 2 == 0 { b"T" } else { b"C" },
                )
            })
            .collect();
        let analysed = [region(1, 600)];
        let per_sample: [&[SampleLocusObservations]; 1] = [&sample];
        let expected = one_at_a_time(&analysed, &per_sample, width(100));
        assert_eq!(expected.loci.len(), 6, "the fixture holds six loci");
        for regions in [1, 2, 6] {
            assert_eq!(
                expected,
                streamed(&analysed, keeping(&per_sample), width(100), regions),
                "rounds of {regions}",
            );
        }
    }

    /// A keeping source whose draw fails at its record `fail_at`, as a damaged file does.
    struct FailingAt {
        inner: KeepingSource,
        drawn: usize,
        fail_at: usize,
    }

    impl ObservationSource for FailingAt {
        type Error = SourceFailed;

        fn next_observation(
            &mut self,
            _spare: Option<SampleLocusObservations>,
        ) -> Option<Result<SampleLocusObservations, SourceFailed>> {
            unreachable!("a keeping source is drawn through `next_drawn`")
        }

        fn next_drawn(
            &mut self,
            spare: Option<SampleLocusObservations>,
        ) -> Option<Result<super::super::observation_cache::Drawn, SourceFailed>> {
            if self.drawn == self.fail_at {
                return Some(Err(SourceFailed("the file is damaged here")));
            }
            self.drawn += 1;
            self.inner.next_drawn(spare)
        }

        fn build(
            &self,
            body: core::ops::Range<usize>,
        ) -> Result<SampleLocusObservations, SourceFailed> {
            self.inner.build(body)
        }

        fn seal(
            &mut self,
            from_body: usize,
        ) -> Option<Arc<dyn super::super::observation_cache::BuildsKeptEvidence<SourceFailed>>>
        {
            self.inner.seal(from_body)
        }
    }

    /// **A read that fails ends the run with its error, after the loci the round driver hands
    /// over before it** — the same loci, in the same order, whatever the round size and wherever
    /// in the file the failure lies. Reading ahead must not change what was handed over first.
    #[test]
    fn a_failed_read_is_reported_after_the_loci_the_round_driver_hands_over() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 600)];
        let failing = |fail_at: usize| -> ObservationCache<FailingAt> {
            ObservationCache::over_fixture(
                per_sample
                    .iter()
                    .enumerate()
                    .map(|(at, sample)| FailingAt {
                        inner: KeepingSource::over(sample),
                        drawn: 0,
                        fail_at: if at == 0 { fail_at } else { usize::MAX },
                    })
                    .collect(),
            )
        };
        let mut some_handed_over = false;
        for fail_at in [5, 30, 55] {
            for regions in [1, 3, 7] {
                let mut by_rounds = Vec::new();
                let mut refused = Vec::new();
                let rounds_answer = super::super::rounds::merge_and_call_in_rounds(
                    &analysed,
                    &mut failing(fail_at),
                    width(20),
                    in_flight(regions),
                    MaxCohortLocusSpan::DEFAULT,
                    MinAltReads::DEFAULT,
                    || (),
                    |(), built, _| format!("{built:?}"),
                    &mut |rendered| by_rounds.push(rendered),
                    &mut refused,
                );
                let mut streamed_loci = Vec::new();
                let mut refused = Vec::new();
                let answer = merge_and_call_streaming(
                    &analysed,
                    &mut failing(fail_at),
                    width(20),
                    in_flight(regions),
                    MaxCohortLocusSpan::DEFAULT,
                    MinAltReads::DEFAULT,
                    || (),
                    |(), built, _| format!("{built:?}"),
                    &mut |rendered| streamed_loci.push(rendered),
                    &mut refused,
                );
                assert!(
                    matches!(rounds_answer, Err(SourceFailed("the file is damaged here"))),
                    "the round driver met the failure",
                );
                assert!(
                    matches!(answer, Err(SourceFailed("the file is damaged here"))),
                    "record {fail_at}, rounds of {regions}: the read's failure is the answer",
                );
                assert_eq!(
                    streamed_loci, by_rounds,
                    "record {fail_at}, rounds of {regions}: handed over before the failure",
                );
                some_handed_over |= !streamed_loci.is_empty();
            }
        }
        assert!(
            some_handed_over,
            "no case handed anything over before failing"
        );
    }

    /// **Several analysed intervals are one run**, rounds never crossing one.
    #[test]
    fn several_analysed_intervals_stream_as_one_run() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 300), region(301, 600)];
        for regions in [1, 3, 8] {
            assert_eq!(
                one_at_a_time(&analysed, &per_sample, width(50)),
                streamed(&analysed, keeping(&per_sample), width(50), regions),
                "rounds of {regions}",
            );
        }
    }
}
