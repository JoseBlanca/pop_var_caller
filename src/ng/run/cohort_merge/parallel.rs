//! The parallel arrangement: rounds of building regions, worked concurrently.
//!
//! [`merge_cohort_in_parallel`] is [`merge_cohort_through_cache`](super::serial::merge_cohort_through_cache)
//! with the builders running at the same time as each other, and its whole claim is that the
//! output is unchanged — which is the oracle's, byte for byte
//! (`doc/devel/ng/spec/cohort_merge.md` §15).
//!
//! **A round is what the design costs and what it buys** (spec §6.2, amended by the owner on
//! 2026-08-18). A builder reads the observation cache while the organiser draws the readers
//! forward, so the two cannot run at once: the organiser covers the ground for a round of
//! regions, the round's builders run, and only when every one of them has finished does it
//! resolve, release and evict. Every builder therefore waits for the slowest in its own round.
//! That is a bounded frontier — one round spans
//! `cohort_locus_builder_regions_in_flight × cohort_locus_builder_regions_len` bases, 3,200 at
//! 16 regions in flight on the default 200-base regions — where production's watermark spanned
//! the whole run.
//! What it buys is one reader per sample instead of one per sample per region in flight: 3,000
//! cursors at a cohort of 3,000 samples and 16 regions in flight, against 48,000.
//!
//! **Nothing in this module's tests can tell a parallel merge from a sequential one**, and that
//! is not a gap that can be closed by a fixture: the whole claim of the step is that the answer
//! is the same either way. What stands in for it is the type of
//! [`in_region_order`] — replacing `par_iter` with `iter` stops compiling — and the memory the
//! round holds, which `tests::the_cache_holds_a_whole_round_and_not_one_region` pins to an exact
//! table. Whether the builders genuinely occupy several threads is untested here.

use rayon::prelude::*;

use super::build::{RegionOutcome, build_region};
use super::observation_cache::{ObservationCache, building_regions_of};
use super::organise::{Organiser, RegionIndex};
use super::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan,
    MinAltReads, refuse_malformed_analysed_regions,
};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion};

/// Merge the cohort over `analysed` with `regions_in_flight` building regions worked at once,
/// reading each sample through the observation cache.
///
/// **The output is the serial drivers' and that is the whole claim of this step**: same cohort
/// observations, same failed spans, at any number of regions in flight and any region width
/// (spec §15). Nothing here decides anything about the answer — the builders are
/// [`build_region`], the ordering is the organiser's, and both are already the serial path's.
///
/// **`regions_in_flight` is how many regions are worked at once, not how many threads work
/// them** ([`CohortLocusBuilderRegionsInFlight`]). Threads come from rayon's pool; what this
/// number sets is the ground the cache must hold — `regions_in_flight ×
/// cohort_locus_builder_regions_len` bases plus the tail of observations reaching past it
/// (spec §6.4). One region in flight gives the serial driver's memory and its answer.
///
/// **A round never crosses an analysed region.** The run's intervals are disjoint and may sit on
/// different contigs, and the eviction that opens a round is at its first region's first base —
/// a round spanning two intervals would evict on the strength of ground the second interval does
/// not contain. The cost is one short round per interval, and it is small: an interval holds
/// `its bases ÷ the region width` regions, and the last round leaves idle whatever part of
/// `regions_in_flight` that count does not fill. A 600-base interval on 20-base regions is 30
/// regions, which at 16 in flight is one full round and one of 14 — **two builders idle, once
/// per interval**. The worst any interval can do is `regions_in_flight - 1` idle, once.
///
/// **`S: Sync` is what sharing the cache costs, and it is broader than what the builders touch.**
/// They read only each sample's held observations, and never the readers themselves; the bound
/// falls out of `&ObservationCache<S>` crossing threads. So a source built on `Rc` or `RefCell`
/// cannot be merged in parallel, which `run_streaming.md` arch §2's `ObservationSource` has to
/// honour. Splitting the held windows out of the cache into a view type would drop the bound,
/// and is not worth a type until some source needs it.
///
/// **The failed spans are gathered here and the organiser's count is the cross-check.** The
/// organiser resolves overlaps and counts the failures that survived (spec §6.1, §3.3) but does
/// not hand the spans back, so this keeps the spans of every outcome it submitted. The two agree
/// unless a failed locus was displaced, which cannot happen under [`build_region`]'s input
/// contract — the argument is on [`Organiser`].
///
/// **The two assertions at the end are safety nets that no test makes fire**, which is worth
/// saying rather than leaving a reader to look for the fixture. Every merge that returns `Ok`
/// evaluates them; what no fixture does is trip either, because both can only fire if the
/// organiser displaced something and nothing a builder produces can make it. Mutation testing
/// says the same from the other side: disabling either leaves the whole suite green, because
/// every defect that would trip them — spans not gathered, outcomes submitted out of order — is
/// already caught by the output being compared against the serial drivers'.
///
/// The cache comes in from the caller for [`merge_cohort_through_cache`](super::serial::merge_cohort_through_cache)'s
/// reason: it holds one reader per sample for the whole run. **A failure ends the merge and
/// leaves the cache advanced**, so the same cache cannot be used to try the same ground again —
/// also that driver's caveat, and still owed.
///
/// **A panicking builder is the same caveat, and it is reachable.** [`build_region`]'s
/// producer-guarantee assertions fire on a rayon worker: two overlapping records of one sample
/// pass the cache, which checks only that starts do not go backwards, and trip `build.rs`'s
/// disjointness check. Rayon joins every builder in the round before it re-raises the panic, so
/// the cache is unborrowed when the caller regains it — what it is *not* is rewound, and
/// [`Organiser::finish`] never runs, so the run-ended-short guard is skipped too. A caller that
/// catches the unwind must abandon the cache for the same reason a failure's caller must.
pub fn merge_cohort_in_parallel<S, E>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    regions_in_flight: CohortLocusBuilderRegionsInFlight,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
) -> Result<RegionOutcome, E>
where
    S: Iterator<Item = Result<SampleLocusObservations, E>> + Sync,
{
    let mut merged = RegionOutcome::default();
    refuse_malformed_analysed_regions(analysed);

    let mut organiser = Organiser::new();
    let mut regions_handed_out = 0u64;
    // Grown rather than reserved: `regions_in_flight` is a knob a caller sets, and reserving on
    // it would let a merge over thirty regions try to allocate for `usize::MAX` of them.
    let mut regions_in_round: Vec<GenomeRegion> = Vec::new();

    for analysed_region in analysed {
        let mut building_regions =
            building_regions_of(*analysed_region, cohort_locus_builder_regions_len);
        loop {
            // The clear is what ends the loop as well as what keeps the buffer honest: without
            // it the round never empties and `first()` is always `Some`.
            regions_in_round.clear();
            regions_in_round.extend(building_regions.by_ref().take(regions_in_flight.get()));
            let Some(first) = regions_in_round.first().copied() else {
                break;
            };

            // **Evict, then cover the whole round, then let it run.** What ends before the
            // round's first base can reach no locus any builder in the round can own, and
            // covering the ground before any builder starts is what lets them all read a
            // cache nobody is drawing forward.
            //
            // **One cover for the round, not one per region.** A cover draws every sample
            // forward until the chain from its region's last base closes, and the chain from
            // the round's last base contains every one of those — so covering the round's
            // whole span draws the same observations and leaves the same window, in one pass
            // over the cohort instead of `regions_in_flight` of them. Measured on 8 threads
            // with 8 regions in flight: 11–13% off the whole merge at 3,000 samples on both
            // the densities `examples/ng_cohort_merge_parallel_cost.rs` walks.
            cache.evict_before(GenomePosition {
                contig: first.contig,
                position: first.start,
            });
            let last = regions_in_round
                .last()
                .copied()
                .expect("the round is not empty");
            cache.cover(GenomeRegion {
                contig: first.contig,
                start: first.start,
                end: last.end,
            })?;

            let cache = &*cache;
            let outcomes = in_region_order(regions_in_round.par_iter().map(|building_region| {
                cache.with_observations(*building_region, |observations_per_sample| {
                    build_region(
                        *building_region,
                        observations_per_sample,
                        max_cohort_locus_span,
                        min_alt_reads,
                    )
                })
            }));

            for outcome in outcomes {
                // Destructured for the reason `organise.rs` gives at its own two consumers of
                // this type: a field `RegionOutcome` gains has to be answered for here —
                // gathered or dropped deliberately — rather than starting at its `Default` and
                // staying there for the whole run.
                let RegionOutcome {
                    cohort_observations: _,
                    failed_locus_spans,
                } = &outcome;
                merged.failed_locus_spans.extend(failed_locus_spans);
                organiser.submit(RegionIndex(regions_handed_out), outcome);
                regions_handed_out += 1;
            }
            merged.cohort_observations.extend(organiser.drain_ready());
        }
    }

    let tally = organiser
        .finish(regions_handed_out)
        .expect("every region handed out was submitted and every locus drained");
    assert_eq!(
        tally.failed_loci,
        merged.failed_locus_spans.len() as u64,
        "the organiser kept {} of the {} refused spans this merge gathered, so one was \
         displaced by an earlier locus and the two no longer describe the same ground",
        tally.failed_loci,
        merged.failed_locus_spans.len(),
    );
    assert_eq!(
        tally.displaced_loci, 0,
        "a builder produced a locus on ground an earlier locus already owned",
    );

    Ok(merged)
}

/// The round's outcomes, in the order its regions were handed out.
///
/// **The bound is the point, not the body.** `collect` gives back the *input's* order whatever
/// order the builders finished in, but only for an **indexed** parallel iterator — so stating
/// that as a bound turns two silent regressions into compile errors: `par_iter` replaced by
/// `iter`, which would drop the parallelism and change nothing a test can see, and an adaptor
/// that gives up the index (`filter`, `flat_map`) inserted before the collect, which would
/// reorder the round and be caught only by whichever locus happened to move.
fn in_region_order<I>(builders: I) -> Vec<RegionOutcome>
where
    I: IndexedParallelIterator<Item = RegionOutcome>,
{
    builders.collect()
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{
        SourceFailed, in_flight, member, refuse_any_difference, region, region_on, source_of,
        three_samples_over_six_hundred_bases, width,
    };
    use super::super::observation_cache::building_regions_of;
    use super::super::serial::{merge_cohort_serially, merge_cohort_through_cache};
    use super::*;
    use crate::ng::types::Position;

    /// A cache over one reader per sample of `layouts`.
    fn cache_over(layouts: &[Vec<SampleLocusObservations>]) -> ObservationCache<SourceOfFixture> {
        ObservationCache::over(
            layouts
                .iter()
                .map(|sample| source_of(sample))
                .collect::<Vec<_>>(),
        )
    }

    type SourceOfFixture = std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>>;

    /// Merge a fixture in parallel and hand back its outcome.
    fn merge_fixture_in_parallel(
        analysed: &[GenomeRegion],
        layouts: &[Vec<SampleLocusObservations>],
        building_region_width: CohortLocusBuilderRegionsLen,
        regions_in_flight: CohortLocusBuilderRegionsInFlight,
    ) -> RegionOutcome {
        merge_cohort_in_parallel(
            analysed,
            &mut cache_over(layouts),
            building_region_width,
            regions_in_flight,
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("the fixture sources hold")
    }

    /// The oracle's answer over the same fixture.
    fn merge_fixture_serially(
        analysed: &[GenomeRegion],
        layouts: &[Vec<SampleLocusObservations>],
    ) -> RegionOutcome {
        let per_sample: Vec<&[SampleLocusObservations]> =
            layouts.iter().map(Vec::as_slice).collect();
        merge_cohort_serially(
            analysed,
            &per_sample,
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
    }

    /// **One region in flight is the serial driver**, and that is the base case the rest rests
    /// on: with a round of one region, this loop is `merge_cohort_through_cache`'s.
    #[test]
    fn one_region_in_flight_gives_the_serial_drivers_answer() {
        let layouts = three_samples_over_six_hundred_bases();
        let analysed = [region(1, 600)];

        refuse_any_difference(
            "one region in flight",
            &merge_fixture_serially(&analysed, &layouts),
            &merge_fixture_in_parallel(&analysed, &layouts, width(20), in_flight(1)),
        );
    }

    /// **And so is every other count of regions in flight**, which is the milestone's claim.
    /// Sixteen in flight on 20-base regions put a round over 320 bases at once, and the deletion
    /// at 305–330 sits inside one such round.
    #[test]
    fn the_answer_is_the_same_at_every_count_of_regions_in_flight() {
        let layouts = three_samples_over_six_hundred_bases();
        let analysed = [region(1, 600)];
        let oracle = merge_fixture_serially(&analysed, &layouts);

        for regions in [1, 2, 4, 8, 16] {
            refuse_any_difference(
                &format!("{regions} regions in flight"),
                &oracle,
                &merge_fixture_in_parallel(&analysed, &layouts, width(20), in_flight(regions)),
            );
        }
    }

    /// **A round holds a round's worth of ground, and that is what the count in flight buys and
    /// costs.** Eviction opens a round rather than a region, so the cache holds the whole round's
    /// window — the one claim of this step no comparison of outputs can see, since the outputs
    /// are identical at every count.
    ///
    /// **The exact table, not merely that it grows.** E3's review showed an inequality passes for
    /// a round of twice the count, or four times, or the whole analysed region: doubling the
    /// round takes the table from 2/4/4/12/29 to 4/4/12/29/62 and every `eight > one` still
    /// holds. Re-measure and rewrite these five if the fixture changes; they are what the
    /// arrangement's memory bound rests on.
    #[test]
    fn the_cache_holds_a_whole_round_and_not_one_region() {
        let layouts = three_samples_over_six_hundred_bases();
        let records_held_at = |regions: usize| {
            let mut cache = cache_over(&layouts);
            merge_cohort_in_parallel(
                &[region(1, 600)],
                &mut cache,
                width(20),
                in_flight(regions),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            )
            .expect("the fixture sources hold");
            cache.held_observations_len()
        };

        assert_eq!(
            [1, 2, 4, 8, 16].map(records_held_at),
            [2, 4, 4, 12, 29],
            "eviction is no longer following the round",
        );
    }

    /// **A round never crosses an analysed region**, on two contigs: a round of eight against
    /// four regions in the first interval would evict at a base of contig 0 and then cover
    /// ground on contig 1, which the cache refuses.
    #[test]
    fn a_round_stops_at_the_end_of_its_analysed_region() {
        let layouts = vec![vec![
            member(region_on(0, 5, 5), b"G", b"T"),
            member(region_on(0, 61, 61), b"G", b"T"),
            member(region_on(1, 7, 7), b"G", b"T"),
        ]];
        let analysed = [region_on(0, 1, 80), region_on(1, 1, 80)];

        let in_parallel = merge_fixture_in_parallel(&analysed, &layouts, width(20), in_flight(8));

        refuse_any_difference(
            "a round of eight over two intervals on two contigs",
            &merge_fixture_serially(&analysed, &layouts),
            &in_parallel,
        );
        assert_eq!(in_parallel.cohort_observations.len(), 3);
    }

    /// **And on one contig, where nothing refuses a round that ran on.** Two intervals with a
    /// gap between them, the first shorter than a 100-base region: a driver that let a round
    /// cross would evict at the wrong base and answer wrong, quietly. The previous test cannot
    /// catch that class — there the cache refuses the cover outright.
    #[test]
    fn a_round_stops_at_an_analysed_regions_end_on_one_contig_too() {
        let layouts = vec![vec![
            member(region(3, 3), b"G", b"T"),
            member(region(9, 12), b"GGGG", b"G"),
            member(region(41, 41), b"G", b"T"),
        ]];
        let analysed = [region(1, 15), region(30, 45)];
        let oracle = merge_fixture_serially(&analysed, &layouts);

        for regions in [1, 2, 3, 8] {
            for bases in [1, 4, 100] {
                refuse_any_difference(
                    &format!("{regions} regions in flight on {bases}-base regions"),
                    &oracle,
                    &merge_fixture_in_parallel(
                        &analysed,
                        &layouts,
                        width(bases),
                        in_flight(regions),
                    ),
                );
            }
        }
    }

    /// The failed spans come through unchanged, and the organiser's count agrees with them.
    #[test]
    fn the_failed_spans_come_through_the_parallel_driver_unchanged() {
        let layouts = vec![vec![
            member(region(10, 10), b"A", b"C"),
            member(region(30, 120), &[b'A'; 91], b"A"),
            member(region(200, 200), b"A", b"C"),
        ]];
        let analysed = [region(1, 300)];
        let oracle = merge_fixture_serially(&analysed, &layouts);
        assert_eq!(
            oracle.failed_locus_spans,
            vec![region(30, 120)],
            "the fixture must contain a locus the width bound refuses",
        );

        for regions in [1, 4, 16] {
            refuse_any_difference(
                &format!("{regions} regions in flight"),
                &oracle,
                &merge_fixture_in_parallel(&analysed, &layouts, width(7), in_flight(regions)),
            );
        }
    }

    /// A reader that fails ends the merge, and its own error comes back untouched.
    #[test]
    fn a_failing_source_ends_the_parallel_merge() {
        let failing = vec![
            Ok(member(region(5, 5), b"G", b"T")),
            Err(SourceFailed("the psp block would not decode")),
        ];
        let mut cache = ObservationCache::over(vec![failing.into_iter()]);

        let merged = merge_cohort_in_parallel(
            &[region(1, 600)],
            &mut cache,
            width(20),
            in_flight(4),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        assert_eq!(
            merged.err(),
            Some(SourceFailed("the psp block would not decode")),
        );
    }

    /// **A builder's refusal reaches the caller, and the cache is theirs again but not rewound.**
    /// Two overlapping records of one sample pass the cache — which checks only that starts do
    /// not go backwards — and trip `build_region`'s disjointness assertion on a rayon worker.
    /// Rayon joins every builder in the round before re-raising, which is what makes the unwind
    /// sound; what it does not do is put the readers back.
    #[test]
    fn a_builder_that_panics_leaves_the_cache_in_the_callers_hands_and_advanced() {
        let layouts = vec![vec![
            member(region(101, 110), &[b'A'; 10], b"A"),
            member(region(105, 106), &[b'A'; 2], b"T"),
        ]];
        let mut cache = cache_over(&layouts);

        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            merge_cohort_in_parallel(
                &[region(1, 600)],
                &mut cache,
                width(20),
                in_flight(16),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            )
        }));

        assert!(
            panicked.is_err(),
            "a builder's refusal must reach the caller"
        );
        assert_eq!(
            cache.held_observations_len(),
            2,
            "the cache is the caller's again — and advanced, not rewound",
        );
    }

    /// Nothing analysed is not an error, at any count of regions in flight.
    #[test]
    fn merging_no_analysed_regions_in_parallel_yields_nothing() {
        let merged = merge_fixture_in_parallel(
            &[],
            &[vec![member(region(12, 12), b"G", b"T")]],
            width(20),
            in_flight(4),
        );

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
    }

    /// A cohort of no samples yields nothing rather than failing.
    #[test]
    fn merging_no_samples_in_parallel_yields_nothing() {
        let merged = merge_fixture_in_parallel(&[region(1, 100)], &[], width(20), in_flight(4));

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
    }

    /// An analysed region with its ends the wrong way round is refused here as in both serial
    /// drivers — they read such a region differently, so no merge may accept one.
    #[test]
    #[should_panic(expected = "the wrong way round")]
    fn an_analysed_region_with_inverted_ends_is_refused_in_parallel() {
        merge_fixture_in_parallel(
            &[GenomeRegion {
                contig: crate::ng::types::ContigId(0),
                start: Position(50),
                end: Position(1),
            }],
            &[vec![member(region(12, 12), b"G", b"T")]],
            width(20),
            in_flight(4),
        );
    }

    /// Analysed regions that overlap are refused too — the second of the shared guard's two
    /// rules, which the inverted fixture above does not reach.
    #[test]
    #[should_panic(expected = "not disjoint and ascending")]
    fn analysed_regions_that_overlap_are_refused_in_parallel() {
        merge_fixture_in_parallel(
            &[region(1, 50), region(40, 90)],
            &[vec![member(region(12, 12), b"G", b"T")]],
            width(20),
            in_flight(4),
        );
    }

    /// **The cached serial driver and the parallel one agree at every width and every count in
    /// flight**, which is what says the round changed nothing the cache had not already settled.
    #[test]
    fn the_parallel_driver_matches_the_cached_one_at_every_width() {
        let layouts = three_samples_over_six_hundred_bases();
        let analysed = [region(1, 600)];

        for bases in [1, 3, 20, 47, 600] {
            let through_cache = merge_cohort_through_cache(
                &analysed,
                &mut cache_over(&layouts),
                width(bases),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            )
            .expect("the fixture sources hold");

            for regions in [1, 2, 4, 8, 16] {
                refuse_any_difference(
                    &format!("{regions} regions in flight on {bases}-base regions"),
                    &through_cache,
                    &merge_fixture_in_parallel(
                        &analysed,
                        &layouts,
                        width(bases),
                        in_flight(regions),
                    ),
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // Milestone E's claim, stated once: **dividing the genome among builders changes nothing
    // a caller sees.** Everything the milestone added — the ordered release, the overlap
    // resolution, the rounds — is about speed and memory, so anything that changes the answer
    // is a defect rather than a trade (spec §15).
    // ---------------------------------------------------------------

    /// **The whole cross-product, against the oracle**: five counts of regions in flight by
    /// five region widths, on ground carrying both shapes the division can get wrong.
    ///
    /// The comparison is with [`merge_cohort_serially`] and not with the cached driver, because
    /// the oracle is what the milestone undertakes to reproduce — the simplest thing that gives
    /// the right answer, holding every sample's observations at once and dividing nothing. The
    /// cached driver's own agreement is a separate claim, and `serial.rs` keeps it.
    ///
    /// **The fixture has to reach both shapes or the sweep proves nothing.** It carries a
    /// deletion at 305–330, which chains two samples into one locus, and a 91-base record the
    /// 50-base default bound refuses, which puts a failed span across many regions at the
    /// narrow widths. Both are asserted present in the oracle before the sweep, so the
    /// twenty-five comparisons cannot pass by comparing nothing.
    ///
    /// The deletion **reaches across a region boundary at four of the five widths** — it touches
    /// 26, 9, 2 and 2 building regions at widths 1, 3, 20 and 47, and 1 at width 600, where the
    /// whole analysed stretch is a single region and nothing is divided. **That row is the
    /// undivided control the other four are read against**, and it is also where one round
    /// stands fifteen-sixteenths idle at 16 in flight.
    ///
    /// **Those counts are asserted per width, not inferred**, because the two shape assertions
    /// below are made against the oracle — which divides nothing — so neither can fail if the
    /// width list is edited to widths that stop dividing the fixture. Without the per-width
    /// check the sweep could come to compare five undivided merges against an undivided oracle
    /// and still report the milestone proved.
    ///
    /// A width of 600 makes one region of the whole analysed stretch, so at 16 in flight the
    /// round is fifteen-sixteenths idle; a width of 1 makes 600 regions, so the round is always
    /// full. Both ends are in the sweep on purpose.
    #[test]
    fn the_parallel_merge_is_the_oracles_at_every_width_and_count() {
        let mut layouts = three_samples_over_six_hundred_bases();
        layouts.push(vec![member(region(420, 510), &[b'A'; 91], b"A")]);
        let analysed = [region(1, 600)];

        let oracle = merge_fixture_serially(&analysed, &layouts);
        assert!(
            oracle
                .cohort_observations
                .iter()
                .any(|observed| observed.region == region(305, 330)),
            "the fixture must carry a locus reaching across a building-region boundary",
        );
        assert_eq!(
            oracle.failed_locus_spans,
            vec![region(420, 510)],
            "the fixture must carry a locus the width bound refuses",
        );
        assert_eq!(
            oracle.cohort_observations.len(),
            50,
            "60 dotted loci, less the two the deletion swallowed, less the nine the refused \
             91-base record at 420–510 absorbed — pinned so that an edit to the shared fixture \
             cannot leave the sweep agreeing on a nearly empty answer",
        );

        for (bases, regions_over_the_deletion) in [(1, 26), (3, 9), (20, 2), (47, 2), (600, 1)] {
            assert_eq!(
                building_regions_of(analysed[0], width(bases))
                    .filter(
                        |building| building.start <= Position(330) && building.end >= Position(305)
                    )
                    .count(),
                regions_over_the_deletion,
                "at {bases}-base regions the sweep no longer divides the straddling locus the \
                 way this test's doc claims",
            );
            for regions in [1, 2, 4, 8, 16] {
                refuse_any_difference(
                    &format!("{regions} regions in flight on {bases}-base regions"),
                    &oracle,
                    &merge_fixture_in_parallel(
                        &analysed,
                        &layouts,
                        width(bases),
                        in_flight(regions),
                    ),
                );
            }
        }
    }

    /// **A round asked for more regions than the interval holds is the interval's round**, and
    /// the count in flight is never turned into an allocation — the driver grows the round's
    /// buffer rather than reserving on a number a caller sets, and nothing else pinned that.
    #[test]
    fn merging_with_more_regions_in_flight_than_the_interval_holds_gives_the_oracles_answer() {
        let layouts = three_samples_over_six_hundred_bases();
        let analysed = [region(1, 600)];

        refuse_any_difference(
            "usize::MAX regions in flight on one-base regions",
            &merge_fixture_serially(&analysed, &layouts),
            &merge_fixture_in_parallel(&analysed, &layouts, width(1), in_flight(usize::MAX)),
        );
    }

    /// **The parallel driver agrees with the oracle on random layouts, at a random count in
    /// flight** — the shape the five hand-written fixtures cannot enumerate: a locus opening in
    /// a round's last region and reaching into the next, an interval whose region count is one
    /// more than the count in flight, a wide record straddling an analysed edge inside a short
    /// round. The generator is `serial.rs`'s, with the count in flight drawn too.
    #[test]
    fn the_parallel_driver_agrees_with_the_oracle_on_random_layouts() {
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

        let ground_end = 400u64;
        let (mut with_observations, mut with_failed, mut with_a_straddling_record) = (0, 0, 0);

        for seed in 0..200u64 {
            let mut draw = Seeded(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0x5EED);
            let samples = 1 + draw.next(6) as usize;
            let contigs = 1 + draw.next(2) as u32;

            let mut layouts: Vec<Vec<SampleLocusObservations>> = Vec::new();
            for _ in 0..samples {
                let mut records = Vec::new();
                for contig in 0..contigs {
                    let mut at_base = 1 + draw.next(10);
                    while at_base <= ground_end {
                        let bases = match draw.next(10) {
                            0 => 1 + draw.next(60),
                            _ => 1 + draw.next(4),
                        };
                        let end = at_base + bases - 1;
                        let span = usize::try_from(end - at_base + 1).expect("small");
                        records.push(SampleLocusObservations {
                            reference_bases: vec![b'A'; span].into_boxed_slice(),
                            ..member(region_on(contig, at_base, end), b"A", b"T")
                        });
                        at_base = end + 1 + draw.next(6);
                    }
                }
                layouts.push(records);
            }

            let analysed_width = 20 + draw.next(80);
            let mut analysed = Vec::new();
            for contig in 0..contigs {
                let mut at_base = 1u64;
                while at_base <= ground_end {
                    let end = at_base + analysed_width - 1;
                    if draw.next(4) != 0 {
                        analysed.push(region_on(contig, at_base, end));
                    }
                    at_base = end + 1;
                }
            }
            let bases = 1 + u32::try_from(draw.next(30)).expect("a small width");
            let regions = 1 + draw.next(20) as usize;

            if layouts.iter().flatten().any(|record| {
                analysed.iter().any(|ground| {
                    ground.contig == record.region.contig
                        && record.region.start >= ground.start
                        && record.region.start <= ground.end
                        && record.reach() > ground.end
                })
            }) {
                with_a_straddling_record += 1;
            }

            let oracle = merge_fixture_serially(&analysed, &layouts);
            if !oracle.cohort_observations.is_empty() {
                with_observations += 1;
            }
            if !oracle.failed_locus_spans.is_empty() {
                with_failed += 1;
            }

            refuse_any_difference(
                &format!("seed {seed}: {regions} regions in flight on {bases}-base regions"),
                &oracle,
                &merge_fixture_in_parallel(&analysed, &layouts, width(bases), in_flight(regions)),
            );
        }

        assert!(
            with_observations > 100 && with_failed > 50 && with_a_straddling_record > 100,
            "the generator stopped producing the shapes this test is for: {with_observations} \
             layouts with observations, {with_failed} with failed spans, \
             {with_a_straddling_record} with a record straddling an analysed edge",
        );
    }
}
