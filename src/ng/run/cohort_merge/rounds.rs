//! **A round of building regions, each region built *and called* on a worker**, with the
//! results folded back in genome order.
//!
//! This is [`merge_cohort_in_parallel`](super::parallel::merge_cohort_in_parallel)'s round
//! structure with one difference, and the difference is what the arrangement is for: the
//! worker does not stop at a [`CohortObservation`]. It hands
//! each locus straight to a caller-supplied closure while the observation is still on the
//! worker's own stack, so the round never materialises a `Vec<CohortObservation>` and the
//! **genotyping runs on the pool beside the body decode**. In a calling run the genotyping is
//! the larger half of the two, and it was the whole reason the calling thread stayed busy
//! while the pool slept.
//!
//! Measured on 63 tomato accessions over 8 Mb of SL4.0 at about three reads a position, 18
//! logical cores: the calling pass falls from 74.8 s to 13.5 s and the whole command from
//! 83.5 s to 23.0 s, with the VCF body identical — 199,641 records.
//!
//! **What it costs is the ground a round holds**, which is
//! `regions_in_flight × cohort_locus_builder_regions_len` bases where the streaming driver
//! holds one region's worth. So the two knobs are one lever and a caller that wants today's
//! memory divides today's width by the count it asks for
//! (`crate::pop_var_caller_exp::calling_run::round_shape_for`). Sized that way, peak resident rises 4% to
//! 13% across cohorts of 1, 4, 8, 16, 32 and 63 accessions.
//!
//! **Ordering.** `regions_in_round.par_iter()` is indexed, so `collect` gives the round's
//! outcomes in region order whatever order the workers finished in; rounds are taken in genome
//! order; so the calls reach `take` in genome order. The organiser's reorder buffer is not
//! needed here for that reason — what is kept from it is the frontier rule (spec §6.1), which
//! is replicated below.
//!
//! **What stays on the calling thread**, and each for a reason that is about order rather than
//! about cost: the sink, because it writes; the eviction and the cover, because a builder reads
//! the cache while the cover draws the readers forward and spec §6.2 forbids the two at once;
//! and the identity of the first failure, because "the first refusal is the run's answer" is a
//! statement about genome order.

use rayon::prelude::*;

use super::build::{CohortObservation, build_region_handing_over_windowed};
use super::observation_cache::{
    ObservationCache, ObservationSource, ReferenceUnreadable, building_regions_of,
};
// **The organiser's own frontier helpers, not a second pair.** The rule below is
// `Organiser::resolve_and_release`'s, and two spellings of "where does this span begin" is
// exactly the drift this module keeps its shared fixtures in one place to avoid.
use super::organise::{first_base_of, last_base_of};
use super::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan,
    MinAltReads, refuse_malformed_analysed_regions,
};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion};

/// What one region's worker produced: what the closure made of each locus, in genome order,
/// and the spans the width bound refused.
pub struct CalledRegion<T> {
    /// Each locus's span beside what the closure made of it.
    pub calls: Vec<(GenomeRegion, T)>,
    /// The spans this region's builder refused, in genome order.
    pub failed: Vec<GenomeRegion>,
}

/// Merge and call the cohort over `analysed`, `regions_in_flight` regions at a time.
///
/// `make_scratch` builds one worker's mutable state; rayon calls it once per task split, not
/// once per locus. `call_locus` is what a worker does with a locus. `take` is the sink, on the
/// calling thread, in genome order.
///
/// **Returns how many loci the frontier rule dropped**, which is
/// [`Organiser::displaced_locus_count`](super::organise::Organiser::displaced_locus_count)'s
/// number and is kept for that type's reason: spec §6.1's argument says a builder cannot
/// produce a locus on ground an earlier locus already owns — given that eviction happens at the
/// round's first region's first base, which is what this driver does — and the count is how a
/// run would say the argument had failed. **It is a count and not an assertion on purpose**:
/// the organiser answers the same question by counting, and a run that has walked a genome
/// should report a surprise rather than abort on it.
///
/// # Errors
///
/// The cover's, and whatever a builder's decode refuses.
#[expect(
    clippy::too_many_arguments,
    reason = "the round's two knobs, the merge's two bounds, and the three closures; grouping \
              them would put a struct between the caller and its own closures for no reader's \
              benefit"
)]
pub fn merge_and_call_in_rounds<S, E, Sc, T, MK, CL>(
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
    E: Send + From<ReferenceUnreadable>,
    T: Send,
    Sc: Send,
    MK: Fn() -> Sc + Sync + Send,
    CL: Fn(&mut Sc, CohortObservation) -> T + Sync + Send,
{
    refuse_malformed_analysed_regions(analysed);

    let mut regions_in_round: Vec<GenomeRegion> = Vec::new();
    let mut graveyard: Vec<SampleLocusObservations> = Vec::new();
    let mut owned_through: Option<GenomePosition> = None;
    let mut displaced: u64 = 0;

    for analysed_region in analysed {
        let mut building_regions =
            building_regions_of(*analysed_region, cohort_locus_builder_regions_len);
        loop {
            regions_in_round.clear();
            regions_in_round.extend(building_regions.by_ref().take(regions_in_flight.get()));
            let Some(first) = regions_in_round.first().copied() else {
                break;
            };
            let last = regions_in_round
                .last()
                .copied()
                .expect("the round is not empty");

            let beside_the_builders = rayon::current_num_threads() > 1;
            let evicted_base = GenomePosition {
                contig: first.contig,
                position: first.start,
            };
            if beside_the_builders {
                cache.evict_before_in_parallel(evicted_base, &mut graveyard);
            } else {
                cache.evict_before(evicted_base);
            }
            cache.cover_in_parallel(GenomeRegion {
                contig: first.contig,
                start: first.start,
                end: last.end,
            })?;

            let cache = &*cache;
            let make_scratch = &make_scratch;
            let call_locus = &call_locus;
            let regions = &regions_in_round;
            let build_the_round = || -> Result<Vec<CalledRegion<T>>, E> {
                regions
                    .par_iter()
                    .map_init(
                        make_scratch,
                        |scratch, building_region| -> Result<CalledRegion<T>, E> {
                            let mut calls: Vec<(GenomeRegion, T)> = Vec::new();
                            let mut failed: Vec<GenomeRegion> = Vec::new();
                            // **The window's own prefix, added back to every member index.**
                            // The closer numbers members inside the window; `build_at` indexes
                            // the whole held list. A round evicts once, at its first region's
                            // first base, so every later region in the round has a prefix and
                            // the two numberings differ by it.
                            let mut starts: Vec<usize> = Vec::new();
                            cache.window_starts_into(*building_region, &mut starts);
                            let starts = &starts;
                            cache.with_observations(*building_region, |window| {
                                build_region_handing_over_windowed(
                                    *building_region,
                                    window,
                                    max_cohort_locus_span,
                                    min_alt_reads,
                                    &mut |built: CohortObservation| {
                                        let span = built.region;
                                        calls.push((span, call_locus(scratch, built)));
                                    },
                                    &mut failed,
                                    &|sample, index| cache.build_at(sample, starts[sample] + index),
                                )
                            })?;
                            Ok(CalledRegion { calls, failed })
                        },
                    )
                    .collect()
            };
            let outcomes: Result<Vec<CalledRegion<T>>, E> = if beside_the_builders {
                let mut dead = std::mem::take(&mut graveyard);
                let (outcomes, emptied) = rayon::join(build_the_round, move || {
                    dead.clear();
                    dead
                });
                graveyard = emptied;
                outcomes
            } else {
                build_the_round()
            };

            // **The frontier rule, on the calling thread, in region order** — spec §6.1's, the
            // organiser's `resolve_and_release` with the failed spans and the calls interleaved
            // by first base. Nothing a healthy run produces reaches the drop.
            for outcome in outcomes? {
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
            }
        }
    }

    Ok(displaced)
}

/// Claim `span`'s ground, or refuse it to whatever got there first.
fn claim(
    owned_through: &mut Option<GenomePosition>,
    displaced: &mut u64,
    span: GenomeRegion,
) -> bool {
    if owned_through.is_some_and(|owned| first_base_of(span) <= owned) {
        *displaced = displaced.saturating_add(1);
        return false;
    }
    *owned_through = Some(last_base_of(span));
    true
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{
        KeepingSource, in_flight, member, region, source_of, three_samples_over_six_hundred_bases,
        width,
    };
    use super::super::serial::merge_cohort_handing_each_locus_over;
    use super::super::{MaxCohortLocusSpan, MinAltReads};
    use super::*;
    use crate::ng::locus_generation::SampleLocusObservations;

    /// What one driver made of a fixture, in the shape both drivers can be compared in: each
    /// locus's span and its cohort observation rendered, and the spans the width bound refused.
    ///
    /// **The observations are rendered rather than compared**, for
    /// [`super::super::fixtures::render`]'s reason — `CohortObservation` has no `PartialEq`,
    /// and two sums that differ in their last bit render as different strings.
    #[derive(Debug, PartialEq, Eq)]
    struct WhatCameOut {
        loci: Vec<String>,
        refused: Vec<GenomeRegion>,
    }

    /// The streaming driver's answer over `per_sample`, which is what the round driver must
    /// reproduce.
    fn streaming(
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
            &mut |built| loci.push(format!("{built:?}")),
            &mut refused,
        )
        .expect("the fixture sources hold");
        WhatCameOut { loci, refused }
    }

    /// The round driver's answer over the same fixture, read through **keeping** sources — the
    /// psp path's shape, where the evidence is decoded inside the builder.
    fn in_rounds(
        analysed: &[GenomeRegion],
        per_sample: &[&[SampleLocusObservations]],
        region_width: CohortLocusBuilderRegionsLen,
        regions: usize,
    ) -> WhatCameOut {
        let mut cache = ObservationCache::over_fixture(
            per_sample
                .iter()
                .map(|sample| KeepingSource::over(sample))
                .collect(),
        );
        let mut loci = Vec::new();
        let mut refused = Vec::new();
        let displaced = merge_and_call_in_rounds(
            analysed,
            &mut cache,
            region_width,
            in_flight(regions),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
            || (),
            |(), built| format!("{built:?}"),
            &mut |rendered| loci.push(rendered),
            &mut refused,
        )
        .expect("the fixture sources hold");
        assert_eq!(displaced, 0, "a fixture merge displaced a locus");
        WhatCameOut { loci, refused }
    }

    /// **The claim of the whole arrangement**: a round of regions built and called on the pool
    /// produces what one region at a time produces, at every width and every round size — and
    /// **over sources that keep their evidence**, which is the shape a run over stored psps has
    /// and the shape no other comparison in this module uses.
    ///
    /// The widths straddle the fixture's 600 bases: 20 gives 30 regions, 610 gives one, and 305
    /// puts a boundary inside the 305–330 deletion, which is the case the plan named.
    #[test]
    fn a_round_of_regions_calls_what_one_region_at_a_time_calls() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 600)];

        for bases in [20, 100, 305, 610] {
            let expected = streaming(&analysed, &per_sample, width(bases));
            for regions in [2, 3, 7, 16] {
                let actual = in_rounds(&analysed, &per_sample, width(bases), regions);
                assert_eq!(
                    expected, actual,
                    "{regions} regions of {bases} bases in flight disagreed with one at a time",
                );
            }
        }
    }

    /// **The regression test for the defect the round driver was written on top of.**
    ///
    /// A closed locus numbers its members inside the window the cache handed the builder;
    /// `build_at` indexes the sample's whole held list. They differ by the prefix the window
    /// skipped, which is empty only for a round's *first* region — so a driver that forgot the
    /// prefix would build every later region's members from an earlier region's bodies. It is
    /// invisible to a source that hands its records over whole, because such a window carries
    /// the records themselves and `build_at` is never called.
    ///
    /// This fixture makes the prefix non-empty on purpose: one record per region over six
    /// regions, with a round holding all six, so five of them sit behind a prefix of one to
    /// five records. Deleting the `starts[sample] +` in the driver fails it.
    #[test]
    fn a_round_builds_each_regions_members_from_that_regions_own_bodies() {
        // One variant per hundred bases, so each 100-base region owns exactly one and the
        // window's prefix grows by one record a region.
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

        let expected = streaming(&analysed, &per_sample, width(100));
        assert_eq!(expected.loci.len(), 6, "the fixture holds six loci");
        let actual = in_rounds(&analysed, &per_sample, width(100), 6);
        assert_eq!(
            expected, actual,
            "a round of six regions built a locus from another region's evidence",
        );
    }

    /// **One region in flight is the streaming driver**, which is what lets the caller choose
    /// between them on that number alone.
    #[test]
    fn one_region_in_flight_is_one_region_at_a_time() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 600)];
        assert_eq!(
            streaming(&analysed, &per_sample, width(100)),
            in_rounds(&analysed, &per_sample, width(100), 1),
        );
    }

    /// **A round never crosses an analysed region**, so two intervals give at least two rounds
    /// and the answer is still one merge's.
    #[test]
    fn several_analysed_intervals_are_merged_as_one_run() {
        let samples = three_samples_over_six_hundred_bases();
        let per_sample: Vec<&[SampleLocusObservations]> =
            samples.iter().map(Vec::as_slice).collect();
        let analysed = [region(1, 300), region(301, 600)];
        assert_eq!(
            streaming(&analysed, &per_sample, width(50)),
            in_rounds(&analysed, &per_sample, width(50), 8),
        );
    }
}
