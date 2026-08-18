//! What the parallel merge buys over the serial one, at cohort sizes spanning the caller's
//! range.
//!
//! **A one-off probe, not a benchmark.** Nothing in CI runs it and nothing gates on it; it
//! exists because the parallel driver's correctness was proved four ways and its speed had
//! never been measured at all.
//!
//! It times the same merge three ways over the same fabricated ground — the oracle
//! (`merge_cohort_serially`, holding every sample's observations at once), the cached serial
//! driver (`merge_cohort_through_cache`, one forward reader per sample), and the parallel one
//! (`merge_cohort_in_parallel`) at several counts of regions in flight. All three produce the
//! same answer, which the test suite already asserts; what differs is time.
//!
//! **The ground is 2,000 bases with a record every four**, so every 20-base building region
//! holds five records per sample and there is real work in each. That is denser than the
//! measured tomato corner, where about one position in a hundred varies — a merge over sparse
//! ground spends proportionally more of its time on the cache and less inside the builders,
//! so the speed-up here is the optimistic end of the range, not the typical one.
//!
//! **It reports the median of seven repeats and the spread**, because this machine's own swing
//! at 3,000 samples reached 30% between runs of one unchanged binary.
//!
//! **What it cannot measure from here is the cover's own share**, which is the thing the
//! numbers turn on: `ObservationCache::cover` and `evict_before` are `pub(super)`, so only the
//! module's own drivers may drive the cache. Measured once with those two temporarily widened,
//! at 1,000 samples over this ground: evicting and covering every building region and building
//! nothing at all took **101.0 ms**, against **111.5 ms** for the whole cached serial merge —
//! so the builders are 9% of it and the cache is the rest. Evicting once per round of sixteen
//! instead of once per region took the same cover-only loop to **79.9 ms**.
//!
//! Run in release: `cargo run --release --example ng_cohort_merge_parallel_cost`.

use std::num::NonZeroUsize;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use pop_var_caller::ng::run::cohort_merge::observation_cache::ObservationCache;
use pop_var_caller::ng::run::cohort_merge::parallel::merge_cohort_in_parallel;
use pop_var_caller::ng::run::cohort_merge::serial::{
    merge_cohort_serially, merge_cohort_through_cache,
};
use pop_var_caller::ng::run::cohort_merge::{
    CohortLocusBuilderRegionsInFlight, CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltObs,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};

/// This probe cannot fail, so its source never errors.
#[derive(Debug)]
struct Never;

/// One sample's record at one position, showing `observed` to three reads.
fn record_at(position: u64, observed: &[u8]) -> SampleLocusObservations {
    SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(0),
            start: Position(position),
            end: Position(position),
        },
        reference_bases: Box::from(&b"A"[..]),
        observations: vec![SequenceObservation {
            bases: Box::from(observed),
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

/// The median, fastest and slowest of seven repeats, in milliseconds.
fn timed(mut one_merge: impl FnMut()) -> (f64, f64, f64) {
    one_merge();
    let mut each_repeat: Vec<f64> = (0..7)
        .map(|_| {
            let started = Instant::now();
            one_merge();
            started.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    each_repeat.sort_by(f64::total_cmp);
    (
        each_repeat[each_repeat.len() / 2],
        each_repeat[0],
        each_repeat[each_repeat.len() - 1],
    )
}

fn main() {
    let ground_end = 2_000u64;
    let analysed = [GenomeRegion {
        contig: ContigId(0),
        start: Position(1),
        end: Position(ground_end),
    }];
    let width = CohortLocusBuilderRegionsLen::DEFAULT;

    println!("threads available: {}", rayon::current_num_threads());
    println!("cohort, driver, median_ms, min_ms, max_ms");

    for samples in [1usize, 10, 63, 250, 1000] {
        // Every sample carries a record every four bases, one of them off-reference, so every
        // building region has work in it for every sample.
        let cohort: Vec<Vec<SampleLocusObservations>> = (0..samples)
            .map(|sample| {
                (0..ground_end / 4)
                    .map(|record| {
                        let at = record * 4 + 1;
                        let observed: &[u8] = if (record as usize + sample) % 8 == 0 {
                            b"C"
                        } else {
                            b"A"
                        };
                        record_at(at, observed)
                    })
                    .collect()
            })
            .collect();
        let sources = || -> Vec<std::vec::IntoIter<Result<SampleLocusObservations, Never>>> {
            cohort
                .iter()
                .map(|sample| {
                    sample
                        .iter()
                        .cloned()
                        .map(Ok)
                        .collect::<Vec<_>>()
                        .into_iter()
                })
                .collect()
        };

        let slices: Vec<&[SampleLocusObservations]> = cohort.iter().map(Vec::as_slice).collect();
        let (median, fastest, slowest) = timed(|| {
            let outcome = merge_cohort_serially(
                &analysed,
                &slices,
                MaxCohortLocusSpan::DEFAULT,
                MinAltObs::DEFAULT,
            );
            std::hint::black_box(&outcome);
        });
        println!("{samples}, oracle, {median:.2}, {fastest:.2}, {slowest:.2}");

        let (median, fastest, slowest) = timed(|| {
            let mut cache = ObservationCache::over(sources());
            let outcome = merge_cohort_through_cache(
                &analysed,
                &mut cache,
                width,
                MaxCohortLocusSpan::DEFAULT,
                MinAltObs::DEFAULT,
            )
            .expect("the probe's sources cannot fail");
            std::hint::black_box(&outcome);
        });
        println!("{samples}, cached serial, {median:.2}, {fastest:.2}, {slowest:.2}");

        for in_flight in [1usize, 2, 4, 8, 16] {
            let regions =
                CohortLocusBuilderRegionsInFlight(NonZeroUsize::new(in_flight).expect("non-zero"));
            let (median, fastest, slowest) = timed(|| {
                let mut cache = ObservationCache::over(sources());
                let outcome = merge_cohort_in_parallel(
                    &analysed,
                    &mut cache,
                    width,
                    regions,
                    MaxCohortLocusSpan::DEFAULT,
                    MinAltObs::DEFAULT,
                )
                .expect("the probe's sources cannot fail");
                std::hint::black_box(&outcome);
            });
            println!("{samples}, parallel x{in_flight}, {median:.2}, {fastest:.2}, {slowest:.2}");
        }
    }
}
