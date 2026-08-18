//! Where the cohort merge spends its time, at cohort sizes and record densities spanning the
//! caller's range.
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
//! **Two densities, because the answer differs between them.** *Dense* is 2,000 bases with a
//! record every four, so every 20-base building region holds five records per sample: the
//! ground the first measurements were taken on. *Sparse* is 20,000 bases with a record every
//! hundred — roughly the measured tomato corner, where about one position in a hundred varies —
//! so four building regions in five hold no record for any sample at all. Sparse ground is
//! where the cache's fixed per-region cost shows: the builders have nothing to do and the walk
//! over the cohort happens anyway.
//!
//! **The fixture is built outside the clock.** The cache owns the observations it is handed, so
//! every repeat needs its own copy of the cohort; at 3,000 samples that copy is a million
//! records and several hundred milliseconds. Cloning it inside the timed span — which this
//! probe did until 2026-08-18 — charged the cached drivers for work the oracle, which borrows
//! its slices, never did, and made the cached path look between two and three times worse than
//! it is.
//!
//! **It sweeps the building region's width as well as the cohort's size**, because that knob
//! turned out to matter more than anything else measured: the merge pays a fixed cost per
//! building region per sample — a cover, an eviction, a window and a builder's setup, each a
//! walk over the whole cohort — and on sparse ground a 20-base region holds a fifth of a
//! locus. What widening costs is memory, which is the whole reason the regions are short
//! (`doc/devel/ng/spec/cohort_merge.md` §8): the cache holds `regions in flight × width`
//! bases of ground plus the tail reaching past it.
//!
//! **It reports the median of five repeats and the spread**, because this machine's own swing
//! at 3,000 samples reached 30% between runs of one unchanged binary. Comparing two builds
//! means running them alternately in one sitting: the machine drifted 14% across two runs of
//! one unchanged binary an hour apart.
//!
//! **It also sweeps the thread count**, in its own pass: the parallel driver is run inside a
//! rayon pool of 1, 2, 3, 4 and 8 threads with one region in flight per thread — the rule a run
//! with no `--cohort-locus-builder-regions-in-flight` takes — against the one-thread driver
//! that uses no pool at all. That last one is the honest baseline for *what threads buy*: the
//! parallel driver at one thread still pays rayon for every round.
//!
//! Run in release: `./scripts/dev.sh cargo run --release --example
//! ng_cohort_merge_parallel_cost`. `NG_MERGE_COST_COHORTS=1,63`, `NG_MERGE_COST_WIDTHS=20`,
//! `NG_MERGE_COST_GROUND=sparse` and `NG_MERGE_COST_THREADS=1,4` narrow it while iterating —
//! **through `env` inside the container**, as in `./scripts/dev.sh env
//! NG_MERGE_COST_COHORTS=63 cargo run …`, because `scripts/dev.sh` forwards two variables of
//! its own and no others. Set on the host side of it they are silently ignored and the whole
//! matrix runs.

use std::num::{NonZeroU32, NonZeroUsize};
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

/// How many repeats each number is the median of.
const REPEATS: usize = 5;

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

/// One stretch of fabricated ground: how long it is, and how far apart one sample's records sit.
struct Ground {
    name: &'static str,
    last_base: u64,
    bases_between_records: u64,
}

impl Ground {
    /// Every sample's records over this ground, in coordinate order.
    ///
    /// One record in eight is off-reference, and which one moves with the sample, so no
    /// building region is empty of work for the whole cohort at once.
    fn cohort(&self, samples: usize) -> Vec<Vec<SampleLocusObservations>> {
        (0..samples)
            .map(|sample| {
                (0..self.last_base / self.bases_between_records)
                    .map(|record| {
                        let at = record * self.bases_between_records + 1;
                        let observed: &[u8] = if (record as usize + sample) % 8 == 0 {
                            b"C"
                        } else {
                            b"A"
                        };
                        record_at(at, observed)
                    })
                    .collect()
            })
            .collect()
    }

    fn analysed(&self) -> [GenomeRegion; 1] {
        [GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(self.last_base),
        }]
    }
}

/// The median, fastest and slowest of [`REPEATS`] runs, in milliseconds.
///
/// **`prepare` runs outside the clock and `one_merge` inside it.** Every repeat gets its own
/// copy of whatever the merge consumes, and building that copy is fixture work: charging it to
/// the driver was this probe's own longest-standing defect.
fn timed<T>(mut prepare: impl FnMut() -> T, mut one_merge: impl FnMut(T)) -> (f64, f64, f64) {
    one_merge(prepare());
    let mut each_repeat: Vec<f64> = (0..REPEATS)
        .map(|_| {
            let prepared = prepare();
            let started = Instant::now();
            one_merge(prepared);
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

/// One sample's source: its records, owned, in coordinate order.
type SampleSource = std::vec::IntoIter<Result<SampleLocusObservations, Never>>;

/// A fresh set of sources over `cohort` — one copy of every record, which the cache then owns.
fn sources_over(cohort: &[Vec<SampleLocusObservations>]) -> Vec<SampleSource> {
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
}

/// What a pool of `threads` buys, against the driver that uses no pool at all.
///
/// **One region in flight per thread**, which is what a run takes when the operator names no
/// count ([`CohortLocusBuilderRegionsInFlight::one_per_worker_thread`]) — so the ground the
/// cache holds grows with the thread count, and the two are swept together rather than one at
/// a time.
fn threads_sweep(
    ground: &Ground,
    cohort: &[Vec<SampleLocusObservations>],
    samples: usize,
    width: CohortLocusBuilderRegionsLen,
    thread_counts: &[usize],
) {
    let analysed = ground.analysed();
    let name = ground.name;
    let bases = width.get();

    let (median, fastest, slowest) = timed(
        || ObservationCache::over(sources_over(cohort)),
        |mut cache| {
            let outcome = merge_cohort_through_cache(
                &analysed,
                &mut cache,
                width,
                MaxCohortLocusSpan::DEFAULT,
                MinAltObs::DEFAULT,
            )
            .expect("the probe's sources cannot fail");
            std::hint::black_box(&outcome);
        },
    );
    println!(
        "{name}, {samples}, one thread no pool, {bases}, {median:.2}, {fastest:.2}, {slowest:.2}"
    );

    for &threads in thread_counts {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool of the asked-for size");
        let regions =
            CohortLocusBuilderRegionsInFlight(NonZeroUsize::new(threads).expect("non-zero"));
        let (median, fastest, slowest) = pool.install(|| {
            timed(
                || ObservationCache::over(sources_over(cohort)),
                |mut cache| {
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
                },
            )
        });
        println!(
            "{name}, {samples}, pool of {threads}, {bases}, {median:.2}, {fastest:.2}, \
             {slowest:.2}"
        );
    }
}

fn main() {
    let grounds = [
        Ground {
            name: "dense",
            last_base: 2_000,
            bases_between_records: 4,
        },
        Ground {
            name: "sparse",
            last_base: 20_000,
            bases_between_records: 100,
        },
    ];
    let wanted_ground = std::env::var("NG_MERGE_COST_GROUND").ok();
    let cohorts: Vec<usize> = match std::env::var("NG_MERGE_COST_COHORTS") {
        Ok(list) => list
            .split(',')
            .map(|count| count.trim().parse().expect("a cohort size"))
            .collect(),
        Err(_) => vec![1, 63, 250, 1_000, 3_000],
    };
    let widths: Vec<CohortLocusBuilderRegionsLen> = match std::env::var("NG_MERGE_COST_WIDTHS") {
        Ok(list) => list
            .split(',')
            .map(|bases| {
                CohortLocusBuilderRegionsLen(
                    NonZeroU32::new(bases.trim().parse().expect("a width in bases"))
                        .expect("a width of at least one base"),
                )
            })
            .collect(),
        Err(_) => vec![
            CohortLocusBuilderRegionsLen(NonZeroU32::new(20).expect("non-zero")),
            CohortLocusBuilderRegionsLen::DEFAULT,
        ],
    };

    let thread_counts: Vec<usize> = match std::env::var("NG_MERGE_COST_THREADS") {
        Ok(list) => list
            .split(',')
            .map(|count| count.trim().parse().expect("a thread count"))
            .collect(),
        Err(_) => vec![1, 2, 3, 4, 8],
    };

    println!("threads available: {}", rayon::current_num_threads());
    println!("ground, cohort, driver, region_bases, median_ms, min_ms, max_ms");

    for ground in &grounds {
        if wanted_ground
            .as_deref()
            .is_some_and(|wanted| wanted != ground.name)
        {
            continue;
        }
        let analysed = ground.analysed();

        for &samples in &cohorts {
            let cohort = ground.cohort(samples);
            let name = ground.name;

            let slices: Vec<&[SampleLocusObservations]> =
                cohort.iter().map(Vec::as_slice).collect();
            let (median, fastest, slowest) = timed(
                || (),
                |()| {
                    let outcome = merge_cohort_serially(
                        &analysed,
                        &slices,
                        MaxCohortLocusSpan::DEFAULT,
                        MinAltObs::DEFAULT,
                    );
                    std::hint::black_box(&outcome);
                },
            );
            println!("{name}, {samples}, oracle, -, {median:.2}, {fastest:.2}, {slowest:.2}");

            for &width in &widths {
                let bases = width.get();
                let (median, fastest, slowest) = timed(
                    || ObservationCache::over(sources_over(&cohort)),
                    |mut cache| {
                        let outcome = merge_cohort_through_cache(
                            &analysed,
                            &mut cache,
                            width,
                            MaxCohortLocusSpan::DEFAULT,
                            MinAltObs::DEFAULT,
                        )
                        .expect("the probe's sources cannot fail");
                        std::hint::black_box(&outcome);
                    },
                );
                println!(
                    "{name}, {samples}, cached serial, {bases}, {median:.2}, {fastest:.2}, \
                     {slowest:.2}"
                );

                for in_flight in [8usize, 16] {
                    let regions = CohortLocusBuilderRegionsInFlight(
                        NonZeroUsize::new(in_flight).expect("non-zero"),
                    );
                    let (median, fastest, slowest) = timed(
                        || ObservationCache::over(sources_over(&cohort)),
                        |mut cache| {
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
                        },
                    );
                    println!(
                        "{name}, {samples}, parallel x{in_flight}, {bases}, {median:.2}, \
                         {fastest:.2}, {slowest:.2}"
                    );
                }
            }

            threads_sweep(
                ground,
                &cohort,
                samples,
                CohortLocusBuilderRegionsLen::DEFAULT,
                &thread_counts,
            );
        }
    }
}
