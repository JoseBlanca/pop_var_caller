//! **The joint parameters fit, timed without a CRAM file in the way.**
//!
//! The joint fit estimates error rates, heterozygosity, inbreeding, contamination and
//! repeat-tract slippage once over a whole cohort, before a single variant is called. Until
//! this file existed the only way to time it was `examples/ng_joint_records_walk`, which opens
//! CRAMs and a reference genome first — so a wall-clock delta could not tell a regression in
//! the per-tract likelihood loop from one in `noodles_cram::Block::decode`, and the
//! 2026-08-15 performance review named this the single highest-value missing measurement
//! (`doc/devel/reports/reviews/perf_ng-census-joint-fit_2026-08-15.md` §3, item 1).
//!
//! Every benchmark here calls the **production** function — `fit_stratum`, `fit_strata`,
//! `fit_jointly` — on evidence drawn in memory. No line of `src/` changed to make that
//! possible: all three were already `pub` on a `pub mod` path and take plain data.
//!
//! ```text
//! cargo bench --features bench-fixtures --bench ng_joint_fit_perf
//! ```
//!
//! Without `--features bench-fixtures` the target still compiles — so `cargo bench --no-run`
//! catches it when it stops building — and refuses at startup rather than reporting a time for
//! a fixture it could not draw.
//!
//! # What the fixtures are
//!
//! The drawn-evidence generators are the ones the two modules' own positive controls use,
//! moved out of `#[cfg(test)]` and behind the `bench-fixtures` feature:
//! `ssr_fit::bench_fixtures::draw_stratum` (a stratum drawn at a known slippage) and
//! `fit::bench_fixtures::draw_cohort` (a cohort drawn at a known error rate, heterozygosity
//! and inbreeding). **One generator, two callers** — a benchmark drawn differently from the
//! oracle would be timing a workload no test has ever checked.
//!
//! # The two axes, and why both are here
//!
//! `CLAUDE.md` §0 commits this caller to degrading gracefully from **one sample to several
//! thousand** and from **three reads a position to several hundred**. Both halves are
//! therefore swept over sample count as well as over their own size axis (tracts for the
//! repeat-tract half, kept positions for the ordinary-position half), because a change that
//! pays at eight samples and loses at one is a change this project has to be able to see.
//!
//! # Threads
//!
//! Both halves parallelise over the global rayon pool, so a number taken on an 18-core laptop
//! and one taken in the 8-core container are not comparable. Every timed body therefore runs
//! inside a pool of a **fixed** size — four threads, matching every other measurement in the
//! review — overridable with `NG_JOINT_FIT_BENCH_THREADS`. Sweeping that variable over
//! 1, 2, 4, 8, 18 against the same benchmark is the thread-count sweep §3 item 4 asks for.
//!
//! # Sizes: why these, and why they are small
//!
//! One whole run of this file is minutes. That is bought by cutting the *number of times* the
//! estimators evaluate their objective — one starting point instead of three, one round or
//! eight passes instead of five rounds or two hundred — while leaving every inner loop at
//! production's own shape: 256 quadrature points, thirteen allele classes, ninety-one
//! genotypes, sixteen frequency nodes. The report's own law says the repeat-tract fit is
//! **linear in tracts fitted** (0.157, 0.134, 0.151 and 0.133 seconds a tract over four runs
//! spanning 299 to 7,824 tracts), so a 32-tract stratum is a faithful and much cheaper proxy
//! for a 1,000-tract one — and the two tract sizes exist precisely so that a change can be
//! checked for having broken that linearity.

use criterion::{Criterion, criterion_group, criterion_main};

#[cfg(feature = "bench-fixtures")]
use criterion::BenchmarkId;
#[cfg(feature = "bench-fixtures")]
use std::hint::black_box;
#[cfg(feature = "bench-fixtures")]
use std::time::Duration;

#[cfg(feature = "bench-fixtures")]
use pop_var_caller::ng::parameter_estimation::joint::census::Stratum;
#[cfg(feature = "bench-fixtures")]
use pop_var_caller::ng::parameter_estimation::joint::fit::bench_fixtures::{
    DrawnCohort, as_cohort, draw_cohort,
};
#[cfg(feature = "bench-fixtures")]
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFitConfig, StartingPoint as GenericStart, fit_jointly,
};
#[cfg(feature = "bench-fixtures")]
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit::bench_fixtures::{
    draw_stratum, spectrum_of,
};
#[cfg(feature = "bench-fixtures")]
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit::{
    Slippage, SsrFitConfig, StartingPoint as SsrStart, StratumEvidence, StratumOutcome, fit_strata,
    fit_stratum,
};

// ---------------------------------------------------------------------------
// The pool every timed body runs in
// ---------------------------------------------------------------------------

/// Threads the timed bodies get when `NG_JOINT_FIT_BENCH_THREADS` says nothing.
///
/// **Four, because that is what every wall number in the 2026-08-15 review was taken at**, so
/// a bench number and the review's harness numbers describe the same machine width.
#[cfg(feature = "bench-fixtures")]
const DEFAULT_THREADS: usize = 4;

/// A pool of a fixed width, built once per benchmark group.
///
/// **Not the global pool.** `fit_stratum` and `fit_jointly` both parallelise over
/// `rayon::current_num_threads()`, which is the machine's core count unless something says
/// otherwise — so the same commit measured on the 18-core macOS host and in the container
/// would produce two numbers that differ for a reason that is not the code.
#[cfg(feature = "bench-fixtures")]
fn fixed_pool() -> rayon::ThreadPool {
    let threads = std::env::var("NG_JOINT_FIT_BENCH_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|threads| *threads > 0)
        .unwrap_or(DEFAULT_THREADS);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a rayon pool of the requested width")
}

// ---------------------------------------------------------------------------
// The repeat-tract half
// ---------------------------------------------------------------------------

/// The truth every drawn stratum here carries: tomato's own dinucleotide numbers — a read
/// slips about eight times in a hundred, and a slipped read shows a shorter tract five times
/// as often as a longer one.
#[cfg(feature = "bench-fixtures")]
const SLIPPAGE: Slippage = Slippage {
    level: 0.08,
    shorter_share: 0.83,
    fall_off: 0.25,
};

/// How far either side of the reference length the draw spreads allele mass, **and** the read
/// span it records in.
///
/// Six, which is `ALLELE_SPAN` — the production value — so a drawn tract is fitted over
/// thirteen allele classes and the ninety-one genotypes of a diploid over them. That number is
/// the trip count of the innermost loop the profile puts 35.2% of busy CPU in, so it is the
/// one thing here that must not be shrunk to make the bench cheap.
#[cfg(feature = "bench-fixtures")]
const SPAN: i32 = 6;

/// Reads a sample puts on a tract. Six is the middle of the range this caller is for: the
/// tomato cohort sits near three and the GIAB trio near thirty.
#[cfg(feature = "bench-fixtures")]
const DEPTH: u32 = 6;

/// How monomorphic a drawn stratum's tracts are, and how far each drawn sample falls short of
/// the heterozygote proportions. Both are the positive control's values.
#[cfg(feature = "bench-fixtures")]
const CONCENTRATION: f64 = 0.5;
#[cfg(feature = "bench-fixtures")]
const HOM_EXCESS: f64 = 0.4;

/// The repeat-tract fit as the bench asks for it: **one starting point and one round**.
///
/// Production climbs from three starts for up to five rounds, and each round is a few hundred
/// evaluations of the same objective. Cutting to one and one divides the wall time by roughly
/// fifteen and changes nothing about the loop being measured — the objective, its quadrature,
/// its allele classes and its genotypes are all production's.
#[cfg(feature = "bench-fixtures")]
fn ssr_config(borrowing_floor: usize) -> SsrFitConfig {
    SsrFitConfig {
        starting_points: vec![SsrStart {
            slippage_level: 0.10,
            concentration: 3.0,
        }],
        max_rounds: 1,
        borrowing_floor,
        // **Four and not the production fifty**, because these strata are deliberately far
        // thinner than a real one: at the default every drawn stratum here would come back
        // `Refused` and the benchmark would time the refusal.
        refusal_floor: 4,
        ..SsrFitConfig::default()
    }
}

/// One drawn stratum, at a given repeat count so a set of them can borrow from each other.
#[cfg(feature = "bench-fixtures")]
fn stratum(tracts: usize, samples: usize, reference_repeats: u64, seed: u64) -> StratumEvidence {
    let spectrum = spectrum_of((2 * SPAN + 1) as usize);
    let mut evidence = draw_stratum(
        SLIPPAGE,
        &spectrum,
        CONCENTRATION,
        HOM_EXCESS,
        tracts,
        samples,
        DEPTH,
        SPAN,
        seed,
    );
    evidence.stratum = Stratum {
        period: 2,
        reference_repeats,
    };
    evidence
}

/// **The guard against timing a fit that fitted nothing.**
///
/// A `fit_stratum` that refuses — because no group has reads, or because the draw produced no
/// tracts — returns `None` in microseconds, and criterion would report that as a very fast
/// fit. Three things are pinned, each of which can fail on its own: the fit happened, it read
/// every tract it was given, and it returned a number rather than an infinity.
#[cfg(feature = "bench-fixtures")]
fn check_the_stratum_was_fitted(evidence: &StratumEvidence, excess: &[f64], config: &SsrFitConfig) {
    let fitted = fit_stratum(evidence, excess, config)
        .expect("a drawn stratum has reads in its one slippage group");
    assert_eq!(
        fitted.tracts_fitted,
        evidence.tracts.len(),
        "the fit must read every drawn tract, not a prefix of them"
    );
    assert!(
        fitted.log_likelihood_a_tract.is_finite(),
        "the mean log-likelihood a tract came back as {}, so the timed body is not doing the \
         arithmetic this bench claims",
        fitted.log_likelihood_a_tract
    );
    assert!(
        fitted.slippage[0].is_some(),
        "the one slippage group put reads in this stratum and must be fitted"
    );
}

/// **Tract count at a fixed eight samples.** The axis the review's own law is stated on, and
/// the one a change to the per-tract likelihood loop moves.
#[cfg(feature = "bench-fixtures")]
fn stratum_by_tracts(c: &mut Criterion) {
    const SAMPLES: usize = 8;
    let pool = fixed_pool();
    let excess = vec![HOM_EXCESS; SAMPLES];
    let config = ssr_config(1);

    let mut group = c.benchmark_group("ng_joint_fit/stratum_by_tracts");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for &tracts in &[32_usize, 128] {
        let evidence = stratum(tracts, SAMPLES, 10, 0x5EED_0001);
        check_the_stratum_was_fitted(&evidence, &excess, &config);
        group.bench_with_input(BenchmarkId::from_parameter(tracts), &tracts, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    black_box(fit_stratum(
                        black_box(&evidence),
                        black_box(&excess),
                        black_box(&config),
                    ))
                })
            });
        });
    }
    group.finish();
}

/// **Sample count at a fixed thirty-two tracts** — the first of the two axes `CLAUDE.md` §0
/// commits to.
///
/// One sample is included because it is the case a cohort method is most likely to get wrong,
/// and because the per-tract cost the review measured (`0.043 + 0.014 × samples` seconds) has
/// a constant term worth a third of the cost at eight samples: whether that constant is still
/// there after a change is a question only this axis can answer.
#[cfg(feature = "bench-fixtures")]
fn stratum_by_samples(c: &mut Criterion) {
    const TRACTS: usize = 32;
    let pool = fixed_pool();
    let config = ssr_config(1);

    let mut group = c.benchmark_group("ng_joint_fit/stratum_by_samples");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for &samples in &[1_usize, 8, 32] {
        let evidence = stratum(TRACTS, samples, 10, 0x5EED_0002);
        let excess = vec![HOM_EXCESS; samples];
        check_the_stratum_was_fitted(&evidence, &excess, &config);
        group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    black_box(fit_stratum(
                        black_box(&evidence),
                        black_box(&excess),
                        black_box(&config),
                    ))
                })
            });
        });
    }
    group.finish();
}

/// How many strata the borrowing benchmark holds, and how many tracts each holds.
///
/// Six thin strata of sixteen tracts. **The product is the same on both sides of the
/// comparison below** — ninety-six tracts fitted either way — which is what makes the pair
/// readable.
#[cfg(feature = "bench-fixtures")]
const STRATA: usize = 6;
#[cfg(feature = "bench-fixtures")]
const TRACTS_A_STRATUM: usize = 16;

/// **The borrowing and dedup path, as a pair of runs that fit the same number of tracts.**
///
/// `fit_strata` lets a stratum too thin to carry its own answer take in its neighbouring
/// repeat counts, and then fits each *distinct pooled set* once. The two cases here differ
/// only in the floor:
///
/// - `stands_alone` — the floor is 1, so each of the six strata is fitted on its own sixteen
///   tracts: **six fits of sixteen tracts**.
/// - `borrows_and_shares` — the floor is unreachable, so every stratum takes in all five
///   neighbours, all six arrive at the identical pooled set, and the dedup collapses them to
///   **one fit of ninety-six tracts** plus five clones of its answer.
///
/// Both fit ninety-six tracts' worth of likelihood. What separates them is everything else:
/// six quadrature rebuilds against one — the review puts the 256-point rebuild at 21.7% of
/// the repeat-tract fit, and it costs the same whether a stratum holds sixteen tracts or a
/// thousand — against six clones of the pooled evidence. **On tomato 68 of 141 strata borrow**,
/// so which way that trade falls is not an edge case.
#[cfg(feature = "bench-fixtures")]
fn strata_borrowing(c: &mut Criterion) {
    const SAMPLES: usize = 8;
    let pool = fixed_pool();
    let excess = vec![HOM_EXCESS; SAMPLES];

    // Six neighbouring repeat counts of one motif length, so every stratum can borrow from
    // every other. Same period throughout: `fit_strata` never pools across motif lengths.
    let strata: Vec<StratumEvidence> = (0..STRATA)
        .map(|which| {
            stratum(
                TRACTS_A_STRATUM,
                SAMPLES,
                8 + which as u64,
                0x5EED_0100 + which as u64,
            )
        })
        .collect();

    let mut group = c.benchmark_group("ng_joint_fit/strata");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for (name, floor, borrowed_each, tracts_each) in [
        ("stands_alone", 1_usize, 0_usize, TRACTS_A_STRATUM),
        (
            "borrows_and_shares",
            usize::MAX,
            STRATA - 1,
            STRATA * TRACTS_A_STRATUM,
        ),
    ] {
        let config = ssr_config(floor);
        check_the_strata_took_the_route_they_were_named_for(
            &strata,
            &excess,
            &config,
            borrowed_each,
            tracts_each,
        );
        group.bench_function(name, |b| {
            b.iter(|| {
                pool.install(|| {
                    black_box(fit_strata(
                        black_box(&strata),
                        black_box(&excess),
                        black_box(&config),
                    ))
                })
            });
        });
    }
    group.finish();
}

/// **The guard that says which route `fit_strata` actually took.**
///
/// The floor is the only difference between the two cases, and a floor that stopped working
/// would leave both timing the same route while both still returned a plausible number.
/// Nothing else in the bench would notice. So each case asserts the shape it is named for:
/// how many neighbours each stratum borrowed, and how many tracts its answer was fitted from.
#[cfg(feature = "bench-fixtures")]
fn check_the_strata_took_the_route_they_were_named_for(
    strata: &[StratumEvidence],
    excess: &[f64],
    config: &SsrFitConfig,
    borrowed_each: usize,
    tracts_each: usize,
) {
    let outcomes = fit_strata(strata, excess, config);
    assert_eq!(outcomes.len(), STRATA);
    for outcome in &outcomes {
        match outcome {
            StratumOutcome::Fitted(fit) => {
                assert_eq!(
                    fit.borrowed.len(),
                    borrowed_each,
                    "a stratum borrowed from {} neighbours where this case is the one that \
                     borrows from {borrowed_each}",
                    fit.borrowed.len()
                );
                assert_eq!(
                    fit.tracts_fitted, tracts_each,
                    "the answer was fitted from {} tracts, not the {tracts_each} this case names",
                    fit.tracts_fitted
                );
            }
            StratumOutcome::Refused { reason, .. } => {
                panic!("a drawn stratum was refused ({reason:?}), so nothing was timed")
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The ordinary-position half
// ---------------------------------------------------------------------------

/// The population the ordinary positions are drawn from: nine positions in ten invariant, one
/// in a hundred fixed for the non-reference base, the rest segregating at a rare-allele-heavy
/// density. These are the oracle test's own numbers.
#[cfg(feature = "bench-fixtures")]
fn density() -> FrequencyDensity {
    FrequencyDensity {
        p_invariant: 0.90,
        p_fixed_alt: 0.01,
        a: 0.7,
        b: 2.5,
    }
}

/// Reads a sample puts on an ordinary position. Eight — between the tomato cohort's three and
/// the GIAB trio's thirty.
#[cfg(feature = "bench-fixtures")]
const MEAN_DEPTH: f64 = 8.0;

/// The chemistry the cohort is drawn at: a clean position's reads disagree 2 times in 1,000,
/// a mismapped one's 6 in 100, and 2 positions in 100 are mismapped.
#[cfg(feature = "bench-fixtures")]
const NOISE: (f64, f64, f64) = (0.002, 0.06, 0.02);

/// How many passes of the alternation the timed fit runs.
///
/// **Eight, and it is the whole reason one iteration is a second rather than a minute.**
/// Production allows two hundred and stops when nothing moves; the review's harness ran at
/// sixty. Every pass is one sweep of `one_position` over every kept position, so the pass
/// count multiplies the kernel this benchmark exists to watch and divides nothing else —
/// which is also the one caveat to read the number with: `fit_jointly` runs its contamination
/// step **once**, after the alternation, so at eight passes that step holds a larger share of
/// this benchmark than it holds of a production run.
#[cfg(feature = "bench-fixtures")]
const PASSES: u32 = 8;

/// The ordinary-position fit as the bench asks for it: one starting point, eight passes, and
/// **everything else production's own** — sixteen frequency nodes, the duplicated-position
/// class on, the depth code read as the range it stands for.
#[cfg(feature = "bench-fixtures")]
fn generic_config() -> JointFitConfig {
    JointFitConfig {
        starting_points: vec![GenericStart {
            clean: 0.002,
            noisy: 0.05,
            noisy_share: 0.01,
            p_invariant: 0.97,
            p_fixed_alt: 0.005,
            a: 0.5,
            b: 2.0,
            duplicated_share: 0.001,
            carrier_a: 1.2,
            carrier_b: 9.5,
        }],
        max_passes: PASSES,
        ..JointFitConfig::default()
    }
}

/// **The guard against timing a fit that walked no positions.**
///
/// `noisy_posterior` carries one value per kept position, in position order, so its length is
/// direct evidence that the timed body swept the whole census rather than an empty one — the
/// failure a wall number cannot show, because a fit over nothing is fast and finite.
#[cfg(feature = "bench-fixtures")]
fn check_the_fit_swept_every_position(drawn: &DrawnCohort, positions: usize, samples: usize) {
    let mut cohort = as_cohort(&drawn.samples);
    let fit = fit_jointly(&mut cohort, &generic_config()).expect("a drawn cohort pools");
    assert_eq!(
        fit.noisy_posterior.len(),
        positions,
        "the fit produced {} per-position posteriors over a census of {positions} positions",
        fit.noisy_posterior.len()
    );
    assert_eq!(fit.rates.len(), samples);
    assert!(
        fit.log_likelihood.is_finite(),
        "the fit's log-likelihood came back as {}",
        fit.log_likelihood
    );
    assert!(
        fit.rates["s0"].value.positions_with_reads > (positions as u64) / 2,
        "the first sample carried reads at {} of {positions} positions, so the draw is not the \
         {MEAN_DEPTH}× one this bench names",
        fit.rates["s0"].value.positions_with_reads
    );
}

/// **Kept positions at a fixed eight samples**, through `fit_jointly` and therefore through
/// `one_position` — the kernel the review's profile puts 30.4% of busy CPU in, and the largest
/// single thing left after round one's ten fixes.
///
/// The cost is expected to be linear in positions; the pair is here so that a change which
/// breaks that — by turning a per-position constant into a per-position search, say — shows up
/// as a bend rather than as a number nobody can interpret.
#[cfg(feature = "bench-fixtures")]
fn generic_by_positions(c: &mut Criterion) {
    const SAMPLES: usize = 8;
    let pool = fixed_pool();
    let config = generic_config();

    let mut group = c.benchmark_group("ng_joint_fit/generic_by_positions");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for &positions in &[5_000_usize, 20_000] {
        let drawn = draw_cohort(
            SAMPLES,
            positions,
            MEAN_DEPTH,
            NOISE,
            density(),
            0.2,
            0x9E37_79B9_7F4A_7C15,
        );
        check_the_fit_swept_every_position(&drawn, positions, SAMPLES);
        // Built once and re-fitted: a resident census lends its sections and takes them back,
        // so nothing is consumed, and rebuilding it per iteration would put a deep clone of
        // every sample's records inside the measurement.
        let mut cohort = as_cohort(&drawn.samples);
        group.bench_with_input(
            BenchmarkId::from_parameter(positions),
            &positions,
            |b, _| {
                b.iter(|| {
                    pool.install(|| {
                        black_box(fit_jointly(black_box(&mut cohort), black_box(&config)))
                            .expect("a drawn cohort pools")
                    })
                });
            },
        );
    }
    group.finish();
}

/// **Sample count at a fixed five thousand positions** — the second half of the axis
/// `CLAUDE.md` §0 commits to, on the half of the fit that holds every sample at once.
///
/// One sample is the case where the cohort-wide quantities have no cohort to be read from, and
/// thirty-two is where the per-sample arrays start to be the thing being swept rather than a
/// detail. The three points hold `samples × positions` roughly constant against the largest
/// point of the sweep above, so the two groups can be read together.
#[cfg(feature = "bench-fixtures")]
fn generic_by_samples(c: &mut Criterion) {
    const POSITIONS: usize = 5_000;
    let pool = fixed_pool();
    let config = generic_config();

    let mut group = c.benchmark_group("ng_joint_fit/generic_by_samples");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for &samples in &[1_usize, 8, 32] {
        let drawn = draw_cohort(
            samples,
            POSITIONS,
            MEAN_DEPTH,
            NOISE,
            density(),
            0.2,
            0x9E37_79B9_7F4A_7C16,
        );
        check_the_fit_swept_every_position(&drawn, POSITIONS, samples);
        let mut cohort = as_cohort(&drawn.samples);
        group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    black_box(fit_jointly(black_box(&mut cohort), black_box(&config)))
                        .expect("a drawn cohort pools")
                })
            });
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Entry
// ---------------------------------------------------------------------------

/// Every benchmark in this file, behind one target so the feature gate is stated once.
fn every_bench(c: &mut Criterion) {
    #[cfg(not(feature = "bench-fixtures"))]
    {
        let _ = c;
        panic!(
            "ng_joint_fit_perf draws its evidence from the joint fit's own fixture generators, \
             which are compiled only behind the `bench-fixtures` feature. Run it as:\n\n    \
             cargo bench --features bench-fixtures --bench ng_joint_fit_perf\n"
        );
    }
    #[cfg(feature = "bench-fixtures")]
    {
        stratum_by_tracts(c);
        stratum_by_samples(c);
        strata_borrowing(c);
        generic_by_positions(c);
        generic_by_samples(c);
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = every_bench
}

criterion_main!(benches);
