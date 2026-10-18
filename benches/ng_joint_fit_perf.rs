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
//! Every benchmark here calls the **production** function, `fit_jointly`, on evidence drawn in
//! memory. **Only the ordinary-position half is timed now**: the repeat-tract benchmarks timed
//! the climb (`fit_stratum`, `fit_strata`), which was deleted on 2026-10-16
//! (`doc/devel/ng/arch/str_penalised_fit.md` §4.1). The penalised fit that replaced it is timed by
//! the research bench (`examples/ng_str_fit_bench.rs`, `examples/str_fit_cost_probe.rs`).
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
//! The drawn-evidence generator is the one the module's own positive control uses, moved out of
//! `#[cfg(test)]` and behind the `bench-fixtures` feature: `fit::bench_fixtures::draw_cohort` (a
//! cohort drawn at a known error rate, heterozygosity and inbreeding). **One generator, two
//! callers** — a benchmark drawn differently from the oracle would be timing a workload no test
//! has ever checked.
//!
//! # The two axes, and why both are here
//!
//! `CLAUDE.md` §0 commits this caller to degrading gracefully from **one sample to several
//! thousand** and from **three reads a position to several hundred**. The fit is therefore swept
//! over sample count as well as over kept positions, because a change that pays at eight samples
//! and loses at one is a change this project has to be able to see.
//!
//! # Threads
//!
//! The fit parallelises over the global rayon pool, so a number taken on an 18-core laptop
//! and one taken in the 8-core container are not comparable. Every timed body therefore runs
//! inside a pool of a **fixed** size — four threads, matching every other measurement in the
//! review — overridable with `NG_JOINT_FIT_BENCH_THREADS`. Sweeping that variable over
//! 1, 2, 4, 8, 18 against the same benchmark is the thread-count sweep §3 item 4 asks for.
//!
//! # Sizes: why these, and why they are small
//!
//! One whole run of this file is minutes. That is bought by cutting the *number of times* the
//! estimator evaluates its objective — eight passes instead of two hundred — while leaving every
//! inner loop at production's own shape, sixteen frequency nodes included.

use criterion::{Criterion, criterion_group, criterion_main};

#[cfg(feature = "bench-fixtures")]
use criterion::BenchmarkId;
#[cfg(feature = "bench-fixtures")]
use std::hint::black_box;
#[cfg(feature = "bench-fixtures")]
use std::time::Duration;

#[cfg(feature = "bench-fixtures")]
use pop_var_caller::parameter_estimation::joint::fit::bench_fixtures::{
    DrawnCohort, as_cohort, draw_cohort,
};
#[cfg(feature = "bench-fixtures")]
use pop_var_caller::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFit, JointFitConfig, StartEnding, StartingPoint as GenericStart,
    fit_jointly,
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
/// **Not the global pool.** `fit_jointly` parallelises over
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
// The judging passes: the expectation and the information together
// ---------------------------------------------------------------------------

/// The pass limit the judging fit is given, the final pass not counted.
///
/// **Ten, so that most passes sum the information.** The first cycle — three passes — is always
/// plain, because `maximise` starts judging only after a cycle has gained less than
/// `log_likelihood_stillness` a position; from the second cycle on every first pass judges, every
/// Newton try sums the information too, and so does the final pass, whatever happened before it.
/// Production sums it on nearly every pass: judging starts by about pass 13 of 131 at 59 tomato
/// samples, and every Newton try sums it as well.
#[cfg(feature = "bench-fixtures")]
const JUDGING_PASSES: u32 = 10;

/// **How many of the judging fit's passes must sum the information, the final pass included.**
///
/// Every case of the group measured 8 of 11 when it was written (2026-10-18). The check allows one
/// fewer, so that a change which moves one Newton try's outcome does not stop the bench, while a fit
/// that stopped judging — only its final pass summing the information, 1 of 11 — still fails it.
#[cfg(feature = "bench-fixtures")]
const FEWEST_INFORMATION_PASSES: usize = 7;

/// The fit as a production run looks after its first few cycles: **every cycle judges**, no
/// start ever converges, and nothing runs after the alternation.
///
/// - `log_likelihood_stillness` is infinite, so the first cycle's gain — any finite gain — turns
///   judging on, and the second cycle's first pass is the first that sums the information: the
///   earliest `maximise` allows.
/// - `decrement_threshold` is zero, so no judging pass finds the start converged and the fit runs
///   until too few passes remain for another cycle or Newton try. That can leave it one or two
///   short of [`JUDGING_PASSES`]; on these draws it is ten, which the setup checks, and the fit is
///   deterministic, so every iteration times the same passes. A drawn cohort that converged at
///   pass 5 would time a different amount of work in every tree whose numbers moved.
/// - Contamination is off. It runs once, after the alternation, and is not the code this group
///   watches; at ten passes it would hold a far larger share of a timed fit than of a production
///   one (see [`PASSES`] for the same caveat on the group above).
/// - The per-pass trace is on: it is how the setup check finds the first judging pass.
#[cfg(feature = "bench-fixtures")]
fn judging_config() -> JointFitConfig {
    JointFitConfig {
        max_passes: JUDGING_PASSES,
        log_likelihood_stillness: f64::INFINITY,
        decrement_threshold: 0.0,
        estimate_contamination: false,
        pass_trace: true,
        ..generic_config()
    }
}

/// **The population the judging fit's positions are drawn from**, chosen so that the positions
/// carry as many candidate alleles as a real cohort's: nine hundred and ninety in a thousand
/// invariant, one fixed for the non-reference base, the rest segregating at the rare-allele-heavy
/// density the groups above use.
///
/// What the fit's per-position work grows with is how many candidate alleles a position sums over
/// — one when no sample showed a non-reference read, two when they showed one base, three when
/// they showed two or more. On 59 tomato2 samples' census (two million positions, 2026-10-18)
/// **67.5, 27.6 and 4.9% of positions had one, two and three candidates**. This density, at
/// [`JUDGING_DEPTH`] and [`JUDGING_HOM_EXCESS`], draws 67.9, 26.7 and 5.4% at 32 samples; the
/// groups above' density, at eight reads, drew 53.0, 36.8 and 10.3% — 1.57 candidates a position
/// against the real cohort's 1.37, about 15% more work a position-sample. The match is at 32 samples only; fewer samples see fewer
/// non-reference reads a position and so fewer candidates, as a real cohort would.
#[cfg(feature = "bench-fixtures")]
fn judging_density() -> FrequencyDensity {
    FrequencyDensity {
        p_invariant: 0.99,
        p_fixed_alt: 0.001,
        a: 0.7,
        b: 2.5,
    }
}

/// Reads a sample puts on a position in the judging fit: five and a half, tomato2's median sample
/// depth (59 samples, 4.2 to 13.7 reads a position).
#[cfg(feature = "bench-fixtures")]
const JUDGING_DEPTH: f64 = 5.5;

/// Each drawn sample's homozygote excess in the judging fit: 0.9, near the 0.94 to 1.00 the fit
/// returns for the oracle's four tomato samples — a selfing crop's.
#[cfg(feature = "bench-fixtures")]
const JUDGING_HOM_EXCESS: f64 = 0.9;

/// The passes of the fit's only start that summed the information, the final pass included — read
/// from the fit's own count, `StartRecord::judging_passes`, which leaves the final pass out though
/// it always sums it.
#[cfg(feature = "bench-fixtures")]
fn information_pass_count(fit: &JointFit) -> usize {
    fit.starts[0].judging_passes as usize + 1
}

/// **The guard against timing a judging fit that did not judge.**
///
/// The information pass is what this group exists for, and nothing in a wall number says it ran:
/// a fit that never left its first cycle, or one that converged at the first judging pass, is
/// fast and finite. So the setup fits once and checks the claims the group's name makes — one
/// start, which ran to its pass limit without converging; judging began at pass 4 (the second
/// cycle, the earliest `maximise` allows); and most passes summed the information — beside the
/// older check's guards against an empty or shallow draw, and prints the count.
#[cfg(feature = "bench-fixtures")]
fn check_the_fit_judged(drawn: &DrawnCohort, positions: usize, samples: usize, depth: f64) {
    let mut cohort = as_cohort(&drawn.samples);
    let fit = fit_jointly(&mut cohort, &judging_config()).expect("a drawn cohort pools");
    assert_eq!(fit.noisy_posterior.len(), positions);
    assert_eq!(fit.rates.len(), samples);
    assert!(
        fit.log_likelihood.is_finite(),
        "the fit's log-likelihood came back as {}",
        fit.log_likelihood
    );
    assert!(
        fit.rates["s0"].value.positions_with_reads > (positions as u64) / 2,
        "the first sample carried reads at {} of {positions} positions, so the draw is not the \
         {depth}× one this bench names",
        fit.rates["s0"].value.positions_with_reads
    );
    assert_eq!(fit.starts.len(), 1, "the judging fit runs one start");
    assert_eq!(
        fit.starts[0].ended,
        StartEnding::AtTheLimit,
        "the judging fit's start was to stop at its pass limit"
    );
    assert!(
        !fit.converged && fit.passes == JUDGING_PASSES,
        "the judging fit was to run to its limit of {JUDGING_PASSES} passes; it ran {} and \
         converged = {}",
        fit.passes,
        fit.converged
    );
    let first_judged = fit
        .trace
        .iter()
        .find(|summary| summary.decrement.is_some())
        .map(|summary| summary.pass);
    assert_eq!(
        first_judged,
        Some(4),
        "judging was to begin with the second cycle's first pass"
    );
    let summed = information_pass_count(&fit);
    let all = fit.passes as usize + 1;
    assert!(
        summed >= FEWEST_INFORMATION_PASSES,
        "only {summed} of {all} passes summed the information; at least \
         {FEWEST_INFORMATION_PASSES} were to"
    );
    eprintln!(
        "ng_joint_fit/generic_judging/{samples}: {summed} of {all} passes summed the information"
    );
}

/// **The judging fit by sample count, at twenty thousand positions** — the bench that sees the
/// hot path of a production fit since the Newton step landed (2026-10-16), where nearly every pass
/// is a judging pass: about half its CPU in `one_position`, the expectation, and about two fifths
/// in the information scoring (`score_position`'s helpers and `InformationSums::add_position`),
/// measured on 59 tomato samples. The two groups above run one information pass in nine.
///
/// The points: one sample, the cohort-less case; eight, a small cohort; thirty-two, above the
/// twenty at which the information switches from the whole matrix to its blocks
/// (`information::FULL_MATRIX_SAMPLES`), which is the layout every cohort of production size
/// takes. The positions are drawn at tomato2's median depth, five and a half reads, with as many
/// candidate alleles a position as tomato2's 59 samples carry ([`judging_density`]).
///
/// **Twenty thousand positions cut a pass into chunks of 256**, the smallest the pass allows,
/// where production's two million positions make chunks of about 16,000: each chunk's setup — its
/// scratch, its information sums, and adding them into the pass's total — is therefore paid about
/// sixty times as often a position here as in production, and is over-weighted in this group
/// (adding the chunks' totals and the allocator were under 0.2% of the workers' samples at 32
/// samples, 2026-10-18). At one sample the work between passes
/// — the maximisation, the quadrature — is over-weighted too, since a real one-sample fit has a
/// hundred times the positions; that case checks the one-sample code path, not production's
/// proportions.
#[cfg(feature = "bench-fixtures")]
fn generic_judging(c: &mut Criterion) {
    const POSITIONS: usize = 20_000;
    let pool = fixed_pool();
    let config = judging_config();

    let mut group = c.benchmark_group("ng_joint_fit/generic_judging");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    for &samples in &[1_usize, 8, 32] {
        let drawn = draw_cohort(
            samples,
            POSITIONS,
            JUDGING_DEPTH,
            NOISE,
            judging_density(),
            JUDGING_HOM_EXCESS,
            0x9E37_79B9_7F4A_7C17,
        );
        pool.install(|| check_the_fit_judged(&drawn, POSITIONS, samples, JUDGING_DEPTH));
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
        generic_by_positions(c);
        generic_by_samples(c);
        generic_judging(c);
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = every_bench
}

criterion_main!(benches);
