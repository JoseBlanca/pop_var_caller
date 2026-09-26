//! **Where the repeat-tract fit's threads are better spent: on one repeat-length class's
//! tracts, or on several classes at once.**
//!
//! The parameters fit has two halves. The first reads ordinary positions; the second reads
//! repeat tracts, one **repeat-length class** at a time — every tract of one motif period whose
//! reference carries the same number of copies. A class is fitted from its own tracts and no
//! other's, and the curves that reconcile the classes are drawn afterwards from the finished
//! answers, so the classes are independent of each other for the whole of the fit.
//!
//! That leaves two grains a thread pool can be spread over, and this harness measures both:
//!
//! - **one class at a time, its tracts split across the pool** — what the fit did until
//!   2026-09-11;
//! - **`THREADS` classes at once, each on one thread** — `SsrFitConfig::strata_at_once` set to
//!   `THREADS`.
//!
//! **Both arms return the same bits** (`ssr_fit`'s own
//! `the_two_ways_of_spending_the_pool_give_the_same_bits`), so what is measured here is wall
//! time and peak memory and nothing else.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example ng_ssr_strata_threads
//! ```
//!
//! Knobs, all environment variables: `THREADS` (default 4), `SAMPLES` (8), `ARM`
//! (`across_classes` or `across_tracts`), `PROFILE` (`skewed`, `even` or `many_thin`),
//! `ROUNDS` (2) and `STARTS` (1). `ROUNDS` and `STARTS` are cut below production's five and
//! three so that a sweep is minutes rather than an hour; they divide both arms equally and so
//! cannot favour either.
//!
//! # Why the class sizes are unequal by default
//!
//! **Because the real ones are, by three orders of magnitude.** On GIAB HG002 the 8-copy
//! homopolymer class holds 4,194 tracts and the 38-copy class holds 5
//! (`doc/devel/ng/reports/str_slippage_curves_on_both_cohorts_2026-08-21.md` §2). No class is
//! split under the second arm, so a run of it cannot finish sooner than its largest class would
//! have taken on one thread — and a fixture of equal classes is the one shape that hides that
//! ceiling. `PROFILE=even` is there to show the difference.

use std::time::Instant;

use pop_var_caller::parameter_estimation::joint::census::Stratum;
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::{
    draw_stratum, spectrum_of,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    Slippage, SsrFitConfig, StartingPoint, StratumEvidence, StratumOutcome, fit_strata,
};

/// Tomato's own dinucleotide numbers, which every drawn class here is drawn at.
const SLIPPAGE: Slippage = Slippage {
    level: 0.08,
    shorter_share: 0.83,
    fall_off: 0.25,
};

/// How far either side of the reference length the draw spreads allele mass. Six is
/// production's `ALLELE_SPAN`: thirteen allele classes and the ninety-one genotypes of a
/// diploid over them, which is the trip count of the innermost loop.
const SPAN: i32 = 6;

/// Reads a sample puts on a tract. Six sits between tomato's three and the GIAB trio's thirty.
const DEPTH: u32 = 6;

const CONCENTRATION: f64 = 0.5;
const HOM_EXCESS: f64 = 0.4;

fn from_env<T: std::str::FromStr>(name: &str, fallback: T) -> T {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(fallback)
}

/// How many tracts each class holds, under each of the three shapes.
///
/// All three hold about the same number of tracts in total — 256, 253 and 240 — so the three
/// rows of a sweep are the same amount of arithmetic arranged differently.
fn tracts_a_class(profile: &str) -> Vec<usize> {
    match profile {
        // **The shape the real data has**: one class carrying a quarter of the work, a tail of
        // thin ones. The largest is what bounds the second arm.
        "skewed" => vec![64, 48, 32, 24, 20, 16, 14, 12, 10, 8, 8],
        // Equal classes, where the second arm has nothing working against it.
        "even" => vec![23; 11],
        // **Many thin classes**, which is what a cohort at low depth looks like: on tomato 49
        // classes carry a read and 17 hold enough tracts to be fitted at all.
        "many_thin" => vec![10; 24],
        other => panic!("PROFILE must be skewed, even or many_thin, not {other}"),
    }
}

fn peak_rss_mib() -> f64 {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return f64::NAN;
    };
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:")
            && let Some(kib) = rest
                .split_whitespace()
                .next()
                .and_then(|it| it.parse::<f64>().ok())
        {
            return kib / 1024.0;
        }
    }
    f64::NAN
}

fn main() {
    let threads: usize = from_env("THREADS", 4);
    let samples: usize = from_env("SAMPLES", 8);
    let rounds: u32 = from_env("ROUNDS", 2);
    let starts: usize = from_env("STARTS", 1);
    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "skewed".to_string());
    let arm = std::env::var("ARM").unwrap_or_else(|_| "across_classes".to_string());
    let strata_at_once = match arm.as_str() {
        "across_classes" => std::num::NonZeroUsize::new(threads).expect("THREADS is at least one"),
        "across_tracts" => std::num::NonZeroUsize::MIN,
        other => panic!("ARM must be across_classes or across_tracts, not {other}"),
    };

    let sizes = tracts_a_class(&profile);
    let spectrum = spectrum_of((2 * SPAN + 1) as usize);
    let strata: Vec<StratumEvidence> = sizes
        .iter()
        .enumerate()
        .map(|(which, tracts)| {
            let mut evidence = draw_stratum(
                SLIPPAGE,
                &spectrum,
                CONCENTRATION,
                HOM_EXCESS,
                *tracts,
                samples,
                DEPTH,
                SPAN,
                0x5EED_0200 + which as u64,
            );
            evidence.stratum = Stratum {
                period: 2,
                reference_repeats: 8 + which as u64,
            };
            evidence
        })
        .collect();
    let drawn_rss = peak_rss_mib();

    let mut starting_points = StartingPoint::spanning_the_monomorphic_range();
    starting_points.truncate(starts.max(1));
    let config = SsrFitConfig {
        max_rounds: rounds,
        starting_points,
        // Eight is production's floor. Every class in every profile clears it, so nothing here
        // is timing a refusal.
        refusal_floor: 8,
        strata_at_once,
        ..SsrFitConfig::default()
    };

    let excess = vec![HOM_EXCESS; samples];
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a pool of the requested width");

    let started = Instant::now();
    let outcomes = pool.install(|| fit_strata(&strata, &excess, &config));
    let seconds = started.elapsed().as_secs_f64();

    // **The guard against timing a fit that fitted nothing.** A refusal returns in microseconds
    // and would read as a very fast arm.
    let fitted = outcomes
        .iter()
        .filter(|outcome| matches!(outcome, StratumOutcome::Fitted(_)))
        .count();
    assert_eq!(
        fitted,
        sizes.len(),
        "every drawn class clears the floor and must be fitted on its own tracts"
    );

    let tracts: usize = sizes.iter().sum();
    println!(
        "arm\tprofile\tthreads\tsamples\tclasses\ttracts\tseconds\tpeak_rss_mib\tafter_draw_mib"
    );
    println!(
        "{arm}\t{profile}\t{threads}\t{samples}\t{}\t{tracts}\t{seconds:.2}\t{:.0}\t{drawn_rss:.0}",
        sizes.len(),
        peak_rss_mib(),
    );
}
