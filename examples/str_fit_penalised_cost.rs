//! **Where penalised counting's time goes, against the counting fit's**, on one cell of the fitting
//! bench's tomato-like truth (`doc/devel/ng/research/str_fit/track2_penalised_em_2026-10-12.md`,
//! the follow-up section). Every time is wall time on a shared machine: provisional.
//!
//! Draws the cell as the bench does (at most 100 loci a class, the shipped window), then runs:
//!
//! - `pass`: one round of the counting fit's E-step over every class, one class after another, at
//!   the starting numbers: what a round costs either fit;
//! - `B`: the counting fit as the bench runs it (pairs of samples, up to 10 readings, 2,000 rounds,
//!   tolerance 10⁻⁷, classes side by side, a class below the refusal floor not fitted);
//! - `D_every`: penalised counting choosing strengths every round (the first schedule), every
//!   class counted every round;
//! - `D_steady`: holding a curve's strength once it moves by less than 0.1;
//! - `D_reuse`: that, and a settled class keeping its last counts (the defaults);
//! - `D_errors`: `D_reuse` with standard errors.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example str_fit_penalised_cost -- \
//!     <period> <samples> <depth> <seed> [arms, default "pass,B,D_every,D_steady,D_reuse,D_errors"]
//! ```
//!
//! Prints one line an arm: `arm  seconds  rounds  ...`, the columns named in a header line.

use std::time::Instant;

use rayon::prelude::*;

use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::cohort_draw::{
    CohortModel, LociPerClass, StutterRegime, draw_period,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::expected_counts::{
    expected_counts, fit_class_by_counting_with_pairs,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, fit_period_by_penalised_counting,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    ALLELE_SPAN, DEFAULT_REFUSAL_FLOOR, Slippage, SsrFitConfig,
};

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    let argument = |index: usize, name: &str| -> String {
        arguments
            .get(index)
            .unwrap_or_else(|| panic!("missing {name}"))
            .clone()
    };
    let period: u8 = argument(1, "period").parse().expect("a period");
    let samples: usize = argument(2, "samples").parse().expect("samples");
    let depth: f64 = argument(3, "depth").parse().expect("depth");
    let seed: u64 = argument(4, "seed").parse().expect("seed");
    let arms: Vec<String> = arguments
        .get(5)
        .map_or("pass,B,D_every,D_steady,D_reuse,D_errors", String::as_str)
        .split(',')
        .map(str::to_owned)
        .collect();
    let mut model = CohortModel::tomato_like(period, StutterRegime::Pcr);
    model.samples = samples;
    model.depth = depth;
    model.seed = seed;
    model.loci = LociPerClass::Tomato2 { at_most: Some(100) };
    let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
    let classes = &drawn.strata;
    let excess = &drawn.homozygote_excess;
    let config = SsrFitConfig {
        allele_span: ALLELE_SPAN,
        ..SsrFitConfig::default()
    };
    let start = Slippage {
        level: 0.05,
        shorter_share: 0.5,
        fall_off: 0.3,
    };
    let cell = format!("p{period}_n{samples}_d{depth}_s{seed}");
    println!(
        "cell\tarm\tseconds\trounds\trounds_per_counting\texpectation_s\tmaximisation_s\tstrength_search_s\trounds_choosing\tpairs_s\tpair_readings\terrors_s\tclass_rounds_counted"
    );
    for arm in &arms {
        match arm.as_str() {
            "pass" => {
                let shares =
                    vec![1.0 / f64::from(2 * ALLELE_SPAN + 1); (2 * ALLELE_SPAN + 1) as usize];
                let started = Instant::now();
                for evidence in classes {
                    let _ = expected_counts(evidence, start, &shares, 1.0, excess, ALLELE_SPAN);
                }
                println!(
                    "{cell}\tpass\t{:.3}\t1\t\t\t\t\t\t\t\t",
                    started.elapsed().as_secs_f64()
                );
            }
            "B" => {
                let started = Instant::now();
                let rounds: Vec<usize> = classes
                    .par_iter()
                    .map(|evidence| {
                        if evidence.tracts_with_reads() < DEFAULT_REFUSAL_FLOOR {
                            return 0;
                        }
                        fit_class_by_counting_with_pairs(
                            evidence, start, excess, &config, 2000, 1e-7, 10,
                        )
                        .rounds_in_all
                    })
                    .collect();
                println!(
                    "{cell}\tB\t{:.3}\t{}\t{}\t\t\t\t\t\t\t",
                    started.elapsed().as_secs_f64(),
                    rounds.iter().sum::<usize>(),
                    rounds
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            arm @ ("D_every" | "D_steady" | "D_reuse" | "D_errors") => {
                let fit = PenalisedFitConfig {
                    allele_span: ALLELE_SPAN,
                    start,
                    max_rounds: 2000,
                    tolerance: 1e-7,
                    pair_alternations: 10,
                    strength_steady_within: if arm == "D_every" { 0.0 } else { 0.1 },
                    reuse_settled_counts: arm == "D_reuse" || arm == "D_errors",
                    standard_errors: arm == "D_errors",
                    ..PenalisedFitConfig::prototype()
                };
                let started = Instant::now();
                let fitted = fit_period_by_penalised_counting(classes, excess, &fit);
                println!(
                    "{cell}\t{arm}\t{:.3}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.3}\t{}\t{:.3}\t{}",
                    started.elapsed().as_secs_f64(),
                    fitted.rounds,
                    fitted
                        .rounds_per_counting
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                    fitted.expectation_seconds,
                    fitted.maximisation_seconds,
                    fitted.strength_search_seconds,
                    fitted.rounds_choosing_strength,
                    fitted.pair_seconds,
                    fitted.pair_alternations,
                    fitted.error_seconds,
                    fitted.class_rounds_counted
                );
            }
            other => panic!("no arm {other}"),
        }
    }
}
