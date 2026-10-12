//! **Penalised counting against the counting fit, on one drawn motif length** — the measurement
//! behind `doc/devel/ng/research/str_fit/track2_penalised_em_2026-10-12.md`.
//!
//! Draws every class of one motif length at known numbers that change smoothly with the reference
//! repeat count (8 to 20 repeats) and bend as they change: none is a straight line on the logit
//! scale, which is the shape the penalty leaves alone, and the fall-off dips and rises again. One
//! deliberate kink: the slip rate at 9 repeats is drawn at
//! 0.53 of its curve, as HG002's 9-repeat homopolymers sit. Loci per class fall from 2,000 at 8
//! repeats to 5 at 20, and 17 repeats has none. Then fits it three ways:
//!
//! - `counting`: each class on its own by the counting fit;
//! - `counting+blend`: the same, then the slip rate drawn across repeat count by the existing curve
//!   and blended with each class's own (`str_slippage_level_curve.md` §4 and §7); only the slip
//!   rate has this arm;
//! - `penalised`: every class at once by penalised counting.
//!
//! Prints one row per (setting, draw, method, class, number): the estimate and the truth, and one
//! cost row per (setting, draw, method).
//!
//! ```text
//! cargo run --release --features bench-fixtures --example str_fit_penalised_em -- \
//!     [allele span, default 4] [draws, default 3] [settings, default "1x3,50x30"] \
//!     [homozygote excess, default 0.3] [pair alternations, default 0] \
//!     [slip rate at 9 repeats as a share of its curve, default 0.53] \
//!     [cost of a bend beyond the knee per log unit, default 2] \
//!     [1 to weigh allele draws as clustered by locus, default 0] \
//!     [a penalty strength to hold every curve at, default: chosen]
//! ```
//!
//! With pair alternations above zero both fits read the between-locus share off pairs of samples
//! at a locus (`fit_class_by_counting_with_pairs`, and the penalised fit's own pairs route).

use std::time::Instant;

use pop_var_caller::float;
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::draw_stratum;
use pop_var_caller::parameter_estimation::joint::ssr_fit::composite_fit::{
    fit_class_by_counting, fit_class_by_counting_with_pairs,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, Smoothing, fit_period_by_penalised_counting,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    FittedCell, Slippage, SlippageCurveConfig, SsrFitConfig, Stratum, StratumEvidence, blend_level,
    choose_rise_shape,
};

const REPEATS: std::ops::RangeInclusive<u64> = 8..=20;
const KINK_AT: u64 = 9;
/// The slip rate at [`KINK_AT`] as a share of its curve: 0.53 unless the run says otherwise.
static KINK_FACTOR: std::sync::OnceLock<f64> = std::sync::OnceLock::new();

fn logistic(x: f64) -> f64 {
    1.0 / (1.0 + float::exp(-x))
}

/// Loci per class: 2,000 at 8 repeats, falling by 0.55 a repeat, at least 5; none at 17.
fn loci_at(repeats: u64) -> usize {
    match repeats {
        17 => 0,
        18 => 12,
        19 => 8,
        20 => 5,
        _ => (2000.0 * float::powi(0.55, (repeats - 8) as i32)).round() as usize,
    }
}

/// The class's true numbers.
struct Truth {
    slippage: Slippage,
    between_locus_share: f64,
    allele_length_variation: Vec<f64>,
}

fn truth_at(repeats: u64, span: i32) -> Truth {
    let x = (repeats - 8) as f64;
    let mut level = logistic(-5.3 + 0.45 * x - 0.012 * x * x);
    if repeats == KINK_AT {
        level *= KINK_FACTOR.get().copied().unwrap_or(0.53);
    }
    // Alleles spread wider at longer tracts, the shorter side a little more; the edge classes
    // carry everything at or past them, so they sit above a falling tail.
    let ratio = 0.15 + 0.035 * x;
    let mut shares: Vec<f64> = (-span..=span)
        .map(|offset| {
            let base = float::powi(ratio, offset.abs());
            let side = if offset < 0 { 1.3 } else { 1.0 };
            let beyond = if offset.abs() == span {
                1.0 / (1.0 - ratio)
            } else {
                1.0
            };
            base * side * beyond
        })
        .collect();
    let total: f64 = shares.iter().sum();
    for share in &mut shares {
        *share /= total;
    }
    Truth {
        slippage: Slippage {
            level,
            shorter_share: logistic(0.4 + 0.25 * x - 0.015 * x * x),
            fall_off: logistic(-1.2 - 0.12 * x + 0.012 * x * x),
        },
        between_locus_share: logistic(-0.85 - 0.15 * x + 0.004 * x * x),
        allele_length_variation: shares,
    }
}

fn draw_period(
    samples: usize,
    depth: u32,
    span: i32,
    homozygote_excess: f64,
    seed: u64,
) -> Vec<StratumEvidence> {
    REPEATS
        .map(|repeats| {
            let truth = truth_at(repeats, span);
            let concentration = 1.0 / truth.between_locus_share - 1.0;
            let mut evidence = draw_stratum(
                truth.slippage,
                &truth.allele_length_variation,
                concentration,
                homozygote_excess,
                loci_at(repeats),
                samples,
                depth,
                span,
                seed * 1000 + repeats,
            );
            evidence.stratum = Stratum {
                period: 1,
                reference_repeats: repeats,
            };
            evidence
        })
        .collect()
}

/// One class's estimate in the terms every arm reports.
struct Estimate {
    slippage: Slippage,
    /// `None` where the method reports none (penalised counting at one sample).
    between_locus_share: Option<f64>,
    allele_length_variation: Vec<f64>,
}

fn print_rows(
    setting: &str,
    seed: u64,
    method: &str,
    repeats: u64,
    loci: usize,
    estimate: &Estimate,
    span: i32,
) {
    let truth = truth_at(repeats, span);
    let row = |number: &str, value: f64, true_value: f64| {
        println!(
            "row\t{setting}\t{seed}\t{method}\t{repeats}\t{loci}\t{number}\t{value:.6e}\t{true_value:.6e}"
        );
    };
    row("level", estimate.slippage.level, truth.slippage.level);
    row(
        "shorter_share",
        estimate.slippage.shorter_share,
        truth.slippage.shorter_share,
    );
    row(
        "fall_off",
        estimate.slippage.fall_off,
        truth.slippage.fall_off,
    );
    if let Some(share) = estimate.between_locus_share {
        row("between_locus_share", share, truth.between_locus_share);
    }
    let distance: f64 = 0.5
        * estimate
            .allele_length_variation
            .iter()
            .zip(&truth.allele_length_variation)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>();
    row("allele_length_distance", distance, 0.0);
}

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    let span: i32 = arguments
        .get(1)
        .map_or(4, |a| a.parse().expect("an allele span"));
    let draws: u64 = arguments
        .get(2)
        .map_or(3, |a| a.parse().expect("a number of draws"));
    let settings: Vec<(usize, u32)> = arguments
        .get(3)
        .map_or("1x3,50x30", String::as_str)
        .split(',')
        .map(|setting| {
            let (samples, depth) = setting.split_once('x').expect("samples x depth");
            (
                samples.parse().expect("samples"),
                depth.parse().expect("depth"),
            )
        })
        .collect();
    let config = SsrFitConfig {
        allele_span: span,
        ..SsrFitConfig::default()
    };
    let homozygote_excess: f64 = arguments
        .get(4)
        .map_or(0.3, |a| a.parse().expect("a homozygote excess"));
    let pair_alternations: usize = arguments
        .get(5)
        .map_or(0, |a| a.parse().expect("a number of pair alternations"));
    let kink_factor: f64 = arguments
        .get(6)
        .map_or(0.53, |a| a.parse().expect("a kink factor"));
    KINK_FACTOR.set(kink_factor).expect("set once");
    let kink_tail: f64 = arguments
        .get(7)
        .map_or(2.0, |a| a.parse().expect("a tail coefficient"));
    let weigh_allele_draws_by_locus = arguments.get(8).is_some_and(|a| a == "1");
    // A strength held for every curve instead of chosen; for measuring what the choice gives.
    let smoothing = arguments.get(9).map_or(Smoothing::MarginalLikelihood, |a| {
        Smoothing::Fixed(a.parse().expect("a penalty strength"))
    });
    let fit = PenalisedFitConfig {
        pair_alternations,
        kink_tail,
        weigh_allele_draws_by_locus,
        smoothing,
        ..PenalisedFitConfig::default()
    };
    println!("row\tsetting\tdraw\tmethod\trepeats\tloci\tnumber\testimate\ttruth");
    println!(
        "cost\tsetting\tdraw\tmethod\tseconds\trounds\texpectation_seconds\tmaximisation_seconds\tlast_change"
    );
    for (samples, depth) in settings {
        let setting = format!("{samples}x{depth}@{homozygote_excess}p{pair_alternations}");
        let excess = vec![homozygote_excess; samples];
        for seed in 1..=draws {
            let classes = draw_period(samples, depth, span, homozygote_excess, seed);

            // Each class on its own.
            let started = Instant::now();
            let mut rounds = Vec::new();
            let mut cells = Vec::new();
            let mut alone = Vec::new();
            for evidence in &classes {
                if evidence.rows_with_reads() == 0 {
                    alone.push(None);
                    continue;
                }
                let counted = if pair_alternations > 0 {
                    fit_class_by_counting_with_pairs(
                        evidence,
                        fit.start,
                        &excess,
                        &config,
                        fit.max_rounds,
                        fit.tolerance,
                        pair_alternations,
                    )
                } else {
                    fit_class_by_counting(
                        evidence,
                        fit.start,
                        &excess,
                        &config,
                        fit.max_rounds,
                        fit.tolerance,
                    )
                };
                rounds.push(counted.rounds);
                if counted.slippage.level > 0.0 {
                    cells.push(FittedCell {
                        repeats: evidence.stratum.reference_repeats,
                        level: counted.slippage.level,
                        slipped_reads: counted.slippage.level * evidence.spanning_reads() as f64,
                    });
                }
                alone.push(Some(counted));
            }
            let counting_seconds = started.elapsed().as_secs_f64();
            println!(
                "cost\t{setting}\t{seed}\tcounting\t{counting_seconds:.3}\t{}\t\t",
                rounds
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            );
            let curve_config = SlippageCurveConfig::default();
            let curve = choose_rise_shape(&[cells.clone()], &curve_config).ok();
            for (evidence, counted) in classes.iter().zip(&alone) {
                let Some(counted) = counted else { continue };
                let repeats = evidence.stratum.reference_repeats;
                let loci = evidence.tracts.len();
                let estimate = Estimate {
                    slippage: counted.slippage,
                    between_locus_share: Some(1.0 / (1.0 + counted.concentration)),
                    allele_length_variation: counted.allele_length_variation.clone(),
                };
                print_rows(&setting, seed, "counting", repeats, loci, &estimate, span);
                let cell = cells.iter().find(|cell| cell.repeats == repeats).copied();
                let curve_here = curve.as_ref().and_then(|c| c.by_group[0].as_ref());
                if let Some(blended) = blend_level(cell, curve_here, repeats, &curve_config) {
                    println!(
                        "row\t{setting}\t{seed}\tcounting+blend\t{repeats}\t{loci}\tlevel\t{:.6e}\t{:.6e}",
                        blended.level,
                        truth_at(repeats, span).slippage.level
                    );
                }
            }

            // Every class at once.
            let started = Instant::now();
            let together = fit_period_by_penalised_counting(&classes, &excess, &config, &fit);
            let penalised_seconds = started.elapsed().as_secs_f64();
            println!(
                "cost\t{setting}\t{seed}\tpenalised\t{penalised_seconds:.3}\t{}\t{:.3}\t{:.3}\t{:.1e}",
                together.rounds,
                together.expectation_seconds,
                together.maximisation_seconds,
                together.last_change
            );
            for (evidence, fitted) in classes.iter().zip(&together.classes) {
                let Some(fitted) = fitted else { continue };
                let estimate = Estimate {
                    slippage: fitted.slippage,
                    between_locus_share: fitted.between_locus_share,
                    allele_length_variation: fitted.allele_length_variation.clone(),
                };
                print_rows(
                    &setting,
                    seed,
                    "penalised",
                    evidence.stratum.reference_repeats,
                    evidence.tracts.len(),
                    &estimate,
                    span,
                );
            }
            if let Some(curves) = &together.curves {
                // The empty class at 17 repeats, from the curves alone.
                let estimate = Estimate {
                    slippage: Slippage {
                        level: curves.level.at(17),
                        shorter_share: curves.shorter_share.at(17),
                        fall_off: curves.fall_off.at(17),
                    },
                    between_locus_share: curves
                        .between_locus_share
                        .as_ref()
                        .map(|curve| curve.at(17)),
                    allele_length_variation: curves
                        .allele_length_variation
                        .iter()
                        .map(|curve| curve.at(17))
                        .collect(),
                };
                print_rows(&setting, seed, "penalised", 17, 0, &estimate, span);
                println!(
                    "curve\t{setting}\t{seed}\tsmoothing\tlevel {:.3e} shorter {:.3e} fall_off {:.3e} between {:.3e} shares {:.3e}\teased level at {:?}",
                    curves.level.smoothing,
                    curves.shorter_share.smoothing,
                    curves.fall_off.smoothing,
                    curves
                        .between_locus_share
                        .as_ref()
                        .map_or(f64::NAN, |curve| curve.smoothing),
                    curves.allele_length_variation[0].smoothing,
                    curves.level.eased_at
                );
            }
        }
    }
}
