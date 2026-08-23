//! Three ways to predict a panel's allele-frequency spectrum, timed against each other and
//! against the exact answer.
//!
//! Step D of the genotype prior fits two numbers — the chromosomes' worth of prior belief on the
//! reference allele and on the alternatives — by searching for the pair whose predicted spectrum
//! matches the one the parameter pre-pass measured. The prediction is the objective of that
//! search, so a fit pays it about a hundred times over, and the shipped version is cubic in the
//! number of samples. This harness asks what that costs and what buying it back costs in
//! accuracy.
//!
//! The three:
//!
//! - **exact** — `ng::calling::genotype_prior::seed_generic::fill_expected_spectrum`, the shipped
//!   one. One exponential per term.
//! - **recurrence** — the same sum, with the innermost factor stepped by an exact ratio instead of
//!   rebuilt from logarithms. Each term is a beta-binomial weight times a hypergeometric one, both
//!   genuine probabilities, so nothing overflows on the way. Arithmetically the same answer.
//! - **tail-trimmed** — the recurrence, plus: stop summing over how many individuals are inbred
//!   once that split is far enough out in its own binomial tail. Each dropped split can move a
//!   class by no more than its own probability, so the error has a bound that is printed beside
//!   the measured one.
//!
//! Run it with the panel sizes to walk, e.g.
//! `cargo run --release --example ng_spectrum_projection_cost -- 26 63 200 400 800 1600`.

use std::time::Instant;

use pop_var_caller::genetics::lgamma;
use pop_var_caller::ng::calling::genotype_prior::seed_generic::fill_expected_spectrum;
use pop_var_caller::ng::types::InbreedingF;

/// `ln` of the chance that exactly `inbred` of `individuals` have two copies that are one
/// ancestral copy counted twice.
fn log_branch_split(individuals: usize, inbred: usize, inbreeding: f64) -> Option<f64> {
    if inbreeding == 0.0 {
        return (inbred == 0).then_some(0.0);
    }
    if inbreeding == 1.0 {
        return (inbred == individuals).then_some(0.0);
    }
    let outbred = individuals - inbred;
    Some(
        lgamma(individuals as f64 + 1.0)
            - lgamma(inbred as f64 + 1.0)
            - lgamma(outbred as f64 + 1.0)
            + inbred as f64 * inbreeding.ln()
            + outbred as f64 * (1.0 - inbreeding).ln(),
    )
}

/// The same sum the shipped function computes, with two changes that do not alter the model.
///
/// `branch_tail_tolerance` drops branch splits whose probability is below that fraction of the
/// most likely split's; pass `0.0` to keep every one of them, which is the *recurrence* variant.
fn fill_by_recurrence(
    alpha_ref: f64,
    alpha_alt: f64,
    individuals: u32,
    inbreeding: f64,
    branch_tail_tolerance: f64,
    out: &mut [f64],
) {
    let n = individuals as usize;
    out.fill(0.0);
    if alpha_alt == 0.0 {
        out[0] = 1.0;
        return;
    }
    let concentration_total = alpha_ref + alpha_alt;
    let log_pair_constant = lgamma(concentration_total) - lgamma(alpha_alt) - lgamma(alpha_ref);
    let log_factorial: Vec<f64> = (0..=2 * n + 1).map(|k| lgamma(k as f64 + 1.0)).collect();
    let log_binomial = |top: usize, chosen: usize| {
        log_factorial[top] - log_factorial[chosen] - log_factorial[top - chosen]
    };

    // The most likely split, so the tail can be measured against it rather than against zero.
    let splits: Vec<Option<f64>> = (0..=n)
        .map(|inbred| log_branch_split(n, inbred, inbreeding))
        .collect();
    let peak = splits
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let floor = if branch_tail_tolerance > 0.0 {
        peak + branch_tail_tolerance.ln()
    } else {
        f64::NEG_INFINITY
    };

    for (inbred, split) in splits.iter().enumerate() {
        let Some(log_branch_weight) = *split else {
            continue;
        };
        if log_branch_weight < floor || log_branch_weight < (1e-300_f64).ln() {
            continue;
        }
        let branch_weight = log_branch_weight.exp();
        let distinct = 2 * n - inbred;
        let singles = distinct - inbred;
        let log_draw_constant = log_pair_constant - lgamma(concentration_total + distinct as f64);

        for alternative_draws in 0..=distinct {
            // How likely it is that exactly this many of the distinct chromosomes are
            // alternative: a beta-binomial weight, and a probability in its own right.
            let beta_binomial = (log_binomial(distinct, alternative_draws)
                + lgamma(alpha_alt + alternative_draws as f64)
                + lgamma(alpha_ref + (distinct - alternative_draws) as f64)
                + log_draw_constant)
                .exp();
            if beta_binomial == 0.0 {
                continue;
            }
            let lowest = alternative_draws.saturating_sub(singles);
            let highest = inbred.min(alternative_draws);
            // How many of the duplicated chromosomes are among the alternative draws, given how
            // many draws there were: a hypergeometric weight, stepped by its own exact ratio.
            //
            // **The walk starts at the mode and goes out both ways**, which is not a refinement.
            // Started at `lowest`, the first weight underflows to zero long before the ones at the
            // mode become small — the whole row then multiplies zero and vanishes. Measured at
            // 1,600 individuals: one class came back 5.7e-16 against its true 6.1e-7, and the
            // spectrum lost 3 parts in 10,000 of its mass with every entry still finite.
            let mode =
                (((alternative_draws + 1) * (inbred + 1)) / (distinct + 2)).clamp(lowest, highest);
            let at_mode = (log_binomial(inbred, mode)
                + log_binomial(singles, alternative_draws - mode)
                - log_binomial(distinct, alternative_draws))
            .exp();
            out[alternative_draws + mode] += branch_weight * beta_binomial * at_mode;

            let mut climbing = at_mode;
            for doubled in mode..highest {
                let rise = ((inbred - doubled) * (alternative_draws - doubled)) as f64;
                let fall = ((doubled + 1) * ((singles + doubled + 1) - alternative_draws)) as f64;
                climbing *= rise / fall;
                if climbing == 0.0 {
                    break;
                }
                out[alternative_draws + doubled + 1] += branch_weight * beta_binomial * climbing;
            }
            let mut falling = at_mode;
            for doubled in (lowest..mode).rev() {
                let rise = ((inbred - doubled) * (alternative_draws - doubled)) as f64;
                let fall = ((doubled + 1) * ((singles + doubled + 1) - alternative_draws)) as f64;
                falling *= fall / rise;
                if falling == 0.0 {
                    break;
                }
                out[alternative_draws + doubled] += branch_weight * beta_binomial * falling;
            }
        }
    }
}

fn worst_relative_gap(got: &[f64], want: &[f64]) -> f64 {
    got.iter()
        .zip(want)
        .map(|(a, b)| if *b > 0.0 { (a - b).abs() / b } else { a.abs() })
        .fold(0.0_f64, f64::max)
}

fn main() {
    let sizes: Vec<u32> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("panel sizes are whole numbers"))
        .collect();
    let sizes = if sizes.is_empty() {
        vec![26, 63, 200, 400, 800]
    } else {
        sizes
    };

    // Tomato's fitted diversity and an inbreeding coefficient in its fitted range: the corner the
    // projection is actually aimed at, and the one where the branch sum is widest.
    let (alpha_ref, alpha_alt, inbreeding) = (1.0, 6e-4, 0.8);
    let tolerances = [1e-18_f64, 1e-12, 1e-8, 1e-6];

    println!(
        "panel   exact       recurrence  |  {}",
        tolerances
            .iter()
            .map(|t| format!("trim {t:>7.0e}"))
            .collect::<Vec<_>>()
            .join("  ")
    );
    println!("{}", "-".repeat(52 + 14 * tolerances.len()));

    for &individuals in &sizes {
        let classes = 2 * individuals as usize + 1;
        let mut exact = vec![0.0; classes];
        let started = Instant::now();
        fill_expected_spectrum(
            alpha_ref,
            alpha_alt,
            individuals,
            InbreedingF::try_new(inbreeding).unwrap(),
            &mut exact,
        );
        let exact_time = started.elapsed();

        let mut recurrence = vec![0.0; classes];
        let started = Instant::now();
        fill_by_recurrence(
            alpha_ref,
            alpha_alt,
            individuals,
            inbreeding,
            0.0,
            &mut recurrence,
        );
        let recurrence_time = started.elapsed();
        let recurrence_gap = worst_relative_gap(&recurrence, &exact);
        {
            let exact_mass: f64 = exact.iter().sum();
            let recurrence_mass: f64 = recurrence.iter().sum();
            let (worst_class, _) = exact
                .iter()
                .zip(&recurrence)
                .enumerate()
                .map(|(k, (a, b))| (k, if *a > 0.0 { (a - b).abs() / a } else { b.abs() }))
                .fold((0, 0.0_f64), |acc, x| if x.1 > acc.1 { x } else { acc });
            eprintln!(
                "  N={individuals} exact_mass-1={:+.3e} recurrence_mass-1={:+.3e} worst_class={worst_class} exact[{worst_class}]={:.6e} recurrence[{worst_class}]={:.6e} exact_min={:.3e}",
                exact_mass - 1.0,
                recurrence_mass - 1.0,
                exact[worst_class],
                recurrence[worst_class],
                exact.iter().copied().fold(f64::INFINITY, f64::min),
            );
        }

        let mut trimmed_cells = Vec::new();
        for &tolerance in &tolerances {
            let mut trimmed = vec![0.0; classes];
            let started = Instant::now();
            fill_by_recurrence(
                alpha_ref,
                alpha_alt,
                individuals,
                inbreeding,
                tolerance,
                &mut trimmed,
            );
            let elapsed = started.elapsed();
            let gap = worst_relative_gap(&trimmed, &exact);
            let mass: f64 = trimmed.iter().sum();
            trimmed_cells.push(format!(
                "{:>7.1?} {gap:>7.0e} {:>+8.0e}",
                elapsed,
                mass - 1.0
            ));
        }

        println!(
            "{individuals:<7} {exact_time:>10.1?}  {recurrence_time:>10.1?} {recurrence_gap:>7.0e}  |  {}",
            trimmed_cells.join("  ")
        );
    }

    println!();
    println!(
        "columns: wall clock; worst class-by-class relative gap against exact; and how far the \
         trimmed classes depart from summing to one."
    );
    println!("α_ref {alpha_ref}, α_alt {alpha_alt} (tomato's fitted diversity), F {inbreeding}.");
}
