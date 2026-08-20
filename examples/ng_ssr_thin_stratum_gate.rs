//! How thin a repeat-tract stratum may be before its fit stops meaning anything.
//!
//! **The question this answers.** Nothing below 50 tracts is fitted today, so a stratum thinner
//! than that has no answer of its own and therefore nothing to contribute to its motif period's
//! curves. The design wants every stratum to contribute, weighted by how precisely it measured
//! itself — a 3-tract stratum moving a curve by about one part in nine hundred against a
//! well-measured one. **That is only true if a fit that thin returns something noisy rather than
//! something broken**, and a climb that fails to converge, or one that returns a slippage level
//! of effectively zero, is a different failure from a noisy one
//! (`doc/devel/ng/spec/str_slippage_level_curve.md` §11).
//!
//! **What it does.** Draws strata from a known truth at a range of tract counts, fits each one on
//! its own tracts with the curve switched off, and reports, per tract count: how often the climb
//! settled, how often the level collapsed, how far the median fit sits from the truth, and how
//! much weight a stratum that thin would carry against a well-measured one.
//!
//! **Two shapes of cohort, at the two ends of the range this caller works over**: one sample at
//! 30 reads a tract, and 63 samples at 3 reads each. A stratum of the same tract count holds
//! about twenty times more reads in the second, so the two ends answer the question differently
//! and both are reported.
//!
//! Run: `./scripts/dev.sh cargo run --release --features bench-fixtures \
//!      --example ng_ssr_thin_stratum_gate`

use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::joint::slippage_curve::SlippageCurveConfig;
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit::{
    self, Slippage, SsrFitConfig, StratumOutcome, bench_fixtures,
};

/// The truth every stratum here is drawn from — a mid-range slippage level with the direction
/// split and fall-off both cohorts sit near.
const TRUTH: Slippage = Slippage {
    level: 0.02,
    shorter_share: 0.65,
    fall_off: 0.40,
};

/// How many tracts each drawn stratum holds. **3 is a handful and 400 is well measured**; the
/// current refusal floor sits at 50.
const TRACT_COUNTS: [usize; 8] = [3, 5, 8, 12, 20, 30, 50, 400];

/// How many strata are drawn at each tract count, each from its own seed.
const DRAWS: usize = 30;

/// One cell of the sweep: a tract count, at one shape of cohort.
struct Fitted {
    converged: bool,
    level: f64,
    shorter_share: f64,
    fall_off: f64,
    slipped_reads: f64,
}

fn main() {
    let spectrum = bench_fixtures::spectrum_of(3);
    let config = SsrFitConfig {
        allele_span: 1,
        // Fit everything, however thin — the point is to see what a thin fit returns.
        refusal_floor: 1,
        curve: SlippageCurveConfig {
            draw_curves: false,
            ..SlippageCurveConfig::default()
        },
        ..SsrFitConfig::default()
    };

    for (cohort, samples, depth) in [
        ("one sample at 30 reads a tract", 1_usize, 30_u32),
        ("63 samples at 3 reads each", 63, 3),
    ] {
        println!("\n=== {cohort} ===");
        println!(
            "  drawn from a level of {:.3}, a direction split of {:.2} and a fall-off of {:.2}; \
             {DRAWS} strata a row",
            TRUTH.level, TRUTH.shorter_share, TRUTH.fall_off
        );

        // **The widest row first**, so every thinner row's weight can be reported against a
        // well-measured stratum's rather than in units nobody has.
        let mut rows: Vec<(usize, Vec<Fitted>)> = TRACT_COUNTS
            .iter()
            .rev()
            .map(|tracts| {
                let at = Instant::now();
                let fits = draw_and_fit(*tracts, samples, depth, &spectrum, &config);
                eprintln!("    {tracts} tracts: {:.1} s", at.elapsed().as_secs_f64());
                (*tracts, fits)
            })
            .collect();
        rows.reverse();
        let well_measured = rows
            .last()
            .map(|(_, fits)| median_weight(fits))
            .unwrap_or(1.0);

        println!(
            "  {:>7} {:>9} {:>9} {:>10} {:>10} {:>8} {:>10} {:>8} {:>9} {:>11}",
            "tracts",
            "reads",
            "settled",
            "collapsed",
            "level",
            "spread",
            "9 in 10 by",
            "split",
            "fall-off",
            "its weight"
        );
        for (tracts, fits) in &rows {
            report(*tracts, fits, well_measured, samples, depth);
        }
    }
}

/// Draw `DRAWS` strata of `tracts` tracts and fit each on its own tracts.
fn draw_and_fit(
    tracts: usize,
    samples: usize,
    depth: u32,
    spectrum: &[f64],
    config: &SsrFitConfig,
) -> Vec<Fitted> {
    let mut fits = Vec::with_capacity(DRAWS);
    for draw in 0..DRAWS {
        let evidence = bench_fixtures::draw_stratum(
            TRUTH,
            spectrum,
            3.0,
            0.4,
            tracts,
            samples,
            depth,
            1,
            (tracts as u64) * 1_000 + draw as u64,
        );
        let outcomes = ssr_fit::fit_strata(&[evidence], &vec![0.4; samples], config);
        let StratumOutcome::Fitted(fit) = &outcomes[0] else {
            continue;
        };
        let Some(slippage) = fit.slippage[0] else {
            continue;
        };
        fits.push(Fitted {
            converged: fit.converged,
            level: slippage.level,
            shorter_share: slippage.shorter_share,
            fall_off: slippage.fall_off,
            slipped_reads: slippage.level * fit.reads_crossing as f64,
        });
    }
    fits
}

fn quantile(values: &[f64], part: f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let at = ((sorted.len() - 1) as f64 * part).round() as usize;
    sorted[at]
}

/// The weight the share curves would give the median stratum of this row: the inverse variance
/// of its logit, `slipped reads × p × (1 − p)`.
fn median_weight(fits: &[Fitted]) -> f64 {
    let weights: Vec<f64> = fits
        .iter()
        .map(|fit| {
            let share = fit.shorter_share.clamp(1e-4, 1.0 - 1e-4);
            fit.slipped_reads.max(1.0) * share * (1.0 - share)
        })
        .collect();
    quantile(&weights, 0.5).max(1e-12)
}

/// One row of the sweep.
fn report(tracts: usize, fits: &[Fitted], well_measured: f64, samples: usize, depth: u32) {
    if fits.is_empty() {
        println!("  {tracts:>7}  nothing was fitted at all");
        return;
    }
    let share = |count: usize| count as f64 / fits.len() as f64;
    let settled = share(fits.iter().filter(|fit| fit.converged).count());
    // **Collapsed** is a fitted level below a tenth of the truth — a stratum that says slippage
    // barely happens where it happens two reads in a hundred. Every one of these was at 1e-4 or
    // below, so the two thresholds the specification asks about are one column.
    let collapsed = share(
        fits.iter()
            .filter(|fit| fit.level < TRUTH.level / 10.0)
            .count(),
    );

    let levels: Vec<f64> = fits.iter().map(|fit| fit.level).collect();
    let median = quantile(&levels, 0.5);
    let spread = quantile(&levels, 0.75) / quantile(&levels, 0.25).max(1e-12);
    // The over-estimating tail, which is the one that costs: a level fitted high carries a
    // weight in proportion to it.
    let nine_in_ten = quantile(&levels, 0.9) / TRUTH.level;

    let median_split = quantile(
        &fits.iter().map(|fit| fit.shorter_share).collect::<Vec<_>>(),
        0.5,
    );
    let median_fall_off = quantile(
        &fits.iter().map(|fit| fit.fall_off).collect::<Vec<_>>(),
        0.5,
    );
    let weight = format!("1 in {:.0}", well_measured / median_weight(fits));

    let reads = tracts * samples * depth as usize;
    println!(
        "  {tracts:>7} {reads:>9} {:>8.0}% {:>9.0}% {median:>10.5} {spread:>7.2}x \
         {nine_in_ten:>9.2}x {median_split:>8.2} {median_fall_off:>9.2} {weight:>11}",
        settled * 100.0,
        collapsed * 100.0,
    );
}
