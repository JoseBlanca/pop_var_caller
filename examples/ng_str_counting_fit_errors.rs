//! **Are the counting fit's standard errors honest?** Draws many repeat classes from a known truth,
//! fits each by counting (`ssr_fit::composite_fit`), computes each number's error three ways
//! (`ssr_fit::counting_fit_errors`), and counts how often the truth lies within one and two
//! reported errors of the estimate. An honest error holds the truth within one error about 68
//! times in 100 and within two about 95 times in 100.
//!
//! The three ways:
//!
//! - **bootstrap**: the class's tracts drawn again with replacement and refitted;
//! - **sandwich**: the composite likelihood's formula, which allows for samples at one tract sharing
//!   its alleles;
//! - **naive**: the inverse curvature alone, as if every sample at a tract were independent.
//!
//! Each setting varies one thing from a base of 200 tracts, 20 samples, 3 reads a sample, a
//! homozygote excess of 0.4 and three allele classes (−1, 0, +1). Every class is drawn at a
//! concentration of 0.5 (a between-locus share of two thirds, near tomato2's) and the stutter
//! `TRUTH`.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example ng_str_counting_fit_errors -- \
//!   --draws 200 --replicates 100 --output tmp/counting_fit_errors [--setting NAME ...]
//! cargo run --release --features bench-fixtures --example ng_str_counting_fit_errors -- \
//!   --cost --draws 3 --replicates 100 --output tmp/counting_fit_errors_cost
//! ```
//!
//! Writes one row a draw and number to `<output>/<setting>.tsv`, and the coverage tables to
//! standard output. `--cost` instead runs the draws one after another and reports how long a fit,
//! the sandwich and the bootstrap take; other programs running on the machine move those times.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;

use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::ssr_fit::composite_fit::{
    ConcentrationSource, CountedClass, fit_class_by_counting, fit_class_by_counting_with_pairs,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::counting_fit_errors::{
    BootstrapSettings, CountedClassErrors, bootstrap_errors, sandwich_errors,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    Slippage, SsrFitConfig, StratumError, StratumEvidence, bench_fixtures,
};

/// The stutter every class is drawn at.
const TRUTH: Slippage = Slippage {
    level: 0.08,
    shorter_share: 0.8,
    fall_off: 0.25,
};

/// The Dirichlet concentration every class is drawn at: a between-locus share of `1 / 1.5`.
const CONCENTRATION: f64 = 0.5;

/// The counting fit's settings, as `examples/ng_str_counting_fit.rs` runs it on real cohorts.
const MAX_ROUNDS: usize = 2000;
const TOLERANCE: f64 = 1e-7;
const ALTERNATIONS: usize = 10;
const START: Slippage = Slippage {
    level: 0.05,
    shorter_share: 0.5,
    fall_off: 0.3,
};

#[derive(Parser)]
#[command(about = "Check the counting fit's standard errors on drawn classes")]
struct Cli {
    /// Classes drawn a setting.
    #[arg(long, default_value_t = 200)]
    draws: usize,
    /// Bootstrap replicates a class.
    #[arg(long, default_value_t = 100)]
    replicates: usize,
    /// Where the per-draw tables go.
    #[arg(long)]
    output: PathBuf,
    /// Run only these settings (by name); every setting when absent.
    #[arg(long)]
    setting: Vec<String>,
    /// Time each method on draws run one after another instead of measuring coverage.
    #[arg(long)]
    cost: bool,
}

/// The allele classes a setting's tracts carry.
#[derive(Clone, Copy)]
enum Window {
    /// −1 ..= +1, reads recorded over the same range.
    One,
    /// −9 ..= +9, reads recorded to ±8: alleles at ±9 send their exact reads into the end buckets.
    Nine,
    /// Alleles drawn over −12 ..= +12 and fitted at −9 ..= +9, reads recorded to ±8: the fit's
    /// outermost lengths stand for "this far or further".
    NineAllelesBeyond,
}

#[derive(Clone, Copy)]
struct Setting {
    name: &'static str,
    tracts: usize,
    samples: usize,
    depth: u32,
    homozygote_excess: f64,
    window: Window,
    source: ConcentrationSource,
}

const fn setting(
    name: &'static str,
    tracts: usize,
    samples: usize,
    depth: u32,
    homozygote_excess: f64,
    window: Window,
    source: ConcentrationSource,
) -> Setting {
    Setting {
        name,
        tracts,
        samples,
        depth,
        homozygote_excess,
        window,
        source,
    }
}

use ConcentrationSource::{EachSamplesTwoCopies as OwnCopies, PairsOfSamples as Pairs};

const SETTINGS: [Setting; 14] = [
    setting("base", 200, 20, 3, 0.4, Window::One, Pairs),
    setting("tracts_20", 20, 20, 3, 0.4, Window::One, Pairs),
    setting("tracts_2000", 2000, 20, 3, 0.4, Window::One, Pairs),
    setting("samples_500", 200, 500, 3, 0.4, Window::One, Pairs),
    setting("reads_30", 200, 20, 30, 0.4, Window::One, Pairs),
    setting("excess_095", 200, 20, 3, 0.95, Window::One, Pairs),
    setting("window_9", 200, 20, 3, 0.4, Window::Nine, Pairs),
    setting(
        "window_9_alleles_beyond",
        200,
        20,
        3,
        0.4,
        Window::NineAllelesBeyond,
        Pairs,
    ),
    setting(
        "window_9_tracts_2000",
        2000,
        20,
        3,
        0.4,
        Window::Nine,
        Pairs,
    ),
    setting("one_sample_reads_3", 2000, 1, 3, 0.4, Window::One, Pairs),
    setting(
        "one_sample_reads_3_own_copies",
        2000,
        1,
        3,
        0.4,
        Window::One,
        OwnCopies,
    ),
    setting("one_sample_reads_30", 2000, 1, 30, 0.4, Window::One, Pairs),
    setting(
        "one_sample_reads_30_own_copies",
        2000,
        1,
        30,
        0.4,
        Window::One,
        OwnCopies,
    ),
    setting(
        "one_sample_reads_30_inbred",
        2000,
        1,
        30,
        1.0,
        Window::One,
        Pairs,
    ),
];

/// **The allele-length variation a setting is drawn at**, over the lengths it is drawn at
/// ([`drawn_span`]). Over three lengths, most chromosomes at the reference length. Over nineteen,
/// the same fall-off from the reference plus extra mass 8 and 9 repeats short, so that the lengths
/// past the recorded reads carry something to estimate. Over twenty-five, the same, the extra mass
/// moved to 10 to 12 repeats short: alleles past the fit's outermost length.
fn spectrum(window: Window) -> Vec<f64> {
    let with_extra = |lengths: usize, extra: &[(usize, f64)]| {
        let mut spectrum = bench_fixtures::spectrum_of(lengths);
        for (index, mass) in extra {
            spectrum[*index] += mass;
        }
        let total: f64 = spectrum.iter().sum();
        spectrum.iter().map(|share| share / total).collect()
    };
    match window {
        Window::One => bench_fixtures::spectrum_of(3),
        Window::Nine => with_extra(19, &[(0, 0.04), (1, 0.03)]),
        Window::NineAllelesBeyond => with_extra(25, &[(0, 0.03), (1, 0.02), (2, 0.02)]),
    }
}

/// **The shares the fit should return**: the drawn variation with every length past the fit's
/// outermost one counted in that outermost one, which means "this far or further".
fn truth_shares(window: Window) -> Vec<f64> {
    let drawn = spectrum(window);
    let beyond = (drawn_span(window) - allele_span(window)) as usize;
    let mut shares = drawn[beyond..drawn.len() - beyond].to_vec();
    let last = shares.len() - 1;
    shares[0] += drawn[..beyond].iter().sum::<f64>();
    shares[last] += drawn[drawn.len() - beyond..].iter().sum::<f64>();
    shares
}

/// The fit's allele span for a window.
fn allele_span(window: Window) -> i32 {
    match window {
        Window::One => 1,
        Window::Nine | Window::NineAllelesBeyond => 9,
    }
}

/// How far either side of the reference the alleles are drawn.
fn drawn_span(window: Window) -> i32 {
    match window {
        Window::One => 1,
        Window::Nine => 9,
        Window::NineAllelesBeyond => 12,
    }
}

/// **One class drawn at `setting`.** At the ±9 windows the reads are drawn over the drawn span and
/// every bucket past ±8 folded into the end bucket, which is the ±8 recording exactly: an end
/// bucket holds every read at least that far.
fn draw(setting: &Setting, seed: u64) -> StratumEvidence {
    let span = drawn_span(setting.window);
    let mut evidence = bench_fixtures::draw_stratum(
        TRUTH,
        &spectrum(setting.window),
        CONCENTRATION,
        setting.homozygote_excess,
        setting.tracts,
        setting.samples,
        setting.depth,
        span,
        seed,
    );
    let recorded = RECORDED_OFFSET_RANGE;
    if span > recorded {
        let beyond = (span - recorded) as usize;
        for tract in &mut evidence.tracts {
            for sample in &mut tract.samples {
                for (_, counts) in &mut sample.by_library {
                    let last = counts.len() - beyond;
                    let mut folded = counts[beyond..last].to_vec();
                    folded[0] += counts[..beyond].iter().sum::<u32>();
                    *folded.last_mut().expect("seventeen buckets") +=
                        counts[last..].iter().sum::<u32>();
                    *counts = folded;
                }
            }
        }
        evidence.read_span = recorded;
    }
    evidence
}

fn config_for(setting: &Setting) -> SsrFitConfig {
    SsrFitConfig {
        allele_span: allele_span(setting.window),
        ..SsrFitConfig::default()
    }
}

fn fit(setting: &Setting, evidence: &StratumEvidence, excess: &[f64]) -> CountedClass {
    let config = config_for(setting);
    match setting.source {
        Pairs => fit_class_by_counting_with_pairs(
            evidence,
            START,
            excess,
            &config,
            MAX_ROUNDS,
            TOLERANCE,
            ALTERNATIONS,
        ),
        OwnCopies => fit_class_by_counting(evidence, START, excess, &config, MAX_ROUNDS, TOLERANCE),
    }
}

fn bootstrap_settings(replicates: usize, seed: u64) -> BootstrapSettings {
    BootstrapSettings {
        replicates,
        seed,
        max_rounds: MAX_ROUNDS,
        tolerance: TOLERANCE,
        alternations: ALTERNATIONS,
    }
}

/// One number of one draw: its truth, its estimate, and its three errors.
struct Number {
    name: String,
    truth: f64,
    estimate: f64,
    bootstrap: StratumError,
    sandwich: StratumError,
    naive: StratumError,
}

/// Every number of a fit, with its errors by each method, in one list.
fn numbers(
    setting: &Setting,
    fitted: &CountedClass,
    bootstrap: &CountedClassErrors,
    sandwich: &CountedClassErrors,
    naive: &CountedClassErrors,
) -> Vec<Number> {
    let span = allele_span(setting.window);
    let truth_shares = truth_shares(setting.window);
    let mut out = Vec::new();
    let mut push = |name: String,
                    truth: f64,
                    estimate: f64,
                    pick: &dyn Fn(&CountedClassErrors) -> StratumError| {
        out.push(Number {
            name,
            truth,
            estimate,
            bootstrap: pick(bootstrap),
            sandwich: pick(sandwich),
            naive: pick(naive),
        });
    };
    push("level".into(), TRUTH.level, fitted.slippage.level, &|e| {
        e.level
    });
    push(
        "shorter_share".into(),
        TRUTH.shorter_share,
        fitted.slippage.shorter_share,
        &|e| e.shorter_share,
    );
    push(
        "fall_off".into(),
        TRUTH.fall_off,
        fitted.slippage.fall_off,
        &|e| e.fall_off,
    );
    for (class, (truth, estimate)) in truth_shares
        .iter()
        .zip(&fitted.allele_length_variation)
        .enumerate()
    {
        push(
            format!("share_{}", class as i32 - span),
            *truth,
            *estimate,
            &|e| e.allele_length_variation[class],
        );
    }
    push(
        "concentration".into(),
        CONCENTRATION,
        fitted.concentration.unwrap_or(f64::NAN),
        &|e| e.concentration,
    );
    push(
        "between_locus_share".into(),
        1.0 / (1.0 + CONCENTRATION),
        fitted.between_locus_share().unwrap_or(f64::NAN),
        &|e| e.between_locus_share,
    );
    out
}

fn error_text(error: StratumError) -> String {
    match error {
        StratumError::Estimated(value) => format!("{value:.6e}"),
        StratumError::NotIdentified => "not_identified".into(),
        StratumError::NotPlaced => "not_placed".into(),
        StratumError::NoShare => "no_share".into(),
    }
}

/// Coverage of one number by one method over the draws.
#[derive(Default)]
struct Coverage {
    with_error: usize,
    without_error: usize,
    within_one: usize,
    within_two: usize,
    error_sum: f64,
}

impl Coverage {
    fn add(&mut self, truth: f64, estimate: f64, error: StratumError) {
        match error.value() {
            Some(error) => {
                self.with_error += 1;
                self.error_sum += error;
                let distance = (estimate - truth).abs();
                self.within_one += usize::from(distance <= error);
                self.within_two += usize::from(distance <= 2.0 * error);
            }
            None => self.without_error += 1,
        }
    }

    fn cell(&self) -> String {
        if self.with_error == 0 {
            return format!("— (none of {})", self.without_error);
        }
        let share = |count: usize| 100.0 * count as f64 / self.with_error as f64;
        let mut cell = format!(
            "{:.0} / {:.0}",
            share(self.within_one),
            share(self.within_two)
        );
        if self.without_error > 0 {
            let _ = write!(cell, " ({} without)", self.without_error);
        }
        cell
    }

    fn mean_error(&self) -> f64 {
        self.error_sum / self.with_error.max(1) as f64
    }
}

fn standard_deviation(values: &[f64]) -> f64 {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    (values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}

fn run_coverage(cli: &Cli, setting: &Setting) {
    let excess = vec![setting.homozygote_excess; setting.samples];
    let config = config_for(setting);
    let started = Instant::now();
    let done = std::sync::atomic::AtomicUsize::new(0);
    let drawn: Vec<(usize, Vec<Number>)> = (0..cli.draws)
        .into_par_iter()
        .map(|draw_index| {
            let seed = 1_000_003 * (draw_index as u64 + 1);
            let evidence = draw(setting, seed);
            let fitted = fit(setting, &evidence, &excess);
            let curvature = sandwich_errors(&evidence, &fitted, &excess, &config);
            let bootstrap = bootstrap_errors(
                &evidence,
                &fitted,
                &excess,
                &config,
                bootstrap_settings(cli.replicates, seed ^ 0x5eed),
            );
            let finished = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if finished.is_multiple_of(10) {
                eprintln!(
                    "{}: {finished} of {} draws, {:.0} s",
                    setting.name,
                    cli.draws,
                    started.elapsed().as_secs_f64()
                );
            }
            (
                fitted.rounds,
                numbers(
                    setting,
                    &fitted,
                    &bootstrap,
                    &curvature.sandwich,
                    &curvature.naive,
                ),
            )
        })
        .collect();
    let seconds = started.elapsed().as_secs_f64();
    let rounds: Vec<usize> = drawn.iter().map(|(rounds, _)| *rounds).collect();
    let per_draw: Vec<Vec<Number>> = drawn.into_iter().map(|(_, numbers)| numbers).collect();

    std::fs::create_dir_all(&cli.output).expect("the output directory can be made");
    let path = cli.output.join(format!("{}.tsv", setting.name));
    let mut file = std::io::BufWriter::new(
        std::fs::File::create(&path).expect("the per-draw table can be written"),
    );
    writeln!(
        file,
        "draw\tnumber\ttruth\testimate\tbootstrap\tsandwich\tnaive"
    )
    .expect("writes");
    for (draw_index, numbers) in per_draw.iter().enumerate() {
        for number in numbers {
            writeln!(
                file,
                "{draw_index}\t{}\t{}\t{}\t{}\t{}\t{}",
                number.name,
                number.truth,
                number.estimate,
                error_text(number.bootstrap),
                error_text(number.sandwich),
                error_text(number.naive)
            )
            .expect("writes");
        }
    }

    println!(
        "\n### {} — {} tracts, {} samples, {} reads a sample, homozygote excess {}, window ±{}, \
         concentration from {}\n",
        setting.name,
        setting.tracts,
        setting.samples,
        setting.depth,
        setting.homozygote_excess,
        allele_span(setting.window),
        match setting.source {
            Pairs => "pairs of samples",
            OwnCopies => "each sample's two copies",
        }
    );
    println!(
        "Rounds of the last counting: median {}, most {}. {} draws, {} bootstrap replicates \
         each, {:.0} s. Cells: % of draws with the truth within 1 / 2 errors.\n",
        median(&rounds),
        rounds.iter().max().copied().unwrap_or(0),
        cli.draws,
        cli.replicates,
        seconds
    );
    println!(
        "| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | \
         mean error: bootstrap / sandwich / naive |"
    );
    println!("|---|---|---|---|---|---|---|---|");
    let names: Vec<String> = per_draw[0].iter().map(|n| n.name.clone()).collect();
    for (index, name) in names.iter().enumerate() {
        let rows: Vec<&Number> = per_draw.iter().map(|numbers| &numbers[index]).collect();
        let mut coverage = [
            Coverage::default(),
            Coverage::default(),
            Coverage::default(),
        ];
        for row in &rows {
            coverage[0].add(row.truth, row.estimate, row.bootstrap);
            coverage[1].add(row.truth, row.estimate, row.sandwich);
            coverage[2].add(row.truth, row.estimate, row.naive);
        }
        let estimates: Vec<f64> = rows.iter().map(|row| row.estimate).collect();
        println!(
            "| {name} | {:.4} | {:.4} | {:.4} | {} | {} | {} | {:.4} / {:.4} / {:.4} |",
            rows[0].truth,
            estimates.iter().sum::<f64>() / estimates.len() as f64,
            standard_deviation(&estimates),
            coverage[0].cell(),
            coverage[1].cell(),
            coverage[2].cell(),
            coverage[0].mean_error(),
            coverage[1].mean_error(),
            coverage[2].mean_error(),
        );
    }
    std::io::stdout().flush().expect("stdout flushes");
}

fn run_cost(cli: &Cli, setting: &Setting) {
    let excess = vec![setting.homozygote_excess; setting.samples];
    let config = config_for(setting);
    let (mut fit_seconds, mut sandwich_seconds, mut bootstrap_seconds) = (0.0, 0.0, 0.0);
    for draw_index in 0..cli.draws {
        let seed = 1_000_003 * (draw_index as u64 + 1);
        let evidence = draw(setting, seed);
        let clock = Instant::now();
        let fitted = fit(setting, &evidence, &excess);
        fit_seconds += clock.elapsed().as_secs_f64();
        let clock = Instant::now();
        let _ = sandwich_errors(&evidence, &fitted, &excess, &config);
        sandwich_seconds += clock.elapsed().as_secs_f64();
        let clock = Instant::now();
        let _ = bootstrap_errors(
            &evidence,
            &fitted,
            &excess,
            &config,
            bootstrap_settings(cli.replicates, seed ^ 0x5eed),
        );
        bootstrap_seconds += clock.elapsed().as_secs_f64();
    }
    let draws = cli.draws as f64;
    println!(
        "| {} | {:.3} | {:.3} ({:.1}×) | {:.3} ({:.1}×) |",
        setting.name,
        fit_seconds / draws,
        sandwich_seconds / draws,
        sandwich_seconds / fit_seconds,
        bootstrap_seconds / draws,
        bootstrap_seconds / fit_seconds,
    );
    std::io::stdout().flush().expect("stdout flushes");
}

fn main() {
    let cli = Cli::parse();
    let chosen: Vec<&Setting> = SETTINGS
        .iter()
        .filter(|setting| cli.setting.is_empty() || cli.setting.iter().any(|n| n == setting.name))
        .collect();
    assert!(!chosen.is_empty(), "no setting is named {:?}", cli.setting);
    if cli.cost {
        println!(
            "| setting | one fit, s | sandwich, s (× a fit) | bootstrap of {} replicates, s (× a fit) |",
            cli.replicates
        );
        println!("|---|---|---|---|");
        for setting in chosen {
            run_cost(&cli, setting);
        }
    } else {
        for setting in chosen {
            run_coverage(&cli, setting);
        }
    }
}

fn median(values: &[usize]) -> usize {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted.get(sorted.len() / 2).copied().unwrap_or(0)
}
